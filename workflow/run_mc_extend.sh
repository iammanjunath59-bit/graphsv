#!/bin/bash
# run_mc_extend.sh
#
# Extends the Minigraph-Cactus arm from 11 windows to the full chr20 set.
#
# Identical to run_mc_comparison.sh in every step that touches data: same
# simulator invocation and seed, same panel-identity assertion against the
# PGGB sweep, same Cactus invocation, same post-processing (Cactus applies
# vcfbub internally, so only vcfwave is added), same evaluation scripts.
#
# The only change is that the window list is derived rather than hardcoded:
# every window with a seed-1 PGGB replicate that does not yet have a Cactus
# comparison. Completed windows are skipped, so this is resumable and safe
# to rerun.
#
# Rationale for extending: the builder claim currently rests on 11 windows
# against 108 PGGB replicates. That asymmetry is the weakest number in the
# paper and is fixable with background compute alone.
#
# ~12 min per window. Check remaining count before launching.

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg
cd "${PANGENOME_DIR:-$PWD}" || exit 1
# Working directory: set PANGENOME_DIR to the analysis tree holding
# sweep54d/, data/windows/ and the mc/ output directory.

OUT=mc
IMG=quay.io/comparative-genomics-toolkit/cactus:latest
SEED=1
MIN_FREE_GB=${MIN_FREE_GB:-20}
mkdir -p $OUT

# Windows with a seed-1 PGGB replicate but no Cactus comparison yet.
WINDOWS=""
for f in sweep54d/chr20_*_s${SEED}.truth_svs.tsv; do
  [ -e "$f" ] || continue
  B=$(basename "$f" _s${SEED}.truth_svs.tsv)
  [ -f "$OUT/${B}_s${SEED}.mccmp.svs.tsv" ] && continue
  [ -f "data/windows/${B}.fa" ] || { echo "no backbone for $B, skipping"; continue; }
  [ -f "data/windows/${B}_svs.bed" ] || { echo "no SV bed for $B, skipping"; continue; }
  WINDOWS="$WINDOWS $B"
done

N=$(echo $WINDOWS | wc -w)
echo "windows to build: $N   (~$((N * 12)) min ≈ $((N * 12 / 60)) h)"
[ "$N" -eq 0 ] && { echo "nothing to do"; exit 0; }
echo "start $(date '+%F %T')"

DONE=0
for W in $WINDOWS; do
  S=${W#chr20_}
  P=$OUT/${W}_s${SEED}

  FREE=$(df -BG --output=avail . | tail -n 1 | tr -dc '0-9')
  if [ "${FREE:-0}" -lt "$MIN_FREE_GB" ]; then
    echo "STOPPING: only ${FREE}G free, need ${MIN_FREE_GB}G"; break
  fi

  python scripts/simulate_pangenome_truth.py \
    --backbone data/windows/${W}.fa --sv-bed data/windows/${W}_svs.bed \
    --window-start "$S" --n-hap 20 --seed $SEED --ne 20000 \
    --mode stratified --out-prefix "$P" 2>/dev/null

  # the PGGB sweep must have seen the same panel, or the comparison is void
  REF=sweep54d/${W}_s${SEED}.truth_svs.tsv
  if [ -f "$REF" ] && ! diff -q "${P}.truth_svs.tsv" "$REF" >/dev/null; then
    echo "PANEL MISMATCH $W -- skipping"; continue
  fi

  D=$OUT/split_${W}_s${SEED}
  rm -rf "$D"; mkdir -p "$D"
  awk -v d="$D" '/^>/{n=substr($1,2); split(n,a,"#"); f=d"/"a[1]".fa";
                      print ">"a[1] > f; next} {print > f}' \
    "${P}.haplotypes.fa"
  for f in $D/*.fa; do
    b=$(basename "$f" .fa); echo -e "${b}\t/data/${f}"
  done > "${P}.seqfile.txt"

  rm -rf "${P}_js" "${P}_out"
  docker run --rm -v "$PWD":/data $IMG \
    cactus-pangenome "/data/${P}_js" "/data/${P}.seqfile.txt" \
    --outDir "/data/${P}_out" --outName mcw \
    --reference anc --vcf --maxCores 4 >/dev/null 2>&1

  if [ ! -f "${P}_out/mcw.vcf.gz" ]; then
    echo "CACTUS FAILED $W"; rm -rf "$D" "${P}_js" "${P}_out"; continue
  fi

  zcat "${P}_out/mcw.vcf.gz" > "${P}.mc.vcf"
  vcfwave -I 1000 "${P}.mc.vcf" > "${P}.mcwave.vcf" 2>/dev/null
  python scripts/compare_truth_graph.py \
    --prefix "$P" --vcf "${P}.mcwave.vcf" --out "${P}.mccmp" >/dev/null 2>&1

  zcat "${P}_out/mcw.gfa.gz" > "${P}.mc.gfa" 2>/dev/null
  python scripts/graphsv.py --gfa "${P}.mc.gfa" --ref-prefix anc \
    --truth "${P}.truth_svs.tsv" --min-seg 50 \
    --out "${P}.mcgsv.tsv" >/dev/null 2>&1

  rm -rf "$D" "${P}_js" "${P}_out" "${P}.haplotypes.fa"* "${P}.mc.vcf" \
         "${P}.mc.gfa"
  DONE=$((DONE + 1))
  echo "$(date +%H:%M) done $W   [$DONE/$N]"
done

echo "COMPLETE $(date '+%F %T')   built $DONE of $N"
echo
echo "=== Minigraph-Cactus, VCF route (all windows) ==="
cat $OUT/*.mccmp.svs.tsv 2>/dev/null | grep -v "^sv_id" | awk -F'\t' \
  '{t[$4]++; if($14==1)k[$4]++} END{for(b in t) printf "  %6s n=%3d %.3f\n",b,t[b],k[b]/t[b]}' | sort -n
echo "=== Minigraph-Cactus, graph route (all windows) ==="
cat $OUT/*.mcgsv.eval.tsv 2>/dev/null | grep -v "^sv_id" | awk -F'\t' \
  '{t[$4]++; if($8==1)k[$4]++} END{for(b in t) printf "  %6s n=%3d %.3f\n",b,t[b],k[b]/t[b]}' | sort -n
echo
echo "windows with a Cactus comparison: $(ls $OUT/*.mccmp.svs.tsv 2>/dev/null | wc -l)"
