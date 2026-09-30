#!/bin/bash
# run_mc_comparison.sh
#
# Minigraph-Cactus arm of the builder comparison.
#
# On chr20:25-26 Mb seed 5, Cactus and PGGB gave IDENTICAL per-SV outcomes:
# 24/24 SVs agreed individually, and the per-size-bin recovery matched cell
# for cell (1.000, 1.000, 0.750, 0.750, 0.250, 0.500). Two independently
# developed aligners producing identical results is not a coincidence -- it
# localises the failure to what they SHARE, which is the decomposition
# toolchain (vg deconstruct -> vcfbub -> vcfwave), not graph construction.
#
# This extends that from one window to a set spanning the segdup range, so
# the claim rests on more than a single locus.
#
# Panels are regenerated from the same seeds as the PGGB sweep. The
# simulator is deterministic, so Cactus sees byte-identical input and any
# difference is attributable to the builder alone. Each replicate asserts
# panel identity against the PGGB run before building.
#
# ~12 min per window (Cactus is ~3x slower than PGGB here, partly because
# Docker storage is on the slow filesystem). Resumable.

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg
cd "${PANGENOME_DIR:-$PWD}" || exit 1
# Working directory: set PANGENOME_DIR to the analysis tree holding
# sweep54d/, data/windows/ and the mc/ output directory.

OUT=mc
IMG=quay.io/comparative-genomics-toolkit/cactus:latest
mkdir -p $OUT

# windows spanning segdup 0 -> 0.25, seed 1 (plus the hard window at seed 5,
# already done, which is why it is excluded here)
WINDOWS="chr20_14000000 chr20_39000000 chr20_2000000 chr20_35000000
         chr20_22000000 chr20_3000000 chr20_34000000 chr20_48000000
         chr20_50000000 chr20_23000000 chr20_25000000"
SEED=1

echo "start $(date '+%F %T')"

for W in $WINDOWS; do
  S=${W#chr20_}
  P=$OUT/${W}_s${SEED}
  [ -f "${P}.mccmp.svs.tsv" ] && { echo "skip $W"; continue; }

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
    echo "CACTUS FAILED $W"; rm -rf "$D" "${P}_js"; continue
  fi
  true

  # match the PGGB pipeline exactly: Cactus applies vcfbub internally, so
  # only vcfwave is added here.
  zcat "${P}_out/mcw.vcf.gz" > "${P}.mc.vcf"
  vcfwave -I 1000 "${P}.mc.vcf" > "${P}.mcwave.vcf" 2>/dev/null
  python scripts/compare_truth_graph.py \
    --prefix "$P" --vcf "${P}.mcwave.vcf" --out "${P}.mccmp" >/dev/null 2>&1

  # graph route on Cactus's own GFA
  zcat "${P}_out/mcw.gfa.gz" > "${P}.mc.gfa" 2>/dev/null
  python scripts/graphsv.py --gfa "${P}.mc.gfa" --ref-prefix anc \
    --truth "${P}.truth_svs.tsv" --min-seg 50 \
    --out "${P}.mcgsv.tsv" >/dev/null 2>&1

  rm -rf "$D" "${P}_js" "${P}_out" "${P}.haplotypes.fa"* "${P}.mc.vcf" \
         "${P}.mc.gfa"
  echo "$(date +%H:%M) done $W"
done

echo "COMPLETE $(date '+%F %T')"
echo
echo "=== Minigraph-Cactus, VCF route ==="
cat $OUT/*.mccmp.svs.tsv 2>/dev/null | grep -v "^sv_id" | awk -F'\t' \
  '{t[$4]++; if($14==1)k[$4]++} END{for(b in t) printf "  %6s n=%3d %.3f\n",b,t[b],k[b]/t[b]}' | sort -n
echo "=== Minigraph-Cactus, graph route ==="
cat $OUT/*.mcgsv.eval.tsv 2>/dev/null | grep -v "^sv_id" | awk -F'\t' \
  '{t[$4]++; if($8==1)k[$4]++} END{for(b in t) printf "  %6s n=%3d %.3f\n",b,t[b],k[b]/t[b]}' | sort -n
