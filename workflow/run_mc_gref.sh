#!/bin/bash
# run_mc_gref.sh
#
# GRef test through the official route: Minigraph-Cactus with
# --gref, which computes graph reference paths with vg paths -u on graphs whose
# reference path is acyclic by construction, as HPRC release 2 does.
#
# One Cactus run per window builds the graph and writes both VCFs:
#   <name>.vcf.gz       deconstructed against the reference sample (vcfbub applied)
#   <name>.gref.vcf.gz  deconstructed against gref_<reference> (vcfbub applied)
# Both then get vcfwave -I 1000. From the gref VCF only reference contig records
# are scored (CHROM not ending _<N>_alt); GRef contig records are counted, not
# pooled, as the Cactus documentation instructs.
#
# Panels are the main sweep's haplotypes and truth (seed 1), so results are
# comparable with every other arm. Same Cactus image as the
# Minigraph-Cactus comparison (workflow/run_mc_extend.sh).
#
# Usage:  bash run_mc_gref.sh one   |   nohup bash run_mc_gref.sh > mc_gref.log 2>&1 &
# Environment: REBUILD_DIR (~/pangenome_rebuild), REPO (~/graphsv)

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg
R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
REPO=${REPO:-$HOME/graphsv}
IMG=quay.io/comparative-genomics-toolkit/cactus:latest
cd "$R" || exit 2
OUT=mcgref
mkdir -p $OUT
echo "image $(docker inspect --format='{{index .RepoDigests 0}}' $IMG)"

exec 9>/tmp/run_mc_gref.lock
flock -n 9 || { echo "already running"; exit 1; }

WINDOWS="chr20_10000000 chr20_14000000 chr20_39000000 chr20_2000000 chr20_35000000
         chr20_22000000 chr20_3000000 chr20_34000000 chr20_48000000
         chr20_50000000 chr20_23000000 chr20_25000000"
[ "${1:-all}" = "one" ] && WINDOWS="chr20_10000000"
SEED=1

for W in $WINDOWS; do
  P=$OUT/${W}_s${SEED}
  [ -s "$P.mcgref.wave.vcf" ] && { echo "skip $W"; continue; }
  H=sweep54d/${W}_s${SEED}.haplotypes.fa
  [ -s "$H" ] || { echo "NO HAPLOTYPES $H"; continue; }
  cp "sweep54d/${W}_s${SEED}.truth_svs.tsv" "$P.truth_svs.tsv"

  D=$OUT/split_${W}_s${SEED}
  rm -rf "$D"; mkdir -p "$D"
  awk -v d="$D" '/^>/{n=substr($1,2); split(n,a,"#"); f=d"/"a[1]".fa"; print ">"a[1] > f; next} {print > f}' "$H"
  for f in "$D"/*.fa; do b=$(basename "$f" .fa); echo -e "${b}\t/data/${f}"; done > "$P.seqfile.txt"

  docker run --rm -v "$PWD":/data --entrypoint rm $IMG -rf "/data/${P}_js" "/data/${P}_out"
  docker run --rm -v "$PWD":/data $IMG \
    cactus-pangenome "/data/${P}_js" "/data/${P}.seqfile.txt" \
    --outDir "/data/${P}_out" --outName mcg --reference anc --vcf --gref --maxCores 4 \
    > "$P.cactus.log" 2>&1

  if [ ! -s "${P}_out/mcg.vcf.gz" ] || [ ! -s "${P}_out/mcg.gref.vcf.gz" ]; then
    echo "CACTUS OUTPUT MISSING $W; files present:"; ls "${P}_out" 2>/dev/null | head -20
    tail -5 "$P.cactus.log"; rm -rf "$D"
    docker run --rm -v "$PWD":/data --entrypoint rm $IMG -rf "/data/${P}_js"; continue
  fi

  zcat "${P}_out/mcg.vcf.gz" > "$P.mc.bub.vcf"
  zcat "${P}_out/mcg.gref.vcf.gz" > "$P.gref.all.vcf"
  awk -F'\t' '/^#/ || $1 !~ /_[0-9]+_alt$/' "$P.gref.all.vcf" > "$P.mcgref.bub.vcf"
  nalt=$(awk -F'\t' '!/^#/ && $1 ~ /_[0-9]+_alt$/' "$P.gref.all.vcf" | wc -l)
  vcfwave -I 1000 "$P.mc.bub.vcf" > "$P.mc.wave.vcf" 2>/dev/null
  vcfwave -I 1000 "$P.mcgref.bub.vcf" > "$P.mcgref.wave.vcf" 2>/dev/null
  cp "${P}_out/mcg.gref.gref-segs.tsv.gz" "$P.gref_segs.tsv.gz" 2>/dev/null
  # Cactus runs as root inside the container, so its output is root owned;
  # remove it from inside a container too
  rm -rf "$D"
  docker run --rm -v "$PWD":/data --entrypoint rm $IMG -rf "/data/${P}_js" "/data/${P}_out"
  echo "$(date +%H:%M) done $W  mc wave records $(grep -vc '^#' $P.mc.wave.vcf)  gref reference contig wave records $(grep -vc '^#' $P.mcgref.wave.vcf)  gref contig records $nalt"
done

cd "$REPO"
for arm in mc mcgref; do
  echo; echo "===== $arm ====="
  python3 analysis/precision_recall_types.py --sweep "$R/$OUT" --vcf-suffix .$arm.wave.vcf --out out/mcgref_$arm \
    | sed -n '/RECALL under/,/PRECISION/p;/callset alleles of/p;/matching a true/p;/RECALL by size/,/30000/p;/RECALL by true variant type/,/INV/p'
done
