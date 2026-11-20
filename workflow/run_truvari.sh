#!/bin/bash
# run_truvari.sh
#
# Compares truvari collapse with the record join (merge_fix.py) on the same
# 108 simulated wave callsets, scored by the same scorer.
#
# Each callset's 20 haplotypes are paired into 10 phased diploid samples
# (pair_haps.py), split to biallelic records, sorted and indexed, then
# collapsed under two settings, and split back to haplotypes for scoring:
#   default  truvari's own defaults, except --sizemax raised to 100 kb so the
#            30 kb class is not skipped
#   loose    matching made as permissive as the join rule: --refdist 50000,
#            --pctsize 0.8, --pctseq 0.8, --chain
# The split but uncollapsed callset is scored too, as the baseline.
#
# Usage:  bash run_truvari.sh            (all replicates)
#         bash run_truvari.sh one        (one replicate, to check it runs)
#
# Environment: REBUILD_DIR (default ~/pangenome_rebuild), REPO (~/graphsv)

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg
# truvari is kept in its own virtual environment so its dependencies do not
# change the pg environment used for every other analysis
TRUVARI=${TRUVARI:-$HOME/tv/bin/truvari}
[ -x "$TRUVARI" ] || { echo "truvari not found at $TRUVARI"; exit 2; }
echo "truvari $($TRUVARI version 2>&1 | head -1)"

R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
REPO=${REPO:-$HOME/graphsv}
OUT=$R/truvari
mkdir -p "$OUT"

n=0
for T in "$R"/sweep54d/*.truth_svs.tsv; do
  B=$(basename "$T" .truth_svs.tsv)
  W=$R/sweep54d/$B.wave.vcf
  [ -s "$W" ] || continue
  Q=$OUT/$B
  cp "$T" "$Q.truth_svs.tsv"
  # haplotypes paired into diploid samples (truvari treats a one field
  # genotype as unknown and would never consolidate it), split to biallelic
  if [ ! -s "$Q.norm.vcf.gz.tbi" ]; then
    python3 "$REPO/analysis/pair_haps.py" to-diploid "$W" \
      | bcftools norm -m- 2>/dev/null | bcftools sort -Oz -o "$Q.norm.vcf.gz" 2>/dev/null
    tabix -f -p vcf "$Q.norm.vcf.gz"
  fi
  python3 "$REPO/analysis/pair_haps.py" to-haploid "$Q.norm.vcf.gz" > "$Q.split.vcf"
  for set in tdef tloose; do
    if [ $set = tdef ]; then ARGS="--sizemax 100000"
    else ARGS="--sizemax 100000 --refdist 50000 --pctsize 0.8 --pctseq 0.8 --chain"; fi
    "$TRUVARI" collapse -i "$Q.norm.vcf.gz" -o "$Q.$set.dip.vcf" -c "$Q.$set.removed.vcf" $ARGS \
      > "$Q.$set.log" 2>&1 || { echo "COLLAPSE FAILED $set $B"; tail -3 "$Q.$set.log"; continue; }
    python3 "$REPO/analysis/pair_haps.py" to-haploid "$Q.$set.dip.vcf" > "$Q.$set.vcf"
  done
  n=$((n+1))
  [ "${1:-all}" = "one" ] && break
done
echo "replicates processed: $n"

cd "$REPO"
for s in split tdef tloose; do
  echo; echo "===== $s ====="
  python3 analysis/precision_recall_types.py --sweep "$OUT" --vcf-suffix .$s.vcf --out out/truvari_$s \
    | sed -n '/RECALL under/,/PRECISION/p;/callset alleles of/p;/matching a true/p;/RECALL by true variant type/,/INV/p'
done
