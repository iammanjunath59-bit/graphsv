#!/bin/bash
# run_gref.sh
#
# Tests GRef graph reference coordinates (HPRC release 2; vg paths -u, added in
# vg v1.74.0) on the simulated PGGB graphs, to test: does
# deconstructing against the GRef sample keep structural variant counts
# together, or does the outer site stay split?
#
# For each replicate, with the same vg binary for both arms:
#   plain  vg deconstruct -C -a -p <anc path>            (no GRef)
#   gref   vg paths -u -Q anc, then vg deconstruct -C -a -P gref_anc
# These are the commands Cactus uses for --gref (cactus_graphmap_join.py).
# Both then go through vcfbub -l 0 -a 100000 and vcfwave -I 1000, the release
# pipeline. GRef contig records (CHROM ending _<N>_alt) are counted and set
# aside; only reference contig records are scored, because Cactus's own
# documentation says not to pool the two when computing allele counts.
#
# Usage:  bash run_gref.sh one   |   bash run_gref.sh
# Environment: REBUILD_DIR (~/pangenome_rebuild), REPO (~/graphsv),
#              VG (default ~/vg178/vg; vg v1.74.0 or newer is required)

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg
VG=${VG:-$HOME/vg178/vg}
[ -x "$VG" ] || { echo "vg not found at $VG"; exit 2; }
echo "$($VG version | head -1)"

R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
REPO=${REPO:-$HOME/graphsv}
OUT=$R/gref
mkdir -p "$OUT"

n=0
for T in "$R"/sweep54d/*.truth_svs.tsv; do
  B=$(basename "$T" .truth_svs.tsv)
  G=$R/sweep54d/$B.final.gfa
  [ -s "$G" ] || { echo "NO GRAPH $B"; continue; }
  Q=$OUT/$B
  cp "$T" "$Q.truth_svs.tsv"
  if [ ! -s "$Q.vg_gref.wave.vcf" ]; then
    "$VG" convert -g "$G" -p > "$Q.vg" 2> "$Q.convert.log" || { echo "CONVERT FAILED $B"; continue; }
    REF=$("$VG" paths -x "$Q.vg" -L | grep -m1 "^anc#")
    [ -n "$REF" ] || { echo "NO anc PATH $B"; continue; }
    "$VG" paths -x "$Q.vg" -u -Q anc --min-gref-len 50 --gref-segs "$Q.gref_segs.tsv" \
      > "$Q.gref.vg" 2> "$Q.gref.log" || { echo "GREF FAILED $B"; tail -3 "$Q.gref.log"; continue; }
    "$VG" deconstruct "$Q.vg" -C -a -p "$REF" > "$Q.vg_plain.raw.vcf" 2> "$Q.plain.log" \
      || { echo "PLAIN DECONSTRUCT FAILED $B"; tail -3 "$Q.plain.log"; continue; }
    "$VG" deconstruct "$Q.gref.vg" -C -a -P gref_anc > "$Q.vg_gref.all.vcf" 2> "$Q.grefdec.log" \
      || { echo "GREF DECONSTRUCT FAILED $B"; tail -3 "$Q.grefdec.log"; continue; }
    # reference contig records only; GRef contig records counted separately
    awk -F'\t' '/^#/ || $1 !~ /_[0-9]+_alt$/' "$Q.vg_gref.all.vcf" > "$Q.vg_gref.raw.vcf"
    awk -F'\t' '!/^#/ && $1 ~ /_[0-9]+_alt$/' "$Q.vg_gref.all.vcf" | wc -l > "$Q.gref_contig_records"
    for arm in vg_plain vg_gref; do
      vcfbub -l 0 -a 100000 -i "$Q.$arm.raw.vcf" 2>/dev/null > "$Q.$arm.bub.vcf"
      vcfwave -I 1000 "$Q.$arm.bub.vcf" > "$Q.$arm.wave.vcf" 2>/dev/null
    done
    rm -f "$Q.vg" "$Q.gref.vg"
  fi
  n=$((n+1))
  echo "$B  plain wave records $(grep -vc '^#' "$Q.vg_plain.wave.vcf")  gref reference contig wave records $(grep -vc '^#' "$Q.vg_gref.wave.vcf")  gref contig raw records $(cat "$Q.gref_contig_records")"
  [ "${1:-all}" = "one" ] && break
done
echo "replicates: $n"

cd "$REPO"
for arm in vg_plain vg_gref; do
  echo; echo "===== $arm ====="
  python3 analysis/precision_recall_types.py --sweep "$OUT" --vcf-suffix .$arm.wave.vcf --out out/gref_$arm \
    | sed -n '/RECALL under/,/PRECISION/p;/callset alleles of/p;/matching a true/p;/RECALL by size/,/30000/p;/RECALL by true variant type/,/INV/p'
done
