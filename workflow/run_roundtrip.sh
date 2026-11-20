#!/bin/bash
# run_roundtrip.sh
#
# vg construct round trip. For each replicate: write the truth as
# a VCF, build a graph from it with vg construct, index the haplotypes from
# the same VCF, deconstruct, and score the deconstructed VCF against the
# truth with the same scorer used for the PGGB callsets.
#
# Usage, one replicate first:
#   bash run_roundtrip.sh one
# then the twelve comparison windows:
#   bash run_roundtrip.sh all
#
# Environment: REBUILD_DIR (default ~/pangenome_rebuild), REPO (~/graphsv)

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg

MODE=${1:-one}
R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
REPO=${REPO:-$HOME/graphsv}
OUT=$R/roundtrip
mkdir -p "$OUT"
cd "$R" || exit 2

ALL="chr20_10000000 chr20_14000000 chr20_39000000 chr20_2000000 chr20_35000000
     chr20_22000000 chr20_3000000 chr20_34000000 chr20_48000000
     chr20_50000000 chr20_23000000 chr20_25000000"
if [ "$MODE" = "one" ]; then WINDOWS="chr20_10000000"; else WINDOWS="$ALL"; fi
SEED=1

for W in $WINDOWS; do
  P=sweep54d/${W}_s${SEED}
  Q=$OUT/${W}_s${SEED}
  python3 "$REPO/workflow/make_truth_vcf.py" --prefix "$P" \
    --backbone "data/windows/${W}.fa" --out "$Q.truth.vcf" || { echo "VCF FAILED $W"; continue; }
  bgzip -f "$Q.truth.vcf" && tabix -f -p vcf "$Q.truth.vcf.gz"
  CHROM=$(zcat "$Q.truth.vcf.gz" | grep -m1 -v "^#" | cut -f1)

  vg construct -r "data/windows/${W}.fa" -v "$Q.truth.vcf.gz" -a > "$Q.vg" 2> "$Q.construct.log" \
    || { echo "CONSTRUCT FAILED $W"; tail -3 "$Q.construct.log"; continue; }
  vg gbwt -x "$Q.vg" -v "$Q.truth.vcf.gz" -o "$Q.gbwt" 2> "$Q.gbwt.log" \
    || { echo "GBWT FAILED $W"; tail -3 "$Q.gbwt.log"; continue; }
  vg deconstruct -p "$CHROM" -g "$Q.gbwt" "$Q.vg" > "$Q.decon.vcf" 2> "$Q.decon.log" \
    || { echo "DECONSTRUCT FAILED $W"; tail -3 "$Q.decon.log"; continue; }

  # the scorer pairs <prefix>.truth_svs.tsv with <prefix><suffix>
  cp "$P.truth_svs.tsv" "$Q.truth_svs.tsv"
  echo "$W  truth records $(zcat "$Q.truth.vcf.gz" | grep -vc '^#')  deconstructed records $(grep -vc '^#' "$Q.decon.vcf")"
done

echo
python3 "$REPO/analysis/precision_recall_types.py" --sweep "$OUT" --vcf-suffix .decon.vcf --out "$REPO/out/roundtrip"
