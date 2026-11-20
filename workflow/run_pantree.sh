#!/bin/bash
# run_pantree.sh
#
# Runs pantree (Salehi Nowbandegani et al. 2026, Cell Genomics) on the 108
# simulated PGGB graphs and scores its VCFs with the same scorer as the
# vg deconstruct + vcfbub + vcfwave callsets. pantree reports an insertion
# once, with nested point variants as separate records, so it is the
# published graph side answer to the splitting this paper measures.
#
# The reference sample in the simulated graphs is "anc". --no-missingness is
# used as the pantree README advises for graphs whose reference walk may be
# cyclic, which tandem duplications can produce.
#
# Install once, in its own environment (Python 3.10 or newer):
#   git clone https://github.com/oclb/pantree ~/pantree_src
#   python3 -m venv ~/pt && ~/pt/bin/pip install ~/pantree_src
#
# Usage:  bash run_pantree.sh one     (one replicate)
#         bash run_pantree.sh         (all)
#
# Environment: REBUILD_DIR (default ~/pangenome_rebuild), REPO (~/graphsv),
#              PANTREE (default ~/pt/bin/pantree)

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg
PANTREE=${PANTREE:-$HOME/pt/bin/pantree}
# extra pantree options and a tag naming the output, e.g.
#   TAG=pantree_contig PANTREE_ARGS="--dfs-method contiguous --priority-samples anc" bash run_pantree.sh
TAG=${TAG:-pantree}
PANTREE_ARGS=${PANTREE_ARGS:-}
echo "options: --ref-name anc --no-missingness $PANTREE_ARGS   tag: $TAG"
[ -x "$PANTREE" ] || { echo "pantree not found at $PANTREE"; exit 2; }
[ -d "$HOME/pantree_src/.git" ] && echo "pantree commit $(git -C "$HOME/pantree_src" rev-parse --short HEAD)"

R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
REPO=${REPO:-$HOME/graphsv}
OUT=$R/$TAG
mkdir -p "$OUT"

n=0; fail=0
for T in "$R"/sweep54d/*.truth_svs.tsv; do
  B=$(basename "$T" .truth_svs.tsv)
  G=$R/sweep54d/$B.final.gfa
  [ -s "$G" ] || { echo "NO GRAPH $B"; continue; }
  Q=$OUT/$B
  cp "$T" "$Q.truth_svs.tsv"
  if [ ! -s "$Q.$TAG.vcf" ]; then
    /usr/bin/time -f "%e s %M KB" -o "$Q.time" \
      "$PANTREE" gfa2vcf "$G" "$Q.$TAG.vcf" --ref-name anc --no-missingness $PANTREE_ARGS \
      > "$Q.log" 2>&1 || { echo "PANTREE FAILED $B"; tail -3 "$Q.log"; fail=$((fail+1)); rm -f "$Q.$TAG.vcf"; continue; }
  fi
  n=$((n+1))
  [ "${1:-all}" = "one" ] && { cat "$Q.time"; break; }
done
echo "replicates with a pantree VCF: $n   failed: $fail"

cd "$REPO"
python3 analysis/precision_recall_types.py --sweep "$OUT" --vcf-suffix .$TAG.vcf --out out/$TAG
