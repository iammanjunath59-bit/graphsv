#!/bin/bash
# run_join_grid.sh
#
# Sensitivity of the VCF join to its tolerances, for Additional file 1,
# Table S3. Runs merge_fix.py on the 108 chr20 release style callsets over
#   size tolerance  --merge-size-frac  0.05 0.10 0.20 0.30   (with 50 bp minimum)
#   start tolerance --merge-pos-frac   0.5 1.0                (with 200 bp minimum)
#   sequence check  --merge-min-kmer   0 (off) 0.5 0.8 0.9
# 32 settings, plus joins restricted to alleles of different records at the
# default setting (size 0.20, start 1.0, k mer 0.8). The scoring tolerances
# are never changed. One line per setting is written to out/join_grid.tsv.
#
# Usage (from ~/graphsv):
#   nohup bash workflow/run_join_grid.sh > join_grid.log 2>&1 &
#   tail -f join_grid.log
#
# Environment: SWEEP (default ~/pangenome_rebuild/sweep54d)
set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg
SWEEP=${SWEEP:-$HOME/pangenome_rebuild/sweep54d}
mkdir -p out/join_grid
T=out/join_grid.tsv
printf "kmer\tsize_frac\tpos_frac\tacross_only\tcarrier_before\tcarrier_after\tgained\tlost\tfalse_fusion\tclusters\tprecision_after\n" > $T

run () {  # kmer size pos across
  local tag="k$1_s$2_p$3_a$4" extra=""
  [ "$4" = 1 ] && extra="--across-records-only"
  python3 analysis/merge_fix.py --sweep "$SWEEP" --vcf-suffix .wave.vcf --out out/join_grid/$tag \
    --merge-min-kmer $1 --merge-size-frac $2 --merge-pos-frac $3 $extra > out/join_grid/$tag.log 2>&1 \
    || { echo "FAILED $tag"; tail -3 out/join_grid/$tag.log; return; }
  S=out/join_grid/$tag.summary.txt
  car=$(grep "exact carrier set, type matched" $S | awk '{print $(NF-1)"\t"$NF}')
  gl=$(grep "variants gained" $S | awk '{print $3"\t"$6}')
  ff=$(awk '/MERGED CLUSTERS/{f=1} f && $1=="false_fusion"{print $2; exit}' $S); ff=${ff:-0}
  tot=$(awk '/MERGED CLUSTERS/{f=1} f && $1=="total"{print $2; exit}' $S)
  pr=$(grep "matching a true variant" $S | awk '{print $NF}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" $1 $2 $3 $4 "$car" "$gl" "$ff" "$tot" "$pr" >> $T
  echo "$(date +%H:%M) $tag  carrier $car  gained/lost $gl  false fusions $ff  precision $pr"
}

for k in 0 0.5 0.8 0.9; do
  for s in 0.05 0.10 0.20 0.30; do
    for p in 0.5 1.0; do
      run $k $s $p 0
    done
  done
done
run 0.8 0.20 1.0 1
echo "COMPLETE"; column -t $T
