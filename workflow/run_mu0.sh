#!/bin/bash
# run_mu0.sh
#
# Causal test of the splitting mechanism. Rebuilds a set of replicates with
# the point mutation rate set to zero and everything else as in the main
# sweep (same backbone, gnomAD beds, seed, Ne, panel size, PGGB settings,
# vcfbub and vcfwave). If duplications are split because their carriers'
# inserted copies differ by point mutations, then without point mutations
# duplication recall should approach that of deletions at every size.
#
# For each replicate the script also reports whether the zero mutation truth
# table is byte identical to the main sweep's. Identical means the same
# variants, sizes and carriers, so the comparison is perfectly paired and the
# only difference is the absence of point mutations. If it is not identical,
# the comparison is still like for like in window, size and type, but not
# paired variant by variant, and the output says so.
#
# Usage, one replicate first to check that --mu 0 is accepted and yields
# structural variants:
#   bash run_mu0.sh one
# then the twelve comparison windows:
#   nohup bash run_mu0.sh all > run_mu0.log 2>&1 &
# or every eligible window (SEEDS="1 2" for both seeds):
#   nohup bash run_mu0.sh full > run_mu0_full.log 2>&1 &
#
# Environment: REBUILD_DIR (default ~/pangenome_rebuild), REPO (~/graphsv)

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg

MODE=${1:-one}
R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
REPO=${REPO:-$HOME/graphsv}
SIM=$REPO/workflow/simulate_pangenome_truth.py
THREADS=4
OUT=sweep_mu0
REF=sweep_mu0_ref      # links to the matching main sweep files, for pairing
cd "$R" || exit 2
mkdir -p $OUT $REF

exec 9>/tmp/run_mu0.lock
flock -n 9 || { echo "already running"; exit 1; }

# the eleven Minigraph-Cactus comparison windows (segdup 0 to 0.25) plus
# chr20_10000000; seed 1 throughout
ALL="chr20_10000000 chr20_14000000 chr20_39000000 chr20_2000000 chr20_35000000
     chr20_22000000 chr20_3000000 chr20_34000000 chr20_48000000
     chr20_50000000 chr20_23000000 chr20_25000000"
if [ "$MODE" = "one" ]; then WINDOWS="chr20_10000000"
elif [ "$MODE" = "full" ]; then
  # every eligible window; already finished replicates are skipped
  WINDOWS=$(awk 'NR>1{print $1"_"$2}' data/eligible_windows.tsv)
else WINDOWS="$ALL"; fi
SEEDS=${SEEDS:-1}

echo "start $(date '+%F %T')  mode $MODE"
for SEED in $SEEDS; do
for W in $WINDOWS; do
  S=${W#chr20_}
  P=$OUT/${W}_s${SEED}
  [ -s "${P}.wave.vcf" ] && { echo "skip $W"; continue; }

  python "$SIM" --backbone "data/windows/${W}.fa" --sv-bed "data/windows/${W}_svs.bed" \
    --window-start "$S" --n-hap 20 --seed $SEED --ne 20000 --mu 0 \
    --mode stratified --out-prefix "$P" > "${P}.sim.log" 2>&1
  if [ ! -s "${P}.truth_svs.tsv" ]; then
    echo "SIMULATION FAILED $W -- see ${P}.sim.log"; tail -5 "${P}.sim.log"; continue
  fi
  nsv=$(($(wc -l < "${P}.truth_svs.tsv") - 1))
  if cmp -s "${P}.truth_svs.tsv" "sweep54d/${W}_s${SEED}.truth_svs.tsv"; then
    pair="IDENTICAL truth table (paired)"
  else
    pair="truth table DIFFERS from main sweep (not paired)"
  fi
  echo "$W  structural variants $nsv  $pair"

  ln -sf "$R/sweep54d/${W}_s${SEED}.truth_svs.tsv" "$REF/"
  ln -sf "$R/sweep54d/${W}_s${SEED}.wave.vcf" "$REF/"

  bgzip -@$THREADS -k -f "${P}.haplotypes.fa"
  samtools faidx "${P}.haplotypes.fa.gz"
  rm -rf "${P}_pggb"; mkdir -p "${P}_pggb"
  pggb -i "${P}.haplotypes.fa.gz" -o "${P}_pggb" \
       -n 21 -t $THREADS -p 95 -s 5000 -V anc:# > "${P}.pggb.log" 2>&1
  GFA=$(ls ${P}_pggb/*.smooth.final.gfa 2>/dev/null | head -1)
  VCF=$(ls ${P}_pggb/*.smooth.final.anc.vcf 2>/dev/null | head -1)
  if [ -z "$GFA" ] || [ -z "$VCF" ]; then echo "PGGB FAILED $W"; continue; fi
  cp "$GFA" "${P}.final.gfa"
  cp "$VCF" "${P}.raw.vcf"
  vcfbub -l 0 -a 100000 -i "${P}.raw.vcf" 2>/dev/null > "${P}.bub.vcf"
  vcfwave -I 1000 "${P}.bub.vcf" > "${P}.wave.vcf" 2>/dev/null
  rm -rf "${P}_pggb" "${P}.haplotypes.fa.gz" "${P}.haplotypes.fa.gz".*
  echo "$(date +%H:%M) done ${W}_s${SEED}   wave records $(grep -vc '^#' ${P}.wave.vcf)"
done
done
echo "COMPLETE $(date '+%F %T')"
echo
echo "Score both sets with:"
echo "  cd $REPO"
echo "  python3 analysis/precision_recall_types.py --sweep $R/$OUT --out out/mu0"
echo "  python3 analysis/precision_recall_types.py --sweep $R/$REF --out out/mu0_ref"
