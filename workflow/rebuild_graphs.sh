#!/bin/bash
# rebuild_graphs.sh
#
# Rebuilds the 108 PGGB graphs and their callsets from the regenerated panels,
# following the original sweep step for step:
#
#   pggb -n 21 -t 4 -p 95 -s 5000 -V anc:#
#   vcfbub -l 0 -a 100000
#   vcfwave -I 1000
#   compare_truth_graph.py, graphsv.py --min-seg 50
#
# Unlike the original sweep, it keeps the graph and all three callsets (raw,
# after vcfbub, after vcfwave), because the revision analyses score precision
# and carrier sets at every stage.
#
# It then compares the regenerated per variant tables with the deposited ones.
# Identical tables mean the graph builds themselves reproduce, not only the
# panels; any difference is reported rather than hidden.
#
# Run in the background; it resumes from where it stopped:
#   nohup bash rebuild_graphs.sh > rebuild_graphs.log 2>&1 &
#
# Environment: REBUILD_DIR (default ~/pangenome_rebuild), REPO (~/graphsv)

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg

exec 9>/tmp/rebuild_graphs.lock
flock -n 9 || { echo "already running"; exit 1; }

R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
REPO=${REPO:-$HOME/graphsv}
THREADS=4
cd "$R" || exit 2

for t in pggb vcfbub vcfwave bgzip samtools; do
  command -v "$t" >/dev/null || { echo "FATAL: $t not on PATH"; exit 2; }
done

echo "start $(date '+%F %T')"
n_done=0; n_fail=0

while read -r chrom start end sd rep; do
  [ "$chrom" = "chrom" ] && continue
  W=${chrom}_${start}
  for SEED in 1 2; do
    P=sweep54d/${W}_s${SEED}
    [ -s "${P}.gsv.eval.tsv" ] && { n_done=$((n_done+1)); continue; }
    [ -s "${P}.haplotypes.fa" ] || { echo "NO PANEL ${W}_s${SEED}"; n_fail=$((n_fail+1)); continue; }

    bgzip -@$THREADS -k -f "${P}.haplotypes.fa"
    samtools faidx "${P}.haplotypes.fa.gz"

    rm -rf "${P}_pggb"; mkdir -p "${P}_pggb"
    pggb -i "${P}.haplotypes.fa.gz" -o "${P}_pggb" \
         -n 21 -t $THREADS -p 95 -s 5000 -V anc:# >"${P}.pggb.log" 2>&1

    GFA=$(ls ${P}_pggb/*.smooth.final.gfa 2>/dev/null | head -1)
    VCF=$(ls ${P}_pggb/*.smooth.final.anc.vcf 2>/dev/null | head -1)
    if [ -z "$GFA" ] || [ -z "$VCF" ]; then
      echo "PGGB FAILED ${W}_s${SEED}"; n_fail=$((n_fail+1)); continue
    fi

    cp "$GFA" "${P}.final.gfa"
    cp "$VCF" "${P}.raw.vcf"
    vcfbub -l 0 -a 100000 -i "${P}.raw.vcf" 2>/dev/null > "${P}.bub.vcf"
    vcfwave -I 1000 "${P}.bub.vcf" > "${P}.wave.vcf" 2>/dev/null

    python "$REPO/workflow/compare_truth_graph.py" \
      --prefix "$P" --vcf "${P}.wave.vcf" --out "${P}.vcfcmp" >/dev/null 2>&1
    python "$REPO/src/graphsv.py" --gfa "${P}.final.gfa" --ref-prefix anc \
      --truth "${P}.truth_svs.tsv" --min-seg 50 --out "${P}.gsv.tsv" >/dev/null 2>&1

    rm -rf "${P}_pggb" "${P}.haplotypes.fa.gz" "${P}.haplotypes.fa.gz".*
    n_done=$((n_done+1))
    echo "$(date +%H:%M) done ${W}_s${SEED}  [$n_done]"
  done
done < data/eligible_windows.tsv

echo "builds finished $(date '+%F %T')   done $n_done   failed $n_fail"

# ---------------------------------------------------------------- reproduction
echo
echo "=== do the regenerated per variant tables match the deposited ones? ==="
for kind in vcfcmp.svs.tsv gsv.eval.tsv; do
  same=0; diffr=0; absent=0
  for f in sweep54d/*.${kind}; do
    b=$(basename "$f")
    d="$REPO/results/sweep54d/$b"
    if [ ! -s "$d" ]; then absent=$((absent+1))
    elif cmp -s "$f" "$d"; then same=$((same+1))
    else diffr=$((diffr+1))
    fi
  done
  echo "  $kind   identical $same   different $diffr   no deposited copy $absent"
done
echo
echo "COMPLETE $(date '+%F %T')"
