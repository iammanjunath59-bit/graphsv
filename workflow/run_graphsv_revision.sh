#!/bin/bash
# run_graphsv_revision.sh
#
# Part A. graphsv on the 108 simulated graphs again, with --carriers, so that
#         its calls can be scored on the exact carrier set like the VCF
#         routes. Checks first that the first seven columns are identical to
#         the deposited calls in results/sweep54d (graphsv 1.1.0 without the
#         new column must reproduce 1.0.0), then scores with
#         analysis/score_graphsv.py.
# Part B. graphsv on the 62 kb deletion locus of the PGGB HPRC release 2
#         chr20 graph: extract the locus (hprc/gfa_extract_locus.py), run
#         graphsv --carriers on all haplotypes and again with the one
#         haplotype that removes the anchors after the event set aside
#         (EXCLUDE), and compare with the graph dosage and the VCF
#         (hprc/diagnose_locus.py, hprc/compare_graphsv_locus.py).
# Part C. Figures 3 and 4.
#
# Usage (from ~/graphsv):
#   bash workflow/run_graphsv_revision.sh A      # minutes
#   nohup bash workflow/run_graphsv_revision.sh B > graphsv_hprc.log 2>&1 &
#   bash workflow/run_graphsv_revision.sh C      # Figures 3 and 4
#
# Environment: REBUILD_DIR (~/pangenome_rebuild), HPRC_GFA (PGGB chr20 GFA)
set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg
R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
HPRC_GFA=${HPRC_GFA:-/mnt/d/pangenome_data/20251014_hprc25272.p98-k311.chr20.gfa}
python3 src/graphsv.py --version
TIME=""; [ -x /usr/bin/time ] && TIME="/usr/bin/time -v"   # memory and time, if GNU time is installed

if [ "${1:-}" = A ]; then
  mkdir -p out/gsv_car
  n=0; bad=0
  for G in $R/sweep54d/*.final.gfa; do
    B=$(basename $G .final.gfa)
    python3 src/graphsv.py --gfa $G --ref-prefix anc --carriers \
      --out out/gsv_car/$B.gsv.tsv 2>> out/gsv_car/run.log || { echo "FAILED $B"; continue; }
    n=$((n+1))
    cut -f1-7 out/gsv_car/$B.gsv.tsv | cmp -s - results/sweep54d/$B.gsv.tsv || { echo "DIFFERS from deposited: $B"; bad=$((bad+1)); }
  done
  echo "graphs run: $n   differing from deposited calls: $bad"
  python3 analysis/score_graphsv.py --truth-dir results/sweep54d --calls-dir out/gsv_car --out out/gsv_car
fi

if [ "${1:-}" = B ]; then
  mkdir -p out/hprc/graph
  O=out/hprc/graph/del62k_locus
  echo "$(date +%H:%M) extracting locus"
  $TIME python3 hprc/gfa_extract_locus.py --gfa $HPRC_GFA --ref-prefix GRCh38 \
    --event 25783207-25845155 --right-start 25850000 --flank 10000 --out $O 2> $O.time.txt \
    || { echo "EXTRACTION FAILED"; tail -20 $O.time.txt; exit 1; }
  grep -E "Maximum resident|Elapsed" $O.time.txt
  echo "$(date +%H:%M) running graphsv"
  $TIME python3 src/graphsv.py --gfa $O.gfa --ref-prefix GRCh38 --carriers \
    --out $O.gsv.tsv 2> $O.gsv.err || { echo "GRAPHSV FAILED"; tail -20 $O.gsv.err; exit 1; }
  grep -E "^reference|^haplotypes|^merged|^called|Maximum resident|Elapsed" $O.gsv.err
  python3 hprc/diagnose_locus.py --gfa $O.gfa --offset $O.offset.txt --event 25783207-25845155
  python3 hprc/compare_graphsv_locus.py --calls $O.gsv.tsv --locus-gfa $O.gfa --offset $O.offset.txt \
    --dosage out/hprc/graph/del62k.paths.tsv --vcf out/hprc/locus.wave.vcf --event 25783207-25845155

  # On all 464 haplotypes every anchor lies before the event: haplotype
  # HG02735#1 does not cross the segments just after it (diagnose_locus.py
  # and the check in REVISION.md). Run again with it set aside.
  X=${EXCLUDE:-HG02735#1#CM089340.1}
  echo "$(date +%H:%M) running graphsv without $X"
  grep -v -P "^P\t${X//./\\.}\t" $O.gfa > ${O}_x1.gfa
  echo "paths kept: $(grep -c '^P' ${O}_x1.gfa) (reference included)"
  python3 src/graphsv.py --gfa ${O}_x1.gfa --ref-prefix GRCh38 --carriers --out ${O}_x1.gsv.tsv 2>&1 \
    | grep -E "^reference|^haplotypes|^merged|^called"
  python3 hprc/diagnose_locus.py --gfa ${O}_x1.gfa --offset $O.offset.txt --event 25783207-25845155
  python3 hprc/compare_graphsv_locus.py --calls ${O}_x1.gsv.tsv --locus-gfa ${O}_x1.gfa --offset $O.offset.txt \
    --dosage out/hprc/graph/del62k.paths.tsv --vcf out/hprc/locus.wave.vcf --event 25783207-25845155
  echo "COMPLETE $(date +%H:%M)"
fi

if [ "${1:-}" = C ]; then
  # Figures 3 and 4 from the tables written by A and B
  python3 figures/fig3_routes.py --out-dir out --gsv out/gsv_car.perSV.tsv --out figures/output/Figure3
  python3 figures/fig4_hprc.py --scale-glob 'out/hprc/scale_*_s*.events.tsv' \
    --full out/hprc/frag_chr20_strict_site.events.tsv \
    --dosage out/hprc/graph/del62k.paths.tsv --vcf out/hprc/locus.wave.vcf \
    --gsv-calls out/hprc/graph/del62k_locus_x1.gsv.tsv \
    --gsv-offset out/hprc/graph/del62k_locus.offset.txt --out figures/output/Figure4
fi
