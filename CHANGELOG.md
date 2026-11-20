# Changelog

## v1.0.1 (November 2026)

Analyses added during revision at BMC Genomics. graphsv is unchanged (1.0.0). See REVISION.md for the command behind each new result.

- Scoring with six nested criteria and precision (analysis/precision_recall_types.py).
- Mechanism: alleles of split duplications (piece_diversity.py), simulation without point mutations (run_mu0.sh), graphs built from a truth VCF (run_roundtrip.sh, make_truth_vcf.py).
- VCF join of split alleles (merge_fix.py, vcf_partition_merge.py) and its sensitivity grid (run_join_grid.sh).
- Comparisons with Truvari 5.4.0 (run_truvari.sh, pair_haps.py), pantree (run_pantree.sh, repolarize_vcf.py, pantree_dup_carriers.py) and GRef (run_gref.sh, run_mc_gref.sh).
- HPRC release 2: genome wide scan (hprc_fragmentation.py, hprc_genome_summary.py, run_hprc_genome.sh) and the 62 kb locus checked on the PGGB graph (gfa_locus_dosage.py, compare_locus_carriers.py).
- New Figures 1 to 4 and their scripts; supplementary table values (supp_tables.py).
- Result tables in out/ and the per window structural variant beds in data/windows/.

## v1.0.0

Initial deposit: graphsv 1.0.0, benchmark pipeline, analysis, verification.
