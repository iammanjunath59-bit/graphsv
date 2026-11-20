# Revision analyses (repository releases v1.0.1 and v1.1.0)

These releases add the analyses made in response to review at BMC Genomics. v1.0.1 left graphsv unchanged (1.0.0). v1.1.0 adds one option to graphsv, --carriers, which writes the carrier haplotypes of each allele; without it the output is identical to 1.0.0 (graphsv 1.1.0). Every number added to the revised manuscript can be regenerated with the commands below, run from the repository root with the `pg` conda environment. Paths to simulated panels assume the rebuilt sweep at `~/pangenome_rebuild` (set `REBUILD_DIR` otherwise); panels regenerate exactly from the deposited seeds.

## Scoring

`analysis/precision_recall_types.py` scores a callset against the truth under six nested criteria (overlap, size, type, count, count on a type matched allele, exact carrier set) and reports precision. It is used for every callset in the revision.

    python3 analysis/precision_recall_types.py --sweep ~/pangenome_rebuild/sweep54d --vcf-suffix .wave.vcf --out out/pr_wave

## Mechanism

| Analysis | Command | Output |
|---|---|---|
| Alleles of split duplications (Table S2) | `python3 analysis/piece_diversity.py --sweep ~/pangenome_rebuild/sweep54d --vcf-suffix .wave.vcf --out out/pieces` | `out/pieces.tsv` |
| Simulation without point mutations | `SEEDS=1 bash workflow/run_mu0.sh full`, then `python3 analysis/precision_recall_types.py --sweep ~/pangenome_rebuild/sweep_mu0 --out out/mu0` and the same with `sweep_mu0_ref` and `--out out/mu0_ref` (the matching main sweep callsets) | `out/mu0.perSV.tsv`, `out/mu0_ref.perSV.tsv` |
| Split duplications without point mutations | `python3 analysis/piece_diversity.py --sweep ~/pangenome_rebuild/sweep_mu0 --vcf-suffix .wave.vcf --out out/pieces_mu0` | `out/pieces_mu0.tsv` |
| Graphs built from a truth VCF with vg construct | `bash workflow/run_roundtrip.sh all` (uses `workflow/make_truth_vcf.py`) | `out/roundtrip.*` |

## Remedies and existing tools

| Analysis | Command | Output |
|---|---|---|
| VCF join, default setting | `python3 analysis/merge_fix.py --sweep ~/pangenome_rebuild/sweep54d --vcf-suffix .wave.vcf --merge-min-kmer 0.8 --out out/mk_0.8_1.0_0.20` | `out/mk_0.8_1.0_0.20.*` |
| VCF join, all settings (Table S3) | `bash workflow/run_join_grid.sh` | `out/join_grid.tsv`, `out/join_grid/` |
| Standalone join for any VCF | `python3 analysis/vcf_partition_merge.py --help` | |
| Truvari 5.4.0 collapse | `bash workflow/run_truvari.sh` (uses `analysis/pair_haps.py`; Truvari in its own environment, `TRUVARI=~/tv/bin/truvari`) | `out/truvari_*` |
| pantree, as written | `bash workflow/run_pantree.sh` | `out/pantree.*` |
| pantree, other readings | `analysis/repolarize_vcf.py <vcf> --ref-sample anc` and `analysis/pantree_dup_carriers.py <vcf> > <out.vcf>` on each pantree VCF (usage in each script's header), then the scorer | `out/pantree_rp.*`, `out/pantree_fair.*` |
| GRef on PGGB graphs (vg v1.78.0) | `VG=~/vg178/vg bash workflow/run_gref.sh` | `out/gref_*` |
| GRef with Minigraph-Cactus --gref | `bash workflow/run_mc_gref.sh` | `out/mcgref_*` |

## HPRC release 2

| Analysis | Command | Output |
|---|---|---|
| Genome wide scan, 22 autosomes, conservative rule | `bash workflow/run_hprc_genome.sh`, then `python3 analysis/hprc_genome_summary.py --dir out/hprc/genome --out out/hprc/genome_summary` | `out/hprc/genome_summary.*` |
| Haplotype subsampling, chromosome 20 | `analysis/hprc_fragmentation.py` with `--size-min 10 --size-frac 0.01 --same-site-only --subsample N --seed S` | `out/hprc/scale_N_sS.events.tsv` |
| Upper bound, chromosome 20 | the same without `--same-site-only` | |
| Dosage of the 62 kb deletion in the PGGB HPRC chr20 graph | `python3 hprc/gfa_locus_dosage.py --gfa <chr20 GFA> --ref-prefix GRCh38 --event 25783207-25845155 --right-start 25850000 --flank 10000 --out out/hprc/graph/del62k` | `out/hprc/graph/del62k.*` |
| VCF and graph carriers compared | `python3 hprc/compare_locus_carriers.py --vcf out/hprc/locus.wave.vcf --dosage out/hprc/graph/del62k.paths.tsv --cut 0.8` | |

## graphsv 1.1.0 (release v1.1.0)

| Analysis | Command | Output |
|---|---|---|
| graphsv with carriers on the 108 rebuilt graphs; check against the deposited calls; score on count, type matched count and exact carrier set | `bash workflow/run_graphsv_revision.sh A` (uses `analysis/score_graphsv.py`) | `out/gsv_car/`, `out/gsv_car.perSV.tsv` |
| graphsv on the 62 kb deletion locus of the PGGB HPRC chr20 graph: extract, run on all haplotypes, diagnose, run without HG02735#1, compare with dosage and VCF | `nohup bash workflow/run_graphsv_revision.sh B > graphsv_hprc.log 2>&1 &` (uses `hprc/gfa_extract_locus.py`, `hprc/diagnose_locus.py`, `hprc/compare_graphsv_locus.py`; `HPRC_GFA` sets the graph) | `out/hprc/graph/del62k_locus.*`, `del62k_locus_x1.gsv.tsv` |
| Figures 3 and 4 | `bash workflow/run_graphsv_revision.sh C` | `figures/output/Figure3.*`, `Figure4.*` |

Notes. graphsv's built in evaluation (`--truth`) checks the count only; `score_graphsv.py` reproduces it exactly on the deposited calls and adds the two stricter criteria used for the VCF routes. On the rebuilt graphs, 20 of 108 replicates give calls slightly different from those deposited in `results/sweep54d` (built from the original graphs); graphsv 1.0.0 and 1.1.0 give identical calls on the same graph, so the difference comes from PGGB, which is not byte reproducible when run with several threads. The extracted locus graph files (`*.gfa`) are not committed (see `.gitignore`); part B regenerates them.

## Figures and supplementary tables

    python3 figures/fig1_type_effect.py --sim out/pr_wave.perSV.tsv --hprc out/hprc/genome_summary.by_chrom.tsv --out figures/output/Figure1
    python3 figures/fig2_mechanism.py --with out/mu0_ref.perSV.tsv --without out/mu0.perSV.tsv --pieces-with out/pieces.tsv --pieces-without out/pieces_mu0.tsv --roundtrip out/roundtrip.perSV.tsv --out figures/output/Figure2
    bash workflow/run_graphsv_revision.sh C      # Figures 3 and 4
    python3 analysis/supp_tables.py --stages results/stages --pieces out/pieces.tsv

Each figure script prints every value it plots. Figures 1 to 4 replace those of v1.0.0; the earlier figure scripts (`fig1_recovery_by_size.py`, `fig2_route_comparison.py`, `fig3_hprc_locus.py`, `plot_supp_fig1.py`) are kept for the record.

## Inputs added

`data/windows/*_svs.bed`: the gnomAD SV v4.1 intervals for each simulated window, as passed to the simulator.
