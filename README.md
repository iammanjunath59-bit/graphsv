# graphsv

Structural variant allele counts computed from pangenome graph traversals,
and the benchmark showing why they are not reliably preserved when a graph is
converted to VCF.

This repository accompanies the manuscript *[title]* and contains the tool,
the simulation and evaluation pipeline, the per variant tables behind every
reported number, and a script that re-derives those numbers from a fresh
checkout.

---

## Quick start

```bash
git clone https://github.com/iammanjunath59-bit/graphsv.git
cd graphsv
conda env create -f environment.yml && conda activate pg
bash verification/verify_all.sh 2>&1 | tee verification.txt
```

The verification script recomputes every value reported in the paper from the
tables in `results/` and prints each beside the value the manuscript states.
It needs only Python and the packages in `environment.yml`; the section that
measures runtime is skipped unless a graph is supplied, because graphs are
archived rather than committed (see Data below).

---

## Install

The conda environment covers Python, the analysis packages, and the tools
that are conda installable:

```bash
conda env create -f environment.yml
conda activate pg
```

Four components are not installed by that file and are needed to rebuild
graphs or regenerate callsets. They are not needed to verify the reported
numbers or to redraw the figures.

| Tool | Version used | Notes |
|---|---|---|
| vg | v1.63.1 | `deconstruct`, `paths` |
| vcfbub | 0.1.0 | https://github.com/pangenome/vcfbub |
| vcfwave | 1.0.10 | part of vcflib |
| Minigraph-Cactus | see digest below | run from the official container |

The Cactus container tag is mutable, so the image is pinned by digest:

```
quay.io/comparative-genomics-toolkit/cactus@sha256:fb4f1bab95d9c3b003d1036989f977e8cde2ffcf90a27691b992bd7d2d1a253d
```

Other versions used: PGGB 0.7.4, msprime and tskit 1.4.2, BCFtools 1.24 with
HTSlib 1.23.1, BEDTools v2.31.1, statsmodels 0.14.6, pandas 3.0.5, NumPy
2.4.6, SciPy 1.17.1, Matplotlib 3.10.9.

---

## Using graphsv on your own graph

graphsv reads a GFA in which one path is the reference and reports, for each
site at which haplotype paths diverge from and rejoin the reference, how many
haplotypes carry each distinct traversal.

```bash
python3 src/graphsv.py \
    --gfa your_graph.gfa \
    --ref-prefix anc \
    --min-seg 50 \
    --out calls.tsv
```

- `--ref-prefix` matches the start of the reference path name. Path names are
  expected in PanSN form, so a reference path named `anc#1#chr20:1-1000000#0`
  is selected by `--ref-prefix anc`. List the path names in a graph with
  `vg paths -Lx your_graph.gfa`.
- `--min-seg` is the minimum segment length at which two traversals are
  considered to differ, in base pairs. It sets a floor on the size of variant
  the method can resolve. Recovery was insensitive to this threshold between
  20 and 50 bp in our benchmark and degraded above roughly 100 bp as genuine
  small variants began to be discarded.
- `--min-alt-len` is the minimum length difference between traversals for a
  site to be called, in base pairs. Default 50.
- `--max-merge-gap` merges adjacent sites that share a carrier set and lie
  within this distance, on the reasoning that the two breakpoints of one
  large event are carried by the same haplotypes. Default 50000; set it to 0
  to disable merging. The merge audit in the paper reports how often this
  fuses genuinely distinct variants (2.07 % of calls) and what it costs,
  which is the reported span rather than the allele count.
- Passing `--truth` with a table of known variants additionally writes an
  evaluation file; this is what the benchmark uses and is not needed for
  ordinary use.

On a one megabase window with 20 haplotypes, graphsv takes about 0.1 s and
27 MB, against about 6 s and 141 MB for `vg deconstruct` followed by vcfbub
and vcfwave on the same graph.

---

## Reproducing the benchmark

The full pipeline rebuilds graphs and is the expensive part: roughly 12 min
per window for Minigraph-Cactus and rather less for PGGB, across 54 windows.
The steps below assume `PANGENOME_DIR` is a working directory containing
`data/windows/` and enough space for the graphs.

**1. Simulate a panel with known truth.**

```bash
python3 workflow/simulate_pangenome_truth.py \
    --backbone data/windows/chr20_10000000.fa \
    --sv-bed   data/windows/chr20_10000000_svs.bed \
    --window-start 10000000 --n-hap 20 --seed 1 --ne 20000 \
    --mode stratified --out-prefix sweep54d/chr20_10000000_s1
```

`--mode realistic` draws allele counts from the empirical frequency spectrum
instead of stratifying across size classes.

**2. Build a graph and convert it.** PGGB with `-n 21 -t 4 -p 95 -s 5000`,
then

```bash
REF=$(vg paths -Lx graph.gfa | grep '^anc#')
vg deconstruct -p "$REF" -a graph.gfa > raw.vcf
vcfbub -l 0 -a 100000 -i raw.vcf > bub.vcf
vcfwave -I 1000 bub.vcf > wave.vcf
```

**3. Evaluate both routes.**

```bash
python3 workflow/compare_truth_graph.py --prefix PREFIX --vcf wave.vcf --out PREFIX.vcfcmp
python3 src/graphsv.py --gfa graph.gfa --ref-prefix anc \
    --truth PREFIX.truth_svs.tsv --min-seg 50 --out PREFIX.gsv.tsv
```

**4. The Minigraph-Cactus arm.** `workflow/run_mc_extend.sh` regenerates each
panel from the same seed, asserts byte level identity against the PGGB truth
table before building, runs Cactus in the pinned container, and evaluates the
result. Set `PANGENOME_DIR` first.

**5. Statistics.**

```bash
python3 analysis/analyse_sweep.py --glob 'results/sweep54d/*.vcfcmp.svs.tsv' \
    --windows data/windows_chr20.tsv --out out/an20
python3 analysis/robustness_segdup.py --perSV out/an20.perSV.tsv --out out/robustness20
python3 analysis/audit_merges.py --dir results/sweep54d
python3 analysis/preflight_checks.py --sweep results/sweep54d \
    --windows data/windows_chr20.tsv --nuc data/nuc.txt --label chr20
```

**6. The staged comparison** behind Supplementary Figure S1 regenerates the
raw and post vcfbub callsets from the stored graphs and evaluates all three
stages with the same criterion: `analysis/stage_recovery.py`.

---

## Reproducing the figures

```bash
mkdir -p figures/output
python3 figures/fig1_recovery_by_size.py --chr20 out/an20.cells.tsv \
    --chr16 out/an16.cells.tsv --out figures/output/Figure1
python3 figures/fig2_route_comparison.py \
    --pggb-vcf 'results/sweep54d/*.vcfcmp.svs.tsv' \
    --pggb-graph 'results/sweep54d/*.gsv.eval.tsv' \
    --mc-vcf 'results/mc/*.mccmp.svs.tsv' \
    --mc-graph 'results/mc/*.mcgsv.eval.tsv' \
    --out figures/output/Figure2
python3 figures/fig3_hprc_locus.py --locus results/hprc/locus.tsv \
    --out figures/output/Figure3
python3 figures/plot_supp_fig1.py --figs-dir results/stages \
    --sweep-dir results/sweep54d --out figures/output/supp_fig1
```

Each script prints the values it plotted, and each carries the expected values
in its docstring so a mismatch is visible immediately.

---

## The HPRC analyses

These read the released Human Pangenome Reference Consortium callsets over
HTTPS through their tabix indices. Nothing is downloaded in full.

```bash
B="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38"
export RAW="$B/hprc-v2.1-mc-grch38.raw.vcf.gz"
export WAVE="$B/hprc-v2.1-mc-grch38.wave.vcf.gz"

python3 hprc/sweep_hprc_bubbles.py --regions results/hprc/large_events.regions \
    --out chr20_events.tsv
python3 hprc/wave_dispersal.py --pos 25783205 --event 111144578_112737255
python3 figures/fig3_hprc_locus.py --dump results/hprc/locus.tsv
```

graphsv is not applied to HPRC data at any point. The reference partition in
these analyses comes from vg's own bubble decomposition, recorded in the `AT`
field of the released raw callset, so the comparison is between two
representations produced by the same published pipeline.

---

## Data

`results/` holds the per variant and per event tables that every reported
number is computed from: truth tables, VCF route comparisons, graph route
evaluations, the staged comparison, the Minigraph-Cactus arm, and the HPRC
event classification. These are the evidence base and are small enough to
commit.

`data/` holds the window annotations, the segmental duplication and GC tracks,
and the structural variant intervals used for placement.

Not committed, and archived at [Zenodo DOI] instead:

- the 108 PGGB graphs and 54 Minigraph-Cactus graphs,
- the simulated haplotype FASTA files,
- the generated callsets at all three conversion stages.

All of these regenerate from the deposited seeds and the commands above, since
the simulator is deterministic.

Reference sequence is not included. GRCh38 chromosomes 20 and 16 are available
from the UCSC Genome Browser or Ensembl; segmental duplication and repeat
annotations are the UCSC `genomicSuperDups` and `rmsk` tracks; structural
variant intervals are drawn from gnomAD SV v4.1.

---

## Layout

```
src/            graphsv
workflow/       simulation, graph building, evaluation
analysis/       statistics, robustness, audits
hprc/           real data analyses, all remote
figures/        figure scripts and generated output
verification/   verify_all.sh and its output at submission
data/           window annotations and input tracks
results/        per variant tables, the evidence base
```

---

## Citation

[preprint DOI, then journal reference]

## License

MIT, see LICENSE
