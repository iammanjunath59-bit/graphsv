#!/usr/bin/env python3
"""
fig1_recovery_by_size.py

Figure 1. Allele count recovery of the VCF route as a function of structural
variant size, stratified by segmental duplication content, on the discovery
chromosome and on the held out chromosome.

The script plots the cells tables written by analyse_sweep.py rather than
recomputing the stratification, so the figure and the reported statistics come
from the same code path. Generate the inputs first:

    python3 scripts/analyse_sweep.py --glob 'sweep54d/*.vcfcmp.svs.tsv' \\
        --windows data/windows_chr20.tsv --out out/an20
    python3 scripts/analyse_sweep.py --glob 'sweep16/*.vcfcmp.svs.tsv' \\
        --windows data/windows_chr16.tsv --out out/an16

then

    python3 fig1_recovery_by_size.py --chr20 out/an20.cells.tsv \\
        --chr16 out/an16.cells.tsv --out figures/Figure1

Expected values, for checking the figure against the reported results:
    chr20 low  0.987 0.982 0.990 0.906 0.786 0.571   (n = 392 per class)
    chr20 high 1.000 0.925 0.975 0.825 0.750 0.375   (n =  40 per class)
    chr16 low  0.994 0.989 0.994 0.972 0.875 0.562   (n = 176 per class)
    chr16 high 0.956 0.926 0.897 0.853 0.691 0.485   (n =  68 per class)

No legend is drawn inside the axes; a single legend is placed below both
panels.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOW_COLOR = "#1b4965"
HIGH_COLOR = "#c1666b"
BINS = [100, 300, 1000, 3000, 10000, 30000]
TICK = ["100", "300", "1k", "3k", "10k", "30k"]


def read_cells(path):
    """Parse the whitespace aligned cells table written by analyse_sweep.py.

    Expected columns: bin stratum n recov lo95 hi95
    Returns {stratum: {bin: (recov, lo, hi, n)}}
    """
    out = {}
    with open(path) as fh:
        header = None
        for line in fh:
            f = line.split()
            if not f:
                continue
            if f[0] == "bin":
                header = f
                continue
            if header is None or len(f) < 6:
                continue
            try:
                b = int(f[0])
                n = int(f[2])
                recov, lo, hi = float(f[3]), float(f[4]), float(f[5])
            except ValueError:
                continue
            out.setdefault(f[1], {})[b] = (recov, lo, hi, n)
    if not out:
        raise SystemExit(f"no cells parsed from {path}")
    return out


def draw(ax, cells, title):
    x = list(range(len(BINS)))
    handles = []
    for stratum, colour, marker, offset in (("low", LOW_COLOR, "o", -0.06),
                                            ("high", HIGH_COLOR, "s", 0.06)):
        if stratum not in cells:
            continue
        xs, ys, lo, hi = [], [], [], []
        for i, b in enumerate(BINS):
            if b not in cells[stratum]:
                continue
            r, l, h, _ = cells[stratum][b]
            xs.append(i + offset)
            ys.append(r)
            lo.append(r - l)
            hi.append(h - r)
        line = ax.errorbar(xs, ys, yerr=[lo, hi], color=colour, marker=marker,
                           ms=5, lw=1.7, capsize=3, elinewidth=1.0,
                           label=("low segmental duplication"
                                  if stratum == "low"
                                  else "high segmental duplication"))
        handles.append(line)
    ax.set_xticks(x)
    ax.set_xticklabels(TICK)
    ax.set_xlim(-0.5, len(BINS) - 0.5)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("structural variant size (bp)")
    ax.grid(alpha=0.25, lw=0.5)
    ax.set_title(title, loc="left", fontsize=10)
    return handles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chr20", required=True, help="an20.cells.tsv")
    ap.add_argument("--chr16", required=True, help="an16.cells.tsv")
    ap.add_argument("--out", default="Figure1")
    args = ap.parse_args()

    c20 = read_cells(args.chr20)
    c16 = read_cells(args.chr16)

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.9), sharey=True)
    fig.subplots_adjust(left=0.095, right=0.985, top=0.90, bottom=0.30, wspace=0.10)

    handles = draw(axes[0], c20, "a   Chromosome 20 (discovery, 108 replicates)")
    draw(axes[1], c16, "b   Chromosome 16 (held out, 61 replicates)")
    axes[0].set_ylabel("fraction of variants whose true\nallele count is recovered")

    labels = [h.get_label() for h in handles]
    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 0.045))

    out = os.path.splitext(args.out)[0]
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=400)
        print("wrote", f"{out}.{ext}")

    print("\nvalues plotted:")
    for name, cells in (("chr20", c20), ("chr16", c16)):
        for stratum in ("low", "high"):
            if stratum not in cells:
                continue
            row = "  ".join(f"{b if b < 1000 else str(b // 1000) + 'k'}:"
                            f"{cells[stratum][b][0]:.3f}"
                            for b in BINS if b in cells[stratum])
            n = cells[stratum][BINS[0]][3]
            print(f"  {name} {stratum:4} (n={n:3d} per class)  {row}")


if __name__ == "__main__":
    main()
