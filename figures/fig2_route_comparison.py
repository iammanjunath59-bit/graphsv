#!/usr/bin/env python3
"""
fig2_route_comparison.py

Figure 2. Allele count recovery by route and by graph builder.

Panel a: recovery against structural variant size for the VCF route and the
graph route, in graphs built by PGGB and by Minigraph-Cactus from byte
identical haplotype panels.
Panel b: the same quantities aggregated over all size classes.

The script reads the per variant comparison tables directly, so nothing is
transcribed. Column names are resolved from each file's header rather than by
position, because the VCF comparison tables carry both exact_ac_present
(column 14) and exact_ac_anywhere (column 15) and the two are easy to confuse.

    python3 fig2_route_comparison.py \\
        --pggb-vcf  'sweep54d/*.vcfcmp.svs.tsv' \\
        --pggb-graph 'sweep54d/*.gsv.eval.tsv' \\
        --mc-vcf    'mc/*.mccmp.svs.tsv' \\
        --mc-graph  'mc/*.mcgsv.eval.tsv' \\
        --out figures/Figure2

Expected values, for checking the figure against the reported results:
    PGGB  VCF   0.988 0.977 0.988 0.898 0.782 0.553   overall 0.8646  (n 2,592)
    PGGB  graph 0.972 0.998 1.000 1.000 0.998 0.986   overall 0.9923  (n 2,592)
    MC    VCF   0.991 0.986 0.991 0.912 0.801 0.569   overall 0.8750  (n 1,296)
    MC    graph 0.977 0.995 0.995 0.995 0.995 0.968   overall 0.9877  (n 1,296)

No legend is drawn inside the axes; a single legend is placed below both
panels.
"""

import argparse
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BINS = [100, 300, 1000, 3000, 10000, 30000]
TICK = ["100", "300", "1k", "3k", "10k", "30k"]

# route determines colour, builder determines line style and marker
COLOUR = {"vcf": "#c1666b", "graph": "#1b4965"}
STYLE = {"pggb": ("-", "o"), "mc": ("--", "^")}
NAME = {"pggb": "PGGB", "mc": "Minigraph-Cactus"}


def tally(pattern, success_col):
    """-> ({size_bin: fraction}, {size_bin: n}, overall_fraction, total_n)."""
    hit, tot = {}, {}
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"no files match {pattern}")
    for path in files:
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            idx = {name: i for i, name in enumerate(head)}
            for need in ("size_bin", success_col):
                if need not in idx:
                    sys.exit(f"{path} has no column {need}; found {head}")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) <= max(idx[c] for c in ("size_bin", success_col)):
                    continue
                try:
                    b = int(f[idx["size_bin"]])
                    ok = int(f[idx[success_col]])
                except ValueError:
                    continue
                tot[b] = tot.get(b, 0) + 1
                hit[b] = hit.get(b, 0) + ok
    frac = {b: hit[b] / tot[b] for b in tot}
    return frac, tot, sum(hit.values()) / sum(tot.values()), sum(tot.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pggb-vcf", required=True)
    ap.add_argument("--pggb-graph", required=True)
    ap.add_argument("--mc-vcf", required=True)
    ap.add_argument("--mc-graph", required=True)
    ap.add_argument("--vcf-col", default="exact_ac_present",
                    help="success column in the VCF comparison tables")
    ap.add_argument("--graph-col", default="exact",
                    help="success column in the graph evaluation tables")
    ap.add_argument("--out", default="Figure2")
    args = ap.parse_args()

    series = {}
    for key, pattern, col in (("pggb_vcf", args.pggb_vcf, args.vcf_col),
                              ("pggb_graph", args.pggb_graph, args.graph_col),
                              ("mc_vcf", args.mc_vcf, args.vcf_col),
                              ("mc_graph", args.mc_graph, args.graph_col)):
        series[key] = tally(pattern, col)

    fig = plt.figure(figsize=(9.4, 4.0))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0],
                          left=0.085, right=0.985, top=0.90, bottom=0.30,
                          wspace=0.26)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])

    handles, labels = [], []
    x = list(range(len(BINS)))
    for builder in ("pggb", "mc"):
        for route in ("vcf", "graph"):
            frac = series[f"{builder}_{route}"][0]
            ls, mk = STYLE[builder]
            y = [frac.get(b, float("nan")) for b in BINS]
            ln, = axa.plot(x, y, ls, color=COLOUR[route], marker=mk, ms=5, lw=1.7)
            handles.append(ln)
            labels.append(f"{NAME[builder]}, {'graph route' if route == 'graph' else 'VCF route'}")
    axa.set_xticks(x)
    axa.set_xticklabels(TICK)
    axa.set_xlim(-0.4, len(BINS) - 0.6)
    axa.set_ylim(0, 1.05)
    axa.set_xlabel("structural variant size (bp)")
    axa.set_ylabel("fraction of variants whose true\nallele count is recovered")
    axa.grid(alpha=0.25, lw=0.5)
    axa.set_title("a   Recovery by variant size", loc="left", fontsize=10)

    positions, heights, colours, ticks = [], [], [], []
    for i, builder in enumerate(("pggb", "mc")):
        for j, route in enumerate(("vcf", "graph")):
            positions.append(i + (j - 0.5) * 0.36)
            heights.append(series[f"{builder}_{route}"][2])
            colours.append(COLOUR[route])
        ticks.append(i)
    axb.bar(positions, heights, width=0.34, color=colours, edgecolor="white", lw=0.8)
    axb.set_xticks(ticks)
    axb.set_xticklabels([NAME["pggb"], NAME["mc"]], fontsize=9)
    axb.set_xlim(-0.55, 1.55)
    axb.set_ylim(0, 1.05)
    axb.set_ylabel("overall recovery")
    axb.grid(alpha=0.25, lw=0.5, axis="y")
    axb.set_axisbelow(True)
    axb.set_title("b   Overall, all size classes", loc="left", fontsize=10)

    fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 0.02), columnspacing=2.4)

    out = os.path.splitext(args.out)[0]
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=400)
        print("wrote", f"{out}.{ext}")

    print("\nvalues plotted:")
    for key in ("pggb_vcf", "pggb_graph", "mc_vcf", "mc_graph"):
        frac, tot, overall, n = series[key]
        row = "  ".join(f"{b if b < 1000 else str(b // 1000) + 'k'}:{frac[b]:.3f}"
                        for b in BINS if b in frac)
        print(f"  {key:11} {row}   overall {overall:.4f}  (n = {n})")


if __name__ == "__main__":
    main()
