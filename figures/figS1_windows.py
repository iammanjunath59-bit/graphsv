#!/usr/bin/env python3
"""
figS1_windows.py

Additional file 2, Figure S1. Recovery along chromosomes 20 and 16, window
by window, with segmental duplication content.

For each 1 Mb window, the fraction of simulated structural variants whose
true allele count is recovered, for the release style VCF (vg deconstruct,
vcfbub, vcfwave) and for graphsv, and the window's segmental duplication
fraction. Both chromosomes are scored with the same criterion, the one
available in the deposited tables of both: the true count among the counts
at an overlapping site (column 14 of the *.vcfcmp.svs.tsv tables, as in
Table S1 and Table S4; column exact of the *.gsv.eval.tsv tables). Chromosome
20 has two replicates per window (48 variants), chromosome 16 one (24).

Inputs (deposited): results/sweep54d, results/sweep16, data/windows_chr20.tsv
(with header), data/windows_chr16.tsv (no header; chrom start end segdup
repeat Nfrac).

Every plotted value is printed, with the Spearman correlation between
segmental duplication content and recovery and the windows with the lowest
recovery.

Usage:
    python3 figures/figS1_windows.py --out figures/output/FigureS1
"""

import argparse
import csv
import glob
import os
import re
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
VCF, GSV, SD = "#52514e", "#1baf7a", "#eb6834"


def segdup(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "chrom":
                continue
            d[int(f[1])] = float(f[3])
    return d


def per_window(pattern, col):
    k, n = defaultdict(int), defaultdict(int)
    for p in glob.glob(pattern):
        m = re.search(r"_(\d+)_s\d+\.", os.path.basename(p))
        w = int(m.group(1))
        with open(p) as fh:
            for r in csv.reader(fh, delimiter="\t"):
                if r[0] == "sv_id":
                    continue
                n[w] += 1
                k[w] += int(r[col] == "1")
    return {w: (k[w], n[w]) for w in n}


def ranks(x):
    o = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and x[o[j + 1]] == x[o[i]]:
            j += 1
        for t in range(i, j + 1):
            r[o[t]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(a, b):
    ra, rb = ranks(a), ranks(b)
    ma, mb = sum(ra) / len(ra), sum(rb) / len(rb)
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return num / den


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=7.5, length=3, width=0.8)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    chroms = [("chr20", "results/sweep54d", "data/windows_chr20.tsv"),
              ("chr16", "results/sweep16", "data/windows_chr16.tsv")]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    fig, axes = plt.subplots(4, 1, figsize=(7.4, 6.6), constrained_layout=True,
                             gridspec_kw={"height_ratios": [2.2, 0.8, 2.2, 0.8]}, sharex=False)
    report = []
    for c, (chrom, d, wf) in enumerate(chroms):
        sd = segdup(wf)
        v = per_window(f"{d}/*.vcfcmp.svs.tsv", 13)
        g = per_window(f"{d}/*.gsv.eval.tsv", 7)
        ws = sorted(v)
        missing = [w for w in ws if w not in sd or w not in g]
        if missing:
            raise SystemExit(f"{chrom}: windows without segdup or graphsv values: {missing[:5]}")
        x = [w / 1e6 + 0.5 for w in ws]
        yv = [v[w][0] / v[w][1] for w in ws]
        yg = [g[w][0] / g[w][1] for w in ws]
        s = [sd[w] for w in ws]
        ax, axs = axes[2 * c], axes[2 * c + 1]
        ax.plot(x, yv, "o", ms=4, color=VCF, mec="white", mew=0.6, label="Release style VCF")
        ax.plot(x, yg, "v", ms=4, color=GSV, mec="white", mew=0.6, label="graphsv")
        ax.set_ylim(0.4, 1.03)
        ax.yaxis.grid(True, color=GRID, linewidth=0.8)
        ax.set_ylabel("Count recovered", color=INK)
        nvar = v[ws[0]][1]
        ax.set_title(f"{'ab'[c]}  Chromosome {chrom[3:]}: {len(ws)} windows, {nvar} variants per window",
                     loc="left", fontsize=8.5, color=INK)
        if c == 0:
            ax.legend(frameon=False, fontsize=7, loc="lower left", ncol=2)
        style(ax)
        axs.bar(x, s, width=0.8, color=SD)
        axs.set_ylabel("Segmental\nduplication", color=INK, fontsize=7)
        axs.set_ylim(0, max(0.3, max(s) * 1.05))
        axs.set_xlabel(f"Position on chromosome {chrom[3:]} (Mb)", color=INK)
        style(axs)
        for A in (ax, axs):
            A.set_xlim(0, (max(ws) / 1e6) + 1.5)

        rv, rg = spearman(s, yv), spearman(s, yg)
        kv, nv = sum(v[w][0] for w in ws), sum(v[w][1] for w in ws)
        kg, ng = sum(g[w][0] for w in ws), sum(g[w][1] for w in ws)
        report.append(f"{chrom}: {len(ws)} windows; release VCF {kv}/{nv} {kv/nv:.4f}, "
                      f"per window {min(yv):.3f} to {max(yv):.3f}; graphsv {kg}/{ng} {kg/ng:.4f}, "
                      f"per window {min(yg):.3f} to {max(yg):.3f}")
        report.append(f"   Spearman correlation with segmental duplication content: release VCF {rv:+.3f}, graphsv {rg:+.3f}")
        low = sorted(ws, key=lambda w: v[w][0] / v[w][1])[:5]
        report.append("   lowest release VCF windows: " + "; ".join(
            f"{chrom}:{w/1e6:.0f} Mb {v[w][0]}/{v[w][1]} (segdup {sd[w]:.3f}, graphsv {g[w][0]}/{g[w][1]})" for w in low))
        report.append("   every window (start Mb, release VCF k/n, graphsv k/n, segdup): " + "; ".join(
            f"{w/1e6:.0f} {v[w][0]}/{v[w][1]} {g[w][0]}/{g[w][1]} {sd[w]:.3f}" for w in ws))

    for ext in ("pdf", "png"):
        fig.savefig(f"{a.out}.{ext}", dpi=300)
    print(f"wrote {a.out}.pdf and {a.out}.png")
    print("\n".join(report))


if __name__ == "__main__":
    main()
