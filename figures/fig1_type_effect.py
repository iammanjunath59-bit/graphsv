#!/usr/bin/env python3
"""
fig1_type_effect.py

Figure 1. The size effect is an insertion effect, in simulation and in the
human pangenome.

  a  Simulation (chr20, 108 replicates, release style VCF: vg deconstruct,
     vcfbub, vcfwave). Fraction of true variants for which a type matched
     callset allele carries exactly the true carrier set, by size class,
     deletions against duplications, with Wilson 95% intervals. Inversions
     (21 variants) are too few to plot and are reported in the printout only.
     Input: the perSV table written by precision_recall_types.py
     (columns type, size_bin, carrier_exact).

  b  HPRC release 2 (Minigraph-Cactus v2.1 GRCh38 wave VCF, 22 autosomes).
     Share of each event's carriers held by its largest single allele (one
     minus the share a per allele reading misses), pooled genome wide, by size class, insertions against deletions, under the
     conservative rule (strict size, one graph site). Grey dots are single
     chromosomes, shown where a chromosome has at least --min-events events in
     that class. Input: the by_chrom table written by hprc_genome_summary.py.

Every plotted value is printed, so the text can quote the figure exactly.

Usage:
    python3 figures/fig1_type_effect.py --sim out/pr_wave.perSV.tsv \\
        --hprc out/hprc/genome_summary.by_chrom.tsv --out figures/output/Figure1
"""

import argparse
import csv
import math
import random
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE = "#2a78d6", "#eb6834"          # validated categorical slots 1 and 2
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
SIM_BINS = ["100", "300", "1000", "3000", "10000", "30000"]
SIM_LAB = ["100 bp", "300 bp", "1 kb", "3 kb", "10 kb", "30 kb"]
HPRC_BINS = ["50-100", "100-300", "300-1000", "1000-3000", "3000-10000", "10000-30000", ">=30000"]
HPRC_LAB = ["50–100 bp", "100–300 bp", "0.3–1 kb", "1–3 kb", "3–10 kb", "10–30 kb", "≥30 kb"]


def wilson(k, n, z=1.959964):
    if n == 0:
        return float("nan"), float("nan")
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=8, length=3, width=0.8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sim", required=True)
    ap.add_argument("--hprc", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-events", type=int, default=20)
    a_ = ap.parse_args()

    # ---------------- simulation ----------------
    k = defaultdict(int); n = defaultdict(int)
    with open(a_.sim) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["type"], str(r["size_bin"]))
            n[key] += 1
            k[key] += int(r["carrier_exact"])
    sim = {}
    for ty in ("DEL", "DUP", "INV"):
        for b in SIM_BINS:
            if n[(ty, b)]:
                lo, hi = wilson(k[(ty, b)], n[(ty, b)])
                sim[(ty, b)] = (k[(ty, b)], n[(ty, b)], k[(ty, b)] / n[(ty, b)], lo, hi)

    # ---------------- HPRC ----------------
    pooled = defaultdict(lambda: [0, 0, 0])    # (kind, bin) -> events, ac, outside
    perchrom = defaultdict(list)                # (kind, bin) -> [share per chromosome]
    chroms = set()
    with open(a_.hprc) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["kind"], r["size_bin"])
            ev, ac, out = int(r["events"]), int(r["ac"]), int(r["outside"])
            p = pooled[key]
            p[0] += ev; p[1] += ac; p[2] += out
            chroms.add(r["chrom"])
            if ev >= a_.min_events and ac:
                perchrom[key].append(1 - out / ac)

    # ---------------- figure ----------------
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.1), constrained_layout=True)
    fig.patch.set_facecolor("white")

    xs = list(range(len(SIM_BINS)))
    for ty, col, name in (("DEL", BLUE, "Deletions"), ("DUP", ORANGE, "Duplications")):
        pts = [sim.get((ty, b)) for b in SIM_BINS]
        x = [i for i, p in zip(xs, pts) if p]
        y = [p[2] for p in pts if p]
        lo = [p[2] - p[3] for p in pts if p]
        hi = [p[4] - p[2] for p in pts if p]
        ntot = sum(p[1] for p in pts if p)
        ax1.errorbar(x, y, yerr=[lo, hi], color=col, lw=2, marker="o", ms=5,
                     mec="white", mew=1, capsize=0, elinewidth=1.2,
                     label=f"{name} (n = {ntot:,})")
        ax1.annotate(name, (x[-1], y[-1]), xytext=(6, 0), textcoords="offset points",
                     va="center", fontsize=8, color=INK)
    ax1.set_xticks(xs); ax1.set_xticklabels(SIM_LAB)
    ax1.set_xlim(-0.4, len(xs) - 0.4 + 1.1)
    ax1.set_ylim(0, 1.03)
    ax1.set_xlabel("Simulated variant size class", color=INK)
    ax1.set_ylabel("Exact carrier set recovered", color=INK)
    ax1.set_title("a  Simulation, release style VCF", loc="left", fontsize=9, color=INK)
    ax1.legend(frameon=False, fontsize=7, loc="lower left")
    style(ax1)

    xs2 = list(range(len(HPRC_BINS)))
    rnd = random.Random(1)
    for kind, col, name in (("DEL", BLUE, "Deletions"), ("INS", ORANGE, "Insertions")):
        for i, b in enumerate(HPRC_BINS):
            for v in perchrom.get((kind, b), []):
                ax2.plot(i + rnd.uniform(-0.12, 0.12) + (-0.08 if kind == "DEL" else 0.08), v,
                         "o", ms=2.4, color="#b9b8b2", zorder=1)
        x = [i for i, b in enumerate(HPRC_BINS) if pooled.get((kind, b), [0, 0, 0])[1]]
        y = [1 - pooled[(kind, HPRC_BINS[i])][2] / pooled[(kind, HPRC_BINS[i])][1] for i in x]
        nev = sum(pooled[(kind, b)][0] for b in HPRC_BINS if (kind, b) in pooled)
        ax2.plot(x, y, color=col, lw=2, marker="o", ms=5, mec="white", mew=1, zorder=3,
                 label=f"{name} ({nev:,} events)")
        ax2.annotate(name, (x[-1], y[-1]), xytext=(6, 0), textcoords="offset points",
                     va="center", fontsize=8, color=INK)
    ax2.set_xticks(xs2); ax2.set_xticklabels(HPRC_LAB, rotation=30, ha="right")
    ax2.set_xlim(-0.4, len(xs2) - 0.4 + 1.3)
    ax2.set_ylim(0, 1.03)
    ax2.set_xlabel("Event size class", color=INK)
    ax2.set_ylabel("Carriers held by the largest allele", color=INK)
    ax2.set_title(f"b  HPRC release 2, {len(chroms)} autosomes", loc="left", fontsize=9, color=INK)
    ax2.legend(frameon=False, fontsize=7, loc="lower left")
    style(ax2)

    for ext in ("pdf", "png"):
        fig.savefig(f"{a_.out}.{ext}", dpi=300)
    print(f"wrote {a_.out}.pdf and {a_.out}.png\n")

    print("a  simulation, exact carrier set, k/n, fraction, Wilson 95% interval")
    for ty in ("DEL", "DUP", "INV"):
        for b, lab in zip(SIM_BINS, SIM_LAB):
            if (ty, b) in sim:
                kk, nn, f, lo, hi = sim[(ty, b)]
                print(f"   {ty} {lab:>7}  {kk:4d}/{nn:<4d} {f:.3f}  ({lo:.3f} to {hi:.3f})")
    print("   (INV not plotted: too few variants)\n")
    print(f"b  HPRC, share of carriers held by the largest allele (plotted) and missed (1 minus it), "
          f"pooled; per chromosome range of the plotted value over chromosomes with >= {a_.min_events} events")
    for kind in ("DEL", "INS"):
        for b, lab in zip(HPRC_BINS, HPRC_LAB):
            if (kind, b) in pooled:
                ev, ac, out = pooled[(kind, b)]
                pc = perchrom.get((kind, b), [])
                rng = f"{min(pc):.3f} to {max(pc):.3f} (n={len(pc)})" if pc else "none plotted"
                print(f"   {kind} {lab:>11}  events {ev:7d}  carriers {ac:9d}  held {1-out/ac:.3f}  missed {out/ac:.3f}  chromosomes {rng}")


if __name__ == "__main__":
    main()
