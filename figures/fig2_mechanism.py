#!/usr/bin/env python3
"""
fig2_mechanism.py

Figure 2. Point variation inside inserted sequence splits duplication counts.

  a  Exact carrier set recovered, by size class, release style VCF, for the
     54 chromosome 20 windows at seed 1, simulated with point mutations (the
     main sweep) and without (mutation rate 0). Duplications and deletions.
     The two sets are like for like in windows, size classes and design but
     not paired variant by variant (the simulator draws structural variants
     from the same random stream as mutations).
     Inputs: perSV tables from precision_recall_types.py for the main sweep
     seed 1 windows (--with) and the zero mutation sweep (--without).

  b  Share of duplications split into two or more callset alleles (no single
     allele with the exact carrier set, two or more pieces whose carriers are
     subsets of the truth), by size class, with and without point mutations.
     Inputs: tables from piece_diversity.py (--pieces-with, --pieces-without).
     Only seed 1 replicates are used from --pieces-with so both sides cover the
     same windows.

  c  Exact carrier set recovered by variant type for three kinds of graph:
     assembled haplotypes with point mutations, assembled haplotypes without
     point mutations, and a graph built with vg construct from a truth VCF in
     which every carrier's inserted copy is identical (12 windows).
     Inputs: the same two perSV tables plus the round trip perSV table.

Every plotted value is printed.

Usage:
    python3 figures/fig2_mechanism.py --with out/mu0_ref.perSV.tsv --without out/mu0.perSV.tsv \\
        --pieces-with out/pieces.tsv --pieces-without out/pieces_mu0.tsv \\
        --roundtrip out/roundtrip.perSV.tsv --out figures/output/Figure2
"""

import argparse
import csv
import math
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BINS = ["100", "300", "1000", "3000", "10000", "30000"]
LAB = ["100 bp", "300 bp", "1 kb", "3 kb", "10 kb", "30 kb"]


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
    ax.tick_params(colors=INK2, labelsize=7.5, length=3, width=0.8)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def read_exact(path):
    """(type, bin) -> [k, n] and type -> [k, n] from a perSV table."""
    by = defaultdict(lambda: [0, 0]); bt = defaultdict(lambda: [0, 0]); reps = set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            e = int(r["carrier_exact"])
            by[(r["type"], str(r["size_bin"]))][0] += e; by[(r["type"], str(r["size_bin"]))][1] += 1
            bt[r["type"]][0] += e; bt[r["type"]][1] += 1
            reps.add(r["replicate"])
    return by, bt, reps


def read_split(path, seed1_only):
    """size bin -> [split, n] for duplications, from a piece_diversity table."""
    out = defaultdict(lambda: [0, 0]); reps = set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["type"] != "DUP":
                continue
            if seed1_only and not r["replicate"].endswith("_s1"):
                continue
            reps.add(r["replicate"])
            split = (r["exact_single"] == "0" and int(r["n_pieces"]) >= 2)
            out[str(r["size_bin"])][0] += split
            out[str(r["size_bin"])][1] += 1
    return out, reps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--with", dest="w", required=True)
    ap.add_argument("--without", dest="wo", required=True)
    ap.add_argument("--pieces-with", required=True)
    ap.add_argument("--pieces-without", required=True)
    ap.add_argument("--roundtrip", required=True)
    ap.add_argument("--out", required=True)
    a_ = ap.parse_args()

    w_by, w_bt, w_reps = read_exact(a_.w)
    o_by, o_bt, o_reps = read_exact(a_.wo)
    r_by, r_bt, r_reps = read_exact(a_.roundtrip)
    sw, sw_reps = read_split(a_.pieces_with, True)
    so, so_reps = read_split(a_.pieces_without, False)
    if {x.rsplit("_s", 1)[0] for x in sw_reps} != {x.rsplit("_s", 1)[0] for x in so_reps}:
        print("WARNING: panel b window sets differ between the two piece tables")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.8), constrained_layout=True,
                             gridspec_kw={"width_ratios": [1.2, 1.0, 1.1]})
    fig.patch.set_facecolor("white")
    xs = list(range(len(BINS)))

    # a
    ax = axes[0]
    for ty, col in (("DEL", BLUE), ("DUP", ORANGE)):
        for by, ls, mfc, tag in ((o_by, "--", "white", "without"), (w_by, "-", col, "with")):
            y = [by[(ty, b)][0] / by[(ty, b)][1] if by[(ty, b)][1] else float("nan") for b in BINS]
            ax.plot(xs, y, ls=ls, color=col, lw=1.8, marker="o", ms=4.5, mfc=mfc, mec=col, mew=1.2,
                    label=f"{'Deletions' if ty == 'DEL' else 'Duplications'}, {tag} mutations")
    ax.set_xticks(xs); ax.set_xticklabels(LAB, rotation=30, ha="right")
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Exact carrier set recovered", color=INK)
    ax.set_xlabel("Size class", color=INK)
    ax.set_title("a  Point mutations on or off", loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=6.3, loc="lower left")
    style(ax)

    # b
    ax = axes[1]
    width = 0.38
    for off, data, fc, ec, lab in ((-width / 2, sw, ORANGE, ORANGE, "with mutations"),
                                    (width / 2, so, "white", ORANGE, "without mutations")):
        y = [data[b][0] / data[b][1] if data[b][1] else 0 for b in BINS]
        ax.bar([x + off for x in xs], y, width=width - 0.04, color=fc, edgecolor=ec, linewidth=1.2,
               label=f"Duplications, {lab}")
        for x, v, b in zip(xs, y, BINS):
            if data[b][1] and data[b][0] == 0:
                ax.text(x + off, 0.012, "0", ha="center", va="bottom", fontsize=6.5, color=INK2)
    ax.set_xticks(xs); ax.set_xticklabels(LAB, rotation=30, ha="right")
    ax.set_ylim(0, 1.03)
    ax.set_ylabel("Split into two or more alleles", color=INK)
    ax.set_xlabel("Size class", color=INK)
    ax.set_title("b  Duplications split", loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=6.3, loc="upper left")
    style(ax)

    # c
    ax = axes[2]
    groups = [("Assembled\nwith\nmutations", w_bt), ("Assembled\nwithout\nmutations", o_bt), ("Built\nfrom a\nVCF", r_bt)]
    gx = list(range(len(groups)))
    for off, ty, col in ((-0.19, "DEL", BLUE), (0.19, "DUP", ORANGE)):
        y = [bt[ty][0] / bt[ty][1] if bt[ty][1] else 0 for _, bt in groups]
        ax.bar([x + off for x in gx], y, width=0.34, color=col,
               label="Deletions" if ty == "DEL" else "Duplications")
    ax.set_xticks(gx); ax.set_xticklabels([g for g, _ in groups], fontsize=6.6)
    ax.set_ylim(0, 1.22)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_ylabel("Exact carrier set recovered", color=INK)
    ax.set_title("c  Graph source", loc="left", fontsize=8.5, color=INK)
    ax.legend(frameon=False, fontsize=6.3, loc="upper center", ncol=2, handlelength=1.2,
              columnspacing=0.8)
    style(ax)

    for ext in ("pdf", "png"):
        fig.savefig(f"{a_.out}.{ext}", dpi=300)
    print(f"wrote {a_.out}.pdf and {a_.out}.png\n")

    print(f"a  exact carrier set, k/n (Wilson 95%); with mutations {len(w_reps)} replicates, "
          f"without {len(o_reps)} replicates")
    for ty in ("DEL", "DUP", "INV"):
        for b, lab in zip(BINS, LAB):
            for name, by in (("with", w_by), ("without", o_by)):
                k, n = by[(ty, b)]
                if n:
                    lo, hi = wilson(k, n)
                    print(f"   {ty} {lab:>7} {name:7}  {k:4d}/{n:<4d} {k/n:.3f}  ({lo:.3f} to {hi:.3f})")
    print(f"\nb  duplications split, k/n; with mutations {len(sw_reps)} seed 1 replicates, "
          f"without {len(so_reps)} replicates")
    for b, lab in zip(BINS, LAB):
        print(f"   {lab:>7}  with {sw[b][0]:3d}/{sw[b][1]:<3d} {sw[b][0]/max(sw[b][1],1):.3f}"
              f"   without {so[b][0]:3d}/{so[b][1]:<3d} {so[b][0]/max(so[b][1],1):.3f}")
    print(f"\nc  exact carrier set by type, k/n; round trip {len(r_reps)} replicates")
    for name, bt in (("assembled with mutations", w_bt), ("assembled without mutations", o_bt),
                     ("built from VCF (vg construct)", r_bt)):
        print("   " + name + ": " + "  ".join(
            f"{ty} {bt[ty][0]}/{bt[ty][1]} {bt[ty][0]/bt[ty][1]:.3f}" for ty in ("DEL", "DUP", "INV") if bt[ty][1]))


if __name__ == "__main__":
    main()
