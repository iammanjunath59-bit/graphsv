#!/usr/bin/env python3
"""
fig3_routes.py

Figure 3. Routes to structural variant allele counts on the same 108
simulated callsets (chr20, 2,592 variants).

Criterion: exact carrier set. A route recovers a variant when one of its
alleles overlapping the variant is carried by exactly the true haplotypes.
For the VCF routes this is column carrier_exact of precision_recall_types.py
(type matched allele); for the VCF join, after_car; for graphsv, column
carrier of analysis/score_graphsv.py on calls written with graphsv 1.1.0
--carriers from the same rebuilt graphs as the VCF routes. The type matched
count is printed alongside for every route.

  a  All variants, per route.
  b  Duplications by size class, per route.

Routes and inputs:
  release VCF       vg deconstruct + vcfbub + vcfwave      pr_wave.perSV.tsv
  Truvari default   truvari 5.4.0 collapse                  truvari_tdef.perSV.tsv
  Truvari loose     refdist 50000, pctsize/pctseq 0.8, chain truvari_tloose.perSV.tsv
  VCF join          this work, 21-mer >= 0.8                 mk_0.8_1.0_0.20.perSV.tsv (after_ac)
  pantree           as written (GT)                          pantree.perSV.tsv
  pantree, best     best of three readings per variant       pantree, pantree_rp, pantree_fair
                    (chosen with the truth: an upper bound, not a usable result)
  graphsv           graph route, this work                   out/gsv_car.perSV.tsv

Usage:
    python3 figures/fig3_routes.py --out-dir out --gsv out/gsv_car.perSV.tsv \\
        --out figures/output/Figure3
"""

import argparse
import csv
import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
BINS = ["100", "300", "1000", "3000", "10000", "30000"]
LAB = ["100 bp", "300 bp", "1 kb", "3 kb", "10 kb", "30 kb"]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=7.5, length=3, width=0.8)
    ax.set_axisbelow(True)


def load_std(path, col="ac_in_typed"):
    """(replicate, sv_id) -> (type, bin, count_ok, carrier_ok or None)"""
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            car = int(r["carrier_exact"]) if "carrier_exact" in r else None
            d[(r["replicate"], r["sv_id"])] = (r["type"], str(r["size_bin"]), int(r[col]), car)
    return d


def load_join(path):
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[(r["replicate"], r["sv_id"])] = (r["type"], str(r["size_bin"]), int(r["after_ac"]), int(r["after_car"]))
    return d


def load_gsv(path):
    """score_graphsv.py table: (type, bin, typed count, exact carrier set)"""
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["carrier"] == "":
                raise SystemExit(f"{path}: no carrier column; run graphsv with --carriers")
            d[(r["replicate"], r["sv_id"])] = (r["type"], str(r["size_bin"]), int(r["typed"]), int(r["carrier"]))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="out")
    ap.add_argument("--gsv", default="out/gsv_car.perSV.tsv")
    ap.add_argument("--out", required=True)
    a_ = ap.parse_args()
    O = a_.out_dir

    routes = [
        ("Release VCF", load_std(f"{O}/pr_wave.perSV.tsv")),
        ("Truvari, default", load_std(f"{O}/truvari_tdef.perSV.tsv")),
        ("Truvari, relaxed", load_std(f"{O}/truvari_tloose.perSV.tsv")),
        ("VCF join (this work)", load_join(f"{O}/mk_0.8_1.0_0.20.perSV.tsv")),
        ("pantree, as written", load_std(f"{O}/pantree.perSV.tsv")),
    ]
    p1, p2, p3 = (load_std(f"{O}/{x}.perSV.tsv") for x in ("pantree", "pantree_rp", "pantree_fair"))
    best = {k: (v[0], v[1], int(v[2] or p2[k][2] or p3[k][2]),
                int(v[3] or p2[k][3] or p3[k][3])) for k, v in p1.items()}
    routes.append(("pantree, best reading", best))
    routes.append(("graphsv (this work)", load_gsv(a_.gsv)))

    keys = set(routes[0][1])
    for name, d in routes:
        if set(d) != keys:
            raise SystemExit(f"{name}: variant set differs from the release VCF table "
                             f"({len(d)} vs {len(keys)} variants); refusing to plot unmatched sets")
        mism = sum(1 for k in keys if (d[k][0], d[k][1]) != routes[0][1][k][:2])
        if any(d[k][3] is None for k in keys):
            raise SystemExit(f"{name}: no exact carrier set value for some variants")
        if mism:
            raise SystemExit(f"{name}: {mism} variants disagree on type or size class with the release VCF table")

    # ---------------- figure ----------------
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.2), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1.0, 1.15]})
    fig.patch.set_facecolor("white")

    COL = {"Release VCF": "#52514e", "Truvari, default": "#f4b399", "Truvari, relaxed": "#eb6834",
           "VCF join (this work)": "#2a78d6", "pantree, as written": "#b9b8b2",
           "pantree, best reading": "#b9b8b2", "graphsv (this work)": "#1baf7a"}
    names = [n for n, _ in routes]
    vals = [sum(d[k][3] for k in keys) / len(keys) for _, d in routes]
    y = list(range(len(routes)))[::-1]
    for yi, n, v in zip(y, names, vals):
        if n == "pantree, best reading":
            ax1.barh(yi, v, height=0.62, color="white", edgecolor=COL[n], hatch="////", linewidth=1.0)
        else:
            ax1.barh(yi, v, height=0.62, color=COL[n])
        ax1.text(v + 0.01, yi, f"{v:.3f}", va="center", fontsize=7, color=INK)
    ax1.set_yticks(y); ax1.set_yticklabels(names, fontsize=7.3)
    ax1.set_xlim(0, 1.12); ax1.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax1.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax1.set_xlabel("Exact carrier set recovered, all variants", color=INK)
    ax1.set_title(f"a  All variants (n = {len(keys):,})", loc="left", fontsize=8.5, color=INK)
    style(ax1)

    # four routes only: five hues fail the all pairs palette check; pantree is
    # in panel a and its per size values are printed below
    plot = [("Release VCF", "#52514e", "o", "--"), ("Truvari, relaxed", "#eb6834", "s", "-"),
            ("VCF join (this work)", "#2a78d6", "D", "-"), ("graphsv (this work)", "#1baf7a", "v", "-")]
    xs = list(range(len(BINS)))
    lookup = dict(routes)
    for n, col, m, ls in plot:
        d = lookup[n]
        yv = []
        for b in BINS:
            sub = [d[k][3] for k in keys if d[k][0] == "DUP" and d[k][1] == b]
            yv.append(sum(sub) / len(sub))
        ax2.plot(xs, yv, color=col, lw=1.8, ls=ls, marker=m, ms=4.5, mec="white", mew=0.8, label=n)
    ax2.set_xticks(xs); ax2.set_xticklabels(LAB)
    ax2.set_ylim(0, 1.03)
    ax2.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax2.set_xlabel("Duplication size class", color=INK)
    ax2.set_ylabel("Exact carrier set recovered", color=INK)
    ndup = sum(1 for k in keys if routes[0][1][k][0] == "DUP")
    ax2.set_title(f"b  Duplications (n = {ndup:,})", loc="left", fontsize=8.5, color=INK)
    ax2.legend(frameon=False, fontsize=6.6, loc="lower left")
    style(ax2)

    for ext in ("pdf", "png"):
        fig.savefig(f"{a_.out}.{ext}", dpi=300)
    print(f"wrote {a_.out}.pdf and {a_.out}.png\n")

    print("a  all variants: exact carrier set (plotted), type matched count")
    for n, d in routes:
        c = sum(d[k][3] for k in keys)
        t = sum(d[k][2] for k in keys)
        print(f"   {n:24} carrier set {c}/{len(keys)} {c/len(keys):.4f}   typed count {t}/{len(keys)} {t/len(keys):.4f}")
    print("\nb  duplications by size class, exact carrier set, k/n")
    for n in [p[0] for p in plot] + ["pantree, as written", "pantree, best reading"]:
        d = lookup[n]
        cells = []
        for b, lab in zip(BINS, LAB):
            sub = [d[k][3] for k in keys if d[k][0] == "DUP" and d[k][1] == b]
            cells.append(f"{lab} {sum(sub)}/{len(sub)} {sum(sub)/len(sub):.3f}")
        print(f"   {n:24} " + "  ".join(cells))
    print("\n   by type, exact carrier set:")
    for n, d in routes:
        print(f"   {n:24} " + "  ".join(
            f"{ty} {sum(d[k][3] for k in keys if d[k][0]==ty)}/{sum(1 for k in keys if d[k][0]==ty)}"
            for ty in ("DEL", "DUP", "INV")))


if __name__ == "__main__":
    main()
