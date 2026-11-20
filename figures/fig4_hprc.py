#!/usr/bin/env python3
"""
fig4_hprc.py

Figure 4. HPRC release 2: splitting grows with panel size, and the record join
recovers a split deletion that an independently built graph confirms.

  a  chr20 insertions, conservative rule (strict size, one graph site): share
     of each event's carriers held by its largest allele, by size class, for
     50, 100 and 200 randomly drawn haplotypes (mean of the seeds, bars span
     the seed range) and for all haplotypes. Genotypes are subsampled within
     a graph built from all haplotypes; graphs were not rebuilt.
     Inputs: events tables from hprc_fragmentation.py
     (--scale-glob for the subsamples, --full for all haplotypes).

  b  The 61,948 bp deletion at chr20:25,783,207. Dosage of the deleted
     interval relative to GRCh38 for every haplotype path that passes both
     flanks in the PGGB HPRC chr20 graph (no VCF used), coloured by the
     haplotype's status in the joined VCF records (no graph used).
     Inputs: gfa_locus_dosage.py paths table; the locus wave VCF.

  c  The same event counted three ways: the largest single VCF record, the
     joined VCF records, and the graph (dosage below --cut).

Every plotted value is printed.

Usage:
    python3 figures/fig4_hprc.py --scale-glob 'out/hprc/scale_*_s*.events.tsv' \\
        --full out/hprc/frag_chr20_strict_site.events.tsv \\
        --dosage out/hprc/graph/del62k.paths.tsv --vcf out/hprc/locus.wave.vcf \\
        --out figures/output/Figure4
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
BLUE, GREY, LIGHT = "#2a78d6", "#8a8984", "#d9d8d3"
RAMP = ["#f6c3a8", "#f19a6d", "#eb6834", "#a8441a"]       # 50, 100, 200, all
BINS = ["50-100", "100-300", "300-1000", "1000-3000", "3000-10000", "10000-30000"]
LAB = ["50–100 bp", "100–300 bp", "0.3–1 kb", "1–3 kb", "3–10 kb", "10–30 kb"]


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=7.5, length=3, width=0.8)
    ax.set_axisbelow(True)


def held(path, kind="INS"):
    """size bin -> share of carriers held by the largest allele"""
    ac = defaultdict(int); big = defaultdict(int)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] != kind:
                continue
            ac[r["size_bin"]] += int(r["ac"]); big[r["size_bin"]] += int(r["largest_ac"])
    return {b: big[b] / ac[b] for b in ac if ac[b]}


def vcf_carriers(path, lo, hi, min_del):
    samples, car, called, recs = [], set(), set(), []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]; continue
            pos = int(f[1])
            if not lo <= pos <= hi:
                continue
            want = {j for j, alt in enumerate(f[4].split(","), 1) if len(f[3]) - len(alt) >= min_del}
            if not want:
                continue
            here = set()
            for si, g in enumerate(f[9:]):
                for h, x in enumerate(g.split(":")[0].replace("/", "|").split("|")):
                    if x == ".":
                        continue
                    called.add((samples[si], str(h + 1)))
                    if x.isdigit() and int(x) in want:
                        here.add((samples[si], str(h + 1)))
            recs.append((pos, len(here)))
            car |= here
    return car, called, recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale-glob", required=True)
    ap.add_argument("--full", required=True)
    ap.add_argument("--dosage", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--cut", type=float, default=0.8)
    ap.add_argument("--min-del", type=int, default=50000)
    ap.add_argument("--from", dest="frm", type=int, default=25700000)
    ap.add_argument("--to", type=int, default=25900000)
    ap.add_argument("--out", required=True)
    a_ = ap.parse_args()

    # ---------------- a ----------------
    sub = defaultdict(list)
    for p in sorted(glob.glob(a_.scale_glob)):
        m = re.search(r"scale_(\d+)_s(\d+)\.events\.tsv$", os.path.basename(p))
        if m:
            sub[int(m.group(1))].append(held(p))
    if not sub:
        raise SystemExit("no subsample events tables matched --scale-glob")
    full = held(a_.full)

    # ---------------- b, c ----------------
    gpaths = defaultdict(list)
    with open(a_.dosage) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["complete"] == "1":
                gpaths[(r["sample"], r["hap"])].append(float(r["dosage"]))
    if any(len(v) > 1 for v in gpaths.values()):
        raise SystemExit("some haplotypes have several complete paths; resolve before plotting")
    dos = {k: v[0] for k, v in gpaths.items()}
    vcar, called, recs = vcf_carriers(a_.vcf, a_.frm, a_.to, a_.min_del)
    gcar = {k for k, v in dos.items() if v < a_.cut}

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    fig = plt.figure(figsize=(7.4, 5.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.95], width_ratios=[1.0, 1.0])
    ax1 = fig.add_subplot(gs[0, :]); ax2 = fig.add_subplot(gs[1, 0]); ax3 = fig.add_subplot(gs[1, 1])
    fig.patch.set_facecolor("white")

    xs = list(range(len(BINS)))
    for col, n in zip(RAMP, sorted(sub)):
        mean = [sum(d[b] for d in sub[n]) / len(sub[n]) for b in BINS]
        lo = [mean[i] - min(d[b] for d in sub[n]) for i, b in enumerate(BINS)]
        hi = [max(d[b] for d in sub[n]) - mean[i] for i, b in enumerate(BINS)]
        ax1.errorbar(xs, mean, yerr=[lo, hi], color=col, lw=1.8, marker="o", ms=4.5, mec="white",
                     mew=0.8, capsize=0, elinewidth=1.0, label=f"{n} haplotypes ({len(sub[n])} draws)")
    ax1.plot(xs, [full[b] for b in BINS], color=RAMP[-1], lw=2.2, marker="o", ms=5, mec="white",
             mew=0.8, label="all haplotypes")
    ax1.set_xticks(xs); ax1.set_xticklabels(LAB)
    ax1.set_ylim(0, 1.03); ax1.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax1.set_xlabel("Insertion size class", color=INK)
    ax1.set_ylabel("Carriers held by the largest allele", color=INK)
    ax1.set_title("a  chr20 insertions: more haplotypes, more splitting", loc="left", fontsize=8.5, color=INK)
    h, l = ax1.get_legend_handles_labels()
    order = sorted(range(len(l)), key=lambda i: (l[i] == "all haplotypes", l[i]))
    order = sorted(order, key=lambda i: (l[i] == "all haplotypes",
                                          int(l[i].split()[0]) if l[i][0].isdigit() else 0))
    ax1.legend([h[i] for i in order], [l[i] for i in order], frameon=False, fontsize=6.8, loc="lower left")
    style(ax1)

    edges = [round(0.30 + 0.05 * i, 2) for i in range(23)]          # 0.30 to 1.40
    cats = [("carrier in joined VCF", lambda k: k in called and k in vcar, BLUE),
            ("non carrier in joined VCF", lambda k: k in called and k not in vcar, GREY),
            ("not called in VCF", lambda k: k not in called, LIGHT)]
    bottoms = [0] * (len(edges) - 1)
    counts = {}
    for name, test, col in cats:
        h = [0] * (len(edges) - 1)
        for k, v in dos.items():
            if test(k):
                i = min(max(int((v - edges[0]) / 0.05), 0), len(h) - 1)
                h[i] += 1
        counts[name] = h
        ax2.bar([e + 0.025 for e in edges[:-1]], h, width=0.045, bottom=bottoms, color=col, label=name)
        bottoms = [b + x for b, x in zip(bottoms, h)]
    ax2.axvline(a_.cut, color=INK2, lw=0.8, ls=":")
    ax2.set_ylim(0, max(bottoms) * 1.55)
    ax2.text(a_.cut + 0.01, max(bottoms) * 0.04, f"cut {a_.cut}", fontsize=6.5, color=INK2)
    ax2.set_xlabel("Graph dosage of the deleted interval (GRCh38 = 1)", color=INK)
    ax2.set_ylabel("Haplotypes", color=INK)
    ax2.set_title("b  62 kb deletion, PGGB graph", loc="left", fontsize=8.5, color=INK)
    ax2.legend(frameon=False, fontsize=6.3, loc="upper left", ncol=1)
    ax2.yaxis.grid(True, color=GRID, linewidth=0.8)
    style(ax2)

    largest = max(r[1] for r in recs) if recs else 0
    both = len(vcar & gcar); gonly_called = len((gcar - vcar) & called); gonly_unc = len(gcar - called)
    rows = [("Largest VCF\nrecord", [(largest, BLUE)]),
            ("Joined VCF\nrecords", [(n, BLUE if i == 0 else "#86b6ef") for i, (_, n) in
                                     enumerate(sorted(recs, key=lambda r: -r[1]))]),
            ("Graph\n(dosage < cut)", [(both, BLUE), (gonly_called, GREY), (gonly_unc, LIGHT)])]
    yy = [2, 1, 0]
    for y, (lab, segs) in zip(yy, rows):
        left = 0
        for n, col in segs:
            ax3.barh(y, n, left=left, height=0.55, color=col, edgecolor="white", linewidth=1.0)
            left += n
        ax3.text(left + 4, y, f"{left}", va="center", fontsize=7, color=INK)
    ax3.set_yticks(yy); ax3.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax3.set_xlim(0, max(260, len(gcar) + 30))
    ax3.set_xlabel("Haplotypes carrying the deletion", color=INK)
    ax3.set_title("c  One event, three counts", loc="left", fontsize=8.5, color=INK)
    ax3.xaxis.grid(True, color=GRID, linewidth=0.8)
    style(ax3)

    for ext in ("pdf", "png"):
        fig.savefig(f"{a_.out}.{ext}", dpi=300)
    print(f"wrote {a_.out}.pdf and {a_.out}.png\n")

    print("a  insertions, share held by the largest allele (mean, min to max over draws)")
    for n in sorted(sub):
        print(f"   {n:4d} haplotypes, {len(sub[n])} draws: " + "  ".join(
            f"{l} {sum(d[b] for d in sub[n])/len(sub[n]):.3f} ({min(d[b] for d in sub[n]):.3f}-{max(d[b] for d in sub[n]):.3f})"
            for b, l in zip(BINS, LAB)))
    print("   all haplotypes: " + "  ".join(f"{l} {full[b]:.3f}" for b, l in zip(BINS, LAB)))
    print(f"\nb  haplotypes with one complete path: {len(dos)}; dosage < {a_.cut}: {len(gcar)}")
    for name in counts:
        print(f"   {name:28} {sum(counts[name])}")
    print(f"   carrier in joined VCF with dosage >= cut (VCF only): {len(((vcar & called) & set(dos)) - gcar)}")
    print(f"\nc  VCF records (position, carriers): {sorted(recs)}")
    print(f"   largest record {largest}; joined {len(vcar)} of {len(called)} called haplotypes")
    print(f"   graph carriers {len(gcar)}: also VCF carriers {both}, VCF non carriers {gonly_called}, "
          f"not called in VCF {gonly_unc}")


if __name__ == "__main__":
    main()
