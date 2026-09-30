#!/usr/bin/env python3
"""
plot_supp_fig1.py

Supplementary Figure 1: the contribution of vcfbub and vcfwave to structural
variant allele count recovery.

The Methods commit to reporting recovery from the raw deconstruction output
separately, "since these steps exist specifically to resolve nested bubble
structure and their contribution should be measured rather than assumed".
This figure is that measurement.

Both panels use the evaluation criterion of the main analysis, computed by
compare_truth_graph.py, so the numbers are directly comparable to the primary
results rather than to a separate metric defined only for the supplement.

  Panel A   exact_ac_present -- the true allele count appears among the
            alternate allele counts at an overlapping site. Permissive: a
            site carrying many alternates passes if the correct count is
            among them.

  Panel B   clean recovery -- the same, additionally requiring the site to
            carry at most two alleles, so the count is readable without
            already knowing the answer.

The gap between the panels is the point: post-processing changes little in A
and a great deal in B.

Inputs are the per-SV tables written by compare_truth_graph.py:
    raw, bub   figs1/*.{raw,bub}.cmp.svs.tsv
    wave       sweep54d/*.vcfcmp.svs.tsv     (already present from the sweep)

Usage:
    python3 plot_supp_fig1.py --figs-dir figs1 --sweep-dir sweep54d \\
        --out figs1/supp_fig1.png
"""

import argparse
import collections
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STAGES = ("raw", "bub", "wave")
LABEL = {"raw": "vg deconstruct", "bub": "+ vcfbub", "wave": "+ vcfwave"}
STYLE = {"raw":  ("#1b4965", "-",  "o"),
         "bub":  ("#5fa8d3", "--", "s"),
         "wave": ("#c1666b", "-",  "^")}


def read_stage(pattern):
    """-> present{size_bin: [0/1]}, clean{size_bin: [0/1]}, n_files"""
    present = collections.defaultdict(list)
    clean = collections.defaultdict(list)
    files = sorted(glob.glob(pattern))
    for path in files:
        rows = [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]
        if len(rows) < 2:
            continue
        h = {n: i for i, n in enumerate(rows[0])}
        for need in ("size_bin", "exact_ac_present", "n_alleles"):
            if need not in h:
                sys.exit(f"{path} lacks column {need}; found {rows[0]}")
        for r in rows[1:]:
            if len(r) <= max(h.values()):
                continue
            try:
                b = int(r[h["size_bin"]])
            except ValueError:
                continue
            ex = int(r[h["exact_ac_present"]])
            na = r[h["n_alleles"]]
            na = int(na) if na not in ("", "NA") else 99
            present[b].append(ex)
            # matches compare_truth_graph.py: clean = exact AND n_alleles <= 2
            clean[b].append(1 if (ex and na <= 2) else 0)
    return present, clean, len(files)


def frac(d, b):
    v = d.get(b, [])
    return sum(v) / len(v) if v else float("nan")


def wilson(k, n, z=1.96):
    """Wilson score interval. Binomial only -- within-window clustering is not
    modelled here, so these bands are narrower than the cluster-robust
    intervals used for the main results."""
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5) / d
    return max(0.0, c - h), min(1.0, c + h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--figs-dir", default="figs1")
    ap.add_argument("--sweep-dir", default="sweep54d")
    ap.add_argument("--out", default="figs1/supp_fig1.png")
    args = ap.parse_args()

    patterns = {
        "raw": os.path.join(args.figs_dir, "*.raw.cmp.svs.tsv"),
        "bub": os.path.join(args.figs_dir, "*.bub.cmp.svs.tsv"),
        "wave": os.path.join(args.sweep_dir, "*.vcfcmp.svs.tsv"),
    }
    data = {}
    for st, pat in patterns.items():
        pres, cln, nf = read_stage(pat)
        if not pres:
            sys.exit(f"no data for stage {st} from {pat}")
        data[st] = (pres, cln, nf)
        sys.stderr.write(f"{st:5} {nf:3d} replicates, "
                         f"{sum(len(v) for v in pres.values())} SVs\n")

    bins = sorted(set().union(*[set(d[0]) for d in data.values()]))
    x = range(len(bins))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    titles = ["A   True allele count present at an overlapping site",
              "B   ...and site is biallelic, so the count is readable"]
    for ax, which, title in zip(axes, (0, 1), titles):
        for st in STAGES:
            d = data[st][which]
            y = [frac(d, b) for b in bins]
            lo = [wilson(sum(d.get(b, [])), len(d.get(b, [])))[0] for b in bins]
            hi = [wilson(sum(d.get(b, [])), len(d.get(b, [])))[1] for b in bins]
            c, ls, mk = STYLE[st]
            ax.plot(x, y, ls, color=c, marker=mk, ms=5, lw=1.7, label=LABEL[st])
            ax.fill_between(x, lo, hi, color=c, alpha=0.12, lw=0)
        ax.set_xticks(list(x))
        ax.set_xticklabels([f"{b // 1000}k" if b >= 1000 else str(b) for b in bins])
        ax.set_xlabel("simulated SV size (bp)")
        ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.25, lw=0.5)
        ax.set_title(title, loc="left", fontsize=10)
    axes[0].set_ylabel("fraction of simulated SVs")
    axes[0].legend(frameon=False, loc="lower left", fontsize=9)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        p = os.path.splitext(args.out)[0] + "." + ext
        fig.savefig(p, dpi=300)
        print("wrote", p)

    print("\nvalues plotted (fraction of SVs, by size bin):")
    for which, name in ((0, "A ac_present"), (1, "B clean")):
        print(f"\n  {name}")
        for st in STAGES:
            d = data[st][which]
            cells = "  ".join(
                f"{b if b < 1000 else str(b // 1000) + 'k'}:{frac(d, b):.3f}"
                for b in bins)
            tot = sum(sum(d.get(b, [])) for b in bins)
            n = sum(len(d.get(b, [])) for b in bins)
            print(f"    {LABEL[st]:16} {cells}   overall={tot / n:.4f}")


if __name__ == "__main__":
    main()
