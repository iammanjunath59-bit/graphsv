#!/usr/bin/env python3
"""
preflight_checks.py

Three statistical checks that a reviewer will ask for, run before any
corresponding claim is written into the manuscript.

1. COLLINEARITY between variant size and segmental duplication content.
   In real genomes large structural variants are enriched in duplicated
   sequence, because non-allelic homologous recombination requires it, and a
   regression on real data could not separate the two. The stratified design
   used here decouples them by construction: four variants are placed in
   every size class in every window, regardless of that window's duplication
   content. This check confirms that the decoupling worked, by reporting the
   correlation between mean variant size and duplication content across
   windows and the variance inflation factors of the two predictors. It also
   refits with an interaction term to confirm the main effects are stable.

2. GC CONTENT as a confounder. wfmash and minimap2 anchor density depends on
   base composition, and chromosome 16 is GC-rich relative to chromosome 20,
   so a held-out validation on chromosome 16 could in principle be confounded
   by alignment difficulty rather than structural complexity. This refits the
   recovery model with mean GC as an additional covariate.

3. SIZE CLASS COMPOSITION relative to the 50 bp segment threshold used in
   traversal comparison. Variants shorter than the threshold cannot be
   resolved by the graph route, so the actual size distribution within the
   smallest class needs to be reported rather than assumed.

Nothing here is written to the manuscript automatically. Read the numbers,
then write what they say.

Usage
-----
  python preflight_checks.py --sweep sweep54d --windows data/eligible_windows.tsv \\
                             --nuc data/nuc.txt --label chr20
  python preflight_checks.py --sweep sweep16 --windows data/eligible_chr16.tsv \\
                             --nuc data/nuc16.txt --label chr16
"""

import argparse
import glob
import math
import os
import sys
from collections import defaultdict

try:
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
except ImportError:
    sys.exit("needs numpy, pandas, statsmodels")


def load_windows(path):
    """window id -> segdup fraction"""
    sd = {}
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 4 and f[1].isdigit():
                sd[f"{f[0]}_{f[1]}"] = float(f[3])
    return sd


def load_gc(path):
    """window id -> GC fraction, from `bedtools nuc` output (no header)."""
    gc = {}
    if not path or not os.path.exists(path):
        return gc
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or not f[1].isdigit():
                continue
            # bedtools nuc: chrom start end ... pct_at pct_gc ...
            try:
                gc[f"{f[0]}_{f[1]}"] = float(f[4])
            except ValueError:
                continue
    return gc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--nuc", default=None,
                    help="bedtools nuc output for the same windows")
    ap.add_argument("--label", default="")
    ap.add_argument("--min-seg", type=int, default=50)
    args = ap.parse_args()

    sd = load_windows(args.windows)
    gc = load_gc(args.nuc)

    rows = []
    for f in sorted(glob.glob(f"{args.sweep}/*.vcfcmp.svs.tsv")):
        base = os.path.basename(f).replace(".vcfcmp.svs.tsv", "")
        win = base.rsplit("_s", 1)[0]
        if win not in sd:
            continue
        gsv = f.replace(".vcfcmp.svs.tsv", ".gsv.eval.tsv")
        graph = {}
        if os.path.exists(gsv):
            with open(gsv) as gh:
                gh.readline()
                for l in gh:
                    c = l.rstrip("\n").split("\t")
                    graph[c[0]] = int(c[7])
        with open(f) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(hdr)}
            for line in fh:
                c = line.rstrip("\n").split("\t")
                try:
                    L = int(c[ix["len"]])
                except ValueError:
                    continue
                rows.append({
                    "window": win, "sv_id": c[ix["sv_id"]],
                    "segdup10": sd[win] * 10,
                    "gc": gc.get(win, float("nan")),
                    "len": L, "bin": c[ix["size_bin"]],
                    "log10size": math.log10(max(L, 1)),
                    "vcf_ok": int(c[ix["exact_ac_anywhere"]])
                    if "exact_ac_anywhere" in ix else int(c[14]),
                    "graph_ok": graph.get(c[ix["sv_id"]], np.nan),
                })
    if not rows:
        sys.exit(f"no rows loaded from {args.sweep}")
    df = pd.DataFrame(rows)
    df["logsize_c"] = df["log10size"] - df["log10size"].mean()
    tag = f"[{args.label}] " if args.label else ""
    print(f"{tag}{len(df)} variants, {df['window'].nunique()} windows\n")

    # ---- 1. collinearity ---------------------------------------------
    print("=== 1. Collinearity between size and segmental duplication ===")
    w = df.groupby("window").agg(mean_len=("len", "mean"),
                                 mean_log=("log10size", "mean"),
                                 segdup10=("segdup10", "first")).reset_index()
    r_p = w["mean_len"].corr(w["segdup10"])
    r_s = w["mean_len"].corr(w["segdup10"], method="spearman")
    r_log = w["mean_log"].corr(w["segdup10"])
    print(f"  across windows: Pearson r(mean SV size, segdup) = {r_p:+.4f}")
    print(f"                  Spearman rho                    = {r_s:+.4f}")
    print(f"                  Pearson r(mean log10 size)      = {r_log:+.4f}")

    # VIF at the level the model is actually fitted (per variant)
    X = df[["logsize_c", "segdup10"]].copy()
    X["const"] = 1.0
    def vif(col):
        y = X[col]
        others = X.drop(columns=[col])
        fit = sm.OLS(y, others).fit()
        return 1.0 / (1.0 - fit.rsquared) if fit.rsquared < 1 else float("inf")
    print(f"  VIF logsize_c = {vif('logsize_c'):.4f}")
    print(f"  VIF segdup10  = {vif('segdup10'):.4f}")
    print("  (VIF near 1 means the predictors are effectively orthogonal,")
    print("   which the stratified design is intended to guarantee)")

    # ---- interaction ---------------------------------------------------
    m_main = smf.glm("vcf_ok ~ logsize_c + segdup10", data=df,
                     family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": df["window"]})
    m_int = smf.glm("vcf_ok ~ logsize_c * segdup10", data=df,
                    family=sm.families.Binomial()).fit(
        cov_type="cluster", cov_kwds={"groups": df["window"]})
    print("\n  main-effects model:")
    for t in ("logsize_c", "segdup10"):
        print(f"    {t:>22}  beta={m_main.params[t]:+.4f}  "
              f"z={m_main.tvalues[t]:+.2f}  p={m_main.pvalues[t]:.3g}")
    print("  with interaction:")
    for t in m_int.params.index:
        if t == "Intercept":
            continue
        print(f"    {t:>22}  beta={m_int.params[t]:+.4f}  "
              f"z={m_int.tvalues[t]:+.2f}  p={m_int.pvalues[t]:.3g}")

    # ---- 2. GC as a covariate -----------------------------------------
    print("\n=== 2. GC content as a confounder ===")
    if df["gc"].isna().all():
        print("  no GC data supplied (--nuc); skipped")
    else:
        print(f"  GC range across windows: {df['gc'].min():.4f} "
              f"to {df['gc'].max():.4f}, mean {df['gc'].mean():.4f}")
        print(f"  r(GC, segdup) across windows = "
              f"{df.groupby('window')['gc'].first().corr(df.groupby('window')['segdup10'].first()):+.4f}")
        d2 = df.dropna(subset=["gc"]).copy()
        d2["gc_c"] = d2["gc"] - d2["gc"].mean()
        m_gc = smf.glm("vcf_ok ~ logsize_c + segdup10 + gc_c", data=d2,
                       family=sm.families.Binomial()).fit(
            cov_type="cluster", cov_kwds={"groups": d2["window"]})
        for t in ("logsize_c", "segdup10", "gc_c"):
            if t in m_gc.params:
                print(f"    {t:>12}  beta={m_gc.params[t]:+.4f}  "
                      f"z={m_gc.tvalues[t]:+.2f}  p={m_gc.pvalues[t]:.3g}")
        print("  compare segdup10 without GC: "
              f"{m_main.params['segdup10']:+.4f}")

    # ---- 3. size class composition ------------------------------------
    print("\n=== 3. Size composition relative to the "
          f"{args.min_seg} bp threshold ===")
    for b in sorted(df["bin"].unique(),
                    key=lambda x: int(x) if str(x).isdigit() else 0):
        g = df[df["bin"] == b]
        below = (g["len"] < args.min_seg).sum()
        print(f"  class {b:>6}  n={len(g):4d}  "
              f"len min={g['len'].min():6d}  median={int(g['len'].median()):6d}  "
              f"max={g['len'].max():6d}  below threshold={below}")
    tot_below = (df["len"] < args.min_seg).sum()
    print(f"\n  variants shorter than {args.min_seg} bp overall: "
          f"{tot_below}/{len(df)}")
    if tot_below == 0:
        print("  No simulated variant falls below the threshold, so the")
        print("  graph route discards none of them; the threshold sets a")
        print("  floor on what the method can resolve but does not remove")
        print("  anything from this benchmark.")


if __name__ == "__main__":
    main()
