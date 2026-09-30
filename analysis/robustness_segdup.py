#!/usr/bin/env python3
"""
robustness_segdup.py

Leverage diagnostics for the segmental-duplication effect.

Why this is needed
------------------
The segdup slope is estimated BETWEEN windows, and the covariate range is
badly unbalanced: ~5 of 54 windows carry the entire upper half of the
segdup axis (40 SVs vs 392 per size bin). Cluster-robust standard errors
assume many clusters and no small subset dominating the covariate -- the
second assumption is violated here. A significant z is therefore not by
itself evidence that the effect would survive dropping one window.

Three checks:

  1. Leave-one-window-out    -- refit 54 times, one window dropped each
                                time. Report the range of the segdup
                                coefficient. If dropping a single window
                                flips the sign or kills significance, the
                                effect is one window wearing a trend coat.

  2. Leave-top-k-out         -- drop the k highest-segdup windows in turn.
                                The honest stress test, since those are
                                the leverage points.

  3. Window-level permutation -- shuffle segdup labels ACROSS WINDOWS
                                (not across SVs) and refit. This is the
                                correct null: it preserves within-window
                                structure and size effects, and destroys
                                only the segdup association. Compare the
                                observed coefficient to this distribution.
                                Robust to any residual clustering the
                                sandwich estimator mishandles.

Usage
-----
  python robustness_segdup.py --perSV out/an20.perSV.tsv --out out/robustness20
"""

import argparse
import random
import sys

try:
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
except ImportError:
    sys.exit("needs pandas/statsmodels/numpy")


FORMULA = "recovered ~ logsize_c + segdup10"


def fit(df, robust=True):
    """Return (beta_segdup, se, z) for the cluster-robust binomial GLM."""
    if df["segdup10"].nunique() < 2:
        return (np.nan,) * 3
    m = smf.glm(FORMULA, data=df, family=sm.families.Binomial())
    if robust and df["window"].nunique() > 2:
        res = m.fit(cov_type="cluster", cov_kwds={"groups": df["window"]})
    else:
        res = m.fit()
    b = res.params.get("segdup10", np.nan)
    se = res.bse.get("segdup10", np.nan)
    return b, se, (b / se if se else np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perSV", required=True,
                    help="perSV table written by analyse_sweep.py")
    ap.add_argument("--nperm", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="robustness",
                    help="output prefix; writes <out>.txt and <out>.loo.tsv")
    args = ap.parse_args()

    df = pd.read_csv(args.perSV, sep="\t")
    df["logsize_c"] = df["log10size"] - df["log10size"].mean()
    # perSV.tsv carries raw segdup fraction; rescale to "per 10 percentage
    # points" so the coefficient matches the main analysis.
    if "segdup10" not in df.columns:
        df["segdup10"] = df["segdup"] * 10
    wins = sorted(df["window"].unique())
    wsd = df.groupby("window")["segdup10"].first().sort_values(ascending=False)

    b0, se0, z0 = fit(df)
    lines = [f"full model: beta_segdup = {b0:+.4f}  SE = {se0:.4f}  "
             f"z = {z0:+.2f}",
             f"windows = {len(wins)}, SVs = {len(df)}",
             f"segdup range = {df['segdup10'].min():.4f} "
             f"to {df['segdup10'].max():.4f} (per 10%)", ""]

    # ---- 1. leave-one-window-out ---------------------------------------
    loo = []
    for w in wins:
        b, se, z = fit(df[df["window"] != w])
        loo.append((w, b, z))
    bs = np.array([b for _, b, _ in loo])
    zs = np.array([z for _, _, z in loo])
    lines += ["=== 1. Leave-one-window-out (n=54 refits) ===",
              f"beta range: {bs.min():+.4f} to {bs.max():+.4f}",
              f"z    range: {zs.min():+.2f} to {zs.max():+.2f}",
              f"refits with z > -1.96: {(zs > -1.96).sum()} / {len(zs)}",
              f"refits with beta >= 0: {(bs >= 0).sum()} / {len(bs)}", ""]
    infl = sorted(loo, key=lambda t: abs(t[1] - b0), reverse=True)[:5]
    lines.append("most influential windows (largest shift in beta):")
    for w, b, z in infl:
        lines.append(f"  drop {w:>20}  segdup10={wsd.get(w, float('nan')):.4f}"
                     f"  beta={b:+.4f} ({b - b0:+.4f})  z={z:+.2f}")
    lines.append("")

    # ---- 2. leave-top-k-out --------------------------------------------
    lines.append("=== 2. Drop the k highest-segdup windows ===")
    lines.append("This is the real stress test: those windows ARE the "
                 "upper covariate range.")
    for k in range(0, 6):
        drop = set(wsd.index[:k])
        sub = df[~df["window"].isin(drop)]
        b, se, z = fit(sub)
        lines.append(f"  k={k}  max segdup10 remaining="
                     f"{sub['segdup10'].max():.4f}  "
                     f"beta={b:+.4f}  SE={se:.4f}  z={z:+.2f}  "
                     f"n_win={sub['window'].nunique()}")
    lines.append("")

    # ---- 3. window-level permutation -----------------------------------
    rng = random.Random(args.seed)
    obs = b0
    wlist = list(wsd.index)
    vals = list(wsd.values)
    null = []
    for _ in range(args.nperm):
        perm = vals[:]
        rng.shuffle(perm)
        mapping = dict(zip(wlist, perm))
        d = df.copy()
        d["segdup10"] = d["window"].map(mapping)
        b, _, _ = fit(d, robust=False)   # SEs irrelevant under permutation
        null.append(b)
    null = np.array(null)
    p_perm = (np.sum(null <= obs) + 1) / (len(null) + 1)
    lines += ["=== 3. Window-level permutation of segdup labels ===",
              f"permutations = {args.nperm}",
              f"observed beta = {obs:+.4f}",
              f"null mean = {null.mean():+.4f}, "
              f"null sd = {null.std():.4f}",
              f"null 2.5th / 97.5th pct = "
              f"{np.percentile(null, 2.5):+.4f} / "
              f"{np.percentile(null, 97.5):+.4f}",
              f"one-sided p (beta <= observed) = {p_perm:.4g}", ""]

    txt = "\n".join(lines)
    print(txt)
    with open(f"{args.out}.txt", "w") as fh:
        fh.write(txt + "\n")
    pd.DataFrame(loo, columns=["dropped_window", "beta", "z"]).to_csv(
        f"{args.out}.loo.tsv", sep="\t", index=False)
    print(f"wrote {args.out}.txt and {args.out}.loo.tsv")


if __name__ == "__main__":
    main()
