#!/usr/bin/env python3
"""
analyse_sweep.py

Fit SV allele-count recovery against variant size and segmental-duplication
content, with window as a random effect.

Why a mixed model
-----------------
The naive regression treats each simulation replicate as independent. It
is not: two seeds share a window, and all 24 SVs within a replicate share
a genomic context and a graph build. Ignoring that clustering inflates the
apparent precision of the segdup slope -- which matters here, because the
top of the segdup range is carried by only ~5 windows. Window is therefore
a random intercept, and the segdup effect is estimated BETWEEN windows,
which is the level at which it actually varies.

Size enters as log10(bp), because recovery is flat to ~1 kb and then falls;
a linear term in raw bp would be dominated by the 30 kb bin.

Outputs
-------
  <out>.perSV.tsv        one row per SV: window, segdup, size, recovered
  <out>.cells.tsv        size bin x segdup stratum table with Wilson CIs
  <out>.model.txt        fitted models and coefficient tables
"""

import argparse
import glob
import math
import os
import sys


def wilson(k, n, z=1.96):
    """Wilson score interval -- correct near 0 and 1, unlike normal approx."""
    if n == 0:
        return (float("nan"),) * 3
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(0.0, c - h), min(1.0, c + h)


def load_segdup(path):
    sd = {}
    with open(path) as fh:
        for i, line in enumerate(fh):
            f = line.split()
            if len(f) < 4:
                continue
            if i == 0 and not f[1].isdigit():
                continue
            sd[f"{f[0]}_{f[1]}"] = float(f[3])
    return sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default="sweep54/*.wcmp.svs.tsv")
    ap.add_argument("--windows", default="data/eligible_windows.tsv")
    ap.add_argument("--out", default="sweep54/analysis")
    ap.add_argument("--hi-segdup", type=float, default=0.05,
                    help="threshold defining the 'high segdup' stratum")
    args = ap.parse_args()

    sd = load_segdup(args.windows)
    if not sd:
        sys.exit(f"no windows parsed from {args.windows}")

    rows = []
    for f in sorted(glob.glob(args.glob)):
        base = os.path.basename(f).replace(".wcmp.svs.tsv", "")
        win, _, seed = base.rpartition("_s")
        if win not in sd:
            continue
        with open(f) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(hdr)}
            for line in fh:
                c = line.rstrip("\n").split("\t")
                if len(c) < len(hdr):
                    continue
                try:
                    size = int(c[ix["len"]])
                    rec = int(c[ix["exact_ac_present"]])
                except (ValueError, KeyError):
                    continue
                rows.append({
                    "window": win, "seed": seed, "segdup": sd[win],
                    "sv_id": c[ix["sv_id"]], "type": c[ix["type"]],
                    "size_bin": c[ix["size_bin"]], "size": size,
                    "log10size": math.log10(max(size, 1)),
                    "recovered": rec,
                })
    if not rows:
        sys.exit("no SV rows loaded -- check --glob")

    print(f"loaded {len(rows)} SVs from "
          f"{len(set(r['window'] for r in rows))} windows, "
          f"{len(set((r['window'], r['seed']) for r in rows))} replicates")

    with open(f"{args.out}.perSV.tsv", "w") as fh:
        cols = ["window", "seed", "segdup", "sv_id", "type", "size_bin",
                "size", "log10size", "recovered"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    # ---- descriptive cell table with Wilson intervals -------------------
    cells = {}
    for r in rows:
        key = (r["size_bin"], r["segdup"] >= args.hi_segdup)
        n, k = cells.get(key, (0, 0))
        cells[key] = (n + 1, k + r["recovered"])

    lines = [f"{'bin':>7} {'stratum':>10} {'n':>5} {'recov':>7} "
             f"{'lo95':>7} {'hi95':>7}"]
    bins = sorted({r["size_bin"] for r in rows},
                  key=lambda b: int(b) if str(b).isdigit() else 0)
    for b in bins:
        for hi in (False, True):
            if (b, hi) not in cells:
                continue
            n, k = cells[(b, hi)]
            p, lo, up = wilson(k, n)
            lines.append(f"{b:>7} {'high' if hi else 'low':>10} {n:>5} "
                         f"{p:>7.3f} {lo:>7.3f} {up:>7.3f}")
    cell_txt = "\n".join(lines)
    print("\n" + cell_txt)
    with open(f"{args.out}.cells.tsv", "w") as fh:
        fh.write(cell_txt + "\n")

    # ---- models ---------------------------------------------------------
    out = [cell_txt, ""]
    try:
        import pandas as pd
        import statsmodels.formula.api as smf
        import statsmodels.api as sm
    except ImportError:
        sys.exit("\nInstall statsmodels to fit the models:\n"
                 "  mamba install -c conda-forge statsmodels pandas")

    df = pd.DataFrame(rows)
    df["segdup10"] = df["segdup"] * 10      # per 10% segdup, interpretable
    df["logsize_c"] = df["log10size"] - df["log10size"].mean()

    # 1. Naive GLM -- WRONG standard errors, shown only for contrast
    m0 = smf.glm("recovered ~ logsize_c + segdup10",
                 data=df, family=sm.families.Binomial()).fit()
    out += ["=== 1. Naive binomial GLM (ignores window clustering) ===",
            "Standard errors here are anticonservative: SVs within a window",
            "and replicates within a window are not independent.",
            str(m0.summary()), ""]

    # 2. GLM with cluster-robust SEs by window -- the honest headline model
    m1 = smf.glm("recovered ~ logsize_c + segdup10",
                 data=df, family=sm.families.Binomial()
                 ).fit(cov_type="cluster",
                       cov_kwds={"groups": df["window"]})
    out += ["=== 2. Binomial GLM, cluster-robust SEs by window ===",
            "Same point estimates, but SEs account for clustering. The",
            "segdup slope is identified BETWEEN windows, so its effective",
            "sample size is the number of windows, not the number of SVs.",
            str(m1.summary()), ""]

    # 3. Interaction
    m2 = smf.glm("recovered ~ logsize_c * segdup10",
                 data=df, family=sm.families.Binomial()
                 ).fit(cov_type="cluster",
                       cov_kwds={"groups": df["window"]})
    out += ["=== 3. Size x segdup interaction, cluster-robust ===",
            str(m2.summary()), ""]

    # 4. Mixed model on the replicate-level proportion
    agg = (df.groupby(["window", "seed", "size_bin"])
             .agg(segdup10=("segdup10", "first"),
                  logsize_c=("logsize_c", "mean"),
                  p=("recovered", "mean"), n=("recovered", "size"))
             .reset_index())
    try:
        m3 = smf.mixedlm("p ~ logsize_c + segdup10", data=agg,
                         groups=agg["window"]).fit()
        out += ["=== 4. LMM on replicate-level recovery proportion ===",
                "Random intercept per window. Coarser than the GLMs (it",
                "throws away within-cell counts) but makes the clustering",
                "explicit and reports between-window variance directly.",
                str(m3.summary()), ""]
    except Exception as e:                       # noqa: BLE001
        out += [f"mixedlm failed: {e}", ""]

    txt = "\n".join(out)
    with open(f"{args.out}.model.txt", "w") as fh:
        fh.write(txt + "\n")

    print("\n--- key coefficients (cluster-robust model 2) ---")
    for term in ("logsize_c", "segdup10"):
        if term in m1.params:
            b, se = m1.params[term], m1.bse[term]
            print(f"  {term:>10}  beta={b:+.3f}  SE={se:.3f}  "
                  f"z={b/se:+.2f}  p={m1.pvalues[term]:.4g}  "
                  f"OR={math.exp(b):.3f}")
    print(f"\nwrote {args.out}.perSV.tsv, .cells.tsv, .model.txt")


if __name__ == "__main__":
    main()
