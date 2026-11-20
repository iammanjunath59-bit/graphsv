#!/usr/bin/env python3
"""
score_graphsv.py

Scores graphsv calls against simulation truth under the criteria used for
the VCF routes, so that the routes can be compared like for like.

For each true variant of length L, the calls overlapping it within
max(200 bp, L/2) are examined (the overlap rule of graphsv's own
evaluation and of precision_recall_types.py). Three criteria:

  count     the true allele count equals the count of an overlapping ALT
            allele or the call's total ALT count (graphsv's built in
            evaluation; no check of type or carriers)
  typed     an overlapping ALT allele of the matching type and size carries
            the true count: a deletion allele is shorter than the reference
            traversal by L, a duplication allele longer by L, each within
            max(50 bp, 20 % of L), and an inversion allele has the same
            length within 50 bp. Lengths are those graphsv reports, which
            count only segments of at least --min-seg bp, and merged calls
            report the merged span, so this criterion can fail on a correct
            count.
  carrier   an overlapping ALT allele is carried by exactly the true
            haplotypes. Requires calls written with graphsv --carriers;
            haplotype paths are matched to truth carrier indices by the
            digits of the sample field of the PanSN name (sim05#1#... -> 5),
            as in precision_recall_types.py.

Inputs: <truth-dir>/<replicate>.truth_svs.tsv and
        <calls-dir>/<replicate>.gsv.tsv for every replicate with both.
Outputs: <out>.perSV.tsv and a printed summary by type and size class.

Usage:
    python3 analysis/score_graphsv.py --truth-dir results/sweep54d \\
        --calls-dir out/gsv_car --out out/gsv_car
"""

import argparse
import csv
import glob
import os
from collections import defaultdict

BINS = ["100", "300", "1000", "3000", "10000", "30000"]


def sample_index(name):
    d = "".join(ch for ch in name.split("#")[0] if ch.isdigit())
    return int(d) if d else None


def score_replicate(truth_path, calls_path):
    with open(calls_path) as fh:
        calls = list(csv.DictReader(fh, delimiter="\t"))
    has_car = bool(calls) and "allele_carriers" in calls[0]
    rows = []
    with open(truth_path) as fh:
        for t in csv.DictReader(fh, delimiter="\t"):
            s, e, L, ac, ty = int(t["start"]), int(t["end"]), int(t["len"]), int(t["ac"]), t["type"]
            truth_car = {int(x) for x in t["carriers"].split(",") if x != ""}
            slack = max(200, L // 2)
            tol = max(50, 0.2 * L)
            count = typed = carrier = 0
            for c in calls:
                if not (int(c["end"]) > s - slack and int(c["pos"]) < e + slack):
                    continue
                acs = [int(x) for x in c["allele_acs"].split(",") if x]
                lens = [int(x) for x in c["allele_lens"].split(",") if x]
                R = int(c["ref_len"])
                if ac in acs or int(c["ac_total"]) == ac:
                    count = 1
                for a, n in zip(lens, acs):
                    d = a - R
                    ok = {"DEL": abs(-d - L) <= tol, "DUP": abs(d - L) <= tol,
                          "INV": abs(d) < 50}.get(ty, False)
                    if ok and n == ac:
                        typed = 1
                if has_car:
                    for grp in c["allele_carriers"].split(","):
                        got = {sample_index(x) for x in grp.split(";") if x}
                        if got == truth_car:
                            carrier = 1
            rows.append({"sv_id": t["sv_id"], "type": ty, "len": L, "size_bin": t["size_bin"],
                         "ac_true": ac, "count": count, "typed": typed,
                         "carrier": carrier if has_car else ""})
    return rows, has_car


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--calls-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    allrows, reps, car_all = [], 0, True
    for tp in sorted(glob.glob(os.path.join(a.truth_dir, "*.truth_svs.tsv"))):
        rep = os.path.basename(tp)[: -len(".truth_svs.tsv")]
        cp = os.path.join(a.calls_dir, rep + ".gsv.tsv")
        if not os.path.exists(cp):
            continue
        rows, has_car = score_replicate(tp, cp)
        car_all &= has_car
        reps += 1
        for r in rows:
            r["replicate"] = rep
        allrows += rows
    if not allrows:
        raise SystemExit("no replicate had both a truth table and a call file")

    cols = ["replicate", "sv_id", "type", "len", "size_bin", "ac_true", "count", "typed", "carrier"]
    with open(a.out + ".perSV.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in allrows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    crit = ["count", "typed"] + (["carrier"] if car_all else [])
    print(f"replicates {reps}   true variants {len(allrows)}"
          + ("" if car_all else "   (no allele_carriers column: carrier criterion not scored)"))
    for k in crit:
        n = sum(int(r[k]) for r in allrows)
        print(f"  {k:8} {n}/{len(allrows)} {n/len(allrows):.4f}")
    print("by type and size class (" + ", ".join(crit) + ")")
    for ty in ("DEL", "DUP", "INV"):
        cells = []
        for b in BINS:
            sub = [r for r in allrows if r["type"] == ty and r["size_bin"] == b]
            if sub:
                cells.append(f"{b}: " + "/".join(str(sum(int(r[k]) for r in sub)) for k in crit) + f" of {len(sub)}")
        sub = [r for r in allrows if r["type"] == ty]
        print(f"  {ty} " + "  ".join(cells) + "   all: " +
              "/".join(str(sum(int(r[k]) for r in sub)) for k in crit) + f" of {len(sub)}")
    print(f"wrote {a.out}.perSV.tsv")


if __name__ == "__main__":
    main()
