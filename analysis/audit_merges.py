#!/usr/bin/env python3
"""
audit_merges.py

Does carrier-set merging ever fuse two distinct structural variants?

The concern
-----------
graphsv merges adjacent bubbles that share an identical carrier set within
50 kb, on the reasoning that the two breakpoints of one large event are by
definition carried by the same haplotypes. With only 20 haplotypes the
genealogy contains few distinct clades, so two INDEPENDENT variants can land
on the same branch by chance and therefore share a carrier set exactly. If
that happens within the merge window, two real variants are fused into one
call. Because the fused call inherits a carrier set that is correct for both,
the permissive scoring criterion would count it as a successful recovery for
whichever variant it is matched to, while the other is silently absorbed.

That would inflate the graph route's reported recovery.

What this checks
----------------
For every merged call, how many distinct truth variants fall inside its
span. A merged call spanning exactly one truth variant is a correct
reassembly of a fragmented event. A merged call spanning two or more is a
false fusion.

Also reports the converse: truth variants that are split across more than one
call, which merging is supposed to prevent.

Usage
-----
  python audit_merges.py --dir sweep54d
  python audit_merges.py --dir sweep16
"""

import argparse
import glob
import os
from collections import Counter


def load_truth(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            rows.append({
                "id": c[ix["sv_id"]],
                "start": int(c[ix["start"]]),
                "end": int(c[ix["end"]]),
                "len": int(c[ix["len"]]),
                "ac": int(c[ix["ac"]]),
                "carriers": frozenset(c[ix["carriers"]].split(",")),
                "bin": c[ix["size_bin"]] if "size_bin" in ix else "",
            })
    return rows


def load_calls(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            try:
                rows.append({
                    "pos": int(c[ix["pos"]]),
                    "end": int(c[ix["end"]]),
                    "ref_len": int(c[ix["ref_len"]]),
                    "ac": int(c[ix["ac_total"]]),
                    "n_alleles": int(c[ix["n_alleles"]]),
                })
            except (ValueError, KeyError):
                continue
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="sweep54d")
    ap.add_argument("--slack", type=int, default=200,
                    help="bp tolerance when testing whether a truth variant "
                         "lies inside a call")
    args = ap.parse_args()

    n_rep = 0
    calls_total = 0
    span_counts = Counter()          # truth variants per call
    fusion_examples = []
    split_counts = Counter()         # calls per truth variant
    split_examples = []
    same_carrier_pairs = 0           # truth pairs at risk of false fusion
    at_risk_pairs = 0                # ... and within the merge window

    for tf in sorted(glob.glob(f"{args.dir}/*.truth_svs.tsv")):
        base = tf.replace(".truth_svs.tsv", "")
        cf = base + ".gsv.tsv"
        if not os.path.exists(cf):
            continue
        n_rep += 1
        truth = load_truth(tf)
        calls = load_calls(cf)
        calls_total += len(calls)

        # --- how often could a false fusion happen at all? -------------
        for i in range(len(truth)):
            for j in range(i + 1, len(truth)):
                a, b = truth[i], truth[j]
                if a["carriers"] == b["carriers"]:
                    same_carrier_pairs += 1
                    gap = max(a["start"], b["start"]) - min(a["end"], b["end"])
                    if gap <= 50000:
                        at_risk_pairs += 1

        # --- truth variants inside each call ---------------------------
        for c in calls:
            lo, hi = c["pos"] - args.slack, c["end"] + args.slack
            inside = [t for t in truth if t["start"] >= lo and t["end"] <= hi]
            span_counts[len(inside)] += 1
            if len(inside) >= 2:
                fusion_examples.append(
                    (os.path.basename(base), c["pos"], c["end"],
                     c["ref_len"], c["ac"],
                     [(t["id"], t["len"], t["ac"]) for t in inside]))

        # --- calls overlapping each truth variant ----------------------
        for t in truth:
            lo, hi = t["start"] - args.slack, t["end"] + args.slack
            over = [c for c in calls if c["end"] > lo and c["pos"] < hi]
            split_counts[len(over)] += 1
            if len(over) >= 2:
                split_examples.append(
                    (os.path.basename(base), t["id"], t["len"], t["ac"],
                     len(over)))

    print(f"replicates: {n_rep}")
    print(f"graph calls: {calls_total}")
    print()

    print("=== truth variants contained in each call ===")
    tot = sum(span_counts.values())
    for k in sorted(span_counts):
        print(f"  {k} truth variant(s): {span_counts[k]:5d} calls "
              f"({span_counts[k]/max(1,tot):.4f})")
    fused = sum(v for k, v in span_counts.items() if k >= 2)
    print(f"\n  FALSE FUSION RATE: {fused}/{tot} = {fused/max(1,tot):.4f}")
    print("  (calls spanning two or more distinct truth variants)")

    print("\n=== calls overlapping each truth variant ===")
    tot2 = sum(split_counts.values())
    for k in sorted(split_counts):
        print(f"  {k} call(s): {split_counts[k]:5d} variants "
              f"({split_counts[k]/max(1,tot2):.4f})")
    split = sum(v for k, v in split_counts.items() if k >= 2)
    print(f"\n  RESIDUAL SPLIT RATE: {split}/{tot2} = {split/max(1,tot2):.4f}")
    print("  (truth variants still represented by more than one call)")

    print("\n=== opportunity for false fusion in the truth data ===")
    print(f"  truth pairs sharing an identical carrier set: {same_carrier_pairs}")
    print(f"  ... and separated by <= 50 kb (at risk): {at_risk_pairs}")
    print("  If this is zero, merging cannot have fused distinct variants in")
    print("  this design, and the merge criterion is untested rather than safe.")

    if fusion_examples:
        print("\n=== example false fusions (up to 10) ===")
        for e in fusion_examples[:10]:
            print(f"  {e[0]}  call {e[1]}-{e[2]} len={e[3]} AC={e[4]}")
            for sid, L, ac in e[5]:
                print(f"      contains {sid} len={L} AC={ac}")

    if split_examples:
        print("\n=== example residual splits (up to 10) ===")
        for e in split_examples[:10]:
            print(f"  {e[0]}  {e[1]} len={e[2]} AC={e[3]} -> {e[4]} calls")


if __name__ == "__main__":
    main()
