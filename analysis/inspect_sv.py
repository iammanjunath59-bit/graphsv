#!/usr/bin/env python3
"""
inspect_sv.py

Show exactly how one simulated SV is represented in the post-processed VCF
versus in the graph.

Why
---
The aggregate benchmark says the VCF route loses 75% of inversions above
3 kb. That is a number, not a mechanism. To claim WHY, you have to look at
what actually appears at the locus. The expectation is that vcfwave, having
no length signal to anchor a balanced inversion, realigns it base-by-base
and emits a run of substitutions -- so the inversion does not merely lose
its allele count, it stops existing as a variant class and reappears as
spurious SNPs. That would also inflate SNP counts and distort the SNP
frequency spectrum, which is a second-order consequence worth measuring.

This prints the evidence for a single site so the claim can be stated
precisely rather than inferred.

Usage
-----
  python inspect_sv.py --prefix inv/chr20_15000000_s1 \
                       --sv-id SV0006          # or --type INV --min-len 5000
"""

import argparse
import sys


def load_truth(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            rows.append({k: c[i] for k, i in ix.items()})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--sv-id")
    ap.add_argument("--type", default=None)
    ap.add_argument("--min-len", type=int, default=0)
    ap.add_argument("--pad", type=int, default=100,
                    help="bp of flank to include when scanning the VCF")
    args = ap.parse_args()

    truth = load_truth(f"{args.prefix}.truth_svs.tsv")
    if args.sv_id:
        sel = [r for r in truth if r["sv_id"] == args.sv_id]
    else:
        sel = [r for r in truth
               if (args.type is None or r["type"] == args.type)
               and int(r["len"]) >= args.min_len]
    if not sel:
        sys.exit("no matching SV in truth table")
    sv = sel[0]

    s, e, L = int(sv["start"]), int(sv["end"]), int(sv["len"])
    print(f"TRUTH  {sv['sv_id']}  {sv['type']}  {s}-{e}  len={L}  "
          f"AC={sv['ac']}  carriers={sv['carriers']}")
    print()

    # ---- what the VCF contains over this interval ----------------------
    lo, hi = s - args.pad, e + args.pad
    recs = []
    with open(f"{args.prefix}.wave.vcf") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            pos = int(f[1])
            reflen = len(f[3])
            if pos + reflen < lo or pos > hi:
                continue
            alts = f[4].split(",")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            recs.append({
                "pos": pos, "reflen": reflen,
                "altlens": [len(a) for a in alts],
                "n_alt": len(alts),
                "ac": info.get("AC", ""), "an": info.get("AN", ""),
                "lv": info.get("LV", ""),
                "is_snp": reflen == 1 and all(len(a) == 1 for a in alts),
            })

    n_snp = sum(1 for r in recs if r["is_snp"])
    n_big = sum(1 for r in recs if max([r["reflen"]] + r["altlens"]) >= 50)
    print(f"VCF over {lo}-{hi}: {len(recs)} records "
          f"({n_snp} SNP-like, {n_big} with an allele >= 50 bp)")
    if n_snp > 10 and n_big == 0:
        print("  -> the event is present ONLY as a run of substitutions: "
              "it has been decomposed, not represented")
    print(f"  {'pos':>9} {'reflen':>7} {'altlens':>18} {'AC':>10} {'LV':>3}")
    for r in recs[:25]:
        al = ",".join(map(str, r["altlens"]))[:18]
        print(f"  {r['pos']:>9} {r['reflen']:>7} {al:>18} "
              f"{r['ac']:>10} {r['lv']:>3}")
    if len(recs) > 25:
        print(f"  ... {len(recs)-25} more")
    print()

    # ---- what the graph caller produced --------------------------------
    try:
        with open(f"{args.prefix}.gsv.tsv") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            gi = {c: i for i, c in enumerate(hdr)}
            hits = []
            for line in fh:
                c = line.rstrip("\n").split("\t")
                p, en = int(c[gi["pos"]]), int(c[gi["end"]])
                if en >= lo and p <= hi:
                    hits.append(c)
        print(f"GRAPH caller over the same interval: {len(hits)} site(s)")
        print("  " + "\t".join(hdr))
        for c in hits:
            print("  " + "\t".join(c))
    except FileNotFoundError:
        print("(no .gsv.tsv found)")


if __name__ == "__main__":
    main()
