#!/usr/bin/env python3
"""
classify_failures.py

When the VCF route fails to yield the correct allele count, what actually
went wrong?

Why this is needed
------------------
Aggregate recovery rates say how often the route fails but not how. Two
failures inspected by hand turned out to be qualitatively different:

  chr20_45000000_s1 SV0000, a 22,648 bp deletion at AC 6, appears in the
  post-processed VCF only as 31 single nucleotide substitutions. No record in
  the interval carries an allele of 50 bp or more. The variant has ceased to
  exist as a structural variant.

  chr20_21000000_s2 SV0022, a 42,000 bp deletion at AC 10, is delimited
  essentially perfectly: one record with a reference allele of 42,001 bp.
  But it reports AC 4. The coordinates are right and the frequency is wrong,
  with nothing in the record to signal it.

The second is the more consequential in practice, because a record that looks
correct will be used without question. Describing both as "fragmentation",
as we initially did, is inaccurate: fragmentation alone would not even cause
a scoring failure here, since the evaluation criterion accepts the true
allele count at any overlapping record, and fragments of one event share a
carrier set and would each report the correct count.

Classification
--------------
For each failed variant, the interval is examined and assigned to:

  ABSENT        no record in the interval carries an allele >= min_sv_len.
                The structural variant is not represented at all.

  MISCOUNTED    a record spans the variant with a size within tolerance, but
                no allele at it carries the true count.

  SPLIT         several large records overlap the variant, none matching its
                size within tolerance, and none carrying the true count.

  UNDERSIZED    large records are present but all are far smaller than the
                variant, consistent with only part of the event surviving.

Successes are classified for contrast: EXACT_SINGLE where one record of the
right size carries the right count, and EXACT_MULTI where the right count is
present but among several records or several alleles, meaning it could not be
identified without prior knowledge of the answer.

Usage
-----
  python classify_failures.py --dir sweep54d --out sweep54d/failure_modes
"""

import argparse
import glob
import os
from collections import Counter, defaultdict


def read_vcf(path, min_sv_len):
    """Large records only: (pos0, end0, size, [acs])."""
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            pos0 = int(f[1]) - 1
            ref = f[3]
            alts = f[4].split(",")
            size = max([len(ref)] + [len(a) for a in alts])
            if size < min_sv_len:
                continue
            info = {}
            for kv in f[7].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info[k] = v
            acs = []
            for a in info.get("AC", "").split(","):
                try:
                    acs.append(int(a))
                except ValueError:
                    pass
            out.append((pos0, pos0 + len(ref), size, acs))
    out.sort()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="sweep54d")
    ap.add_argument("--min-sv-len", type=int, default=50)
    ap.add_argument("--size-tol", type=float, default=0.25,
                    help="a record 'matches' the variant size within this "
                         "relative tolerance")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    modes = defaultdict(Counter)      # size_bin -> mode counts
    examples = defaultdict(list)
    rows = []

    for tf in sorted(glob.glob(f"{args.dir}/*.truth_svs.tsv")):
        base = tf.replace(".truth_svs.tsv", "")
        wave = base + ".wave.vcf"
        if not os.path.exists(wave):
            continue
        recs = read_vcf(wave, args.min_sv_len)

        with open(tf) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {c: i for i, c in enumerate(hdr)}
            for line in fh:
                c = line.rstrip("\n").split("\t")
                s, e = int(c[ix["start"]]), int(c[ix["end"]])
                L, ac = int(c[ix["len"]]), int(c[ix["ac"]])
                sbin = c[ix["size_bin"]] if "size_bin" in ix else ""
                slack = max(200, int(0.5 * L))
                over = [r for r in recs if r[1] > s - slack and r[0] < e + slack]

                has_true = any(ac in r[3] for r in over)
                size_match = [r for r in over
                              if abs(r[2] - L) <= args.size_tol * L]

                if has_true:
                    n_alleles = sum(len(r[3]) for r in over)
                    mode = "EXACT_SINGLE" if (len(over) == 1 and n_alleles <= 1) \
                        else "EXACT_MULTI"
                elif not over:
                    mode = "ABSENT"
                elif size_match:
                    mode = "MISCOUNTED"
                elif max(r[2] for r in over) < 0.5 * L:
                    mode = "UNDERSIZED"
                else:
                    mode = "SPLIT"

                modes[sbin][mode] += 1
                modes["ALL"][mode] += 1
                rows.append((os.path.basename(base), c[ix["sv_id"]],
                             c[ix["type"]], L, ac, sbin, mode,
                             len(over),
                             max([r[2] for r in over], default=0),
                             ",".join(str(x) for r in over for x in r[3])[:60]))
                if mode in ("ABSENT", "MISCOUNTED") and len(examples[mode]) < 8:
                    examples[mode].append(
                        (os.path.basename(base), c[ix["sv_id"]], c[ix["type"]],
                         L, ac, len(over),
                         max([r[2] for r in over], default=0),
                         [r[3] for r in over][:3]))

    order = ["EXACT_SINGLE", "EXACT_MULTI", "MISCOUNTED", "SPLIT",
             "UNDERSIZED", "ABSENT"]
    bins = [b for b in modes if b != "ALL"]
    bins.sort(key=lambda x: int(x) if str(x).isdigit() else 0)

    hdr = f"{'bin':>7} {'n':>5} " + " ".join(f"{m:>13}" for m in order)
    print(hdr)
    for b in bins + ["ALL"]:
        tot = sum(modes[b].values())
        cells = " ".join(f"{modes[b][m]/tot:>13.3f}" for m in order)
        print(f"{b:>7} {tot:>5} {cells}")

    print("\ncounts")
    print(hdr)
    for b in bins + ["ALL"]:
        tot = sum(modes[b].values())
        cells = " ".join(f"{modes[b][m]:>13d}" for m in order)
        print(f"{b:>7} {tot:>5} {cells}")

    for m in ("MISCOUNTED", "ABSENT"):
        if examples[m]:
            print(f"\n=== {m} examples ===")
            for e in examples[m]:
                print(f"  {e[0]} {e[1]} {e[2]} len={e[3]} AC_true={e[4]}  "
                      f"large_records={e[5]} max_size={e[6]}  acs={e[7]}")

    if args.out:
        with open(f"{args.out}.tsv", "w") as fh:
            fh.write("replicate\tsv_id\ttype\tlen\tac_true\tsize_bin\tmode\t"
                     "n_large_records\tmax_record_size\trecord_acs\n")
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        print(f"\nwrote {args.out}.tsv")


if __name__ == "__main__":
    main()
