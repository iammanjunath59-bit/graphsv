#!/usr/bin/env python3
"""
repolarize_vcf.py

Returns records to reference polarity before scoring. A tool whose
reference tree leaves the reference haplotype (pantree, at sites where the
reference path revisits a node) can report a deletion as an insertion of
the same sequence carried by everyone except the deletion carriers, with the
reference sample itself carrying the "alternate" allele.

Rule: in a biallelic record where the reference sample's genotype is the
alternate allele, swap REF and ALT and invert every called genotype (0 <-> 1).
Missing genotypes stay missing; other FORMAT fields are dropped (only GT is
kept). Multiallelic records and records where the reference sample is not
called are written unchanged and counted.

Usage:
    python3 repolarize_vcf.py --ref-sample anc in.vcf > out.vcf
"""

import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("vcf")
    ap.add_argument("--ref-sample", required=True)
    a_ = ap.parse_args()
    out = sys.stdout
    ri = None
    n = flipped = multi = uncalled = 0
    with open(a_.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                if a_.ref_sample not in f[9:]:
                    sys.exit(f"sample {a_.ref_sample} not in header")
                ri = 9 + f[9:].index(a_.ref_sample)
                out.write(line)
                continue
            n += 1
            if "," in f[4]:
                multi += 1
                out.write(line)
                continue
            rg = f[ri].split(":")[0]
            if rg in (".", "./.", ".|."):
                uncalled += 1
                out.write(line)
                continue
            if "1" not in rg.replace("/", "|").split("|"):
                out.write(line)
                continue
            flipped += 1
            sw = {"0": "1", "1": "0"}
            gts = []
            for g in f[9:]:
                gt = g.split(":")[0]
                sep = "|" if "|" in gt else ("/" if "/" in gt else "")
                parts = gt.split(sep) if sep else [gt]
                gts.append(sep.join(sw.get(p, p) for p in parts))
            f[3], f[4] = f[4], f[3]
            f[8] = "GT"
            out.write("\t".join(f[:9] + gts) + "\n")
    sys.stderr.write(f"records {n}  repolarized {flipped}  multiallelic left {multi}  "
                     f"reference uncalled {uncalled}\n")


if __name__ == "__main__":
    main()
