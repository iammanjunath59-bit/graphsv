#!/usr/bin/env python3
"""
pantree_dup_carriers.py

pantree writes a tandem duplication as a loop edge (INFO VT=DUP). Its GT
field records only whether a haplotype visits the alternate edge at all, and
in a loop every haplotype that passes the duplicated segment can visit it
once; the copy number is in the CA field ("Number of times visiting the ALT
allele"). A duplication carrier visits the edge at least twice.

For records with VT=DUP this script sets each called haplotype's GT to 1 if
CA >= 2, otherwise 0. Other records are copied unchanged. A tally of the CA
values seen in VT=DUP records is printed so the rule can be checked.

Usage:
    python3 pantree_dup_carriers.py in.vcf > out.vcf
"""

import sys
from collections import Counter


def main():
    path = sys.argv[1]
    out = sys.stdout
    n_dup = 0
    ca_tally = Counter()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
            if info.get("VT") != "DUP" or "," in f[4]:
                out.write(line)
                continue
            keys = f[8].split(":")
            if "CA" not in keys or "GT" not in keys:
                out.write(line)
                continue
            ig, ic = keys.index("GT"), keys.index("CA")
            n_dup += 1
            new = []
            for g in f[9:]:
                p = g.split(":")
                if p[ig] == ".":
                    new.append(g)
                    continue
                ca = int(p[ic]) if p[ic] not in (".", "") else 0
                ca_tally[ca] += 1
                p[ig] = "1" if ca >= 2 else "0"
                new.append(":".join(p))
            out.write("\t".join(f[:9] + new) + "\n")
    sys.stderr.write(f"VT=DUP records {n_dup}  CA values over haplotypes: "
                     + " ".join(f"{k}:{v}" for k, v in sorted(ca_tally.items())) + "\n")


if __name__ == "__main__":
    main()
