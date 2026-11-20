#!/usr/bin/env python3
"""
pair_haps.py

Converts a haploid multi-sample VCF to phased diploid samples and back, so
that tools expecting diploid genotypes (truvari collapse treats a genotype
with other than two fields as unknown) can be run on the simulated callsets.

  to-diploid: columns 0,1 -> sample "A__B" with GT a|b; 2,3 -> next; ...
  to-haploid: each "A__B" sample split back into columns A and B

Only GT is kept. Usage:
  python3 pair_haps.py to-diploid in.vcf > out.vcf
  python3 pair_haps.py to-haploid in.vcf > out.vcf
"""

import gzip
import sys


def main():
    mode, path = sys.argv[1], sys.argv[2]
    op = gzip.open if path.endswith(".gz") else open
    out = sys.stdout
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                if not line.startswith("##FORMAT=") or "ID=GT," in line:
                    out.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                s = f[9:]
                if mode == "to-diploid":
                    if len(s) % 2:
                        sys.exit("odd number of haploid samples")
                    names = [f"{s[i]}__{s[i + 1]}" for i in range(0, len(s), 2)]
                else:
                    names = [p for n in s for p in n.split("__")]
                out.write("\t".join(f[:9] + names) + "\n")
                continue
            gts = [g.split(":")[0] for g in f[9:]]
            if mode == "to-diploid":
                new = [f"{gts[i]}|{gts[i + 1]}" for i in range(0, len(gts), 2)]
            else:
                new = [a for g in gts for a in g.replace("/", "|").split("|")]
            out.write("\t".join(f[:8] + ["GT"] + new) + "\n")


if __name__ == "__main__":
    main()
