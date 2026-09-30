#!/usr/bin/env python3
"""
hprc_dosage_correlation.py

Are overlapping large VCF records at a locus fragments of one event, or
independent variants?

The question
------------
At several HPRC loci the published callset reports multiple records of
near-identical size within a few kilobases, with allele counts that decay
(for example 200, 12, 3, 1 at chr20:25.78 Mb). We interpret this as one
structural polymorphism split across records. A reviewer will reasonably
object that these could be genuinely distinct variants: a smaller deletion
nested inside a larger inversion, say, or independent events in a
duplication-rich region.

Genotype dosage settles it. If two records are fragments of the same
underlying allelic partition, the per-sample alternate dosages must track
each other closely across all 232 samples, because the same haplotypes carry
both. Independent variants have no such requirement, and in a region with
extensive haplotype diversity would not be expected to correlate strongly.

This script computes pairwise Pearson correlations of per-sample alternate
allele dosage between every pair of large records within a locus, and reports
the distribution.

Interpretation
--------------
High correlation (say above 0.9) across most pairs supports fragmentation:
the records partition the same haplotypes. Low or mixed correlation means
the records carry distinct information and the fragmentation reading does not
hold for that locus. Both outcomes are reportable; do not run this expecting
one of them.

Usage
-----
  python hprc_dosage_correlation.py \\
      --vcf <url or path to hprc wave.vcf.gz> \\
      --loci hprc_loci/loci.txt \\
      --out hprc_loci/dosage
"""

import argparse
import itertools
import subprocess
import sys
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("needs numpy")


def dosages(vcf, region, min_len):
    """Return [(pos, size, ac, np.array of per-sample alt dosage)]."""
    cmd = ["bcftools", "query", "-r", region,
           "-f", "%POS\t%REF\t%ALT\t%INFO/AC[\t%GT]\n", vcf]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True,
                             timeout=600).stdout
    except subprocess.TimeoutExpired:
        return []
    recs = []
    for line in out.rstrip("\n").split("\n"):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        pos = int(f[0])
        ref, alts, ac = f[1], f[2].split(","), f[3]
        size = max([len(ref)] + [len(a) for a in alts])
        if size < min_len:
            continue
        gts = f[4:]
        dos = []
        for g in gts:
            g = g.replace("|", "/")
            if g in (".", "./.", ".|."):
                dos.append(np.nan)
                continue
            alleles = [a for a in g.split("/") if a != "."]
            if not alleles:
                dos.append(np.nan)
                continue
            # dosage = number of non-reference alleles carried
            dos.append(sum(1 for a in alleles if a != "0"))
        recs.append((pos, size, ac, np.array(dos, dtype=float)))
    return recs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--loci", required=True,
                    help="one per line: chrom:start:end:label")
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--out", default="dosage")
    args = ap.parse_args()

    summary = []
    detail_lines = ["locus\tpos_a\tsize_a\tpos_b\tsize_b\tn_shared\tpearson_r"]

    for raw in open(args.loci):
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split(":")
        if len(parts) < 3:
            continue
        chrom, s, e = parts[0], parts[1], parts[2]
        label = parts[3] if len(parts) > 3 else f"{chrom}_{s}"
        region = f"{chrom}:{s}-{e}"

        recs = dosages(args.vcf, region, args.min_len)
        if len(recs) < 2:
            summary.append((label, len(recs), 0, float("nan"), float("nan")))
            continue

        rs = []
        for (pa, sa, aca, da), (pb, sb, acb, db) in itertools.combinations(recs, 2):
            ok = ~(np.isnan(da) | np.isnan(db))
            if ok.sum() < 10:
                continue
            x, y = da[ok], db[ok]
            if x.std() == 0 or y.std() == 0:
                continue
            r = float(np.corrcoef(x, y)[0, 1])
            rs.append(r)
            detail_lines.append(
                f"{label}\t{pa}\t{sa}\t{pb}\t{sb}\t{int(ok.sum())}\t{r:.4f}")
        if rs:
            summary.append((label, len(recs), len(rs),
                            float(np.median(rs)), float(np.max(rs))))
        else:
            summary.append((label, len(recs), 0, float("nan"), float("nan")))

    print(f"{'locus':<24} {'records':>8} {'pairs':>7} {'median r':>10} {'max r':>8}")
    for lab, nrec, npair, med, mx in summary:
        print(f"{lab:<24} {nrec:>8} {npair:>7} "
              f"{med:>10.3f} {mx:>8.3f}" if npair else
              f"{lab:<24} {nrec:>8} {npair:>7} {'--':>10} {'--':>8}")

    with open(f"{args.out}.pairs.tsv", "w") as fh:
        fh.write("\n".join(detail_lines) + "\n")
    print(f"\nwrote {args.out}.pairs.tsv")
    print("\nHigh median r indicates the records partition the same haplotypes")
    print("and are therefore fragments of one event. Low r indicates they")
    print("carry distinct information. Report whichever you observe.")


if __name__ == "__main__":
    main()
