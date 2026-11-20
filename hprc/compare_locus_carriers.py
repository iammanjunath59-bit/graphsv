#!/usr/bin/env python3
"""
compare_locus_carriers.py

Compares, haplotype by haplotype, the carriers of a large deletion as found
by joining VCF records (no graph used) with the carriers implied by graph
dosage on an independently built graph (no VCF used).

VCF side: every deletion allele of at least --min-del bp between --from and
--to is collected; a haplotype carries the event if it carries any of them.
This is the same set the record joining tool produced.

Graph side: gfa_locus_dosage.py output; a complete path with dosage below
--cut is a carrier. Haplotypes with several complete paths, or none, are
reported separately and left out of the comparison.

Haplotypes are matched by sample name and haplotype number. In the VCF the
first genotype field is taken as haplotype 1. Because phase labels need not
agree between builders, agreement is also reported per sample as dosage
0, 1 or 2, which does not depend on phase.

Usage:
    python3 compare_locus_carriers.py --vcf out/hprc/locus.wave.vcf \\
        --dosage out/hprc/graph/del62k.paths.tsv --cut 0.5
"""

import argparse
import csv
from collections import Counter, defaultdict


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--dosage", required=True)
    ap.add_argument("--cut", type=float, required=True,
                    help="graph dosage below this marks a carrier; choose from the histogram")
    ap.add_argument("--min-del", type=int, default=50000)
    ap.add_argument("--from", dest="frm", type=int, default=25700000)
    ap.add_argument("--to", type=int, default=25900000)
    a_ = ap.parse_args()

    samples, vcar, called = [], set(), set()
    recs = 0
    with open(a_.vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            pos = int(f[1])
            if not (a_.frm <= pos <= a_.to):
                continue
            ref = f[3]
            want = {j for j, alt in enumerate(f[4].split(","), start=1)
                    if len(ref) - len(alt) >= a_.min_del}
            if not want:
                continue
            recs += 1
            for si, g in enumerate(f[9:]):
                for hi, x in enumerate(g.split(":")[0].replace("/", "|").split("|")):
                    if x == ".":
                        continue
                    called.add((samples[si], str(hi + 1)))
                    if x.isdigit() and int(x) in want:
                        vcar.add((samples[si], str(hi + 1)))

    gpaths = defaultdict(list)
    with open(a_.dosage) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["complete"] == "1":
                gpaths[(r["sample"], r["hap"])].append(float(r["dosage"]))
    multi = {k for k, v in gpaths.items() if len(v) > 1}
    gcar = {k for k, v in gpaths.items() if len(v) == 1 and v[0] < a_.cut}
    gscored = {k for k, v in gpaths.items() if len(v) == 1}

    both = called & gscored
    tab = Counter((k in vcar, k in gcar) for k in both)
    print(f"VCF: {recs} deletion records of at least {a_.min_del} bp; carriers {len(vcar)} "
          f"of {len(called)} called haplotypes")
    print(f"graph: {len(gscored)} haplotypes with one complete path, carriers {len(gcar)} "
          f"(dosage < {a_.cut}); {len(multi)} with several complete paths left out")
    print(f"\nhaplotypes scored on both sides: {len(both)}")
    print(f"  carrier in both            {tab[(True, True)]}")
    print(f"  carrier in VCF only        {tab[(True, False)]}")
    print(f"  carrier in graph only      {tab[(False, True)]}")
    print(f"  carrier in neither         {tab[(False, False)]}")
    if both:
        print(f"  haplotype agreement        {(tab[(True, True)] + tab[(False, False)]) / len(both):.4f}")

    sv, sg = Counter(), Counter()
    sboth = defaultdict(int)
    for k in both:
        sboth[k[0]] += 1
        sv[k[0]] += k in vcar
        sg[k[0]] += k in gcar
    full = [s for s, n in sboth.items() if n == 2]
    agree = sum(1 for s in full if sv[s] == sg[s])
    print(f"\nsamples with both haplotypes scored on both sides: {len(full)}")
    if full:
        print(f"  dosage agreement (0, 1 or 2 copies, phase free) {agree}/{len(full)} = {agree/len(full):.4f}")
    dis = sorted(k for k in both if (k in vcar) != (k in gcar))
    if dis:
        print("\ndisagreeing haplotypes (sample, hap, VCF carrier, graph dosage):")
        for k in dis[:40]:
            print(f"  {k[0]}#{k[1]}  {int(k in vcar)}  {gpaths[k][0]:.3f}")


if __name__ == "__main__":
    main()
