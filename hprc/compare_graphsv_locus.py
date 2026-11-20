#!/usr/bin/env python3
"""
compare_graphsv_locus.py

Compares three counts of one large deletion, haplotype by haplotype:
graphsv run on the extracted locus (gfa_extract_locus.py, then graphsv
--carriers), the graph dosage of the deleted interval (gfa_locus_dosage.py),
and the VCF deletion records at the locus.

graphsv carriers are those of the ALT allele, among calls overlapping the
event, with the largest length deficit against the reference traversal (as
in Figure 4c), provided it is at least --min-graph-del bp. graphsv lengths
count only segments of at least its --min-seg, so a 62 kb deletion can show
a smaller deficit. Every graphsv call overlapping the event is printed, so
the choice can be checked. Haplotypes are matched on
(sample, haplotype) from PanSN names; VCF haplotypes are 1 and 2 in GT order.

Usage:
    python3 hprc/compare_graphsv_locus.py --calls out/hprc/graph/del62k_locus.gsv.tsv \\
        --locus-gfa out/hprc/graph/del62k_locus.gfa \\
        --offset out/hprc/graph/del62k_locus.offset.txt \\
        --dosage out/hprc/graph/del62k.paths.tsv --vcf out/hprc/locus.wave.vcf \\
        --event 25783207-25845155
"""

import argparse
import csv


def vcf_carriers(path, lo, hi, min_del):
    samples, car, called, recs = [], set(), set(), []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            pos = int(f[1])
            if not lo <= pos <= hi:
                continue
            want = {j for j, alt in enumerate(f[4].split(","), 1) if len(f[3]) - len(alt) >= min_del}
            if not want:
                continue
            here = set()
            for si, g in enumerate(f[9:]):
                for h, x in enumerate(g.split(":")[0].replace("/", "|").split("|")):
                    if x == ".":
                        continue
                    called.add((samples[si], str(h + 1)))
                    if x.isdigit() and int(x) in want:
                        here.add((samples[si], str(h + 1)))
            recs.append((pos, len(here)))
            car |= here
    return car, called, recs


def hap(name):
    p = name.split("#")
    return (p[0], p[1] if len(p) > 2 else "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--locus-gfa", required=True, help="the GFA graphsv was run on")
    ap.add_argument("--offset", required=True)
    ap.add_argument("--dosage", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--event", required=True)
    ap.add_argument("--cut", type=float, default=0.8)
    ap.add_argument("--min-del", type=int, default=50000, help="VCF deletion alleles of at least this size")
    ap.add_argument("--min-graph-del", type=int, default=10000,
                    help="the graphsv allele taken is the one with the largest length deficit against the "
                         "reference traversal, if at least this; graphsv lengths count only segments >= --min-seg")
    ap.add_argument("--from", dest="frm", type=int, default=25700000)
    ap.add_argument("--to", type=int, default=25900000)
    a = ap.parse_args()
    s, e = (int(x) for x in a.event.split("-"))
    off = int(open(a.offset).read().split()[0])

    gsv, haps_in_graph, best = set(), set(), None
    with open(a.locus_gfa) as fh:
        for line in fh:
            if line[0] == "P":
                n = line.split("\t", 2)[1]
                if not n.startswith("GRCh38"):
                    haps_in_graph.add(hap(n))
    print("graphsv calls overlapping the event (reference coordinates):")
    with open(a.calls) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if rows and "allele_carriers" not in rows[0]:
        raise SystemExit("calls were written without --carriers")
    for r in rows:
        p, q, R = int(r["pos"]) + off, int(r["end"]) + off, int(r["ref_len"])
        lens = [int(x) for x in r["allele_lens"].split(",") if x]
        groups = r["allele_carriers"].split(",")
        if not (q > s and p < e):
            continue
        print(f"  {p:,}-{q:,}  ref_len {R:,}  alleles {r['allele_lens']}  counts {r['allele_acs']}")
        for L, g in zip(lens, groups):
            if best is None or R - L > best[0]:
                best = (R - L, L, R, g)
    if best is not None and best[0] >= a.min_graph_del:
        gsv = {hap(x) for x in best[3].split(";") if x}
        print(f"  taken: allele of {best[1]:,} bp against a reference traversal of {best[2]:,} bp "
              f"({best[0]:,} bp shorter; lengths count segments >= graphsv --min-seg only)")
    else:
        print(f"  no allele at least {a.min_graph_del:,} bp shorter than the reference traversal")

    dos = {}
    with open(a.dosage) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["complete"] == "1":
                dos[(r["sample"], r["hap"])] = float(r["dosage"])
    dcar = {k for k, v in dos.items() if v < a.cut}
    vcar, called, recs = vcf_carriers(a.vcf, a.frm, a.to, a.min_del)

    print(f"\ngraphsv carriers {len(gsv)}   graph dosage < {a.cut}: {len(dcar)}   "
          f"VCF joined records: {len(vcar)} of {len(called)} called (records {sorted(recs)})")
    both = haps_in_graph & set(dos)
    print(f"haplotypes in the graphsv locus graph {len(haps_in_graph)}; also scored by dosage {len(both)}")
    g, d = gsv & both, dcar & both
    print(f"graphsv vs dosage, on {len(both)} shared haplotypes: both {len(g & d)}, "
          f"graphsv only {len(g - d)}, dosage only {len(d - g)}")
    sc = both & called
    print(f"graphsv vs VCF, on {len(sc)} haplotypes in both: both {len(gsv & vcar & sc)}, "
          f"graphsv only {len((gsv - vcar) & sc)}, VCF only {len((vcar - gsv) & sc)}")


if __name__ == "__main__":
    main()
