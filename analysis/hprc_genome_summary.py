#!/usr/bin/env python3
"""
hprc_genome_summary.py

Combines the per chromosome event tables written by hprc_fragmentation.py.

Reports, for deletions and insertions in each size class:
  - the genome wide share of carriers outside each event's largest allele
  - the same share per chromosome, with its range across chromosomes
and writes a per megabase table for a genome wide plot.

Usage:
    python3 analysis/hprc_genome_summary.py --dir out/hprc/genome \\
        --out out/hprc/genome_summary

Outputs <out>.txt, <out>.by_chrom.tsv and <out>.per_mb.tsv
"""

import argparse
import csv
import glob
import os
import re

ORDER = ["50-100", "100-300", "300-1000", "1000-3000", "3000-10000",
         "10000-30000", ">=30000"]


def chrom_key(c):
    m = re.match(r"chr(\d+)$", c)
    return int(m.group(1)) if m else 99


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mb", type=int, default=1_000_000, help="bin width for the per Mb table")
    a_ = ap.parse_args()

    files = sorted(glob.glob(os.path.join(a_.dir, "*.events.tsv")),
                   key=lambda p: chrom_key(os.path.basename(p).split(".")[0]))
    if not files:
        raise SystemExit(f"no *.events.tsv in {a_.dir}")

    agg = {}      # (chrom, kind, bin) -> [events, fragmented, ac, outside]
    per_mb = {}   # (chrom, mb, kind) -> [ac, outside, events]
    for p in files:
        with open(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["kind"] not in ("DEL", "INS"):
                    continue
                ac, big, na = int(r["ac"]), int(r["largest_ac"]), int(r["n_alleles"])
                key = (r["chrom"], r["kind"], r["size_bin"])
                a = agg.setdefault(key, [0, 0, 0, 0])
                a[0] += 1; a[1] += na >= 2; a[2] += ac; a[3] += ac - big
                mk = (r["chrom"], int(r["pos"]) // a_.mb, r["kind"])
                m = per_mb.setdefault(mk, [0, 0, 0])
                m[0] += ac; m[1] += ac - big; m[2] += 1

    chroms = sorted({k[0] for k in agg}, key=chrom_key)
    L = [f"chromosomes: {len(chroms)}  ({', '.join(chroms)})", "",
         "Genome wide share of carriers outside each event's largest allele,",
         "with the range of the per chromosome values (classes with at least",
         "20 events on that chromosome).", "",
         f"{'kind':4} {'size':>11} {'events':>8} {'fragmented':>10} "
         f"{'carriers missed':>15}   {'per chromosome min to max':>25}"]
    with open(a_.out + ".by_chrom.tsv", "w") as fh:
        fh.write("chrom\tkind\tsize_bin\tevents\tfragmented\tac\toutside\tshare_missed\n")
        for c in chroms:
            for kind in ("DEL", "INS"):
                for b in ORDER:
                    v = agg.get((c, kind, b))
                    if v:
                        fh.write(f"{c}\t{kind}\t{b}\t{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t"
                                 f"{(v[3]/v[2] if v[2] else 0):.4f}\n")
    for kind in ("DEL", "INS"):
        for b in ORDER:
            vs = [agg[(c, kind, b)] for c in chroms if (c, kind, b) in agg]
            if not vs:
                continue
            ev = sum(v[0] for v in vs); fr = sum(v[1] for v in vs)
            ac = sum(v[2] for v in vs); out = sum(v[3] for v in vs)
            per = [v[3] / v[2] for v in vs if v[0] >= 20 and v[2]]
            rng = f"{min(per):.3f} to {max(per):.3f} (n={len(per)})" if per else "too few events"
            L.append(f"{kind:4} {b:>11} {ev:8d} {fr:10d} {(out/ac if ac else 0):15.3f}   {rng:>25}")
    with open(a_.out + ".per_mb.tsv", "w") as fh:
        fh.write("chrom\tmb\tkind\tevents\tac\toutside\tshare_missed\n")
        for (c, mb, kind) in sorted(per_mb, key=lambda k: (chrom_key(k[0]), k[1], k[2])):
            ac, out, ev = per_mb[(c, mb, kind)]
            fh.write(f"{c}\t{mb}\t{kind}\t{ev}\t{ac}\t{out}\t{(out/ac if ac else 0):.4f}\n")
    rep = "\n".join(L)
    open(a_.out + ".txt", "w").write(rep + "\n")
    print(rep)
    print(f"\nwrote {a_.out}.txt, .by_chrom.tsv, .per_mb.tsv")


if __name__ == "__main__":
    main()
