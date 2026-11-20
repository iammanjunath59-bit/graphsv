#!/usr/bin/env python3
"""
supp_tables.py

Prints the values for Additional file 1, Tables S1 and S2, from the files the
paper's numbers come from. Nothing is computed that the main text does not
already report; the tables only break those numbers down by size class.

Table S1  Recovery at each stage of the conversion to VCF, by size class.
          Input: results/stages/*.<stage>.cmp.svs.tsv for stage raw, bub and,
          if present, wave. Columns are read by position, exactly as
          verify_all.sh reads them: 4 size bin, 8 number of alleles at the
          matched site, 14 true count present (1/0).
            count present   column 14 == 1
            at most two     column 14 == 1 and column 8 <= 2

Table S2  Duplications split across alleles, by size class (all 108 chr20
          replicates). Input: out/pieces.tsv from piece_diversity.py.
            split           exact_single == 0 and n_pieces >= 2
          For split duplications: median and range of n_pieces, number with
          same_length == 1, number with pieces_cover_truth == 1, median and
          range of min_subs (substitutions between the two most similar
          pieces), and number at a single position (n_positions == 1).

Usage (from ~/graphsv):
    python3 analysis/supp_tables.py --stages results/stages --pieces out/pieces.tsv
"""

import argparse
import csv
import glob
import os
import statistics as st
from collections import defaultdict

BINS = ["100", "300", "1000", "3000", "10000", "30000"]
LAB = ["100 bp", "300 bp", "1 kb", "3 kb", "10 kb", "30 kb"]


def stage_table(d):
    print("TABLE S1  recovery by stage and size class")
    print(f"  {'stage':6} {'files':>5}  " + "  ".join(f"{l:>13}" for l in LAB) + "  overall")
    for stage in ("raw", "bub", "wave"):
        files = sorted(glob.glob(os.path.join(d, f"*.{stage}.cmp.svs.tsv")))
        if not files:
            print(f"  {stage:6} no files")
            continue
        n = defaultdict(int); p = defaultdict(int); c = defaultdict(int)
        for f in files:
            with open(f) as fh:
                for line in fh:
                    if line.startswith("sv_id"):
                        continue
                    x = line.rstrip("\n").split("\t")
                    b = x[3]
                    n[b] += 1
                    if x[13] == "1":
                        p[b] += 1
                        if int(x[7]) <= 2:
                            c[b] += 1
        N = sum(n.values())
        for name, k in (("present", p), ("at most two", c)):
            cells = "  ".join(f"{k[b]:4d}/{n[b]:<4d}{k[b]/n[b]:.3f}" if n[b] else f"{'none':>13}" for b in BINS)
            print(f"  {stage:6} {len(files):5d}  {cells}  {sum(k.values())}/{N} {sum(k.values())/N:.4f}  [{name}]")


def pieces_table(path):
    rows = defaultdict(list); reps = set()
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["type"] == "DUP":
                rows[r["size_bin"]].append(r); reps.add(r["replicate"])
    print(f"\nTABLE S2  duplications split across alleles ({len(reps)} replicates)")
    print("  size     n  split  pieces median (range)  same length  cover truth  one position  min subs median (range)")
    tot = 0
    for b, l in zip(BINS, LAB):
        rs = rows.get(b, [])
        sp = [r for r in rs if r["exact_single"] == "0" and int(r["n_pieces"]) >= 2]
        tot += len(sp)
        if not sp:
            print(f"  {l:>6} {len(rs):4d}  {0:5d}"); continue
        npc = [int(r["n_pieces"]) for r in sp]
        ms = [float(r["min_subs"]) for r in sp if r["min_subs"] not in ("", "NA", "nan")]
        print(f"  {l:>6} {len(rs):4d}  {len(sp):5d}  {st.median(npc):6g} ({min(npc)} to {max(npc)})"
              f"  {sum(r['same_length'] == '1' for r in sp):11d}  {sum(r['pieces_cover_truth'] == '1' for r in sp):11d}"
              f"  {sum(r['n_positions'] == '1' for r in sp):12d}"
              f"  {(f'{st.median(ms):g} ({min(ms):g} to {max(ms):g})') if ms else 'none'}  [min_subs values {len(ms)}]")
    print(f"  total split {tot}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages", default="results/stages")
    ap.add_argument("--pieces", default="out/pieces.tsv")
    a = ap.parse_args()
    stage_table(a.stages)
    pieces_table(a.pieces)


if __name__ == "__main__":
    main()
