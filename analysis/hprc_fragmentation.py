#!/usr/bin/env python3
"""
hprc_fragmentation.py

Measures, across a chromosome of a released pangenome VCF, how often one
structural variant is spread over several alternate alleles whose sequences
are nearly identical. In simulation this splitting is caused by point
variation inside inserted sequence and disappears when point mutations are
removed; this script asks how common the same pattern is in real data.

Alleles are joined with the rule of vcf_partition_merge.py (same kind,
similar size, nearby, disjoint haplotype carriers, shared 21-mers), applied
in a single streaming pass so a whole chromosome fits in memory.

For each joined event it records the number of alleles and records, the
summed carrier count, and the carrier count of the largest single allele.
An event is fragmented when it has two or more alleles. The share of
carriers outside the largest allele measures how much of the event's count
a per allele reading of the VCF misses.

Run two settings and report both:
  strict   pieces of equal length up to max(10 bp, 1%)  --size-min 10 --size-frac 0.01
  default  the simulation rule, max(50 bp, 20%)

Usage, a test region first:
    bcftools view -r chr20:25700000-25900000 RELEASE.wave.vcf.gz \\
      | python3 analysis/hprc_fragmentation.py --vcf - --out out/hprc/frag_test

Outputs <out>.events.tsv and <out>.summary.txt
"""

import argparse
import gzip
import os
import random
import re
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from vcf_partition_merge import identity, revcomp, kmers, kmer_sim  # noqa: E402

BINS = [(50, 100), (100, 300), (300, 1000), (1000, 3000), (3000, 10000),
        (10000, 30000), (30000, 10**12)]


def bin_label(size):
    for lo, hi in BINS:
        if lo <= size < hi:
            return f"{lo}-{hi}" if hi < 10**12 else f">={lo}"
    return "<50"


def alleles_from(fh, min_size, k, subsample=0, seed=1, info=None):
    """Yield callset alleles in file order, skipping records with no
    structural allele before parsing genotypes. With subsample > 0, only a
    fixed random set of that many haplotypes is counted, chosen once from
    the genotype layout of the first structural record."""
    rec = -1
    chosen = None
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t", 9)
        if len(f) < 10:
            continue
        rec += 1
        ref = f[3]
        alts = f[4].split(",")
        keep = []
        for j, alt in enumerate(alts, start=1):
            if alt.startswith("<") or alt == "*":
                continue
            d = len(alt) - len(ref)
            if d <= -min_size:
                keep.append((j, "DEL", -d, ref))
            elif d >= min_size:
                keep.append((j, "INS", d, alt))
            elif len(ref) >= min_size and identity(alt, revcomp(ref)) >= 0.90:
                keep.append((j, "INV", len(ref), ref))
        if not keep:
            continue
        cols = f[9].split("\t")
        if subsample and chosen is None:
            allh = [(si, hi) for si, g in enumerate(cols)
                    for hi in range(len(g.split(":", 1)[0].replace("/", "|").split("|")))]
            if subsample > len(allh):
                sys.exit(f"--subsample {subsample} exceeds the {len(allh)} haplotypes")
            chosen = set(random.Random(seed).sample(allh, subsample))
            if info is not None:
                info["total_haplotypes"] = len(allh)
        carriers, an = {}, 0
        for si, g in enumerate(cols):
            gt = g.split(":", 1)[0].replace("/", "|").split("|")
            for hi, a in enumerate(gt):
                if a == ".":
                    continue
                if chosen is not None and (si, hi) not in chosen:
                    continue
                an += 1
                if a != "0":
                    carriers.setdefault(int(a), set()).add((si, hi))
        pos = int(f[1])
        # graph site: the snarl ID with any vcfwave _N suffix removed
        site = re.sub(r"_\d+$", "", f[2]) if f[2] != "." else f"{f[0]}:{pos}"
        for j, kind, size, sq in keep:
            c = carriers.get(j)
            if not c:
                continue
            yield {"chrom": f[0], "pos": pos, "end": pos + len(ref), "rec": rec,
                   "kind": kind, "size": size, "an": an, "carriers": c,
                   "site": site,
                   "kmers": kmers(sq, k)}


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--vcf", required=True, help="VCF path, .gz, or - for stdin")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-size", type=int, default=50)
    ap.add_argument("--size-frac", type=float, default=0.20)
    ap.add_argument("--size-min", type=int, default=50)
    ap.add_argument("--pos-frac", type=float, default=1.0)
    ap.add_argument("--pos-min", type=int, default=200)
    ap.add_argument("--min-kmer", type=float, default=0.8)
    ap.add_argument("--k", type=int, default=21)
    ap.add_argument("--subsample", type=int, default=0,
                    help="count only this many randomly chosen haplotypes (0 = all)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--same-site-only", action="store_true",
                    help="join only alleles from the same graph site (snarl ID); "
                         "gives a conservative lower bound")
    ap.add_argument("--lookback", type=int, default=200,
                    help="open groups compared against each new allele")
    a_ = ap.parse_args()

    if a_.vcf == "-":
        fh = sys.stdin
    elif a_.vcf.endswith(".gz"):
        fh = gzip.open(a_.vcf, "rt")
    else:
        fh = open(a_.vcf)

    os.makedirs(os.path.dirname(a_.out) or ".", exist_ok=True)
    ev = open(a_.out + ".events.tsv", "w")
    cols = ["chrom", "pos", "end", "kind", "size", "size_bin", "n_alleles",
            "n_records", "n_sites", "ac", "largest_ac", "max_an"]
    ev.write("\t".join(cols) + "\n")
    events = []

    def close(g):
        m = g["members"]
        row = {"chrom": g["chrom"], "pos": g["pos"], "end": g["end"],
               "kind": g["kind"], "size": g["anchor_size"],
               "size_bin": bin_label(g["anchor_size"]),
               "n_alleles": len(m), "n_records": len(g["recs"]),
               "n_sites": len(g["sites"]),
               "ac": len(g["carriers"]), "largest_ac": max(len(x) for x in m),
               "max_an": g["max_an"]}
        ev.write("\t".join(str(row[c]) for c in cols) + "\n")
        events.append({k: row[k] for k in ("kind", "size_bin", "n_alleles",
                                           "n_records", "n_sites", "ac", "largest_ac")})

    open_groups = []
    n_alleles = 0
    info = {}
    for x in alleles_from(fh, a_.min_size, a_.k, a_.subsample, a_.seed, info):
        n_alleles += 1
        # retire groups that no later allele can reach
        still = []
        for g in open_groups:
            reach = g["end"] + max(a_.pos_min, a_.pos_frac * g["anchor_size"])
            if g["chrom"] != x["chrom"] or x["pos"] > reach:
                close(g)
            else:
                still.append(g)
        open_groups = still

        best, best_d = None, None
        for g in open_groups[-a_.lookback:]:
            if g["kind"] != x["kind"]:
                continue
            s0 = g["anchor_size"]
            if abs(x["size"] - s0) > max(a_.size_min, a_.size_frac * max(x["size"], s0)):
                continue
            ptol = max(a_.pos_min, a_.pos_frac * s0)
            if x["pos"] > g["end"] + ptol or x["end"] < g["pos"] - ptol:
                continue
            if x["carriers"] & g["carriers"]:
                continue
            if a_.same_site_only and x["site"] not in g["sites"]:
                continue
            if kmer_sim(x["kmers"], g["anchor_kmers"]) < a_.min_kmer:
                continue
            d = abs(x["size"] - s0)
            if best is None or d < best_d:
                best, best_d = g, d
        if best is None:
            open_groups.append({
                "chrom": x["chrom"], "kind": x["kind"], "pos": x["pos"],
                "end": x["end"], "anchor_size": x["size"],
                "anchor_kmers": x["kmers"], "carriers": set(x["carriers"]),
                "members": [x["carriers"]], "recs": {x["rec"]},
                "sites": {x["site"]},
                "max_an": x["an"]})
        else:
            best["pos"] = min(best["pos"], x["pos"])
            best["end"] = max(best["end"], x["end"])
            best["carriers"] |= x["carriers"]
            best["members"].append(x["carriers"])
            best["recs"].add(x["rec"])
            best["sites"].add(x["site"])
            best["max_an"] = max(best["max_an"], x["an"])
        if n_alleles % 20000 == 0:
            sys.stderr.write(f"  {n_alleles} alleles, at {x['chrom']}:{x['pos']}\n")
    for g in open_groups:
        close(g)
    ev.close()

    L = [(f"haplotypes counted: {a_.subsample} of {info.get('total_haplotypes')} "
          f"(seed {a_.seed})" if a_.subsample else "haplotypes counted: all"),
         f"alleles of at least {a_.min_size} bp with a carrier: {n_alleles}",
         f"events after joining: {len(events)}",
         f"rule: size within max({a_.size_min} bp, {a_.size_frac:.0%}); start within "
         f"max({a_.pos_min} bp, {a_.pos_frac} x size); disjoint haplotypes; "
         f"shared {a_.k}-mers >= {a_.min_kmer}; same site only: {a_.same_site_only}",
         "",
         "kind  size        events  fragmented  frac   alleles med  "
         "multi record  of which one site  carriers outside largest allele",
         ]
    order = [bin_label(lo) for lo, _ in BINS]
    for kind in ("DEL", "INS", "INV"):
        for b in order:
            s = [e for e in events if e["kind"] == kind and e["size_bin"] == b]
            if not s:
                continue
            fr = [e for e in s if e["n_alleles"] >= 2]
            tot = sum(e["ac"] for e in s)
            outside = sum(e["ac"] - e["largest_ac"] for e in s)
            L.append(f"{kind:4}  {b:>11} {len(s):7d} {len(fr):10d}  {len(fr)/len(s):.3f}  "
                     f"{(st.median(e['n_alleles'] for e in fr) if fr else 0):11.1f}  "
                     f"{sum(e['n_records'] >= 2 for e in fr):12d}  "
                     f"{sum(e['n_records'] >= 2 and e['n_sites'] == 1 for e in fr):15d}  "
                     f"{(outside/tot if tot else 0):.3f}")
    L += ["",
          "fragmented        events joined from two or more alleles",
          "alleles med       median alleles per fragmented event",
          "multi record      fragmented events whose alleles span two or more records",
          "one site          of those, events whose records all come from one graph",
          "                  site (same snarl ID); the rule never reads the ID",
          "carriers outside  share of all carriers in the class that are not in their",
          "                  event's largest allele, i.e. missed by a per allele count",
          "",
          "No truth set exists for real data, so these are joins by the rule, not",
          "verified events. In simulation the rule made no false joins at 21-mer",
          "thresholds of 0.5 to 0.9."]
    rep = "\n".join(L)
    open(a_.out + ".summary.txt", "w").write(rep + "\n")
    print(rep)
    print(f"\nwrote {a_.out}.events.tsv and {a_.out}.summary.txt")


if __name__ == "__main__":
    main()
