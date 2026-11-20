#!/usr/bin/env python3
"""
vcf_partition_merge.py

Rejoins a structural variant whose carriers a graph to VCF conversion has
split across several records. Uses the VCF alone; the AT field is not needed.

Two callset alleles are joined when all of these hold:

  same kind        DEL with DEL, INS with INS, INV with INV
  similar size     within max(--size-min bp, --size-frac x size)
  nearby           start within max(--pos-min bp, --pos-frac x size) of the
                   group span
  disjoint         no haplotype carries both (a haplotype carries one allele)
  same sequence    the variant sequences share at least --min-kmer of their
                   21-mers (inserted sequence for INS, deleted reference for
                   DEL, reference span for INV); pieces of one tandem event
                   carry the same unit, possibly rotated

Carriers are counted per haplotype (sample, phase index), so diploid and
haploid genotypes are both handled and a homozygous sample counts twice.

Input: a VCF, plain or gzipped, already restricted to the region of interest
(for example with bcftools view -r). Output: one row per merged event.

    bcftools view -r chr20:25700000-25900000 RELEASE.wave.vcf.gz > locus.vcf
    python3 vcf_partition_merge.py --vcf locus.vcf --out locus.merged.tsv
"""

import argparse
import gzip
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def identity(a, b):
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] == b[i]) / n if n else 0.0


def kmers(seq, k=21):
    seq = seq.upper()
    return {hash(seq[i:i + k]) for i in range(len(seq) - k + 1)}


def kmer_sim(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def read_alleles(path, min_size):
    op = gzip.open if path.endswith(".gz") else open
    alleles, an_by_pos = [], {}
    rec = -1
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, ref = f[0], int(f[1]), f[3]
            rec += 1
            carriers, an = {}, 0
            for si, g in enumerate(f[9:]):
                gt = g.split(":")[0].replace("/", "|").split("|")
                for hi, a in enumerate(gt):
                    if a == ".":
                        continue
                    an += 1
                    if a.isdigit() and int(a) > 0:
                        carriers.setdefault(int(a), set()).add((si, hi))
            an_by_pos[(chrom, pos)] = an
            for j, alt in enumerate(f[4].split(","), start=1):
                if alt.startswith("<") or alt == "*":
                    continue
                d = len(alt) - len(ref)
                if d <= -min_size:
                    kind, size, sq = "DEL", -d, ref
                elif d >= min_size:
                    kind, size, sq = "INS", d, alt
                elif len(ref) >= min_size and identity(alt, revcomp(ref)) >= 0.90:
                    kind, size, sq = "INV", len(ref), ref
                else:
                    continue
                if not carriers.get(j):
                    continue  # allele present in the record but carried by no haplotype
                alleles.append({
                    "chrom": chrom, "pos": pos, "end": pos + len(ref),
                    "kind": kind, "size": size, "an": an, "rec": rec,
                    "carriers": carriers.get(j, set()), "kmers": kmers(sq),
                })
    return alleles


def merge(alleles, a_):
    order = sorted(range(len(alleles)),
                   key=lambda i: (alleles[i]["chrom"], alleles[i]["pos"]))
    groups = []
    for i in order:
        x = alleles[i]
        best, best_d = None, None
        for g in groups[-a_.lookback:]:
            if g["chrom"] != x["chrom"] or g["kind"] != x["kind"]:
                continue
            s0 = g["anchor"]["size"]
            if abs(x["size"] - s0) > max(a_.size_min, a_.size_frac * max(x["size"], s0)):
                continue
            ptol = max(a_.pos_min, a_.pos_frac * s0)
            if x["pos"] > g["end"] + ptol or x["end"] < g["pos"] - ptol:
                continue
            if x["carriers"] & g["carriers"]:
                continue
            if a_.across_records_only and x["rec"] in g["recs"]:
                continue
            if kmer_sim(x["kmers"], g["anchor"]["kmers"]) < a_.min_kmer:
                continue
            d = abs(x["size"] - s0)
            if best is None or d < best_d:
                best, best_d = g, d
        if best is None:
            groups.append({"chrom": x["chrom"], "kind": x["kind"], "anchor": x,
                           "pos": x["pos"], "end": x["end"],
                           "carriers": set(x["carriers"]), "members": [x],
                           "recs": {x["rec"]}})
        else:
            best["pos"] = min(best["pos"], x["pos"])
            best["end"] = max(best["end"], x["end"])
            best["carriers"] |= x["carriers"]
            best["members"].append(x)
            best["recs"].add(x["rec"])
    return groups


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-size", type=int, default=50)
    ap.add_argument("--size-frac", type=float, default=0.20)
    ap.add_argument("--size-min", type=int, default=50)
    ap.add_argument("--pos-frac", type=float, default=1.0)
    ap.add_argument("--pos-min", type=int, default=200)
    ap.add_argument("--min-kmer", type=float, default=0.8)
    ap.add_argument("--lookback", type=int, default=50)
    ap.add_argument("--across-records-only", action="store_true",
                    help="never join two alternates of the same VCF record; "
                         "those are already reported as distinct alleles of one site")
    a_ = ap.parse_args()

    alleles = read_alleles(a_.vcf, a_.min_size)
    groups = merge(alleles, a_)
    cols = ["chrom", "pos", "end", "kind", "size", "n_records", "ac",
            "max_an", "n_distinct_records", "member_pos", "member_size", "member_ac"]
    with open(a_.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in groups:
            m = sorted(g["members"], key=lambda r: r["pos"])
            fh.write("\t".join(str(v) for v in (
                g["chrom"], g["pos"], g["end"], g["kind"], g["anchor"]["size"],
                len(m), len(g["carriers"]), max(r["an"] for r in m),
                len(g["recs"]),
                ",".join(str(r["pos"]) for r in m),
                ",".join(str(r["size"]) for r in m),
                ",".join(str(len(r["carriers"])) for r in m))) + "\n")
    multi = [g for g in groups if len(g["members"]) > 1]
    print(f"alleles of at least {a_.min_size} bp: {len(alleles)}")
    print(f"events after merging: {len(groups)}")
    print(f"events joined from two or more records: {len(multi)}")
    for g in sorted(multi, key=lambda g: -len(g["carriers"]))[:20]:
        m = sorted(g["members"], key=lambda r: r["pos"])
        print(f"  {g['chrom']}:{g['pos']} {g['kind']} {g['anchor']['size']} bp  "
              f"alleles {len(m)} from {len(g['recs'])} records  ac {len(g['carriers'])}  "
              f"pieces {'+'.join(str(len(r['carriers'])) for r in m)}")
    print(f"wrote {a_.out}")


if __name__ == "__main__":
    main()
