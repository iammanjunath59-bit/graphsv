#!/usr/bin/env python3
"""
make_truth_vcf.py

Writes a phased, haploid VCF of a simulated panel's truth, for the vg
construct round trip: build a graph from the VCF, deconstruct
it, and check whether every structural variant comes back with its count.

Each structural variant becomes one record with one alternate allele:
  DEL  REF = anchor + deleted segment, ALT = anchor
  DUP  tandem: inserted after the segment, ALT = anchor + copy of the segment
  INV  REF = anchor + segment, ALT = anchor + reverse complement
Point mutations from truth_snps.tsv become ordinary SNP records.

A VCF cannot place a point mutation inside an inserted copy, because the
copy has no reference coordinates. In this VCF a point mutation inside a
duplicated segment therefore lies on the reference copy only, and every
carrier's inserted copy is identical. That is the point of the test: it asks
whether deconstruction is exact when the graph carries no variation nested
inside insertions.

Coordinates are checked, not assumed: the script tests whether SNP positions
are 0 or 1 based by comparing the recorded reference base with the backbone,
and stops if neither convention matches.

Point mutations that fall inside a deletion or inversion interval, or at a
structural variant's anchor base, overlap that record. By default they are
dropped (vg construct handles overlapping records poorly); the number dropped
is reported. Structural variant truth is unaffected.

Usage:
    python3 make_truth_vcf.py --prefix sweep54d/chr20_10000000_s1 \\
        --backbone data/windows/chr20_10000000.fa --out rt/chr20_10000000_s1.truth.vcf
"""

import argparse
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path):
    name, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    sys.exit("backbone FASTA must hold a single record")
                name = line[1:].split()[0]
            else:
                seq.append(line.strip())
    return name, "".join(seq).upper()


def read_tsv(path):
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True, help="replicate prefix holding .truth_svs.tsv and .truth_snps.tsv")
    ap.add_argument("--backbone", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-hap", type=int, default=20)
    ap.add_argument("--no-snps", action="store_true", help="structural variants only")
    ap.add_argument("--keep-overlapping-snps", action="store_true")
    a_ = ap.parse_args()

    chrom, bb = read_fasta(a_.backbone)
    svs = read_tsv(a_.prefix + ".truth_svs.tsv")
    snps = [] if a_.no_snps else read_tsv(a_.prefix + ".truth_snps.tsv")

    # SNP coordinate convention
    off = 0
    if snps:
        m0 = sum(bb[int(s["pos"])] == s["ref"].upper() for s in snps if int(s["pos"]) < len(bb))
        m1 = sum(bb[int(s["pos"]) - 1] == s["ref"].upper() for s in snps if 0 < int(s["pos"]) <= len(bb))
        n = len(snps)
        sys.stderr.write(f"SNP ref base matches backbone: 0 based {m0}/{n}, 1 based {m1}/{n}\n")
        if m0 == n:
            off = 0
        elif m1 == n:
            off = -1
        else:
            sys.exit("SNP positions match neither convention fully; stopping rather than guessing")

    recs = []      # (pos1, ref, alt, carriers, id)
    occupied = []  # 0 based [lo, hi) intervals an SNP must not overlap
    for s in svs:
        st, en, ty = int(s["start"]), int(s["end"]), s["type"].upper()
        car = {int(x) for x in s["carriers"].split(",") if x != ""}
        seg = bb[st:en]
        if len(seg) != en - st or st < 1:
            sys.exit(f"{s['sv_id']}: interval {st}-{en} outside the backbone")
        if ty == "DEL":
            recs.append((st, bb[st - 1:en], bb[st - 1], car, s["sv_id"]))
            occupied.append((st - 1, en))
        elif ty == "DUP":
            a = bb[en - 1]
            recs.append((en, a, a + seg, car, s["sv_id"]))
            occupied.append((en - 1, en))
        elif ty == "INV":
            recs.append((st, bb[st - 1:en], bb[st - 1] + seg.translate(COMP)[::-1], car, s["sv_id"]))
            occupied.append((st - 1, en))
        else:
            sys.exit(f"{s['sv_id']}: unknown type {ty}")

    occupied.sort()
    dropped = 0
    for s in snps:
        p0 = int(s["pos"]) + off
        if not a_.keep_overlapping_snps and any(lo <= p0 < hi for lo, hi in occupied):
            dropped += 1
            continue
        car = {int(x) for x in s["carriers"].split(",") if x != ""}
        recs.append((p0 + 1, bb[p0], s["alt"].upper(), car, "."))

    recs.sort(key=lambda r: r[0])
    seen = set()
    with open(a_.out, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##contig=<ID={chrom},length={len(bb)}>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(f"h{i}" for i in range(a_.n_hap)) + "\n")
        n_written = 0
        for pos, ref, alt, car, vid in recs:
            if (pos, ref, alt) in seen:
                continue
            seen.add((pos, ref, alt))
            gts = "\t".join("1" if i in car else "0" for i in range(a_.n_hap))
            fh.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alt}\t60\tPASS\t.\tGT\t{gts}\n")
            n_written += 1
    sys.stderr.write(f"contig {chrom}  structural variants {len(svs)}  point mutations "
                     f"{len(snps) - dropped} kept, {dropped} dropped as overlapping  "
                     f"records written {n_written}\n")


if __name__ == "__main__":
    main()
