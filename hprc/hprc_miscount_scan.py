#!/usr/bin/env python3
"""
hprc_miscount_scan.py

Look for the documented failure mode in the published human pangenome.

What the simulation established
-------------------------------
Two mechanisms account for VCF-route failures at large structural variants,
and one dominates. In the 30 kb size class, 39.4% of variants on chromosome
20 and 40.6% on chromosome 16 were MISCOUNTED: a record of essentially the
correct size sits at essentially the correct position, but its allele count
is wrong, because the carriers of a single event are partitioned across
several records. Truth AC 6 appearing as records of 1, 5 and 9 is typical.
Outright absence of the variant is comparatively rare (1.2% and 2.5%).

The distinguishing signature is therefore NOT that a large event is missing,
nor that records are scattered. It is that several records of near-identical
size occur in the same interval, each carrying part of the true carrier set,
so that no single record reports the number of haplotypes that actually
carry the event.

What this script does
---------------------
For a set of loci it extracts, from the graph, the grouping of haplotypes
into distinct structural alleles, and from the published VCF, the large
records present. It then asks whether any single VCF record reports an allele
count matching the graph's haplotype support for a corresponding allele.

Two earlier mistakes this avoids
--------------------------------
The extraction window must exceed the largest allele at the locus. A 62 kb
allele cannot be represented inside a 40 kb window, and comparing a truncated
graph traversal against a full-length VCF record produces a spurious
mismatch. Windows are therefore sized from the largest record present, with a
wide margin, and the script refuses to report a locus where the reference
walk is shorter than the largest allele.

And a multi-allelic locus correctly decomposed into biallelic records is not
a failure. Several records with uncorrelated per-sample dosages are distinct
alleles carried by distinct haplotypes, which is the correct representation.
The signature of interest is the opposite: records whose dosages ARE
correlated, meaning they partition one carrier set. Dosage correlation is
therefore reported alongside, and a locus is only flagged when correlated
records fail to reproduce the graph's haplotype support.

Usage
-----
  python hprc_miscount_scan.py --db <gbz.db> --vcf <wave.vcf.gz> \\
      --loci loci.txt --margin 100000 --out hprc_loci/miscount
"""

import argparse
import itertools
import os
import subprocess
import sys
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("needs numpy")


def run(cmd, timeout=900):
    try:
        return subprocess.run(cmd, capture_output=True, text=True,
                              timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        return ""


def vcf_records(vcf, region, min_len):
    """(pos, size, [acs], dosage array) for records with a large allele."""
    out = run(["bcftools", "query", "-r", region,
               "-f", "%POS\t%REF\t%ALT\t%INFO/AC[\t%GT]\n", vcf])
    recs = []
    for line in out.rstrip("\n").split("\n"):
        if not line.strip():
            continue
        f = line.split("\t")
        if len(f) < 5:
            continue
        ref, alts = f[1], f[2].split(",")
        size = max([len(ref)] + [len(a) for a in alts])
        if size < min_len:
            continue
        acs = []
        for a in f[3].split(","):
            try:
                acs.append(int(a))
            except ValueError:
                pass
        dos = []
        for g in f[4:]:
            g = g.replace("|", "/")
            al = [x for x in g.split("/") if x != "."]
            dos.append(sum(1 for x in al if x != "0") if al else np.nan)
        recs.append({"pos": int(f[0]), "size": size, "acs": acs,
                     "dos": np.array(dos, dtype=float)})
    return recs


def graph_alleles(db, chrom, start, end, ref_prefix, min_seg, workdir, label):
    """Distinct traversals and their haplotype counts, via gbz-base."""
    gfa = os.path.join(workdir, f"{label}.gfa")
    if not os.path.exists(gfa) or os.path.getsize(gfa) == 0:
        out = run(["gbz-base", "query", db, "--sample", ref_prefix,
                   "--contig", chrom, "-i", f"{start}..{end}",
                   "--context", "5000", "--format", "gfa"])
        if not out.strip():
            return None, 0
        with open(gfa, "w") as fh:
            fh.write(out)

    seg_len, walks = {}, {}
    for line in open(gfa):
        if line.startswith("S\t"):
            f = line.rstrip("\n").split("\t")
            seg_len[f[1]] = len(f[2]) if len(f) > 2 and f[2] != "*" else 0
        elif line.startswith("W\t"):
            f = line.rstrip("\n").split("\t")
            name = f"{f[1]}#{f[2]}#{f[3]}"
            w, i, steps = f[6], 0, []
            while i < len(w):
                j = i + 1
                while j < len(w) and w[j] not in "><":
                    j += 1
                steps.append((w[i + 1:j], "+" if w[i] == ">" else "-"))
                i = j
            walks[name] = steps
    if not walks:
        return None, 0

    ref = next((n for n in walks if n.startswith(ref_prefix)), None)
    if ref is None:
        return None, 0
    ref_span = sum(seg_len.get(s, 0) for s, _ in walks[ref])

    groups = defaultdict(list)
    for name, steps in walks.items():
        if name == ref:
            continue
        key = tuple(t for t in steps if seg_len.get(t[0], 0) >= min_seg)
        groups[key].append(name)
    alleles = sorted(
        ((len(v), sum(seg_len.get(s, 0) for s, _ in k)) for k, v in groups.items()),
        reverse=True)
    return alleles, ref_span


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--loci", required=True,
                    help="chrom:pos:label per line, centred on a large record")
    ap.add_argument("--margin", type=int, default=100000,
                    help="half-width of the extraction window; must exceed "
                         "the largest allele at the locus")
    ap.add_argument("--min-len", type=int, default=5000)
    ap.add_argument("--min-seg", type=int, default=50)
    ap.add_argument("--ref-prefix", default="GRCh38")
    ap.add_argument("--workdir", default="hprc_loci")
    ap.add_argument("--out", default="hprc_loci/miscount")
    args = ap.parse_args()

    os.makedirs(args.workdir, exist_ok=True)
    rows = []

    print(f"{'locus':<22} {'graph alleles':>28}   {'VCF records':>26}  {'dosage r':>9}  flag")
    for raw in open(args.loci):
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        p = raw.split(":")
        if len(p) < 2:
            continue
        chrom, pos = p[0], int(p[1])
        label = p[2] if len(p) > 2 else f"{chrom}_{pos}"
        s, e = max(0, pos - args.margin), pos + args.margin

        recs = vcf_records(args.vcf, f"{chrom}:{s}-{e}", args.min_len)
        if not recs:
            print(f"{label:<22} {'--':>28}   {'no large records':>26}")
            continue
        biggest = max(r["size"] for r in recs)

        alleles, ref_span = graph_alleles(
            args.db, chrom, s, e, args.ref_prefix, args.min_seg,
            args.workdir, label)
        if alleles is None:
            print(f"{label:<22} {'query failed':>28}")
            continue
        if ref_span < biggest:
            print(f"{label:<22} WINDOW TOO SMALL: ref walk {ref_span:,} bp "
                  f"< largest allele {biggest:,} bp; increase --margin")
            continue

        # dosage correlation among the large records
        rs = []
        for a, b in itertools.combinations(recs, 2):
            ok = ~(np.isnan(a["dos"]) | np.isnan(b["dos"]))
            if ok.sum() < 20:
                continue
            x, y = a["dos"][ok], b["dos"][ok]
            if x.std() == 0 or y.std() == 0:
                continue
            rs.append(float(np.corrcoef(x, y)[0, 1]))
        med_r = float(np.median(rs)) if rs else float("nan")

        top = alleles[:3]
        gstr = "; ".join(f"{n}hap/{L//1000}kb" for n, L in top)
        allacs = sorted({a for r in recs for a in r["acs"]}, reverse=True)[:5]
        vstr = f"{len(recs)} recs, AC {allacs}"

        # flag: correlated records (one partition) whose counts do not
        # reproduce the graph's haplotype support for any allele
        graph_counts = {n for n, _ in alleles}
        matched = any(a in graph_counts for r in recs for a in r["acs"])
        flag = ""
        if len(recs) > 1 and not np.isnan(med_r) and med_r > 0.5 and not matched:
            flag = "MISCOUNT"
        elif len(recs) > 1 and not np.isnan(med_r) and med_r < 0.2:
            flag = "distinct alleles"
        elif matched:
            flag = "consistent"

        print(f"{label:<22} {gstr:>28}   {vstr:>26}  {med_r:>9.3f}  {flag}")
        rows.append((label, chrom, pos, ref_span, len(recs), biggest,
                     med_r, gstr, str(allacs), flag))

    with open(f"{args.out}.tsv", "w") as fh:
        fh.write("label\tchrom\tpos\tref_span\tn_records\tlargest_allele\t"
                 "median_dosage_r\tgraph_alleles\tvcf_acs\tflag\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"\nwrote {args.out}.tsv")
    print("\nMISCOUNT: records are correlated (one carrier set) yet no record")
    print("reproduces the graph's haplotype support. distinct alleles:")
    print("uncorrelated records, i.e. a genuinely multi-allelic locus, which")
    print("is correct behaviour and not a failure.")


if __name__ == "__main__":
    main()
