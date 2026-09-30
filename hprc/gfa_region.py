#!/usr/bin/env python3
"""
gfa_region.py -- extract the bubble structure of a target reference
interval from a large pangenome GFA, by streaming.

Why streaming
-------------
An HPRC per-chromosome GFA is ~2.4 GB compressed, 15-25 GB as text, with
hundreds of haplotype paths. Loading every path into memory (as graphsv.py
does, which is fine at 20 haplotypes x 1 Mb) would need hundreds of GB.
But the question here is local: what does the graph do over a few kb? So we
make ONE pass, keep only what that interval needs, and discard the rest.

Pass structure
--------------
The file is read twice because segment lengths are needed before path
steps can be converted to reference coordinates, and S lines are not
guaranteed to precede P/W lines.

  pass 1: S lines only -> segment lengths (one int per segment; for a human
          chromosome graph this is tens of millions of entries, a few GB in
          a plain dict, so lengths are stored in an array keyed by a compact
          id map to keep it manageable).
  pass 2: locate the reference path, walk it accumulating coordinates, and
          record the set of segments falling inside the target interval.
          Then re-scan paths, keeping for each haplotype only its steps
          that touch those segments.

What it answers
---------------
The motivating observation: HPRC v2.1's vcfwaved VCF reports four records
near chr20:25,783,000-25,787,500, all ~61,950 bp. Either that is ONE large
event shattered into four records by VCF conversion (the failure mode this
project quantifies in simulation), or it is four genuine near-identical
variants in a duplicated array (which is what segmental duplications
actually contain). The coordinates alone cannot distinguish these.

The graph can: if all four positions fall inside a single bubble whose
alternate traversal is carried by one consistent set of haplotypes, it is
one event. If there are four separate bubbles with different carrier sets,
they are four events and the VCF is right.

Usage
-----
  zstdcat chr20.gfa.zst | python gfa_region.py \
      --region 25780000-25790000 --ref-prefix GRCh38 --out chr20_region

  # or from an uncompressed file
  python gfa_region.py --gfa chr20.gfa --region ... --ref-prefix GRCh38
"""

import argparse
import sys
from collections import defaultdict


def open_input(path):
    if path is None or path == "-":
        return sys.stdin
    return open(path)


def parse_walk(w):
    """W-line walk string '>12<13>14' -> [(id, orient), ...]"""
    out, i = [], 0
    while i < len(w):
        o = w[i]
        j = i + 1
        while j < len(w) and w[j] not in "><":
            j += 1
        out.append((w[i + 1:j], "+" if o == ">" else "-"))
        i = j
    return out


def parse_path(p):
    """P-line path string '12+,13-,14+' -> [(id, orient), ...]"""
    return [(t[:-1], t[-1]) for t in p.split(",") if t]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", default=None,
                    help="GFA file; omit or '-' to read stdin (allows "
                         "`zstdcat x.gfa.zst | ...`). Requires --two-pass-file "
                         "if reading stdin is not possible twice.")
    ap.add_argument("--region", required=True,
                    help="START-END in reference-path coordinates, 0-based")
    ap.add_argument("--ref-prefix", required=True,
                    help="prefix of the reference path name, e.g. GRCh38")
    ap.add_argument("--pad", type=int, default=5000)
    ap.add_argument("--min-seg", type=int, default=50)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lo, hi = (int(x) for x in args.region.split("-"))
    lo -= args.pad
    hi += args.pad

    if args.gfa in (None, "-"):
        sys.exit("two passes are needed, so a seekable file is required.\n"
                 "Decompress first, or write to a temp file:\n"
                 "  zstdcat chr20.gfa.zst > /tmp/chr20.gfa")

    # ---------------- pass 1: segment lengths -------------------------
    seg_len = {}
    n = 0
    with open(args.gfa) as fh:
        for line in fh:
            if line[0] != "S":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            seg_len[f[1]] = len(f[2]) if f[2] != "*" else 0
            n += 1
            if n % 5_000_000 == 0:
                print(f"  pass1: {n:,} segments", file=sys.stderr)
    print(f"pass1 done: {len(seg_len):,} segments", file=sys.stderr)

    # ---------------- pass 2a: reference walk -> coordinates ----------
    ref_name, ref_walk = None, None
    with open(args.gfa) as fh:
        for line in fh:
            c = line[0]
            if c not in "PW":
                continue
            f = line.rstrip("\n").split("\t")
            if c == "P":
                name = f[1]
                if not name.startswith(args.ref_prefix):
                    continue
                ref_name, ref_walk = name, parse_path(f[2])
                break
            else:
                name = f"{f[1]}#{f[2]}#{f[3]}"
                if not name.startswith(args.ref_prefix):
                    continue
                ref_name, ref_walk = name, parse_walk(f[6])
                break
    if ref_walk is None:
        sys.exit(f"no path starting with '{args.ref_prefix}' found")
    print(f"reference path: {ref_name} ({len(ref_walk):,} steps)",
          file=sys.stderr)
    pos = 0
    target = []
    for sid, orient in ref_walk:
        L = seg_len.get(sid, 0)
        if pos + L > lo and pos < hi:
            target.append((sid, orient, pos, pos + L))
        pos += L
        if pos > hi:
            break
    if not target:
        sys.exit(f"interval {lo}-{hi} not covered by the reference path "
                 f"(path spans 0-{pos:,})")
    tset = {sid for sid, _, _, _ in target}
    print(f"reference passes through {len(target)} segments in the interval",
          file=sys.stderr)

    # ---------------- pass 2b: every haplotype's steps there ----------
    hits = defaultdict(list)
    npath = 0
    with open(args.gfa) as fh:
        for line in fh:
            c = line[0]
            if c not in "PW":
                continue
            f = line.rstrip("\n").split("\t")
            if c == "P":
                name, walk = f[1], parse_path(f[2])
            else:
                name, walk = f"{f[1]}#{f[2]}#{f[3]}", parse_walk(f[6])
            npath += 1
            if npath % 100 == 0:
                print(f"  pass2: {npath} paths", file=sys.stderr)
            # Compare only segments >= min_seg bp. Including 1 bp SNP
            # bubbles makes every haplotype's traversal unique and the
            # allele count collapses to n_haplotypes -- the artefact that
            # made 489 haplotypes look like 424 distinct alleles.
            keep = [(i, sid, o) for i, (sid, o) in enumerate(walk)
                    if sid in tset and seg_len.get(sid, 0) >= args.min_seg]
            if keep:
                hits[name] = keep

    print(f"pass2 done: {npath} paths, {len(hits)} touch the interval",
          file=sys.stderr)

    # ---------------- report ------------------------------------------
    # Group haplotypes by the traversal they take through the interval.
    # Identical traversal == same allele. If the four VCF records are one
    # shattered event, the haplotypes will fall into few groups spanning
    # the whole interval; if they are four real variants, the grouping
    # will be finer and the carrier sets will differ between sub-intervals.
    sig = defaultdict(list)
    for name, keep in hits.items():
        sig[tuple((sid, o) for _, sid, o in keep)].append(name)

    with open(f"{args.out}.alleles.tsv", "w") as fh:
        fh.write("allele\tn_haps\tn_steps\thaps\n")
        for k, (trav, names) in enumerate(
                sorted(sig.items(), key=lambda kv: -len(kv[1]))):
            fh.write(f"A{k}\t{len(names)}\t{len(trav)}\t"
                     f"{','.join(sorted(names)[:40])}\n")

    print(f"\ndistinct traversals through the interval: {len(sig)}")
    for k, (trav, names) in enumerate(
            sorted(sig.items(), key=lambda kv: -len(kv[1]))[:10]):
        span = sum(seg_len.get(s, 0) for s, _ in trav)
        rev = sum(1 for _, o in trav if o == "-")
        print(f"  A{k}: {len(names):>4} haplotypes, {len(trav):>5} steps, "
              f"{span:>8,} bp, {rev} reverse")

    with open(f"{args.out}.refsegs.tsv", "w") as fh:
        fh.write("segment\torient\tref_start\tref_end\tlen\n")
        for sid, o, a, b in target:
            fh.write(f"{sid}\t{o}\t{a}\t{b}\t{b-a}\n")
    print(f"\nwrote {args.out}.alleles.tsv and {args.out}.refsegs.tsv")


if __name__ == "__main__":
    main()
