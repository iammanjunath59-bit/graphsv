#!/usr/bin/env python3
"""
graph_ac.py

Recover structural-variant allele counts directly from a pangenome graph
by comparing haplotype path traversals, bypassing VCF conversion entirely.

The claim being tested
----------------------
Deconstructing a graph to VCF flattens nested bubble structure into linear
records. For large SVs this is lossy: the PGGB pipeline (vg deconstruct +
vcfbub + vcfwave) fails to represent the true allele count for ~45% of
30 kb SVs. Minigraph's native traversal-based genotyping recovers 23/23 of
the same SVs exactly. That suggests the information survives in the graph
and is destroyed at the VCF step -- but minigraph and PGGB build DIFFERENT
graphs, so that comparison confounds builder with representation.

This script closes the confound: it computes allele counts from PGGB's own
GFA by traversal identity. If PGGB's graph yields correct counts while
PGGB's VCF does not, the loss is unambiguously in the conversion.

Method
------
A GFA path is a sequence of oriented segments. For an interval on the
reference path, each haplotype's traversal is the segment list it takes
between the two flanking segments that all haplotypes share (the bubble
boundary). Two haplotypes carry the same allele iff their traversals are
identical. Allele count = number of haplotypes whose traversal differs
from the reference haplotype's.

This needs no variant calling, no realignment, and no allele
normalisation -- the three places VCF conversion loses information.

Boundary finding: walk outward from the SV's reference-coordinate span
until reaching segments visited by every haplotype exactly once. Those
anchors define the bubble. If no such anchor exists within --max-flank,
the site is reported unresolved rather than guessed at.

Usage
-----
  python graph_ac.py \
      --gfa runs/pggb_x/*.smooth.final.gfa \
      --truth runs/x.truth_svs.tsv \
      --ref-prefix anc \
      --out runs/x.graph_ac.tsv
"""

import argparse
import sys
from collections import defaultdict


def parse_gfa(path):
    """Return (seg_len, paths) where paths maps name -> list of segment ids.

    Handles both P lines (explicit paths) and W lines (walks, which is what
    modern vg/PGGB emit). Orientation is kept in the segment token so that
    an inversion is not silently treated as the reference allele.
    """
    seg_len, paths = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                f = line.rstrip("\n").split("\t")
                seg_len[f[1]] = len(f[2]) if len(f) > 2 and f[2] != "*" else 0
            elif line.startswith("P\t"):
                f = line.rstrip("\n").split("\t")
                name = f[1]
                toks = [t for t in f[2].split(",") if t]
                paths[name] = [(t[:-1], t[-1]) for t in toks]
            elif line.startswith("W\t"):
                f = line.rstrip("\n").split("\t")
                # W <sample> <hap> <seq> <start> <end> <walk>
                name = f"{f[1]}#{f[2]}#{f[3]}"
                walk, cur = [], f[6]
                i = 0
                while i < len(cur):
                    o = cur[i]
                    j = i + 1
                    while j < len(cur) and cur[j] not in "><":
                        j += 1
                    walk.append((cur[i + 1:j], "+" if o == ">" else "-"))
                    i = j
                paths[name] = walk
    return seg_len, paths


def ref_offsets(walk, seg_len):
    """Cumulative reference coordinate at the START of each step."""
    offs, pos = [], 0
    for sid, _ in walk:
        offs.append(pos)
        pos += seg_len.get(sid, 0)
    return offs, pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--ref-prefix", default="anc")
    ap.add_argument("--max-flank", type=int, default=200000,
                    help="how far to search for a shared anchor segment")
    ap.add_argument("--min-seg", type=int, default=50,
                    help="ignore segments shorter than this "
                         "when comparing traversals")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    MINSEG = args.min_seg
    seg_len, paths = parse_gfa(args.gfa)
    if not paths:
        sys.exit("no P or W lines found in GFA")

    ref_name = next((n for n in paths if n.startswith(args.ref_prefix)), None)
    if ref_name is None:
        sys.exit(f"no path starting with '{args.ref_prefix}'. "
                 f"Found: {sorted(paths)[:5]}")
    hap_names = [n for n in paths if n != ref_name]
    n_hap = len(hap_names)
    print(f"reference path: {ref_name}", file=sys.stderr)
    print(f"haplotypes: {n_hap}", file=sys.stderr)

    ref_walk = paths[ref_name]
    offs, total = ref_offsets(ref_walk, seg_len)
    print(f"reference path length: {total:,} bp", file=sys.stderr)

    # segments visited exactly once by every haplotype = safe anchors
    visits = defaultdict(lambda: defaultdict(int))
    for name, walk in paths.items():
        for sid, _ in walk:
            visits[sid][name] += 1
    anchors = {sid for sid, v in visits.items()
               if len(v) == len(paths) and all(c == 1 for c in v.values())}
    print(f"anchor segments: {len(anchors):,} / {len(seg_len):,}",
          file=sys.stderr)

    # index each haplotype's walk by segment for fast slicing
    index = {name: {sid: i for i, (sid, _) in enumerate(walk)}
             for name, walk in paths.items()}

    def find_anchor(start_step, direction):
        i = start_step
        while 0 <= i < len(ref_walk):
            sid = ref_walk[i][0]
            if sid in anchors:
                return i, sid
            i += direction
            if abs(offs[min(max(i, 0), len(offs) - 1)]
                   - offs[min(max(start_step, 0), len(offs) - 1)]) \
                    > args.max_flank:
                return None, None
        return None, None

    def step_at(pos):
        lo, hi = 0, len(offs) - 1
        while lo < hi:
            mid = (lo + hi + 1) // 2
            if offs[mid] <= pos:
                lo = mid
            else:
                hi = mid - 1
        return lo

    rows = []
    with open(args.truth) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            s, e = int(c[ix["start"]]), int(c[ix["end"]])
            ac_true = int(c[ix["ac"]])
            rec = {"sv_id": c[ix["sv_id"]], "type": c[ix["type"]],
                   "len": int(c[ix["len"]]),
                   "size_bin": c[ix.get("size_bin", 0)] if "size_bin" in ix else "",
                   "ac_true": ac_true}

            li, lsid = find_anchor(step_at(s), -1)
            ri, rsid = find_anchor(step_at(e), +1)
            if lsid is None or rsid is None or li >= ri:
                rec.update({"ac_graph": "", "n_alleles": "",
                            "status": "no_anchor", "exact": 0})
                rows.append(rec)
                continue

            def traversal(name):
                # Compare only segments above MINSEG bp. Anchors generally
                # enclose many linked SNP bubbles, each contributing a 1 bp
                # segment; including those makes every haplotype's traversal
                # unique and the allele count collapses to n_hap. Structural
                # differences live in the large segments.
                idx = index[name]
                if lsid not in idx or rsid not in idx:
                    return None
                a, b = idx[lsid], idx[rsid]
                if a > b:
                    a, b = b, a
                return tuple(t for t in paths[name][a + 1:b]
                             if seg_len.get(t[0], 0) >= MINSEG)

            ref_trav = traversal(ref_name)
            if ref_trav is None:
                rec.update({"ac_graph": "", "n_alleles": "",
                            "status": "ref_missing", "exact": 0})
                rows.append(rec)
                continue

            alleles, alt = defaultdict(int), 0
            missing = 0
            for h in hap_names:
                t = traversal(h)
                if t is None:
                    missing += 1
                    continue
                alleles[t] += 1
                if t != ref_trav:
                    alt += 1
            rec.update({
                "ac_graph": alt,
                "n_alleles": len(alleles),
                "status": "ok" if missing == 0 else f"missing{missing}",
                "exact": int(alt == ac_true),
            })
            rows.append(rec)

    cols = ["sv_id", "type", "len", "size_bin", "ac_true", "ac_graph",
            "n_alleles", "status", "exact"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    ok = [r for r in rows if r["status"] == "ok"]
    ex = sum(r["exact"] for r in rows)
    print(f"\nresolved {len(ok)}/{len(rows)} sites; "
          f"exact AC {ex}/{len(rows)} = {ex/max(1,len(rows)):.3f}",
          file=sys.stderr)
    bins = defaultdict(lambda: [0, 0])
    for r in rows:
        b = bins[r["size_bin"] or "NA"]
        b[0] += 1
        b[1] += r["exact"]
    for b in sorted(bins, key=lambda x: int(x) if str(x).isdigit() else 0):
        n, k = bins[b]
        print(f"  bin {b:>6}  n={n:3d}  exact={k/n:.3f}", file=sys.stderr)
    print(f"wrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
