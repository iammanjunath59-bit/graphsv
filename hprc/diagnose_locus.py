#!/usr/bin/env python3
"""
diagnose_locus.py

Explains what graphsv sees at one locus of a large graph. Reads the GFA
written by gfa_extract_locus.py and reports, along the reference path:

  * segment lengths of the reference inside the event interval (how much of
    it lies in segments of at least --min-seg bp, which graphsv compares)
  * graphsv anchors (segments visited exactly once by every path) and where
    they fall relative to the event
  * for each pair of consecutive anchors that overlaps the event: the number
    of distinct haplotype traversals, with and without the --min-seg filter,
    and the largest length difference from the reference traversal

The graphsv functions are imported from src/graphsv.py, so the anchors and
traversals are exactly those graphsv uses.

Usage:
    python3 hprc/diagnose_locus.py --gfa out/hprc/graph/del62k_locus.gfa \\
        --offset out/hprc/graph/del62k_locus.offset.txt --event 25783207-25845155
"""

import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
import graphsv  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--offset", required=True)
    ap.add_argument("--event", required=True)
    ap.add_argument("--ref-prefix", default="GRCh38")
    ap.add_argument("--min-seg", type=int, default=50)
    a = ap.parse_args()
    s, e = (int(x) for x in a.event.split("-"))
    off = int(open(a.offset).read().split()[0])

    seg_len, paths = graphsv.parse_gfa(a.gfa)
    ref = next(n for n in paths if n.startswith(a.ref_prefix))
    rw = paths[ref]
    haps = [n for n in paths if n != ref]
    print(f"paths {len(paths)} (haplotypes {len(haps)}), segments {len(seg_len):,}")

    pos, coords = 0, []
    for sid, _ in rw:
        coords.append(pos + off)
        pos += seg_len.get(sid, 0)

    inside = [(sid, seg_len.get(sid, 0)) for (sid, _), c in zip(rw, coords) if s <= c < e]
    tot = sum(L for _, L in inside)
    print(f"\nreference inside the event: {len(inside):,} segments, {tot:,} bp")
    for t in (1, 10, 20, 50, 100):
        k = [L for _, L in inside if L >= t]
        print(f"  segments >= {t:3d} bp: {len(k):6,}  covering {sum(k):8,} bp ({sum(k)/max(tot,1):.3f})")

    anchors = graphsv.find_anchors(paths)
    steps = [i for i, (sid, _) in enumerate(rw) if sid in anchors]
    ac = [coords[i] for i in steps]
    print(f"\nanchors {len(anchors)}: before the event {sum(c < s for c in ac)}, "
          f"inside {sum(s <= c < e for c in ac)}, after {sum(c >= e for c in ac)}")
    if any(s <= c < e for c in ac):
        print("  anchors inside the event at " + ", ".join(f"{c:,}" for c in ac if s <= c < e)[:600])

    index = {n: {sid: i for i, (sid, _) in enumerate(w)} for n, w in paths.items()}

    def trav(n, l, r, ms):
        ix = index[n]
        if l not in ix or r not in ix:
            return None
        x, y = sorted((ix[l], ix[r]))
        return tuple(t for t in paths[n][x + 1:y] if seg_len.get(t[0], 0) >= ms)

    def tl(t):
        return sum(seg_len.get(x, 0) for x, _ in t)

    print(f"\nbubbles (consecutive anchors) overlapping the event")
    print("  ref start - end            ref bp  distinct traversals: all segs / >= min-seg   "
          "largest length difference (>= min-seg)")
    for k in range(len(steps) - 1):
        i, j = steps[k], steps[k + 1]
        lo, hi = coords[i], coords[j]
        if not (hi > s and lo < e):
            continue
        l, r = rw[i][0], rw[j][0]
        full = Counter(trav(n, l, r, 0) for n in haps)
        filt = Counter(trav(n, l, r, a.min_seg) for n in haps)
        rt = trav(ref, l, r, a.min_seg)
        diff = max((abs(tl(t) - tl(rt)) for t in filt if t is not None), default=0)
        print(f"  {lo:,} - {hi:,}  {hi - lo:9,}  {len(full):5d} / {len(filt):5d}   {diff:,}"
              f"   (haplotypes without both anchors: {full.get(None, 0)})")


if __name__ == "__main__":
    main()
