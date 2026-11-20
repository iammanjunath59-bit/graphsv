#!/usr/bin/env python3
"""
gfa_locus_dosage.py

Measures, for every haplotype in a pangenome GFA, how much of a reference
interval it carries, as a dosage relative to the reference itself. Built to
count carriers of the 61,948 bp deletion at chr20:25,783,207 independently of
any VCF, on a graph from a different builder.

Reference segments are split into three sets along the reference path:
  L  the flank before the event        [start - flank, start)
  D  the event interval                [start, end)
  R  the flank after the event         [end, end + flank)
Segments shorter than --min-seg are ignored, and any segment that the
reference visits in more than one set is dropped from all three, so that
collapsed copies cannot be counted as both flank and event.

For each haplotype path, steps are taken between its first and last visit to
a flank segment (so visits to the same segments elsewhere on the chromosome
are not counted). The path is complete if it visits both L and R. Its event
dosage is the summed length times visits of D segments, divided by the same
sum for the reference path. A haplotype lacking the event should sit well
below the reference; in a segmental duplication collapsed by the builder it
need not reach zero, which is why the dosage distribution is printed rather
than a fixed cut being assumed.

Usage (reads the GFA twice, about as long as gfa_region.py):
    python3 gfa_locus_dosage.py --gfa chr20.gfa --ref-prefix GRCh38 \\
        --event 25783207-25845155 --out out/hprc/graph/del62k

Outputs <out>.paths.tsv (one row per path touching the locus) and
<out>.summary.txt (dosage histogram).
"""

import argparse
import sys
from collections import Counter, defaultdict


def parse_walk(w):
    out, i = [], 0
    while i < len(w):
        j = i + 1
        while j < len(w) and w[j] not in "><":
            j += 1
        out.append(w[i + 1:j])
        i = j
    return out


def parse_path(p):
    return [t[:-1] for t in p.split(",") if t]


def paths(gfa):
    with open(gfa) as fh:
        for line in fh:
            c = line[0]
            if c not in "PW":
                continue
            f = line.rstrip("\n").split("\t")
            if c == "P":
                yield f[1], parse_path(f[2])
            else:
                yield f"{f[1]}#{f[2]}#{f[3]}", parse_walk(f[6])


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--ref-prefix", required=True)
    ap.add_argument("--event", required=True, help="START-END, 0 based, reference path coordinates")
    ap.add_argument("--flank", type=int, default=20000)
    ap.add_argument("--right-start", type=int, default=None,
                    help="start of the right flank window (default: event end)")
    ap.add_argument("--min-seg", type=int, default=50)
    ap.add_argument("--out", required=True)
    a_ = ap.parse_args()
    s, e = (int(x) for x in a_.event.split("-"))

    seg_len = {}
    with open(a_.gfa) as fh:
        for line in fh:
            if line[0] == "S":
                f = line.split("\t", 3)
                seg_len[f[1]] = len(f[2].rstrip("\n")) if f[2] != "*" else 0
    sys.stderr.write(f"pass1: {len(seg_len):,} segments\n")

    ref_name, ref_walk = None, None
    for name, walk in paths(a_.gfa):
        if name.startswith(a_.ref_prefix):
            ref_name, ref_walk = name, walk
            break
    if ref_walk is None:
        sys.exit(f"no path starting with {a_.ref_prefix}")

    # flank windows: left [s - flank, s), right [rs, rs + flank) where rs
    # defaults to the event end but can be moved past a repeat with --right-start
    rs = a_.right_start if a_.right_start is not None else e
    lo, hi = s - a_.flank, rs + a_.flank
    cats = defaultdict(set)
    visits = Counter()
    pos = 0
    for sid in ref_walk:
        L = seg_len.get(sid, 0)
        if pos + L > lo and pos < hi and L >= a_.min_seg:
            mid = pos + L / 2
            cats[sid].add("L" if mid < s else ("D" if mid < e else ("R" if mid >= rs else "G")))
            visits[sid] += 1
        pos += L
        if pos > hi:
            break
    # flank segments must be visited once by the reference, so they mark a
    # single place; event segments may repeat (a collapsed duplication)
    Lset = {k for k, v in cats.items() if v == {"L"} and visits[k] == 1}
    Dset = {k for k, v in cats.items() if v == {"D"}}
    Rset = {k for k, v in cats.items() if v == {"R"} and visits[k] == 1}
    mixed = sum(1 for v in cats.values() if len(v) > 1)
    rep_flank = sum(1 for k, v in cats.items() if v in ({"L"}, {"R"}) and visits[k] > 1)
    sys.stderr.write(f"reference {ref_name}: L {len(Lset)}  D {len(Dset)}  R {len(Rset)} "
                     f"segments >= {a_.min_seg} bp; {mixed} dropped as visited in more than one "
                     f"set, {rep_flank} flank segments dropped as visited more than once\n")
    sys.stderr.write(f"windows: left {s - a_.flank}-{s}  event {s}-{e}  right {rs}-{rs + a_.flank}\n")
    if not Lset or not Rset:
        sys.exit("a flank has no unique segments; widen --flank or move --right-start")
    flank = Lset | Rset

    def measure(walk):
        idx = [i for i, sid in enumerate(walk) if sid in flank]
        if not idx:
            return None
        sub = walk[idx[0]: idx[-1] + 1]
        hasL = any(sid in Lset for sid in sub)
        hasR = any(sid in Rset for sid in sub)
        c = Counter(sid for sid in sub if sid in Dset)
        dbp = sum(seg_len[sid] * n for sid, n in c.items())
        # full sequence length between the outermost flank visits, every
        # segment counted including novel sequence: the haplotype's own span
        span = sum(seg_len.get(sid, 0) for sid in sub)
        return hasL, hasR, dbp, span

    ref_m = measure(ref_walk)
    ref_d = ref_m[2] if ref_m else 0
    if ref_m:
        sys.stderr.write(f"reference span between flanks: {ref_m[3]:,} bp\n")
    if not ref_d:
        sys.exit("reference has no event segments between its flanks; check --event")
    sys.stderr.write(f"reference event dosage basis: {ref_d:,} bp\n")

    rows = []
    for name, walk in paths(a_.gfa):
        m = measure(walk)
        if m is None:
            continue
        hasL, hasR, dbp, span = m
        parts = name.split("#")
        rows.append({"path": name, "sample": parts[0],
                     "hap": parts[1] if len(parts) > 2 else "",
                     "complete": int(hasL and hasR), "event_bp": dbp,
                     "dosage": round(dbp / ref_d, 4), "span_bp": span})

    with open(a_.out + ".paths.tsv", "w") as fh:
        cols = ["path", "sample", "hap", "complete", "event_bp", "dosage", "span_bp"]
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    comp = [r for r in rows if r["complete"]]
    haps = Counter((r["sample"], r["hap"]) for r in comp)
    L = [f"paths touching the locus: {len(rows)}   complete (both flanks): {len(comp)}",
         f"haplotypes with a complete path: {len(haps)}   "
         f"with more than one complete path: {sum(1 for v in haps.values() if v > 1)}",
         "", "dosage of the event interval relative to the reference, complete paths",
         "  bin           paths"]
    bins = Counter(min(int(r["dosage"] * 10), 20) for r in comp)
    for b in range(0, 21):
        if bins.get(b):
            lab = f"{b/10:.1f} to {(b+1)/10:.1f}" if b < 20 else ">= 2.0"
            L.append(f"  {lab:12} {bins[b]:6d}  " + "#" * min(bins[b], 60))
    rep = "\n".join(L)
    open(a_.out + ".summary.txt", "w").write(rep + "\n")
    print(rep)
    print(f"\nwrote {a_.out}.paths.tsv and {a_.out}.summary.txt")


if __name__ == "__main__":
    main()
