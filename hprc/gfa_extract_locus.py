#!/usr/bin/env python3
"""
gfa_extract_locus.py

Cuts one locus out of a large GFA so that graphsv can be run on it.

graphsv anchors bubbles on segments that every path visits exactly once.
In a graph of hundreds of haplotypes, a haplotype that ends inside the
locus, or is split into two contigs there, removes every anchor. This
script therefore keeps only paths that pass through both flanks of the
locus, one per haplotype, and trims each to the stretch between its
outermost flank visits.

Flanks are chosen as in gfa_locus_dosage.py: reference segments of at least
--min-seg bp in [start - flank, start) and [right_start, right_start +
flank) that the reference visits once within the window. Each path is cut
from the outermost left flank segment it contains to the outermost right
flank segment, using the closest pair of occurrences: flank segments can
recur elsewhere on the chromosome, and cutting at the first and last flank
visit would then carry the locus far beyond its end. A path that reaches
the right flank before the left runs against the reference and is reversed
(order reversed, orientation flipped), so that graphsv does not read its
alleles as inversions. Haplotypes with more than one complete path are left out and
counted.

The file is read three times: segment lengths, paths, then the sequences of
the segments kept. Memory holds the length of every segment and the kept
sub walks only.

Outputs
  <out>.gfa          S lines (real sequence) and one P line per kept path
  <out>.offset.txt   reference coordinate of the first base of the trimmed
                     reference path, to convert graphsv positions
  <out>.log.txt      counts of kept and excluded paths

Usage:
    python3 hprc/gfa_extract_locus.py --gfa chr20.gfa --ref-prefix GRCh38 \\
        --event 25783207-25845155 --right-start 25850000 --flank 10000 \\
        --out out/hprc/graph/del62k_locus
"""

import argparse
import sys
from collections import Counter, defaultdict


def walk_tokens(line):
    """(name, [(segment, '+'|'-'), ...]) from a P or W line."""
    f = line.rstrip("\n").split("\t")
    if line[0] == "P":
        return f[1], [(t[:-1], t[-1]) for t in f[2].split(",") if t]
    w, out, i = f[6], [], 0
    while i < len(w):
        j = i + 1
        while j < len(w) and w[j] not in "><":
            j += 1
        out.append((w[i + 1:j], "+" if w[i] == ">" else "-"))
        i = j
    return f"{f[1]}#{f[2]}#{f[3]}", out


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--ref-prefix", required=True)
    ap.add_argument("--event", required=True, help="START-END, reference path coordinates")
    ap.add_argument("--flank", type=int, default=10000)
    ap.add_argument("--right-start", type=int, default=None)
    ap.add_argument("--min-seg", type=int, default=50)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    s, e = (int(x) for x in a.event.split("-"))
    rs = a.right_start if a.right_start is not None else e
    lo, hi = s - a.flank, rs + a.flank

    seg_len = {}
    with open(a.gfa) as fh:
        for line in fh:
            if line[0] == "S":
                f = line.split("\t", 3)
                seg_len[f[1]] = len(f[2].rstrip("\n")) if f[2] != "*" else 0
    sys.stderr.write(f"pass 1: {len(seg_len):,} segments\n")

    # pass 2: reference first, then every path
    ref = None
    with open(a.gfa) as fh:
        for line in fh:
            if line[0] in "PW":
                name, walk = walk_tokens(line)
                if name.startswith(a.ref_prefix):
                    ref = (name, walk)
                    break
    if ref is None:
        sys.exit(f"no path starting with {a.ref_prefix}")
    ref_name, ref_walk = ref

    cats, visits, pos, starts = defaultdict(set), Counter(), 0, {}
    for k, (sid, _) in enumerate(ref_walk):
        L = seg_len.get(sid, 0)
        if pos + L > lo and pos < hi and L >= a.min_seg:
            mid = pos + L / 2
            cats[sid].add("L" if mid < s else ("D" if mid < e else ("R" if mid >= rs else "G")))
            visits[sid] += 1
        starts.setdefault(k, pos)
        pos += L
        if pos > hi:
            break
    Lset = {k for k, v in cats.items() if v == {"L"} and visits[k] == 1}
    Rset = {k for k, v in cats.items() if v == {"R"} and visits[k] == 1}
    if not Lset or not Rset:
        sys.exit("a flank has no segment the reference visits once; widen --flank")
    flank = Lset | Rset
    sys.stderr.write(f"reference {ref_name}: left flank {len(Lset)}, right flank {len(Rset)} segments\n")

    # flank segments in reference order, outermost first on each side
    ref_order = []
    for sid, _ in ref_walk:
        if sid in flank and sid not in ref_order:
            ref_order.append(sid)
    L_order = [x for x in ref_order if x in Lset]            # left to right
    R_order = [x for x in ref_order if x in Rset][::-1]      # right to left

    def trim_tokens(toks):
        """toks: list of 'id+'/'id-' strings. The kept stretch runs from the
        outermost left flank segment the path contains to the outermost
        right flank segment, choosing the closest pair of occurrences, so a
        flank segment that the path also visits elsewhere on the chromosome
        does not stretch the locus. Only that stretch is turned into
        (segment, orientation) pairs."""
        pos = defaultdict(list)
        for i, t in enumerate(toks):
            sid = t[:-1]
            if sid in flank:
                pos[sid].append(i)
        if not pos:
            return None, "no flank"
        a = next((x for x in L_order if x in pos), None)
        b = next((x for x in R_order if x in pos), None)
        if a is None or b is None:
            return None, "one flank only"
        best = None
        for i in pos[a]:
            for j in pos[b]:
                if best is None or abs(j - i) < abs(best[1] - best[0]):
                    best = (i, j)
        i, j = best
        if i <= j:
            return [(t[:-1], t[-1]) for t in toks[i: j + 1]], "forward"
        sub = [(t[:-1], t[-1]) for t in toks[j: i + 1]]
        return [(sid, "-" if o == "+" else "+") for sid, o in reversed(sub)], "reversed"

    def line_tokens(line):
        f = line.rstrip("\n").split("\t")
        if line[0] == "P":
            return f[1], [t for t in f[2].split(",") if t]
        w, out, i = f[6], [], 0
        while i < len(w):
            j = i + 1
            while j < len(w) and w[j] not in "><":
                j += 1
            out.append(w[i + 1:j] + ("+" if w[i] == ">" else "-"))
            i = j
        return f"{f[1]}#{f[2]}#{f[3]}", out

    kept, why, per_hap = {}, Counter(), defaultdict(list)
    with open(a.gfa) as fh:
        for line in fh:
            if line[0] not in "PW":
                continue
            name, toks = line_tokens(line)
            sub, how = trim_tokens(toks)
            del toks
            why[how] += 1
            if sub is None:
                continue
            parts = name.split("#")
            hap = (parts[0], parts[1] if len(parts) > 2 else "")
            per_hap[hap].append(name)
            kept[name] = sub
    multi = {h for h, v in per_hap.items() if len(v) > 1}
    for h in multi:
        for n in per_hap[h]:
            if not n.startswith(a.ref_prefix):
                kept.pop(n, None)
    if ref_name not in kept:
        sys.exit("the reference path was not kept; check --event and flanks")

    rep = Counter(sid for sid, _ in kept[ref_name] if sid in flank)
    rep_n = sum(1 for v in rep.values() if v > 1)

    # offset of the trimmed reference path
    first = kept[ref_name][0][0]
    off, p = None, 0
    for sid, _ in ref_walk:
        if sid == first:
            off = p
            break
        p += seg_len.get(sid, 0)

    need = {sid for w in kept.values() for sid, _ in w}
    seqs = {}
    with open(a.gfa) as fh:
        for line in fh:
            if line[0] == "S":
                f = line.rstrip("\n").split("\t")
                if f[1] in need:
                    seqs[f[1]] = f[2]
    with open(a.out + ".gfa", "w") as fh:
        fh.write("H\tVN:Z:1.0\n")
        for sid in sorted(need, key=lambda x: (0, int(x), "") if x.isdigit() else (1, 0, x)):
            fh.write(f"S\t{sid}\t{seqs.get(sid, '*')}\n")
        for name, w in kept.items():
            fh.write(f"P\t{name}\t{','.join(sid + o for sid, o in w)}\t*\n")
    open(a.out + ".offset.txt", "w").write(f"{off}\n")
    ref_bp = sum(seg_len.get(sid, 0) for sid, _ in kept[ref_name])
    log = [f"reference {ref_name}; trimmed reference starts at {off}, "
           f"{ref_bp:,} bp, {rep_n} flank segments visited more than once",
           f"paths seen: {sum(why.values())}",
           *[f"  {k}: {v}" for k, v in sorted(why.items())],
           f"haplotypes with more than one complete path (left out): {len(multi)}",
           f"paths written: {len(kept)} (including the reference); segments: {len(need):,}"]
    open(a.out + ".log.txt", "w").write("\n".join(log) + "\n")
    print("\n".join(log))
    print(f"wrote {a.out}.gfa, .offset.txt, .log.txt")


if __name__ == "__main__":
    main()
