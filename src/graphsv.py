#!/usr/bin/env python3
"""
graphsv.py -- call structural variant alleles and allele counts directly
from a pangenome graph, without deconstructing to VCF.

Motivation
----------
`vg deconstruct` (+ vcfbub + vcfwave) projects a graph's nested bubble
structure onto linear VCF records. That projection is lossy, and the loss
scales with variant size: in a 54-window, 2,592-SV benchmark, the true
allele count was unrecoverable from the VCF for ~45% of 30 kb SVs, versus
~1% of SVs below 1 kb. The same information is intact in the graph -- read
by traversal identity, the same PGGB graphs return the correct allele
count for 100% of SVs from 300 bp to 30 kb.

This tool reads the graph.

Algorithm
---------
1. ANCHORS. A segment traversed exactly once by every haplotype is an
   anchor. Anchors are points where all haplotypes agree, so consecutive
   anchors along the reference path bracket exactly one top-level bubble.
   This is a path-based equivalent of a top-level snarl decomposition and
   needs no external snarl computation.

   A large deletion removes the anchors inside it for carriers, so those
   segments stop being global anchors and the flanking anchors bracket the
   whole event. That is what makes large SVs come out as single sites
   rather than shattering.

2. TRAVERSALS. Within a bubble, each haplotype's allele is the sequence of
   oriented segments it takes between the flanking anchors. Two haplotypes
   carry the same allele iff their traversals are identical. Orientation is
   retained, so an inversion is a distinct allele rather than a match.

3. SEGMENT FILTER (--min-seg). Bubbles bracketed by distant anchors enclose
   many linked SNP bubbles, each contributing a ~1 bp segment. Comparing
   raw traversals then makes every haplotype unique and the allele count
   collapses to n_haplotypes -- this was the dominant failure mode during
   development. Restricting the comparison to segments of at least
   --min-seg bp removes incidental SNP variation and leaves structural
   differences. Results are stable across 20-50 bp; above ~200 bp genuine
   small SVs start being discarded.

4. CALLING. A bubble with more than one distinct filtered traversal is a
   site. The reference haplotype's traversal defines the REF allele; every
   other distinct traversal is an ALT with its own allele count.

Output
------
TSV, one row per polymorphic site:
  chrom  pos  end  ref_len  n_alleles  ac_total  allele_lens  allele_acs

`pos`/`end` are reference-path coordinates spanned by the bubble interior.
`ac_total` is the number of haplotypes not carrying the reference allele.

Usage
-----
  # call
  python graphsv.py --gfa graph.gfa --ref-prefix anc --out calls.tsv

  # call and evaluate against simulation truth
  python graphsv.py --gfa graph.gfa --ref-prefix anc --out calls.tsv \
                    --truth run.truth_svs.tsv
"""

__version__ = "1.0.0"

import argparse
import sys
from collections import defaultdict


# ----------------------------------------------------------------------
# GFA
# ----------------------------------------------------------------------

def parse_gfa(path):
    """Return (seg_len, paths). Handles both P lines and W lines (walks)."""
    seg_len, paths = {}, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("S\t"):
                f = line.rstrip("\n").split("\t")
                seg_len[f[1]] = len(f[2]) if len(f) > 2 and f[2] != "*" else 0
            elif line.startswith("P\t"):
                f = line.rstrip("\n").split("\t")
                toks = [t for t in f[2].split(",") if t]
                paths[f[1]] = [(t[:-1], t[-1]) for t in toks]
            elif line.startswith("W\t"):
                f = line.rstrip("\n").split("\t")
                name = f"{f[1]}#{f[2]}#{f[3]}"
                walk, cur, i = [], f[6], 0
                while i < len(cur):
                    o = cur[i]
                    j = i + 1
                    while j < len(cur) and cur[j] not in "><":
                        j += 1
                    walk.append((cur[i + 1:j], "+" if o == ">" else "-"))
                    i = j
                paths[name] = walk
    return seg_len, paths


# ----------------------------------------------------------------------
# calling
# ----------------------------------------------------------------------

def find_anchors(paths):
    """Segments traversed exactly once by every path."""
    visits = defaultdict(lambda: defaultdict(int))
    for name, walk in paths.items():
        for sid, _ in walk:
            visits[sid][name] += 1
    n = len(paths)
    return {sid for sid, v in visits.items()
            if len(v) == n and all(c == 1 for c in v.values())}


def call_sites(seg_len, paths, ref_name, min_seg, min_alt_len):
    ref_walk = paths[ref_name]
    haps = [n for n in paths if n != ref_name]

    anchors = find_anchors(paths)
    index = {name: {sid: i for i, (sid, _) in enumerate(walk)}
             for name, walk in paths.items()}

    # reference coordinate at the start of each reference step
    offs, pos = [], 0
    for sid, _ in ref_walk:
        offs.append(pos)
        pos += seg_len.get(sid, 0)
    ref_len_total = pos

    anchor_steps = [i for i, (sid, _) in enumerate(ref_walk)
                    if sid in anchors]

    def traversal(name, lsid, rsid):
        idx = index[name]
        if lsid not in idx or rsid not in idx:
            return None
        a, b = idx[lsid], idx[rsid]
        if a > b:
            a, b = b, a
        return tuple(t for t in paths[name][a + 1:b]
                     if seg_len.get(t[0], 0) >= min_seg)

    def tlen(t):
        return sum(seg_len.get(s, 0) for s, _ in t)

    sites = []
    for k in range(len(anchor_steps) - 1):
        i, j = anchor_steps[k], anchor_steps[k + 1]
        lsid, rsid = ref_walk[i][0], ref_walk[j][0]

        ref_trav = traversal(ref_name, lsid, rsid)
        if ref_trav is None:
            continue

        groups = defaultdict(list)
        for h in haps:
            t = traversal(h, lsid, rsid)
            if t is None:
                continue
            groups[t].append(h)
        if len(groups) <= 1 and (not groups or ref_trav in groups):
            continue                                    # monomorphic

        ref_l = tlen(ref_trav)
        alts = [(t, hs) for t, hs in groups.items() if t != ref_trav]
        if not alts:
            continue
        # Event size must NOT be length difference alone: a balanced
        # inversion has zero length difference, so a length-based filter
        # discards every inversion before it can be called. Size is the
        # larger of the length change and the amount of sequence that
        # differs in CONTENT between the two traversals.
        def event_size(t):
            ref_segs = {sid for sid, _ in ref_trav}
            alt_segs = {sid for sid, _ in t}
            uniq_alt = sum(seg_len.get(x, 0) for x in alt_segs - ref_segs)
            uniq_ref = sum(seg_len.get(x, 0) for x in ref_segs - alt_segs)
            # orientation-only difference (inversion): same segments, flipped
            flipped = sum(seg_len.get(sid, 0) for sid, o in t
                          if (sid, o) not in set(ref_trav) and sid in ref_segs)
            return max(abs(tlen(t) - ref_l), uniq_alt, uniq_ref, flipped)

        if max(event_size(t) for t, _ in alts) < min_alt_len:
            continue                                    # too small to be an SV

        start = offs[i] + seg_len.get(lsid, 0)
        end = offs[j]
        sites.append({
            "pos": start,
            "end": end,
            "ref_len": ref_l,
            "n_alleles": len(groups) + (0 if ref_trav in groups else 1),
            "ac_total": sum(len(hs) for _, hs in alts),
            "allele_lens": ",".join(str(tlen(t)) for t, _ in alts),
            "allele_acs": ",".join(str(len(hs)) for _, hs in alts),
            # carrier sets drive breakpoint merging (see merge_breakpoints)
            "_carriers": [frozenset(hs) for _, hs in alts],
            "_altlens": [tlen(t) for t, _ in alts],
        })
    return sites, ref_len_total, len(haps), len(anchors)


def merge_breakpoints(sites, max_gap):
    """Merge adjacent sites that are the two ends of one large event.

    Consecutive-anchor decomposition finds ONE bubble per event only when
    the event's interior contains no global anchors. For very large
    deletions in SV-dense windows that fails: enough haplotypes still
    traverse interior segments that those segments stay anchors, so the
    event is emitted as two small breakpoint sites (often with ref_len 0)
    with nothing spanning the interval between them.

    The two breakpoints of a single event are carried by exactly the same
    haplotypes, because they are the same mutation. Carrier-set identity is
    therefore the correct grouping key -- not proximity, which would merge
    genuinely distinct neighbouring SVs that happen to be close.

    Sites are merged when they share an identical carrier set and lie
    within max_gap bp. The merged record spans breakpoint to breakpoint,
    so its length is the true event length rather than the fragment
    lengths.
    """
    if not sites:
        return sites
    sites = sorted(sites, key=lambda s: s["pos"])
    used = [False] * len(sites)
    out = []
    for a in range(len(sites)):
        if used[a]:
            continue
        cur = sites[a]
        used[a] = True
        last = a
        for b in range(a + 1, len(sites)):
            if used[b]:
                continue
            if sites[b]["pos"] - cur["end"] > max_gap:
                break
            shared = set(cur["_carriers"]) & set(sites[b]["_carriers"])
            if not shared:
                continue
            cs = max(shared, key=len)

            # Carrier-set identity alone is NOT sufficient. With 20
            # haplotypes the genealogy has few distinct clades, so two
            # INDEPENDENT variants can arise on the same branch and share a
            # carrier set exactly; merging them would fuse two real events
            # into one call. Audited against ground truth, identity alone
            # gave a false fusion rate of ~2%.
            #
            # The two breakpoints of one deletion are distinguished by a
            # physical property: carriers have NO sequence between them, so
            # they cannot appear at any intervening polymorphic site. Two
            # independent variants do retain sequence in between and their
            # shared carriers will surface there. Require that.
            blocked = False
            for k in range(last + 1, b):
                if used[k]:
                    continue
                for car in sites[k]["_carriers"]:
                    if car & cs:
                        blocked = True
                        break
                if blocked:
                    break
            if blocked:
                continue
            last = b
            cur = {
                "pos": cur["pos"],
                "end": sites[b]["end"],
                "ref_len": sites[b]["end"] - cur["pos"],
                "n_alleles": max(cur["n_alleles"], sites[b]["n_alleles"]),
                "ac_total": len(cs),
                "allele_lens": str(sites[b]["end"] - cur["pos"]),
                "allele_acs": str(len(cs)),
                "_carriers": [cs],
                "_altlens": [sites[b]["end"] - cur["pos"]],
                "merged": True,
            }
            used[b] = True
        out.append(cur)
    return out


# ----------------------------------------------------------------------
# evaluation (optional)
# ----------------------------------------------------------------------

def evaluate(sites, truth_path, slack_frac=0.5, slack_min=200):
    rows = []
    with open(truth_path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        for line in fh:
            c = line.rstrip("\n").split("\t")
            s, e = int(c[ix["start"]]), int(c[ix["end"]])
            L, ac = int(c[ix["len"]]), int(c[ix["ac"]])
            slack = max(slack_min, int(slack_frac * L))
            over = [st for st in sites
                    if st["end"] > s - slack and st["pos"] < e + slack]
            hit = 0
            best = ""
            for st in over:
                accs = [int(x) for x in st["allele_acs"].split(",") if x]
                if ac in accs or st["ac_total"] == ac:
                    hit = 1
                    best = st["allele_acs"]
                    break
            rows.append({
                "sv_id": c[ix["sv_id"]], "type": c[ix["type"]], "len": L,
                "size_bin": c[ix["size_bin"]] if "size_bin" in ix else "",
                "ac_true": ac, "n_sites": len(over),
                "acs_found": best, "exact": hit,
            })
    return rows


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,)
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--version", action="version",
                    version=f"graphsv {__version__}")
    ap.add_argument("--ref-prefix", default="anc")
    ap.add_argument("--min-seg", type=int, default=50,
                    help="ignore segments below this length when comparing "
                         "traversals (removes linked SNP bubbles)")
    ap.add_argument("--min-alt-len", type=int, default=50,
                    help="minimum length difference to call a site")
    ap.add_argument("--max-merge-gap", type=int, default=50000,
                    help="merge adjacent sites sharing a carrier set within "
                         "this distance; 0 disables merging")
    ap.add_argument("--truth", help="optional truth_svs.tsv to evaluate against")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seg_len, paths = parse_gfa(args.gfa)
    if not paths:
        sys.exit("no P or W lines in GFA")
    ref_name = next((n for n in paths if n.startswith(args.ref_prefix)), None)
    if ref_name is None:
        sys.exit(f"no path with prefix '{args.ref_prefix}'; "
                 f"have: {sorted(paths)[:5]}")

    sites, ref_total, n_hap, n_anchor = call_sites(
        seg_len, paths, ref_name, args.min_seg, args.min_alt_len)
    n_raw = len(sites)
    if args.max_merge_gap > 0:
        sites = merge_breakpoints(sites, args.max_merge_gap)
        print(f"merged {n_raw} raw sites -> {len(sites)} events",
              file=sys.stderr)

    print(f"reference: {ref_name} ({ref_total:,} bp)", file=sys.stderr)
    print(f"haplotypes: {n_hap}   anchors: {n_anchor:,}   "
          f"segments: {len(seg_len):,}", file=sys.stderr)
    print(f"called {len(sites)} polymorphic sites", file=sys.stderr)

    cols = ["pos", "end", "ref_len", "n_alleles", "ac_total",
            "allele_lens", "allele_acs"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for st in sites:
            fh.write("\t".join(str(st[c]) for c in cols) + "\n")
    print(f"wrote {args.out}", file=sys.stderr)

    if args.truth:
        rows = evaluate(sites, args.truth)
        ev = args.out.replace(".tsv", "") + ".eval.tsv"
        ecols = ["sv_id", "type", "len", "size_bin", "ac_true",
                 "n_sites", "acs_found", "exact"]
        with open(ev, "w") as fh:
            fh.write("\t".join(ecols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in ecols) + "\n")
        tot = len(rows)
        ex = sum(r["exact"] for r in rows)
        print(f"\nexact AC: {ex}/{tot} = {ex/max(1,tot):.3f}", file=sys.stderr)
        bins = defaultdict(lambda: [0, 0])
        for r in rows:
            b = bins[r["size_bin"] or "NA"]
            b[0] += 1
            b[1] += r["exact"]
        for b in sorted(bins, key=lambda x: int(x) if str(x).isdigit() else 0):
            n, k = bins[b]
            print(f"  bin {b:>6}  n={n:3d}  exact={k/n:.3f}", file=sys.stderr)
        print(f"wrote {ev}", file=sys.stderr)


if __name__ == "__main__":
    main()
