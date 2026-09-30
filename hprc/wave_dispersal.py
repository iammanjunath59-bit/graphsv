#!/usr/bin/env python3
"""
wave_dispersal.py

An event can be missing from the wave VCF in two very different senses:

  (a) the sequence difference is gone, or
  (b) the sequence difference is present but scattered across many small
      records, so no single record represents the structural event.

These have different implications. (b) means SV callers with a size filter
miss the event while sequence-level analyses still see it; (a) would mean
outright loss. This script distinguishes them.

Method: take one deletion event at a raw bubble, identify its carrier
haplotypes from the AT traversal field, then sum deleted base pairs per
haplotype across every wave record in the bubble span. If carriers sum to
roughly the event size and non-carriers do not, the event is dispersed,
not lost.

Usage:
    python3 wave_dispersal.py --pos 28948271 --event 111822236_111820427
    python3 wave_dispersal.py --pos 25783205 --event 111144578_112737255
"""

import argparse
import os
import re
import statistics
import subprocess
import sys

NODE_RE = re.compile(r"[<>](\d+)")

B = ("https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/"
     "freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/")
RAW = os.environ.get("RAW") or B + "hprc-v2.1-mc-grch38.raw.vcf.gz"
WAVE = os.environ.get("WAVE") or B + "hprc-v2.1-mc-grch38.wave.vcf.gz"


def bcf(fmt, region, vcf):
    p = subprocess.run(["bcftools", "query", "-r", region, "-f", fmt, vcf],
                       capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"bcftools failed on {region}: {p.stderr.strip().splitlines()[-1:]}")
    return [ln for ln in p.stdout.rstrip("\n").split("\n") if ln]


def nodes(s):
    return NODE_RE.findall(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr20")
    ap.add_argument("--pos", type=int, required=True)
    ap.add_argument("--event", required=True, help="left_right node pair, as in the sweep table")
    ap.add_argument("--min-event", type=int, default=1000)
    args = ap.parse_args()

    left, right = args.event.split("_")

    # ---- raw bubble: find the event's alleles and its carrier haplotypes ----
    rows = bcf("%POS\\t%REF\\t%ALT\\t%INFO/AT\\t%INFO/AC\\t%INFO/AN\\t[%SAMPLE=%GT;]\\n",
               f"{args.chrom}:{args.pos}-{args.pos}", RAW)
    row = next((r for r in rows if r.split("\t")[0] == str(args.pos)), None)
    if row is None:
        sys.exit(f"no raw record at {args.chrom}:{args.pos}")
    _, ref, alt, at, ac, an, gts = row.split("\t")[:7]
    alts, ats = alt.split(","), at.split(",")
    acs = [int(x) for x in ac.split(",")]

    refn = nodes(ats[0])
    ref_index = {}
    for i, n in enumerate(refn):
        ref_index.setdefault(n, i)

    ev_alleles, sizes = [], []
    for i, a in enumerate(alts):
        deleted = len(ref) - len(a)
        if deleted < args.min_event or i + 1 >= len(ats):
            continue
        kept = sorted({ref_index[n] for n in nodes(ats[i + 1]) if n in ref_index})
        if len(kept) < 2:
            continue
        gap, l, r = 0, None, None
        for p1, p2 in zip(kept, kept[1:]):
            if p2 - p1 > gap:
                gap, l, r = p2 - p1, refn[p1], refn[p2]
        if (l, r) == (left, right):
            ev_alleles.append(i + 1)
            sizes.append(deleted)

    if not ev_alleles:
        sys.exit(f"event {args.event} not found at {args.chrom}:{args.pos}")
    event_size = int(statistics.median(sizes))

    carriers = set()
    for chunk in gts.strip(";").split(";"):
        if "=" not in chunk:
            continue
        smp, gt = chunk.split("=", 1)
        for slot, a in enumerate(gt.replace("/", "|").split("|")):
            if a.isdigit() and int(a) in ev_alleles:
                carriers.add((smp, slot))

    span_end = args.pos + len(ref)
    print(f"event {args.event} at {args.chrom}:{args.pos}")
    print(f"  size {event_size:,} bp | {len(ev_alleles)} raw alleles | "
          f"{len(carriers)} carrier haplotypes | AN {an}")
    print(f"  bubble span {args.chrom}:{args.pos}-{span_end} ({len(ref):,} bp)")

    # ---- wave: sum deleted bp per haplotype across all records in the span ----
    wrows = bcf("%POS\\t%REF\\t%ALT\\t[%SAMPLE=%GT;]\\n",
                f"{args.chrom}:{args.pos}-{span_end}", WAVE)
    deleted_bp, present, biggest = {}, set(), 0
    for line in wrows:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        _, wref, walt, wgts = parts[:4]
        alt_del = {j + 1: len(wref) - len(a) for j, a in enumerate(walt.split(","))}
        biggest = max(biggest, max(alt_del.values(), default=0))
        for chunk in wgts.strip(";").split(";"):
            if "=" not in chunk:
                continue
            smp, gt = chunk.split("=", 1)
            for slot, a in enumerate(gt.replace("/", "|").split("|")):
                present.add((smp, slot))
                if a.isdigit() and int(a) in alt_del and alt_del[int(a)] > 0:
                    deleted_bp[(smp, slot)] = deleted_bp.get((smp, slot), 0) + alt_del[int(a)]

    print(f"  wave: {len(wrows):,} records in span, largest single deletion {biggest:,} bp")

    car = [deleted_bp.get(h, 0) for h in carriers if h in present]
    non = [deleted_bp.get(h, 0) for h in present - carriers]
    if not car or not non:
        sys.exit("  not enough haplotypes on both sides to compare")

    def q(v):
        v = sorted(v)
        return v[len(v) // 2], v[0], v[-1]

    cm, clo, chi = q(car)
    nm, nlo, nhi = q(non)
    print(f"\n  summed deleted bp per haplotype across all wave records:")
    print(f"    carriers     (n={len(car):3d})  median {cm:>8,}  range {clo:,}-{chi:,}")
    print(f"    non-carriers (n={len(non):3d})  median {nm:>8,}  range {nlo:,}-{nhi:,}")
    print(f"    event size                    {event_size:>8,}")
    print(f"    carrier median / event size   {cm / event_size:.3f}")
    print(f"    excess over non-carriers      {cm - nm:>8,} bp "
          f"({(cm - nm) / event_size:.3f} of event)")
    print("\n  Interpretation: a carrier excess near the event size means the")
    print("  sequence survives, dispersed across small records; near zero means")
    print("  the difference is absent from the wave VCF entirely.")


if __name__ == "__main__":
    main()
