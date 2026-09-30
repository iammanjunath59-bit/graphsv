#!/usr/bin/env python3
"""
fig3_hprc_locus.py

Figure 3. A common deletion in the HPRC Release 2 pangenome, as represented in
the raw bubble decomposition and in the released callset.

Panel a: alleles at the top level bubble at chr20:25,783,205, grouped by
alternate allele length. For each group the figure shows the number of
haplotype copies carried by the group as a whole and the largest count carried
by any single allele in it. The gap between the two is the quantity that a
reader of the allele level fields cannot recover.
Panel b: the haplotype copies of the deletion as partitioned across the
records of the released callset that carry it.

The script takes a small tab separated summary rather than querying the
pangenome, so the figure is reproducible offline and the numbers plotted are
exactly those reported. Produce that summary with:

    python3 fig3_hprc_locus.py --dump locus.tsv          # requires RAW and WAVE
    python3 fig3_hprc_locus.py --locus locus.tsv --out figures/Figure3

With --dump the script reads the released VCFs over HTTPS through their tabix
indices, using the environment variables RAW and WAVE, and writes the summary;
it performs no other analysis.

Expected contents of the summary, for checking:
    class   deletion   alleles  84   copies 216   max_single_allele_count 25
    class   reflike    alleles 132   copies 132   max_single_allele_count  1
    class   insertion  alleles 112   copies 113   max_single_allele_count  2
    record  25783207   copies 200
    record  25783899   copies  12
    record  25785257   copies   3
    record  25787417   copies   1
    AN 461, bubble REF length 68,466 bp, event size 61,963 bp

No legend is drawn inside the axes; a single legend is placed below both
panels.
"""

import argparse
import os
import re
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

POS = 25783205
CHROM = "chr20"
GROUP_COLOUR = "#1b4965"
MAX_COLOUR = "#c1666b"
RECORD_COLOUR = "#5fa8d3"


# ---------------------------------------------------------------- dump mode
def bcf(fmt, region, vcf):
    p = subprocess.run(["bcftools", "query", "-r", region, "-f", fmt, vcf],
                       capture_output=True, text=True)
    if p.returncode != 0:
        tail = p.stderr.strip().splitlines()[-1:] or ["unknown"]
        sys.exit(f"bcftools failed on {region}: {tail[0]}")
    return [l for l in p.stdout.rstrip("\n").split("\n") if l]


def dump(path):
    raw = os.environ.get("RAW")
    wave = os.environ.get("WAVE")
    if not raw or not wave:
        sys.exit("set RAW and WAVE to the released raw and wave VCF URLs")

    row = bcf("%POS\\t%REF\\t%ALT\\t%INFO/AT\\t%INFO/AC\\t%INFO/AN\\t[%SAMPLE=%GT;]\\n",
              f"{CHROM}:{POS}-{POS}", raw)
    row = next((r for r in row if r.split("\t")[0] == str(POS)), None)
    if row is None:
        sys.exit(f"no raw record at {CHROM}:{POS}")
    _, ref, alt, at, ac, an, gts = row.split("\t")[:7]
    alts = alt.split(",")
    ats = at.split(",")
    acs = [int(x) for x in ac.split(",")]

    # length classes, as used in the text: the three modes of the ALT length
    # distribution at this bubble
    def cls(n):
        if n < 20000:
            return "deletion"
        if n < 73000:
            return "reflike"
        return "insertion"

    groups = {}
    for i, a in enumerate(alts):
        g = groups.setdefault(cls(len(a)), {"alleles": 0, "copies": 0, "max": 0,
                                            "idx": []})
        g["alleles"] += 1
        g["copies"] += acs[i]
        g["max"] = max(g["max"], acs[i])
        g["idx"].append(i + 1)

    # carriers of the deletion class, as haplotype copies
    want = set(groups["deletion"]["idx"])
    carriers = set()
    for chunk in gts.strip(";").split(";"):
        if "=" not in chunk:
            continue
        s, gt = chunk.split("=", 1)
        for slot, a in enumerate(gt.replace("/", "|").split("|")):
            if a.isdigit() and int(a) in want:
                carriers.add((s, slot))

    import statistics
    event = int(statistics.median(len(ref) - len(alts[i - 1])
                                  for i in groups["deletion"]["idx"]))

    # wave records in the bubble span that carry those haplotypes
    end = POS + len(ref)
    recs = []
    for line in bcf("%POS\\t%REF\\t%ALT\\t[%SAMPLE=%GT;]\\n",
                    f"{CHROM}:{POS}-{end}", wave):
        p_, wref, walt, wgts = line.split("\t")[:4]
        keep = {j + 1 for j, a in enumerate(walt.split(","))
                if abs((len(wref) - len(a)) - event) <= max(200, 0.2 * event)}
        if not keep:
            continue
        hit = set()
        for chunk in wgts.strip(";").split(";"):
            if "=" not in chunk:
                continue
            s, gt = chunk.split("=", 1)
            for slot, a in enumerate(gt.replace("/", "|").split("|")):
                if a.isdigit() and int(a) in keep and (s, slot) in carriers:
                    hit.add((s, slot))
        if hit:
            recs.append((int(p_), len(hit)))
    recs.sort(key=lambda t: -t[1])

    with open(path, "w") as fh:
        fh.write("kind\tlabel\talleles\tcopies\tmax_single\n")
        for name in ("deletion", "reflike", "insertion"):
            g = groups.get(name)
            if g:
                fh.write(f"class\t{name}\t{g['alleles']}\t{g['copies']}\t{g['max']}\n")
        for p_, n in recs:
            fh.write(f"record\t{p_}\t\t{n}\t\n")
        fh.write(f"meta\tAN\t\t{an}\t\n")
        fh.write(f"meta\tref_len\t\t{len(ref)}\t\n")
        fh.write(f"meta\tevent_size\t\t{event}\t\n")
    print("wrote", path)


# ---------------------------------------------------------------- plot mode
def read_locus(path):
    classes, records, meta = [], [], {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "class":
                classes.append((f[1], int(f[2]), int(f[3]), int(f[4])))
            elif f[0] == "record":
                records.append((int(f[1]), int(f[3])))
            elif f[0] == "meta":
                meta[f[1]] = int(f[3])
    return classes, records, meta


LABEL = {"deletion": "deletion allele\n(~6.5 kb)",
         "reflike": "reference-like\n(~68.5 kb)",
         "insertion": "insertion allele\n(~79 kb)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dump", help="write the locus summary and exit")
    ap.add_argument("--locus", help="locus summary written by --dump")
    ap.add_argument("--out", default="Figure3")
    args = ap.parse_args()

    if args.dump:
        dump(args.dump)
        return
    if not args.locus:
        sys.exit("give --locus (or --dump to produce it)")

    classes, records, meta = read_locus(args.locus)
    an = meta.get("AN")

    fig = plt.figure(figsize=(9.4, 4.1))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0],
                          left=0.095, right=0.985, top=0.90, bottom=0.30,
                          wspace=0.28)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])

    x = list(range(len(classes)))
    copies = [c[2] for c in classes]
    maxima = [c[3] for c in classes]
    b1 = axa.bar([i - 0.19 for i in x], copies, width=0.36,
                 color=GROUP_COLOUR, edgecolor="white", lw=0.8)
    b2 = axa.bar([i + 0.19 for i in x], maxima, width=0.36,
                 color=MAX_COLOUR, edgecolor="white", lw=0.8)
    axa.set_xticks(x)
    axa.set_xticklabels([LABEL.get(c[0], c[0]) for c in classes], fontsize=9)
    axa.set_ylabel("haplotype copies")
    axa.grid(alpha=0.25, lw=0.5, axis="y")
    axa.set_axisbelow(True)
    if an:
        axa.set_ylim(0, an * 0.55)
    axa.set_title("a   Alleles at the bubble, grouped by length",
                  loc="left", fontsize=10)

    xs = list(range(len(records)))
    b3 = axb.bar(xs, [r[1] for r in records], width=0.55,
                 color=RECORD_COLOUR, edgecolor="white", lw=0.8)
    total = sum(r[1] for r in records)
    l1 = axb.axhline(total, color=MAX_COLOUR, lw=1.4, ls="--")
    axb.set_xticks(xs)
    axb.set_xticklabels([f"{r[0]:,}" for r in records], fontsize=8, rotation=30,
                        ha="right")
    axb.set_xlabel("position of record in the released callset")
    axb.set_ylabel("haplotype copies")
    axb.set_ylim(0, total * 1.18)
    axb.grid(alpha=0.25, lw=0.5, axis="y")
    axb.set_axisbelow(True)
    axb.set_title("b   The same deletion in the released callset",
                  loc="left", fontsize=10)

    fig.legend([b1, b2, b3, l1],
               ["copies carried by the allele group",
                "largest count on any single allele",
                "copies carried by one record",
                f"total copies of the deletion ({total})"],
               loc="lower center", ncol=2, frameon=False, fontsize=9,
               bbox_to_anchor=(0.5, 0.02), columnspacing=2.4)

    out = os.path.splitext(args.out)[0]
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out}.{ext}", dpi=400)
        print("wrote", f"{out}.{ext}")

    print("\nvalues plotted:")
    for name, nall, cop, mx in classes:
        print(f"  {name:10} alleles {nall:4d}  copies {cop:4d}  "
              f"largest single allele count {mx}")
    for p_, n in records:
        print(f"  record {p_:,}  copies {n}")
    print(f"  records sum to {total}; AN at the bubble is {an}")


if __name__ == "__main__":
    main()
