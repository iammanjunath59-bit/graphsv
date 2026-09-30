#!/usr/bin/env python3
"""
stage_recovery.py

Supplementary Fig. 1 analysis: where in the VCF conversion chain does the
allele count of a structural variant stop being recoverable?

For every simulated SV with known truth (size, allele count, carrier
haplotypes) this scores three VCFs produced from the same graph:

    raw    vg deconstruct -a          (bubble traversals, one allele per path)
    bub    + vcfbub -l 0 -a 100000    (bubble filtering / flattening)
    wave   + vcfwave -I 1000          (realignment into atomic records)

Scoring is by carrier set, not by record count, because a record of the right
size that names the wrong haplotypes has not recovered the variant:

    EXACT_SINGLE   one record's carrier set equals the truth carrier set
    FRAGMENTED     several records together cover it, none alone
    PARTIAL        records found, union is a strict subset of truth
    INFLATED       union includes haplotypes that do not carry the SV
    ABSENT         no size-compatible record anywhere in the interval

Coordinate systems are reconciled automatically: the offset between truth
coordinates and VCF POS is estimated once per file from large SVs with
unambiguous sizes, then applied throughout, and reported so it can be checked.

Usage:
    python3 stage_recovery.py --dir ~/pangenome/figs1 \\
        --truth-dir ~/pangenome/sweep54d --out stage_recovery.tsv
"""

import argparse
import collections
import glob
import os
import sys

STAGES = ("raw", "bub", "wave")


def read_truth(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(hdr):
                continue
            carriers = {int(x) for x in f[ix["carriers"]].split(",") if x != ""}
            rows.append({
                "sv_id": f[ix["sv_id"]],
                "type": f[ix["type"]],
                "start": int(f[ix["start"]]),
                "end": int(f[ix["end"]]),
                "len": int(f[ix["len"]]),
                "ac": int(f[ix["ac"]]),
                "size_bin": f[ix["size_bin"]],
                "synthetic": f[ix["synthetic"]],
                "carriers": carriers,
            })
    return rows


def read_vcf(path):
    """-> (records, sample_names). Each record: pos, size_change, carrier index set."""
    recs, samples = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            pos, ref, alt = int(f[1]), f[3], f[4]
            gts = f[9:]
            alts = alt.split(",")
            # carriers per ALT index, and that ALT's size change vs REF
            per_alt = collections.defaultdict(set)
            for si, g in enumerate(gts):
                call = g.split(":")[0]
                for a in call.replace("/", "|").split("|"):
                    if a.isdigit() and int(a) > 0:
                        per_alt[int(a)].add(si)
            for j, a in enumerate(alts, start=1):
                if j not in per_alt:
                    continue
                recs.append({
                    "pos": pos,
                    "size": abs(len(ref) - len(a)),
                    "carriers": per_alt[j],
                })
    return recs, samples


def estimate_offset(truth, recs, min_len=5000, tol=0.02):
    """
    Offset between truth 'start' and VCF POS, from large SVs whose size is
    distinctive. Returns (offset, n_supporting). Falls back to 0.
    """
    votes = collections.Counter()
    big = [t for t in truth if t["len"] >= min_len]
    by_size = collections.defaultdict(list)
    for r in recs:
        by_size[r["size"]].append(r)
    sizes = sorted(by_size)
    for t in big:
        lo, hi = t["len"] * (1 - tol), t["len"] * (1 + tol)
        for s in sizes:
            if s < lo:
                continue
            if s > hi:
                break
            for r in by_size[s]:
                votes[r["pos"] - t["start"]] += 1
    if not votes:
        return 0, 0
    off, n = votes.most_common(1)[0]
    return off, n


def score(t, recs_by_pos, offset, pos_tol, size_tol_frac, size_tol_min,
          name_to_idx):
    want = t["start"] + offset
    stol = max(size_tol_min, size_tol_frac * t["len"])
    hits = []
    for r in recs_by_pos:
        if abs(r["pos"] - want) > pos_tol:
            continue
        if abs(r["size"] - t["len"]) > stol:
            continue
        hits.append(r)

    truth_c = {name_to_idx[i] for i in t["carriers"] if i in name_to_idx}
    if not hits:
        return "ABSENT", 0, 0.0, 0.0, 0

    exact = [r for r in hits if r["carriers"] == truth_c]
    union = set().union(*[r["carriers"] for r in hits])
    recall = len(union & truth_c) / len(truth_c) if truth_c else 0.0
    prec = len(union & truth_c) / len(union) if union else 0.0
    best = max(len(r["carriers"] & truth_c) for r in hits)

    if exact:
        verdict = "EXACT_SINGLE"
    elif union >= truth_c and len(hits) > 1:
        verdict = "FRAGMENTED"
    elif union - truth_c:
        verdict = "INFLATED"
    else:
        verdict = "PARTIAL"
    return verdict, len(hits), round(recall, 4), round(prec, 4), best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, help="directory with *.raw/bub/wave.vcf")
    ap.add_argument("--truth-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pos-tol", type=int, default=200)
    ap.add_argument("--size-tol-frac", type=float, default=0.05)
    ap.add_argument("--size-tol-min", type=int, default=20)
    ap.add_argument("--min-len", type=int, default=50, help="skip SVs below this")
    ap.add_argument("--types", default="DEL",
                    help="comma-separated SV types to score. Size matching uses "
                         "|len(REF)-len(ALT)|, which is meaningful for DEL but "
                         "not for INV (no length change) or DUP (length change "
                         "does not equal the duplicated span), so those are "
                         "excluded by default rather than scored as absent. "
                         "Use ALL to include everything.")
    args = ap.parse_args()

    reps = sorted({os.path.basename(p).rsplit(".", 2)[0]
                   for p in glob.glob(os.path.join(args.dir, "*.raw.vcf"))})
    if not reps:
        sys.exit(f"no *.raw.vcf in {args.dir}")
    sys.stderr.write(f"{len(reps)} replicates\n")

    out = open(args.out, "w")
    out.write("\t".join(["replicate", "sv_id", "type", "size_bin", "len", "ac",
                         "stage", "verdict", "n_records", "recall",
                         "precision", "best_overlap"]) + "\n")

    tally = collections.Counter()
    for rep in reps:
        tpath = os.path.join(args.truth_dir, rep + ".truth_svs.tsv")
        if not os.path.exists(tpath):
            sys.stderr.write(f"  {rep}: no truth file, skipped\n")
            continue
        truth = [t for t in read_truth(tpath) if t["len"] >= args.min_len]
        if args.types.upper() != "ALL":
            keep = {x.strip().upper() for x in args.types.split(",")}
            truth = [t for t in truth if t["type"].upper() in keep]
        for stage in STAGES:
            vpath = os.path.join(args.dir, f"{rep}.{stage}.vcf")
            if not os.path.exists(vpath):
                continue
            recs, samples = read_vcf(vpath)
            # truth carriers are haplotype indices; VCF columns are sim00..simNN
            name_to_idx = {}
            for si, s in enumerate(samples):
                digits = "".join(ch for ch in s.split("#")[0] if ch.isdigit())
                if digits:
                    name_to_idx[int(digits)] = si
            off, nsup = estimate_offset(truth, recs)
            if stage == "raw":
                sys.stderr.write(f"  {rep}: offset {off} (n={nsup}), "
                                 f"{len(samples)} samples, {len(recs)} alt-records\n")
            for t in truth:
                v, n, rec, pre, best = score(
                    t, recs, off, args.pos_tol, args.size_tol_frac,
                    args.size_tol_min, name_to_idx)
                out.write("\t".join(str(x) for x in [
                    rep, t["sv_id"], t["type"], t["size_bin"], t["len"], t["ac"],
                    stage, v, n, rec, pre, best]) + "\n")
                tally[(stage, v)] += 1
    out.close()

    print(f"\nwrote {args.out}\n")
    verdicts = sorted({v for _, v in tally})
    print(f"{'stage':6} " + " ".join(f"{v:>13}" for v in verdicts))
    for stage in STAGES:
        tot = sum(tally[(stage, v)] for v in verdicts)
        if not tot:
            continue
        cells = " ".join(f"{tally[(stage, v)]:>13}" for v in verdicts)
        print(f"{stage:6} {cells}   (n={tot})")
    print()
    for stage in STAGES:
        tot = sum(tally[(stage, v)] for v in verdicts)
        if tot:
            ex = tally[(stage, "EXACT_SINGLE")]
            print(f"  {stage:5} exact single-record recovery: {ex}/{tot} = {ex/tot:.3f}")


if __name__ == "__main__":
    main()
