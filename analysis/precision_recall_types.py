#!/usr/bin/env python3
"""
precision_recall_types.py

Scores a callset against the truth in three ways:

  1. Report precision, not only recall. The v1.0.0 evaluation scored each
     true variant against the callset and never asked the reverse question,
     namely what fraction of the structural variant records in the callset
     correspond to a true variant.

  2. Tighten the matching criterion. The v1.0.0 criterion scored a variant
     as recovered when its true allele count appeared among the alternate
     allele counts at any overlapping site, without requiring the matched
     allele to be the same event or even the same type. This script adds
     successively stricter criteria and reports all of them side by side, so
     the effect of the criterion is itself visible.

  3. Identify which graph structures are associated with failure, by
     recording the properties of the matched site (allele count, span,
     allele number) alongside the outcome.

Matching proceeds in nested stages, each a subset of the previous one:

  overlap      a record overlaps the true interval within max(200 bp, L/2)
  size         one of its alternate alleles changes length by about L
  type         that alternate's implied type agrees with the true type
  count        the true allele count appears among the alternates
  count_typed  the true allele count appears among type compatible alternates
  carriers     a type compatible alternate's carrier set equals the truth set

Type is inferred from the record itself, not from any INFO field: an alternate
shorter than the reference by at least the minimum size is a deletion, one
longer is an insertion, and a same length alternate that matches the reverse
complement of the reference above a sequence identity threshold is an
inversion. Simulated duplications are tandem, so they present as insertions
and are matched against the insertion class; this is stated in the output
rather than hidden, because it is a limitation of what a VCF record can
express.

Usage, one replicate first:

    python3 precision_recall_types.py \\
        --truth results/sweep54d/chr20_10000000_s1.truth_svs.tsv \\
        --vcf   results/sweep54d/chr20_10000000_s1.wave.vcf \\
        --out   out/pr_one

then across a sweep:

    python3 precision_recall_types.py \\
        --sweep results/sweep54d --vcf-suffix .wave.vcf --out out/pr_wave

Outputs <out>.perSV.tsv, <out>.percall.tsv and <out>.summary.txt
"""

import argparse
import glob
import os
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def identity(a, b):
    """Fraction of matching positions over the shorter sequence. Cheap and
    sufficient to tell a reverse complement from an unrelated allele."""
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    same = sum(1 for i in range(n) if a[i] == b[i])
    return same / n


def read_truth(path):
    rows = []
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(hdr)}
        for need in ("sv_id", "type", "start", "end", "len", "ac", "carriers"):
            if need not in ix:
                sys.exit(f"{path} lacks column {need}; found {hdr}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ix.values()):
                continue
            rows.append({
                "sv_id": f[ix["sv_id"]],
                "type": f[ix["type"]].upper(),
                "start": int(f[ix["start"]]),
                "end": int(f[ix["end"]]),
                "len": int(f[ix["len"]]),
                "ac": int(f[ix["ac"]]),
                "size_bin": f[ix["size_bin"]] if "size_bin" in ix else "",
                "carriers": {int(x) for x in f[ix["carriers"]].split(",") if x != ""},
            })
    return rows


def read_vcf(path, min_size):
    """-> (alts, samples). One entry per ALT allele that changes length by at
    least min_size or looks like an inversion."""
    alts, samples = [], []
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
            pos, ref = int(f[1]), f[3]
            gts = f[9:]
            per_alt = {}
            for si, g in enumerate(gts):
                for a in g.split(":")[0].replace("/", "|").split("|"):
                    if a.isdigit() and int(a) > 0:
                        per_alt.setdefault(int(a), set()).add(si)
            for j, alt in enumerate(f[4].split(","), start=1):
                d = len(alt) - len(ref)
                if d <= -min_size:
                    kind, size = "DEL", -d
                elif d >= min_size:
                    kind, size = "INS", d
                elif len(ref) >= min_size and identity(alt, revcomp(ref)) >= 0.90:
                    kind, size = "INV", len(ref)
                else:
                    continue
                alts.append({
                    "pos": pos, "end": pos + len(ref), "alt_index": j,
                    "kind": kind, "size": size,
                    "ac": len(per_alt.get(j, ())),
                    "carriers": per_alt.get(j, set()),
                })
    return alts, samples


# truth type -> the record types that can legitimately express it
TYPE_MAP = {"DEL": {"DEL"}, "DUP": {"INS"}, "INS": {"INS"}, "INV": {"INV"}}


def score(truth, alts, name_to_idx, size_tol_frac, size_tol_min, pos_tol_min):
    out = []
    for t in truth:
        ptol = max(pos_tol_min, t["len"] // 2)
        stol = max(size_tol_min, size_tol_frac * t["len"])
        want = TYPE_MAP.get(t["type"], set())
        truth_c = {name_to_idx[i] for i in t["carriers"] if i in name_to_idx}

        overlap = [a for a in alts
                   if a["end"] > t["start"] - ptol and a["pos"] < t["end"] + ptol]
        sized = [a for a in overlap if abs(a["size"] - t["len"]) <= stol]
        typed = [a for a in sized if a["kind"] in want]

        best_j = 0.0
        exact_carrier = 0
        for a in typed:
            inter = len(a["carriers"] & truth_c)
            union = len(a["carriers"] | truth_c)
            if union:
                best_j = max(best_j, inter / union)
            if a["carriers"] == truth_c and truth_c:
                exact_carrier = 1

        out.append({
            "sv_id": t["sv_id"], "type": t["type"], "len": t["len"],
            "size_bin": t["size_bin"], "ac_true": t["ac"],
            "n_overlap": len(overlap), "n_sized": len(sized), "n_typed": len(typed),
            "alleles_at_site": max((a["ac"] for a in overlap), default=0),
            "hit_overlap": int(bool(overlap)),
            "hit_size": int(bool(sized)),
            "hit_type": int(bool(typed)),
            "ac_in_overlap": int(any(a["ac"] == t["ac"] for a in overlap)),
            "ac_in_typed": int(any(a["ac"] == t["ac"] for a in typed)),
            "carrier_exact": exact_carrier,
            "carrier_best_jaccard": round(best_j, 4),
        })
    return out


def precision(truth, alts, size_tol_frac, size_tol_min, pos_tol_min):
    """For each ALT in the callset, does a true variant of compatible type and
    size overlap it?"""
    rows = []
    for a in alts:
        ok = 0
        for t in truth:
            if a["kind"] not in TYPE_MAP.get(t["type"], set()):
                continue
            ptol = max(pos_tol_min, t["len"] // 2)
            stol = max(size_tol_min, size_tol_frac * t["len"])
            if (a["end"] > t["start"] - ptol and a["pos"] < t["end"] + ptol
                    and abs(a["size"] - t["len"]) <= stol):
                ok = 1
                break
        rows.append({"pos": a["pos"], "alt_index": a["alt_index"],
                     "kind": a["kind"], "size": a["size"], "ac": a["ac"],
                     "matches_truth": ok})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth")
    ap.add_argument("--vcf")
    ap.add_argument("--sweep", help="directory of replicates, used with --vcf-suffix")
    ap.add_argument("--vcf-suffix", default=".wave.vcf")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-size", type=int, default=50)
    ap.add_argument("--size-tol-frac", type=float, default=0.20)
    ap.add_argument("--size-tol-min", type=int, default=50)
    ap.add_argument("--pos-tol-min", type=int, default=200)
    args = ap.parse_args()

    pairs = []
    if args.sweep:
        for t in sorted(glob.glob(os.path.join(args.sweep, "*.truth_svs.tsv"))):
            base = t[: -len(".truth_svs.tsv")]
            v = base + args.vcf_suffix
            if os.path.exists(v):
                pairs.append((t, v))
        if not pairs:
            sys.exit(f"no truth/{args.vcf_suffix} pairs in {args.sweep}")
    else:
        if not (args.truth and args.vcf):
            sys.exit("give --truth and --vcf, or --sweep")
        pairs = [(args.truth, args.vcf)]
    sys.stderr.write(f"{len(pairs)} replicate(s)\n")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    sv_fields = ["replicate", "sv_id", "type", "len", "size_bin", "ac_true",
                 "n_overlap", "n_sized", "n_typed", "alleles_at_site",
                 "hit_overlap", "hit_size", "hit_type",
                 "ac_in_overlap", "ac_in_typed",
                 "carrier_exact", "carrier_best_jaccard"]
    call_fields = ["replicate", "pos", "alt_index", "kind", "size", "ac",
                   "matches_truth"]

    all_sv, all_call = [], []
    with open(args.out + ".perSV.tsv", "w") as fsv, \
         open(args.out + ".percall.tsv", "w") as fc:
        fsv.write("\t".join(sv_fields) + "\n")
        fc.write("\t".join(call_fields) + "\n")
        for tpath, vpath in pairs:
            rep = os.path.basename(tpath)[: -len(".truth_svs.tsv")]
            truth = read_truth(tpath)
            alts, samples = read_vcf(vpath, args.min_size)
            name_to_idx = {}
            for si, s in enumerate(samples):
                digits = "".join(ch for ch in s.split("#")[0] if ch.isdigit())
                if digits:
                    name_to_idx[int(digits)] = si
            rows = score(truth, alts, name_to_idx, args.size_tol_frac,
                         args.size_tol_min, args.pos_tol_min)
            calls = precision(truth, alts, args.size_tol_frac,
                              args.size_tol_min, args.pos_tol_min)
            for r in rows:
                r["replicate"] = rep
                fsv.write("\t".join(str(r[k]) for k in sv_fields) + "\n")
            for r in calls:
                r["replicate"] = rep
                fc.write("\t".join(str(r[k]) for k in call_fields) + "\n")
            all_sv.extend(rows)
            all_call.extend(calls)

    # ---------------------------------------------------------------- report
    def frac(rows, key):
        return (sum(r[key] for r in rows) / len(rows)) if rows else float("nan")

    bins = sorted({r["size_bin"] for r in all_sv if str(r["size_bin"]).isdigit()},
                  key=int)
    lines = []
    lines.append(f"replicates: {len(pairs)}   true variants: {len(all_sv)}   "
                 f"callset alleles scored: {len(all_call)}")
    lines.append("")
    lines.append("RECALL under successively stricter criteria, all variants")
    for key, label in (("hit_overlap", "a record overlaps"),
                       ("hit_size", "and its size matches"),
                       ("hit_type", "and its type matches"),
                       ("ac_in_overlap", "true count among overlapping alternates"),
                       ("ac_in_typed", "true count among type matched alternates"),
                       ("carrier_exact", "a type matched alternate has the exact carrier set")):
        lines.append(f"  {label:58} {frac(all_sv, key):.4f}")
    lines.append("")
    lines.append("PRECISION of the callset")
    lines.append(f"  callset alleles of at least {args.min_size} bp: {len(all_call)}")
    lines.append(f"  matching a true variant of compatible type and size: "
                 f"{frac(all_call, 'matches_truth'):.4f}")
    by_kind = {}
    for r in all_call:
        by_kind.setdefault(r["kind"], []).append(r)
    for k in sorted(by_kind):
        lines.append(f"    {k:4} n={len(by_kind[k]):6d}  "
                     f"precision {frac(by_kind[k], 'matches_truth'):.4f}")
    lines.append("")
    lines.append("RECALL by size class: overlap / size / type / count / count typed / carriers")
    for b in bins:
        s = [r for r in all_sv if str(r["size_bin"]) == b]
        lines.append(f"  {b:>7} n={len(s):5d}  " + "  ".join(
            f"{frac(s, k):.3f}" for k in ("hit_overlap", "hit_size", "hit_type",
                                          "ac_in_overlap", "ac_in_typed",
                                          "carrier_exact")))
    lines.append("")
    lines.append("RECALL by true variant type")
    for ty in sorted({r["type"] for r in all_sv}):
        s = [r for r in all_sv if r["type"] == ty]
        lines.append(f"  {ty:4} n={len(s):5d}  " + "  ".join(
            f"{frac(s, k):.3f}" for k in ("hit_overlap", "hit_size", "hit_type",
                                          "ac_in_overlap", "ac_in_typed",
                                          "carrier_exact")))
    lines.append("")
    lines.append("FAILURE against the number of alleles at the matched site")
    buckets = [(1, 2), (3, 4), (5, 9), (10, 29), (30, 10**9)]
    for lo, hi in buckets:
        s = [r for r in all_sv if lo <= r["alleles_at_site"] <= hi]
        if not s:
            continue
        label = f"{lo} to {hi}" if hi < 10**9 else f"{lo} or more"
        lines.append(f"  alleles {label:>12}  n={len(s):5d}  "
                     f"count {frac(s, 'ac_in_overlap'):.3f}  "
                     f"carriers {frac(s, 'carrier_exact'):.3f}")
    lines.append("")
    lines.append("Note: simulated duplications are tandem and appear in the callset as "
                 "insertions, so DUP is matched against the insertion class. A VCF "
                 "record cannot distinguish a tandem duplication from any other "
                 "insertion of the same length without inspecting the sequence.")

    report = "\n".join(lines)
    with open(args.out + ".summary.txt", "w") as fh:
        fh.write(report + "\n")
    print(report)
    print(f"\nwrote {args.out}.perSV.tsv, .percall.tsv, .summary.txt")


if __name__ == "__main__":
    main()
