#!/usr/bin/env python3
"""
merge_fix.py

Tests a correction that uses the VCF alone, for tools that cannot read the
AT field.

The partition test showed that most failed duplications are not lost: their
carriers are split across several biallelic records whose carrier sets are
disjoint and whose counts sum to the true count. This script asks whether
merging such records restores the correct count, and how often merging
wrongly joins two genuinely distinct variants.

Merge rule, applied to every callset alternate of at least --min-size bp:

  same kind      DEL with DEL, INS with INS, INV with INV
  similar size   within max(--size-tol-min, --size-tol-frac x size) of the
                 size of the first record in the cluster
  nearby         start within max(--pos-tol-min, size) of the cluster span
  disjoint       no haplotype carries both; a haplotype carries one allele

A merged allele takes the union of carriers, so its count is the sum of the
member counts. Nothing about the truth is used to decide a merge; the truth
is used only to score the result.

Scoring reuses precision_recall_types.py so that the before numbers are the
same numbers already reported, and the after numbers are directly comparable.

False fusion: a cluster whose members are each attributable to a true
variant (compatible type and size, carriers a subset of that variant's
carriers) but not all to the same one.

Usage:

    python3 analysis/merge_fix.py --sweep ~/pangenome_rebuild/sweep54d \\
        --vcf-suffix .wave.vcf --out out/merge_wave

Outputs <out>.perSV.tsv, <out>.clusters.tsv and <out>.summary.txt
"""

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import precision_recall_types as P  # noqa: E402


def kmers(seq, k=21):
    seq = seq.upper()
    return {hash(seq[i:i + k]) for i in range(len(seq) - k + 1)}


def kmer_sim(a, b):
    """Shared k-mers over the smaller set. Two pieces of one tandem
    duplication carry the same unit, possibly rotated, so they share nearly
    all k-mers; two unrelated events share almost none."""
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def attach_kmers(alts, vcf_path, min_size):
    """Add the k-mer set of each allele's variant sequence (inserted sequence
    for INS, deleted reference for DEL, reference span for INV). Records are
    read in the same order and with the same filter as P.read_vcf."""
    seqs = []
    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            ref = f[3]
            for alt in f[4].split(","):
                d = len(alt) - len(ref)
                if d <= -min_size:
                    seqs.append(ref)
                elif d >= min_size:
                    seqs.append(alt)
                elif len(ref) >= min_size and P.identity(alt, P.revcomp(ref)) >= 0.90:
                    seqs.append(ref)
    if len(seqs) != len(alts):
        sys.exit(f"{vcf_path}: allele count mismatch {len(seqs)} vs {len(alts)}")
    for a, sq in zip(alts, seqs):
        a["kmers"] = kmers(sq)


def size_ok(a, b, frac, tmin):
    return abs(a - b) <= max(tmin, frac * max(a, b))


def compatible(alt, t, frac, tmin, pmin):
    """Same test precision_recall_types uses: kind, size and position."""
    if alt["kind"] not in P.TYPE_MAP.get(t["type"], set()):
        return False
    ptol = max(pmin, t["len"] // 2)
    stol = max(tmin, frac * t["len"])
    return (alt["end"] > t["start"] - ptol and alt["pos"] < t["end"] + ptol
            and abs(alt["size"] - t["len"]) <= stol)


def merge(alts, frac, tmin, pmin, pfrac=1.0, min_kmer=0.0, across=False, lookback=12):
    """Greedy clustering in position order. Returns a list of clusters, each a
    dict with the merged allele and the indices of its members."""
    order = sorted(range(len(alts)), key=lambda i: (alts[i]["pos"], alts[i]["end"]))
    clusters = []
    for i in order:
        a = alts[i]
        best, best_d = None, None
        for c in clusters[-lookback:]:
            if c["kind"] != a["kind"]:
                continue
            if not size_ok(a["size"], c["anchor_size"], frac, tmin):
                continue
            ptol = max(pmin, pfrac * c["anchor_size"])
            if a["pos"] > c["end"] + ptol or a["end"] < c["pos"] - ptol:
                continue
            if a["carriers"] & c["carriers"]:
                continue
            if across and a["pos"] in c["recs"]:
                continue
            if min_kmer > 0 and kmer_sim(a["kmers"], c["anchor_kmers"]) < min_kmer:
                continue
            d = abs(a["size"] - c["anchor_size"])
            if best is None or d < best_d:
                best, best_d = c, d
        if best is None:
            clusters.append({
                "kind": a["kind"], "anchor_size": a["size"],
                "pos": a["pos"], "end": a["end"],
                "carriers": set(a["carriers"]), "members": [i],
                "lead_size": a["size"], "lead_ac": a["ac"],
                "anchor_kmers": a.get("kmers"), "recs": {a["pos"]},
            })
        else:
            best["pos"] = min(best["pos"], a["pos"])
            best["end"] = max(best["end"], a["end"])
            best["carriers"] |= a["carriers"]
            best["members"].append(i)
            best["recs"].add(a["pos"])
            if a["ac"] > best["lead_ac"]:
                best["lead_ac"], best["lead_size"] = a["ac"], a["size"]
    merged = []
    for k, c in enumerate(clusters):
        merged.append({
            "pos": c["pos"], "end": c["end"], "alt_index": k,
            "kind": c["kind"], "size": c["lead_size"],
            "ac": len(c["carriers"]), "carriers": c["carriers"],
            "members": c["members"],
        })
    return merged


def attribute(alt, truth_c, truth, frac, tmin, pmin):
    """Index of the true variant this callset allele belongs to, or None."""
    hits = [k for k, t in enumerate(truth)
            if compatible(alt, t, frac, tmin, pmin)
            and alt["carriers"] and alt["carriers"] <= truth_c[k]]
    return hits[0] if len(hits) == 1 else (hits if hits else None)


def classify_failure(t, tc, alts, frac, tmin, pmin):
    """Why a true variant still lacks an exact carrier match, judged on the
    unmerged callset."""
    typed = [a for a in alts if compatible(a, t, frac, tmin, pmin)]
    if not typed:
        return "no_matching_record"
    if any(a["carriers"] == tc for a in typed):
        return "correct_before_fused_by_merge"
    subs = [a for a in typed if a["carriers"] and a["carriers"] <= tc]
    union = set().union(*[a["carriers"] for a in subs]) if subs else set()
    disjoint = sum(len(a["carriers"]) for a in subs) == len(union)
    if union == tc and disjoint:
        return "partition_not_merged"
    if union == tc:
        return "pieces_overlap"
    if any(a["carriers"] > tc for a in typed):
        return "extra_carriers"
    if union < tc:
        return "carriers_missing"
    return "other"


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--vcf-suffix", default=".wave.vcf")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-size", type=int, default=50)
    ap.add_argument("--size-tol-frac", type=float, default=0.20)
    ap.add_argument("--size-tol-min", type=int, default=50)
    ap.add_argument("--pos-tol-min", type=int, default=200)
    # merge parameters, kept separate from the scoring tolerances above so a
    # sensitivity analysis changes only the correction, never the yardstick
    ap.add_argument("--merge-size-frac", type=float, default=0.20)
    ap.add_argument("--merge-size-min", type=int, default=50)
    ap.add_argument("--merge-pos-min", type=int, default=200)
    ap.add_argument("--merge-pos-frac", type=float, default=1.0)
    ap.add_argument("--merge-min-kmer", type=float, default=0.0,
                    help="minimum shared 21-mer fraction between variant "
                         "sequences; 0 disables the sequence check")
    ap.add_argument("--across-records-only", action="store_true",
                    help="never join two alternates of the same record")
    args = ap.parse_args()
    F, T, M = args.size_tol_frac, args.size_tol_min, args.pos_tol_min

    pairs = []
    for t in sorted(glob.glob(os.path.join(args.sweep, "*.truth_svs.tsv"))):
        v = t[: -len(".truth_svs.tsv")] + args.vcf_suffix
        if os.path.exists(v):
            pairs.append((t, v))
    if not pairs:
        sys.exit(f"no truth/{args.vcf_suffix} pairs in {args.sweep}")
    sys.stderr.write(f"{len(pairs)} replicate(s)\n")
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    per_sv, n_before, n_after = [], 0, 0
    prec_before, prec_after = [], []
    cl_rows = []
    for tpath, vpath in pairs:
        rep = os.path.basename(tpath)[: -len(".truth_svs.tsv")]
        truth = P.read_truth(tpath)
        alts, samples = P.read_vcf(vpath, args.min_size)
        name_to_idx = {}
        for si, s in enumerate(samples):
            digits = "".join(ch for ch in s.split("#")[0] if ch.isdigit())
            if digits:
                name_to_idx[int(digits)] = si
        truth_c = [{name_to_idx[i] for i in t["carriers"] if i in name_to_idx}
                   for t in truth]

        if args.merge_min_kmer > 0:
            attach_kmers(alts, vpath, args.min_size)
        merged = merge(alts, args.merge_size_frac, args.merge_size_min,
                       args.merge_pos_min, args.merge_pos_frac, args.merge_min_kmer,
                       args.across_records_only)
        n_before += len(alts)
        n_after += len(merged)

        before = P.score(truth, alts, name_to_idx, F, T, M)
        after = P.score(truth, merged, name_to_idx, F, T, M)
        prec_before += P.precision(truth, alts, F, T, M)
        prec_after += P.precision(truth, merged, F, T, M)

        for k, (b, a) in enumerate(zip(before, after)):
            why = ""
            if not a["carrier_exact"]:
                why = classify_failure(truth[k], truth_c[k], alts, F, T, M)
            per_sv.append({
                "replicate": rep, "sv_id": b["sv_id"], "type": b["type"],
                "len": b["len"], "size_bin": b["size_bin"],
                "ac_true": b["ac_true"],
                "before_ac": b["ac_in_typed"], "before_car": b["carrier_exact"],
                "after_ac": a["ac_in_typed"], "after_car": a["carrier_exact"],
                "failure_after": why,
            })

        for c in merged:
            if len(c["members"]) < 2:
                continue
            owners = []
            for i in c["members"]:
                o = attribute(alts[i], truth_c, truth, F, T, M)
                owners.append(o)
            ints = [o for o in owners if isinstance(o, int)]
            if any(o is None for o in owners):
                status = "unattributed_member"
            elif any(isinstance(o, list) for o in owners):
                status = "ambiguous_member"
            elif len(set(ints)) == 1:
                status = "correct"
            else:
                status = "false_fusion"
            cl_rows.append({
                "replicate": rep, "pos": c["pos"], "end": c["end"],
                "kind": c["kind"], "size": c["size"],
                "n_members": len(c["members"]), "ac": c["ac"],
                "status": status,
                "truth_ids": ",".join(sorted({truth[o]["sv_id"] for o in ints})),
            })

    # ---------------------------------------------------------------- write
    f1 = ["replicate", "sv_id", "type", "len", "size_bin", "ac_true",
          "before_ac", "before_car", "after_ac", "after_car", "failure_after"]
    with open(args.out + ".perSV.tsv", "w") as fh:
        fh.write("\t".join(f1) + "\n")
        for r in per_sv:
            fh.write("\t".join(str(r[k]) for k in f1) + "\n")
    f2 = ["replicate", "pos", "end", "kind", "size", "n_members", "ac",
          "status", "truth_ids"]
    with open(args.out + ".clusters.tsv", "w") as fh:
        fh.write("\t".join(f2) + "\n")
        for r in cl_rows:
            fh.write("\t".join(str(r[k]) for k in f2) + "\n")

    def frac(rows, key):
        return sum(r[key] for r in rows) / len(rows) if rows else float("nan")

    L = []
    L.append(f"replicates {len(pairs)}   true variants {len(per_sv)}")
    L.append(f"merge rule: same kind; size within max({args.merge_size_min} bp, "
             f"{args.merge_size_frac:.0%}); start within max({args.merge_pos_min} bp, "
             f"{args.merge_pos_frac} x size); disjoint carriers; "
             f"shared 21-mers >= {args.merge_min_kmer}; "
             f"across records only: {args.across_records_only}")
    L.append("")
    L.append("RECALL                                 before    after")
    L.append(f"  exact count, type matched           {frac(per_sv,'before_ac'):.4f}   "
             f"{frac(per_sv,'after_ac'):.4f}")
    L.append(f"  exact carrier set, type matched     {frac(per_sv,'before_car'):.4f}   "
             f"{frac(per_sv,'after_car'):.4f}")
    gained = sum(1 for r in per_sv if r["after_car"] and not r["before_car"])
    lost = sum(1 for r in per_sv if r["before_car"] and not r["after_car"])
    L.append(f"  variants gained {gained}   variants lost {lost}")
    L.append("")
    L.append("PRECISION                              before    after")
    L.append(f"  alleles reported                    {len(prec_before):6d}   {len(prec_after):6d}")
    L.append(f"  matching a true variant             {frac(prec_before,'matches_truth'):.4f}   "
             f"{frac(prec_after,'matches_truth'):.4f}")
    L.append("")
    L.append("MERGED CLUSTERS (two or more records)")
    st = {}
    for r in cl_rows:
        st[r["status"]] = st.get(r["status"], 0) + 1
    L.append(f"  total {len(cl_rows)}")
    for k in sorted(st):
        L.append(f"    {k:22} {st[k]}")
    L.append("")
    L.append("RECALL by type and size, exact carrier set: before -> after")
    bins = sorted({r["size_bin"] for r in per_sv if str(r["size_bin"]).isdigit()}, key=int)
    for ty in sorted({r["type"] for r in per_sv}):
        cells = []
        for b in bins:
            s = [r for r in per_sv if r["type"] == ty and r["size_bin"] == b]
            if s:
                cells.append(f"{b}:{frac(s,'before_car'):.3f}->{frac(s,'after_car'):.3f}(n={len(s)})")
        L.append(f"  {ty:4} " + "  ".join(cells))
    L.append("")
    L.append("REMAINING FAILURES after merge, judged on the unmerged callset")
    for ty in sorted({r["type"] for r in per_sv}):
        s = [r for r in per_sv if r["type"] == ty and not r["after_car"]]
        cnt = {}
        for r in s:
            cnt[r["failure_after"]] = cnt.get(r["failure_after"], 0) + 1
        L.append(f"  {ty:4} n={len(s):4d}  " +
                 "  ".join(f"{k} {v}" for k, v in sorted(cnt.items())))
    report = "\n".join(L)
    with open(args.out + ".summary.txt", "w") as fh:
        fh.write(report + "\n")
    print(report)
    print(f"\nwrote {args.out}.perSV.tsv, .clusters.tsv, .summary.txt")


if __name__ == "__main__":
    main()
