#!/usr/bin/env python3
"""
piece_diversity.py

Tests why a true structural variant is split across several callset alleles.

For every true variant, the "pieces" are the callset alternates that match it
in type, size and position and whose carriers are a subset of its carriers.
For each variant the script records how many pieces there are, at how many
distinct positions, whether the pieces have identical length, and how many
substitutions separate the pieces' variant sequences.

If pieces of one event sit at one position, have one length and differ by a
few substitutions, the split comes from sequence variation inside the
variant, not from placement of the event.

Usage:
    python3 analysis/piece_diversity.py --sweep ~/pangenome_rebuild/sweep54d \\
        --vcf-suffix .wave.vcf --out out/pieces
"""

import argparse
import glob
import itertools
import os
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import precision_recall_types as P  # noqa: E402
from merge_fix import compatible  # noqa: E402


def read_seqs(path, min_size):
    """Variant sequence of each alternate, in the order P.read_vcf emits them."""
    out = []
    with open(path) as fh:
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
                    out.append(ref)
                elif d >= min_size:
                    out.append(alt)
                elif len(ref) >= min_size and P.identity(alt, P.revcomp(ref)) >= 0.90:
                    out.append(ref)
    return out


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", required=True)
    ap.add_argument("--vcf-suffix", default=".wave.vcf")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-size", type=int, default=50)
    a_ = ap.parse_args()
    F, T, M = 0.20, 50, 200

    rows = []
    for tpath in sorted(glob.glob(os.path.join(a_.sweep, "*.truth_svs.tsv"))):
        vpath = tpath[: -len(".truth_svs.tsv")] + a_.vcf_suffix
        if not os.path.exists(vpath):
            continue
        rep = os.path.basename(tpath)[: -len(".truth_svs.tsv")]
        truth = P.read_truth(tpath)
        alts, samples = P.read_vcf(vpath, a_.min_size)
        seqs = read_seqs(vpath, a_.min_size)
        if len(seqs) != len(alts):
            sys.exit(f"{vpath}: allele count mismatch")
        for a, s in zip(alts, seqs):
            a["seq"] = s.upper()
        n2i = {}
        for si, s in enumerate(samples):
            dg = "".join(ch for ch in s.split("#")[0] if ch.isdigit())
            if dg:
                n2i[int(dg)] = si
        for t in truth:
            tc = {n2i[i] for i in t["carriers"] if i in n2i}
            pieces = [a for a in alts if compatible(a, t, F, T, M)
                      and a["carriers"] and a["carriers"] <= tc]
            exact = any(a["carriers"] == tc for a in alts
                        if compatible(a, t, F, T, M))
            lens = {len(a["seq"]) for a in pieces}
            hd = []
            if len(lens) == 1 and len(pieces) > 1:
                for x, y in itertools.combinations(pieces, 2):
                    hd.append(hamming(x["seq"], y["seq"]))
            union = set().union(*[a["carriers"] for a in pieces]) if pieces else set()
            rows.append({
                "replicate": rep, "sv_id": t["sv_id"], "type": t["type"],
                "len": t["len"], "size_bin": t["size_bin"], "ac_true": t["ac"],
                "exact_single": int(exact),
                "n_pieces": len(pieces),
                "n_positions": len({a["pos"] for a in pieces}),
                "same_length": int(len(lens) == 1) if pieces else "",
                "pieces_cover_truth": int(union == tc) if pieces else 0,
                "min_subs": min(hd) if hd else "",
                "max_subs": max(hd) if hd else "",
            })

    cols = list(rows[0].keys())
    with open(a_.out + ".tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    L = ["Split variants: no single allele carries the exact carrier set, but",
         "two or more pieces exist.", ""]
    hdr = (f"{'type':4} {'bin':>6} {'n':>4} {'split':>5} {'pieces med':>10} "
           f"{'one pos':>8} {'same len':>8} {'cover':>6} {'min subs med':>12}")
    L.append(hdr)
    bins = sorted({r["size_bin"] for r in rows if str(r["size_bin"]).isdigit()}, key=int)
    for ty in sorted({r["type"] for r in rows}):
        for b in bins:
            s = [r for r in rows if r["type"] == ty and r["size_bin"] == b]
            if not s:
                continue
            sp = [r for r in s if not r["exact_single"] and r["n_pieces"] >= 2]
            if not sp:
                L.append(f"{ty:4} {b:>6} {len(s):4d} {0:5d}")
                continue
            subs = [r["min_subs"] for r in sp if r["min_subs"] != ""]
            L.append(f"{ty:4} {b:>6} {len(s):4d} {len(sp):5d} "
                     f"{st.median(r['n_pieces'] for r in sp):10.1f} "
                     f"{sum(r['n_positions'] == 1 for r in sp):8d} "
                     f"{sum(r['same_length'] == 1 for r in sp):8d} "
                     f"{sum(r['pieces_cover_truth'] for r in sp):6d} "
                     f"{(st.median(subs) if subs else float('nan')):12.1f}")
    L += ["",
          "split       variants with no exact single allele and two or more pieces",
          "pieces med  median number of pieces per split variant",
          "one pos     split variants whose pieces all share one position",
          "same len    split variants whose pieces all have one length",
          "cover       split variants whose pieces' carriers together equal the truth",
          "min subs    median, over split variants, of the smallest substitution",
          "            count between two pieces (same length pieces only)"]
    rep = "\n".join(L)
    open(a_.out + ".summary.txt", "w").write(rep + "\n")
    print(rep)
    print(f"\nwrote {a_.out}.tsv and {a_.out}.summary.txt")


if __name__ == "__main__":
    main()
