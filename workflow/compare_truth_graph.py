#!/usr/bin/env python3
"""
compare_truth_graph.py

Join simulated ground truth to a pangenome-graph VCF and quantify what the
graph got right.

Two different questions, deliberately kept separate:

  SNPs  -- exact position match. Was the site called at all, and is the
           allele count right? Cheap, unambiguous.

  SVs   -- interval overlap, NOT position match. Large SVs in a PGGB graph
           do not appear as one clean biallelic record: they land inside a
           multiallelic bubble that may carry a dozen alleles and whose
           reported POS can sit well away from the true breakpoint. So for
           each truth SV we collect every graph record whose span overlaps
           the truth interval, and ask the question that actually matters:

             is the true allele count recoverable from this site
             WITHOUT already knowing the answer?

           `exact_ac_present` = some ALT allele at an overlapping site has
           AC equal to the truth AC. `n_alleles` = how many candidates you
           would have to choose between. A site with n_alleles=14 that
           happens to contain the right AC is not a successful call -- it is
           a lottery ticket. Report both.

Usage
-----
  python compare_truth_graph.py \
      --prefix runs/easy_strat01 \
      --vcf runs/pggb_easy_strat01/*.smooth.final.anc.vcf \
      --out runs/easy_strat01.comparison

Outputs
-------
  <out>.snps.tsv   per-SNP recovery
  <out>.svs.tsv    per-SV recovery
  <out>.summary.tsv  one row -- append these across runs for the analysis
"""

import argparse
import csv
import json
import os
import sys


def read_vcf(path):
    """Return list of records: dict(pos0, ref, alts, acs, an, n_alleles)."""
    recs = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            pos0 = int(f[1]) - 1          # VCF is 1-based
            ref = f[3]
            alts = f[4].split(",")
            info = {}
            for kv in f[7].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    info[k] = v
            acs = []
            if "AC" in info:
                for a in info["AC"].split(","):
                    try:
                        acs.append(int(a))
                    except ValueError:
                        acs.append(-1)
            an = int(info.get("AN", -1))
            recs.append({
                "pos0": pos0,
                "end0": pos0 + len(ref),   # half-open span on the ref path
                "ref": ref,
                "alts": alts,
                "acs": acs,
                "an": an,
                "n_alleles": len(alts),
                "lv": int(info.get("LV", 0)),
            })
    recs.sort(key=lambda r: r["pos0"])
    return recs


def compare_snps(truth_path, recs):
    by_pos = {}
    for r in recs:
        by_pos.setdefault(r["pos0"], []).append(r)

    rows = []
    with open(truth_path) as fh:
        for t in csv.DictReader(fh, delimiter="\t"):
            pos = int(t["pos"])
            ac_true = int(t["ac"])
            hit = by_pos.get(pos, [])
            if not hit:
                rows.append({"pos": pos, "ac_true": ac_true, "found": 0,
                             "ac_graph": "", "an": "", "ac_error": ""})
                continue
            r = hit[0]
            # pick the ALT whose AC is closest to truth (multiallelic sites)
            best = min(r["acs"], key=lambda a: abs(a - ac_true)) if r["acs"] else -1
            rows.append({"pos": pos, "ac_true": ac_true, "found": 1,
                         "ac_graph": best, "an": r["an"],
                         "ac_error": best - ac_true})
    return rows


def compare_svs(truth_path, recs, slack_frac=0.5, slack_min=200):
    rows = []
    with open(truth_path) as fh:
        for t in csv.DictReader(fh, delimiter="\t"):
            s, e = int(t["start"]), int(t["end"])
            L = int(t["len"])
            ac_true = int(t["ac"])
            slack = max(slack_min, int(slack_frac * L))
            lo, hi = s - slack, e + slack

            over = [r for r in recs
                    if r["end0"] > lo and r["pos0"] < hi
                    and (len(r["ref"]) > 50 or any(len(a) > 50 for a in r["alts"]))]

            if not over:
                rows.append({
                    "sv_id": t["sv_id"], "type": t["type"], "len": L,
                    "size_bin": t.get("size_bin", ""),
                    "synthetic": t.get("synthetic", ""),
                    "ac_true": ac_true, "n_sites": 0, "n_alleles": 0,
                    "n_alleles_all": 0,
                    "best_pos": "", "best_reflen": "", "pos_offset": "",
                    "len_ratio": "", "exact_ac_present": 0,
                    "exact_ac_anywhere": 0, "min_an": "",
                })
                continue

            # "best" = the record whose REF span is closest in length to truth
            # An SV may be written REF-long (deletion-style) or ALT-long
            # (insertion/duplication-style). Comparing only REF length to
            # truth mismeasures the second class entirely.
            def var_size(r):
                return max(len(r["ref"]),
                           max((len(a) for a in r["alts"]), default=0))
            best = min(over, key=lambda r: abs(var_size(r) - L))
            # Alleles AT THE BEST BUBBLE, not summed over the interval: a
            # 30 kb SV overlaps ~300x more reference than a 100 bp one, so
            # a summed count would rise with size for trivial reasons and
            # manufacture the very trend we are trying to measure.
            n_alleles = best["n_alleles"]
            n_alleles_all = sum(r["n_alleles"] for r in over)
            exact = int(ac_true in best["acs"])
            exact_any = int(any(ac_true in r["acs"] for r in over))
            rows.append({
                "sv_id": t["sv_id"], "type": t["type"], "len": L,
                "size_bin": t.get("size_bin", ""),
                "synthetic": t.get("synthetic", ""),
                "ac_true": ac_true,
                "n_sites": len(over),
                "n_alleles": n_alleles,
                "n_alleles_all": n_alleles_all,
                "best_pos": best["pos0"],
                "best_reflen": var_size(best),
                "pos_offset": best["pos0"] - s,
                "len_ratio": round(var_size(best) / L, 4),
                "exact_ac_present": exact,
                "exact_ac_anywhere": exact_any,
                "min_an": min(r["an"] for r in over),
            })
    return rows


def write_tsv(path, rows, cols):
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prefix", required=True,
                    help="simulation out-prefix (finds .truth_*.tsv/.json)")
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    recs = read_vcf(args.vcf)
    stats = json.load(open(f"{args.prefix}.truth_stats.json"))

    snp_rows = compare_snps(f"{args.prefix}.truth_snps.tsv", recs)
    sv_rows = compare_svs(f"{args.prefix}.truth_svs.tsv", recs)

    write_tsv(f"{args.out}.snps.tsv", snp_rows,
              ["pos", "ac_true", "found", "ac_graph", "an", "ac_error"])
    write_tsv(f"{args.out}.svs.tsv", sv_rows,
              ["sv_id", "type", "len", "size_bin", "synthetic", "ac_true",
               "n_sites", "n_alleles", "n_alleles_all", "best_pos",
               "best_reflen", "pos_offset", "len_ratio",
               "exact_ac_present", "exact_ac_anywhere", "min_an"])

    n_snp = len(snp_rows)
    found = sum(r["found"] for r in snp_rows)
    ac_ok = sum(1 for r in snp_rows if r["found"] and r["ac_error"] == 0)
    n_sv = len(sv_rows)
    sv_found = sum(1 for r in sv_rows if r["n_sites"] > 0)
    sv_exact = sum(r["exact_ac_present"] for r in sv_rows)
    sv_clean = sum(1 for r in sv_rows
                   if r["exact_ac_present"] and r["n_alleles"] <= 2)

    summary = {
        "prefix": os.path.basename(args.prefix),
        "mode": stats.get("mode", ""),
        "contig": stats.get("contig", ""),
        "seed": stats.get("seed", ""),
        "n_snp_truth": n_snp,
        "snp_found_frac": round(found / n_snp, 5) if n_snp else "",
        "snp_ac_exact_frac": round(ac_ok / n_snp, 5) if n_snp else "",
        "n_sv_truth": n_sv,
        "sv_found_frac": round(sv_found / n_sv, 4) if n_sv else "",
        "sv_ac_present_frac": round(sv_exact / n_sv, 4) if n_sv else "",
        "sv_clean_frac": round(sv_clean / n_sv, 4) if n_sv else "",
        "n_vcf_records": len(recs),
        "frac_sites_an_lt_max": round(
            sum(1 for r in recs if 0 < r["an"] < max(
                (x["an"] for x in recs), default=0)) / len(recs), 5)
        if recs else "",
    }
    write_tsv(f"{args.out}.summary.tsv", [summary], list(summary))

    print(json.dumps(summary, indent=2))
    print("\n-- SV recovery by size bin --", file=sys.stderr)
    bins = {}
    for r in sv_rows:
        b = r["size_bin"] or "NA"
        bins.setdefault(b, []).append(r)
    for b in sorted(bins, key=lambda x: (x == "NA", int(x) if x.isdigit() else 0)):
        g = bins[b]
        print(f"  bin {b:>6}  n={len(g):2d}  "
              f"found={sum(1 for r in g if r['n_sites']>0)}  "
              f"ac_present={sum(r['exact_ac_present'] for r in g)}  "
              f"alleles_bubble={sum(r['n_alleles'] for r in g)/len(g):5.1f}  "
              f"alleles_all={sum(r['n_alleles_all'] for r in g)/len(g):6.1f}  "
              f"mean|offset|="
              f"{sum(abs(r['pos_offset']) for r in g if r['pos_offset']!='')/max(1,len(g)):8.0f}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
