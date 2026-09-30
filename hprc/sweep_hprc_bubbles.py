#!/usr/bin/env python3
"""
sweep_hprc_bubbles.py

For each top-level bubble in the HPRC raw VCF that contains a deletion-type
structural event, reconstruct the event's true haplotype carrier set from the
AT traversal field, then ask how the post-processed (wave) VCF represents it.

Truth here is not simulated. It is the graph's own traversal encoding, grouped
by shared breakpoint. The comparison is therefore between two representations
of the same underlying object, with no ground-truth assumption and no code of
ours inside the graph construction.

Per event it reports:
    true_copies      haplotype copies carrying the event (from AT grouping)
    true_af          true_copies / AN
    event_max_ac     largest single ALT allele count among the event's own
                     alleles, i.e. the frequency a naive reader would take
    wave_n_records   how many wave records carry part of the event
    wave_best_af     AF of the single wave record holding the most copies
    wave_recall      union of wave copies / true_copies
    wave_dup_haps    haplotype copies claimed by more than one wave record

Requires: bcftools on PATH, network access to the HPRC S3 bucket.
Nothing is downloaded in full; bcftools fetches only the BGZF blocks covering
each requested region, so cost scales with total span queried, not with the
fields printed.

Usage:
    python3 sweep_hprc_bubbles.py --out chr20_sweep.tsv
    python3 sweep_hprc_bubbles.py --regions my_windows.txt --out sweep.tsv
    python3 sweep_hprc_bubbles.py --limit 3 --out test.tsv      # smoke test
"""

import argparse
import collections
import os
import re
import statistics
import shutil
import subprocess
import sys
import time

BASE = ("https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/"
        "freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38/")
RAW_DEFAULT = BASE + "hprc-v2.1-mc-grch38.raw.vcf.gz"
WAVE_DEFAULT = BASE + "hprc-v2.1-mc-grch38.wave.vcf.gz"

NODE_RE = re.compile(r"[<>](\d+)")

COLUMNS = [
    "chrom", "bubble_pos", "event_id", "event_size_bp", "n_ref_nodes",
    "gap_ref_nodes", "left_node", "right_node", "an", "n_alt_alleles",
    "n_event_alleles", "n_exact_subgroups", "size_spread_bp",
    "true_copies", "true_af", "event_max_ac", "event_max_af",
    "raw_af_ratio", "wave_n_records", "wave_best_copies", "wave_best_af",
    "wave_best_recall", "wave_union_copies", "wave_recall",
    "wave_records_for_90pct", "wave_dup_haps", "wave_extra_haps",
    "carrier_del_bp", "other_del_bp", "dispersal_ratio", "n_uncalled",
    "verdict",
]


# ----------------------------------------------------------------------------
# shell helpers
# ----------------------------------------------------------------------------

class ToolMissing(RuntimeError):
    """bcftools is not on PATH. Fatal: every query would fail identically."""


def preflight(raw, wave):
    """
    Verify bcftools exists and both remote VCFs are reachable BEFORE the sweep.

    Without this a missing bcftools makes every window fail with rc=255, the
    run writes zero rows, and still prints DONE -- an empty result that looks
    like a completed sweep. Fail loudly at the start instead.
    """
    if shutil.which("bcftools") is None:
        raise ToolMissing(
            "bcftools not found on PATH.\n"
            "  If you use conda, you are probably in the wrong environment:\n"
            "      conda activate <env>\n"
            "  Note that 'nohup ... &' inherits the launching shell's env.")
    for label, url in (("raw", raw), ("wave", wave)):
        probe = subprocess.run(
            ["bcftools", "query", "-r", "chr20:25783205-25783205", "-f", "%POS\\n", url],
            capture_output=True, text=True, timeout=300)
        if probe.returncode != 0:
            tail = probe.stderr.strip().splitlines()[-1] if probe.stderr.strip() else "unknown"
            raise RuntimeError(f"cannot read {label} VCF: {tail}")
    sys.stderr.write("preflight OK: bcftools present, both VCFs reachable\n")


def bcf(fmt, region, vcf, timeout=1800, retries=3, backoff=5):
    """
    Run bcftools query, returning rstripped lines.

    Retries transient remote failures with linear backoff. Raises ToolMissing
    if bcftools itself is absent, so the caller aborts rather than silently
    recording every window as empty.
    """
    cmd = ["bcftools", "query", "-r", region, "-f", fmt, vcf]
    last = "unknown"
    for attempt in range(1, retries + 1):
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        except FileNotFoundError:
            raise ToolMissing("bcftools disappeared from PATH mid-run")
        except subprocess.TimeoutExpired:
            last = "timeout"
        else:
            if p.returncode == 0:
                return [ln for ln in p.stdout.rstrip("\n").split("\n") if ln]
            last = p.stderr.strip().splitlines()[-1] if p.stderr.strip() else "unknown"
            if "not found" in last and "bcftools" in last:
                raise ToolMissing(last)
        if attempt < retries:
            sys.stderr.write(f"  [retry {attempt}/{retries - 1}] {region}: {last}\n")
            time.sleep(backoff * attempt)
    sys.stderr.write(f"  [FAILED after {retries}] {region}: {last}\n")
    return None


def nodes(traversal):
    """Node IDs of a traversal string, orientation stripped."""
    return NODE_RE.findall(traversal)


def hap_copies(gt_blob, wanted):
    """
    Parse a '[%SAMPLE=%GT;]' blob into the set of (sample, slot) pairs whose
    allele index is in `wanted`. Slot is the position within the genotype, so
    each element is one haplotype copy rather than one diploid sample.
    """
    out = set()
    for chunk in gt_blob.strip(";").split(";"):
        if not chunk or "=" not in chunk:
            continue
        sample, gt = chunk.split("=", 1)
        for slot, allele in enumerate(gt.replace("/", "|").split("|")):
            if allele.isdigit() and int(allele) in wanted:
                out.add((sample, slot))
    return out


# ----------------------------------------------------------------------------
# event reconstruction from the raw VCF
# ----------------------------------------------------------------------------

def _breakpoints(ref, alts, ats, min_event):
    """
    Per ALT allele, locate its largest gap in reference-path index.

    Returns a list of dicts with the flanking reference-path indices (p1, p2),
    the flanking node IDs, the gap in nodes and the deleted bp. Alleles with no
    usable projection onto the reference path are dropped.
    """
    refn = nodes(ats[0])
    if len(refn) < 3:
        return [], refn
    ref_index = {}
    for i, n in enumerate(refn):
        ref_index.setdefault(n, i)

    out = []
    for i, alt in enumerate(alts):
        deleted = len(ref) - len(alt)
        if deleted < min_event or i + 1 >= len(ats):
            continue
        kept = sorted({ref_index[n] for n in nodes(ats[i + 1]) if n in ref_index})
        if len(kept) < 2:
            continue
        gap = p1 = p2 = 0
        for a, b in zip(kept, kept[1:]):
            if b - a > gap:
                gap, p1, p2 = b - a, a, b
        if gap < 2:
            continue
        out.append({"idx": i + 1, "p1": p1, "p2": p2, "gap": gap,
                    "deleted": deleted, "left": refn[p1], "right": refn[p2]})
    return out, refn


def find_events(ref, alts, ats, acs, min_event,
                tol_nodes=10, tol_frac=0.02, exact=False):
    """
    Group ALT alleles into deletion events by breakpoint.

    Requiring identical flanking node IDs is too strict: haplotypes carrying
    one deletion often rejoin the reference a few nodes apart, which splits a
    single event into many near-identical "events" and inflates any per-event
    count. Alleles are therefore clustered by single linkage on breakpoint
    position, with tolerance

        tol = max(tol_nodes, tol_frac * gap)

    in reference-path node units. `exact=True` restores the old identical-node
    behaviour so the effect of the threshold can be reported.

    Each returned event records how many distinct exact node pairs it absorbed
    (n_exact_subgroups) and the spread of deleted lengths among its alleles
    (size_spread_bp), both of which say how heterogeneous the grouping was.
    """
    if not ats or len(ats) < 2:
        return {}

    per, refn = _breakpoints(ref, alts, ats, min_event)
    if not per:
        return {}

    n = len(per)
    parent = list(range(n))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = root(a), root(b)
        if ra != rb:
            parent[rb] = ra

    for a in range(n):
        for b in range(a + 1, n):
            A, B = per[a], per[b]
            if exact:
                if A["left"] == B["left"] and A["right"] == B["right"]:
                    union(a, b)
                continue
            tol = max(tol_nodes, tol_frac * max(A["gap"], B["gap"]))
            if abs(A["p1"] - B["p1"]) <= tol and abs(A["p2"] - B["p2"]) <= tol:
                union(a, b)

    clusters = collections.defaultdict(list)
    for i in range(n):
        clusters[root(i)].append(per[i])

    events = {}
    for members in clusters.values():
        copies = collections.Counter()
        for m in members:
            copies[(m["left"], m["right"])] += (
                acs[m["idx"] - 1] if m["idx"] - 1 < len(acs) else 0)
        rep = max(copies, key=copies.get)          # most-supported breakpoint pair
        sizes = [m["deleted"] for m in members]
        events[rep] = {
            "alt_idx": [m["idx"] for m in members],
            "copies": sum(copies.values()),
            "sizes": sizes,
            "gap": max(m["gap"] for m in members),
            "n_exact_subgroups": len({(m["left"], m["right"]) for m in members}),
            "size_spread": max(sizes) - min(sizes),
        }
    return events


# ----------------------------------------------------------------------------
# wave-VCF side
# ----------------------------------------------------------------------------

def wave_span(chrom, start, end, wave_vcf):
    """
    Every wave record in the bubble span, with per-haplotype deleted bp.

    Returns (records, deleted_bp, called) where

        records    [(pos, {alt_idx: deleted_bp}, {(sample, slot): alt_idx})]
        deleted_bp {(sample, slot): total deleted bp across the whole span}
        called     haplotypes with at least one non-missing call in the span

    Size matching is deliberately not applied here. Asking whether any record
    has the same length as the event assumes both representations agree on
    where the event begins and ends, which is exactly what is in question at
    high-diversity bubbles: the wave VCF may cover the interval with one much
    larger record, or with many smaller ones, and the earlier size filter
    scored both as "no record found".
    """
    fmt = "%POS\\t%REF\\t%ALT\\t[%SAMPLE=%GT;]\\n"
    lines = bcf(fmt, f"{chrom}:{start}-{end}", wave_vcf)
    if lines is None:
        return None
    records, deleted_bp, called = [], collections.Counter(), set()
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        pos, ref, alt, gts = parts[0], parts[1], parts[2], parts[3]
        alt_del = {j + 1: len(ref) - len(a) for j, a in enumerate(alt.split(","))}
        carriers = {}
        for chunk in gts.strip(";").split(";"):
            if not chunk or "=" not in chunk:
                continue
            sample, gt = chunk.split("=", 1)
            for slot, a in enumerate(gt.replace("/", "|").split("|")):
                h = (sample, slot)
                if not a.isdigit():
                    continue
                called.add(h)
                j = int(a)
                if j in alt_del and alt_del[j] > 0:
                    deleted_bp[h] += alt_del[j]
                    carriers[h] = j
        records.append((int(pos), alt_del, carriers))
    return records, deleted_bp, called


def score_event(truth, event_size, records, deleted_bp, called, tol_frac, tol_min):
    """
    Score how the wave VCF represents one event.

    Two independent questions, kept separate because they have different
    downstream consequences:

      dispersal  do the event's carrier haplotypes carry the deleted sequence
                 at all, anywhere in the span? (ratio of carrier excess in
                 summed deleted bp to the event size; ~1 means represented,
                 ~0 means absent, >1 means the wave record spans more than the
                 event does)

      records    how many records the carriers are spread over, and whether a
                 single record of compatible size reproduces the carrier set

    A haplotype with no call in the span contributes nothing to the dispersal
    sum, so `n_uncalled` is reported alongside rather than folded in.
    """
    tol = max(tol_min, tol_frac * event_size)
    matched = []
    for pos, alt_del, carriers in records:
        if not any(abs(d - event_size) <= tol for d in alt_del.values()):
            continue
        hit = {h for h, j in carriers.items()
               if abs(alt_del[j] - event_size) <= tol} & truth
        if hit:
            matched.append((pos, {h for h, j in carriers.items()
                                  if abs(alt_del[j] - event_size) <= tol}, len(hit)))
    matched.sort(key=lambda x: -x[2])

    carried = [deleted_bp.get(h, 0) for h in truth if h in called]
    others = [deleted_bp.get(h, 0) for h in called - truth]
    med = lambda v: sorted(v)[len(v) // 2] if v else 0
    excess = med(carried) - med(others)

    counted = collections.Counter()
    for _, carriers, _ in matched:
        for h in carriers:
            counted[h] += 1
    union = set().union(*[c for _, c, _ in matched]) if matched else set()

    cum, n90 = 0, 0
    for _, _, hit in matched:
        cum += hit
        n90 += 1
        if cum >= 0.9 * len(truth):
            break

    return {
        "n_records": len(matched),
        "best_copies": matched[0][2] if matched else 0,
        "best_total": len(matched[0][1]) if matched else 0,
        "union_hit": len(union & truth),
        "n90": n90 if matched else 0,
        "dup": sum(1 for h, c in counted.items() if c > 1 and h in truth),
        "extra": len(union - truth),
        "carrier_del_bp": med(carried),
        "other_del_bp": med(others),
        "dispersal_ratio": round(excess / event_size, 3) if event_size else 0,
        "n_uncalled": sum(1 for h in truth if h not in called),
    }


# ----------------------------------------------------------------------------
# per-window driver
# ----------------------------------------------------------------------------

def process_window(chrom, wstart, wend, args, writer, seen):
    region = f"{chrom}:{wstart}-{wend}"
    fmt = "%POS\\t%REF\\t%ALT\\t%INFO/LV\\t%INFO/AC\\t%INFO/AN\\n"
    lines = bcf(fmt, region, args.raw)
    if lines is None:
        sys.stderr.write(f"{region}: SKIPPED (query failed)\n")
        return 0, 1
    if not lines:
        sys.stderr.write(f"{region}: 0 records returned -- empty region or bad contig name\n")

    candidates = []
    for line in lines:
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        pos, ref, alt, lv, ac, an = parts[:6]
        if args.top_level_only and lv != "0":
            continue
        if not any(len(ref) - len(a) >= args.min_event for a in alt.split(",")):
            continue
        candidates.append(int(pos))

    sys.stderr.write(f"{region}: {len(lines)} records, {len(candidates)} candidate bubbles\n")
    n_events = 0
    failed_bubbles = 0

    for pos in candidates:
        if (chrom, pos) in seen:
            continue
        fmt2 = "%POS\\t%REF\\t%ALT\\t%INFO/AT\\t%INFO/AC\\t%INFO/AN\\t[%SAMPLE=%GT;]\\n"
        rows = bcf(fmt2, f"{chrom}:{pos}-{pos}", args.raw)
        if not rows:
            failed_bubbles += 1
            continue
        row = next((r for r in rows if r.split("\t")[0] == str(pos)), None)
        if row is None:
            continue
        parts = row.split("\t")
        if len(parts) < 7:
            continue
        _, ref, alt, at, ac, an, gts = parts[:7]

        alts = alt.split(",")
        ats = at.split(",")
        acs = [int(x) if x.lstrip("-").isdigit() else 0 for x in ac.split(",")]
        try:
            an_i = int(an)
        except ValueError:
            an_i = 0
        if an_i <= 0:
            continue

        seen.add((chrom, pos))          # prevents duplicate rows when the
                                        # bubble spans several query windows
        events = find_events(ref, alts, ats, acs, args.min_event,
                             tol_nodes=args.cluster_tol_nodes,
                             tol_frac=args.cluster_tol_frac,
                             exact=args.exact_grouping)
        if not events:
            continue

        refn_count = len(nodes(ats[0])) if ats else 0
        span_end_pre = pos + len(ref)
        span = wave_span(chrom, pos, span_end_pre, args.wave)
        if span is None:
            failed_bubbles += 1
            continue
        wrecords, wdeleted, wcalled = span
        span_end = pos + len(ref)

        for (left, right), ev in events.items():
            if ev["copies"] < args.min_copies:
                continue
            size = int(statistics.median(ev["sizes"]))
            # largest count among THIS event's alleles, not the whole bubble:
            # a bubble hosts several events, and the bubble-wide maximum can
            # belong to an unrelated allele, which inverts the ratio.
            event_max_ac = max((acs[i - 1] for i in ev["alt_idx"]
                                if 0 < i <= len(acs)), default=0)
            truth = hap_copies(gts, set(ev["alt_idx"]))
            if not truth:
                continue

            sc = score_event(truth, size, wrecords, wdeleted, wcalled,
                             args.tol_frac, args.tol_min)

            # A single verdict, so the ambiguous middle cannot be silently read
            # as loss. RECOVERED: one compatible record reproduces the carrier
            # set. FRAGMENTED: several do. RESCALED: the sequence is present
            # (dispersal ~1 or more) but no record matches the event's size --
            # the two representations disagree on the unit of variation, which
            # is not the same as the event being absent. ABSENT: carriers show
            # no excess deleted sequence at all.
            if sc["n_records"] == 1 and sc["union_hit"] >= len(truth):
                verdict = "RECOVERED"
            elif sc["n_records"] > 1:
                verdict = "FRAGMENTED"
            elif sc["n_records"] == 1:
                verdict = "PARTIAL"
            elif sc["dispersal_ratio"] >= 0.5:
                verdict = "RESCALED"
            elif sc["n_uncalled"] > 0.25 * len(truth):
                verdict = "UNCALLED"
            else:
                verdict = "ABSENT"

            writer.write("\t".join(str(x) for x in [
                chrom, pos, f"{left}_{right}", size, refn_count, ev["gap"],
                left, right, an_i, len(alts), len(ev["alt_idx"]),
                ev["n_exact_subgroups"], ev["size_spread"],
                len(truth), round(len(truth) / an_i, 4),
                event_max_ac, round(event_max_ac / an_i, 4),
                round(len(truth) / event_max_ac, 3) if event_max_ac else "NA",
                sc["n_records"], sc["best_copies"],
                round(sc["best_total"] / an_i, 4),
                round(sc["best_copies"] / len(truth), 4),
                sc["union_hit"], round(sc["union_hit"] / len(truth), 4),
                sc["n90"], sc["dup"], sc["extra"],
                sc["carrier_del_bp"], sc["other_del_bp"],
                sc["dispersal_ratio"], sc["n_uncalled"], verdict,
            ]) + "\n")
            writer.flush()
            n_events += 1

    note = f", {failed_bubbles} bubbles skipped after query failure" if failed_bubbles else ""
    sys.stderr.write(f"{region}: wrote {n_events} events{note}\n")
    return n_events, failed_bubbles


# ----------------------------------------------------------------------------

def default_windows(chrom, start, stop, width, step):
    w, out = start, []
    while w + width <= stop:
        out.append((chrom, w, w + width))
        w += step
    return out


def load_regions(path):
    out = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" in line:                      # BED-like
                f = line.split("\t")
                out.append((f[0], int(f[1]), int(f[2])))
            else:                                 # chr:start-end
                c, rng = line.split(":")
                s, e = rng.split("-")
                out.append((c, int(s), int(e)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", default=RAW_DEFAULT)
    ap.add_argument("--wave", default=WAVE_DEFAULT)
    ap.add_argument("--out", required=True)
    ap.add_argument("--regions", help="file of chr:start-end or BED lines")
    ap.add_argument("--chrom", default="chr20")
    ap.add_argument("--start", type=int, default=2_000_000)
    ap.add_argument("--stop", type=int, default=62_000_000)
    ap.add_argument("--width", type=int, default=400_000)
    ap.add_argument("--step", type=int, default=2_500_000)
    ap.add_argument("--limit", type=int, help="process only the first N windows")
    ap.add_argument("--min-event", type=int, default=1000,
                    help="minimum deleted bp for an ALT to count (default 1000)")
    ap.add_argument("--min-copies", type=int, default=2,
                    help="skip events below this many haplotype copies")
    ap.add_argument("--cluster-tol-nodes", type=int, default=10,
                    help="breakpoint clustering tolerance, reference-path nodes")
    ap.add_argument("--cluster-tol-frac", type=float, default=0.02,
                    help="clustering tolerance as a fraction of gap size")
    ap.add_argument("--exact-grouping", action="store_true",
                    help="require identical flanking nodes (old behaviour; "
                         "use to report sensitivity to the threshold)")
    ap.add_argument("--tol-frac", type=float, default=0.20,
                    help="wave size match tolerance as fraction of event size")
    ap.add_argument("--tol-min", type=int, default=200,
                    help="wave size match tolerance floor in bp")
    ap.add_argument("--top-level-only", action="store_true", default=True)
    ap.add_argument("--include-nested", dest="top_level_only", action="store_false")
    args = ap.parse_args()

    windows = load_regions(args.regions) if args.regions else default_windows(
        args.chrom, args.start, args.stop, args.width, args.step)
    if args.limit:
        windows = windows[:args.limit]

    seen = set()
    resume = os.path.exists(args.out) and os.path.getsize(args.out) > 0
    if resume:
        with open(args.out) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) > 2 and f[0] != "chrom":
                    seen.add((f[0], int(f[1])))
        sys.stderr.write(f"resuming: {len(seen)} bubbles already in {args.out}\n")

    try:
        preflight(args.raw, args.wave)
    except ToolMissing as exc:
        sys.stderr.write(f"\nFATAL: {exc}\n")
        sys.exit(2)
    except Exception as exc:
        sys.stderr.write(f"\nFATAL: {exc}\n")
        sys.exit(3)

    t0 = time.time()
    total = 0
    failed = 0
    try:
        with open(args.out, "a" if resume else "w") as writer:
            if not resume:
                writer.write("\t".join(COLUMNS) + "\n")
            for i, (chrom, ws, we) in enumerate(windows, 1):
                sys.stderr.write(f"[{i}/{len(windows)}] ")
                n_ev, n_fail = process_window(chrom, ws, we, args, writer, seen)
                total += n_ev
                failed += n_fail
    except ToolMissing as exc:
        sys.stderr.write(f"\nFATAL: {exc}\nPartial results kept in {args.out}\n")
        sys.exit(2)

    status = "DONE" if failed == 0 else "INCOMPLETE"
    sys.stderr.write(f"\n{status}: {total} events, {len(windows)} windows, "
                     f"{time.time() - t0:.0f}s -> {args.out}\n")
    if failed:
        sys.stderr.write(f"WARNING: {failed} queries failed after retries. "
                         f"Rerun the same command to fill the gaps (it resumes).\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
