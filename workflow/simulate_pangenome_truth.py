#!/usr/bin/env python3
"""
simulate_pangenome_truth.py

Simulate a haplotype panel with EXACT ground truth, for benchmarking
pangenome graph construction (PGGB / Minigraph-Cactus).

Design
------
We do NOT simulate sequence. We simulate the evolutionary PROCESS on top of
a real reference window, so that repeat content, segmental duplications and
TEs are authentic -- these are exactly what break graph builders, and a
neutrally-simulated sequence would hand you a false negative.

  1. msprime  -> genealogy of n haplotypes (with recombination)
  2. SNPs     -> dropped by msprime's mutation model
  3. SVs      -> sampled from a REAL callset (gnomAD-SV / HPRC) that falls in
                 the window, each assigned to one branch of the tree at its
                 start position; all descendants of that branch carry it
  4. Apply    -> edits spliced into the backbone per haplotype
  5. Truth    -> every variant's position, type, carrier set, allele frequency

Outputs
-------
  <prefix>.haplotypes.fa   all haplotypes, PanSN-named (bgzip + faidx this)
  <prefix>.truth_snps.tsv
  <prefix>.truth_svs.tsv
  <prefix>.truth_stats.json

Usage
-----
  python simulate_pangenome_truth.py \
      --backbone chr20_low_repeat_1Mb.fa \
      --sv-bed   gnomad_svs_in_window.bed \
      --window-start 30000000 \
      --n-hap 20 --n-svs 60 --seed 1 \
      --out-prefix runs/lowrep_rep01

Then:
  bgzip -@4 <prefix>.haplotypes.fa
  samtools faidx <prefix>.haplotypes.fa.gz
  pggb -i <prefix>.haplotypes.fa.gz -o out_pggb -n 21 -t 4 -p 95 -s 5000

Notes / known simplifications (v1 -- document these in the paper)
----------------------------------------------------------------
* An SV is placed using the marginal tree at its START position. With
  recombination inside the SV span the true carrier set could differ. Real
  recombination is suppressed inside inversions, so this is defensible, but
  it is an approximation.
* SVs are forced non-overlapping. Overlapping/nested SVs make coordinate
  bookkeeping much harder and are deferred to v2.
* Haploid samples. Real pangenomes are haplotype-resolved diploids; for
  graph construction each haplotype is a separate sequence anyway.
"""

import argparse
import json
import random
import sys

try:
    import msprime
    import tskit
except ImportError:
    sys.exit("Need msprime and tskit:  pip install msprime tskit")


# ----------------------------------------------------------------------
# sequence helpers
# ----------------------------------------------------------------------

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def read_single_fasta(path):
    """Read a one-record FASTA. Returns (name, uppercased sequence)."""
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    raise ValueError("backbone FASTA must contain exactly one record")
                name = line[1:].split()[0]
            else:
                chunks.append(line)
    if name is None:
        raise ValueError("no FASTA record found in backbone")
    return name, "".join(chunks).upper()


def write_fasta(fh, name, seq, width=60):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


# ----------------------------------------------------------------------
# SV input
# ----------------------------------------------------------------------

def load_svs(bed_path, window_start, window_len, min_len, max_len, min_af):
    """
    BED: chrom  start  end  TYPE  AF     (AF = real population allele freq)

    The AF column is essential. Without it, SV frequency is drawn from the
    genealogy alone and there is NO size-frequency coupling -- large
    deletions end up common, which is biologically false (they are
    deleterious and stay rare) and inflates per-haplotype deletion burden
    by ~10x. We instead use the real AF as a target and pick a clade of the
    matching size, which preserves genealogical structure AND the real
    size-frequency relationship.

    SVs below min_af are dropped: they would essentially never be seen in a
    panel of this size, so including them is unrealistic.
    """
    keep = []
    with open(bed_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split()
            if len(f) < 5:
                raise ValueError(
                    "SV BED needs 5 columns: chrom start end TYPE AF. "
                    "Re-extract with %INFO/AF in the bcftools query."
                )
            start = int(f[1]) - window_start
            end = int(f[2]) - window_start
            svtype = f[3].upper()
            if svtype not in ("DEL", "INV", "DUP"):
                continue
            try:
                af = float(f[4])
            except ValueError:
                continue  # '.' or missing
            if af < min_af or af >= 1.0:
                continue
            if start < 1 or end > window_len - 1 or end <= start:
                continue
            length = end - start
            if length < min_len or length > max_len:
                continue
            keep.append({"start": start, "end": end, "type": svtype,
                         "len": length, "af": af})
    return keep


def choose_nonoverlapping(svs, n_want, rng, pad=100):
    """Greedy: shuffle, accept an SV only if it clears all accepted ones."""
    rng.shuffle(svs)
    accepted = []
    for sv in svs:
        if all(sv["end"] + pad <= a["start"] or sv["start"] >= a["end"] + pad
               for a in accepted):
            accepted.append(sv)
            if len(accepted) >= n_want:
                break
    accepted.sort(key=lambda s: s["start"])
    return accepted


# ----------------------------------------------------------------------
# branch sampling
# ----------------------------------------------------------------------

def sample_branch_carriers(ts, position, target_k, rng):
    """
    Pick a clade in the marginal tree at `position` whose size is as close as
    possible to target_af * n_samples, breaking ties by branch length.

    This keeps carriers a genuine clade (so LD and haplotype structure are
    real) while forcing the panel frequency to match the SV's real
    population frequency -- which is what couples SV size to SV frequency.

    Returns the sorted tuple of sample IDs descending from that branch.
    """
    n = ts.num_samples
    target_k = max(1, min(n - 1, int(target_k)))
    tree = ts.at(position)
    by_count = {}
    for u in tree.nodes():
        if u == tree.root:
            continue
        bl = tree.branch_length(u)
        if bl <= 0:
            continue
        k = tree.num_samples(u)
        if k == 0 or k == n:
            continue
        by_count.setdefault(k, []).append((u, bl))
    if not by_count:
        return tuple()
    best_k = min(by_count, key=lambda k: (abs(k - target_k), k))
    cand = by_count[best_k]
    chosen = rng.choices([c[0] for c in cand],
                         weights=[c[1] for c in cand], k=1)[0]
    return tuple(sorted(tree.samples(chosen)))



def select_stratified(pool, args, rng, rs, window_len, backbone):
    """
    Sample SVs across log-spaced size bins, frequency drawn INDEPENDENTLY
    of size.

    A realistic panel is not a benchmark. Conditioning on presence in a
    20-haplotype panel selects almost only for small SVs, because large SVs
    are deleterious and therefore rare -- so you never observe how graph
    construction behaves at 10 kb or 30 kb. To measure error as a FUNCTION
    of size, size has to be a controlled variable.

    Real SVs from the window's callset are preferred (real breakpoints, real
    flanking context -- the point of using a real backbone). Only when a bin
    cannot be filled from the pool is an interval synthesised.
    """
    bins = [int(b) for b in args.bins.split(",")]
    chosen, n_synth = [], 0

    def clear(a, b, pad=200):
        if a < 1 or b > window_len - 1:
            return False
        if "N" in backbone[a:b]:
            return False
        return all(b + pad <= c["start"] or a >= c["end"] + pad
                   for c in chosen)

    # Force a minimum number of inversions per size bin. gnomAD SV is
    # short-read based and severely undercalls inversions -- the pool for a
    # 1 Mb window typically contains none at common frequency, so a
    # pool-driven sample can say nothing about the SV class where graph
    # representation is most likely to break. These are synthesised at real
    # genomic positions, so flanking sequence context is still real.
    for target in bins:
        placed_inv = 0
        tries = 0
        while placed_inv < args.min_inv_per_bin and tries < 5000:
            tries += 1
            start = int(rs.integers(1, max(2, window_len - target - 1)))
            end = start + target
            if not clear(start, end):
                continue
            chosen.append({
                "start": start, "end": end, "len": target, "type": "INV",
                "af": float("nan"),
                "target_k": int(rs.integers(2, 11)),
                "synth": True, "size_bin": target,
            })
            placed_inv += 1
            n_synth += 1

    for target in bins:
        lo, hi = target / 2.0, target * 2.0
        cand = [sv for sv in pool if lo <= sv["len"] <= hi]
        rng.shuffle(cand)
        placed = 0
        for sv in cand:
            if placed >= args.per_bin:
                break
            if not clear(sv["start"], sv["end"]):
                continue
            sv = dict(sv)
            sv["target_k"] = int(rs.integers(2, 11))
            sv["synth"] = False
            sv["size_bin"] = target
            chosen.append(sv)
            placed += 1
        tries = 0
        while placed < args.per_bin and tries < 5000:
            tries += 1
            start = int(rs.integers(1, max(2, window_len - target - 1)))
            end = start + target
            if not clear(start, end):
                continue
            chosen.append({
                "start": start, "end": end, "len": target,
                "type": rng.choice(["DEL", "DEL", "DUP", "INV"]),
                "af": float("nan"),
                "target_k": int(rs.integers(2, 11)),
                "synth": True, "size_bin": target,
            })
            placed += 1
            n_synth += 1
    chosen.sort(key=lambda s: s["start"])
    return chosen, n_synth


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backbone", required=True,
                    help="single-record FASTA of the reference window")
    ap.add_argument("--sv-bed", required=True,
                    help="BED of real SVs (chrom start end TYPE)")
    ap.add_argument("--window-start", type=int, default=0,
                    help="genome coord of backbone base 0, to offset the BED")
    ap.add_argument("--n-hap", type=int, default=20)
    ap.add_argument("--n-svs", type=int, default=60)
    ap.add_argument("--ne", type=float, default=10000)
    ap.add_argument("--mu", type=float, default=1.25e-8)
    ap.add_argument("--rec", type=float, default=1e-8)
    ap.add_argument("--min-sv-len", type=int, default=50)
    ap.add_argument("--max-sv-len", type=int, default=50000)
    ap.add_argument("--mode", choices=["realistic", "stratified"],
                    default="realistic",
                    help="realistic: panel counts from real gnomAD AF "
                         "(reproduces real pangenome composition, but yields "
                         "almost only small SVs). stratified: SVs across "
                         "log-spaced size bins with frequency decoupled from "
                         "size, to map error as a function of size.")
    ap.add_argument("--bins", default="100,300,1000,3000,10000,30000")
    ap.add_argument("--per-bin", type=int, default=4)
    ap.add_argument("--min-inv-per-bin", type=int, default=0,
                    help="synthesise at least this many inversions per\n"
                         "size bin; gnomAD supplies almost none")
    ap.add_argument("--min-af", type=float, default=0.0,
                    help="drop SVs below this real population AF; they would\n"
                         "essentially never appear in a panel this size")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--no-ancestor", action="store_true",
                    help="omit the ancestral backbone from the FASTA. Default "
                         "is to INCLUDE it as anc#1#<contig>: use it as the "
                         "`vg deconstruct` reference path so graph VCF "
                         "coordinates land in backbone space and match truth "
                         "directly. Costs realism, saves enormous pain.")
    args = ap.parse_args()

    rng = random.Random(args.seed)

    contig, backbone = read_single_fasta(args.backbone)
    L = len(backbone)
    print(f"[backbone] {contig}  {L:,} bp  "
          f"N-frac={backbone.count('N') / L:.4f}", file=sys.stderr)

    # ---- 1. genealogy --------------------------------------------------
    ts = msprime.sim_ancestry(
        samples=args.n_hap,
        ploidy=1,
        sequence_length=L,
        recombination_rate=args.rec,
        population_size=args.ne,
        random_seed=args.seed,
    )
    ts = msprime.sim_mutations(ts, rate=args.mu, random_seed=args.seed + 1)
    print(f"[msprime] {ts.num_trees} trees, {ts.num_sites} SNP sites",
          file=sys.stderr)

    # ---- 2. SNP table --------------------------------------------------
    # per-haplotype: list of (pos, derived_base); plus a truth row per site
    snp_edits = {i: [] for i in range(args.n_hap)}
    snp_truth = []
    for var in ts.variants():
        pos = int(var.site.position)
        if pos <= 0 or pos >= L:
            continue
        if backbone[pos] == "N":
            continue  # don't place SNPs in assembly gaps
        anc = backbone[pos]
        # map msprime's 0/1 alleles onto real bases: ancestral = backbone base
        alts = [b for b in "ACGT" if b != anc]
        derived = rng.choice(alts)
        carriers = [i for i, g in enumerate(var.genotypes) if g != 0]
        if not carriers or len(carriers) == args.n_hap:
            continue  # monomorphic in the panel: invisible to the graph
        for i in carriers:
            snp_edits[i].append((pos, derived))
        snp_truth.append({
            "pos": pos, "ref": anc, "alt": derived,
            "ac": len(carriers), "af": len(carriers) / args.n_hap,
            "carriers": ",".join(map(str, carriers)),
        })
    print(f"[snps] {len(snp_truth):,} polymorphic sites retained",
          file=sys.stderr)

    # ---- 3. SVs onto branches -----------------------------------------
    pool = load_svs(args.sv_bed, args.window_start, L,
                    args.min_sv_len, args.max_sv_len, args.min_af)
    print(f"[svs] {len(pool)} candidates in window (AF >= {args.min_af})",
          file=sys.stderr)
    import numpy as np
    rs = np.random.default_rng(args.seed + 2)

    if args.mode == "realistic":
        # Realise panel counts FIRST from the real AF, then resolve overlaps
        # among survivors -- the other order discards common SVs before they
        # can be drawn.
        survivors, n_zero = [], 0
        for sv in pool:
            k = int(rs.binomial(args.n_hap, sv["af"]))
            if k == 0:
                n_zero += 1
                continue
            sv = dict(sv)
            sv["target_k"] = k
            survivors.append(sv)
        print(f"[svs] {len(survivors)} present in panel; {n_zero} absent",
              file=sys.stderr)
        svs = choose_nonoverlapping(survivors, args.n_svs, rng)
    else:
        svs, n_synth = select_stratified(pool, args, rng, rs, L, backbone)
        print(f"[svs] stratified: {len(svs)} placed, {n_synth} synthesised",
              file=sys.stderr)

    sv_edits = {i: [] for i in range(args.n_hap)}
    sv_truth = []
    for k, sv in enumerate(svs):
        carriers = sample_branch_carriers(ts, sv["start"], sv["target_k"], rng)
        if not carriers or len(carriers) == args.n_hap:
            continue
        for i in carriers:
            sv_edits[i].append(sv)
        sv_truth.append({
            "sv_id": f"SV{k:04d}", "type": sv["type"],
            "start": sv["start"], "end": sv["end"], "len": sv["len"],
            "ac": len(carriers), "af": len(carriers) / args.n_hap,
            "gnomad_af": sv["af"],
            "size_bin": sv.get("size_bin", ""),
            "synthetic": sv.get("synth", False),
            "carriers": ",".join(map(str, carriers)),
        })
    print(f"[svs] {len(sv_truth)} polymorphic in panel", file=sys.stderr)

    # ---- 4. build haplotypes -------------------------------------------
    # CRITICAL: apply edits in DESCENDING coordinate order so that upstream
    # coordinates are never shifted by an earlier edit. SNPs falling inside
    # a deletion this haplotype carries are dropped (they are physically gone).
    fa_path = f"{args.out_prefix}.haplotypes.fa"
    surviving = {}          # pos -> haplotypes that actually carry the SNP
    with open(fa_path, "w") as fh:
        if not args.no_ancestor:
            write_fasta(fh, f"anc#1#{contig}", backbone)
        for i in range(args.n_hap):
            my_svs = sv_edits[i]
            dels = [(s["start"], s["end"]) for s in my_svs if s["type"] == "DEL"]

            def in_deletion(p):
                return any(a <= p < b for a, b in dels)

            edits = []
            for pos, base in snp_edits[i]:
                if not in_deletion(pos):
                    edits.append((pos, "SNP", base, None))
                    surviving.setdefault(pos, []).append(i)
            for s in my_svs:
                edits.append((s["start"], s["type"], None, s["end"]))
            edits.sort(key=lambda e: e[0], reverse=True)

            seq = list(backbone)
            for start, kind, base, end in edits:
                if kind == "SNP":
                    seq[start] = base
                elif kind == "DEL":
                    del seq[start:end]
                elif kind == "INV":
                    seq[start:end] = list(revcomp("".join(seq[start:end])))
                elif kind == "DUP":
                    seq[end:end] = seq[start:end]  # tandem, inserted after end
            write_fasta(fh, f"sim{i:02d}#1#{contig}", "".join(seq))
    print(f"[out] {fa_path}", file=sys.stderr)

    # A SNP whose carriers all sit inside a deletion is physically absent
    # from every haplotype that had it: the panel is monomorphic there and
    # the graph is CORRECT to emit nothing. Counting it as truth would
    # charge the graph for our own bookkeeping. Rebuild truth from what
    # actually made it into sequence.
    n_before = len(snp_truth)
    snp_truth = [r for r in snp_truth if 0 < len(surviving.get(r["pos"], [])) < args.n_hap]
    for r in snp_truth:
        car = sorted(surviving[r["pos"]])
        r["ac"] = len(car)
        r["af"] = len(car) / args.n_hap
        r["carriers"] = ",".join(map(str, car))
    print(f"[snps] {len(snp_truth)} still polymorphic after deletions "
          f"({n_before - len(snp_truth)} lost to deleted sequence)",
          file=sys.stderr)

    # ---- 5. truth tables -----------------------------------------------
    with open(f"{args.out_prefix}.truth_snps.tsv", "w") as fh:
        cols = ["pos", "ref", "alt", "ac", "af", "carriers"]
        fh.write("\t".join(cols) + "\n")
        for r in snp_truth:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    with open(f"{args.out_prefix}.truth_svs.tsv", "w") as fh:
        cols = ["sv_id", "type", "start", "end", "len", "ac", "af",
                "gnomad_af", "size_bin", "synthetic", "carriers"]
        fh.write("\t".join(cols) + "\n")
        for r in sv_truth:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    stats = {
        "contig": contig,
        "window_len": L,
        "n_hap": args.n_hap,
        "seed": args.seed,
        "mode": args.mode,
        "n_trees": int(ts.num_trees),
        "n_snps_polymorphic": len(snp_truth),
        "n_snps_lost_to_deletions": n_before - len(snp_truth),
        "n_svs_polymorphic": len(sv_truth),
        "sv_type_counts": {t: sum(1 for r in sv_truth if r["type"] == t)
                           for t in ("DEL", "INV", "DUP")},
        # truth summary stats from the tree sequence (SNP-based, exact)
        "true_diversity_pi": float(ts.diversity()),
        "true_tajimas_D": float(ts.Tajimas_D()),
        "n_backbone": backbone.count("N"),
    }
    with open(f"{args.out_prefix}.truth_stats.json", "w") as fh:
        json.dump(stats, fh, indent=2)
    print(json.dumps(stats, indent=2), file=sys.stderr)


if __name__ == "__main__":
    main()
