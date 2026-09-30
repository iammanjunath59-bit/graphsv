#!/usr/bin/env bash
# verify_all.sh
#
# Re-derives every numeric claim in Results_draft.md from the files on disk,
# and prints the expected value beside each so they can be compared by eye.
#
# Purpose: the draft numbers were produced during an interactive session.
# This script exists so they can be checked independently, in one run, without
# trusting any earlier output. Anything that prints MISMATCH or a blank value
# should be treated as unverified until resolved.
#
#   bash verify_all.sh 2>&1 | tee verification.txt
#
# Run from the repository root, or set PANGENOME_DIR to a tree laid out
# the same way. Graphs are archived rather than committed, so the runtime
# section is skipped unless GFA points at one.

set -u
cd "${PANGENOME_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}" || exit 1

hdr () { printf '\n\033[1m=== %s ===\033[0m\n' "$1"; }
exp () { printf '    expected: %s\n' "$1"; }

hdr "0. environment"
which vg bcftools vcfbub vcfwave python3 2>&1 | sed 's/^/  /'
python3 -c "import statsmodels, pandas, numpy, matplotlib; print('  statsmodels', statsmodels.__version__)"

# ---------------------------------------------------------------------------
hdr "1. benchmark size"
echo -n "  chr20 stratified replicates: "; ls results/sweep54d/*.truth_svs.tsv | wc -l
exp "108"
echo -n "  chr20 SVs total:             "
cat results/sweep54d/*.truth_svs.tsv | grep -vc '^sv_id'
exp "2592"
echo -n "  chr20 distinct windows:      "
ls results/sweep54d/*.truth_svs.tsv | sed 's/_s[0-9]*\.truth_svs\.tsv//' | sort -u | wc -l
exp "54"
echo -n "  chr16 replicates:            "; ls results/sweep16/*.truth_svs.tsv | wc -l
exp "61"
echo -n "  chr16 SVs total:             "
cat results/sweep16/*.truth_svs.tsv | grep -vc '^sv_id'
exp "1464  (2592 + 1464 = 4056 total)"
echo -n "  SVs per replicate (chr20):   "
grep -vc '^sv_id' results/sweep54d/chr20_10000000_s1.truth_svs.tsv
exp "24"
echo -n "  size bins present:           "
cat results/sweep54d/*.truth_svs.tsv | grep -v '^sv_id' | cut -f9 | sort -u | tr '\n' ' '; echo
echo -n "  SV types present:            "
cat results/sweep54d/*.truth_svs.tsv | grep -v '^sv_id' | cut -f2 | sort | uniq -c | tr '\n' ' '; echo

# ---------------------------------------------------------------------------
hdr "2. PGGB VCF route (wave) -- overall and by size bin"
for f in results/sweep54d/*.vcfcmp.summary.tsv; do tail -n 1 "$f"; done | \
  awk -F'\t' '{n++; p+=$10; c+=$11} END {printf "  windows=%d  ac_present=%.4f  clean=%.4f\n", n, p/n, c/n}'
exp "108  0.8646  0.8646"
cat results/sweep54d/*.vcfcmp.svs.tsv | grep -v '^sv_id' | \
  awk -F'\t' '{t[$4]++; if($14==1)k[$4]++} END {for(b in t) printf "  %7s n=%4d %.3f\n", b, t[b], k[b]/t[b]}' | sort -n
exp "100:0.988  300:0.977  1000:0.988  3000:0.898  10000:0.782  30000:0.553"

# ---------------------------------------------------------------------------
hdr "3. Stage comparison (Supplementary Fig. 1)"
echo -n "  raw comparison files:  "; ls results/stages/*.raw.cmp.svs.tsv 2>/dev/null | wc -l
echo -n "  bub comparison files:  "; ls results/stages/*.bub.cmp.svs.tsv 2>/dev/null | wc -l
exp "108 each"
for st in raw bub; do
  echo "  --- $st ---"
  cat results/stages/*.$st.cmp.svs.tsv 2>/dev/null | grep -v '^sv_id' | \
    awk -F'\t' '{n++; if($14==1){e++; if($8<=2)c++}}
      END {printf "    ac_present=%.4f  clean=%.4f  (n=%d)\n", e/n, c/n, n}'
done
exp "raw 0.8777 / 0.4664   bub 0.8773 / 0.4653"
echo "  clean by size bin:"
for st in raw bub; do
  printf "    %-4s " "$st"
  cat results/stages/*.$st.cmp.svs.tsv 2>/dev/null | grep -v '^sv_id' | \
    awk -F'\t' '{t[$4]++; if($14==1 && $8<=2)k[$4]++}
      END {n=split("100 300 1000 3000 10000 30000",o," ");
           for(i=1;i<=n;i++){b=o[i]; printf "%s:%.3f ", b, k[b]/t[b]} print ""}'
done
exp "raw 0.975 0.903 0.597 0.227 0.062 0.035"

# ---------------------------------------------------------------------------
hdr "4. PGGB graph route (graphsv)"
awk -F'\t' 'FNR>1 {n++; e+=$8} END {printf "  overall exact=%.4f  (n=%d SVs)\n", e/n, n}' \
  results/sweep54d/*.gsv.eval.tsv
exp "0.9923  n=2592"
cat results/sweep54d/*.gsv.eval.tsv | grep -v '^sv_id' | \
  awk -F'\t' '{t[$4]++; if($8==1)k[$4]++} END {for(b in t) printf "  %7s n=%4d %.4f\n", b, t[b], k[b]/t[b]}' | sort -n
exp "100:0.9722  300:0.9977  1000:1.0000  3000:1.0000  10000:0.9977  30000:0.9861"
echo -n "  100bp misses that are NO-CALL (empty acs_found): "
awk -F'\t' 'FNR>1 && $4==100 && $8==0 {n++; if($7=="")z++} END {print z"/"n}' results/sweep54d/*.gsv.eval.tsv
exp "12/12  -- all misses are no-calls, not miscounts"

# ---------------------------------------------------------------------------
hdr "5. Minigraph-Cactus"
echo -n "  windows with comparison: "; ls results/mc/*.mccmp.svs.tsv | wc -l
exp "54"
for f in results/mc/*.mccmp.summary.tsv; do tail -n 1 "$f"; done | \
  awk -F'\t' '{n++; p+=$10; c+=$11} END {printf "  VCF route: n=%d  ac_present=%.4f  clean=%.4f\n", n, p/n, c/n}'
exp "54  0.8750  0.8750"
awk -F'\t' 'FNR>1 {n++; e+=$8} END {printf "  graph route: exact=%.4f  (n=%d SVs)\n", e/n, n}' results/mc/*.mcgsv.eval.tsv
exp "0.9877  n=1296"
echo "  VCF route by bin:"
cat results/mc/*.mccmp.svs.tsv | grep -v '^sv_id' | \
  awk -F'\t' '{t[$4]++; if($14==1)k[$4]++} END {for(b in t) printf "    %7s n=%4d %.3f\n", b, t[b], k[b]/t[b]}' | sort -n
exp "100:0.991 300:0.986 1000:0.991 3000:0.912 10000:0.801 30000:0.569"
echo "  graph route by bin:"
cat results/mc/*.mcgsv.eval.tsv | grep -v '^sv_id' | \
  awk -F'\t' '{t[$4]++; if($8==1)k[$4]++} END {for(b in t) printf "    %7s n=%4d %.3f\n", b, t[b], k[b]/t[b]}' | sort -n
exp "100:0.977 300:0.995 1000:0.995 3000:0.995 10000:0.995 30000:0.968"
echo "  panel identity check (MC vs PGGB truth tables):"
mis=0; tot=0
for f in results/mc/*.truth_svs.tsv; do
  b=$(basename "$f" .truth_svs.tsv); r=results/sweep54d/${b}.truth_svs.tsv
  [ -f "$r" ] || continue
  tot=$((tot+1)); cmp -s "$f" "$r" || { mis=$((mis+1)); echo "    MISMATCH $b"; }
done
echo "    identical panels: $((tot-mis))/$tot"
exp "all identical -- builder comparison is only valid if this holds"

# ---------------------------------------------------------------------------
hdr "6. Regression -- chr20 and chr16"
python3 analysis/analyse_sweep.py --glob 'results/sweep54d/*.vcfcmp.svs.tsv' \
  --windows data/windows_chr20.tsv --out /tmp/vf20 2>/dev/null | \
  sed -n '/key coefficients/,$p'
exp "logsize_c beta=-2.173 SE=0.139 ; segdup10 beta=-0.419 SE=0.063"
python3 analysis/analyse_sweep.py --glob 'results/sweep16/*.vcfcmp.svs.tsv' \
  --windows data/windows_chr16.tsv --out /tmp/vf16 2>/dev/null | \
  sed -n '/key coefficients/,$p'
exp "logsize_c beta=-2.183 SE=0.252 ; segdup10 beta=-0.281 SE=0.061"
echo "  95% CIs (cluster-robust model 2, chr20):"
awk '/2\. Binomial GLM, cluster-robust/,/^====.*$/' /tmp/vf20.model.txt | grep -E "logsize_c|segdup10"
exp "logsize_c CI -2.446 to -1.900 ; segdup10 CI -0.544 to -0.295"
echo "  interaction term (chr20 model 3):"
grep "logsize_c:segdup10" /tmp/vf20.model.txt
exp "-0.1229  SE 0.206  p=0.551  (no interaction)"
echo "  MIXED MODEL CONVERGENCE -- must read 'No':"
grep -i "Converged" /tmp/vf20.model.txt
exp "Converged: No  -- so the LMM is NOT reportable; use the cluster-robust GLM"
echo "  segdup-stratified cells (chr20):"
head -14 /tmp/vf20.model.txt
exp "30000 low 0.571 / high 0.375 ; 100 low 0.987 / high 1.000"

# ---------------------------------------------------------------------------
hdr "7. Segdup robustness (permutation)"
python3 analysis/robustness_segdup.py --perSV /tmp/vf20.perSV.tsv --out /tmp/vf20_robustness 2>&1 | \
  sed -n '/permutation of segdup/,$p'
exp "observed -0.4193 ; null mean +0.013 sd 0.179 ; one-sided p = 5.0e-04"

# ---------------------------------------------------------------------------
hdr "8. Merge audit"
python3 analysis/audit_merges.py --dir results/sweep54d 2>&1 | sed -n '1,25p'
exp "FALSE FUSION 52/2512 = 0.0207 ; RESIDUAL SPLIT 2/2592 = 0.0008"
exp "79 truth pairs share a carrier set, 63 within 50 kb  (opportunity exists)"
echo "  NOTE: the value 0.017 quoted in older notes does NOT appear here."

# ---------------------------------------------------------------------------
hdr "9. Realistic mode (sweep54r)"
for f in results/sweep54r/*.vcfcmp.summary.tsv; do tail -n 1 "$f"; done | \
  awk -F'\t' '{n++; p+=$10; c+=$11} END {printf "  n=%d ac_present=%.4f clean=%.4f\n", n, p/n, c/n}'
exp "108  clean=0.9847"
echo -n "  mode field: "; tail -n 1 results/sweep54r/chr20_10000000_s1.vcfcmp.summary.tsv | cut -f2
exp "realistic   (sweep54d should read 'stratified')"
echo -n "  sweep54d mode field: "; tail -n 1 results/sweep54d/chr20_10000000_s1.vcfcmp.summary.tsv | cut -f2

# ---------------------------------------------------------------------------
hdr "10. Runtime"
G=${GFA:-results/sweep54d/chr20_10000000_s1.final.gfa}
if [ ! -s "$G" ]; then
  echo "  SKIPPED: no GFA at $G."
  echo "  Graphs are archived rather than committed; set GFA to one,"
  echo "  or point PANGENOME_DIR at a tree that contains them."
else
R=$(vg paths -Lx "$G" | grep '^anc#')
echo "  graphsv:"
/usr/bin/time -f "    %e s  %M KB" python3 src/graphsv.py --gfa "$G" \
  --ref-prefix anc --truth results/sweep54d/chr20_10000000_s1.truth_svs.tsv \
  --min-seg 50 --out /tmp/vt.tsv >/dev/null
echo "  vg deconstruct + vcfbub + vcfwave:"
/usr/bin/time -f "    %e s  %M KB" bash -c \
  "vg deconstruct -p '$R' -a '$G' > /tmp/vr.vcf 2>/dev/null; \
   vcfbub -l 0 -a 100000 -i /tmp/vr.vcf > /tmp/vb.vcf 2>/dev/null; \
   vcfwave -I 1000 /tmp/vb.vcf > /tmp/vw.vcf 2>/dev/null"
echo -n "  wave records produced: "; grep -vc '^#' /tmp/vw.vcf
fi
exp "0.08-0.13 s vs 6.07-6.13 s ; 1751 records -- if records=0 the chain did not run"

# ---------------------------------------------------------------------------
hdr "11. UNVERIFIED -- these two are NOT re-derived above"
cat <<'TXT'
  The Results draft states VIF = 1.0000 and GC p = 0.60 as evidence that the
  size and segdup effects are separable. Neither was reproduced in the session
  that produced the draft. Locate the script that computed them and run it, or
  cut both clauses from the text. Candidates:
TXT
grep -rln "vif\|variance_inflation\|gc_frac\|gc_content" analysis/*.py src/*.py workflow/*.py hprc/*.py 2>/dev/null | sed 's/^/    /'

# ---------------------------------------------------------------------------
hdr "12. HPRC case study (requires network)"
cat <<'TXT'
  Run only if re-checking the real-data section. Needs RAW/WAVE exported:

    B="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38"
    export RAW="$B/hprc-v2.1-mc-grch38.raw.vcf.gz"
    export WAVE="$B/hprc-v2.1-mc-grch38.wave.vcf.gz"
    python3 wave_dispersal.py --pos 25783205 --event 111144578_112737255

  expected: size 61,963 bp | 84 raw alleles | 216 carrier haplotypes | AN 461
            carrier median deleted bp 61,964 -> dispersal ratio 1.000
TXT

hdr "done"
echo "Compare each printed value against its 'expected' line."
echo "Any mismatch, blank, or error means that number is unverified."
