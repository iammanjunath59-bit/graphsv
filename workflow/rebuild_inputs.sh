#!/bin/bash
# rebuild_inputs.sh
#
# Rebuilds the simulation inputs lost with ~/pangenome, using the exact
# commands recorded when they were first made, and checks every step against
# values recorded at the time.
#
#   stage 1  reference      chr20.fa from UCSC, indexed
#   stage 2  windows        54 backbones and 54 gnomAD SV beds, built with the
#                           original bcftools query; each bed's line count is
#                           compared with the count recorded on 16 August
#   stage 3  one panel      simulate chr20_10000000 seed 1 and compare its truth
#                           table byte for byte with the deposited one
#   stage 4  all panels     the same for all 108 replicates
#
# Nothing here builds graphs. If stage 4 passes, the regenerated panels are
# identical to the originals and the graph builds can follow.
#
# Usage:
#   bash rebuild_inputs.sh            run every stage, stopping at the first failure
#   bash rebuild_inputs.sh 3          run up to and including stage 3
#
# Environment:
#   REBUILD_DIR  working directory, default ~/pangenome_rebuild
#   REPO         repository, default ~/graphsv

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg

UNTIL=${1:-4}
R=${REBUILD_DIR:-$HOME/pangenome_rebuild}
REPO=${REPO:-$HOME/graphsv}
URL=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz
SIM=$REPO/workflow/simulate_pangenome_truth.py
TRUTH=$REPO/results/sweep54d

for need in "$SIM" "$REPO/data/eligible_windows.tsv"; do
  [ -e "$need" ] || { echo "FATAL: missing $need"; exit 2; }
done

mkdir -p "$R/data/windows" "$R/sweep54d"
cd "$R" || exit 2
cp -f "$REPO/data/eligible_windows.tsv" data/

# --------------------------------------------------------------- stage 1
echo "=== stage 1: reference ==="
if [ ! -s data/chr20.fa.fai ]; then
  wget -q -O data/chr20.fa.gz \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz \
    || { echo "FATAL: download of chr20.fa.gz failed"; exit 2; }
  gunzip -f data/chr20.fa.gz
  samtools faidx data/chr20.fa
fi
echo "  chr20 length: $(cut -f2 data/chr20.fa.fai | head -1)"
[ "$UNTIL" -ge 2 ] || exit 0

# --------------------------------------------------------------- stage 2
echo "=== stage 2: windows ==="
# Line counts recorded on 16 August when the per window beds were first built.
cat > data/expected_counts.tsv <<'EOF'
chr20_14000000	761
chr20_15000000	428
chr20_16000000	251
chr20_19000000	242
chr20_20000000	296
chr20_24000000	357
chr20_37000000	304
chr20_39000000	234
chr20_40000000	217
chr20_42000000	361
chr20_55000000	355
chr20_58000000	265
chr20_60000000	285
chr20_8000000	252
chr20_9000000	229
chr20_2000000	320
chr20_51000000	384
chr20_38000000	322
chr20_33000000	356
chr20_52000000	292
chr20_18000000	304
chr20_7000000	273
chr20_21000000	224
chr20_35000000	337
chr20_44000000	324
chr20_49000000	422
chr20_57000000	330
chr20_45000000	282
chr20_6000000	246
chr20_11000000	260
chr20_56000000	371
chr20_22000000	260
chr20_10000000	206
chr20_4000000	340
chr20_59000000	303
chr20_12000000	303
chr20_54000000	223
chr20_43000000	298
chr20_17000000	254
chr20_3000000	513
chr20_41000000	244
chr20_36000000	388
chr20_13000000	260
chr20_5000000	408
chr20_32000000	349
chr20_34000000	380
chr20_46000000	232
chr20_53000000	365
chr20_1000000	290
chr20_48000000	454
chr20_47000000	399
chr20_50000000	367
chr20_23000000	326
chr20_25000000	278
EOF

mismatch=0; n=0
while read -r chrom start end sd rep; do
  [ "$chrom" = "chrom" ] && continue
  W=${chrom}_${start}
  n=$((n+1))
  if [ ! -s "data/windows/${W}.fa" ]; then
    samtools faidx data/chr20.fa "${chrom}:$((start+1))-${end}" > "data/windows/${W}.fa"
  fi
  if [ ! -s "data/windows/${W}_svs.bed" ]; then
    bcftools query -r "${chrom}:${start}-${end}" \
      -i 'FILTER="PASS" && (SVTYPE="DEL" || SVTYPE="INV" || SVTYPE="DUP")' \
      -f '%CHROM\t%POS\t%INFO/END\t%INFO/SVTYPE\t%INFO/AF\n' "$URL" \
      > "data/windows/${W}_svs.bed"
  fi
  got=$(wc -l < "data/windows/${W}_svs.bed")
  want=$(awk -v w="$W" '$1==w{print $2}' data/expected_counts.tsv)
  if [ "$got" != "$want" ]; then
    echo "  MISMATCH $W  rebuilt $got  recorded $want"
    mismatch=$((mismatch+1))
  fi
done < data/eligible_windows.tsv
echo "  windows: $n   count mismatches: $mismatch"
if [ "$mismatch" -gt 0 ]; then
  echo "STOP: rebuilt SV pools differ from the originals; panels will not reproduce."
  exit 1
fi
[ "$UNTIL" -ge 3 ] || exit 0

# --------------------------------------------------------------- stage 3
simulate () {
  local W=$1 start=$2 SEED=$3
  local P=sweep54d/${W}_s${SEED}
  python "$SIM" \
    --backbone "data/windows/${W}.fa" --sv-bed "data/windows/${W}_svs.bed" \
    --window-start "$start" --n-hap 20 --seed "$SEED" --ne 20000 \
    --mode stratified --out-prefix "$P" >>sweep54d/sim.log 2>&1
}

echo "=== stage 3: one panel ==="
simulate chr20_10000000 10000000 1
if cmp -s sweep54d/chr20_10000000_s1.truth_svs.tsv "$TRUTH/chr20_10000000_s1.truth_svs.tsv"; then
  echo "  chr20_10000000_s1 truth table: IDENTICAL to the deposited one"
else
  echo "  chr20_10000000_s1 truth table: DIFFERS from the deposited one"
  diff <(head -6 sweep54d/chr20_10000000_s1.truth_svs.tsv) \
       <(head -6 "$TRUTH/chr20_10000000_s1.truth_svs.tsv")
  echo "STOP: the simulator does not reproduce the original panel."
  exit 1
fi
[ "$UNTIL" -ge 4 ] || exit 0

# --------------------------------------------------------------- stage 4
echo "=== stage 4: all 108 panels ==="
same=0; diffr=0; missing=0
while read -r chrom start end sd rep; do
  [ "$chrom" = "chrom" ] && continue
  W=${chrom}_${start}
  for SEED in 1 2; do
    P=sweep54d/${W}_s${SEED}
    [ -s "${P}.truth_svs.tsv" ] || simulate "$W" "$start" "$SEED"
    if [ ! -s "$TRUTH/${W}_s${SEED}.truth_svs.tsv" ]; then
      missing=$((missing+1))
    elif cmp -s "${P}.truth_svs.tsv" "$TRUTH/${W}_s${SEED}.truth_svs.tsv"; then
      same=$((same+1))
    else
      diffr=$((diffr+1)); echo "  DIFFERS ${W}_s${SEED}"
    fi
  done
done < data/eligible_windows.tsv
echo "  identical: $same   different: $diffr   no deposited table: $missing"
[ "$diffr" -eq 0 ] && [ "$same" -eq 108 ] \
  && echo "ALL 108 PANELS REPRODUCE EXACTLY. Graph builds can proceed." \
  || echo "Not all panels reproduce; do not build graphs until this is understood."
