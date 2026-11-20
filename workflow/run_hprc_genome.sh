#!/bin/bash
# run_hprc_genome.sh
#
# Runs the fragmentation scan on every autosome of the HPRC release 2
# Minigraph-Cactus GRCh38 wave VCF, under the conservative rule (strict size,
# one graph site). Each chromosome is fetched once through the tabix index and
# kept, so the run resumes where it stopped.
#
#   nohup bash run_hprc_genome.sh > hprc_genome.log 2>&1 &
#
# Then combine:
#   python3 analysis/hprc_genome_summary.py --dir out/hprc/genome --out out/hprc/genome_summary
#
# Environment: DATA (default /mnt/d/pangenome_data/hprc_wave_by_chrom),
#              REPO (default ~/graphsv), CHROMS (default chr1..chr22)

set -u
source ~/miniforge3/etc/profile.d/conda.sh && conda activate pg

B="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/release2/minigraph-cactus/v2.1/hprc-v2.1-mc-grch38"
WAVE="$B/hprc-v2.1-mc-grch38.wave.vcf.gz"
DATA=${DATA:-/mnt/d/pangenome_data/hprc_wave_by_chrom}
REPO=${REPO:-$HOME/graphsv}
CHROMS=${CHROMS:-$(seq -f "chr%g" 1 22)}
OUT=$REPO/out/hprc/genome
mkdir -p "$DATA" "$OUT"

exec 9>/tmp/run_hprc_genome.lock
flock -n 9 || { echo "already running"; exit 1; }

echo "start $(date '+%F %T')"
for C in $CHROMS; do
  V=$DATA/hprc-v2.1-mc-grch38.wave.${C}.vcf.gz
  if [ ! -s "$V" ]; then
    bcftools view -r "$C" -Oz -o "$V.part" "$WAVE" && mv "$V.part" "$V" \
      || { echo "FETCH FAILED $C"; rm -f "$V.part"; continue; }
  fi
  [ -s "$OUT/${C}.summary.txt" ] && { echo "skip $C"; continue; }
  python3 "$REPO/analysis/hprc_fragmentation.py" --vcf "$V" \
    --size-min 10 --size-frac 0.01 --same-site-only \
    --out "$OUT/${C}" > "$OUT/${C}.log" 2>&1 \
    || { echo "SCAN FAILED $C"; continue; }
  echo "$(date +%H:%M) done $C  $(du -h "$V" | cut -f1)  $(sed -n 3p "$OUT/${C}.summary.txt")"
done
echo "COMPLETE $(date '+%F %T')"
