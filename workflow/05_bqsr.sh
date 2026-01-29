#!/usr/bin/env bash
set -e

IN_DIR="results/markduplicates"
OUT_DIR="results/bqsr"
LOG_DIR="logs"

REF="ref/GRCh38.fa"
DBSNP="ref/known_sites/dbsnp_146.hg38.vcf.gz"
MILLS="ref/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

mkdir -p "$OUT_DIR" "$LOG_DIR"

for SAMPLE in father mother proband; do
  echo "=== BQSR for $SAMPLE ==="

  IN_BAM="$IN_DIR/$SAMPLE.sorted.rg.dedup.bam"
  RECAL_TABLE="$OUT_DIR/$SAMPLE.recal_data.table"
  OUT_BAM="$OUT_DIR/$SAMPLE.sorted.rg.dedup.recal.bam"
  LOG="$LOG_DIR/$SAMPLE.05_bqsr.log"

  # 1) Build recalibration model
  gatk BaseRecalibrator \
    -R "$REF" \
    -I "$IN_BAM" \
    --known-sites "$DBSNP" \
    --known-sites "$MILLS" \
    -O "$RECAL_TABLE" \
    > "$LOG" 2>&1

  # 2) Apply recalibration
  gatk ApplyBQSR \
    -R "$REF" \
    -I "$IN_BAM" \
    --bqsr-recal-file "$RECAL_TABLE" \
    -O "$OUT_BAM" \
    >> "$LOG" 2>&1

  # Index output BAM 
  samtools index "$OUT_BAM"
done

echo "Step [05] Is Done."
