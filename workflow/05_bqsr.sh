#!/usr/bin/env bash
set -euo pipefail

IN_DIR="results/markduplicates"
OUT_DIR="results/bqsr"
LOG_DIR="logs"

REF="ref/GRCh38.fa"
DBSNP="ref/known_sites/dbsnp_146.hg38.vcf.gz"
MILLS="ref/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# Basic reference checks (common failure points)
[[ -f "${REF}" ]] || { echo "ERROR: missing reference FASTA: ${REF}" >&2; exit 1; }
[[ -f "${REF}.fai" ]] || { echo "ERROR: missing reference index (.fai): ${REF}.fai" >&2; exit 1; }
[[ -f "${DBSNP}" ]] || { echo "ERROR: missing known-sites VCF: ${DBSNP}" >&2; exit 1; }
[[ -f "${MILLS}" ]] || { echo "ERROR: missing known-sites VCF: ${MILLS}" >&2; exit 1; }

# Known-sites VCFs are typically bgzip-compressed and tabix-indexed
[[ -f "${DBSNP}.tbi" ]] || echo "WARN: missing VCF index: ${DBSNP}.tbi (GATK may fail)"
[[ -f "${MILLS}.tbi" ]] || echo "WARN: missing VCF index: ${MILLS}.tbi (GATK may fail)"


for SAMPLE in father mother proband; do

  IN_BAM="$IN_DIR/$SAMPLE.sorted.rg.dedup.bam"
  RECAL_TABLE="$OUT_DIR/$SAMPLE.recal_data.table"
  OUT_BAM="$OUT_DIR/$SAMPLE.sorted.rg.dedup.recal.bam"
  LOG="$LOG_DIR/$SAMPLE.05_bqsr.log"

  echo "[05] BQSR: ${SAMPLE}"
  [[ -f "${IN_BAM}" ]] || { echo "ERROR: missing input BAM: ${IN_BAM}" >&2; exit 1; }

  # 1) Build recalibration model
  gatk BaseRecalibrator \
    -R "$REF" \
    -I "$IN_BAM" \
    --known-sites "$DBSNP" \
    --known-sites "$MILLS" \
    -O "$RECAL_TABLE" \
    > "$LOG" 2>&1

  [[ -s "${RECAL_TABLE}" ]] || { echo "ERROR: recal table not created or empty: ${RECAL_TABLE}" >&2; exit 1; }

  # 2) Apply recalibration
  gatk ApplyBQSR \
    -R "$REF" \
    -I "$IN_BAM" \
    --bqsr-recal-file "$RECAL_TABLE" \
    -O "$OUT_BAM" \
    >> "$LOG" 2>&1

  [[ -f "${OUT_BAM}" ]] || { echo "ERROR: output BAM missing: ${OUT_BAM}" >&2; exit 1; }


  # Index output BAM (needed for downstream GATK steps)
  samtools index "${OUT_BAM}"
  [[ -f "${OUT_BAM}.bai" ]] || { echo "ERROR: output BAM index missing: ${OUT_BAM}.bai" >&2; exit 1; }
done

echo "[05] BQSR — completed."
