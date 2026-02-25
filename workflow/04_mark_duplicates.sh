#!/usr/bin/env bash
set -euo pipefail

IN_DIR="results/readgroups"
OUT_DIR="results/markduplicates"
LOG_DIR="logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

for SAMPLE in father mother proband; do
  
  IN_BAM="${IN_DIR}/${SAMPLE}.sorted.rg.bam"
  OUT_BAM="${OUT_DIR}/${SAMPLE}.sorted.rg.dedup.bam"
  METRICS="${OUT_DIR}/${SAMPLE}.dedup.metrics.txt"
  LOG="${LOG_DIR}/${SAMPLE}.04_markdup.log"

  echo "[04] MarkDuplicates: ${SAMPLE}"
  [[ -f "${IN_BAM}" ]] || { echo "ERROR: missing ${IN_BAM}" >&2; exit 1; }

  picard MarkDuplicates \
    I="${IN_BAM}" \
    O="${OUT_BAM}" \
    M="${METRICS}" \
    CREATE_INDEX=true \
    > "${LOG}" 2>&1

  [[ -f "${OUT_BAM}" && -f "${OUT_BAM}.bai" ]] || { echo "ERROR: output/index missing for ${SAMPLE}" >&2; exit 1; }

  samtools flagstat "${OUT_BAM}" > "${OUT_DIR}/${SAMPLE}.dedup.flagstat.txt"
done

echo "Step [04] completed."
