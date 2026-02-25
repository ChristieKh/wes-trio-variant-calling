#!/usr/bin/env bash
set -euo pipefail

REF="ref/GRCh38.fa"
IN_DIR="results/bqsr"
OUT_DIR="results/vcf"
LOG_DIR="logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

# Reference checks
[[ -f "${REF}" ]] || { echo "ERROR: missing reference FASTA: ${REF}" >&2; exit 1; }
[[ -f "${REF}.fai" ]] || { echo "ERROR: missing reference index (.fai): ${REF}.fai" >&2; exit 1; }

for SAMPLE in father mother proband; do
 
  IN_BAM="${IN_DIR}/${SAMPLE}.sorted.rg.dedup.recal.bam"
  OUT_GVCF="${OUT_DIR}/${SAMPLE}.g.vcf.gz"
  LOG="${LOG_DIR}/${SAMPLE}.06_haplotypecaller_gvcf.log"

  echo "[06] HaplotypeCaller (gVCF): ${SAMPLE}"
  [[ -f "${IN_BAM}" ]] || { echo "ERROR: missing input BAM: ${IN_BAM}" >&2; exit 1; }

  gatk HaplotypeCaller \
    -R "$REF" \
    -I "$IN_BAM" \
    -O "$OUT_GVCF" \
    -ERC GVCF \
    > "$LOG" 2>&1

  # Output sanity checks
  [[ -s "${OUT_GVCF}" ]] || { echo "ERROR: gVCF not created or empty: ${OUT_GVCF}" >&2; exit 1; }
done

echo "Step [06] HaplotypeCaller (gVCF) completed."
