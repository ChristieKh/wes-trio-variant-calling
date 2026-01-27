#!/usr/bin/env bash
set -euo pipefail

# EDIT THESE PATHS
IN_DIR="results/bam"  
OUT_DIR="results/readgroups"
LOG_DIR="logs"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"


SAMPLES=("father" "mother" "proband")

for SAMPLE in "${SAMPLES[@]}"; do
  IN_BAM="${IN_DIR}/${SAMPLE}.sorted.bam"
  OUT_BAM="${OUT_DIR}/${SAMPLE}.sorted.rg.bam"
  LOG="${LOG_DIR}/${SAMPLE}.03_add_rg.log"

  echo "[03] Add RG: ${SAMPLE}"
  [[ -f "${IN_BAM}" ]] || { echo "ERROR: missing ${IN_BAM}" >&2; exit 1; }

  picard AddOrReplaceReadGroups \
    I="${IN_BAM}" \
    O="${OUT_BAM}" \
    RGID="${SAMPLE}_L001" \
    RGLB="lib1" \
    RGPL="ILLUMINA" \
    RGPU="${SAMPLE}_L001" \
    RGSM="${SAMPLE}" \
    CREATE_INDEX=true \
    > "${LOG}" 2>&1

  # quick sanity check: RG exists
  samtools view -H "${OUT_BAM}" | grep -q "^@RG" || { echo "ERROR: no @RG in ${OUT_BAM}" >&2; exit 1; }
done

echo "Step [03] Is Done."
