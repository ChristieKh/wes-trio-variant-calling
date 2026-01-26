#!/bin/bash
set -euo pipefail

REF="ref/GRCh38.fa"
FASTQ_DIR="data/fastq"
BAM_DIR="results/bam"
BAM_QC_DIR="results/bam/qc"
THREADS=8

mkdir -p "${BAM_DIR}" "${QC_DIR}"


tail -n +2 samples.tsv | while IFS=$'\t' read -r SAMPLE_ID ROLE SRA; do
  echo "==> Aligning ${SAMPLE_ID} (${ROLE}) [${SRA}]"

  bwa mem -t "${THREADS}" "${REF}" \
    "${FASTQ_DIR}/${SRA}_1.fastq.gz" \
    "${FASTQ_DIR}/${SRA}_2.fastq.gz" | \
  samtools sort -o "${BAM_DIR}/${SAMPLE_ID}.sorted.bam"

  samtools index "${BAM_DIR}/${SAMPLE_ID}.sorted.bam"
  samtools flagstat "${BAM_DIR}/${SAMPLE_ID}.sorted.bam" > "${BAM_QC_DIR}/${SAMPLE_ID}.flagstat.txt"
done

