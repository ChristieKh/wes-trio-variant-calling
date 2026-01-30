#!/usr/bin/env bash
set -e

REF="ref/GRCh38.fa"
IN_DIR="results/bqsr"
OUT_DIR="results/vcf"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

for SAMPLE in father mother proband; do
  echo "=== HaplotypeCaller (gVCF) for $SAMPLE ==="

  gatk HaplotypeCaller \
    -R "$REF" \
    -I "$IN_DIR/${SAMPLE}.sorted.rg.dedup.recal.bam" \
    -O "$OUT_DIR/${SAMPLE}.g.vcf.gz" \
    -ERC GVCF \
    > "$LOG_DIR/${SAMPLE}.06_hc.log" 2>&1
done

echo "Step [06] HaplotypeCaller (gVCF) done."
