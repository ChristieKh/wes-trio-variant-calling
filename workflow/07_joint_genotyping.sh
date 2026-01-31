#!/usr/bin/env bash
set -e

REF="ref/GRCh38.fa"
VCF_DIR="results/vcf"
LOG_DIR="logs"

mkdir -p "$LOG_DIR"

echo "=== Step 07: Joint genotyping (CombineGVCFs) ==="

gatk CombineGVCFs \
  -R "$REF" \
  -V "$VCF_DIR/father.g.vcf.gz" \
  -V "$VCF_DIR/mother.g.vcf.gz" \
  -V "$VCF_DIR/proband.g.vcf.gz" \
  -O "$VCF_DIR/trio.combined.g.vcf.gz" \
  > "$LOG_DIR/07_combine_gvcfs.log" 2>&1

echo "=== Step 07: Joint genotyping (GenotypeGVCFs) ==="

gatk GenotypeGVCFs \
  -R "$REF" \
  -V "$VCF_DIR/trio.combined.g.vcf.gz" \
  -O "$VCF_DIR/trio.raw.vcf.gz" \
  > "$LOG_DIR/07_genotype_gvcfs.log" 2>&1

echo "Step [07] Joint genotyping done."
