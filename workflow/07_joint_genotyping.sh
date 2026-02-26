#!/usr/bin/env bash
set -euo pipefail

REF="ref/GRCh38.fa"
VCF_DIR="results/vcf"
LOG_DIR="logs"

mkdir -p "${LOG_DIR}"

# --- Reference checks ---
[[ -f "${REF}" ]] || { echo "ERROR: missing reference FASTA: ${REF}" >&2; exit 1; }
[[ -f "${REF}.fai" ]] || { echo "ERROR: missing reference index (.fai)" >&2; exit 1; }

# --- Input gVCF checks ---
for SAMPLE in father mother proband; do
  GVCF="${VCF_DIR}/${SAMPLE}.g.vcf.gz"
  [[ -f "${GVCF}" ]] || { echo "ERROR: missing gVCF: ${GVCF}" >&2; exit 1; }
  [[ -f "${GVCF}.tbi" ]] || { echo "WARN: missing index for ${GVCF}" >&2; }
done

echo "Step 07: Joint genotyping (CombineGVCFs)"

gatk CombineGVCFs \
  -R "$REF" \
  -V "$VCF_DIR/father.g.vcf.gz" \
  -V "$VCF_DIR/mother.g.vcf.gz" \
  -V "$VCF_DIR/proband.g.vcf.gz" \
  -O "$VCF_DIR/trio.combined.g.vcf.gz" \
  > "$LOG_DIR/07_combine_gvcfs.log" 2>&1
  
[[ -s "${VCF_DIR}/trio.combined.g.vcf.gz" ]] || { echo "ERROR: combined gVCF not created" >&2; exit 1; }

echo "Step 07: Joint genotyping (GenotypeGVCFs)"

gatk GenotypeGVCFs \
  -R "$REF" \
  -V "$VCF_DIR/trio.combined.g.vcf.gz" \
  -O "$VCF_DIR/trio.raw.vcf.gz" \
  > "$LOG_DIR/07_genotype_gvcfs.log" 2>&1

  [[ -s "${VCF_DIR}/trio.raw.vcf.gz" ]] || { echo "ERROR: raw VCF not created" >&2; exit 1; }


echo "Step [07] Joint genotyping completed."
