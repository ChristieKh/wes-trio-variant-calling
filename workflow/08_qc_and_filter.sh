#!/usr/bin/env bash
set -euo pipefail

VCF_IN="results/06_vcf/trio.raw.vcf.gz"
VCF_OUT="results/06_vcf/trio.filtered.vcf.gz"
LOG_DIR="logs"
QC_DIR="results/08_qc"

LOG="${LOG_DIR}/08_qc_filter.log"

mkdir -p "${QC_DIR}" "${LOG_DIR}"

# --- input checks ---
[[ -f "${VCF_IN}" ]] || { echo "ERROR: missing input VCF: ${VCF_IN}" >&2; exit 1; }
[[ -f "${VCF_IN}.tbi" || -f "${VCF_IN}.csi" ]] || echo "WARN: missing VCF index for ${VCF_IN} (bcftools may be slower)" >&2


echo "Step 08: QC summary" | tee -a "${LOG}"

echo "[1] Samples in VCF:" | tee -a "${LOG}"
bcftools query -l "$VCF_IN" | tee "$QC_DIR/trio.samples.txt" | tee -a "${LOG}"

# verify proband sample exists
grep -qx "proband" "${QC_DIR}/trio.samples.txt" || { echo "ERROR: sample 'proband' not found in VCF samples" >&2; exit 1; }

echo "[2] Total variant records:" | tee -a "${LOG}"
bcftools view -H "$VCF_IN" | wc -l | tee "$QC_DIR/trio.raw.count.txt" | tee -a "${LOG}"

echo "[3] Basic bcftools stats:" | tee -a "${LOG}"
bcftools stats "$VCF_IN" > "$QC_DIR/trio.raw.bcftools_stats.txt"

echo "Step 08: Filtering (proband DP/GQ)" | tee -a "${LOG}"

# 1) Build regions list (BED) from proband passing genotype-level thresholds
# Keep non-missing genotypes in proband with DP>=10 and GQ>=20
bcftools view -s proband -Ou "$VCF_IN" \
| bcftools filter -Ou -i 'GT!="mis" && FORMAT/DP>=10 && FORMAT/GQ>=20' \
| bcftools query -f '%CHROM\t%POS\n' \
| awk '{print $1"\t"($2-1)"\t"$2}' \
> "$QC_DIR/proband.pass.bed"

echo "[4] Passing positions in proband:" | tee -a "${LOG}"
wc -l "$QC_DIR/proband.pass.bed" | tee -a "${LOG}"

# 2) Apply to full trio VCF by region list
bcftools view -R "$QC_DIR/proband.pass.bed" -Oz -o "$VCF_OUT" "$VCF_IN"
bcftools index -t "$VCF_OUT"

echo "[5] Filtered variant records:" | tee -a "${LOG}"
bcftools view -H "$VCF_OUT" | wc -l | tee "$QC_DIR/trio.filtered.count.txt" | tee -a "${LOG}"

echo "Step [08] QC + filtering completed"
echo "Output: $VCF_OUT"
