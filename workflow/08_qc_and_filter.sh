#!/usr/bin/env bash
set -e

VCF_IN="results/vcf/trio.raw.vcf.gz"
VCF_OUT="results/vcf/trio.filtered.vcf.gz"
LOG_DIR="logs"
QC_DIR="results/qc"

mkdir -p "$LOG_DIR" "$QC_DIR"

echo "=== Step 08: QC summary ===" | tee "$LOG_DIR/08_qc_filter.log"

echo "[1] Samples in VCF:" | tee -a "$LOG_DIR/08_qc_filter.log"
bcftools query -l "$VCF_IN" | tee "$QC_DIR/trio.samples.txt" | tee -a "$LOG_DIR/08_qc_filter.log"

echo "[2] Total variant records:" | tee -a "$LOG_DIR/08_qc_filter.log"
bcftools view -H "$VCF_IN" | wc -l | tee "$QC_DIR/trio.raw.count.txt" | tee -a "$LOG_DIR/08_qc_filter.log"

echo "[3] Basic bcftools stats:" | tee -a "$LOG_DIR/08_qc_filter.log"
bcftools stats "$VCF_IN" > "$QC_DIR/trio.raw.bcftools_stats.txt"

echo "=== Step 08: Filtering (proband DP/GQ) ===" | tee -a "$LOG_DIR/08_qc_filter.log"

# 1) build regions list from proband only
bcftools view -s proband -Ou "$VCF_IN" \
| bcftools filter -Ou -i 'GT!="mis" && FORMAT/DP>=10 && FORMAT/GQ>=20' \
| bcftools query -f '%CHROM\t%POS\n' \
| awk '{print $1"\t"($2-1)"\t"$2}' \
> "$QC_DIR/proband.pass.bed"

echo "[4] Passing positions in proband:" | tee -a "$LOG_DIR/08_qc_filter.log"
wc -l "$QC_DIR/proband.pass.bed" | tee -a "$LOG_DIR/08_qc_filter.log"

# 2) apply to trio VCF
bcftools view -R "$QC_DIR/proband.pass.bed" -Oz -o "$VCF_OUT" "$VCF_IN"
bcftools index -t "$VCF_OUT"

echo "[5] Filtered variant records:" | tee -a "$LOG_DIR/08_qc_filter.log"
bcftools view -H "$VCF_OUT" | wc -l | tee "$QC_DIR/trio.filtered.count.txt" | tee -a "$LOG_DIR/08_qc_filter.log"

echo "Step [08] QC + filtering done."
echo "Output: $VCF_OUT"
