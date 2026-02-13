#!/usr/bin/env bash
set -euo pipefail

IN_DIR="results/13B_gnomad_filtered"
CLINVAR_VCF="data/clinvar/clinvar.vcf.gz"
OUT_DIR="results/14_clinvar"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"
LOG="$LOG_DIR/14_clinvar_join.log"
: > "$LOG"

echo "=== Step 14: ClinVar join (no heredocs) ===" | tee -a "$LOG"
echo "Input rare dir: $IN_DIR" | tee -a "$LOG"
echo "ClinVar VCF:    $CLINVAR_VCF" | tee -a "$LOG"
echo "Output dir:     $OUT_DIR" | tee -a "$LOG"

if [[ ! -f "$CLINVAR_VCF" ]]; then
  echo "ERROR: ClinVar VCF not found: $CLINVAR_VCF" | tee -a "$LOG"
  echo "You need to download ClinVar GRCh38 VCF into results/13_external/clinvar.vcf.gz" | tee -a "$LOG"
  exit 1
fi

if [[ ! -f "${CLINVAR_VCF}.tbi" ]]; then
  echo "[0] Index ClinVar VCF" | tee -a "$LOG"
  tabix -p vcf "$CLINVAR_VCF"
fi

CHR_PREFIX=""
if bcftools view -h "$CLINVAR_VCF" | grep -q '##contig=<ID=chr'; then
  CHR_PREFIX="chr"
fi
echo "[1] ClinVar contig naming: prefix='${CHR_PREFIX}'" | tee -a "$LOG"

REGIONS="$OUT_DIR/clinvar.regions.tsv"
SUBSET_VCF="$OUT_DIR/clinvar.subset.vcf.gz"
LOOKUP_TSV="$OUT_DIR/clinvar.lookup.tsv"

echo "[2] Build regions from rare TSVs" | tee -a "$LOG"
python3 workflow/14A_build_regions.py "$IN_DIR" "$REGIONS" "$CHR_PREFIX" | tee -a "$LOG"

echo "[3] Subset ClinVar by regions" | tee -a "$LOG"
bcftools view -R "$REGIONS" -Oz -o "$SUBSET_VCF" "$CLINVAR_VCF"
tabix -p vcf "$SUBSET_VCF"
echo "    subset variants: $(bcftools view -H "$SUBSET_VCF" | wc -l)" | tee -a "$LOG"

echo "[4] Build lookup TSV" | tee -a "$LOG"
{
  echo -e "CHROM\tPOS\tREF\tALT\tCLNSIG\tCLNREVSTAT\tCLNDN\tCLNDISDB\tCLNVC\tRS"
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\t%INFO/CLNDISDB\t%INFO/CLNVC\t%ID\n' "$SUBSET_VCF"
} > "$LOOKUP_TSV"
echo "    lookup lines (incl header): $(wc -l < "$LOOKUP_TSV")" | tee -a "$LOG"

echo "[5] Join ClinVar into rare TSVs" | tee -a "$LOG"
python3 workflow/14B_join_clinvar.py "$IN_DIR" "$LOOKUP_TSV" "$OUT_DIR" | tee -a "$LOG"

echo "Step 14 complete. Outputs: $OUT_DIR/*_clinvar.tsv" | tee -a "$LOG"
