#!/usr/bin/env bash
set -euo pipefail

IN_DIR="results/13B_gnomad_filtered"
CLINVAR_VCF="data/clinvar/clinvar.vcf.gz"
OUT_DIR="results/14_clinvar"
LOG_DIR="logs"
PATTERN="${PATTERN:-*rare*.tsv}"   # keeps old behavior if nothing is passed

# Parse optional CLI args (override defaults)
while [[ $# -gt 0 ]]; do
  case "$1" in
    --in-dir) IN_DIR="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --clinvar-vcf) CLINVAR_VCF="$2"; shift 2 ;;
    --log-dir) LOG_DIR="$2"; shift 2 ;;
    --pattern) PATTERN="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 [--in-dir DIR] [--out-dir DIR] [--clinvar-vcf FILE] [--pattern GLOB] [--log-dir DIR]"
      exit 0
      ;;
    *) echo "ERROR: unknown argument: $1" >&2; exit 1 ;;
  esac
done

mkdir -p "$OUT_DIR" "$LOG_DIR"
LOG="$LOG_DIR/14_clinvar_join.log"
: > "$LOG"

echo "=== Step 14: ClinVar join ===" | tee -a "$LOG"
echo "Input dir:       $IN_DIR" | tee -a "$LOG"
echo "Input pattern:   $PATTERN" | tee -a "$LOG"
echo "ClinVar VCF:     $CLINVAR_VCF" | tee -a "$LOG"
echo "Output dir:      $OUT_DIR" | tee -a "$LOG"

if [[ ! -f "$CLINVAR_VCF" ]]; then
  echo "ERROR: ClinVar VCF not found: $CLINVAR_VCF" | tee -a "$LOG"
  exit 1
fi

if [[ ! -f "${CLINVAR_VCF}.tbi" ]]; then
  echo "[0] Index ClinVar VCF" | tee -a "$LOG"
  tabix -p vcf "$CLINVAR_VCF"
fi

# Sanity: check that pattern matches at least one file
shopt -s nullglob
TSVS=( "$IN_DIR"/$PATTERN )
shopt -u nullglob
if [[ ${#TSVS[@]} -eq 0 ]]; then
  echo "ERROR: no TSVs matched: ${IN_DIR}/${PATTERN}" | tee -a "$LOG"
  exit 1
fi

CHR_PREFIX=""
if bcftools view -h "$CLINVAR_VCF" | grep -q '##contig=<ID=chr'; then
  CHR_PREFIX="chr"
fi
echo "[1] ClinVar contig naming: prefix='${CHR_PREFIX}'" | tee -a "$LOG"

REGIONS="$OUT_DIR/clinvar.regions.tsv"
SUBSET_VCF="$OUT_DIR/clinvar.subset.vcf.gz"
LOOKUP_TSV="$OUT_DIR/clinvar.lookup.tsv"

echo "[2] Build regions from TSVs" | tee -a "$LOG"
python3 workflow/14A_build_regions.py "$IN_DIR" "$REGIONS" "$CHR_PREFIX" --pattern "$PATTERN" | tee -a "$LOG"

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

echo "[5] Join ClinVar into TSVs" | tee -a "$LOG"
python3 workflow/14B_join_clinvar.py "$IN_DIR" "$LOOKUP_TSV" "$OUT_DIR" --pattern "$PATTERN" | tee -a "$LOG"

echo "Step 14 complete. Outputs: $OUT_DIR/*_clinvar.tsv" | tee -a "$LOG"
