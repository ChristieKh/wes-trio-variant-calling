#!/usr/bin/env bash
set -euo pipefail

CAND_TSV="results/final_candidates/final_candidates.tsv"

RES_DIR="results/resources"
OUT_DIR="results/13_external"
LOG_DIR="logs"

mkdir -p "$RES_DIR" "$OUT_DIR" "$LOG_DIR"
LOG="$LOG_DIR/13A_clinvar.log"
: > "$LOG"

CLINVAR_VCF_GZ="$RES_DIR/clinvar_GRCh38.vcf.gz"
CLINVAR_TBI="$RES_DIR/clinvar_GRCh38.vcf.gz.tbi"

REGIONS="$OUT_DIR/candidates.clinvar.regions.tsv"
CLINVAR_SUB_VCF="$OUT_DIR/clinvar.subset.vcf.gz"
CLINVAR_LOOKUP="$OUT_DIR/clinvar.lookup.tsv"

OUT_TSV="$OUT_DIR/13_candidates_with_clinvar.tsv"

echo "=== Step 13A: ClinVar (GRCh38) annotation ===" | tee -a "$LOG"
echo "Candidates: $CAND_TSV" | tee -a "$LOG"
echo "Output:     $OUT_TSV" | tee -a "$LOG"

if [[ ! -f "$CAND_TSV" ]]; then
  echo "ERROR: Candidate TSV not found: $CAND_TSV" | tee -a "$LOG"
  exit 1
fi

if [[ ! -f "$CLINVAR_VCF_GZ" ]]; then
  echo "[1] Download ClinVar GRCh38 VCF" | tee -a "$LOG"
  wget -O "$CLINVAR_VCF_GZ" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
fi

if [[ ! -f "$CLINVAR_TBI" ]]; then
  echo "[2] Index ClinVar VCF with tabix" | tee -a "$LOG"
  tabix -p vcf "$CLINVAR_VCF_GZ"
fi

echo "[3] Build ClinVar-style regions list (no chr prefix; M -> MT)" | tee -a "$LOG"
awk -F'\t' '
  function norm(c){
    sub(/^chr/, "", c)
    if (c=="M") c="MT"
    return c
  }
  NR==1{
    for(i=1;i<=NF;i++){ if($i=="CHROM") c=i; if($i=="POS") p=i }
    next
  }
  {
    chr = norm($c)
    pos = $p
    # skip contigs like GL000..., KI..., etc. ClinVar typically won’t have them
    if (chr ~ /^(GL|KI)/) next
    print chr "\t" pos "\t" pos
  }
' "$CAND_TSV" | sort -k1,1 -k2,2n > "$REGIONS"

echo "    regions: $(wc -l < "$REGIONS")" | tee -a "$LOG"

echo "[4] Subset ClinVar to candidate positions" | tee -a "$LOG"
bcftools view -R "$REGIONS" -Oz -o "$CLINVAR_SUB_VCF" "$CLINVAR_VCF_GZ"
tabix -p vcf "$CLINVAR_SUB_VCF"
echo "    subset variants: $(bcftools view -H "$CLINVAR_SUB_VCF" | wc -l)" | tee -a "$LOG"

echo "[5] Build ClinVar lookup TSV" | tee -a "$LOG"
echo -e "CHROM\tPOS\tREF\tALT\tCLNSIG\tCLNREVSTAT\tCLNDN" > "$CLINVAR_LOOKUP"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/CLNDN\n' \
  "$CLINVAR_SUB_VCF" >> "$CLINVAR_LOOKUP"
echo "    lookup lines (incl header): $(wc -l < "$CLINVAR_LOOKUP")" | tee -a "$LOG"

echo "[6] Join ClinVar fields into candidate table (chrom normalized)" | tee -a "$LOG"
python3 workflow/13A_join_clinvar.py \
  --candidates "$CAND_TSV" \
  --clinvar "$CLINVAR_LOOKUP" \
  --out "$OUT_TSV" \
  | tee -a "$LOG"

echo "Step 13A done: $OUT_TSV" | tee -a "$LOG"
