#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Step 12: Create final prioritized candidate table
# Combines snpEff annotations with evidence tables
# -----------------------------

JAVA21="/usr/lib/jvm/java-21-openjdk-amd64/bin/java"

ANN_DIR="results/annotation"
EVID_DIR="results/candidates_evidence_simple"
OUT_DIR="results/12_final_candidates"
LOG_DIR="logs"

SNPSIFT_JAR="tools/snpeff/snpEff/SnpSift.jar"

mkdir -p "$OUT_DIR" "$LOG_DIR"
LOG="$LOG_DIR/12_final_table.log"
: > "$LOG"

echo "=== Step 12: Create final candidate table ===" | tee -a "$LOG"
echo "Annotation dir: $ANN_DIR" | tee -a "$LOG"
echo "Evidence dir: $EVID_DIR" | tee -a "$LOG"
echo "Output dir: $OUT_DIR" | tee -a "$LOG"

# --- sanity checks ---
if [[ ! -f "$SNPSIFT_JAR" ]]; then
  echo "ERROR: SnpSift.jar not found at: $SNPSIFT_JAR" | tee -a "$LOG"
  exit 1
fi

# Function: extract annotations from VCF to TSV
extract_annotations () {
  local label="$1"
  local vcf_ann="$ANN_DIR/${label}.snpeff.vcf.gz"
  local tsv_out="$OUT_DIR/${label}_annotated.tsv"

  if [[ ! -f "$vcf_ann" ]]; then
    echo "SKIP [$label]: missing $vcf_ann" | tee -a "$LOG"
    return 0
  fi

  echo "--- [$label] ---" | tee -a "$LOG"
  echo "[1] Extracting fields from $vcf_ann" | tee -a "$LOG"

  # Extract key fields using SnpSift
  "$JAVA21" -jar "$SNPSIFT_JAR" extractFields \
    -s "," -e "." \
    "$vcf_ann" \
    CHROM POS ID REF ALT QUAL FILTER \
    "ANN[0].GENE" \
    "ANN[0].GENEID" \
    "ANN[0].FEATUREID" \
    "ANN[0].EFFECT" \
    "ANN[0].IMPACT" \
    "ANN[0].HGVS_C" \
    "ANN[0].HGVS_P" \
    "ANN[0].CDNA_POS" \
    "ANN[0].CDS_POS" \
    "ANN[0].AA_POS" \
    "GEN[0].GT" "GEN[0].AD" "GEN[0].DP" "GEN[0].GQ" \
    "GEN[1].GT" "GEN[1].AD" "GEN[1].DP" "GEN[1].GQ" \
    "GEN[2].GT" "GEN[2].AD" "GEN[2].DP" "GEN[2].GQ" \
    > "$tsv_out"

  echo "    Variants annotated: $(tail -n +2 "$tsv_out" | wc -l)" | tee -a "$LOG"
  echo "DONE [$label]: $tsv_out" | tee -a "$LOG"
}

# Extract annotations for each inheritance mode
extract_annotations "de_novo"
extract_annotations "inherited_het"
extract_annotations "recessive"

echo "" | tee -a "$LOG"
echo "[2] Merging with evidence tables (optional)" | tee -a "$LOG"

# Optional: merge with original evidence tables to add inheritance info

echo "" | tee -a "$LOG"
echo "=== Step 12 complete ===" | tee -a "$LOG"
echo "Annotated tables in: $OUT_DIR" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo "Next steps:" | tee -a "$LOG"
echo "  1. Review *_annotated.tsv files" | tee -a "$LOG"
echo "  2. Prioritize based on IMPACT, GENE, HGVS" | tee -a "$LOG"
echo "  3. Look up genes in OMIM/ClinVar" | tee -a "$LOG"


echo "[3] Merging and creating Excel reports" | tee -a "$LOG"
python3 workflow/12B_merge_annotations.py