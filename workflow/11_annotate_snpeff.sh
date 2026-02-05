#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Step 11: Functional annotation with snpEff (no Docker)
# Uses Java 21 explicitly (your system path)
# -----------------------------

JAVA21="/usr/lib/jvm/java-21-openjdk-amd64/bin/java"

VCF_IN="results/vcf/trio.filtered.vcf.gz"
EVID_DIR="results/candidates_evidence_simple"
OUT_DIR="results/annotation"
LOG_DIR="logs"

SNPEFF_JAR="tools/snpeff/snpEff/snpEff.jar"
SNPEFF_DB="GRCh38.mane.1.2.ensembl"   # GRCh38 + MANE transcripts (cleaner transcript set)

mkdir -p "$OUT_DIR" "$LOG_DIR"
LOG="$LOG_DIR/11_annotation_snpeff.log"
: > "$LOG"

echo "=== Step 11: snpEff annotation (GRCh38) ===" | tee -a "$LOG"
echo "Java: $JAVA21" | tee -a "$LOG"
echo "Input VCF: $VCF_IN" | tee -a "$LOG"
echo "Evidence dir: $EVID_DIR" | tee -a "$LOG"
echo "snpEff jar: $SNPEFF_JAR" | tee -a "$LOG"
echo "snpEff DB:  $SNPEFF_DB" | tee -a "$LOG"
echo "Output dir: $OUT_DIR" | tee -a "$LOG"

# --- sanity checks ---
if [[ ! -x "$JAVA21" ]]; then
  echo "ERROR: Java 21 not found/executable at: $JAVA21" | tee -a "$LOG"
  exit 1
fi

if [[ ! -f "$SNPEFF_JAR" ]]; then
  echo "ERROR: snpEff jar not found at: $SNPEFF_JAR" | tee -a "$LOG"
  exit 1
fi

if [[ ! -f "$VCF_IN" ]]; then
  echo "ERROR: Input VCF not found: $VCF_IN" | tee -a "$LOG"
  exit 1
fi

# Ensure VCF index exists
if [[ ! -f "${VCF_IN}.tbi" ]]; then
  echo "[0] Indexing input VCF with tabix" | tee -a "$LOG"
  tabix -p vcf "$VCF_IN"
fi

# Build regions file (CHROM POS POS) from TSV (header-aware)
make_regions () {
  local tsv="$1"
  local out="$2"

  awk -F'\t' '
    NR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="CHROM") c=i;
        if ($i=="POS")   p=i;
      }
      next
    }
    { print $c "\t" $p "\t" $p }
  ' "$tsv" | sort -k1,1 -k2,2n > "$out"
}

annotate_one () {
  local label="$1"
  local tsv_pass="$EVID_DIR/${label}_evidence_pass.tsv"

  if [[ ! -f "$tsv_pass" ]]; then
    echo "SKIP [$label]: missing $tsv_pass" | tee -a "$LOG"
    return 0
  fi

  echo "--- [$label] ---" | tee -a "$LOG"

  local regions="$OUT_DIR/${label}.regions.tsv"
  local vcf_sub="$OUT_DIR/${label}.subset.vcf.gz"
  local vcf_ann="$OUT_DIR/${label}.snpeff.vcf.gz"

  echo "[1] Build regions from: $tsv_pass" | tee -a "$LOG"
  make_regions "$tsv_pass" "$regions"
  echo "    regions lines: $(wc -l < "$regions")" | tee -a "$LOG"

  echo "[2] Subset VCF to candidate regions" | tee -a "$LOG"
  bcftools view -R "$regions" -Oz -o "$vcf_sub" "$VCF_IN"
  tabix -p vcf "$vcf_sub"
  echo "    subset variants: $(bcftools view -H "$vcf_sub" | wc -l)" | tee -a "$LOG"

  echo "[3] Annotate with snpEff -> $vcf_ann" | tee -a "$LOG"
  # snpEff reads VCF from stdin and writes annotated VCF to stdout
  bcftools view -Ov "$vcf_sub" \
    | "$JAVA21" -Xmx4g -jar "$SNPEFF_JAR" "$SNPEFF_DB" \
    | bgzip -c > "$vcf_ann"
  tabix -p vcf "$vcf_ann"

  echo "DONE [$label]: $vcf_ann" | tee -a "$LOG"
}

annotate_one "de_novo"
annotate_one "inherited_het"
annotate_one "recessive"

echo "Step 11 done. Annotated VCFs in: $OUT_DIR" | tee -a "$LOG"
