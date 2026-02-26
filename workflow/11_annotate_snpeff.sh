#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Step 11: Functional annotation with snpEff (one-pass, then subset)
# -----------------------------

JAVA21="/usr/lib/jvm/java-21-openjdk-amd64/bin/java"

VCF_IN="results/vcf/trio.filtered.vcf.gz"
EVID_DIR="results/10_candidates_evidence"
OUT_DIR="results/11_annotation"
LOG_DIR="logs"

SNPEFF_JAR="tools/snpeff/snpEff/snpEff.jar"
SNPEFF_DB="GRCh38.mane.1.2.ensembl"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"
LOG="$LOG_DIR/11_annotation_snpeff.log"
: > "$LOG"

echo "Step 11: snpEff annotation (one-pass)" | tee -a "$LOG"
echo "Java: $JAVA21" | tee -a "$LOG"
echo "Input VCF: $VCF_IN" | tee -a "$LOG"
echo "Evidence dir: $EVID_DIR" | tee -a "$LOG"
echo "snpEff jar: $SNPEFF_JAR" | tee -a "$LOG"
echo "snpEff DB:  $SNPEFF_DB" | tee -a "$LOG"
echo "Output dir: $OUT_DIR" | tee -a "$LOG"

# --- sanity checks ---
[[ -x "${JAVA21}" ]] || { echo "ERROR: Java not found/executable at: ${JAVA21}" | tee -a "${LOG}"; exit 1; }
[[ -f "${SNPEFF_JAR}" ]] || { echo "ERROR: snpEff jar not found at: ${SNPEFF_JAR}" | tee -a "${LOG}"; exit 1; }
[[ -f "${VCF_IN}" ]] || { echo "ERROR: Input VCF not found: ${VCF_IN}" | tee -a "${LOG}"; exit 1; }

if [[ ! -d "${EVID_DIR}" ]]; then
  echo "WARN: Evidence dir not found: ${EVID_DIR} (Step 10 may not have been run)" | tee -a "${LOG}"
fi

# Ensure VCF index exists
if [[ ! -f "${VCF_IN}.tbi" ]]; then
  echo "[0] Indexing input VCF with tabix" | tee -a "${LOG}"
  tabix -p vcf "${VCF_IN}"
fi

# Output: full annotated VCF (one-pass)
VCF_ANN="$OUT_DIR/trio.filtered.snpeff.vcf.gz"

has_ann_header () {
  bcftools view -h "$1" | grep -q '##INFO=<ID=ANN'
}

echo "[1] Check if input VCF already has ANN header" | tee -a "${LOG}"
if has_ann_header "${VCF_IN}"; then
  echo "ANN already present in input VCF -> will reuse as annotated source" | tee -a "${LOG}"
  VCF_SOURCE="${VCF_IN}"
else
  echo "[2] Annotate FULL VCF with snpEff -> ${VCF_ANN}" | tee -a "${LOG}"

  # Stream VCF -> snpEff -> bgzip
  bcftools view -Ov "${VCF_IN}" \
    | "${JAVA21}" -Xmx6g -jar "${SNPEFF_JAR}" "${SNPEFF_DB}" \
    | bgzip -c > "${VCF_ANN}"

  tabix -p vcf "${VCF_ANN}"

  echo "annotated variants: $(bcftools view -H "${VCF_ANN}" | wc -l)" | tee -a "${LOG}"
  VCF_SOURCE="${VCF_ANN}"
fi

# Build regions file (CHROM POS POS) from TSV pass list (header-aware)
make_regions () {
  local tsv="$1"
  local out="$2"

  [[ -f "${tsv}" ]] || { echo "ERROR: TSV not found: ${tsv}" | tee -a "${LOG}"; return 2; }

  awk -F'\t' '
    NR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="CHROM") c=i;
        if ($i=="POS")   p=i;
      }
      if (!c || !p) {
        print "ERROR: CHROM/POS columns not found in header" > "/dev/stderr";
        exit 2
      }
      next
    }
    { print $c "\t" $p "\t" $p }
  ' "${tsv}" | sort -k1,1 -k2,2n > "${out}"
}

subset_one () {
  local label="$1"
  local tsv_pass="${EVID_DIR}/${label}_evidence_pass.tsv"

  if [[ ! -f "${tsv_pass}" ]]; then
    echo "SKIP [${label}]: missing ${tsv_pass}" | tee -a "${LOG}"
    return 0
  fi

  echo "--- [${label}] ---" | tee -a "${LOG}"

  local regions="${OUT_DIR}/${label}.regions.txt"
  local vcf_sub="${OUT_DIR}/${label}.snpeff.subset.vcf.gz"

  echo "[A] Build regions from: ${tsv_pass}" | tee -a "${LOG}"
  make_regions "${tsv_pass}" "${regions}"
  echo "regions lines: $(wc -l < "${regions}")" | tee -a "${LOG}"

  echo "[B] Subset annotated VCF -> ${vcf_sub}" | tee -a "${LOG}"
  bcftools view -R "${regions}" -Oz -o "${vcf_sub}" "${VCF_SOURCE}"
  tabix -p vcf "${vcf_sub}"
  echo "subset variants: $(bcftools view -H "${vcf_sub}" | wc -l)" | tee -a "${LOG}"
}

# Models
subset_one "de_novo"
subset_one "ar_homo"
subset_one "ar_het"
subset_one "x_linked"

echo "=== Step 11 output stats ===" | tee -a "${LOG}"
for f in "${OUT_DIR}"/*.subset.vcf.gz; do
  [[ -f "${f}" ]] || continue
  echo "$(basename "${f}") variants=$(bcftools view -H "${f}" | wc -l)" | tee -a "${LOG}"
done

echo "Step 11 done. Annotated source: ${VCF_SOURCE}" | tee -a "${LOG}"
echo "Subsets in: ${OUT_DIR}" | tee -a "${LOG}"

subset_one () {
  local label="$1"
  local tsv_pass="${EVID_DIR}/${label}_evidence_pass.tsv"

  if [[ ! -f "${tsv_pass}" ]]; then
    echo "SKIP [${label}]: missing ${tsv_pass}" | tee -a "${LOG}"
    return 0
  fi

  echo "--- [${label}] ---" | tee -a "${LOG}"

  local regions="${OUT_DIR}/${label}.regions.txt"
  local vcf_sub="${OUT_DIR}/${label}.snpeff.subset.vcf.gz"

  echo "[A] Build regions from: ${tsv_pass}" | tee -a "${LOG}"
  make_regions "${tsv_pass}" "${regions}"
  echo "    regions lines: $(wc -l < "${regions}")" | tee -a "${LOG}"

  echo "[B] Subset annotated VCF -> ${vcf_sub}" | tee -a "${LOG}"
  bcftools view -R "${regions}" -Oz -o "${vcf_sub}" "${VCF_SOURCE}"
  tabix -p vcf "${vcf_sub}"
  echo "    subset variants: $(bcftools view -H "${vcf_sub}" | wc -l)" | tee -a "${LOG}"
}

# Models
subset_one "de_novo"
subset_one "ar_homo"
subset_one "ar_het"
subset_one "x_linked"

echo "Step 11 output stats" | tee -a "${LOG}"
for f in "${OUT_DIR}"/*.subset.vcf.gz; do
  [[ -f "${f}" ]] || continue
  echo "$(basename "${f}") variants=$(bcftools view -H "${f}" | wc -l)" | tee -a "${LOG}"
done

echo "Step [11] conpleted. Annotated source: ${VCF_SOURCE}" | tee -a "${LOG}"
echo "Subsets in: ${OUT_DIR}" | tee -a "${LOG}"
