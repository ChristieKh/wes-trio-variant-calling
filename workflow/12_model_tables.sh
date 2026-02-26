#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Step 12: Build per-model annotated TSV tables
# - Uses annotated source VCF produced by Step 11 (INFO/ANN present)
# - Subsets by regions derived from Step 10 evidence-pass TSVs (position-based)
# - Exports unified RAW TSV schema (includes ANN_RAW + trio genotypes)
# - Parses ANN_RAW into structured columns via Python parser
# NOTE: Subsetting is position-based (CHROM:POS). 
# ------------------------------------------------------------

ANN_VCF="results/11_annotation/trio.filtered.snpeff.vcf.gz"
EVID_DIR="results/10_candidates_evidence"
OUT_DIR="results/12_model_tables"
LOG_DIR="logs"
ANN_PARSER="workflow/12_split_snpeff_ann.py"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"
LOG="${LOG_DIR}/12_model_tables.log"
: > "${LOG}"

echo "Step 12: Build per-model annotated tables" | tee -a "${LOG}"
echo "Annotated VCF: ${ANN_VCF}" | tee -a "${LOG}"
echo "Evidence dir:  ${EVID_DIR}" | tee -a "${LOG}"
echo "Output dir:    ${OUT_DIR}" | tee -a "${LOG}"
echo "ANN parser:    ${ANN_PARSER}" | tee -a "${LOG}"

# --- sanity checks ---
[[ -f "${ANN_VCF}" ]] || { echo "ERROR: Annotated VCF not found: ${ANN_VCF}" | tee -a "${LOG}"; exit 1; }
[[ -f "${ANN_PARSER}" ]] || { echo "ERROR: ANN parser not found: ${ANN_PARSER}" | tee -a "${LOG}"; exit 1; }

[[ -d "${EVID_DIR}" ]] || { echo "WARN: Evidence dir not found: ${EVID_DIR} (Step 10 may not have been run)" | tee -a "${LOG}"; }

if [[ ! -f "${ANN_VCF}.tbi" ]]; then
  echo "[0] Index annotated VCF" | tee -a "${LOG}"
  tabix -p vcf "${ANN_VCF}"
fi

# Check ANN exists
if ! bcftools view -h "${ANN_VCF}" | grep -q '##INFO=<ID=ANN'; then
  echo "ERROR: ANN field not found in VCF header. Step 11 annotation may have failed." | tee -a "${LOG}"
  exit 1
fi

echo "Samples in annotated VCF (order matters for GEN[0..2]):" | tee -a "${LOG}"
bcftools query -l "${ANN_VCF}" | nl -ba | tee -a "${LOG}"

MODELS=(de_novo ar_homo ar_het x_linked)

# --- bcftools query format ---
# Export raw ANN string in one column ANN_RAW, then parse in Python
QFMT='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER'\
'\t%INFO/ANN'\
'[\t%GT\t%AD\t%DP\t%GQ]\n'

# Raw header (ANN_RAW)
HEADER_RAW="inheritance\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tANN_RAW\tGEN[0].GT\tGEN[0].AD\tGEN[0].DP\tGEN[0].GQ\tGEN[1].GT\tGEN[1].AD\tGEN[1].DP\tGEN[1].GQ\tGEN[2].GT\tGEN[2].AD\tGEN[2].DP\tGEN[2].GQ"

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
        print "ERROR: CHROM/POS not found in header" > "/dev/stderr";
        exit 2
      }
      next
    }
    { print $c "\t" $p "\t" $p }
  ' "${tsv}" | sort -k1,1 -k2,2n > "${out}"
}

# Check that parsed gene column is not still raw (contains '|') for ANY data row
check_ann_parsed () {
  local tsv="$1"
  awk -F'\t' '
    NR==1 {
      for(i=1;i<=NF;i++) if($i=="ANN[0].GENE") gi=i;
      if(!gi){ exit 2 }
      next
    }
    NR>1 {
      if($gi ~ /\|/){ exit 1 }
    }
    END { exit 0 }
  ' "${tsv}"
}

for m in "${MODELS[@]}"; do
  pass_tsv="${EVID_DIR}/${m}_evidence_pass.tsv"
  regions="${OUT_DIR}/${m}.regions.tsv"
  sub_vcf="${OUT_DIR}/${m}.subset.snpeff.vcf.gz"

  out_raw="${OUT_DIR}/${m}_raw.tsv"
  out_tsv="${OUT_DIR}/${m}_annotated.tsv"

  if [[ ! -f "${pass_tsv}" ]]; then
    echo "SKIP [${m}]: missing ${pass_tsv}" | tee -a "${LOG}"
    continue
  fi

  echo "--- [${m}] ---" | tee -a "${LOG}"
  echo "Pass TSV: ${pass_tsv}" | tee -a "${LOG}"

  echo "[1] Build regions" | tee -a "${LOG}"
  make_regions "${pass_tsv}" "${regions}"
  n_regions=$(wc -l < "${regions}")
  echo "    regions: ${n_regions}" | tee -a "${LOG}"

  if [[ "${n_regions}" -eq 0 ]]; then
    echo "WARN: no regions for model ${m} (empty pass TSV) -> skipping" | tee -a "${LOG}"
    continue
  fi

  echo "[2] Subset annotated VCF" | tee -a "${LOG}"
  bcftools view -R "${regions}" -Oz -o "${sub_vcf}" "${ANN_VCF}"
  tabix -p vcf "${sub_vcf}"
  n_sub=$(bcftools view -H "${sub_vcf}" | wc -l)
  echo "    subset variants: ${n_sub}" | tee -a "${LOG}"

  if [[ "${n_sub}" -eq 0 ]]; then
    echo "WARN: subset contains 0 variants for model ${m} -> skipping export/parse" | tee -a "${LOG}"
    continue
  fi

  echo "[3] Export RAW TSV (ANN_RAW)" | tee -a "${LOG}"
  {
    echo -e "${HEADER_RAW}"
    bcftools query -f "${QFMT}" "${sub_vcf}" | awk -v m="${m}" 'BEGIN{FS=OFS="\t"} {print m, $0}'
  } > "${out_raw}"
  echo "    wrote: ${out_raw} (rows=$(($(wc -l < "${out_raw}")-1)))" | tee -a "${LOG}"

  echo "[4] Parse ANN_RAW -> annotated TSV" | tee -a "${LOG}"
  python3 "${ANN_PARSER}" "${out_raw}" "${out_tsv}"

  if check_ann_parsed "${out_tsv}"; then
    echo "    annotated OK: ${out_tsv}" | tee -a "${LOG}"
  else
    rc=$?
    if [[ "${rc}" -eq 2 ]]; then
      echo "ERROR: ANN[0].GENE column not found in parsed TSV header: ${out_tsv}" | tee -a "${LOG}"
    else
      echo "ERROR: Parsed gene field still contains '|' in at least one row (ANN parsing failed). File: ${out_tsv}" | tee -a "${LOG}"
    fi
    exit 1
  fi

  echo "wrote: ${out_tsv} (rows=$(($(wc -l < "${out_tsv}")-1)))" | tee -a "${LOG}"
done

echo "Step [12] completed" | tee -a "${LOG}"
echo "Output tables: ${OUT_DIR}/*_annotated.tsv" | tee -a "${LOG}"