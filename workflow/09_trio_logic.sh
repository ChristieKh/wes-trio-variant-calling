#!/usr/bin/env bash
set -euo pipefail

VCF_IN="results/vcf/trio.filtered.vcf.gz"
OUT_DIR="results/09_candidates"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

LOG="$LOG_DIR/09_trio_logic.log"
: > "$LOG"

echo "=== Step 09: Trio logic (model-based) ===" | tee -a "$LOG"
echo "Input: $VCF_IN" | tee -a "$LOG"

bcftools query -l "$VCF_IN" | nl -ba | tee -a "$LOG"

FATHER="father"
MOTHER="mother"
PROBAND="proband"

COMMON_QC='GT[0]!="mis" && GT[1]!="mis" && GT[2]!="mis" && FORMAT/DP[0]>=10 && FORMAT/DP[1]>=10 && FORMAT/DP[2]>=10 && FORMAT/GQ[0]>=20 && FORMAT/GQ[1]>=20 && FORMAT/GQ[2]>=20'

HEADER="CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tGT_father\tDP_father\tGQ_father\tAD_father\tGT_mother\tDP_mother\tGQ_mother\tAD_mother\tGT_proband\tDP_proband\tGQ_proband\tAD_proband"

QUERY_FMT='%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO[\t%GT\t%DP\t%GQ\t%AD]\n'

TRIO_STREAM_CMD=(bcftools view -Ou -s "${FATHER},${MOTHER},${PROBAND}" "$VCF_IN")

# 1 AD de novo
echo -e "$HEADER" > "$OUT_DIR/de_novo.tsv"

"${TRIO_STREAM_CMD[@]}" \
| bcftools filter -Ou -i "$COMMON_QC" \
| bcftools filter -Ou -i '((GT[0]="0/0" || GT[0]="0|0") && (GT[1]="0/0" || GT[1]="0|0") && (GT[2]!="0/0" && GT[2]!="0|0"))' \
| bcftools query -f "$QUERY_FMT" \
>> "$OUT_DIR/de_novo.tsv"

# 2 AR homozygous
echo -e "$HEADER" > "$OUT_DIR/ar_homo.tsv"

"${TRIO_STREAM_CMD[@]}" \
| bcftools filter -Ou -i "$COMMON_QC" \
| bcftools filter -Ou -i '((GT[2]="1/1" || GT[2]="1|1") &&
                          (GT[0]="0/1" || GT[0]="1/0" || GT[0]="0|1" || GT[0]="1|0") &&
                          (GT[1]="0/1" || GT[1]="1/0" || GT[1]="0|1" || GT[1]="1|0"))' \
| bcftools query -f "$QUERY_FMT" \
>> "$OUT_DIR/ar_homo.tsv"

# 3️ AR heterozygous (raw data for comp-het)
echo -e "$HEADER" > "$OUT_DIR/ar_het.tsv"

"${TRIO_STREAM_CMD[@]}" \
| bcftools filter -Ou -i "$COMMON_QC" \
| bcftools filter -Ou -i '((GT[2]="0/1" || GT[2]="1/0" || GT[2]="0|1" || GT[2]="1|0") &&
                          ((GT[0]="0/1" || GT[0]="1/0" || GT[0]="0|1" || GT[0]="1|0") ||
                           (GT[1]="0/1" || GT[1]="1/0" || GT[1]="0|1" || GT[1]="1|0")))' \
| bcftools query -f "$QUERY_FMT" \
>> "$OUT_DIR/ar_het.tsv"

# 4️ X-linked (proban is a boy)
# X-linked recessive (male proband): mother carries ALT, proband has ALT
X_QC="$COMMON_QC && (CHROM=\"X\" || CHROM=\"chrX\")"

X_MOTHER_CARRIER='(GT[1]="0/1" || GT[1]="1/0" || GT[1]="0|1" || GT[1]="1|0")'

X_PROBAND_HAS_ALT='(GT[2]!="0/0" && GT[2]!="0|0")'

echo -e "$HEADER" > "$OUT_DIR/x_linked.tsv"
"${TRIO_STREAM_CMD[@]}" \
| bcftools filter -Ou -i "$X_QC" \
| bcftools filter -Ou -i "($X_MOTHER_CARRIER) && ($X_PROBAND_HAS_ALT)" \
| bcftools query -f "$QUERY_FMT" \
>> "$OUT_DIR/x_linked.tsv"


# counter for statistics
count_lines () {
  local f="$1"
  if [[ -f "$f" ]]; then
    local n
    n=$(wc -l < "$f")
    local data=$(( n > 0 ? n-1 : 0 ))
    echo "$(basename "$f")  lines_total=$n  data_rows=$data" | tee -a "$LOG"
  else
    echo "$(basename "$f")  MISSING" | tee -a "$LOG"
  fi
}

echo "=== Step 09 output stats ===" | tee -a "$LOG"
count_lines "$OUT_DIR/de_novo.tsv"
count_lines "$OUT_DIR/ar_homo.tsv"
count_lines "$OUT_DIR/ar_het.tsv"
count_lines "$OUT_DIR/x_linked.tsv"


echo "Step 09 complete. Model-based candidate lists created." | tee -a "$LOG"
