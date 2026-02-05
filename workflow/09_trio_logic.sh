#!/usr/bin/env bash
set -euo pipefail

VCF_IN="results/vcf/trio.filtered.vcf.gz"
OUT_DIR="results/candidates"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

LOG="$LOG_DIR/09_trio_logic.log"
: > "$LOG"

echo "=== Step 09: Trio logic ===" | tee -a "$LOG"
echo "Input: $VCF_IN" | tee -a "$LOG"

echo "[0] Samples in VCF (current order):" | tee -a "$LOG"
bcftools query -l "$VCF_IN" | nl -ba | tee -a "$LOG"

# Sample names 
FATHER="father"
MOTHER="mother"
PROBAND="proband"

echo "[1] Enforcing trio order: ${FATHER},${MOTHER},${PROBAND}" | tee -a "$LOG"

# After we reorder samples with -s:
#   sample[0]=father, sample[1]=mother, sample[2]=proband
# QC must be applied per-sample explicitly using indices.
COMMON_QC='GT[0]!="mis" && GT[1]!="mis" && GT[2]!="mis" && FORMAT/DP[0]>=10 && FORMAT/DP[1]>=10 && FORMAT/DP[2]>=10 && FORMAT/GQ[0]>=20 && FORMAT/GQ[1]>=20 && FORMAT/GQ[2]>=20'

HEADER="CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tGT_father\tDP_father\tGQ_father\tAD_father\tGT_mother\tDP_mother\tGQ_mother\tAD_mother\tGT_proband\tDP_proband\tGQ_proband\tAD_proband"

# IMPORTANT FIX:
# Put a leading \t inside the repeated block [ ... ]
# so that each sample group starts with a tab and never "sticks" to the previous field.
QUERY_FMT='%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO[\t%GT\t%DP\t%GQ\t%AD]\n'

# Stream with enforced trio order
TRIO_STREAM_CMD=(bcftools view -Ou -s "${FATHER},${MOTHER},${PROBAND}" "$VCF_IN")

# 1) de novo candidates:
# parents 0/0, proband 0/1 or 1/1
echo -e "$HEADER" > "$OUT_DIR/de_novo.tsv"
"${TRIO_STREAM_CMD[@]}" \
| bcftools filter -Ou -i "$COMMON_QC" \
| bcftools filter -Ou -i '((GT[0]="0/0" || GT[0]="0|0") && (GT[1]="0/0" || GT[1]="0|0") && (GT[2]="0/1" || GT[2]="1/0" || GT[2]="0|1" || GT[2]="1|0" || GT[2]="1/1" || GT[2]="1|1"))' \
| bcftools query -f "$QUERY_FMT" \
>> "$OUT_DIR/de_novo.tsv"

echo "[de_novo] lines (including header):" | tee -a "$LOG"
wc -l "$OUT_DIR/de_novo.tsv" | tee -a "$LOG"

# 2) autosomal recessive candidates:
# parents 0/1, proband 1/1
echo -e "$HEADER" > "$OUT_DIR/recessive.tsv"
"${TRIO_STREAM_CMD[@]}" \
| bcftools filter -Ou -i "$COMMON_QC" \
| bcftools filter -Ou -i '((GT[0]="0/1" || GT[0]="1/0" || GT[0]="0|1" || GT[0]="1|0") && (GT[1]="0/1" || GT[1]="1/0" || GT[1]="0|1" || GT[1]="1|0") && (GT[2]="1/1" || GT[2]="1|1"))' \
| bcftools query -f "$QUERY_FMT" \
>> "$OUT_DIR/recessive.tsv"

echo "[recessive] lines (including header):" | tee -a "$LOG"
wc -l "$OUT_DIR/recessive.tsv" | tee -a "$LOG"

# 3) inherited heterozygous (simple transmission):
# proband 0/1 and exactly one parent 0/1, other 0/0
echo -e "$HEADER" > "$OUT_DIR/inherited_het.tsv"
"${TRIO_STREAM_CMD[@]}" \
| bcftools filter -Ou -i "$COMMON_QC" \
| bcftools filter -Ou -i '((GT[2]="0/1" || GT[2]="1/0" || GT[2]="0|1" || GT[2]="1|0") && (((GT[0]="0/1" || GT[0]="1/0" || GT[0]="0|1" || GT[0]="1|0") && (GT[1]="0/0" || GT[1]="0|0")) || ((GT[1]="0/1" || GT[1]="1/0" || GT[1]="0|1" || GT[1]="1|0") && (GT[0]="0/0" || GT[0]="0|0"))))' \
| bcftools query -f "$QUERY_FMT" \
>> "$OUT_DIR/inherited_het.tsv"

echo "[inherited_het] lines (including header):" | tee -a "$LOG"
wc -l "$OUT_DIR/inherited_het.tsv" | tee -a "$LOG"

echo "Step [09] Trio candidate lists created in: $OUT_DIR" | tee -a "$LOG"

