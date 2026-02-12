#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------------------------------
# Step 12: Build per-model annotated TSV tables 
# - Uses one annotated source VCF produced by Step 11
# - Subsets by regions derived from Step 10 evidence-pass TSVs
# - Exports unified TSV schema with snpEff ANN[0] + trio genotypes
# Sample order is enforced earlier: sample[0]=father, sample[1]=mother, sample[2]=proband
# ------------------------------------------------------------

ANN_VCF="results/11_annotation/trio.filtered.snpeff.vcf.gz"
EVID_DIR="results/10_candidates_evidence"
OUT_DIR="results/12_model_tables"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"
LOG="$LOG_DIR/12_model_tables.log"
: > "$LOG"

echo "=== Step 12: Build per-model annotated tables ===" | tee -a "$LOG"
echo "Annotated VCF: $ANN_VCF" | tee -a "$LOG"
echo "Evidence dir:  $EVID_DIR" | tee -a "$LOG"
echo "Output dir:    $OUT_DIR" | tee -a "$LOG"

if [[ ! -f "$ANN_VCF" ]]; then
  echo "ERROR: Annotated VCF not found: $ANN_VCF" | tee -a "$LOG"
  exit 1
fi
if [[ ! -f "${ANN_VCF}.tbi" ]]; then
  echo "[0] Index annotated VCF" | tee -a "$LOG"
  tabix -p vcf "$ANN_VCF"
fi

# Check ANN exists
if ! bcftools view -h "$ANN_VCF" | grep -q '##INFO=<ID=ANN'; then
  echo "ERROR: ANN field not found in VCF header. Step 11 annotation may have failed." | tee -a "$LOG"
  exit 1
fi

MODELS=(de_novo ar_homo ar_het x_linked)

# Header + bcftools query format
QFMT='%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER'\
'\t%ANN[0].GENE\t%ANN[0].GENEID\t%ANN[0].FEATUREID\t%ANN[0].EFFECT\t%ANN[0].IMPACT\t%ANN[0].HGVS_C\t%ANN[0].HGVS_P\t%ANN[0].CDNA_POS\t%ANN[0].CDS_POS\t%ANN[0].AA_POS'\
'[\t%GT\t%AD\t%DP\t%GQ]\n'

HEADER="inheritance\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tANN[0].GENE\tANN[0].GENEID\tANN[0].FEATUREID\tANN[0].EFFECT\tANN[0].IMPACT\tANN[0].HGVS_C\tANN[0].HGVS_P\tANN[0].CDNA_POS\tANN[0].CDS_POS\tANN[0].AA_POS\tGEN[0].GT\tGEN[0].AD\tGEN[0].DP\tGEN[0].GQ\tGEN[1].GT\tGEN[1].AD\tGEN[1].DP\tGEN[1].GQ\tGEN[2].GT\tGEN[2].AD\tGEN[2].DP\tGEN[2].GQ"

make_regions () {
  local tsv="$1"
  local out="$2"

  awk -F'\t' '
    NR==1 {
      for (i=1; i<=NF; i++) {
        if ($i=="CHROM") c=i;
        if ($i=="POS")   p=i;
        if ($i=="REF")   r=i;
        if ($i=="ALT")   a=i;
      }
      if (!c || !p) { print "ERROR: CHROM/POS not found in header" > "/dev/stderr"; exit 2 }
      next
    }
    { print $c "\t" $p "\t" $p }
  ' "$tsv" | sort -k1,1 -k2,2n > "$out"
}

for m in "${MODELS[@]}"; do
  pass_tsv="$EVID_DIR/${m}_evidence_pass.tsv"
  regions="$OUT_DIR/${m}.regions.tsv"
  sub_vcf="$OUT_DIR/${m}.subset.snpeff.vcf.gz"
  out_tsv="$OUT_DIR/${m}_annotated.tsv"

  if [[ ! -f "$pass_tsv" ]]; then
    echo "SKIP [$m]: missing $pass_tsv" | tee -a "$LOG"
    continue
  fi

  echo "--- [$m] ---" | tee -a "$LOG"
  echo "Pass TSV: $pass_tsv" | tee -a "$LOG"

  echo "[1] Build regions" | tee -a "$LOG"
  make_regions "$pass_tsv" "$regions"
  echo "    regions: $(wc -l < "$regions")" | tee -a "$LOG"

  echo "[2] Subset annotated VCF" | tee -a "$LOG"
  bcftools view -R "$regions" -Oz -o "$sub_vcf" "$ANN_VCF"
  tabix -p vcf "$sub_vcf"
  echo "    subset variants: $(bcftools view -H "$sub_vcf" | wc -l)" | tee -a "$LOG"

  echo "[3] Export TSV" | tee -a "$LOG"
  {
    echo -e "$HEADER"
    bcftools query -f "$QFMT" "$sub_vcf" | awk -v m="$m" 'BEGIN{FS=OFS="\t"} {print m, $0}'
  } > "$out_tsv"

  echo "    wrote: $out_tsv (rows=$(($(wc -l < "$out_tsv")-1)))" | tee -a "$LOG"
done

echo "=== Step 12 complete ===" | tee -a "$LOG"
echo "Output tables: $OUT_DIR/*_annotated.tsv" | tee -a "$LOG"
