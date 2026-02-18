#!/usr/bin/env bash
set -euo pipefail

REGIONS="results/13_gnomad/master_variants.regions.tsv"
OUTDIR="results/13_gnomad"
TMPDIR="${OUTDIR}/tmp_gnomad_master"
LOGDIR="logs"
LOG="${LOGDIR}/13_gnomad_ucsc_master.log"
LOOKUP="${OUTDIR}/gnomad_exomes_v4_lookup.tsv"
JOBS=2

mkdir -p "$OUTDIR" "$TMPDIR" "$LOGDIR"
: > "$LOG"

echo "[$(date)] Start gnomAD UCSC MASTER lookup" | tee -a "$LOG"
echo "Regions: $REGIONS" | tee -a "$LOG"
echo "Lookup : $LOOKUP" | tee -a "$LOG"
echo "Jobs   : $JOBS" | tee -a "$LOG"

echo -e "CHROM\tPOS\tREF\tALT\tAC\tAN\tAF" > "$LOOKUP"

cut -f1 "$REGIONS" | sort -u > "${TMPDIR}/chroms.txt"

run_chr () {
  chr="$1"
  reg_chr="${TMPDIR}/${chr}.regions.tsv"
  out_chr="${TMPDIR}/${chr}.lookup.tsv"
  done_flag="${TMPDIR}/${chr}.DONE"

  if [[ -f "$done_flag" && -s "$out_chr" ]]; then
    echo "[$(date)] [${chr}] SKIP (resume)" | tee -a "$LOG"
    return 0
  fi

  grep -P "^${chr}\t" "$REGIONS" > "$reg_chr"
  nreg=$(wc -l < "$reg_chr" || echo 0)
  if [[ "$nreg" -eq 0 ]]; then
    : > "$out_chr"
    touch "$done_flag"
    echo "[$(date)] [${chr}] SKIP (0 regions)" | tee -a "$LOG"
    return 0
  fi

  url="http://hgdownload.soe.ucsc.edu/gbdb/hg38/gnomAD/v4/exomes/gnomad.exomes.v4.0.sites.${chr}.vcf.bgz"

  echo "[$(date)] [${chr}] START regions=${nreg}" | tee -a "$LOG"
  t0=$(date +%s)

  tmp_out="${out_chr}.tmp"
  rm -f "$tmp_out"

  bcftools view -R "$reg_chr" "$url" -Ou \
    | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\n' \
    > "$tmp_out"

  mv "$tmp_out" "$out_chr"
  touch "$done_flag"

  t1=$(date +%s)
  dt=$((t1 - t0))
  nout=$(wc -l < "$out_chr" || echo 0)
  echo "[$(date)] [${chr}] DONE secs=${dt} rows=${nout}" | tee -a "$LOG"
}

export -f run_chr
export REGIONS TMPDIR LOG

cat "${TMPDIR}/chroms.txt" | xargs -n 1 -P "$JOBS" -I {} bash -lc 'run_chr "$@"' _ {}

# merge
echo -e "CHROM\tPOS\tREF\tALT\tAC\tAN\tAF" > "$LOOKUP"
cat "${TMPDIR}"/*.lookup.tsv >> "$LOOKUP"

echo "[$(date)] DONE. Lookup rows: $(wc -l < "$LOOKUP")" | tee -a "$LOG"
echo "Log: $LOG" | tee -a "$LOG"
