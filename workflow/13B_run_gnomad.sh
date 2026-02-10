#!/usr/bin/env bash
set -euo pipefail

IN_DIR="results/final_candidates"
OUT_DIR="results/13_external"
mkdir -p "$OUT_DIR"


python3 workflow/13B_gnomad_api.py \
  --in  "$IN_DIR/shortlist_de_novo.tsv" \
  --out "$OUT_DIR/shortlist_de_novo.gnomad.tsv" \
  --dataset gnomad_r4 \
  --retry-errors

python3 workflow/13B_gnomad_api.py \
  --in  "$IN_DIR/shortlist_high_impact.tsv" \
  --out "$OUT_DIR/shortlist_high_impact.gnomad.tsv" \
  --dataset gnomad_r4 \
  --retry-errors

python3 workflow/13B_gnomad_api.py \
  --in  "$IN_DIR/shortlist_ptv.tsv" \
  --out "$OUT_DIR/shortlist_ptv.gnomad.tsv" \
  --dataset gnomad_r4\
  --retry-errors
