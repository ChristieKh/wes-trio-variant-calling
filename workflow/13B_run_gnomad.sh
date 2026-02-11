#!/usr/bin/env bash
set -euo pipefail

IN_DIR="results/final_candidates"
OUT_DIR="results/13_external"
mkdir -p "$OUT_DIR"


# python3 workflow/13B_gnomad_api.py \
#   --in  "$IN_DIR/shortlist_de_novo.tsv" \
#   --out "$OUT_DIR/shortlist_de_novo.gnomad.tsv" \
#   --dataset gnomad_r4 \
#   --sleep 1.2 \
#   --retry-errors

# python3 workflow/13B_gnomad_api.py \
#   --in  "$IN_DIR/shortlist_high_impact.tsv" \
#   --out "$OUT_DIR/shortlist_high_impact.gnomad.tsv" \
#   --dataset gnomad_r4 \
#   --sleep 1.2 \
#   --retry-errors

# python3 workflow/13B_gnomad_api.py \
#   --in  "$IN_DIR/shortlist_ptv.tsv" \
#   --out "$OUT_DIR/shortlist_ptv.gnomad.tsv" \
#   --dataset gnomad_r4\
#   --sleep 1.2 \
#   --retry-errors

  python3 workflow/13B_gnomad_api.py \
  --in  "$IN_DIR/shortlist_recessive.tsv" \
  --out "$OUT_DIR/shortlist_recessive_gnomad.tsv" \
  --dataset gnomad_r4 \
  --cache results/13_external/gnomad_cache_v2_faf95.jsonl \
  --sleep 1.2 \
  --retry-errors
