#!/usr/bin/env bash
set -e

IN_DIR="results/readgroups"
OUT_DIR="results/markduplicates"
LOG_DIR="logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"

for SAMPLE in father mother proband; do
  echo "=== MarkDuplicates for $SAMPLE ==="

  picard MarkDuplicates \
    I="$IN_DIR/$SAMPLE.sorted.rg.bam" \
    O="$OUT_DIR/$SAMPLE.sorted.rg.dedup.bam" \
    M="$OUT_DIR/$SAMPLE.dedup.metrics.txt" \
    CREATE_INDEX=true \
    > "$LOG_DIR/$SAMPLE.markdup.log" 2>&1

  samtools flagstat "$OUT_DIR/$SAMPLE.sorted.rg.dedup.bam" > "$OUT_DIR/$SAMPLE.dedup.flagstat.txt"
done

echo "Step [04] Is Done."
