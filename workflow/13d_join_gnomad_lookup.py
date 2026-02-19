#!/usr/bin/env python3
"""
Join gnomAD lookup TSV (CHROM POS REF ALT AC AN AF) into an input TSV by
exact key: CHROM+POS+REF+ALT (chr prefix is normalized away on both sides).

- Input TSV must have: CHROM, POS, REF, ALT
- Lookup TSV must have: CHROM, POS, REF, ALT, AC, AN, AF
- Output appends: gnomAD_AC, gnomAD_AN, gnomAD_AF
- Prints matched/unmatched stats.

Usage:
  python3 workflow/13d_join_gnomad_lookup.py \
    --in results/12_model_tables/x_linked_functional.tsv \
    --lookup results/13_gnomad/gnomad_exomes_v4_lookup.tsv \
    --out results/13_gnomad/final_results/x_linked_functional.gnomad.tsv
"""
import argparse
import csv
from typing import Dict, Tuple, Optional


def norm_chrom(chrom: str) -> str:
    c = (chrom or "").strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c


def make_key(chrom: str, pos: str, ref: str, alt: str) -> Tuple[str, str, str, str]:
    return (norm_chrom(chrom), (pos or "").strip(), (ref or "").strip(), (alt or "").strip())


def safe_get(row: Dict[str, str], k: str) -> str:
    v = row.get(k, "")
    return (v or "").strip()


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Join gnomAD lookup (CHROM POS REF ALT AC AN AF) into TSV by CHROM+POS+REF+ALT."
    )
    ap.add_argument("--in", dest="in_tsv", required=True, help="Input TSV (must contain CHROM POS REF ALT)")
    ap.add_argument("--lookup", required=True, help="Lookup TSV (must contain CHROM POS REF ALT AC AN AF)")
    ap.add_argument("--out", dest="out_tsv", required=True, help="Output TSV with gnomAD_* columns")
    args = ap.parse_args()

    # ---- Load lookup into dict ----
    lookup: Dict[Tuple[str, str, str, str], Dict[str, str]] = {}
    dup_count = 0

    with open(args.lookup, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit("ERROR: lookup TSV has no header")

        required_lk = {"CHROM", "POS", "REF", "ALT"}
        missing = sorted(list(required_lk - set(r.fieldnames)))
        if missing:
            raise SystemExit(f"ERROR: lookup missing required columns: {missing}")

        # AC/AN/AF can be missing in some sources; we still join and fill '.'
        for row in r:
            k = make_key(safe_get(row, "CHROM"), safe_get(row, "POS"), safe_get(row, "REF"), safe_get(row, "ALT"))
            if any(x == "" for x in k):
                continue

            if k in lookup:
                dup_count += 1  # keep first; duplicates typically come from multiple lines but same key
                continue

            lookup[k] = {
                "gnomAD_AC": safe_get(row, "AC") or ".",
                "gnomAD_AN": safe_get(row, "AN") or ".",
                "gnomAD_AF": safe_get(row, "AF") or ".",
            }

    # ---- Stream input -> output with joined fields ----
    matched = 0
    unmatched = 0

    with open(args.in_tsv, encoding="utf-8") as fin, open(args.out_tsv, "w", newline="", encoding="utf-8") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: input TSV has no header")

        required_in = {"CHROM", "POS", "REF", "ALT"}
        missing = sorted(list(required_in - set(reader.fieldnames)))
        if missing:
            raise SystemExit(f"ERROR: input missing required columns: {missing}")

        out_fields = list(reader.fieldnames)
        for col in ["gnomAD_AC", "gnomAD_AN", "gnomAD_AF"]:
            if col not in out_fields:
                out_fields.append(col)

        writer = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            k = make_key(safe_get(row, "CHROM"), safe_get(row, "POS"), safe_get(row, "REF"), safe_get(row, "ALT"))
            hit = lookup.get(k)
            if hit is None:
                row["gnomAD_AC"] = "."
                row["gnomAD_AN"] = "."
                row["gnomAD_AF"] = "."
                unmatched += 1
            else:
                row.update(hit)
                matched += 1
            writer.writerow(row)

    print(f"Lookup entries : {len(lookup)}")
    if dup_count:
        print(f"Lookup dups    : {dup_count} (ignored; kept first occurrence)")
    print(f"Matched rows   : {matched}")
    print(f"Unmatched rows : {unmatched}")
    print(f"Wrote          : {args.out_tsv}")


if __name__ == "__main__":
    main()
