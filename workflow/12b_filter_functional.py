#!/usr/bin/env python3
import argparse
import csv
import sys
from typing import Set

DEFAULT_IMPACTS: Set[str] = {"HIGH", "MODERATE"}

def main():
    ap = argparse.ArgumentParser(
        description="Step 12b: Filter annotated TSV to functional variants (by snpEff IMPACT, optionally EFFECT)."
    )
    ap.add_argument("--in", dest="in_tsv", required=True, help="Input annotated TSV (Step 12 output)")
    ap.add_argument("--out", dest="out_tsv", required=True, help="Output filtered TSV")
    ap.add_argument(
        "--impact",
        nargs="+",
        default=sorted(DEFAULT_IMPACTS),
        help="Allowed IMPACT values (default: HIGH MODERATE)",
    )
    ap.add_argument(
        "--exclude-effect",
        nargs="*",
        default=[],
        help="EFFECT values to exclude explicitly (optional)",
    )
    args = ap.parse_args()

    allowed_impacts = set(args.impact)
    excluded_effects = set(args.exclude_effect)

    with open(args.in_tsv, encoding="utf-8") as fin, open(args.out_tsv, "w", newline="", encoding="utf-8") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            sys.exit("ERROR: input TSV has no header")

        required = {"CHROM", "POS", "REF", "ALT", "ANN[0].IMPACT"}
        missing = sorted(list(required - set(reader.fieldnames)))
        if missing:
            sys.exit(f"ERROR: missing required columns: {missing}")

        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        total = kept = drop_impact = drop_effect = 0

        for row in reader:
            total += 1
            impact = (row.get("ANN[0].IMPACT") or "").strip()
            effect = (row.get("ANN[0].EFFECT") or "").strip()

            if impact not in allowed_impacts:
                drop_impact += 1
                continue

            if effect and effect in excluded_effects:
                drop_effect += 1
                continue

            writer.writerow(row)
            kept += 1

    print(f"Input variants : {total}")
    print(f"Kept variants  : {kept}")
    print(f"Dropped impact : {drop_impact}")
    print(f"Dropped effect : {drop_effect}")
    print(f"Output written : {args.out_tsv}")

if __name__ == "__main__":
    main()
