#!/usr/bin/env python3
"""
Usage:
python3 workflow/13a_build_master_variants.py \
  --out results/13_gnomad/master_variants.tsv \
  --inputs \
    results/12_model_tables/de_novo_functional.tsv \
    results/12_model_tables/ar_homo_functional.tsv \
    results/16_comphet/comphet_variants.tsv \
    results/12_model_tables/x_linked_functional.tsv
"""
import argparse
import csv
import os
from typing import Dict, Tuple, List

def norm_chrom(c: str) -> str:
    c = (c or "").strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c

def key(row: Dict[str, str]) -> Tuple[str, str, str, str]:
    return (
        norm_chrom(row.get("CHROM","")),
        (row.get("POS","") or "").strip(),
        (row.get("REF","") or "").strip(),
        (row.get("ALT","") or "").strip(),
    )

def read_tsv(path: str) -> Tuple[List[str], List[Dict[str,str]]]:
    with open(path, encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"ERROR: no header in {path}")
        rows = list(r)
        return r.fieldnames, rows

def main():
    ap = argparse.ArgumentParser(description="Build master unique variants list across inheritance models (by CHROM+POS+REF+ALT).")
    ap.add_argument("--out", required=True, help="Output master TSV")
    ap.add_argument("--inputs", nargs="+", required=True, help="Input TSV files (models), must contain CHROM POS REF ALT")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    required = {"CHROM","POS","REF","ALT"}

    # We'll output a minimal schema + sources
    out_fields = ["CHROM","POS","REF","ALT","sources"]

    uniq: Dict[Tuple[str,str,str,str], Dict[str,str]] = {}

    for path in args.inputs:
        _, rows = read_tsv(path)
        for row in rows:
            if not required.issubset(row.keys()):
                raise SystemExit(f"ERROR: {path} missing one of {sorted(required)}")
            k = key(row)
            if any(x == "" for x in k):
                continue
            if k not in uniq:
                # keep CHROM as chr-prefixed for regions convenience later
                chrom = row.get("CHROM","").strip()
                if not chrom.startswith("chr"):
                    chrom = "chr" + norm_chrom(chrom)
                uniq[k] = {
                    "CHROM": chrom,
                    "POS": (row.get("POS","") or "").strip(),
                    "REF": (row.get("REF","") or "").strip(),
                    "ALT": (row.get("ALT","") or "").strip(),
                    "sources": os.path.basename(path),
                }
            else:
                s = uniq[k]["sources"]
                b = os.path.basename(path)
                if b not in s.split(","):
                    uniq[k]["sources"] = s + "," + b

    with open(args.out, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=out_fields, delimiter="\t")
        w.writeheader()
        for _, row in uniq.items():
            w.writerow(row)

    print(f"Inputs: {len(args.inputs)}")
    print(f"Unique variants written: {len(uniq)}")
    print(f"Output: {args.out}")

if __name__ == "__main__":
    main()
