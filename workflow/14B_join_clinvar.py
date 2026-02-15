#!/usr/bin/env python3
import argparse
import csv
import glob
import os
from typing import Dict, Tuple

def norm_chrom(c: str) -> str:
    c = (c or "").strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c

ADD_COLS = [
    "ClinVar_CLNSIG",
    "ClinVar_CLNREVSTAT",
    "ClinVar_CLNDN",
    "ClinVar_CLNDISDB",
    "ClinVar_CLNVC",
    "ClinVar_ID",
]

def load_lookup(lookup_path: str) -> Dict[Tuple[str, str, str, str], Dict[str, str]]:
    lookup: Dict[Tuple[str, str, str, str], Dict[str, str]] = {}
    with open(lookup_path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            chrom = norm_chrom(row.get("CHROM", ""))
            key = (chrom, row.get("POS", ""), row.get("REF", ""), row.get("ALT", ""))
            lookup[key] = row
    return lookup

def join_one(path_in: str, path_out: str, lookup: Dict[Tuple[str, str, str, str], Dict[str, str]]):
    with open(path_in, newline="", encoding="utf-8") as fin, open(path_out, "w", newline="", encoding="utf-8") as fout:
        r = csv.DictReader(fin, delimiter="\t")
        if r.fieldnames is None:
            raise RuntimeError(f"Missing header: {path_in}")

        required = {"CHROM", "POS", "REF", "ALT"}
        missing = sorted(list(required - set(r.fieldnames)))
        if missing:
            raise RuntimeError(f"Missing required columns in {path_in}: {missing}")

        out_fields = list(r.fieldnames) + ADD_COLS
        w = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        total = 0
        hits = 0

        for row in r:
            total += 1
            chrom = norm_chrom(row.get("CHROM", ""))
            key = (chrom, row.get("POS", ""), row.get("REF", ""), row.get("ALT", ""))
            cv = lookup.get(key)

            if cv:
                hits += 1
                row["ClinVar_CLNSIG"] = cv.get("CLNSIG", ".") or "."
                row["ClinVar_CLNREVSTAT"] = cv.get("CLNREVSTAT", ".") or "."
                row["ClinVar_CLNDN"] = cv.get("CLNDN", ".") or "."
                row["ClinVar_CLNDISDB"] = cv.get("CLNDISDB", ".") or "."
                row["ClinVar_CLNVC"] = cv.get("CLNVC", ".") or "."
                row["ClinVar_ID"] = cv.get("RS", ".") or "."
            else:
                row["ClinVar_CLNSIG"] = "."
                row["ClinVar_CLNREVSTAT"] = "."
                row["ClinVar_CLNDN"] = "."
                row["ClinVar_CLNDISDB"] = "."
                row["ClinVar_CLNVC"] = "."
                row["ClinVar_ID"] = "."

            w.writerow(row)

    print(f"{os.path.basename(path_in)} -> {os.path.basename(path_out)} | clinvar_hits={hits}/{total}")

def main():
    ap = argparse.ArgumentParser(
        description="Join ClinVar lookup TSV into candidate TSV files by (CHROM,POS,REF,ALT)."
    )
    ap.add_argument("in_dir", help="Input directory with TSV files")
    ap.add_argument("lookup_tsv", help="ClinVar lookup TSV (CHROM POS REF ALT ...)")
    ap.add_argument("out_dir", help="Output directory")
    ap.add_argument("--pattern", default="*rare*.tsv", help="Glob pattern for input TSVs (default: *rare*.tsv)")
    ap.add_argument("--suffix", default="_clinvar.tsv", help="Output suffix (default: _clinvar.tsv)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.in_dir, args.pattern)))
    if not files:
        raise SystemExit(f"ERROR: no input TSV files matched: {os.path.join(args.in_dir, args.pattern)}")

    lookup = load_lookup(args.lookup_tsv)
    print(f"lookup_records={len(lookup)} from {args.lookup_tsv}")

    for path_in in files:
        base = os.path.basename(path_in)
        if base.endswith(".tsv"):
            base = base[:-4]
        path_out = os.path.join(args.out_dir, base + args.suffix)
        join_one(path_in, path_out, lookup)

if __name__ == "__main__":
    main()
