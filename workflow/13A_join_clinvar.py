#!/usr/bin/env python3
import argparse
import csv
from typing import Dict, Tuple

Key = Tuple[str, str, str, str]  # CHROM, POS, REF, ALT


def norm_chrom(c: str) -> str:
    c = c or ""
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c


def load_clinvar(path: str) -> Dict[Key, Dict[str, str]]:
    m: Dict[Key, Dict[str, str]] = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            chrom = norm_chrom(row["CHROM"])
            key: Key = (chrom, row["POS"], row["REF"], row["ALT"])
            m[key] = {
                "CLNSIG": row.get("CLNSIG", ".") or ".",
                "CLNREVSTAT": row.get("CLNREVSTAT", ".") or ".",
                "CLNDN": row.get("CLNDN", ".") or ".",
            }
    return m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--clinvar", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    clin = load_clinvar(args.clinvar)

    kept = 0
    hit = 0

    with open(args.candidates, newline="") as inp, open(args.out, "w", newline="") as out:
        r = csv.DictReader(inp, delimiter="\t")
        fieldnames = list(r.fieldnames or [])
        add_cols = ["CLNSIG", "CLNREVSTAT", "CLNDN"]
        out_fields = fieldnames + add_cols

        w = csv.DictWriter(out, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        for row in r:
            kept += 1
            chrom = norm_chrom(row.get("CHROM", "."))
            key: Key = (chrom, row.get("POS", "."), row.get("REF", "."), row.get("ALT", "."))
            ann = clin.get(key)
            if ann:
                hit += 1
                row.update(ann)
            else:
                row.update({c: "." for c in add_cols})
            w.writerow(row)

    print(f"Candidates: {kept}")
    print(f"ClinVar hits (normalized CHROM + exact POS:REF:ALT): {hit}")


if __name__ == "__main__":
    main()
