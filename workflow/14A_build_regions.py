#!/usr/bin/env python3
import argparse
import csv
import glob
import os
from typing import Set, Tuple

def norm_chrom(c: str) -> str:
    c = (c or "").strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c

def main():
    ap = argparse.ArgumentParser(
        description="Build a bcftools -R regions TSV from candidate TSV files (CHROM, POS)."
    )
    ap.add_argument("in_dir", help="Input directory with TSV files")
    ap.add_argument("out_regions", help="Output regions TSV (CHROM\\tSTART\\tEND)")
    ap.add_argument("chr_prefix", help="Prefix to add to contigs ('' or 'chr')")
    ap.add_argument("--pattern", default="*rare*.tsv", help="Glob pattern for input TSVs (default: *rare*.tsv)")
    args = ap.parse_args()

    in_dir = args.in_dir
    out_regions = args.out_regions
    chr_prefix = args.chr_prefix or ""

    if not os.path.isdir(in_dir):
        raise SystemExit(f"ERROR: input directory not found: {in_dir}")

    files = sorted(glob.glob(os.path.join(in_dir, args.pattern)))
    if not files:
        raise SystemExit(f"ERROR: no TSV files matched: {os.path.join(in_dir, args.pattern)}")

    regions: Set[Tuple[str, int]] = set()
    scanned = 0

    for path in files:
        scanned += 1
        with open(path, newline="", encoding="utf-8") as f:
            r = csv.DictReader(f, delimiter="\t")
            if not r.fieldnames or "CHROM" not in r.fieldnames or "POS" not in r.fieldnames:
                continue
            for row in r:
                chrom = norm_chrom(row.get("CHROM", ""))
                pos = row.get("POS", "")
                if not chrom or not pos:
                    continue
                try:
                    pos_i = int(pos)
                except ValueError:
                    continue
                chrom2 = (chr_prefix + chrom) if chr_prefix else chrom
                regions.add((chrom2, pos_i))

    with open(out_regions, "w", encoding="utf-8") as out:
        for chrom, pos in sorted(regions, key=lambda x: (x[0], x[1])):
            out.write(f"{chrom}\t{pos}\t{pos}\n")

    print(f"files_scanned={scanned} pattern={args.pattern} regions={len(regions)} -> {out_regions}")

if __name__ == "__main__":
    main()
