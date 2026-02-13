#!/usr/bin/env python3
import os
import sys
import csv

# args: IN_DIR OUT_REGIONS CHR_PREFIX
in_dir, out_regions, chr_prefix = sys.argv[1], sys.argv[2], sys.argv[3]

def norm_chrom(c: str) -> str:
    c = (c or "").strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c

regions = set()

for fn in os.listdir(in_dir):
    if not fn.endswith(".tsv"):
        continue
    if "rare" not in fn:
        continue
    path = os.path.join(in_dir, fn)
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames or "CHROM" not in r.fieldnames or "POS" not in r.fieldnames:
            continue
        for row in r:
            chrom = norm_chrom(row.get("CHROM", ""))
            pos = row.get("POS", "")
            if not chrom or not pos:
                continue
            chrom2 = (chr_prefix + chrom) if chr_prefix else chrom
            regions.add((chrom2, int(pos)))

with open(out_regions, "w", encoding="utf-8") as out:
    for chrom, pos in sorted(regions, key=lambda x: (x[0], x[1])):
        out.write(f"{chrom}\t{pos}\t{pos}\n")

print(f"regions={len(regions)} -> {out_regions}")
