#!/usr/bin/env python3
import csv
import os
from typing import Optional

IN_DIR = "results/13_gnomad"
OUT_DIR = "results/13B_gnomad_filtered"

THRESHOLDS = {
    "de_novo": 1e-4,
    "x_linked": 1e-4,
    "ar_homo": 1e-3,
    "ar_het": 1e-3,
    "comphet": 1e-3,
}

def parse_float(x: str) -> Optional[float]:
    if not x or x == ".":
        return None
    try:
        return float(x)
    except:
        return None

def max_af(row) -> Optional[float]:
    vals = [
        parse_float(row.get("gnomAD_exome_af")),
        parse_float(row.get("gnomAD_genome_af")),
        parse_float(row.get("gnomAD_faf95_popmax")),
    ]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None

def model_from_filename(fn: str) -> str:
    for m in THRESHOLDS.keys():
        if fn.startswith(m + "_"):
            return m
    return "unknown"

os.makedirs(OUT_DIR, exist_ok=True)

for fn in os.listdir(IN_DIR):
    if not fn.endswith("_annotated_gnomad.tsv"):
        continue

    model = model_from_filename(fn)
    thr = THRESHOLDS.get(model, 1e-3)

    in_path = os.path.join(IN_DIR, fn)
    out_path = os.path.join(
        OUT_DIR,
        fn.replace("_annotated_gnomad.tsv", f"_rare_maxaf{thr}.tsv")
    )

    with open(in_path, newline="", encoding="utf-8") as fin, \
         open(out_path, "w", newline="", encoding="utf-8") as fout:

        r = csv.DictReader(fin, delimiter="\t")
        fields = r.fieldnames + ["gnomAD_max_af", "gnomAD_filter"]

        w = csv.DictWriter(fout, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        total = kept = 0

        for row in r:
            total += 1
            status = row.get("gnomAD_status", "")
            maf = max_af(row)

            row["gnomAD_max_af"] = "." if maf is None else f"{maf:.6g}"

            if status == "OK":
                keep = (maf is None) or (maf <= thr)
                row["gnomAD_filter"] = "PASS" if keep else "FAIL"

            elif status == "NOT_FOUND":
                keep = True
                row["gnomAD_filter"] = "PASS"

            else:
                # API errors: keep but flag
                keep = True
                row["gnomAD_filter"] = "API_ERROR"

            if keep:
                kept += 1
                w.writerow(row)

    print(f"[{model}] {fn}: kept {kept}/{total} -> {out_path}")

print("Step 13B complete.")
