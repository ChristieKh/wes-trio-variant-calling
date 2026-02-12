#!/usr/bin/env python3
import csv
import os
from typing import Optional

IN_DIR = "results/13_gnomad"
OUT_DIR = "results/13B_gnomad_filtered"

# Настройки порогов (можешь менять)
# Обычно для rare Mendelian:
# - de novo / X-linked: очень редкие (<= 1e-4 или 1e-3)
# - AR (homo/comp-het): можно чуть мягче (<= 1e-3 или 1e-2), но лучше начинать строго
THRESHOLDS = {
    "de_novo": 1e-4,
    "x_linked": 1e-4,
    "ar_homo": 1e-3,
    "ar_het": 1e-3,  # это только сырьё для pairing; потом станет ещё меньше
}

def parse_float(x: str) -> Optional[float]:
    if x is None:
        return None
    x = x.strip()
    if x == "" or x == ".":
        return None
    try:
        return float(x)
    except Exception:
        return None

def max_af(row) -> Optional[float]:
    vals = [
        parse_float(row.get("gnomAD_exome_af", ".")),
        parse_float(row.get("gnomAD_genome_af", ".")),
        parse_float(row.get("gnomAD_faf95_popmax", ".")),
    ]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None

def model_from_filename(fn: str) -> str:
    # ожидаем: de_novo_annotated_gnomad.tsv, ar_homo_..., ar_het_..., x_linked_...
    for m in ("de_novo", "ar_homo", "ar_het", "x_linked"):
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
    out_path = os.path.join(OUT_DIR, fn.replace("_annotated_gnomad.tsv", f"_rare_maxaf{thr}.tsv"))

    with open(in_path, newline="", encoding="utf-8") as fin, \
         open(out_path, "w", newline="", encoding="utf-8") as fout:

        r = csv.DictReader(fin, delimiter="\t")
        fieldnames = (r.fieldnames or []) + ["gnomAD_max_af", "gnomAD_filter"]

        w = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        total = kept = 0
        for row in r:
            total += 1
            maf = max_af(row)
            row["gnomAD_max_af"] = "." if maf is None else f"{maf:.6g}"

            # NOT_FOUND (maf None) обычно трактуем как "очень редкий" → оставляем
            keep = (maf is None) or (maf <= thr)

            row["gnomAD_filter"] = "PASS" if keep else "FAIL"
            if keep:
                kept += 1
                w.writerow(row)

    print(f"[{model}] {fn}: kept {kept}/{total} (max_af <= {thr} OR NOT_FOUND) -> {out_path}")

print("Step 13B complete.")
