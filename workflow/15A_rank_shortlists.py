#!/usr/bin/env python3
import csv
import os
from typing import List, Dict, Tuple

IN_DIR = "results/14_clinvar"
OUT_DIR = "results/15_shortlists"
os.makedirs(OUT_DIR, exist_ok=True)

MODELS = {"de_novo", "ar_homo"}  # keep this tight for now

def model_from_inheritance(row: Dict[str, str]) -> str:
    return (row.get("inheritance") or "").strip()

def parse_gq(x: str) -> int:
    try:
        return int(float(x))
    except Exception:
        return 0

def gq_sum(row: Dict[str, str]) -> int:
    return parse_gq(row.get("GEN[0].GQ", "0")) + parse_gq(row.get("GEN[1].GQ", "0")) + parse_gq(row.get("GEN[2].GQ", "0"))

def impact_rank(impact: str) -> int:
    imp = (impact or "").upper()
    return {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}.get(imp, 0)

def clinvar_rank(clnsig: str) -> int:
    s = (clnsig or "").lower()
    # simple, conservative ordering
    if "pathogenic" in s and "benign" not in s:
        if "likely" in s:
            return 3
        return 4
    if "uncertain" in s or "vus" in s:
        return 2
    if "benign" in s:
        return 0
    return 1  # includes "." or not provided

def sort_key(row: Dict[str, str]) -> Tuple[int, int, int, int]:
    # higher is better (we sort reverse=True)
    return (
        clinvar_rank(row.get("ClinVar_CLNSIG", ".")),
        impact_rank(row.get("ANN[0].IMPACT", "")),
        gq_sum(row),
        # push explicit FAIL down even if present
        1 if (row.get("gnomAD_filter", "") == "PASS") else 0,
    )

def load_rows(path: str) -> Tuple[List[Dict[str, str]], List[str]]:
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        header = r.fieldnames or []
        rows = list(r)
    return rows, header

def write_rows(path: str, header: List[str], rows: List[Dict[str, str]]):
    with open(path, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    # Collect rows per model across any *_clinvar.tsv files
    per_model: Dict[str, List[Dict[str, str]]] = {m: [] for m in MODELS}
    header_ref: List[str] = []

    for fn in os.listdir(IN_DIR):
        if not fn.endswith("_clinvar.tsv"):
            continue
        path = os.path.join(IN_DIR, fn)
        rows, header = load_rows(path)
        if not header_ref:
            header_ref = header

        for row in rows:
            m = model_from_inheritance(row)
            if m not in MODELS:
                continue
            # rely on your existing 13B result
            if row.get("gnomAD_filter", "") != "PASS":
                continue
            per_model[m].append(row)

    for m, rows in per_model.items():
        rows.sort(key=sort_key, reverse=True)
        top200 = rows[:200]
        out_path = os.path.join(OUT_DIR, f"{m}_top200.tsv")
        write_rows(out_path, header_ref, top200)
        print(f"{m}: input_pass={len(rows)} wrote={len(top200)} -> {out_path}")

if __name__ == "__main__":
    main()
