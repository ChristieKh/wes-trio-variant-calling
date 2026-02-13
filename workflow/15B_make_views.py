#!/usr/bin/env python3
import csv
import os
from typing import List, Dict

IN_DIR = "results/15_shortlists"
OUT_DIR = "results/15_shortlists"
os.makedirs(OUT_DIR, exist_ok=True)

FILES = [
    "de_novo_top200.tsv",
    "ar_homo_top200.tsv",
]

def load_tsv(path: str) -> (List[Dict[str, str]], List[str]):
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        return list(r), (r.fieldnames or [])

def write_tsv(path: str, header: List[str], rows: List[Dict[str, str]]):
    with open(path, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def is_clinvar_pathogenic(sig: str) -> bool:
    s = (sig or "").lower()
    # include Likely_pathogenic and Pathogenic; exclude benign-only
    return ("pathogenic" in s) and ("benign" not in s)

def is_high_impact(impact: str) -> bool:
    return (impact or "").upper() == "HIGH"

def main():
    summary_lines = []
    for fn in FILES:
        in_path = os.path.join(IN_DIR, fn)
        if not os.path.exists(in_path):
            continue

        rows, header = load_tsv(in_path)
        if not rows:
            continue

        model = fn.replace("_top200.tsv", "")

        clinvar_rows = [r for r in rows if is_clinvar_pathogenic(r.get("ClinVar_CLNSIG", "."))]
        high_rows = [r for r in rows if is_high_impact(r.get("ANN[0].IMPACT", ""))]

        out_clinvar = os.path.join(OUT_DIR, f"{model}_clinvar_path.tsv")
        out_high = os.path.join(OUT_DIR, f"{model}_high_impact.tsv")

        write_tsv(out_clinvar, header, clinvar_rows)
        write_tsv(out_high, header, high_rows)

        summary_lines.append(f"{model}: total={len(rows)} clinvar_path={len(clinvar_rows)} high_impact={len(high_rows)}")
        print(f"{model}: wrote {len(clinvar_rows)} -> {out_clinvar}")
        print(f"{model}: wrote {len(high_rows)} -> {out_high}")

    summary_path = os.path.join(OUT_DIR, "summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")
    print(f"Wrote -> {summary_path}")

if __name__ == "__main__":
    main()
