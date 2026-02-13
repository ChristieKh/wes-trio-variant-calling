#!/usr/bin/env python3
import csv
import os
from typing import Dict, List, Optional, Tuple

IN_DIR = "results/14_clinvar"
OUT_DIR = "results/15_shortlists"
os.makedirs(OUT_DIR, exist_ok=True)

TOP_N = 200

DROP_EFFECT_SUBSTR = (
    "intergenic_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intron_variant",
    "5_prime_utr_variant",
    "3_prime_utr_variant",
    "non_coding_transcript",
)

def parse_float(x: str) -> Optional[float]:
    if x is None:
        return None
    x = x.strip()
    if x in ("", "."):
        return None
    try:
        return float(x)
    except Exception:
        return None

def gt_norm(gt: str) -> str:
    return (gt or "").replace("|", "/")

def gt_is_hom_alt(gt: str) -> bool:
    return gt_norm(gt) == "1/1"

def gq_int(x: str) -> int:
    try:
        return int(float(x))
    except Exception:
        return 0

def gq_sum(row: Dict[str, str]) -> int:
    return gq_int(row.get("GEN[0].GQ", "0")) + gq_int(row.get("GEN[1].GQ", "0")) + gq_int(row.get("GEN[2].GQ", "0"))

def should_drop_by_effect(text: str) -> bool:
    t = (text or "").lower()
    return any(k in t for k in DROP_EFFECT_SUBSTR)

def extract_ann_impact(row: Dict[str, str]) -> str:
    # Prefer explicit column if it's clean; otherwise parse from ANN-like text
    imp = (row.get("ANN[0].IMPACT") or "").strip().upper()
    if imp in {"HIGH", "MODERATE", "LOW", "MODIFIER"}:
        return imp

    ann_text = row.get("ANN[0].EFFECT", "") or row.get("ANN", "")
    if not ann_text or ann_text == ".":
        return ""
    first = ann_text.split(",")[0]
    parts = first.split("|")
    if len(parts) >= 3:
        return parts[2].strip().upper()
    return ""

def impact_rank(imp: str) -> int:
    return {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}.get((imp or "").upper(), 0)

def clinvar_bucket(sig: str) -> str:
    s = (sig or "").lower()
    if "conflicting" in s:
        return "conflicting"
    if "benign" in s and "pathogenic" not in s:
        return "benign"
    if "pathogenic" in s:
        if "likely" in s:
            return "likely_pathogenic"
        return "pathogenic"
    if "uncertain" in s or "vus" in s:
        return "vus"
    return "none"

def clinvar_rank(sig: str) -> int:
    b = clinvar_bucket(sig)
    return {
        "pathogenic": 4,
        "likely_pathogenic": 3,
        "vus": 1,
        "conflicting": 1,
        "none": 0,
        "benign": -2,
    }.get(b, 0)

def sort_key(row: Dict[str, str]) -> Tuple[int, int, float, int]:
    maf = parse_float(row.get("gnomAD_max_af", "."))
    maf_key = -(maf if maf is not None else 1.0)  # smaller AF is better
    imp = extract_ann_impact(row)
    return (
        clinvar_rank(row.get("ClinVar_CLNSIG", ".")),
        impact_rank(imp),
        maf_key,
        gq_sum(row),
    )

def write_tsv(path: str, header: List[str], rows: List[Dict[str, str]]):
    with open(path, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    summary = []

    for fn in sorted(os.listdir(IN_DIR)):
        if not fn.endswith("_clinvar.tsv"):
            continue

        in_path = os.path.join(IN_DIR, fn)
        with open(in_path, newline="", encoding="utf-8") as f:
            r = csv.DictReader(f, delimiter="\t")
            header = r.fieldnames or []
            rows = list(r)

        if not rows:
            continue

        model = (rows[0].get("inheritance") or "").strip() or fn.replace("_clinvar.tsv", "")

        total = len(rows)
        pass_rows = [x for x in rows if x.get("gnomAD_filter", "") == "PASS"]
        after_pass = len(pass_rows)

        # Views (no phenotype)
        cleaned = []
        for row in pass_rows:
            ann_text = row.get("ANN[0].EFFECT", "") or row.get("ANN", "")
            if should_drop_by_effect(ann_text):
                continue
            cleaned.append(row)

        imp_high = [x for x in cleaned if extract_ann_impact(x) == "HIGH"]
        imp_mod  = [x for x in cleaned if extract_ann_impact(x) == "MODERATE"]

        cv_path = [x for x in cleaned if clinvar_bucket(x.get("ClinVar_CLNSIG", ".")) in {"pathogenic", "likely_pathogenic"}]
        cv_conf = [x for x in cleaned if clinvar_bucket(x.get("ClinVar_CLNSIG", ".")) == "conflicting"]

        # TopN
        ranked = cleaned[:]
        ranked.sort(key=sort_key, reverse=True)
        topn = ranked[:TOP_N]

        stem = fn.replace("_clinvar.tsv", "")
        write_tsv(os.path.join(OUT_DIR, f"{stem}_top{TOP_N}.tsv"), header, topn)
        write_tsv(os.path.join(OUT_DIR, f"{stem}_clinvar_path.tsv"), header, cv_path)
        write_tsv(os.path.join(OUT_DIR, f"{stem}_clinvar_conflicting.tsv"), header, cv_conf)
        write_tsv(os.path.join(OUT_DIR, f"{stem}_high_impact.tsv"), header, imp_high)
        write_tsv(os.path.join(OUT_DIR, f"{stem}_moderate.tsv"), header, imp_mod)

        summary.append(
            f"{stem} | model={model} total={total} pass={after_pass} cleaned={len(cleaned)} "
            f"clinvar_path={len(cv_path)} clinvar_conflicting={len(cv_conf)} high={len(imp_high)} moderate={len(imp_mod)} top={len(topn)}"
        )

    summary_path = os.path.join(OUT_DIR, "summary_step15.txt")
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("\n".join(summary) + "\n")

    print(f"Wrote summary -> {summary_path}")
    for line in summary:
        print(line)

if __name__ == "__main__":
    main()
