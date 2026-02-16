#!/usr/bin/env python3
import csv
import os
import logging
from typing import Dict, List, Optional, Tuple

IN_DIR = "results/14_clinvar"
OUT_DIR = "results/15_shortlists"
LOG_FILE = "logs/15_shortlists.log"

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs("logs", exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, encoding="utf-8"), logging.StreamHandler()],
)
log = logging.getLogger("step15")

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

API_ERROR_STATUSES = {"GRAPHQL_ERROR", "HTTP_ERROR", "ERROR", "TIMEOUT"}

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
    return (
        gq_int(row.get("GEN[0].GQ", "0")) +
        gq_int(row.get("GEN[1].GQ", "0")) +
        gq_int(row.get("GEN[2].GQ", "0"))
    )

def should_drop_by_effect(effect: str) -> bool:
    t = (effect or "").lower()
    return any(k in t for k in DROP_EFFECT_SUBSTR)

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

def max_af(row) -> Optional[float]:
    vals = [
        parse_float(row.get("gnomAD_exome_af", ".")),
        parse_float(row.get("gnomAD_genome_af", ".")),
        parse_float(row.get("gnomAD_faf95_popmax", ".")),
        parse_float(row.get("gnomAD_max_af", ".")),  # if already computed
    ]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None

def api_penalty(row: Dict[str, str]) -> int:
    # Slight penalty so API_ERROR doesn't dominate TopN
    st = (row.get("gnomAD_status") or "").strip()
    return -1 if st in API_ERROR_STATUSES else 0

def sort_key(row: Dict[str, str]) -> Tuple[int, int, float, int, int]:
    maf = max_af(row)
    maf_sort = 0.0 if maf is None else maf          # None treated as rare (best)
    imp = (row.get("ANN[0].IMPACT") or "").strip().upper()
    return (
        clinvar_rank(row.get("ClinVar_CLNSIG", ".")),
        impact_rank(imp),
        -maf_sort,                                  # smaller AF => larger key
        gq_sum(row),
        api_penalty(row),
    )

def write_tsv(path: str, header: List[str], rows: List[Dict[str, str]]):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    summary_lines = []

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
        model_dir = os.path.join(OUT_DIR, model)
        os.makedirs(model_dir, exist_ok=True)

        total = len(rows)

        # Keep PASS and API_ERROR, but drop explicit FAIL if present
        pass_rows = []
        for x in rows:
            gf = (x.get("gnomAD_filter") or "").strip()
            if gf == "FAIL":
                continue
            pass_rows.append(x)

        # Optional model-specific sanity: ar_homo must be hom-alt in proband
        if model == "ar_homo":
            pass_rows = [x for x in pass_rows if gt_is_hom_alt(x.get("GEN[2].GT", ""))]

        cleaned = []
        for row in pass_rows:
            effect = row.get("ANN[0].EFFECT", "") or ""
            if should_drop_by_effect(effect):
                continue
            cleaned.append(row)

        imp_high = [x for x in cleaned if (x.get("ANN[0].IMPACT") or "").upper() == "HIGH"]
        imp_mod  = [x for x in cleaned if (x.get("ANN[0].IMPACT") or "").upper() == "MODERATE"]

        cv_path = [x for x in cleaned if clinvar_bucket(x.get("ClinVar_CLNSIG", ".")) in {"pathogenic", "likely_pathogenic"}]
        cv_conf = [x for x in cleaned if clinvar_bucket(x.get("ClinVar_CLNSIG", ".")) == "conflicting"]

        ranked = cleaned[:]
        ranked.sort(key=sort_key, reverse=True)
        topn = ranked[:TOP_N]

        stem = fn.replace("_clinvar.tsv", "")

        write_tsv(os.path.join(model_dir, f"{stem}_top{TOP_N}.tsv"), header, topn)
        write_tsv(os.path.join(model_dir, f"{stem}_clinvar_path.tsv"), header, cv_path)
        write_tsv(os.path.join(model_dir, f"{stem}_clinvar_conflicting.tsv"), header, cv_conf)
        write_tsv(os.path.join(model_dir, f"{stem}_high_impact.tsv"), header, imp_high)
        write_tsv(os.path.join(model_dir, f"{stem}_moderate.tsv"), header, imp_mod)

        summary_lines.append(
            f"{stem} | model={model} total={total} kept_after_gnomad={len(pass_rows)} cleaned={len(cleaned)} "
            f"clinvar_path={len(cv_path)} clinvar_conflicting={len(cv_conf)} high={len(imp_high)} moderate={len(imp_mod)} top={len(topn)}"
        )

        log.info(summary_lines[-1])

    summary_path = os.path.join(OUT_DIR, "summary_step15.txt")
    with open(summary_path, "w", encoding="utf-8") as out:
        out.write("\n".join(summary_lines) + "\n")

    log.info(f"Wrote summary -> {summary_path}")
    log.info(f"Log -> {LOG_FILE}")

if __name__ == "__main__":
    main()
