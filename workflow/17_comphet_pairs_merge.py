#!/usr/bin/env python3
"""
Step 17 — Merge evidence back to compound-het pairs.

Goal:
  Build a review-ready table where each row is a (gene, paternal_variant, maternal_variant) pair,
  enriched with annotations and evidence for BOTH variants.

Inputs:
  1) results/16_comphet/comphet_pairs.tsv
     - columns include: GENE, paternal_vid, maternal_vid
  2) results/14_clinvar/comphet_rare_maxaf0.001_clinvar.tsv
     - variant-level table with _vid and evidence columns (gnomAD_max_af, ClinVar_*, snpEff ANN fields, etc.)

Output:
  results/17_comphet/comphet_pairs_with_evidence.tsv
"""

import csv
import os
from typing import Dict, Optional, Tuple
import logging

PAIRS_TSV = "results/16_comphet/comphet_pairs.tsv"
VARIANTS_TSV = "results/14_clinvar/comphet_rare_maxaf0.001_clinvar.tsv"
OUT_DIR = "results/17_comphet"
OUT_TSV = os.path.join(OUT_DIR, "comphet_pairs_with_evidence.tsv")

os.makedirs(OUT_DIR, exist_ok=True)

# Severity helpers (pair-level summaries)
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
DAMAGING_MARKERS = ("splice", "frameshift", "stop_gained", "start_lost", "stop_lost")

LOG_DIR = "logs"
LOG_FILE = os.path.join(LOG_DIR, "17_comphet_pairs.log")
os.makedirs(LOG_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)

def parse_float(x: str) -> Optional[float]:
    x = (x or "").strip()
    if x == "" or x == ".":
        return None
    try:
        return float(x)
    except Exception:
        return None

def max2(a: Optional[float], b: Optional[float]) -> Optional[float]:
    if a is None:
        return b
    if b is None:
        return a
    return a if a >= b else b

def impact_rank(impact: str) -> int:
    return IMPACT_RANK.get((impact or "").upper(), -1)

def is_damaging(effect: str, impact: str) -> bool:
    if (impact or "").upper() == "HIGH":
        return True
    e = (effect or "").lower()
    return any(m in e for m in DAMAGING_MARKERS)

def load_variants_by_vid(path: str) -> Dict[str, Dict[str, str]]:
    """Load variant table keyed by _vid for O(1) joins."""
    variants: Dict[str, Dict[str, str]] = {}
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            raise SystemExit("ERROR: variants TSV has no header")
        if "_vid" not in r.fieldnames:
            raise SystemExit("ERROR: variants TSV must contain column: _vid")
        for row in r:
            variants[row["_vid"]] = row
    return variants

def pick(v: Optional[Dict[str, str]], key: str) -> str:
    """Safe getter for missing variants."""
    if not v:
        return "."
    return v.get(key, ".") or "."

def main():
    # 1) Load variant evidence table into memory
    variants = load_variants_by_vid(VARIANTS_TSV)
    logger.info(f"Loaded variants: {len(variants)} from {VARIANTS_TSV}")
    
    # 2) Read pairs and write enriched output
    with open(PAIRS_TSV, newline="", encoding="utf-8") as fin:
        pr = csv.DictReader(fin, delimiter="\t")
        if not pr.fieldnames:
            raise SystemExit("ERROR: pairs TSV has no header")

        required = {"GENE", "paternal_vid", "maternal_vid"}
        missing = sorted(list(required - set(pr.fieldnames)))
        if missing:
            raise SystemExit(f"ERROR: pairs TSV missing columns: {missing}")

        # Keep original pair columns + add pair-level + selected paternal/maternal evidence
        out_fields = list(pr.fieldnames) + [
            "pair_max_af",
            "pair_top_impact",
            "pair_has_damaging",

            "pat_effect", "pat_impact", "pat_hgvs_c", "pat_hgvs_p",
            "pat_gnomAD_max_af", "pat_gnomAD_status",
            "pat_ClinVar_CLNSIG", "pat_ClinVar_CLNREVSTAT", "pat_ClinVar_ID",

            "mat_effect", "mat_impact", "mat_hgvs_c", "mat_hgvs_p",
            "mat_gnomAD_max_af", "mat_gnomAD_status",
            "mat_ClinVar_CLNSIG", "mat_ClinVar_CLNREVSTAT", "mat_ClinVar_ID",
        ]

        with open(OUT_TSV, "w", newline="", encoding="utf-8") as fout:
            w = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()

            total = 0
            ok_both = 0

            for pair in pr:
                total += 1
                pvid = pair["paternal_vid"]
                mvid = pair["maternal_vid"]

                pv = variants.get(pvid)
                mv = variants.get(mvid)

                if pv and mv:
                    ok_both += 1

                # Variant-level fields (from your snpEff+evidence schema)
                pat_effect = pick(pv, "ANN[0].EFFECT")
                pat_impact = pick(pv, "ANN[0].IMPACT")
                pat_hgvs_c = pick(pv, "ANN[0].HGVS_C")
                pat_hgvs_p = pick(pv, "ANN[0].HGVS_P")

                mat_effect = pick(mv, "ANN[0].EFFECT")
                mat_impact = pick(mv, "ANN[0].IMPACT")
                mat_hgvs_c = pick(mv, "ANN[0].HGVS_C")
                mat_hgvs_p = pick(mv, "ANN[0].HGVS_P")

                # AF fields (computed earlier in 13B)
                pat_af = parse_float(pick(pv, "gnomAD_max_af"))
                mat_af = parse_float(pick(mv, "gnomAD_max_af"))
                pair_af = max2(pat_af, mat_af)

                # Pair-level summaries
                top_imp = pat_impact if impact_rank(pat_impact) >= impact_rank(mat_impact) else mat_impact
                has_dmg = is_damaging(pat_effect, pat_impact) or is_damaging(mat_effect, mat_impact)

                pair["pair_max_af"] = "." if pair_af is None else f"{pair_af:.6g}"
                pair["pair_top_impact"] = top_imp or "."
                pair["pair_has_damaging"] = "YES" if has_dmg else "NO"

                # Attach paternal evidence
                pair["pat_effect"] = pat_effect
                pair["pat_impact"] = pat_impact
                pair["pat_hgvs_c"] = pat_hgvs_c
                pair["pat_hgvs_p"] = pat_hgvs_p
                pair["pat_gnomAD_max_af"] = pick(pv, "gnomAD_max_af")
                pair["pat_gnomAD_status"] = pick(pv, "gnomAD_status")
                pair["pat_ClinVar_CLNSIG"] = pick(pv, "ClinVar_CLNSIG")
                pair["pat_ClinVar_CLNREVSTAT"] = pick(pv, "ClinVar_CLNREVSTAT")
                pair["pat_ClinVar_ID"] = pick(pv, "ClinVar_ID")

                # Attach maternal evidence
                pair["mat_effect"] = mat_effect
                pair["mat_impact"] = mat_impact
                pair["mat_hgvs_c"] = mat_hgvs_c
                pair["mat_hgvs_p"] = mat_hgvs_p
                pair["mat_gnomAD_max_af"] = pick(mv, "gnomAD_max_af")
                pair["mat_gnomAD_status"] = pick(mv, "gnomAD_status")
                pair["mat_ClinVar_CLNSIG"] = pick(mv, "ClinVar_CLNSIG")
                pair["mat_ClinVar_CLNREVSTAT"] = pick(mv, "ClinVar_CLNREVSTAT")
                pair["mat_ClinVar_ID"] = pick(mv, "ClinVar_ID")

                w.writerow(pair)

    logger.info(f"Pairs processed: {total}")
    logger.info(f"Pairs with both variants found: {ok_both}")
    logger.info(f"Wrote -> {OUT_TSV}")

if __name__ == "__main__":
    main()
