#!/usr/bin/env python3
"""
Step 18 — Build automatic shortlist for comphet pairs.

Input:
  results/17_comphet/comphet_pairs_with_evidence.tsv

Outputs:
  results/18_shortlists/comphet_pairs_priority.tsv        (all pairs, sorted, with priority columns)
  results/18_shortlists/comphet_pairs_shortlist_top50.tsv (top N pairs)

Ranking logic (simple + robust):
  Tier 1: pair_has_damaging=YES AND pair_top_impact=HIGH
  Tier 2: pair_has_damaging=YES (impact not HIGH)
  Tier 3: pair_top_impact=HIGH (damaging=NO)
  Tier 4: everything else

Within a tier, sort by:
  - ClinVar pathogenic signal present (any variant) first
  - then by pair_max_af ascending (missing '.' treated as very rare -> ranks higher)
  - then by gene name
"""

import csv
import os
import logging
from typing import Optional, Tuple

IN_TSV = "results/17_comphet/comphet_pairs_with_evidence.tsv"
OUT_DIR = "results/18_shortlists"
OUT_ALL = os.path.join(OUT_DIR, "comphet_pairs_priority.tsv")
TOP_N = 50
OUT_TOP = os.path.join(OUT_DIR, f"comphet_pairs_shortlist_top{TOP_N}.tsv")

LOG_DIR = "logs"
LOG_FILE = os.path.join(LOG_DIR, "18_comphet_shortlist.log")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, encoding="utf-8"), logging.StreamHandler()],
)
logger = logging.getLogger("step18_comphet_shortlist")

IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}

PATHOGENIC_TOKENS = (
    "pathogenic",          # includes Likely_pathogenic if normalized to lower
    "likely_pathogenic",
)
BENIGN_TOKENS = (
    "benign",
    "likely_benign",
)

def parse_float(x: str) -> Optional[float]:
    x = (x or "").strip()
    if x == "" or x == ".":
        return None
    try:
        return float(x)
    except Exception:
        return None

def impact_rank(x: str) -> int:
    return IMPACT_RANK.get((x or "").upper(), -1)

def clinvar_signal(clnsig: str) -> str:
    """
    Returns one of: PATHOGENIC, BENIGN, VUS_OR_OTHER, NONE
    """
    s = (clnsig or "").strip()
    if s == "" or s == ".":
        return "NONE"
    low = s.lower()
    # quick normalization: ClinVar may contain pipe/semicolon-separated values
    # We just look for key tokens.
    if any(tok in low for tok in PATHOGENIC_TOKENS):
        return "PATHOGENIC"
    if any(tok in low for tok in BENIGN_TOKENS):
        return "BENIGN"
    return "VUS_OR_OTHER"

def compute_tier(pair_has_damaging: str, pair_top_impact: str) -> str:
    dmg = (pair_has_damaging or "").upper()
    imp = (pair_top_impact or "").upper()
    if dmg == "YES" and imp == "HIGH":
        return "T1"
    if dmg == "YES":
        return "T2"
    if imp == "HIGH":
        return "T3"
    return "T4"

def tier_rank(tier: str) -> int:
    return {"T1": 1, "T2": 2, "T3": 3, "T4": 4}.get(tier, 9)

def sort_key(row: dict) -> Tuple:
    tier = row.get("pair_tier", "T4")
    # ClinVar: any pathogenic on either variant -> best
    cv_any = row.get("pair_clinvar_any", "NONE")
    cv_rank = {"PATHOGENIC": 0, "VUS_OR_OTHER": 1, "BENIGN": 2, "NONE": 3}.get(cv_any, 9)

    # AF: smaller is better; missing treated as ultra-rare (rank first)
    af = parse_float(row.get("pair_max_af", "."))
    af_rank = -1.0 if af is None else af

    gene = row.get("GENE", ".")
    # Keep stable ordering within same scores
    return (tier_rank(tier), cv_rank, af_rank, gene)

def main() -> None:
    if not os.path.exists(IN_TSV):
        raise SystemExit(f"ERROR: missing input: {IN_TSV}")

    with open(IN_TSV, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            raise SystemExit("ERROR: input TSV has no header")

        needed = {"GENE", "pair_max_af", "pair_top_impact", "pair_has_damaging"}
        missing = sorted(list(needed - set(r.fieldnames)))
        if missing:
            raise SystemExit(f"ERROR: input missing required columns: {missing}")

        rows = list(r)

    logger.info(f"Loaded pairs: {len(rows)} from {IN_TSV}")

    # Add priority columns
    for row in rows:
        tier = compute_tier(row.get("pair_has_damaging"), row.get("pair_top_impact"))
        row["pair_tier"] = tier

        # ClinVar aggregation across both variants
        pat_sig = clinvar_signal(row.get("pat_ClinVar_CLNSIG", "."))
        mat_sig = clinvar_signal(row.get("mat_ClinVar_CLNSIG", "."))

        # "any" signal: PATHOGENIC wins, then VUS, then BENIGN, then NONE
        if "PATHOGENIC" in (pat_sig, mat_sig):
            any_sig = "PATHOGENIC"
        elif "VUS_OR_OTHER" in (pat_sig, mat_sig):
            any_sig = "VUS_OR_OTHER"
        elif "BENIGN" in (pat_sig, mat_sig):
            any_sig = "BENIGN"
        else:
            any_sig = "NONE"

        row["pair_clinvar_pat"] = pat_sig
        row["pair_clinvar_mat"] = mat_sig
        row["pair_clinvar_any"] = any_sig

    # Sort all
    rows_sorted = sorted(rows, key=sort_key)

    # Write all (priority)
    out_fields = list(rows_sorted[0].keys())
    # Ensure priority fields appear near the front
    for k in ["pair_tier", "pair_clinvar_any", "pair_clinvar_pat", "pair_clinvar_mat"]:
        if k in out_fields:
            out_fields.remove(k)
    out_fields = ["pair_tier", "pair_clinvar_any", "pair_clinvar_pat", "pair_clinvar_mat"] + out_fields

    with open(OUT_ALL, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows_sorted:
            w.writerow(row)

    # Write top N
    top = rows_sorted[:TOP_N]
    with open(OUT_TOP, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in top:
            w.writerow(row)

    # Quick stats
    from collections import Counter
    tier_counts = Counter([x.get("pair_tier","?") for x in rows_sorted])
    cv_counts = Counter([x.get("pair_clinvar_any","?") for x in rows_sorted])

    logger.info(f"Wrote ALL -> {OUT_ALL}")
    logger.info(f"Wrote TOP -> {OUT_TOP}")
    logger.info(f"Tier counts: {dict(tier_counts)}")
    logger.info(f"ClinVar(any) counts: {dict(cv_counts)}")
    logger.info(f"Log -> {LOG_FILE}")

if __name__ == "__main__":
    main()
