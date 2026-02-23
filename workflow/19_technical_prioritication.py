#!/usr/bin/env python3
import csv
import os
import glob
import logging
from typing import Optional, Dict, Tuple, List

IN_CLINVAR = "results/14_clinvar"
# Берём пары из Step17 (у тебя он уже 353/353), но:
# - tier считаем на лету
# - clinvar_any считаем на лету
# - QC делаем по AB (если есть), потому что DP/GQ в pairs часто нет
IN_COMP = "results/17_comphet/comphet_pairs_with_evidence.tsv"

OUT_DIR = "results/19_master"
TOP_N = 100

LOG_DIR = "logs"
LOG_FILE = os.path.join(LOG_DIR, "19_master_scoring.log")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, encoding="utf-8"), logging.StreamHandler()],
)
logger = logging.getLogger("step19")

# ---------------------------
# Helpers
# ---------------------------

def parse_float(x) -> Optional[float]:
    if x in (None, "", "."):
        return None
    try:
        return float(x)
    except Exception:
        return None

def impact_score(imp: str) -> int:
    return {"HIGH": 4, "MODERATE": 2}.get((imp or "").upper(), 0)

def clinvar_score(sig: str) -> int:
    s = (sig or "").lower()

    if "conflicting" in s:
        return 1   # treat as VUS-like
    if "uncertain" in s or "vus" in s:
        return 1
    if "likely benign" in s or "likely_benign" in s:
        return -2
    if "benign" in s:
        return -4
    if "likely pathogenic" in s or "likely_pathogenic" in s:
        return 4
    if "pathogenic" in s:
        return 6

    return 0

def clinvar_bucket(sig: str) -> str:
    s = (sig or "").lower()

    if "conflicting" in s:
        return "VUS_OR_OTHER"
    if "uncertain" in s or "vus" in s:
        return "VUS_OR_OTHER"
    if "likely_benign" in s or "likely benign" in s:
        return "BENIGN"
    if "benign" in s:
        return "BENIGN"
    if "likely_pathogenic" in s:
        return "LIKELY_PATHOGENIC"
    if "pathogenic" in s:
        return "PATHOGENIC"

    return "NONE"

def af_score(af: Optional[float]) -> int:
    if af is None:
        return 0
    if af <= 1e-4:
        return 3
    if af <= 1e-3:
        return 2
    if af <= 1e-2:
        return 1
    return -3

def qc_variant(row: Dict[str, str]) -> bool:
    # Variant QC: proband DP/GQ
    try:
        dp = int(float(row.get("GEN[2].DP", "0") or 0))
        gq = int(float(row.get("GEN[2].GQ", "0") or 0))
    except Exception:
        return False
    return dp >= 10 and gq >= 20

def compute_pair_tier(pair_has_damaging: str, pair_top_impact: str) -> str:
    dmg = (pair_has_damaging or "").upper()
    imp = (pair_top_impact or "").upper()
    if dmg == "YES" and imp == "HIGH":
        return "T1"
    if dmg == "YES":
        return "T2"
    if imp == "HIGH":
        return "T3"
    return "T4"

def qc_comphet_pair(row: Dict[str, str]) -> bool:
    # Pair QC: use AB if available in pairs file.
    # Many pair files have: proband_ab_paternal / proband_ab_maternal
    abp = parse_float(row.get("proband_ab_paternal", row.get("proband_ab_paternal", ".")))
    abm = parse_float(row.get("proband_ab_maternal", row.get("proband_ab_maternal", ".")))

    # If AB is missing, don't fail hard (keep it, but it will rank by other signals).
    if abp is None or abm is None:
        return True
    return (0.25 <= abp <= 0.75) and (0.25 <= abm <= 0.75)

def safe_get(row: Dict[str, str], k: str, default=".") -> str:
    v = row.get(k, default)
    return (v if v is not None else default)

# ---------------------------
# Main
# ---------------------------

def main() -> None:
    candidates: List[Dict[str, object]] = []
    excluded: List[Tuple[str, str]] = []

    # ---- Variant-level inputs (from step14 clinvar outputs)
    variant_files = sorted(glob.glob(os.path.join(IN_CLINVAR, "*_variants_rare_*_clinvar.tsv")))
    logger.info("Variant inputs matched: %d files", len(variant_files))
    for p in variant_files:
        logger.info("  variant_in: %s", p)

    v_total = v_qc_fail = 0
    for fn in variant_files:
        with open(fn, newline="", encoding="utf-8") as f:
            r = csv.DictReader(f, delimiter="\t")
            for row in r:
                v_total += 1
                vid = safe_get(row, "_vid", ".")

                if vid in ("", "."):
                    chrom = safe_get(row, "CHROM", ".")
                    pos   = safe_get(row, "POS", ".")
                    ref   = safe_get(row, "REF", ".")
                    alt   = safe_get(row, "ALT", ".")
                    vid = f"{chrom}-{pos}-{ref}-{alt}"

                if not qc_variant(row):
                    v_qc_fail += 1
                    excluded.append(("variant_qc_fail", vid))
                    continue

                af = parse_float(row.get("gnomAD_max_af") or row.get("gnomAD_AF"))
                score = (
                    impact_score(row.get("ANN[0].IMPACT"))
                    + clinvar_score(row.get("ClinVar_CLNSIG"))
                    + af_score(af)
                    + 1  # qc bonus
                )

                candidates.append({
                    "kind": "variant",
                    "model": safe_get(row, "inheritance", "."),
                    "gene": safe_get(row, "ANN[0].GENE", "."),
                    "key": vid,
                    "score": int(score),
                    "impact": safe_get(row, "ANN[0].IMPACT", "."),
                    "clinvar": safe_get(row, "ClinVar_CLNSIG", "."),
                    "af": "." if af is None else af,
                    "source": os.path.basename(fn),
                })

    logger.info("Variants: total=%d qc_fail=%d kept=%d", v_total, v_qc_fail, len([c for c in candidates if c["kind"]=="variant"]))

    # ---- Comphet pairs (from step17)
    if not os.path.exists(IN_COMP):
        logger.warning("Comphet input missing: %s (skipping pairs)", IN_COMP)
    else:
        with open(IN_COMP, newline="", encoding="utf-8") as f:
            r = csv.DictReader(f, delimiter="\t")
            header = r.fieldnames or []
            logger.info("Comphet input: %s", IN_COMP)
            logger.info("Comphet columns present (sample): %s", ", ".join(header[:20]))

            p_total = p_qc_fail = p_tier_skip = 0
            kept_pairs = 0

            for row in r:
                p_total += 1

                # Compute tier on the fly (because pair_tier may not exist in step17 output)
                tier = compute_pair_tier(row.get("pair_has_damaging", ""), row.get("pair_top_impact", ""))

                # Keep only T1/T2 for master TOP (you can later keep all tiers in ALL if you want)
                if tier not in ("T1", "T2"):
                    p_tier_skip += 1
                    continue

                if not qc_comphet_pair(row):
                    p_qc_fail += 1
                    excluded.append(("comphet_qc_fail", safe_get(row, "GENE", ".")))
                    continue

                # AF for pair
                af = parse_float(row.get("pair_max_af", "."))
                # ClinVar any: compute from pat/mat CLNSIG if available
                pat_sig = row.get("pat_ClinVar_CLNSIG", row.get("pat_ClinVar_CLNSIG", "."))
                mat_sig = row.get("mat_ClinVar_CLNSIG", row.get("mat_ClinVar_CLNSIG", "."))
                pat_b = clinvar_bucket(pat_sig)
                mat_b = clinvar_bucket(mat_sig)
                if "PATHOGENIC" in (pat_b, mat_b):
                    cv_any = "PATHOGENIC"
                    cv_sig = "PATHOGENIC"
                    cv_points = 5
                elif "LIKELY_PATHOGENIC" in (pat_b, mat_b):
                    cv_any = "LIKELY_PATHOGENIC"
                    cv_sig = "LIKELY_PATHOGENIC"
                    cv_points = 3
                elif "BENIGN" in (pat_b, mat_b):
                    cv_any = "BENIGN"
                    cv_sig = "BENIGN"
                    cv_points = -3
                elif "VUS_OR_OTHER" in (pat_b, mat_b):
                    cv_any = "VUS_OR_OTHER"
                    cv_sig = "VUS_OR_OTHER"
                    cv_points = 1
                else:
                    cv_any = "NONE"
                    cv_sig = "NONE"
                    cv_points = 0

                tier_score = {"T1": 8, "T2": 5}.get(tier, 0)

                score = tier_score + cv_points + af_score(af) + 1  # small QC bonus

                gene = safe_get(row, "GENE", ".")
                pvid = safe_get(row, "paternal_vid", safe_get(row, "pat_vid", "."))
                mvid = safe_get(row, "maternal_vid", safe_get(row, "mat_vid", "."))

                # CRITICAL FIX: unique key per PAIR, not per gene
                pair_key = f"{gene}|{pvid}|{mvid}"

                candidates.append({
                    "kind": "pair",
                    "model": "comphet",
                    "gene": gene,
                    "key": pair_key,
                    "score": int(score),
                    "impact": safe_get(row, "pair_top_impact", "."),
                    "clinvar": cv_sig,
                    "af": "." if af is None else af,
                    "tier": tier,
                    "source": os.path.basename(IN_COMP),
                })
                kept_pairs += 1

            logger.info(
                "Comphet pairs: total=%d tier_skip=%d qc_fail=%d kept=%d",
                p_total, p_tier_skip, p_qc_fail, kept_pairs
            )

    # ---- Dedup (keep best score per (kind, model, key))
    before = len(candidates)
    unique: Dict[Tuple[str, str, str], Dict[str, object]] = {}
    for c in candidates:
        k = (str(c["kind"]), str(c["model"]), str(c["key"]))
        if k not in unique or int(c["score"]) > int(unique[k]["score"]):
            unique[k] = c
    candidates = list(unique.values())
    logger.info("Candidates before dedup: %d | after dedup: %d", before, len(candidates))

    if not candidates:
        raise SystemExit("ERROR: no candidates left after filtering/dedup. Check inputs/QC.")

    # ---- Sort (score desc, then AF asc; unknown AF treated as 'later')
    def af_sort(v):
        if v in (None, ".", ""):
            return 9.0
        try:
            return float(v)
        except:
            return 9.0

    candidates_sorted = sorted(candidates, key=lambda x: (-int(x["score"]), af_sort(x.get("af"))))

    # ---- Write outputs
    out_all = os.path.join(OUT_DIR, "master_scored_all.tsv")
    out_top = os.path.join(OUT_DIR, f"master_scored_top{TOP_N}.tsv")
    out_exc = os.path.join(OUT_DIR, "master_excluded_qc.tsv")

    header = list(candidates_sorted[0].keys())
    with open(out_all, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(candidates_sorted)

    with open(out_top, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(candidates_sorted[:TOP_N])

    with open(out_exc, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["reason", "key"])
        w.writerows(excluded)

    logger.info("Wrote ALL -> %s (rows=%d)", out_all, len(candidates_sorted))
    logger.info("Wrote TOP -> %s (rows=%d)", out_top, min(TOP_N, len(candidates_sorted)))
    logger.info("Wrote EXC -> %s (rows=%d)", out_exc, len(excluded))
    logger.info("Log -> %s", LOG_FILE)

if __name__ == "__main__":
    main()
