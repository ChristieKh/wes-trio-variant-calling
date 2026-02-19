#!/usr/bin/env python3
"""
Step 17 — Merge evidence back to compound-het pairs (schema-aligned).

Inputs:
  1) results/16_comphet/comphet_pairs_strict.tsv
     - must contain: GENE, paternal_vid, maternal_vid
  2) results/14_clinvar/comphet_variants_rare_af0.001_clinvar.tsv
     - must contain: _vid plus evidence columns:
       gnomAD_max_af, gnomAD_filter, ClinVar_*, ANN[0].*, GEN[*].*, _ab_p, _origin, etc.

Output:
  results/17_comphet/comphet_pairs_with_evidence.tsv
"""

import csv
import os
import logging
from typing import Dict, Optional

PAIRS_TSV = "results/16_comphet/comphet_pairs_strict.tsv"
VARIANTS_TSV = "results/14_clinvar/comphet_variants_rare_af0.001_clinvar.tsv"
OUT_DIR = "results/17_comphet"
OUT_TSV = os.path.join(OUT_DIR, "comphet_pairs_with_evidence.tsv")

LOG_DIR = "logs"
LOG_FILE = os.path.join(LOG_DIR, "17_comphet_pairs.log")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, encoding="utf-8"), logging.StreamHandler()],
)
logger = logging.getLogger("step17_comphet")


IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
DAMAGING_MARKERS = ("splice", "frameshift", "stop_gained", "start_lost", "stop_lost")


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


def pick(v: Optional[Dict[str, str]], key: str) -> str:
    if not v:
        return "."
    val = v.get(key, ".")
    return (val if val not in (None, "") else ".")


def load_variants_by_vid(path: str) -> Dict[str, Dict[str, str]]:
    variants: Dict[str, Dict[str, str]] = {}
    dups = 0
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            raise SystemExit("ERROR: variants TSV has no header")
        if "_vid" not in r.fieldnames:
            raise SystemExit("ERROR: variants TSV must contain column: _vid")

        for row in r:
            vid = (row.get("_vid") or "").strip()
            if not vid:
                continue
            if vid in variants:
                dups += 1
                # keep first occurrence for stability
                continue
            variants[vid] = row

    logger.info(f"Loaded variants: {len(variants)} (dups skipped: {dups}) from {path}")
    return variants


def main() -> None:
    variants = load_variants_by_vid(VARIANTS_TSV)

    with open(PAIRS_TSV, newline="", encoding="utf-8") as fin:
        pr = csv.DictReader(fin, delimiter="\t")
        if not pr.fieldnames:
            raise SystemExit("ERROR: pairs TSV has no header")

        required = {"GENE", "paternal_vid", "maternal_vid"}
        missing = sorted(list(required - set(pr.fieldnames)))
        if missing:
            raise SystemExit(f"ERROR: pairs TSV missing columns: {missing}")

        out_fields = list(pr.fieldnames) + [
            # pair-level
            "pair_max_af",
            "pair_top_impact",
            "pair_has_damaging",

            # paternal selected evidence
            "pat_vid",
            "pat_CHROM", "pat_POS", "pat_REF", "pat_ALT",
            "pat_effect", "pat_impact", "pat_hgvs_c", "pat_hgvs_p",
            "pat_gt_proband", "pat_ad_proband", "pat_dp_proband", "pat_gq_proband", "pat_ab_proband",
            "pat_gnomAD_max_af", "pat_gnomAD_filter",
            "pat_ClinVar_CLNSIG", "pat_ClinVar_CLNREVSTAT", "pat_ClinVar_ID",

            # maternal selected evidence
            "mat_vid",
            "mat_CHROM", "mat_POS", "mat_REF", "mat_ALT",
            "mat_effect", "mat_impact", "mat_hgvs_c", "mat_hgvs_p",
            "mat_gt_proband", "mat_ad_proband", "mat_dp_proband", "mat_gq_proband", "mat_ab_proband",
            "mat_gnomAD_max_af", "mat_gnomAD_filter",
            "mat_ClinVar_CLNSIG", "mat_ClinVar_CLNREVSTAT", "mat_ClinVar_ID",
        ]

        with open(OUT_TSV, "w", newline="", encoding="utf-8") as fout:
            w = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()

            total = 0
            ok_both = 0
            miss_pat = 0
            miss_mat = 0

            for pair in pr:
                total += 1
                pvid = (pair.get("paternal_vid") or "").strip()
                mvid = (pair.get("maternal_vid") or "").strip()

                pv = variants.get(pvid)
                mv = variants.get(mvid)

                if pv is None:
                    miss_pat += 1
                if mv is None:
                    miss_mat += 1
                if pv is not None and mv is not None:
                    ok_both += 1

                # Variant-level fields
                pat_effect = pick(pv, "ANN[0].EFFECT")
                pat_impact = pick(pv, "ANN[0].IMPACT")
                pat_hgvs_c = pick(pv, "ANN[0].HGVS_C")
                pat_hgvs_p = pick(pv, "ANN[0].HGVS_P")

                mat_effect = pick(mv, "ANN[0].EFFECT")
                mat_impact = pick(mv, "ANN[0].IMPACT")
                mat_hgvs_c = pick(mv, "ANN[0].HGVS_C")
                mat_hgvs_p = pick(mv, "ANN[0].HGVS_P")

                # AF fields (computed in 13B; present in your schema)
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
                pair["pat_vid"] = pvid or "."
                pair["pat_CHROM"] = pick(pv, "CHROM")
                pair["pat_POS"] = pick(pv, "POS")
                pair["pat_REF"] = pick(pv, "REF")
                pair["pat_ALT"] = pick(pv, "ALT")
                pair["pat_effect"] = pat_effect
                pair["pat_impact"] = pat_impact
                pair["pat_hgvs_c"] = pat_hgvs_c
                pair["pat_hgvs_p"] = pat_hgvs_p

                # proband is GEN[2] in your schema
                pair["pat_gt_proband"] = pick(pv, "GEN[2].GT")
                pair["pat_ad_proband"] = pick(pv, "GEN[2].AD")
                pair["pat_dp_proband"] = pick(pv, "GEN[2].DP")
                pair["pat_gq_proband"] = pick(pv, "GEN[2].GQ")
                pair["pat_ab_proband"] = pick(pv, "_ab_p")

                pair["pat_gnomAD_max_af"] = pick(pv, "gnomAD_max_af")
                pair["pat_gnomAD_filter"] = pick(pv, "gnomAD_filter")
                pair["pat_ClinVar_CLNSIG"] = pick(pv, "ClinVar_CLNSIG")
                pair["pat_ClinVar_CLNREVSTAT"] = pick(pv, "ClinVar_CLNREVSTAT")
                pair["pat_ClinVar_ID"] = pick(pv, "ClinVar_ID")

                # Attach maternal evidence
                pair["mat_vid"] = mvid or "."
                pair["mat_CHROM"] = pick(mv, "CHROM")
                pair["mat_POS"] = pick(mv, "POS")
                pair["mat_REF"] = pick(mv, "REF")
                pair["mat_ALT"] = pick(mv, "ALT")
                pair["mat_effect"] = mat_effect
                pair["mat_impact"] = mat_impact
                pair["mat_hgvs_c"] = mat_hgvs_c
                pair["mat_hgvs_p"] = mat_hgvs_p

                pair["mat_gt_proband"] = pick(mv, "GEN[2].GT")
                pair["mat_ad_proband"] = pick(mv, "GEN[2].AD")
                pair["mat_dp_proband"] = pick(mv, "GEN[2].DP")
                pair["mat_gq_proband"] = pick(mv, "GEN[2].GQ")
                pair["mat_ab_proband"] = pick(mv, "_ab_p")

                pair["mat_gnomAD_max_af"] = pick(mv, "gnomAD_max_af")
                pair["mat_gnomAD_filter"] = pick(mv, "gnomAD_filter")
                pair["mat_ClinVar_CLNSIG"] = pick(mv, "ClinVar_CLNSIG")
                pair["mat_ClinVar_CLNREVSTAT"] = pick(mv, "ClinVar_CLNREVSTAT")
                pair["mat_ClinVar_ID"] = pick(mv, "ClinVar_ID")

                w.writerow(pair)

    logger.info(f"Pairs processed: {total}")
    logger.info(f"Pairs with both variants found: {ok_both}")
    logger.info(f"Missing paternal variant rows: {miss_pat}")
    logger.info(f"Missing maternal variant rows: {miss_mat}")
    logger.info(f"Wrote -> {OUT_TSV}")
    logger.info(f"Log -> {LOG_FILE}")


if __name__ == "__main__":
    main()



head -n 1 results/16_comphet/comphet_pairs_strict.tsv
cut -f2-3 results/16_comphet/comphet_pairs_strict.tsv | head
# или точнее:
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="paternal_vid") p=i; if($i=="maternal_vid") m=i; next}
NR<=6{print $p"\n"$m}' results/16_comphet/comphet_pairs_strict.tsv


head -n 1 results/14_clinvar/comphet_variants_rare_af0.001_clinvar.tsv
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="_vid") v=i; next}
NR<=10{print $v}' results/14_clinvar/comphet_variants_rare_af0.001_clinvar.tsv


python3 - <<'PY'
import csv

pairs="results/16_comphet/comphet_pairs_strict.tsv"
vars="results/14_clinvar/comphet_variants_rare_af0.001_clinvar.tsv"

def read_vids_pairs(path):
    vids=set()
    with open(path, newline='', encoding='utf-8') as f:
        r=csv.DictReader(f, delimiter='\t')
        for row in r:
            vids.add((row.get("paternal_vid") or "").strip())
            vids.add((row.get("maternal_vid") or "").strip())
    vids.discard("")
    return vids

def read_vids_vars(path):
    vids=set()
    with open(path, newline='', encoding='utf-8') as f:
        r=csv.DictReader(f, delimiter='\t')
        for row in r:
            vids.add((row.get("_vid") or "").strip())
    vids.discard("")
    return vids

p=read_vids_pairs(pairs)
v=read_vids_vars(vars)
print("pairs unique vids:", len(p))
print("variants vids:", len(v))
print("intersection:", len(p & v))
print("example only-in-pairs:", list(sorted(p - v))[:5])
print("example only-in-variants:", list(sorted(v - p))[:5])
PY
