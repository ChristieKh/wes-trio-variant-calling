#!/usr/bin/env python3
import csv
import os
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Tuple

IN_TSV = "results/12_model_tables/ar_het_annotated.tsv"
OUT_DIR = "results/16_comphet"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PAIRS = os.path.join(OUT_DIR, "comphet_pairs.tsv")
OUT_VARIANTS = os.path.join(OUT_DIR, "comphet_variants.tsv")

# --- QC thresholds  ---
MIN_PROBAND_DP = 10
MIN_PROBAND_GQ = 20
MIN_PARENT_DP = 8
MIN_AB = 0.25
MAX_AB = 0.75

# --- column names (from SnpSift extractFields) ---
GENE_KEY = "ANN[0].GENE"
EFF_KEY = "ANN[0].EFFECT"
IMPACT_KEY = "ANN[0].IMPACT"

F_GT, M_GT, P_GT = "GEN[0].GT", "GEN[1].GT", "GEN[2].GT"
F_DP, M_DP, P_DP = "GEN[0].DP", "GEN[1].DP", "GEN[2].DP"
P_GQ = "GEN[2].GQ"
P_AD = "GEN[2].AD"


def norm_gt(gt: str) -> str:
    return (gt or "").strip().replace("|", "/")


def gt_is_het(gt: str) -> bool:
    return norm_gt(gt) in {"0/1", "1/0"}


def gt_is_hom_ref(gt: str) -> bool:
    return norm_gt(gt) == "0/0"


def parse_int(x: str) -> Optional[int]:
    x = (x or "").strip()
    if not x or x == ".":
        return None
    try:
        return int(float(x))
    except ValueError:
        return None


def parse_ad(ad: str) -> Tuple[Optional[int], Optional[int]]:
    ad = (ad or "").strip()
    if not ad or ad == ".":
        return None, None
    parts = ad.split(",")
    if len(parts) < 2:
        return None, None
    try:
        return int(parts[0]), int(parts[1])
    except ValueError:
        return None, None


def calc_ab(ad: str) -> Optional[float]:
    refc, altc = parse_ad(ad)
    if refc is None or altc is None:
        return None
    denom = refc + altc
    if denom <= 0:
        return None
    return altc / denom


def make_vid(row: Dict[str, str]) -> str:
    return f"{row.get('CHROM','')}-{row.get('POS','')}-{row.get('REF','')}-{row.get('ALT','')}"


def origin(f_gt: str, m_gt: str, p_gt: str) -> str:
    if not gt_is_het(p_gt):
        return "other"
    if gt_is_het(f_gt) and gt_is_hom_ref(m_gt):
        return "paternal"
    if gt_is_het(m_gt) and gt_is_hom_ref(f_gt):
        return "maternal"
    if gt_is_het(f_gt) and gt_is_het(m_gt):
        return "both_parents_het"
    return "other"


def passes_qc(row: Dict[str, str]) -> bool:
    p_dp = parse_int(row.get(P_DP, ""))
    p_gq = parse_int(row.get(P_GQ, ""))
    f_dp = parse_int(row.get(F_DP, ""))
    m_dp = parse_int(row.get(M_DP, ""))

    if p_dp is None or p_dp < MIN_PROBAND_DP:
        return False
    if p_gq is None or p_gq < MIN_PROBAND_GQ:
        return False
    if f_dp is None or f_dp < MIN_PARENT_DP:
        return False
    if m_dp is None or m_dp < MIN_PARENT_DP:
        return False

    ab = calc_ab(row.get(P_AD, ""))
    if ab is None or ab < MIN_AB or ab > MAX_AB:
        return False
    row["_ab_p"] = f"{ab:.4f}"
    return True


def main():
    with open(IN_TSV, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        header = r.fieldnames or []
        rows = list(r)

    required = [GENE_KEY, F_GT, M_GT, P_GT, P_DP, P_GQ, P_AD, F_DP, M_DP]
    missing = [c for c in required if c not in header]
    if missing:
        raise SystemExit(f"ERROR: missing columns: {missing}")

    by_gene = defaultdict(list)
    qc_kept = 0

    for row in rows:
        gene = (row.get(GENE_KEY) or "").strip()
        if not gene:
            continue

        o = origin(row.get(F_GT, ""), row.get(M_GT, ""), row.get(P_GT, ""))
        if o not in {"paternal", "maternal"}:
            continue
        if not passes_qc(row):
            continue

        row["_origin"] = o
        row["_vid"] = make_vid(row)
        by_gene[gene].append(row)
        qc_kept += 1

    pair_rows: List[Dict[str, str]] = []
    variant_rows: Dict[str, Dict[str, str]] = {}

    for gene, vars_ in by_gene.items():
        paternal = [v for v in vars_ if v["_origin"] == "paternal"]
        maternal = [v for v in vars_ if v["_origin"] == "maternal"]
        if not paternal or not maternal:
            continue

        for pv in paternal:
            for mv in maternal:
                if pv["_vid"] == mv["_vid"]:
                    continue
                pair_rows.append({
                    "GENE": gene,
                    "paternal_vid": pv["_vid"],
                    "maternal_vid": mv["_vid"],
                    "paternal_pos": pv.get("POS", ""),
                    "maternal_pos": mv.get("POS", ""),
                    "paternal_effect": pv.get(EFF_KEY, ""),
                    "maternal_effect": mv.get(EFF_KEY, ""),
                    "paternal_impact": pv.get(IMPACT_KEY, ""),
                    "maternal_impact": mv.get(IMPACT_KEY, ""),
                    "proband_ab_paternal": pv.get("_ab_p", "."),
                    "proband_ab_maternal": mv.get("_ab_p", "."),
                })
                variant_rows[pv["_vid"]] = pv
                variant_rows[mv["_vid"]] = mv

    pair_header = [
        "GENE",
        "paternal_vid", "maternal_vid",
        "paternal_pos", "maternal_pos",
        "paternal_effect", "maternal_effect",
        "paternal_impact", "maternal_impact",
        "proband_ab_paternal", "proband_ab_maternal",
    ]
    with open(OUT_PAIRS, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=pair_header, delimiter="\t")
        w.writeheader()
        w.writerows(pair_rows)

    out_header = header + ["_origin", "_vid", "_ab_p"]
    with open(OUT_VARIANTS, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=out_header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for _, row in variant_rows.items():
            w.writerow(row)

    genes_with_pairs = len({p["GENE"] for p in pair_rows})
    print(f"Input ar_het rows: {len(rows)}")
    print(f"QC-kept paternal/maternal variants: {qc_kept}")
    print(f"Genes with paternal+maternal het: {genes_with_pairs}")
    print(f"Pairs: {len(pair_rows)} -> {OUT_PAIRS}")
    print(f"Unique variants in pairs: {len(variant_rows)} -> {OUT_VARIANTS}")
    print("Thresholds:",
          f"proband DP>={MIN_PROBAND_DP}, proband GQ>={MIN_PROBAND_GQ}, parent DP>={MIN_PARENT_DP}, AB in [{MIN_AB},{MAX_AB}]")

if __name__ == "__main__":
    main()
