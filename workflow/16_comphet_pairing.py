#!/usr/bin/env python3
import csv
import os
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Tuple

IN_TSV = "results/12_model_tables/ar_het_annotated.comphet.tsv"
OUT_DIR = "results/16_comphet"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PAIRS_STRICT = os.path.join(OUT_DIR, "comphet_pairs_strict.tsv")
OUT_PAIRS_CAND   = os.path.join(OUT_DIR, "comphet_pairs_candidate.tsv")
OUT_VARIANTS     = os.path.join(OUT_DIR, "comphet_variants.tsv")

# --- QC thresholds  ---
MIN_PROBAND_DP = 10
MIN_PROBAND_GQ = 20
MIN_PARENT_DP = 8
MIN_AB = 0.25
MAX_AB = 0.75

# --- column names ---
GENE_KEY = "ANN[0].GENE"
EFF_KEY = "ANN[0].EFFECT"
IMPACT_KEY = "ANN[0].IMPACT"

F_GT, M_GT, P_GT = "GEN[0].GT", "GEN[1].GT", "GEN[2].GT"
F_DP, M_DP, P_DP = "GEN[0].DP", "GEN[1].DP", "GEN[2].DP"
P_GQ = "GEN[2].GQ"
P_AD = "GEN[2].AD"

ALLOW_PARENT_NOCALL_AS_REF = True  # <- ключевая "мягкость"


def norm_gt(gt: str) -> str:
    return (gt or "").strip().replace("|", "/")

def gt_is_het(gt: str) -> bool:
    return norm_gt(gt) in {"0/1", "1/0"}

def gt_is_hom_ref(gt: str) -> bool:
    return norm_gt(gt) == "0/0"

def gt_is_nocall(gt: str) -> bool:
    g = norm_gt(gt)
    return g in {".", "./.", ".|."}

def parent_ref_like(gt: str) -> bool:
    # treat no-call as "ref-like" if option enabled
    if gt_is_hom_ref(gt):
        return True
    if ALLOW_PARENT_NOCALL_AS_REF and gt_is_nocall(gt):
        return True
    return False

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
        return int(float(parts[0])), int(float(parts[1]))
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

def passes_qc(row: Dict[str, str]) -> bool:
    p_dp = parse_int(row.get(P_DP, ""))
    p_gq = parse_int(row.get(P_GQ, ""))
    f_dp = parse_int(row.get(F_DP, ""))
    m_dp = parse_int(row.get(M_DP, ""))

    # For candidate mode: if we allow parent nocall, parent DP might be low / missing.
    # But we still require some minimal DP when present.
    if p_dp is None or p_dp < MIN_PROBAND_DP:
        return False
    if p_gq is None or p_gq < MIN_PROBAND_GQ:
        return False

    if f_dp is not None and f_dp < MIN_PARENT_DP:
        return False
    if m_dp is not None and m_dp < MIN_PARENT_DP:
        return False

    ab = calc_ab(row.get(P_AD, ""))
    if ab is None or ab < MIN_AB or ab > MAX_AB:
        return False
    row["_ab_p"] = f"{ab:.4f}"
    return True

def origin_strict(f_gt: str, m_gt: str, p_gt: str) -> str:
    if not gt_is_het(p_gt):
        return "other"
    if gt_is_het(f_gt) and gt_is_hom_ref(m_gt):
        return "paternal"
    if gt_is_het(m_gt) and gt_is_hom_ref(f_gt):
        return "maternal"
    if gt_is_het(f_gt) and gt_is_het(m_gt):
        return "both_parents_het"
    return "other"

def origin_candidate(f_gt: str, m_gt: str, p_gt: str) -> str:
    # same as strict, but allow other parent to be "ref-like" (0/0 or nocall)
    if not gt_is_het(p_gt):
        return "other"
    if gt_is_het(f_gt) and parent_ref_like(m_gt):
        return "paternal"
    if gt_is_het(m_gt) and parent_ref_like(f_gt):
        return "maternal"
    if gt_is_het(f_gt) and gt_is_het(m_gt):
        return "both_parents_het"
    return "other"

def build_pairs(by_gene: Dict[str, List[Dict[str, str]]], mode_label: str) -> Tuple[List[Dict[str, str]], Dict[str, Dict[str, str]]]:
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
                    "MODE": mode_label,
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

    return pair_rows, variant_rows

def write_tsv(path: str, header: List[str], rows: List[Dict[str, str]]):
    with open(path, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

def main():
    with open(IN_TSV, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        header = r.fieldnames or []
        rows = list(r)

    required = [GENE_KEY, F_GT, M_GT, P_GT, P_DP, P_GQ, P_AD]
    missing = [c for c in required if c not in header]
    if missing:
        raise SystemExit(f"ERROR: missing columns: {missing}")

    # ---- stats before QC ----
    origin_counts_strict_pre = Counter()
    origin_counts_cand_pre = Counter()

    # ---- collect after QC ----
    by_gene_strict = defaultdict(list)
    by_gene_cand = defaultdict(list)

    qc_kept_strict = 0
    qc_kept_cand = 0

    for row in rows:
        gene = (row.get(GENE_KEY) or "").strip()
        if not gene:
            continue

        o_strict = origin_strict(row.get(F_GT, ""), row.get(M_GT, ""), row.get(P_GT, ""))
        o_cand = origin_candidate(row.get(F_GT, ""), row.get(M_GT, ""), row.get(P_GT, ""))
        origin_counts_strict_pre[o_strict] += 1
        origin_counts_cand_pre[o_cand] += 1

        if not passes_qc(row):
            continue

        # STRICT
        if o_strict in {"paternal", "maternal"}:
            r1 = dict(row)
            r1["_origin"] = o_strict
            r1["_vid"] = make_vid(r1)
            by_gene_strict[gene].append(r1)
            qc_kept_strict += 1

        # CANDIDATE
        if o_cand in {"paternal", "maternal"}:
            r2 = dict(row)
            r2["_origin"] = o_cand
            r2["_vid"] = make_vid(r2)
            by_gene_cand[gene].append(r2)
            qc_kept_cand += 1

    # ---- build pairs ----
    strict_pairs, strict_variants = build_pairs(by_gene_strict, "strict")
    cand_pairs, cand_variants = build_pairs(by_gene_cand, "candidate")

    # Merge unique variants (candidate superset usually)
    all_variants = {}
    all_variants.update(strict_variants)
    all_variants.update(cand_variants)

    pair_header = [
        "MODE", "GENE",
        "paternal_vid", "maternal_vid",
        "paternal_pos", "maternal_pos",
        "paternal_effect", "maternal_effect",
        "paternal_impact", "maternal_impact",
        "proband_ab_paternal", "proband_ab_maternal",
    ]
    write_tsv(OUT_PAIRS_STRICT, pair_header, strict_pairs)
    write_tsv(OUT_PAIRS_CAND, pair_header, cand_pairs)

    out_header = header + ["_origin", "_vid", "_ab_p"]
    with open(OUT_VARIANTS, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=out_header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for _, row in all_variants.items():
            w.writerow(row)

    # ---- stats ----
    strict_gene_p = {g for g, vs in by_gene_strict.items() if any(v["_origin"]=="paternal" for v in vs)}
    strict_gene_m = {g for g, vs in by_gene_strict.items() if any(v["_origin"]=="maternal" for v in vs)}

    cand_gene_p = {g for g, vs in by_gene_cand.items() if any(v["_origin"]=="paternal" for v in vs)}
    cand_gene_m = {g for g, vs in by_gene_cand.items() if any(v["_origin"]=="maternal" for v in vs)}

    print(f"Input ar_het rows: {len(rows)}")
    print("")
    print("Origin counts (pre-QC):")
    print("  strict   :", dict(origin_counts_strict_pre))
    print("  candidate:", dict(origin_counts_cand_pre))
    print("")
    print(f"QC-kept paternal/maternal variants: strict={qc_kept_strict} candidate={qc_kept_cand}")
    print(f"Genes paternal: strict={len(strict_gene_p)} maternal={len(strict_gene_m)} intersection={len(strict_gene_p & strict_gene_m)}")
    print(f"Genes paternal: cand  ={len(cand_gene_p)} maternal={len(cand_gene_m)} intersection={len(cand_gene_p & cand_gene_m)}")
    print("")
    print(f"Pairs strict   : {len(strict_pairs)} -> {OUT_PAIRS_STRICT}")
    print(f"Pairs candidate: {len(cand_pairs)} -> {OUT_PAIRS_CAND}")
    print(f"Unique variants in any pairs: {len(all_variants)} -> {OUT_VARIANTS}")
    print("Thresholds:",
          f"proband DP>={MIN_PROBAND_DP}, proband GQ>={MIN_PROBAND_GQ}, parent DP>={MIN_PARENT_DP} (if present), AB in [{MIN_AB},{MAX_AB}]",
          f"allow_parent_nocall_as_ref={ALLOW_PARENT_NOCALL_AS_REF}")

if __name__ == "__main__":
    main()
