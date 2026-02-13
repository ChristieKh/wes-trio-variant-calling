#!/usr/bin/env python3
import csv
import os
from collections import defaultdict
from typing import Dict, List, Tuple

IN_TSV = "results/12_model_tables/ar_het_annotated.tsv"
OUT_DIR = "results/16_comphet"
os.makedirs(OUT_DIR, exist_ok=True)

OUT_PAIRS = os.path.join(OUT_DIR, "comphet_pairs.tsv")
OUT_VARIANTS = os.path.join(OUT_DIR, "comphet_variants.tsv")

def gt_is_het(gt: str) -> bool:
    gt = (gt or "").replace("|", "/")
    return gt in {"0/1", "1/0"}

def gt_is_hom_ref(gt: str) -> bool:
    gt = (gt or "").replace("|", "/")
    return gt == "0/0"

def make_vid(row: Dict[str, str]) -> str:
    return f"{row.get('CHROM','')}-{row.get('POS','')}-{row.get('REF','')}-{row.get('ALT','')}"

def inheritance_origin(f_gt: str, m_gt: str, p_gt: str) -> str:
    if not gt_is_het(p_gt):
        return "other"
    if gt_is_het(f_gt) and gt_is_hom_ref(m_gt):
        return "paternal"
    if gt_is_het(m_gt) and gt_is_hom_ref(f_gt):
        return "maternal"
    if gt_is_het(f_gt) and gt_is_het(m_gt):
        return "both_parents_het"
    return "other"

def main():
    with open(IN_TSV, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        header = r.fieldnames or []
        rows = list(r)

    gene_key = "ANN[0].GENE"
    f_gt_key, m_gt_key, p_gt_key = "GEN[0].GT", "GEN[1].GT", "GEN[2].GT"

    by_gene = defaultdict(list)
    for row in rows:
        gene = (row.get(gene_key) or "").strip()
        if not gene:
            continue
        origin = inheritance_origin(row.get(f_gt_key, ""), row.get(m_gt_key, ""), row.get(p_gt_key, ""))
        if origin in {"paternal", "maternal"}:
            row["_origin"] = origin
            row["_vid"] = make_vid(row)
            by_gene[gene].append(row)

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
                    "paternal_effect": pv.get("ANN[0].EFFECT", ""),
                    "maternal_effect": mv.get("ANN[0].EFFECT", ""),
                    "paternal_impact": pv.get("ANN[0].IMPACT", ""),
                    "maternal_impact": mv.get("ANN[0].IMPACT", ""),
                })
                variant_rows[pv["_vid"]] = pv
                variant_rows[mv["_vid"]] = mv

    pair_header = [
        "GENE",
        "paternal_vid", "maternal_vid",
        "paternal_pos", "maternal_pos",
        "paternal_effect", "maternal_effect",
        "paternal_impact", "maternal_impact",
    ]

    with open(OUT_PAIRS, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=pair_header, delimiter="\t")
        w.writeheader()
        w.writerows(pair_rows)

    # Write unique variant table (same schema as input + origin + vid)
    out_header = header + ["_origin", "_vid"]
    with open(OUT_VARIANTS, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=out_header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for vid, row in variant_rows.items():
            w.writerow(row)

    print(f"Input ar_het rows: {len(rows)}")
    print(f"Genes with paternal+maternal het: {len({p['GENE'] for p in pair_rows})}")
    print(f"Pairs: {len(pair_rows)} -> {OUT_PAIRS}")
    print(f"Unique variants in pairs: {len(variant_rows)} -> {OUT_VARIANTS}")

if __name__ == "__main__":
    main()
