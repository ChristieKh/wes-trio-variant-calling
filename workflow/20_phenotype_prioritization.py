#!/usr/bin/env python3
import csv
import os
import argparse
import logging
from typing import Dict, List, Optional, Tuple

DEFAULT_MASTER = "results/19_master/master_scored_all.tsv"
DEFAULT_PHENO_DIR = "results/20_phenotype"
DEFAULT_OUT_DIR = "results/20_interpretation"
DEFAULT_LOG = "logs/20_phenotype_prioritization.log"


# ----------------------------
# Helpers
# ----------------------------

def setup_logger(log_path: str) -> logging.Logger:
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path, encoding="utf-8"), logging.StreamHandler()],
    )
    return logging.getLogger("step20")


def norm_gene(g: str) -> str:
    return (g or "").strip().upper()


def parse_float(x) -> Optional[float]:
    if x is None:
        return None
    s = str(x).strip()
    if s in ("", ".", "NA", "nan", "None"):
        return None
    try:
        return float(s)
    except Exception:
        return None


def load_gene_panels(path: str, logger: logging.Logger) -> Dict[str, Dict[str, str]]:
    """
    TSV columns (required): gene, category, weight
    optional: comment
    """
    if not os.path.exists(path):
        raise SystemExit(f"ERROR: gene_panels.tsv not found: {path}")

    panels: Dict[str, Dict[str, str]] = {}
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            raise SystemExit("ERROR: gene_panels.tsv has no header")
        header_l = [h.strip().lower() for h in r.fieldnames]
        required = {"gene", "category", "weight"}
        if not required.issubset(set(header_l)):
            raise SystemExit(
                f"ERROR: gene_panels.tsv must have columns: gene, category, weight (got: {r.fieldnames})"
            )

        for row in r:
            g = norm_gene(row.get("gene", ""))
            if not g:
                continue
            panels[g] = {
                "category": (row.get("category", "") or "").strip(),
                "weight": (row.get("weight", "") or "").strip(),
                "comment": (row.get("comment", "") or "").strip(),
            }

    logger.info("Loaded gene_panels: %d genes", len(panels))
    return panels


def detect_tech_score_col(header: List[str]) -> str:
    """
    Choose technical score column from master.
    Prefer: score_final -> score_tech -> score
    """
    for c in ("score_final", "score_tech", "score"):
        if c in header:
            return c
    return ""


# ----------------------------
# Noise penalty (strong)
# ----------------------------

NOISY_PREFIXES = ("KRTAP", "KRT", "MUC", "HLA", "OR")
NOISY_GENES = {"CCHCR1", "PSORS1C1"}  # extend as needed


def is_noisy_gene(gene: str) -> bool:
    g = norm_gene(gene)
    if not g or g == ".":
        return False
    if g in NOISY_GENES:
        return True
    return any(g.startswith(p) for p in NOISY_PREFIXES)


def noise_penalty(gene: str) -> int:
    # strong penalty (not exclusion)
    return -10 if is_noisy_gene(gene) else 0


# ----------------------------
# Main
# ----------------------------

def main():
    ap = argparse.ArgumentParser(description="Step20: phenotype-driven final shortlist from Step19 master")
    ap.add_argument("--master", default=DEFAULT_MASTER, help="Input master_scored_all.tsv (Step19)")
    ap.add_argument("--pheno_dir", default=DEFAULT_PHENO_DIR, help="Directory with phenotype files")
    ap.add_argument("--out_dir", default=DEFAULT_OUT_DIR, help="Output directory")
    ap.add_argument("--top_n", type=int, default=20, help="How many rows to keep in final shortlist")
    ap.add_argument("--log", default=DEFAULT_LOG, help="Log file path")
    args = ap.parse_args()

    logger = setup_logger(args.log)

    master_path = args.master
    pheno_dir = args.pheno_dir
    out_dir = args.out_dir

    os.makedirs(out_dir, exist_ok=True)

    panels_path = os.path.join(pheno_dir, "gene_panels.tsv")
    summary_path = os.path.join(pheno_dir, "phenotype_summary.md")

    if not os.path.exists(summary_path):
        logger.warning("phenotype_summary.md missing: %s (recommended for documentation)", summary_path)
    else:
        logger.info("Phenotype summary: %s", summary_path)

    panels = load_gene_panels(panels_path, logger)

    if not os.path.exists(master_path):
        raise SystemExit(f"ERROR: master file not found: {master_path}")

    # Read master
    with open(master_path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            raise SystemExit("ERROR: master TSV has no header.")
        header = r.fieldnames

        tech_col = detect_tech_score_col(header)
        if not tech_col:
            raise SystemExit("ERROR: master TSV must contain score_final or score_tech or score.")
        logger.info("Using technical score column: %s", tech_col)

        rows = list(r)

    logger.info("Master rows: %d", len(rows))

    enriched: List[Dict[str, str]] = []

    for row in rows:
        gene_raw = row.get("gene", row.get("GENE", "."))  # defensive
        gene = norm_gene(gene_raw)
        model = (row.get("model", "") or row.get("inheritance", "") or ".").strip()
        kind = (row.get("kind", ".") or ".").strip()

        tech = parse_float(row.get(tech_col, "0")) or 0.0

        # phenotype bonus from gene_panels
        panel = panels.get(gene, None)
        panel_w = int(parse_float(panel["weight"]) or 0) if panel else 0
        panel_cat = panel["category"] if panel else "."
        panel_comment = panel["comment"] if panel else "."

        # inheritance nudge for AR (tiny)
        inh = 0
        ml = model.lower()
        if ml in ("ar_homo", "comphet"):
            inh = 1
        elif ml == "de_novo":
            inh = -1

        pen = noise_penalty(gene)

        final = tech + panel_w + inh + pen

        out = dict(row)
        out["tech_score_used"] = tech_col
        out["tech_score"] = str(int(tech)) if float(tech).is_integer() else str(tech)

        out["panel_category"] = panel_cat
        out["panel_weight"] = str(panel_w)
        out["inheritance_bonus"] = str(inh)
        out["noise_penalty"] = str(pen)

        out["final_score"] = str(int(final)) if float(final).is_integer() else str(final)
        out["phenotype_note"] = (panel_comment or panel_cat or ".")[:180]
        out["is_noisy_gene"] = "1" if is_noisy_gene(gene) else "0"

        enriched.append(out)

    # Sort by final_score desc, then AF asc (unknown AF last), then tech desc
    def af_sort(v) -> float:
        x = parse_float(v)
        return x if x is not None else 9.0

    enriched_sorted = sorted(
        enriched,
        key=lambda x: (
            -float(parse_float(x.get("final_score", "0")) or 0.0),
            af_sort(x.get("af", ".")),
            -float(parse_float(x.get("tech_score", "0")) or 0.0),
        ),
    )

    final_rows = enriched_sorted[: args.top_n]

    # Write final shortlist (single file)
    out_final = os.path.join(out_dir, "final_shortlist.tsv")
    out_md = os.path.join(out_dir, "top_candidates.md")

    # stable output header
    base_cols = []
    for c in ["kind", "model", "gene", "key", "impact", "clinvar", "af", "tier", "source"]:
        if c in header:
            base_cols.append(c)
    # keep whatever exists but ensure our appended columns are present
    appended = [
        "tech_score_used", "tech_score",
        "panel_category", "panel_weight",
        "inheritance_bonus", "noise_penalty",
        "final_score", "is_noisy_gene", "phenotype_note"
    ]
    out_header = list(dict.fromkeys((header or []) + appended))

    with open(out_final, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(final_rows)

    # Write MD scaffold
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("# Final shortlist (phenotype-driven)\n\n")
        f.write(f"- Master: `{master_path}`\n")
        f.write(f"- Phenotype: `{summary_path}`\n")
        f.write(f"- Gene panels: `{panels_path}`\n\n")
        f.write("## Candidates\n\n")
        for i, row in enumerate(final_rows, start=1):
            f.write(f"### {i}. {row.get('gene','.')}\n\n")
            f.write(f"- model/kind: {row.get('model','.')}, {row.get('kind','.')}\n")
            f.write(f"- key: `{row.get('key','.')}`\n")
            f.write(f"- impact: `{row.get('impact','.')}` | ClinVar: `{row.get('clinvar','.')}` | AF: `{row.get('af','.')}`\n")
            f.write(f"- scores: tech={row.get('tech_score','.')}, panel={row.get('panel_weight','.')}, inh={row.get('inheritance_bonus','.')}, noise={row.get('noise_penalty','.')}, final={row.get('final_score','.')}\n")
            f.write(f"- phenotype note: {row.get('phenotype_note','.')}\n")
            f.write("- IGV: [ ] proband genotype  [ ] parents segregation  [ ] coverage/artefacts\n")
            f.write("- ACMG draft: \n")
            f.write("- Phenotype fit (why/why not): \n\n")

    logger.info("Wrote FINAL shortlist -> %s (rows=%d)", out_final, len(final_rows))
    logger.info("Wrote MD scaffold -> %s", out_md)
    logger.info("Log -> %s", args.log)


if __name__ == "__main__":
    main()
