#!/usr/bin/env python3
import csv
import os
from typing import Optional, Tuple, List

IN_DIR = "results/final_candidates"
OUT_DIR = "results/final_candidates"
OUT_TSV = os.path.join(OUT_DIR, "final_candidates.tsv")

FILES = [
    ("de_novo",       os.path.join(IN_DIR, "de_novo_annotated.tsv")),
    ("inherited_het", os.path.join(IN_DIR, "inherited_het_annotated.tsv")),
    ("recessive",     os.path.join(IN_DIR, "recessive_annotated.tsv")),
]

# Keep only meaningful effects (first-pass)
KEEP_EFFECT_SUBSTR = (
    "missense_variant",
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
)

DROP_EFFECT_SUBSTR = (
    "intergenic_region",
)

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3, ".": 9, "": 9}


def parse_ad(ad: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse AD like 'REF,ALT' -> (ref, alt).
    Returns (None, None) if missing.
    """
    if not ad or ad == ".":
        return None, None
    parts = ad.split(",")
    if len(parts) < 2:
        return None, None
    try:
        ref = int(float(parts[0]))
        alt = int(float(parts[1]))
        return ref, alt
    except ValueError:
        return None, None


def calc_ab(ref: Optional[int], alt: Optional[int]) -> Optional[float]:
    if ref is None or alt is None:
        return None
    total = ref + alt
    if total <= 0:
        return None
    return alt / total


def keep_row(effect: str, impact: str) -> bool:
    e = (effect or "").lower()
    imp = (impact or "").upper()

    # Drop obvious non-coding/noise
    for bad in DROP_EFFECT_SUBSTR:
        if bad in e:
            return False
    if "modifier" in (imp or "").lower():
        return False

    # Keep only key coding-like effects (first-pass)
    if not any(k in e for k in KEEP_EFFECT_SUBSTR):
        return False

    # Usually keep HIGH/MODERATE; keep LOW only if you want broader list
    if imp not in ("HIGH", "MODERATE"):
        return False

    return True


def sort_key(r: dict) -> Tuple[int, int, int]:
    # de_novo first, then recessive, then inherited
    inh_rank = {"de_novo": 0, "recessive": 1, "inherited_het": 2}.get(r["inheritance"], 9)
    impact_rank = IMPACT_RANK.get((r.get("ANN[0].IMPACT") or "").upper(), 9)
    # higher AB_proband is better -> sort ascending by negative (i.e., descending AB)
    ab = r.get("AB_proband")
    ab_rank = 9 if ab in ("", ".") else -float(ab)
    return (inh_rank, impact_rank, ab_rank)


def main() -> None:
    os.makedirs(OUT_DIR, exist_ok=True)

    rows: List[dict] = []
    base_header: Optional[List[str]] = None

    for inh, path in FILES:
        if not os.path.exists(path):
            continue
        with open(path, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
            if base_header is None:
                base_header = header

            # map column index
            idx = {name: i for i, name in enumerate(header)}

            def get(col: str, rec: List[str]) -> str:
                i = idx.get(col)
                return rec[i] if i is not None and i < len(rec) else "."

            for rec in reader:
                row = {h: (rec[i] if i < len(rec) else ".") for i, h in enumerate(header)}
                row["inheritance"] = inh

                # Evidence from proband AD (GEN[2].AD in your extract)
                ad_p = get("GEN[2].AD", rec)
                ref_p, alt_p = parse_ad(ad_p)
                ab_p = calc_ab(ref_p, alt_p)

                row["ALT_reads_proband"] = "." if alt_p is None else str(alt_p)
                row["AB_proband"] = "." if ab_p is None else f"{ab_p:.3f}"

                effect = row.get("ANN[0].EFFECT", ".")
                impact = row.get("ANN[0].IMPACT", ".")

                if keep_row(effect, impact):
                    rows.append(row)

    if base_header is None:
        raise SystemExit("No input TSVs found. Run Step 12A first.")

    # Output header: add our computed columns + inheritance
    out_header = ["inheritance"] + base_header + ["ALT_reads_proband", "AB_proband"]

    # Sort
    rows.sort(key=sort_key)

    with open(OUT_TSV, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_header, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Step 12B done: {OUT_TSV}")
    print(f"Kept variants (after filters): {len(rows)}")


if __name__ == "__main__":
    main()
