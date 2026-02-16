#!/usr/bin/env python3
import csv
import os
from typing import Optional, Tuple, Dict, Callable

IN_DIR = "results/09_candidates"
OUT_DIR = "results/10_candidates_evidence"

FILES = {
    "de_novo": "de_novo.tsv",
    "ar_homo": "ar_homo.tsv",
    "ar_het": "ar_het.tsv",
    "x_linked": "x_linked.tsv",
}

os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------------
# Canonical contigs
# -----------------------------
def norm_chrom(chrom: str) -> str:
    c = (chrom or "").strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c

def is_canonical_chrom(chrom: str) -> bool:
    c = norm_chrom(chrom)
    if c in {"X", "Y", "MT"}:
        return True
    if c.isdigit():
        n = int(c)
        return 1 <= n <= 22
    return False

# -----------------------------
# Parsing helpers
# -----------------------------
def parse_ad(ad_str: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse AD like "ref,alt" (may contain floats depending on upstream tools).
    Returns (ref, alt) as ints or (None, None).
    """
    if not ad_str or ad_str == ".":
        return None, None
    parts = ad_str.split(",")
    if len(parts) < 2:
        return None, None
    try:
        ref = int(float(parts[0]))
        alt = int(float(parts[1]))
        return ref, alt
    except Exception:
        return None, None

def allele_balance(ref: Optional[int], alt: Optional[int]) -> Optional[float]:
    if ref is None or alt is None:
        return None
    tot = ref + alt
    if tot <= 0:
        return None
    return alt / tot

def fmt_float(x: Optional[float]) -> str:
    return "." if x is None else f"{x:.4f}"

def safe_int(x: str) -> Optional[int]:
    x = (x or "").strip()
    if x == "" or x == ".":
        return None
    try:
        return int(float(x))
    except Exception:
        return None

# -----------------------------
# Thresholds / rules config
# -----------------------------
THR = {
    "de_novo": {"min_alt_p": 5, "min_ab_p": 0.25, "max_alt_parent": 0},
    "ar_homo": {"min_alt_p": 8, "min_ab_p": 0.85},
    "ar_het":  {"min_alt_p": 3, "min_ab_p": 0.20, "max_ab_p": 0.80},
    "x_linked":{"min_alt_p": 5, "min_ab_p": 0.25, "min_alt_m": 2},
}

def rule_de_novo(row: Dict[str, str]) -> Tuple[bool, str]:
    rf, af = parse_ad(row.get("AD_father", "."))
    rm, am = parse_ad(row.get("AD_mother", "."))
    rp, ap = parse_ad(row.get("AD_proband", "."))
    ab_p = allele_balance(rp, ap)

    if ap is None or ap < THR["de_novo"]["min_alt_p"]:
        return False, "ALT_PROBAND_TOO_LOW"
    if ab_p is None or ab_p < THR["de_novo"]["min_ab_p"]:
        return False, "AB_PROBAND_TOO_LOW"
    if af is None or am is None:
        return False, "PARENT_AD_MISSING"
    if af > THR["de_novo"]["max_alt_parent"] or am > THR["de_novo"]["max_alt_parent"]:
        return False, "ALT_IN_PARENT"
    return True, "PASS"

def rule_ar_homo(row: Dict[str, str]) -> Tuple[bool, str]:
    rp, ap = parse_ad(row.get("AD_proband", "."))
    ab_p = allele_balance(rp, ap)

    if ap is None or ap < THR["ar_homo"]["min_alt_p"]:
        return False, "ALT_PROBAND_TOO_LOW"
    if ab_p is None or ab_p < THR["ar_homo"]["min_ab_p"]:
        return False, "AB_PROBAND_TOO_LOW_FOR_HOM"
    return True, "PASS"

def rule_ar_het(row: Dict[str, str]) -> Tuple[bool, str]:
    rp, ap = parse_ad(row.get("AD_proband", "."))
    ab_p = allele_balance(rp, ap)

    if ap is None or ap < THR["ar_het"]["min_alt_p"]:
        return False, "ALT_PROBAND_TOO_LOW"
    if ab_p is None:
        return False, "AB_PROBAND_MISSING"
    if not (THR["ar_het"]["min_ab_p"] <= ab_p <= THR["ar_het"]["max_ab_p"]):
        return False, "AB_PROBAND_OUT_OF_RANGE"
    return True, "PASS"

def rule_x_linked(row: Dict[str, str]) -> Tuple[bool, str]:
    rm, am = parse_ad(row.get("AD_mother", "."))
    rp, ap = parse_ad(row.get("AD_proband", "."))
    ab_p = allele_balance(rp, ap)

    if ap is None or ap < THR["x_linked"]["min_alt_p"]:
        return False, "ALT_PROBAND_TOO_LOW"
    if ab_p is None or ab_p < THR["x_linked"]["min_ab_p"]:
        return False, "AB_PROBAND_TOO_LOW"
    if am is None or am < THR["x_linked"]["min_alt_m"]:
        return False, "ALT_MOTHER_TOO_LOW"
    return True, "PASS"

RULES: Dict[str, Callable[[Dict[str, str]], Tuple[bool, str]]] = {
    "de_novo": rule_de_novo,
    "ar_homo": rule_ar_homo,
    "ar_het": rule_ar_het,
    "x_linked": rule_x_linked,
}

# -----------------------------
# Main
# -----------------------------
for mode, filename in FILES.items():
    in_path = os.path.join(IN_DIR, filename)
    if not os.path.exists(in_path):
        continue

    pass_path = os.path.join(OUT_DIR, f"{mode}_evidence_pass.tsv")
    fail_path = os.path.join(OUT_DIR, f"{mode}_evidence_fail.tsv")

    with open(in_path, "r", encoding="utf-8") as fin, \
         open(pass_path, "w", encoding="utf-8", newline="") as fpass, \
         open(fail_path, "w", encoding="utf-8", newline="") as ffail:

        reader = csv.DictReader(fin, delimiter="\t")
        base_fields = reader.fieldnames or []

        # Add helpful computed columns (does not break downstream; they can ignore)
        extra_fields = [
            "AB_father", "AB_mother", "AB_proband",
            "ALT_father", "ALT_mother", "ALT_proband",
            "canonical_chrom",
            "fail_reason",
        ]
        out_fields = base_fields + [c for c in extra_fields if c not in base_fields]

        writer_pass = csv.DictWriter(fpass, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        writer_fail = csv.DictWriter(ffail, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        writer_pass.writeheader()
        writer_fail.writeheader()

        total = 0
        okc = 0
        noncanon = 0
        rule = RULES[mode]

        for row in reader:
            total += 1

            chrom = row.get("CHROM", "")
            canonical = is_canonical_chrom(chrom)
            row["canonical_chrom"] = "YES" if canonical else "NO"
            if not canonical:
                noncanon += 1
                row["fail_reason"] = "NONCANONICAL_CONTIG"
                writer_fail.writerow(row)
                continue

            # compute AB/ALT for convenience
            rf, af = parse_ad(row.get("AD_father", "."))
            rm, am = parse_ad(row.get("AD_mother", "."))
            rp, ap = parse_ad(row.get("AD_proband", "."))

            row["ALT_father"] = "." if af is None else str(af)
            row["ALT_mother"] = "." if am is None else str(am)
            row["ALT_proband"] = "." if ap is None else str(ap)

            row["AB_father"] = fmt_float(allele_balance(rf, af))
            row["AB_mother"] = fmt_float(allele_balance(rm, am))
            row["AB_proband"] = fmt_float(allele_balance(rp, ap))

            ok, reason = rule(row)
            row["fail_reason"] = reason

            if ok:
                okc += 1
                writer_pass.writerow(row)
            else:
                writer_fail.writerow(row)

    print(f"[{mode}] input={total} pass={okc} fail={total-okc} noncanonical={noncanon} -> {pass_path}")

print("Step [10] Evidence filtering per inheritance model completed.")
