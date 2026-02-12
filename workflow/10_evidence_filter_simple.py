#!/usr/bin/env python3
import csv
import os

IN_DIR = "results/09_candidates"
OUT_DIR = "results/10_candidates_evidence"

FILES = {
    "de_novo": "de_novo.tsv",
    "ar_homo": "ar_homo.tsv",
    "ar_het": "ar_het.tsv",
    "x_linked": "x_linked.tsv",
}

os.makedirs(OUT_DIR, exist_ok=True)

def parse_ad(ad_str):
    if not ad_str or ad_str == ".":
        return None, None
    parts = ad_str.split(",")
    if len(parts) < 2:
        return None, None
    ref, alt = parts[0], parts[1]
    return int(float(ref)), int(float(alt))

def allele_balance(ref, alt):
    if ref is None or alt is None:
        return None
    total = ref + alt
    if total <= 0:
        return None
    return alt / total

def gt_is_ref(gt: str) -> bool:
    return gt in ("0/0", "0|0")

def gt_has_alt(gt: str) -> bool:
    return gt not in ("0/0", "0|0", ".", "./.", ".|.")

# -----------------------------
# Evidence rules per model
# -----------------------------

def pass_de_novo(row):
    rf, af = parse_ad(row["AD_father"])
    rm, am = parse_ad(row["AD_mother"])
    rp, ap = parse_ad(row["AD_proband"])
    ab_p = allele_balance(rp, ap)

    if ap is None or ap < 5:
        return False
    if ab_p is None or ab_p < 0.25:
        return False
    # parents should have zero ALT reads
    if af is None or am is None:
        return False
    if af > 0 or am > 0:
        return False
    return True

def pass_ar_homo(row):
    rp, ap = parse_ad(row["AD_proband"])
    ab_p = allele_balance(rp, ap)

    if ap is None or ap < 8:
        return False
    if ab_p is None or ab_p < 0.85:
        return False
    return True

def pass_ar_het(row):
    rp, ap = parse_ad(row["AD_proband"])
    ab_p = allele_balance(rp, ap)

    if ap is None or ap < 3:
        return False
    if ab_p is None or not (0.20 <= ab_p <= 0.80):
        return False
    return True

def pass_x_linked(row):
    # Male proband on chrX: caller may encode as 0/1 or 1/1; we just need solid ALT support.
    rm, am = parse_ad(row["AD_mother"])
    rp, ap = parse_ad(row["AD_proband"])
    ab_p = allele_balance(rp, ap)

    if ap is None or ap < 5:
        return False
    if ab_p is None or ab_p < 0.25:
        return False
    # mother should have evidence of being a carrier
    if am is None or am < 2:
        return False
    return True

RULES = {
    "de_novo": pass_de_novo,
    "ar_homo": pass_ar_homo,
    "ar_het": pass_ar_het,
    "x_linked": pass_x_linked,
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
        writer_pass = csv.DictWriter(fpass, fieldnames=reader.fieldnames, delimiter="\t", extrasaction="ignore")
        writer_fail = csv.DictWriter(ffail, fieldnames=reader.fieldnames, delimiter="\t", extrasaction="ignore")
        writer_pass.writeheader()
        writer_fail.writeheader()

        total = okc = 0
        rule = RULES.get(mode)

        for row in reader:
            total += 1
            ok = rule(row) if rule else False
            if ok:
                okc += 1
                writer_pass.writerow(row)
            else:
                writer_fail.writerow(row)

    print(f"[{mode}] input={total} pass={okc} fail={total-okc} -> {pass_path}")

print("Step [10] Evidence filtering per inheritance model completed.")
