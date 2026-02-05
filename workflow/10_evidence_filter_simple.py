#!/usr/bin/env python3

import csv
import os

# -----------------------------
# Input / Output
# -----------------------------
IN_DIR = "results/candidates"
OUT_DIR = "results/candidates_evidence_simple"

FILES = {
    "de_novo": "de_novo.tsv",
    "inherited_het": "inherited_het.tsv",
    "recessive": "recessive.tsv",
}

os.makedirs(OUT_DIR, exist_ok=True)

# -----------------------------
# Helper functions
# -----------------------------

def parse_ad(ad_str):
    """
    Parse AD field of the form "REF,ALT"
    Returns (ref_reads, alt_reads)
    """
    if ad_str is None or ad_str == ".":
        return None, None
    ref, alt = ad_str.split(",")
    return int(float(ref)), int(float(alt))


def allele_balance(ref, alt):
    """
    AB = ALT / (REF + ALT)
    """
    if ref is None or alt is None:
        return None
    total = ref + alt
    if total == 0:
        return None
    return alt / total


# -----------------------------
# Evidence rules
# -----------------------------

def pass_de_novo(row):
    rf, af = parse_ad(row["AD_father"])
    rm, am = parse_ad(row["AD_mother"])
    rp, ap = parse_ad(row["AD_proband"])

    ab_p = allele_balance(rp, ap)

    if ap < 5:
        return False
    if ab_p is None or ab_p < 0.25:
        return False
    if af > 0 or am > 0:
        return False

    return True


def pass_inherited_het(row):
    rp, ap = parse_ad(row["AD_proband"])
    ab_p = allele_balance(rp, ap)

    if ap < 3:
        return False
    if ab_p is None or not (0.25 <= ab_p <= 0.75):
        return False

    return True


def pass_recessive(row):
    rp, ap = parse_ad(row["AD_proband"])
    ab_p = allele_balance(rp, ap)

    if ap < 8:
        return False
    if ab_p is None or ab_p < 0.85:
        return False

    return True


# -----------------------------
# Main processing
# -----------------------------

for mode, filename in FILES.items():
    in_path = os.path.join(IN_DIR, filename)
    pass_path = os.path.join(OUT_DIR, f"{mode}_evidence_pass.tsv")
    fail_path = os.path.join(OUT_DIR, f"{mode}_evidence_fail.tsv")

    if not os.path.exists(in_path):
        continue

    with open(in_path, "r") as fin, \
         open(pass_path, "w") as fpass, \
         open(fail_path, "w") as ffail:

        reader = csv.DictReader(fin, delimiter="\t")
        writer_pass = csv.DictWriter(fpass, fieldnames=reader.fieldnames, delimiter="\t")
        writer_fail = csv.DictWriter(ffail, fieldnames=reader.fieldnames, delimiter="\t")

        writer_pass.writeheader()
        writer_fail.writeheader()

        for row in reader:
            if mode == "de_novo":
                ok = pass_de_novo(row)
            elif mode == "inherited_het":
                ok = pass_inherited_het(row)
            elif mode == "recessive":
                ok = pass_recessive(row)
            else:
                ok = False

            if ok:
                writer_pass.writerow(row)
            else:
                writer_fail.writerow(row)

print("Step [10] (simple evidence filtering with pass/fail) completed.")
