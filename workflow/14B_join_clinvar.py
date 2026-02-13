#!/usr/bin/env python3
import os
import sys
import csv

# args: IN_DIR LOOKUP_TSV OUT_DIR
in_dir, lookup_path, out_dir = sys.argv[1], sys.argv[2], sys.argv[3]

def norm_chrom(c: str) -> str:
    c = (c or "").strip()
    if c.startswith("chr"):
        c = c[3:]
    if c == "M":
        c = "MT"
    return c

# Load lookup keyed by normalized (chrom,pos,ref,alt)
lookup = {}
with open(lookup_path, newline="", encoding="utf-8") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        chrom = norm_chrom(row.get("CHROM", ""))
        key = (chrom, row.get("POS", ""), row.get("REF", ""), row.get("ALT", ""))
        lookup[key] = row

add_cols = [
    "ClinVar_CLNSIG",
    "ClinVar_CLNREVSTAT",
    "ClinVar_CLNDN",
    "ClinVar_CLNDISDB",
    "ClinVar_CLNVC",
    "ClinVar_ID",
]

def join_one(path_in: str, path_out: str):
    with open(path_in, newline="", encoding="utf-8") as fin, \
         open(path_out, "w", newline="", encoding="utf-8") as fout:
        r = csv.DictReader(fin, delimiter="\t")
        out_fields = (r.fieldnames or []) + add_cols
        w = csv.DictWriter(fout, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        total = hits = 0
        for row in r:
            total += 1
            chrom = norm_chrom(row.get("CHROM", ""))
            pos = row.get("POS", "")
            ref = row.get("REF", "")
            alt = row.get("ALT", "")
            key = (chrom, pos, ref, alt)
            cv = lookup.get(key)

            if cv:
                hits += 1
                row["ClinVar_CLNSIG"] = cv.get("CLNSIG", ".") or "."
                row["ClinVar_CLNREVSTAT"] = cv.get("CLNREVSTAT", ".") or "."
                row["ClinVar_CLNDN"] = cv.get("CLNDN", ".") or "."
                row["ClinVar_CLNDISDB"] = cv.get("CLNDISDB", ".") or "."
                row["ClinVar_CLNVC"] = cv.get("CLNVC", ".") or "."
                # ClinVar ID может быть в RS или в ID — мы в lookup кладём RS
                row["ClinVar_ID"] = cv.get("RS", ".") or "."
            else:
                row["ClinVar_CLNSIG"] = "."
                row["ClinVar_CLNREVSTAT"] = "."
                row["ClinVar_CLNDN"] = "."
                row["ClinVar_CLNDISDB"] = "."
                row["ClinVar_CLNVC"] = "."
                row["ClinVar_ID"] = "."

            w.writerow(row)

    print(f"{os.path.basename(path_in)} -> {os.path.basename(path_out)} | clinvar_hits={hits}/{total}")

os.makedirs(out_dir, exist_ok=True)

for fn in os.listdir(in_dir):
    if not fn.endswith(".tsv"):
        continue
    if "rare" not in fn:
        continue
    path_in = os.path.join(in_dir, fn)
    path_out = os.path.join(out_dir, fn.replace(".tsv", "_clinvar.tsv"))
    join_one(path_in, path_out)
