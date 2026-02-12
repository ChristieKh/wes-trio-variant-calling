#!/usr/bin/env python3
import csv
import json
import os
import time
import requests
from typing import Dict

API_URL = "https://gnomad.broadinstitute.org/api"

QUERY = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    exome { ac an }
    genome { ac an }
    joint {
      faf95 { popmax popmax_population }
    }
  }
}
"""

INPUT_DIR = "results/12_model_tables"
OUTPUT_DIR = "results/13_gnomad"
CACHE_FILE = "results/13_gnomad/gnomad_cache.jsonl"
DATASET = "gnomad_r4"

os.makedirs(OUTPUT_DIR, exist_ok=True)

def norm_chrom(chrom: str) -> str:
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom

def make_variant_id(row):
    return f"{norm_chrom(row['CHROM'])}-{row['POS']}-{row['REF']}-{row['ALT']}"

def ac_an_to_af(ac, an):
    if not ac or not an or an == 0:
        return None
    return ac / an

# Load cache
cache: Dict[str, dict] = {}
if os.path.exists(CACHE_FILE):
    with open(CACHE_FILE) as f:
        for line in f:
            obj = json.loads(line.strip())
            cache[obj["variant_id"]] = obj

def cache_put(obj):
    with open(CACHE_FILE, "a") as f:
        f.write(json.dumps(obj) + "\n")

session = requests.Session()

def fetch_variant(variant_id):
    payload = {
        "query": QUERY,
        "variables": {
            "variantId": variant_id,
            "datasetId": DATASET
        }
    }

    try:
        r = session.post(API_URL, json=payload, timeout=30)
        r.raise_for_status()
        j = r.json()

        if "errors" in j:
            return {"status": "GRAPHQL_ERROR", "data": {}}

        variant = j.get("data", {}).get("variant")

        if variant is None:
            return {"status": "NOT_FOUND", "data": {}}

        return {"status": "OK", "data": variant}

    except requests.exceptions.HTTPError:
        return {"status": "HTTP_ERROR", "data": {}}
    except Exception:
        return {"status": "ERROR", "data": {}}

# Process each model table
for fname in os.listdir(INPUT_DIR):
    if not fname.endswith("_annotated.tsv"):
        continue

    in_path = os.path.join(INPUT_DIR, fname)
    out_path = os.path.join(OUTPUT_DIR, fname.replace("_annotated.tsv", "_annotated_gnomad.tsv"))

    print(f"\nProcessing {fname}")

    with open(in_path) as fin, open(out_path, "w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = reader.fieldnames + [
            "gnomAD_exome_af",
            "gnomAD_genome_af",
            "gnomAD_faf95_popmax",
            "gnomAD_faf95_popmax_population",
            "gnomAD_status"
        ]

        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            vid = make_variant_id(row)

            if vid in cache:
                result = cache[vid]
            else:
                result = fetch_variant(vid)
                result["variant_id"] = vid
                cache[vid] = result
                cache_put(result)
                time.sleep(0.6)

            status = result["status"]
            data = result.get("data", {})

            ex_af = "."
            ge_af = "."
            faf = "."
            faf_pop = "."

            if status == "OK":
                ex = data.get("exome") or {}
                ge = data.get("genome") or {}
                joint = data.get("joint") or {}

                ex_af_val = ac_an_to_af(ex.get("ac"), ex.get("an"))
                ge_af_val = ac_an_to_af(ge.get("ac"), ge.get("an"))

                if ex_af_val is not None:
                    ex_af = f"{ex_af_val:.6g}"
                if ge_af_val is not None:
                    ge_af = f"{ge_af_val:.6g}"

                faf_block = joint.get("faf95") or {}
                if faf_block:
                    faf = faf_block.get("popmax", ".")
                    faf_pop = faf_block.get("popmax_population", ".")

            row["gnomAD_exome_af"] = ex_af
            row["gnomAD_genome_af"] = ge_af
            row["gnomAD_faf95_popmax"] = faf
            row["gnomAD_faf95_popmax_population"] = faf_pop
            row["gnomAD_status"] = status

            writer.writerow(row)

    print(f"Written -> {out_path}")

print("\nStep 13 gnomAD complete.")
