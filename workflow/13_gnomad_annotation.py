#!/usr/bin/env python3
import argparse
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

DEFAULT_DATASET = "gnomad_r4"
DEFAULT_CACHE_FILE = "results/13_gnomad/gnomad_cache.jsonl"


def norm_chrom(chrom: str) -> str:
    chrom = (chrom or "").strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom


def make_variant_id(row: Dict[str, str]) -> str:
    return f"{norm_chrom(row['CHROM'])}-{row['POS']}-{row['REF']}-{row['ALT']}"


def ac_an_to_af(ac, an):
    if ac is None or an is None or an == 0:
        return None
    return ac / an


def load_cache(cache_file: str) -> Dict[str, dict]:
    cache: Dict[str, dict] = {}
    if os.path.exists(cache_file):
        with open(cache_file, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                obj = json.loads(line)
                vid = obj.get("variant_id")
                if vid:
                    cache[vid] = obj
    return cache


def cache_put(cache_file: str, obj: dict):
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    with open(cache_file, "a", encoding="utf-8") as f:
        f.write(json.dumps(obj) + "\n")


def fetch_variant(session: requests.Session, variant_id: str, dataset: str) -> dict:
    payload = {
        "query": QUERY,
        "variables": {"variantId": variant_id, "datasetId": dataset},
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


def main():
    ap = argparse.ArgumentParser(
        description="Step 13: Annotate TSV variants with gnomAD (GraphQL API). Requires CHROM POS REF ALT columns."
    )
    ap.add_argument("--in", dest="in_tsv", required=True, help="Input TSV (must contain CHROM POS REF ALT)")
    ap.add_argument("--out", dest="out_tsv", required=True, help="Output TSV path")
    ap.add_argument("--dataset", default=DEFAULT_DATASET, help=f"gnomAD datasetId (default: {DEFAULT_DATASET})")
    ap.add_argument("--cache", default=DEFAULT_CACHE_FILE, help=f"Cache file jsonl (default: {DEFAULT_CACHE_FILE})")
    ap.add_argument("--sleep", type=float, default=0.6, help="Sleep seconds between API calls (default: 0.6)")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)

    cache = load_cache(args.cache)
    session = requests.Session()

    print(f"Input : {args.in_tsv}")
    print(f"Output: {args.out_tsv}")
    print(f"Cache : {args.cache} (loaded {len(cache)} entries)")
    print(f"Dataset: {args.dataset}")
    print("")

    with open(args.in_tsv, encoding="utf-8") as fin, open(args.out_tsv, "w", newline="", encoding="utf-8") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: input TSV has no header")

        required = {"CHROM", "POS", "REF", "ALT"}
        missing = sorted(list(required - set(reader.fieldnames)))
        if missing:
            raise SystemExit(f"ERROR: missing required columns: {missing}")

        fieldnames = reader.fieldnames + [
            "gnomAD_exome_af",
            "gnomAD_genome_af",
            "gnomAD_faf95_popmax",
            "gnomAD_faf95_popmax_population",
            "gnomAD_status",
        ]
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        n = 0
        n_api = 0

        for row in reader:
            n += 1
            vid = make_variant_id(row)

            if vid in cache:
                result = cache[vid]
            else:
                result = fetch_variant(session, vid, args.dataset)
                result["variant_id"] = vid
                cache[vid] = result
                cache_put(args.cache, result)
                n_api += 1
                time.sleep(args.sleep)

            status = result["status"]
            data = result.get("data", {}) or {}

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

    print(f"Done. Rows: {n}, API calls: {n_api} (rest from cache)")


if __name__ == "__main__":
    main()
