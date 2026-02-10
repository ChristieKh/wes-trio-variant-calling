#!/usr/bin/env python3
import csv
import json
import os
import time
import argparse
from typing import Dict, Optional, Tuple

import requests

API_URL = "https://gnomad.broadinstitute.org/api"

QUERY = """
query GnomadVariant($variantId: String!, $datasetId: DatasetId!) {
  variant(variantId: $variantId, dataset: $datasetId) {
    variantId
    exome { ac an }
    genome { ac an }
    joint {
      ac
      an
      fafmax {
        faf95_max
        faf95_max_gen_anc
      }
    }
  }
}
"""


def norm_chrom(chrom: str) -> str:
    chrom = chrom or ""
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom


def make_variant_id(chrom: str, pos: str, ref: str, alt: str) -> str:
    # gnomAD variant id format: "1-55051215-G-GA" (no "chr")
    return f"{norm_chrom(chrom)}-{pos}-{ref}-{alt}"


def ac_an_to_af(ac, an) -> Optional[float]:
    try:
        ac = int(ac) if ac is not None else None
        an = int(an) if an is not None else None
        if ac is None or an is None or an == 0:
            return None
        return ac / an
    except Exception:
        return None


def fetch_variant(
    session: requests.Session,
    variant_id: str,
    dataset_id: str,
    sleep_s: float,
    log_path: str,
) -> Tuple[dict, str]:
    """
    Returns: (data_dict, status_str)
    status_str: OK | NOT_FOUND | ERROR:HTTP_429 | ERROR:HTTP_403 | ERROR:HTTP_400 | ERROR:HTTP_5xx | ERROR:GRAPHQL | ERROR:REQUEST_EXCEPTION
    """
    time.sleep(sleep_s)

    payload = {"query": QUERY, "variables": {"variantId": variant_id, "datasetId": dataset_id}}

    try:
        r = session.post(API_URL, json=payload, timeout=30)

        if r.status_code >= 400:
            body = (r.text or "")[:300].replace("\n", " ")
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(f"{variant_id}\tHTTP_{r.status_code}\t{body}\n")
            if 500 <= r.status_code <= 599:
                return {}, "ERROR:HTTP_5xx"
            return {}, f"ERROR:HTTP_{r.status_code}"

        j = r.json()

        # GraphQL-level errors can come with HTTP 200
        if isinstance(j, dict) and j.get("errors"):
            msg = str(j["errors"])[:300].replace("\n", " ")
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(f"{variant_id}\tGRAPHQL_ERRORS\t{msg}\n")
            return j.get("data", {}) or {}, "ERROR:GRAPHQL"

        data = j.get("data", {}) or {}
        v = (data.get("variant") or {}) if isinstance(data, dict) else {}

        # "variant" can be null => not in gnomAD (not an HTTP error)
        if not v:
            return data, "NOT_FOUND"

        return data, "OK"

    except requests.exceptions.RequestException as e:
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(f"{variant_id}\tREQUEST_EXCEPTION\t{repr(e)}\n")
        return {}, "ERROR:REQUEST_EXCEPTION"
    except ValueError as e:
        # JSON decode errors
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(f"{variant_id}\tJSON_DECODE_ERROR\t{repr(e)}\n")
        return {}, "ERROR:JSON_DECODE"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="in_tsv", required=True, help="Input TSV (e.g. shortlist_ptv.tsv)")
    ap.add_argument("--out", dest="out_tsv", required=True, help="Output TSV with gnomAD columns")
    ap.add_argument("--dataset", default="gnomad_r4", help="DatasetId to use (try gnomad_r4, else gnomad_r3)")
    ap.add_argument("--sleep", type=float, default=1.2, help="Seconds to sleep between requests")
    ap.add_argument("--cache", default="results/13_external/gnomad_cache.jsonl", help="JSONL cache path")
    ap.add_argument(
        "--retry-errors",
        action="store_true",
        help="If set, will re-query variants cached with ERROR:* statuses.",
    )
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)
    os.makedirs(os.path.dirname(args.cache), exist_ok=True)

    log_path = os.path.join(os.path.dirname(args.out_tsv), "gnomad_api_errors.log")

    # Load cache (variant_id -> result dict)
    cache: Dict[str, dict] = {}
    if os.path.exists(args.cache):
        with open(args.cache, encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                obj = json.loads(line)
                if "variant_id" in obj:
                    cache[obj["variant_id"]] = obj

    def cache_put(variant_id: str, obj: dict):
        cache[variant_id] = obj
        with open(args.cache, "a", encoding="utf-8") as f:
            f.write(json.dumps(obj) + "\n")

    with open(args.in_tsv, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        header = r.fieldnames or []

    add_cols = [
        "gnomAD_variant_id",
        "gnomAD_dataset",
        "gnomAD_exome_af",
        "gnomAD_genome_af",
        "gnomAD_grpmax_faf95",
        "gnomAD_grpmax_faf95_gen_anc",
        "gnomAD_status",
    ]
    out_header = header + add_cols

    session = requests.Session()
    session.headers.update(
        {
            "Content-Type": "application/json",
            "Accept": "application/json",
            "User-Agent": "wes-trio-gnomad-af/1.0",
        }
    )

    hits = 0
    misses = 0

    with open(args.out_tsv, "w", newline="", encoding="utf-8") as out:
        w = csv.DictWriter(out, fieldnames=out_header, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        for row in rows:
            chrom = row.get("CHROM", ".")
            pos = row.get("POS", ".")
            ref = row.get("REF", ".")
            alt = row.get("ALT", ".")
            vid = make_variant_id(chrom, pos, ref, alt)

            row["gnomAD_variant_id"] = vid
            row["gnomAD_dataset"] = args.dataset

            use_cache = False
            cached_status = None
            cached_data = {}

            if vid in cache:
                cached_status = cache[vid].get("status", "CACHED")
                cached_data = cache[vid].get("data", {}) or {}
                # re-run only if requested and cached status is an ERROR
                if not (args.retry_errors and isinstance(cached_status, str) and cached_status.startswith("ERROR:")):
                    use_cache = True

            if use_cache:
                data = cached_data
                status = cached_status or "CACHED"
            else:
                data = {}
                status = "ERROR:UNKNOWN"

                # Retry with backoff only for transient-ish statuses
                for attempt in range(1, 6):
                    data, status = fetch_variant(
                        session=session,
                        variant_id=vid,
                        dataset_id=args.dataset,
                        sleep_s=args.sleep,
                        log_path=log_path,
                    )

                    if status in ("OK", "NOT_FOUND"):
                        break

                    # Only retry for throttling / server-side / request exceptions
                    if status in ("ERROR:HTTP_429", "ERROR:HTTP_5xx", "ERROR:REQUEST_EXCEPTION", "ERROR:JSON_DECODE"):
                        time.sleep(min(30, 2 ** (attempt - 1)))
                        continue

                    # For 400/403/etc - usually not worth retrying
                    break

                cache_put(vid, {"variant_id": vid, "status": status, "data": data})

            v = (data.get("variant") or {}) if isinstance(data, dict) else {}
            ex = v.get("exome") or {}
            ge = v.get("genome") or {}
            joint = v.get("joint") or {}
            fafmax = joint.get("fafmax") or {}

            ex_af = ac_an_to_af(ex.get("ac"), ex.get("an"))
            ge_af = ac_an_to_af(ge.get("ac"), ge.get("an"))

            row["gnomAD_exome_af"] = "." if ex_af is None else f"{ex_af:.6g}"
            row["gnomAD_genome_af"] = "." if ge_af is None else f"{ge_af:.6g}"
            row["gnomAD_joint_faf95_popmax"] = fafmax.get("faf95_max", ".") if fafmax else "."
            row["gnomAD_joint_faf95_popmax_population"] = fafmax.get("faf95_max_gen_anc", ".") if fafmax else "."
            row["gnomAD_status"] = status

            if v:
                hits += 1
            else:
                misses += 1

            w.writerow(row)

    print(f"Input rows: {len(rows)}")
    print(f"gnomAD variant hits: {hits}")
    print(f"gnomAD misses: {misses}")
    print(f"Output: {args.out_tsv}")
    print(f"Cache: {args.cache}")
    print(f"API error log: {log_path}")


if __name__ == "__main__":
    main()
