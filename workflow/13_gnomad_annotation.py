#!/usr/bin/env python3
import argparse
import csv
import json
import os
import time
import logging
import requests
from typing import Dict, Optional, Tuple

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
DEFAULT_LOG_FILE = "logs/13_gnomad.log"


def setup_logger(log_file: str) -> logging.Logger:
    os.makedirs(os.path.dirname(log_file), exist_ok=True)
    logger = logging.getLogger("step13_gnomad")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setFormatter(fmt)
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)

    logger.addHandler(fh)
    logger.addHandler(sh)
    return logger


def norm_chrom(chrom: str) -> str:
    chrom = (chrom or "").strip()
    if chrom.startswith("chr"):
        chrom = chrom[3:]
    if chrom == "M":
        chrom = "MT"
    return chrom


def is_canonical_chrom(chrom: str) -> bool:
    c = norm_chrom(chrom)
    if c in {"X", "Y", "MT"}:
        return True
    if c.isdigit():
        n = int(c)
        return 1 <= n <= 22
    return False


def make_variant_id(row: Dict[str, str]) -> str:
    return f"{norm_chrom(row['CHROM'])}-{row['POS']}-{row['REF']}-{row['ALT']}"


def ac_an_to_af(ac, an) -> Optional[float]:
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


def fetch_variant_once(session: requests.Session, variant_id: str, dataset: str) -> dict:
    payload = {
        "query": QUERY,
        "variables": {"variantId": variant_id, "datasetId": dataset},
    }
    r = session.post(API_URL, json=payload, timeout=30)
    r.raise_for_status()
    j = r.json()

    if "errors" in j:
        msg = ""
        try:
            msg = j["errors"][0].get("message", "")
        except Exception:
            msg = ""
        return {"status": "GRAPHQL_ERROR", "error": msg, "data": {}}

    variant = j.get("data", {}).get("variant")
    if variant is None:
        return {"status": "NOT_FOUND", "error": "", "data": {}}

    return {"status": "OK", "error": "", "data": variant}


def fetch_variant_with_retry(
    session: requests.Session,
    variant_id: str,
    dataset: str,
    retries: int,
    backoff: float,
    logger: logging.Logger,
) -> dict:
    """
    Retries only for transient failures (HTTP_ERROR/ERROR/GRAPHQL_ERROR).
    NOT_FOUND/OK are terminal.
    """
    last: dict = {"status": "ERROR", "error": "uninitialized", "data": {}}

    for attempt in range(retries + 1):
        try:
            res = fetch_variant_once(session, variant_id, dataset)
            status = res.get("status", "ERROR")

            if status in {"OK", "NOT_FOUND"}:
                return res

            # GRAPHQL_ERROR could be transient; retry
            last = res
            if attempt < retries:
                sleep_s = backoff * (2 ** attempt)
                logger.warning(f"{variant_id}: {status} (retry {attempt+1}/{retries}, sleep {sleep_s:.2f}s) {res.get('error','')}".strip())
                time.sleep(sleep_s)
                continue
            return res

        except requests.exceptions.HTTPError as e:
            last = {"status": "HTTP_ERROR", "error": str(e), "data": {}}
        except Exception as e:
            last = {"status": "ERROR", "error": str(e), "data": {}}

        if attempt < retries:
            sleep_s = backoff * (2 ** attempt)
            logger.warning(f"{variant_id}: {last['status']} (retry {attempt+1}/{retries}, sleep {sleep_s:.2f}s) {last.get('error','')}".strip())
            time.sleep(sleep_s)

    return last


def main():
    ap = argparse.ArgumentParser(
        description="Step 13: Annotate TSV variants with gnomAD (GraphQL API). Requires CHROM POS REF ALT columns."
    )
    ap.add_argument("--in", dest="in_tsv", required=True, help="Input TSV (must contain CHROM POS REF ALT)")
    ap.add_argument("--out", dest="out_tsv", required=True, help="Output TSV path")
    ap.add_argument("--dataset", default=DEFAULT_DATASET, help=f"gnomAD datasetId (default: {DEFAULT_DATASET})")
    ap.add_argument("--cache", default=DEFAULT_CACHE_FILE, help=f"Cache file jsonl (default: {DEFAULT_CACHE_FILE})")
    ap.add_argument("--sleep", type=float, default=0.6, help="Sleep seconds between successful API calls (default: 0.6)")
    ap.add_argument("--retries", type=int, default=3, help="Retries for transient errors (default: 3)")
    ap.add_argument("--backoff", type=float, default=0.5, help="Backoff base seconds for retries (default: 0.5)")
    ap.add_argument("--cache-errors", action="store_true", help="If set, cache ERROR/HTTP/GRAPHQL statuses (default: OFF)")
    ap.add_argument("--skip-noncanonical", action="store_true", help="If set, skip non-canonical contigs with status SKIPPED_NONCANONICAL")
    ap.add_argument("--max-rows", type=int, default=0, help="If >0, process only first N rows (for quick tests)")
    ap.add_argument("--log", default=DEFAULT_LOG_FILE, help=f"Log file path (default: {DEFAULT_LOG_FILE})")
    args = ap.parse_args()

    logger = setup_logger(args.log)

    os.makedirs(os.path.dirname(args.out_tsv), exist_ok=True)

    cache = load_cache(args.cache)
    session = requests.Session()

    logger.info(f"Input : {args.in_tsv}")
    logger.info(f"Output: {args.out_tsv}")
    logger.info(f"Cache : {args.cache} (loaded {len(cache)} entries)")
    logger.info(f"Dataset: {args.dataset}")
    logger.info(f"Sleep between API calls: {args.sleep}s")
    logger.info(f"Retries: {args.retries}, backoff: {args.backoff}s, cache_errors={args.cache_errors}")
    logger.info(f"skip_noncanonical={args.skip_noncanonical}, max_rows={args.max_rows if args.max_rows>0 else 'ALL'}")

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
            "gnomAD_error",
        ]
        writer = csv.DictWriter(fout, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        n = 0
        n_api = 0
        n_cache = 0
        n_skipped = 0

        for row in reader:
            n += 1
            if args.max_rows > 0 and n > args.max_rows:
                break

            if args.skip_noncanonical and not is_canonical_chrom(row.get("CHROM", "")):
                row["gnomAD_exome_af"] = "."
                row["gnomAD_genome_af"] = "."
                row["gnomAD_faf95_popmax"] = "."
                row["gnomAD_faf95_popmax_population"] = "."
                row["gnomAD_status"] = "SKIPPED_NONCANONICAL"
                row["gnomAD_error"] = ""
                writer.writerow(row)
                n_skipped += 1
                continue

            vid = make_variant_id(row)

            if vid in cache:
                result = cache[vid]
                n_cache += 1
            else:
                result = fetch_variant_with_retry(
                    session=session,
                    variant_id=vid,
                    dataset=args.dataset,
                    retries=args.retries,
                    backoff=args.backoff,
                    logger=logger,
                )
                result["variant_id"] = vid

                status = result.get("status", "ERROR")
                if status in {"OK", "NOT_FOUND"} or args.cache_errors:
                    cache[vid] = result
                    cache_put(args.cache, result)

                n_api += 1
                time.sleep(args.sleep)

            status = result.get("status", "ERROR")
            err = result.get("error", "") or ""
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
            row["gnomAD_error"] = err

            writer.writerow(row)

    logger.info(f"Done. Rows: {n if args.max_rows<=0 else min(n, args.max_rows)} | API calls: {n_api} | cache hits: {n_cache} | skipped_noncanonical: {n_skipped}")
    logger.info(f"Wrote: {args.out_tsv}")
    logger.info(f"Log: {args.log}")


if __name__ == "__main__":
    main()
