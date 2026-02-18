#!/usr/bin/env python3
import csv
import os
import logging
from typing import Optional

IN_DIR = "results/13_gnomad"
OUT_DIR = "results/13B_gnomad_filtered"
LOG_FILE = "logs/13B_gnomad_filter.log"

THRESHOLDS = {
    "de_novo": 1e-4,
    "x_linked": 1e-4,
    "ar_homo": 1e-3,
    "ar_het": 1e-3,
    "comphet": 1e-3,
}

# -----------------------
# Logging setup
# -----------------------
os.makedirs("logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE, encoding="utf-8"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger("step13B")

# -----------------------
# Helpers
# -----------------------

def parse_float(x: str) -> Optional[float]:
    if not x or x == ".":
        return None
    try:
        return float(x)
    except:
        return None

def max_af(row) -> Optional[float]:
    vals = [
        parse_float(row.get("gnomAD_exome_af")),
        parse_float(row.get("gnomAD_genome_af")),
        parse_float(row.get("gnomAD_faf95_popmax")),
    ]
    vals = [v for v in vals if v is not None]
    return max(vals) if vals else None

def model_from_filename(fn: str) -> str:
    for m in THRESHOLDS.keys():
        if fn.startswith(m + "_"):
            return m
    return "unknown"

# -----------------------
# Main
# -----------------------

os.makedirs(OUT_DIR, exist_ok=True)

logger.info("=== Step 13B: gnomAD rare filtering ===")
logger.info(f"Input dir : {IN_DIR}")
logger.info(f"Output dir: {OUT_DIR}")

for fn in os.listdir(IN_DIR):
    if not fn.endswith("_annotated_gnomad.tsv"):
        continue

    model = model_from_filename(fn)
    thr = THRESHOLDS.get(model, 1e-3)

    in_path = os.path.join(IN_DIR, fn)
    out_path = os.path.join(
        OUT_DIR,
        fn.replace("_annotated_gnomad.tsv", f"_rare_maxaf{thr}.tsv")
    )

    total = kept = api_errors = 0

    with open(in_path, newline="", encoding="utf-8") as fin, \
         open(out_path, "w", newline="", encoding="utf-8") as fout:

        r = csv.DictReader(fin, delimiter="\t")
        fields = r.fieldnames + ["gnomAD_max_af", "gnomAD_filter"]

        w = csv.DictWriter(fout, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        for row in r:
            total += 1
            status = row.get("gnomAD_status", "")
            maf = max_af(row)

            row["gnomAD_max_af"] = "." if maf is None else f"{maf:.6g}"

            if status == "OK":
                keep = (maf is None) or (maf <= thr)
                row["gnomAD_filter"] = "PASS" if keep else "FAIL"

            elif status == "NOT_FOUND":
                keep = True
                row["gnomAD_filter"] = "PASS"

            else:
                keep = True
                row["gnomAD_filter"] = "API_ERROR"
                api_errors += 1

            if keep:
                kept += 1
                w.writerow(row)

    logger.info(
        f"[{model}] {fn}: kept {kept}/{total} "
        f"(thr={thr}, api_errors={api_errors}) -> {out_path}"
    )

logger.info("Step 13B complete.")
logger.info(f"Log written to: {LOG_FILE}")
