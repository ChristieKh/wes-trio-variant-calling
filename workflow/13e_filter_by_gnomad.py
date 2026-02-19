#!/usr/bin/env python3
import csv
import os
import logging
from typing import Optional

IN_DIR = "results/13_gnomad/final_results"
OUT_DIR = "results/13B_gnomad_filtered"
LOG_FILE = "logs/13B_gnomad_filter.log"

THRESHOLDS = {
    "de_novo": 1e-4,
    "x_linked": 1e-4,
    "ar_homo": 1e-3,
    "ar_het": 1e-3,
    "comphet": 1e-3,
}

os.makedirs("logs", exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, encoding="utf-8"), logging.StreamHandler()],
)
logger = logging.getLogger("step13B")


def parse_float(x: str) -> Optional[float]:
    if x is None:
        return None
    x = str(x).strip()
    if x == "" or x == ".":
        return None
    try:
        return float(x)
    except ValueError:
        return None


def model_from_filename(fn: str) -> str:
    for m in THRESHOLDS.keys():
        if fn.startswith(m + "_"):
            return m
    return "unknown"


os.makedirs(OUT_DIR, exist_ok=True)

logger.info("=== Step 13B: gnomAD rare filtering (using gnomAD_AF) ===")
logger.info(f"Input dir : {IN_DIR}")
logger.info(f"Output dir: {OUT_DIR}")

for fn in os.listdir(IN_DIR):
    if not fn.endswith("_functional.gnomad.tsv"):
        continue

    model = model_from_filename(fn)
    thr = THRESHOLDS.get(model, 1e-3)

    in_path = os.path.join(IN_DIR, fn)
    out_path = os.path.join(
        OUT_DIR,
        fn.replace("_functional.gnomad.tsv", f"_rare_af{thr}.tsv"),
    )

    total = kept = af_missing = 0

    with open(in_path, newline="", encoding="utf-8") as fin, open(
        out_path, "w", newline="", encoding="utf-8"
    ) as fout:
        r = csv.DictReader(fin, delimiter="\t")
        if r.fieldnames is None:
            raise SystemExit(f"ERROR: no header in {in_path}")

        # sanity check
        if "gnomAD_AF" not in r.fieldnames:
            raise SystemExit(
                f"ERROR: {in_path} missing gnomAD_AF column. Found columns: {r.fieldnames}"
            )

        fields = list(r.fieldnames)
        if "gnomAD_max_af" not in fields:
            fields.append("gnomAD_max_af")
        if "gnomAD_filter" not in fields:
            fields.append("gnomAD_filter")

        w = csv.DictWriter(fout, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()

        for row in r:
            total += 1
            af = parse_float(row.get("gnomAD_AF"))
            if af is None:
                af_missing += 1
                row["gnomAD_max_af"] = "."
                row["gnomAD_filter"] = "PASS_NO_GNOMAD"
                keep = True
            else:
                row["gnomAD_max_af"] = f"{af:.6g}"
                keep = af <= thr
                row["gnomAD_filter"] = "PASS" if keep else "FAIL"

            if keep:
                kept += 1
                w.writerow(row)

    logger.info(
        f"[{model}] {fn}: kept {kept}/{total} (thr={thr}, no_gnomad={af_missing}) -> {out_path}"
    )

logger.info("Step 13B complete.")
logger.info(f"Log written to: {LOG_FILE}")
