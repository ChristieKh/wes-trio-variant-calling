#!/usr/bin/env python3
import csv
import os

IN_TSV = "results/final_candidates/final_candidates.tsv"
OUT_DIR = "results/final_candidates"
os.makedirs(OUT_DIR, exist_ok=True)

OUTS = {
    "de_novo": os.path.join(OUT_DIR, "shortlist_de_novo.tsv"),
    "high":    os.path.join(OUT_DIR, "shortlist_high_impact.tsv"),
    "ptv":     os.path.join(OUT_DIR, "shortlist_ptv.tsv"),
}

PTV_EFFECTS = (
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
)

def contains_any(text: str, keys) -> bool:
    t = (text or "").lower()
    return any(k in t for k in keys)

def main():
    with open(IN_TSV, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        header = r.fieldnames or []

    def write(path, filt, limit=200):
        out_rows = [x for x in rows if filt(x)]
        out_rows = out_rows[:limit]
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=header, delimiter="\t")
            w.writeheader()
            w.writerows(out_rows)
        print(f"Wrote {len(out_rows)} rows -> {path}")

    write(OUTS["de_novo"],
          lambda x: x.get("inheritance") == "de_novo",
          limit=300)

    write(OUTS["high"],
          lambda x: (x.get("ANN[0].IMPACT", "").upper() == "HIGH"),
          limit=300)

    write(OUTS["ptv"],
          lambda x: contains_any(x.get("ANN[0].EFFECT", ""), PTV_EFFECTS),
          limit=300)

if __name__ == "__main__":
    main()
