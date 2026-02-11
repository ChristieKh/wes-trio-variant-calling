#!/usr/bin/env python3
import csv
import os

IN_TSV = "results/final_candidates/final_candidates.tsv"
OUT_DIR = "results/final_candidates"
os.makedirs(OUT_DIR, exist_ok=True)

OUTS = {
    "recessive": os.path.join(OUT_DIR, "shortlist_recessive.tsv"),
}


LOW_VALUE_EFFECTS = (
    "synonymous_variant",
    "intron_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "intergenic_variant",
    "non_coding_transcript_variant",
    "non_coding_transcript_exon_variant",
    "5_prime_utr_variant",
    "3_prime_utr_variant",
)

def contains_any(text: str, keys) -> bool:
    t = (text or "").lower()
    return any(k.lower() in t for k in keys)

def main():
    with open(IN_TSV, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        rows = list(r)
        header = r.fieldnames or []

    def write(path, filt, limit=None):
        out_rows = [x for x in rows if filt(x)]
        if limit is not None:
            out_rows = out_rows[:limit]
        with open(path, "w", newline="", encoding="utf-8") as out:
            w = csv.DictWriter(out, fieldnames=header, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(out_rows)
        print(f"Wrote {len(out_rows)} rows -> {path}")


    def recessive_filter(x):
        if x.get("inheritance") != "recessive":
            return False
        eff = x.get("ANN[0].EFFECT", "")
        if eff and contains_any(eff, LOW_VALUE_EFFECTS):
            return False
        return True

    write(OUTS["recessive"], recessive_filter, limit=None)

if __name__ == "__main__":
    main()
