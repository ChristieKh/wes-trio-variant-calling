#!/usr/bin/env python3
import csv
import sys

IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}

def pick_best_ann(ann_raw: str):
    """
    ann_raw: comma-separated list of ANN entries.
    ANN entry: Allele|Effect|Impact|Gene|GeneID|FeatureType|FeatureID|Biotype|Rank|HGVS.c|HGVS.p|cDNA.pos/cDNA.len|CDS.pos/CDS.len|AA.pos/AA.len|Distance|Errors
    We pick best by impact rank, then prefer protein_coding.
    """
    if not ann_raw or ann_raw == ".":
        return None

    best = None
    best_key = (-1, -1)  # (impact_rank, protein_coding_flag)
    for entry in ann_raw.split(","):
        parts = entry.split("|")
        # need at least Effect/Impact/Gene
        if len(parts) < 8:
            continue
        effect = parts[1].strip()
        impact = parts[2].strip()
        gene = parts[3].strip()
        geneid = parts[4].strip() if len(parts) > 4 else ""
        featureid = parts[6].strip() if len(parts) > 6 else ""
        biotype = parts[7].strip() if len(parts) > 7 else ""

        impact_rank = IMPACT_RANK.get(impact, -1)
        pc = 1 if biotype == "protein_coding" else 0
        key = (impact_rank, pc)

        if key > best_key and gene:
            best_key = key
            best = parts

    return best

def main(in_tsv: str, out_tsv: str):
    with open(in_tsv, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        fields = r.fieldnames or []
        if "ANN_RAW" not in fields:
            raise SystemExit("ERROR: missing ANN_RAW column")

        out_fields = []
        # Replace ANN_RAW with parsed ANN[0].* fields
        for c in fields:
            if c == "ANN_RAW":
                out_fields += [
                    "ANN[0].GENE","ANN[0].GENEID","ANN[0].FEATUREID",
                    "ANN[0].EFFECT","ANN[0].IMPACT",
                    "ANN[0].HGVS_C","ANN[0].HGVS_P",
                    "ANN[0].CDNA_POS","ANN[0].CDS_POS","ANN[0].AA_POS",
                ]
            else:
                out_fields.append(c)

        with open(out_tsv, "w", newline="", encoding="utf-8") as out:
            w = csv.DictWriter(out, fieldnames=out_fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()

            for row in r:
                best = pick_best_ann(row.get("ANN_RAW",""))
                # defaults
                gene=geneid=featureid=effect=impact=hgvsc=hgvsp=cdna=cds=aa=""

                if best:
                    # indices per snpEff ANN format
                    effect = best[1].strip() if len(best) > 1 else ""
                    impact = best[2].strip() if len(best) > 2 else ""
                    gene   = best[3].strip() if len(best) > 3 else ""
                    geneid = best[4].strip() if len(best) > 4 else ""
                    featureid = best[6].strip() if len(best) > 6 else ""
                    hgvsc = best[9].strip() if len(best) > 9 else ""
                    hgvsp = best[10].strip() if len(best) > 10 else ""
                    cdna  = best[11].strip() if len(best) > 11 else ""
                    cds   = best[12].strip() if len(best) > 12 else ""
                    aa    = best[13].strip() if len(best) > 13 else ""

                out_row = {}
                for c in fields:
                    if c == "ANN_RAW":
                        out_row.update({
                            "ANN[0].GENE": gene,
                            "ANN[0].GENEID": geneid,
                            "ANN[0].FEATUREID": featureid,
                            "ANN[0].EFFECT": effect,
                            "ANN[0].IMPACT": impact,
                            "ANN[0].HGVS_C": hgvsc,
                            "ANN[0].HGVS_P": hgvsp,
                            "ANN[0].CDNA_POS": cdna,
                            "ANN[0].CDS_POS": cds,
                            "ANN[0].AA_POS": aa,
                        })
                    else:
                        out_row[c] = row.get(c, "")
                w.writerow(out_row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("USAGE: split_snpeff_ann.py <in_raw.tsv> <out_annotated.tsv>")
    main(sys.argv[1], sys.argv[2])
