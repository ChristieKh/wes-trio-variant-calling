# Quality & Thresholds (SOP)

This document defines the quality acceptance criteria and filtering thresholds used in the WES trio analysis.
Thresholds are aligned with the implemented workflow scripts (Step 08/09/10/16/13e).

---

## 1. Sample-level QC (high-level)

Sample-level QC is assessed prior to variant interpretation to ensure sufficient sequencing and alignment quality.
Key outputs are summarized in `reports/qc_summary.tsv`.

*(Detailed metrics and cutoffs are documented in the QC report; this SOP focuses on variant-level decision thresholds.)*

---

## 2. Variant-level baseline QC (applied before model logic)

Baseline genotype QC criteria (Step 08/09):

- Genotype must be non-missing: `GT != "mis"`
- Minimum depth: `DP >= 10`
- Minimum genotype quality: `GQ >= 20`

Implementation:
- `workflow/08_qc_and_filter.sh` (proband-focused initial filtering)
- `workflow/09_trio_logic.sh` (`COMMON_QC` applied to all trio members)

---

## 3. Evidence-based thresholds (allele balance + ALT read support)

Evidence thresholds are applied to reduce false positives and enforce model-consistent read support.

Definitions:
- `alt_count` = number of ALT-supporting reads (from AD)
- `AB` (allele balance) = ALT / (REF + ALT), computed where AD is available

### 3.1 De novo (secondary)

Criteria (Step 10):
- Proband: `ALT reads >= 5`, `AB >= 0.25`
- Parents: `ALT reads == 0` (max ALT in either parent = 0)

Source: `workflow/10_evidence_filter_simple.py` (`THR["de_novo"]`)

### 3.2 Autosomal recessive homozygous (primary)

Criteria (Step 10):
- Proband: `ALT reads >= 8`, `AB >= 0.85`

Source: `workflow/10_evidence_filter_simple.py` (`THR["ar_homo"]`)

### 3.3 Autosomal recessive heterozygous / carrier-like (used for comp-het)

Criteria (Step 10):
- Proband (and/or relevant carrier genotypes): `ALT reads >= 3`
- `0.20 <= AB <= 0.80`

Source: `workflow/10_evidence_filter_simple.py` (`THR["ar_het"]`)

### 3.4 X-linked (secondary)

Criteria (Step 10):
- Proband: `ALT reads >= 5`, `AB >= 0.25`
- Mother: `ALT reads >= 2` (carrier support)

Source: `workflow/10_evidence_filter_simple.py` (`THR["x_linked"]`)

---

## 4. Compound Heterozygous Pairing QC (Step 16)

Additional QC constraints applied during comp-het pairing:

- Proband:
  - DP >= 10
  - GQ >= 20

- Parents:
  - DP >= 8 (if genotype present)

- Allele balance (AB):
  - 0.25 <= AB <= 0.75

Source: `workflow/16_comphet_pairing.py`

---

## 5. Population Frequency Thresholds (gnomAD)

Population frequency annotation is joined as `gnomAD_AF` (Step 13d) and filtered in Step 13e.

Model-specific allele frequency thresholds:

- De novo & X-linked: AF < 1e-4
- AR homozygous & heterozygous: AF < 1e-3

Source: `workflow/13e_filter_by_gnomad.py`

---

## 6. Notes

- Baseline DP/GQ filtering is enforced before inheritance logic.
- AR models are primary due to the reported autosomal recessive inheritance pattern.
- De novo and X-linked are screened for completeness and deprioritized unless strongly phenotype-concordant.