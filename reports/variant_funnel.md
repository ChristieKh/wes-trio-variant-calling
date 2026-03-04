# Variant Filtering Funnel

At the VCF level, baseline genotype QC reduced the joint callset from **696,572** to **164,170** variants.
This step removes low-confidence genotypes by enforcing non-missing calls and minimum support (DP>=10, GQ>=20), providing a higher-quality starting point for inheritance-based prioritization.

Next stages quantify how many candidates remain after inheritance modeling (AR homozygous / comp-het), functional constraints, population frequency filtering (gnomAD), and clinical annotation (ClinVar).

## AR Homozygous Branch

Evidence-based filtering resulted in 4,989 candidate AR-homozygous variants.

Functional filtering reduced this to 473 coding/splice-impact variants.

Population frequency filtering (gnomAD AF < 0.001) further reduced the set to 19 rare candidates.

ClinVar annotation did not alter the candidate count but provided clinical context.

The major reduction steps were:
- Functional impact filtering (~90% reduction)
- Population rarity filtering (~96% reduction from functional set)

The final set of **19 AR-homozygous rare variants** represents a manageable and biologically plausible candidate space for interpretation.


## Compound Heterozygous Branch

Initial pairing identified 353 candidate gene-level pairs.

Strict segregation and evidence validation did not substantially reduce this set, indicating consistent genotype support.

Population frequency filtering (AF < 0.001) reduced the candidate set to 80 rare pairs.

Technical prioritization further narrowed the set to **42 high-confidence candidates.**

## De Novo Branch (Secondary Model)

Evidence-based filtering identified 121 candidate de novo variants.

Application of a stringent population frequency threshold (gnomAD AF < 1e-4) reduced this set to **1 ultra-rare candidate.**

The substantial reduction after frequency filtering indicates that most initial de novo calls represent background variation rather than plausible causal candidates.

Given the autosomal recessive inheritance pattern described in the case, the de novo model was retained for completeness but did not represent the primary explanatory pathway.


## X-Linked Branch (Secondary Model)

Initial X-linked screening yielded 388 candidate variants passing genotype and evidence thresholds.

After applying the ultra-rare frequency filter (gnomAD AF < 1e-4), **5 candidates** remained.

Although X-linked inheritance was evaluated to ensure comprehensive model coverage, the limited number of rare candidates and the documented autosomal recessive inheritance pattern suggest that this model is less likely to explain the phenotype.

X-linked variants were therefore deprioritized relative to autosomal recessive models.


## Secondary Models Summary

Screening of de novo and X-linked models was performed to ensure comprehensive inheritance evaluation.

Both secondary models produced substantially fewer rare candidates compared to autosomal recessive branches, supporting the initial hypothesis of recessive inheritance as the most plausible explanation in this case.



## Model Aggregation and Final Candidate Space

After applying inheritance-specific evidence thresholds and population frequency filtering across all evaluated models (autosomal recessive homozygous, compound heterozygous, de novo, and X-linked), a total of 67 rare, technically supported candidates remained.

This number represents the union of all inheritance-consistent variants prior to phenotype-based refinement.

Technical prioritization was applied to rank these 67 candidates based on quality metrics, rarity, and inheritance consistency. No additional filtering was performed at this stage; the candidate space remained constant.

Subsequent phenotype-guided refinement reduced the candidate set from 67 to 20 final variants. This final reduction reflects biological plausibility and gene–phenotype concordance rather than purely technical criteria.

The clear separation between model-driven filtering (**67 candidates**) and phenotype-driven prioritization (**20 candidates**) demonstrates a structured and interpretable analytical workflow.