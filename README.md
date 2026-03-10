# WES Trio Variant Discovery (Case Study)

This project reconstructs a **real clinical whole-exome trio analysis ([PMC11748688](https://pmc.ncbi.nlm.nih.gov/articles/PMC11748688/))** to identify rare, high-confidence germline variants consistent with the reported phenotype using an inheritance-driven prioritization workflow applied to WES data from a proband and both parents.

## Data
- Study: [SRP490127](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP490127&o=acc_s%3Aa) (NCBI SRA)
- Samples: Trio (proband + parents)
- Reference: GRCh38
- Variant calling: GATK HaplotypeCaller → joint genotyping

## Project Highlights

- **696,572 → 20 variants** through inheritance-driven filtering
- Evaluation of **four inheritance models**: AR-homo, compound heterozygous, de novo, X-linked
- Identification of **GBE1** as the most plausible causal candidate for Glycogen Storage Disease type IV

## Workflow
The pipeline performs:
```
FASTQ
 ↓
Alignment (BWA)
 ↓
Variant Calling (GATK)
 ↓
Joint VCF
 ↓
Trio Inheritance Models
 ├─ AR-homo (main priority)
 ├─ Comp-het (main priority)
 ├─ De novo
 └─ X-linked
 ↓
Functional (SnpEff) + Population Filtering (gnomAD)
 ↓
Clinical annotation (ClinVar)
 ↓
Technical Prioritization 
 ↓
Phenotype Refinement
 ↓
Final shortlist 
```

Full workflow description:  
`docs/02_pipeline_map.md`

## Variant Filtering Funnel

| Stage | Variants |
|------|---------|
| Joint VCF | 696,572 |
| After QC | 164,170 |
| AR-homo rare | 19 |
| Comp-het rare | 80 |
| De novo rare | 1 |
| X-linked rare | 5 |
| Technical candidates | 67 |
| Final shortlist | 20 |

Detailed funnel:  
`reports/variant_funnel.md`

## Key Finding

A homozygous missense variant in **GBE1** was identified as the most plausible causal candidate under the autosomal recessive model.  
GBE1 is associated with **Glycogen Storage Disease type IV (Andersen disease)** and shows strong concordance with the reported phenotype.

Full interpretation:  
`docs/03_candidate_interpretation.md`

## Repository Structure
```
docs/
00_case_question.md
01_quality_and_thresholds.md
02_pipeline_map.md
03_candidate_interpretation.md

reports/
variant_funnel.tsv
variant_funnel.md

workflow/
analysis scripts

results/
intermediate and final outputs
```

## What this project demonstrates

- Trio-based inheritance analysis
- Variant prioritization using population databases
- Integration of functional and clinical annotations
- Candidate-level interpretation using phenotype evidence


