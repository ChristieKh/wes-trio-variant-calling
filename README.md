# WES Trio Analysis — Germline Variant Calling & Interpretation

This project implements a reproducible **whole-exome sequencing (WES) trio workflow**
following **GATK Best Practices**, with a focus on **germline variant discovery and interpretation**.

The analysis is inspired by published trio-WES studies and aims not only to build a technical pipeline,
but also to explore the **biological and clinical interpretation** of identified variants.

🔗 **Key reference paper:**  
“Clinical phenotype and trio whole exome sequencing data from a patient with glycogen storage disease IV in Indonesia”  
https://pmc.ncbi.nlm.nih.gov/articles/PMC11748688

## 🔄 Pipeline overview

1. [x] Read alignment to the reference genome (BWA MEM)
2. [x] BAM sorting and indexing
3. [x] Read group assignment (Picard)
4. [ ] PCR duplicate marking (Picard MarkDuplicates)
5. [ ] Base Quality Score Recalibration (BQSR)
6. [ ] Variant calling (GATK HaplotypeCaller, gVCF mode)
7. [ ] Joint genotyping of the trio
8. [ ] Variant filtering and interpretation

## 🧠 Interpretation focus

This project emphasizes:
- **de novo variants**
- recessive inheritance patterns
- genotype–phenotype reasoning
- comparison with diagnostic yields reported in trio-WES studies

Variant interpretation will be performed using public databases
and clinical genetics principles.

