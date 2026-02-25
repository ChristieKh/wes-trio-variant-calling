# WES Trio Analysis — End-to-End Germline Variant Discovery Pipeline

This repository implements a reproducible **whole-exome sequencing (WES) trio workflow**
following **GATK Best Practices**, covering the full path from raw alignment to
candidate variant prioritization.

The project combines:

- technical pipeline engineering
- inheritance-based filtering (trio logic)
- population database integration
- clinical-oriented variant prioritization

🔗 Key reference inspiration:  
“Clinical phenotype and trio whole exome sequencing data from a patient with glycogen storage disease IV in Indonesia”  
https://pmc.ncbi.nlm.nih.gov/articles/PMC11748688

---

## 🔄 Pipeline overview

### Preprocessing (per sample)

- Read alignment (BWA-MEM)
- BAM sorting and indexing
- Read group assignment (Picard)
- Duplicate marking (Picard MarkDuplicates)
- Base Quality Score Recalibration (GATK BQSR)
- Variant calling in gVCF mode (GATK HaplotypeCaller)

### Joint genotyping (trio-level)

- gVCF aggregation
- Joint genotyping (GATK)
- Raw multi-sample VCF generation

### Variant filtering & annotation

- Quality-based filtering
- Functional annotation (SnpEff)
- Population frequency integration (gnomAD)
- Clinical database integration (ClinVar)

### Inheritance modeling

- De novo detection
- Autosomal recessive (homozygous)
- Compound heterozygous pairing
- X-linked logic

### Candidate prioritization

- Technical scoring
- Population frequency thresholds
- Phenotype-driven prioritization

---

## 📁 Repository structure

- `workflow/` — step-by-step pipeline scripts (01–20)
- `ref/` — reference genome and known-sites resources
- `results/` — generated outputs (not version-controlled)
- `logs/` — tool logs (not version-controlled)
- `samples.tsv` — trio sample definitions
- `environment.yml` — reproducible software environment

---

## 🧬 Interpretation focus

The project emphasizes:

- detection of **rare, potentially pathogenic variants**
- correct modeling of trio inheritance patterns
- integration of population and clinical evidence
- structured candidate shortlist generation

---

## ⚙️ Reproducibility

Each pipeline step is implemented as an individual script.
Steps can be executed sequentially from `workflow/`.

Generated results and logs are excluded from version control
and can be reproduced by running the pipeline.

