# Pipeline Map

This document summarizes the end-to-end analytical workflow of the WES trio project.
It provides a stepwise map of how raw sequencing data were transformed into ranked candidate variants and final phenotype-guided shortlists.

## Workflow Overview

| Step | Script | Main input | Main output | Purpose |
|------|--------|------------|-------------|---------|
| 01 | `01_qc.sh` | FASTQ files | FastQC reports | Verify read-level sequencing quality before alignment |
| 02 | `02_alignment.sh` | FASTQ, GRCh38 reference | aligned BAM | Align reads to the reference genome |
| 03 | `03_add_read_group.sh` | aligned BAM | BAM with read groups | Add read group information required for downstream GATK processing |
| 04 | `04_mark_duplicates.sh` | BAM with read groups | deduplicated BAM, duplication metrics | Mark PCR/optical duplicates and retain QC metrics |
| 05 | `05_bqsr.sh` | deduplicated BAM, known sites | recalibrated BAM | Recalibrate base quality scores |
| 06 | `06_haplotypecaller.sh` | recalibrated BAM | per-sample gVCF | Call germline variants in gVCF mode |
| 07 | `07_joint_genotyping.sh` | trio gVCFs | joint VCF | Perform joint genotyping across the trio |
| 08 | `08_qc_and_filter.sh` | joint VCF | baseline-filtered VCF, bcftools QC outputs | Apply baseline genotype QC (non-missing GT, DP, GQ) |
| 09 | `09_trio_logic.sh` | baseline-filtered VCF | inheritance-specific candidate sets | Classify variants by inheritance model (AR-homo, de novo, X-linked, heterozygous carriers) |
| 10 | `10_evidence_filter_simple.py` | inheritance candidate tables | evidence-pass candidate tables | Enforce model-specific ALT read and allele balance thresholds |
| 11 | `11_annotate_snpeff.sh` | evidence-pass variant sets | SnpEff-annotated VCFs | Add functional consequence annotations |
| 12 | `12_model_tables.sh` | annotated variant sets | model tables (`raw`, `annotated`, `functional`) | Convert annotated variants into structured tables and retain functional candidates |
| 12a | `12_split_snpeff_ann.py` | annotated tables | normalized annotation fields | Split SnpEff annotation into structured columns |
| 12b | `12b_filter_functional.py` | annotated tables | functional subsets | Retain coding and splice-relevant variants |
| 13 | `13a_build_master_variants.py` | model tables | master variant table | Build a unified analytical dataset across models |
| 13c | `13c_gnomad_ucsc_lookup_master.sh` | master / model variant tables | gnomAD lookup table | Retrieve population frequency annotations |
| 13d | `13d_join_gnomad_lookup.py` | variant tables + gnomAD lookup | gnomAD-annotated tables | Join gnomAD AF/AC/AN fields to candidate variants |
| 13e | `13e_filter_by_gnomad.py` | gnomAD-annotated tables | rare-variant subsets | Apply model-specific rarity thresholds |
| 14 | `14_clinvar.sh` | rare-variant subsets | ClinVar lookup files | Retrieve ClinVar annotations |
| 14b | `14B_join_clinvar.py` | rare-variant tables + ClinVar | ClinVar-annotated candidate tables | Add clinical database context |
| 15 | `15_prioritize_and_make_shortlists.py` | annotated candidate tables | ranked shortlists | Rank variants using technical and annotation-based criteria |
| 16 | `16_comphet_pairing.py` | heterozygous candidate variants | comp-het pair tables | Construct candidate compound heterozygous pairs |
| 17 | `17_comphet_pairs_merge.py` | comp-het pairs + annotations | merged comp-het candidate table | Merge evidence and annotations for compound heterozygous pairs |
| 18 | `18_comphet_shortlist.py` | merged comp-het pairs | comp-het priority tables / top shortlists | Prioritize compound heterozygous candidates |
| 19 | `19_technical_prioritication.py` | all rare model-specific candidates | `master_scored_all.tsv`, technical shortlist | Aggregate and rank technically supported candidates across models |
| 20 | `20_phenotype_prioritization.py` | technical shortlist, phenotype/gene panel information | phenotype-guided shortlist, final interpretation tables | Refine candidates based on phenotype concordance and produce final shortlist |

## Final Analytical Layers

The workflow can be viewed as four major analytical layers:

1. **Read and variant generation**
   - Steps 01–07  
   - FASTQ → BAM → gVCF → joint VCF

2. **Quality control and inheritance modeling**
   - Steps 08–10  
   - Baseline QC, trio logic, evidence filtering

3. **Functional and population-based prioritization**
   - Steps 11–18  
   - Functional annotation, model tables, gnomAD, ClinVar, comp-het construction

4. **Candidate ranking and biological refinement**
   - Steps 19–20  
   - Technical prioritization followed by phenotype-guided interpretation

## Key Final Outputs

- `results/19_master/master_scored_all.tsv`
- `results/20_interpretation/final_shortlist.tsv`
- `reports/variant_funnel.tsv`
- `reports/variant_funnel.md`