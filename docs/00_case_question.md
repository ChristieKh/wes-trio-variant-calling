# Case Origin

This analytical case is based on a published real-world WES trio study (PMCID: PMC11748688).
The phenotype and clinical context are derived from the reported case (see Table 1 in the publication).


This case study independently reconstructs the variant discovery and prioritization strategy using the reported inheritance context.

### 1. Objective

To identify rare, high-confidence germline variants in a WES trio (proband + parents) that are consistent with:

- Autosomal recessive inheritance
- The reported phenotype
- Mendelian inheritance models
- Sequencing quality constraints
- Population allele frequency expectations

The primary analytical goal is to produce a transparent, reproducible candidate prioritization workflow.

### 2. Data Specification

- Data type: Whole-Exome Sequencing (WES)
- Samples: Trio (proband, mother, father) from NCBI SRA study **SRP490127**.
- Variant types: SNVs and short indels
- Reference genome: GRCh38
- Variant calling strategy:
    - GATK HaplotypeCaller (gVCF mode)
    - Joint genotyping

####  Input FASTQ Status

According to the source publication, FASTQ files were already adapter-trimmed and cleaned prior to deposition.

No additional trimming was applied in this workflow.
Read-level quality was verified independently using FastQC prior to alignment

### 3. Analytical Hypothesis

Given the documented autosomal recessive inheritance pattern, the working hypothesis is:

> The causal variant is consistent with autosomal recessive inheritance.

Therefore, prioritization focuses on:

#### Tier 1 (Primary Models)

1. Autosomal recessive homozygous variants
    -  Homozygous alternative in proband
    -  Heterozygous in both parents

 2. Compound heterozygous variants
    -  Two heterozygous variants in the same gene
    -  Each inherited from a different parent

#### Tier 2 (Secondary Model)

3. De novo variants
    - Heterozygous in proband
    - Homozygous reference in both parents
    - Considered only if AR models do not yield convincing candidates

4. X-linked
    - Male proband: hemizygous alternative genotype
    - Mother: heterozygous carrier
    - Father: hemizygous reference
    - Screened for completeness but deprioritized unless strong phenotypic concordance is observed

###  4. Analytical Strategy

The analytical workflow follows a structured, inheritance-driven prioritization framework:

**1. Read-level validation**
- Verification of sequencing quality (no additional trimming applied as input FASTQ files were pre-cleaned).

**2. Alignment and variant calling**
- Alignment to GRCh38.
- Post-processing (duplicate marking, base quality recalibration).
- Variant calling in gVCF mode followed by joint genotyping.

**3. Variant-level technical filtering**
- Filtering based on depth (DP), genotype quality (GQ), allele balance, and Mendelian consistency.

**4. Functional annotation**
- Annotation of coding and splice-region variants.
- Retention of functionally relevant consequences.

**5. Population frequency filtering**
- Annotation with gnomAD.
- Filtering aligned with recessive inheritance expectations.

**6. Inheritance model classification**
- Primary focus on autosomal recessive models (homozygous and compound heterozygous).
- Secondary screening for de novo and X-linked variants.

**7. Prioritization and scoring**
- Integration of technical quality, population rarity, inheritance consistency, and gene-level aggregation.

**8. Phenotype-guided refinement**
- Alignment of candidate genes with phenotype-specific gene panels.
- Shortlist generation.

**9. Candidate-level interpretation**
- Manual review of top-ranked variants.
- Evaluation of biological plausibility.
- Assessment of gene–phenotype concordance.
- Cross-check with ClinVar and literature evidence.
- Selection of the most plausible candidate variant.

###  5. Analytical Boundaries

This case study does not include:
- Structural variant detection
- CNV analysis
- Deep intronic or regulatory variant analysis
- Experimental validation
- Clinical diagnosis or medical reporting

This is an analytical prioritization exercise.

###  6. Analytical Outputs

The analysis is expected to produce the following structured outputs:
- Sample-level sequencing quality assessment
- Documented variant filtering funnel with counts per stage
- Inheritance-classified candidate tables:
    - Autosomal recessive homozygous variants
    - Compound heterozygous variants
    - Secondary models (de novo / X-linked)
- Fully annotated and scored master variant table
- Phenotype-aligned ranked shortlist
- Candidate-level interpretation summaries