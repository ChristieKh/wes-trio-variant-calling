WES trio analysis (father, mother, affected proband) on GRCh38.
Variant calling and prioritization under Mendelian inheritance models.

## Quality Control (FastQC & MultiQC)

Raw paired-end WES FASTQ files for all trio members (father, mother, proband) were assessed using FastQC and summarized with MultiQC.

**Key observations:**
- Overall per-base sequence quality was high across all samples.
- Adapter contamination was not detected.
- GC content profiles were consistent between all trio members.
- A FAIL flag was observed for *Per Base Sequence Content*, which is expected for targeted exome sequencing due to capture bias.
- All samples showed identical sequence content patterns, indicating a protocol-level effect rather than sample-specific issues.

**Conclusion:**
No trimming or additional filtering was applied. Data quality was deemed sufficient for downstream alignment and variant calling.

