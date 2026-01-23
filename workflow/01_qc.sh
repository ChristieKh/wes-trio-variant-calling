#!/bin/bash

fastqc -t 4 data/fastq/*.fastq.gz -o results/qc
multiqc results/qc -o results/qc
