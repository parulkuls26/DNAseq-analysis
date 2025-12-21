#!/bin/bash

## create the working directory and transfer files
# STEP 2: Quality Control with FastQC and MultiQC

# Load required module
module load fastqc

# Create directories
mkdir -p QC
mkdir -p FASTQ

# Move FASTQ files to dedicated folder
mv DNAseq/*.fq.gz FASTQ/

# Verify files
ls -lh FASTQ/

# Run FastQC on all files
fastqc -o QC/ FASTQ/*.fq.gz

# Load Miniforge / Conda for MultiQC
module load miniforge

# Create conda environment (only needed once)
mamba create --name multiqc -c bioconda multiqc

# Activate and run MultiQC
mamba activate multiqc
multiqc QC/ -o QC/
mamba deactivate
