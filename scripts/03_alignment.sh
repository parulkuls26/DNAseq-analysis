#!/bin/bash

# STEP 3: Align reads to chromosome 17 using Bowtie2

mkdir -p Alignment

module unload miniforge
module load bowtie2
module load samtools

# Align tumour sample
time bowtie2 -p 8 \
    --rg ID:tumour \
    --rg SM:tumour \
    --rg PL:ILLUMINA \
    --rg LB:tumour \
    -x Reference/Bowtie2Idx/GRCh38.108.chr17 \
    -1 FASTQ/tumour_R1.fq.gz \
    -2 FASTQ/tumour_R2.fq.gz \
    2> Alignment/tumour_align.log | \
samtools sort -@ 8 -o Alignment/tumour.bam

samtools index Alignment/tumour.bam

# Align germline sample
time bowtie2 -p 8 \
    --rg ID:germline \
    --rg SM:germline \
    --rg PL:ILLUMINA \
    --rg LB:germline \
    -x Reference/Bowtie2Idx/GRCh38.108.chr17 \
    -1 FASTQ/germline_R1.fq.gz \
    -2 FASTQ/germline_R2.fq.gz \
    2> Alignment/germline_align.log | \
samtools sort -@ 8 -o Alignment/germline.bam

samtools index Alignment/germline.bam

# Verify BAM files
samtools view Alignment/germline.bam | head -n 3
samtools view Alignment/tumour.bam | head -n 3
