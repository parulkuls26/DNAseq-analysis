#!/bin/bash

# STEP 4: Mark duplicate reads using GATK

module load gatk

# Mark duplicates in germline
gatk MarkDuplicates \
    -I Alignment/germline.bam \
    -O Alignment/germline.marked.bam \
    -M Alignment/germline.metrics.txt \
    --CREATE_INDEX true

# Mark duplicates in tumour
gatk MarkDuplicates \
    -I Alignment/tumour.bam \
    -O Alignment/tumour.marked.bam \
    -M Alignment/tumour.metrics.txt \
    --CREATE_INDEX true

# Check duplication statistics (line 8 contains summary)
sed -n '8p' Alignment/germline.metrics.txt
sed -n '8p' Alignment/tumour.metrics.txt
