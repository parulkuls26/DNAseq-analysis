#!/bin/bash

# STEP 5: Base Quality Score Recalibration (GATK BQSR)


module load gatk
module load samtools

REF="/data/home/ha251948/PK_OMICS/Reference/Homo_sapiens.GRCh38.108.dna.chromosome.17.fa"
KNOWN_SITES="/data/home/ha251948/PK_OMICS/Reference/gatkResources/resources_broad_hg38_v0_1000G_omni2.5.hg38.noCHR.vcf"
ALIGNMENT_DIR="/data/home/ha251948/PK_OMICS/Alignment"

# --- Germline ---
# BaseRecalibrator
gatk --java-options "-Xmx4G" BaseRecalibrator \
    -I ${ALIGNMENT_DIR}/germline.marked.bam \
    -R ${REF} \
    --known-sites ${KNOWN_SITES} \
    -O ${ALIGNMENT_DIR}/germline.table

# ApplyBQSR
gatk --java-options "-Xmx4G" ApplyBQSR \
    -R ${REF} \
    -I ${ALIGNMENT_DIR}/germline.marked.bam \
    --bqsr-recal-file ${ALIGNMENT_DIR}/germline.table \
    -O ${ALIGNMENT_DIR}/germline.recal.bam

samtools index ${ALIGNMENT_DIR}/germline.recal.bam

# --- Tumour ---
# BaseRecalibrator
gatk --java-options "-Xmx4G" BaseRecalibrator \
    -I ${ALIGNMENT_DIR}/tumour.marked.bam \
    -R ${REF} \
    --known-sites ${KNOWN_SITES} \
    -O ${ALIGNMENT_DIR}/tumour.table

# ApplyBQSR
gatk --java-options "-Xmx4G" ApplyBQSR \
    -R ${REF} \
    -I ${ALIGNMENT_DIR}/tumour.marked.bam \
    --bqsr-recal-file ${ALIGNMENT_DIR}/tumour.table \
    -O ${ALIGNMENT_DIR}/tumour.recal.bam

samtools index ${ALIGNMENT_DIR}/tumour.recal.bam

# Generate flagstat reports
mkdir -p QC
samtools flagstat ${ALIGNMENT_DIR}/tumour.recal.bam > QC/tumour.flagstat.txt
samtools flagstat ${ALIGNMENT_DIR}/germline.recal.bam > QC/germline.flagstat.txt
