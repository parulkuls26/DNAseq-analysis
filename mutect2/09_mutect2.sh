#!/bin/bash

# STEP 9: Somatic Variant Calling with Mutect2

# Mutect2 compares tumour against matched germline to identify somatic mutations
# This is different from VarScan tumour-only analysis

# Create output directory
mkdir -p Mutect2_VCF

# Load required modules
module load gatk
module load samtools

# Run Mutect2 (tumour vs germline)
gatk Mutect2 \
    -R Reference/Homo_sapiens.GRCh38.108.dna.chromosome.17.fa \
    -I Alignment/tumour.recal.bam \
    -I Alignment/germline.recal.bam \
    -normal germline \
    -O Mutect2_VCF/somatic_raw.vcf.gz


# Filter somatic variants

gatk FilterMutectCalls \
    -R Reference/Homo_sapiens.GRCh38.108.dna.chromosome.17.fa \
    -V Mutect2_VCF/somatic_raw.vcf.gz \
    -O Mutect2_VCF/somatic_filtered.vcf.gz

# Decompress filtered VCF
gunzip Mutect2_VCF/somatic_filtered.vcf.gz

# -----------------------------------------------------------------------------
# Annotate with ANNOVAR
# -----------------------------------------------------------------------------
module unload perl
module load annovar

# Convert VCF to ANNOVAR format
convert2annovar.pl --format vcf4 \
    Mutect2_VCF/somatic_filtered.vcf \
    --includeinfo \
    -allsample \
    -withfreq \
    --outfile Annovar/tumour_mutect2_somatic.pass.vcf

# Annotate somatic variants
table_annovar.pl \
    -buildver hg38 \
    -out Annovar/tumour_mutect2_somatic_final_annotated \
    Annovar/tumour_mutect2_somatic.pass.vcf \
    /data/teaching/bci_teaching/DNAseq/Reference/humandb/ \
    -remove \
    -otherinfo \
    -protocol refGene,avsnp150,cosmic92_coding,cytoband \
    -operation g,f,f,r \
    -nastring .
