#!/bin/bash

# STEP 7: Variant Annotation with ANNOVAR

mkdir -p Annovar

module load annovar/2019Oct24

# Convert VCF to ANNOVAR input format
convert2annovar.pl --format vcf4 \
    VCF/tumour.vcf \
    --includeinfo \
    --filter PASS \
    --outfile Annovar/tumour.pass.vcf

# Filter 1: Remove variants with >1% frequency in 1000 Genomes
annotate_variation.pl -filter \
    -dbtype 1000g2015aug_all \
    -buildver hg38 \
    -out VCF/tumour_filtered_temp \
    Annovar/tumour.pass.vcf \
    Reference/humandb/ \
    -maf 0.01

# Filter 2: Remove variants with >1% frequency in ExAC
annotate_variation.pl -filter \
    -dbtype exac03 \
    -buildver hg38 \
    -out VCF/tumour_doubly_filtered \
    VCF/tumour_filtered_temp.hg38_ALL.sites.2015_08_filtered \
    Reference/humandb/ \
    -score_threshold 0.01

# Annotate with gene names, dbSNP, COSMIC
table_annovar.pl \
    -buildver hg38 \
    -out Annovar/tumour_final_annotated \
    VCF/tumour_doubly_filtered.hg38_exac03_filtered \
    Reference/humandb/ \
    -remove \
    -otherinfo \
    -protocol refGene,avsnp150,cosmic92_coding,cytoband \
    -operation g,f,f,r \
    -nastring .

# Check output
head -n 4 Annovar/tumour_final_annotated.hg38_multianno.txt
