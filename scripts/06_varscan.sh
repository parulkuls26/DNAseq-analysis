#!/bin/bash

# STEP 6: Variant Calling with VarScan (tumour-only)


mkdir -p VCF

module load openjdk/17.0.8.1_1-gcc-12.2.0
module load samtools

REF="/data/home/ha251948/PK_OMICS/Reference/Homo_sapiens.GRCh38.108.dna.chromosome.17.fa"
TUMOUR_BAM="/data/home/ha251948/PK_OMICS/Alignment/tumour.recal.bam"
VARSCAN_JAR="/data/home/ha251948/PK_OMICS/VarScan.v2.4.3.jar"

# Call SNPs
samtools mpileup -q 20 -f ${REF} ${TUMOUR_BAM} | \
java -jar ${VARSCAN_JAR} mpileup2snp \
    --min-coverage 20 \
    --min-avg-qual 20 \
    --min-read2 4 \
    --p-value 0.2 \
    --min-var-freq 0.01 \
    --strand-filter 1 \
    --output-vcf 1 > VCF/tumour_snps_filtered.vcf

# Call Indels
samtools mpileup -q 20 -f ${REF} ${TUMOUR_BAM} | \
java -jar ${VARSCAN_JAR} mpileup2indel \
    --min-coverage 20 \
    --min-avg-qual 20 \
    --min-read2 4 \
    --p-value 0.2 \
    --min-var-freq 0.01 \
    --strand-filter 1 \
    --output-vcf 1 > VCF/tumour_indels_filtered.vcf

# Combine SNPs and Indels into one VCF
cat VCF/tumour_snps_filtered.vcf \
    <(grep -v '^#' VCF/tumour_indels_filtered.vcf) \
    > VCF/tumour.vcf

# Verify combined variant count
grep -v '^#' VCF/tumour.vcf | wc -l
# Expected: ~94,258 variants
