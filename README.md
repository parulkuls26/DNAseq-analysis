# Somatic Variant Calling Pipeline

A bioinformatics pipeline for identifying somatic mutations in matched tumour-normal paired-end sequencing data. This workflow processes raw FASTQ files through alignment, quality control, variant calling, and annotation to identify cancer-specific variants.

## Table of Contents

- [Overview](#overview)
- [Pipeline Structure](#pipeline-structure)
- [Input Data](#input-data)
- [Tools Used](#tools-used)
- [Pipeline Steps](#pipeline-steps)
  - [Preprocessing](#preprocessing)
  - [VarScan Pipeline](#varscan-pipeline-tumour-only)
  - [Mutect2 Pipeline](#mutect2-pipeline-somatic-paired)
- [Key Findings](#key-findings)
- [Output Files](#output-files)
- [Usage](#usage)
- [References](#references)

## Overview

This pipeline implements two complementary approaches for variant detection:

1. **VarScan (tumour-only)**: Identifies all variants in the tumour sample, then filters against population databases to remove common polymorphisms
2. **Mutect2 (somatic, paired)**: Directly compares tumour against matched germline to identify true somatic mutations

Both approaches identified the same clinically significant TP53 mutation (p.Y234S), providing cross-validation of results.

## Pipeline Structure
```
├── README.md
├── scripts/
│   ├── preprocessing/
│   │   ├── 02_qc.sh                 # Quality control (FastQC, MultiQC)
│   │   ├── 03_alignment.sh          # Bowtie2 alignment to chr17
│   │   ├── 04_mark_duplicates.sh    # GATK MarkDuplicates
│   │   └── 05_bqsr.sh               # Base quality score recalibration
│   ├── varscan/
│   │   ├── 06_varscan.sh            # VarScan variant calling
│   │   ├── 07_annotation.sh         # ANNOVAR annotation & filtering
│   │   └── 08_vaf_filtering.R       # VAF filtering & TP53 analysis
│   └── mutect2/
│       ├── 09_mutect2.sh            # Mutect2 somatic calling
│       └── 10_mutect2_analysis.R    # TP53 variant analysis
├── docs/
│   └── analysis_report.md           # Detailed methods and results
└── results/
    ├── varscan/
    └── mutect2/
```

## Input Data

| Sample | Files | Description |
|--------|-------|-------------|
| Tumour | `tumour_R1.fq.gz`, `tumour_R2.fq.gz` | Paired-end tumour sequencing |
| Germline | `germline_R1.fq.gz`, `germline_R2.fq.gz` | Matched normal control |

## Tools Used

| Tool | Version | Purpose |
|------|---------|---------|
| FastQC | - | Quality control |
| MultiQC | - | QC report aggregation |
| Bowtie2 | - | Read alignment |
| SAMtools | - | BAM manipulation |
| GATK | 4.x | Duplicate marking, BQSR, Mutect2 |
| VarScan | 2.4.3 | Variant calling |
| ANNOVAR | 2019Oct24 | Variant annotation |
| R/tidyverse | - | Data analysis and filtering |

## Pipeline Steps

### Preprocessing

These steps are shared by both variant calling approaches.

**Step 2: Quality Control**
- FastQC analysis on all FASTQ files
- MultiQC aggregation of QC metrics
- Verified Phred scores ~35, GC content 43-44%

**Step 3: Alignment**
- Aligned to human chromosome 17 (GRCh38) using Bowtie2
- Germline: 99.03% alignment rate, 15.2M read pairs
- Tumour: 98.97% alignment rate, 12.9M read pairs

**Step 4: Mark Duplicates**
- GATK MarkDuplicates to flag PCR duplicates
- Germline duplication rate: 4.95%
- Tumour duplication rate: 4.07%

**Step 5: Base Quality Score Recalibration**
- GATK BaseRecalibrator and ApplyBQSR
- Known sites: 1000 Genomes Omni 2.5 variants

### VarScan Pipeline (Tumour-only)

**Step 6: Variant Calling**
- Called SNPs and indels using VarScan mpileup2snp/mpileup2indel
- Parameters: min-coverage 20, min-avg-qual 20, min-var-freq 0.01
- Total variants: 94,258 (77,876 SNPs + 16,382 indels)

**Step 7: Annotation & Filtering**
- Filtered against 1000 Genomes (MAF >1%)
- Filtered against ExAC (frequency >1%)
- Annotated with refGene, dbSNP (avsnp150), COSMIC

**Step 8: VAF Filtering**
- Applied VAF ≥10% threshold
- Remaining variants: 12,470
- Extracted TP53 variants for analysis

### Mutect2 Pipeline (Somatic, Paired)

**Step 9: Somatic Variant Calling**
- GATK Mutect2 with germline as matched normal
- FilterMutectCalls to remove technical artifacts
- ANNOVAR annotation (refGene, dbSNP, COSMIC)

**Step 10: TP53 Analysis**
- Extracted TP53 variants with filter status and VAF
- Compared results with VarScan findings

## Key Findings

### TP53 p.Y234S Mutation

Both pipelines independently identified the same somatic mutation in TP53:

| Attribute | Value |
|-----------|-------|
| Position | chr17:7,674,262 |
| Reference | T |
| Alternate | G |
| cDNA change | c.A701C |
| Protein change | p.Y234S (Tyr → Ser) |
| Exon | 7 |
| dbSNP | rs587780073 |
| COSMIC | COSV52686167 |

**Variant Allele Frequency:**
| Method | VAF | Somatic Confirmation |
|--------|-----|----------------------|
| VarScan | 36% | Not assessed (tumour-only) |
| Mutect2 | 35.9% | Yes (PASS, germline AF 3.8%) |

The concordant VAF across both methods and Mutect2's somatic classification confirms this as a bona fide tumour-specific mutation in the TP53 tumour suppressor gene.

## Output Files

### VarScan Results
- `tumour.vcf` — Combined SNP and indel calls
- `tumour_final_annotated.hg38_multianno.txt` — Annotated variants
- `tumour_final_annotated_vaf10.txt` — VAF-filtered variants
- `TP53_variants_report.csv` — TP53 variant summary

### Mutect2 Results
- `somatic_filtered.vcf` — Filtered somatic calls
- `tumour_mutect2_somatic_final_annotated.hg38_multianno.txt` — Annotated somatic variants
- `TP53_Mutect2_all_variants_with_filters.csv` — TP53 variants with filter status

## Usage

### Prerequisites
- Access to HPC with required modules (FastQC, Bowtie2, GATK, SAMtools, ANNOVAR)
- R with tidyverse package
- Reference genome (GRCh38 chromosome 17)
- ANNOVAR databases (refGene, avsnp150, cosmic92_coding, 1000g2015aug_all, exac03)

### Running the Pipeline
```bash
# Preprocessing
bash scripts/preprocessing/02_qc.sh
bash scripts/preprocessing/03_alignment.sh
bash scripts/preprocessing/04_mark_duplicates.sh
bash scripts/preprocessing/05_bqsr.sh

# VarScan pipeline
bash scripts/varscan/06_varscan.sh
bash scripts/varscan/07_annotation.sh
Rscript scripts/varscan/08_vaf_filtering.R

# Mutect2 pipeline
bash scripts/mutect2/09_mutect2.sh
Rscript scripts/mutect2/10_mutect2_analysis.R
```

## Limitations

This pipeline is designed for SNPs and small indels. It does not detect:
- Large structural variants (SVs)
- Copy number variants (CNVs)
- Complex rearrangements

For these variant types, specialised tools such as CNVkit, Manta, or DELLY would be required.

## References

- [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)
- [VarScan Documentation](http://varscan.sourceforge.net/)
- [ANNOVAR Documentation](https://annovar.openbioinformatics.org/)
- [Bowtie2 Manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
