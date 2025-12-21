# Analysis Report

## 1. Quality Control

MultiQC was run to assess sequencing quality before proceeding with alignment.

### Key Findings

- **Phred scores**: Peak around 35 (~5M reads), indicating 99.97% base call accuracy
- **GC content**: 43-44% for both samples, within normal range for human genome
- **Duplication levels**: 25.5M (germline) and 30.2M (tumour) — expected for PCR-amplified libraries
- **Overall assessment**: All FastQC parameters passed (green), no significant quality issues

---

## 2. Alignment Statistics

| Metric | Germline | Tumour |
|--------|----------|--------|
| Read pairs | 15,233,202 | 12,904,869 |
| Overall alignment rate | 99.03% | 98.97% |
| Concordant unique mappings | 72.64% | 70.40% |
| Failed to align | 1.97% | 2.03% |

Both samples show excellent alignment rates (>98%) with high concordant mapping, suitable for downstream variant calling.

---

## 3. Duplicate Marking

| Metric | Germline | Tumour |
|--------|----------|--------|
| Read pairs examined | 12,719,219 | 15,029,782 |
| Duplicate pairs | 621,166 | 602,134 |
| Duplication rate | 4.95% | 4.07% |

Duplication rates <5% indicate excellent library preparation with minimal PCR over-amplification.

---

## 4. Base Quality Score Recalibration

**Known sites used**: 1000 Genomes Project Omni 2.5 variants (`resources_broad_hg38_v0_1000G_omni2.5.hg38.noCHR.vcf`)

While GATK best practices recommend multiple known-sites files (dbSNP, Mills indels, 1000G variants), only the 1000G Omni file was available in the teaching dataset. This resource contains high-confidence SNPs sufficient for BQSR on chromosome 17 data.

---

## 5. Variant Calling with VarScan

### Rationale

The input data (paired-end Illumina short reads) and chosen tools (VarScan → ANNOVAR) are designed for detecting SNPs and small indels — the most clinically relevant variants in cancer. This pipeline is not suited for:
- Large structural variants (SVs)
- Copy number variants (CNVs)
- Complex rearrangements

### Results

| Variant type | Count |
|--------------|-------|
| SNPs | 77,876 |
| Indels | 16,382 |
| **Total** | **94,258** |

---

## 6. Variant Annotation and Filtering

### Population Database Filtering

We used **ExAC** rather than ESP6500 for filtering common variants:
- ExAC contains ~60,000 exomes vs ESP6500's ~6,500 samples
- ExAC includes ESP6500 data plus additional large-scale projects
- ExAC is the current standard for exome variant filtering

### VAF Filtering

Variants with <10% variant allele frequency were removed.

**Result**: 12,470 variants passed the VAF ≥10% filter

---

## 7. TP53 Variants (VarScan)

### Variant 1: TP53 p.Y234S (Clinically Significant)

| Attribute | Value |
|-----------|-------|
| Chromosome | 17 |
| Position | 7,674,262 |
| Reference allele | T |
| Mutated allele | G |
| cDNA change | c.A701C (NM_000546) |
| Protein change | p.Y234S (Tyr → Ser) |
| Exon | 7 |
| VAF | 36% |
| dbSNP | rs587780073 |
| COSMIC | COSV52686167 |

This is a nonsynonymous missense mutation in the TP53 tumour suppressor gene, reported in ovary, skin, prostate, lung, and upper aerodigestive tract cancers.

### Variant 2: TP53I13 (Not TP53)

| Attribute | Value |
|-----------|-------|
| Position | 29,571,116 |
| Gene | TP53I13 (TP53 Inducible protein 13) |
| Location | Intronic |
| VAF | 46.15% |
| dbSNP | rs117301138 |

**Note**: This is NOT in the TP53 tumour suppressor gene itself.

---

## 8. Somatic Variant Calling with Mutect2

Mutect2 was run with the germline sample as matched normal control (`-normal germline`) to distinguish somatic mutations from germline variants. FilterMutectCalls was applied to remove technical artifacts.

### TP53 Variants Identified

| Position | Type | Change | VAF | Filter Status |
|----------|------|--------|-----|---------------|
| 7,674,262 | Exonic (missense) | p.Y234S | 35.9% | **PASS** |
| 3' UTR | UTR variant | — | 18.5% | base_qual;weak_evidence |
| 3' UTR | UTR variant | — | 18.5% | base_qual;weak_evidence |

The exonic p.Y234S mutation passed all quality filters and represents a high-confidence somatic mutation.

---

## 9. Concordance Between Methods

The TP53 p.Y234S mutation was independently identified by both pipelines:

| Method | VAF | Somatic Confirmation |
|--------|-----|----------------------|
| VarScan (tumour-only) | 36% | Not assessed |
| Mutect2 (paired) | 35.9% | Yes (PASS, germline AF 3.8%) |

The concordant VAF across two distinct algorithms and Mutect2's somatic classification (minimal germline presence) confirms this as a bona fide tumour-specific mutation.

---

## 10. Output Files

### Flagstat Reports

#### Germline (`germline.flagstat.txt`)
```
30438404 + 0 in total (QC-passed reads + QC-failed reads)
30438404 + 0 primary
0 + 0 secondary
0 + 0 supplementary
1242332 + 0 duplicates
1242332 + 0 primary duplicates
30150884 + 0 mapped (99.06% : N/A)
30150884 + 0 primary mapped (99.06% : N/A)
30438404 + 0 paired in sequencing
15219202 + 0 read1
15219202 + 0 read2
29470286 + 0 properly paired (96.82% : N/A)
30002056 + 0 with itself and mate mapped
148828 + 0 singletons (0.49% : N/A)
417038 + 0 with mate mapped to a different chr
294946 + 0 with mate mapped to a different chr (mapQ>=5)
```

#### Tumour (`tumour.flagstat.txt`)
```
25809738 + 0 in total (QC-passed reads + QC-failed reads)
25809738 + 0 primary
0 + 0 secondary
0 + 0 supplementary
1204268 + 0 duplicates
1204268 + 0 primary duplicates
25542254 + 0 mapped (98.96% : N/A)
25542254 + 0 primary mapped (98.96% : N/A)
25809738 + 0 paired in sequencing
12904869 + 0 read1
12904869 + 0 read2
24674024 + 0 properly paired (95.60% : N/A)
25381832 + 0 with itself and mate mapped
160422 + 0 singletons (0.62% : N/A)
557SEQ + 0 with mate mapped to a different chr
395862 + 0 with mate mapped to a different chr (mapQ>=5)
```

### VarScan Results

#### TP53 Variants Report (`TP53_variants_report.csv`)

| chr | start | end | ref | alt | gene.refgene | exonicfunc.refgene | aachange.refgene | avsnp150 | cosmic92_coding | freq |
|-----|-------|-----|-----|-----|--------------|-------------------|------------------|----------|-----------------|------|
| 17 | 7674262 | 7674262 | T | G | TP53 | nonsynonymous SNV | NM_000546:exon7:c.A701C:p.Y234S | rs587780073 | COSV52686167 | 36 |
| 17 | 29571116 | 29571116 | G | A | TP53I13 | . | . | rs117301138 | . | 46.15 |

### Mutect2 Results

#### TP53 Variants with Filter Status (`TP53_Mutect2_all_variants_with_filters.csv`)

| Chr | Start | End | Ref | Alt | Func.refGene | Gene.refGene | ExonicFunc.refGene | AAChange.refGene | avsnp150 | cosmic92_coding | Filter_Status | VAF |
|-----|-------|-----|-----|-----|--------------|--------------|-------------------|------------------|----------|-----------------|---------------|-----|
| 17 | 7674262 | 7674262 | T | G | exonic | TP53 | nonsynonymous SNV | NM_000546:exon7:c.A701C:p.Y234S | rs587780073 | COSV52686167 | PASS | 35.9 |
| 17 | 7673700 | 7673700 | C | T | UTR3 | TP53 | . | . | . | . | base_qual;weak_evidence | 18.5 |
| 17 | 7673767 | 7673767 | G | A | UTR3 | TP53 | . | . | . | . | base_qual;weak_evidence | 18.5 |

---

## 11. Conclusions

1. **High-quality data**: Both samples showed excellent sequencing quality (Phred ~35) and alignment rates (>98%)

2. **Low duplication**: Duplication rates <5% indicate good library preparation

3. **TP53 somatic mutation confirmed**: The p.Y234S mutation was identified by both VarScan and Mutect2 with concordant VAFs (~36%), and Mutect2 confirmed it as somatic (minimal presence in germline)

4. **Clinical relevance**: TP53 p.Y234S is a known cancer-associated mutation (COSMIC: COSV52686167) affecting the DNA-binding domain of the p53 tumour suppressor protein
```

---

