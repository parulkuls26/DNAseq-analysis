# =============================================================================
# STEP 8: VAF Filtering and TP53 Analysis (R Script)
# =============================================================================

library(tidyverse)

# Load ANNOVAR output
file_path <- "~/Downloads/tumour_final_annotated.hg38_multianno.txt"
variants <- read.delim(file_path, header = FALSE)

# Define column headings
headings <- c("id", "ref", "alt", "info", "format", "sample")
Annotated_variants <- setNames(variants[-1,], 
                               c(variants[1, 1:(ncol(variants)-length(headings))] %>% unlist(), headings))

# Extract allele counts from sample column
headings_ac <- str_split(Annotated_variants$format[1], ":") %>% unlist()
AlleleCounts <- str_split(Annotated_variants$sample, ":") %>% 
    do.call("rbind", .) %>% 
    as.data.frame() %>% 
    setNames(headings_ac)

# Clean VAF column
AlleleCounts <- mutate(AlleleCounts, FREQ = gsub("%", "", FREQ) %>% as.numeric())
Annotated_variants <- cbind(Annotated_variants, AlleleCounts)
colnames(Annotated_variants) <- tolower(colnames(Annotated_variants))

# Filter for VAF >= 10%
Annotated_variants_vaf10 <- subset(Annotated_variants, freq >= 10)
cat(paste("Total variants passing VAF >= 10% filter:", nrow(Annotated_variants_vaf10), "\n"))

# Save filtered variants
write_tsv(Annotated_variants_vaf10, "~/Downloads/tumour_final_annotated_vaf10.txt")

# Extract TP53 variants
TP53_variants <- Annotated_variants_vaf10 %>% 
  filter(grepl("TP53", gene.refgene, ignore.case = TRUE))

TP53_report <- TP53_variants %>% 
  select(chr, start, end, ref, alt, gene.refgene, exonicfunc.refgene, 
         aachange.refgene, avsnp150, cosmic92_coding, freq)

print(TP53_report)
write_csv(TP53_report, "~/Downloads/TP53_variants_report.csv")
