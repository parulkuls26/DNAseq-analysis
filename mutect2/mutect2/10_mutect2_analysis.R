
# STEP 10: Mutect2 Results Analysis - TP53 Variants (R Script)

library(tidyverse)

# Load Mutect2 annotated file
file_path <- "~/Downloads/tumour_mutect2_somatic_Dec5_final_annotated.hg38_multianno.txt"
mutect2_data <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)

# Filter for TP53 gene
TP53_somatic <- mutect2_data %>%
  filter(Gene.refGene == "TP53")

print(paste("Total TP53 variants found:", nrow(TP53_somatic)))

# Extract VAF and filter status
TP53_with_VAF <- TP53_somatic %>%
  mutate(
    Filter_Status = Otherinfo10,
    VAF = str_split(Otherinfo14, ":") %>%
          map_chr(~.x[3]) %>%
          as.numeric() * 100
  )

# Create report
TP53_report_all <- TP53_with_VAF %>%
  select(
    Chr, Start, End, Ref, Alt,
    Func.refGene, Gene.refGene, ExonicFunc.refGene, AAChange.refGene,
    avsnp150, cosmic92_coding,
    Filter_Status, VAF
  )

# Display results
print("All TP53 Variants with Filter Status and VAF:")
print(TP53_report_all)

print(paste("PASS variants:", sum(grepl("PASS", TP53_report_all$Filter_Status))))
print(paste("Filtered variants:", sum(!grepl("PASS", TP53_report_all$Filter_Status))))

# Save results
write_csv(TP53_report_all, "~/Downloads/TP53_Mutect2_all_variants_with_filters.csv")
```

