# Load required packages
library(tidyverse)

# Load logCPM matrix
logCPM <- read.csv("data/processed/TCGA_LUAD_logCPM.csv", row.names = 1, check.names = FALSE)

# Prepare for CIBERSORTx
# Reset rownames as a column named "Gene"
logCPM_cibersortx <- logCPM %>%
  rownames_to_column(var = "Gene")

# Save as TSV (preferred by CIBERSORTx)
write.table(logCPM_cibersortx,
            "data/processed/TCGA_LUAD_logCPM_CIBERSORTx.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("✅ CIBERSORTx input file created: TCGA_LUAD_logCPM_CIBERSORTx.tsv\n")
