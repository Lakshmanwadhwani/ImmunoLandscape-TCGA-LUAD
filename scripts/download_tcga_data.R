# Load required libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(tidyverse)

# Create data directories if they don't exist
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# Query RNA-seq STAR - Counts
query_exp <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Download RNA-seq STAR - Counts
# Note: GDCdownload() will skip already-downloaded files in cache.
GDCdownload(query_exp)
LUAD_data <- GDCprepare(query = query_exp)

# Save raw counts
counts_matrix <- assay(LUAD_data)
write.csv(counts_matrix, "data/raw/TCGA_LUAD_STAR_Counts.csv")

# Normalize to logCPM
dge <- DGEList(counts = counts_matrix)
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
write.csv(logCPM, "data/processed/TCGA_LUAD_logCPM.csv")

# Download Clinical Data (Updated for TCGAbiolinks >= 2.26+)
# Note: This replaces the old file.type = "bcr xml" usage.
# Download Clinical Data (Updated for TCGAbiolinks >= 2.26+)
clinical_data <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

# Flatten list columns
clinical_data_clean <- clinical_data %>%
  mutate(across(where(is.list), ~ sapply(., function(x) paste(unlist(x), collapse = ";"))))

# Save to CSV
write.csv(clinical_data_clean, "data/processed/TCGA_LUAD_clinical_metadata.csv", row.names = FALSE)


# Final message
cat("✅ Data download and processing complete.\n")

