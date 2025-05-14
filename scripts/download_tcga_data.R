#!/usr/bin/env Rscript

# ------------------------------
# Download and Preprocess TCGA-LUAD RNA-seq Data
# ------------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(edgeR)
library(tidyverse)

# ------------------------------
# Step 1: Set output directories
# ------------------------------
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ------------------------------
# Step 2: Query TCGA-LUAD expression data
# ------------------------------
query_exp <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)

# ------------------------------
# Step 3: Download data
# ------------------------------
GDCdownload(query_exp)

# ------------------------------
# Step 4: Prepare SummarizedExperiment
# ------------------------------
LUAD_data <- GDCprepare(query = query_exp)

# ------------------------------
# Step 5: Extract raw count matrix
# ------------------------------
counts_matrix <- assay(LUAD_data)
write.csv(counts_matrix, "data/raw/TCGA_LUAD_HTSeq_Counts.csv")

# ------------------------------
# Step 6: Normalize to log2 CPM using edgeR
# ------------------------------
dge <- DGEList(counts = counts_matrix)
dge <- calcNormFactors(dge, method = "TMM")
logCPM <- cpm(dge, log = TRUE, prior.count = 1)
write.csv(logCPM, "data/processed/TCGA_LUAD_logCPM.csv")

# ------------------------------
# Step 7: Save clinical metadata
# ------------------------------
clinical <- colData(LUAD_data) %>% as.data.frame()
write.csv(clinical, "data/processed/TCGA_LUAD_clinical_metadata.csv")

cat("✅ Data download and processing complete.\n")
