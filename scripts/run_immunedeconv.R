# run_immunedeconv.R
# Author: [Your Name]
# Project: ImmunoLandscape-TCGA-LUAD
# Purpose: Run immune deconvolution locally using immunedeconv package

# -------------------------------------
# Install immunedeconv if not installed
# -------------------------------------
if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("icbi-lab/immunedeconv")
}

# -------------------------------------
# Load required libraries
# -------------------------------------
library(immunedeconv)
library(tidyverse)

# -------------------------------------
# Create output directory
# -------------------------------------
dir.create("data/processed", showWarnings = FALSE)

# -------------------------------------
# Load logCPM matrix
# -------------------------------------
logCPM <- read.csv("data/processed/TCGA_LUAD_logCPM.csv", row.names = 1, check.names = FALSE)

# -------------------------------------
# Run CIBERSORT ABS deconvolution
# -------------------------------------
cat("Running CIBERSORT ABS deconvolution locally...\n")

res_cibersort_abs <- deconvolute(logCPM, method = "cibersort_abs")

# -------------------------------------
# Save results
# -------------------------------------
write.csv(res_cibersort_abs,
          "data/processed/Immunedeconv_CIBERSORT_ABS.csv",
          row.names = FALSE)

cat("✅ Immunedeconv results saved: data/processed/Immunedeconv_CIBERSORT_ABS.csv\n")