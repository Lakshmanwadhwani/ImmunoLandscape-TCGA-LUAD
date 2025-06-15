# process_cibersortx_output.R
# Author: [Your Name]
# Project: ImmunoLandscape-TCGA-LUAD
# Purpose: Process CIBERSORTx output and generate immune landscape plots

# Load required packages
library(tidyverse)
library(reshape2)
library(ggplot2)

# Create figures directory if not exists
dir.create("figures", showWarnings = FALSE)

# Load CIBERSORTx results
# IMPORTANT: make sure this path matches where your CIBERSORTx_Results.csv is saved

ciber_results <- read.csv("data/processed/Immunedeconv_CIBERSORT_ABS.csv", check.names = FALSE)

# Clean: remove metrics columns (keep only immune fractions)
# If your CIBERSORTx file contains these columns, they will be removed:vg
#   P-value, Correlation, RMSE
# If your file does not have them → this line will still work safely
immune_cols <- ciber_results %>% select(-matches("P.value|Correlation|RMSE"))

# OPTIONAL: rename first column to "SampleID" if needed
colnames(immune_cols)[1] <- "SampleID"

# Reshape to long format for plotting
immune_long <- melt(immune_cols, id.vars = "SampleID",
                    variable.name = "CellType",
                    value.name = "Fraction")

# Plot: Stacked barplot of immune landscape per tumor
p <- ggplot(immune_long, aes(x = SampleID, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "Immune Landscape of TCGA-LUAD (CIBERSORTx)",
       x = "Sample",
       y = "Estimated Fraction",
       fill = "Immune Cell Type") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

# Save the plot
ggsave("figures/TCGA_LUAD_ImmuneLandscape_Barplot.png", p,
       width = 12, height = 6)

cat("✅ Immune landscape plot saved: figures/TCGA_LUAD_ImmuneLandscape_Barplot.png\n")
