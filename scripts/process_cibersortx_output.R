# process_cibersortx_output.R
# Author: [Your Name]
# Project: ImmunoLandscape-TCGA-LUAD
# Purpose: Process EPIC deconvolution results and generate immune landscape plot

# -------------------------------------
# Load required packages
# -------------------------------------
library(tidyverse)  # includes dplyr, tidyr, ggplot2, etc.

# -------------------------------------
# Create figures directory if it doesn't exist
# -------------------------------------
dir.create("figures", showWarnings = FALSE)

# -------------------------------------
# Load EPIC immunedeconv results (wide format)
# -------------------------------------
epic_raw <- read.csv(
  "data/processed/Immunedeconv_EPIC.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# -------------------------------------
# Pivot from wide → long format
# -------------------------------------
immune_long <- epic_raw %>%
  pivot_longer(
    cols      = -cell_type,
    names_to  = "SampleID",
    values_to = "Fraction"
  ) %>%
  rename(CellType = cell_type)

# -------------------------------------
# Plot: stacked barplot of immune fractions per tumor
# -------------------------------------
p <- ggplot(immune_long, aes(x = SampleID, y = Fraction, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(
    title    = "Immune Landscape of TCGA-LUAD (EPIC Deconvolution)",
    x        = "Sample",
    y        = "Estimated Fraction",
    fill     = "Cell Type"
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  )

# -------------------------------------
# Save the plot
# -------------------------------------
ggsave(
  filename = "figures/TCGA_LUAD_ImmuneLandscape_Barplot.png",
  plot     = p,
  width    = 12,
  height   = 6
)

cat("✅ Immune landscape plot saved: figures/TCGA_LUAD_ImmuneLandscape_Barplot.png\n")

