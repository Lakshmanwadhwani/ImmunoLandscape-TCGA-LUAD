# ImmunoLandscape-TCGA-LUAD
## ImmunoLandscape: TCGA-LUAD Immune Microenvironment Pipeline

This project analyzes the immune microenvironment of lung adenocarcinoma (LUAD) using transcriptional profiling and immune deconvolution.

### Overview

- **Cohort:** TCGA-LUAD (The Cancer Genome Atlas)
- **Data:** RNA-seq (STAR - Counts), Clinical metadata
- **Goal:** Infer immune cell type proportions in tumor samples

---

### Pipeline Steps

✅ **1. Download & preprocess RNA-seq data**  
[`scripts/download_tcga_data.R`](scripts/download_tcga_data.R)

- Download TCGA-LUAD RNA-seq raw counts (STAR - Counts)
- Normalize to log2 CPM (logCPM)
- Save `TCGA_LUAD_logCPM.csv`

✅ **2. Prepare input for immune deconvolution**  
[`scripts/format_for_cibersortx.R`](scripts/format_for_cibersortx.R)

- Format logCPM matrix for deconvolution
- Save `TCGA_LUAD_logCPM_CIBERSORTx.tsv`

✅ **3. Run immune deconvolution (local option)**  
[`scripts/run_immunedeconv.R`](scripts/run_immunedeconv.R)

- Run immune deconvolution locally using [`immunedeconv`](https://github.com/icbi-lab/immunedeconv)
- Methods supported: CIBERSORT ABS, EPIC, quanTIseq, xCell, MCP-counter
- Save `Immunedeconv_CIBERSORT_ABS.csv`

✅ **4. Visualize immune landscape**  
[`scripts/process_cibersortx_output.R`](scripts/process_cibersortx_output.R)

- Generate stacked barplot of immune cell proportions per tumor
- Save `figures/TCGA_LUAD_ImmuneLandscape_Barplot.png`

---

### Notes

- This project currently supports **fully local analysis** (no CIBERSORTx registration required).
- You can optionally run the pipeline using CIBERSORTx online if an academic email is available.
- The `immunedeconv` package provides several state-of-the-art deconvolution methods.

---

### Example Output

![Immune Landscape Barplot](figures/TCGA_LUAD_ImmuneLandscape_Barplot.png)

---

### References

- [CIBERSORTx](https://cibersortx.stanford.edu/)
- [TCGAbiolinks R package](https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html)
- [immunedeconv R package](https://github.com/icbi-lab/immunedeconv)

---

### License
## Future Work

- **Multi‐omics Integration**  
  Extend the pipeline to incorporate somatic mutation (MAF) and copy-number (CNV) data from TCGA for combined analysis of genetic alterations and immune profiles.

- **Pan-Cancer Comparison**  
  Apply the same workflow across multiple TCGA cohorts (e.g., LUAD vs. LUSC or BRCA) to identify shared and tumor-type-specific immune-microenvironment signatures.

- **Single-Cell Validation**  
  Integrate publicly available single-cell RNA-seq datasets (e.g., Human Tumor Atlas) to validate bulk deconvolution estimates and refine cell-type signature matrices.

- **Survival and Clinical Correlation**  
  Perform Kaplan–Meier and Cox regression analyses to link inferred immune cell proportions with overall and progression-free survival, as well as treatment response metadata.

- **Machine Learning Classifiers**  
  Train predictive models (e.g., random forests, neural networks) using immune-deconvolution features to classify tumors into “immune hot” vs. “immune cold” and predict immunotherapy response.

- **Interactive Dashboard**  
  Develop a Shiny or Dash application that lets users explore immune profiles, clinical annotations, and survival curves in an interactive web interface.

- **Benchmarking & Method Comparison**  
  Systematically compare multiple deconvolution methods (CIBERSORT, EPIC, quanTIseq, xCell, MCP-counter) on the same dataset to assess robustness and concordance.

- **Extension to Spatial Data**  
  Adapt the workflow to integrate spatial transcriptomics data, enabling mapping of immune cell spatial distributions within the tumor microenvironment.

> _These enhancements will further deepen biological insights, improve clinical utility, and broaden the impact of the ImmunoLandscape pipeline._



For non-commercial academic use only.

---

