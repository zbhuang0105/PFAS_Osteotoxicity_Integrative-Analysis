# PFAS Toxicity in Osteoporosis and Osteosarcoma: Integrative Analysis

This repository contains the source code for the computational analysis presented in the manuscript:
**"From Network Topology to Atomic Dynamics: An Integrative Computational Strategy Decodes the Structural Basis of PFAS Toxicity in Osteoporosis and Osteosarcoma"**

## 📂 Repository Contents

* **scRNA_seq_analysis.R**: Main R script for single-cell RNA sequencing analysis.
    * Quality Control & Normalization (Seurat)
    * Batch Effect Removal (Harmony)
    * Dimensionality Reduction (PCA, t-SNE, UMAP)
    * Cell Type Annotation & Marker Identification
    * Visualization (Violin plots, Dot plots, Feature plots)

## 🔧 Prerequisites & Dependencies

The analysis was performed using **R**. Please ensure the following packages are installed:

```r
install.packages(c("Seurat", "tidyverse", "Matrix", "dplyr", "patchwork", "ggplot2", "cowplot", "harmony", "devtools"))
# Note: Some packages (e.g., SingleR, monocle, irGSEA) may require installation via BiocManager or GitHub.
