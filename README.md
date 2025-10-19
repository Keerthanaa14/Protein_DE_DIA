# Protein_DE_DIA
R script for data independent MS differential analysis
# Protein Differential Expression Analysis (limma) suitable for data independent MS

This repository provides a full reproducible pipeline for differential protein expression analysis using **limma** and **imputeLCMD**.

## 📋 Features
- Handles missing values using MinDet imputation
- Performs log2 transformation
- Runs pairwise group comparisons with `limma`
- Generates:
  - Differential expression tables (CSV)
  - Volcano plots
  - Bar plots of top up/down proteins
  - Combined heatmaps of significant proteins

## 🧰 Requirements
R (≥ 4.2) with the following packages:
```r
tidyverse
limma
imputeLCMD
pheatmap
RColorBrewer
