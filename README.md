# Integrated-Transcriptomics-and-Survival-Analysis-of-TCGA-LUAD-using-DESeq2-and-GSEA
RNA-seq analysis pipeline of TCGA-LUAD using DESeq2 for differential expression, clusterProfiler for GO/KEGG enrichment, GSEA for pathway analysis, and survival analysis for prognostic genes. Identifies key dysregulated pathways and potential biomarkers in lung adenocarcinoma.
# 🧬 TCGA-LUAD RNA-Seq Analysis Pipeline

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![TCGA](https://img.shields.io/badge/TCGA-LUAD-blue?style=for-the-badge)
![DESeq2](https://img.shields.io/badge/DESeq2-Differential%20Expression-orange?style=for-the-badge)
![License](https://img.shields.io/badge/License-MIT-green?style=for-the-badge)

---

## 📌 Overview

This project presents an end-to-end RNA-seq bioinformatics pipeline for analyzing **lung adenocarcinoma (TCGA-LUAD)** data. It integrates differential expression analysis, functional enrichment, gene set enrichment analysis (GSEA), and survival analysis to identify key molecular pathways and potential prognostic biomarkers.

---

## 🎯 Objectives

- Identify differentially expressed genes (DEGs) in LUAD vs normal tissue  
- Perform functional enrichment analysis (GO & KEGG)  
- Conduct Gene Set Enrichment Analysis (GSEA)  
- Evaluate prognostic significance using survival analysis  
- Visualize transcriptomic patterns and pathway alterations  

---

## 🧪 Methodology

The workflow is implemented in **R** using Bioconductor packages:

- `TCGAbiolinks` → Data acquisition  
- `DESeq2` → Differential expression analysis  
- `clusterProfiler` → GO & KEGG enrichment  
- `org.Hs.eg.db` → Gene annotation  
- `enrichplot` → Visualization  
- `survival`, `survminer` → Survival analysis  

---

## 📊 Workflow Summary

1. Download TCGA-LUAD RNA-seq count data  
2. Preprocess and define tumor vs normal groups  
3. Perform differential expression analysis (DESeq2)  
4. Identify significant genes and visualize (volcano plot)  
5. Conduct GO and KEGG enrichment analysis  
6. Perform GSEA using ranked gene lists  
7. Run Kaplan–Meier survival analysis for selected genes  
8. Visualize results (dot plots, ridge plots, enrichment plots)  

---

## 📈 Key Outputs

- Differentially expressed genes (DEGs)  
- Volcano plot of gene expression changes  
- GO and KEGG enriched pathways  
- GSEA enrichment and ridge plots  
- Kaplan–Meier survival curves  

---

## 🧠 Biological Insights

This analysis identifies dysregulated biological processes in lung adenocarcinoma, including cell cycle regulation, immune response, and metabolic reprogramming. Survival analysis highlights potential prognostic biomarkers linked to patient outcomes.

---

## ⚙️ Requirements

- R (≥ 4.0)
- Bioconductor packages:
  - `TCGAbiolinks`
  - `DESeq2`
  - `clusterProfiler`
  - `org.Hs.eg.db`
  - `enrichplot`
  - `survival`
  - `survminer`
  - `ggplot2`

---

## 🚀 Future Work

- Multi-gene survival risk modeling  
- Integration with clinical variables (stage, age, gender)  
- Machine learning-based biomarker discovery  
- Validation using external GEO datasets  

---

## 📂 Reproducibility

This pipeline is fully reproducible and can be adapted to other TCGA cancer types such as BRCA, COAD, and others with minimal modifications.

---

## 📜 License

This project is licensed under the MIT License.

---

## 👨‍🔬 Author

Bioinformatics RNA-seq analysis project for TCGA-LUAD.
