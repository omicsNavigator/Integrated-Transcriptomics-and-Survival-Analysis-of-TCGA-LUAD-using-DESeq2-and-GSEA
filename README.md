# Integrated-Transcriptomics-and-Survival-Analysis-of-TCGA-LUAD-using-DESeq2-and-GSEA
RNA-seq analysis pipeline of TCGA-LUAD using DESeq2 for differential expression, clusterProfiler for GO/KEGG enrichment, GSEA for pathway analysis, and survival analysis for prognostic genes. Identifies key dysregulated pathways and potential biomarkers in lung adenocarcinoma.

# TCGA-LUAD RNA-Seq Analysis Pipeline
## Overview

This project focuses on the analysis of RNA-seq data from TCGA-LUAD (lung adenocarcinoma). The workflow covers differential gene expression analysis, functional enrichment, pathway analysis, and survival analysis to explore molecular changes associated with tumor progression.

The pipeline is designed to process raw count data and extract biologically meaningful insights using standard bioinformatics tools in R.

---

## Objectives

- Identify differentially expressed genes between tumor and normal samples  
- Perform GO and KEGG enrichment analysis to identify affected biological pathways  
- Conduct Gene Set Enrichment Analysis (GSEA) using ranked gene lists  
- Evaluate the association of gene expression with patient survival  
- Visualize key results using standard plots  

---

## Workflow

1. Download TCGA-LUAD RNA-seq count data using TCGAbiolinks  
2. Preprocess data and define tumor vs normal groups  
3. Perform differential expression analysis using DESeq2  
4. Filter and extract significant genes  
5. Conduct GO and KEGG enrichment analysis using clusterProfiler  
6. Generate ranked gene list for GSEA  
7. Perform GSEA for GO and KEGG pathways  
8. Perform Kaplan–Meier survival analysis using clinical data  
9. Visualize results using ggplot2 and enrichplot  

---

## Tools and Packages

- TCGAbiolinks  
- DESeq2  
- clusterProfiler  
- org.Hs.eg.db  
- enrichplot  
- survival  
- survminer  
- ggplot2  

---

## Outputs

- Differentially expressed gene table  
- Volcano plot of gene expression changes  
- GO and KEGG enriched pathways  
- GSEA enrichment plots and ridge plots  
- Kaplan–Meier survival curves for selected genes  

---

## Notes

The pipeline can be applied to other TCGA cancer datasets with minor modifications to the project ID. It is structured for reproducibility and further extension into multi-omics or clinical integration studies.

---

## License

This project is intended for academic and research use.
