# Fabry_hUREC_Transcriptomics
This repository contains scripts for the analysis of bulk RNA sequencing (RNAseq) data from human urine-derived renal epithelial cells (hURECs) in Fabry disease. The scripts cover key steps in transcriptomic analysis, including quality control, normalization, differential expression analysis, and visualization.


## 1. RNASeq_Fabry_Differential_Expression_Analysis.R

Performs differential gene expression (DEG) analysis using RNA sequencing (RNA-Seq) data from Fabry disease samples. Involves normalization, statistical testing (DESeq2), and identification of significantly differentially expressed genes.

## 2. Gene_Expression_QC_and_DEG_Visualization.R

Handles quality control (QC) checks on RNA-Seq data. Generates visualizations of PCA plots, celltype gene marker expression scatter plots, and volcano plots.

## 3. gene_set_enrichment_analysis.R

Performs Gene Set Enrichment Analysis (GSEA) to determine whether predefined gene sets (e.g., KEGG pathways, GO terms) are enriched in Fabry disease samples.
Uses tools like GSEA, fgsea, or clusterProfiler for pathway analysis.

## 4. Heatmaps.R

Creates heatmaps of gene expression data for visualization of DEGs.
Uses clustering techniques to group genes or samples based on expression patterns.

## 5. UPC_lnregression.R

Performs linear regression analysis on Urine Protein Creatinine (UPC) ratio over time.

## 6. eGFR_lnregression.R

Similar to the UPC analysis but focused on estimated Glomerular Filtration Rate (eGFR) trends over time.
Uses regression modeling to study kidney function decline and the impact of treatments.
