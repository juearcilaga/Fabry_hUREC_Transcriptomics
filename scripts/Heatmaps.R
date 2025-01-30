# -------------------------------------------------------------------------
# R Script: Heatmaps.R
# Description:
# This script generates heatmaps to visualize the expression patterns of significant genes across different samples in a Fabry disease study. 
# It uses DESeq2 results and filtered gene expression data (LCPM, CPM) to select genes and plot heatmaps based on their clusters and conditions.
#
# Author: Dr. Juliana E. Arcila-Galvis
# Date: 2024-12-12
#
# Inputs:
# - Deseq_results.csv: DESeq2 results file containing differential expression analysis.
# - Genes_hurecs.rds: A dataset containing gene annotations with gene stable IDs, gene names, and Entrez gene IDs.
# - lcpm_counts.txt: Raw gene expression data in log-transformed CPM (LCPM) format.
# - cpm_counts.txt: Raw gene expression data in CPM (counts per million) format.
# - heatmap_genes.rds: A pre-defined list of genes with associated GO cluster terms to be plotted in the heatmap.
#
# Outputs:
# - Heatmap visualizations: Heatmaps of gene expression across samples with annotations for different conditions and clusters.
# -------------------------------------------------------------------------

# Load necessary libraries
library(pheatmap)
library(dplyr)

# Set working directory (adjust the path as necessary)
setwd("C:/Users/nja147/OneDrive - Newcastle University/POSDOC/Sayer/projects/NPHP1/Marco_handout_data/Fabry_project/Juliana")

# Read the DESeq2 results and filter for significant results
Deseq_results <- read.delim("Deseq_results.csv", sep = ",")
Deseq_fresults <- Deseq_results %>% filter(padj < 0.05 | is.na(padj) & abs(log2FoldChange) > 1)

# Read gene expression data in LCPM and CPM format
lcpm <- read.delim("lcpm_counts.txt")
cpm <- read.delim("cpm_counts.txt")

# Load gene annotations
Genes_hurecs <- readRDS("Genes_hurecs.rds")

# Prepare annotations for heatmap
annotations <- unique(Genes_hurecs[, c("Gene.stable.ID", "Gene.name", "NCBI.gene..formerly.Entrezgene..ID")])
colnames(annotations) <- c("ensembl", "gene_name", "entrez")

# Define sample names and treatments
samples <- c("Fabry1_t0", "Fabry2_t0", "Fabry3_t0", "Control1_rep1", "Control1_rep2", "Control2")
patient_treatment <- c("Fabry_t1", "Fabry_t2", "Fabry_t3", "Fabry_t4", "Fabry_t5")
culture_treatment <- c("Fabry_t1_invitro", "Fabry_t2_invitro", "Fabry_t3_invitro", "Fabry_t4_invitro")

# Filter expression data for significant genes
pathway1_expression_data <- data.frame(lcpm[rownames(lcpm) %in% Deseq_fresults$ensembl, 
                                            colnames(lcpm) %in% c(samples[1:3], patient_treatment, culture_treatment, samples[4:6])])
pathway1_expression_data$ensembl <- rownames(pathway1_expression_data)

# Order the expression data by the DESeq2 results
pathway1_expression_data <- pathway1_expression_data[match(Deseq_fresults$ensembl, pathway1_expression_data$ensembl), ]
pathway1_expression_data <- left_join(pathway1_expression_data, annotations)

# Filter genes for heatmap generation
heatmap_genes <- readRDS("heatmap_genes.rds")
heatmap_genes[heatmap_genes$Gene.name %in% c("FASN", "LIPA", "PTK7", "UCHL1")]$cluster_GO <- "carboxylic acid metabolic process"
heatmap_genes[heatmap_genes$Gene.name %in% c("ADAMTS7", "COL6A1", "ENG", "NDST1", "PDLIM7")]$cluster_GO <- "regulation of transmembrane receptor protein serine/threonine kinase signaling pathway"

# Remove genes with "axon guidance" cluster
heatmap_genes <- heatmap_genes[heatmap_genes$cluster_GO != "axon guidance", ]
heatmap_genes <- unique(heatmap_genes)

# Merge expression data with heatmap genes
data <- pathway1_expression_data[pathway1_expression_data$gene_name %in% heatmap_genes$Gene.name, ]
rownames(data) <- data$gene_name
data$gene_name <- as.factor(data$gene_name)

# Prepare data for heatmap (selecting specific samples)
matrix_data <- data[, c("Fabry1_t0", "Fabry2_t0", "Fabry3_t0", "Control1_rep1", "Control1_rep2", "Control2")]
colnames(matrix_data) <- c("Fabry_rep1", "Fabry_rep2", "Fabry_rep3", "Control1_rep1", "Control1_rep2", "Control2_rep1")

# Hierarchical clustering of rows
hClustering <- matrix_data %>% dist %>% hclust
order <- hClustering$order
order <- rownames(matrix_data)[order]
order <- order[c(1:12, 21:25, 35:40, 13:20, 26:34)]
matrix_data <- matrix_data[match(order, rownames(matrix_data)), ]

# Define annotations for rows and columns
annotations_rows_type <- data.frame(cluster_GO = heatmap_genes[,"cluster_GO"])
rownames(annotations_rows_type) <- heatmap_genes$Gene.name

annotations_col_type <- data.frame(sample_type = as.factor(c(rep("Fabry", 3), rep("Control", 3))))
rownames(annotations_col_type) <- colnames(matrix_data)[1:6]

# Color palette for clustering
color_clusters <- c("#9E9D24", "#374c80", "#7a5195", "#bc5090", "#ef5675", "#ff764a", "#ffa600")
names(color_clusters) <- unique(annotations_rows_type$cluster_GO)

ann_colors <- list(cluster_GO = color_clusters, sample_type = c(Control = "#367156", Fabry = "#b87800"))

# Generate the heatmap
pheatmap(
  matrix_data,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "row",
  annotation_row = annotations_rows_type,
  annotation_col = annotations_col_type,
  annotation_colors = ann_colors,
  border_color = "white",
  cellwidth = 10,
  cellheight = 12,
  treeheight_col = 12,
  treeheight_row = 12,
  color = colorRampPalette(c("#395A88", "white", "#C0272C"))(1000),
  fontsize = 10
)

# Additional heatmap generation with updated configurations for temporal conditions and treatments
matrix_data <- data[, c("Fabry1_t0", "Fabry2_t0", "Fabry3_t0", "Fabry_t1", "Fabry_t2", "Fabry_t3", "Fabry_t4", "Fabry_t5", "Control1_rep1", "Control1_rep2", "Control2")]
colnames(matrix_data) <- c("t0_rep1", "t0_rep2", "t0_rep3", "t1", "t2", "t3", "t4", "t5")

matrix_data <- matrix_data[match(order, rownames(matrix_data)), ]

annotations_col_type <- data.frame(sample_type = as.factor(c("Untreated", "Untreated", "Untreated", "Galafold", "Galafold", "Galafold", "Galafold", "Chaperone", "Galafold_miglastat", "Galafold_miglastat", "Galafold_miglastat", "Galafold_miglastat", "Control", "Control", "Control")))
rownames(annotations_col_type) <- colnames(matrix_data)

ann_colors <- list(cluster_GO = color_clusters, sample_type = c(Untreated = "#b87800", Galafold = "#155c76", Chaperone = "#76151f", Control = "#367156"))

pheatmap(
  matrix_data,
  clustering_method = "median",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = annotations_col_type,
  annotation_colors = ann_colors,
  cutree_cols = 4,
  border_color = "white",
  cellwidth = 10,
  cellheight = 12,
  treeheight_col = 12,
  treeheight_row = 12,
  color = colorRampPalette(c("#395A88", "white", "#C0272C"))(1000),
  fontsize = 10
)
