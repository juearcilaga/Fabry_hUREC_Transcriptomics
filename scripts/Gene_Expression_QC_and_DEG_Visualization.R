
# Description:
# This script performs quality control (QC) and visualizes differential gene 
# expression (DEG) results for kidney gene expression data. The DEG analysis has 
# already been conducted, and this script focuses on generating PCA plots, 
# volcano plots, and visualizing the expression of cell-type marker genes.
# It uses previously obtained DESeq2 results, normalized counts from EdgeR, and 
# gene annotations for visualizing and interpreting DEG results.
#
# Author: Dr. Juliana E. Arcila-Galvis
# Date: 2024-11-12
#
# Inputs:
# - Genes_hurecs.rds: Gene annotation data with Ensembl IDs and gene names.
# - Deseq_results.csv: Differential expression results from DESeq2 analysis.
# - lcpm_TMM_counts.txt: Log-transformed normalized counts from EdgeR.
# - cpm_counts.txt: CPM counts from EdgeR, normalized gene expression data.
#
# Outputs:
# - PCA Plot: A plot showing the distribution of samples in PCA space.
# - Volcano Plot: A plot displaying differentially expressed genes.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Load Required Libraries
#-------------------------------------------------------------------------------
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(ggtext)
library(ggrepel)

#-------------------------------------------------------------------------------
# Import Ensembl Biomart Annotations for the Genes
#-------------------------------------------------------------------------------
Genes_hurecs <- readRDS("Genes_hurecs.rds")
annotations <- unique(Genes_hurecs[, c("Gene.stable.ID", "Gene.name")])
colnames(annotations) <- c("ensembl", "gene_name")

#-------------------------------------------------------------------------------
# Import Differential Expression Results from Deseq2
#-------------------------------------------------------------------------------
Deseq_results <- read.delim("Deseq_results.csv", sep = ",")
Deseq_results <- left_join(Deseq_results, annotations)

# Filter for significant results (padj < 0.05 & log2FoldChange > 1)
Deseq_fresults <- Deseq_results %>%
  filter(padj < 0.05 | is.na(padj)) %>%
  filter(abs(log2FoldChange) > 1)
# Filter for more stringent results (padj < 0.05)
Deseq_fresults_strict <- Deseq_fresults %>%
  filter(padj < 0.05)

#-------------------------------------------------------------------------------
# Import Counts Normalized with EdgeR
#-------------------------------------------------------------------------------
lcpm_TMM <- read.delim("lcpm_TMM_counts.txt")
lcpm_TMM$ensembl <- rownames(lcpm_TMM)
lcpm_annotated <- left_join(lcpm_TMM, annotations)

# Import cpm counts (normalized with EdgeR)
count_data <- read.delim("cpm_counts.txt")
count_data$ensembl <- rownames(count_data)
count_data_annotated <- left_join(count_data, annotations)

#-------------------------------------------------------------------------------
# Define Gene Markers for Various Cell Types
#-------------------------------------------------------------------------------

# Vasculature genes
vasculature_genes <- c("NRP1", "CDH5", "ELN", "PLAT", "EMCN", "TSAPN7", "MAPT", 
                       "KDR", "SMAD6", "EHD3", "LPL", "FLT1", "FLNB2", "MGP", "TRPV4", 
                       "BMX", "SOX17", "GJA5")

# Smooth Muscle Cells / Juxtaglomerular Cells genes
smcs_jgs_genes <- c("SERPINE2", "FHL2", "DES", "PRKCA", "ART3", "NT5E", "PDGFRB", 
                    "TAGLN", "MYH11", "ACTA2", "GATA3", "RERG1", "MAP3K7CL", "REN1", "AKR1B7", "RGS5")

# Podocytes genes
podocytes_genes <- c("NPHS1", "NPHS2", "SYNPO", "CDKN1C", "FOXC2", "MAFB", "EMX2", "FOXI")

# Proximal Tubule genes
proximal_tubule_genes <- c("SLC34A1", "LRP2", "HXRD2", "HSPA12", "ACSM1", "ACSM2", "CPT1A", 
                           "ACO3", "GLUD1", "PCK1", "AQP8", "HNF4A", "PPARA", "CUBN", "GCLC", 
                           "GGT1", "LAP3")

# Loop of Henle genes
loop_of_henle_genes <- c("FST", "AQP1", "SLC14A2", "BST1", "EPHA7", "CRYAB", "TSHZ2", 
                         "CLAD1", "BST1", "LYPD2", "CLDN1", "PROM1", "AKR1B1", "CLCNKA", 
                         "CLDN10", "NR3C2", "SLC19A7")

# Distal Convoluted Tubule / Collecting Duct genes
dct_cnt_genes <- c("SLC12A3", "TRPM7", "WNK1", "KLHL3", "SCNN1B", "CYP11A1", "SLC8A1", "TRPM6")
collecting_duct_genes <- c("SCNN1B", "AQP2", "AVPR2", "HSD11B2", "RHBG", "FYS", "AQP3", "ATP6V0D2")

# Immune Cells genes
immune_cells_genes <- c("CIQA", "CLYB", "ITGAM", "APOE", "CIQB", "CSF1R", "CTSS", "FCGR1", 
                        "ITGAX", "CXCL10")

# Combine all the gene markers into a single vector
genes <- c(vasculature_genes, smcs_jgs_genes, podocytes_genes, proximal_tubule_genes, 
           loop_of_henle_genes, dct_cnt_genes, collecting_duct_genes, immune_cells_genes)

# Corresponding categories for the genes
categories <- c(rep("Vasculature", length(vasculature_genes)),
                rep("SMCs/JGs", length(smcs_jgs_genes)),
                rep("Podocytes", length(podocytes_genes)),
                rep("Proximal Tubule", length(proximal_tubule_genes)),
                rep("Loop of Henle", length(loop_of_henle_genes)),
                rep("DCT/CNT", length(dct_cnt_genes)),
                rep("Collecting Duct", length(collecting_duct_genes)),
                rep("Immune Cells", length(immune_cells_genes)))

# Create a data frame for the marker genes
marker_df <- data.frame(gene_name = genes, categories = categories)

#-------------------------------------------------------------------------------
# Annotate Count Data with Gene Markers
#-------------------------------------------------------------------------------
df <- left_join(count_data_annotated[, c(2:7, 20)], marker_df)
df <- df %>% filter(!is.na(df$categories)) %>% unique()

# Reshape data for plotting
df_long <- df %>%
  gather(key = "Sample", value = "Expression", -gene_name, -categories)

# Add column for sample type (Control or Fabry)
df_long$Type <- ifelse(grepl("Fabry", df_long$Sample), "Fabry", "Control")
df_long$categories <- factor(df_long$categories, levels = unique(df_long$categories))

# Update gene names for disregarded markers (red for disregarded)
disreg_markers <- Deseq_markers %>% filter(padj < 0.05)
df_long <- df_long %>%
  mutate(gene_name = ifelse(gene_name %in% disreg_markers$gene_name, 
                            glue::glue("<i style='color:red'>{gene_name}</i>"), 
                            glue::glue("<i>{gene_name}</i>")))

#-------------------------------------------------------------------------------
# PCA Analysis of Count Data
#-------------------------------------------------------------------------------
log_counts <- log2(count_data[, c("Control1_rep1", "Control1_rep2", "Control2", "Fabry0_rep1", "Fabry0_rep2", "Fabry0_rep3")] + 1)
log_counts <- log_counts[apply(log_counts, 1, var) > 0, ]  # Remove zero variance rows

# Perform PCA
pca <- prcomp(t(log_counts), scale. = TRUE)  # Transpose so samples are rows
pca_data <- as.data.frame(pca$x)
pca_data$Sample <- rownames(pca_data)  # Add sample names

# Calculate variance explained by each PC
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_var <- round(percentVar[1], 2)
pc2_var <- round(percentVar[2], 2)

# Create PCA plot
pca_data$Type <- c("Control", "Control", "Control", "Fabry", "Fabry", "Fabry")
PCA_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(aes(color = Type), size = 2) + 
  geom_text(vjust = -1.5, size = 3) +
  labs(title = "PCA of Count Data", x = paste0("PC1: ", pc1_var, "% variance"), 
       y = paste0("PC2: ", pc2_var, "% variance")) +
  theme_minimal() + scale_color_manual(values = c(Control = "#367156", Fabry = "#b87800"))

# Print PCA plot
print(PCA_plot)

#-------------------------------------------------------------------------------
# Volcano Plot for Differential Expression
#-------------------------------------------------------------------------------
Deseq_results$Significance <- ifelse(Deseq_results$padj < 0.05 & abs(Deseq_results$log2FoldChange) > 1, 
                                     ifelse(Deseq_results$log2FoldChange > 1, "Upregulated", "Downregulated"),
                                     "Not Significant")

# Adjust for y-axis values
Deseq_results <- Deseq_results %>%
  mutate(neg_log10_padj = -log10(padj),
         adjusted_y = if_else(neg_log10_padj > 10, 10 + log2(neg_log10_padj - 10), neg_log10_padj))

# Identify most significant genes
top_genes <- Deseq_results %>%
  filter(Significance != "Not Significant") %>%
  arrange(padj) %>%
  slice_head(n = 10) %>%
  bind_rows(
    Deseq_results %>%
      filter(Significance != "Not Significant") %>%
      arrange(padj) %>%
      filter(log2FoldChange < 0) %>%
      slice_head(n = 10)
  )

# Create Volcano plot
volcano_plot <- ggplot(Deseq_results, aes(x = log2FoldChange, y = adjusted_y)) +
  geom_point(aes(fill = Significance, color = Significance), size = 3, alpha = 0.4, shape = 21, stroke = 0.9) +
  scale_fill_manual(values = c("Upregulated" = "#C0272C", "Downregulated" = "#395A88", "Not Significant" = "grey")) +
  scale_color_manual(values = c("Upregulated" = "#C0272C", "Downregulated" = "#395A88", "Not Significant" = "grey")) +
  geom_text_repel(data = top_genes, aes(label = gene_name), size = 3, max.overlaps = Inf, nudge_y = 1) +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 padj (Compressed above 10)") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))

# Print Volcano plot
print(volcano_plot)
