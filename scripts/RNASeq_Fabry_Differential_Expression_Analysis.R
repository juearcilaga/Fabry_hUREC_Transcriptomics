# ------------------------------------------------------------------------------ 
# RNA-Seq Data Analysis Script
# Author: Dr. Juliana E. Arcila-Galvis
# Date: 29/10/2024
# 
# Description:
# This script processes RNA-Seq data for differential gene expression analysis.
# It loads raw count data, performs quality control, normalizes the data, and 
# prepares it for downstream differential expression analysis.
#
# Inputs:
# - Fabry_experiments_batch1.RData (raw counts and metadata from first batch)
# - Fabry_experiments_batch2.RData (raw counts and metadata from second batch)
# - mart_export_GRCh38.p14.txt (gene annotations from Ensembl Biomart)
#
# Outputs:
# - raw_counts.rds (processed raw count matrix)
# - Genes_hurecs.rds (gene annotations)
# - raw_counts_noMito_widelyexp.rds (filtered gene expression data excluding mitochondrial genes)
# - sampleInfo.rds (sample metadata)
# ------------------------------------------------------------------------------

# Load required libraries for RNA-Seq analysis
library(DESeq2)    # For differential expression analysis
library(edgeR)     # For working with count data
library(dplyr)     # For data manipulation
library(RColorBrewer)  # For color schemes
library(ggplot2)   # For plotting
library(ggfortify) # For PCA visualization
library(stats)     # For statistical tests
library(viridis)   # For color palettes

# Set the working directory (change this path as needed)
setwd("Fabry_project/")

# ------------------------------------------------------------------------------ 
# Step 1: Load Raw Count Data
# ------------------------------------------------------------------------------ 

# This section loads two RNA-Seq datasets: one from Fabry experiments and another from hURECs of different patients.
# The raw counts data is extracted, and sample metadata is stored.

# Load the first dataset (Batch 1)
load("Fabry_experiments_batch1.RData")
raw_counts <- assay(se)  # Extract raw counts from the first dataset
sampleInfo1 <- sampleInfo  # Store sample metadata for Batch 1

# Load the second dataset (Batch 2)
load("Fabry_experiments_batch2.RData")
raw_counts_newsample <- assay(se)  # Extract raw counts from the second dataset
colnames(raw_counts_newsample) <- c("Fabry_t5", "other_hUREC_4")  # Rename columns for clarity
sampleInfo2 <- sampleInfo  # Store sample metadata for Batch 2

# Combine datasets from both batches
raw_counts2 <- cbind(raw_counts, raw_counts_newsample)

# Reorder columns by donor name to maintain consistency
columnames <- c(
  "Control0", "Control1_rep1", "Control1_rep2", "Control2", "AS_1", "AS_2", "AS_3", "Fabry_t1",
  "Fabry_t2", "Fabry_t3", "Fabry_t4", "Fabry_t5", "Fabry_t1_treat", "Fabry_t2_invitro",
  "Fabry_t3_invitro", "Fabry_t4_invitro", "other_hUREC_1", "other_hUREC_2", "other_hUREC_3", "other_hUREC_4"
)
raw_counts2 <- raw_counts2[, match(columnames, colnames(raw_counts2))]
colnames(raw_counts2)[2] <- "Control1_rep1"  # Rename the second column

# Save processed raw count data
saveRDS(raw_counts2, "raw_counts.rds")

# ------------------------------------------------------------------------------ 
# Step 2: Load Gene Annotations and Filter Mitochondrial Genes
# ------------------------------------------------------------------------------ 

# Load gene annotations from Ensembl Biomart
Biomart <- read.delim("mart_export_GRCh38.p14.txt")

# Assign Ensembl ID as gene name where gene symbol is missing
Genes_hurecs <- Biomart[Biomart$Gene.stable.ID %in% rownames(raw_counts), ]
Genes_hurecs$Gene.name[Genes_hurecs$Gene.name == ""] <- Genes_hurecs$Gene.stable.ID[Genes_hurecs$Gene.name == ""]

# Save processed gene annotation data
saveRDS(Genes_hurecs, "Genes_hurecs.rds")

# Filter out mitochondrial genes based on gene annotations
Genes_hurecs <- readRDS("Genes_hurecs.rds")
raw_counts_noMito <- raw_counts2[rownames(raw_counts2) %in% Genes_hurecs$Gene.stable.ID, ]

# ------------------------------------------------------------------------------ 
# Step 3: Merge Metadata and Define Batches
# ------------------------------------------------------------------------------ 

# Combine metadata from both batches
sampleInfo <- rbind(sampleInfo1, sampleInfo2)

# Assign batch information based on sequencing runs (this might vary depending on your experimental design)
batch_labels <- c(
  "X204SC24083061-Z01-F002", "X204SC23032810-Z01-F002", 
  "X204SC23073497-Z01-F003_01", "X204SC23063379-Z01-F001",
  "X204SC24010986-Z01-F004_01", "X204SC23073497-Z01-F006"
)

# Add batch information to the sample metadata based on the path of each sample
for (batch in batch_labels) {
  sampleInfo$BATCH[grepl(batch, sampleInfo$path)] <- batch
}

# Select relevant samples for downstream analysis
selected_columns <- c(
  "Control1_rep1", "Control1_rep2", "Control2", "AS_1", "AS_2", "AS_3", 
  "Fabry_t1", "Fabry_t2", "Fabry_t3", "Fabry_t4", "Fabry_t5", "Fabry_t1_treat",
  "Fabry_t2_invitro", "Fabry_t3_invitro", "Fabry_t4_invitro"
)
raw_counts_noMito <- raw_counts_noMito[, colnames(raw_counts_noMito) %in% selected_columns]

# ------------------------------------------------------------------------------ 
# Step 4: Filter and Process Data
# ------------------------------------------------------------------------------ 

# Remove genes that have zero counts across all samples
nsamples <- length(colnames(raw_counts_noMito))
zeros <- table(rowSums(raw_counts_noMito == 0) == nsamples)

# Print the number of genes with zero counts across all samples
print("Features with zero counts across all samples:")
print(zeros)

# Remove genes with zero counts across all samples
raw_counts_noMito_noZeros <- raw_counts_noMito[(rowSums(raw_counts_noMito == 0) == nsamples) == FALSE, ]
dim(raw_counts_noMito_noZeros)

# Prepare counts for edgeR analysis
counts <- raw_counts_noMito_noZeros
dim(counts)  # Check the dimensions of the counts matrix

# Initialize DGEList object for edgeR
counts <- DGEList(counts = counts, lib.size = colSums(counts), 
                  norm.factors = rep(1, ncol(counts)), 
                  samples = NULL, group = NULL, 
                  genes = NULL, remove.zeros = FALSE)

# Organize sample information
metadata <- sampleInfo
colnames(metadata)[1] <- "Sample"
metadata$condition <- paste(metadata$type, metadata$treatment, sep="_")

# Merge metadata with counts object
counts$samples$Sample <- rownames(counts$samples)
counts$samples <- left_join(counts$samples, metadata, by = "Sample")

# Save merged sample information
saveRDS(counts$samples, "sampleInfo.rds")

# ------------------------------------------------------------------------------ 
# Step 5: Transform Data to Counts Per Million (CPM)
# ------------------------------------------------------------------------------ 

# Calculate counts per million (CPM)
cpm <- cpm(counts)

# Filter low-expressed genes by setting a cutoff threshold
cutoff <- 1
drop <- which(apply(cpm, 1, max) < cutoff)
counts <- counts[-drop, ]

# Log-transformed counts (log-CPM)
lcpm <- cpm(counts, log = TRUE)

# Export the normalized counts for visualizations (e.g., PCA, heatmaps)
write.table(cpm, "cpm_counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)
write.table(lcpm, "lcpm_counts.txt", sep = "\t", row.names = TRUE, quote = FALSE)

# ------------------------------------------------------------------------------ 
# Step 6: Differential Expression Analysis
# ------------------------------------------------------------------------------ 

# Select samples for the baseline comparison (control vs. Fabry)
selected_samples <- c("AS_1", "AS_2", "AS_3", "Control1_rep1", "Control1_rep2", "Control2")

# Uncomment the next line for temporal comparison instead
# selected_samples <- c("AS_1", "AS_2", "AS_3", "Fabry_t1", "Fabry_t2", "Fabry_t3", "Fabry_t4", "Fabry_t5")

# Filter counts and sample metadata for the selected samples
countsdds <- counts
countsdds$samples <- counts$samples[counts$samples$Sample %in% selected_samples, ]
countsdds$counts <- countsdds$counts[, colnames(countsdds$counts) %in% selected_samples]

# Create a DESeqDataSet for differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = countsdds$counts,
                              colData = countsdds$samples,
                              design = ~ type)

# Perform pre-filtering to retain genes with at least 10 counts in a minimal number of samples
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep, ]

# Set the factor levels for the comparison (control vs. Fabry)
dds$type <- factor(dds$type, levels = c("control", "fabry"))

# Normalize the data and run DESeq2 analysis
dds <- estimateSizeFactors(dds)
sizeFactors(dds)  # Check size factors
normalized_counts <- counts(dds, normalized = TRUE)

# Save the normalized counts for further analysis
write.table(normalized_counts, file = "Deseq_norm_counts.txt", sep = "\t", quote = FALSE, col.names = NA)

# Perform differential expression analysis
dds <- DESeq(dds)
resultsNames(dds)  # List the coefficients for differential analysis
res <- results(dds, name = "type_fabry_vs_control")
res$ensembl = rownames(res)

# Order the results by p-value
resOrdered <- res[order(res$pvalue), ]
summary(res)

# Count the number of significant results (adjusted p-value < 0.05)
sum(res$padj < 0.05, na.rm = TRUE)

# Filter results with adjusted p-value < 0.05
res05 <- results(dds, alpha = 0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm = TRUE)

# Merge gene annotations with differential expression results
colnames(Genes_hurecs)[c(1, 7)] <- c("ensembl", "Entrez")
dfres <- as.data.frame(resOrdered)
dfres <- left_join(dfres, unique(Genes_hurecs[c("ensembl", "Gene.name", "Entrez")]), by = "ensembl") %>% unique

# Save the results as CSV
write.csv(dfres, file = "Deseq_results.csv")

# Filter significant results (padj < 0.05) and calculate fold change (FC)
resSig <- subset(resOrdered, padj < 0.05)
dfressig <- as.data.frame(resSig)
dfressig <- left_join(dfressig, unique(Genes_hurecs[c("ensembl", "Gene.name", "Entrez")]), by = "ensembl") %>% unique
dfressig$FC <- 2^dfressig$log2FoldChange

# Further filter based on fold change thresholds (|FC| > 2)
dfressig <- dfressig[dfressig$FC > 2 | dfressig$FC < 0.5, ]

# Apply additional filter based on log fold change standard error
dfressig <- dfressig[dfressig$lfcSE < 0.5, ]

# Save significant results
write.csv(dfressig, file = "Deseq_results_significant.csv")

# Identify upregulated and downregulated genes
down <- unique(dfressig[dfressig$FC < 0.5, ]$Gene.name)
up <- unique(dfressig[dfressig$FC > 2, ]$Gene.name)

# Save the upregulated and downregulated genes to text files
write.table(down, "results_down_contrl_fabry_juli_deseq.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(up, "results_up_contrl_fabry_juli_deseq.txt", sep = "\t", quote = FALSE, row.names = FALSE)
