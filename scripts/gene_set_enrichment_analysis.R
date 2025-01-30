# ---------------------------
# Script Name: gene_set_enrichment_analysis.R
# Description: This script performs gene set enrichment analysis (GSEA) and pathway clustering 
# for differential gene expression (DEG) results using the DESeq2 output. 
# The script generates pathway enrichment visualizations,  pathway clustering results, 
# and prepares data for downstream analysis, such as heatmap creation.
# Author: Dr. Juliana E. Arcila Galvis
# Date: 2024-11-1
# Inputs: 
# - DESeq2 results: "Deseq_results_significant.csv", "Deseq_results.csv"
# - GO term annotations: "mart_export_GRCh38.p14.txt"
# Outputs: 
# - Enrichment results: "enrichBP.rds", "df_enrichBP.rds"
# - Pathway clusters: "clusters.rds"
# - Gene data for heatmap: "heatmap_genes.rds"
# ---------------------------

# Loading required libraries
library(aPEAR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(dplyr)

# Reading in DESeq2 results and filtering for significant genes
Deseq_results_significant2 <- read.delim("Deseq_results_significant.csv", sep = ",")
Deseq_results <- read.delim("Deseq_results.csv", sep = ",")
Deseq_fresults <- Deseq_results %>% filter(padj < 0.05 | is.na(padj) & abs(log2FoldChange) > 1)

# Reading GO terms file
Goterms <- read.delim("mart_export_GRCh38.p14.txt")

# Filtering DESeq results for significant genes and sorting them
Deseq_results_significant <- Deseq_fresults %>% arrange(-stat)
Deseq_no_NA <- Deseq_results_significant %>% filter(Entrez != "NA")

# Preparing gene list for GSEA
gene_list <- Deseq_no_NA$stat
names(gene_list) <- Deseq_no_NA$Entrez

# Perform gene set enrichment analysis (GSEA) for Biological Process (BP) terms
set.seed(892323)
enrichBP <- gseGO(gene_list, OrgDb = org.Hs.eg.db, ont = 'BP', minGSSize = 10, maxGSSize = 200)

# Convert results to a dataframe for easy manipulation
df_enrichBP <- data.frame(enrichBP@result)

# Saving results as RDS files for later use
saveRDS(enrichBP, "enrichBP.rds")
saveRDS(df_enrichBP, "df_enrichBP.rds")

# Performing pathway clustering using aPEAR
clusters <- aPEAR::findPathClusters(enrichBP@result, methods = aPEAR.methods, minClusterSize = 3)

# Saving pathway clusters as an RDS file
saveRDS(clusters, "clusters.rds")

# Plotting the pathway clusters
aPEAR::plotPathClusters(
  enrichment = enrichBP@result,
  sim = clusters$similarity,
  clusters = clusters$clusters,
  fontSize = 4,
  outerCutoff = 1, # Decrease cutoff to show connections between clusters
  drawEllipses = TRUE
)

# Extracting genes associated with specific pathways for further analysis
get_genes <- function(GOtermname) {
  listofinterest <- clusters$clusters[clusters$clusters$Cluster == GOtermname, ]$Pathway
  Goofinterest <- Goterms[Goterms$GO.term.name %in% listofinterest, "GO.term.accession"] %>% unique
  genesofinterest <- Goterms[Goterms$GO.term.name %in% listofinterest, "NCBI.gene..formerly.Entrezgene..ID"] %>% unique
  genesforheatmap <- Deseq_no_NA[Deseq_no_NA$Entrez %in% names(gene_list)[names(gene_list) %in% genesofinterest], ]
  return(genesforheatmap)
}

# Example of extracting genes for various biological processes
genesforheatmap <- get_genes("carboxylic acid metabolic process")
carboxilic_acid_genes <- genesforheatmap[genesforheatmap$lfcSE < 0.5, ]
carboxilic_acid_genes$cluster_GO <- "carboxylic acid metabolic process"

genesforheatmap <- get_genes("inorganic ion homeostasis")
ion_genes <- genesforheatmap[genesforheatmap$lfcSE < 0.5, ]
ion_genes$cluster_GO <- "inorganic ion homeostasis"

genesforheatmap <- get_genes("cellular response to growth factor stimulus")
Path1_genes <- genesforheatmap[genesforheatmap$lfcSE < 0.5, ]
Path1_genes$cluster_GO <- "cellular response to growth factor stimulus"

# Combine all the significant genes across selected pathways
heatmap_genes <- rbind(ion_genes, carboxilic_acid_genes, Path1_genes)
saveRDS(heatmap_genes, "heatmap_genes.rds")

# You can add more pathways and genes as needed here

# End of script
