#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(enrichplot)
  library(clusterProfiler)
  library(DOSE)
})

cat("=== Step 10: Generating Gene-Pathway Network (Cnetplot) ===\n")

# 1. Load Data
de_res <- read.csv("results/colon_cancer_annotated.csv") %>%
  filter(!is.na(symbol)) %>%
  distinct(symbol, .keep_all = TRUE)

# Named vector for coloring
genelist <- de_res$log2FoldChange
names(genelist) <- de_res$symbol

# 2. Load GSEA Results
gsea_csv <- read.csv("results/gsea_hallmark_results.csv")

# 3. Format Pathway-Gene list
# We pick the top 3 pathways to avoid the "spaghetti ball" effect
top_pathways <- gsea_csv %>%
  filter(padj < 0.05) %>%
  arrange(padj) %>%
  head(3)

# Split the genes into a list format cnetplot likes
pathway_genes <- setNames(strsplit(top_pathways$leadingEdge, ";"), top_pathways$pathway)

# 4. Generate the Plot (Simplified Arguments)
cat("-> Rendering network...\n")

# Using the most stable argument set for cnetplot
plot <- cnetplot(pathway_genes, 
                 foldChange = genelist,
                 layout = "circle",    # Replacement for circular = TRUE
                 node_label = "all") + 
  labs(title = "Gene-Concept Network: Top Hallmark Pathways",
       subtitle = "TCGA-COAD: High-Impact Genes connecting Pathways",
       color = "Log2 Fold Change")

# 5. Save the plot
if(!dir.exists("plots")) dir.create("plots")
# Use ggsave for better compatibility with ggplot objects
ggsave("plots/gene_pathway_network.png", plot, width = 12, height = 10, bg = "white")

cat("[OK] Network plot saved to plots/gene_pathway_network.png\n")
