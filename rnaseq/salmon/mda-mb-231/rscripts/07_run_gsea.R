#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(fs)
})

message("\n--- Gene Set Enrichment Analysis (GSEA) ---")

out_dir <- Sys.getenv("OUT_DIR")
res_dir <- path(out_dir, "DE_results")
plot_dir <- path(out_dir, "plots")

# Dynamically find the full genome list
all_csv <- dir_ls(res_dir, glob = "*_all_genes.csv")[1]
if(is.na(all_csv)) stop("Error: Cannot find any *all_genes.csv file!")

all_genes_df <- read_csv(all_csv, show_col_types = FALSE)

ranked_df <- all_genes_df %>% drop_na(log2FoldChange, GeneID) %>% arrange(desc(log2FoldChange))
gene_list <- ranked_df$log2FoldChange

# Strip version numbers from the Ensembl IDs before assigning them as names
names(gene_list) <- sub("\\..*", "", ranked_df$GeneID)

set.seed(42) 
# Note: Relaxed pvalueCutoff to 0.20 to catch more exploratory pathways
gsea_results <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = "ENSEMBLTRANS", ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.20, verbose = FALSE)

# Convert results to a dataframe to save and check
gsea_df <- as.data.frame(gsea_results)
out_csv <- path(res_dir, "GSEA_Biological_Processes.csv")
write_csv(gsea_df, out_csv)

# Graceful check: Only plot if there are actual results!
if (nrow(gsea_df) == 0) {
  message("⚠️ No significant GSEA terms found at the current p-value cutoff. Skipping plots to prevent crash.")
} else {
  message("Generating GSEA Plots...")
  png(path(plot_dir, "gsea_split_dotplot.png"), width = 900, height = 700, res = 120)
  print(dotplot(gsea_results, showCategory = 10, split = ".sign", title = "GSEA: Activated vs Suppressed") + facet_grid(.~.sign))
  dev.off()

  png(path(plot_dir, "gsea_ridgeplot.png"), width = 900, height = 700, res = 120)
  print(ridgeplot(gsea_results, showCategory = 10) + ggtitle("GSEA Ridgeplot"))
  dev.off()
}
