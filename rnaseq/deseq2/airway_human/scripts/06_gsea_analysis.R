#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(msigdbr)
  library(fs)
})

cat("=== Step 6: Running Gene Set Enrichment Analysis (GSEA) ===\n")

# 1. Load annotated results
res <- read_csv("results/airway_de_results_annotated.csv", show_col_types = FALSE)

# 2. Prepare the Ranked List using SYMBOLS
# GSEA needs a unique list. We'll use the 'symbol' column we created in Step 4.
gene_list_df <- res %>%
  filter(!is.na(symbol)) %>%
  group_by(symbol) %>%
  slice_max(order_by = abs(log2FoldChange), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(desc(log2FoldChange))

gene_list <- gene_list_df$log2FoldChange
names(gene_list) <- gene_list_df$symbol

# 3. Pull Hallmark Gene Sets (specifically requesting Gene Symbols)
# 'gs_name' is the pathway, 'gene_symbol' is the ID
h_df <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  select(gs_name, gene_symbol)

cat(paste0("[INFO] Ranking ", length(gene_list), " unique gene symbols.\n"))

# 4. Run GSEA
cat("[INFO] Calculating enrichment scores...\n")
set.seed(42)
gsea_res <- GSEA(gene_list, 
                 TERM2GENE = h_df, 
                 pvalueCutoff = 0.05, 
                 verbose = FALSE)

# 5. Save results
dir_create("results")
write_csv(as.data.frame(gsea_res), "results/gsea_hallmark_results.csv")

# 6. Plotting
library(enrichplot)
if(nrow(as.data.frame(gsea_res)) > 0) {
  plot <- dotplot(gsea_res, showCategory=10, split=".sign") + 
    facet_grid(.~.sign) +
    theme_minimal()
  
  dir_create("figures")
  ggsave("figures/gsea_dotplot.png", plot, width = 10, height = 7)
  cat("[OK] GSEA complete! Results saved in results/ and figures/\n")
} else {
  cat("[WARN] No significant pathways found at p < 0.05. Try a higher pvalueCutoff?\n")
}
