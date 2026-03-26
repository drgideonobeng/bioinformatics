#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
})

cat("=== Step 5: Gene Set Enrichment Analysis (GSEA) ===\n")

# 1. Load annotated results
input_file <- "results/colon_cancer_annotated.csv"
if (!file.exists(input_file)) stop("Missing input! Run Step 3 first.")

# Use 'stat' for ranking (the Wald statistic)
res_df <- read.csv(input_file) %>% 
  filter(!is.na(symbol) & !is.na(stat)) %>%
  distinct(symbol, .keep_all = TRUE)

# 2. Prepare the Ranked List
ranks <- res_df$stat
names(ranks) <- res_df$symbol
ranks <- sort(ranks, decreasing = TRUE)

# 3. Load Hallmark Gene Sets (Using the stable 'H' category)
cat("-> Loading MSigDB Hallmark gene sets...\n")
# Hallmark sets are always category 'H'
h_df <- msigdbr(species = "human", category = "H") 
h_list <- split(x = h_df$gene_symbol, f = h_df$gs_name)

# 4. Run fgsea
cat("-> Running fgsea...\n")
set.seed(42) # For reproducibility
fgsea_res <- fgsea(pathways = h_list, 
                  stats = ranks,
                  minSize = 15,
                  maxSize = 500)

# 5. Format results for CSV (Flattening leadingEdge list)
cat("-> Formatting results for CSV export...\n")
fgsea_res_txt <- fgsea_res %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

# 6. Save Results
if(!dir.exists("results")) dir.create("results")
write.csv(fgsea_res_txt, "results/gsea_hallmark_results.csv", row.names = FALSE)

# 7. Quick Summary Output
cat("\nTop 5 Upregulated Pathways (Positive NES):\n")
print(head(fgsea_res_txt %>% filter(NES > 0) %>% select(pathway, NES, padj), 5))

cat("\nTop 5 Downregulated Pathways (Negative NES):\n")
print(head(fgsea_res_txt %>% filter(NES < 0) %>% select(pathway, NES, padj), 5))

cat("\n[OK] GSEA complete! Results in results/gsea_hallmark_results.csv\n")
