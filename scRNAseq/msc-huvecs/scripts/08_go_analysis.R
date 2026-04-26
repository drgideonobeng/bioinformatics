#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
library(fs)
library(clusterProfiler)
library(org.Mm.eg.db)
library(readr)
library(glue)

# 1. Setup paths
table_dir <- Sys.getenv("TABLE_DIR", unset = "results/tables")
plot_dir  <- Sys.getenv("PLOT_DIR", unset = "results/plots")

# 2. Load the FILTERED DE Results
message("Loading filtered DE results...")
de_results <- read_csv(fs::path(table_dir, "06_Cluster_Endothelial_E105_vs_E80_DE.csv"))

# 3. Isolate the Significant Genes for each timepoint
# We use a strict p-value and a Log2FC cutoff of 0.5 (approx 1.5x fold change)
up_in_E105 <- de_results %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 0.5) %>% 
  pull(gene)

up_in_E80 <- de_results %>% 
  filter(p_val_adj < 0.05 & avg_log2FC < -0.5) %>% 
  pull(gene)

message(sprintf("Found %d genes UP in E10.5 and %d genes UP in E8.0", length(up_in_E105), length(up_in_E80)))

# 4. Run Gene Ontology (Biological Process)
message("Running GO Analysis for E10.5...")
go_E105 <- enrichGO(gene          = up_in_E105,
                    OrgDb         = org.Mm.eg.db,
                    keyType       = "SYMBOL", # We are using standard gene names
                    ont           = "BP",     # BP = Biological Process
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05)

message("Running GO Analysis for E8.0...")
go_E80 <- enrichGO(gene          = up_in_E80,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "SYMBOL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

# 5. Save the raw pathway tables
dir_create(table_dir)
write_csv(as.data.frame(go_E105), fs::path(table_dir, "08_GO_Pathways_UP_E105.csv"))
write_csv(as.data.frame(go_E80), fs::path(table_dir, "08_GO_Pathways_UP_E80.csv"))

# 6. Generate DotPlots
message("Generating DotPlots...")
dir_create(plot_dir)

# We use a tryCatch just in case one of the lists doesn't have enough genes to form a pathway
if (nrow(as.data.frame(go_E105)) > 0) {
  p1 <- dotplot(go_E105, showCategory = 10) + 
        ggtitle("Top 10 Biological Processes: UP in E10.5") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
  pdf(fs::path(plot_dir, "08_GO_DotPlot_E105.pdf"), width = 8, height = 6)
  print(p1)
  dev.off()
}

if (nrow(as.data.frame(go_E80)) > 0) {
  p2 <- dotplot(go_E80, showCategory = 10) + 
        ggtitle("Top 10 Biological Processes: UP in E8.0") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        
  pdf(fs::path(plot_dir, "08_GO_DotPlot_E80.pdf"), width = 8, height = 6)
  print(p2)
  dev.off()
}

message("Success! GO tables and plots saved.")
