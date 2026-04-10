#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(fs)
})

message("\n--- Gene Ontology (GO) Enrichment ---")

out_dir <- Sys.getenv("OUT_DIR")
res_dir <- path(out_dir, "DE_results")
plot_dir <- path(out_dir, "plots")

# Dynamically find the significant DEGs file
sig_csv <- dir_ls(res_dir, glob = "*_sig_DEGs.csv")[1]
if(is.na(sig_csv)) stop("Error: Cannot find any *sig_DEGs.csv file!")

sig_degs <- read_csv(sig_csv, show_col_types = FALSE)
significant_gene_ids <- sig_degs$GeneID

message("Running GO Enrichment...")
go_results <- enrichGO(gene = significant_gene_ids, OrgDb = org.At.tair.db, keyType = "TAIR", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

out_csv <- path(res_dir, "GO_Biological_Processes.csv")
write_csv(as.data.frame(go_results), out_csv)

message("Generating GO Dot Plot...")
png(path(plot_dir, "go_enrichment_dotplot.png"), width = 900, height = 700, res = 120)
print(dotplot(go_results, showCategory = 15, title = "Top 15 Enriched GO Terms"))
dev.off()
