#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(fs)
})

message("\n--- Gene Ontology (GO) Enrichment ---")

# 1. Load the significant DEGs from the previous step
res_dir <- path("results", "DE_results")
plot_dir <- path("results", "plots")

sig_csv <- path(res_dir, "CLE9_vs_Water_sig_DEGs.csv")
if(!file_exists(sig_csv)) stop("Error: Cannot find the sig_DEGs file!")

# Read the file into the 'sig_degs' variable
message("Loading significant genes...")
sig_degs <- read_csv(sig_csv, show_col_types = FALSE)

# 2. Pull out just the Gene IDs
significant_gene_ids <- sig_degs$GeneID

# 3. Run the GO Enrichment
message("Running GO Enrichment against the Arabidopsis database...")
go_results <- enrichGO(gene          = significant_gene_ids,
                       OrgDb         = org.At.tair.db,
                       keyType       = "TAIR",  # Arabidopsis standard ID format
                       ont           = "BP",    # BP = Biological Process
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05)

# 4. Save the table of pathways
out_csv <- path(res_dir, "CLE9_GO_Biological_Processes.csv")
write_csv(as.data.frame(go_results), out_csv)
message("=> Saved GO pathway data to: ", out_csv)

# 5. Plot the top 15 most significantly enriched pathways
message("Generating GO Dot Plot...")
png(path(plot_dir, "go_enrichment_dotplot.png"), width = 900, height = 700, res = 120)

# Note: The print() wrapper is required for terminal-executed scripts!
print(dotplot(go_results, showCategory = 15, title = "Top 15 Enriched GO Terms: CLE9 vs Water"))

dev.off()

message("=> GO Enrichment complete! Plot saved to results/plots/")
