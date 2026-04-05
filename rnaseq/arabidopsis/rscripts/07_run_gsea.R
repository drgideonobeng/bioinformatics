#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.At.tair.db)
  library(enrichplot) # Required for advanced GSEA visualizations
  library(fs)
})

message("\n--- Gene Set Enrichment Analysis (GSEA) ---")

# 1. Load the UNFILTERED genome-wide results
res_dir <- path("results", "DE_results")
plot_dir <- path("results", "plots")

all_csv <- path(res_dir, "CLE9_vs_Water_all_genes.csv")
if(!file_exists(all_csv)) stop("Error: Cannot find the all_genes.csv file!")

message("Loading unfiltered genome-wide results...")
all_genes_df <- read_csv(all_csv, show_col_types = FALSE)

# 2. Prepare the Ranked List (The most critical step for GSEA!)
message("Ranking genes by Log2 Fold Change...")

# Drop any rows where DESeq2 couldn't calculate an LFC, then sort descending
ranked_df <- all_genes_df %>%
  drop_na(log2FoldChange, GeneID) %>%
  arrange(desc(log2FoldChange))

# Extract the Log2 Fold Changes into a standalone list
gene_list <- ranked_df$log2FoldChange

# Attach the TAIR Gene IDs as the "names" of that list
names(gene_list) <- ranked_df$GeneID

# 3. Run GSEA
message("Running GSEA against the Arabidopsis GO database...")
# We set a seed because GSEA uses random permutations. This ensures 
# you get the exact same P-values if you run the script twice.
set.seed(42) 

gsea_results <- gseGO(geneList     = gene_list,
                      OrgDb        = org.At.tair.db,
                      keyType      = "TAIR",
                      ont          = "BP", # Biological Process
                      minGSSize    = 10,   # Ignore tiny pathways
                      maxGSSize    = 500,  # Ignore overly broad pathways
                      pvalueCutoff = 0.05,
                      verbose      = FALSE)

# 4. Save the table of pathways
out_csv <- path(res_dir, "CLE9_GSEA_Biological_Processes.csv")
write_csv(as.data.frame(gsea_results), out_csv)
message("=> Saved GSEA pathway data to: ", out_csv)

# 5. Visualization: Split Dot Plot
message("Generating GSEA Plots...")

png(path(plot_dir, "gsea_split_dotplot.png"), width = 900, height = 700, res = 120)
# This splits the plot into "Activated" (upregulated) and "Suppressed" (downregulated) pathways
print(dotplot(gsea_results, showCategory = 10, split = ".sign", title = "GSEA: Activated vs Suppressed Pathways") + facet_grid(.~.sign))
dev.off()

# 6. Visualization: Ridgeplot
png(path(plot_dir, "gsea_ridgeplot.png"), width = 900, height = 700, res = 120)
print(ridgeplot(gsea_results, showCategory = 10) + ggtitle("GSEA Ridgeplot: CLE9 vs Water"))
dev.off()

message("=> GSEA complete! Plots saved to results/plots/")
