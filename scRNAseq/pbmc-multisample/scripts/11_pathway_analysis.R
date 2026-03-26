# scripts/11_pathway_analysis.R
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)
library(fs)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 11...")
table_dir <- Sys.getenv("TABLE_DIR")
plot_dir  <- Sys.getenv("PLOT_DIR")

# ========== 2. LOAD DIFFERENTIAL GENES ==========
message("Loading the volcano plot data...")
de_results <- read.csv(fs::path(table_dir, "09_treated_vs_untreated_Bcells.csv"), row.names = 1)
de_results$gene <- rownames(de_results)

# ========== 3. ISOLATE UPREGULATED GENES ==========
message("Filtering for genes turned ON by the treated environment...")
# Grab genes with a log2FC > 1 and significant p-value
upregulated_genes <- de_results %>% 
  filter(avg_log2FC > 1 & p_val_adj < 0.05) %>% 
  pull(gene)

message(paste("Found", length(upregulated_genes), "upregulated genes. Running pathway analysis..."))

# ========== 4. RUN GENE ONTOLOGY (GO) ENRICHMENT ==========
message("Querying the Human Gene Ontology Database (This takes a moment!)...")
go_results <- enrichGO(
  gene          = upregulated_genes,
  OrgDb         = org.Hs.eg.db,      # The human gene database
  keyType       = "SYMBOL",          # We are using standard gene names (like ISG15)
  ont           = "BP",              # BP = Biological Process
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# ========== 5. PLOT THE PATHWAYS ==========
message("Generating Pathway DotPlot...")
pdf(fs::path(plot_dir, "11_Bcell_Upregulated_Pathways.pdf"), width = 10, height = 7)

# We plot the top 15 most significantly activated biological pathways
pathway_plot <- dotplot(go_results, showCategory = 15) + 
  ggtitle("Top Biological Pathways Activated in Treated B-cells") +
  theme(plot.title = element_text(face = "bold", size = 14))

print(pathway_plot)
dev.off()

# Save the raw pathway data to a CSV for your records
write.csv(as.data.frame(go_results), fs::path(table_dir, "11_Bcell_Pathway_Results.csv"))

message("Step 11 Complete! Pipeline officially finished.")

