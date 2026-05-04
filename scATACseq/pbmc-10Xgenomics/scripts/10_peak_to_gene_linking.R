#!/usr/bin/env Rscript
library(Seurat)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(fs)
library(glue)
library(dplyr)
library(ggplot2)
library(future)
library(readr)

# Parallel Processing (Increase threads to 8)
plan(multicore, workers = 8)
options(future.globals.maxSize = 8 * 1024^3)

message("=========Parallel processing enabled=========")

# Fetch configurations
obj_dir    <- Sys.getenv("OBJ_DIR", unset = "results/objects")
plot_dir   <- Sys.getenv("PLOT_DIR", unset = "results/plots")
tables_dir <- Sys.getenv("TABLES_DIR", unset = "results/tables")

dir_create(tables_dir)

message("===============================================")
message(" Step 10: Unbiased Peak-to-Gene Linking (GRN)  ")
message("===============================================")

# 1. Load Object
message("Loading ChromVAR-annotated ATAC object...")
atac_obj <- readRDS(path(obj_dir, "09_atac_chromvar.rds"))

# 2. Prepare the Peaks (Calculate GC Content)
message("Calculating RegionStats (GC content & peak length)...")
DefaultAssay(atac_obj) <- 'peaks'
atac_obj <- RegionStats(atac_obj, genome = BSgenome.Hsapiens.UCSC.hg38)

# 3. Create Gene Activity Matrix (Inferred RNA)
if (!"RNA" %in% Assays(atac_obj)) {
  message("No RNA assay found. Inferring Gene Activity from ATAC fragments...")
  gene_activities <- GeneActivity(atac_obj)
  
  # Add the inferred gene activity as a new 'RNA' assay
  atac_obj[['RNA']] <- CreateAssayObject(counts = gene_activities)
  
  # Normalize the new RNA assay
  atac_obj <- NormalizeData(
    object = atac_obj,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(atac_obj$nCount_RNA)
  )
}

# 4. Dynamically Discover Top Monocyte Target Genes
message("Finding unbiased top marker genes for CD14+ Monocytes...")

DefaultAssay(atac_obj) <- 'RNA'
Idents(atac_obj) <- "predicted.id" 

mono_markers <- FindMarkers(
  object = atac_obj,
  ident.1 = "CD14+ Monocytes",
  only.pos = TRUE,
  min.pct = 0.25 # Gene must be active in at least 25% of monocytes
)

# Sort by the most statistically significant and extract the top 6 names
top_discovered_genes <- mono_markers %>%
  arrange(p_val_adj, desc(avg_log2FC)) %>%
  head(6) %>%
  rownames()

message(glue("-> Top discovered targets: {paste(top_discovered_genes, collapse=', ')}"))

# Save these discovered markers for your records
write_csv(
  x = mono_markers %>% tibble::rownames_to_column("gene"), 
  file = path(tables_dir, "10_CD14_Mono_discovered_targets.csv")
)

# 5. Perform Peak-to-Gene Linking
message(glue("Linking distal enhancers to these unbiased target genes..."))
message("(Calculating Pearson correlations across all cells - this may take a moment)")

atac_obj <- LinkPeaks(
  object = atac_obj,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = top_discovered_genes
)

# 6. Visualize the 3D Enhancer Loops (CoveragePlot)
top_gene <- top_discovered_genes[1]
message(glue("Plotting 3D chromatin loops for #1 discovered gene: {top_gene}..."))

# Switch back to the DNA layer for plotting
DefaultAssay(atac_obj) <- 'peaks' 

# Dynamically grab a valid negative control cluster from your specific dataset
all_clusters <- levels(Idents(atac_obj))
neg_control <- all_clusters[all_clusters != "CD14+ Monocytes"][1]

p1 <- CoveragePlot(
  object = atac_obj,
  region = top_gene,
  features = top_gene,
  expression.assay = "RNA",
  idents = c("CD14+ Monocytes", neg_control), # Dynamically compares against an existing cluster
  extend.upstream = 50000,
  extend.downstream = 50000
) + 
  ggtitle(glue("Unbiased Peak-to-Gene Links: {top_gene} Regulatory Landscape")) +
  theme_minimal()

plot_path <- path(plot_dir, glue("10_{top_gene}_Enhancer_Links.pdf"))
ggsave(filename = plot_path, plot = p1, width = 10, height = 8)
message(glue("-> Saved CoveragePlot to: {plot_path}"))

# 7. Save Final Object
saveRDS(atac_obj, path(obj_dir, "10_atac_linked.rds"))
message("===============================================")
message("Step 10: Unbiased GRN Linking Completed")


