#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)
  library(org.At.tair.db) # The specific database for Arabidopsis
  library(fs)
})

message("--- Functional Gene Annotation ---")

# Capture command-line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: No input CSV file provided. Usage: Rscript 08_annotate_genes.R <path_to_csv>", call. = FALSE)
}

# 1. Load the significant results dynamically from the command line
results_file <- args[1]
if(!file_exists(results_file)) stop(paste("Error: Cannot find the DE results CSV at:", results_file))

message(paste("Loading DE results from:", results_file))
res_tidy <- read_csv(results_file, show_col_types = FALSE)

# 2. Map the IDs to biological names
message("Querying org.At.tair.db for gene symbols and descriptions...")

res_annotated <- res_tidy %>%
  mutate(
    # Get Gene Symbols (Short names)
    symbol = mapIds(org.At.tair.db,
                    keys = GeneID,          
                    column = "SYMBOL",
                    keytype = "TAIR",
                    multiVals = "first"),
    
    # Get Gene Names (Full biological descriptions)
    description = mapIds(org.At.tair.db,
                         keys = GeneID,    
                         column = "GENENAME",
                         keytype = "TAIR",
                         multiVals = "first")
  ) %>%
  # Reorder columns
  relocate(symbol, description, .after = GeneID)

# Extract the base name of the input file to dynamically name outputs
# e.g., "results/DE_results_ConditionX.csv" becomes "DE_results_ConditionX"
base_name <- path_ext_remove(path_file(results_file))

# 3. Save the final, biologist-friendly CSV
out_file <- path("results", paste0(base_name, "_Annotated.csv"))
message("Saving full annotated results table...")
write_csv(res_annotated, out_file)
message("=> Saved full annotated results to: ", out_file)

# 4. Extract and save the Top 50 most significant genes
message("\nExtracting the Top 50 significant DE genes for quick reference...")

top_50_genes <- res_annotated %>%
  filter(padj < 0.05) %>%       # Keep only statistically significant results
  arrange(padj) %>%             # Sort by the smallest p-value first
  slice_head(n = 50)            # Snip off exactly the top 50 rows

# Save the Top 50 to its own dynamic file
top50_file <- path("results", paste0(base_name, "_Top50_Annotated.csv"))
write_csv(top_50_genes, top50_file)
message("=> Saved Top 50 snapshot to: ", top50_file)

# Print a quick terminal preview of the top 10 (just the most important columns)
message("\n--- Top 10 DE Genes ---")
top_50_genes %>% 
  dplyr::select(GeneID, symbol, log2FoldChange, padj) %>% 
  head(10) %>% 
  print()
message("----------------------------------\n")

message("======================================================")
message("Pipeline Expansion Complete! Data is ready for biological interpretation.")
message("======================================================")
