#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(tidyverse)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(fs)
  library(glue)
})

message("--- Functional Gene Annotation ---")

out_dir <- Sys.getenv("OUT_DIR")
res_dir <- path(out_dir, "DE_results")

# Automatically find only the DE results files (ignoring GO/GSEA and already annotated files)
csv_files <- dir_ls(res_dir, regexp = ".*(all_genes|sig_DEGs)\\.csv$")

if(length(csv_files) == 0) stop("No DE result CSV files found to annotate!")

for (results_file in csv_files) {
    # Extra safety check: Skip if it's already an annotated file
    if(grepl("Annotated", results_file)) next
    
    message(paste("Annotating:", path_file(results_file)))
    res_tidy <- read_csv(results_file, show_col_types = FALSE)

    res_annotated <- res_tidy %>%
      mutate(GeneID_clean = sub("\\..*", "", GeneID)) %>%
      mutate(
        symbol = mapIds(org.Hs.eg.db, keys = GeneID_clean, column = "SYMBOL", keytype = "ENSEMBLTRANS", multiVals = "first"),
        description = mapIds(org.Hs.eg.db, keys = GeneID_clean, column = "GENENAME", keytype = "ENSEMBLTRANS", multiVals = "first")
      ) %>%
      dplyr::select(-GeneID_clean) %>%
      relocate(symbol, description, .after = GeneID)

    base_name <- path_ext_remove(path_file(results_file))
    out_file <- path(res_dir, paste0(base_name, "_Annotated.csv"))
    write_csv(res_annotated, out_file)
}

