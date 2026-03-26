# scripts/00_fetch_benchmark_data.R
library(Seurat)
library(Matrix)
library(SeuratData)
library(fs)      # Modern file system operations
library(glue)    # Elegant string interpolation
library(readr)   # Fast, clean TSV writing

message(glue("Loading the Kang 'ifnb' dataset..."))
data("ifnb")

message(glue("Updating object to match your modern Seurat version..."))
ifnb <- UpdateSeuratObject(ifnb)

message(glue("Extracting raw count matrices (the Seurat v5 way!)..."))
counts <- GetAssayData(ifnb, assay = "RNA", layer = "counts")

# Split the cell barcodes into Untreated (CTRL) and Stimulated (STIM)
ctrl_barcodes <- rownames(ifnb@meta.data[ifnb@meta.data$stim == "CTRL", ])
stim_barcodes <- rownames(ifnb@meta.data[ifnb@meta.data$stim == "STIM", ])

ctrl_counts <- counts[, ctrl_barcodes]
stim_counts <- counts[, stim_barcodes]

message(glue("Formatting and writing raw 10x matrices to your folders..."))

write_10x <- function(mat, out_dir) {
  # 1. Safely create the directory using fs
  fs::dir_create(out_dir)
  
  # 2. Define clean paths using fs::path
  mat_file  <- fs::path(out_dir, "matrix.mtx")
  feat_file <- fs::path(out_dir, "features.tsv")
  bar_file  <- fs::path(out_dir, "barcodes.tsv")
  
  # 3. Write files (using readr for the dataframes to avoid base R clutter)
  Matrix::writeMM(mat, file = mat_file)
  
  features_df <- data.frame(
    V1 = rownames(mat), 
    V2 = rownames(mat), 
    V3 = "Gene Expression"
  )
  readr::write_tsv(features_df, file = feat_file, col_names = FALSE)
  
  barcodes_df <- data.frame(V1 = colnames(mat))
  readr::write_tsv(barcodes_df, file = bar_file, col_names = FALSE)
              
  # 4. Automatically compress the files
  message(glue("Compressing files in {out_dir} to .gz format..."))
  system2("gzip", args = c("-f", mat_file))
  system2("gzip", args = c("-f", feat_file))
  system2("gzip", args = c("-f", bar_file))
}

write_10x(ctrl_counts, "data/untreated")
write_10x(stim_counts, "data/treated")

message(glue("Success! Data folders are primed, compressed, and ready."))
