#!/usr/bin/env Rscript

# scripts/01_data_download.R
library(fs)
library(glue)

# Increase timeout for large GEO files on Mac
options(timeout = max(1000, getOption("timeout")))

# Pull variables from config
geo_acc <- Sys.getenv("GEO_ACCESSION")
geo_pref <- Sys.getenv("GEO_PREFIX")
data_matrix_dir <- Sys.getenv("DATA_MATRIX_DIR")

if (geo_acc == "") stop("GEO_ACCESSION is empty. Check config.sh")

# Ensure the exact matrix directory exists
dir_create(data_matrix_dir)

# Define the three files we need from GEO and what Seurat expects them to be called
files_to_download <- list(
  "barcodes.tsv.gz" = glue("{geo_pref}_barcodes.tsv.gz"),
  "features.tsv.gz" = glue("{geo_pref}_genes.tsv.gz"), # GEO often uses 'genes', Seurat likes 'features'
  "matrix.mtx.gz"   = glue("{geo_pref}_matrix.mtx.gz")
)

# Base GEO download URL structure
base_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc={geo_acc}&format=file&file={geo_file}"

message(glue("Downloading dataset {geo_acc} directly from GEO..."))

# Loop through and download each file, renaming it on the fly
for (seurat_name in names(files_to_download)) {
  
  geo_file <- files_to_download[[seurat_name]]
  download_url <- glue(base_url)
  dest_file <- path(data_matrix_dir, seurat_name)
  
  message(glue(" -> Fetching {seurat_name}..."))
  
  download.file(
    url = download_url, 
    destfile = dest_file, 
    mode = "wb", # Important for downloading compressed files on Mac/Windows
    quiet = TRUE
  )
}

message("Download complete! Files are named correctly and ready for Seurat.")
