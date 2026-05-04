#!/usr/bin/env Rscript

# scripts/01_data_download.R
library(fs)
library(glue)

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
  "features.tsv.gz" = glue("{geo_pref}_features.tsv.gz"),#GEO may use 'genes';NB:Seurat prefers 'features'
  "matrix.mtx.gz"   = glue("{geo_pref}_matrix.mtx.gz")
)

# Base GEO download URL structure
base_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc={geo_acc}&format=file&file={geo_file}"

message(glue("Downloading dataset {geo_acc} directly from GEO..."))

# Loop through and download each file, renaming it on the fly
for (name in names(files_to_download)) {
  
  geo_file <- files_to_download[[name]]
  download_url <- glue(base_url)
  dest_file <- path(data_matrix_dir, name)
  
  message(glue(" -> Fetching {name}..."))
  
 curl:: curl_download(
    url = download_url, 
    destfile = dest_file, 
    quiet = TRUE
  )
}

message("Download complete! Files are named correctly and ready for Seurat.")
