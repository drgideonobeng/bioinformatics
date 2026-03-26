# scripts/01_data_download.R

# Get ROOT_DIR from environment
root_dir <- Sys.getenv("ROOT_DIR")
if (root_dir == "") root_dir <- getwd() # Fallback

data_dir <- file.path(root_dir, "data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# URL for 10x Genomics PBMC 3k dataset
url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
dest_file <- file.path(data_dir, "pbmc3k_filtered_gene_bc_matrices.tar.gz")

message("Downloading PBMC 3k dataset...")
download.file(url, destfile = dest_file)

message("Extracting dataset...")
untar(dest_file, exdir = data_dir)

message("Download and extraction complete. Data is in: ", data_dir)
