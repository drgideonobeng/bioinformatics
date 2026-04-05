#!/usr/bin/env Rscript

cat("=== RNA-seq Project: Environment Check ===\n\n")

# Print R version info
cat("[INFO] R version:\n")
cat(R.version.string, "\n\n")

# Packages we want to test
packages_to_check <- c(
  "DESeq2",
  "apeglm",
  "tidyverse",
  "readr",
  "pheatmap",
  "optparse",
  "rmarkdown",
  "knitr"
)

# Function to test one package
check_package <- function(pkg_name) {
  cat("[INFO] Checking package:", pkg_name, "\n")
  
  # Try to load the package quietly
  loaded_ok <- suppressPackageStartupMessages(
    require(pkg_name, character.only = TRUE, quietly = TRUE)
  )
  
  if (loaded_ok) {
    pkg_version <- as.character(packageVersion(pkg_name))
    cat("[OK]   Loaded", pkg_name, "version", pkg_version, "\n\n")
  } else {
    cat("[FAIL] Could not load", pkg_name, "\n\n")
  }
}

# Run checks
for (pkg in packages_to_check) {
  check_package(pkg)
}

cat("[OK] Environment check complete.\n")
