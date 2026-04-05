#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(survival)
  library(survminer)
  library(tidyverse)
})

cat("=== Step 15: Survival Analysis (Deep Metadata Search) ===\n")

# 1. Load Data
dds <- readRDS("results/dds_output.rds")
vsd <- vst(dds, blind = FALSE)
pc1_drivers <- read.csv("results/pc1_driver_genes.csv")
meta <- as.data.frame(colData(dds))

# 2. DEBUG: What are the clinical columns actually called?
# We'll look for anything containing 'status' or 'vital'
potential_status = grep("status|vital", colnames(meta), ignore.case = TRUE, value = TRUE)
potential_time = grep("days|time|follow", colnames(meta), ignore.case = TRUE, value = TRUE)

cat("Metadata check:\n")
cat("Found potential status columns:", paste(potential_status, collapse=", "), "\n")
cat("Found potential time columns:", paste(potential_time, collapse=", "), "\n\n")

# 3. Assign Columns (Prioritizing TCGA standard names)
# We pick the first match that isn't empty
status_col <- potential_status[1]
time_col <- potential_time[grep("death", potential_time)][1] 
if(is.na(time_col)) time_col <- potential_time[1]

if (is.na(status_col) || is.na(time_col)) {
  stop("CRITICAL ERROR: Could not find survival columns. Please check the 'Found potential...' list above.")
}

# 4. Prepare DataFrame
# Filter for Tumor patients with non-NA survival data
surv_df <- meta %>% 
  filter(condition == "Tumor") %>%
  mutate(
    event = ifelse(tolower(get(status_col)) == "dead", 1, 0),
    time = as.numeric(get(time_col))
  ) %>%
  filter(!is.na(time) & time > 0)

# 5. Link with Top PC1 Gene
top_gene_id <- pc1_drivers$ensembl_id[1]
top_gene_symbol <- pc1_drivers$symbol[1]

# Get expression for the valid patients
gene_expr <- assay(vsd)[top_gene_id, rownames(surv_df)]
surv_df$group <- ifelse(gene_expr > median(gene_expr), "High", "Low")

# 6. Fit & Plot
cat("Fitting survival model for", top_gene_symbol, "...\n")
fit <- survfit(Surv(time, event) ~ group, data = surv_df)

# 
km_plot <- ggsurvplot(fit, data = surv_df,
                      pval = TRUE, 
                      risk.table = TRUE,
                      conf.int = TRUE,
                      palette = c("#E41A1C", "#377EB8"),
                      title = paste("TCGA-COAD Survival:", top_gene_symbol),
                      xlab = "Time (Days)",
                      legend.labs = c("High Expr", "Low Expr"))

# 7. Save
pdf("plots/survival_km_plot.pdf", width = 8, height = 9, onefile = FALSE)
print(km_plot)
dev.off()

cat("[OK] Survival plot saved to plots/survival_km_plot.pdf\n")
