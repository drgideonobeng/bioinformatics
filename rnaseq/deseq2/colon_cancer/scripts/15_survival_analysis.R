#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(survival)
  library(survminer)
  library(tidyverse)
})

cat("=== Step 15: Survival Analysis (PC1 Top Gene) ===\n")

# 1. Load Data
dds <- readRDS("results/dds_output.rds")
vsd <- vst(dds, blind = FALSE)
pc1_drivers <- read.csv("results/pc1_driver_genes.csv")

# 2. Pick the #1 Gene from your PC1 analysis
top_gene_symbol <- pc1_drivers$symbol[1]
top_gene_id <- pc1_drivers$ensembl_id[1]

# 3. Prepare Clinical Data from Metadata
meta <- as.data.frame(colData(dds))

# Standardize Vital Status (1 = Event/Death, 0 = Censored/Alive)
meta$event <- ifelse(tolower(meta$tcga.gdc_cases.demographic.vital_status) == "dead", 1, 0)

# Calculate Overall Survival Time
# If dead, use days_to_death; if alive, use days_to_last_follow_up
meta$time <- as.numeric(ifelse(meta$event == 1, 
                               meta$tcga.gdc_cases.demographic.days_to_death, 
                               meta$tcga.gdc_cases.diagnoses.days_to_last_follow_up))

# Filter out samples with missing time or status (common in TCGA)
surv_df <- meta %>% 
  filter(!is.na(time) & time > 0) %>%
  filter(condition == "Tumor") # Only perform survival on tumor patients

# 4. Add Gene Expression Groups (High vs Low)
# Pull expression for our top gene and split at the median
gene_expr <- assay(vsd)[top_gene_id, rownames(surv_df)]
surv_df$expression_group <- ifelse(gene_expr > median(gene_expr), "High", "Low")

# 5. Fit Survival Model
fit <- survfit(Surv(time, event) ~ expression_group, data = surv_df)

# 6. Generate Kaplan-Meier Plot
km_plot <- ggsurvplot(fit, data = surv_df,
                      pval = TRUE,             # Show p-value
                      risk.table = TRUE,       # Show number at risk
                      conf.int = TRUE,         # Show confidence intervals
                      palette = c("#E41A1C", "#377EB8"),
                      title = paste("Overall Survival:", top_gene_symbol),
                      xlab = "Time (Days)",
                      legend.labs = c("High Expression", "Low Expression"))

# 7. Save
pdf("plots/survival_km_plot.pdf", width = 8, height = 9, onefile = FALSE)
print(km_plot)
dev.off()

cat(paste("[OK] Survival plot for", top_gene_symbol, "saved to plots/survival_km_plot.pdf\n"))
