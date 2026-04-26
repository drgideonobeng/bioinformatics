#!/usr/bin/env Rscript

# ========== 0. LOAD PACKAGES ==========
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db) # Mouse database
  library(tidyverse)
  library(fs)
  library(glue)
})

# ========== 1. SETUP PATHS ==========
table_dir <- Sys.getenv("TABLE_DIR", unset = "results/tables")
plot_dir  <- Sys.getenv("PLOT_DIR", unset = "results/plots")

input_file <- path(table_dir, "07_pseudotime_correlated_genes.csv")

if (!file_exists(input_file)) {
  stop(glue("Could not find {input_file}. Did you run script 09?"))
}

# ========== 2. LOAD & SPLIT GENES ==========
message("Loading pseudotime-correlated genes...")
cor_df <- read_csv(input_file, show_col_types = FALSE)

# Positive Rho = Genes increasing over time (Maturation)
# Negative Rho = Genes decreasing over time (Progenitor state)
genes_up <- cor_df %>% filter(Spearman_Rho > 0.3) %>% pull(Gene)
genes_dn <- cor_df %>% filter(Spearman_Rho < -0.3) %>% pull(Gene)

message(glue("Found {length(genes_up)} genes increasing and {length(genes_dn)} genes decreasing."))

# ========== 3. TRANSLATE TO ENTREZ IDs FOR KEGG ==========
message("Translating Gene Symbols to Entrez IDs for KEGG...")
entrez_up <- bitr(genes_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID
entrez_dn <- bitr(genes_dn, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")$ENTREZID

# ========== 4. RUN ENRICHMENT (GO & KEGG) ==========
message("Running GO and KEGG enrichments...")

# --- GO Biological Process ---
ego_up <- enrichGO(gene = genes_up, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

ego_dn <- enrichGO(gene = genes_dn, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', 
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)

# --- KEGG Pathways ---
kegg_up <- enrichKEGG(gene = entrez_up, organism = 'mmu', pvalueCutoff = 0.05)
kegg_dn <- enrichKEGG(gene = entrez_dn, organism = 'mmu', pvalueCutoff = 0.05)

# Convert KEGG IDs back to readable gene symbols for the final output tables
if (!is.null(kegg_up)) kegg_up <- setReadable(kegg_up, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
if (!is.null(kegg_dn)) kegg_dn <- setReadable(kegg_dn, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

# ========== 5. VISUALIZATION ==========
dir_create(plot_dir)

# --- Plot GO ---
pdf_go <- path(plot_dir, "10_GO_Enrichment.pdf")
message(glue("Saving GO dotplots to {pdf_go}"))
pdf(pdf_go, width = 10, height = 7)
if(!is.null(ego_up)) print(dotplot(ego_up, showCategory = 15) + ggtitle("GO: Biological Processes Increasing (Maturation)"))
if(!is.null(ego_dn)) print(dotplot(ego_dn, showCategory = 15) + ggtitle("GO: Biological Processes Decreasing (Progenitors)"))
invisible(dev.off())

# --- Plot KEGG ---
pdf_kegg <- path(plot_dir, "10_KEGG_Enrichment.pdf")
message(glue("Saving KEGG dotplots to {pdf_kegg}"))
pdf(pdf_kegg, width = 10, height = 7)
if(!is.null(kegg_up)) print(dotplot(kegg_up, showCategory = 15) + ggtitle("KEGG: Pathways Increasing (Maturation)"))
if(!is.null(kegg_dn)) print(dotplot(kegg_dn, showCategory = 15) + ggtitle("KEGG: Pathways Decreasing (Progenitors)"))
invisible(dev.off())

# ========== 6. EXPORT TABLES ==========
dir_create(table_dir)
if(!is.null(ego_up)) write_csv(as.data.frame(ego_up), path(table_dir, "10_GO_Maturation_UP.csv"))
if(!is.null(ego_dn)) write_csv(as.data.frame(ego_dn), path(table_dir, "10_GO_Progenitor_DN.csv"))
if(!is.null(kegg_up)) write_csv(as.data.frame(kegg_up), path(table_dir, "10_KEGG_Maturation_UP.csv"))
if(!is.null(kegg_dn)) write_csv(as.data.frame(kegg_dn), path(table_dir, "10_KEGG_Progenitor_DN.csv"))

message("Step 10 Complete! GO and KEGG analyses finished.")
