# scripts/09_differential_expression.R
library(Seurat)
library(fs)
library(glue)

# ========== 1. PULL VARIABLES FROM CONFIG ==========
message("Loading configuration variables for Step 09...")
obj_dir   <- Sys.getenv("OBJ_DIR")
table_dir <- Sys.getenv("TABLE_DIR")

# ========== 2. LOAD ANNOTATED OBJECT ==========
message("Loading automated annotated dataset (SingleR)...")
umap_obj <- readRDS(fs::path(obj_dir, "08b_annotated.rds"))

# Fix for Seurat v5 to ensure layers can be compared mathematically
umap_obj <- JoinLayers(umap_obj)

# ========== 3. SUBSET THE CELLS ==========
message("Isolating the B cells...")
# First, we tell Seurat to look at the AI-generated labels
Idents(umap_obj) <- "automated_cell_type"

# MonacoImmuneData labels them as "B cells"
b_cells <- subset(umap_obj, idents = "B cells")

# Now switch the identity to the condition (Untreated vs Treated folders)
# From Step 07, we know your folder names are stored in orig.ident
Idents(b_cells) <- "orig.ident"

# ========== 4. RUN DIFFERENTIAL EXPRESSION ==========
message("Comparing Treated B-cells vs Untreated B-cells...")

# ident.1 = Treated
# ident.2 = Untreated
treated_vs_untreated_markers <- FindMarkers(
  b_cells, 
  ident.1 = "Treated",   
  ident.2 = "Untreated", 
  logfc.threshold = 0.25
)

# ========== 5. SAVE RESULTS ==========
message("Saving results...")
write.csv(treated_vs_untreated_markers, fs::path(table_dir, "09_treated_vs_untreated_Bcells.csv"))

message("Step 09 Complete! Differential expression test finished.")
message(glue("Total significant genes altered by the treated environment: {nrow(treated_vs_untreated_markers)}"))
