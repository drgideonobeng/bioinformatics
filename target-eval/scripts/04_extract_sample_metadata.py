#!/usr/bin/env python3

import pandas as pd

print("=== Step 4C: Extract Sample Metadata from Column Names ===")

# Load expression matrix
expr_file = "data/raw/mock_expression_matrix.csv"
print(f"[INFO] Loading expression matrix: {expr_file}")
df = pd.read_csv(expr_file)

# Identify sample columns (everything except gene_id in this simple example)
gene_col = "gene_id"
sample_cols = [col for col in df.columns if col != gene_col]

print("\n[INFO] Sample columns found:")
for col in sample_cols:
    print(f" - {col}")

# Parse sample metadata from column names
# Expected format: sample_001_Tumor
metadata_rows = []

for col in sample_cols:
    parts = col.split("_")
    
    # For this mock format:
    # parts = ["sample", "001", "Tumor"]
    sample_prefix = parts[0]      # "sample"
    sample_number = parts[1]      # "001"
    group = parts[2]              # "Tumor" or "Normal"
    
    sample_id = f"{sample_prefix}_{sample_number}"
    
    metadata_rows.append({
        "sample_column": col,
        "sample_id": sample_id,
        "group": group
    })

metadata_df = pd.DataFrame(metadata_rows)

# Save metadata table
out_file = "metadata/mock_sample_metadata.csv"
metadata_df.to_csv(out_file, index=False)

print(f"\n[OK] Sample metadata saved to: {out_file}")
print("\n[INFO] Sample metadata preview:")
print(metadata_df)

print("\n[OK] Sample metadata extraction complete.")
