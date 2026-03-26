#!/usr/bin/env python3

import pandas as pd

print("=== Step 7A: Check Sample Order Alignment ===")

expr_file = "data/processed/mock_expression_matrix_messy_clean.csv"
meta_file = "metadata/mock_sample_metadata_from_helper.csv"

print(f"[INFO] Loading cleaned expression matrix: {expr_file}")
expr_df = pd.read_csv(expr_file)

print(f"[INFO] Loading metadata:                 {meta_file}")
meta_df = pd.read_csv(meta_file)

# Extract sample columns from expression matrix
expr_sample_cols = [c for c in expr_df.columns if c != "gene_id"]

# Sample order from metadata
meta_sample_cols = meta_df["sample_column"].tolist()

print("\n[INFO] Expression sample order:")
for i, s in enumerate(expr_sample_cols, start=1):
    print(f" {i}. {s}")

print("\n[INFO] Metadata sample order:")
for i, s in enumerate(meta_sample_cols, start=1):
    print(f" {i}. {s}")

# Exact order check
if expr_sample_cols == meta_sample_cols:
    print("\n[OK] Sample order is aligned (exact match).")
else:
    print("\n[WARN] Sample order is NOT aligned.")
    print("[INFO] The sample sets may still match, but the order differs.")
    
    # Show position-by-position mismatches
    max_len = max(len(expr_sample_cols), len(meta_sample_cols))
    print("\n[INFO] Position-by-position comparison:")
    for i in range(max_len):
        expr_val = expr_sample_cols[i] if i < len(expr_sample_cols) else "<missing>"
        meta_val = meta_sample_cols[i] if i < len(meta_sample_cols) else "<missing>"
        status = "OK" if expr_val == meta_val else "MISMATCH"
        print(f" {i+1}. expr={expr_val} | meta={meta_val} [{status}]")

print("\n[OK] Alignment check complete.")
