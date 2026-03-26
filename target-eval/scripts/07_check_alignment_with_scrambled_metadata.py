#!/usr/bin/env python3

import pandas as pd

print("=== Step 7C: Check Alignment with Scrambled Metadata ===")

expr_file = "data/processed/mock_expression_matrix_messy_clean.csv"
meta_file = "metadata/mock_sample_metadata_scrambled.csv"

print(f"[INFO] Loading cleaned expression matrix: {expr_file}")
expr_df = pd.read_csv(expr_file)

print(f"[INFO] Loading scrambled metadata:       {meta_file}")
meta_df = pd.read_csv(meta_file)

# Get sample order from each file
expr_sample_cols = [c for c in expr_df.columns if c != "gene_id"]
meta_sample_cols = meta_df["sample_column"].tolist()

print("\n[INFO] Expression sample order:")
for i, sample_name in enumerate(expr_sample_cols, start=1):
    print(f" {i}. {sample_name}")

print("\n[INFO] Scrambled metadata order:")
for i, sample_name in enumerate(meta_sample_cols, start=1):
    print(f" {i}. {sample_name}")

# Check exact order match
print("\n[INFO] Comparing order...")
if expr_sample_cols == meta_sample_cols:
    print("[OK] Sample order is aligned (exact match).")
else:
    print("[WARN] Sample order is NOT aligned.")
    print("[INFO] Showing position-by-position comparison:")

    max_len = max(len(expr_sample_cols), len(meta_sample_cols))
    for i in range(max_len):
        expr_name = expr_sample_cols[i] if i < len(expr_sample_cols) else "<missing>"
        meta_name = meta_sample_cols[i] if i < len(meta_sample_cols) else "<missing>"

        if expr_name == meta_name:
            status = "OK"
        else:
            status = "MISMATCH"

        print(f" {i+1}. expr={expr_name} | meta={meta_name} [{status}]")

print("\n[OK] Alignment check complete.")
