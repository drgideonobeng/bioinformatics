#!/usr/bin/env python3

import pandas as pd

print("=== Step 6C: Validate Helper-Generated Metadata ===")

expr_file = "data/processed/mock_expression_matrix_messy_clean.csv"
meta_file = "metadata/mock_sample_metadata_from_helper.csv"

print(f"[INFO] Loading cleaned expression matrix: {expr_file}")
expr_df = pd.read_csv(expr_file)

print(f"[INFO] Loading helper metadata:         {meta_file}")
meta_df = pd.read_csv(meta_file)

# Sample columns from expression matrix
expr_sample_cols = [c for c in expr_df.columns if c != "gene_id"]

# Sample columns listed in metadata
meta_sample_cols = meta_df["sample_column"].tolist()

print("\n[INFO] Sample count comparison:")
print(f" - Expression sample columns: {len(expr_sample_cols)}")
print(f" - Metadata rows:            {len(meta_sample_cols)}")

expr_set = set(expr_sample_cols)
meta_set = set(meta_sample_cols)

missing_in_meta = sorted(expr_set - meta_set)
extra_in_meta = sorted(meta_set - expr_set)

print("\n[INFO] Validation results:")
if not missing_in_meta:
    print("[OK] No cleaned expression columns are missing from helper metadata.")
else:
    print("[FAIL] Missing in metadata:")
    for x in missing_in_meta:
        print(f" - {x}")

if not extra_in_meta:
    print("[OK] No extra metadata sample columns.")
else:
    print("[FAIL] Extra in metadata:")
    for x in extra_in_meta:
        print(f" - {x}")

if meta_df["sample_id"].is_unique:
    print("[OK] sample_id values are unique.")
else:
    print("[FAIL] Duplicate sample_id values found.")

print("\n[INFO] Metadata preview:")
print(meta_df)

print("\n[OK] Helper metadata validation complete.")
