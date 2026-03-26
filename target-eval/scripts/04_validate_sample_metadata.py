#!/usr/bin/env python3

import pandas as pd

print("=== Step 4D: Validate Sample Metadata Against Expression Matrix ===")

expr_file = "data/raw/mock_expression_matrix.csv"
meta_file = "metadata/mock_sample_metadata.csv"

print(f"[INFO] Loading expression matrix: {expr_file}")
expr_df = pd.read_csv(expr_file)

print(f"[INFO] Loading sample metadata:   {meta_file}")
meta_df = pd.read_csv(meta_file)

# Extract sample columns from expression matrix
gene_col = "gene_id"
expr_sample_cols = [col for col in expr_df.columns if col != gene_col]

# Extract sample columns listed in metadata
meta_sample_cols = meta_df["sample_column"].tolist()

print("\n[INFO] Sample count comparison:")
print(f" - Expression matrix sample columns: {len(expr_sample_cols)}")
print(f" - Metadata rows:                    {len(meta_sample_cols)}")

# Compare sets
expr_set = set(expr_sample_cols)
meta_set = set(meta_sample_cols)

missing_in_meta = sorted(expr_set - meta_set)
extra_in_meta = sorted(meta_set - expr_set)

print("\n[INFO] Validation results:")

if not missing_in_meta:
    print("[OK] No expression sample columns are missing from metadata.")
else:
    print("[FAIL] These expression sample columns are missing from metadata:")
    for s in missing_in_meta:
        print(f"  - {s}")

if not extra_in_meta:
    print("[OK] No extra sample columns found in metadata.")
else:
    print("[FAIL] These metadata sample columns are not in the expression matrix:")
    for s in extra_in_meta:
        print(f"  - {s}")

# Optional check: sample IDs unique?
if meta_df["sample_id"].is_unique:
    print("[OK] sample_id values are unique.")
else:
    print("[FAIL] Duplicate sample_id values found.")

print("\n[OK] Metadata validation complete.")
