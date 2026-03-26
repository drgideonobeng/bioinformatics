#!/usr/bin/env python3

import pandas as pd

print("=== Step 5B: Inspect Messy Expression Matrix ===")

file_path = "data/raw/mock_expression_matrix_messy.csv"
print(f"[INFO] Loading file: {file_path}")
df = pd.read_csv(file_path)

print("\n[INFO] Shape (rows, columns):", df.shape)

print("\n[INFO] Column names (repr):")
for col in df.columns:
    print(f" - {repr(col)}")

# Gene ID checks
print("\n[INFO] Gene IDs (repr, first 5):")
for g in df["gene_id"].head():
    print(f" - {repr(g)}")

# Missing values per column
print("\n[INFO] Missing values per column:")
print(df.isna().sum())

# Identify sample columns (everything except gene_id)
gene_col = "gene_id"
sample_cols = [c for c in df.columns if c != gene_col]

print("\n[INFO] Sample columns detected:")
for c in sample_cols:
    print(f" - {repr(c)}")

# Check sample label format (expected pattern for this exercise)
# Expected: sample_###_Group
print("\n[INFO] Sample label format check (expected: sample_###_Group)")
for c in sample_cols:
    parts = c.strip().split("_")  # strip first so trailing spaces don't hide format
    if len(parts) != 3:
        print(f"[WARN] Bad format: {repr(c)} -> split parts = {parts}")
    else:
        print(f"[OK]   Format looks valid: {repr(c)} -> {parts}")

# Check group naming consistency (Tumor/Normal capitalization)
print("\n[INFO] Group label consistency check")
for c in sample_cols:
    parts = c.strip().split("_")
    if len(parts) == 3:
        group = parts[2]
        if group not in {"Tumor", "Normal"}:
            print(f"[WARN] Non-standard group label in {repr(c)}: {repr(group)}")
        else:
            print(f"[OK]   Group label in {repr(c)}: {group}")

print("\n[OK] Messy matrix inspection complete.")
