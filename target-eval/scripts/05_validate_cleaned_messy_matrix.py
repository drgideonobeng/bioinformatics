#!/usr/bin/env python3

import pandas as pd

print("=== Step 5D: Validate Cleaned Messy Expression Matrix ===")

file_path = "data/processed/mock_expression_matrix_messy_clean.csv"
print(f"[INFO] Loading cleaned file: {file_path}")
df = pd.read_csv(file_path)

print("\n[INFO] Shape:", df.shape)

print("\n[INFO] Column names (repr):")
for c in df.columns:
    print(f" - {repr(c)}")

# Check sample column naming format
sample_cols = [c for c in df.columns if c != "gene_id"]

print("\n[INFO] Sample label format check (expected: sample_###_Group)")
all_format_ok = True
for c in sample_cols:
    parts = c.split("_")
    if len(parts) != 3:
        all_format_ok = False
        print(f"[FAIL] Bad format: {repr(c)} -> {parts}")
        continue

    prefix, number, group = parts
    if prefix != "sample":
        all_format_ok = False
        print(f"[FAIL] Bad prefix in {repr(c)}: {prefix}")
        continue

    if not number.isdigit() or len(number) != 3:
        all_format_ok = False
        print(f"[FAIL] Bad sample number in {repr(c)}: {number}")
        continue

    if group not in {"Tumor", "Normal"}:
        all_format_ok = False
        print(f"[FAIL] Bad group label in {repr(c)}: {group}")
        continue

    print(f"[OK]   {repr(c)}")

if all_format_ok:
    print("[OK] All sample labels are standardized.")

# Missing value report
print("\n[INFO] Missing values per column:")
print(df.isna().sum())

# Locate rows with any missing values
rows_with_na = df[df.isna().any(axis=1)]

print("\n[INFO] Rows containing missing values:")
if rows_with_na.empty:
    print("[OK] No rows with missing values.")
else:
    print(rows_with_na)

print("\n[OK] Validation complete.")
