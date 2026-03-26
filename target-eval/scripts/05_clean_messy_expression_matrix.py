#!/usr/bin/env python3

import pandas as pd
import re

print("=== Step 5C: Clean Messy Expression Matrix ===")

input_file = "data/raw/mock_expression_matrix_messy.csv"
output_file = "data/processed/mock_expression_matrix_messy_clean.csv"

print(f"[INFO] Loading raw file: {input_file}")
df = pd.read_csv(input_file)

# -----------------------------
# 1) Clean gene IDs (strip whitespace)
# -----------------------------
df["gene_id"] = df["gene_id"].str.strip()

# -----------------------------
# 2) Clean column names (strip whitespace)
# -----------------------------
original_cols = df.columns.tolist()
cleaned_cols = [c.strip() for c in original_cols]
df.columns = cleaned_cols

# -----------------------------
# 3) Standardize sample column names
#    Goal format: sample_###_Group
# -----------------------------
new_columns = []

for col in df.columns:
    if col == "gene_id":
        new_columns.append(col)
        continue

    fixed = col

    # Fix malformed "sample003_Normal" -> "sample_003_Normal"
    fixed = re.sub(r"^sample(\d{3})_", r"sample_\1_", fixed)

    # Standardize group labels: tumor -> Tumor, normal -> Normal
    parts = fixed.split("_")
    if len(parts) == 3:
        prefix, num, group = parts
        group = group.capitalize()   # tumor -> Tumor
        fixed = f"{prefix}_{num}_{group}"

    new_columns.append(fixed)

df.columns = new_columns

# -----------------------------
# 4) Report missing values (do not fill yet)
# -----------------------------
missing_counts = df.isna().sum()

print("\n[INFO] Missing values after structural cleaning:")
print(missing_counts)

# Save cleaned file
df.to_csv(output_file, index=False)

print(f"\n[OK] Cleaned matrix saved to: {output_file}")

print("\n[INFO] Cleaned column names (repr):")
for c in df.columns:
    print(f" - {repr(c)}")

print("\n[INFO] Cleaned gene IDs (repr, first 5):")
for g in df["gene_id"].head():
    print(f" - {repr(g)}")

print("\n[INFO] Preview of cleaned matrix:")
print(df)

print("\n[OK] Cleaning complete.")
