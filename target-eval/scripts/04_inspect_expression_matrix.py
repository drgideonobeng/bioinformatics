#!/usr/bin/env python3

import pandas as pd

print("=== Step 4B: Inspect Mock Expression Matrix ===")

file_path = "data/raw/mock_expression_matrix.csv"
print(f"[INFO] Loading file: {file_path}")

df = pd.read_csv(file_path)

print("\n[INFO] Shape (rows, columns):", df.shape)

print("\n[INFO] First 5 rows:")
print(df.head())

# Identify the gene ID column (first column in this simple example)
gene_col = df.columns[0]
sample_cols = df.columns[1:]

print(f"\n[INFO] Gene ID column: {gene_col}")

print("\n[INFO] Sample columns:")
for col in sample_cols:
    print(f" - {col}")

print("\n[INFO] Data types:")
print(df.dtypes)

# Check whether sample columns are numeric
print("\n[INFO] Numeric check for sample columns:")
for col in sample_cols:
    is_numeric = pd.api.types.is_numeric_dtype(df[col])
    print(f" - {col}: {'numeric' if is_numeric else 'NOT numeric'}")

# Basic summary statistics for expression values
print("\n[INFO] Summary statistics (sample columns):")
print(df[sample_cols].describe())

print("\n[OK] Expression matrix inspection complete.")
