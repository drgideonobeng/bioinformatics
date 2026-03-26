#!/usr/bin/env python3

import pandas as pd

print("=== Step 3B: Inspect Practice Data ===")

# File path
file_path = "data/raw/countries_practice.csv"
print(f"[INFO] Loading file: {file_path}")

# Load CSV into a pandas DataFrame
df = pd.read_csv(file_path)

print("\n[INFO] Data loaded successfully.")
print(f"[INFO] Shape (rows, columns): {df.shape}")

print("\n[INFO] First 5 rows:")
print(df.head())

print("\n[INFO] Column names:")
for col in df.columns:
    print(f" - {col}")

print("\n[INFO] Data types:")
print(df.dtypes)

print("\n[INFO] Missing values per column:")
print(df.isna().sum())

print("\n[OK] Inspection complete.")
