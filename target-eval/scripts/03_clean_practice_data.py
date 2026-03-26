#!/usr/bin/env python3

import pandas as pd

print("=== Step 3C: Clean Practice Data ===")

# Input and output paths
input_file = "data/raw/countries_practice.csv"
output_file = "data/processed/countries_practice_clean.csv"

print(f"[INFO] Loading raw file: {input_file}")
df = pd.read_csv(input_file)

# Basic cleaning (good habits, even if data already looks clean)
# 1) Remove leading/trailing whitespace from text columns
df["Country"] = df["Country"].str.strip()
df["Region"] = df["Region"].str.strip()

# 2) Standardize region names to Title Case (e.g., AFRICA -> Africa)
df["Region"] = df["Region"].str.title()

# 3) Drop duplicate rows if any
before_rows = len(df)
df = df.drop_duplicates()
after_rows = len(df)

print(f"[INFO] Rows before deduplication: {before_rows}")
print(f"[INFO] Rows after deduplication:  {after_rows}")

# Save cleaned file
df.to_csv(output_file, index=False)
print(f"[OK] Cleaned file saved to: {output_file}")

print("\n[INFO] Preview of cleaned data:")
print(df.head())

print("\n[OK] Cleaning complete.")
