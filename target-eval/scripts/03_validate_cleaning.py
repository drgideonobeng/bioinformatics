#!/usr/bin/env python3

import pandas as pd

print("=== Step 3D: Validate Cleaning ===")

raw_file = "data/raw/countries_practice.csv"
clean_file = "data/processed/countries_practice_clean.csv"

print(f"[INFO] Loading raw file:   {raw_file}")
raw_df = pd.read_csv(raw_file)

print(f"[INFO] Loading clean file: {clean_file}")
clean_df = pd.read_csv(clean_file)

print("\n[INFO] Shapes")
print(f" - Raw:   {raw_df.shape}")
print(f" - Clean: {clean_df.shape}")

print("\n[INFO] First 5 rows (raw):")
print(raw_df.head())

print("\n[INFO] First 5 rows (clean):")
print(clean_df.head())

# Compare unique Region values
raw_regions = sorted(raw_df["Region"].dropna().unique().tolist())
clean_regions = sorted(clean_df["Region"].dropna().unique().tolist())

print("\n[INFO] Unique Region values (raw):")
print(raw_regions)

print("\n[INFO] Unique Region values (clean):")
print(clean_regions)

# Basic validation checks
if raw_df.shape[0] != clean_df.shape[0]:
    print("\n[WARN] Row count changed during cleaning.")
else:
    print("\n[OK] Row count unchanged.")

if raw_regions != clean_regions:
    print("[OK] Region labels changed as expected (standardization applied).")
else:
    print("[WARN] Region labels did not change.")

print("\n[OK] Validation complete.")
