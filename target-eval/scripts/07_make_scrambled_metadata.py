#!/usr/bin/env python3

import pandas as pd

print("=== Step 7B: Create Scrambled Metadata (for practice) ===")

input_file = "metadata/mock_sample_metadata_from_helper.csv"
output_file = "metadata/mock_sample_metadata_scrambled.csv"

print(f"[INFO] Loading metadata: {input_file}")
meta_df = pd.read_csv(input_file)

print("\n[INFO] Original metadata order:")
print(meta_df["sample_column"].tolist())

# Scramble row order on purpose (fixed order for reproducibility)
scrambled_order = [2, 0, 3, 1]
scrambled_df = meta_df.iloc[scrambled_order].reset_index(drop=True)

print("\n[INFO] Scrambled metadata order:")
print(scrambled_df["sample_column"].tolist())

# Save scrambled file
scrambled_df.to_csv(output_file, index=False)
print(f"\n[OK] Scrambled metadata saved to: {output_file}")

print("\n[INFO] Scrambled metadata preview:")
print(scrambled_df)

print("\n[OK] Done.")
