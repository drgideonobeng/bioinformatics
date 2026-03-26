#!/usr/bin/env python3

import pandas as pd
import numpy as np

print("=== Step 5A: Create Messy Mock Expression Matrix ===")

# Create a mock expression matrix with intentional issues:
# - whitespace in gene IDs
# - inconsistent sample labels
# - missing values
df = pd.DataFrame({
    "gene_id": [" TP53", "EGFR ", "MYC", "GAPDH", "ACTB"],
    "sample_001_Tumor":   [12.5, 45.2, 30.1, 100.0, 95.4],
    "sample_002_tumor":   [10.8, np.nan, 28.9, 102.1, 97.0],   # lowercase 'tumor' + missing value
    "sample003_Normal":   [25.4, 5.2, 8.1, 98.7, 96.5],        # missing underscore after 'sample'
    "sample_004_Normal ": [24.8, 4.9, 7.5, 101.2, 98.1],       # trailing space in column name
})

out_file = "data/raw/mock_expression_matrix_messy.csv"
df.to_csv(out_file, index=False)

print(f"[OK] Messy mock expression matrix saved to: {out_file}")
print("\n[INFO] Preview:")
print(df)

print("\n[INFO] Column names (repr, so spaces are visible):")
for c in df.columns:
    print(repr(c))
