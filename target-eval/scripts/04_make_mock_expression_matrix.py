#!/usr/bin/env python3

import pandas as pd

print("=== Step 4A: Create Mock Expression Matrix ===")

# Create a small gene-expression-style table
# Rows = genes
# Columns = samples
df = pd.DataFrame({
    "gene_id": ["TP53", "EGFR", "MYC", "GAPDH", "ACTB"],
    "sample_001_Tumor":  [12.5, 45.2, 30.1, 100.0, 95.4],
    "sample_002_Tumor":  [10.8, 50.3, 28.9, 102.1, 97.0],
    "sample_003_Normal": [25.4,  5.2,  8.1,  98.7, 96.5],
    "sample_004_Normal": [24.8,  4.9,  7.5, 101.2, 98.1],
})

out_file = "data/raw/mock_expression_matrix.csv"
df.to_csv(out_file, index=False)

print(f"[OK] Mock expression matrix saved to: {out_file}")
print("\n[INFO] Preview:")
print(df)
