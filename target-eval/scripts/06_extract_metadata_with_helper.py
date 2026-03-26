#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import sys

# Allow importing helper script from the same folder
sys.path.append(str(Path(__file__).parent))

from sample_label_helpers import parse_sample_label  # noqa: E402

print("=== Step 6B: Extract Metadata Using Helper Functions ===")

expr_file = "data/raw/mock_expression_matrix_messy.csv"
out_file = "metadata/mock_sample_metadata_from_helper.csv"

print(f"[INFO] Loading expression matrix: {expr_file}")
df = pd.read_csv(expr_file)

# Identify sample columns
sample_cols = [c for c in df.columns if c != "gene_id"]

print("\n[INFO] Parsing sample columns with helper...")
metadata_rows = []
invalid_labels = []

for col in sample_cols:
    parsed = parse_sample_label(col)
    if parsed is None:
        invalid_labels.append(col)
    else:
        metadata_rows.append(parsed)

metadata_df = pd.DataFrame(metadata_rows)

print(f"[INFO] Parsed valid sample columns: {len(metadata_rows)}")
print(f"[INFO] Invalid sample columns: {len(invalid_labels)}")

if invalid_labels:
    print("\n[WARN] Invalid labels detected:")
    for x in invalid_labels:
        print(f" - {repr(x)}")

# Save valid metadata rows
metadata_df.to_csv(out_file, index=False)
print(f"\n[OK] Saved metadata table to: {out_file}")

print("\n[INFO] Metadata preview:")
print(metadata_df)

print("\n[OK] Helper-based metadata extraction complete.")
