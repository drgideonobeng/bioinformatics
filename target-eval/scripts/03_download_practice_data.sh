#!/usr/bin/env bash
set -euo pipefail

echo "=== Step 3: Download Practice Data ==="

# Create raw data folder if it does not exist
mkdir -p data/raw

# Download a small public CSV for table-inspection practice
URL="https://raw.githubusercontent.com/cs109/2014_data/master/countries.csv"
OUT="data/raw/countries_practice.csv"

echo "[INFO] Downloading practice dataset..."
curl -L "$URL" -o "$OUT"

echo "[OK] Saved to: $OUT"

echo
echo "[INFO] File size:"
wc -c "$OUT"

echo
echo "[INFO] File preview (first 5 lines):"
head -n 5 "$OUT"
