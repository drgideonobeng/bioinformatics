#!/usr/bin/env bash
set -euo pipefail

echo "=== Step 2: Python Environment Setup ==="

if [ ! -d ".venv" ]; then
  echo "[INFO] Creating virtual environment (.venv)..."
  python3 -m venv .venv
else
  echo "[OK] .venv already exists"
fi

echo "[INFO] Activating virtual environment..."
# shellcheck disable=SC1091
source .venv/bin/activate

echo "[INFO] Upgrading pip..."
pip install --upgrade pip

echo "[INFO] Installing packages from env/requirements.txt..."
pip install -r env/requirements.txt

echo
echo "[OK] Python environment is ready."
echo "To activate it later, run: source .venv/bin/activate"
