#!/usr/bin/env bash
set -euo pipefail

echo "=== Target Eval Project Setup Check ==="
echo "Current directory: $(pwd)"
echo

echo "Checking required folders..."
for d in data/raw data/processed metadata scripts notebooks results figures docs env; do
  if [ -d "$d" ]; then
    echo "  [OK] $d"
  else
    echo "  [MISSING] $d"
  fi
done
echo

echo "Checking key files..."
for f in README.md .gitignore; do
  if [ -f "$f" ]; then
    echo "  [OK] $f"
  else
    echo "  [MISSING] $f"
  fi
done
echo

echo "Git repo root (if inside a repo):"
git rev-parse --show-toplevel 2>/dev/null || echo "  Not inside a git repo"
echo

echo "Git status (from current folder):"
git status --short || true
echo
echo "Setup check complete."
