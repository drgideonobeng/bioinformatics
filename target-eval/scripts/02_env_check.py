#!/usr/bin/env python3

import sys

print("=== Step 2: Environment Check ===")
print(f"Python executable: {sys.executable}")
print(f"Python version: {sys.version.split()[0]}")
print()

# Import and print package versions
packages = {}

try:
    import pandas as pd
    packages["pandas"] = pd.__version__
except Exception as e:
    packages["pandas"] = f"IMPORT FAILED ({e})"

try:
    import numpy as np
    packages["numpy"] = np.__version__
except Exception as e:
    packages["numpy"] = f"IMPORT FAILED ({e})"

try:
    import matplotlib
    packages["matplotlib"] = matplotlib.__version__
except Exception as e:
    packages["matplotlib"] = f"IMPORT FAILED ({e})"

try:
    import openpyxl
    packages["openpyxl"] = openpyxl.__version__
except Exception as e:
    packages["openpyxl"] = f"IMPORT FAILED ({e})"

try:
    import IPython
    packages["IPython"] = IPython.__version__
except Exception as e:
    packages["IPython"] = f"IMPORT FAILED ({e})"

print("Package check:")
for name, version in packages.items():
    status = "[OK]" if not str(version).startswith("IMPORT FAILED") else "[FAIL]"
    print(f"  {status} {name}: {version}")

print()

# Simple pandas test table (important for later metadata work)
try:
    import pandas as pd
    df = pd.DataFrame({
        "sample_id": ["S1", "S2", "S3"],
        "group": ["Tumor", "Normal", "Tumor"],
        "expression": [12.4, 3.8, 9.1]
    })
    print("Test DataFrame:")
    print(df)
    print("\n[OK] pandas DataFrame test passed.")
except Exception as e:
    print(f"[FAIL] pandas DataFrame test failed: {e}")
