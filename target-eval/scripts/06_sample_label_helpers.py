#!/usr/bin/env python3

import re
from typing import Optional, Dict


def clean_sample_label(label: str) -> str:
    """
    Clean and standardize a sample label.

    Rules (for this training project):
    - strip leading/trailing whitespace
    - fix 'sample003_Normal' -> 'sample_003_Normal'
    - standardize group capitalization (tumor -> Tumor)
    """
    cleaned = label.strip()

    # Fix missing underscore after 'sample'
    cleaned = re.sub(r"^sample(\d{3})_", r"sample_\1_", cleaned)

    # Standardize group capitalization if format is close to expected
    parts = cleaned.split("_")
    if len(parts) == 3:
        prefix, number, group = parts
        cleaned = f"{prefix}_{number}_{group.capitalize()}"

    return cleaned


def parse_sample_label(label: str) -> Optional[Dict[str, str]]:
    """
    Parse a standardized sample label.

    Expected format: sample_###_Group
    Example: sample_001_Tumor

    Returns:
        dict with keys sample_column, sample_id, group
        OR None if parsing fails
    """
    cleaned = clean_sample_label(label)
    parts = cleaned.split("_")

    if len(parts) != 3:
        return None

    prefix, number, group = parts

    if prefix != "sample":
        return None

    if not number.isdigit() or len(number) != 3:
        return None

    if group not in {"Tumor", "Normal"}:
        return None

    return {
        "sample_column": cleaned,
        "sample_id": f"{prefix}_{number}",
        "group": group,
    }


if __name__ == "__main__":
    print("=== Step 6A: Sample Label Helper Demo ===")

    test_labels = [
        "sample_001_Tumor",
        "sample_002_tumor",
        "sample003_Normal",
        "sample_004_Normal ",
        "bad_label",
    ]

    for label in test_labels:
        cleaned = clean_sample_label(label)
        parsed = parse_sample_label(label)

        print(f"\nOriginal: {repr(label)}")
        print(f"Cleaned:  {repr(cleaned)}")
        print(f"Parsed:   {parsed}")
