#!/usr/bin/env python3
"""Validate the portable Figure 2 / Table 4 display runner outputs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


PACKAGE = Path(__file__).resolve().parent
WORK = PACKAGE / "work" / "revised_figures"
REFERENCE_TABLE4 = PACKAGE / "source_data" / "Table4_reference.csv"
PACKAGE_TABLE4 = WORK / "tables" / "Table4_revised_primary_SNCA_autophagy_estimates.csv"
PRIMARY = PACKAGE / "source_data" / "primary86.csv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    primary = pd.read_csv(PRIMARY)
    table4 = pd.read_csv(PACKAGE_TABLE4)
    reference = pd.read_csv(REFERENCE_TABLE4)
    if primary.shape[0] != 86:
        raise AssertionError(f"Expected 86 primary rows, found {primary.shape[0]}")
    if primary["gene"].nunique() != 43:
        raise AssertionError(f"Expected 43 primary genes, found {primary['gene'].nunique()}")
    pd.testing.assert_frame_equal(table4, reference, check_dtype=False, check_like=False)
    expected_files = [
        WORK / "figures" / "Figure_2_revised.svg",
        WORK / "figures" / "Figure_2_revised.pdf",
        WORK / "figures" / "Figure_2_revised.png",
        WORK / "figures" / "Figure_2_revised.tiff",
        WORK / "figures" / "Figure2_revised.png",
        PACKAGE_TABLE4,
    ]
    missing = [str(p.relative_to(PACKAGE)) for p in expected_files if not p.exists()]
    if missing:
        raise AssertionError(f"Missing output files: {missing}")
    result = {
        "status": "passed",
        "scope": "Figure 2 and Table 4 display runner only",
        "primary_rows": int(primary.shape[0]),
        "primary_gene_count": int(primary["gene"].nunique()),
        "table4_rows": int(table4.shape[0]),
        "table4_equals_06_reference": True,
        "package_table4_sha256": sha256(PACKAGE_TABLE4),
        "reference_table4_sha256": sha256(REFERENCE_TABLE4),
        "outputs": [str(p.relative_to(PACKAGE)) for p in expected_files],
    }
    out = WORK / "fig2_table4_validation.json"
    out.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
