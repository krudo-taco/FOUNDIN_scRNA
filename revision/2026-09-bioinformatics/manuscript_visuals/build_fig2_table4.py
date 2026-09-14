#!/usr/bin/env python3
"""Portable runner for revised Figure 2 and Table 4 display assets.

This script intentionally consumes only the small display CSVs bundled in
``source_data``. It does not rerun the upstream bioinformatics analyses and does
not require the full differential-expression, CAMERA, or donor-level input files.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
from pathlib import Path


PACKAGE = Path(__file__).resolve().parent
ROOT = PACKAGE
OUT = PACKAGE / "work" / "revised_figures"
FIG_DIR = OUT / "figures"
TAB_DIR = OUT / "tables"
FIG_SOURCE = FIG_DIR / "source_data"
TAB_SOURCE = TAB_DIR / "source_data"
SOURCE = PACKAGE / "source_data"
PROVENANCE_SCRIPT = PACKAGE / "provenance" / "build_figures_tables_final06.py"


def load_original_module():
    spec = importlib.util.spec_from_file_location("final06_figures", PROVENANCE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load provenance script: {PROVENANCE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare_work_inputs(module) -> None:
    for path in [FIG_DIR, TAB_DIR, FIG_SOURCE, TAB_SOURCE]:
        path.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(SOURCE / "primary86.csv", FIG_SOURCE / "Figure2_AB_primary_86_estimates.csv")
    shutil.copyfile(
        SOURCE / "Figure2_C_coverage_source.csv",
        FIG_SOURCE / "Figure2_C_donor_cell_coverage_source.csv",
    )
    shutil.copyfile(
        SOURCE / "Figure2_D_profile_source.csv",
        FIG_SOURCE / "Figure2_D_fixed_direction_profile_source.csv",
    )
    module.ROOT = ROOT
    module.OUT = OUT
    module.FIG_DIR = FIG_DIR
    module.TAB_DIR = TAB_DIR
    module.FIG_SOURCE = FIG_SOURCE
    module.TAB_SOURCE = TAB_SOURCE
    module.PRIMARY = SOURCE / "primary86.csv"


def main() -> None:
    module = load_original_module()
    prepare_work_inputs(module)
    primary = module.load_primary()
    table4 = module.build_table4(primary)
    profile = module.pd.read_csv(FIG_SOURCE / "Figure2_D_fixed_direction_profile_source.csv")
    coverage = module.pd.read_csv(FIG_SOURCE / "Figure2_C_donor_cell_coverage_source.csv")
    fig2_outputs = module.build_figure2(primary, profile, coverage)
    shutil.copyfile(FIG_DIR / "Figure_2_revised.png", FIG_DIR / "Figure2_revised.png")
    summary = {
        "status": "built",
        "scope": "Figure 2 and Table 4 display assets only",
        "primary_rows": int(primary.shape[0]),
        "primary_genes": int(primary["gene"].nunique()),
        "table4_rows": int(table4.shape[0]),
        "fig2_outputs": fig2_outputs,
        "table4_csv": str((TAB_DIR / "Table4_revised_primary_SNCA_autophagy_estimates.csv").relative_to(PACKAGE)),
    }
    (OUT / "fig2_table4_build_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
