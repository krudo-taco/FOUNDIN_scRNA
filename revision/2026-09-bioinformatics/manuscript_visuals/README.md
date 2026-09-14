# Manuscript Figure 2 / Table 4 Display Runner

This directory is a small public-facing display package for the revised manuscript
Figure 2 and Table 4. It does not rerun the full bioinformatics analyses.

## Contents

- `source_data/primary86.csv` - final 86 primary SNCA-autophagy estimates.
- `source_data/Figure2_C_coverage_source.csv` - final donor coverage panel data.
- `source_data/Figure2_D_profile_source.csv` - final fixed-direction profile panel data.
- `provenance/build_figures_tables_final06.py` - exact final build script copied from `06-manuscript-revision`.
- `build_fig2_table4.py` - portable wrapper that imports the provenance script and runs only `load_primary`, `build_table4`, and `build_figure2`.
- `validate_fig2_table4.py` - checks source row counts, required outputs, and Table 4 equality to the bundled reference table.
- `captions.json`, `Table2_manuscript.md`, `Table3_manuscript.md`, `Table4_manuscript.md` - final manuscript-facing display text.

Complete differential-expression, CAMERA/ROAST and supplementary source tables
accompany the manuscript supplement; their analysis scripts are retained in the
code package. This runner is intentionally
limited to Figure 2 and Table 4 display regeneration from small CSV files.

## Tested Environment

The current local run used:

- matplotlib 3.1.2
- Pillow 7.0.0
- tabulate 0.9.0

`requirements.txt` records these display dependencies. Existing project
environments may already provide pandas and numpy.

## Run

From this directory:

```bash
python build_fig2_table4.py
python validate_fig2_table4.py
```

Outputs are written under `work/revised_figures/`.

The validator is expected to report:

- 86 primary estimate rows.
- 43 primary genes.
- 43 Table 4 rows.
- regenerated Table 4 equals the final `source_data/Table4_reference.csv`.
- Figure 2 SVG/PDF/PNG/TIFF and parent-compatible `Figure2_revised.png` exist.
