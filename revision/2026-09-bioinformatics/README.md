# 2026-09 Bioinformatics Reproducibility Package

This directory packages the verified 04-reviewer-reanalysis and
05-claim-support analyses for public repository deposition. It is designed to be
copied into a public repository as:

```text
revision/2026-09-bioinformatics/
```

The package separates code portability from data access. Public raw data can be
downloaded and rerun. Restricted FOUNDIN/PPMI donor-level files are not included
and must only be used on machines with authorized local access.

## What Is Included

- `src/repro_package/`: portable statistical kernels used by tests and wrappers.
- `src/original/`: source analysis scripts copied from the verified local
  analyses for provenance. These files intentionally preserve legacy source
  text, including machine-specific paths. The runnable wrappers copy them into
  `work/` and patch the staged copy only.
- `run_public_gse243639.py`: full rerun of the public GSE243639 pilot from raw
  GEO files.
- `run_authorized_foundin.py`: authorized cached-input rerun of the primary
  FOUNDIN donor-level pipeline, followed by all-86-edge comparison.
- `run_manual_kamath_heavy.py`: local-existing public Kamath rerun from staged
  raw/processed inputs, followed by all-1008-edge and all-48-test comparison.
- `tools/make_summary_tables.py`: regenerates a compact summary table and SVG
  figure from the public pilot outputs.
- `tools/compare_results.py`: exact row-family numerical comparisons.
- `tests/test_synthetic_kernels.py`: non-identifying synthetic tests for
  statistical kernels.
- `manifests/data_access_manifest.csv`: access route and redistribution status
  for each input class.
- `ci/github-actions/revision-2026-09-bioinformatics.yml`: GitHub Actions
  workflow template for the public repository. At deposition, copy this file to
  the repository-root `.github/workflows/` directory. A duplicate is also kept
  under this package's `.github/workflows/` for direct inspection.

## Verification Levels

Synthetic/kernel checks validate code portability and selected numerical
properties. They do not reproduce the biological analyses.

The public GSE243639 rerun downloads public raw files, regenerates processed
inputs, reruns the pilot disease-effect analysis, and rebuilds a table and SVG
figure. This is the main external-machine validation target.

The Kamath/GSE178265 path is public but heavier. The wrapper validates local
public raw/processed input hashes, runs metadata preprocessing plus the frozen
outcome scripts `13` and `14`, and compares all 1008 edge rows plus 48
panel/profile tests. The default local-existing mode does not re-stream the
5GB matrix; that heavier extraction code is preserved under `src/original/`.

The FOUNDIN rerun requires restricted local inputs. The included
verified-cached mode reruns the current primary donor-level pipeline from
checked cached 04 inputs and compares all 86 primary edges; it does not reload
the original large Seurat object.

## Quick Local Checks

Run commands from this package directory unless an absolute path is supplied.
Config paths are relative to this directory; edit `config/*.json` for another
machine's local input layout.

```bash
python -m pip install -r requirements.txt
python tests/test_synthetic_kernels.py
python tests/test_staging_helpers.py
python tools/verify_package.py
```

## Full Public GSE243639 Rerun

To download from GEO and rerun the public pilot:

```bash
python run_public_gse243639.py --config config/public_gse243639.json
python tools/make_summary_tables.py --work-root work/public_gse243639
python tools/compare_results.py \
  --profile gse243639_edges264 \
  --observed work/public_gse243639/results/gse243639_pilot_all_edges.csv \
  --reference results_reference/public_gse243639/gse243639_pilot_all_edges.csv \
  --out work/public_gse243639/gse243639_edges264_comparison.json
python tools/compare_results.py \
  --profile gse243639_panel16 \
  --observed work/public_gse243639/results/gse243639_pilot_panel_tests.csv \
  --reference results_reference/public_gse243639/gse243639_pilot_panel_tests.csv \
  --out work/public_gse243639/gse243639_panel16_comparison.json
```

If the public files are already present elsewhere:

```bash
python run_public_gse243639.py \
  --config config/public_gse243639.json \
  --external-root /path/to/GSE243639/files \
  --skip-download
python tools/make_summary_tables.py --work-root work/public_gse243639
```

Expected public input files:

```text
GSE243639_Clinical_data.csv.gz
GSE243639_UMAP_coordinates.xlsx
GSE243639_Filtered_count_table.csv.gz
```

The rerun writes outputs under `work/public_gse243639/`.

## Authorized FOUNDIN Rerun

The public package intentionally excludes restricted donor-level inputs. On an
authorized local machine with the verified cached 04 inputs available:

```bash
python run_authorized_foundin.py --config config/authorized_foundin_cached.json
```

This stages a private working copy under `work/authorized_foundin/`, patches the
gene-list path inside the staged script only, verifies cached input hashes,
reruns `03_within_donor_coupling.py`, extracts the 86 primary edges, and compares
effect sizes, HC3 statistics, P values and q values against the reference table.

## Kamath Local-Existing Rerun

For the public Kamath analysis using existing local public inputs and processed
count/annotation files:

```bash
python run_manual_kamath_heavy.py --config config/kamath_local_existing.json
```

This stages all declared public inputs, verifies SHA256 hashes, runs
`04_prepare_kamath_metadata.py`, `13_test_kamath.py` and
`14_verify_kamath_results.py`, then compares all 1008 planned edge rows and all
48 panel/profile tests. It is a validated processed-input rerun, not a full
matrix re-extraction.

For a different machine, edit `config/kamath_local_existing.json` so
`external_root` and `processed_input_root` point to the local public files before
running the command.

## Code And Data Availability Draft

Code and non-restricted aggregate outputs are prepared for the authorized public
code deposit after the local checks reported in `manifests/local_validation.json`.
Public GSE243639 files are retrieved from GEO by script; Kamath public
local-existing inputs are hash-checked before rerun.
Restricted FOUNDIN/PPMI donor-level files are not redistributed; full rerun of
those analyses requires authorized local access. The Kamath qs annotation was
used from a public access route during local analysis, but redistribution and
license status should be confirmed before bundling that file in any public
repository.
