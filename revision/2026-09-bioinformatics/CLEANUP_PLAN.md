# Cleanup Plan

## Scope

This package is a portable, public-deposition-ready wrapper for the verified
`04-reviewer-reanalysis` and `05-claim-support` computational analyses. Edits
are limited to files under `07-reproducibility-package/`.

The package will not include restricted donor-level FOUNDIN/PPMI clinical
workbooks, private manuscript or reviewer text, credentials, or large raw data
that require separate repository access. It may include source code, safe
aggregate result tables, provenance manifests, synthetic non-identifying test
data, and instructions for authorized users to place raw inputs locally.

## Minimal Necessary Changes

1. Copy analysis code into `src/` and preserve original file names where
   possible.
2. Add portable configuration through explicit input and output roots instead of
   hard-coded machine paths.
3. Add wrapper scripts for three verification levels:
   - synthetic/kernel tests that run without restricted data;
   - public-data checks for Kamath/GSE inputs when users download public files;
   - authorized FOUNDIN raw-data reruns for users with local access.
4. Add data-access manifests that classify each input as public, restricted, or
   derived aggregate.
5. Add dependency/version records for Python and R.
6. Add GitHub Actions workflow that validates the public, synthetic layer on a
   fresh runner.
7. Add regression checks for statistical kernels and published aggregate outputs
   without redistributing restricted raw donor-level data.
8. Add functional staging wrappers that copy original scripts into `work/`,
   patch machine-specific paths in the staged copy only, verify declared input
   SHA256 hashes before computation, run the relevant pipeline, and compare
   complete result families against package references.

## Baseline Tests Before Packaging

Run these existing validators before making copied-code adaptations:

1. `python 04-reviewer-reanalysis/check_statistics.py`
2. `python 04-reviewer-reanalysis/check_actual_design.py`
3. `python 05-claim-support/05_verify_results.py`
4. `python 05-claim-support/14_verify_kamath_results.py`

Passing these tests confirms the current source directories are internally
consistent before the packaging layer copies or wraps them.

## Regression Tests After Packaging

1. `python 07-reproducibility-package/tests/test_synthetic_kernels.py`
2. `python 07-reproducibility-package/tools/verify_package.py`
3. Public GSE rerun:
   `python 07-reproducibility-package/run_public_gse243639.py --config 07-reproducibility-package/config/public_gse243639.json`
4. Authorized cached FOUNDIN rerun:
   `python 07-reproducibility-package/run_authorized_foundin.py --config 07-reproducibility-package/config/authorized_foundin_cached.json`
5. Local-existing Kamath rerun:
   `python 07-reproducibility-package/run_manual_kamath_heavy.py --config 07-reproducibility-package/config/kamath_local_existing.json`

Synthetic/kernel tests prove code portability and key numerical behavior only.
They do not constitute full external reproduction of the paper analyses.

## Preservation Rules

- Preserve random seeds, family sizes, multiple-testing corrections, and
  statistical choices from the verified analyses.
- Preserve the Kamath frozen plan hash behavior by recording and checking the
  source plan hash, rather than editing the original `05-claim-support`
  `analysis_plan.md`.
- Do not silently relabel gene activity as transcriptomics or upgrade
  exploratory associations to confirmatory claims.
- Flag licensing and data-access uncertainty explicitly.
