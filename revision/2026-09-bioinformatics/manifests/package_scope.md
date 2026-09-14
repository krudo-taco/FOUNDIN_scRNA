# Package Scope And Verification Levels

## Included

- Portable statistical kernel code copied from the verified analyses.
- Original analysis scripts from `04-reviewer-reanalysis` and `05-claim-support`
  under `src/original/` for provenance.
- A public GSE243639 rerun path that downloads public raw files, rebuilds
  processed inputs, reruns the pilot disease-effect analysis, and regenerates a
  summary table plus SVG figure.
- Synthetic non-identifying tests for HC3/wild-bootstrap behavior, multiple
  testing, and partial-correlation constant-vector handling.
- Aggregate reference outputs that do not expose restricted donor-level raw
  inputs.
- Functional staging wrappers for public GSE243639, local-existing Kamath, and
  authorized cached FOUNDIN reruns.
- Exact reference comparisons for 86 FOUNDIN primary edges, 264 GSE243639 edges,
  16 GSE243639 panel tests, 1008 Kamath edges and 48 Kamath panel/profile tests.

## Excluded

- Restricted FOUNDIN/PPMI donor-level clinical workbooks and raw inputs.
- Credentials. Legacy provenance scripts may contain historical local path
  strings; portable runners patch staged copies and do not require those paths
  for cached-input or public reruns.
- Private manuscript drafts, reviewer letters, or coauthor response documents.
- Large raw public files; these are downloaded or placed by the user.

## Verification Levels

1. **Synthetic/kernel**: runs without restricted data and validates numerical
   primitives. This is CI-friendly but is not full biological reproduction.
2. **Public GSE243639 rerun**: downloads public raw data and reruns the complete
   low-coverage pilot analysis plus table/figure generation. This is the main
   fresh-run external validation target for the public package.
3. **Kamath local-existing rerun**: verifies staged public raw/processed input
   hashes, reruns metadata preprocessing plus scripts 13 and 14, and compares all
   planned outcome rows. It does not re-stream the 5GB matrix by default.
4. **Authorized FOUNDIN cached rerun**: requires restricted local cached inputs,
   reruns the primary donor-level pipeline, and compares all 86 primary edge
   rows. Public users can inspect code and aggregate outputs but cannot fully
   rerun restricted donor-level analyses without authorized data access.
