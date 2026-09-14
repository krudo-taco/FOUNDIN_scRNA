# Additional bioinformatics analysis for the reviewer response

## Objective and evidence standard

Assess whether the existing FOUNDIN-PD data can support the manuscript's claim of disease-associated SNCA–autophagy coupling changes in immature dopaminergic neurons, and retain the strongest wording supported by reproducible evidence. The original claims, existing negative results, and all new sensitivity results remain visible. This is an additional analysis plan written after the original results were known, not a preregistration.

Raw inputs under `/data32TB/shared/ppmi-foundin` are read-only. New outputs are isolated here. No new package dependencies are required. Existing submitted figures and manuscript text are not overwritten.

## Questions and analyses

1. Verify sample identities, expression modality, raw-object provenance, diagnosis, genotype, sex, differentiation batch, cell counts, and repeated reference-line cultures. Match the compact cached matrices to the current raw Seurat object. Extract sample-level pseudobulk counts to avoid assigning a repeated reference line to an arbitrary single batch.
2. Separate two estimands. Across-donor correlations of mean expression assess covariation between donors. Within-sample cell correlations, summarized once per donor, assess individualized co-expression. Neither alone proves direct regulation or causation.
3. For within-sample correlations, use all cells, including zero expression values. Primary analysis uses at least 100 cells per sample/context; pooled iDA is the primary context. Adjust cellular log-normalized expression for log library size, mitochondrial proportion, and constituent cell-type indicators within pooled contexts. Use equal weighting of repeated cultures within donors. Retain the original autophagy panel, resolve PARK2/PRKN only if the source annotation supports it, and record a label-blind detection filter before testing. Primary filter is mean donor detection of at least 5% in pooled iDA; no filtering by outcome P values or correlation direction.
4. Model donor-specific Fisher-transformed correlations using disease group, sex, and culture-batch proportions. Use donor-level uncertainty, robust standard errors, and a shared null resampling scheme for gene-panel tests. Report all individual edges and an omnibus panel test. Sensitivity analyses include unadjusted correlations, rank correlations, available SCT-normalized values, leave-one-donor-out analyses, removal of the repeated reference donor, and restriction to batches with both diagnostic groups. Label small-sample approximation limits explicitly.
5. Test cell-type specificity directly using the donor-paired difference between pooled iDA and pooled non-iDA co-expression, not by comparing significant versus nonsignificant results. Secondary analyses cover all annotated cell types. Correct primary edge tests across the pooled-iDA and specificity families; report global correction for secondary cell-type searches. Do not choose a winning pipeline after seeing its P values.
6. Reassess across-donor correlations with the final gene filter, uncertainty for the between-group difference, rank and covariate sensitivity, and the actual power to detect a correlation difference at the observed donor counts. Report power as design sensitivity under stated assumptions, not observed-effect post-hoc proof.
7. Reanalyse donor/sample pseudobulk differential expression with appropriate normalization, sex and batch terms, and donor blocking where relevant. Provide complete gene-level results. Replace arbitrary top-300 pathway interpretation with a continuous-statistic gene-set analysis that accounts for inter-gene correlation. Recompute the legacy over-representation results with all eligible gene sets in the testing family, including zero-overlap sets.
8. Assess bulk RNA-seq at day 65 as a separate assay. Quantify donor overlap first. Any donor-disjoint subset must be reported separately; mixed cultures cannot establish dopaminergic cell specificity. Days 0 and 25 can contextualize differentiation, but are not independent validation replicates. Audit protein and scATAC coverage; use them only for questions supported by their real independent sample size and assay identity.

## Deliverables and completion checks

- Input manifest and provenance/cache validation, including exact donor overlap and batch balance.
- Re-runnable scripts, fixed seeds, full result tables, exclusions, test-family definitions, and software versions.
- Figures showing estimates, uncertainty, donor influence, and supported versus unsupported aspects of the claim.
- An evidence-based Chinese analysis report and updated English bioinformatics reviewer Q&A DOCX using the new results. Unperformed manuscript edits and unavailable independent validation remain explicit.
- Validate statistical primitives on synthetic and null examples; check numerical reproduction, matrix/sample alignment, correction scope, and generated-document rendering.
- Finish by stating which claim can be retained, which wording requires qualification, and which part cannot be supported. A failure to reach significance is not proof of no biological effect.

## Initial design findings, before additional outcome testing

- Local cache has 80 scRNA-seq donors. The primary HC/iPD comparison has 8 HC and 30 iPD donors, represented by 41 cultures.
- HC reference donor 3966 has four cultures across batches 1, 3 and 5; these must not be treated as four independent controls.
- Batch 4 has ten iPD cultures and no HC cultures. Batch handling and shared-batch sensitivity are necessary.
- The available bulk RNA count table has 301 assay columns. The accompanying source workbook contains diagnosis, genetic status, age bins and culture batches.
- Original donor-mean correlation, DEG and enrichment outputs have already been inspected; new analyses are explicitly additional and exploratory.

## Metadata amendment before new outcome testing

The source workbook explicitly labels donors 3954 and 4106 as genetic PD (GBA+ and LRRK2+, respectively), while missing subtype fields in the Seurat object previously placed them in `na_iPD`. The primary additional analysis will therefore use a genotype-consistent subset excluding those two donors (8 HC and 28 iPD). The original 8-versus-30 grouping will remain as a named sensitivity analysis. The discrepancy and source rows will be delivered for author confirmation; no genotype is inferred from a missing annotation.

Bulk day-65 has no donor-disjoint healthy controls, although it contains four donor-disjoint iPD donors. It cannot provide an independent HC-versus-iPD replication. It can only provide explicitly labelled cross-assay support using overlapping donors.

## Additional diagnostics after the first new results

The rank-based within-donor sensitivity showed a panel signal while the prespecified RNA/Pearson analysis did not. This prompted diagnostics, not a replacement primary analysis: correction over all eight method/context omnibus tests, conditioning on donor gene-detection fractions and cell counts, and exact permutation within batch-by-sex strata after jointly restricting to shared batches and excluding the reference donor. Reference exclusion alone leaves a saturated batch nuisance stratum; those HC3 estimates are explicitly marked non-estimable rather than assigned a P value.

Bulk day-65 showed an unadjusted SNCA–GABARAPL1 correlation difference. Follow-up retained the complete 43-gene panel and matched bulk and scRNA cultures by source identities. It examined sex/batch adjustment, the scRNA-derived iDA fraction as a composition proxy, leave-one-donor-out effects, and donor-bootstrap intervals with nuisance models refitted in every resample. Conditional slope-interaction tests are labelled separately from correlation-difference tests. These diagnostics are exploratory and do not create an independent validation cohort.
