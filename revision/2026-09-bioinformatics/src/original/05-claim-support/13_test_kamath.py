#!/usr/bin/env python3
"""Frozen donor-level independent analysis of the public Kamath GEO release."""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
INPUT, OUT = HERE / 'inputs', HERE / 'results'
sys.path.insert(0, str(HERE.parent / '04-reviewer-reanalysis'))
from stats_core import fisher, independent_columns, wild_test
from numerics import partial_correlations

SEED, DRAWS = 20260913, 19999
SCENARIOS = ['published_primary', 'published_PMI', 'GEO_clinical_fields',
             'GEO_diagnosis', 'exclude_4775', 'Sepulveda_only']
METHODS = ['rna_qc_library', 'spearman_qc_library']
CONTEXTS = ['DA', 'DA_minus_NonDA']


def donor_design(md, scenario):
    cols = [np.ones(len(md)), md.sex.eq('F').to_numpy(float), (md.age.to_numpy() - 75) / 10]
    names = ['intercept', 'female', 'age_per_10_years']
    bank = pd.get_dummies(md.bank, prefix='bank', drop_first=True)
    cols.extend(bank.to_numpy().T)
    names.extend(bank.columns.tolist())
    if scenario in ['published_PMI', 'GEO_clinical_fields']:
        cols.append((md.PMI.to_numpy() - 20) / 10)
        names.append('PMI_per_10_hours')
    x0, names = independent_columns(np.column_stack(cols), names)
    return np.column_stack([x0, md.group.eq('PD').to_numpy(float)]), names + ['PD']


def scenario_metadata(clinical, scenario):
    md = clinical.copy()
    if scenario == 'GEO_diagnosis':
        md['group'] = md.GEO_group
        md['age'] = md.GEO_age
    elif scenario == 'GEO_clinical_fields':
        md['age'], md['PMI'] = md.GEO_age, md.GEO_PMI
    elif scenario == 'exclude_4775':
        md = md.loc[md.index != '4775']
    elif scenario == 'Sepulveda_only':
        md = md[md.bank.eq('Sepulveda')]
    return md[md.group.isin(['HC', 'PD'])]


def fit(values, md, scenario, context, method, test, forced_failure=None):
    rows, observed, null = [], [], []
    weights = np.random.default_rng(SEED).choice([-1., 1.], size=(DRAWS, len(md)))
    masks = {}
    for gene in values.columns:
        masks.setdefault(tuple(np.isfinite(values[gene])), []).append(gene)
    for mask, genes in masks.items():
        mask = np.array(mask)
        m, y = md.loc[mask], values.loc[mask, genes].to_numpy(float)
        info = dict(scenario=scenario, context=context, method=method, test=test,
                    n_HC=int(m.group.eq('HC').sum()), n_PD=int(m.group.eq('PD').sum()))
        try:
            if forced_failure:
                raise ValueError(forced_failure)
            if min(info['n_HC'], info['n_PD']) < 5:
                raise ValueError('Fewer than five independent donors in an arm')
            x, terms = donor_design(m, scenario)
            info['design'] = ' + '.join(terms)
            f = wild_test(y, x, weights[:, mask])
            for j, gene in enumerate(genes):
                if not np.isfinite(f['t'][j]):
                    rows.append(dict(info, gene=gene, status='Undefined HC3 t statistic', p_wild=np.nan))
                    continue
                assert np.isfinite(f['tstar'][:, j]).all()
                rows.append(dict(info, gene=gene, status='tested', effect_z=f['effect'][j],
                                 ci_low=f['ci_low'][j], ci_high=f['ci_high'][j],
                                 se_hc3=f['se'][j], t_hc3=f['t'][j], df=f['df'],
                                 max_leverage=max(f['leverage']), p_wild=f['p_wild'][j],
                                 p_hc3_t=f['p_hc3_t'][j],
                                 mean_value_HC=np.mean(y[m.group.eq('HC'), j]),
                                 mean_value_PD=np.mean(y[m.group.eq('PD'), j])))
                observed.append(f['t'][j])
                null.append(f['tstar'][:, j])
        except ValueError as exc:
            for gene in genes:
                rows.append(dict(info, gene=gene, status=str(exc), p_wild=np.nan))
    panel = dict(scenario=scenario, context=context, method=method, test=test,
                 genes_tested=len(observed), p_omnibus=np.nan,
                 status='tested' if null else 'non-estimable')
    if null:
        statistic = np.mean(np.square(observed))
        distribution = np.mean(np.column_stack(null) ** 2, axis=1)
        panel.update(statistic_mean_t2=statistic,
                     p_omnibus=(1 + (distribution >= statistic).sum()) / (DRAWS + 1))
    return pd.DataFrame(rows), panel


def adjusted_column(values, method):
    adjusted = multipletests(values.fillna(1), method=method)[1]
    return np.where(values.notna(), adjusted, np.nan)


def main():
    plan_hash = hashlib.sha256((HERE / 'analysis_plan.md').read_bytes()).hexdigest()
    assert plan_hash == (INPUT / 'kamath_pre_outcome_plan_sha256.txt').read_text().split()[0]
    audit = json.loads((INPUT / 'kamath_annotation_mapping_audit.json').read_text())
    assert audit['mapped_cells'] == 387558 and audit['unmatched_cells'] == 0
    assert audit['author_and_matrix_QC_discrepancies'] == 0
    settings = json.loads((OUT / 'foundin_regrouping_settings.json').read_text())
    signs = settings['frozen_direction_signs']
    meta = pd.read_csv(INPUT / 'kamath_cell_annotation_mapping.csv.gz', dtype={'donor': str})
    meta = meta[meta.clusters.isin(['DA', 'NonDA'])].reset_index(drop=True)
    clinical = pd.read_csv(INPUT / 'kamath_published_donor_metadata.csv', dtype={'donor': str}).set_index('donor')
    clinical = clinical.rename(columns={'clinical_group': 'group', 'Sex': 'sex',
                                        'Age (yrs)': 'age', 'PMI (hrs)': 'PMI', 'Brain Bank': 'bank'})
    libraries = pd.read_csv(INPUT / 'kamath_GEO_library_manifest.csv', dtype={'donor_from_unique_library_name_match': str})
    known = libraries[libraries.SNpc_library & libraries.donor_from_unique_library_name_match.notna()]
    for source_field, output in [('status', 'GEO_group'), ('age', 'GEO_age'), ('pmi', 'GEO_PMI')]:
        grouped = known.groupby('donor_from_unique_library_name_match')[source_field]
        assert grouped.nunique().eq(1).all()
        clinical[output] = grouped.first().reindex(clinical.index)
    clinical['GEO_group'] = clinical.GEO_group.replace({'Ctrl': 'HC'})
    clinical.to_csv(OUT / 'kamath_analysis_clinical_metadata.csv')
    coverage = meta.groupby(['donor', 'clusters']).size().rename('n_cells').reset_index()
    coverage['clinical_group'] = coverage.donor.map(clinical.group)
    coverage['passes_50_nuclei'] = coverage.n_cells.ge(50)
    coverage.to_csv(OUT / 'kamath_cell_coverage.csv', index=False)
    primary_donors = coverage.loc[coverage.clusters.eq('DA') & coverage.passes_50_nuclei &
                                 coverage.clinical_group.isin(['HC', 'PD']), 'donor']
    primary_counts = clinical.loc[primary_donors].group.value_counts().to_dict()
    assert primary_counts == {'HC': 8, 'PD': 5}
    with np.load(INPUT / 'kamath_panel_counts.npz', allow_pickle=False) as saved:
        source_genes = saved['genes'].tolist()
        indices = pd.Index(saved['cells']).get_indexer(meta.cell)
        assert (indices >= 0).all() and len(set(indices)) == len(indices)
        wanted = ['SNCA'] + settings['genes']
        raw = saved['counts'][[source_genes.index(g) for g in wanted]][:, indices].T
    counts = pd.DataFrame(raw, columns=wanted)
    selected = meta.clusters.eq('DA') & meta.donor.isin(primary_donors)
    detection = counts.loc[selected, settings['genes']].gt(0).groupby(meta.loc[selected, 'donor']).mean().mean()
    eligibility = detection.rename('mean_primary_donor_detection_DA').rename_axis('gene').reset_index()
    eligibility['eligible'] = eligibility.mean_primary_donor_detection_DA.ge(.05)
    eligibility.to_csv(OUT / 'kamath_gene_eligibility.csv', index=False)
    genes = eligibility.loc[eligibility.eligible, 'gene'].tolist()
    profile_genes = [g for g in signs if g in genes]
    frozen = dict(plan_sha256=plan_hash, primary_donor_counts=primary_counts, genes=genes,
                  profile_genes=profile_genes, profile_signs={g: signs[g] for g in profile_genes},
                  minimum_nuclei=50, minimum_donors_per_arm=5, seed=SEED, draws=DRAWS,
                  scenarios=SCENARIOS, methods=METHODS, contexts=CONTEXTS,
                  planned_edges=len(genes) * 24, planned_panel_and_profile_tests=48,
                  public_release_DA_cells=audit['mapped_DA_nuclei'], published_DA_cells=22048,
                  final_published_cell_set_reproduced=False,
                  normalization='log1p(10000 * raw count / full-matrix library total)',
                  cellular_adjustment='log1p library total, mitochondrial fraction, source library indicators')
    (OUT / 'kamath_analysis_settings.json').write_text(json.dumps(frozen, indent=2) + '\n')
    print('Frozen primary counts', primary_counts, 'eligible panel', len(genes), 'profile', len(profile_genes), flush=True)
    expr = np.log1p(counts[['SNCA'] + genes].to_numpy() / meta.nCount_RNA.to_numpy()[:, None] * 10000)
    records = []
    for (context, donor), m in meta.groupby(['clusters', 'donor']):
        if len(m) < 50:
            continue
        a = expr[m.index]
        nuisance = np.column_stack([np.ones(len(m)), np.log1p(m.nCount_RNA), m.percent_mt / 100,
                                    pd.get_dummies(m.library, drop_first=True).to_numpy()])
        ranked = np.column_stack([stats.rankdata(a[:, j]) for j in range(a.shape[1])])
        rank_nuisance = nuisance.copy()
        rank_nuisance[:, 1] = stats.rankdata(nuisance[:, 1])
        rank_nuisance[:, 2] = stats.rankdata(nuisance[:, 2])
        for method, matrix, covariates in [(METHODS[0], a, nuisance), (METHODS[1], ranked, rank_nuisance)]:
            r = partial_correlations(matrix, covariates)
            for j, gene in enumerate(genes, 1):
                records.append(dict(context=context, donor=donor, clinical_group=clinical.loc[donor, 'group'],
                                    method=method, gene=gene, n_cells=len(m), n_libraries=m.library.nunique(),
                                    r=r[j], z=fisher(r[j])))
    correlations = pd.DataFrame(records)
    correlations.to_csv(OUT / 'kamath_donor_correlations.csv', index=False)
    edges, profiles, tests, scores, donor_values = [], [], [], [], []
    for method in METHODS:
        d = correlations[(correlations.context == 'DA') & correlations.method.eq(method)].pivot(index='donor', columns='gene', values='z')
        n = correlations[(correlations.context == 'NonDA') & correlations.method.eq(method)].pivot(index='donor', columns='gene', values='z')
        common = d.index.intersection(n.index)
        for context, values in [('DA', d), ('DA_minus_NonDA', d.loc[common] - n.loc[common])]:
            for scenario in SCENARIOS:
                md = scenario_metadata(clinical, scenario)
                donors = values.index.intersection(md.index).sort_values()
                md, v = md.loc[donors], values.loc[donors]
                donor_values.append(v.assign(scenario=scenario, context=context, method=method).reset_index())
                result, panel = fit(v, md, scenario, context, method, 'panel')
                edges.append(result)
                tests.append(panel)
                score = v[profile_genes].mul(pd.Series(signs).reindex(profile_genes), axis=1).mean(axis=1, skipna=False)
                scores.append(pd.DataFrame({'donor': donors, 'group': md.group, 'scenario': scenario,
                                            'context': context, 'method': method, 'profile_score': score}))
                failure = 'Fewer than 12 frozen profile genes pass eligibility' if len(profile_genes) < 12 else None
                result, panel = fit(score.to_frame('fixed_direction_profile'), md, scenario, context, method, 'profile', failure)
                profiles.append(result)
                tests.append(panel)
    edges = pd.concat(edges, ignore_index=True)
    assert len(edges) == len(genes) * 24
    edges['q_BH_planned_edges'] = adjusted_column(edges.p_wild, 'fdr_bh')
    edges.to_csv(OUT / 'kamath_all_edges.csv', index=False)
    pd.concat(profiles, ignore_index=True).to_csv(OUT / 'kamath_profile_effects.csv', index=False)
    pd.concat(scores, ignore_index=True).to_csv(OUT / 'kamath_donor_profile_scores.csv', index=False)
    pd.concat(donor_values, ignore_index=True).to_csv(OUT / 'kamath_model_donor_values.csv', index=False)
    tests = pd.DataFrame(tests)
    assert len(tests) == 48
    tests['p_Holm_Kamath_48'] = adjusted_column(tests.p_omnibus, 'holm')
    primary = tests.scenario.eq('published_primary') & tests.method.eq(METHODS[0])
    assert primary.sum() == 4
    tests['p_Holm_primary_4'] = np.nan
    tests.loc[primary, 'p_Holm_primary_4'] = adjusted_column(tests.loc[primary, 'p_omnibus'], 'holm')
    previous = pd.read_csv(OUT / 'gse243639_pilot_panel_tests.csv')
    combined_p = pd.concat([tests.p_omnibus, previous.p_omnibus], ignore_index=True)
    combined_adjusted = adjusted_column(combined_p, 'holm')
    tests['p_Holm_external_64'] = combined_adjusted[:48]
    tests.to_csv(OUT / 'kamath_panel_profile_tests.csv', index=False)
    combined = pd.concat([tests.assign(dataset='GSE178265'), previous.assign(dataset='GSE243639')], ignore_index=True, sort=False)
    combined['p_Holm_external_64'] = combined_adjusted
    combined.to_csv(OUT / 'external_combined_64_tests.csv', index=False)
    print(tests.to_string(index=False), flush=True)
    print('Edges planned / estimable / BH < .05:', len(edges), edges.status.eq('tested').sum(),
          edges.q_BH_planned_edges.lt(.05).sum(), flush=True)


if __name__ == '__main__':
    main()
