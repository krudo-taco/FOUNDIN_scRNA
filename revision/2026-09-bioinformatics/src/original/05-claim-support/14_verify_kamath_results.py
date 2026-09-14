#!/usr/bin/env python3
"""Independent estimation checks and finite-support bootstrap diagnostics."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

HERE = Path(__file__).resolve().parent
INPUT, OUT = HERE / 'inputs', HERE / 'results'


def adjust(p, method):
    values = np.asarray(p.fillna(1), float)
    order = np.argsort(values)
    if method == 'holm':
        q = np.maximum.accumulate(values[order] * np.arange(len(values), 0, -1))
    else:
        q = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty(len(values))
    result[order] = np.minimum(q, 1)
    return np.where(p.notna(), result, np.nan)


def reference_metadata(clinical, scenario, donors):
    md = clinical.copy()
    if scenario == 'GEO_diagnosis':
        md['group'], md['age'] = md.GEO_group, md.GEO_age
    elif scenario == 'GEO_clinical_fields':
        md['age'], md['PMI'] = md.GEO_age, md.GEO_PMI
    elif scenario == 'exclude_4775':
        md = md.drop('4775')
    elif scenario == 'Sepulveda_only':
        md = md[md.bank.eq('Sepulveda')]
    md = md[md.group.isin(['HC', 'PD'])].loc[donors]
    return md


def reference_design(md, scenario):
    columns = [np.ones(len(md)), np.array([float(s == 'F') for s in md.sex]), (md.age.to_numpy() - 75) / 10]
    for b in sorted(md.bank.unique())[1:]:
        columns.append(np.array([float(s == b) for s in md.bank]))
    if scenario in ['published_PMI', 'GEO_clinical_fields']:
        columns.append((md.PMI.to_numpy() - 20) / 10)
    columns.append(np.array([float(s == 'PD') for s in md.group]))
    return np.asarray(columns).T


def bootstrap_reference(y, x, weights):
    null_fit = sm.OLS(y, x[:, :-1]).fit()
    h0 = null_fit.get_influence().hat_matrix_diag
    responses = null_fit.fittedvalues[None, :] + weights * (null_fit.resid / np.sqrt(1 - h0))[None, :]
    inverse = np.linalg.pinv(x)
    beta = inverse @ responses.T
    residual = responses.T - x @ beta
    leverage = np.diag(x @ inverse)
    se = np.sqrt(((inverse[-1, :, None] * residual / (1 - leverage[:, None])) ** 2).sum(axis=0))
    return beta[-1] / se


def main():
    settings = json.loads((OUT / 'kamath_analysis_settings.json').read_text())
    assert settings['plan_sha256'] == hashlib.sha256((HERE / 'analysis_plan.md').read_bytes()).hexdigest()
    edges = pd.read_csv(OUT / 'kamath_all_edges.csv')
    profiles = pd.read_csv(OUT / 'kamath_profile_effects.csv')
    tests = pd.read_csv(OUT / 'kamath_panel_profile_tests.csv')
    values = pd.read_csv(OUT / 'kamath_model_donor_values.csv', dtype={'donor': str})
    scores = pd.read_csv(OUT / 'kamath_donor_profile_scores.csv', dtype={'donor': str})
    clinical = pd.read_csv(OUT / 'kamath_analysis_clinical_metadata.csv', dtype={'donor': str}).set_index('donor')
    assert len(edges) == 42 * 24 == settings['planned_edges']
    assert not edges[['scenario', 'context', 'method', 'gene']].duplicated().any()
    assert edges.status.eq('tested').sum() == 672
    assert len(profiles) == 24 and profiles.status.eq('tested').sum() == 16
    assert len(tests) == 48 and tests.status.eq('tested').sum() == 32
    failed = edges[edges.scenario.isin(['GEO_diagnosis', 'exclude_4775'])]
    assert failed.n_PD.eq(4).all() and failed.p_wild.isna().all()
    assert np.allclose(adjust(edges.p_wild, 'bh'), edges.q_BH_planned_edges, equal_nan=True)
    assert np.allclose(adjust(tests.p_omnibus, 'holm'), tests.p_Holm_Kamath_48, equal_nan=True)
    primary = tests.scenario.eq('published_primary') & tests.method.eq('rna_qc_library')
    assert np.allclose(adjust(tests.loc[primary, 'p_omnibus'], 'holm'), tests.loc[primary, 'p_Holm_primary_4'])
    combined = pd.read_csv(OUT / 'external_combined_64_tests.csv')
    assert len(combined) == 64
    assert np.allclose(adjust(combined.p_omnibus, 'holm'), combined.p_Holm_external_64, equal_nan=True)

    checked, snapshots = 0, {}
    for key, d in values.groupby(['scenario', 'context', 'method']):
        scenario, context, method = key
        d = d.set_index('donor').sort_index()
        assert d.index.is_unique
        md = reference_metadata(clinical, scenario, d.index)
        s = scores[scores.scenario.eq(scenario) & scores.context.eq(context) & scores.method.eq(method)].set_index('donor').loc[d.index]
        expected_score = sum(d[g] * sign for g, sign in settings['profile_signs'].items()) / len(settings['profile_signs'])
        assert np.allclose(expected_score, s.profile_score, equal_nan=True)
        selected = edges[edges.scenario.eq(scenario) & edges.context.eq(context) & edges.method.eq(method) & edges.status.eq('tested')]
        selected_profiles = profiles[profiles.scenario.eq(scenario) & profiles.context.eq(context) & profiles.method.eq(method) & profiles.status.eq('tested')]
        for _, row in pd.concat([selected, selected_profiles], ignore_index=True).iterrows():
            y = s.profile_score if row.test == 'profile' else d[row.gene]
            valid = y.notna()
            m, y = md[valid], y[valid]
            x = reference_design(m, scenario)
            fitted = sm.OLS(y.to_numpy(), x).fit(cov_type='HC3', use_t=True)
            assert row.n_HC == m.group.eq('HC').sum() and row.n_PD == m.group.eq('PD').sum()
            assert np.isclose(fitted.params[-1], row.effect_z, rtol=0, atol=1e-12)
            assert np.isclose(fitted.bse[-1], row.se_hc3, rtol=0, atol=1e-12)
            assert np.allclose(fitted.conf_int()[-1], [row.ci_low, row.ci_high], rtol=0, atol=1e-12)
            assert np.isclose(fitted.pvalues[-1], row.p_hc3_t, rtol=0, atol=1e-12)
            checked += 1
            if (scenario == 'Sepulveda_only' and context == 'DA_minus_NonDA' and
                ((row.gene == 'TOLLIP' and method == 'rna_qc_library') or
                 (row.gene == 'ATG4B' and method == 'spearman_qc_library'))):
                snapshots[(row.gene, method)] = (row, y.to_numpy(), x)

    # Independently rebuild selected cellular partial correlations from raw counts.
    meta = pd.read_csv(INPUT / 'kamath_cell_annotation_mapping.csv.gz', dtype={'donor': str})
    correlations = pd.read_csv(OUT / 'kamath_donor_correlations.csv', dtype={'donor': str})
    with np.load(INPUT / 'kamath_panel_counts.npz', allow_pickle=False) as saved:
        all_genes, all_cells, raw = saved['genes'].tolist(), pd.Index(saved['cells']), saved['counts']
    partial_checks = 0
    for donor in ['3298', '3887', '4956']:
        for context in ['DA', 'NonDA']:
            m = meta[meta.donor.eq(donor) & meta.clusters.eq(context)]
            cell_ids = all_cells.get_indexer(m.cell)
            for gene in ['SPG11', 'UVRAG']:
                a = np.log1p(raw[[all_genes.index('SNCA'), all_genes.index(gene)]][:, cell_ids].T / m.nCount_RNA.to_numpy()[:, None] * 10000)
                cov = np.column_stack([np.ones(len(m)), np.log1p(m.nCount_RNA), m.percent_mt / 100,
                                       pd.get_dummies(m.library, drop_first=True).to_numpy()])
                for method in ['rna_qc_library', 'spearman_qc_library']:
                    matrix, nuisance = a.copy(), cov.copy()
                    if method.startswith('spearman'):
                        matrix = np.column_stack([stats.rankdata(matrix[:, j]) for j in range(2)])
                        nuisance[:, 1] = stats.rankdata(nuisance[:, 1])
                        nuisance[:, 2] = stats.rankdata(nuisance[:, 2])
                    residuals = [sm.OLS(matrix[:, j], nuisance).fit().resid for j in range(2)]
                    expected_r = stats.pearsonr(*residuals)[0]
                    row = correlations[correlations.donor.eq(donor) & correlations.context.eq(context) & correlations.gene.eq(gene) & correlations.method.eq(method)].iloc[0]
                    assert np.isclose(expected_r, row.r, rtol=0, atol=2e-11)
                    partial_checks += 1

    diagnostics = []
    for (gene, method), (row, y, x) in snapshots.items():
        random_weights = np.random.default_rng(settings['seed']).choice([-1., 1.], size=(settings['draws'], len(y)))
        reference_t = bootstrap_reference(y, x, random_weights)
        observed_t = sm.OLS(y, x).fit(cov_type='HC3', use_t=True).tvalues[-1]
        exceedances = int((np.abs(reference_t) >= abs(observed_t)).sum())
        assert np.isclose((1 + exceedances) / (settings['draws'] + 1), row.p_wild, atol=1e-12)
        patterns = 2 ** len(y)
        signs = 2 * ((np.arange(patterns)[:, None] >> np.arange(len(y))) & 1) - 1
        enumerated_t = bootstrap_reference(y, x, signs)
        enumerated_exceedances = int((np.abs(enumerated_t) >= abs(observed_t)).sum())
        diagnostics.append(dict(gene=gene, scenario=row.scenario, context=row.context, method=method,
                                n_donors=len(y), observed_t=observed_t, MC_draws=settings['draws'],
                                MC_tail_exceedances=exceedances, MC_P=row.p_wild,
                                finite_sign_patterns=patterns, enumerated_tail_exceedances=enumerated_exceedances,
                                largest_abs_enumerated_t=float(np.abs(enumerated_t).max()),
                                enumerated_tail_fraction=enumerated_exceedances / patterns,
                                enumerated_add_one_P=(1 + enumerated_exceedances) / (patterns + 1),
                                HC3_t_P=row.p_hc3_t, HC3_ci_low=row.ci_low, HC3_ci_high=row.ci_high))
    # Preserve all original P values. This bounded arithmetic sensitivity changes
    # only the two zero-tail floors, not the prespecified primary analysis.
    annotated = edges.copy()
    annotated['p_wild_floor_sensitivity'] = annotated.p_wild
    annotated['bootstrap_caution'] = ''
    for item in diagnostics:
        matched = (annotated.gene.eq(item['gene']) & annotated.method.eq(item['method']) &
                   annotated.context.eq(item['context']) & annotated.scenario.eq(item['scenario']))
        assert matched.sum() == 1
        annotated.loc[matched, 'p_wild_floor_sensitivity'] = item['enumerated_add_one_P']
        annotated.loc[matched, 'bootstrap_caution'] = 'Zero tail exceedances in MC and all 2048 Rademacher patterns; reported MC P is an add-one floor'
    annotated['q_BH_floor_sensitivity'] = adjust(annotated.p_wild_floor_sensitivity, 'bh')
    annotated.to_csv(OUT / 'kamath_all_edges_with_QA_flags.csv', index=False)
    for item in diagnostics:
        row = annotated[annotated.gene.eq(item['gene']) & annotated.method.eq(item['method']) &
                        annotated.context.eq(item['context']) & annotated.scenario.eq(item['scenario'])].iloc[0]
        item['original_MC_BH_q'] = row.q_BH_planned_edges
        item['BH_q_under_two_floor_sensitivity'] = row.q_BH_floor_sensitivity
    diagnostic = pd.DataFrame(diagnostics)
    diagnostic.to_csv(OUT / 'kamath_bootstrap_floor_diagnostic.csv', index=False)
    report = {
        'status': 'PASS', 'scope': 'Count mapping, independent coefficient/interval calculations, multiple-testing families and selected bootstrap tails; not proof of biological replication',
        'independent_statsmodels_HC3_effect_and_interval_checks': checked,
        'independent_cellular_partial_correlation_checks': partial_checks,
        'planned_edge_tests': len(edges), 'estimable_edge_tests': int(edges.status.eq('tested').sum()),
        'planned_panel_profile_tests': len(tests), 'estimable_panel_profile_tests': int(tests.status.eq('tested').sum()),
        'signed_scores_match_frozen_directions': True, 'corrections_independently_recomputed': True,
        'primary_RNA_test_P': tests.loc[primary, ['context', 'test', 'p_omnibus', 'p_Holm_primary_4']].to_dict('records'),
        'BH_significant_secondary_edges_under_original_MC_rule': int(edges.q_BH_planned_edges.lt(.05).sum()),
        'BH_significant_edges_under_two_floor_sensitivity': int(annotated.q_BH_floor_sensitivity.lt(.05).sum()),
        'bootstrap_floor_diagnostics': diagnostics,
        'warning': 'The two smallest Monte Carlo P values have zero exceedances. Exhaustive Rademacher evaluation is a finite-support diagnostic, not an exact randomized disease-label test. Do not present these secondary results as robust independent replication.'
    }
    (INPUT / 'kamath_statistical_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
