#!/usr/bin/env python3
"""GSE243639 fixed low-coverage pilot. The primary replication coverage gate fails."""
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
from stats_core import anchor_correlations, fisher, independent_columns, residualize, wild_test
from numerics import partial_correlations

SETTINGS = json.loads((OUT / 'foundin_regrouping_settings.json').read_text())
SEED, DRAWS = 20260913, 19999
SIGNS = SETTINGS['frozen_direction_signs']


def donor_design(md, model):
    columns = [np.ones(len(md)), (md.Sex.str.strip() == 'female').astype(float).values,
               (md.Age.values - 75) / 10]
    names = ['intercept', 'female', 'age_per_10_years']
    if model == 'age_sex_PMI_RIN':
        columns.extend([(md['PMI hours'].values - 20) / 10, md['RIN measure'].values - 7])
        names.extend(['PMI_per_10_hours', 'RIN'])
    x0, names = independent_columns(np.column_stack(columns), names)
    return np.column_stack([x0, md.group.eq('iPD').astype(float)]), names + ['iPD']


def fit(values, md, context, method, model, test):
    rows, null_t, observed_t = [], [], []
    weights = np.random.default_rng(SEED).choice([-1., 1.], size=(DRAWS, len(md)))
    masks = {}
    for g in values.columns:
        masks.setdefault(tuple(np.isfinite(values[g])), []).append(g)
    for mask, genes in masks.items():
        mask = np.array(mask)
        m = md.loc[mask]
        y = values.loc[mask, genes].to_numpy(float)
        x, terms = donor_design(m, model)
        info = dict(context=context, method=method, model=model, test=test,
                    n_HC=int(m.group.eq('HC').sum()), n_iPD=int(m.group.eq('iPD').sum()),
                    design=' + '.join(terms))
        try:
            if min(info['n_HC'], info['n_iPD']) < 5:
                raise ValueError('Fewer than five independent donors in an arm')
            f = wild_test(y, x, weights[:, mask])
            for j, gene in enumerate(genes):
                rows.append(dict(info, gene=gene, status='tested', effect_z=f['effect'][j],
                                 ci_low=f['ci_low'][j], ci_high=f['ci_high'][j],
                                 se_hc3=f['se'][j], t_hc3=f['t'][j], df=f['df'],
                                 max_leverage=max(f['leverage']), p_wild=f['p_wild'][j],
                                 p_hc3_t=f['p_hc3_t'][j]))
                observed_t.append(f['t'][j])
                null_t.append(f['tstar'][:, j])
        except ValueError as exc:
            for gene in genes:
                rows.append(dict(info, gene=gene, status=str(exc), p_wild=np.nan))
    panel = dict(context=context, method=method, model=model, test=test,
                 genes_tested=len(observed_t), p_omnibus=np.nan)
    if null_t:
        obs = np.mean(np.square(observed_t))
        null = np.mean(np.column_stack(null_t) ** 2, axis=1)
        panel.update(statistic_mean_t2=obs, p_omnibus=(1 + (null >= obs).sum()) / (DRAWS + 1))
    return pd.DataFrame(rows), panel


def main():
    meta = pd.read_csv(INPUT / 'gse243639_cell_metadata.csv.gz').set_index('CELL_ID')
    qc = pd.read_csv(INPUT / 'gse243639_qc.csv.gz', index_col='cell')
    meta = meta.join(qc)
    counts = pd.read_csv(INPUT / 'gse243639_panel_counts.tsv.gz', sep='\t', index_col='cell')
    neurons = meta[meta.neuron_type.notna()].copy()
    counts = counts.loc[neurons.index]
    expr = np.log1p(counts.div(neurons.nCount_RNA, axis=0) * 10000)
    coverage = neurons.assign(context=np.where(neurons.neuron_type.eq('neurons00'), 'DA', 'other_neurons')).groupby(
        ['donor', 'group', 'context']).size().rename('n_cells').reset_index()
    coverage.to_csv(OUT / 'gse243639_pilot_coverage.csv', index=False)
    primary = coverage[(coverage.context == 'DA') & (coverage.n_cells >= 50)]
    primary_n = primary.groupby('group').donor.nunique().to_dict()
    assert primary_n == {'HC': 7, 'iPD': 2}, primary_n
    eligible_donors = coverage[(coverage.context == 'DA') & (coverage.n_cells >= 20)].donor
    da = neurons[neurons.neuron_type.eq('neurons00') & neurons.donor.isin(eligible_donors)]
    detections = counts.loc[da.index, SETTINGS['genes']].gt(0).groupby(da.donor).mean().mean()
    eligibility = detections.rename('mean_donor_detection_DA').reset_index().rename(columns={'index': 'gene'})
    eligibility['eligible'] = eligibility.mean_donor_detection_DA >= .05
    eligibility.to_csv(OUT / 'gse243639_pilot_gene_eligibility.csv', index=False)
    genes = eligibility.loc[eligibility.eligible, 'gene'].tolist()
    profile_genes = [g for g in SIGNS if g in genes]
    print('Primary coverage gate fails:', primary_n, flush=True)
    print('20-cell pilot panel/profile genes', len(genes), len(profile_genes), flush=True)
    print('Profile genes', profile_genes, flush=True)
    markers = ['TH', 'SLC6A3', 'SLC18A2', 'ALDH1A1', 'GAD1', 'GAD2', 'SNCA']
    counts[markers].gt(0).groupby(neurons.neuron_type).mean().to_csv(OUT / 'gse243639_neuron_marker_detection.csv')
    all_corr = []
    for context, selected in [('DA', neurons.neuron_type.eq('neurons00')),
                              ('other_neurons', ~neurons.neuron_type.eq('neurons00'))]:
        for donor, m in neurons[selected].groupby('donor'):
            n = len(m)
            if n < 20:
                continue
            a = expr.loc[m.index, ['SNCA'] + genes].to_numpy(float)
            cov = np.column_stack([np.ones(n), np.log1p(m.nCount_RNA), m.percent_mt / 100,
                                   pd.get_dummies(m.neuron_type, drop_first=True).values])
            ranks = np.column_stack([stats.rankdata(a[:, j]) for j in range(a.shape[1])])
            rcov = cov.copy()
            rcov[:, 1] = stats.rankdata(cov[:, 1])
            rcov[:, 2] = stats.rankdata(cov[:, 2])
            for method, matrix, nuisance in [('rna_qc', a, cov), ('spearman_qc', ranks, rcov)]:
                r = partial_correlations(matrix, nuisance)
                for j, gene in enumerate(genes, 1):
                    all_corr.append(dict(context=context, donor=donor, group=m.iloc[0].group,
                                         method=method, gene=gene, n_cells=n, r=r[j], z=fisher(r[j])))
    corr = pd.DataFrame(all_corr)
    corr.to_csv(OUT / 'gse243639_pilot_donor_correlations.csv', index=False)
    clinical = pd.read_csv(INPUT / 'gse243639_donors.csv').set_index('donor')
    edges, tests, profiles, donor_scores = [], [], [], []
    for method in ['rna_qc', 'spearman_qc']:
        d = corr[(corr.context == 'DA') & (corr.method == method)].pivot(index='donor', columns='gene', values='z')
        n = corr[(corr.context == 'other_neurons') & (corr.method == method)].pivot(index='donor', columns='gene', values='z')
        common = d.index.intersection(n.index)
        for context, values in [('DA', d), ('DA_minus_other_neurons', d.loc[common] - n.loc[common])]:
            md = clinical.loc[values.index]
            score = values[profile_genes].mul(pd.Series(SIGNS).reindex(profile_genes), axis=1).mean(axis=1, skipna=False)
            if len(profile_genes) < 12:
                score[:] = np.nan
            named = pd.DataFrame({'manuscript_direction_score': score})
            ds = md[['group', 'Age', 'Sex', 'PMI hours', 'RIN measure']].copy()
            ds['score'] = score
            ds['context'], ds['method'] = context, method
            donor_scores.append(ds.reset_index())
            for model in ['age_sex', 'age_sex_PMI_RIN']:
                r, panel = fit(values, md, context, method, model, 'full_panel')
                edges.append(r)
                tests.append(panel)
                r, panel = fit(named, md, context, method, model, 'frozen_direction_profile')
                r['profile_genes'] = len(profile_genes)
                profiles.append(r)
                tests.append(panel)
    edges = pd.concat(edges, ignore_index=True)
    mask = edges.p_wild.notna()
    adjusted = multipletests(edges.p_wild.fillna(1), method='fdr_bh')[1]
    edges['q_joint_pilot_family'] = np.where(mask, adjusted, np.nan)
    edges.to_csv(OUT / 'gse243639_pilot_all_edges.csv', index=False)
    pd.concat(profiles, ignore_index=True).to_csv(OUT / 'gse243639_pilot_profile_effects.csv', index=False)
    pd.concat(donor_scores, ignore_index=True).to_csv(OUT / 'gse243639_pilot_donor_scores.csv', index=False)
    tests = pd.DataFrame(tests)
    tests['p_holm_all_pilot_tests'] = multipletests(tests.p_omnibus.fillna(1), method='holm')[1]
    tests.to_csv(OUT / 'gse243639_pilot_panel_tests.csv', index=False)
    (OUT / 'gse243639_pilot_settings.json').write_text(json.dumps({
        'primary_minimum_nuclei': 50, 'primary_donor_counts': primary_n, 'primary_replication_gate': 'FAIL',
        'pilot_minimum_nuclei': 20, 'eligible_genes': genes,
        'profile_genes': profile_genes, 'profile_signs': {g: SIGNS[g] for g in profile_genes},
        'seed': SEED, 'draws': DRAWS,
        'interpretation': 'Exploratory low-coverage pilot; insufficient to establish independent replication alone',
        'normalization': 'log(1 + 10000 * raw gene count / total count from supplied matrix)',
        'nuclei_used': 'Author-provided neuronal subclusters only; no post hoc reassignment of excluded neurons'
    }, indent=2) + '\n')
    print(tests.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
