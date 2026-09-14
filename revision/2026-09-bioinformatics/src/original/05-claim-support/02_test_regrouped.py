#!/usr/bin/env python3
"""Frozen post hoc biological regrouping of FOUNDIN-PD with complete multiplicity accounting."""
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / '04-reviewer-reanalysis'
OUT = HERE / 'results'
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(PREVIOUS))
from stats_core import anchor_correlations, bh, fisher, residualize
from numerics import partial_correlations

spec = importlib.util.spec_from_file_location('previous_coupling', PREVIOUS / '03_within_donor_coupling.py')
previous = importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)
SETTINGS = json.loads((PREVIOUS / 'results/analysis_settings.json').read_text())
GENES = SETTINGS['eligible_genes']
POSITIVE = ['ATG5', 'BCL2', 'UVRAG', 'ATG10', 'ATG7', 'MAP1LC3C', 'MAP1LC3A']
NEGATIVE = ['OPTN', 'EPG5', 'SPG11', 'RB1CC1', 'PINK1', 'SNX4', 'SNX30', 'SNX14', 'CALCOCO2', 'WIPI2', 'ATG16L2']
SIGNS = {g: (1 if g in POSITIVE else -1) for g in POSITIVE + NEGATIVE if g in GENES}
METHODS = ['rna_qc', 'spearman_qc']


def calculate():
    meta = pd.read_csv(PREVIOUS / 'inputs/cell_metadata_full.tsv.gz', sep='\t', index_col='cell')
    meta = meta[meta.subtype.isin(['na_HC', 'na_iPD']) & ~meta.PPMI_ID.isin([3954, 4106])].copy()
    meta['group'] = meta.subtype.map({'na_HC': 'HC', 'na_iPD': 'iPD'})
    expr = pd.read_csv(PREVIOUS / 'inputs/rna_panel_logexpr.tsv.gz', sep='\t', index_col='cell').loc[meta.index, ['SNCA'] + GENES]
    assert len(SIGNS) == 16
    rows, coverage = [], []
    for context in ['DA123', 'immature4', 'remaining_cells']:
        is_da = meta.CellType.isin(['iDA1', 'iDA2', 'iDA3'])
        selected = is_da if context == 'DA123' else meta.CellType.eq('iDA4') if context == 'immature4' else ~is_da
        for sample, sm in meta[selected].groupby('SampleID'):
            first, n = sm.iloc[0], len(sm)
            identity = dict(context=context, sample=sample, PPMI_ID=int(first.PPMI_ID), group=first.group,
                            genetic_sex=int(first.genetic_sex), BATCH=int(first.BATCH), n_cells=n)
            assert sm.PPMI_ID.nunique() == sm.genetic_sex.nunique() == sm.BATCH.nunique() == sm.group.nunique() == 1
            coverage.append(dict(identity, included=n >= 100))
            if n < 100:
                continue
            a = expr.loc[sm.index].to_numpy(float)
            cov = np.column_stack([np.ones(n), np.log1p(sm.nCount_RNA), sm['percent.mt'] / 100,
                                   pd.get_dummies(sm.CellType, drop_first=True).values])
            ranked = np.column_stack([stats.rankdata(a[:, j]) for j in range(a.shape[1])])
            rcov = cov.copy()
            rcov[:, 1] = stats.rankdata(cov[:, 1])
            rcov[:, 2] = stats.rankdata(cov[:, 2])
            for method, matrix, nuisance in [('rna_qc', a, cov), ('spearman_qc', ranked, rcov)]:
                correlations = partial_correlations(matrix, nuisance)
                for j, gene in enumerate(GENES, 1):
                    rows.append(dict(identity, method=method, gene=gene, r=float(correlations[j]),
                                     z=float(fisher(correlations[j]))))
        print(time.strftime('%H:%M:%S'), 'Calculated', context, flush=True)
    corr, cov = pd.DataFrame(rows), pd.DataFrame(coverage)
    corr.to_csv(OUT / 'foundin_regrouped_culture_correlations.csv.gz', index=False)
    cov.to_csv(OUT / 'foundin_regrouped_culture_coverage.csv', index=False)
    return corr, cov


def paired(corr, coverage, other, label):
    a = corr[corr.context == 'DA123'].copy()
    b = corr[corr.context == other]
    joined = a.merge(b[['sample', 'method', 'gene', 'z']], on=['sample', 'method', 'gene'], suffixes=('', '_other'), validate='one_to_one')
    joined['z'] = joined.z - joined.z_other
    joined['r'] = np.nan
    joined['context'] = label
    used = set(joined['sample'])
    cov = coverage[(coverage.context == 'DA123') & coverage['sample'].isin(used)].copy()
    cov['context'] = label
    return joined.drop(columns='z_other'), cov


def main():
    corr, coverage = calculate()
    contrasts = ['DA123', 'DA123_minus_immature4', 'DA123_minus_remaining_cells']
    for other, label in [('immature4', contrasts[1]), ('remaining_cells', contrasts[2])]:
        c, v = paired(corr, coverage, other, label)
        corr = pd.concat([corr, c], ignore_index=True)
        coverage = pd.concat([coverage, v], ignore_index=True)
    all_edges, panels, profiles, all_values = [], [], [], []
    for method in METHODS:
        for context in contrasts:
            values, md = previous.aggregate(corr, coverage, context, method)
            estimate, panel = previous.estimate(values, md, context, method, GENES)
            estimate = estimate.drop(columns=['mean_r_HC', 'mean_r_iPD'], errors='ignore')
            all_edges.append(estimate)
            panels.append(dict(panel, test='full_panel'))
            score = values[list(SIGNS)].mul(pd.Series(SIGNS), axis=1).mean(axis=1, skipna=False)
            named = pd.DataFrame({'manuscript_direction_score': score})
            res, unused = previous.estimate(named, md, context, method, list(named.columns))
            res = res.drop(columns=['mean_r_HC', 'mean_r_iPD'], errors='ignore')
            res['profile_genes'] = len(SIGNS)
            profiles.append(res)
            row = res.iloc[0]
            panels.append(dict(context=context, method=method, test='frozen_direction_profile',
                               genes_tested=len(SIGNS), p_omnibus=row.p_wild,
                               statistic_mean_t2=row.t_hc3 ** 2 if row.status == 'tested' else np.nan))
            dv = values.copy()
            dv['profile_score'] = score
            dv = dv.join(md)
            dv['context'] = context
            dv['method'] = method
            all_values.append(dv.reset_index())
            print(time.strftime('%H:%M:%S'), 'Fit', context, method, 'panel P', panel['p_omnibus'],
                  'profile effect/P', row.get('effect_z'), row.p_wild, flush=True)
    edges = pd.concat(all_edges, ignore_index=True)
    # Count all 258 planned tests, assigning P=1 only inside the adjustment for
    # non-estimable tests; keep their displayed estimates and q values missing.
    planned_q = multipletests(edges.p_wild.fillna(1), method='fdr_bh')[1]
    edges['q_joint_258_edges'] = np.where(edges.p_wild.notna(), planned_q, np.nan)
    edges.to_csv(OUT / 'foundin_regrouped_all_edges.csv', index=False)
    pd.concat(profiles, ignore_index=True).to_csv(OUT / 'foundin_regrouped_profile_effects.csv', index=False)
    pd.concat(all_values, ignore_index=True).to_csv(OUT / 'foundin_regrouped_donor_values.csv', index=False)
    panel = pd.DataFrame(panels)
    assert len(panel) == 12
    panel['p_holm_new_12'] = multipletests(panel.p_omnibus.fillna(1), method='holm')[1]
    old = pd.read_csv(PREVIOUS / 'results/omnibus_multipipeline_sensitivity.csv')
    old = old[old.context.isin(['iDA_pooled', 'iDA_specificity'])]
    assert len(old) == 8
    combined = np.r_[panel.p_omnibus.fillna(1), old.p_omnibus]
    panel['p_holm_including_previous_8'] = multipletests(combined, method='holm')[1][:12]
    panel.to_csv(OUT / 'foundin_regrouped_panel_tests.csv', index=False)
    (OUT / 'foundin_regrouping_settings.json').write_text(json.dumps({
        'genes': GENES, 'frozen_direction_signs': SIGNS, 'contrasts': contrasts,
        'seed': previous.SEED, 'wild_draws': previous.N_WILD, 'cell_threshold': 100,
        'methods': METHODS, 'post_hoc': True, 'independent_validation': False,
        'cell_group_mapping_status': 'Published group inferred from exact cell totals, supported by archived earlier annotations',
        'paired_contrast_rule': 'Compute within identical cultures before equal culture weighting within donor'
    }, indent=2) + '\n')
    print(panel.to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
