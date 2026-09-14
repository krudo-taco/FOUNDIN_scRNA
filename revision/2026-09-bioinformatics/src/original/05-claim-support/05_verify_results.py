#!/usr/bin/env python3
"""Verify analysis families, donor units, independent HC3 fits, and input provenance."""
import ast
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'
sys.path.insert(0, str(HERE.parent / '04-reviewer-reanalysis'))
from stats_core import design


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 ** 2), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    for path in HERE.glob('*.py'):
        ast.parse(path.read_text(), filename=str(path))
    local = pd.read_csv(OUT / 'foundin_regrouped_all_edges.csv')
    assert len(local) == 258
    assert not local[['context', 'method', 'gene']].duplicated().any()
    tested = local.status.eq('tested')
    assert tested.sum() == 172
    assert np.allclose(multipletests(local.p_wild.fillna(1), method='fdr_bh')[1][tested], local.loc[tested, 'q_joint_258_edges'])
    assert local.loc[~tested, 'p_wild'].isna().all()
    donor = pd.read_csv(OUT / 'foundin_regrouped_donor_values.csv')
    checked = []
    for context in ['DA123', 'DA123_minus_remaining_cells']:
        for method in ['rna_qc', 'spearman_qc']:
            d = donor[(donor.context == context) & (donor.method == method)].set_index('PPMI_ID').sort_index()
            assert d.index.is_unique
            assert d.group.value_counts().to_dict() == {'iPD': 28, 'HC': 8}
            for gene in ['SPG11', 'UVRAG', 'GABARAPL1']:
                r = local[(local.context == context) & (local.method == method) & (local.gene == gene)].iloc[0]
                m = d[d[gene].notna()]
                x, names = design(m)
                reference = sm.OLS(m[gene].values, x).fit(cov_type='HC3', use_t=True)
                assert np.isclose(reference.params[-1], r.effect_z, atol=1e-12)
                assert np.isclose(reference.bse[-1], r.se_hc3, atol=1e-12)
                assert np.allclose(reference.conf_int()[-1], [r.ci_low, r.ci_high], atol=1e-12)
                checked.append([context, method, gene])
    settings = json.loads((OUT / 'foundin_regrouping_settings.json').read_text())
    assert len(settings['genes']) == 43 and len(settings['frozen_direction_signs']) == 16
    for _, d in donor.groupby(['context', 'method']):
        expected = d[list(settings['frozen_direction_signs'])].mul(pd.Series(settings['frozen_direction_signs']), axis=1).mean(axis=1, skipna=False)
        assert np.allclose(expected, d.profile_score, equal_nan=True)
    panel = pd.read_csv(OUT / 'foundin_regrouped_panel_tests.csv')
    assert len(panel) == 12
    assert np.allclose(multipletests(panel.p_omnibus.fillna(1), method='holm')[1], panel.p_holm_new_12)
    external = pd.read_csv(OUT / 'gse243639_pilot_all_edges.csv')
    assert len(external) == 264 and external.status.eq('tested').sum() == 92
    assert not external[['context', 'method', 'model', 'gene']].duplicated().any()
    tested_ext = external.status.eq('tested')
    assert np.allclose(multipletests(external.p_wild.fillna(1), method='fdr_bh')[1][tested_ext], external.loc[tested_ext, 'q_joint_pilot_family'])
    assert external.loc[tested_ext, 'n_iPD'].eq(5).all()
    ext_panel = pd.read_csv(OUT / 'gse243639_pilot_panel_tests.csv')
    assert len(ext_panel) == 16
    assert np.allclose(multipletests(ext_panel.p_omnibus.fillna(1), method='holm')[1], ext_panel.p_holm_all_pilot_tests)
    extraction = json.loads((HERE / 'inputs/gse243639_independent_extract_validation.json').read_text())
    assert extraction['matrix_rows'] == 33537 and extraction['genes_with_nonzero_counts'] == 30194
    assert len(extraction['independently_verified_genes']) == 3
    crosswalk = pd.read_csv(HERE / 'sources/foundin_annotation_crosswalk.csv', index_col=0)
    assert crosswalk.values.sum() == 416216
    assert crosswalk.loc[['iDA1', 'iDA2', 'iDA3']].values.sum() == 96623
    assert crosswalk.loc['iDA4'].sum() == 41267
    old = pd.read_csv(HERE.parent / '04-reviewer-reanalysis/results/culture_level_correlations.csv.gz')
    artifacts = old[old.detection.eq(0) & old.r.notna()]
    relevant = artifacts[artifacts.context.isin(['iDA_pooled', 'non_iDA_pooled']) & artifacts.gene.isin(settings['genes'])]
    assert len(relevant) == 0
    kamath_path = HERE / 'inputs/kamath_independent_input_validation.json'
    kamath = None
    if kamath_path.is_file():
        kamath = json.loads(kamath_path.read_text())
        assert kamath['status'] == 'PASS'
        assert kamath['cell_ids_and_order_checked'] == 434340
        assert kamath['independent_line_parser_cells'] == 1000
        for relative_path, expected_hash in kamath['verified_file_sha256'].items():
            assert digest(HERE / relative_path) == expected_hash, relative_path
    kamath_statistics = None
    statistics_path = HERE / 'inputs/kamath_statistical_validation.json'
    if statistics_path.is_file():
        kamath_statistics = json.loads(statistics_path.read_text())
        assert kamath_statistics['status'] == 'PASS'
        assert kamath_statistics['independent_statsmodels_HC3_effect_and_interval_checks'] == 688
    document_validation = None
    document_path = HERE / 'deliverables/document_validation.json'
    if document_path.is_file():
        document_validation = json.loads(document_path.read_text())
        assert document_validation['visual_validation'] == 'PASS'
        assert document_validation['rendered_text_complete']
        assert digest(HERE.parent / document_validation['docx']) == document_validation['sha256']
    report = {
        'status': 'PASS', 'local_planned_edges': len(local), 'local_estimable_edges': int(tested.sum()),
        'local_significant_edges': int((local.q_joint_258_edges < .05).sum()),
        'local_panel_and_profile_tests': len(panel),
        'independent_statsmodels_HC3_checks': checked,
        'profile_scores_match_fixed_gene_directions': True,
        'external_pilot_planned_edges': len(external), 'external_pilot_estimable_edges': int(tested_ext.sum()),
        'external_pilot_significant_edges': int((external.q_joint_pilot_family < .05).sum()),
        'external_panel_and_profile_tests': len(ext_panel),
        'external_input_validation': extraction,
        'kamath_count_and_metadata_validation': kamath,
        'kamath_statistical_validation': kamath_statistics,
        'bioinformatics_QA_document_validation': document_validation,
        'prior_primary_eligible_panel_affected_by_constant_rank_artifact': False,
        'constant_rank_artifact_note': 'Additional analyses guard against numerical residuals from constant rank vectors. In earlier pooled source correlations, affected records were MAP1LC3C, which was excluded from the frozen 43-gene primary panel.',
        'all_python_scripts_parse': True,
        'goal_status': ('Both external datasets and the biologically regrouped FOUNDIN-PD analysis completed; '
                        'strong disease-specific replication not established' if kamath_statistics else
                        'Additional routes tested; Kamath expression validation remains outstanding')
    }
    (HERE / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    paths = (sorted(HERE.glob('*.py')) + sorted(HERE.glob('*.R')) + sorted(HERE.glob('*.md'))
             + [HERE / 'validation.json'] + sorted(OUT.glob('*'))
             + sorted((HERE / 'inputs').glob('*')) + sorted((HERE / 'deliverables').rglob('*')))
    manifest = [{'path': str(p.relative_to(HERE)), 'bytes': p.stat().st_size,
                 'sha256': digest(p)} for p in paths if p.is_file()]
    (HERE / 'artifact_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
