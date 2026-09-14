#!/usr/bin/env python3
"""Independently check extracted counts; do not infer cell types or test disease effects."""
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXT, OUT = HERE / 'external', HERE / 'inputs'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 ** 2), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    download = json.loads((OUT / 'kamath_matrix_download_audit.json').read_text())
    extraction = json.loads((OUT / 'kamath_count_extraction_audit.json').read_text())
    matrix = EXT / 'GSE178265_Homo_matrix.mtx.gz'
    assert matrix.stat().st_size == download['bytes']
    matrix_hash = digest(matrix)
    assert matrix_hash == download['sha256']
    features = pd.read_csv(EXT / 'GSE178265_Homo_features.tsv.gz', sep='\t', header=None,
                           names=['feature_id', 'gene', 'type'], dtype=str)
    cells = pd.read_csv(EXT / 'GSE178265_Homo_bcd.tsv.gz', header=None, dtype=str)[0].to_numpy()
    with np.load(OUT / 'kamath_panel_counts.npz', allow_pickle=False) as saved:
        genes, panel, panel_cells = saved['genes'], saved['counts'], saved['cells']
    qc = pd.read_csv(OUT / 'kamath_cell_qc.csv.gz')
    totals = pd.read_csv(OUT / 'kamath_feature_audit.csv')
    assert len(set(cells)) == len(cells) == 434340
    assert np.array_equal(cells, panel_cells)
    assert np.array_equal(cells, qc.cell.to_numpy())
    assert panel.shape == (70, len(cells)) and np.issubdtype(panel.dtype, np.integer)
    assert (panel >= 0).all() and len(set(genes)) == len(genes)
    assert np.array_equal(features.to_numpy(), totals[['feature_id', 'gene', 'type']].to_numpy())
    canonical = features.gene.replace({'PARK2': 'PRKN'})
    assert np.array_equal(canonical.to_numpy(), totals.canonical_panel_symbol.to_numpy())
    for i, gene in enumerate(genes):
        assert panel[i].sum() == totals.loc[canonical.eq(gene), 'total_UMIs'].sum()
    assert qc.nCount_RNA.sum() == totals.total_UMIs.sum() == extraction['total_UMIs'] == 5738238594
    assert qc.nFeature_RNA.sum() == extraction['matrix_nonzero_entries'] == 1654688953
    assert qc.nCount_RNA.gt(0).all() and qc.nFeature_RNA.between(1, len(features)).all()
    assert qc.percent_mt.between(0, 100).all()
    mito = features.gene.str.upper().str.startswith('MT-').to_numpy()
    recovered_mt = qc.percent_mt.to_numpy() * qc.nCount_RNA.to_numpy() / 100
    assert np.allclose(recovered_mt, np.rint(recovered_mt), rtol=0, atol=1e-7)
    assert int(np.rint(recovered_mt).sum()) == int(totals.loc[mito, 'total_UMIs'].sum())

    # This reference parser uses Python integers and one text line at a time;
    # it shares neither the numpy chunk parser nor the extraction aggregators.
    n_reference_cells = 1000
    reference = np.zeros((len(features), n_reference_cells), dtype=np.int64)
    entries = 0
    with gzip.open(matrix, 'rt') as stream:
        assert stream.readline().strip().lower() == '%%matrixmarket matrix coordinate integer general'
        line = stream.readline()
        while line.startswith('%'):
            line = stream.readline()
        assert tuple(map(int, line.split())) == (len(features), len(cells), extraction['matrix_nonzero_entries'])
        for line in stream:
            fields = line.split()
            assert len(fields) == 3
            row, col, value = map(int, fields)
            if col > n_reference_cells:
                break
            assert reference[row - 1, col - 1] == 0 and value > 0
            reference[row - 1, col - 1] = value
            entries += 1
    assert np.array_equal(reference.sum(axis=0), qc.nCount_RNA.to_numpy()[:n_reference_cells])
    assert np.array_equal((reference > 0).sum(axis=0), qc.nFeature_RNA.to_numpy()[:n_reference_cells])
    expected_mt = reference[mito].sum(axis=0) / reference.sum(axis=0) * 100
    assert np.allclose(expected_mt, qc.percent_mt.to_numpy()[:n_reference_cells], rtol=0, atol=1e-12)
    for i, gene in enumerate(genes):
        assert np.array_equal(reference[canonical.eq(gene)].sum(axis=0), panel[i, :n_reference_cells])

    spec = importlib.util.spec_from_file_location('kamath_extractor', HERE / '07_extract_kamath_counts.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    expected = np.array([[1, 1, 23], [35, 1, 1], [6, 2, 789], [345, 2, 2], [3, 3, 5]])
    encoded = '\n'.join(' '.join(map(str, row)) for row in expected).encode()
    for ending in [b'', b'\n']:
        for block_bytes in [1, 2, 7, 17, 1024]:
            result = np.concatenate(list(module.triplet_chunks(io.BytesIO(encoded + ending), block_bytes)))
            assert np.array_equal(result, expected)

    clinical = pd.read_csv(OUT / 'kamath_published_donor_metadata.csv', dtype={'donor': str})
    libraries = pd.read_csv(OUT / 'kamath_GEO_library_manifest.csv', dtype={'donor_from_unique_library_name_match': str})
    actual_library_counts = pd.Series(cells).str.rsplit('_', n=1).str[0].value_counts().sort_index()
    saved_library_counts = libraries.set_index('library').n_barcodes.sort_index()
    assert np.array_equal(actual_library_counts.index, saved_library_counts.index)
    assert np.array_equal(actual_library_counts.to_numpy(), saved_library_counts.to_numpy())
    known = libraries[libraries.SNpc_library & libraries.donor_from_unique_library_name_match.notna()]
    assert set(known.donor_from_unique_library_name_match) == set(clinical.donor)
    assert len(known) == 89 and known.n_barcodes.sum() == 387558
    assert clinical.published_DA_nuclei.sum() == 22048
    # These are coverage scenarios from published counts, not measured DA-cell
    # membership in the downloaded matrix and not biological hypothesis tests.
    scenarios = []
    for name in ['published_table_primary', 'GEO_diagnosis_sensitivity', 'exclude_4775', 'Sepulveda_only']:
        d = clinical.copy()
        if name == 'GEO_diagnosis_sensitivity':
            diagnosis = known.groupby('donor_from_unique_library_name_match').status.agg(lambda x: sorted(set(x)))
            assert diagnosis.map(len).eq(1).all()
            d['clinical_group'] = d.donor.map(diagnosis.map(lambda x: x[0])).replace({'Ctrl': 'HC'})
            assert d.clinical_group.notna().all()
        elif name == 'exclude_4775':
            d = d[d.donor.ne('4775')]
        elif name == 'Sepulveda_only':
            d = d[d['Brain Bank'].eq('Sepulveda')]
        d = d[d.published_DA_nuclei.ge(50)]
        counts = d.clinical_group.value_counts()
        scenarios.append({'scenario': name, 'n_HC': int(counts.get('HC', 0)),
                          'n_PD': int(counts.get('PD', 0)), 'n_LBD': int(counts.get('LBD', 0)),
                          'PD_vs_HC_passes_published_coverage_gate': bool(counts.get('HC', 0) >= 5 and counts.get('PD', 0) >= 5)})
    pd.DataFrame(scenarios).to_csv(OUT / 'kamath_published_coverage_scenarios.csv', index=False)
    report = {
        'status': 'PASS', 'matrix_sha256_matches_download': True,
        'cell_ids_and_order_checked': len(cells), 'panel_gene_totals_checked': len(genes),
        'whole_matrix_QC_and_gene_totals_agree': True,
        'independent_line_parser_cells': n_reference_cells,
        'independent_line_parser_nonzero_entries': entries,
        'independent_line_parser_all_panel_counts_and_QC_match': True,
        'chunk_boundary_and_missing_final_newline_cases': 10,
        'known_clinical_SNpc_libraries': len(known), 'known_clinical_SNpc_matrix_cells': int(known.n_barcodes.sum()),
        'published_SNpc_cells': 387483, 'unreconciled_cell_count_difference': int(known.n_barcodes.sum() - 387483),
        'published_DA_nuclei': int(clinical.published_DA_nuclei.sum()),
        'coverage_scenarios_from_published_counts': scenarios,
        'limitation': 'Independent direct count comparison covers the first 1000 cells; global totals cover all cells. Count verification does not establish author cell types or final published cell membership; consult the separate annotation and outcome audits.',
        'disease_effect_testing': 'Not assessed by count verification; see separate annotation and statistical audits'
    }
    report['verified_file_sha256'] = {
        str(p.relative_to(HERE)): digest(p) for p in [
            EXT / 'GSE178265_Homo_features.tsv.gz', EXT / 'GSE178265_Homo_bcd.tsv.gz',
            OUT / 'kamath_panel_counts.npz', OUT / 'kamath_cell_qc.csv.gz',
            OUT / 'kamath_feature_audit.csv', OUT / 'kamath_GEO_library_manifest.csv',
            OUT / 'kamath_published_donor_metadata.csv', OUT / 'kamath_published_coverage_scenarios.csv'
        ]
    }
    report['verified_matrix_sha256'] = matrix_hash
    (OUT / 'kamath_independent_input_validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
