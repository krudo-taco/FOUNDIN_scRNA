#!/usr/bin/env python3
"""Stream the full public matrix into a compact panel and QC summaries, without disease testing."""
import gzip
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXT, OUT = HERE / 'external', HERE / 'inputs'


def triplet_chunks(stream, block_bytes=8 * 1024 ** 2):
    carry = b''
    while True:
        block = stream.read(block_bytes)
        if not block:
            if carry.strip():
                values = np.fromstring(carry, sep=' ', dtype=np.int64)
                assert len(values) % 3 == 0
                yield values.reshape(-1, 3)
            return
        data = carry + block
        end = data.rfind(b'\n')
        if end < 0:
            carry = data
            continue
        body, carry = data[:end + 1], data[end + 1:]
        values = np.fromstring(body, sep=' ', dtype=np.int64)
        assert len(values) % 3 == 0
        if len(values):
            yield values.reshape(-1, 3)


def main():
    features = pd.read_csv(EXT / 'GSE178265_Homo_features.tsv.gz', sep='\t', header=None,
                           names=['feature_id', 'gene', 'type'], dtype=str)
    cells = pd.read_csv(EXT / 'GSE178265_Homo_bcd.tsv.gz', header=None, names=['cell'], dtype=str).cell
    assert cells.is_unique
    previous = json.loads((HERE.parent / '04-reviewer-reanalysis/results/analysis_settings.json').read_text())
    markers = ['TH', 'SLC6A3', 'SLC18A2', 'ALDH1A1', 'KCNJ6', 'NR4A2', 'TPH1', 'TPH2', 'SLC6A4',
               'SNAP25', 'GAD1', 'GAD2', 'MAP2', 'DCX', 'TUBB3', 'RBFOX3', 'SOX2', 'MKI67',
               'GFAP', 'AQP4', 'MBP', 'PLP1', 'P2RY12']
    requested = list(dict.fromkeys(['SNCA'] + previous['original_panel_genes'] + markers))
    canonical = features.gene.replace({'PARK2': 'PRKN'})
    genes = [g for g in requested if canonical.eq(g).any()]
    assert 'SNCA' in genes
    lookup = canonical.map({g: i for i, g in enumerate(genes)}).fillna(-1).to_numpy(np.int64)
    mito = features.gene.str.upper().str.startswith('MT-').to_numpy()
    n, g = len(cells), len(features)
    panel = np.zeros((len(genes), n), dtype=np.int64)
    library = np.zeros(n, dtype=np.float64)
    detected = np.zeros(n, dtype=np.int64)
    mt = np.zeros(n, dtype=np.float64)
    gene_total = np.zeros(g, dtype=np.float64)
    seen, last_report = 0, time.time()
    last_pair = (-1, -1)
    path = EXT / 'GSE178265_Homo_matrix.mtx.gz'
    with gzip.open(path, 'rb') as stream:
        assert stream.readline().strip().lower() == b'%%matrixmarket matrix coordinate integer general'
        line = stream.readline()
        while line.startswith(b'%'):
            line = stream.readline()
        ng, nc, expected = map(int, line.split())
        assert (ng, nc) == (g, n)
        for block in triplet_chunks(stream):
            row, col, value = block[:, 0] - 1, block[:, 1] - 1, block[:, 2]
            assert row.min() >= 0 and row.max() < g and col.min() >= 0 and col.max() < n
            assert value.min() > 0
            # A strictly increasing (column, row) order proves unique coordinates
            # and makes the nonzero count a valid detected-feature count.
            assert (int(col[0]), int(row[0])) > last_pair
            assert ((col[1:] > col[:-1]) | ((col[1:] == col[:-1]) & (row[1:] > row[:-1]))).all()
            last_pair = (int(col[-1]), int(row[-1]))
            lo, hi = int(col[0]), int(col[-1]) + 1
            library[lo:hi] += np.bincount(col - lo, weights=value, minlength=hi - lo)
            detected[lo:hi] += np.bincount(col - lo, minlength=hi - lo)
            is_mt = mito[row]
            mt[lo:hi] += np.bincount(col[is_mt] - lo, weights=value[is_mt], minlength=hi - lo)
            gene_total += np.bincount(row, weights=value, minlength=g)
            selected = lookup[row] >= 0
            np.add.at(panel, (lookup[row[selected]], col[selected]), value[selected])
            seen += len(block)
            if time.time() - last_report > 20:
                print('Processed entries', seen, 'of', expected, 'percent', round(100 * seen / expected, 1), flush=True)
                last_report = time.time()
    assert seen == expected and np.all(library > 0)
    assert np.sum(library) == np.sum(gene_total)
    np.savez_compressed(OUT / 'kamath_panel_counts.npz', counts=panel, genes=np.array(genes), cells=cells.to_numpy(str))
    qc = pd.DataFrame({'cell': cells, 'nCount_RNA': library.astype(np.int64),
                       'nFeature_RNA': detected, 'percent_mt': 100 * mt / library})
    qc.to_csv(OUT / 'kamath_cell_qc.csv.gz', index=False)
    features['canonical_panel_symbol'] = canonical
    features['total_UMIs'] = gene_total.astype(np.int64)
    features.to_csv(OUT / 'kamath_feature_audit.csv', index=False)
    audit = {'matrix_genes': g, 'matrix_cells': n, 'matrix_nonzero_entries': int(seen),
             'full_gzip_stream_and_CRC': 'PASS', 'strictly_sorted_unique_coordinates': True,
             'positive_integer_counts': True, 'library_total_equals_gene_total': True,
             'total_UMIs': int(library.sum()), 'retained_panel_genes': genes,
             'missing_requested_genes': [x for x in requested if x not in genes],
             'symbol_mapping': 'PARK2 to PRKN; records sharing a canonical panel symbol are summed',
             'no_donor_or_cell_type_assignments_inferred': True,
             'disease_effect_testing': 'Not assessed by count extraction; see separate annotation and statistical audits'}
    (OUT / 'kamath_count_extraction_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == '__main__':
    main()
