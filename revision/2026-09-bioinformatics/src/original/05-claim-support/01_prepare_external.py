#!/usr/bin/env python3
"""Audit published annotations and stream public GSE243639 counts without testing disease effects."""
import csv
import gzip
import hashlib
import io
import json
import re
import time
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXT = HERE / 'external'
OUT = HERE / 'inputs'
OUT.mkdir(exist_ok=True)
NS = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}


def xlsx_sheet(z, sheet, shared):
    root = ET.fromstring(z.read('xl/worksheets/sheet%d.xml' % sheet))
    ns = {'m': root.tag.split('}')[0][1:]}
    records = []
    for row in root.findall('m:sheetData/m:row', ns):
        values = {}
        for cell in row.findall('m:c', ns):
            column = re.match('[A-Z]+', cell.attrib['r']).group()
            value = cell.find('m:v', ns)
            if value is None:
                continue
            value = value.text
            values[column] = shared[int(value)] if cell.attrib.get('t') == 's' else value
        records.append(values)
    header = records[0]
    return pd.DataFrame([{header[k]: v for k, v in row.items() if k in header}
                         for row in records[1:]])


def metadata():
    text = gzip.open(EXT / 'GSE243639_Clinical_data.csv.gz', 'rt').read()
    start = next(i for i, line in enumerate(text.splitlines()) if line.startswith('N;Brain Bank ID;'))
    clinical = pd.read_csv(io.StringIO('\n'.join(text.splitlines()[start:])), sep=';')
    clinical.columns = clinical.columns.str.strip()
    clinical = clinical[clinical['Sample ID'].notna()].copy()
    clinical['donor'] = clinical['Sample ID'].str.strip()
    clinical['group'] = clinical['Clinical diagnosis'].str.strip().map({"Parkinson's": 'iPD', 'Control': 'HC'})
    assert clinical.donor.is_unique and clinical.group.notna().all(), clinical.to_string()
    clinical.to_csv(OUT / 'gse243639_donors.csv', index=False)
    with ZipFile(EXT / 'GSE243639_UMAP_coordinates.xlsx') as z:
        shared = [''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
        broad = xlsx_sheet(z, 1, shared).rename(columns={'IDENT': 'broad_type', 'CLUSTER': 'broad_cluster'})
        neurons = xlsx_sheet(z, 4, shared).rename(columns={'IDENT': 'neuron_type', 'CLUSTER': 'neuron_cluster'})
    assert broad.CELL_ID.is_unique and neurons.CELL_ID.is_unique
    assert set(neurons.CELL_ID) <= set(broad.CELL_ID)
    broad = broad.merge(neurons[['CELL_ID', 'neuron_type', 'neuron_cluster']], how='left', on='CELL_ID', validate='one_to_one')
    broad['donor'] = broad.CELL_ID.str.split('_').str[0]
    broad = broad.merge(clinical, on='donor', how='left', validate='many_to_one')
    assert broad.group.notna().all()
    broad.to_csv(OUT / 'gse243639_cell_metadata.csv.gz', index=False)
    coverage = broad.groupby(['donor', 'group', 'broad_type'], observed=True).size().rename('n_cells').reset_index()
    coverage.to_csv(OUT / 'gse243639_broad_coverage.csv', index=False)
    nc = broad.dropna(subset=['neuron_type']).groupby(['donor', 'group', 'neuron_type'], observed=True).size().unstack(fill_value=0)
    nc.to_csv(OUT / 'gse243639_neuron_coverage.csv')
    print('Donors', clinical.groupby('group').size().to_dict(), flush=True)
    print('Broad annotations', broad.broad_type.value_counts().to_dict(), flush=True)
    print('Neuronal coverage by donor\n' + nc.to_string(), flush=True)
    return broad


def counts(meta):
    settings = json.loads((HERE.parent / '04-reviewer-reanalysis/results/analysis_settings.json').read_text())
    markers = ['TH', 'SLC6A3', 'SLC18A2', 'ALDH1A1', 'KCNJ6', 'NR4A2', 'TPH1', 'TPH2', 'SLC6A4',
               'SNAP25', 'GAD1', 'GAD2', 'MAP2', 'DCX', 'TUBB3', 'RBFOX3', 'SOX2', 'MKI67',
               'GFAP', 'AQP4', 'MBP', 'PLP1', 'P2RY12']
    requested = list(dict.fromkeys(['SNCA'] + settings['original_panel_genes'] + markers))
    all_genes = []
    selected = {}
    path = EXT / 'GSE243639_Filtered_count_table.csv.gz'
    t = time.time()
    with gzip.open(path, 'rt') as f:
        header = next(csv.reader([f.readline()]))
        cells = [re.sub(r'\.1$', '-1', s) for s in header[1:]]
        assert len(cells) == len(set(cells)) == len(meta)
        assert set(cells) == set(meta.CELL_ID)
        total = np.zeros(len(cells), dtype=np.float64)
        mitochondrial = np.zeros_like(total)
        detected = np.zeros(len(cells), dtype=np.int64)
        for i, line in enumerate(f, 1):
            label, values = line.split(',', 1)
            gene = next(csv.reader([label]))[0]
            a = np.fromstring(values, sep=',', dtype=np.float64)
            assert len(a) == len(cells) and np.isfinite(a).all() and (a >= 0).all(), gene
            assert np.equal(a, np.floor(a)).all(), 'Count table is not raw integer counts: ' + gene
            total += a
            detected += a > 0
            if gene.upper().startswith('MT-'):
                mitochondrial += a
            if gene in requested:
                assert gene not in selected, 'Duplicate gene symbol ' + gene
                selected[gene] = a.astype(np.float32)
            all_genes.append(gene)
            if i % 3000 == 0:
                print('Processed genes', i, 'elapsed seconds', round(time.time() - t, 1), flush=True)
    assert len(all_genes) == len(set(all_genes))
    assert (total > 0).all()
    kept = [g for g in requested if g in selected]
    matrix = pd.DataFrame({g: selected[g] for g in kept}, index=pd.Index(cells, name='cell'))
    matrix.to_csv(OUT / 'gse243639_panel_counts.tsv.gz', sep='\t', float_format='%.0f')
    qc = pd.DataFrame({'cell': cells, 'nCount_RNA': total, 'nFeature_RNA': detected,
                       'percent_mt': 100 * mitochondrial / total})
    qc.to_csv(OUT / 'gse243639_qc.csv.gz', index=False)
    (OUT / 'gse243639_all_genes.txt').write_text('\n'.join(all_genes) + '\n')
    audit = {'public_accession': 'GSE243639', 'n_cells': len(cells), 'n_genes': len(all_genes),
             'source_counts_are_nonnegative_integers': True, 'cell_ids_exactly_match_published_metadata': True,
             'cell_id_conversion': 'Final .1 in R-exported CSV converted to -1 in XLSX only',
             'selected_genes': kept, 'missing_requested_genes': [g for g in requested if g not in kept],
             'no_disease_effects_tested_during_preparation': True,
             'sources': {'counts': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243639/suppl/' + path.name,
                         'paper': 'https://doi.org/10.1186/s13024-023-00699-0'},
             'sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in EXT.glob('GSE243639*')}}
    (OUT / 'gse243639_input_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    print(json.dumps(audit, indent=2), flush=True)


if __name__ == '__main__':
    counts(metadata())
