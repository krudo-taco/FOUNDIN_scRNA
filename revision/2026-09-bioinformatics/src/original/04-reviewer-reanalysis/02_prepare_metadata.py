#!/usr/bin/env python3
"""Build assay manifests and bulk inputs without inspecting outcome significance."""
import json
import re
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = Path('/data32TB/shared/ppmi-foundin')
OUT = HERE / 'inputs'
OUT.mkdir(exist_ok=True, parents=True)


def workbook_table(path):
    ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
    with ZipFile(path) as z:
        strings = [''.join(t.text or '' for t in si.findall('.//m:t', ns))
                   for si in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', ns)]
        rows = []
        sheet = ET.fromstring(z.read('xl/worksheets/sheet1.xml'))
        for row in sheet.findall('m:sheetData/m:row', ns):
            vals = {}
            for cell in row.findall('m:c', ns):
                v = cell.find('m:v', ns)
                val = v.text if v is not None else ''
                if cell.attrib.get('t') == 's':
                    val = strings[int(val)]
                vals[re.sub(r'\d', '', cell.attrib['r'])] = val
            rows.append(vals)
    cols = list('ABCDEFGHIJKLMNOPQRS')
    # Repeated MAP2 heading is intentionally disambiguated by worksheet column.
    header = [f'{rows[1][c]}__{c}' if c in ['K', 'M'] else rows[1][c] for c in cols]
    return pd.DataFrame([[r.get(c, '') for c in cols] for r in rows[2:]], columns=header)


def main():
    clinical = workbook_table(DATA / 'mmc2-2.xlsx')
    clinical = clinical[clinical.PATNO.str.startswith('PPMI')].copy()
    clinical['PPMI_ID'] = clinical.PATNO.str.extract(r'PPMI(\d+)')[0].astype(int)
    clinical['BATCH'] = pd.to_numeric(clinical.BATCH).astype(int)
    clinical['culture_key'] = clinical.PATNO.str.replace('_', '', regex=False)
    clinical.to_csv(OUT / 'source_workbook_metadata.tsv', sep='\t', index=False)

    cell = pd.read_csv(ROOT / '02-reanalysis/cell_meta.tsv.gz', sep='\t')
    atac = pd.read_csv(ROOT / '02-reanalysis/scatac_cell_meta.tsv.gz', sep='\t')
    sex = pd.concat([cell[['PPMI_ID', 'genetic_sex']], atac[['PPMI_ID', 'genetic_sex']]]).drop_duplicates()
    assert sex.groupby('PPMI_ID').genetic_sex.nunique().max() == 1
    sex = sex.set_index('PPMI_ID').genetic_sex
    sample = cell.groupby(['SampleID', 'PPMI_ID', 'subtype', 'genetic_sex', 'BATCH']).size().reset_index(name='n_cells')
    sample.to_csv(OUT / 'scrna_sample_manifest.tsv', sep='\t', index=False)

    h5 = sorted((DATA / 'ibm-aspera-gene-data/processed/SCRN/cellRanger').glob('*/outs/filtered_feature_bc_matrix.h5'))[0]
    with h5py.File(h5, 'r') as f:
        ids = f['matrix/features/id'][:].astype(str)
        symbols = f['matrix/features/name'][:].astype(str)
    annotation = pd.DataFrame({'gene_id': ids, 'symbol': symbols})
    assert annotation.gene_id.is_unique
    annotation.to_csv(OUT / 'source_gene_annotation.tsv', sep='\t', index=False)

    bulk_path = DATA / 'ibm-aspera-gene-data/processed/RNAB/aggregated_expression/countTable.tsv'
    bulk = pd.read_csv(bulk_path, sep='\t', index_col=0)
    assert bulk.index.is_unique and bulk.columns.is_unique
    gene_map = annotation.set_index('gene_id').symbol
    mapped = gene_map.reindex(bulk.index)
    gene_labels = pd.Series(bulk.index, index=bulk.index).where(mapped.isna(), mapped)
    bulk = bulk.groupby(gene_labels, sort=False).sum()
    bulk.index.name = 'gene'
    bulk.to_csv(OUT / 'bulk_rna_counts_by_symbol.tsv.gz', sep='\t')

    records = []
    for name in bulk.columns:
        match = re.match(r'RNAB_(PPMI(\d+)[^_]*)_([^_]+)_da?(0|25|65)_v(\d+)$', name)
        if not match:
            raise ValueError(f'Unrecognized assay name {name}')
        patno, donor, line, day, version = match.groups()
        donor = int(donor)
        culture_key = 'PPMI3966' if patno == 'PPMI3966B1' else patno
        candidates = clinical[clinical.culture_key == culture_key]
        if len(candidates) != 1:
            raise ValueError(f'Ambiguous workbook match {name}: {len(candidates)}')
        c = candidates.iloc[0]
        group = {'Healthy Control': 'HC', 'Idiopathic PD': 'iPD'}.get(c.Disease_Status, 'other')
        # Only mutation-negative healthy controls are eligible for the primary comparison.
        if group == 'HC' and c.Genetic_Status != 'LRRK2-/SNCA-/GBA-':
            group = 'other'
        records.append(dict(sample=name, PATNO=patno, PPMI_ID=donor, day=int(day), version=int(version),
                            workbook_PATNO=c.PATNO,
                            BATCH=int(c.BATCH), group=group, disease_status=c.Disease_Status,
                            genetic_status=c.Genetic_Status, age_bin=c.ageBin,
                            genetic_sex=sex.get(donor, np.nan),
                            overlaps_scrna=donor in set(cell.PPMI_ID),
                            library_counts=int(bulk[name].sum())))
    manifest = pd.DataFrame(records)
    manifest.to_csv(OUT / 'bulk_rna_sample_manifest.tsv', sep='\t', index=False)
    overlap = manifest[manifest.group.isin(['HC', 'iPD'])].groupby(['day', 'group']).agg(
        assays=('sample', 'size'), donors=('PPMI_ID', 'nunique'))
    holdout = manifest[(manifest.day == 65) & ~manifest.overlaps_scrna & manifest.group.isin(['HC', 'iPD'])]
    summary = {'bulk_shape': list(bulk.shape), 'source_gene_mapping': str(h5),
               'bulk_comparison_counts': overlap.reset_index().to_dict('records'),
               'day65_donor_disjoint_counts': holdout.groupby('group').PPMI_ID.nunique().to_dict(),
               'missing_bulk_sex_donors': int(manifest[manifest.genetic_sex.isna()].PPMI_ID.nunique()),
               'source_workbook_rows': len(clinical), 'source_workbook_donors': clinical.PPMI_ID.nunique(),
               'reference_culture_rows': clinical[clinical.PPMI_ID == 3966][['PATNO', 'BATCH']].to_dict('records')}
    (OUT / 'metadata_audit.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
