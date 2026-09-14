#!/usr/bin/env python3
"""Match every matrix library to GEO metadata and compare donor attributes with final paper tables."""
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
INPUT, EXT = HERE / 'inputs', HERE / 'external'
source = EXT / 'GSE178265_family.soft.gz'
text = gzip.open(source, 'rt').read()
records = []
for block in text.split('^SAMPLE = ')[1:]:
    row = {'GSM': block.splitlines()[0].strip()}
    for line in block.splitlines():
        if line.startswith('!Sample_title = '):
            row['library'] = line.split(' = ', 1)[1]
        elif line.startswith('!Sample_organism_ch1 = '):
            row['organism'] = line.split(' = ', 1)[1]
        elif line.startswith('!Sample_characteristics_ch1 = '):
            key, value = line.split(' = ', 1)[1].split(':', 1)
            assert key.strip().lower() not in row
            row[key.strip().lower()] = value.strip()
    records.append(row)
geo = pd.DataFrame(records)
assert geo.GSM.is_unique and geo.library.is_unique
coverage = pd.read_csv(INPUT / 'kamath_GEO_library_coverage.csv')
joined = coverage.merge(geo, on='library', how='left', validate='one_to_one')
assert len(joined) == 97 and joined.GSM.notna().all() and joined.organism.eq('Homo sapiens').all()
clinical = pd.read_csv(INPUT / 'kamath_published_donor_metadata.csv', dtype={'donor': str}).set_index('donor')
conflicts = []
donors = []
for i, row in joined.iterrows():
    matches = [d for d in clinical.index if d in row.library]
    assert len(matches) <= 1, row.library
    donor = matches[0] if matches else None
    donors.append(donor)
    if donor is None or row.tissue != 'Substantia nigra pars compacta':
        continue
    published = clinical.loc[donor]
    comparisons = {
        'diagnosis': (published.clinical_group, {'Ctrl': 'HC', 'PD': 'PD', 'LBD': 'LBD'}.get(row.status, row.status)),
        'sex': (published['Sex'], {'Female': 'F', 'Male': 'M'}.get(row.sex, row.sex)),
        'age': (float(published['Age (yrs)']), float(row.age)),
        'PMI': (float(published['PMI (hrs)']), float(row.pmi)),
    }
    for key, (paper, deposited) in comparisons.items():
        same = paper == deposited if isinstance(paper, str) else np.isclose(paper, deposited)
        if not same:
            conflicts.append(dict(donor=donor, library=row.library, GSM=row.GSM,
                                  field=key, published_table=paper, GEO_characteristic=deposited))
joined['donor_from_unique_library_name_match'] = donors
joined['SNpc_library'] = joined.tissue.eq('Substantia nigra pars compacta')
joined.to_csv(INPUT / 'kamath_GEO_library_manifest.csv', index=False)
conflicts = pd.DataFrame(conflicts, columns=['donor', 'library', 'GSM', 'field', 'published_table', 'GEO_characteristic'])
conflicts.to_csv(INPUT / 'kamath_clinical_source_conflicts.csv', index=False)
summary = joined.groupby(['tissue', 'status'], dropna=False).agg(libraries=('GSM', 'nunique'), matrix_cells=('n_barcodes', 'sum'))
summary.to_csv(INPUT / 'kamath_GEO_tissue_summary.csv')
audit = {
    'all_matrix_libraries_match_GEO_sample_titles': True,
    'GEO_samples': len(geo), 'human_matrix_libraries': len(joined),
    'matched_matrix_barcodes': int(joined.n_barcodes.sum()),
    'unmapped_SNpc_libraries': joined[joined.SNpc_library & joined.donor_from_unique_library_name_match.isna()].library.tolist(),
    'clinical_conflicts': conflicts[['donor', 'field', 'published_table', 'GEO_characteristic']].drop_duplicates().to_dict('records'),
    'donor_mapping_basis': 'Unique published donor ID in exact GEO library title, checked against tissue and clinical fields',
    'clinical_primary_source': 'Final publication supplementary clinical tables; discrepancies preserved for sensitivity planning',
    'cell_type_annotation_status': 'Not established by library identity or FACS enrichment',
    'source': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178265/soft/GSE178265_family.soft.gz',
    'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
}
(INPUT / 'kamath_GEO_library_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
print(summary.to_string())
print(json.dumps(audit, indent=2))
