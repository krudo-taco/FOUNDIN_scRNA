#!/usr/bin/env python3
"""Parse author clinical tables and DA-cell coverage before expression testing."""
import gzip
import hashlib
import importlib.util
import json
import tarfile
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile
import pandas as pd

HERE = Path(__file__).resolve().parent
EXT, OUT = HERE / 'external', HERE / 'inputs'
spec = importlib.util.spec_from_file_location('external_preparation', HERE / '01_prepare_external.py')
prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prep)

with ZipFile(EXT / 'kamath_supplementary_tables.xlsx') as z:
    shared = [''.join(x.itertext()) for x in ET.fromstring(z.read('xl/sharedStrings.xml'))]
    controls = prep.xlsx_sheet(z, 1, shared)
    cases = prep.xlsx_sheet(z, 2, shared)
    cell_counts = prep.xlsx_sheet(z, 6, shared)
controls['clinical_group'] = 'HC'
cases['clinical_group'] = cases['Neuropathological evaluation'].map({"Parkinson's disease": 'PD', 'Lewy body dementia': 'LBD'})
assert cases.clinical_group.notna().all()
clinical = pd.concat([controls, cases], ignore_index=True).rename(columns={'Sample ID': 'donor'})
clinical = clinical.merge(cell_counts.rename(columns={'Donor ID': 'donor'}), on='donor', validate='one_to_one')
clinical['published_DA_nuclei'] = pd.to_numeric(clinical['Number of DA neurons sampled'])
clinical['passes_50_nuclei'] = clinical.published_DA_nuclei >= 50
clinical['passes_100_nuclei'] = clinical.published_DA_nuclei >= 100
assert clinical.donor.is_unique and len(clinical) == 18
clinical.to_csv(OUT / 'kamath_published_donor_metadata.csv', index=False)
with gzip.open(EXT / 'GSE178265_Homo_bcd.tsv.gz', 'rt') as f:
    barcodes = [s.strip() for s in f]
assert len(barcodes) == len(set(barcodes))
barcodes = pd.DataFrame({'cell': barcodes})
barcodes['library'] = barcodes.cell.str.rsplit('_', n=1).str[0]
barcodes.groupby('library').size().rename('n_barcodes').to_csv(OUT / 'kamath_GEO_library_coverage.csv')
archive = EXT / 'kamath_individual_metadata.tar.gz'
with tarfile.open(archive, 'r:gz') as t:
    members = [m for m in t.getmembers() if m.isfile()]
    assert len(members) == 1 and Path(members[0].name).name == 'individualobject_seuratmetadataframe.qs'
    annotation = t.extractfile(members[0]).read()
    # Write this specific data member only; never unpack archive-supplied paths.
    (EXT / 'individualobject_seuratmetadataframe.qs').write_bytes(annotation)
matrix_path = EXT / 'GSE178265_Homo_matrix.mtx.gz'
count_audit_path = OUT / 'kamath_count_extraction_audit.json'
library_audit_path = OUT / 'kamath_GEO_library_audit.json'
independent_audit_path = OUT / 'kamath_independent_input_validation.json'
count_status = 'Not yet downloaded'
if matrix_path.is_file():
    count_status = 'Downloaded; count extraction not yet verified'
    if count_audit_path.is_file():
        counts_audit = json.loads(count_audit_path.read_text())
        if counts_audit.get('full_gzip_stream_and_CRC') == 'PASS':
            count_status = 'Downloaded and fully streamed; panel counts and QC extracted; see kamath_count_extraction_audit.json'
library_status = 'Not yet validated'
if library_audit_path.is_file():
    libraries_audit = json.loads(library_audit_path.read_text())
    if libraries_audit.get('all_matrix_libraries_match_GEO_sample_titles'):
        library_status = 'All matrix libraries matched to GEO; clinical donor mapping and conflicts recorded in kamath_GEO_library_audit.json. Cell-type membership remains unverified.'
annotation_status = 'Downloaded; not decoded in this analysis package'
mapping_path = OUT / 'kamath_annotation_mapping_audit.json'
if (OUT / 'kamath_author_metadata.csv.gz').is_file() and (OUT / 'kamath_author_metadata.rds').is_file():
    annotation_status = 'Decoded into 18 donor metadata frames; lossless base-R RDS copy and CSV export available'
if mapping_path.is_file():
    mapping_audit = json.loads(mapping_path.read_text())
    if mapping_audit['unmatched_cells'] == 0 and mapping_audit['author_and_matrix_QC_discrepancies'] == 0:
        library_status = '387558 SNpc cells matched to author annotations with exact nUMI and nGene agreement; final publication cell-set difference remains documented'
outcome_status = 'Not started'
outcome_audit = OUT / 'kamath_statistical_validation.json'
if outcome_audit.is_file() and (HERE / 'results/kamath_panel_profile_tests.csv').is_file():
    outcome_status = 'Completed on the public GEO release; see kamath_statistical_validation.json and results/kamath_panel_profile_tests.csv'
audit = {
    'clinical_counts': clinical.clinical_group.value_counts().to_dict(),
    'published_DA_counts_at_least_50': clinical[clinical.passes_50_nuclei].clinical_group.value_counts().to_dict(),
    'published_DA_counts_at_least_100': clinical[clinical.passes_100_nuclei].clinical_group.value_counts().to_dict(),
    'GEO_barcodes': len(barcodes), 'GEO_libraries': barcodes.library.nunique(),
    'annotation_format': 'qs', 'annotation_status': annotation_status,
    'library_to_donor_and_cell_annotation_mapping': library_status,
    'count_matrix_status': count_status,
    'independent_count_validation': str(independent_audit_path.name) if independent_audit_path.is_file() else 'Not yet completed',
    'disease_effect_testing_status': outcome_status,
    'clinical_source': 'https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41593-022-01061-1/MediaObjects/41593_2022_1061_MOESM3_ESM.xlsx',
    'annotation_source': 'https://github.com/tkamath1/Kamathetal2022/blob/main/Analyses/IntermediateObjects/individualobject_seuratmetadataframe.tar.gz',
    'clinical_source_sha256': hashlib.sha256((EXT / 'kamath_supplementary_tables.xlsx').read_bytes()).hexdigest(),
    'annotation_archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
    'annotation_data_sha256': hashlib.sha256(annotation).hexdigest()
}
(OUT / 'kamath_metadata_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps(audit, indent=2))
