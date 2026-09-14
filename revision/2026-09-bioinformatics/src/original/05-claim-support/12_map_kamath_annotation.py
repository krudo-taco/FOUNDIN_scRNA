#!/usr/bin/env python3
"""Map author library aliases using barcode identity and audit final cell membership."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
INPUT = HERE / 'inputs'


def main():
    author = pd.read_csv(INPUT / 'kamath_author_metadata.csv.gz')
    assert author.source_cell_id.is_unique
    parts = author.source_cell_id.str.rsplit('_', n=1, expand=True)
    assert parts[0].eq(author['orig.ident']).all()
    author['barcode'] = parts[1]
    author['donor'] = author['orig.ident'].str.extract(r'(\d+)', expand=False)
    assert author.groupby('source_list_index').donor.nunique().eq(1).all()
    clinical = pd.read_csv(INPUT / 'kamath_published_donor_metadata.csv', dtype={'donor': str})
    assert set(author.donor) == set(clinical.donor)
    libraries = pd.read_csv(INPUT / 'kamath_GEO_library_manifest.csv', dtype={'donor_from_unique_library_name_match': str})
    qc = pd.read_csv(INPUT / 'kamath_cell_qc.csv.gz')
    parts = qc.cell.str.rsplit('_', n=1, expand=True)
    qc['library'], qc['barcode'] = parts[0], parts[1]
    qc = qc.merge(libraries[['library', 'donor_from_unique_library_name_match', 'SNpc_library']], on='library', validate='many_to_one')
    known = qc[qc.SNpc_library & qc.donor_from_unique_library_name_match.notna()].copy()
    candidates, chosen = [], []
    for library, cells in known.groupby('library'):
        donor = cells.donor_from_unique_library_name_match.iloc[0]
        target = set(cells.barcode)
        scores = []
        for alias, sub in author[author.donor.eq(donor)].groupby('orig.ident'):
            overlap = len(target.intersection(sub.barcode))
            record = dict(donor=donor, library=library, author_alias=alias,
                          matrix_cells=len(target), author_cells=len(sub), barcode_overlap=overlap)
            candidates.append(record)
            scores.append(record)
        scores = sorted(scores, key=lambda x: x['barcode_overlap'], reverse=True)
        assert scores[0]['barcode_overlap'] > scores[1]['barcode_overlap']
        assert scores[0]['barcode_overlap'] / len(target) >= .95
        chosen.append(dict(scores[0], second_best_overlap=scores[1]['barcode_overlap']))
    mapping = pd.DataFrame(chosen)
    assert mapping.library.is_unique and mapping.author_alias.is_unique and len(mapping) == 89
    pd.DataFrame(candidates).to_csv(INPUT / 'kamath_library_alias_candidates.csv', index=False)
    mapping.to_csv(INPUT / 'kamath_library_alias_map.csv', index=False)
    author = author.merge(mapping[['author_alias', 'library']], left_on='orig.ident', right_on='author_alias', validate='many_to_one')
    author['cell'] = author.library + '_' + author.barcode
    assert author.cell.is_unique
    selected = ['cell', 'source_cell_id', 'source_list_index', 'donor', 'orig.ident', 'clusters', 'nGene', 'nUMI']
    joined = known.merge(author[selected], on='cell', how='left', validate='one_to_one', indicator=True)
    joined.to_csv(INPUT / 'kamath_cell_annotation_mapping.csv.gz', index=False)
    missing = joined[joined['_merge'].ne('both')]
    missing.to_csv(INPUT / 'kamath_annotation_unmatched_cells.csv', index=False)
    matched = joined[joined['_merge'].eq('both')].copy()
    assert matched.donor.eq(matched.donor_from_unique_library_name_match).all()
    differences = matched[matched.nGene.ne(matched.nFeature_RNA) | matched.nUMI.ne(matched.nCount_RNA)]
    differences.to_csv(INPUT / 'kamath_annotation_QC_discrepancies.csv', index=False)
    class_counts = matched.groupby(['donor', 'clusters']).size().unstack(fill_value=0)
    class_counts.to_csv(INPUT / 'kamath_mapped_author_celltype_counts.csv')
    da = class_counts['DA'].rename('mapped_author_DA_nuclei').reset_index()
    da = clinical[['donor', 'clinical_group', 'published_DA_nuclei']].merge(da, on='donor', validate='one_to_one')
    da['difference'] = da.mapped_author_DA_nuclei - da.published_DA_nuclei
    da.to_csv(INPUT / 'kamath_DA_publication_reconciliation.csv', index=False)
    retained = matched[matched.clusters.ne('REMOVE')]
    out = {
        'author_rows': len(author), 'author_library_aliases': len(mapping),
        'known_clinical_SNpc_matrix_cells': len(known), 'mapped_cells': len(matched),
        'unmatched_cells': len(missing), 'author_and_matrix_QC_discrepancies': len(differences),
        'mapped_celltype_counts': matched.clusters.value_counts().to_dict(),
        'retained_non_REMOVE_cells': len(retained), 'published_SNpc_cells': 387483,
        'DA_counts_match_each_published_donor': bool(da.difference.eq(0).all()),
        'mapped_DA_nuclei': int(da.mapped_author_DA_nuclei.sum()),
        'published_DA_nuclei': int(da.published_DA_nuclei.sum()),
        'mapping_basis': 'Unique maximum barcode overlap within donor; one-to-one library alias assignment. Cell IDs and source QC values checked separately.',
        'disease_effects_tested': False,
        'source_qs_sha256': hashlib.sha256((HERE / 'external/individualobject_seuratmetadataframe.qs').read_bytes()).hexdigest()
    }
    (INPUT / 'kamath_annotation_mapping_audit.json').write_text(json.dumps(out, indent=2) + '\n')
    print(da.to_string(index=False))
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
