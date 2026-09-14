#!/usr/bin/env python3
"""Recompute the full legacy ORA family and audit orthogonal assay coverage."""
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from stats_core import bh

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=HERE/'results'
DATA=Path('/data32TB/shared/ppmi-foundin')


def main():
    legacy=pd.read_csv(ROOT/'03-deliverables/TableS_pathway_enrichment.csv')
    deg=pd.read_csv(ROOT/'03-deliverables/TableS_DEG_HC_vs_iPD_all_celltypes.csv')
    deg['gene']=deg.gene.str.upper()
    libraries={}
    for path in sorted((ROOT/'02-reanalysis/genesets').glob('*.gmt')):
        terms={}
        for line in path.read_text().splitlines():
            parts=line.split('\t')
            if len(parts)<3:continue
            terms[parts[0].strip()]={g.split(',')[0].strip().upper() for g in parts[2:] if g.strip()}
        libraries[path.stem]=terms
    rows=[]
    for ct,d in deg.groupby('CellType'):
        universe=set(d.gene)
        significant=set(d.loc[d['adj.P.Val']<.05,'gene'])
        query=significant if len(significant)>=10 else set(d.nsmallest(300,'P.Value').gene)
        for library,terms in libraries.items():
            for name,genes in terms.items():
                g=genes&universe
                if not 5<=len(g)<=500:continue
                k=len(query&g)
                p=stats.hypergeom.sf(k-1,len(universe),len(g),len(query))
                rows.append(dict(CellType=ct,library=library,term=name,universe_size=len(universe),
                                 query_size=len(query),set_size=len(g),overlap=k,p_value=p,
                                 overlap_genes=';'.join(sorted(query&g))))
    result=pd.DataFrame(rows);result['q_full_family']=np.nan
    for _,idx in result.groupby('CellType').groups.items():
        result.loc[idx,'q_full_family']=bh(result.loc[idx,'p_value'])
    joined=legacy.merge(result,on=['CellType','library','term'],suffixes=('_legacy','_recomputed'),validate='one_to_one')
    assert len(joined)==len(legacy)
    assert np.allclose(joined.p_value_legacy,joined.p_value_recomputed,rtol=1e-10,atol=1e-14)
    assert np.all(joined.overlap_legacy==joined.overlap_recomputed)
    result.to_csv(OUT/'legacy_ORA_full_testing_family.csv.gz',index=False)
    summary=dict(legacy_rows=len(legacy),full_family_rows=len(result),
                 original_significant=int((legacy.q_value<.05).sum()),
                 corrected_significant=int((result.q_full_family<.05).sum()),
                 original_pvalues_reproduced=True)
    (OUT/'legacy_ORA_validation.json').write_text(json.dumps(summary,indent=2)+'\n')

    clinical=pd.read_csv(HERE/'inputs/source_workbook_metadata.tsv',sep='\t').drop_duplicates('PPMI_ID')
    clinical=clinical.set_index('PPMI_ID')
    protein_rows=[];assay_rows=[]
    for replicate in [1,2,3]:
        path=DATA/f'ibm-aspera-gene-data/processed/PROT/PROT_Replicate{replicate}.csv'
        d=pd.read_csv(path)
        ratios=[c for c in d if c.startswith('Abundance Ratio:') and 'PPMI' in c]
        for column in ratios:
            donor=int(re.search(r'PPMI(\d+)',column).group(1))
            c=clinical.loc[donor]
            group={'Healthy Control':'HC','Idiopathic PD':'iPD'}.get(c.Disease_Status,'other')
            assay_rows.append(dict(assay='proteomics',replicate=replicate,PPMI_ID=donor,group=group,
                                   disease_status=c.Disease_Status,genetic_status=c.Genetic_Status))
            if 'Gene Symbol' in d:
                for _,row in d[d['Gene Symbol'].isin(['SNCA','PRKN','PARK2','PINK1','ATG5','ATG7','LAMP2','MAP1LC3B'])].iterrows():
                    protein_rows.append(dict(replicate=replicate,PPMI_ID=donor,group=group,gene=row['Gene Symbol'],
                                             abundance_ratio=row[column]))
    assays=pd.DataFrame(assay_rows)
    assays.to_csv(OUT/'proteomics_assay_manifest.csv',index=False)
    pd.DataFrame(protein_rows).to_csv(OUT/'proteomics_targeted_coverage.csv',index=False)
    atac=pd.read_csv(ROOT/'02-reanalysis/scatac_cell_meta.tsv.gz',sep='\t')
    am=atac.drop_duplicates('PPMI_ID')[['PPMI_ID','subtype']].copy()
    am['source_disease_status']=am.PPMI_ID.map(clinical.Disease_Status)
    am['source_genetic_status']=am.PPMI_ID.map(clinical.Genetic_Status)
    am.to_csv(OUT/'scatac_donor_manifest_audit.csv',index=False)
    audit=dict(protein_donors_by_group=assays.groupby('group').PPMI_ID.nunique().to_dict(),
               protein_assay_replicates=3,
               protein_note='Three assay replicates do not create three independent donor cohorts. Small HC/iPD groups cannot validate differential correlations.',
               scatac_cached_groups=am.groupby('subtype').PPMI_ID.nunique().to_dict(),
               scatac_note='The archived gene-level RNA assay is accessibility-derived activity; it is not measured RNA expression. Only two HC donors are present.',
               bulk_note='Day-65 bulk includes zero HC donors disjoint from scRNA. This is cross-assay triangulation, not independent replication.')
    (OUT/'orthogonal_assay_feasibility.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(summary,indent=2));print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
