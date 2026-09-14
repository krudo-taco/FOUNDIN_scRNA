#!/usr/bin/env python3
"""Between-donor correlation uncertainty, design sensitivity and bulk triangulation."""
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from stats_core import anchor_correlations, bh, design, fisher, residualize, rank_axis1

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
OUT=HERE/'results';OUT.mkdir(exist_ok=True)
INPUT=HERE/'inputs'
SEED=20260913
N_BOOT=10000
EXCLUDED=[3954,4106]


def batch_correlations(a):
    a=a-a.mean(axis=1,keepdims=True)
    var=np.sum(a*a,axis=1)
    numerator=np.sum(a[:,:,0,None]*a,axis=1)
    denom=np.sqrt(var[:,0,None]*var)
    return np.divide(numerator,denom,out=np.full_like(numerator,np.nan),where=denom>1e-12)


def correlation_analysis(values,md,label,genes,scenario='primary'):
    values=values.reindex(index=md.index,columns=['SNCA']+genes)
    keep=values.notna().all(axis=1)
    values=values.loc[keep];md=md.loc[keep]
    x=values.values
    masks=[md.group.values=='HC',md.group.values=='iPD']
    ns=[int(m.sum()) for m in masks]
    if min(ns)<5:
        return pd.DataFrame(),dict(dataset=label,scenario=scenario,n_HC=ns[0],n_iPD=ns[1],status='insufficient donors')
    nuisance,names=design(md)
    nuisance=nuisance[:,:-1]
    versions={'Pearson':x,'Spearman':np.column_stack([stats.rankdata(x[:,j]) for j in range(x.shape[1])]),
              'pooled_nuisance_residualized':residualize(x,nuisance)}
    rows=[]
    for method,a in versions.items():
        r=[anchor_correlations(a[m]) for m in masks]
        delta=r[1]-r[0]
        # Spearman correlations are ranked separately within diagnostic groups.
        if method=='Spearman':
            r=[anchor_correlations(np.column_stack([stats.rankdata(x[m,j]) for j in range(x.shape[1])])) for m in masks]
            delta=r[1]-r[0]
        rng=np.random.default_rng(SEED)
        boot=[]
        for m in masks:
            arm=a[m]
            picks=rng.integers(0,len(arm),size=(N_BOOT,len(arm)))
            if method=='Spearman':
                # Ranking within bootstrap samples handles the additional bootstrap ties.
                raw=x[m][picks]
                br=np.empty((N_BOOT,x.shape[1]))
                for start in range(0,N_BOOT,200):
                    chunk=raw[start:start+200]
                    ranked=rank_axis1(chunk)
                    br[start:start+len(chunk)]=batch_correlations(ranked)
            else:
                br=batch_correlations(arm[picks])
            boot.append(br)
        dboot=boot[1]-boot[0]
        lo,hi=np.nanpercentile(dboot,[2.5,97.5],axis=0)
        for j,g in enumerate(genes,1):
            z=(fisher(r[1][j])-fisher(r[0][j]))/np.sqrt(1/(ns[0]-3)+1/(ns[1]-3))
            p=float(2*stats.norm.sf(abs(z))) if method=='Pearson' else np.nan
            rows.append(dict(dataset=label,scenario=scenario,method=method,gene=g,n_HC=ns[0],n_iPD=ns[1],
                             r_HC=r[0][j],r_iPD=r[1][j],delta_r=delta[j],ci_delta_low=lo[j],ci_delta_high=hi[j],
                             p_fisher=p,bootstrap_valid=int(np.isfinite(dboot[:,j]).sum()),
                             uncertainty_note='Percentile donor bootstrap. Pooled-residual sensitivity conditions on the fitted nuisance adjustment; not a fully refitted adjusted CI.' if method=='pooled_nuisance_residualized' else 'Percentile donor bootstrap; small-control-cohort limitations apply.'))
    result=pd.DataFrame(rows);result['q_fisher_panel']=np.nan
    idx=result.method=='Pearson';result.loc[idx,'q_fisher_panel']=bh(result.loc[idx,'p_fisher'])
    return result,dict(dataset=label,scenario=scenario,n_HC=ns[0],n_iPD=ns[1],status='tested',nuisance_columns=names[:-1])


def metadata_from_cultures(cultures,weights=None):
    info=cultures.groupby('PPMI_ID').agg(group=('group','first'),genetic_sex=('genetic_sex','first'))
    b=pd.crosstab(cultures.PPMI_ID,cultures.BATCH).reindex(columns=range(1,6),fill_value=0)
    b=b.div(b.sum(axis=1),axis=0)
    for i in range(1,6):info[f'batch_fraction_{i}']=b[i]
    return info


def main():
    settings=json.loads((OUT/'analysis_settings.json').read_text())
    genes=settings['eligible_genes']
    results=[];audit=[]
    meta=pd.read_csv(ROOT/'02-reanalysis/cell_meta.tsv.gz',sep='\t',index_col=0)
    meta=meta[meta.subtype.isin(['na_HC','na_iPD'])]
    meta['group']=meta.subtype.map({'na_HC':'HC','na_iPD':'iPD'})
    expr=pd.read_csv(INPUT/'rna_panel_logexpr.tsv.gz',sep='\t',index_col=0).loc[meta.index,['SNCA']+genes]
    for context in ['iDA_pooled']+sorted(meta.CellType.unique()):
        mask=meta.CellType.str.startswith('iDA') if context=='iDA_pooled' else meta.CellType==context
        cm=meta.loc[mask]
        sizes=cm.SampleID.value_counts()
        cm=cm[cm.SampleID.isin(sizes[sizes>=20].index)]
        samples=cm.drop_duplicates('SampleID').set_index('SampleID')
        means=expr.loc[cm.index].groupby(cm.SampleID).mean()
        values=means.groupby(samples.PPMI_ID).mean()
        md=metadata_from_cultures(samples)
        keep=~md.index.isin(EXCLUDED)
        r,a=correlation_analysis(values.loc[keep],md.loc[keep],f'scRNA_{context}',genes)
        results.append(r);audit.append(a)
        if context=='iDA_pooled':
            for scenario,k in [('original_30_iPD_grouping',np.ones(len(md),bool)),
                               ('exclude_reference_donor',~md.index.isin(EXCLUDED+[3966]))]:
                r,a=correlation_analysis(values.loc[k],md.loc[k],f'scRNA_{context}',genes,scenario)
                results.append(r);audit.append(a)
        print(time.strftime('%H:%M:%S'),context,'between-donor complete',flush=True)

    bulk=pd.read_csv(INPUT/'bulk_rna_counts_by_symbol.tsv.gz',sep='\t',index_col=0)
    bm=pd.read_csv(INPUT/'bulk_rna_sample_manifest.tsv',sep='\t')
    # Version denotes a rerun in the release README. Only the latest version is
    # used for each culture/day; cultures still collapse to one donor value.
    bm['selected_latest_version']=bm.version==bm.groupby(['PATNO','day']).version.transform('max')
    bm.to_csv(OUT/'bulk_assay_selection.csv',index=False)
    bm=bm[bm.selected_latest_version & bm.group.isin(['HC','iPD'])].copy()
    for day in [65,0,25]:
        dm=bm[bm.day==day].set_index('sample').copy()
        counts=bulk.loc[:,dm.index]
        positive=(counts>0).all(axis=1)
        geom=np.exp(np.log(counts.loc[positive]).mean(axis=1))
        sf=counts.loc[positive].div(geom,axis=0).median(axis=0)
        sf=sf/np.exp(np.log(sf).mean())
        effective_lib=sf*np.exp(np.log(counts.sum(axis=0)).mean())
        lcpm=np.log2(counts.reindex(['SNCA']+genes).div(effective_lib,axis=1)*1e6+.5)
        values=lcpm.T.groupby(dm.PPMI_ID).mean()
        # Unknown sex is not imputed as male. The complete-assay analysis omits
        # sex, followed by a separate analysis restricted to reported-sex donors.
        md=metadata_from_cultures(dm)
        md['sex_reported']=md.genetic_sex.notna()
        original_sex=md.genetic_sex.copy();md['genetic_sex']=np.nan
        r,a=correlation_analysis(values,md,f'bulk_day{day}',genes,'all_donors_batch_adjustment')
        results.append(r);audit.append(a)
        md['genetic_sex']=original_sex
        keep=md.sex_reported
        if not keep.all():
            r,a=correlation_analysis(values.loc[keep],md.loc[keep],f'bulk_day{day}',genes,'reported_sex_subset')
            results.append(r);audit.append(a)
        if day==65:
            values.to_csv(OUT/'bulk_day65_donor_logcpm.csv')
            md.to_csv(OUT/'bulk_day65_donor_metadata.csv')
    res=pd.concat(results,ignore_index=True)
    res['q_fisher_global_primary_scrna']=np.nan
    idx=res.dataset.str.startswith('scRNA_') & (res.scenario=='primary') & (res.method=='Pearson')
    res.loc[idx,'q_fisher_global_primary_scrna']=bh(res.loc[idx,'p_fisher'])
    res.to_csv(OUT/'between_donor_and_bulk_correlations.csv',index=False)
    (OUT/'between_donor_analysis_audit.json').write_text(json.dumps(audit,indent=2)+'\n')

    power=[]
    for nh,npd in [(8,28),(8,30),(7,18),(8,32)]:
        for rh in [0.,.3,.6]:
            for rp in [-.6,-.3,0.,.3,.6,.8]:
                for alpha in [.05,.05/len(genes)]:
                    shift=(fisher(rp)-fisher(rh))/np.sqrt(1/(nh-3)+1/(npd-3))
                    critical=stats.norm.ppf(1-alpha/2)
                    pwr=stats.norm.sf(critical-shift)+stats.norm.cdf(-critical-shift)
                    power.append(dict(n_HC=nh,n_iPD=npd,r_HC=rh,r_iPD=rp,delta_r=rp-rh,
                                      alpha=alpha,approximate_power=float(pwr)))
    pd.DataFrame(power).to_csv(OUT/'correlation_difference_design_power.csv',index=False)
    print('Between-donor and bulk analyses complete',flush=True)


if __name__=='__main__':main()
