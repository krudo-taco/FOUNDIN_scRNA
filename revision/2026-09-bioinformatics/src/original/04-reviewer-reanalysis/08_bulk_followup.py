#!/usr/bin/env python3
"""Follow up bulk association with matched cultures and fully refitted uncertainty."""
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from stats_core import anchor_correlations,bh,design,independent_columns,residualize,wild_test

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent;OUT=HERE/'results';INPUT=HERE/'inputs'


def culture_key(sample,batch):
    key=re.search(r'PPMI\d+[^_]*',sample).group(0)
    if key=='PPMI3966':key='PPMI3966B1'
    elif key in ['PPMI3966E6','PPMI3966E8']:key=key.replace('3966E','3966B5E')
    return key


def get_values():
    genes=json.loads((OUT/'analysis_settings.json').read_text())['eligible_genes']
    meta=pd.read_csv(ROOT/'02-reanalysis/cell_meta.tsv.gz',sep='\t')
    meta['is_iDA']=meta.CellType.str.startswith('iDA')
    fractions=meta.groupby('SampleID').agg(PPMI_ID=('PPMI_ID','first'),BATCH=('BATCH','first'),
        genetic_sex=('genetic_sex','first'),iDA_fraction=('is_iDA','mean'))
    fractions['culture_key']=[culture_key(s,b) for s,b in zip(fractions.index,fractions.BATCH)]
    assert fractions.culture_key.is_unique
    bulk=pd.read_csv(INPUT/'bulk_rna_counts_by_symbol.tsv.gz',sep='\t',index_col=0)
    manifest=pd.read_csv(OUT/'bulk_assay_selection.csv')
    m=manifest[(manifest.day==65)&manifest.selected_latest_version&manifest.group.isin(['HC','iPD'])]
    matched=m.merge(fractions[['culture_key','iDA_fraction']],left_on='PATNO',right_on='culture_key',validate='one_to_one')
    matched.to_csv(OUT/'bulk_scrna_matched_cultures.csv',index=False)
    assert not matched.genetic_sex.isna().any()
    counts=bulk.loc[:,matched['sample']]
    positive=(counts>0).all(axis=1)
    geo=np.exp(np.log(counts.loc[positive]).mean(axis=1))
    sf=counts.loc[positive].div(geo,axis=0).median(axis=0);sf/=np.exp(np.log(sf).mean())
    lib=sf*np.exp(np.log(counts.sum(axis=0)).mean())
    expr=np.log2(counts.loc[['SNCA']+genes].div(lib,axis=1)*1e6+.5).T
    mm=matched.set_index('sample').loc[expr.index]
    values=expr.groupby(mm.PPMI_ID).mean()
    md=mm.groupby('PPMI_ID').agg(group=('group','first'),genetic_sex=('genetic_sex','first'),
                                  iDA_fraction=('iDA_fraction','mean'))
    batches=pd.crosstab(mm.PPMI_ID,mm.BATCH).reindex(columns=range(1,6),fill_value=0)
    batches=batches.div(batches.sum(axis=1),axis=0)
    for b in range(1,6):md[f'batch_fraction_{b}']=batches[b]
    return values,md,genes


def refitted_correlations(values,md,genes,scenario,composition):
    a=values.loc[md.index,['SNCA']+genes].values
    x,_=design(md)
    # Include diagnosis when estimating nuisance effects so group mean shifts
    # are not inadvertently absorbed into correlated batch covariates.
    if composition:x=np.column_stack([x,md.iDA_fraction.values])
    x,_=independent_columns(x,list(range(x.shape[1])))
    residual=residualize(a,x)
    hc=md.group.values=='HC';pd_mask=~hc
    r0=anchor_correlations(residual[hc]);r1=anchor_correlations(residual[pd_mask])
    rng=np.random.default_rng(20260913)
    boot=np.empty((5000,a.shape[1]))
    hc_idx=np.where(hc)[0];pd_idx=np.where(pd_mask)[0]
    for i in range(len(boot)):
        idx=np.r_[rng.choice(hc_idx,len(hc_idx),replace=True),rng.choice(pd_idx,len(pd_idx),replace=True)]
        residual=residualize(a[idx],x[idx])
        boot[i]=anchor_correlations(residual[len(hc_idx):])-anchor_correlations(residual[:len(hc_idx)])
    lo,hi=np.nanpercentile(boot,[2.5,97.5],axis=0)
    return [dict(scenario=scenario,adjust_cell_fraction=composition,gene=g,n_HC=int(hc.sum()),n_iPD=int(pd_mask.sum()),
                 r_HC=r0[j],r_iPD=r1[j],delta_r=r1[j]-r0[j],ci_low=lo[j],ci_high=hi[j],bootstrap_draws=5000,
                 valid_draws=int(np.isfinite(boot[:,j]).sum())) for j,g in enumerate(genes,1)]


def main():
    values,md,genes=get_values()
    all_interactions=[];all_correlations=[];omnibus=[];influence=[]
    # This model tests group-dependent slopes in pooled-SD units. Its coefficient
    # is not reported as a Pearson correlation difference.
    scenarios=[('matched_cultures',np.ones(len(md),bool)),
               ('matched_without_reference',md.index!=3966),
               ('shared_batches_without_reference',(md.index!=3966)&(md.batch_fraction_4==0)&(md.batch_fraction_5==0))]
    for scenario,mask in scenarios:
        m=md.loc[mask];v=values.loc[mask]
        for composition in [False,True]:
            a=v.loc[:,['SNCA']+genes].values
            z=(a-a.mean(axis=0))/a.std(axis=0,ddof=1)
            x,names=design(m)
            cov=[x,z[:,0,None]]
            if composition:cov.append(m.iDA_fraction.values[:,None])
            nuisance=np.column_stack(cov)
            nuisance,_=independent_columns(nuisance,list(range(nuisance.shape[1])))
            full=np.column_stack([nuisance,z[:,0]*(m.group.values=='iPD')])
            weights=np.random.default_rng(20260913).choice([-1.,1.],size=(19999,len(m)))
            try:
                fit=wild_test(z[:,1:],full,weights)
                for j,g in enumerate(genes):
                    all_interactions.append(dict(scenario=scenario,adjust_cell_fraction=composition,gene=g,
                        n_HC=int((m.group=='HC').sum()),n_iPD=int((m.group=='iPD').sum()),
                        interaction_slope=fit['effect'][j],ci_low=fit['ci_low'][j],ci_high=fit['ci_high'][j],
                        p_wild=fit['p_wild'][j],status='tested',max_leverage=max(fit['leverage'])))
                obs=np.mean(fit['t']**2);null=np.mean(fit['tstar']**2,axis=1)
                omnibus.append(dict(scenario=scenario,adjust_cell_fraction=composition,
                                    p_omnibus=(1+(null>=obs).sum())/(len(null)+1)))
            except ValueError as e:
                for g in genes:all_interactions.append(dict(scenario=scenario,adjust_cell_fraction=composition,
                    gene=g,status=str(e),p_wild=np.nan))
            all_correlations.extend(refitted_correlations(v,m,genes,scenario,composition))
            print(scenario,'composition',composition,'complete',flush=True)
        for donor in m.index:
            k=m.index!=donor
            a=v.loc[k,['SNCA']+genes].values;mm=m.loc[k]
            x,_=design(mm);res=residualize(a,x)
            r0=anchor_correlations(res[mm.group.values=='HC']);r1=anchor_correlations(res[mm.group.values=='iPD'])
            for j,g in enumerate(genes,1):
                influence.append(dict(scenario=scenario,excluded_donor=int(donor),gene=g,delta_r=r1[j]-r0[j]))
    ir=pd.DataFrame(all_interactions);ir['q_wild_scenario_family']=np.nan
    for _,idx in ir.groupby('scenario').groups.items():ir.loc[idx,'q_wild_scenario_family']=bh(ir.loc[idx,'p_wild'])
    ir.to_csv(OUT/'bulk_conditional_interactions.csv',index=False)
    pd.DataFrame(all_correlations).to_csv(OUT/'bulk_refitted_adjusted_correlations.csv',index=False)
    pd.DataFrame(omnibus).to_csv(OUT/'bulk_interaction_omnibus.csv',index=False)
    pd.DataFrame(influence).to_csv(OUT/'bulk_adjusted_leave_one_donor_out.csv',index=False)
    values.to_csv(OUT/'bulk_matched_donor_logcpm.csv');md.to_csv(OUT/'bulk_matched_donor_metadata.csv')
    print('Bulk follow-up complete',flush=True)


if __name__=='__main__':main()
