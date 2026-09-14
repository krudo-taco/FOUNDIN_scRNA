#!/usr/bin/env python3
"""Robustness diagnostics motivated by the initial method-dependent panel signal."""
import importlib.util
import itertools
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from stats_core import bh, design, hc3_fit, independent_columns, wild_test

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'
spec=importlib.util.spec_from_file_location('coupling',HERE/'03_within_donor_coupling.py')
coupling=importlib.util.module_from_spec(spec);spec.loader.exec_module(coupling)


def holm(p):
    p=np.asarray(p,float);o=np.argsort(p)
    v=np.maximum.accumulate(p[o]*(len(p)-np.arange(len(p))))
    r=np.empty(len(p));r[o]=np.minimum(v,1);return r


def main():
    genes=json.loads((OUT/'analysis_settings.json').read_text())['eligible_genes']
    corr=pd.read_csv(OUT/'culture_level_correlations.csv.gz')
    cov=pd.read_csv(OUT/'culture_context_coverage.csv')
    om=pd.read_csv(OUT/'within_donor_panel_omnibus.csv')
    selected=om.context.isin(['iDA_pooled','iDA_specificity'])
    om['p_holm_across_all_8_method_context_tests']=np.nan
    om.loc[selected,'p_holm_across_all_8_method_context_tests']=holm(om.loc[selected,'p_omnibus'])
    om.to_csv(OUT/'omnibus_multipipeline_sensitivity.csv',index=False)

    # A combined common-batch/reference-exclusion analysis avoids saturated
    # nuisance strata left when the reference donor is removed on its own.
    excluded=[3954,4106,3966]
    cohort=cov[~cov.PPMI_ID.isin(excluded)]
    shared=set(cohort.groupby('BATCH').group.nunique().loc[lambda x:x==2].index)
    sc=corr[~corr.PPMI_ID.isin(excluded)&corr.BATCH.isin(shared)]
    cv=cov[~cov.PPMI_ID.isin(excluded)&cov.BATCH.isin(shared)]
    exact_rows=[];exact_om=[];loo=[]
    for method in coupling.METHODS:
        da,md=coupling.aggregate(sc,cv,'iDA_pooled',method)
        non,nm=coupling.aggregate(sc,cv,'non_iDA_pooled',method)
        common=da.index.intersection(non.index)
        md=md.loc[common];da=da.loc[common];non=non.loc[common]
        assert (md.n_cultures==1).all() and (nm.loc[common].n_cultures==1).all()
        md['batch']=md[[f'batch_fraction_{b}' for b in range(1,6)]].idxmax(axis=1)
        base,design_names=design(md)
        strata=[]
        for _,idx in md.reset_index().groupby(['batch','genetic_sex']).groups.items():
            idx=list(idx);cases=int((md.iloc[idx].group=='iPD').sum())
            strata.append(list(itertools.combinations(idx,cases)))
        total=math.prod(len(s) for s in strata)
        assert total<=200000,'Exact enumeration unexpectedly large; revise with a documented Monte Carlo fallback.'
        labels=np.zeros((total,len(md)))
        for k,choice in enumerate(itertools.product(*strata)):
            for indices in choice:labels[k,list(indices)]=1
        for context,values in [('iDA_pooled',da),('iDA_specificity',da-non)]:
            vals=values.reindex(columns=genes)
            complete=vals.notna().all(axis=0)
            g=vals.columns[complete].tolist();y=vals.loc[:,complete].values
            fit=hc3_fit(y,base)
            null=np.empty((total,len(g)))
            for k,label in enumerate(labels):
                x=base.copy();x[:,-1]=label
                null[k]=hc3_fit(y,x)['t']
            p=(np.abs(null)>=np.abs(fit['t'])[None,:]-1e-12).mean(axis=0)
            panel_stat=np.mean(fit['t']**2)
            ppanel=float((np.mean(null**2,axis=1)>=panel_stat-1e-12).mean())
            exact_om.append(dict(method=method,context=context,n_HC=int((md.group=='HC').sum()),
                                 n_iPD=int((md.group=='iPD').sum()),permutation_states=total,
                                 shared_batches=','.join(str(b) for b in sorted(shared)),
                                 p_exact_omnibus=ppanel,statistic_mean_t2=panel_stat))
            for j,gene in enumerate(g):
                exact_rows.append(dict(method=method,context=context,gene=gene,effect_z=fit['effect'][j],
                                       p_exact=p[j],permutation_states=total))
            if method in ['rna_qc','spearman_qc']:
                for donor in md.index:
                    keep=md.index!=donor
                    rows,_=coupling.estimate(values.loc[keep],md.loc[keep],context,method,genes,resample=False)
                    rows['excluded_donor']=int(donor);rows['scenario']='common_batches_no_reference';loo.append(rows)
    exact=pd.DataFrame(exact_rows)
    exact['q_exact_method_family']=np.nan
    for _,idx in exact.groupby('method').groups.items():exact.loc[idx,'q_exact_method_family']=bh(exact.loc[idx,'p_exact'])
    exact.to_csv(OUT/'common_batch_exact_permutation_edges.csv',index=False)
    eo=pd.DataFrame(exact_om);eo['p_holm_all_methods_contexts']=holm(eo.p_exact_omnibus)
    eo.to_csv(OUT/'common_batch_exact_permutation_omnibus.csv',index=False)
    pd.concat(loo,ignore_index=True).to_csv(OUT/'common_batch_leave_one_donor_out.csv',index=False)

    # Detection conditioning is a diagnostic, not a newly selected primary model.
    detection_rows=[];detection_om=[]
    for method in coupling.METHODS:
        da,md=coupling.aggregate(corr,cov,'iDA_pooled',method)
        keep=~md.index.isin([3954,4106]);md=md.loc[keep];da=da.loc[keep]
        c=corr[(corr.context=='iDA_pooled')&(corr.method==method)&~corr.PPMI_ID.isin([3954,4106])]
        det=c.groupby(['PPMI_ID','gene'])[['detection','snca_detection']].mean()
        weights=np.random.default_rng(20260913).choice([-1.,1.],size=(19999,len(md)))
        nulls=[];observed=[]
        for gene in genes:
            y=da[gene].values
            gene_det=det.xs(gene,level='gene').reindex(md.index)
            x,names=design(md)
            technical=np.column_stack([np.log1p(md.n_cells.values),
                np.log((gene_det.snca_detection.values+1e-4)/(1-gene_det.snca_detection.values+1e-4)),
                np.log((gene_det.detection.values+1e-4)/(1-gene_det.detection.values+1e-4))])
            nuisance,_=independent_columns(np.column_stack([x[:,:-1],technical]),list(range(x.shape[1]+2)))
            x=np.column_stack([nuisance,x[:,-1]])
            valid=np.isfinite(y)&np.isfinite(x).all(axis=1)
            try:
                f=wild_test(y[valid],x[valid],weights[:,valid])
                detection_rows.append(dict(method=method,gene=gene,status='tested',effect_z=f['effect'][0],
                    p_wild=f['p_wild'][0],ci_low=f['ci_low'][0],ci_high=f['ci_high'][0]))
                nulls.append(f['tstar'][:,0]);observed.append(f['t'][0])
            except ValueError as e:
                detection_rows.append(dict(method=method,gene=gene,status=str(e),p_wild=np.nan))
        if nulls:
            null=np.mean(np.column_stack(nulls)**2,axis=1);obs=np.mean(np.array(observed)**2)
            detection_om.append(dict(method=method,genes_tested=len(observed),p_omnibus=(1+(null>=obs).sum())/(len(null)+1)))
    dr=pd.DataFrame(detection_rows);dr['q_wild_panel']=np.nan
    for _,idx in dr.groupby('method').groups.items():dr.loc[idx,'q_wild_panel']=bh(dr.loc[idx,'p_wild'])
    dr.to_csv(OUT/'detection_conditioned_diagnostic.csv',index=False)
    do=pd.DataFrame(detection_om);do['p_holm_methods']=holm(do.p_omnibus)
    do.to_csv(OUT/'detection_conditioned_omnibus.csv',index=False)
    print(eo.to_string(index=False));print(do.to_string(index=False))


if __name__=='__main__':main()
