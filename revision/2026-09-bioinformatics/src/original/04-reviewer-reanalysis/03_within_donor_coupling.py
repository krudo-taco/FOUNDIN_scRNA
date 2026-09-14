#!/usr/bin/env python3
"""Individualized co-expression and direct DA-versus-non-DA interaction tests."""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from stats_core import anchor_correlations, bh, design, fisher, hc3_fit, residualize, wild_test

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
INPUT=HERE/'inputs'
OUT=HERE/'results'
OUT.mkdir(exist_ok=True)
SEED=20260913
N_WILD=19999
MIN_CELLS=100
METHODS=['rna_qc','rna_raw','spearman_qc','sct_qc']
COHORT_EXCLUSIONS=[3954,4106]


def log(message):
    print(time.strftime('%H:%M:%S'),message,flush=True)


def get_panel(columns):
    original=pd.read_csv('/data32TB/shared/ppmi-foundin/analyse/RNA-processed/shaohua-gene-list-2-10-25.csv').gene.tolist()
    mapped=['PRKN' if g=='PARK2' and 'PRKN' in columns else g for g in original]
    assert len(mapped)==len(set(mapped))
    return [g for g in mapped if g in columns]


def calculate_sample_correlations():
    meta=pd.read_csv(INPUT/'cell_metadata_full.tsv.gz',sep='\t',index_col='cell')
    meta=meta[meta.subtype.isin(['na_HC','na_iPD'])].copy()
    meta['group']=meta.subtype.map({'na_HC':'HC','na_iPD':'iPD'})
    expr=pd.read_csv(INPUT/'rna_panel_logexpr.tsv.gz',sep='\t',index_col='cell').loc[meta.index]
    sct=pd.read_csv(INPUT/'sct_panel_logexpr.tsv.gz',sep='\t',index_col='cell').loc[meta.index]
    panel=get_panel(expr.columns)
    genes=['SNCA']+panel
    primary_da=meta.CellType.str.startswith('iDA') & ~meta.PPMI_ID.isin(COHORT_EXCLUSIONS)
    donor_detection=(expr.loc[primary_da,panel]>0).groupby(meta.loc[primary_da,'PPMI_ID']).mean()
    filter_table=pd.DataFrame({'gene':panel,'mean_donor_detection_pooled_iDA':donor_detection.mean().reindex(panel).values})
    filter_table['eligible_primary']=filter_table.mean_donor_detection_pooled_iDA>=.05
    filter_table.to_csv(OUT/'gene_panel_eligibility.csv',index=False)
    eligible=filter_table.loc[filter_table.eligible_primary,'gene'].tolist()
    (OUT/'analysis_settings.json').write_text(json.dumps(dict(seed=SEED,wild_draws=N_WILD,
        min_cells_per_culture_context=MIN_CELLS,original_panel_genes=panel,eligible_genes=eligible,
        primary_cohort_excluded_donors=COHORT_EXCLUSIONS,
        cohort_exclusion_reason='Explicit genetic-PD labels in source workbook conflict with missing Seurat subtype fields; original grouping retained as sensitivity.',
        eligibility='At least 5% mean donor detection in pooled iDA, with HC/iPD labels ignored.',
        primary='RNA log expression; partial cellular correlation adjusted for log library size, mitochondrial fraction and constituent cell type.',
        donor_model='Equal-culture mean Fisher z per donor, disease + sex + equal-culture batch fractions.',
        main_edge_family='All eligible genes in pooled iDA and donor-paired iDA-minus-non-iDA contrast.',
        secondary_edge_family='All eligible genes across the eleven annotated cell types.',
        omnibus='Mean squared HC3 t, calibrated with common donor wild-bootstrap weights; Holm over the two primary omnibus tests.'),indent=2)+'\n')
    log(f'Locked label-blind panel: {len(eligible)}/{len(panel)} genes; {len(meta)} HC/iPD cells')

    rows=[]; coverage=[]
    contexts=['iDA_pooled','non_iDA_pooled']+sorted(meta.CellType.unique())
    expr=expr[genes]
    # Rare genes absent from SCT are retained as explicitly missing sensitivity values.
    sct=sct.reindex(columns=genes)
    for context in contexts:
        select=(meta.CellType.str.startswith('iDA') if context=='iDA_pooled' else
                ~meta.CellType.str.startswith('iDA') if context=='non_iDA_pooled' else meta.CellType==context)
        cm=meta.loc[select]
        for sample,sm in cm.groupby('SampleID',sort=True):
            first=sm.iloc[0]; n=len(sm)
            identity=dict(context=context,sample=sample,PPMI_ID=int(first.PPMI_ID),group=first.group,
                          genetic_sex=int(first.genetic_sex),BATCH=int(first.BATCH),n_cells=n)
            assert sm.PPMI_ID.nunique()==sm.genetic_sex.nunique()==sm.BATCH.nunique()==sm.group.nunique()==1
            coverage.append(dict(identity,included=n>=MIN_CELLS))
            if n<MIN_CELLS:
                continue
            a=expr.loc[sm.index].to_numpy(float)
            st=sct.loc[sm.index].to_numpy(float)
            technical=np.column_stack([np.ones(n),np.log1p(sm.nCount_RNA.values),sm['percent.mt'].values/100])
            if context in ['iDA_pooled','non_iDA_pooled']:
                technical=np.column_stack([technical,pd.get_dummies(sm.CellType,drop_first=True).values])
            ranked=np.column_stack([stats.rankdata(a[:,j]) for j in range(a.shape[1])])
            ranked_cov=technical.copy()
            ranked_cov[:,1]=stats.rankdata(technical[:,1])
            ranked_cov[:,2]=stats.rankdata(technical[:,2])
            sct_ok=np.all(np.isfinite(st),axis=0)
            st_res=np.full(st.shape,np.nan)
            st_res[:,sct_ok]=residualize(st[:,sct_ok],technical)
            matrices={'rna_raw':a,'rna_qc':residualize(a,technical),
                      'spearman_qc':residualize(ranked,ranked_cov),'sct_qc':st_res}
            for method,values in matrices.items():
                r=anchor_correlations(values)
                for j,g in enumerate(panel,1):
                    rows.append(dict(identity,method=method,gene=g,r=float(r[j]),z=float(fisher(r[j])),
                                     detection=float(np.mean(a[:,j]>0)),snca_detection=float(np.mean(a[:,0]>0))))
        log(f'Calculated culture correlations for {context}')
    corr=pd.DataFrame(rows)
    cov=pd.DataFrame(coverage)
    corr.to_csv(OUT/'culture_level_correlations.csv.gz',index=False)
    cov.to_csv(OUT/'culture_context_coverage.csv',index=False)
    return corr,cov,eligible


def aggregate(corr,cov,context,method):
    records=corr[(corr.context==context)&(corr.method==method)].copy()
    values=records.groupby(['PPMI_ID','gene']).z.mean().unstack('gene')
    md=cov[(cov.context==context)&cov.included].copy()
    info=md.groupby('PPMI_ID').agg(group=('group','first'),genetic_sex=('genetic_sex','first'),
                                    n_cultures=('sample','nunique'),n_cells=('n_cells','sum'))
    batch=pd.crosstab(md.PPMI_ID,md.BATCH).reindex(columns=[1,2,3,4,5],fill_value=0)
    batch=batch.div(batch.sum(axis=1),axis=0)
    for b in range(1,6):
        info[f'batch_fraction_{b}']=batch[b]
    return values,info


def estimate(values,metadata,label,method,genes,resample=True,adjusted=True):
    values=values.reindex(index=metadata.index,columns=genes)
    rng=np.random.default_rng(SEED)
    weights=rng.choice([-1.,1.],size=(N_WILD,len(metadata))) if resample else None
    groups={}
    rows=[]; nulls=[]; observed=[]
    for g in genes:
        valid=np.isfinite(values[g].values)
        key=tuple(valid)
        groups.setdefault(key,[]).append(g)
    for key,subset in groups.items():
        valid=np.array(key)
        md=metadata.loc[valid]
        nh=int((md.group=='HC').sum()); npd=int((md.group=='iPD').sum())
        y=values.loc[valid,subset].values
        x,names=design(md,adjusted=adjusted)
        common=dict(context=label,method=method,n_HC=nh,n_iPD=npd,model=' + '.join(names))
        try:
            if min(nh,npd)<5:
                raise ValueError('Fewer than five independent donors in an arm')
            fit=wild_test(y,x,weights[:,valid]) if resample else hc3_fit(y,x)
            if not resample:
                fit['p_hc3_t']=2*stats.t.sf(abs(fit['t']),fit['df'])
                crit=stats.t.ppf(.975,fit['df'])
                fit['ci_low']=fit['effect']-crit*fit['se']; fit['ci_high']=fit['effect']+crit*fit['se']
            for j,g in enumerate(subset):
                row=dict(common,gene=g,status='tested',effect_z=float(fit['effect'][j]),
                         se_hc3=float(fit['se'][j]),t_hc3=float(fit['t'][j]),df=fit['df'],
                         ci_low=float(fit['ci_low'][j]),ci_high=float(fit['ci_high'][j]),
                         p_hc3_t=float(fit['p_hc3_t'][j]),max_leverage=float(max(fit['leverage'])),
                         p_wild=float(fit['p_wild'][j]) if resample else np.nan)
                if label!='iDA_specificity':
                    row['mean_r_HC']=float(np.mean(np.tanh(y[md.group.values=='HC',j])))
                    row['mean_r_iPD']=float(np.mean(np.tanh(y[md.group.values=='iPD',j])))
                rows.append(row)
                if resample:
                    observed.append(fit['t'][j]); nulls.append(fit['tstar'][:,j])
        except ValueError as exc:
            for g in subset:
                rows.append(dict(common,gene=g,status=str(exc),p_wild=np.nan,p_hc3_t=np.nan))
    omnibus=dict(context=label,method=method,genes_tested=len(observed),p_omnibus=np.nan)
    if nulls:
        t0=np.column_stack(nulls)
        obs=float(np.mean(np.square(observed)))
        null=np.mean(t0**2,axis=1)
        omnibus.update(statistic_mean_t2=obs,p_omnibus=float((1+(null>=obs).sum())/(len(null)+1)))
    return pd.DataFrame(rows),omnibus


def analyse(corr,cov,genes):
    results=[]; omnibus=[]; donor_values=[]; donor_meta=[]
    contexts=['iDA_pooled','non_iDA_pooled']+sorted(cov[~cov.context.isin(['iDA_pooled','non_iDA_pooled'])].context.unique())
    cache={}; legacy={}
    for method in METHODS:
        for context in contexts:
            values,md=aggregate(corr,cov,context,method)
            legacy[(context,method)]=(values,md)
            keep=~md.index.isin(COHORT_EXCLUSIONS)
            values=values.loc[keep];md=md.loc[keep]
            cache[(context,method)]=(values,md)
            if method=='rna_qc':
                v=values.copy();v['PPMI_ID']=v.index;v['context']=context
                donor_values.append(v)
                m=md.copy();m['PPMI_ID']=m.index;m['context']=context;donor_meta.append(m)
            if method!='rna_qc' and context not in ['iDA_pooled','non_iDA_pooled']:
                continue
            run_wild=context!='non_iDA_pooled'
            r,o=estimate(values,md,context,method,genes,resample=run_wild)
            results.append(r)
            if run_wild: omnibus.append(o)
            log(f'Fit {context}, {method}: {sum(r.status=="tested")} edges')
        da,md=cache[('iDA_pooled',method)]
        non,md_non=cache[('non_iDA_pooled',method)]
        common=da.index.intersection(non.index)
        # Both summaries must represent the same culture mix for the paired estimand.
        consistent=(md.loc[common,[f'batch_fraction_{b}' for b in range(1,6)]].values ==
                    md_non.loc[common,[f'batch_fraction_{b}' for b in range(1,6)]].values).all(axis=1)
        common=common[consistent]
        diff=da.loc[common]-non.loc[common]
        cache[('iDA_specificity',method)]=(diff,md.loc[common])
        r,o=estimate(diff,md.loc[common],'iDA_specificity',method,genes)
        results.append(r);omnibus.append(o)
    res=pd.concat(results,ignore_index=True)
    res['q_wild_family']=np.nan
    primary=res.context.isin(['iDA_pooled','iDA_specificity'])
    for method in METHODS:
        idx=primary & (res.method==method)
        res.loc[idx,'q_wild_family']=bh(res.loc[idx,'p_wild'])
    secondary=(res.method=='rna_qc') & ~res.context.isin(['iDA_pooled','iDA_specificity','non_iDA_pooled'])
    res.loc[secondary,'q_wild_family']=bh(res.loc[secondary,'p_wild'])
    res.to_csv(OUT/'within_donor_coupling_all_results.csv',index=False)
    om=pd.DataFrame(omnibus)
    om['p_holm_primary']=np.nan
    for method in METHODS:
        idx=(om.method==method)&om.context.isin(['iDA_pooled','iDA_specificity'])
        p=om.loc[idx,'p_omnibus'].values
        order=np.argsort(p); adj=np.maximum.accumulate(p[order]*(len(p)-np.arange(len(p))))
        restored=np.empty(len(p)); restored[order]=np.minimum(adj,1)
        om.loc[idx,'p_holm_primary']=restored
    secondary=(om.method=='rna_qc')&~om.context.isin(['iDA_pooled','iDA_specificity'])
    om['q_secondary']=np.nan;om.loc[secondary,'q_secondary']=bh(om.loc[secondary,'p_omnibus'])
    om.to_csv(OUT/'within_donor_panel_omnibus.csv',index=False)
    pd.concat(donor_values).to_csv(OUT/'donor_individualized_fisher_z.csv',index=False)
    pd.concat(donor_meta).to_csv(OUT/'donor_context_metadata.csv',index=False)

    sensitivity=[];influence=[]
    for context in ['iDA_pooled','iDA_specificity']:
        values,md=cache[(context,'rna_qc')]
        for scenario,mask in [('exclude_reference_donor',md.index!=3966),
                              ('unadjusted_donor_model',np.ones(len(md),bool))]:
            r,o=estimate(values.loc[mask],md.loc[mask],context,'rna_qc',genes,
                         adjusted=scenario!='unadjusted_donor_model')
            r['scenario']=scenario;sensitivity.append(r)
            o['scenario']=scenario
            log(f'Sensitivity {context}: {scenario}, omnibus P={o["p_omnibus"]:.4g}')
            with (OUT/'sensitivity_omnibus.jsonl').open('a') as f: f.write(json.dumps(o)+'\n')
        primary_cultures=cov[~cov.PPMI_ID.isin(COHORT_EXCLUSIONS)]
        shared=set(primary_cultures.groupby('BATCH').group.nunique().loc[lambda x:x==2].index)
        shared_corr=corr[corr.BATCH.isin(shared) & ~corr.PPMI_ID.isin(COHORT_EXCLUSIONS)]
        shared_cov=cov[cov.BATCH.isin(shared) & ~cov.PPMI_ID.isin(COHORT_EXCLUSIONS)]
        shared_values,shared_md=aggregate(shared_corr,shared_cov,'iDA_pooled','rna_qc')
        if context=='iDA_specificity':
            shared_non,shared_non_md=aggregate(shared_corr,shared_cov,'non_iDA_pooled','rna_qc')
            common=shared_values.index.intersection(shared_non.index)
            consistent=(shared_md.loc[common,[f'batch_fraction_{b}' for b in range(1,6)]].values ==
                        shared_non_md.loc[common,[f'batch_fraction_{b}' for b in range(1,6)]].values).all(axis=1)
            common=common[consistent]
            shared_values=shared_values.loc[common]-shared_non.loc[common]
            shared_md=shared_md.loc[common]
        r,o=estimate(shared_values,shared_md,context,'rna_qc',genes)
        r['scenario']='shared_batches_only';sensitivity.append(r)
        o.update(scenario='shared_batches_only',included_batches=sorted(shared))
        with (OUT/'sensitivity_omnibus.jsonl').open('a') as f: f.write(json.dumps(o)+'\n')
        legacy_da,legacy_md=legacy[('iDA_pooled','rna_qc')]
        if context=='iDA_specificity':
            legacy_non,legacy_non_md=legacy[('non_iDA_pooled','rna_qc')]
            common=legacy_da.index.intersection(legacy_non.index)
            consistent=(legacy_md.loc[common,[f'batch_fraction_{b}' for b in range(1,6)]].values ==
                        legacy_non_md.loc[common,[f'batch_fraction_{b}' for b in range(1,6)]].values).all(axis=1)
            common=common[consistent]
            legacy_da=legacy_da.loc[common]-legacy_non.loc[common]
            legacy_md=legacy_md.loc[common]
        r,o=estimate(legacy_da,legacy_md,context,'rna_qc',genes)
        r['scenario']='original_30_iPD_grouping';sensitivity.append(r)
        o['scenario']='original_30_iPD_grouping'
        with (OUT/'sensitivity_omnibus.jsonl').open('a') as f: f.write(json.dumps(o)+'\n')
        for donor in md.index:
            mask=md.index!=donor
            r,_=estimate(values.loc[mask],md.loc[mask],context,'rna_qc',genes,resample=False)
            r['excluded_donor']=int(donor);influence.append(r)
    sens=pd.concat(sensitivity,ignore_index=True)
    sens['q_wild_family']=np.nan
    for _,idx in sens.groupby('scenario').groups.items():
        sens.loc[idx,'q_wild_family']=bh(sens.loc[idx,'p_wild'])
    sens.to_csv(OUT/'within_donor_sensitivity.csv',index=False)
    pd.concat(influence,ignore_index=True).to_csv(OUT/'leave_one_donor_out.csv',index=False)
    log('Within-donor analyses complete')


def main():
    if (OUT/'sensitivity_omnibus.jsonl').exists():
        (OUT/'sensitivity_omnibus.jsonl').unlink()
    if (OUT/'culture_level_correlations.csv.gz').exists():
        corr=pd.read_csv(OUT/'culture_level_correlations.csv.gz')
        cov=pd.read_csv(OUT/'culture_context_coverage.csv')
        genes=json.loads((OUT/'analysis_settings.json').read_text())['eligible_genes']
    else:
        corr,cov,genes=calculate_sample_correlations()
    analyse(corr,cov,genes)


if __name__=='__main__':
    main()
