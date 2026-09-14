#!/usr/bin/env python3
"""Independent numerical checks for the new scientific calculations."""
import json
from pathlib import Path
import numpy as np
import statsmodels.api as sm
from scipy import stats
from stats_core import bh, hc3_fit, wild_test, anchor_correlations, rank_axis1

rng=np.random.default_rng(20260913)
n=38
group=np.r_[np.zeros(8),np.ones(30)]
x=np.column_stack([np.ones(n),rng.normal(size=n),group])
y=rng.normal(size=(n,4))+group[:,None]*np.array([0,0,2,3])
fit=hc3_fit(y,x)
errors=[]
for j in range(y.shape[1]):
    reference=sm.OLS(y[:,j],x).fit(cov_type='HC3',use_t=True)
    errors.append(float(abs(reference.bse[-1]-fit['se'][j])))
    assert np.allclose(reference.params[-1],fit['effect'][j],rtol=1e-10,atol=1e-10)
    assert np.allclose(reference.bse[-1],fit['se'][j],rtol=1e-10,atol=1e-10)
assert np.allclose(bh([.01,.04,.03,.002]),[.02,.04,.04,.008])
a=rng.normal(size=(300,5)); r=anchor_correlations(a)
assert np.allclose(r,np.corrcoef(a,rowvar=False)[0])
ties=rng.integers(0,6,size=(50,28,12)).astype(float)
assert np.array_equal(rank_axis1(ties),stats.rankdata(ties,axis=1))
weights=rng.choice([-1.,1.],size=(1999,n))
wild=wild_test(y,x,weights)
assert wild['p_wild'][3]<.01
assert np.all(np.isfinite(wild['tstar']))
# A bounded null calibration check reports Monte Carlo variability rather than
# interpreting one random false positive as a failure of the implementation.
null_y=rng.normal(size=(n,300))
null_result=wild_test(null_y,x,weights[:999])
rate=float(np.mean(null_result['p_wild']<.05))
assert rate<.10, f'Unexpected null rejection rate {rate}'
report={'hc3_matches_statsmodels_max_se_error':max(errors),'bh_reference':'PASS',
        'bootstrap_average_ranks_match_scipy_with_ties':'PASS',
        'correlation_matches_numpy':'PASS','strong_signal_p':float(wild['p_wild'][3]),
        'null_tests':300,'null_rejection_rate_alpha_05':rate,
        'null_note':'One Gaussian design calibration only; does not establish universal finite-sample validity.'}
out=Path(__file__).resolve().parent/'statistics_validation.json'
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
