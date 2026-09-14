"""Donor-level inference primitives with explicit small-sample assumptions."""
import numpy as np
from scipy import stats


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    order = np.argsort(p[valid])
    if not len(order):
        return out
    q = np.minimum.accumulate((p[valid][order] * len(order) / np.arange(1, len(order)+1))[::-1])[::-1]
    restored = np.empty(len(order))
    restored[order] = np.minimum(q, 1)
    out[valid] = restored
    return out


def independent_columns(x, names):
    kept = []
    for j in range(x.shape[1]):
        if np.linalg.matrix_rank(x[:, kept + [j]]) > len(kept):
            kept.append(j)
    return x[:, kept], [names[j] for j in kept]


def design(metadata, adjusted=True):
    cols = [np.ones(len(metadata))]
    names = ['intercept']
    if adjusted:
        cols.append((metadata.genetic_sex.values == 2).astype(float))
        names.append('female')
        for b in [2, 3, 4, 5]:
            key = f'batch_fraction_{b}'
            if key in metadata:
                cols.append(metadata[key].values.astype(float)); names.append(key)
    x0, names = independent_columns(np.column_stack(cols), names)
    group = (metadata.group.values == 'iPD').astype(float)
    x = np.column_stack([x0, group])
    return x, names + ['iPD']


def hc3_fit(y, x):
    """Return last-coefficient estimate, HC3 SE and t; rows are independent donors."""
    y = np.asarray(y, float)
    if y.ndim == 1:
        y = y[:, None]
    n, k = x.shape
    if np.linalg.matrix_rank(x) != k or n-k < 4:
        raise ValueError('Non-estimable design or fewer than four residual degrees of freedom')
    pinv = np.linalg.pinv(x)
    h = np.sum(x * pinv.T, axis=1)
    if np.max(h) >= 1-1e-8:
        raise ValueError('Saturated donor leverage prevents HC3 inference')
    beta = pinv @ y
    resid = y - x @ beta
    a = pinv[-1]
    se = np.sqrt(np.sum((a[:, None] * resid / (1-h[:, None]))**2, axis=0))
    t = np.divide(beta[-1], se, out=np.full(y.shape[1], np.nan), where=se>0)
    return dict(effect=beta[-1], se=se, t=t, df=n-k, leverage=h, residual=resid, pinv=pinv)


def wild_test(y, x, weights):
    """Restricted-null Rademacher wild bootstrap with HC2 residual scaling and HC3 t.

    A common donor weight vector is used for every gene, preserving cross-gene
    dependence. This is an approximate small-sample procedure, not an exact
    randomized treatment test. Fixed nuisance covariates are held fixed.
    """
    fit = hc3_fit(y, x)
    y = np.asarray(y, float)
    if y.ndim == 1:
        y = y[:, None]
    x0 = x[:, :-1]
    p0 = np.linalg.pinv(x0)
    h0 = np.sum(x0*p0.T, axis=1)
    residual0 = (y-x0@(p0@y))/np.sqrt(1-h0[:, None])
    a = fit['pinv'][-1]
    projection = np.eye(len(x))-x@fit['pinv']
    denominator_weights = (a/(1-fit['leverage']))**2
    b = len(weights)
    tstar = np.empty((b, y.shape[1]))
    for start in range(0, b, 256):
        w = weights[start:start+256]
        draws = w[:, :, None] * residual0[None, :, :]
        effect = np.einsum('n,bng->bg', a, draws)
        residual = np.einsum('ij,bjg->big', projection, draws, optimize=True)
        se = np.sqrt(np.einsum('n,bng->bg', denominator_weights, residual**2))
        tstar[start:start+len(w)] = np.divide(effect,se,out=np.zeros_like(effect),where=se>0)
    fit['p_wild'] = (1+np.sum(np.abs(tstar)>=np.abs(fit['t'])[None, :],axis=0))/(b+1)
    fit['p_hc3_t'] = 2*stats.t.sf(np.abs(fit['t']),fit['df'])
    critical = stats.t.ppf(.975,fit['df'])
    fit['ci_low'] = fit['effect']-critical*fit['se']
    fit['ci_high'] = fit['effect']+critical*fit['se']
    fit['tstar'] = tstar
    return fit


def residualize(a, cov):
    cov, _ = independent_columns(np.asarray(cov,float),list(range(cov.shape[1])))
    return a-cov@np.linalg.lstsq(cov,a,rcond=None)[0]


def anchor_correlations(a, anchor=0):
    centered = a-np.mean(a,axis=0)
    scale = np.sqrt(np.sum(centered**2,axis=0))
    denom = scale[anchor]*scale
    return np.divide(centered[:,anchor]@centered,denom,
                     out=np.full(a.shape[1],np.nan),where=denom>1e-12)


def fisher(r):
    return np.arctanh(np.clip(r,-1+1e-7,1-1e-7))


def rank_axis1(a):
    """Average ranks over donors for a bootstrap x donor x gene array."""
    shape=a.shape
    flat=a.transpose(0,2,1).reshape(-1,shape[1])
    order=np.argsort(flat,axis=1,kind='stable')
    sorted_values=np.take_along_axis(flat,order,axis=1)
    first=np.column_stack([np.ones(len(flat),bool),sorted_values[:,1:]!=sorted_values[:,:-1]])
    last=np.column_stack([sorted_values[:,:-1]!=sorted_values[:,1:],np.ones(len(flat),bool)])
    positions=np.arange(1,shape[1]+1)[None,:]
    lo=np.maximum.accumulate(np.where(first,positions,0),axis=1)
    hi=np.minimum.accumulate(np.where(last,positions,shape[1]+1)[:,::-1],axis=1)[:,::-1]
    ranked=np.empty(flat.shape,float)
    np.put_along_axis(ranked,order,(lo+hi)/2,axis=1)
    return ranked.reshape(shape[0],shape[2],shape[1]).transpose(0,2,1)
