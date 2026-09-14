"""Guard partial correlations against constant-input and projection roundoff artifacts."""
import numpy as np
from stats_core import anchor_correlations, residualize


def partial_correlations(values, covariates):
    values = np.asarray(values, float)
    residual = residualize(values, covariates)
    original_norm = np.linalg.norm(values - values.mean(axis=0), axis=0)
    residual_norm = np.linalg.norm(residual - residual.mean(axis=0), axis=0)
    tolerance = 100 * np.finfo(float).eps * max(values.shape[0], covariates.shape[1])
    valid = (original_norm > 0) & (residual_norm > tolerance * np.maximum(original_norm, 1))
    result = anchor_correlations(residual)
    result[~valid] = np.nan
    if not valid[0]:
        result[:] = np.nan
    return result
