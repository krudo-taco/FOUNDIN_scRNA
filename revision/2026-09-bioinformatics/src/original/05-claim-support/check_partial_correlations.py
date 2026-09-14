"""Independent regression checks for the observed constant-rank correlation artifact."""
import sys
from pathlib import Path
import numpy as np
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / '04-reviewer-reanalysis'))
from numerics import partial_correlations

rng = np.random.default_rng(1893)
n = 65
technical = rng.normal(size=n)
anchor = 2 * technical + rng.normal(size=n)
gene = -.7 * technical + .8 * anchor + rng.normal(size=n)
x = np.column_stack([np.ones(n), technical])
a = np.column_stack([anchor, gene, np.zeros(n), np.repeat((n + 1) / 2, n), 4 + 3 * technical])
r = partial_correlations(a, x)
assert np.isnan(r[2:]).all(), 'Constant and completely confounded variables must have undefined correlation'
ra = stats.linregress(technical, anchor)
rg = stats.linregress(technical, gene)
expected = stats.pearsonr(anchor - ra.intercept - ra.slope * technical,
                        gene - rg.intercept - rg.slope * technical)[0]
assert np.isclose(r[1], expected, atol=1e-12)
assert np.isclose(partial_correlations(a * 100, x)[1], expected, atol=1e-12)
a[:, 0] = 0
assert np.isnan(partial_correlations(a, x)).all()
print('PASS: constants, constant average ranks, complete confounding, independent Pearson check and scale invariance')
