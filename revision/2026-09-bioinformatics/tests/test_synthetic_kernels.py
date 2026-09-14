#!/usr/bin/env python3
"""Synthetic non-identifying checks for statistical kernels."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "repro_package"))
from numerics import partial_correlations
from stats_core import bh, hc3_fit, independent_columns, wild_test


def assert_close(actual, expected, tol=1e-10):
    if not np.allclose(actual, expected, atol=tol, rtol=tol, equal_nan=True):
        raise AssertionError(f"{actual!r} != {expected!r}")


def test_bh_reference():
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.20])
    assert_close(bh(p), np.array([0.04, 0.05333333333333334, np.nan, 0.05333333333333334, 0.2]))


def test_independent_columns_drops_redundant_column():
    x = np.column_stack([np.ones(6), np.arange(6), 2 * np.arange(6)])
    kept, names = independent_columns(x, ["intercept", "x", "two_x"])
    if kept.shape[1] != 2 or names != ["intercept", "x"]:
        raise AssertionError(names)


def test_hc3_and_wild_detect_synthetic_signal():
    rng = np.random.default_rng(20260913)
    n = 18
    group = np.r_[np.zeros(9), np.ones(9)]
    x = np.column_stack([np.ones(n), group])
    y = 0.2 + 1.1 * group + rng.normal(0, 0.25, n)
    fit = hc3_fit(y, x)
    if fit["effect"][0] <= 0.8 or fit["se"][0] <= 0:
        raise AssertionError(fit)
    weights = rng.choice([-1.0, 1.0], size=(999, n))
    wt = wild_test(y, x, weights)
    if wt["p_wild"][0] >= 0.05:
        raise AssertionError(wt["p_wild"])


def test_partial_correlation_constant_guard():
    values = np.column_stack(
        [
            np.arange(10, dtype=float),
            np.arange(10, dtype=float) + np.linspace(0, 0.2, 10),
            np.ones(10),
        ]
    )
    cov = np.ones((10, 1))
    result = partial_correlations(values, cov)
    if not np.isfinite(result[1]):
        raise AssertionError(result)
    if not np.isnan(result[2]):
        raise AssertionError(result)


def main():
    tests = [
        test_bh_reference,
        test_independent_columns_drops_redundant_column,
        test_hc3_and_wild_detect_synthetic_signal,
        test_partial_correlation_constant_guard,
    ]
    for test in tests:
        test()
    print(json.dumps({"status": "PASS", "synthetic_tests": len(tests)}, indent=2))


if __name__ == "__main__":
    main()

