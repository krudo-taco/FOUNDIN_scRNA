#!/usr/bin/env python3
"""Bounded null calibration using the actual primary donor design."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests
from stats_core import design, wild_test

HERE = Path(__file__).resolve().parent
meta = pd.read_csv(HERE / "results/donor_context_metadata.csv")
meta = meta[meta.context.eq("iDA_pooled") & ~meta.PPMI_ID.isin([3954, 4106])].copy()
assert meta.PPMI_ID.is_unique and len(meta) == 36
x, names = design(meta)
rng = np.random.default_rng(20260913)
weights = rng.choice([-1., 1.], size=(1999, len(meta)))
nulls = {
    "Gaussian_equal_variance": rng.normal(size=(len(meta), 500)),
    "Gaussian_HC_twice_SD": rng.normal(size=(len(meta), 500)) *
                            np.where(meta.group.eq("HC"), 2., 1.)[:, None],
    "heavy_tailed_df4": rng.standard_t(4, size=(len(meta), 500)),
}
calibration = {}
for label, y in nulls.items():
    result = wild_test(y, x, weights)
    rate = float(np.mean(result["p_wild"] < .05))
    assert rate < .10, (label, rate)
    calibration[label] = {"null_tests":500, "rejection_rate_alpha05":rate}
d = pd.read_csv(HERE / "results/within_donor_coupling_all_results.csv")
d = d[d.method.eq("rna_qc") & d.context.isin(["iDA_pooled", "iDA_specificity"])]
assert len(d) == 86 and d.status.eq("tested").all()
assert np.allclose(multipletests(d.p_wild, method="fdr_bh")[1], d.q_wild_family)
report = {
    "seed":20260913, "wild_draws":1999, "actual_design_columns":names,
    "actual_design_donors":len(meta), "actual_design_null_calibrations":calibration,
    "primary_edges":len(d), "primary_BH_matches_statsmodels":True,
    "limitation":"Bounded simulation checks, not a proof of universal finite-sample validity.",
}
(HERE / "actual_design_validation.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps(report, indent=2))
