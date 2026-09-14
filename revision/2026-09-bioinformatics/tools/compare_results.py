#!/usr/bin/env python3
"""Numerically compare rerun CSV outputs against staged reference outputs."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PROFILES = {
    "foundin_primary86": {
        "keys": ["context", "method", "gene"],
        "numeric": ["effect_z", "se_hc3", "t_hc3", "ci_low", "ci_high", "p_wild", "p_hc3_t", "q_wild_family"],
        "required_rows": 86,
    },
    "gse243639_edges264": {
        "keys": ["context", "method", "model", "test", "gene"],
        "numeric": ["effect_z", "ci_low", "ci_high", "se_hc3", "t_hc3", "p_wild", "p_hc3_t", "q_joint_pilot_family"],
        "required_rows": 264,
    },
    "gse243639_panel16": {
        "keys": ["context", "method", "model", "test"],
        "numeric": ["genes_tested", "statistic_mean_t2", "p_omnibus", "p_holm_all_pilot_tests"],
        "required_rows": 16,
    },
    "kamath_edges1008": {
        "keys": ["scenario", "context", "method", "test", "gene"],
        "numeric": ["effect_z", "ci_low", "ci_high", "se_hc3", "t_hc3", "p_wild", "p_hc3_t", "q_BH_planned_edges"],
        "required_rows": 1008,
    },
    "kamath_panel48": {
        "keys": ["scenario", "context", "method", "test"],
        "numeric": ["genes_tested", "statistic_mean_t2", "p_omnibus", "p_Holm_Kamath_48", "p_Holm_primary_4", "p_Holm_external_64"],
        "required_rows": 48,
    },
}


def compare(profile, observed_path, reference_path, tolerance):
    spec = PROFILES[profile]
    observed = pd.read_csv(observed_path)
    reference = pd.read_csv(reference_path)
    if len(observed) != spec["required_rows"] or len(reference) != spec["required_rows"]:
        raise AssertionError((len(observed), len(reference), spec["required_rows"]))
    merged = observed.merge(reference, on=spec["keys"], suffixes=("_observed", "_reference"), validate="one_to_one")
    if len(merged) != spec["required_rows"]:
        raise AssertionError(f"Key mismatch: {len(merged)} != {spec['required_rows']}")
    checks = []
    for column in spec["numeric"]:
        left_name, right_name = f"{column}_observed", f"{column}_reference"
        if left_name not in merged or right_name not in merged:
            raise AssertionError(f"Missing numeric column {column}")
        left = merged[left_name].fillna(-999999999.0).astype(float)
        right = merged[right_name].fillna(-999999999.0).astype(float)
        diff = (left - right).abs()
        ok = bool((diff <= tolerance).all())
        checks.append({"column": column, "match": ok, "max_abs_diff": float(diff.max())})
    status = "PASS" if all(item["match"] for item in checks) else "FAIL"
    return {"status": status, "profile": profile, "rows": int(len(merged)), "tolerance": tolerance, "checks": checks}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", required=True, choices=sorted(PROFILES))
    parser.add_argument("--observed", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--out")
    parser.add_argument("--tolerance", type=float, default=1e-12)
    args = parser.parse_args()
    result = compare(args.profile, Path(args.observed), Path(args.reference), args.tolerance)
    print(json.dumps(result, indent=2))
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(result, indent=2) + "\n")
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()

