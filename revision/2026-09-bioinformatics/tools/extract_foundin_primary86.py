#!/usr/bin/env python3
"""Extract the 86 primary FOUNDIN edges from a staged all-results table."""
import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-results", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    all_results = pd.read_csv(args.all_results)
    reference = pd.read_csv(args.reference)
    keys = ["context", "method", "gene"]
    selected = reference[keys].merge(all_results, on=keys, how="left", validate="one_to_one")
    if len(selected) != 86 or selected.status.isna().any():
        raise SystemExit("Failed to extract all 86 primary edges")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    selected[reference.columns].to_csv(args.out, index=False)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()

