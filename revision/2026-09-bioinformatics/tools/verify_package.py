#!/usr/bin/env python3
"""Static package checks that do not require restricted data."""
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def require(path):
    if not (ROOT / path).exists():
        raise SystemExit(f"Missing required package file: {path}")


def sha256(path):
    h = hashlib.sha256()
    with open(ROOT / path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    required = [
        "CLEANUP_PLAN.md",
        "README.md",
        "requirements.txt",
        "run_public_gse243639.py",
        "src/repro_package/stats_core.py",
        "src/repro_package/numerics.py",
        "src/original/05-claim-support/analysis_plan.md",
        "src/original/05-claim-support/inputs/kamath_pre_outcome_plan_sha256.txt",
        "manifests/data_access_manifest.csv",
        "manifests/package_file_manifest.csv",
        "manifests/required_input_hashes.json",
        "ci/github-actions/revision-2026-09-bioinformatics.yml",
        "reference_settings/analysis_settings.json",
        "reference_settings/foundin_regrouping_settings.json",
    ]
    for path in required:
        require(path)
    routes = {row["access_route"] for row in csv.DictReader(open(ROOT / "manifests/data_access_manifest.csv"))}
    expected_routes = {"restricted", "reused public source", "derived aggregate", "private"}
    if not expected_routes <= routes:
        raise SystemExit(f"Missing data-access routes: {sorted(expected_routes - routes)}")
    plan_hash = sha256("src/original/05-claim-support/analysis_plan.md")
    recorded = (ROOT / "src/original/05-claim-support/inputs/kamath_pre_outcome_plan_sha256.txt").read_text().split()[0]
    if plan_hash != recorded:
        raise SystemExit(f"Kamath frozen plan hash mismatch: {plan_hash} != {recorded}")
    manifest_rows = list(csv.DictReader(open(ROOT / "manifests/package_file_manifest.csv")))
    checked_manifest = 0
    for row in manifest_rows:
        if row["path"] == "manifests/package_file_manifest.csv":
            raise SystemExit("package_file_manifest.csv must not include itself")
        path = ROOT / row["path"]
        if not path.exists():
            raise SystemExit(f"Manifest file missing: {row['path']}")
        if str(path.stat().st_size) != row["bytes"]:
            raise SystemExit(f"Manifest size mismatch: {row['path']}")
        observed = sha256(row["path"])
        if observed != row["sha256"]:
            raise SystemExit(f"Manifest hash mismatch for {row['path']}: {observed} != {row['sha256']}")
        checked_manifest += 1
    result = {
        "status": "PASS",
        "required_files": len(required),
        "data_access_routes": sorted(routes),
        "copied_05_analysis_plan_sha256": plan_hash,
        "kamath_frozen_plan_recorded_sha256": recorded,
        "manifest_files_checked": checked_manifest,
        "verification_level": "static package integrity; not a full biological rerun",
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
