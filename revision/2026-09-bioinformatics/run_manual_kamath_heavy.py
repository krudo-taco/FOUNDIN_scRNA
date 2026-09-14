#!/usr/bin/env python3
"""Stage public Kamath inputs and run preprocessing plus outcome verification."""
import argparse
import json
import sys
from pathlib import Path


PACKAGE = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE / "src" / "repro_package"))
from staging import copy_file, copy_tree, ensure_clean_dir, load_json, run_command, verify_hash_records


INPUT_NAMES = [
    "kamath_panel_counts.npz",
    "kamath_cell_qc.csv.gz",
    "kamath_cell_annotation_mapping.csv.gz",
    "kamath_author_metadata.csv.gz",
    "kamath_published_donor_metadata.csv",
    "kamath_GEO_library_manifest.csv",
    "kamath_library_alias_map.csv",
    "kamath_metadata_audit.json",
    "kamath_independent_input_validation.json",
    "kamath_count_extraction_audit.json",
    "kamath_annotation_mapping_audit.json",
    "kamath_pre_outcome_plan_sha256.txt",
]


def stage_sources(config):
    work_root = (PACKAGE / config["work_root"]).resolve()
    stage_root = ensure_clean_dir(work_root / "stage")
    logs = work_root / "logs"
    stage_04 = stage_root / "04-reviewer-reanalysis"
    stage_05 = stage_root / "05-claim-support"
    copy_tree(PACKAGE / "src" / "original" / "04-reviewer-reanalysis", stage_04)
    copy_tree(PACKAGE / "src" / "original" / "05-claim-support", stage_05)
    (stage_04 / "results").mkdir(exist_ok=True)
    copy_file(PACKAGE / "reference_settings" / "analysis_settings.json", stage_04 / "results" / "analysis_settings.json")
    for subdir in ["external", "inputs", "results", "r-library"]:
        (stage_05 / subdir).mkdir(exist_ok=True)
    external_root = Path(config["external_root"]).resolve()
    input_root = Path(config["processed_input_root"]).resolve()
    records = load_json(PACKAGE / "manifests" / "required_input_hashes.json")["kamath_local_existing"]
    checked = verify_hash_records(records, search_roots=[external_root, input_root])
    for record in records:
        basename = record["basename"]
        src = external_root / basename
        if not src.exists():
            src = input_root / basename
        dest_dir = stage_05 / ("external" if (external_root / basename).exists() else "inputs")
        copy_file(src, dest_dir / basename)
    for name in INPUT_NAMES:
        src = input_root / name
        if src.exists():
            copy_file(src, stage_05 / "inputs" / name)
    copy_file(PACKAGE / "reference_settings" / "foundin_regrouping_settings.json", stage_05 / "results" / "foundin_regrouping_settings.json")
    copy_file(
        PACKAGE / "results_reference" / "public_gse243639" / "gse243639_pilot_panel_tests.csv",
        stage_05 / "results" / "gse243639_pilot_panel_tests.csv",
    )
    return work_root, stage_root, stage_05, logs, checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PACKAGE / "config" / "kamath_local_existing.json"))
    parser.add_argument("--timeout", type=int, default=5400)
    args = parser.parse_args()
    config = load_json(args.config)
    work_root, stage_root, stage_05, logs, checked = stage_sources(config)
    commands = []
    for script in ["04_prepare_kamath_metadata.py", "13_test_kamath.py", "14_verify_kamath_results.py"]:
        commands.append(run_command(["python", script], cwd=stage_05, log_path=logs / f"{script}.log", timeout=args.timeout))
    commands.append(
        run_command(
            [
                "python",
                str(PACKAGE / "tools" / "compare_results.py"),
                "--profile",
                "kamath_edges1008",
                "--observed",
                str(stage_05 / "results" / "kamath_all_edges.csv"),
                "--reference",
                str(PACKAGE / "results_reference" / "aggregate_summaries" / "kamath_all_edges.csv"),
                "--out",
                str(work_root / "kamath_edges1008_comparison.json"),
            ],
            cwd=PACKAGE,
            log_path=logs / "compare_kamath_edges1008.log",
            timeout=600,
        )
    )
    commands.append(
        run_command(
            [
                "python",
                str(PACKAGE / "tools" / "compare_results.py"),
                "--profile",
                "kamath_panel48",
                "--observed",
                str(stage_05 / "results" / "kamath_panel_profile_tests.csv"),
                "--reference",
                str(PACKAGE / "results_reference" / "aggregate_summaries" / "kamath_panel_profile_tests.csv"),
                "--out",
                str(work_root / "kamath_panel48_comparison.json"),
            ],
            cwd=PACKAGE,
            log_path=logs / "compare_kamath_panel48.log",
            timeout=600,
        )
    )
    report = {
        "status": "PASS",
        "mode": "local_existing_public_inputs",
        "work_root": str(work_root),
        "stage_root": str(stage_root),
        "checked_inputs": checked,
        "commands": commands,
        "edge_comparison": json.loads((work_root / "kamath_edges1008_comparison.json").read_text()),
        "panel_comparison": json.loads((work_root / "kamath_panel48_comparison.json").read_text()),
        "note": "This mode stages public Kamath raw/processed inputs from local files, runs metadata preprocessing plus 13 and 14, and verifies all planned edge and panel/profile outputs.",
    }
    (work_root / "kamath_rerun_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
