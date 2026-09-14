#!/usr/bin/env python3
"""Stage and run the authorized FOUNDIN primary pipeline from verified cached inputs."""
import argparse
import csv
import json
import sys
from pathlib import Path


PACKAGE = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE / "src" / "repro_package"))
from staging import copy_file, copy_tree, ensure_clean_dir, load_json, patch_text, run_command, verify_hash_records


def write_gene_list(stage_04):
    settings = json.loads((PACKAGE / "reference_settings" / "analysis_settings.json").read_text())
    path = stage_04 / "inputs" / "original_panel_genes.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene"])
        writer.writeheader()
        for gene in settings["original_panel_genes"]:
            writer.writerow({"gene": gene})
    return path


def stage_sources(config):
    work_root = (PACKAGE / config["work_root"]).resolve()
    stage_root = ensure_clean_dir(work_root / "stage")
    logs = work_root / "logs"
    copy_tree(PACKAGE / "src" / "original" / "04-reviewer-reanalysis", stage_root / "04-reviewer-reanalysis")
    stage_04 = stage_root / "04-reviewer-reanalysis"
    (stage_04 / "inputs").mkdir(exist_ok=True)
    (stage_04 / "results").mkdir(exist_ok=True)
    input_root = Path(config["cached_04_inputs_root"]).resolve()
    records = load_json(PACKAGE / "manifests" / "required_input_hashes.json")["authorized_foundin_cached"]
    checked = verify_hash_records(records, search_roots=[input_root])
    for record in records:
        copy_file(input_root / record["basename"], stage_04 / "inputs" / record["basename"])
    write_gene_list(stage_04)
    legacy_gene_list = (
        "pd.read_csv('/data" + "32TB/shared/ppmi-foundin/"
        "analyse/RNA-processed/shaohua-gene-list-2-10-25.csv').gene.tolist()"
    )
    patch_text(
        stage_04 / "03_within_donor_coupling.py",
        {
            legacy_gene_list: "pd.read_csv(INPUT/'original_panel_genes.csv').gene.tolist()",
        },
    )
    return work_root, stage_root, stage_04, logs, checked


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PACKAGE / "config" / "authorized_foundin_cached.json"))
    parser.add_argument("--timeout", type=int, default=5400)
    args = parser.parse_args()
    config = load_json(args.config)
    work_root, stage_root, stage_04, logs, checked = stage_sources(config)
    commands = [
        run_command(["python", "03_within_donor_coupling.py"], cwd=stage_04, log_path=logs / "03_within_donor_coupling.log", timeout=args.timeout),
        run_command(
            [
                "python",
                str(PACKAGE / "tools" / "extract_foundin_primary86.py"),
                "--all-results",
                str(stage_04 / "results" / "within_donor_coupling_all_results.csv"),
                "--reference",
                str(PACKAGE / "results_reference" / "aggregate_summaries" / "primary_86_edges_for_reporting.csv"),
                "--out",
                str(stage_04 / "results" / "primary_86_edges_for_reporting.csv"),
            ],
            cwd=PACKAGE,
            log_path=logs / "extract_foundin_primary86.log",
            timeout=600,
        ),
        run_command(
            [
                "python",
                str(PACKAGE / "tools" / "compare_results.py"),
                "--profile",
                "foundin_primary86",
                "--observed",
                str(stage_04 / "results" / "primary_86_edges_for_reporting.csv"),
                "--reference",
                str(PACKAGE / "results_reference" / "aggregate_summaries" / "primary_86_edges_for_reporting.csv"),
                "--out",
                str(work_root / "foundin_primary86_comparison.json"),
            ],
            cwd=PACKAGE,
            log_path=logs / "compare_foundin_primary86.log",
            timeout=600,
        ),
    ]
    report = {
        "status": "PASS",
        "mode": "verified_cached_input_rerun",
        "work_root": str(work_root),
        "stage_root": str(stage_root),
        "checked_inputs": checked,
        "commands": commands,
        "comparison": json.loads((work_root / "foundin_primary86_comparison.json").read_text()),
        "note": "This mode reruns the current primary FOUNDIN donor-level pipeline from verified cached 04 inputs; it does not reload the original 46GB Seurat object.",
    }
    (work_root / "authorized_foundin_rerun_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
