#!/usr/bin/env python3
"""Download and rerun the public GSE243639 pilot analysis from raw public files."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import sys
import time
import urllib.request
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

PACKAGE = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGE / "src" / "repro_package"))
from numerics import partial_correlations
from staging import load_json, verify_hash_records
from stats_core import fisher, independent_columns, wild_test

SEED = 20260913
DRAWS = 19999
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def read_config(path):
    cfg = json.loads(Path(path).read_text())
    root = Path(path).resolve().parent
    for key in ["external_root", "work_root"]:
        cfg[key] = Path(cfg[key])
        if not cfg[key].is_absolute():
            cfg[key] = (root.parent / cfg[key]).resolve()
    return cfg


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def download_sources(cfg):
    cfg["external_root"].mkdir(parents=True, exist_ok=True)
    for name, url in cfg["urls"].items():
        dest = cfg["external_root"] / name
        if dest.exists():
            continue
        print(f"Downloading {name}", flush=True)
        with urllib.request.urlopen(url) as response, open(dest, "wb") as out:
            while True:
                block = response.read(1024 * 1024)
                if not block:
                    break
                out.write(block)


def xlsx_sheet(zipped, sheet, shared):
    root = ET.fromstring(zipped.read(f"xl/worksheets/sheet{sheet}.xml"))
    ns = {"m": root.tag.split("}")[0][1:]}
    records = []
    for row in root.findall("m:sheetData/m:row", ns):
        values = {}
        for cell in row.findall("m:c", ns):
            column = re.match("[A-Z]+", cell.attrib["r"]).group()
            value = cell.find("m:v", ns)
            if value is None:
                continue
            text = value.text
            values[column] = shared[int(text)] if cell.attrib.get("t") == "s" else text
        records.append(values)
    header = records[0]
    return pd.DataFrame(
        [{header[k]: v for k, v in row.items() if k in header} for row in records[1:]]
    )


def prepare_metadata(external_root, input_root):
    text = gzip.open(external_root / "GSE243639_Clinical_data.csv.gz", "rt").read()
    start = next(i for i, line in enumerate(text.splitlines()) if line.startswith("N;Brain Bank ID;"))
    clinical = pd.read_csv(io.StringIO("\n".join(text.splitlines()[start:])), sep=";")
    clinical.columns = clinical.columns.str.strip()
    clinical = clinical[clinical["Sample ID"].notna()].copy()
    clinical["donor"] = clinical["Sample ID"].str.strip()
    clinical["group"] = clinical["Clinical diagnosis"].str.strip().map(
        {"Parkinson's": "iPD", "Control": "HC"}
    )
    assert clinical.donor.is_unique and clinical.group.notna().all()
    clinical.to_csv(input_root / "gse243639_donors.csv", index=False)

    with ZipFile(external_root / "GSE243639_UMAP_coordinates.xlsx") as zipped:
        shared = ["".join(x.itertext()) for x in ET.fromstring(zipped.read("xl/sharedStrings.xml"))]
        broad = xlsx_sheet(zipped, 1, shared).rename(
            columns={"IDENT": "broad_type", "CLUSTER": "broad_cluster"}
        )
        neurons = xlsx_sheet(zipped, 4, shared).rename(
            columns={"IDENT": "neuron_type", "CLUSTER": "neuron_cluster"}
        )
    assert broad.CELL_ID.is_unique and neurons.CELL_ID.is_unique
    assert set(neurons.CELL_ID) <= set(broad.CELL_ID)
    broad = broad.merge(
        neurons[["CELL_ID", "neuron_type", "neuron_cluster"]],
        how="left",
        on="CELL_ID",
        validate="one_to_one",
    )
    broad["donor"] = broad.CELL_ID.str.split("_").str[0]
    broad = broad.merge(clinical, on="donor", how="left", validate="many_to_one")
    assert broad.group.notna().all()
    broad.to_csv(input_root / "gse243639_cell_metadata.csv.gz", index=False)
    broad.groupby(["donor", "group", "broad_type"], observed=True).size().rename(
        "n_cells"
    ).reset_index().to_csv(input_root / "gse243639_broad_coverage.csv", index=False)
    return broad


def prepare_counts(external_root, input_root, meta):
    settings = json.loads((PACKAGE / "reference_settings" / "analysis_settings.json").read_text())
    markers = [
        "TH", "SLC6A3", "SLC18A2", "ALDH1A1", "KCNJ6", "NR4A2", "TPH1", "TPH2",
        "SLC6A4", "SNAP25", "GAD1", "GAD2", "MAP2", "DCX", "TUBB3", "RBFOX3",
        "SOX2", "MKI67", "GFAP", "AQP4", "MBP", "PLP1", "P2RY12",
    ]
    requested = list(dict.fromkeys(["SNCA"] + settings["original_panel_genes"] + markers))
    all_genes = []
    selected = {}
    path = external_root / "GSE243639_Filtered_count_table.csv.gz"
    t0 = time.time()
    with gzip.open(path, "rt") as handle:
        header = next(csv.reader([handle.readline()]))
        cells = [re.sub(r"\.1$", "-1", s) for s in header[1:]]
        assert len(cells) == len(set(cells)) == len(meta)
        assert set(cells) == set(meta.CELL_ID)
        total = np.zeros(len(cells), dtype=np.float64)
        mitochondrial = np.zeros_like(total)
        detected = np.zeros(len(cells), dtype=np.int64)
        for i, line in enumerate(handle, 1):
            label, values = line.split(",", 1)
            gene = next(csv.reader([label]))[0]
            counts = np.fromstring(values, sep=",", dtype=np.float64)
            assert len(counts) == len(cells) and np.isfinite(counts).all() and (counts >= 0).all()
            assert np.equal(counts, np.floor(counts)).all(), gene
            total += counts
            detected += counts > 0
            if gene.upper().startswith("MT-"):
                mitochondrial += counts
            if gene in requested:
                selected[gene] = counts.astype(np.float32)
            all_genes.append(gene)
            if i % 3000 == 0:
                print(f"Processed genes {i}; elapsed {time.time() - t0:.1f}s", flush=True)
    kept = [gene for gene in requested if gene in selected]
    pd.DataFrame({gene: selected[gene] for gene in kept}, index=pd.Index(cells, name="cell")).to_csv(
        input_root / "gse243639_panel_counts.tsv.gz", sep="\t", float_format="%.0f"
    )
    pd.DataFrame(
        {
            "cell": cells,
            "nCount_RNA": total,
            "nFeature_RNA": detected,
            "percent_mt": 100 * mitochondrial / total,
        }
    ).to_csv(input_root / "gse243639_qc.csv.gz", index=False)
    (input_root / "gse243639_all_genes.txt").write_text("\n".join(all_genes) + "\n")
    audit = {
        "public_accession": "GSE243639",
        "n_cells": len(cells),
        "n_genes": len(all_genes),
        "source_counts_are_nonnegative_integers": True,
        "cell_ids_exactly_match_published_metadata": True,
        "selected_genes": kept,
        "missing_requested_genes": [gene for gene in requested if gene not in kept],
        "no_disease_effects_tested_during_preparation": True,
        "sha256": {p.name: sha256(p) for p in external_root.glob("GSE243639*")},
    }
    (input_root / "gse243639_input_audit.json").write_text(json.dumps(audit, indent=2) + "\n")


def donor_design(md, model):
    columns = [
        np.ones(len(md)),
        (md.Sex.str.strip() == "female").astype(float).values,
        (md.Age.values - 75) / 10,
    ]
    names = ["intercept", "female", "age_per_10_years"]
    if model == "age_sex_PMI_RIN":
        columns.extend([(md["PMI hours"].values - 20) / 10, md["RIN measure"].values - 7])
        names.extend(["PMI_per_10_hours", "RIN"])
    x0, names = independent_columns(np.column_stack(columns), names)
    return np.column_stack([x0, md.group.eq("iPD").astype(float)]), names + ["iPD"]


def fit(values, md, context, method, model, test):
    rows, null_t, observed_t = [], [], []
    weights = np.random.default_rng(SEED).choice([-1.0, 1.0], size=(DRAWS, len(md)))
    masks = {}
    for gene in values.columns:
        masks.setdefault(tuple(np.isfinite(values[gene])), []).append(gene)
    for mask, genes in masks.items():
        mask = np.array(mask)
        m = md.loc[mask]
        y = values.loc[mask, genes].to_numpy(float)
        info = dict(
            context=context,
            method=method,
            model=model,
            test=test,
            n_HC=int(m.group.eq("HC").sum()),
            n_iPD=int(m.group.eq("iPD").sum()),
        )
        try:
            if min(info["n_HC"], info["n_iPD"]) < 5:
                raise ValueError("Fewer than five independent donors in an arm")
            x, terms = donor_design(m, model)
            info["design"] = " + ".join(terms)
            f = wild_test(y, x, weights[:, mask])
            for j, gene in enumerate(genes):
                rows.append(
                    dict(
                        info,
                        gene=gene,
                        status="tested",
                        effect_z=f["effect"][j],
                        ci_low=f["ci_low"][j],
                        ci_high=f["ci_high"][j],
                        se_hc3=f["se"][j],
                        t_hc3=f["t"][j],
                        df=f["df"],
                        max_leverage=max(f["leverage"]),
                        p_wild=f["p_wild"][j],
                        p_hc3_t=f["p_hc3_t"][j],
                    )
                )
                observed_t.append(f["t"][j])
                null_t.append(f["tstar"][:, j])
        except ValueError as exc:
            for gene in genes:
                rows.append(dict(info, gene=gene, status=str(exc), p_wild=np.nan))
    panel = dict(context=context, method=method, model=model, test=test, genes_tested=len(observed_t))
    if null_t:
        obs = np.mean(np.square(observed_t))
        null = np.mean(np.column_stack(null_t) ** 2, axis=1)
        panel.update(statistic_mean_t2=obs, p_omnibus=(1 + (null >= obs).sum()) / (DRAWS + 1))
    else:
        panel["p_omnibus"] = np.nan
    return pd.DataFrame(rows), panel


def run_pilot(input_root, result_root):
    settings = json.loads((PACKAGE / "reference_settings" / "foundin_regrouping_settings.json").read_text())
    signs = settings["frozen_direction_signs"]
    meta = pd.read_csv(input_root / "gse243639_cell_metadata.csv.gz").set_index("CELL_ID")
    qc = pd.read_csv(input_root / "gse243639_qc.csv.gz", index_col="cell")
    meta = meta.join(qc)
    counts = pd.read_csv(input_root / "gse243639_panel_counts.tsv.gz", sep="\t", index_col="cell")
    neurons = meta[meta.neuron_type.notna()].copy()
    counts = counts.loc[neurons.index]
    expr = np.log1p(counts.div(neurons.nCount_RNA, axis=0) * 10000)
    coverage = neurons.assign(
        context=np.where(neurons.neuron_type.eq("neurons00"), "DA", "other_neurons")
    ).groupby(["donor", "group", "context"]).size().rename("n_cells").reset_index()
    coverage.to_csv(result_root / "gse243639_pilot_coverage.csv", index=False)
    primary = coverage[(coverage.context == "DA") & (coverage.n_cells >= 50)]
    primary_n = primary.groupby("group").donor.nunique().to_dict()
    assert primary_n == {"HC": 7, "iPD": 2}, primary_n

    eligible_donors = coverage[(coverage.context == "DA") & (coverage.n_cells >= 20)].donor
    da = neurons[neurons.neuron_type.eq("neurons00") & neurons.donor.isin(eligible_donors)]
    detections = counts.loc[da.index, settings["genes"]].gt(0).groupby(da.donor).mean().mean()
    eligibility = detections.rename("mean_donor_detection_DA").reset_index().rename(columns={"index": "gene"})
    eligibility["eligible"] = eligibility.mean_donor_detection_DA >= 0.05
    eligibility.to_csv(result_root / "gse243639_pilot_gene_eligibility.csv", index=False)
    genes = eligibility.loc[eligibility.eligible, "gene"].tolist()
    profile_genes = [gene for gene in signs if gene in genes]

    markers = ["TH", "SLC6A3", "SLC18A2", "ALDH1A1", "GAD1", "GAD2", "SNCA"]
    counts[markers].gt(0).groupby(neurons.neuron_type).mean().to_csv(
        result_root / "gse243639_neuron_marker_detection.csv"
    )
    all_corr = []
    for context, selected in [
        ("DA", neurons.neuron_type.eq("neurons00")),
        ("other_neurons", ~neurons.neuron_type.eq("neurons00")),
    ]:
        for donor, donor_meta in neurons[selected].groupby("donor"):
            n = len(donor_meta)
            if n < 20:
                continue
            a = expr.loc[donor_meta.index, ["SNCA"] + genes].to_numpy(float)
            cov = np.column_stack(
                [
                    np.ones(n),
                    np.log1p(donor_meta.nCount_RNA),
                    donor_meta.percent_mt / 100,
                    pd.get_dummies(donor_meta.neuron_type, drop_first=True).values,
                ]
            )
            ranks = np.column_stack([stats.rankdata(a[:, j]) for j in range(a.shape[1])])
            rcov = cov.copy()
            rcov[:, 1] = stats.rankdata(cov[:, 1])
            rcov[:, 2] = stats.rankdata(cov[:, 2])
            for method, matrix, nuisance in [("rna_qc", a, cov), ("spearman_qc", ranks, rcov)]:
                result = partial_correlations(matrix, nuisance)
                for j, gene in enumerate(genes, 1):
                    all_corr.append(
                        dict(
                            context=context,
                            donor=donor,
                            group=donor_meta.iloc[0].group,
                            method=method,
                            gene=gene,
                            n_cells=n,
                            r=result[j],
                            z=fisher(result[j]),
                        )
                    )
    corr = pd.DataFrame(all_corr)
    corr.to_csv(result_root / "gse243639_pilot_donor_correlations.csv", index=False)
    clinical = pd.read_csv(input_root / "gse243639_donors.csv").set_index("donor")
    edges, tests, profiles, donor_scores = [], [], [], []
    for method in ["rna_qc", "spearman_qc"]:
        da_values = corr[(corr.context == "DA") & (corr.method == method)].pivot(
            index="donor", columns="gene", values="z"
        )
        other_values = corr[(corr.context == "other_neurons") & (corr.method == method)].pivot(
            index="donor", columns="gene", values="z"
        )
        common = da_values.index.intersection(other_values.index)
        for context, values in [
            ("DA", da_values),
            ("DA_minus_other_neurons", da_values.loc[common] - other_values.loc[common]),
        ]:
            md = clinical.loc[values.index]
            score = values[profile_genes].mul(pd.Series(signs).reindex(profile_genes), axis=1).mean(
                axis=1, skipna=False
            )
            if len(profile_genes) < 12:
                score[:] = np.nan
            named = pd.DataFrame({"manuscript_direction_score": score})
            ds = md[["group", "Age", "Sex", "PMI hours", "RIN measure"]].copy()
            ds["score"] = score
            ds["context"], ds["method"] = context, method
            donor_scores.append(ds.reset_index())
            for model in ["age_sex", "age_sex_PMI_RIN"]:
                edge_result, panel = fit(values, md, context, method, model, "full_panel")
                edges.append(edge_result)
                tests.append(panel)
                profile_result, panel = fit(named, md, context, method, model, "frozen_direction_profile")
                profile_result["profile_genes"] = len(profile_genes)
                profiles.append(profile_result)
                tests.append(panel)
    edges = pd.concat(edges, ignore_index=True)
    adjusted = multipletests(edges.p_wild.fillna(1), method="fdr_bh")[1]
    edges["q_joint_pilot_family"] = np.where(edges.p_wild.notna(), adjusted, np.nan)
    edges.to_csv(result_root / "gse243639_pilot_all_edges.csv", index=False)
    pd.concat(profiles, ignore_index=True).to_csv(result_root / "gse243639_pilot_profile_effects.csv", index=False)
    pd.concat(donor_scores, ignore_index=True).to_csv(result_root / "gse243639_pilot_donor_scores.csv", index=False)
    tests = pd.DataFrame(tests)
    tests["p_holm_all_pilot_tests"] = multipletests(tests.p_omnibus.fillna(1), method="holm")[1]
    tests.to_csv(result_root / "gse243639_pilot_panel_tests.csv", index=False)
    (result_root / "gse243639_pilot_settings.json").write_text(
        json.dumps(
            {
                "primary_minimum_nuclei": 50,
                "primary_donor_counts": primary_n,
                "primary_replication_gate": "FAIL",
                "pilot_minimum_nuclei": 20,
                "eligible_genes": genes,
                "profile_genes": profile_genes,
                "profile_signs": {gene: signs[gene] for gene in profile_genes},
                "seed": SEED,
                "draws": DRAWS,
                "interpretation": "Exploratory low-coverage pilot; insufficient to establish independent replication alone",
                "normalization": "log(1 + 10000 * raw gene count / total count from supplied matrix)",
            },
            indent=2,
        )
        + "\n"
    )
    print(tests.to_string(index=False), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PACKAGE / "config" / "public_gse243639.json"))
    parser.add_argument("--external-root")
    parser.add_argument("--work-root")
    parser.add_argument("--skip-download", action="store_true")
    args = parser.parse_args()
    cfg = read_config(args.config)
    if args.external_root:
        cfg["external_root"] = Path(args.external_root).resolve()
    if args.work_root:
        cfg["work_root"] = Path(args.work_root).resolve()
    if args.skip_download:
        cfg["download"] = False
    input_root = cfg["work_root"] / "inputs"
    result_root = cfg["work_root"] / "results"
    input_root.mkdir(parents=True, exist_ok=True)
    result_root.mkdir(parents=True, exist_ok=True)
    if cfg.get("download", True):
        download_sources(cfg)
    required = load_json(PACKAGE / "manifests" / "required_input_hashes.json")["public_gse243639"]
    verify_hash_records(required, search_roots=[cfg["external_root"]])
    meta = prepare_metadata(cfg["external_root"], input_root)
    prepare_counts(cfg["external_root"], input_root, meta)
    run_pilot(input_root, result_root)


if __name__ == "__main__":
    main()
