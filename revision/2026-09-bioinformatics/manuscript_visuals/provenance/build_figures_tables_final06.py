#!/usr/bin/env python3
"""Build revised bioinformatics figures and tables for manuscript integration.

All statistics are read from the verified 04/05 analysis outputs. This script
formats source-backed display files only; it does not fit new statistical
models or change the upstream result sets.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import shutil
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "06-manuscript-revision"
FIG_DIR = OUT / "figures"
TAB_DIR = OUT / "tables"
FIG_SOURCE = FIG_DIR / "source_data"
TAB_SOURCE = TAB_DIR / "source_data"

PRIMARY = ROOT / "04-reviewer-reanalysis/results/primary_86_edges_for_reporting.csv"
COVERAGE = ROOT / "04-reviewer-reanalysis/results/culture_context_coverage.csv"
PSEUDOBULK_SUMMARY = ROOT / "04-reviewer-reanalysis/results/pseudobulk_DE_summary.csv"
PSEUDOBULK_ALL = ROOT / "04-reviewer-reanalysis/results/pseudobulk_DE_all_results.csv.gz"
ROAST = ROOT / "04-reviewer-reanalysis/results/autophagy_panel_roast.csv"
CAMERA = ROOT / "04-reviewer-reanalysis/results/camera_all_pathways.csv.gz"
FOUNDIN_PROFILE = ROOT / "05-claim-support/results/foundin_regrouped_profile_effects.csv"
FOUNDIN_PANEL = ROOT / "05-claim-support/results/foundin_regrouped_panel_tests.csv"
KAMATH_PROFILE = ROOT / "05-claim-support/results/kamath_profile_effects.csv"
KAMATH_PANEL = ROOT / "05-claim-support/results/kamath_panel_profile_tests.csv"
KAMATH_COVERAGE = ROOT / "05-claim-support/results/kamath_cell_coverage.csv"

CONTEXT_ORDER = [
    "iDA_pooled",
    "iDA1",
    "iDA2",
    "iDA3",
    "iDA4",
    "eProg1",
    "eProg2",
    "lProg1",
    "lProg2",
    "PFPP",
    "NE",
    "Ependymal",
]

CONTEXT_LABELS = {
    "iDA_pooled": "Pooled iDA1-4",
    "iDA_specificity": "Direct specificity",
    "iDA1": "iDA1",
    "iDA2": "iDA2",
    "iDA3": "iDA3",
    "iDA4": "iDA4",
    "eProg1": "eProg1",
    "eProg2": "eProg2",
    "lProg1": "lProg1",
    "lProg2": "lProg2",
    "PFPP": "PFPP",
    "NE": "NE",
    "Ependymal": "Ependymal",
}

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7.2,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "legend.frameon": False,
    }
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def ensure_dirs() -> None:
    for path in [FIG_DIR, TAB_DIR, FIG_SOURCE, TAB_SOURCE]:
        path.mkdir(parents=True, exist_ok=True)


def write_csv(df: pd.DataFrame, path: Path) -> None:
    df.to_csv(path, index=False)


def write_gzip_csv(df: pd.DataFrame, path: Path) -> None:
    with gzip.open(path, "wt", newline="") as f:
        df.to_csv(f, index=False)


def fmt_p(x: float | int | None) -> str:
    if pd.isna(x):
        return "NA"
    x = float(x)
    if x < 0.0001:
        return f"{x:.2e}"
    if x < 0.01:
        return f"{x:.4f}"
    return f"{x:.3f}"


def fmt_est_ci(est: float, low: float, high: float) -> str:
    if pd.isna(est):
        return "NA"
    return f"{est:.4f} ({low:.4f}, {high:.4f})"


def save_figure(fig: mpl.figure.Figure, stem: Path) -> list[str]:
    return save_figure_with_options(fig, stem, tight=True)


def save_figure_with_options(fig: mpl.figure.Figure, stem: Path, tight: bool) -> list[str]:
    outputs = []
    save_kwargs = {"bbox_inches": "tight"} if tight else {}
    for ext, kwargs in [
        (".svg", {}),
        (".pdf", {}),
        (".png", {"dpi": 600}),
    ]:
        path = stem.with_suffix(ext)
        fig.savefig(path, **save_kwargs, **kwargs)
        outputs.append(str(path.relative_to(ROOT)))
    png_path = stem.with_suffix(".png")
    tiff_path = stem.with_suffix(".tiff")
    with Image.open(png_path) as image:
        image.save(tiff_path, dpi=(600, 600))
    outputs.append(str(tiff_path.relative_to(ROOT)))
    return outputs


def load_primary() -> pd.DataFrame:
    df = pd.read_csv(PRIMARY)
    expected_contexts = {"iDA_pooled", "iDA_specificity"}
    observed_contexts = set(df["context"])
    if observed_contexts != expected_contexts:
        raise ValueError(f"Primary contexts mismatch: {observed_contexts}")
    if len(df) != 86:
        raise ValueError(f"Primary estimate count must be 86, found {len(df)}")
    by_context = df.groupby("context")["gene"].nunique().to_dict()
    if by_context != {"iDA_pooled": 43, "iDA_specificity": 43}:
        raise ValueError(f"Primary gene counts mismatch: {by_context}")
    if not (df["status"] == "tested").all():
        raise ValueError("Primary estimate table contains untested rows")
    return df.sort_values(["gene", "context"]).reset_index(drop=True)


def build_table4(primary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for gene in sorted(primary["gene"].unique()):
        out = {"gene": gene}
        for context, prefix in [
            ("iDA_pooled", "pooled_iDA1_4"),
            ("iDA_specificity", "direct_specificity"),
        ]:
            row = primary[(primary["gene"] == gene) & (primary["context"] == context)].iloc[0]
            out.update(
                {
                    f"{prefix}_n_HC": int(row["n_HC"]),
                    f"{prefix}_n_iPD": int(row["n_iPD"]),
                    f"{prefix}_disease_coefficient_fisher_z": row["effect_z"],
                    f"{prefix}_HC3_95CI_low": row["ci_low"],
                    f"{prefix}_HC3_95CI_high": row["ci_high"],
                    f"{prefix}_raw_wild_P": row["p_wild"],
                    f"{prefix}_BH_q_joint86": row["q_wild_family"],
                    f"{prefix}_mean_r_HC": row["mean_r_HC"],
                    f"{prefix}_mean_r_iPD": row["mean_r_iPD"],
                    f"{prefix}_display": fmt_est_ci(row["effect_z"], row["ci_low"], row["ci_high"]),
                    f"{prefix}_raw_wild_P_display": fmt_p(row["p_wild"]),
                    f"{prefix}_BH_q_display": fmt_p(row["q_wild_family"]),
                }
            )
        rows.append(out)
    table = pd.DataFrame(rows)
    if table.shape[0] != 43:
        raise ValueError(f"Table 4 must have 43 gene rows, found {table.shape[0]}")
    write_csv(primary, TAB_SOURCE / "Table4_source_primary_86_long.csv")
    write_csv(table, TAB_DIR / "Table4_revised_primary_SNCA_autophagy_estimates.csv")
    return table


def write_markdown_table(path: Path, title: str, df: pd.DataFrame, footnotes: list[str]) -> None:
    lines = [f"# {title}", "", df.to_markdown(index=False), ""]
    if footnotes:
        lines.append("## Footnotes")
        lines.append("")
        for note in footnotes:
            lines.append(f"- {note}")
    path.write_text("\n".join(lines) + "\n")


def build_table2() -> pd.DataFrame:
    de = pd.read_csv(PSEUDOBULK_SUMMARY)
    de_primary = de[de["scenario"] == "primary"].copy()
    rows = []
    for context in CONTEXT_ORDER:
        d = de_primary[de_primary["context"] == context]
        if d.empty:
            raise ValueError(f"Missing primary pseudobulk DE summary for {context}")
        d = d.iloc[0]
        rows.append(
            {
                "context": context,
                "display_context": CONTEXT_LABELS[context],
                "n_HC_donors": int(d["n_HC"]),
                "n_iPD_donors": int(d["n_iPD"]),
                "genes_tested_DE": int(d["genes_tested"]),
                "DEG_FDR05": int(d["DEG_FDR05"]),
            }
        )
    table = pd.DataFrame(rows)
    if table.shape[0] != 12:
        raise ValueError(f"Table 2 must cover 11 cell types plus pooled iDA, found {table.shape[0]}")
    write_csv(table, TAB_DIR / "Table2_revised_pseudobulk_DE_scope_summary.csv")
    md = table[
        [
            "display_context",
            "n_HC_donors",
            "n_iPD_donors",
            "genes_tested_DE",
            "DEG_FDR05",
        ]
    ].copy()
    md.columns = [
        "Context",
        "HC donors",
        "iPD donors",
        "Genes tested",
        "DEG at FDR < 0.05",
    ]
    write_markdown_table(
        TAB_DIR / "Table2_manuscript.md",
        "Table 2. Scope of primary pseudobulk differential-expression analyses",
        md,
        [
            "Rows cover pooled iDA1-4 plus the 11 individual annotated cell-type contexts.",
            "Donor counts are taken directly from the strict primary pseudobulk differential-expression model summaries.",
            "DEG counts use the context-level limma-voom FDR threshold of 0.05.",
            "ROAST and CAMERA pathway statistics are reported in Table 3.",
            "The full tested-gene companion table is Table2_companion_full_primary_DE_results.csv.gz.",
        ],
    )

    all_de = pd.read_csv(PSEUDOBULK_ALL)
    all_de_primary = all_de[all_de["scenario"] == "primary"].copy()
    write_gzip_csv(all_de_primary, TAB_SOURCE / "Table2_companion_full_primary_DE_results.csv.gz")
    return table


def build_table3() -> pd.DataFrame:
    roast = pd.read_csv(ROAST)
    roast_primary = roast[roast["scenario"] == "primary"].copy()
    camera = pd.read_csv(CAMERA)
    rows = []
    for context in CONTEXT_ORDER:
        r = roast_primary[roast_primary["context"] == context]
        if r.empty:
            raise ValueError(f"Missing ROAST row for {context}")
        r = r.iloc[0]
        c = camera[camera["context"] == context].copy()
        if c.empty:
            raise ValueError(f"Missing CAMERA rows for {context}")
        best = c.sort_values(["q_global_contexts", "PValue", "term"]).iloc[0]
        rows.append(
            {
                "context": context,
                "display_context": CONTEXT_LABELS[context],
                "n_HC": int(r["n_HC"]),
                "n_iPD": int(r["n_iPD"]),
                "ROAST_autophagy_panel_NGenes": int(r["NGenes"]),
                "ROAST_directional_P": r["PValue"],
                "ROAST_directional_FDR": r["FDR"],
                "ROAST_directional_q_across_primary_contexts": r[
                    "q_across_primary_contexts_directional"
                ],
                "ROAST_mixed_P": r["PValue.Mixed"],
                "ROAST_mixed_FDR": r["FDR.Mixed"],
                "ROAST_mixed_q_across_primary_contexts": r["q_across_primary_contexts_mixed"],
                "CAMERA_pathways_tested": int(c["term"].nunique()),
                "CAMERA_min_raw_P": float(c["PValue"].min()),
                "CAMERA_min_q_within_context": float(c["q_within_context"].min()),
                "CAMERA_min_q_global_contexts": float(c["q_global_contexts"].min()),
                "CAMERA_n_q_global_lt_0_05": int((c["q_global_contexts"] < 0.05).sum()),
                "CAMERA_best_global_q_term": best["term"],
            }
        )
    table = pd.DataFrame(rows)
    write_csv(table, TAB_DIR / "Table3_revised_measured_expression_pathway_tests.csv")
    md = table[
        [
            "display_context",
            "ROAST_autophagy_panel_NGenes",
            "ROAST_directional_P",
            "ROAST_directional_q_across_primary_contexts",
            "ROAST_mixed_P",
            "ROAST_mixed_q_across_primary_contexts",
            "CAMERA_pathways_tested",
            "CAMERA_min_q_global_contexts",
            "CAMERA_n_q_global_lt_0_05",
        ]
    ].copy()
    md.columns = [
        "Context",
        "ROAST genes",
        "ROAST directional P",
        "ROAST directional q",
        "ROAST mixed P",
        "ROAST mixed q",
        "CAMERA pathways tested",
        "Minimum CAMERA global q",
        "CAMERA pathways q < 0.05",
    ]
    write_markdown_table(
        TAB_DIR / "Table3_manuscript.md",
        "Table 3. Measured expression pathway tests",
        md,
        [
            "ROAST values test the measured autophagy-panel expression signature in each context.",
            "CAMERA values summarize the measured-expression pathway-family screen; pathway names are not interpreted as activation states.",
            "Across-context q values are reported for the primary-context family.",
        ],
    )
    camera_summary = (
        camera.groupby("context")
        .agg(
            CAMERA_pathways_tested=("term", "nunique"),
            CAMERA_min_raw_P=("PValue", "min"),
            CAMERA_min_q_within_context=("q_within_context", "min"),
            CAMERA_min_q_global_contexts=("q_global_contexts", "min"),
            CAMERA_n_q_global_lt_0_05=("q_global_contexts", lambda s: int((s < 0.05).sum())),
        )
        .reset_index()
    )
    write_csv(camera_summary, TAB_SOURCE / "Table3_CAMERA_all_pathway_summary_by_context.csv")
    return table


def build_fig2_source(primary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    fig2_source = primary.copy()
    fig2_source["display_context"] = fig2_source["context"].map(CONTEXT_LABELS)
    fig2_source["effect_definition"] = np.where(
        fig2_source["context"] == "iDA_pooled",
        "Adjusted iPD-minus-HC difference in within-culture SNCA-gene Fisher z for pooled iDA1-4.",
        "Adjusted iPD-minus-HC difference in the donor-level target-minus-comparison Fisher z contrast.",
    )
    write_csv(fig2_source, FIG_SOURCE / "Figure2_AB_primary_86_estimates.csv")

    foundin_prof = pd.read_csv(FOUNDIN_PROFILE)
    kamath_prof = pd.read_csv(KAMATH_PROFILE)
    foundin_panel = pd.read_csv(FOUNDIN_PANEL)
    kamath_panel = pd.read_csv(KAMATH_PANEL)

    f = foundin_prof[
        (foundin_prof["context"] == "DA123")
        & (foundin_prof["method"] == "rna_qc")
        & (foundin_prof["gene"] == "manuscript_direction_score")
    ].iloc[0]
    fp = foundin_panel[
        (foundin_panel["context"] == "DA123")
        & (foundin_panel["method"] == "rna_qc")
        & (foundin_panel["test"] == "frozen_direction_profile")
    ].iloc[0]
    k = kamath_prof[
        (kamath_prof["scenario"] == "published_primary")
        & (kamath_prof["context"] == "DA")
        & (kamath_prof["method"] == "rna_qc_library")
        & (kamath_prof["test"] == "profile")
    ].iloc[0]
    kp = kamath_panel[
        (kamath_panel["scenario"] == "published_primary")
        & (kamath_panel["context"] == "DA")
        & (kamath_panel["method"] == "rna_qc_library")
        & (kamath_panel["test"] == "profile")
    ].iloc[0]
    profile = pd.DataFrame(
        [
            {
                "dataset": "FOUNDIN DA123 sensitivity",
                "display_dataset": "FOUNDIN DA123",
                "effect_z": f["effect_z"],
                "ci_low": f["ci_low"],
                "ci_high": f["ci_high"],
                "n_HC": int(f["n_HC"]),
                "n_PD_or_iPD": int(f["n_iPD"]),
                "profile_genes": int(f["profile_genes"]),
                "raw_profile_P": f["p_wild"],
                "Holm_P": fp["p_holm_new_12"],
                "status": f["status"],
            },
            {
                "dataset": "Kamath primary DA",
                "display_dataset": "Kamath DA",
                "effect_z": k["effect_z"],
                "ci_low": k["ci_low"],
                "ci_high": k["ci_high"],
                "n_HC": int(k["n_HC"]),
                "n_PD_or_iPD": int(k["n_PD"]),
                "profile_genes": int(kp["genes_tested"]),
                "raw_profile_P": k["p_wild"],
                "Holm_P": kp["p_Holm_primary_4"],
                "status": k["status"],
            },
        ]
    )
    write_csv(profile, FIG_SOURCE / "Figure2_D_fixed_direction_profile_source.csv")

    kamath_cov = pd.read_csv(KAMATH_COVERAGE)
    primary_n = primary[primary["context"] == "iDA_pooled"][["n_HC", "n_iPD"]].drop_duplicates()
    if primary_n.shape[0] != 1:
        raise ValueError("Primary Figure 2 n_HC/n_iPD are not unique")
    coverage_rows = [
        {
            "dataset": "FOUNDIN pooled iDA1-4 primary",
            "HC_donors": int(primary_n.iloc[0]["n_HC"]),
            "PD_or_iPD_donors": int(primary_n.iloc[0]["n_iPD"]),
            "context": "pooled iDA1-4",
        },
        {
            "dataset": "FOUNDIN DA123 sensitivity",
            "HC_donors": int(f["n_HC"]),
            "PD_or_iPD_donors": int(f["n_iPD"]),
            "context": "DA123",
        },
        {
            "dataset": "Kamath primary DA",
            "HC_donors": int(
                kamath_cov[
                    (kamath_cov["clinical_group"] == "HC")
                    & (kamath_cov["clusters"] == "DA")
                    & (kamath_cov["passes_50_nuclei"])
                ]["donor"].nunique()
            ),
            "PD_or_iPD_donors": int(
                kamath_cov[
                    (kamath_cov["clinical_group"] == "PD")
                    & (kamath_cov["clusters"] == "DA")
                    & (kamath_cov["passes_50_nuclei"])
                ]["donor"].nunique()
            ),
            "context": "author DA nuclei",
        },
    ]
    coverage = pd.DataFrame(coverage_rows)
    write_csv(coverage, FIG_SOURCE / "Figure2_C_donor_cell_coverage_source.csv")
    return profile, coverage


def plot_forest_axis(ax, data: pd.DataFrame, title: str, color: str, show_y: bool) -> None:
    genes = sorted(data["gene"].unique())
    pos = np.arange(len(genes))
    data = data.set_index("gene").loc[genes].reset_index()
    ax.axvline(0, color="#444444", lw=0.7, zorder=1)
    for y, row in zip(pos, data.itertuples(index=False)):
        ax.plot([row.ci_low, row.ci_high], [y, y], color="#6b7280", lw=0.65, zorder=2)
        ax.scatter(row.effect_z, y, s=11, color=color, edgecolor="white", linewidth=0.25, zorder=3)
    ax.set_ylim(-0.8, len(genes) - 0.2)
    ax.invert_yaxis()
    ax.set_title(title, loc="left", fontsize=7.8, fontweight="bold", pad=3)
    ax.set_xlabel("Disease coefficient in Fisher z (95% HC3 CI)", fontsize=6.8)
    ax.grid(axis="x", color="#e5e7eb", lw=0.4)
    ax.tick_params(axis="x", labelsize=6.2, length=2)
    ax.tick_params(axis="y", labelsize=6.2, length=0)
    if show_y:
        ax.set_yticks(pos)
        ax.set_yticklabels(genes)
    else:
        ax.set_yticks(pos)
        ax.set_yticklabels([])
    lo = min(data["ci_low"].min(), -0.075)
    hi = max(data["ci_high"].max(), 0.075)
    pad = (hi - lo) * 0.05
    ax.set_xlim(lo - pad, hi + pad)


def build_figure2(primary: pd.DataFrame, profile: pd.DataFrame, coverage: pd.DataFrame) -> list[str]:
    width = 180 / 25.4
    height = 238 / 25.4
    fig = plt.figure(figsize=(width, height), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[5.6, 1.45],
        width_ratios=[1, 1],
        left=0.13,
        right=0.965,
        top=0.948,
        bottom=0.145,
        hspace=0.44,
        wspace=0.28,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    plot_forest_axis(
        ax_a,
        primary[primary["context"] == "iDA_pooled"],
        "a  Pooled iDA1-4 primary context",
        "#2563eb",
        True,
    )
    plot_forest_axis(
        ax_b,
        primary[primary["context"] == "iDA_specificity"],
        "b  Direct specificity contrast",
        "#7c3aed",
        False,
    )

    coverage_plot = coverage.iloc[::-1].reset_index(drop=True)
    coverage_plot["display_dataset"] = coverage_plot["dataset"].replace(
        {
            "FOUNDIN pooled iDA1-4 primary": "FOUNDIN iDA1-4",
            "FOUNDIN DA123 sensitivity": "FOUNDIN iDA1-3",
            "Kamath primary DA": "Kamath DA",
        }
    )
    y = np.arange(len(coverage_plot))
    ax_c.barh(y - 0.17, coverage_plot["HC_donors"], height=0.28, color="#94a3b8", label="HC donors")
    ax_c.barh(
        y + 0.17,
        coverage_plot["PD_or_iPD_donors"],
        height=0.28,
        color="#fb923c",
        label="Disease donors",
    )
    for yi, row in zip(y, coverage_plot.itertuples(index=False)):
        ax_c.text(
            max(row.HC_donors - 0.35, 0.35),
            yi - 0.17,
            str(row.HC_donors),
            va="center",
            ha="right",
            fontsize=6.2,
            color="white",
        )
        ax_c.text(
            max(row.PD_or_iPD_donors - 0.35, 0.35),
            yi + 0.17,
            str(row.PD_or_iPD_donors),
            va="center",
            ha="right",
            fontsize=6.2,
            color="white",
        )
    ax_c.set_yticks(y)
    ax_c.set_yticklabels(coverage_plot["display_dataset"], fontsize=6.2)
    ax_c.set_xlabel("Donors passing analysis gate", fontsize=6.8)
    ax_c.set_title("c  Donor coverage", loc="left", fontsize=7.8, fontweight="bold", pad=3)
    ax_c.tick_params(axis="x", labelsize=6.2)
    ax_c.grid(axis="x", color="#e5e7eb", lw=0.4)
    ax_c.set_xlim(0, float(coverage_plot[["HC_donors", "PD_or_iPD_donors"]].max().max()) + 4)
    ax_c.legend(loc="lower right", fontsize=6.2, ncol=2, handlelength=1.2, columnspacing=0.8)

    prof_plot = profile.iloc[::-1].reset_index(drop=True)
    y2 = np.arange(len(prof_plot))
    ax_d.axvline(0, color="#444444", lw=0.7)
    for yi, row in zip(y2, prof_plot.itertuples(index=False)):
        ax_d.plot([row.ci_low, row.ci_high], [yi, yi], color="#6b7280", lw=0.8)
        ax_d.scatter(row.effect_z, yi, s=18, color="#059669", edgecolor="white", linewidth=0.3)
        label_y = yi - 0.16 if yi > 0 else yi + 0.16
        ax_d.text(
            0.018,
            label_y,
            f"P={fmt_p(row.raw_profile_P)}; Holm={fmt_p(row.Holm_P)}",
            va="center",
            ha="left",
            fontsize=6.2,
        )
    ax_d.set_yticks(y2)
    ax_d.set_yticklabels(prof_plot["display_dataset"], fontsize=6.2)
    ax_d.set_xlabel("Fixed-direction profile coefficient (95% HC3 CI)", fontsize=6.8)
    ax_d.set_title("d  Frozen 16-gene profile tests", loc="left", fontsize=7.8, fontweight="bold", pad=3)
    ax_d.tick_params(axis="x", labelsize=6.2)
    ax_d.grid(axis="x", color="#e5e7eb", lw=0.4)
    xmin = min(prof_plot["ci_low"].min() - 0.01, -0.075)
    xmax = max(prof_plot["ci_high"].max() + 0.035, 0.065)
    ax_d.set_xlim(xmin, xmax)

    fig.text(
        0.13,
        0.985,
        "Figure 2 | Donor-level SNCA-autophagy coupling estimates",
        fontsize=8.6,
        fontweight="bold",
        ha="left",
        va="top",
    )
    fig.text(
        0.13,
        0.022,
        "Intervals are pointwise HC3 95% CIs. Direct specificity is a target-minus-comparison donor contrast,\nnot a raw delta-r.",
        fontsize=6.2,
        ha="left",
    )
    outputs = save_figure_with_options(fig, FIG_DIR / "Figure_2_revised", tight=False)
    plt.close(fig)
    return outputs


def build_supplementary_figure2() -> tuple[list[str], pd.DataFrame]:
    de = pd.read_csv(PSEUDOBULK_ALL)
    primary = de[de["scenario"] == "primary"].copy()
    primary["minus_log10_p"] = -np.log10(primary["P.Value"].clip(lower=np.finfo(float).tiny))
    primary["significant_FDR05"] = primary["adj.P.Val"] < 0.05
    write_gzip_csv(primary, FIG_SOURCE / "Supplementary_Figure2_all_primary_pseudobulk_volcano_source.csv.gz")

    width = 180 / 25.4
    height = 214 / 25.4
    fig, axes = plt.subplots(4, 3, figsize=(width, height), sharex=True, sharey=True)
    axes = axes.flatten()
    finite_logfc = primary.loc[np.isfinite(primary["logFC"]), "logFC"]
    finite_mlogp = primary.loc[np.isfinite(primary["minus_log10_p"]), "minus_log10_p"]
    xlim = max(1.0, float(np.abs(finite_logfc).max()) * 1.04)
    ymax = max(2.0, float(finite_mlogp.max()) * 1.04)
    summary_rows = []
    for ax, context in zip(axes, CONTEXT_ORDER):
        d = primary[primary["context"] == context]
        sig = d[d["significant_FDR05"]]
        nonsig = d[~d["significant_FDR05"]]
        ax.scatter(
            nonsig["logFC"],
            nonsig["minus_log10_p"],
            s=1.2,
            color="#94a3b8",
            alpha=0.35,
            linewidths=0,
            rasterized=True,
        )
        if len(sig):
            ax.scatter(sig["logFC"], sig["minus_log10_p"], s=3, color="#dc2626", alpha=0.7, linewidths=0)
        ax.axvline(0, color="#4b5563", lw=0.4)
        ax.set_title(CONTEXT_LABELS[context], fontsize=6.7, pad=2)
        ax.tick_params(axis="both", labelsize=5.6, length=2)
        ax.grid(color="#eef2f7", lw=0.35)
        ax.set_xlim(-xlim, xlim)
        ax.set_ylim(0, ymax)
        summary_rows.append(
            {
                "context": context,
                "display_context": CONTEXT_LABELS[context],
                "genes_tested": int(d["gene"].nunique()),
                "DEG_FDR05": int(d["significant_FDR05"].sum()),
                "min_raw_P": float(d["P.Value"].min()),
                "min_adj_P": float(d["adj.P.Val"].min()),
            }
        )
    fig.text(0.5, 0.018, "log2 fold-change, iPD versus HC", ha="center", fontsize=7)
    fig.text(
        0.012,
        0.5,
        "-log10 raw P",
        va="center",
        rotation="vertical",
        rotation_mode="anchor",
        fontsize=7,
    )
    fig.suptitle("Supplementary Figure 2 | Pseudobulk differential-expression screens across primary contexts", fontsize=8, fontweight="bold", y=0.985)
    fig.tight_layout(rect=[0.04, 0.04, 1, 0.965])
    outputs = save_figure(fig, FIG_DIR / "Supplementary_Figure_2_revised")
    plt.close(fig)
    summary = pd.DataFrame(summary_rows)
    write_csv(summary, FIG_SOURCE / "Supplementary_Figure2_context_summary.csv")
    return outputs, summary


def build_captions() -> None:
    primary = pd.read_csv(FIG_SOURCE / "Figure2_AB_primary_86_estimates.csv")
    profile = pd.read_csv(FIG_SOURCE / "Figure2_D_fixed_direction_profile_source.csv")
    coverage = pd.read_csv(FIG_SOURCE / "Figure2_C_donor_cell_coverage_source.csv")
    supp = pd.read_csv(FIG_SOURCE / "Supplementary_Figure2_context_summary.csv")
    pooled_n = primary[primary["context"] == "iDA_pooled"][["n_HC", "n_iPD"]].drop_duplicates().iloc[0]
    spec_n = primary[primary["context"] == "iDA_specificity"][["n_HC", "n_iPD"]].drop_duplicates().iloc[0]
    foundin_prof = profile[profile["dataset"] == "FOUNDIN DA123 sensitivity"].iloc[0]
    kamath_prof = profile[profile["dataset"] == "Kamath primary DA"].iloc[0]
    coverage_text_parts = []
    for row in coverage.itertuples(index=False):
        disease_label = "PD" if "Kamath" in row.dataset else "iPD"
        coverage_text_parts.append(
            f"{row.dataset}: {int(row.HC_donors)} HC/{int(row.PD_or_iPD_donors)} {disease_label} donors"
        )
    coverage_text = "; ".join(coverage_text_parts)
    supp_context_text = "; ".join(
        f"{row.display_context}: {int(row.genes_tested)} genes, {int(row.DEG_FDR05)} FDR<0.05"
        for row in supp.itertuples(index=False)
    )
    captions_json = {
        "Figure2": {
            "title": "Donor-level SNCA-autophagy coupling estimates do not establish reproducible disease-specific rewiring",
            "panel_map": {
                "a": "All 43 eligible SNCA-autophagy gene correlations in the strict primary FOUNDIN pooled iDA1-4 donor analysis.",
                "b": "All 43 direct specificity coefficients from the same strict primary donor set; the coefficient is computed on the pooled iDA1-4 minus other-cells donor contrast.",
                "c": "Donor coverage for the primary and support analyses, using HC/iPD for FOUNDIN and HC/PD for Kamath.",
                "d": "Frozen 16-gene fixed-direction profile estimates for the FOUNDIN DA123 sensitivity and independent Kamath primary DA analyses.",
            },
            "caption": (
                "Figure 2 | Donor-level SNCA-autophagy coupling estimates do not establish reproducible "
                "disease-specific rewiring. a, Disease coefficients for all 43 eligible SNCA-autophagy "
                f"genes in the strict primary FOUNDIN pooled iDA1-4 context (n={int(pooled_n.n_HC)} HC and "
                f"n={int(pooled_n.n_iPD)} iPD donors). The model was intercept + female + batch-fraction "
                "covariates + iPD, fitted to donor-level within-culture SNCA-gene Fisher z values. Points "
                "show adjusted iPD-minus-HC coefficients and horizontal bars show pointwise HC3 95% confidence "
                "intervals. b, Direct specificity coefficients for the same 43 genes "
                f"(n={int(spec_n.n_HC)} HC and n={int(spec_n.n_iPD)} iPD donors), defined as the disease "
                "effect on the donor-level pooled iDA1-4 minus other-cells Fisher z contrast. These coefficients are "
                "not raw differences between pooled-cell correlations. Raw P values are wild-bootstrap P values; "
                "q values are Benjamini-Hochberg adjusted jointly across the 86 primary tests. c, Donor coverage "
                f"for the displayed analyses ({coverage_text}); panel c uses HC/iPD for FOUNDIN analyses and HC/PD "
                "for the Kamath analysis. d, Frozen 16-gene fixed-direction profile tests "
                f"for FOUNDIN DA123 (n={int(foundin_prof.n_HC)} HC and n={int(foundin_prof.n_PD_or_iPD)} iPD donors; "
                f"16 genes; raw wild-bootstrap P={fmt_p(foundin_prof.raw_profile_P)}, "
                f"Holm-adjusted P={fmt_p(foundin_prof.Holm_P)} across the 12 FOUNDIN regrouped tests) and Kamath primary DA "
                f"(n={int(kamath_prof.n_HC)} HC and n={int(kamath_prof.n_PD_or_iPD)} PD donors; "
                "age-, sex- and brain-bank-adjusted model; 16 genes summarized as the same frozen profile; "
                f"raw wild-bootstrap P={fmt_p(kamath_prof.raw_profile_P)}, "
                f"Holm-adjusted P={fmt_p(kamath_prof.Holm_P)} across the 4 primary RNA Kamath tests)."
            ),
            "statistics": {
                "effect_unit": "donor-level Fisher z coefficient",
                "error_bars": "pointwise HC3 95% confidence intervals",
                "primary_multiplicity": "Benjamini-Hochberg adjustment across 86 primary tests",
                "raw_p_values": "wild-bootstrap P values",
                "direct_specificity_definition": "disease effect on target-minus-comparison donor contrasts, not raw delta-r",
            },
            "source_data": [
                "figures/source_data/Figure2_AB_primary_86_estimates.csv",
                "figures/source_data/Figure2_C_donor_cell_coverage_source.csv",
                "figures/source_data/Figure2_D_fixed_direction_profile_source.csv",
            ],
        },
        "SupplementaryFigure2": {
            "title": "Pseudobulk differential-expression screens across primary contexts",
            "panel_map": {
                "all_panels": "Twelve small-multiple volcano plots for pooled iDA1-4 and 11 primary annotated cell-type contexts."
            },
            "caption": (
                "Supplementary Figure 2 | Pseudobulk differential-expression screens across primary contexts. "
                "Volcano plots show the existing sex- and batch-adjusted limma-voom primary pseudobulk "
                "differential-expression results for pooled iDA1-4 and 11 annotated cell-type contexts. "
                "Pseudobulk samples were aggregated at the culture/context level using the >=20 cells per "
                "culture/context gate; Table 2 reports the exact donor n for each context. Each point is one tested gene. The x axis "
                "shows log2 fold-change for iPD versus HC and the y axis shows -log10 raw P. Axis limits are set "
                "from the full finite data range, so points are not percentile-clipped. Context-level Benjamini-Hochberg "
                f"FDR was used for DEG calls. Summary by panel: {supp_context_text}."
            ),
            "statistics": {
                "model": "sex- and batch-adjusted limma-voom pseudobulk differential-expression primary models",
                "aggregation": "culture/context pseudobulk aggregation with >=20 cells per culture/context; exact donor n in Table 2",
                "x_axis": "log2 fold-change, iPD versus HC",
                "y_axis": "-log10 raw P",
                "multiplicity": "context-level Benjamini-Hochberg FDR for DEG calls",
                "axis_policy": "full finite data range with padding; no unreported percentile clipping",
            },
            "source_data": [
                "figures/source_data/Supplementary_Figure2_all_primary_pseudobulk_volcano_source.csv.gz",
                "figures/source_data/Supplementary_Figure2_context_summary.csv",
            ],
        },
    }
    (OUT / "figures" / "captions.json").write_text(json.dumps(captions_json, indent=2))
    text = """# Revised Bioinformatics Figure and Table Captions

## Figure 2 | Donor-level SNCA-autophagy coupling estimates do not establish reproducible disease-specific rewiring

**a,** Adjusted iPD-minus-HC disease coefficients for all 43 eligible SNCA-autophagy gene correlations in the strict primary pooled iDA1-4 context (n=8 HC and n=28 iPD donors). Points show the donor-level coefficient in Fisher z units and horizontal bars show pointwise HC3 95% confidence intervals. **b,** Direct specificity coefficients for the same 43 genes and strict n=8/28 donor set, computed as the adjusted disease effect on the donor-level pooled iDA1-4 minus other-cells contrast. This panel is not a raw difference between pooled-cell correlations. Raw P values are wild-bootstrap P values and q values are Benjamini-Hochberg adjusted jointly across the 86 primary tests. **c,** Donor coverage for the primary FOUNDIN-PD analysis (HC/iPD), the DA123 sensitivity analysis (HC/iPD), and the independent Kamath DA analysis (HC/PD). **d,** Frozen 16-gene fixed-direction profile tests in the biologically motivated FOUNDIN DA123 sensitivity analysis and in the independent Kamath primary DA analysis. Holm families are the 12 FOUNDIN regrouped tests and the 4 primary RNA Kamath tests, respectively.

## Supplementary Figure 2 | Pseudobulk differential-expression screens across primary contexts

Volcano plots show existing sex- and batch-adjusted limma-voom pseudobulk differential-expression results for all 12 primary contexts. Pseudobulk samples were aggregated at the culture/context level using the >=20 cells per culture/context gate; Table 2 gives exact donor n for each context. Points represent tested genes; red points would indicate genes passing context-level FDR < 0.05. Axis limits cover all finite points. No primary context contained FDR-significant differentially expressed genes, and these panels replace the earlier selected-volcano/correlation display.

## Table 2 | Scope of primary pseudobulk differential-expression analyses

Table 2 reports the strict primary DE donor n, tested-gene count and FDR-significant DEG count for pooled iDA1-4 and the 11 individual annotated cell-type contexts.

## Table 3 | Measured expression pathway tests

Table 3 reports measured-expression ROAST autophagy-panel tests and CAMERA pathway-family summaries. The table does not report pathway activation states from descriptive enrichment labels.

## Table 4 | Complete primary donor-level SNCA-autophagy estimates

Table 4 reports all 43 eligible genes in fixed alphabetical order, with disease coefficients, pointwise HC3 95% confidence intervals, raw wild-bootstrap P values, and Benjamini-Hochberg q values for the primary pooled iDA1-4 context and the direct specificity contrast.
"""
    (OUT / "figures" / "captions_and_panel_notes.md").write_text(text)


def validation_record(outputs: dict[str, object]) -> None:
    source_files = [
        PRIMARY,
        COVERAGE,
        PSEUDOBULK_SUMMARY,
        PSEUDOBULK_ALL,
        ROAST,
        CAMERA,
        FOUNDIN_PROFILE,
        FOUNDIN_PANEL,
        KAMATH_PROFILE,
        KAMATH_PANEL,
        KAMATH_COVERAGE,
    ]
    output_files = []
    for path in [FIG_DIR, TAB_DIR]:
        output_files.extend(
            sorted(
                p
                for p in path.rglob("*")
                if p.is_file() and p != (OUT / "figures" / "figure_table_validation.json")
            )
        )
    record = {
        "status": "built",
        "script": str((OUT / "build_figures_tables.py").relative_to(ROOT)),
        "figure_contract": {
            "backend": "python/matplotlib",
            "figure2_width_mm": 180,
            "figure2_height_mm": 238,
            "supplementary_figure2_width_mm": 180,
            "supplementary_figure2_height_mm": 214,
            "glyph_floor_pt": 6.2,
            "primary_estimates_required": 86,
            "primary_genes_required_per_context": 43,
            "interval_definition": "pointwise HC3 95% confidence intervals",
        },
        "checks": outputs,
        "source_files": [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in source_files],
        "output_files": [{"path": str(p.relative_to(ROOT)), "sha256": sha256(p)} for p in output_files],
    }
    (OUT / "figures" / "figure_table_validation.json").write_text(json.dumps(record, indent=2))


def main() -> None:
    ensure_dirs()
    primary = load_primary()
    table4 = build_table4(primary)
    table2 = build_table2()
    table3 = build_table3()
    profile, coverage = build_fig2_source(primary)
    fig2_outputs = build_figure2(primary, profile, coverage)
    supp_outputs, supp_summary = build_supplementary_figure2()
    build_captions()
    shutil.copyfile(FIG_DIR / "Figure_2_revised.png", FIG_DIR / "Figure2_revised.png")
    shutil.copyfile(FIG_DIR / "Supplementary_Figure_2_revised.png", FIG_DIR / "Supplementary_Figure2_revised.png")
    table4_md = table4[
        [
            "gene",
            "pooled_iDA1_4_display",
            "pooled_iDA1_4_raw_wild_P_display",
            "pooled_iDA1_4_BH_q_display",
            "direct_specificity_display",
            "direct_specificity_raw_wild_P_display",
            "direct_specificity_BH_q_display",
        ]
    ].copy()
    table4_md.columns = [
        "Gene",
        "Pooled coefficient (95% CI)",
        "Pooled raw P",
        "Pooled BH q",
        "Specificity coefficient (95% CI)",
        "Specificity raw P",
        "Specificity BH q",
    ]
    write_markdown_table(
        TAB_DIR / "Table4_manuscript.md",
        "Table 4. Complete primary donor-level SNCA-autophagy estimates",
        table4_md,
        [
            "All 43 eligible genes are shown in fixed alphabetical order.",
            "Both primary families used n = 8 HC and n = 28 iPD donors.",
            "Coefficients are adjusted iPD-minus-HC differences in donor-level Fisher z values.",
            "The specificity coefficient is the disease effect on the target-minus-comparison donor contrast and should not be interpreted as a raw delta-r.",
            "P values are wild-bootstrap P values; q values are Benjamini-Hochberg adjusted jointly across the 86 primary tests.",
            "The complete source CSV preserves repeated n, group mean correlations and model fields.",
        ],
    )
    outputs = {
        "table4_rows": int(table4.shape[0]),
        "table2_rows": int(table2.shape[0]),
        "table3_rows": int(table3.shape[0]),
        "table2_primary_DE_counts": table2[
            ["context", "n_HC_donors", "n_iPD_donors", "genes_tested_DE", "DEG_FDR05"]
        ].to_dict(orient="records"),
        "figure2_donor_coverage": coverage.to_dict(orient="records"),
        "fig2_outputs": fig2_outputs,
        "supplementary_figure2_outputs": supp_outputs,
        "supplementary_figure2_contexts": int(supp_summary.shape[0]),
        "primary_estimate_rows": int(primary.shape[0]),
        "primary_gene_count": int(primary["gene"].nunique()),
        "table2_full_DE_companion_rows": int(pd.read_csv(TAB_SOURCE / "Table2_companion_full_primary_DE_results.csv.gz", usecols=["gene"]).shape[0]),
    }
    validation_record(outputs)
    print(json.dumps(outputs, indent=2))


if __name__ == "__main__":
    main()
