#!/usr/bin/env python3
"""Export source-backed figures and descriptive correlation confidence intervals."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
FIGURES.mkdir(exist_ok=True)
TEAL, ORANGE, INK, GREY = "#17647A", "#BA5C35", "#263746", "#6A7882"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.edgecolor": "#AAB7BD", "pdf.fonttype": 42, "ps.fonttype": 42,
})


def save(fig, stem):
    fig.savefig(FIGURES / (stem + ".png"), dpi=210, facecolor="white")
    fig.savefig(FIGURES / (stem + ".pdf"), facecolor="white")
    plt.close(fig)


def main():
    all_edges = pd.read_csv(RESULTS / "within_donor_coupling_all_results.csv")
    contexts = ["iDA_pooled", "iDA_specificity"]
    primary = all_edges[all_edges.method.eq("rna_qc") & all_edges.context.isin(contexts)]
    assert len(primary) == 86 and primary.status.eq("tested").all()
    primary.to_csv(RESULTS / "primary_86_edges_for_reporting.csv", index=False)
    order = sorted(primary.gene.unique())
    fig, axes = plt.subplots(1, 2, figsize=(8.3, 11.3), sharey=True)
    fig.subplots_adjust(left=.15, right=.97, top=.9, bottom=.17, wspace=.13)
    titles = ["A  Pooled immature dopaminergic cells",
              "B  Direct iDA specificity contrast"]
    for ax, context, title in zip(axes, contexts, titles):
        d = primary[primary.context.eq(context)].set_index("gene").loc[order]
        yy = np.arange(len(d))
        ax.errorbar(d.effect_z, yy,
                    xerr=np.vstack([d.effect_z-d.ci_low, d.ci_high-d.effect_z]),
                    fmt="o", color=TEAL, ecolor="#92B0BC", ms=3, lw=.8, capsize=1.5)
        ax.axvline(0, color=GREY, ls="--", lw=.8)
        ax.set_yticks(yy)
        ax.set_yticklabels(order, fontsize=8)
        ax.set_title(title, loc="left", fontsize=9, pad=12)
        ax.set_xlabel("Adjusted disease effect in Fisher z\n(iPD minus HC)", fontsize=9)
        ax.grid(axis="x", color="#E9EDF0", lw=.5)
        ax.set_axisbelow(True)
    axes[0].invert_yaxis()
    fig.suptitle("Primary SNCA–autophagy co-expression estimates", x=.15, ha="left",
                 fontsize=14, fontweight="bold", y=.965)
    fig.text(.15, .928, "8 HC and 28 iPD donors · 43 eligible genes · 0/86 edges with BH q < 0.05",
             fontsize=9)
    fig.text(.15, .038,
             "Points are donor-model disease coefficients; bars are pointwise 95% HC3 t intervals.\n"
             "Panel B models each donor's iDA-minus-non-iDA Fisher z difference. These are not Δr intervals.\n"
             "Wild-bootstrap P values use 19,999 draws; BH correction covers both panels together.\n"
             "Minimum primary q = 0.6081. Pointwise intervals and adjusted P values answer different questions.",
             fontsize=8, linespacing=1.45)
    save(fig, "01_primary_coupling_forest")

    correlations = pd.read_csv(RESULTS / "between_donor_and_bulk_correlations.csv")
    pearson = correlations[correlations.method.eq("Pearson")].copy()
    for arm in ["HC", "iPD"]:
        r = pearson["r_" + arm].to_numpy()
        z = np.arctanh(np.clip(r, -1 + 1e-7, 1 - 1e-7))
        half = stats.norm.ppf(.975) / np.sqrt(pearson["n_" + arm].to_numpy()-3)
        pearson["ci_r_" + arm + "_low"] = np.tanh(z-half)
        pearson["ci_r_" + arm + "_high"] = np.tanh(z+half)
    pearson["group_interval_note"] = (
        "Pointwise Fisher-normal 95% interval for each group's Pearson r; "
        "assumes independent donors and approximate bivariate normality. "
        "Difference intervals are the existing donor-bootstrap intervals."
    )
    pearson.to_csv(RESULTS / "between_donor_Pearson_with_group_CIs.csv", index=False)

    # Repeat only the deterministic two-gene aggregation used in script 05.
    meta = pd.read_csv(ROOT / "02-reanalysis/cell_meta.tsv.gz", sep="\t", index_col=0)
    cm = meta[meta.subtype.isin(["na_HC", "na_iPD"]) &
              meta.CellType.str.startswith("iDA") & ~meta.PPMI_ID.isin([3954, 4106])].copy()
    sizes = cm.SampleID.value_counts()
    cm = cm[cm.SampleID.isin(sizes[sizes >= 20].index)]
    expression = pd.read_csv(HERE / "inputs/rna_panel_logexpr.tsv.gz", sep="\t",
                             usecols=["cell", "SNCA", "GABARAPL1"], index_col="cell")
    samples = cm.drop_duplicates("SampleID").set_index("SampleID")
    single = expression.loc[cm.index].groupby(cm.SampleID).mean().groupby(samples.PPMI_ID).mean()
    single["group"] = samples.groupby("PPMI_ID").subtype.first().map({"na_HC":"HC", "na_iPD":"iPD"})
    single.index.name = "PPMI_ID"
    single.to_csv(RESULTS / "GABARAPL1_scrna_scatter_source.csv")
    bulk = pd.read_csv(RESULTS / "bulk_day65_donor_logcpm.csv", index_col="PPMI_ID")
    bulk["group"] = pd.read_csv(RESULTS / "bulk_day65_donor_metadata.csv",
                                index_col="PPMI_ID").group
    bulk[["SNCA", "GABARAPL1", "group"]].to_csv(RESULTS / "GABARAPL1_bulk_scatter_source.csv")
    single_row = pearson[pearson.dataset.eq("scRNA_iDA_pooled") &
                         pearson.scenario.eq("primary") & pearson.gene.eq("GABARAPL1")].iloc[0]
    bulk_row = pearson[pearson.dataset.eq("bulk_day65") &
                       pearson.scenario.eq("all_donors_batch_adjustment") &
                       pearson.gene.eq("GABARAPL1")].iloc[0]
    for values, row in [(single, single_row), (bulk, bulk_row)]:
        for arm in ["HC", "iPD"]:
            arm_values = values[values.group.eq(arm)]
            r = arm_values[["SNCA", "GABARAPL1"]].corr().iloc[0, 1]
            assert abs(r - row["r_" + arm]) < 1e-10
            assert len(arm_values) == row["n_" + arm]

    fig = plt.figure(figsize=(9, 8.7))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, .72], hspace=.53,
                           left=.12, right=.97, top=.88, bottom=.19, wspace=.29)
    for i, (values, row, title, unit) in enumerate([
        (single, single_row, "A  scRNA pooled iDA donor means", "Mean RNA log-normalized expression"),
        (bulk, bulk_row, "B  Bulk day-65 donor means", "Mean log2(normalized CPM + 0.5)"),
    ]):
        ax = fig.add_subplot(grid[0, i])
        for arm, color in [("HC", TEAL), ("iPD", ORANGE)]:
            a = values[values.group.eq(arm)]
            ax.scatter(a.SNCA, a.GABARAPL1, color=color, s=27, alpha=.85,
                       edgecolor="white", lw=.3,
                       label=f"{arm} (n={len(a)}), r={row['r_' + arm]:.2f}")
            coef = np.polyfit(a.SNCA, a.GABARAPL1, 1)
            xx = np.linspace(a.SNCA.min(), a.SNCA.max(), 50)
            ax.plot(xx, np.polyval(coef, xx), color=color, lw=1)
        ax.set_title(title, loc="left", fontsize=10)
        ax.set_xlabel("SNCA\n" + unit, fontsize=8)
        ax.set_ylabel("GABARAPL1", fontsize=9)
        ax.legend(frameon=False, fontsize=7.5, loc="best")
    adjusted = pd.read_csv(RESULTS / "bulk_refitted_adjusted_correlations.csv")
    a = adjusted[adjusted.scenario.eq("matched_cultures") & adjusted.gene.eq("GABARAPL1")]
    uncomposed, composed = a[a.adjust_cell_fraction.eq(False)].iloc[0], a[a.adjust_cell_fraction.eq(True)].iloc[0]
    rows = [
        ("scRNA iDA, unadjusted (8/28)", single_row.delta_r, single_row.ci_delta_low, single_row.ci_delta_high, TEAL),
        ("Bulk, unadjusted (8/32)", bulk_row.delta_r, bulk_row.ci_delta_low, bulk_row.ci_delta_high, ORANGE),
        ("Bulk, sex + batch (8/28)", uncomposed.delta_r, uncomposed.ci_low, uncomposed.ci_high, ORANGE),
        ("Bulk, + iDA fraction (8/28)", composed.delta_r, composed.ci_low, composed.ci_high, ORANGE),
    ]
    ax = fig.add_subplot(grid[1, :])
    ax.set_position([.43, .22, .53, .21])
    for yy, (_, effect, lo, hi, color) in enumerate(rows):
        ax.errorbar(effect, yy, xerr=[[effect-lo], [hi-effect]], color=color,
                    fmt="o", ms=4, lw=1.2, capsize=3)
    ax.set_yticks(range(4))
    ax.set_yticklabels([r[0] for r in rows], fontsize=8)
    ax.invert_yaxis()
    ax.axvline(0, ls="--", color=GREY, lw=.8)
    ax.set_xlabel("Correlation difference Δr (iPD minus HC), with 95% bootstrap CI")
    ax.set_title("C  Different assay directions and uncertain adjusted effects", loc="right",
                 fontsize=10, pad=12)
    fig.suptitle("SNCA–GABARAPL1 is an exploratory cross-assay lead",
                 x=.12, ha="left", y=.97, fontsize=14, fontweight="bold")
    fig.text(.12, .925, "All bulk controls overlap with scRNA; this is not independent replication.", fontsize=9)
    fig.text(.12, .045,
             "Each point in A/B is one donor; culture means are weighted equally. Lines are descriptive unadjusted fits.\n"
             "Raw bulk q = 0.00579 applies to 43 Pearson tests at day 65, not the whole sensitivity search.\n"
             "Adjusted bulk intervals refit nuisance models in each of 5,000 donor resamples; both cross zero.\n"
             "The scRNA iDA fraction is a composition proxy. Mixed-culture bulk cannot establish iDA specificity.",
             fontsize=8, linespacing=1.5)
    save(fig, "02_GABARAPL1_cross_assay")

    omnibus = pd.read_csv(RESULTS / "omnibus_multipipeline_sensitivity.csv")
    om = omnibus[omnibus.context.isin(contexts)].copy()
    method_order = ["rna_qc", "rna_raw", "spearman_qc", "sct_qc"]
    method_names = {"rna_qc":"RNA + QC (primary)", "rna_raw":"RNA unadjusted",
                    "spearman_qc":"Ranks + QC", "sct_qc":"SCT + QC"}
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.5))
    fig.subplots_adjust(left=.205, right=.97, top=.79, bottom=.23, wspace=.6)
    for i, context in enumerate(contexts):
        d = om[om.context.eq(context)].set_index("method").loc[method_order]
        for yy, method in enumerate(method_order):
            raw = d.loc[method, "p_omnibus"]
            adj = d.loc[method, "p_holm_across_all_8_method_context_tests"]
            axes[i].plot([raw, adj], [yy, yy], color="#A6B5BD", lw=1)
        axes[i].scatter(d.p_omnibus, range(4), color=TEAL, s=28, label="Nominal P")
        axes[i].scatter(d.p_holm_across_all_8_method_context_tests, range(4),
                        color=ORANGE, s=28, marker="s", label="Holm, 8 tests")
        axes[i].axvline(.05, color=GREY, ls="--", lw=.8)
        axes[i].set_xscale("log")
        axes[i].set_xlim(.009, .75)
        axes[i].set_xticks([.01, .05, .1, .5])
        axes[i].set_xticklabels(["0.01", "0.05", "0.1", "0.5"])
        axes[i].set_yticks(range(4))
        axes[i].set_yticklabels([method_names[x] for x in method_order] if i == 0 else [])
        axes[i].set_ylim(3.5, -.5)
        axes[i].set_xlabel("Omnibus P value")
        axes[i].set_title(["A  Pooled iDA", "B  Direct iDA specificity"][i], loc="left", fontsize=11)
    axes[1].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0, -.2), ncol=2)
    fig.suptitle("The panel signal depends on the analysis method", x=.08, ha="left",
                 y=.95, fontsize=14, fontweight="bold")
    fig.text(.08, .86, "Primary two-test Holm P values are 0.2438 (pooled iDA) and 0.4213 (specificity).", fontsize=9)
    fig.text(.08, .048,
             "The additional eight-test correction covers four methods and both contexts. Its minimum is 0.0996.\n"
             "Rank results remain exploratory. In common batches without the reference donor (6 HC / 18 iPD),\n"
             "rank-based exact permutation gives P = 0.0260 for pooled iDA and P = 0.1484 for specificity.",
             fontsize=8, linespacing=1.5)
    save(fig, "03_method_sensitivity")

    # Display all 43 primary genes, avoiding selection of a favorable donor-influence example.
    loo = pd.read_csv(RESULTS / "leave_one_donor_out.csv")
    loo = loo[loo.context.eq("iDA_pooled") & loo.status.eq("tested")]
    main = primary[primary.context.eq("iDA_pooled")].set_index("gene").loc[order]
    bounds = loo.groupby("gene").effect_z.agg(["min", "max", "count"]).loc[order]
    bounds["full_sample_effect_z"] = main.effect_z
    bounds.to_csv(RESULTS / "pooled_iDA_donor_influence_summary.csv")
    fig, ax = plt.subplots(figsize=(7, 10.7))
    fig.subplots_adjust(left=.18, right=.95, top=.9, bottom=.125)
    yy = np.arange(len(order))
    ax.hlines(yy, bounds["min"], bounds["max"], color="#90ADB9", lw=2)
    ax.scatter(main.effect_z, yy, color=TEAL, s=13, zorder=3)
    ax.axvline(0, color=GREY, ls="--", lw=.8)
    ax.set_yticks(yy)
    ax.set_yticklabels(order, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Adjusted disease effect in Fisher z")
    fig.suptitle("Donor influence across the complete primary panel",
                 x=.18, ha="left", y=.965, fontsize=12, fontweight="bold")
    fig.text(.18, .925, "Dots are full-cohort estimates; bars span estimable leave-one-donor-out fits.", fontsize=8)
    fig.text(.18, .038,
             "Bars are influence ranges, not confidence intervals. Each refit removes one independent donor.\n"
             "Removing the reference donor alone saturates a batch stratum; that HC3 fit is non-estimable\n"
             "and is explicitly excluded from these ranges. It is retained with its status in the full table.",
             fontsize=8, linespacing=1.5)
    save(fig, "04_donor_influence")

    (FIGURES / "figure_validation.json").write_text(json.dumps({
        "primary_edges": len(primary), "primary_significant_edges": int((primary.q_wild_family < .05).sum()),
        "scatter_correlations_match_result_tables": True,
        "scatter_independent_donors": {"scRNA":len(single), "bulk_day65":len(bulk)},
        "figures": 4, "formats":["png", "pdf"],
        "group_CI_method":"Pointwise Fisher-normal, descriptive",
        "visual_inspection":"pending",
    }, indent=2) + "\n")
    print("Four source-backed figure pairs and descriptive group-r intervals written.")


if __name__ == "__main__":
    main()
