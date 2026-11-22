#!/usr/bin/env python3
"""
ANALYSIS 2 (combined modalities): Genetic Intervention Type Patterns
====================================================================

Pools mouse RNA-seq + human RNA-seq + human DNAm into a single set of
genetic interventions, classifies them as Knockout / Knockdown /
Overexpression / Mutation, and produces:

  figures/analysis_2_genetic_intervention_types_combined.png
  figures/analysis_2_genetic_intervention_types_combined.pdf

Score-floor per modality (mirrors combined volcano):
    mouse RNA-seq : 100
    human RNA-seq : 100
    human DNAm    :  80
"""

import re
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import binomtest, chi2_contingency, fisher_exact, gaussian_kde
from statsmodels.stats.multitest import fdrcorrection

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")

SCORE_FLOOR = {"mouse RNA-seq": 100, "human RNA-seq": 100, "human DNAm": 80}


def load_data():
    """Load the pre-built combined 3-modality table."""
    print("Loading combined 3-modality data...")
    df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)
    print(f"  raw rows: {len(df):,}")
    # Already pre-filtered to outlier_status == 'Normal' in builder, but be safe.
    if "outlier_status" in df.columns:
        df = df[df["outlier_status"] == "Normal"].copy()

    # No SCORE_FLOOR applied here: the canonical atlas is outlier_status=='Normal'
    # only (104,243 rows). Applying a score floor of 100 would drop 98% of the data
    # and contradict the combined volcano, which also uses all Normal rows.
    print(f"  rows after outlier filter: {len(df):,}")
    print(df["modality"].value_counts().to_string())
    return df


def classify_genetic_interventions(main_df):
    print("\n" + "=" * 80)
    print("ANALYSIS 2 (combined): GENETIC INTERVENTION TYPE PATTERNS")
    print("=" * 80)

    genetic_df = main_df[
        (main_df["condition_category"] == "Genetic") & (main_df["fdr"] < 0.05)
    ].copy()
    print(f"\nAnalyzing {len(genetic_df):,} significant genetic interventions (combined)")

    intervention_types = {
        "Knockout": ["knockout", "knock-out", "ko", "null", "deletion", "deleted", "-/-"],
        "Knockdown": ["knockdown", "knock-down", "kd", "shrna", "sirna", "rnai"],
        "Overexpression": ["overexpression", "overexpress", "transgenic", "tg", "oe"],
        "Mutation": ["mutant", "mutation", "mut", "variant"],
    }

    genetic_df["intervention_type"] = "Other"
    for int_type, keywords in intervention_types.items():
        pattern = "|".join(keywords)
        matches = genetic_df["condition_name"].astype(str).str.contains(
            pattern, case=False, na=False
        )
        genetic_df.loc[matches, "intervention_type"] = int_type

    results = []
    for int_type in ["Knockout", "Knockdown", "Overexpression", "Mutation"]:
        type_df = genetic_df[genetic_df["intervention_type"] == int_type]
        if len(type_df) == 0:
            continue
        results.append({
            "Type": int_type,
            "N_Interventions": len(type_df),
            "Mean_Effect": type_df["log2FoldChange"].mean(),
            "Median_Effect": type_df["log2FoldChange"].median(),
            "Prop_Anti_Aging": (type_df["log2FoldChange"] < 0).mean(),
            "N_Anti_Aging": int((type_df["log2FoldChange"] < 0).sum()),
            "N_Pro_Aging": int((type_df["log2FoldChange"] > 0).sum()),
        })

    results_df = pd.DataFrame(results).sort_values("Prop_Anti_Aging", ascending=False)
    print("\nGenetic Intervention Type Analysis (combined):")
    print(results_df.to_string(index=False))
    results_df.to_csv("result/analysis_2_genetic_types_combined.csv", index=False)
    print("Saved: result/analysis_2_genetic_types_combined.csv")
    return genetic_df, results_df


def perform_statistical_tests(genetic_df, results_df):
    print("\n" + "-" * 60)
    print("STATISTICAL TESTS (combined)")
    print("-" * 60)
    if len(results_df) < 2:
        return np.nan, 1.0, pd.DataFrame()
    contingency = [[r["N_Anti_Aging"], r["N_Pro_Aging"]] for _, r in results_df.iterrows()]
    chi2, p_chi2, dof, _ = chi2_contingency(contingency)
    print(f"Chi-square: chi2={chi2:.2f}, p={p_chi2:.2e}, dof={dof}")

    types = results_df["Type"].tolist()
    rows = []
    for i in range(len(types)):
        for j in range(i + 1, len(types)):
            t1 = genetic_df[genetic_df["intervention_type"] == types[i]]
            t2 = genetic_df[genetic_df["intervention_type"] == types[j]]
            if len(t1) and len(t2):
                stat, p = stats.mannwhitneyu(
                    t1["log2FoldChange"], t2["log2FoldChange"],
                    alternative="two-sided",
                )
                rows.append({
                    "Comparison": f"{types[i]} vs {types[j]}",
                    "Type1": types[i], "Type2": types[j],
                    "P_value": p,
                    "Median1": t1["log2FoldChange"].median(),
                    "Median2": t2["log2FoldChange"].median(),
                })
    comp_df = pd.DataFrame(rows)
    if len(comp_df):
        _, comp_df["FDR"] = fdrcorrection(comp_df["P_value"])
        comp_df = comp_df.sort_values("FDR")
        comp_df.to_csv("result/analysis_2_genetic_pairwise_tests_combined.csv", index=False)
        print("Saved: result/analysis_2_genetic_pairwise_tests_combined.csv")
    return chi2, p_chi2, comp_df


def create_figure(genetic_df, results_df, main_df):
    """Combined-modality figure: bar plot (FC vs baseline) + density plot."""
    print("\nCreating combined bar + density figure...")

    # Baseline = proportion of anti-aging among significant + classified interventions
    # (same population as the per-type bars, so FC=1 is the meaningful reference)
    if len(results_df) == 0:
        print("  no data to plot"); return
    total_n = int(results_df["N_Interventions"].sum())
    total_anti = int(results_df["N_Anti_Aging"].sum())
    baseline_prop = total_anti / total_n if total_n else 0.5
    print(f"  baseline (sig + classified, n={total_n:,}): {baseline_prop:.3f}")

    results_df = results_df.copy()
    results_df["Fold_Change"] = results_df["Prop_Anti_Aging"] / baseline_prop
    results_df = results_df.sort_values("Fold_Change", ascending=True)

    p_values = []
    for _, row in results_df.iterrows():
        n_success = int(row["N_Anti_Aging"])
        n_total = int(row["N_Interventions"])
        if n_total == 0:
            p_values.append(1.0)
            continue
        result = binomtest(n_success, n_total, baseline_prop, alternative="two-sided")
        p_values.append(result.pvalue)
    results_df["P_value"] = p_values
    results_df["Significant"] = results_df["P_value"] < 0.05

    types = results_df["Type"].tolist()
    fold_changes = results_df["Fold_Change"].tolist()

    color_map = {}
    colors = []
    for typ, fc in zip(types, fold_changes):
        if fc > 1.0:
            intensity = min((fc - 1.0) / 0.5, 1.0)
            color = plt.cm.Blues(0.5 + 0.5 * intensity)
        else:
            intensity = min((1.0 - fc) / 0.5, 1.0)
            color = plt.cm.YlOrRd(0.5 + 0.5 * intensity)
        colors.append(color)
        color_map[typ] = color

    # Two-panel figure (bar + density)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 5))
    fig.suptitle(
        "Genetic intervention types (mouse RNA-seq + human RNA-seq + human DNAm pooled)",
        fontsize=14, fontweight="bold",
    )

    # Bar plot
    y_pos = np.arange(len(types))
    bars = ax1.barh(y_pos, fold_changes, color=colors, alpha=0.6,
                    edgecolor="black", linewidth=1.2)
    ax1.axvline(1.0, color="black", linestyle="--", linewidth=2, alpha=0.7,
                label=f"Baseline (overall mean: {baseline_prop*100:.1f}%)")
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(types, fontsize=11, fontweight="bold")
    ax1.set_xlabel(
        f"Enrichment of age-lowering signatures (baseline = {baseline_prop*100:.1f}%)",
        fontsize=12, fontweight="bold",
    )
    for i, (bar, row) in enumerate(zip(bars, results_df.iterrows())):
        fc = row[1]["Fold_Change"]
        p_val = row[1]["P_value"]
        is_sig = row[1]["Significant"]
        n = int(row[1]["N_Interventions"])
        sig_marker = "*" if is_sig else ""
        text_color = "black" if is_sig else "grey"
        ax1.text(fc + 0.05, i,
                 f"FC={fc:.2f} (n={n}), p={p_val:.2e}{sig_marker}",
                 va="center", ha="left", fontsize=9,
                 fontweight="bold", color=text_color)
    ax1.legend(loc="lower right", fontsize=10, frameon=True, fancybox=True, shadow=True)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(False)
    ax1.set_xlim(0, max(max(fold_changes) * 1.5 if fold_changes else 1.5, 1.6))

    # Density plot
    for typ in types:
        type_df = genetic_df[genetic_df["intervention_type"] == typ]
        log2fc_values = type_df["log2FoldChange"].dropna().values
        if len(log2fc_values) > 1:
            kde = gaussian_kde(log2fc_values)
            x_range = np.linspace(log2fc_values.min(), log2fc_values.max(), 300)
            density = kde(x_range)
            ax2.fill_between(x_range, 0, density,
                             color=color_map[typ], alpha=0.2,
                             label=f"{typ} (n={len(log2fc_values)})")
            ax2.plot(x_range, density, color=color_map[typ],
                     linewidth=2.5, alpha=0.8)
    ax2.axvline(0, color="black", linestyle="--", linewidth=2.5,
                alpha=0.8, label="No change (log2FC=0)")
    ax2.set_xlabel("Effect Size", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Density", fontsize=12, fontweight="bold")
    ax2.legend(loc="upper right", fontsize=10, frameon=True, fancybox=True, shadow=True)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(False)

    plt.tight_layout()
    plt.savefig(
        "figures/analysis_2_genetic_intervention_types_combined.pdf",
        dpi=300, bbox_inches="tight",
    )
    plt.savefig(
        "figures/analysis_2_genetic_intervention_types_combined.png",
        dpi=300, bbox_inches="tight",
    )
    plt.close()
    print("Saved: figures/analysis_2_genetic_intervention_types_combined.{pdf,png}")

    results_df.to_csv("result/analysis_2_genetic_types_fold_change_combined.csv", index=False)
    print("Saved: result/analysis_2_genetic_types_fold_change_combined.csv")


def main():
    print("\n" + "=" * 80)
    print("GENETIC INTERVENTION TYPE ANALYSIS (combined modalities)")
    print("=" * 80)
    main_df = load_data()
    genetic_df, results_df = classify_genetic_interventions(main_df)
    if len(results_df) == 0:
        print("No significant genetic interventions to analyse.")
        return
    perform_statistical_tests(genetic_df, results_df)
    create_figure(genetic_df, results_df, main_df)
    print("\nDone.")


if __name__ == "__main__":
    main()
