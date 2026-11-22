#!/usr/bin/env python3
"""
Combined-modality ridge plot of effect-size distributions per intervention
category, pooled across mouse RNA-seq + human RNA-seq + human DNAm. Single
KDE per category (no per-modality stratification).

Outputs:
  figures/create_ridge_plot_combined.png
  figures/create_ridge_plot_combined.pdf
"""

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

# Load combined 3-modality data
df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)

# Analytic set: exclude extreme-tail-flagged rows (mirrors analysis_2/3/4_combined)
if "outlier_status" in df.columns:
    df = df[df["outlier_status"] == "Normal"].copy()

# Canonicalise category to canonical 5
CANONICAL = {"Disease", "Drug", "Environment", "Genetic", "Other"}
df["condition_category"] = df["condition_category"].where(
    df["condition_category"].isin(CANONICAL), "Other"
)

# Significant rows only (FDR<0.05, matches the figure caption), pooled across modalities
df_sig = df[df["fdr"] < 0.05].copy()
df_sig = df_sig.dropna(subset=["log2FoldChange"])

# Order categories by median effect size (log2FoldChange)
category_order = (
    df_sig.groupby("condition_category")["log2FoldChange"]
    .median()
    .sort_values()
    .index.tolist()
)
n_categories = len(category_order)

# Statistical tests (pooled)
groups = [
    df_sig[df_sig["condition_category"] == cat]["log2FoldChange"].values
    for cat in category_order
]
kruskal_stat, kruskal_p = stats.kruskal(*groups)

positive_effects = df_sig[df_sig["log2FoldChange"] > 0]["log2FoldChange"].values
negative_effects = df_sig[df_sig["log2FoldChange"] < 0]["log2FoldChange"].values
ttest_stat, ttest_p = stats.ttest_ind(
    np.abs(positive_effects), np.abs(negative_effects)
)

COLORS = {
    "Disease": "#A23B72",
    "Drug": "#2E86AB",
    "Environment": "#C73E1D",
    "Genetic": "#F18F01",
    "Other": "#6B7A8F",
}

sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
fig, axs = plt.subplots(nrows=n_categories, ncols=1, figsize=(9, 6), sharex=True)
axs = axs.flatten() if n_categories > 1 else [axs]

bandwidth = 0.5

for i, category in enumerate(category_order):
    cat_data = df_sig[df_sig["condition_category"] == category]["log2FoldChange"].values

    sns.kdeplot(cat_data, ax=axs[i], bw_adjust=bandwidth,
                fill=True, alpha=0.7, color=COLORS[category],
                linewidth=2, edgecolor="white")
    sns.kdeplot(cat_data, ax=axs[i], bw_adjust=bandwidth,
                fill=False, color=COLORS[category],
                linewidth=2.5, alpha=0.9)

    median = np.median(cat_data)
    n_samples = len(cat_data)
    # Wilcoxon vs zero (only if non-zero entries)
    nonzero = cat_data[cat_data != 0]
    if len(nonzero) > 0:
        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(nonzero)
    else:
        wilcoxon_stat, wilcoxon_p = (np.nan, 1.0)

    axs[i].axvline(0, color="black", linestyle="--", alpha=0.5, linewidth=1)
    axs[i].axvline(median, color=COLORS[category], linestyle="-",
                   alpha=0.8, linewidth=2)

    axs[i].text(-0.02, 0.5, f"{category}",
                transform=axs[i].transAxes,
                fontsize=12, fontweight="bold",
                verticalalignment="center", horizontalalignment="right")
    axs[i].text(-0.02, 0.2, f"(n={n_samples:,})",
                transform=axs[i].transAxes,
                fontsize=9, verticalalignment="center",
                horizontalalignment="right", style="italic")

    if wilcoxon_p < 0.0001:
        p_text = f"p={wilcoxon_p:.2e}"
    else:
        p_text = f"p={wilcoxon_p:.4f}"

    axs[i].text(1.02, 0.5, p_text,
                transform=axs[i].transAxes,
                ha="left", va="center",
                fontsize=10, fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=0.4",
                          facecolor=COLORS[category],
                          edgecolor=COLORS[category], alpha=0.9))

    axs[i].set_yticks([])
    axs[i].set_ylabel("")
    axs[i].set_xlim(-25, 25)
    axs[i].spines["left"].set_visible(False)
    axs[i].spines["right"].set_visible(False)
    axs[i].spines["top"].set_visible(False)

    if i < n_categories - 1:
        axs[i].spines["bottom"].set_visible(False)
        axs[i].set_xticks([])
    else:
        axs[i].spines["bottom"].set_visible(True)
        axs[i].set_xlabel(
            "log2 Fold Change (Age Acceleration)",
            fontsize=12, fontweight="bold",
        )
        axs[i].set_xticks([-20, -10, 0, 10, 20])

plt.subplots_adjust(hspace=-0.5)

total_sig = len(df_sig)
fig.text(0, 1,
         f"Effect size distribution of {total_sig:,} significant interventions "
         f"(mouse RNA-seq + human RNA-seq + human DNAm pooled)",
         ha="left", fontsize=11, style="italic")

legend_elements = [
    mpatches.Patch(color="none", label="Reference lines:"),
    plt.Line2D([0], [0], color="black", linestyle="--", alpha=0.5,
               label="No effect (log2FC = 0)"),
    plt.Line2D([0], [0], color="gray", linestyle="-", linewidth=2,
               label="Category median"),
]
axs[0].legend(handles=legend_elements, loc="upper right",
              frameon=False, fontsize=10)

fig.text(0.25, 0.1, "Decelerating Aging (left)",
         fontsize=11, color="darkgray", fontweight="bold", ha="center")
fig.text(0.9, 0.1, "Accelerating Aging (right)",
         fontsize=11, color="darkgray", fontweight="bold", ha="center")

plt.tight_layout()
plt.savefig("figures/create_ridge_plot_combined.pdf", dpi=300, bbox_inches="tight")
plt.savefig("figures/create_ridge_plot_combined.png", dpi=600, bbox_inches="tight")
print("Saved figures/create_ridge_plot_combined.{pdf,png}")

print("\nRidge Plot Statistics (combined modalities)")
print("=" * 60)
print(f"\nTotal significant interventions: {total_sig:,}")
print(f"\nKruskal-Wallis test (categories): H={kruskal_stat:.4f}, p={kruskal_p:.2e}")
print(f"\nT-test |+effects| vs |-effects|: t={ttest_stat:.4f}, p={ttest_p:.2e}")
print(f"  positive: n={len(positive_effects):,}, mean(abs)={np.mean(np.abs(positive_effects)):.3f}")
print(f"  negative: n={len(negative_effects):,}, mean(abs)={np.mean(np.abs(negative_effects)):.3f}")

print("\nCategory-wise statistics:")
print("-" * 60)
for category in category_order:
    cat_data = df_sig[df_sig["condition_category"] == category]["log2FoldChange"].values
    nonzero = cat_data[cat_data != 0]
    if len(nonzero) > 0:
        w_stat, w_p = stats.wilcoxon(nonzero)
    else:
        w_stat, w_p = (np.nan, 1.0)
    print(f"\n{category}: n={len(cat_data):,}, median={np.median(cat_data):.3f}, mean={np.mean(cat_data):.3f}")
    print(f"  Wilcoxon vs 0: W={w_stat:.1f}, p={w_p:.6e}")
