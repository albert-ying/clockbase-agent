#!/usr/bin/env python3
"""
Create Nature publication-quality ridge plot showing effect size distributions
across intervention categories with statistical test annotations
"""

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle
from scipy import stats

# Load data
df = pd.read_csv("figures/combined_volcano_plot_data.csv")

# Filter for significant interventions only
df_sig = df[df["is_significant"] == True].copy()

# Order categories by median effect size (log2FoldChange)
category_order = (
    df_sig.groupby("condition_category")["log2FoldChange"]
    .median()
    .sort_values()
    .index.tolist()
)
n_categories = len(category_order)

# Statistical tests
# 1. Kruskal-Wallis test (non-parametric ANOVA)
groups = [
    df_sig[df_sig["condition_category"] == cat]["log2FoldChange"].values
    for cat in category_order
]
kruskal_stat, kruskal_p = stats.kruskal(*groups)

# 2. T-test comparing positive vs negative effect sizes
positive_effects = df_sig[df_sig["log2FoldChange"] > 0]["log2FoldChange"].values
negative_effects = df_sig[df_sig["log2FoldChange"] < 0]["log2FoldChange"].values
ttest_stat, ttest_p = stats.ttest_ind(
    np.abs(positive_effects), np.abs(negative_effects)
)

# Color scheme
COLORS = {
    "Disease": "#A23B72",
    "Drug": "#2E86AB",
    "Environment": "#C73E1D",
    "Genetic": "#F18F01",
    "Other": "#6B7A8F",
}

# Set style
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

# Create figure
fig, axs = plt.subplots(nrows=n_categories, ncols=1, figsize=(9, 6), sharex=True)
axs = axs.flatten()

# Plot parameters
bandwidth = 0.5  # Smoothness of KDE

# Create ridge plots
for i, category in enumerate(category_order):
    # Get data for this category
    cat_data = df_sig[df_sig["condition_category"] == category]["log2FoldChange"].values

    # Plot KDE
    sns.kdeplot(
        cat_data,
        ax=axs[i],
        bw_adjust=bandwidth,
        fill=True,
        alpha=0.7,
        color=COLORS[category],
        linewidth=2,
        edgecolor="white",
    )

    # Add darker outline
    sns.kdeplot(
        cat_data,
        ax=axs[i],
        bw_adjust=bandwidth,
        fill=False,
        color=COLORS[category],
        linewidth=2.5,
        alpha=0.9,
    )

    # Calculate statistics
    median = np.median(cat_data)
    mean = np.mean(cat_data)
    n_samples = len(cat_data)

    # Wilcoxon signed-rank test (test if distribution differs from zero)
    wilcoxon_stat, wilcoxon_p = stats.wilcoxon(cat_data)

    # Add vertical reference lines
    axs[i].axvline(0, color="black", linestyle="--", alpha=0.5, linewidth=1)
    axs[i].axvline(
        median, color=COLORS[category], linestyle="-", alpha=0.8, linewidth=2
    )

    # Add category label and sample size
    axs[i].text(
        -0.02,
        0.5,
        f"{category}",
        transform=axs[i].transAxes,
        fontsize=12,
        fontweight="bold",
        verticalalignment="center",
        horizontalalignment="right",
    )

    # Add sample size
    axs[i].text(
        -0.02,
        0.2,
        f"(n={n_samples:,})",
        transform=axs[i].transAxes,
        fontsize=9,
        verticalalignment="center",
        horizontalalignment="right",
        style="italic",
    )

    # Format p-value text
    if wilcoxon_p < 0.0001:
        p_text = f"p={wilcoxon_p:.2e}"
    else:
        p_text = f"p={wilcoxon_p:.4f}"

    # Add p-value on the right side with category theme color
    # Position it to the right of the plot with padding
    axs[i].text(
        1.02,
        0.5,
        p_text,
        transform=axs[i].transAxes,
        ha="left",
        va="center",
        fontsize=10,
        fontweight="bold",
        color="white",
        bbox=dict(
            boxstyle="round,pad=0.4",
            facecolor=COLORS[category],
            edgecolor=COLORS[category],
            alpha=0.9,
        ),
    )

    # Remove y-axis
    axs[i].set_yticks([])
    axs[i].set_ylabel("")

    # Set x-limits
    axs[i].set_xlim(-25, 25)

    # Clean up spines
    axs[i].spines["left"].set_visible(False)
    axs[i].spines["right"].set_visible(False)
    axs[i].spines["top"].set_visible(False)

    # Only show bottom spine for last plot
    if i < n_categories - 1:
        axs[i].spines["bottom"].set_visible(False)
        axs[i].set_xticks([])
    else:
        axs[i].spines["bottom"].set_visible(True)
        axs[i].set_xlabel(
            "log₂ Fold Change (Age Acceleration)", fontsize=12, fontweight="bold"
        )
        axs[i].set_xticks([-20, -10, 0, 10, 20])

# Adjust spacing
plt.subplots_adjust(hspace=-0.5)

# Add subtitle with sample size
total_sig = len(df_sig)
fig.text(
    0,
    1,
    f"Effect size distribution of {total_sig:,} significant interventions",
    ha="left",
    fontsize=12,
    style="italic",
)


# Add legend for reference lines
legend_elements = [
    mpatches.Patch(color="none", label="Reference lines:"),
    plt.Line2D(
        [0],
        [0],
        color="black",
        linestyle="--",
        alpha=0.5,
        label="No effect (log2FC = 0)",
    ),
    plt.Line2D(
        [0], [0], color="gray", linestyle="-", linewidth=2, label="Category median"
    ),
]
axs[0].legend(handles=legend_elements, loc="upper right", frameon=False, fontsize=10)

# Add arrows and labels for aging direction
fig.text(
    0.25,
    0.1,
    "← Decelerating Aging",
    fontsize=11,
    color="darkgray",
    fontweight="bold",
    ha="center",
)
fig.text(
    0.9,
    0.1,
    "Accelerating Aging →",
    fontsize=11,
    color="darkgray",
    fontweight="bold",
    ha="center",
)

# Save figure
plt.tight_layout()
plt.savefig("figures/ridge_plot_effect_sizes.pdf", dpi=300, bbox_inches="tight")
plt.savefig("figures/ridge_plot_effect_sizes.png", dpi=600, bbox_inches="tight")
plt.show()

# Print detailed statistics
print("\nRidge Plot Statistics")
print("=" * 60)
print(f"\nTotal significant interventions: {total_sig:,}")
print(f"\nKruskal-Wallis test (comparing distributions across categories):")
print(f"  H-statistic: {kruskal_stat:.4f}")
print(f"  p-value: {kruskal_p:.2e}")

print(f"\nT-test (|positive effects| vs |negative effects|):")
print(
    f"  Positive effects (accelerating): n={len(positive_effects):,}, mean(abs)={np.mean(np.abs(positive_effects)):.3f}"
)
print(
    f"  Negative effects (decelerating): n={len(negative_effects):,}, mean(abs)={np.mean(np.abs(negative_effects)):.3f}"
)
print(f"  t-statistic: {ttest_stat:.4f}")
print(f"  p-value: {ttest_p:.2e}")

print("\nCategory-wise statistics:")
print("-" * 60)
for category in category_order:
    cat_data = df_sig[df_sig["condition_category"] == category]["log2FoldChange"].values
    wilcoxon_stat, wilcoxon_p = stats.wilcoxon(cat_data)
    print(f"\n{category}:")
    print(f"  n = {len(cat_data):,}")
    print(f"  Median: {np.median(cat_data):.3f}")
    print(f"  Mean: {np.mean(cat_data):.3f}")
    print(f"  Std: {np.std(cat_data):.3f}")
    print(f"  Min: {np.min(cat_data):.3f}")
    print(f"  Max: {np.max(cat_data):.3f}")
    print(
        f"  % Positive (accelerating): {(cat_data > 0).sum() / len(cat_data) * 100:.1f}%"
    )
    print(
        f"  % Negative (decelerating): {(cat_data < 0).sum() / len(cat_data) * 100:.1f}%"
    )
    print(f"  Wilcoxon test vs 0: W={wilcoxon_stat:.1f}, p={wilcoxon_p:.6e}")

# Perform pairwise comparisons
print("\nPairwise comparisons (Mann-Whitney U test):")
print("-" * 60)
for i, cat1 in enumerate(category_order):
    for cat2 in category_order[i + 1 :]:
        data1 = df_sig[df_sig["condition_category"] == cat1]["log2FoldChange"].values
        data2 = df_sig[df_sig["condition_category"] == cat2]["log2FoldChange"].values
        u_stat, p_val = stats.mannwhitneyu(data1, data2, alternative="two-sided")
        if p_val < 0.05:
            print(f"{cat1} vs {cat2}: U={u_stat:.1f}, p={p_val:.2e} ***")
        else:
            print(f"{cat1} vs {cat2}: U={u_stat:.1f}, p={p_val:.3f}")
