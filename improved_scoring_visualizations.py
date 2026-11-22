import warnings

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from aquarel import load_theme
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore")

# Apply consistent theme
theme = (
    load_theme("scientific")
    .set_grid(draw=False)
    .set_font(size=14)
    .set_ticks(direction="out")
    .set_axis_labels(pad=10)
)
theme.apply()
plt.rcParams["savefig.transparent"] = True
plt.rcParams["figure.facecolor"] = "none"
plt.rcParams["axes.facecolor"] = "none"
plt.rcParams["savefig.facecolor"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

# Load the scoring data
df_scores = pd.read_csv("data/scores_summary.csv")

# Load the detailed results for FDR filtering
df_results = pd.read_csv("data/detailed_scores_result_mm_rna.csv")

# Filter based on FDR < 0.05
df_results_sig = df_results[df_results["fdr"] < 0.05].copy()

# Merge with scores
df_merged = pd.merge(
    df_results_sig,
    df_scores,
    left_on="case_id",
    right_on="dataset",
    how="inner",
    suffixes=("_results", "_scores"),
)

print(f"Total scored interventions: {len(df_scores):,}")
print(f"FDR < 0.05 interventions: {len(df_merged):,}")


# ====================
# Figure 1: Final Score Distribution
# ====================
fig1 = plt.figure(figsize=(4, 3))
ax1 = plt.subplot(111)

scores = df_scores["final_score"]
bins = np.arange(-20, 145, 5)
counts, edges, patches = ax1.hist(
    scores, bins=bins, color="#2E86AB", alpha=0.7, edgecolor="black", linewidth=1.2
)

ax1.axvline(
    scores.median(),
    color="#A23B72",
    linestyle="--",
    linewidth=2.5,
    label=f"Median = {scores.median():.0f}",
)
ax1.axvline(
    100,
    color="#F18F01",
    linestyle="--",
    linewidth=2.5,
    label="High-score threshold (100)",
)

ax1.set_xlabel("Final Score", fontsize=16, fontweight="bold")
ax1.set_ylabel("Number of Interventions", fontsize=16, fontweight="bold")
ax1.set_title("Distribution of Final Scores", fontsize=18, fontweight="bold", pad=20)
ax1.legend(frameon=False, fontsize=12, loc="upper right")
ax1.set_xlim(-20, 145)

# Add statistics text
stats_text = f"N = {len(scores):,}\nHigh-scoring (≥100): {len(scores[scores >= 100]):,} ({len(scores[scores >= 100])/len(scores)*100:.1f}%)"
ax1.text(
    0.02,
    0.98,
    stats_text,
    transform=ax1.transAxes,
    fontsize=11,
    verticalalignment="top",
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
)

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig1_score_distribution.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig1_score_distribution.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 2: Score Components Box Plot
# ====================
fig2 = plt.figure(figsize=(4, 3))
ax2 = plt.subplot(111)

score_components = [
    "q1_experimental_model",
    "q2_1_aging_pathways",
    "q2_2_age_diseases",
    "q3_1_experimental_rigor",
    "q4_1_translational_potential",
    "q5_1_literature_saturation",
    "q5_2_mechanistic_novelty",
]
component_labels = [
    "Model\nQuality",
    "Aging\nPathways",
    "Age\nDiseases",
    "Experimental\nRigor",
    "Translational\nPotential",
    "Literature\nSaturation",
    "Mechanistic\nNovelty",
]
component_data = [df_scores[col].values for col in score_components]

bp = ax2.boxplot(
    component_data,
    labels=component_labels,
    patch_artist=True,
    medianprops=dict(color="black", linewidth=2.5),
    boxprops=dict(facecolor="#6FAADB", alpha=0.7),
    flierprops=dict(marker="o", markersize=3, alpha=0.3),
    widths=0.6,
)

ax2.set_ylabel("Score Points", fontsize=16, fontweight="bold")
ax2.set_title(
    "Distribution of Scoring Components", fontsize=18, fontweight="bold", pad=20
)
ax2.tick_params(axis="x", rotation=0, labelsize=12)
ax2.grid(axis="y", alpha=0.3, linestyle="--")

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig2_score_components.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig2_score_components.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 3: Plausibility Filter Impact
# ====================
fig3 = plt.figure(figsize=(10, 7))
ax3 = plt.subplot(111)

plausibility_counts = df_scores["q0_plausibility"].value_counts()
colors_map = {-50: "#C73E1D", 0: "#6FAADB", 30: "#50C878"}
bars = []
labels = []
colors_list = []

for val in [-50, 0, 30]:
    if val in plausibility_counts.index:
        bars.append(plausibility_counts[val])
        if val == -50:
            labels.append("Penalty\n(-50 points)")
        elif val == 0:
            labels.append("Neutral\n(0 points)")
        else:
            labels.append("Bonus\n(+30 points)")
        colors_list.append(colors_map[val])

bar_plot = ax3.bar(
    range(len(bars)),
    bars,
    color=colors_list,
    alpha=0.8,
    edgecolor="black",
    linewidth=1.5,
)
ax3.set_xticks(range(len(labels)))
ax3.set_xticklabels(labels, fontsize=14)
ax3.set_ylabel("Number of Interventions", fontsize=16, fontweight="bold")
ax3.set_title(
    "Biological Plausibility Assessment", fontsize=18, fontweight="bold", pad=20
)

# Add percentage labels on bars
for i, (bar, count) in enumerate(zip(bar_plot, bars)):
    pct = count / len(df_scores) * 100
    ax3.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 100,
        f"{count:,}\n({pct:.1f}%)",
        ha="center",
        va="bottom",
        fontweight="bold",
        fontsize=12,
    )

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig3_plausibility_filter.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig3_plausibility_filter.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 4: Model Quality Distribution (Pie Chart)
# ====================
fig4 = plt.figure(figsize=(5, 4))
ax4 = plt.subplot(111)

model_mapping = {
    25: "In vivo",
    15: "Primary cells",
    0: "Cell lines",
}
model_counts = (
    df_scores["q1_experimental_model"].value_counts().sort_index(ascending=False)
)
model_labels = []
model_values = []
for score in [25, 15, 0]:
    if score in model_counts.index:
        model_labels.append(f"{model_mapping[score]}\n({model_counts[score]:,})")
        model_values.append(model_counts[score])

colors = ["#2E86AB", "#6FAADB", "#B5D4E3"]
wedges, texts, autotexts = ax4.pie(
    model_values,
    labels=model_labels,
    colors=colors,
    autopct="%1.1f%%",
    startangle=90,
    textprops=dict(fontsize=12),
)
for autotext in autotexts:
    autotext.set_color("white")
    autotext.set_fontweight("bold")
    autotext.set_fontsize(14)

ax4.set_title("Experimental Model Distribution", fontsize=18, fontweight="bold", pad=20)

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig4_model_distribution.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig4_model_distribution.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 5: Top Genetic Interventions (FDR < 0.05)
# ====================
fig5 = plt.figure(figsize=(10, 6))
ax5 = plt.subplot(111)

# Filter for genetic interventions with FDR < 0.05
df_genetic = df_merged[df_merged["condition_category"] == "Genetic"].copy()
df_genetic = df_genetic.sort_values("final_score", ascending=False)
# Remove duplicates based on intervention name
df_genetic = df_genetic.drop_duplicates(subset=["intervention_scores"], keep="first")

top_n = min(20, len(df_genetic))
if top_n > 0:
    top_genetic = df_genetic.head(top_n)
    top_genetic_names = top_genetic["intervention_scores"].values
    # remove any ~ in the names
    top_genetic_names = [name.replace("∼", "") for name in top_genetic_names]

    y_pos = np.arange(len(top_genetic))

    # Color bars based on direction column
    # direction = 1.0: pro-aging (red), direction = -1.0: anti-aging (blue)
    bar_colors = [
        "#D55E00" if direction > 0 else "#0072B2"
        for direction in top_genetic["direction"].values
    ]

    bars = ax5.barh(
        y_pos,
        top_genetic["final_score"].values,
        color=bar_colors,
        alpha=0.8,
        edgecolor="#333333",
        linewidth=1,
    )

    ax5.set_yticks(y_pos)
    ax5.set_yticklabels(top_genetic_names, fontsize=10)
    ax5.set_xlabel("Final Score", fontsize=16, fontweight="bold")
    ax5.set_title(
        f"Top {top_n} Genetic Interventions (FDR < 0.05)",
        fontsize=18,
        fontweight="bold",
        pad=20,
    )
    ax5.invert_yaxis()

    # Add effect size as text (absolute value, with direction indicated by color)
    for i, (idx, row) in enumerate(top_genetic.iterrows()):
        effect = row["effect_size"]
        ax5.text(
            row["final_score"] + 1,
            i,
            f"|β|={abs(effect):.2f}",
            va="center",
            fontsize=9,
            color="#333333",
        )

    # Add legend for color coding
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#D55E00", edgecolor="#333333", label="Pro-aging"),
        Patch(facecolor="#0072B2", edgecolor="#333333", label="Anti-aging"),
    ]
    ax5.legend(
        handles=legend_elements,
        loc="lower right",
        frameon=True,
        fontsize=10,
        bbox_to_anchor=(1.3, 0.02),
        framealpha=0.9,
    )

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig5_top_genetic_interventions.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig5_top_genetic_interventions.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 6: Top Drug Interventions (FDR < 0.05)
# ====================
fig6 = plt.figure(figsize=(10, 6))
ax6 = plt.subplot(111)

# Filter for drug interventions with FDR < 0.05
df_drugs = df_merged[df_merged["condition_category"] == "Drug"].copy()
df_drugs = df_drugs.sort_values("final_score", ascending=False)
# Remove duplicates based on intervention name
df_drugs = df_drugs.drop_duplicates(subset=["intervention_scores"], keep="first")

top_n = min(20, len(df_drugs))
if top_n > 0:
    top_drugs = df_drugs.head(40)
    top_drugs_names = top_drugs["intervention_scores"].values
    # Manual dictionary to exclude non-drug interventions that were misclassified
    non_drug_exclusions = {}

    # Filter out non-drug interventions
    top_drugs = top_drugs[~top_drugs["intervention_scores"].isin(non_drug_exclusions)]
    top_drugs = top_drugs.head(top_n)

    y_pos = np.arange(len(top_drugs))

    # Color bars based on direction column
    # direction = 1.0: pro-aging (red), direction = -1.0: anti-aging (blue)
    bar_colors = [
        "#D55E00" if direction > 0 else "#0072B2"
        for direction in top_drugs["direction"].values
    ]

    bars = ax6.barh(
        y_pos,
        top_drugs["final_score"].values,
        color=bar_colors,
        alpha=0.8,
        edgecolor="#333333",
        linewidth=1,
    )

    ax6.set_yticks(y_pos)
    ax6.set_yticklabels(top_drugs["intervention_scores"].values, fontsize=10)
    ax6.set_xlabel("Final Score", fontsize=16, fontweight="bold")
    ax6.set_title(
        f"Top {top_n} Drug Interventions (FDR < 0.05)",
        fontsize=18,
        fontweight="bold",
        pad=20,
    )
    ax6.invert_yaxis()

    # Add effect size as text (absolute value, with direction indicated by color)
    for i, (idx, row) in enumerate(top_drugs.iterrows()):
        effect = row["effect_size"]
        ax6.text(
            row["final_score"] + 1,
            i,
            f"|β|={abs(effect):.2f}",
            va="center",
            fontsize=9,
            color="#333333",
        )

    # Add legend for color coding
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#D55E00", edgecolor="#333333", label="Pro-aging"),
        Patch(facecolor="#0072B2", edgecolor="#333333", label="Anti-aging"),
    ]
    ax6.legend(
        handles=legend_elements,
        loc="lower right",
        frameon=True,
        fontsize=10,
        bbox_to_anchor=(1.3, 0.02),
        framealpha=0.9,
    )

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig6_top_drug_interventions.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig6_top_drug_interventions.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 6b: Top Environmental Interventions (FDR < 0.05)
# ====================
fig6b = plt.figure(figsize=(10, 6))
ax6b = plt.subplot(111)

# Filter for environmental interventions with FDR < 0.05
df_environment = df_merged[df_merged["condition_category"] == "Environment"].copy()
df_environment = df_environment.sort_values("final_score", ascending=False)
# Remove duplicates based on intervention name
df_environment = df_environment.drop_duplicates(
    subset=["intervention_scores"], keep="first"
)

top_n_env = min(20, len(df_environment))
if top_n_env > 0:
    top_environment = df_environment.head(40)

    # Manual dictionary to exclude non-environmental interventions that were misclassified
    non_environment_exclusions = {}

    # Filter out non-environmental interventions
    top_environment = top_environment[
        ~top_environment["intervention_scores"].isin(non_environment_exclusions)
    ]
    top_environment = top_environment.head(top_n_env)

    y_pos = np.arange(len(top_environment))

    # Color bars based on direction column
    # direction = 1.0: pro-aging (red), direction = -1.0: anti-aging (blue)
    bar_colors = [
        "#D55E00" if direction > 0 else "#0072B2"
        for direction in top_environment["direction"].values
    ]

    bars = ax6b.barh(
        y_pos,
        top_environment["final_score"].values,
        color=bar_colors,
        alpha=0.8,
        edgecolor="#333333",
        linewidth=1,
    )

    ax6b.set_yticks(y_pos)
    ax6b.set_yticklabels(top_environment["intervention_scores"].values, fontsize=10)
    ax6b.set_xlabel("Final Score", fontsize=16, fontweight="bold")
    ax6b.set_title(
        f"Top {top_n_env} Environmental Interventions (FDR < 0.05)",
        fontsize=18,
        fontweight="bold",
        pad=20,
    )
    ax6b.invert_yaxis()

    # Add effect size as text (absolute value, with direction indicated by color)
    for i, (idx, row) in enumerate(top_environment.iterrows()):
        effect = row["effect_size"]
        ax6b.text(
            row["final_score"] + 1,
            i,
            f"|β|={abs(effect):.2f}",
            va="center",
            fontsize=9,
            color="#333333",
        )

    # Add legend for color coding
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#D55E00", edgecolor="#333333", label="Pro-aging"),
        Patch(facecolor="#0072B2", edgecolor="#333333", label="Anti-aging"),
    ]
    ax6b.legend(
        handles=legend_elements,
        loc="lower right",
        frameon=True,
        fontsize=10,
        bbox_to_anchor=(1.3, 0.02),
        framealpha=0.9,
    )

theme.apply_transforms()
plt.tight_layout()
plt.savefig(
    "figures/fig6b_top_environmental_interventions.pdf", bbox_inches="tight", dpi=300
)
plt.savefig(
    "figures/fig6b_top_environmental_interventions.png", dpi=600, bbox_inches="tight"
)
plt.show()

# ====================
# Figure 6c: Top Disease Interventions (FDR < 0.05)
# ====================
fig6c = plt.figure(figsize=(10, 6))
ax6c = plt.subplot(111)

# Filter for disease interventions with FDR < 0.05
df_disease = df_merged[df_merged["condition_category"] == "Disease"].copy()
# IMPORTANT: Only show diseases with pro-aging effect (direction = 1.0)
df_disease = df_disease[df_disease["direction"] > 0].copy()
df_disease = df_disease.sort_values("final_score", ascending=False)
# Remove duplicates based on intervention name
df_disease = df_disease.drop_duplicates(subset=["intervention_scores"], keep="first")

top_n_dis = min(20, len(df_disease))
if top_n_dis > 0:
    top_disease = df_disease.head(40)

    # Manual dictionary to exclude non-disease interventions that were misclassified
    non_disease_exclusions = {}

    # Filter out non-disease interventions
    top_disease = top_disease[
        ~top_disease["intervention_scores"].isin(non_disease_exclusions)
    ]
    top_disease = top_disease.head(top_n_dis)

    y_pos = np.arange(len(top_disease))

    # Color bars based on direction column (all should be red since we filtered for direction > 0)
    # direction = 1.0: pro-aging (red)
    bar_colors = [
        "#D55E00" if direction > 0 else "#0072B2"
        for direction in top_disease["direction"].values
    ]

    bars = ax6c.barh(
        y_pos,
        top_disease["final_score"].values,
        color=bar_colors,
        alpha=0.8,
        edgecolor="#333333",
        linewidth=1,
    )

    ax6c.set_yticks(y_pos)
    ax6c.set_yticklabels(top_disease["intervention_scores"].values, fontsize=10)
    ax6c.set_xlabel("Final Score", fontsize=16, fontweight="bold")
    ax6c.set_title(
        f"Top {top_n_dis} Disease Interventions (FDR < 0.05)",
        fontsize=18,
        fontweight="bold",
        pad=20,
    )
    ax6c.invert_yaxis()

    # Add effect size as text (absolute value, with direction indicated by color)
    for i, (idx, row) in enumerate(top_disease.iterrows()):
        effect = row["effect_size"]
        ax6c.text(
            row["final_score"] + 1,
            i,
            f"|β|={abs(effect):.2f}",
            va="center",
            fontsize=9,
            color="#333333",
        )

    # Add legend for color coding (all should be pro-aging)
    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor="#D55E00", edgecolor="#333333", label="Pro-aging"),
    ]
    ax6c.legend(
        handles=legend_elements,
        loc="lower right",
        frameon=True,
        fontsize=10,
        bbox_to_anchor=(1.3, 0.02),
        framealpha=0.9,
    )

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig6c_top_disease_interventions.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig6c_top_disease_interventions.png", dpi=600, bbox_inches="tight")
plt.show()

import os

# ====================
# Figure 6d: Combined Top Interventions by Category (FDR < 0.05) with Coolwarm Colormap
# ====================
from matplotlib import cm
from matplotlib.colors import Normalize

# ============================================================================
# ANNOTATION MANAGEMENT FOR FIG6D - Use Same Manual Adjustments as Volcano Plot
# ============================================================================
print("\n" + "=" * 80)
print("FIG6D ANNOTATION MANAGEMENT - LOADING VOLCANO PLOT ADJUSTMENTS")
print("=" * 80)

# Use the same annotation file as the volcano plot
volcano_annotation_file = (
    "figures/volcano_plot_annotations_FOR_MANUAL_ADJUSTMENT_result_mm_rna.csv"
)

# Prepare data for all four categories (reuse existing filtered dataframes)
# Genetic
df_genetic_combined = df_merged[df_merged["condition_category"] == "Genetic"].copy()
df_genetic_combined = df_genetic_combined.sort_values("final_score", ascending=False)
df_genetic_combined = df_genetic_combined.drop_duplicates(
    subset=["intervention_scores"], keep="first"
)
top_genetic_combined = df_genetic_combined.head(min(20, len(df_genetic_combined)))

# Drug
df_drugs_combined = df_merged[df_merged["condition_category"] == "Drug"].copy()
df_drugs_combined = df_drugs_combined.sort_values("final_score", ascending=False)
df_drugs_combined = df_drugs_combined.drop_duplicates(
    subset=["intervention_scores"], keep="first"
)
non_drug_exclusions = {}
top_drugs_temp = df_drugs_combined.head(40)
top_drugs_temp = top_drugs_temp[
    ~top_drugs_temp["intervention_scores"].isin(non_drug_exclusions)
]
top_drugs_combined = top_drugs_temp.head(min(20, len(top_drugs_temp)))

# Environment
df_env_combined = df_merged[df_merged["condition_category"] == "Environment"].copy()
df_env_combined = df_env_combined.sort_values("final_score", ascending=False)
df_env_combined = df_env_combined.drop_duplicates(
    subset=["intervention_scores"], keep="first"
)
non_environment_exclusions = {}
top_env_temp = df_env_combined.head(40)
top_env_temp = top_env_temp[
    ~top_env_temp["intervention_scores"].isin(non_environment_exclusions)
]
top_env_combined = top_env_temp.head(min(20, len(top_env_temp)))

# Disease (only pro-aging)
df_dis_combined = df_merged[df_merged["condition_category"] == "Disease"].copy()
df_dis_combined = df_dis_combined[df_dis_combined["direction"] > 0].copy()
df_dis_combined = df_dis_combined.sort_values("final_score", ascending=False)
df_dis_combined = df_dis_combined.drop_duplicates(
    subset=["intervention_scores"], keep="first"
)
non_disease_exclusions = {}
top_dis_temp = df_dis_combined.head(40)
top_dis_temp = top_dis_temp[
    ~top_dis_temp["intervention_scores"].isin(non_disease_exclusions)
]
top_dis_combined = top_dis_temp.head(min(20, len(top_dis_temp)))

# Prepare categories data
categories_data = [
    ("Genetic", top_genetic_combined),
    ("Drug", top_drugs_combined),
    ("Environmental", top_env_combined),
    ("Disease", top_dis_combined),
]

# Check if volcano plot manual adjustment file exists
manual_adjustments = {}
if os.path.exists(volcano_annotation_file):
    print(f"\n✓ Found volcano plot manual adjustment file!")
    print(f"  Loading: {volcano_annotation_file}")

    try:
        volcano_manual_df = pd.read_csv(volcano_annotation_file)
        print(f"  Volcano annotation file has {len(volcano_manual_df)} entries")

        # Create mapping from gse_id to adjusted values
        # The volcano file uses gse_id, intervention, and adjusted_label
        for _, row in volcano_manual_df.iterrows():
            gse_id = row.get("gse_id", "")
            intervention = row.get("intervention", row.get("full_intervention", ""))

            # Store by both gse_id and intervention for flexible matching
            key = gse_id
            manual_adjustments[key] = {
                "adjusted_label": row.get("adjusted_label", intervention),
                "intervention": intervention,
            }

        print(
            f"  ✓ Loaded {len(manual_adjustments)} manual adjustments from volcano plot"
        )

        # Apply adjustments to dataframes based on dataset matching
        adjustment_count = 0
        for cat_name, df_cat in categories_data:
            if len(df_cat) > 0:
                # Create adjusted label column by matching gse_id or intervention
                def get_adjusted_label(row):
                    # Try matching by case_id (which contains gse_id)
                    case_id = row.get("case_id", row.get("dataset", ""))
                    gse_id = row.get("gse_id", "")
                    intervention = row.get("intervention_scores", "")

                    # Check if we have manual adjustment for this gse_id
                    if gse_id in manual_adjustments:
                        return manual_adjustments[gse_id]["adjusted_label"]

                    # Check by intervention name match
                    for key, adj in manual_adjustments.items():
                        if adj["intervention"] == intervention:
                            return adj["adjusted_label"]

                    # No match found, use original
                    return intervention

                df_cat["adjusted_intervention_label"] = df_cat.apply(
                    get_adjusted_label, axis=1
                )

                # Count how many were adjusted
                adjusted_mask = (
                    df_cat["adjusted_intervention_label"]
                    != df_cat["intervention_scores"]
                )
                adjustment_count += adjusted_mask.sum()

        if adjustment_count > 0:
            print(
                f"  ✓ Applied {adjustment_count} label adjustments from volcano plot to fig6d"
            )
            print(f"  → These adjusted labels will be used in the fig6d plot")
        else:
            print(
                f"  Note: No matching interventions found between volcano plot and fig6d"
            )

    except Exception as e:
        print(f"\n⚠ Error loading volcano plot adjustments: {e}")
        print("  Using original labels...")
        manual_adjustments = {}

else:
    print(f"\n📝 No volcano plot manual adjustment file found.")
    print(f"   File: {volcano_annotation_file}")
    print(
        f"   Run the volcano plot script (s01_visualize_results.py) first to create manual adjustments."
    )
    print(f"   Using original intervention labels for fig6d...")
    print(f"\nNote: Fig6d shares the same manual adjustment file as the volcano plot.")
    print(
        f"      Edit labels in the volcano plot annotation file to update both figures."
    )

# Ensure all dataframes have adjusted_intervention_label column (use original if no adjustments)
for cat_name, df_cat in categories_data:
    if len(df_cat) > 0 and "adjusted_intervention_label" not in df_cat.columns:
        df_cat["adjusted_intervention_label"] = df_cat["intervention_scores"]

print("=" * 80 + "\n")

# Create combined figure
fig6d = plt.figure(figsize=(32, 6))

# Find global min/max for consistent color scaling across all categories
all_signed_effects = []
for cat_name, df_cat in categories_data:
    if len(df_cat) > 0:
        signed_effects = (df_cat["direction"] * np.abs(df_cat["effect_size"])).values
        all_signed_effects.extend(signed_effects)

if len(all_signed_effects) > 0:
    vmin, vmax = min(all_signed_effects), max(all_signed_effects)
    # Make symmetric around 0 for better visualization
    vmax_abs = max(abs(vmin), abs(vmax))
    vmin, vmax = -vmax_abs, vmax_abs

    # Create colormap - RdBu_r provides better contrast for aging effects
    # Red for pro-aging (positive), Blue for anti-aging (negative)
    cmap = cm.get_cmap("RdBu_r")
    norm = Normalize(vmin=vmin, vmax=vmax)

    # Plot each category
    for idx, (cat_name, df_cat) in enumerate(categories_data):
        ax = plt.subplot(1, 4, idx + 1)

        if len(df_cat) > 0:
            # Use adjusted labels if available, otherwise use original
            if "adjusted_intervention_label" in df_cat.columns:
                cat_names = df_cat["adjusted_intervention_label"].values
            else:
                cat_names = df_cat["intervention_scores"].values
            cat_names = [name.replace("∼", "") for name in cat_names]
            y_pos = np.arange(len(df_cat))

            # Calculate signed effect for coloring
            signed_effects = (df_cat["direction"] * abs(df_cat["effect_size"])).values
            colors = [cmap(norm(se)) for se in signed_effects]

            bars = ax.barh(
                y_pos,
                df_cat["final_score"].values,
                color=colors,
                edgecolor="#333333",
                linewidth=0.8,
            )

            ax.set_yticks(y_pos)
            ax.set_yticklabels(cat_names, fontsize=9)
            ax.set_xlabel("Final Score", fontsize=12, fontweight="bold")
            ax.set_title(
                f"{cat_name}\n(n={len(df_cat)})", fontsize=13, fontweight="bold", pad=10
            )
            ax.invert_yaxis()
            ax.grid(axis="x", alpha=0.3, linestyle="--")

    # Add shared colorbar
    cbar_ax = fig6d.add_axes([0.92, 0.15, 0.01, 0.7])
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=cbar_ax)
    cbar.set_label(
        "Aging Effect\n← Pro-aging | Anti-aging →",
        fontsize=12,
        fontweight="bold",
        rotation=270,
        labelpad=25,
    )

theme.apply_transforms()
plt.tight_layout(rect=[0, 0, 0.91, 1])
plt.savefig("figures/fig6d_combined_interventions.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig6d_combined_interventions.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 7: Score vs Model Quality (Improved)
# ====================
fig7 = plt.figure(figsize=(7, 5))
ax7 = plt.subplot(111)

# Create categorical x-axis
model_categories = ["Cell lines", "Primary cells", "In vivo"]
model_scores_map = {0: 0, 15: 1, 25: 2}

# Map scores to categories
x_categorical = [
    model_scores_map.get(score, -1) for score in df_scores["q1_experimental_model"]
]
x_jittered = np.array(x_categorical) + np.random.normal(0, 0.15, len(x_categorical))
y = df_scores["final_score"].values

# Calculate density for coloring
valid_indices = [i for i, x in enumerate(x_categorical) if x != -1]
x_valid = np.array([x_jittered[i] for i in valid_indices])
y_valid = np.array([y[i] for i in valid_indices])

if len(x_valid) > 0:
    xy = np.vstack([x_valid, y_valid])
    z = gaussian_kde(xy)(xy)

    # Plot with outline
    ax7.scatter(x_valid, y_valid, s=30 * 1.5, c="black", alpha=0.2, zorder=1)
    scatter = ax7.scatter(
        x_valid, y_valid, s=30, c=z, cmap="viridis", alpha=0.6, zorder=2
    )

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax7)
    cbar.set_label("Density", fontsize=12)

ax7.set_xticks(range(3))
ax7.set_xticklabels(model_categories, fontsize=12)
ax7.set_xlabel("Experimental Model Quality", fontsize=16, fontweight="bold")
ax7.set_ylabel("Final Score", fontsize=16, fontweight="bold")
ax7.set_title("Model Quality vs Final Score", fontsize=18, fontweight="bold", pad=20)
ax7.grid(axis="y", alpha=0.3, linestyle="--")

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig7_score_vs_model_quality.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig7_score_vs_model_quality.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 8: Aging Pathway Relevance
# ====================
fig8 = plt.figure(figsize=(6, 4))
ax8 = plt.subplot(111)

pathway_mapping = {
    15: "Multiple pathways",
    10: "Single pathway",
    5: "Indirect effect",
    0: "No connection",
}
pathway_counts = (
    df_scores["q2_1_aging_pathways"].value_counts().sort_index(ascending=False)
)
pathway_data = []
pathway_labels = []

for score, label in pathway_mapping.items():
    if score in pathway_counts.index:
        count = pathway_counts[score]
        pct = count / len(df_scores) * 100
        pathway_data.append(count)
        pathway_labels.append(f"{label}\n({score} pts)")

bars = ax8.barh(
    range(len(pathway_data)),
    pathway_data,
    color="#F18F01",
    alpha=0.8,
    edgecolor="black",
    linewidth=1.5,
)
ax8.set_yticks(range(len(pathway_labels)))
ax8.set_yticklabels(pathway_labels, fontsize=12)
ax8.set_xlabel("Number of Interventions", fontsize=16, fontweight="bold")
ax8.set_title(
    "Aging Pathway Relevance Distribution", fontsize=18, fontweight="bold", pad=20
)
ax8.invert_yaxis()

# Add count labels
for i, (bar, count) in enumerate(zip(bars, pathway_data)):
    pct = count / len(df_scores) * 100
    ax8.text(
        count + 50,
        i,
        f"{count:,} ({pct:.1f}%)",
        va="center",
        fontsize=11,
        fontweight="bold",
    )

theme.apply_transforms()
plt.tight_layout()
plt.savefig("figures/fig8_pathway_relevance.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig8_pathway_relevance.png", dpi=600, bbox_inches="tight")
plt.show()

# ====================
# Figure 9: Score Component Correlations
# ====================
fig9 = plt.figure(figsize=(10, 9))
ax9 = plt.subplot(111)

correlation_cols = ["q0_plausibility"] + score_components + ["final_score"]
corr_matrix = df_scores[correlation_cols].corr()

# Rename columns for display
display_names = {
    "q0_plausibility": "Plausibility",
    "q1_experimental_model": "Model",
    "q2_1_aging_pathways": "Pathways",
    "q2_2_age_diseases": "Diseases",
    "q3_1_experimental_rigor": "Rigor",
    "q4_1_translational_potential": "Translational",
    "q5_1_literature_saturation": "Literature",
    "q5_2_mechanistic_novelty": "Novelty",
    "final_score": "Final Score",
}
corr_matrix.rename(index=display_names, columns=display_names, inplace=True)

sns.heatmap(
    corr_matrix,
    annot=True,
    fmt=".2f",
    cmap="coolwarm",
    center=0,
    square=True,
    ax=ax9,
    cbar_kws={"shrink": 0.8},
    vmin=-0.5,
    vmax=1,
    linewidths=0.5,
    linecolor="gray",
)
ax9.set_title("Score Component Correlations", fontsize=18, fontweight="bold", pad=20)

theme.apply_transforms()
plt.setp(ax9.get_xticklabels(), rotation=45, ha="right")
plt.setp(ax9.get_yticklabels(), rotation=0)
plt.tight_layout()
plt.savefig("figures/fig9_correlations.pdf", bbox_inches="tight", dpi=300)
plt.savefig("figures/fig9_correlations.png", dpi=600, bbox_inches="tight")
plt.show()

# Print summary statistics
print("\n" + "=" * 60)
print("SUMMARY STATISTICS")
print("=" * 60)
print(f"Total interventions scored: {len(df_scores):,}")
print(f"FDR < 0.05 interventions: {len(df_merged):,}")
print(
    f"High-scoring (≥100): {len(df_scores[df_scores['final_score'] >= 100]):,} ({len(df_scores[df_scores['final_score'] >= 100])/len(df_scores)*100:.1f}%)"
)
print(
    f"Very high-scoring (≥120): {len(df_scores[df_scores['final_score'] >= 120]):,} ({len(df_scores[df_scores['final_score'] >= 120])/len(df_scores)*100:.1f}%)"
)
print(f"\nIntervention types (FDR < 0.05):")
print(f"  Genetic: {len(df_genetic):,}")
print(f"  Drug: {len(df_drugs):,}")
print(f"  Environmental: {len(df_environment):,}")
print(f"  Disease: {len(df_disease):,}")
print(f"  Other: {len(df_merged[df_merged['condition_category'] == 'Other']):,}")
print(f"\nMedian score: {df_scores['final_score'].median():.0f}")
print(
    f"Score range: {df_scores['final_score'].min()} to {df_scores['final_score'].max()}"
)
print(f"\nPlausibility filter:")
print(
    f"  Penalties (-50): {len(df_scores[df_scores['q0_plausibility'] == -50]):,} ({len(df_scores[df_scores['q0_plausibility'] == -50])/len(df_scores)*100:.1f}%)"
)
print(
    f"  Bonuses (+30): {len(df_scores[df_scores['q0_plausibility'] == 30]):,} ({len(df_scores[df_scores['q0_plausibility'] == 30])/len(df_scores)*100:.1f}%)"
)

# Find ouabain
ouabain_entries = df_scores[
    df_scores["intervention"].str.contains("ouabain", case=False, na=False)
]
if len(ouabain_entries) > 0:
    print(f"\nOuabain score: {ouabain_entries.iloc[0]['final_score']}")

print("\nTop 5 overall interventions:")
top5 = df_scores.nlargest(5, "final_score")[["intervention", "final_score"]]
for idx, row in top5.iterrows():
    print(f"  {row['intervention']}: {row['final_score']}")

print("\nAll figures saved in 'figures/' directory")
print("=" * 60)
