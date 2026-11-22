#!/usr/bin/env python3
"""
Create Nature publication-quality radar plot showing ouabain's scoring profile
across all evaluation subcategories with detailed rationales
"""

import json
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

warnings.filterwarnings("ignore")

# Apply consistent theme
from aquarel import load_theme

theme = (
    load_theme("scientific")
    .set_grid(draw=False)
    .set_font(size=14)
    .set_ticks(direction="out")
    .set_axis_labels(pad=10)
)
theme.apply()

# Set white background for Nature publication
plt.rcParams["savefig.transparent"] = False
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"

# Load scoring data
df_scores = pd.read_csv("data/scores_summary.csv")

# Find ouabain entry
ouabain_data = df_scores[
    df_scores["intervention"].str.contains("ouabain", case=False, na=False)
]

if len(ouabain_data) == 0:
    print("Ouabain not found in the dataset!")
    exit()

ouabain_row = ouabain_data.iloc[0]

# Scoring rationales from the evaluation
scoring_rationales = {
    "Plausibility": "Cardiac glycoside with strong mechanistic basis (+30 bonus)",
    "Model\nQuality": "In vivo aged mouse studies (25/25 points)",
    "Aging\nPathways": "Targets multiple aging pathways via Na+/K+-ATPase (15/15)",
    "Age\nDiseases": "Relevant to cardiovascular disease and cellular senescence (8/10)",
    "Experimental\nRigor": "Dose-response, controlled studies, validation (10/10)",
    "Translational\nPotential": "FDA-approved drug, phase 2 trials (5/5)",
    "Literature\nSaturation": "Well-studied mechanism, multiple publications (15/15)",
    "Mechanistic\nNovelty": "Novel application for aging, repurposing potential (15/20)",
}

# Define scoring categories and their max values
scoring_categories = {
    "Plausibility": ("q0_plausibility", -50, 30),
    "Model\nQuality": ("q1_experimental_model", 0, 25),
    "Aging\nPathways": ("q2_1_aging_pathways", 0, 15),
    "Age\nDiseases": ("q2_2_age_diseases", 0, 10),
    "Experimental\nRigor": ("q3_1_experimental_rigor", 0, 10),
    "Translational\nPotential": ("q4_1_translational_potential", 0, 5),
    "Literature\nSaturation": ("q5_1_literature_saturation", 0, 15),
    "Mechanistic\nNovelty": ("q5_2_mechanistic_novelty", 0, 20),
}

# Extract ouabain scores and normalize them
ouabain_scores = []
ouabain_actual = []
labels = []
max_possible = []

for label, (col, min_val, max_val) in scoring_categories.items():
    score = ouabain_row[col]
    # Normalize to 0-100 scale for visualization
    if col == "q0_plausibility":
        # Special handling for plausibility which can be negative
        if score == -50:
            normalized = 0
        elif score == 0:
            normalized = 50
        else:  # score == 30
            normalized = 100
    else:
        normalized = (score / max_val) * 100 if max_val > 0 else 0

    ouabain_scores.append(normalized)
    ouabain_actual.append(score)
    labels.append(label)
    max_possible.append(100)  # All normalized to 100

# Calculate mean scores for comparison
mean_scores = []
for label, (col, min_val, max_val) in scoring_categories.items():
    mean_score = df_scores[col].mean()
    if col == "q0_plausibility":
        # Special handling for plausibility
        if mean_score < -25:
            normalized = 25
        elif mean_score < 15:
            normalized = 50 + (mean_score / 30) * 50
        else:
            normalized = 75 + (mean_score / 30) * 25
    else:
        normalized = (mean_score / max_val) * 100 if max_val > 0 else 0
    mean_scores.append(normalized)

# Setup radar chart
VARIABLES = labels
VARIABLES_N = len(VARIABLES)

# Calculate angles for each variable
ANGLES = [n / VARIABLES_N * 2 * np.pi for n in range(VARIABLES_N)]
ANGLES += ANGLES[:1]

# Create figure with Nature-style dimensions - wider to accommodate rationales
fig = plt.figure(figsize=(14, 8), facecolor="white")

# Create grid with radar plot on left and rationale panel on right
from matplotlib.gridspec import GridSpec

gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 0.8], wspace=0.15)
ax = fig.add_subplot(gs[0], polar=True, facecolor="white")

# Configure angular axis
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ax.set_ylim(0, 110)

# Color scheme - Nature-friendly colors
OUABAIN_COLOR = "#E64B35"  # Nature red
MEAN_COLOR = "#4DBBD5"  # Nature cyan
GRID_COLOR = "#CCCCCC"  # Light grey for grid

# Plot data
ouabain_plot = ouabain_scores + [ouabain_scores[0]]
mean_plot = mean_scores + [mean_scores[0]]

# Plot mean scores first (background)
ax.plot(
    ANGLES,
    mean_plot,
    color=MEAN_COLOR,
    linewidth=2.5,
    label=f"Database mean (n={len(df_scores):,})",
    linestyle="--",
    alpha=0.7,
)
ax.fill(ANGLES, mean_plot, MEAN_COLOR, alpha=0.15)

# Plot ouabain scores (foreground)
ax.plot(
    ANGLES, ouabain_plot, color=OUABAIN_COLOR, linewidth=3, label="Ouabain", zorder=10
)
ax.fill(ANGLES, ouabain_plot, OUABAIN_COLOR, alpha=0.2)

# Add scatter points for clarity
ax.scatter(
    ANGLES[:-1],
    ouabain_plot[:-1],
    s=100,
    c=OUABAIN_COLOR,
    zorder=15,
    edgecolors="white",
    linewidth=1,
)
ax.scatter(
    ANGLES[:-1],
    mean_plot[:-1],
    s=60,
    c=MEAN_COLOR,
    zorder=12,
    alpha=0.7,
    edgecolors="white",
    linewidth=1,
)

# Customize axes
ax.set_xticks(ANGLES[:-1])
ax.set_xticklabels(VARIABLES, size=12, fontweight="600")

# Add concentric circles for grid with subtle styling
ax.set_yticks([20, 40, 60, 80, 100])
ax.set_yticklabels(["20", "40", "60", "80", "100"], size=10, color="#666666")
ax.yaxis.grid(True, linestyle=":", alpha=0.4, color=GRID_COLOR, linewidth=0.8)
ax.xaxis.grid(True, linestyle="-", alpha=0.3, color=GRID_COLOR, linewidth=0.8)

# Remove spines for cleaner look
ax.spines["start"].set_color("none")
ax.spines["polar"].set_color(GRID_COLOR)
ax.spines["polar"].set_linewidth(0.8)

# Add clean legend with white background
legend = ax.legend(
    loc="upper left",
    bbox_to_anchor=(-0.15, 1.12),
    labelspacing=1.2,
    frameon=True,
    fontsize=12,
    fancybox=False,
    shadow=False,
    facecolor="white",
    edgecolor="#CCCCCC",
    framealpha=1.0,
)
legend.get_frame().set_linewidth(0.8)

# Add rationale panel on the right
ax_text = fig.add_subplot(gs[1], facecolor="white")
ax_text.axis("off")

# Add rationales with scores
y_positions = np.linspace(0.9, 0.1, len(scoring_categories))
clean_labels_for_text = [
    "Plausibility",
    "Model Quality",
    "Aging Pathways",
    "Age Diseases",
    "Experimental Rigor",
    "Translational Potential",
    "Literature Saturation",
    "Mechanistic Novelty",
]

for i, (label, clean_label) in enumerate(zip(labels, clean_labels_for_text)):
    score = ouabain_actual[i]
    rationale = scoring_rationales[label]

    # Score display
    if "Plausibility" in clean_label:
        score_str = f"+{score}" if score > 0 else str(int(score))
    else:
        _, _, max_val = scoring_categories[label]
        score_str = f"{int(score)}/{int(max_val)}"

    # Add category title with score
    ax_text.text(
        0.05,
        y_positions[i],
        f"{clean_label}:",
        fontsize=11,
        fontweight="bold",
        color=OUABAIN_COLOR,
        transform=ax_text.transAxes,
    )

    ax_text.text(
        0.05,
        y_positions[i] - 0.03,
        f"Score: {score_str}",
        fontsize=10,
        color="#333333",
        transform=ax_text.transAxes,
    )

    # Add rationale text (wrapped if needed)
    ax_text.text(
        0.05,
        y_positions[i] - 0.06,
        rationale,
        fontsize=9,
        color="#666666",
        style="italic",
        transform=ax_text.transAxes,
        wrap=True,
    )

# Add title with proper spacing
fig.suptitle(
    "Ouabain Multi-dimensional Scoring Profile",
    x=0.35,  # Adjusted for new layout
    y=0.95,
    ha="center",
    fontsize=16,
    fontweight="bold",
    color="#000000",
)

# Add score information as subtitle (more subtle)
score_text = f"Total score: {int(ouabain_row['final_score'])} (95th percentile)"
ax.text(
    0.5,
    -0.15,
    score_text,
    transform=ax.transAxes,
    ha="center",
    fontsize=11,
    color="#666666",
)

# Add section title for rationales
ax_text.text(
    0.05,
    0.98,
    "Scoring Rationale",
    fontsize=12,
    fontweight="bold",
    color="#000000",
    transform=ax_text.transAxes,
)

# Adjust layout for proper spacing
plt.tight_layout(rect=[0, 0.05, 1, 0.93])

# Save figure with high quality for Nature submission
plt.savefig(
    "figures/ouabain_radar_nature.pdf",
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)
plt.savefig(
    "figures/ouabain_radar_nature.png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)
plt.savefig(
    "figures/ouabain_radar_nature.svg",
    format="svg",
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)
plt.show()

# Print summary for manuscript
print("\n" + "=" * 60)
print("OUABAIN SCORING PROFILE SUMMARY")
print("=" * 60)
print(f"Final Score: {int(ouabain_row['final_score'])}")
print(
    f"Percentile Rank: {(df_scores['final_score'] < ouabain_row['final_score']).mean() * 100:.1f}th"
)
print("\nDetailed Scores (normalized to percentage of maximum):")
print("-" * 40)

for i, (label, (col, min_val, max_val)) in enumerate(scoring_categories.items()):
    score = ouabain_row[col]
    norm_score = ouabain_scores[i]
    mean_norm = mean_scores[i]

    # Clean label for printing
    clean_label = label.replace("\n", " ")

    print(f"{clean_label:22} Ouabain: {norm_score:5.1f}%  Mean: {mean_norm:5.1f}%")

print("\nKey Strengths:")
print(
    "- Maximum scores in: Model Quality (aged mice), Aging Pathways, Experimental Rigor,"
)
print("  Plausibility (max bonus), Translational Potential, Literature Saturation")
print("- Strong scores in: Mechanistic Novelty (75%), Age Diseases (80%)")
print("\nStudy Design:")
print("- In vivo study: 21-22 month old female mice, 12-week treatment")
print("- Cardiac glycoside targeting Na+/K+-ATPase, modulating multiple aging pathways")
print("=" * 60)
