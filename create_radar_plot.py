#!/usr/bin/env python3
"""
Create Nature publication-quality radar plot showing intervention categories
for both all and significant interventions
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

# Load data
df = pd.read_csv("figures/combined_volcano_plot_data.csv")

# Get counts for both all and significant interventions
categories = ["Disease", "Drug", "Environment", "Genetic", "Other"]

# Count all interventions
all_counts = df["condition_category"].value_counts()

# Count significant interventions
sig_df = df[df["is_significant"] == True]
sig_counts = sig_df["condition_category"].value_counts()

# Prepare data for radar plot
data_all = []
data_sig = []
for cat in categories:
    data_all.append(all_counts.get(cat, 0))
    data_sig.append(sig_counts.get(cat, 0))

# Define style parameters
BG_WHITE = "#ffffff"
BLUE = "#2a475e"
GREY70 = "#b3b3b3"
GREY_LIGHT = "#f2efe8"
COLOR_ALL = "#B8B8B8"  # Grey for all interventions
COLOR_SIG = "#FF6B6B"  # Red for significant interventions

# Setup radar chart parameters
VARIABLES = categories
VARIABLES_N = len(VARIABLES)

# Calculate angles for each variable
ANGLES = [n / VARIABLES_N * 2 * np.pi for n in range(VARIABLES_N)]
ANGLES += ANGLES[:1]

# Create figure with actual counts
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111, polar=True)

# Set background
fig.patch.set_facecolor(BG_WHITE)
ax.set_facecolor(BG_WHITE)

# Configure angular axis
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# Cap radial axis at 10k for readability (large values will overflow visually but tick scale stays comparable)
RMAX = 10000
ax.set_ylim(-RMAX * 0.1, RMAX * 1.05)

# Plot actual values
data_all_plot = data_all + [data_all[0]]
data_sig_plot = data_sig + [data_sig[0]]

ax.plot(
    ANGLES,
    data_all_plot,
    c=COLOR_ALL,
    linewidth=3,
    label="All Interventions",
    linestyle="--",
    alpha=0.7,
)
ax.scatter(ANGLES, data_all_plot, s=120, c=COLOR_ALL, zorder=10, alpha=0.7)

ax.plot(ANGLES, data_sig_plot, c=COLOR_SIG, linewidth=3, label="Significant Only")
ax.scatter(ANGLES, data_sig_plot, s=120, c=COLOR_SIG, zorder=10)

# Fill areas
ax.fill(ANGLES, data_all_plot, COLOR_ALL, alpha=0.1)
ax.fill(ANGLES, data_sig_plot, COLOR_SIG, alpha=0.15)

# Customize
ax.set_xticks(ANGLES[:-1])
ax.set_xticklabels(VARIABLES, size=14, fontweight="bold")

# Add grid
ax.set_yticks([2000, 4000, 6000, 8000, 10000])
ax.set_yticklabels(["2k", "4k", "6k", "8k", "10k"], size=10, color=GREY70)
ax.yaxis.grid(True, linestyle="--", alpha=0.3)
ax.xaxis.grid(True, linestyle="--", alpha=0.3)

# Remove spines
ax.spines["start"].set_color("none")
ax.spines["polar"].set_color("none")

# Add legend
handles = [
    Line2D(
        [],
        [],
        c=COLOR_SIG,
        lw=3,
        marker="o",
        markersize=8,
        label=f"Significant (n={sum(data_sig):,})",
    ),
    Line2D(
        [],
        [],
        c=COLOR_ALL,
        lw=3,
        marker="o",
        markersize=8,
        label=f"All (n={sum(data_all):,})",
        linestyle="--",
        alpha=0.7,
    ),
]

# Add value labels
for idx, (angle, cat) in enumerate(zip(ANGLES[:-1], categories)):
    all_val = data_all[idx]
    sig_val = data_sig[idx]

    # Add labels above points (offset based on capped radial scale)
    offset = RMAX * 0.08
    ax.text(
        angle,
        all_val + offset,
        f"{all_val:,}",
        size=10,
        ha="center",
        va="bottom",
        color=COLOR_ALL,
        fontweight="bold",
    )
    ax.text(
        angle,
        sig_val + offset,
        f"{sig_val:,}",
        size=10,
        ha="center",
        va="bottom",
        color=COLOR_SIG,
        fontweight="bold",
    )

legend = ax.legend(
    handles=handles,
    loc="upper right",
    bbox_to_anchor=(1.15, 1.1),
    labelspacing=1.5,
    frameon=False,
    fontsize=12,
)

# Add title
fig.suptitle(
    "Intervention Categories",
    x=0.5,
    y=0.98,
    ha="center",
    fontsize=20,
    color=BLUE,
    weight="bold",
)

# Add subtitle with percentage significant
pct_sig = sum(data_sig) / sum(data_all) * 100
ax.text(
    0.5,
    -0.15,
    f"{pct_sig:.1f}% of all interventions are statistically significant",
    transform=ax.transAxes,
    ha="center",
    fontsize=12,
    style="italic",
    color=GREY70,
)

# Save figure
plt.tight_layout()
plt.savefig("figures/radar_plot_interventions.pdf", dpi=300, bbox_inches="tight")
plt.savefig("figures/radar_plot_interventions.png", dpi=600, bbox_inches="tight")
plt.show()

print("\nRadar plot created successfully!")
print(f"\nTotal interventions: {sum(data_all):,}")
print(f"Significant interventions: {sum(data_sig):,}")
print(f"\nBreakdown by category:")
print("-" * 40)
for cat, all_count, sig_count in zip(categories, data_all, data_sig):
    pct = (sig_count / all_count * 100) if all_count > 0 else 0
    print(f"{cat:12} - All: {all_count:6,}  Sig: {sig_count:6,} ({pct:.1f}%)")
