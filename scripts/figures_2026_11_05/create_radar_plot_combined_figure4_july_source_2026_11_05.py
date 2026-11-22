#!/usr/bin/env python3
"""
Combined-modality radar plot showing intervention categories pooled across all
three modalities (mouse RNA-seq + human RNA-seq + human DNAm). Outputs a single
radar with one polygon for "All interventions" and one polygon for
"Significant only", with category counts pooled across modalities.

Outputs:
  figures/create_radar_plot_combined.png
  figures/create_radar_plot_combined.pdf
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# Load combined 3-modality plot data
df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)

# Analytic set: exclude extreme-tail-flagged rows (mirrors analysis_2/3/4_combined)
if "outlier_status" in df.columns:
    df = df[df["outlier_status"] == "Normal"].copy()

# Canonicalise condition_category to the canonical 5 buckets (others become Other)
CANONICAL = {"Disease", "Drug", "Environment", "Genetic", "Other"}
df["condition_category"] = df["condition_category"].where(
    df["condition_category"].isin(CANONICAL), "Other"
)

# Order kept consistent with create_radar_plot.py
categories = ["Disease", "Drug", "Environment", "Genetic", "Other"]

# Pool counts across all modalities (no stratification)
# "Significant" is FDR<0.05 (matches the figure caption); the is_significant
# column additionally applies the score floor and is not used here.
all_counts = df["condition_category"].value_counts()
sig_df = df[df["fdr"] < 0.05]
sig_counts = sig_df["condition_category"].value_counts()

data_all = [all_counts.get(cat, 0) for cat in categories]
data_sig = [sig_counts.get(cat, 0) for cat in categories]

# Style
BG_WHITE = "#ffffff"
BLUE = "#2a475e"
GREY70 = "#b3b3b3"
COLOR_ALL = "#B8B8B8"
COLOR_SIG = "#FF6B6B"

VARIABLES = categories
VARIABLES_N = len(VARIABLES)
ANGLES = [n / VARIABLES_N * 2 * np.pi for n in range(VARIABLES_N)]
ANGLES += ANGLES[:1]

fig = plt.figure(figsize=(8.5, 7))
ax = fig.add_subplot(111, polar=True)
fig.subplots_adjust(left=0.05, right=0.78)
fig.patch.set_facecolor(BG_WHITE)
ax.set_facecolor(BG_WHITE)

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# Hard-cap radial axis at exactly 35k.
RMAX = 35000
ax.set_ylim(-RMAX * 0.1, RMAX)

data_all_plot = data_all + [data_all[0]]
data_sig_plot = data_sig + [data_sig[0]]

ax.plot(
    ANGLES, data_all_plot, c=COLOR_ALL, linewidth=3,
    label="All Interventions", linestyle="--", alpha=0.7,
)
ax.scatter(ANGLES, data_all_plot, s=120, c=COLOR_ALL, zorder=10, alpha=0.7)

ax.plot(ANGLES, data_sig_plot, c=COLOR_SIG, linewidth=3, label="Significant Only")
ax.scatter(ANGLES, data_sig_plot, s=120, c=COLOR_SIG, zorder=10)

ax.fill(ANGLES, data_all_plot, COLOR_ALL, alpha=0.1)
ax.fill(ANGLES, data_sig_plot, COLOR_SIG, alpha=0.15)

ax.set_xticks(ANGLES[:-1])
ax.set_xticklabels(VARIABLES, size=14, fontweight="bold")

ax.set_yticks([5000, 10000, 15000, 20000, 25000, 30000, 35000])
ax.set_yticklabels(["5k", "10k", "15k", "20k", "25k", "30k", "35k"], size=10, color=GREY70)
ax.yaxis.grid(True, linestyle="--", alpha=0.3)
ax.xaxis.grid(True, linestyle="--", alpha=0.3)

ax.spines["start"].set_color("none")
ax.spines["polar"].set_color("none")

handles = [
    Line2D([], [], c=COLOR_SIG, lw=3, marker="o", markersize=8,
           label=f"Significant (n={sum(data_sig):,})"),
    Line2D([], [], c=COLOR_ALL, lw=3, marker="o", markersize=8,
           label=f"All (n={sum(data_all):,})", linestyle="--", alpha=0.7),
]

# Value labels (category counts)
for idx, (angle, cat) in enumerate(zip(ANGLES[:-1], categories)):
    all_val = data_all[idx]
    sig_val = data_sig[idx]
    offset = RMAX * 0.08
    ax.text(angle, all_val + offset, f"{all_val:,}",
            size=10, ha="center", va="bottom",
            color=COLOR_ALL, fontweight="bold")
    ax.text(angle, sig_val + offset, f"{sig_val:,}",
            size=10, ha="center", va="bottom",
            color=COLOR_SIG, fontweight="bold")

ax.legend(handles=handles, loc="upper left",
          bbox_to_anchor=(1.02, 1.0), labelspacing=1.0,
          frameon=False, fontsize=11)

# Title (combined-modality)
fig.suptitle(
    "Intervention Categories (all modalities combined)",
    x=0.5, y=0.98, ha="center",
    fontsize=18, color=BLUE, weight="bold",
)

# Subtitle: percentage significant
pct_sig = (sum(data_sig) / sum(data_all) * 100) if sum(data_all) else 0
ax.text(
    0.5, -0.15,
    f"{pct_sig:.1f}% of interventions are statistically significant "
    f"(mouse RNA-seq + human RNA-seq + human DNAm pooled)",
    transform=ax.transAxes, ha="center",
    fontsize=11, style="italic", color=GREY70,
)

plt.tight_layout()
plt.savefig("figures/create_radar_plot_combined.pdf", dpi=300, bbox_inches="tight")
plt.savefig("figures/create_radar_plot_combined.png", dpi=600, bbox_inches="tight")
print("Saved figures/create_radar_plot_combined.{pdf,png}")

print(f"\nTotal interventions: {sum(data_all):,}")
print(f"Significant interventions: {sum(data_sig):,}")
print("\nBreakdown by category (pooled):")
print("-" * 40)
for cat, all_count, sig_count in zip(categories, data_all, data_sig):
    pct = (sig_count / all_count * 100) if all_count > 0 else 0
    print(f"{cat:12} - All: {all_count:6,}  Sig: {sig_count:6,} ({pct:.1f}%)")
