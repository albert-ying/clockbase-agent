"""
Visualize validation results from both evaluators.
Creates publication-quality figures for manuscript.
"""

import sys
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Set style for publication
plt.style.use("seaborn-v0_8-whitegrid")
sns.set_palette("Set2")
plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 10
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]

# Load data from both evaluators
eval1 = pd.read_csv(
    paths.REPO / "validation_results_100_traces.csv"
)
eval2 = pd.read_csv(
    paths.REPO / "validation_results_100_traces_evaluator2.csv"
)

# Criteria columns
criteria_cols = [
    "hypothesis",
    "statistical_test",
    "p_value_validity",
    "effect_size",
    "comparison",
    "completeness",
    "data_integration",
    "protocol",
]

criteria_labels = [
    "Hypothesis",
    "Statistical\nTest",
    "P-value\nValidity",
    "Effect\nSize",
    "Comparison\nValidity",
    "Completeness",
    "Data\nIntegration",
    "Protocol",
]


def calculate_metrics(df, criteria_cols):
    """Calculate correctness metrics for each criterion."""
    metrics = []
    for col in criteria_cols:
        # Convert column to numeric, treating non-numeric values as NaN
        col_numeric = pd.to_numeric(df[col], errors="coerce")

        # 1 = correct, 0 = data issue, anything else (including NaN from text) = error
        correct = (col_numeric == 1).sum()
        data_issues = (col_numeric == 0).sum()
        errors = len(df) - correct - data_issues  # Everything that's not 1 or 0
        total = len(df)

        # Correctness rate excluding data issues
        analyzable = total - data_issues
        correctness = (correct / analyzable * 100) if analyzable > 0 else 0

        metrics.append(
            {
                "criterion": col,
                "correct": correct,
                "errors": errors,
                "data_issues": data_issues,
                "total": total,
                "correctness_rate": correctness,
            }
        )

    return pd.DataFrame(metrics)


# Calculate metrics for both evaluators
metrics_eval1 = calculate_metrics(eval1, criteria_cols)
metrics_eval2 = calculate_metrics(eval2, criteria_cols)

# Add evaluator labels
metrics_eval1["evaluator"] = "Evaluator 1"
metrics_eval2["evaluator"] = "Evaluator 2"

# Combine
combined_metrics = pd.concat([metrics_eval1, metrics_eval2])

# Overall statistics
eval1_total_correct = metrics_eval1["correct"].sum()
eval1_total_errors = metrics_eval1["errors"].sum()
eval1_total_data = metrics_eval1["data_issues"].sum()
eval1_analyzable = eval1_total_correct + eval1_total_errors
eval1_overall = (
    (eval1_total_correct / eval1_analyzable * 100) if eval1_analyzable > 0 else 0
)

eval2_total_correct = metrics_eval2["correct"].sum()
eval2_total_errors = metrics_eval2["errors"].sum()
eval2_total_data = metrics_eval2["data_issues"].sum()
eval2_analyzable = eval2_total_correct + eval2_total_errors
eval2_overall = (
    (eval2_total_correct / eval2_analyzable * 100) if eval2_analyzable > 0 else 0
)

print(f"Evaluator 1: {eval1_overall:.1f}% ({eval1_total_correct}/{eval1_analyzable})")
print(f"Evaluator 2: {eval2_overall:.1f}% ({eval2_total_correct}/{eval2_analyzable})")

# Create figure with multiple subplots
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

# ============================================================================
# Panel A: Correctness Rate by Criterion (Grouped Bar Chart)
# ============================================================================
ax1 = fig.add_subplot(gs[0, :])

x = np.arange(len(criteria_labels))
width = 0.35

eval1_rates = metrics_eval1["correctness_rate"].values
eval2_rates = metrics_eval2["correctness_rate"].values

bars1 = ax1.bar(
    x - width / 2,
    eval1_rates,
    width,
    label="Evaluator 1",
    color="#4C72B0",
    alpha=0.8,
    edgecolor="black",
    linewidth=0.5,
)
bars2 = ax1.bar(
    x + width / 2,
    eval2_rates,
    width,
    label="Evaluator 2",
    color="#DD8452",
    alpha=0.8,
    edgecolor="black",
    linewidth=0.5,
)

ax1.set_ylabel("Correctness Rate (%)", fontsize=11, fontweight="bold")
ax1.set_xlabel("Evaluation Criterion", fontsize=11, fontweight="bold")
ax1.set_xticks(x)
ax1.set_xticklabels(criteria_labels, fontsize=9)
ax1.legend(frameon=True, loc="lower right", fontsize=10)
ax1.set_ylim([90, 102])
ax1.axhline(
    y=95, color="gray", linestyle="--", linewidth=1, alpha=0.5, label="95% threshold"
)
ax1.grid(axis="y", alpha=0.3)

# Add Nature-style panel label
ax1.text(
    -0.08,
    1.05,
    "a",
    transform=ax1.transAxes,
    fontsize=18,
    fontweight="bold",
    va="top",
    ha="right",
    family="Arial",
)

# Add value labels on bars
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax1.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + 0.5,
            f"{height:.1f}%",
            ha="center",
            va="bottom",
            fontsize=7,
        )

# ============================================================================
# Panel B: Error Distribution (Stacked Bar)
# ============================================================================
ax2 = fig.add_subplot(gs[1, 0])

error_types = ["Correct", "Agent Errors", "Data Issues"]
eval1_counts = [eval1_total_correct, eval1_total_errors, eval1_total_data]
eval2_counts = [eval2_total_correct, eval2_total_errors, eval2_total_data]

x_pos = [0, 1]
width = 0.6

colors = ["#55A868", "#C44E52", "#8172B3"]

bottom1 = [0, 0]
bottom2 = [0, 0]

for i, (error_type, color) in enumerate(zip(error_types, colors)):
    if i == 0:
        ax2.bar(
            x_pos,
            [eval1_counts[i], eval2_counts[i]],
            width,
            label=error_type,
            color=color,
            alpha=0.8,
            edgecolor="black",
            linewidth=0.5,
        )
        bottom1[0] = eval1_counts[i]
        bottom2[1] = eval2_counts[i]
    else:
        ax2.bar(
            x_pos,
            [eval1_counts[i], eval2_counts[i]],
            width,
            bottom=[bottom1[0], bottom2[1]],
            label=error_type,
            color=color,
            alpha=0.8,
            edgecolor="black",
            linewidth=0.5,
        )
        bottom1[0] += eval1_counts[i]
        bottom2[1] += eval2_counts[i]

ax2.set_ylabel("Number of Criteria", fontsize=11, fontweight="bold")
ax2.set_xticks(x_pos)
ax2.set_xticklabels(["Evaluator 1", "Evaluator 2"], fontsize=10)
ax2.legend(frameon=True, fontsize=9)
ax2.set_ylim([0, 850])

# Add Nature-style panel label
ax2.text(
    -0.18,
    1.05,
    "b",
    transform=ax2.transAxes,
    fontsize=18,
    fontweight="bold",
    va="top",
    ha="right",
    family="Arial",
)

# Add percentage labels
for i, (e1, e2) in enumerate(zip(eval1_counts, eval2_counts)):
    pct1 = e1 / 800 * 100
    pct2 = e2 / 800 * 100
    if i == 0:
        y1, y2 = e1 / 2, e2 / 2
    elif i == 1:
        y1 = eval1_counts[0] + e1 / 2
        y2 = eval2_counts[0] + e2 / 2
    else:
        y1 = eval1_counts[0] + eval1_counts[1] + e1 / 2
        y2 = eval2_counts[0] + eval2_counts[1] + e2 / 2

    if pct1 > 1:
        ax2.text(
            0,
            y1,
            f"{pct1:.1f}%",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color="white",
        )
    if pct2 > 1:
        ax2.text(
            1,
            y2,
            f"{pct2:.1f}%",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color="white",
        )

# ============================================================================
# Panel C: Overall Correctness (Pie Charts)
# ============================================================================
ax3 = fig.add_subplot(gs[1, 1])

# Combined results
combined_correct = eval1_total_correct + eval2_total_correct
combined_errors = eval1_total_errors + eval2_total_errors
combined_total = combined_correct + combined_errors
combined_rate = combined_correct / combined_total * 100

sizes = [combined_correct, combined_errors]
labels = [
    f"Correct\n({combined_correct}/{combined_total})",
    f"Errors\n({combined_errors}/{combined_total})",
]
colors_pie = ["#55A868", "#C44E52"]
explode = (0.05, 0)

wedges, texts, autotexts = ax3.pie(
    sizes,
    explode=explode,
    labels=labels,
    colors=colors_pie,
    autopct="%1.1f%%",
    startangle=90,
    textprops={"fontsize": 10, "weight": "bold"},
)

# Add Nature-style panel label
ax3.text(
    -0.18,
    1.05,
    "c",
    transform=ax3.transAxes,
    fontsize=18,
    fontweight="bold",
    va="top",
    ha="right",
    family="Arial",
)

# ============================================================================
# Panel D: Comparison Statistics Table
# ============================================================================
ax4 = fig.add_subplot(gs[2, :])
ax4.axis("tight")
ax4.axis("off")

table_data = [
    ["Metric", "Evaluator 1", "Evaluator 2", "Combined"],
    ["Traces Evaluated", "100", "100", "200"],
    ["Total Criteria", "800", "800", "1,600"],
    [
        "Correct",
        f"{eval1_total_correct}",
        f"{eval2_total_correct}",
        f"{combined_correct}",
    ],
    [
        "Agent Errors",
        f"{eval1_total_errors}",
        f"{eval2_total_errors}",
        f"{combined_errors}",
    ],
    [
        "Data Issues",
        f"{eval1_total_data}",
        f"{eval2_total_data}",
        f"{eval1_total_data + eval2_total_data}",
    ],
    [
        "Correctness Rate",
        f"{eval1_overall:.1f}%",
        f"{eval2_overall:.1f}%",
        f"{combined_rate:.1f}%",
    ],
]

table = ax4.table(
    cellText=table_data, cellLoc="center", loc="center", colWidths=[0.25, 0.2, 0.2, 0.2]
)

table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1, 2)

# Style header row
for i in range(4):
    cell = table[(0, i)]
    cell.set_facecolor("#4C72B0")
    cell.set_text_props(weight="bold", color="white")

# Alternate row colors
for i in range(1, len(table_data)):
    for j in range(4):
        cell = table[(i, j)]
        if i % 2 == 0:
            cell.set_facecolor("#F0F0F0")
        else:
            cell.set_facecolor("white")

        # Bold the metric names
        if j == 0:
            cell.set_text_props(weight="bold")

# Add Nature-style panel label
ax4.text(
    -0.04,
    1.05,
    "d",
    transform=ax4.transAxes,
    fontsize=18,
    fontweight="bold",
    va="top",
    ha="right",
    family="Arial",
)

# ============================================================================
# Save figure
# ============================================================================
plt.tight_layout()
output_path = paths.REPO / "figures/validation_results.png"
plt.savefig(output_path, dpi=300, bbox_inches="tight")
print(f"\nFigure saved to: {output_path}")

# Also save as PDF for publication
output_path_pdf = output_path.replace(".png", ".pdf")
plt.savefig(output_path_pdf, bbox_inches="tight")
print(f"PDF saved to: {output_path_pdf}")

plt.show()

# ============================================================================
# Create supplementary figure: Detailed error breakdown
# ============================================================================
fig2, axes = plt.subplots(1, 2, figsize=(12, 5))

# Error types breakdown for each evaluator
error_categories = {
    "Evaluator 1": {
        "Forbidden age\ncomparison": 3,
        "No clear\ncomparison": 1,
        "Unrecognized\ntest": 1,
    },
    "Evaluator 2": {"No clear\ncomparison": 3, "Unrecognized\ntest": 3},
}

for idx, (evaluator, errors) in enumerate(error_categories.items()):
    ax = axes[idx]
    categories = list(errors.keys())
    counts = list(errors.values())
    colors_bar = sns.color_palette("Set2", len(categories))

    bars = ax.bar(
        categories,
        counts,
        color=colors_bar,
        alpha=0.8,
        edgecolor="black",
        linewidth=0.5,
    )
    ax.set_ylabel("Number of Errors", fontsize=11, fontweight="bold")
    ax.set_title(
        f"{evaluator}\n(Total: {sum(counts)} errors)", fontsize=12, fontweight="bold"
    )
    ax.set_ylim([0, max(counts) + 1])

    # Add value labels
    for bar in bars:
        height = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            height + 0.1,
            f"{int(height)}",
            ha="center",
            va="bottom",
            fontsize=10,
            fontweight="bold",
        )

plt.tight_layout()
output_path_supp = paths.REPO / "figures/validation_error_types.png"
plt.savefig(output_path_supp, dpi=300, bbox_inches="tight")
print(f"Supplementary figure saved to: {output_path_supp}")

output_path_supp_pdf = output_path_supp.replace(".png", ".pdf")
plt.savefig(output_path_supp_pdf, bbox_inches="tight")
print(f"PDF saved to: {output_path_supp_pdf}")

plt.show()

print("\n" + "=" * 80)
print("VISUALIZATION COMPLETE")
print("=" * 80)
print(f"\nMain figure: validation_results.png/.pdf")
print(f"Supplementary: validation_error_types.png/.pdf")
