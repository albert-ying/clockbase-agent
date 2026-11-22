#!/usr/bin/env python3
"""
Combined-modality ridge plot of effect-size distributions per intervention
category, pooled across mouse RNA-seq + human RNA-seq + human DNAm. Single
KDE per category (no per-modality stratification).

Drawn at its print size (124 x 50 mm, Figure 4b) and placed at 100% scale.
The data, categories, order, colours and encodings are those of the July 2026
generator; only size, type, line weights and label placement were changed.

Outputs:
  figures/create_ridge_plot_combined.png
  figures/create_ridge_plot_combined.pdf
"""

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import text_to_path
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

# Per-category statistics (Wilcoxon signed-rank vs zero on non-zero entries)
cat_stats = {}
for category in category_order:
    cat_data = df_sig[df_sig["condition_category"] == category]["log2FoldChange"].values
    nonzero = cat_data[cat_data != 0]
    if len(nonzero) > 0:
        wilcoxon_stat, wilcoxon_p = stats.wilcoxon(nonzero)
    else:
        wilcoxon_stat, wilcoxon_p = (np.nan, 1.0)
    cat_stats[category] = dict(data=cat_data, median=np.median(cat_data),
                               n=len(cat_data), p=wilcoxon_p)

total_sig = len(df_sig)

# ---------------------------------------------------------------------------
# Print-size typography and geometry (all lengths in mm on a 124 x 50 mm page)
# ---------------------------------------------------------------------------
sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
plt.rcParams.update({
    "font.family": "Arial",
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.linewidth": 0.6,
    "axes.unicode_minus": True,
})

TEXT_INK = plt.rcParams["text.color"]       # seaborn '.15', as in July
HINT_INK = "darkgray"
MM_PER_IN = 25.4
PT_PER_MM = 72 / 25.4
W_MM, H_MM = 124.0, 50.0

FS_CATEGORY = 7.0     # category name, bold
FS_COUNT = 6.0        # (n = ...) under the category name
FS_PBOX = 6.0         # P value box text
FS_PBOX_SUP = 5.0     # exponent in the P value box (>= 4.8 pt)
FS_LEGEND = 6.0       # reference-line key
FS_TICK = 6.5
FS_AXIS = 7.0
FS_HINT = 6.5
FS_TITLE = 7.0

LW_FILL_EDGE = 0.8    # July: 2.0 (white edge of the filled KDE)
LW_OUTLINE = 1.0      # July: 2.5 (coloured KDE outline)
LW_MEDIAN = 0.8       # July: 2.0
LW_ZERO = 0.5         # July: 1.0
ZERO_DASHES = (4.0, 2.2)   # in units of the line width


def text_width_pt(s, size, weight="normal", style="normal"):
    """Advance width of a plain Arial string in points (same metric as the PDF backend)."""
    prop = FontProperties(family="Arial", size=size, weight=weight, style=style)
    w, _, _ = text_to_path.get_text_width_height_descent(s, prop, ismath=False)
    return w


def space_pt(size):
    return 0.278 * size  # Arial word space, 569/2048 em


def run_width_pt(pieces):
    """pieces: list of dicts(text, size, weight, style, dy, gap_before)."""
    total = 0.0
    for pc in pieces:
        total += pc.get("gap_before", 0.0) + text_width_pt(
            pc["text"], pc["size"], pc.get("weight", "normal"), pc.get("style", "normal"))
    return total


def draw_run(fig, x_mm, y_mm, pieces, color, ha="left", zorder=4):
    """Draw a run of plain-text pieces on one baseline, left to right, in mm.

    Used where a bold string needs an italic statistic symbol or a superscript
    (mathtext cannot set bold italic, and its superscripts shrink to 0.7 em).
    Returns the run's left and right x in mm.
    """
    width_mm = run_width_pt(pieces) / PT_PER_MM
    if ha == "right":
        x_mm -= width_mm
    elif ha == "center":
        x_mm -= width_mm / 2
    x_left = x_mm
    for pc in pieces:
        x_mm += pc.get("gap_before", 0.0) / PT_PER_MM
        fig.text(x_mm / MM_PER_IN, (y_mm + pc.get("dy", 0.0) / PT_PER_MM) / MM_PER_IN,
                 pc["text"], transform=fig.dpi_scale_trans,
                 fontsize=pc["size"], fontweight=pc.get("weight", "normal"),
                 fontstyle=pc.get("style", "normal"), color=color,
                 ha="left", va="baseline", zorder=zorder)
        x_mm += text_width_pt(pc["text"], pc["size"], pc.get("weight", "normal"),
                              pc.get("style", "normal")) / PT_PER_MM
    return x_left, x_mm


def p_box_pieces(p):
    """P value at two significant figures: 6.67e-166 -> P = 6.7 x 10^-166."""
    if p < 0.0001:
        mant, exp = f"{p:.1e}".split("e")
        exp = str(int(exp)).replace("-", "−")
        return [
            dict(text="P", size=FS_PBOX, weight="bold", style="italic"),
            dict(text=f"= {mant}×10", size=FS_PBOX, weight="bold",
                 gap_before=space_pt(FS_PBOX)),
            dict(text=exp, size=FS_PBOX_SUP, weight="bold", dy=0.42 * FS_PBOX,
                 gap_before=0.15),
        ]
    return [
        dict(text="P", size=FS_PBOX, weight="bold", style="italic"),
        dict(text=f"= {p:.2g}", size=FS_PBOX, weight="bold",
             gap_before=space_pt(FS_PBOX)),
    ]


# Horizontal geometry from the widest labels
label_w_mm = max(text_width_pt(c, FS_CATEGORY, "bold") for c in category_order) / PT_PER_MM
PBOX_PAD_PT = 0.4 * FS_PBOX                  # July: round,pad=0.4 (em)
pbox_w_mm = max(run_width_pt(p_box_pieces(cat_stats[c]["p"])) for c in category_order) \
    / PT_PER_MM + 2 * PBOX_PAD_PT / PT_PER_MM
LABEL_GAP_MM = 1.6                           # July: 0.02 of the axes width
PBOX_GAP_MM = 1.6
EDGE_MM = 0.5
L_MM = EDGE_MM + label_w_mm + LABEL_GAP_MM
R_MM = W_MM - EDGE_MM - pbox_w_mm - PBOX_GAP_MM

# Vertical geometry: x-title row, tick labels, five rows, title row
Y_HINT_BASE = 1.3          # baseline of the x-title / hint row
Y0_MM = 6.5                # bottom of the lowest row (x axis)
YTOP_MM = 45.6             # top of the highest row
ROW_GAP_MM = 0.5
row_h_mm = (YTOP_MM - Y0_MM - (n_categories - 1) * ROW_GAP_MM) / n_categories

fig, axs = plt.subplots(nrows=n_categories, ncols=1, sharex=True,
                        figsize=(W_MM / MM_PER_IN, H_MM / MM_PER_IN))
axs = axs.flatten() if n_categories > 1 else [axs]
fig.subplots_adjust(left=L_MM / W_MM, right=R_MM / W_MM,
                    bottom=Y0_MM / H_MM, top=YTOP_MM / H_MM,
                    hspace=ROW_GAP_MM / row_h_mm)

bandwidth = 0.5

for i, category in enumerate(category_order):
    st = cat_stats[category]
    cat_data = st["data"]

    sns.kdeplot(cat_data, ax=axs[i], bw_adjust=bandwidth,
                fill=True, alpha=0.7, color=COLORS[category],
                linewidth=LW_FILL_EDGE, edgecolor="white")
    sns.kdeplot(cat_data, ax=axs[i], bw_adjust=bandwidth,
                fill=False, color=COLORS[category],
                linewidth=LW_OUTLINE, alpha=0.9)

    axs[i].axvline(0, color="black", linestyle=(0, ZERO_DASHES), alpha=0.5,
                   linewidth=LW_ZERO)
    axs[i].axvline(st["median"], color=COLORS[category], linestyle="-",
                   alpha=0.8, linewidth=LW_MEDIAN)

    row_bottom = Y0_MM + (n_categories - 1 - i) * (row_h_mm + ROW_GAP_MM)

    # Category name (7 pt bold) with its count (6 pt) underneath, right-aligned
    fig.text((L_MM - LABEL_GAP_MM) / W_MM, (row_bottom + 0.42 * row_h_mm) / H_MM,
             category, fontsize=FS_CATEGORY, fontweight="bold", color=TEXT_INK,
             ha="right", va="baseline")
    fig.text((L_MM - LABEL_GAP_MM) / W_MM, (row_bottom + 0.10 * row_h_mm) / H_MM,
             f"($n$ = {st['n']:,})", fontsize=FS_COUNT, color=TEXT_INK,
             ha="right", va="baseline")

    # P value box: white bold text on the category colour, as in July
    pieces = p_box_pieces(st["p"])
    base_y = row_bottom + 0.5 * row_h_mm - 0.36 * FS_PBOX / PT_PER_MM
    x_text = R_MM + PBOX_GAP_MM + PBOX_PAD_PT / PT_PER_MM
    x_l, x_r = draw_run(fig, x_text, base_y, pieces, "white", zorder=4)
    asc_mm, desc_mm = 0.905 * FS_PBOX / PT_PER_MM, 0.212 * FS_PBOX / PT_PER_MM
    box = mpatches.FancyBboxPatch(
        (x_l / MM_PER_IN, (base_y - desc_mm) / MM_PER_IN),
        (x_r - x_l) / MM_PER_IN, (asc_mm + desc_mm) / MM_PER_IN,
        boxstyle=f"round,pad={PBOX_PAD_PT / 72:.4f}",
        transform=fig.dpi_scale_trans, facecolor=COLORS[category],
        edgecolor=COLORS[category], alpha=0.9, linewidth=0.5, zorder=3)
    fig.add_artist(box)

    axs[i].set_yticks([])
    axs[i].set_ylabel("")
    axs[i].set_xlim(-25, 25)
    axs[i].spines["left"].set_visible(False)
    axs[i].spines["right"].set_visible(False)
    axs[i].spines["top"].set_visible(False)

    if i < n_categories - 1:
        axs[i].spines["bottom"].set_visible(False)
        axs[i].tick_params(axis="x", length=0, labelbottom=False)
    else:
        axs[i].spines["bottom"].set_visible(True)
        axs[i].spines["bottom"].set_linewidth(0.6)
        axs[i].set_xticks([-20, -10, 0, 10, 20])
        axs[i].tick_params(axis="x", length=0, pad=1.5, labelsize=FS_TICK,
                           labelbottom=True)

# x-axis title and the two direction hints share one row under the tick labels
axs[-1].set_xlabel("Cohen's $d$", fontsize=FS_AXIS, color=TEXT_INK)
axs[-1].xaxis.set_label_coords(((L_MM + R_MM) / 2) / W_MM, Y_HINT_BASE / H_MM,
                               transform=fig.transFigure)
axs[-1].xaxis.label.set_va("baseline")
fig.text(L_MM / W_MM, Y_HINT_BASE / H_MM, "Lower clock output (left)",
         fontsize=FS_HINT, color=HINT_INK, fontweight="bold",
         ha="left", va="baseline")
fig.text(R_MM / W_MM, Y_HINT_BASE / H_MM, "Higher clock output (right)",
         fontsize=FS_HINT, color=HINT_INK, fontweight="bold",
         ha="right", va="baseline")

# Short panel title: what is pooled and how many rows
TITLE_BASE = H_MM - 0.8 - 0.716 * FS_TITLE / PT_PER_MM
draw_run(fig, L_MM, TITLE_BASE, [
    dict(text="Significant rows (", size=FS_TITLE, weight="bold"),
    dict(text="n", size=FS_TITLE, weight="bold", style="italic"),
    dict(text=f"= {total_sig:,})", size=FS_TITLE, weight="bold",
         gap_before=space_pt(FS_TITLE)),
], TEXT_INK)

legend_elements = [
    mpatches.Patch(color="none", label="Reference lines:"),
    plt.Line2D([0], [0], color="black", linestyle=(0, ZERO_DASHES), alpha=0.5,
               linewidth=LW_ZERO, label="No effect (Cohen's $d$ = 0)"),
    plt.Line2D([0], [0], color="gray", linestyle="-", linewidth=LW_MEDIAN,
               label="Category median"),
]
axs[0].legend(handles=legend_elements, loc="upper right",
              bbox_to_anchor=(R_MM / W_MM, (H_MM - 0.5) / H_MM),
              bbox_transform=fig.transFigure, frameon=False, fontsize=FS_LEGEND,
              handlelength=2.4, handletextpad=0.5, labelspacing=0.25,
              borderpad=0.0, borderaxespad=0.0)

plt.savefig("figures/create_ridge_plot_combined.pdf", dpi=300)
plt.savefig("figures/create_ridge_plot_combined.png", dpi=600)
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
summary_rows = []
for category in category_order:
    cat_data = df_sig[df_sig["condition_category"] == category]["log2FoldChange"].values
    nonzero = cat_data[cat_data != 0]
    if len(nonzero) > 0:
        w_stat, w_p = stats.wilcoxon(nonzero)
    else:
        w_stat, w_p = (np.nan, 1.0)
    print(f"\n{category}: n={len(cat_data):,}, median={np.median(cat_data):.3f}, mean={np.mean(cat_data):.3f}")
    print(f"  Wilcoxon vs 0: W={w_stat:.1f}, p={w_p:.6e}")
    summary_rows.append({'category': category, 'n': len(cat_data),
                         'median_d': np.median(cat_data), 'wilcoxon_p': w_p,
                         'outside_display_n': int(((cat_data < axs[-1].get_xlim()[0]) |
                                                   (cat_data > axs[-1].get_xlim()[1])).sum())})
pd.DataFrame(summary_rows).to_csv('result/panel_b_ridge_statistics.csv', index=False)
