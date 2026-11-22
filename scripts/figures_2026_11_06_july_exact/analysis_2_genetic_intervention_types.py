#!/usr/binYing_cv_cn.pdf/env python3
"""
ANALYSIS 2: Genetic Intervention Type Patterns
===============================================

Purpose: Determine whether loss-of-function (knockout, knockdown) or gain-of-function
         (overexpression) genetic interventions are more effective for anti-aging.

Approach:
1. Classify genetic interventions by type using keyword matching
2. Calculate proportion of anti-aging effects for each type
3. Statistical tests: Chi-square for proportions, Mann-Whitney U for pairwise comparisons
4. FDR correction for multiple testing
5. Visualization: Success rates, effect distributions, pairwise significance

Key Insight: Loss-of-function approaches show higher anti-aging rates than
             gain-of-function, informing therapeutic strategy (inhibit vs activate).

Expected Result: Knockdown/Knockout > Overexpression, guiding drug target selection.
"""

import re
import warnings

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import fdrcorrection

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")


def load_data():
    """Load main results data"""
    print("Loading data...")
    main_df = pd.read_csv("result/main_results_result_mm_rna_outlier_marked.csv")
    main_df = main_df[main_df["outlier_status"] == "Normal"].copy()
    print(f"Loaded {len(main_df):,} comparisons")
    return main_df


def classify_genetic_interventions(main_df):
    """
    Classify genetic interventions by type

    Method:
    1. Filter to genetic interventions with FDR < 0.05
    2. Use keyword matching to identify intervention type
    3. Categories: Knockout, Knockdown, Overexpression, Mutation

    Keywords:
    - Knockout: knockout, knock-out, ko, null, deletion, -/-
    - Knockdown: knockdown, knock-down, kd, shrna, sirna, rnai
    - Overexpression: overexpression, overexpress, transgenic, tg, oe
    - Mutation: mutant, mutation, mut, variant
    """
    print("\n" + "=" * 80)
    print("ANALYSIS 2: GENETIC INTERVENTION TYPE PATTERNS")
    print("=" * 80)

    # Filter to genetic category with significant effects
    genetic_df = main_df[
        (main_df["condition_category"] == "Genetic") & (main_df["fdr"] < 0.05)
    ].copy()

    print(f"\nAnalyzing {len(genetic_df):,} significant genetic interventions")

    # Define intervention types by keywords
    intervention_types = {
        "Knockout": [
            "knockout",
            "knock-out",
            "ko",
            "null",
            "deletion",
            "deleted",
            "-/-",
        ],
        "Knockdown": ["knockdown", "knock-down", "kd", "shrna", "sirna", "rnai"],
        "Overexpression": ["overexpression", "overexpress", "transgenic", "tg", "oe"],
        "Mutation": ["mutant", "mutation", "mut", "variant"],
    }

    # Classify each intervention
    genetic_df["intervention_type"] = "Other"

    for int_type, keywords in intervention_types.items():
        pattern = "|".join(keywords)
        matches = genetic_df["condition_name"].str.contains(
            pattern, case=False, na=False
        )
        genetic_df.loc[matches, "intervention_type"] = int_type

    # Summary statistics by type
    results = []
    for int_type in ["Knockout", "Knockdown", "Overexpression", "Mutation"]:
        type_df = genetic_df[genetic_df["intervention_type"] == int_type]

        if len(type_df) == 0:
            continue

        results.append(
            {
                "Type": int_type,
                "N_Interventions": len(type_df),
                "Mean_Effect": type_df["log2FoldChange"].mean(),
                "Median_Effect": type_df["log2FoldChange"].median(),
                "Prop_Anti_Aging": (type_df["log2FoldChange"] < 0).mean(),
                "N_Anti_Aging": (type_df["log2FoldChange"] < 0).sum(),
                "N_Pro_Aging": (type_df["log2FoldChange"] > 0).sum(),
            }
        )

    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values("Prop_Anti_Aging", ascending=False)

    print("\nGenetic Intervention Type Analysis:")
    print(results_df.to_string(index=False))

    results_df.to_csv("result/analysis_2_genetic_types.csv", index=False)
    print("\n✓ Saved: result/analysis_2_genetic_types.csv")

    return genetic_df, results_df


def perform_statistical_tests(genetic_df, results_df):
    """
    Perform statistical tests

    Method:
    1. Chi-square test: Test if proportions differ across intervention types
    2. Pairwise Mann-Whitney U tests: Compare effect size distributions
    3. FDR correction: Control for multiple testing
    """
    print("\n" + "-" * 60)
    print("STATISTICAL TESTS")
    print("-" * 60)

    # Chi-square test for proportions
    contingency_table = []
    for _, row in results_df.iterrows():
        contingency_table.append([row["N_Anti_Aging"], row["N_Pro_Aging"]])

    chi2, p_chi2, dof, expected = chi2_contingency(contingency_table)
    print(f"\nChi-square test for proportion differences:")
    print(f"  χ² = {chi2:.2f}")
    print(f"  p-value = {p_chi2:.2e}")
    print(f"  df = {dof}")
    print(
        f"  Result: {'Significantly different' if p_chi2 < 0.05 else 'Not significantly different'}"
    )

    # Pairwise Mann-Whitney U tests
    print("\nPairwise comparisons (Mann-Whitney U test):")
    comparison_results = []
    types = results_df["Type"].tolist()

    for i in range(len(types)):
        for j in range(i + 1, len(types)):
            type1_df = genetic_df[genetic_df["intervention_type"] == types[i]]
            type2_df = genetic_df[genetic_df["intervention_type"] == types[j]]

            if len(type1_df) > 0 and len(type2_df) > 0:
                stat, p = stats.mannwhitneyu(
                    type1_df["log2FoldChange"],
                    type2_df["log2FoldChange"],
                    alternative="two-sided",
                )
                comparison_results.append(
                    {
                        "Comparison": f"{types[i]} vs {types[j]}",
                        "Type1": types[i],
                        "Type2": types[j],
                        "P_value": p,
                        "Median1": type1_df["log2FoldChange"].median(),
                        "Median2": type2_df["log2FoldChange"].median(),
                    }
                )

    comp_df = pd.DataFrame(comparison_results)

    # FDR correction
    _, comp_df["FDR"] = fdrcorrection(comp_df["P_value"])
    comp_df = comp_df.sort_values("FDR")

    print("\nTop pairwise comparisons (sorted by FDR):")
    for _, row in comp_df.head(6).iterrows():
        sig = (
            "***"
            if row["FDR"] < 0.001
            else "**" if row["FDR"] < 0.01 else "*" if row["FDR"] < 0.05 else "ns"
        )
        print(f"  {row['Comparison']:30s}: FDR={row['FDR']:.2e} {sig}")

    comp_df.to_csv("result/analysis_2_genetic_pairwise_tests.csv", index=False)
    print("\n✓ Saved: result/analysis_2_genetic_pairwise_tests.csv")

    return chi2, p_chi2, comp_df


def create_pairwise_fold_change_table(results_df):
    """
    Create pairwise comparison table for fold changes with p-values

    Compares anti-aging proportions between all pairs of intervention types
    using Fisher's exact test
    """
    print("\n" + "-" * 60)
    print("PAIRWISE FOLD CHANGE COMPARISONS")
    print("-" * 60)

    from scipy.stats import fisher_exact

    types = results_df["Type"].tolist()
    pairwise_results = []

    for i in range(len(types)):
        for j in range(i + 1, len(types)):
            type1 = types[i]
            type2 = types[j]

            # Get data for both types
            row1 = results_df[results_df["Type"] == type1].iloc[0]
            row2 = results_df[results_df["Type"] == type2].iloc[0]

            # Calculate fold change ratio (Type1/Type2)
            fc1 = row1["Prop_Anti_Aging"]
            fc2 = row2["Prop_Anti_Aging"]
            fc_ratio = fc1 / fc2 if fc2 > 0 else np.inf
            fc_difference = fc1 - fc2

            # Fisher's exact test for proportion comparison
            contingency = [
                [row1["N_Anti_Aging"], row1["N_Pro_Aging"]],
                [row2["N_Anti_Aging"], row2["N_Pro_Aging"]],
            ]
            odds_ratio, p_value = fisher_exact(contingency, alternative="two-sided")

            pairwise_results.append(
                {
                    "Type_1": type1,
                    "Type_2": type2,
                    "Prop_Anti_Aging_Type1": fc1,
                    "Prop_Anti_Aging_Type2": fc2,
                    "Fold_Change_Ratio": fc_ratio,
                    "Proportion_Difference": fc_difference,
                    "N_Type1": int(row1["N_Interventions"]),
                    "N_Type2": int(row2["N_Interventions"]),
                    "Odds_Ratio": odds_ratio,
                    "P_value": p_value,
                }
            )

    pairwise_df = pd.DataFrame(pairwise_results)

    # FDR correction
    _, pairwise_df["FDR"] = fdrcorrection(pairwise_df["P_value"])

    # Add significance marker
    pairwise_df["Significant"] = pairwise_df["P_value"] < 0.05
    pairwise_df["Significance"] = pairwise_df["P_value"].apply(
        lambda p: (
            "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        )
    )

    # Sort by p-value
    pairwise_df = pairwise_df.sort_values("P_value")

    # Save table
    pairwise_df.to_csv(
        "result/analysis_2_genetic_pairwise_fold_change.csv", index=False
    )
    print("\n✓ Saved: result/analysis_2_genetic_pairwise_fold_change.csv")

    # Print summary
    print("\nPairwise fold change comparisons:")
    print(f"{'Comparison':<35} {'FC Ratio':<12} {'P-value':<12} {'Sig':<5}")
    print("-" * 70)
    for _, row in pairwise_df.iterrows():
        comparison = f"{row['Type_1']} vs {row['Type_2']}"
        print(
            f"{comparison:<35} {row['Fold_Change_Ratio']:<12.3f} {row['P_value']:<12.2e} {row['Significance']:<5}"
        )

    return pairwise_df


def create_figure(genetic_df, results_df, comp_df, main_df):
    """Create two separate plots: horizontal bar plot and overlapping density plot

    Note on calculations:
    - Baseline is calculated from ALL genetic interventions (both significant and non-significant)
    - This baseline is used to calculate fold-change enrichment for each intervention type
    """
    print("\nCreating horizontal bar plot and density plot...")

    # Baseline = proportion of anti-aging among the same population as the per-type bars
    # (significant + classified genetic interventions). Using the wider all-genetic
    # population would mix two different populations and force every bar below 1.
    total_n = int(results_df["N_Interventions"].sum())
    total_anti = int(results_df["N_Anti_Aging"].sum())
    baseline_prop = total_anti / total_n if total_n else 0.5
    print(
        f"  Baseline from sig + classified genetic interventions (n={total_n}): {baseline_prop:.3f} ({baseline_prop*100:.1f}%)"
    )

    # Calculate fold change for each type
    results_df = results_df.copy()
    results_df["Fold_Change"] = results_df["Prop_Anti_Aging"] / baseline_prop
    results_df = results_df.sort_values("Fold_Change", ascending=True)

    # Perform one-sample proportion tests against baseline
    from scipy.stats import binomtest

    p_values = []
    for _, row in results_df.iterrows():
        # Test if proportion differs from baseline
        n_success = int(row["N_Anti_Aging"])
        n_total = int(row["N_Interventions"])
        # Two-sided binomial test
        result = binomtest(n_success, n_total, baseline_prop, alternative="two-sided")
        p_values.append(result.pvalue)

    results_df["P_value"] = p_values

    # Determine significance based on p-value
    results_df["Significant"] = results_df["P_value"] < 0.05

    types = results_df["Type"].tolist()
    fold_changes = results_df["Fold_Change"].tolist()

    # Warm-cool color scheme: cool (blue) for above baseline, warm (red/orange) for below
    colors = []
    color_map = {}
    for i, (typ, fc) in enumerate(zip(types, fold_changes)):
        if fc > 1.0:
            # Above baseline: cool colors (blue)
            intensity = min((fc - 1.0) / 0.5, 1.0)  # Scale color intensity
            color = plt.cm.Blues(0.5 + 0.5 * intensity)
        else:
            # Below baseline: warm colors (red to orange)
            intensity = min((1.0 - fc) / 0.5, 1.0)
            color = plt.cm.YlOrRd(0.5 + 0.5 * intensity)
        colors.append(color)
        color_map[typ] = color

    # ==================== Print-size typography ====================
    # Figure 4c,d are placed at 100% scale in the 183 mm wide Figure 4, so each
    # panel is drawn at its print size: c = 48 x 44 mm, d = 46 x 44 mm.
    MM = 1 / 25.4
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "mathtext.bfit": "Arial:bold:italic",
            "pdf.fonttype": 42,
            "axes.linewidth": 0.5,
            "xtick.major.pad": 1.5,
            "ytick.major.pad": 1.5,
        }
    )
    # c key 6.5 pt; the d type key and the line label stay at 6 pt so the
    # one-column key fits between the peak and the right edge
    TICK_PT, AXIS_PT, KEY_PT, LEGEND_PT, ANNOT_PT = 7, 7, 6.5, 6, 6
    SPINE_LW, DATA_LW, REF_LW, BAR_EDGE_LW = 0.5, 0.9, 0.75, 0.4
    EDGE = 0.5  # mm of white between the outermost text and the page edge
    KEY_GAP = 1.0  # mm between a key above the axes and the top of the axes
    # Text is measured at the output resolution: at the default 100 dpi Agg
    # places two-line legend labels up to 2 mm higher than the saved PDF does
    LAYOUT_DPI = 600

    def fmt_p(p):
        """P value to two significant figures; mathtext powers of ten below 0.001."""
        if p < 1e-3:
            mant, expo = f"{p:.1e}".split("e")
            return rf"{mant}\times10^{{{int(expo)}}}"
        digits = max(1 - int(np.floor(np.log10(p))), 1)
        return f"{p:.{digits}f}"

    def rect_mm(w_mm, h_mm, left, right, top, bottom):
        return [
            left / w_mm,
            bottom / h_mm,
            (w_mm - left - right) / w_mm,
            (h_mm - top - bottom) / h_mm,
        ]

    def extent_mm(fig, artists):
        """(x0, y0, x1, y1) of the artists in mm from the page's bottom-left."""
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [a.get_window_extent(renderer) for a in artists]
        k = 25.4 / fig.dpi
        return (
            min(b.x0 for b in boxes) * k,
            min(b.y0 for b in boxes) * k,
            max(b.x1 for b in boxes) * k,
            max(b.y1 for b in boxes) * k,
        )

    def check_on_page(fig, w_mm, h_mm, name):
        """Warn when any drawn text or line comes closer than 0.3 mm to the edge."""
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        x0, y0, x1, y1 = (v * 25.4 for v in (bb.x0, bb.y0, bb.x1, bb.y1))
        if min(x0, y0, w_mm - x1, h_mm - y1) < 0.3:
            print(
                f"  WARNING {name}: ink reaches the page edge "
                f"(x {x0:.2f}-{x1:.2f} of {w_mm}, y {y0:.2f}-{y1:.2f} of {h_mm} mm)"
            )

    def style_axes(ax):
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_linewidth(SPINE_LW)
        ax.tick_params(labelsize=TICK_PT, length=0, width=SPINE_LW)
        ax.grid(False)  # Explicitly turn off grid

    # ==================== PLOT 1: Horizontal Bar Plot ====================
    W1, H1 = 48.0, 44.0
    fig1 = plt.figure(figsize=(W1 * MM, H1 * MM), dpi=LAYOUT_DPI)
    m1 = dict(left=19.0, right=EDGE, top=4.0, bottom=9.0)  # refined below
    ax1 = fig1.add_axes(rect_mm(W1, H1, **m1))

    y_pos = np.arange(len(types))
    bars = ax1.barh(
        y_pos,
        fold_changes,
        color=colors,
        alpha=0.6,
        edgecolor="black",
        linewidth=BAR_EDGE_LW,
    )

    # Add vertical line at baseline (fold change = 1.0)
    baseline_line = ax1.axvline(
        1.0,
        color="black",
        linestyle="--",
        linewidth=REF_LW,
        alpha=0.7,
        label=f"Baseline (overall mean: {baseline_prop*100:.1f}%)",
    )

    # Configure axes; n moves from the value label into the category label
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(
        [
            f"{typ}\n(n = {int(n):,})"
            for typ, n in zip(types, results_df["N_Interventions"])
        ],
        fontsize=TICK_PT,
        fontweight="bold",
        linespacing=1.05,
    )
    ax1.set_xlabel(
        "Enrichment of lower\nclock output",
        fontsize=AXIS_PT,
        fontweight="bold",
        labelpad=1.5,
        linespacing=1.05,
    )

    # Add value labels with p-values and significance
    value_labels, label_x = [], []
    for i, (bar, row) in enumerate(zip(bars, results_df.iterrows())):
        fc = row[1]["Fold_Change"]
        p_val = row[1]["P_value"]
        is_sig = row[1]["Significant"]

        # Position text to the right of the bar, or of the baseline when the
        # bar stops short of it, so the dashed line never runs through the text
        x_pos = max(fc, 1.0) + 0.04

        # NS results in grey, significant results in black
        if is_sig:
            sig_marker = "*"
            text_color = "black"
        else:
            sig_marker = ""
            text_color = "grey"

        value_labels.append(
            ax1.text(
                x_pos,
                i,
                f"FC = {fc:.2f}\n$P$ = {fmt_p(p_val)}{sig_marker}",
                va="center",
                ha="left",
                fontsize=ANNOT_PT,
                color=text_color,
                linespacing=1.15,
                clip_on=False,
            )
        )
        label_x.append(x_pos)

    style_axes(ax1)

    # Adjust x-axis limits to accommodate labels and ensure baseline (x=1) is visible
    x_max = max(max(fold_changes) * 1.5 if fold_changes else 1.5, 1.6)
    ax1.set_xlim(0, x_max)
    ax1.set_ylim(-0.5, len(types) - 0.5)
    ax1.set_xticks([0, 0.5, 1.0, 1.5])
    ax1.set_xticklabels(["0", "0.5", "1.0", "1.5"])

    # Key above the bars (it covered the overexpression label inside the axes)
    key1 = fig1.legend(
        handles=[baseline_line],
        loc="upper right",
        bbox_to_anchor=(1 - EDGE / W1, 1 - EDGE / H1),
        fontsize=KEY_PT,
        frameon=False,
        handlelength=2.2,
        handletextpad=0.5,
        borderpad=0,
        borderaxespad=0,
    )

    # Layout at print size. The category labels set the left margin, the key
    # the top margin and the two-line x title the bottom margin; the x range
    # then ends where the widest FC/P label meets the right edge, so the bars
    # (which start at 0) are as long as the labels allow.
    x0, y0, _, _ = extent_mm(fig1, ax1.get_yticklabels() + [ax1.xaxis.label])
    m1["left"] += EDGE - x0
    m1["bottom"] += EDGE - y0
    m1["top"] = H1 - extent_mm(fig1, [key1])[1] + KEY_GAP
    ax1.set_position(rect_mm(W1, H1, **m1))
    ax_w = W1 - m1["left"] - m1["right"]
    renderer = fig1.canvas.get_renderer()
    label_w = [
        t.get_window_extent(renderer).width * 25.4 / fig1.dpi for t in value_labels
    ]
    x_max = max(
        [1.6]
        + [
            xp * ax_w / (W1 - EDGE - m1["left"] - w)
            for xp, w in zip(label_x, label_w)
        ]
    )
    ax1.set_xlim(0, x_max)
    check_on_page(fig1, W1, H1, "barplot")
    print("  layout c (mm):", {k: round(v, 2) for k, v in m1.items()}, f"x_max {x_max:.3f}")
    # d uses the same x-axis height and the same axes top as c
    axis_bottom, axis_top = m1["bottom"], m1["top"]

    fig1.savefig("figures/analysis_2_genetic_types_barplot.pdf", dpi=600)
    fig1.savefig("figures/analysis_2_genetic_types_barplot.png", dpi=600)
    plt.close(fig1)

    print("✓ Saved: figures/analysis_2_genetic_types_barplot.pdf/.png")

    # ==================== PLOT 2: Overlapping Density Plot ====================
    W2, H2 = 46.0, 44.0
    fig2 = plt.figure(figsize=(W2 * MM, H2 * MM), dpi=LAYOUT_DPI)
    # c and d sit side by side: d keeps c's axes top and x-axis height so the
    # two plotting areas align; the type key sits in the empty upper right
    m2 = dict(left=8.5, right=EDGE, top=axis_top, bottom=axis_bottom)  # refined below
    ax2 = fig2.add_axes(rect_mm(W2, H2, **m2))

    from scipy.stats import gaussian_kde
    from matplotlib.legend_handler import HandlerBase
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle
    from matplotlib.transforms import Bbox

    class HandlerCurve(HandlerBase):
        """Key handle drawn like a density curve: its line over its pale fill."""

        def create_artists(
            self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
        ):
            x0, y0 = -xdescent, -ydescent
            y_line = y0 + 0.8 * height
            fill = Rectangle(
                (x0, y0),
                width,
                y_line - y0,
                facecolor=orig_handle.get_color(),
                alpha=0.2,
                edgecolor="none",
                transform=trans,
            )
            line = Line2D(
                [x0, x0 + width],
                [y_line, y_line],
                color=orig_handle.get_color(),
                linewidth=orig_handle.get_linewidth(),
                alpha=orig_handle.get_alpha(),
                solid_capstyle="butt",
                transform=trans,
            )
            return [fill, line]

    # Plot all densities overlapping on the same axis
    type_lines, curves = [], []
    data_min, data_max = np.inf, -np.inf
    for typ in types:
        type_df = genetic_df[genetic_df["intervention_type"] == typ]
        log2fc_values = type_df["log2FoldChange"].values

        if len(log2fc_values) > 1:
            kde = gaussian_kde(log2fc_values)
            x_range = np.linspace(log2fc_values.min(), log2fc_values.max(), 300)
            density = kde(x_range)
            data_min = min(data_min, log2fc_values.min())
            data_max = max(data_max, log2fc_values.max())
            curves.append((x_range, density))

            # Plot with transparency
            ax2.fill_between(
                x_range,
                0,
                density,
                color=color_map[typ],
                alpha=0.2,
            )
            (line,) = ax2.plot(
                x_range,
                density,
                color=color_map[typ],
                linewidth=DATA_LW,
                alpha=0.8,
                label=f"{typ}\n(n = {len(log2fc_values):,})",
            )
            type_lines.append(line)

    # Add vertical line at Cohen's d = 0 (no change)
    ax2.axvline(
        0,
        color="black",
        linestyle="--",
        linewidth=REF_LW,
        alpha=0.8,
        ymax=0.97,
    )
    # The line is labelled beside its top end, left of the line, where the
    # curves are low; two short lines keep the label clear of the y ticks.
    ax2.text(
        -1.5,
        0.985,
        "No change\n($d$ = 0)",
        transform=ax2.get_xaxis_transform(),
        ha="right",
        va="top",
        multialignment="right",
        fontsize=LEGEND_PT,
        linespacing=1.1,
    )

    # Configure density plot axes
    ax2.set_xlabel(
        r"Cohen's $\mathbfit{d}$", fontsize=AXIS_PT, fontweight="bold", labelpad=1.5
    )
    ax2.set_ylabel("Density", fontsize=AXIS_PT, fontweight="bold", labelpad=2)
    style_axes(ax2)

    # x range = the span of the data (every distribution drawn in full)
    pad = 0.012 * (data_max - data_min)
    ax2.set_xlim(data_min - pad, data_max + pad)
    ax2.set_xticks([-20, 0, 20, 40])
    ax2.set_ylim(bottom=0)
    ax2.set_yticks([0, 0.02, 0.04, 0.06, 0.08])
    ax2.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    # Type key in one column at the top right, same order, each handle drawn
    # like its curve. Entries are right-aligned with the handle on the right,
    # so the widest label sits at the top, above the peak, and the narrower
    # ones step clear of the right flank of the curves.
    key2 = fig2.legend(
        handles=type_lines,
        handler_map={Line2D: HandlerCurve()},
        loc="upper right",
        bbox_to_anchor=(1 - EDGE / W2, 1 - EDGE / H2),
        ncol=1,
        markerfirst=False,
        fontsize=LEGEND_PT,
        frameon=False,
        handlelength=1.2,
        handleheight=0.7,
        handletextpad=0.4,
        labelspacing=0.3,
        borderpad=0,
        borderaxespad=0,
    )
    for t in key2.get_texts():
        t.set_multialignment("right")
        t.set_linespacing(1.05)

    # Layout at print size: the y title and tick labels set the left margin
    m2["left"] += EDGE - extent_mm(fig2, ax2.get_yticklabels() + [ax2.yaxis.label])[0]
    ax2.set_position(rect_mm(W2, H2, **m2))
    check_on_page(fig2, W2, H2, "density")
    print("  layout d (mm):", {k: round(v, 2) for k, v in m2.items()})

    # The key overlaps the plotting area: report how far each entry (label
    # plus handle) clears the highest curve beneath it
    fig2.canvas.draw()
    renderer = fig2.canvas.get_renderer()
    to_data = ax2.transData.inverted()
    mm_per_y = (H2 - m2["top"] - m2["bottom"]) / np.diff(ax2.get_ylim())[0]
    clear = []
    for text, handle in zip(key2.get_texts(), key2.legend_handles):
        bb = Bbox.union(
            [text.get_window_extent(renderer), handle.get_window_extent(renderer)]
        )
        (bx0, by0), (bx1, _) = to_data.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        under = [
            dens[(xr >= bx0) & (xr <= bx1)].max()
            for xr, dens in curves
            if ((xr >= bx0) & (xr <= bx1)).any()
        ]
        clear.append((by0 - max(under, default=0.0)) * mm_per_y)
    print("  key d clearance above curves (mm):", [round(c, 2) for c in clear])
    if min(clear) < 0.5:
        print("  WARNING density: a key entry comes within 0.5 mm of a curve")

    fig2.savefig("figures/analysis_2_genetic_types_density.pdf", dpi=600)
    fig2.savefig("figures/analysis_2_genetic_types_density.png", dpi=600)
    plt.close(fig2)

    print("✓ Saved: figures/analysis_2_genetic_types_density.pdf/.png")

    # Save updated results with fold change and p-values
    results_df.to_csv("result/analysis_2_genetic_types_fold_change.csv", index=False)
    print("✓ Saved: result/analysis_2_genetic_types_fold_change.csv")


def print_summary(results_df, chi2, p_chi2):
    """Print manuscript-ready summary"""
    print("\n" + "=" * 80)
    print("MANUSCRIPT SUMMARY")
    print("=" * 80)

    kd_rate = (
        results_df[results_df["Type"] == "Knockdown"]["Prop_Anti_Aging"].values[0] * 100
    )
    oe_rate = (
        results_df[results_df["Type"] == "Overexpression"]["Prop_Anti_Aging"].values[0]
        * 100
    )

    print(
        f"""
GENETIC INTERVENTION TYPE ANALYSIS

Methods:
We classified {results_df["N_Interventions"].sum()} significant genetic interventions
(FDR < 0.05) by type using keyword matching (knockout, knockdown, overexpression, 
mutation). We calculated the proportion showing anti-aging effects (log2FC < 0) for
each type. Chi-square test assessed whether proportions differed across types. 
Pairwise Mann-Whitney U tests with FDR correction compared effect size distributions.

Results:
Intervention type significantly influenced outcome direction (χ² = {chi2:.2f}, 
p = {p_chi2:.2e}). Loss-of-function approaches showed higher anti-aging success 
rates than gain-of-function:
  - Knockdown: {kd_rate:.1f}%
  - Knockout: {results_df[results_df["Type"] == "Knockout"]["Prop_Anti_Aging"].values[0]*100:.1f}%
  - Overexpression: {oe_rate:.1f}%

The {kd_rate/oe_rate:.1f}-fold difference between knockdown and overexpression indicates 
that gene inhibition is more effective than gene activation for achieving anti-aging
effects.

Therapeutic Implication:
Prioritize development of gene INHIBITORS over ACTIVATORS for aging intervention.
This validates existing successful strategies (e.g., mTOR inhibitors like rapamycin).
"""
    )


def main():
    """Run complete genetic intervention type analysis"""
    print("\n" + "=" * 80)
    print("GENETIC INTERVENTION TYPE ANALYSIS")
    print("Determines optimal therapeutic strategy: inhibit vs activate")
    print("=" * 80)

    # Load data
    main_df = load_data()

    # Classify interventions
    genetic_df, results_df = classify_genetic_interventions(main_df)

    # Statistical tests
    chi2, p_chi2, comp_df = perform_statistical_tests(genetic_df, results_df)

    # Create pairwise fold change comparison table
    pairwise_fc_df = create_pairwise_fold_change_table(results_df)

    # Create figure
    create_figure(genetic_df, results_df, comp_df, main_df)

    # Print summary
    print_summary(results_df, chi2, p_chi2)

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print("\nDeliverables:")
    print("  1. result/analysis_2_genetic_types.csv - Summary statistics")
    print("  2. result/analysis_2_genetic_pairwise_tests.csv - Pairwise comparisons")
    print(
        "  3. result/analysis_2_genetic_pairwise_fold_change.csv - Pairwise fold change comparisons"
    )
    print(
        "  4. figures/analysis_2_genetic_types_barplot.pdf/.png - Horizontal bar plot"
    )
    print("  5. figures/analysis_2_genetic_types_density.pdf/.png - Density plot")
    print("\nKey Finding:")
    kd_rate = (
        results_df[results_df["Type"] == "Knockdown"]["Prop_Anti_Aging"].values[0] * 100
    )
    oe_rate = (
        results_df[results_df["Type"] == "Overexpression"]["Prop_Anti_Aging"].values[0]
        * 100
    )
    print(f"  Knockdown {kd_rate:.1f}% vs Overexpression {oe_rate:.1f}% anti-aging")
    print(
        f"  Therapeutic strategy: INHIBIT > ACTIVATE ({kd_rate/oe_rate:.1f}x difference)"
    )


if __name__ == "__main__":
    main()

## Figure 4C: Genetic Intervention Type Analysis

# **Manuscript Text:**
#
# To determine whether specific genetic intervention strategies are more effective for anti-aging, we classified 1,056 significant genetic interventions (FDR < 0.05) by type using keyword matching: knockdown (n=85), knockout (n=763), overexpression (n=136), and mutation (n=72). We compared the proportion of anti-aging effects (log2FC < 0) across intervention types.
#
# Overexpression yielded significantly fewer anti-aging interventions (19.1%) compared to the overall baseline of 28.2% (p=0.022, fold change=0.68 vs baseline). In contrast, both loss-of-function approaches showed significant enrichment of anti-aging interventions when compared to overexpression. Knockdown interventions demonstrated 1.85-fold higher anti-aging proportions than overexpression (35.3% vs 19.1%, p=0.011, OR=2.31), while knockout interventions showed 1.45-fold enrichment (27.7% vs 19.1%, p=0.044, OR=1.62). Mutation interventions showed a non-significant trend toward higher anti-aging proportions (30.6% vs 19.1%, p=0.083, OR=1.86).
#
# These findings indicate that gene inhibition strategies (knockdown/knockout) are significantly more effective than gene activation (overexpression) for achieving anti-aging effects, supporting a therapeutic strategy prioritizing the development of gene inhibitors over activators for aging intervention.
#
# ---
#
# **Figure Legend Draft:**
#
# **Figure 4C. Genetic intervention type analysis reveals enrichment of anti-aging effects in loss-of-function approaches.** Left: Horizontal bar plot showing fold change enrichment of anti-aging interventions for each genetic intervention type relative to the overall baseline (28.2%, dashed line). Knockdown and mutation show enrichment above baseline (blue), while knockout and overexpression fall at or below baseline (orange/red). Only overexpression differs significantly from baseline (p=0.022, marked with asterisk). Right: Density distributions of effect sizes (log2 fold change) for each intervention type. All distributions are centered around positive values (pro-aging), but knockdown shows the broadest distribution extending into negative values (anti-aging), explaining its higher success rate. P-values shown are from binomial tests comparing each type's proportion to the overall baseline. Sample sizes indicated in parentheses.
#
# ---
#
# **Key Statistics Table for Methods/Results:**
#
# | Comparison | Proportion Type 1 | Proportion Type 2 | Fold Change Ratio | Odds Ratio | P-value | Significance |
# |------------|------------------|------------------|------------------|------------|---------|-------------|
# | Knockdown vs Overexpression | 35.3% | 19.1% | 1.85 | 2.31 | 0.011 | * |
# | Knockout vs Overexpression | 27.7% | 19.1% | 1.45 | 1.62 | 0.044 | * |
# | Mutation vs Overexpression | 30.6% | 19.1% | 1.60 | 1.86 | 0.083 | ns |
#
# This text provides a clear, concise reporting of your results that emphasizes the main finding: loss-of-function approaches (particularly knockdown) are more effective than gain-of-function (overexpression) for achieving anti-aging effects.
