#!/usr/bin/env python3
"""
Figure 3b-e analog: top-N interventions per condition category, ranked by final_score.

Reads result/main_results_hs_{MODALITY}_with_scores.csv (output of
s03_merge_scores_human.py) and produces a 4-panel horizontal bar chart —
one panel per category (Genetic, Drug, Disease, Environment) — with bars
colored by log2FoldChange (red = pro-aging, blue = anti-aging).

Flip MODALITY and rerun for the other modality.
"""

# ---- HUMAN CONFIG --------------------------------------------------
MODALITY = "rna"  # "rna" or "dnam"
METHOD_NAME = f"hs_{MODALITY}"
SUFFIX = f"_{METHOD_NAME}"
TOP_N = 15  # bars per panel
CATEGORIES = ["Genetic", "Drug", "Disease", "Environment"]
# --------------------------------------------------------------------

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"


def make_label(row):
    """Display label: prefer intervention (LLM-generated), fallback to condition_name."""
    iv = row.get("intervention")
    if isinstance(iv, str) and iv.strip():
        return iv.strip()
    return str(row.get("condition_name", ""))[:60]


def plot_panel(ax, sub, category, cmap):
    if sub.empty:
        ax.set_title(f"{category}  (no scored cases)")
        ax.axis("off")
        return
    sub = sub.sort_values("final_score", ascending=True).tail(TOP_N)
    labels = [make_label(r) for _, r in sub.iterrows()]
    values = sub["final_score"].values
    # Color by log2FoldChange sign + magnitude
    log2fc = sub["log2FoldChange"].values
    max_abs = max(np.abs(log2fc).max(), 1e-6)
    norm = plt.Normalize(vmin=-max_abs, vmax=max_abs)
    colors = [cmap(norm(v)) for v in log2fc]

    bars = ax.barh(labels, values, color=colors, edgecolor="#333", linewidth=0.3)

    # Annotate log2FC inside each bar
    for bar, v in zip(bars, log2fc):
        ax.text(
            bar.get_width() * 0.5, bar.get_y() + bar.get_height() / 2,
            f"log2FC={v:+.2f}",
            ha="center", va="center", fontsize=7, color="white", alpha=0.9,
        )

    ax.set_title(f"{category}  ({len(sub)} top of {len(_raw_by_cat.get(category, []))})", fontsize=11)
    ax.set_xlabel("Final score (0–100 + Q0 bonus)")
    ax.grid(axis="x", alpha=0.3, linestyle=":")
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    in_path = f"result/main_results_{METHOD_NAME}_with_scores.csv"
    print(f"loading {in_path}")
    df = pd.read_csv(in_path, low_memory=False)
    print(f"  {len(df)} scored rows, {df['condition_category'].nunique()} categories")
    print(df["condition_category"].value_counts().head(8))

    global _raw_by_cat
    _raw_by_cat = {c: df[df["condition_category"] == c] for c in CATEGORIES}

    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    cmap = plt.cm.RdBu_r  # red for +, blue for − (so anti-aging = blue)
    for ax, cat in zip(axes.flat, CATEGORIES):
        plot_panel(ax, _raw_by_cat[cat], cat, cmap)

    fig.suptitle(
        f"Top {TOP_N} scored interventions per category  —  {METHOD_NAME.upper()}",
        fontsize=13, fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    Path("figures").mkdir(exist_ok=True)
    pdf = f"figures/top_hits_barplot_by_category{SUFFIX}.pdf"
    png = f"figures/top_hits_barplot_by_category{SUFFIX}.png"
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    print(f"\n✓ saved {pdf} and {png}")

    # Also save the top-N tables as CSV for inspection
    for cat in CATEGORIES:
        sub = _raw_by_cat[cat].sort_values("final_score", ascending=False).head(TOP_N).copy()
        sub["display_label"] = sub.apply(make_label, axis=1)
        cols = [c for c in ["case_id", "gse_id", "display_label", "intervention",
                            "condition_name", "log2FoldChange", "fdr",
                            "final_score"] if c in sub.columns]
        out_csv = f"result/top_hits_{cat.lower()}{SUFFIX}.csv"
        sub[cols].to_csv(out_csv, index=False)
        print(f"  {cat:<12}: {out_csv}")


if __name__ == "__main__":
    main()
