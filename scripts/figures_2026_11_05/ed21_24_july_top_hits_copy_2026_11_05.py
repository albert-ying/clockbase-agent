#!/usr/bin/env python3
"""
Figure 3b-e v2 — combined-modality top hits ranked by COMPOSITE SCORE:
   composite = -log10(pvalue) × scaled_final_score
   scaled_final_score = (final_score − min) / (max − min)  ∈ [0, 1]

For each of the 4 scoreable categories (Genetic, Drug, Disease, Environment),
rank all scored interventions from (mouse RNA, human RNA, human DNAm) together
by composite and render the top-N as horizontal bars.

Outputs:
  figures/top_hits_barplot_combined_modalities_v2.{pdf,png}
  result/top_hits_combined_modalities_v2.csv
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"

CATEGORIES = ["Genetic", "Drug", "Disease", "Environment"]
TOP_N = 15
# Cap on -log10(p) when forming the composite. Extreme-p outliers (e.g. sorted
# immune-cell comparisons with p < 1e-100) are retained but their bar length is
# bounded so they do not swamp the panel. The 95th percentile of -log10(p) across
# qualifying hits is ~20, so a cap of 30 leaves nearly all hits untouched.
CAP_NEGLOG10P = 30
MODALITY_MARKERS = {
    "mouse RNA-seq": "o",
    "human RNA-seq": "s",
    "human DNAm":    "^",
}
MODALITY_HATCH = {   # bar hatch patterns for print-safe secondary cue
    "mouse RNA-seq": "",
    "human RNA-seq": "////",
    "human DNAm":    "....",
}
COLOR_PRO  = "#D55E00"  # orange-red for pro-aging
COLOR_ANTI = "#0072B2"  # blue for anti-aging


def load_modality_scored():
    """Scored rows for all three modalities, including the extreme-effect-size /
    p-value flagged outliers, read from the same combined table that feeds
    Tables 1-2 (figures/combined_3modality_plot_data.csv) so the bar plot and the
    tables draw from the same data. Rows with a null outlier_status are unscored
    and dropped."""
    df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)
    df = df[df["outlier_status"].notna()].copy()
    keep = ["modality", "gse_id", "condition_category", "condition_name",
            "intervention", "log2FoldChange", "fdr", "final_score", "direction",
            "pvalue", "outlier_status"]
    return df[[c for c in keep if c in df.columns]]


def make_label(row):
    iv = row.get("intervention")
    if isinstance(iv, str) and iv.strip():
        return iv.strip()[:45]
    return str(row.get("condition_name", ""))[:45]


def main():
    all_ = load_modality_scored()
    print(f"loaded {len(all_):,} scored rows across modalities")

    # Filter: significant + final_score ≥ modality-specific floor
    thresholds = {"mouse RNA-seq": 100, "human RNA-seq": 100, "human DNAm": 80}
    qual = all_[all_["fdr"] < 0.05].copy()
    qual = qual[qual.apply(lambda r: r["final_score"] >= thresholds.get(r["modality"], 100), axis=1)]
    print(f"after FDR<0.05 + modality floor: {len(qual):,}")

    # ---- COMPOSITE SCORE = -log10(pvalue) × scaled_final_score [0,1] ----
    # Use pvalue if present; fall back to fdr if not.
    if "pvalue" in qual.columns and qual["pvalue"].notna().any():
        p = qual["pvalue"]
    else:
        p = qual["fdr"]
    qual["_neg_log10_p"] = np.minimum(-np.log10(p.clip(lower=1e-300)), CAP_NEGLOG10P)
    fs_min, fs_max = qual["final_score"].min(), qual["final_score"].max()
    if fs_max > fs_min:
        qual["_scaled_score"] = (qual["final_score"] - fs_min) / (fs_max - fs_min)
    else:
        qual["_scaled_score"] = 1.0
    qual["composite"] = qual["_neg_log10_p"] * qual["_scaled_score"]

    # Collect top-N per category, ranked by COMPOSITE
    picks = {}
    for cat in CATEGORIES:
        sub = qual[qual["condition_category"] == cat].copy()
        if sub.empty:
            picks[cat] = sub
            continue
        sub["display"] = sub.apply(make_label, axis=1)
        sub = sub.sort_values("composite", ascending=False).drop_duplicates(subset=["modality", "display"])
        picks[cat] = sub.head(TOP_N)

    # Figure: 2×2 of categories
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    for ax, cat in zip(axes.flat, CATEGORIES):
        sub = picks[cat]
        if sub.empty:
            ax.set_title(f"{cat}  (no scored cases)", fontsize=14)
            ax.axis("off")
            continue
        sub = sub.sort_values("composite", ascending=True)
        labels = [f"{r['display']}" for _, r in sub.iterrows()]
        values = sub["composite"].values
        colors = [COLOR_PRO if v > 0 else COLOR_ANTI for v in sub["log2FoldChange"].values]
        hatches = [MODALITY_HATCH[m] for m in sub["modality"].values]

        y_pos = np.arange(len(sub))
        bars = ax.barh(y_pos, values, color=colors, edgecolor="#111", linewidth=0.5)
        for bar, h in zip(bars, hatches):
            bar.set_hatch(h)

        # Modality marker at end of each bar
        x_pad = (values.max() - values.min()) * 0.02 if values.max() > values.min() else 0.5
        for y, (_, r) in zip(y_pos, sub.iterrows()):
            ax.plot(r["composite"] + x_pad, y,
                    marker=MODALITY_MARKERS[r["modality"]],
                    markersize=7, color="#222", markerfacecolor="#f7f7f7",
                    markeredgewidth=0.8, zorder=5)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=12)
        ax.set_xlabel("Score", fontsize=12)
        ax.set_title(f"{cat}", fontsize=16, fontweight="bold")
        ax.grid(axis="x", alpha=0.3, linestyle=":")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", labelsize=11)
        ax.tick_params(axis="x", labelsize=11)

    # Shared legend on the figure
    direction_handles = [
        Patch(facecolor=COLOR_ANTI, edgecolor="#111", label="anti-aging (log2FC < 0)"),
        Patch(facecolor=COLOR_PRO,  edgecolor="#111", label="pro-aging (log2FC > 0)"),
    ]
    modality_handles = [
        Line2D([0], [0], marker=MODALITY_MARKERS[m], color="#111", markerfacecolor="#eee",
               markersize=10, linestyle="None", label=m)
        for m in MODALITY_MARKERS
    ]
    fig.legend(handles=direction_handles + modality_handles,
               loc="upper center", ncol=5, frameon=False, fontsize=12,
               bbox_to_anchor=(0.5, 0.99))

    fig.suptitle(
        f"Top {TOP_N} interventions per category, ranked by composite (−log10 p × scaled score)",
        fontsize=16, fontweight="bold", y=0.945,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    Path("figures").mkdir(exist_ok=True)
    Path("result").mkdir(exist_ok=True)
    pdf = "figures/top_hits_barplot_combined_modalities_v2.pdf"
    png = "figures/top_hits_barplot_combined_modalities_v2.png"
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    # Save the picked table
    all_picks = pd.concat(
        [p.assign(category=c) for c, p in picks.items() if not p.empty],
        ignore_index=True,
    )
    all_picks.to_csv("result/top_hits_combined_modalities_v2.csv", index=False)
    print(f"\n✓ {pdf}")
    print(f"✓ {png}")
    print(f"✓ result/top_hits_combined_modalities_v2.csv")

    # Quick modality-coverage summary
    print("\n--- modality coverage across picked top-hits ---")
    print(all_picks.groupby(["category", "modality"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
