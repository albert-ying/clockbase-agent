#!/usr/bin/env python3
"""
Figure 3b-e analog across ALL three modalities.

For each of the 4 scoreable categories (Genetic, Drug, Disease, Environment), rank
all scored interventions from (mouse RNA, human RNA, human DNAm) together by
final_score and render the top-N as horizontal bars.

Styling:
  - bar color        = direction (red = pro-aging, blue = anti-aging)
  - bar edge pattern = modality  (solid = mouse RNA, dashed = human RNA, dotted = human DNAm)
  - circle marker at bar end = modality (●/■/▲)  — redundant shape cue for print clarity

Outputs:
  figures/top_hits_barplot_combined_modalities.{pdf,png}
  result/top_hits_combined_modalities.csv   (table of what's plotted)
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
    """Pull one DataFrame per modality with [modality, intervention, condition_category,
    log2FoldChange, fdr, final_score]."""
    out = []

    # Human RNA
    df = pd.read_csv("result/main_results_hs_rna_with_scores.csv", low_memory=False)
    df["modality"] = "human RNA-seq"
    out.append(df)

    # Human DNAm — keep GrimAgeV2 only (primary clock) to match volcano convention
    df = pd.read_csv("result/main_results_hs_dnam_with_scores.csv", low_memory=False)
    if "clock" in df.columns:
        df = df[df["clock"] == "GrimAgeV2"].copy()
    df["modality"] = "human DNAm"
    out.append(df)

    # Mouse RNA — merge detailed_scores with scores_summary
    det = pd.read_csv("data/detailed_scores_result_mm_rna.csv", low_memory=False)
    scs = pd.read_csv("data/scores_summary.csv")[["dataset", "final_score"]]
    mm  = det.merge(scs, left_on="case_id", right_on="dataset", how="inner")
    mm["modality"] = "mouse RNA-seq"
    out.append(mm)

    # Align columns
    keep = ["modality", "gse_id", "condition_category", "condition_name",
            "intervention", "log2FoldChange", "fdr", "final_score", "direction"]
    out = [d[[c for c in keep if c in d.columns]] for d in out]
    all_ = pd.concat(out, ignore_index=True)
    return all_


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

    # Collect top-N per category, ranked by final_score
    picks = {}
    for cat in CATEGORIES:
        sub = qual[qual["condition_category"] == cat].copy()
        if sub.empty:
            picks[cat] = sub
            continue
        sub["display"] = sub.apply(make_label, axis=1)
        # Dedup by (modality, intervention) to avoid near-duplicates
        sub = sub.sort_values("final_score", ascending=False).drop_duplicates(subset=["modality", "display"])
        picks[cat] = sub.head(TOP_N)

    # Figure: 2×2 of categories
    fig, axes = plt.subplots(2, 2, figsize=(17, 11))
    for ax, cat in zip(axes.flat, CATEGORIES):
        sub = picks[cat]
        if sub.empty:
            ax.set_title(f"{cat}  (no scored cases)", fontsize=12)
            ax.axis("off")
            continue
        sub = sub.sort_values("final_score", ascending=True)
        labels = [f"{r['display']}" for _, r in sub.iterrows()]
        values = sub["final_score"].values
        colors = [COLOR_PRO if v > 0 else COLOR_ANTI for v in sub["log2FoldChange"].values]
        hatches = [MODALITY_HATCH[m] for m in sub["modality"].values]

        y_pos = np.arange(len(sub))
        bars = ax.barh(y_pos, values, color=colors, edgecolor="#111", linewidth=0.5)
        for bar, h in zip(bars, hatches):
            bar.set_hatch(h)

        # Modality marker at end of each bar
        for y, (_, r) in zip(y_pos, sub.iterrows()):
            ax.plot(r["final_score"] + 2, y,
                    marker=MODALITY_MARKERS[r["modality"]],
                    markersize=7, color="#222", markerfacecolor="#f7f7f7",
                    markeredgewidth=0.8, zorder=5)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Final score")
        ax.set_title(f"{cat}", fontsize=13, fontweight="bold")
        ax.grid(axis="x", alpha=0.3, linestyle=":")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", labelsize=8)

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
               loc="upper center", ncol=5, frameon=False, fontsize=10,
               bbox_to_anchor=(0.5, 0.99))

    fig.suptitle(
        f"Top {TOP_N} scored interventions per category, across three modalities",
        fontsize=14, fontweight="bold", y=0.945,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])

    Path("figures").mkdir(exist_ok=True)
    Path("result").mkdir(exist_ok=True)
    pdf = "figures/top_hits_barplot_combined_modalities.pdf"
    png = "figures/top_hits_barplot_combined_modalities.png"
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")

    # Save the picked table
    all_picks = pd.concat(
        [p.assign(category=c) for c, p in picks.items() if not p.empty],
        ignore_index=True,
    )
    all_picks.to_csv("result/top_hits_combined_modalities.csv", index=False)
    print(f"\n✓ {pdf}")
    print(f"✓ {png}")
    print(f"✓ result/top_hits_combined_modalities.csv")

    # Quick modality-coverage summary
    print("\n--- modality coverage across picked top-hits ---")
    print(all_picks.groupby(["category", "modality"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
