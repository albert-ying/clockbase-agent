#!/usr/bin/env python3
"""
Figure 3a equivalent for human: effect-size vs. −log10(p) scatter, ALL comparisons,
points sized by composite final_score, colored by intervention category, with text
annotations on top-scoring interventions at extreme positions.

Mirrors the mouse Fig 3 caption from clockbase_agent_v3.qmd:
  - x = log2FoldChange
  - y = −log10(p-value) (capped for readability)
  - point size ∝ final_score (default minimum when unscored)
  - color by category (Genetic, Drug, Environment, Disease, Other)
  - horizontal dashed line at p=0.05
  - text annotations labeling top-N per category by final_score * |log2FC|

Flip MODALITY and rerun for each modality.
"""

# ---- HUMAN CONFIG --------------------------------------------------
MODALITY = "rna"  # "rna" or "dnam"
METHOD_NAME = f"hs_{MODALITY}"
SUFFIX = f"_{METHOD_NAME}"

TOP_ANNOT_PER_SIDE = 12   # # labels for pro-aging side, # for anti-aging
TOP_ANNOT_PER_CATEGORY = 3  # within each category, try this many labels
# --------------------------------------------------------------------

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
from adjustText import adjust_text

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"

CATEGORY_ORDER = ["Genetic", "Drug", "Environment", "Disease", "Other"]
# Nature-style category palette (matches mouse Fig 3 conventions)
COLORS = {
    "Genetic":     "#D62728",  # red
    "Drug":        "#1F77B4",  # blue
    "Environment": "#FF7F0E",  # orange
    "Disease":     "#9467BD",  # purple
    "Other":       "#7F7F7F",  # gray
}


def _canonicalize(df):
    _CANON = set(CATEGORY_ORDER)
    df = df.copy()
    df["condition_category"] = df["condition_category"].where(
        df["condition_category"].isin(_CANON), "Other"
    )
    return df


def _safe_label(row):
    iv = row.get("intervention")
    if isinstance(iv, str) and iv.strip():
        return iv.strip()[:45]
    return str(row.get("condition_name", ""))[:45]


def main():
    master_path = Path(f"result/main_results_{METHOD_NAME}_outlier_marked.csv")
    scores_path = Path(f"result/main_results_{METHOD_NAME}_with_scores.csv")

    print(f"loading master: {master_path}")
    master = pd.read_csv(master_path, low_memory=False)
    master = master[master["outlier_status"] == "Normal"].copy()
    # DNAm master has 5 clock rows per comparison; collapse to the primary clock
    # (GrimAgeV2) so each comparison is ONE dot, matching the mouse Fig 3 convention.
    if MODALITY == "dnam":
        master = master[master["clock"] == "GrimAgeV2"].copy()
    master = _canonicalize(master)
    print(f"  {len(master):,} non-outlier rows")

    print(f"loading scores: {scores_path}")
    scored = pd.read_csv(scores_path, low_memory=False)
    print(f"  {len(scored):,} scored rows")

    # Merge final_score + intervention onto all comparisons, keyed on
    # (gse_id, comparison_name, condition_name). Unscored rows get NaN.
    # Master already has final_score/intervention columns (defaulted to 50 by s01);
    # drop those before merge so the real scored values land cleanly.
    key = ["gse_id", "comparison_name", "condition_name"]
    master = master.drop(columns=[c for c in ["final_score", "intervention"] if c in master.columns])
    scored_slim = scored[key + ["final_score", "intervention"]].drop_duplicates(subset=key)
    df = master.merge(scored_slim, on=key, how="left")
    scored_n = df["final_score"].notna().sum()
    print(f"  merged: {scored_n:,} rows with final_score, {len(df) - scored_n:,} without")

    # y-axis transform: −log10(p). Cap for readability.
    df["_nlp"] = -np.log10(df["pvalue"].clip(lower=1e-300))

    # Point size mapping: floor at 4, ceil at 80, scaled by final_score
    fs = df["final_score"].fillna(40.0).clip(lower=-30, upper=140)
    df["_size"] = 4 + ((fs + 30) / 170) ** 2 * 76   # cubic-ish scaling

    # Alpha for unscored vs scored
    df["_alpha"] = np.where(df["final_score"].notna(), 0.8, 0.15)

    # --- figure ---
    fig, ax = plt.subplots(figsize=(12, 9))

    # Plot per-category so zorder stacks reasonably (big cats first, Genetic on top)
    for cat in ["Other", "Disease", "Environment", "Drug", "Genetic"]:
        sub = df[df["condition_category"] == cat]
        if sub.empty:
            continue
        ax.scatter(
            sub["log2FoldChange"], sub["_nlp"],
            s=sub["_size"], c=COLORS[cat],
            alpha=sub["_alpha"].values,
            edgecolors="none", rasterized=True,
            label=f"{cat} (n={len(sub):,}, {100*sub.get('fdr', pd.Series(dtype=float)).lt(0.05).sum()/max(1,len(sub)):.1f}% sig)",
        )

    # Significance line at p=0.05
    ax.axhline(-np.log10(0.05), ls="--", color="#666", lw=0.8, zorder=0)
    ax.text(
        ax.get_xlim()[1], -np.log10(0.05),
        "p = 0.05",
        fontsize=9, color="#666", ha="right", va="bottom",
    )

    # --- annotations: top-scoring interventions per side (pro vs anti-aging) ---
    texts = []
    ann_source = df[df["final_score"].notna() & (df["fdr"] < 0.05)].copy()
    ann_source["_rank_score"] = ann_source["final_score"] * np.abs(ann_source["log2FoldChange"])

    # Pick top-N per category per direction
    picks = []
    for cat in CATEGORY_ORDER:
        cat_rows = ann_source[ann_source["condition_category"] == cat]
        for side_mask, n in [
            (cat_rows["log2FoldChange"] < 0, TOP_ANNOT_PER_CATEGORY),
            (cat_rows["log2FoldChange"] > 0, TOP_ANNOT_PER_CATEGORY),
        ]:
            picks.append(cat_rows[side_mask].nlargest(n, "_rank_score"))
    ann_df = pd.concat(picks).drop_duplicates(subset=["gse_id", "condition_name"])
    # Also cap overall at 2 × TOP_ANNOT_PER_SIDE for clutter
    ann_df = ann_df.nlargest(2 * TOP_ANNOT_PER_SIDE, "_rank_score")
    print(f"  annotating {len(ann_df)} interventions")

    for _, row in ann_df.iterrows():
        label = _safe_label(row)
        t = ax.text(
            row["log2FoldChange"], row["_nlp"], label,
            fontsize=7.5, color="#111",
            ha="center", va="center",
            path_effects=[pe.withStroke(linewidth=2.2, foreground="white")],
        )
        texts.append(t)

    # adjust_text to avoid overlap
    try:
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="#888", lw=0.4, alpha=0.6),
            expand_points=(1.2, 1.4),
            force_points=0.5,
            force_text=0.8,
        )
    except Exception as e:
        print(f"  (adjust_text skipped: {e})")

    # Axis cosmetics
    modality_label = {"rna": "Human RNA-seq", "dnam": "Human DNA methylation (GrimAgeV2 + 4 supporting clocks)"}[MODALITY]
    total_comparisons = len(df)
    total_datasets = df["gse_id"].nunique()
    sig_count = (df["fdr"] < 0.05).sum()

    ax.set_xlabel("log2 effect size (β)\n← decelerating aging     accelerating aging →", fontsize=11)
    ax.set_ylabel("−log10(p-value)", fontsize=11)
    ax.set_title(
        f"{modality_label}: AI-agent-identified age-modifying interventions\n"
        f"{total_comparisons:,} intervention–control comparisons from {total_datasets:,} datasets  |  "
        f"significant (FDR<0.05): {sig_count:,} ({100*sig_count/total_comparisons:.1f}%)",
        fontsize=11,
    )
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.2, linestyle=":")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    Path("figures").mkdir(exist_ok=True)
    pdf = f"figures/fig3a_volcano{SUFFIX}.pdf"
    png = f"figures/fig3a_volcano{SUFFIX}.png"
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    print(f"\n✓ saved {pdf} and {png}")


if __name__ == "__main__":
    main()
