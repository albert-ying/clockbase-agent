#!/usr/bin/env python3
"""
Regenerate ED Fig 13 (cross-clock consistency of intervention rankings,
human DNA methylation) as a 5x5 heatmap pair including YingDamAge.

The atlas scoring pipeline evaluates each human DNAm comparison against
five clocks (GrimAgeV2 as primary + PhenoAge, Horvathv1, YingDamAge,
DunedinPACE as cross-clock supporters; see
Aging-AI-Agent/s02_prepare_scoring_input_human.py). The previous 4x4
ED Fig 13 excluded YingDamAge; this regen restores the missing fifth
clock so the figure matches the scoring pipeline.

Method (mirrors the manuscript composite formula):
  1. Load the canonical outlier-removed human DNAm analytic table.
  2. Pivot to wide-by-clock per (gse_id, comparison_name,
     condition_name); each comparison carries five (pvalue,
     effect_size) pairs.
  3. Restrict to comparisons that have non-null p-values across all five
     clocks (so all clock pairs are computed on the same set of rows).
  4. For each clock, composite_k = -log10(p_k) * normalised(final_score),
     where final_score is the agent's rubric score (shared across clocks
     per comparison; the score is rubric-derived and clock-agnostic) and
     the normalisation is min-max within the comparison set.
  5. Pairwise Spearman rank-correlation of composite_k across clocks.
  6. Pairwise Jaccard of the top-100 intervention candidates by composite.
  7. Render two side-by-side heatmaps and write the 5x5 CSVs.

Outputs (overwrite the previous 4x4 versions):
  figure_rev1/sensitivity_analyses/multiclock_spearman_dnam.csv
  figure_rev1/sensitivity_analyses/multiclock_jaccard_dnam.csv
  figure_rev1/ExtendedDataFigure13_multiclock_dnam.png
  figure_rev1/ExtendedDataFigure13_multiclock_dnam.pdf
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
DNAM_MASTER = REPO_ROOT / "Aging-AI-Agent/result/main_results_hs_dnam_outlier_marked.csv"
OUT_DIR = REPO_ROOT / "figure_rev1/sensitivity_analyses"
FIG_DIR = REPO_ROOT / "figure_rev1"

CLOCKS = ["GrimAgeV2", "PhenoAge", "Horvathv1", "YingDamAge", "DunedinPACE"]
# Display labels for the axes. "Horvathv1" is the BioLearn identifier for the
# Horvath 2013 multi-tissue clock; show it as "Horvath" to match the manuscript
# caption ("Horvath multi-tissue") and avoid the confusing "v1" suffix. The data
# filter still keys on the raw `clock` column value (CLOCKS), so this is a
# label-only change with no effect on the computed correlations.
CLOCK_LABELS = {
    "GrimAgeV2": "GrimAgeV2",
    "PhenoAge": "PhenoAge",
    "Horvathv1": "Horvath",
    "YingDamAge": "YingDamAge",
    "DunedinPACE": "DunedinPACE",
}
TICK_LABELS = [CLOCK_LABELS[c] for c in CLOCKS]
TOP_N = 25  # candidate set size for Jaccard (qualifying subset is small;
            # top-25 ~ top-half so the overlap test has real power)
SCORE_FLOOR = 80  # manuscript scoring floor for human DNAm
FDR_THRESHOLD = 0.05  # primary-clock FDR gate from s02_prepare_scoring_input_human.py


def load_dnam_wide() -> pd.DataFrame:
    df = pd.read_csv(DNAM_MASTER, low_memory=False)
    df = df[df["outlier_status"] == "Normal"]
    # final_score is rubric-derived and per-comparison (clock-agnostic).
    # Take the first non-null final_score per comparison key.
    key = ["gse_id", "comparison_name", "condition_name"]
    fs = (df.dropna(subset=["final_score"])
            .drop_duplicates(subset=key, keep="first")[key + ["final_score"]])
    sub = df[df["clock"].isin(CLOCKS)].copy()
    wide = sub.pivot_table(index=key,
                           columns="clock",
                           values=["pvalue", "effect_size"],
                           aggfunc="first")
    wide.columns = [f"{m}_{c}" for m, c in wide.columns]
    wide = wide.merge(fs, on=key, how="left")
    return wide


def composites(wide: pd.DataFrame) -> pd.DataFrame:
    """One column per clock holding -log10(p_clock) * normalised(final_score)."""
    fs = wide["final_score"].astype(float)
    # min-max normalise final_score among rows with data
    lo, hi = fs.min(), fs.max()
    fs_norm = (fs - lo) / (hi - lo) if hi > lo else fs * 0
    out = pd.DataFrame(index=wide.index)
    for clk in CLOCKS:
        p = wide.get(f"pvalue_{clk}").astype(float)
        nlp = -np.log10(p.clip(lower=1e-300))
        out[clk] = nlp * fs_norm
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"loading {DNAM_MASTER.name} ...")
    wide = load_dnam_wide()
    print(f"  {len(wide):,} unique comparisons across the 5 clocks")

    # Restrict to the qualifying scored subset that the manuscript uses for
    # the multi-clock comparison: comparisons that received a final_score at
    # or above the DNAm score floor (80), are FDR-significant on the primary
    # clock (GrimAgeV2 FDR<0.05; this is the gate from
    # s02_prepare_scoring_input_human.py), and have non-null pvalue across
    # ALL five clocks (so every pairwise comparison is computed on the same
    # set of rows).
    pv_cols = [f"pvalue_{c}" for c in CLOCKS]
    keep = wide.dropna(subset=pv_cols + ["final_score"])
    # Primary-clock FDR gate via -log10(p) > -log10(0.05) on GrimAgeV2 fdr column.
    # The master table also carries an fdr column per (comparison, clock); we use
    # the GrimAgeV2 row's fdr by joining it in.
    fdr_g = (pd.read_csv(DNAM_MASTER, low_memory=False)
                .query("outlier_status == 'Normal' and clock == 'GrimAgeV2'")
                [["gse_id", "comparison_name", "condition_name", "fdr"]]
                .rename(columns={"fdr": "fdr_GrimAgeV2"}))
    keep = keep.merge(fdr_g, on=["gse_id", "comparison_name", "condition_name"], how="left")
    keep = keep[(keep["fdr_GrimAgeV2"] < FDR_THRESHOLD) &
                (keep["final_score"] >= SCORE_FLOOR)]
    print(f"  {len(keep):,} comparisons in qualifying subset "
          f"(GrimAgeV2 FDR<{FDR_THRESHOLD}, final_score>={SCORE_FLOOR}, all-5-clocks non-null)")

    comp = composites(keep)
    # Drop any rows where composite is undefined for at least one clock
    comp = comp.dropna()
    print(f"  {len(comp):,} comparisons with well-defined composites in every clock")

    # Spearman pairwise
    spear = pd.DataFrame(index=CLOCKS, columns=CLOCKS, dtype=float)
    for i, a in enumerate(CLOCKS):
        for b in CLOCKS:
            spear.loc[a, b] = spearmanr(comp[a], comp[b]).correlation

    # Jaccard top-N pairwise
    tops = {c: set(comp[c].nlargest(TOP_N).index) for c in CLOCKS}
    jacc = pd.DataFrame(index=CLOCKS, columns=CLOCKS, dtype=float)
    for a in CLOCKS:
        for b in CLOCKS:
            inter = len(tops[a] & tops[b])
            union = len(tops[a] | tops[b])
            jacc.loc[a, b] = inter / union if union else float("nan")

    spear.to_csv(OUT_DIR / "multiclock_spearman_dnam.csv")
    jacc.to_csv(OUT_DIR / "multiclock_jaccard_dnam.csv")
    print(f"saved {OUT_DIR.name}/multiclock_spearman_dnam.csv  (5x5)")
    print(f"saved {OUT_DIR.name}/multiclock_jaccard_dnam.csv  (5x5)")

    # Summaries used by the figure caption (off-diagonal stats)
    off = lambda m: m.values[np.triu_indices(len(CLOCKS), k=1)]
    print(f"\nSpearman off-diagonal: median {np.median(off(spear)):.2f} "
          f"(range {off(spear).min():.2f}--{off(spear).max():.2f})")
    print(f"Jaccard  off-diagonal: median {np.median(off(jacc)):.2f} "
          f"(range {off(jacc).min():.2f}--{off(jacc).max():.2f})")

    # Render
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    for ax, m, title in [
        (axes[0], spear, "a. Spearman rank-correlation\nof composite scores"),
        (axes[1], jacc,  f"b. Jaccard overlap of top-{TOP_N}\nintervention candidates"),  # top-N adjusted to qualifying-subset size
    ]:
        im = ax.imshow(m.values, cmap="RdYlBu_r", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(CLOCKS))); ax.set_yticks(range(len(CLOCKS)))
        ax.set_xticklabels(TICK_LABELS, rotation=35, ha="right")
        ax.set_yticklabels(TICK_LABELS)
        ax.set_title(title)
        for i in range(len(CLOCKS)):
            for j in range(len(CLOCKS)):
                v = m.values[i, j]
                txt_color = "white" if v < 0.55 else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=10, color=txt_color)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    plt.tight_layout()
    for ext in ("png", "pdf"):
        out = FIG_DIR / f"ExtendedDataFigure13_multiclock_dnam.{ext}"
        plt.savefig(out, dpi=220, bbox_inches="tight")
        print(f"saved {out.name}")


if __name__ == "__main__":
    main()
