#!/usr/bin/env python3
"""
Pooled 3-modality scoring-weight sensitivity analysis.

Re-runs the per-component weight perturbation analysis from
sensitivity_results.csv (mouse RNA-seq only) on the POOLED atlas
(mouse RNA-seq + human RNA-seq + human DNAm) used in Figure 3.

Method:
  - Load the pooled 3-modality table (combined_3modality_plot_data.csv).
  - Apply FDR<0.05 and modality-specific final_score floors
    (mouse=100, hs_rna=100, hs_dnam=80) to match the manuscript filter.
  - Join q-component scores per modality from:
      * mouse RNA-seq:   data/scores_summary.csv (key: case_id == dataset)
      * human RNA-seq:   result/main_results_hs_rna_with_scores.csv (key: gse_id+comparison_name)
      * human DNAm:      result/main_results_hs_dnam_with_scores.csv (key: gse_id+comparison_name)
  - Baseline: composite = -log10(pvalue) * scaled_final_score, where
              scaled_final_score = (qsum - min)/(max - min) over the
              pooled qualifying set and qsum is the sum of the 8 q-components.
              (qsum is essentially identical to final_score; using qsum keeps the
              perturbation algebra clean.)
  - For each q-component c in {q0_plausibility, q1_experimental_model,
       q2_1_aging_pathways, q2_2_age_diseases, q3_1_experimental_rigor,
       q4_1_translational_potential, q5_1_literature_saturation,
       q5_2_mechanistic_novelty} and each factor f in {0.5, 0.75, 1.5, 2.0}:
       perturbed final = qsum - c + f*c; re-scale; recompute composite.
  - Report Spearman rho vs baseline (over the qualifying pool) and Jaccard
    top-100 overlap for anti-aging (log2FC<0) and pro-aging (log2FC>0).

Outputs:
  figure_rev1/sensitivity_analyses/sensitivity_results_pooled_3modality.csv
  figure_rev1/ExtendedDataFigure9_scoring_sensitivity_pooled.{png,pdf}
  figure_rev1/sensitivity_analyses/sensitivity_pooled_summary.md
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"

ROOT = Path("/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base")
POOLED_CSV = ROOT / "Aging-AI-Agent/figures/combined_3modality_plot_data.csv"
MM_SCORES  = ROOT / "Aging-AI-Agent/data/scores_summary.csv"
HR_SCORES  = ROOT / "Aging-AI-Agent/result/main_results_hs_rna_with_scores.csv"
HD_SCORES  = ROOT / "Aging-AI-Agent/result/main_results_hs_dnam_with_scores.csv"

OUT_DIR    = ROOT / "figure_rev1/sensitivity_analyses"
RESULTS_OUT = OUT_DIR / "sensitivity_results_pooled_3modality.csv"
SUMMARY_OUT = OUT_DIR / "sensitivity_pooled_summary.md"
FIG_PNG     = ROOT / "figure_rev1/ExtendedDataFigure9_scoring_sensitivity_pooled.png"
FIG_PDF     = ROOT / "figure_rev1/ExtendedDataFigure9_scoring_sensitivity_pooled.pdf"

Q_COLS = [
    "q0_plausibility",
    "q1_experimental_model",
    "q2_1_aging_pathways",
    "q2_2_age_diseases",
    "q3_1_experimental_rigor",
    "q4_1_translational_potential",
    "q5_1_literature_saturation",
    "q5_2_mechanistic_novelty",
]
FACTORS = [0.5, 0.75, 1.5, 2.0]
FLOORS = {"mouse RNA-seq": 100, "human RNA-seq": 100, "human DNAm": 80}
TOP_N = 100


def load_pooled_with_qcomponents():
    pooled = pd.read_csv(POOLED_CSV, low_memory=False)
    pooled = pooled[pooled["fdr"] < 0.05].copy()
    pooled = pooled[pooled.apply(
        lambda r: r["final_score"] >= FLOORS.get(r["modality"], 100), axis=1
    )].reset_index(drop=True)

    mm_q = pd.read_csv(MM_SCORES)[["dataset"] + Q_COLS].rename(columns={"dataset": "case_id"})
    hr_q = pd.read_csv(HR_SCORES, low_memory=False)[
        ["gse_id", "comparison_name"] + Q_COLS
    ].drop_duplicates(subset=["gse_id", "comparison_name"], keep="first")
    hd_q = pd.read_csv(HD_SCORES, low_memory=False)[
        ["gse_id", "comparison_name"] + Q_COLS
    ].drop_duplicates(subset=["gse_id", "comparison_name"], keep="first")

    mm = pooled[pooled["modality"] == "mouse RNA-seq"].merge(mm_q, on="case_id", how="left")
    hr = pooled[pooled["modality"] == "human RNA-seq"].merge(
        hr_q, on=["gse_id", "comparison_name"], how="left"
    )
    hd = pooled[pooled["modality"] == "human DNAm"].merge(
        hd_q, on=["gse_id", "comparison_name"], how="left"
    )

    merged = pd.concat([mm, hr, hd], ignore_index=True)
    n_total = len(merged)
    n_q_complete = merged[Q_COLS].notna().all(axis=1).sum()
    print(f"pooled qualifying rows: {n_total}, with all 8 q-components: {n_q_complete}")
    for mod, sub in merged.groupby("modality"):
        n_ok = sub[Q_COLS].notna().all(axis=1).sum()
        print(f"  {mod}: {len(sub)} rows ({n_ok} with q-components)")
    merged = merged[merged[Q_COLS].notna().all(axis=1)].reset_index(drop=True)
    print(f"after dropping rows missing q-components: {len(merged)}")
    return merged


def composite_score(df, qsum_col="qsum"):
    """composite = -log10(p) * scaled qsum (min-max [0,1])."""
    p = df["pvalue"].clip(lower=1e-300)
    neg_log10 = -np.log10(p)
    fs = df[qsum_col]
    fs_min, fs_max = fs.min(), fs.max()
    if fs_max > fs_min:
        scaled = (fs - fs_min) / (fs_max - fs_min)
    else:
        scaled = pd.Series(1.0, index=df.index)
    return neg_log10 * scaled


def jaccard(a, b):
    a, b = set(a), set(b)
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def top_ids_by_direction(df, comp_col, direction, top_n=TOP_N):
    if direction == "anti":
        sub = df[df["log2FoldChange"] < 0]
    else:
        sub = df[df["log2FoldChange"] > 0]
    return sub.sort_values(comp_col, ascending=False).head(top_n).index.tolist()


def run_sensitivity(df):
    df = df.copy()
    df["qsum"] = df[Q_COLS].sum(axis=1)
    df["composite_baseline"] = composite_score(df, "qsum")

    base_anti = top_ids_by_direction(df, "composite_baseline", "anti")
    base_pro  = top_ids_by_direction(df, "composite_baseline", "pro")

    rows = []
    for c in Q_COLS:
        for f in FACTORS:
            qsum_p = df["qsum"] - df[c] + f * df[c]
            df["_qsum_p"] = qsum_p
            df["_comp_p"] = composite_score(df, "_qsum_p")

            rho, _ = spearmanr(df["composite_baseline"], df["_comp_p"])
            anti_p = top_ids_by_direction(df, "_comp_p", "anti")
            pro_p  = top_ids_by_direction(df, "_comp_p", "pro")
            rows.append({
                "component": c,
                "perturbation_factor": f,
                "spearman_rho": rho,
                "jaccard_top100_anti": jaccard(base_anti, anti_p),
                "jaccard_top100_pro":  jaccard(base_pro,  pro_p),
            })
    return pd.DataFrame(rows), df


def plot_heatmap(res, n_rows, out_png, out_pdf):
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    metrics = [
        ("spearman_rho", "Spearman rho (rank stability)", "RdBu_r", (0.97, 1.0)),
        ("jaccard_top100_anti", "Jaccard top-100 anti-aging", "RdBu_r", (0.0, 1.0)),
        ("jaccard_top100_pro",  "Jaccard top-100 pro-aging",  "RdBu_r", (0.0, 1.0)),
    ]
    for ax, (col, title, cmap, vrng) in zip(axes, metrics):
        mat = res.pivot(index="component", columns="perturbation_factor", values=col)
        mat = mat.reindex(Q_COLS)
        mat = mat[FACTORS]
        im = ax.imshow(mat.values, aspect="auto", cmap=cmap, vmin=vrng[0], vmax=vrng[1])
        ax.set_xticks(range(len(FACTORS)))
        ax.set_xticklabels([f"x{f}" for f in FACTORS], fontsize=11)
        ax.set_yticks(range(len(Q_COLS)))
        ax.set_yticklabels(Q_COLS, fontsize=11)
        ax.set_xlabel("perturbation factor", fontsize=12)
        ax.set_title(title, fontsize=13, fontweight="bold")
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat.values[i, j]
                text = f"{v:.4f}" if col == "spearman_rho" else f"{v:.2f}"
                if col == "spearman_rho":
                    # midpoint at 0.985; use white text outside the band
                    color = "white" if (v < 0.985 or v > 0.998) else "black"
                else:
                    color = "white" if (v < 0.25 or v > 0.75) else "black"
                ax.text(j, i, text, ha="center", va="center", fontsize=10, color=color)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    # suptitle removed: title text lives in the docx caption only
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out_png}")
    print(f"saved {out_pdf}")


def write_summary(res, n_rows, modality_counts, out_md):
    med_rho = res["spearman_rho"].median()
    min_rho = res["spearman_rho"].min()
    med_anti = res["jaccard_top100_anti"].median()
    med_pro  = res["jaccard_top100_pro"].median()
    min_anti = res["jaccard_top100_anti"].min()
    min_pro  = res["jaccard_top100_pro"].min()

    comp_rank = (
        res.groupby("component")
           .agg(median_rho=("spearman_rho", "median"),
                min_jacc_anti=("jaccard_top100_anti", "min"),
                min_jacc_pro=("jaccard_top100_pro", "min"))
           .sort_values("min_jacc_anti")
    )
    most_sens = comp_rank.index[0]
    most_sens_jacc = comp_rank.iloc[0]["min_jacc_anti"]

    mod_str = ", ".join(f"{m}: {n}" for m, n in sorted(modality_counts.items()))

    table_lines = ["| Component | Median rho | Min Jaccard anti | Min Jaccard pro |",
                   "|---|---|---|---|"]
    for c, row in comp_rank.iterrows():
        table_lines.append(
            f"| {c} | {row['median_rho']:.4f} | {row['min_jacc_anti']:.3f} | {row['min_jacc_pro']:.3f} |"
        )
    table = "\n".join(table_lines)

    md = f"""# Sensitivity analysis: ClockBase Agent scoring weights (pooled 3-modality)

Atlas-wide ranking robustness. We re-ran the
mouse-only sensitivity analysis (ED Fig 9) on the POOLED three-modality
candidate set used in Figure 3, perturbing each of 8 q-component weights by
factors {{0.5, 0.75, 1.5, 2.0}} (one component at a time), recomputing the
composite ranking (-log10(p) x min-max scaled qsum, where qsum is the sum of
the 8 q-components) over the qualifying pool (FDR<0.05 with modality-specific
final_score floors: mouse RNA=100, human RNA=100, human DNAm=80), and
quantifying deviation from baseline via Spearman rank correlation and Jaccard
top-100 overlap for anti-aging (log2FC<0) and pro-aging (log2FC>0) lists.

## Scope

- Pooled qualifying interventions: **n = {n_rows}** ({mod_str})
- Perturbations: 8 q-components x 4 factors = 32 conditions

## Headline

Rankings are stable across the perturbation envelope at atlas scale:
median Spearman rho = **{med_rho:.4f}** (min {min_rho:.4f}),
median top-100 Jaccard = **{med_anti:.3f}** (anti) / **{med_pro:.3f}** (pro),
minima {min_anti:.3f} (anti) / {min_pro:.3f} (pro). The most sensitive
component is **{most_sens}** (min Jaccard anti = {most_sens_jacc:.3f} at the
strongest perturbation tested).

## Components ranked by impact (most sensitive first; lowest min Jaccard anti)

{table}

## Manuscript-ready paragraph

To verify that the prioritised interventions in the pooled three-modality atlas
do not depend on the specific q-component weights chosen by the ClockBase Agent
scoring scheme, we perturbed each of the eight component weights by factors of
0.5, 0.75, 1.5 and 2.0 (one component at a time) and recomputed the composite
ranking (-log10(p) multiplied by the min-max scaled q-component sum) over the
qualifying pool of n = {n_rows} significant interventions (FDR<0.05 with
modality-specific final_score floors). Across the 32 perturbations, the median
Spearman rank correlation against the baseline ranking was {med_rho:.4f}
(minimum {min_rho:.4f}), and the median Jaccard overlap among the top 100
prioritised anti-aging and pro-aging hits was {med_anti:.3f} and {med_pro:.3f}
respectively (minima {min_anti:.3f} and {min_pro:.3f}). The most sensitive
component was {most_sens}; even at its strongest perturbation, the top-100
anti-aging set retained a Jaccard of {most_sens_jacc:.3f} with the baseline.
These results indicate that the headline atlas-wide rankings are not artefacts
of the specific weighting choices.
"""
    out_md.write_text(md)
    print(f"saved {out_md}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    merged = load_pooled_with_qcomponents()
    modality_counts = merged["modality"].value_counts().to_dict()
    res, _ = run_sensitivity(merged)
    res.to_csv(RESULTS_OUT, index=False)
    print(f"saved {RESULTS_OUT}")
    plot_heatmap(res, len(merged), FIG_PNG, FIG_PDF)
    write_summary(res, len(merged), modality_counts, SUMMARY_OUT)
    print("\nSpearman rho summary:")
    print(res["spearman_rho"].describe())
    print("\nJaccard anti summary:")
    print(res["jaccard_top100_anti"].describe())
    print("\nJaccard pro summary:")
    print(res["jaccard_top100_pro"].describe())


if __name__ == "__main__":
    main()
