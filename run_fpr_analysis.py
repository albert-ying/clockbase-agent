#!/usr/bin/env python3
"""
Computational false-positive-rate (FPR) analysis for the ClockBase Agent.

Procedure (Option A, parametric Monte Carlo permutation):
    1. Load the agent's existing per-comparison results across the three
       modalities (mouse RNA-seq, human RNA-seq, human DNAm).
    2. For each Monte-Carlo iteration k = 1..K, simulate the global null
       (label-permutation null): under H0 of no true intervention effect,
       the t/Welch p-value of each comparison is Uniform(0,1).
    3. Apply the same BH-FDR correction the agent applies (within modality)
       and the same final_score >= SCORE_FLOOR plausibility filter the
       agent applies for "high-confidence" hits.
    4. Compute empirical FPR = mean fraction of comparisons that "pass"
       across permutations, with 95% CIs from K iterations.
    5. Repeat per intervention category and with/without q0 plausibility
       filter to show the filter's effect on empirical FPR.

Procedure (Option A2, real-sample permutation, validation subset):
    For the subset of mouse RNA-seq studies that have per-sample clock
    predictions saved to disk (n=25), run a true within-study label
    permutation: keep the same clock matrix, shuffle the test/control
    labels, recompute the Welch t-test p-value, apply BH-FDR, count
    rejections. Compare to the parametric estimate to confirm the
    parametric approximation is consistent.

Outputs:
    - figure_rev1/sensitivity_analyses/fpr_results.csv
    - figure_rev1/ExtendedDataFigure12_FPR_calibration.{png,pdf}

Constraints:
    - MPLBACKEND=Agg, pdf.fonttype=42 (Nature compliance).
    - np.random.seed(42) for reproducibility.
    - Only real data, no fabricated numbers.
"""

from __future__ import annotations

import json
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ----- Configuration -----
SEED = 42
K_PERMUTATIONS = 1000
ALPHAS = np.array([0.01, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20])
ALPHA_HEADLINE = 0.05

REPO = Path("/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base")
AGENT = REPO / "Aging-AI-Agent"
COMBINED_DATA = AGENT / "figures" / "combined_3modality_plot_data.csv"
OUT_TABLE = REPO / "figure_rev1" / "sensitivity_analyses" / "fpr_results.csv"
OUT_FIG_PNG = REPO / "figure_rev1" / "ExtendedDataFigure12_FPR_calibration.png"
OUT_FIG_PDF = REPO / "figure_rev1" / "ExtendedDataFigure12_FPR_calibration.pdf"

# Score floors used by the agent in create_combined_modality_volcano.py
SCORE_FLOOR = {
    "mouse RNA-seq": 100,
    "human RNA-seq": 100,
    "human DNAm": 80,
}
MODALITY_ORDER = ["mouse RNA-seq", "human RNA-seq", "human DNAm"]
MODALITY_SHORT = {
    "mouse RNA-seq": "Mouse RNA",
    "human RNA-seq": "Human RNA",
    "human DNAm": "Human DNAm",
}

# Categories the agent reports in main figures
MAIN_CATEGORIES = ["Genetic", "Drug", "Disease", "Environment"]

CATEGORY_COLORS = {
    "Genetic": "#228B22",
    "Drug": "#4B0082",
    "Disease": "#DC143C",
    "Environment": "#FF8C00",
}

MOD_COLORS = {
    "mouse RNA-seq": "#1f77b4",
    "human RNA-seq": "#d62728",
    "human DNAm": "#2ca02c",
}


def load_combined() -> pd.DataFrame:
    if not COMBINED_DATA.exists():
        sys.exit(f"Cannot proceed - required data file not found: {COMBINED_DATA}")
    df = pd.read_csv(COMBINED_DATA, low_memory=False)
    keep = ["modality", "gse_id", "condition_category", "fdr", "pvalue", "final_score"]
    df = df[keep].copy()
    df = df.dropna(subset=["modality", "gse_id"])
    df["final_score"] = df["final_score"].fillna(50.0)
    print(
        f"[load] combined: shape={df.shape}, modalities={df['modality'].value_counts().to_dict()}"
    )
    return df


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (statsmodels)."""
    pvals = np.asarray(pvals, dtype=float)
    if pvals.size == 0:
        return pvals
    _, fdr, _, _ = multipletests(pvals, alpha=0.05, method="fdr_bh")
    return fdr


def simulate_one_permutation(
    n: int, alphas: np.ndarray, score_pass: np.ndarray, rng: np.random.Generator
) -> dict:
    """
    Simulate one global-null permutation for n comparisons.

    Under H0 (random labels), the p-value of every t/Welch test is
    Uniform(0,1). We draw n iid Uniform(0,1) p-values, apply BH-FDR,
    and count rejections under each (alpha, q0_filter) configuration.

    Returns a dict mapping (alpha, q0_filter_on) -> rejection count.
    """
    p = rng.uniform(0.0, 1.0, size=n)
    fdr = bh_fdr(p)
    out = {}
    for alpha in alphas:
        rej = fdr < alpha
        out[(alpha, False)] = int(rej.sum())
        out[(alpha, True)] = int((rej & score_pass).sum())
    return out


def empirical_fpr(
    n: int,
    alphas: np.ndarray,
    score_pass: np.ndarray,
    K: int,
    rng: np.random.Generator,
) -> dict:
    """
    Run K permutations, return per-(alpha,q0_filter) array of FPRs and
    per-permutation rejection counts.

    Per-comparison FPR = # rejections / denom, where denom is:
        - q0_filter off : n  (all comparisons in the modality)
        - q0_filter on  : max(score_pass.sum(), 1)
    The function also returns the per-permutation rejection-count array
    (useful for FWER and for per-permutation interpretation since FPR
    fractions can be vanishingly small under the global null).
    """
    n_pass = int(score_pass.sum())
    fpr_results = {(a, q): np.zeros(K) for a in alphas for q in (False, True)}
    count_results = {(a, q): np.zeros(K, dtype=int) for a in alphas for q in (False, True)}
    fwer_results = {(a, q): np.zeros(K, dtype=int) for a in alphas for q in (False, True)}
    for k in range(K):
        out = simulate_one_permutation(n, alphas, score_pass, rng)
        for (a, q), count in out.items():
            denom = n if not q else max(n_pass, 1)
            fpr_results[(a, q)][k] = count / denom
            count_results[(a, q)][k] = count
            fwer_results[(a, q)][k] = int(count > 0)
    return fpr_results, count_results, fwer_results


def empirical_fpr_per_category(
    df_mod: pd.DataFrame,
    alphas: np.ndarray,
    K: int,
    rng: np.random.Generator,
) -> dict:
    """
    For each category, compute empirical FPR among comparisons in that
    category, where BH-FDR is applied across the entire modality (the
    agent's procedure). The category subset is selected after BH.

    Returns dict: (category, alpha, q0_filter) -> array length K.
    """
    n = len(df_mod)
    cat = df_mod["condition_category"].values
    score_pass = df_mod["score_pass"].values
    cats = MAIN_CATEGORIES
    out = {(c, a, q): np.zeros(K) for c in cats for a in alphas for q in (False, True)}
    for k in range(K):
        p = rng.uniform(0.0, 1.0, size=n)
        fdr = bh_fdr(p)
        for a in alphas:
            rej = fdr < a
            for c in cats:
                m_cat = cat == c
                denom_no = max(int(m_cat.sum()), 1)
                denom_q0 = max(int((m_cat & score_pass).sum()), 1)
                out[(c, a, False)][k] = (rej & m_cat).sum() / denom_no
                out[(c, a, True)][k] = (rej & m_cat & score_pass).sum() / denom_q0
    return out


def confidence_interval(arr: np.ndarray) -> tuple[float, float, float]:
    """Mean and 95% CI from K Monte-Carlo iterations."""
    return float(arr.mean()), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def real_sample_validation_mouse(K: int, rng: np.random.Generator) -> pd.DataFrame:
    """
    Real-sample within-study permutation on mouse studies that saved
    per-sample clock predictions during the original agent run. For each
    study, identify a binary group column, shuffle group assignments,
    Welch-t-test on a primary clock column, BH across studies, count
    rejections at FDR<0.05 across K shuffles.
    """
    candidates = [
        AGENT / "result_mm_rna" / d / "merged_data.csv"
        for d in os.listdir(AGENT / "result_mm_rna")
        if (AGENT / "result_mm_rna" / d / "merged_data.csv").exists()
    ]
    candidates += [
        AGENT / "result_mm_rna" / "GSE65358" / "GSE65383_merged_data.csv",
        AGENT / "result_mm_rna" / "GSE47835" / "gse47835_merged_data.csv",
        AGENT / "result_mm_rna" / "GSE180951" / "merged_data_temp.csv",
    ]
    print(f"[validation] candidate per-sample CSVs: {len(candidates)}")

    # We pick the standard mouse hazard clock as the primary clock test
    primary_clock_candidates = [
        "Hazard_Mouse_All_All_scaleddiff_EN_tAge",
        "Hazard_Rodents_All_All_scaleddiff_EN_tAge",
        "Hazard_All_All_All_scaleddiff_EN_tAge",
    ]

    # Collect (clock_vector, group_vector_binary) per usable study
    studies = []
    for path in candidates:
        if not path.exists():
            continue
        try:
            sub = pd.read_csv(path)
        except Exception:
            continue
        clock_col = None
        for c in primary_clock_candidates:
            if c in sub.columns:
                clock_col = c
                break
        if clock_col is None:
            continue
        # Find a binary grouping column. Prefer 'is_control' if present.
        if "is_control" in sub.columns and sub["is_control"].nunique() == 2:
            group_col = "is_control"
        else:
            # Pick a column with exactly 2 unique values, balanced enough to
            # give >=2 in each group, that is NOT the clock column itself
            group_col = None
            for c in sub.columns:
                if c == clock_col:
                    continue
                vals = sub[c].dropna().unique()
                if len(vals) == 2:
                    counts = sub[c].value_counts()
                    if counts.min() >= 2:
                        group_col = c
                        break
            if group_col is None:
                continue
        v = sub[[clock_col, group_col]].dropna()
        if len(v) < 4:
            continue
        v["__grp"] = (v[group_col] == v[group_col].iloc[0]).astype(int)
        if v["__grp"].sum() < 2 or (1 - v["__grp"]).sum() < 2:
            continue
        studies.append({"path": str(path), "clock": v[clock_col].values, "grp": v["__grp"].values})

    print(f"[validation] usable mouse RNA studies with binary groups: {len(studies)}")
    if len(studies) < 5:
        return pd.DataFrame()

    rejections = np.zeros(K, dtype=int)
    n_studies = len(studies)
    for k in range(K):
        pvals = np.empty(n_studies)
        for i, s in enumerate(studies):
            grp = s["grp"].copy()
            rng.shuffle(grp)
            x1 = s["clock"][grp == 1]
            x0 = s["clock"][grp == 0]
            if len(x1) < 2 or len(x0) < 2:
                pvals[i] = 1.0
                continue
            try:
                t, p = stats.ttest_ind(x1, x0, equal_var=False)
            except Exception:
                p = 1.0
            pvals[i] = p if np.isfinite(p) else 1.0
        fdr = bh_fdr(pvals)
        rejections[k] = int((fdr < 0.05).sum())
    rate = rejections / n_studies
    rec = {
        "n_studies": n_studies,
        "K": K,
        "alpha": 0.05,
        "empirical_fpr": float(rate.mean()),
        "ci_low": float(np.percentile(rate, 2.5)),
        "ci_high": float(np.percentile(rate, 97.5)),
    }
    print(f"[validation] real-sample permutation FPR (mouse RNA, n_studies={n_studies}, K={K}): {rec}")
    return pd.DataFrame([rec])


def main():
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FIG_PNG.parent.mkdir(parents=True, exist_ok=True)

    np.random.seed(SEED)
    rng = np.random.default_rng(SEED)

    df = load_combined()

    # Annotate plausibility filter
    df["score_floor"] = df["modality"].map(SCORE_FLOOR).fillna(80)
    df["score_pass"] = df["final_score"] >= df["score_floor"]

    # ------ Per-modality permutation ------
    rows = []
    per_mod_curves = {}  # for calibration plot
    cat_curves = {}      # for category bars

    mod_seeds = {m: SEED + 1000 * (i + 1) for i, m in enumerate(MODALITY_ORDER)}
    for mod in MODALITY_ORDER:
        sub = df[df["modality"] == mod].reset_index(drop=True)
        n = len(sub)
        score_pass = sub["score_pass"].values
        n_pass = int(score_pass.sum())
        print(f"[perm] {mod}: n={n}, n_score_pass={n_pass}")

        # Subseed per modality for reproducibility (deterministic)
        rng_m = np.random.default_rng(mod_seeds[mod])

        # Empirical FPR across all alphas, with/without q0
        fpr_dict, count_dict, fwer_dict = empirical_fpr(n, ALPHAS, score_pass, K_PERMUTATIONS, rng_m)

        # Calibration curves: track FPR, mean rejection count, and FWER per alpha
        for q0_on in (False, True):
            curve = []
            for a in ALPHAS:
                arr = fpr_dict[(a, q0_on)]
                m, lo, hi = confidence_interval(arr)
                fwer_arr = fwer_dict[(a, q0_on)]
                fwer_mean, fwer_lo, fwer_hi = confidence_interval(fwer_arr)
                count_arr = count_dict[(a, q0_on)]
                count_mean, count_lo, count_hi = confidence_interval(count_arr)
                curve.append((a, m, lo, hi, fwer_mean, fwer_lo, fwer_hi, count_mean, count_lo, count_hi))
            per_mod_curves[(mod, q0_on)] = curve

        # Headline (alpha=0.05) per-modality summary
        for q0_on in (False, True):
            arr = fpr_dict[(ALPHA_HEADLINE, q0_on)]
            counts = count_dict[(ALPHA_HEADLINE, q0_on)]
            m, lo, hi = confidence_interval(arr)
            # FWER: P(at least one rejection per permutation)
            fwer = float((counts > 0).mean())
            mean_counts = float(counts.mean())
            rows.append(
                {
                    "modality": mod,
                    "category": "ALL",
                    "q0_filter": "on" if q0_on else "off",
                    "alpha": ALPHA_HEADLINE,
                    "K": K_PERMUTATIONS,
                    "n_comparisons": n,
                    "n_score_pass": n_pass,
                    "empirical_fpr": m,
                    "ci_low": lo,
                    "ci_high": hi,
                    "mean_rejections_per_perm": mean_counts,
                    "fwer_at_least_one_rejection": fwer,
                }
            )

        # Per-category empirical FPR
        cat_dict = empirical_fpr_per_category(sub, np.array([ALPHA_HEADLINE]), K_PERMUTATIONS, rng_m)
        for cat in MAIN_CATEGORIES:
            for q0_on in (False, True):
                arr = cat_dict[(cat, ALPHA_HEADLINE, q0_on)]
                m, lo, hi = confidence_interval(arr)
                n_cat = int((sub["condition_category"] == cat).sum())
                n_cat_pass = int(((sub["condition_category"] == cat) & (sub["score_pass"])).sum())
                rows.append(
                    {
                        "modality": mod,
                        "category": cat,
                        "q0_filter": "on" if q0_on else "off",
                        "alpha": ALPHA_HEADLINE,
                        "K": K_PERMUTATIONS,
                        "n_comparisons": n_cat,
                        "n_score_pass": n_cat_pass,
                        "empirical_fpr": m,
                        "ci_low": lo,
                        "ci_high": hi,
                        "mean_rejections_per_perm": np.nan,
                        "fwer_at_least_one_rejection": np.nan,
                    }
                )
                cat_curves[(mod, cat, q0_on)] = (m, lo, hi)

    fpr_table = pd.DataFrame(rows)
    fpr_table.to_csv(OUT_TABLE, index=False)
    print(f"[write] {OUT_TABLE}")
    print(fpr_table.to_string(index=False))

    # ------ Real-sample validation on mouse subset ------
    val_rng = np.random.default_rng(SEED + 7)
    val_df = real_sample_validation_mouse(K=200, rng=val_rng)
    if not val_df.empty:
        val_path = OUT_TABLE.parent / "fpr_real_sample_validation_mouse.csv"
        val_df.to_csv(val_path, index=False)
        print(f"[write] {val_path}")

    # ------ Figure ------
    make_figure(per_mod_curves, cat_curves, fpr_table, val_df)

    print("\n[done] FPR analysis complete.")


def make_figure(per_mod_curves, cat_curves, fpr_table, val_df):
    fig = plt.figure(figsize=(13.5, 10.5))
    gs = fig.add_gridspec(
        2, 3, hspace=0.40, wspace=0.30, left=0.07, right=0.985, top=0.93, bottom=0.07
    )

    # Panel a/b/c: FWER calibration per modality (the most interpretable metric
    # under the global null with thousands of comparisons; per-comparison FPR
    # is essentially zero because BH-FDR is conservative under H0_global).
    for i, mod in enumerate(MODALITY_ORDER):
        ax = fig.add_subplot(gs[0, i])
        for q0_on, ls, label_suffix in [(False, "-", "no q0 filter"), (True, "--", "q0 filter on")]:
            curve = per_mod_curves[(mod, q0_on)]
            xs = np.array([c[0] for c in curve])
            fwer_m = np.array([c[4] for c in curve])
            fwer_lo = np.array([c[5] for c in curve])
            fwer_hi = np.array([c[6] for c in curve])
            ax.plot(xs, fwer_m, ls, color=MOD_COLORS[mod], lw=2, label=label_suffix, marker="o", ms=4)
            ax.fill_between(xs, fwer_lo, fwer_hi, color=MOD_COLORS[mod], alpha=0.15)
        # y=x reference (perfect calibration of FWER ~ alpha)
        ax.plot([0, ALPHAS.max()], [0, ALPHAS.max()], color="gray", lw=1, ls=":", label="y = x (well-calibrated)")
        ax.set_xlabel("Nominal alpha (BH-FDR threshold)")
        ax.set_ylabel("Empirical FWER\n(prob. of >=1 false rejection)")
        ax.set_title(f"{MODALITY_SHORT[mod]} calibration", fontsize=12)
        ax.set_xlim(0, ALPHAS.max() * 1.05)
        ax.set_ylim(0, max(0.30, ALPHAS.max() * 1.4))
        ax.legend(fontsize=8, loc="upper left", frameon=True)

    # Panel d: per-category FPR bars (alpha=0.05, no q0)
    # Per-category empirical FPRs under the global null are ~1e-6 (orders of
    # magnitude below alpha), so we autoscale the y-axis to make the bars
    # visible and annotate each bar with its numeric value in scientific
    # notation. The alpha=0.05 reference line is shown in panel f only,
    # where its scale is appropriate.
    def _plot_category_panel(ax, q0_on: bool, panel_title: str, ylabel: str):
        width = 0.22
        x = np.arange(len(MAIN_CATEGORIES))
        bar_values = []  # for autoscale
        bar_handles = []  # (rect, value)
        for j, mod in enumerate(MODALITY_ORDER):
            ms = []
            errs_lo = []
            errs_hi = []
            for cat in MAIN_CATEGORIES:
                m, lo, hi = cat_curves[(mod, cat, q0_on)]
                ms.append(m)
                errs_lo.append(max(m - lo, 0.0))
                errs_hi.append(max(hi - m, 0.0))
            rects = ax.bar(
                x + (j - 1) * width,
                ms,
                width,
                yerr=[errs_lo, errs_hi],
                color=MOD_COLORS[mod],
                label=MODALITY_SHORT[mod],
                capsize=2,
                edgecolor="black",
                linewidth=0.5,
            )
            for r, v, eh in zip(rects, ms, errs_hi):
                bar_values.append(v + eh)
                bar_handles.append((r, v))

        # Autoscale: floor at 1e-6 so a panel of all-zeros still has visible axis
        ymax_data = max(bar_values) if bar_values else 0.0
        ymax = max(ymax_data * 1.6, 1.0e-6)
        ax.set_ylim(0, ymax)

        # Annotate each bar with its value. Use scientific notation for
        # nonzero tiny values; explicit "0" for true zeros so the reader
        # sees the bar rendered.
        for r, v in bar_handles:
            cx = r.get_x() + r.get_width() / 2.0
            if v == 0.0:
                ax.text(cx, ymax * 0.02, "0", ha="center", va="bottom",
                        fontsize=6, rotation=90, color="black")
            else:
                # value formatted like "1.0e-6"
                if v >= 1e-3:
                    txt = f"{v:.2e}"
                else:
                    # express as a x 10^-6 style compact
                    exp = int(np.floor(np.log10(v)))
                    coef = v / (10.0 ** exp)
                    txt = f"{coef:.1f}e{exp}"
                # Place label just above the bar (or above error cap)
                y_label = max(v, 0.0) + ymax * 0.04
                ax.text(cx, y_label, txt, ha="center", va="bottom",
                        fontsize=6, rotation=90, color="black")

        ax.set_xticks(x)
        ax.set_xticklabels(MAIN_CATEGORIES, fontsize=9)
        ax.set_ylabel(ylabel)
        ax.set_title(panel_title, fontsize=11)
        # Use scientific notation on tick labels
        ax.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        ax.legend(fontsize=7, frameon=True, loc="upper right")
        # Note that alpha=0.05 is far above the visible y-range
        ax.text(
            0.02,
            0.98,
            f"alpha=0.05 ref line\nis far above y-range\n(values ~ 1e-6)",
            transform=ax.transAxes,
            fontsize=7,
            va="top",
            ha="left",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                      edgecolor="gray", alpha=0.85),
        )

    ax_d = fig.add_subplot(gs[1, 0])
    _plot_category_panel(
        ax_d,
        q0_on=False,
        panel_title="Per-category FPR, alpha=0.05 (no q0 filter)",
        ylabel="Empirical per-comparison FPR\n(no q0 filter)",
    )

    ax_e = fig.add_subplot(gs[1, 1])
    _plot_category_panel(
        ax_e,
        q0_on=True,
        panel_title="Per-category FPR, alpha=0.05 (q0 filter on)",
        ylabel="Empirical per-comparison FPR\n(q0 filter on)",
    )

    # Panel f: headline FWER per modality (q0 off vs on) + real-sample validation
    # Previously this panel plotted empirical_fpr (~1e-6), which was invisible at
    # a 0-0.05 scale producing an essentially empty box. Now we plot
    # fwer_at_least_one_rejection (probability >= 1 false rejection per modality),
    # which is on the order of 0.01-0.06 and directly comparable to alpha=0.05.
    ax = fig.add_subplot(gs[1, 2])
    head = fpr_table[(fpr_table["category"] == "ALL")].copy()
    head["x"] = head["modality"].map({m: i for i, m in enumerate(MODALITY_ORDER)})
    width = 0.35
    bar_rects_with_vals = []  # collect for annotation
    for j, q0 in enumerate(["off", "on"]):
        sub = head[head["q0_filter"] == q0].sort_values("x")
        fwer_vals = sub["fwer_at_least_one_rejection"].values
        offset = (j - 0.5) * width
        rects = ax.bar(
            sub["x"] + offset,
            fwer_vals,
            width,
            label=f"q0 filter {q0}",
            color="#888888" if q0 == "off" else "#222222",
            capsize=3,
            edgecolor="black",
            linewidth=0.5,
        )
        for r, v in zip(rects, fwer_vals):
            bar_rects_with_vals.append((r, v))
    ax.axhline(0.05, color="red", lw=1, ls="--", label="alpha = 0.05")
    # Annotate each bar with percentage label (e.g. "4.4%")
    ymax_bars = max((v for _, v in bar_rects_with_vals), default=0.0)
    ymax_f = max(ymax_bars * 1.5, 0.08)
    ax.set_ylim(0, ymax_f)
    for r, v in bar_rects_with_vals:
        cx = r.get_x() + r.get_width() / 2.0
        label_txt = f"{v*100:.1f}%" if v > 0 else "0%"
        y_label = v + ymax_f * 0.02
        ax.text(cx, y_label, label_txt, ha="center", va="bottom",
                fontsize=7, color="black")
    # Overlay mouse real-sample validation as a distinct marker with error bar
    if not val_df.empty:
        v = val_df.iloc[0]
        mouse_x = MODALITY_ORDER.index("mouse RNA-seq")
        val_fpr = float(v["empirical_fpr"])
        val_lo = float(v["ci_low"])
        val_hi = float(v["ci_high"])
        ax.errorbar(
            mouse_x,
            val_fpr,
            yerr=[[val_fpr - val_lo], [val_hi - val_fpr]],
            fmt="D",
            color="#e6550d",
            markersize=7,
            capsize=4,
            lw=1.5,
            label=f"Mouse real-sample (within-study)",
            zorder=5,
        )
    ax.set_xticks(np.arange(len(MODALITY_ORDER)))
    ax.set_xticklabels([MODALITY_SHORT[m] for m in MODALITY_ORDER], fontsize=9)
    ax.set_ylabel("Empirical FWER\n(prob. of >=1 false rejection)")
    title = "Headline FWER per modality (alpha=0.05)"
    if not val_df.empty:
        v = val_df.iloc[0]
        title += f"\nMouse real-sample perm: {v['empirical_fpr']:.3f} [{v['ci_low']:.3f}, {v['ci_high']:.3f}], n={int(v['n_studies'])}"
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, frameon=True)

    # suptitle removed: title text lives in the docx caption only
    fig.savefig(OUT_FIG_PNG, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_FIG_PDF, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {OUT_FIG_PNG}")
    print(f"[write] {OUT_FIG_PDF}")


if __name__ == "__main__":
    main()
