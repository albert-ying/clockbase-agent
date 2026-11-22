#!/usr/bin/env python3
"""
Combined-modality volcano, rendered in the same style as
volcano_plot_{method}_all_categories_new.pdf in the paper.

Data: all three modalities pooled (mouse RNA-seq, human RNA-seq, human DNAm
primary-clock GrimAgeV2).

Visual encodings (matches s01_visualize_results_human.py Figure-3 block):
  - x            = log2 effect size
  - y            = −log10(FDR)
  - color        = log2FoldChange, diverging map of the Figure 3 palette
                   (fig3_palette_2026_11_08.py, FIG3_PALETTE; RdBu_r in July) (via colorbar)
  - shape        = MODALITY  (● mouse RNA, ■ human RNA, ▲ human DNAm)
  - size         = July unscored-point default (constant)
  - black outline layer behind every coloured point
  - non-significant points rendered as faint light-gray background
  - KDE density contour overlay (palette density map, rocket_r in July,
    α = 0.35), clipped at −log10 FDR = 0
  - bracket-connector annotations on the periphery of the plot: exactly the
    interventions of panels b-e (higher outputs on the right, lower on the
    left, with the b-e display text) plus the three comparisons followed up
    in aged mice (marked †); each labelled point is redrawn with a black ring
  - Effect Size colorbar on the right; category and modality keys in the
    empty corners below the FDR line
  - display caps marked where they are drawn: the x ticks end at the cap,
    "≤−30" and "≥30" (as "≥15" on y), and so do the colour-bar ends; the
    dashed line is labelled "FDR = 0.05" at its right end
  - y title "−log10 FDR" set from three pieces so the subscript is 6 pt
    (mathtext would set it at 70% of 7 pt, 4.9 pt)
  - a b-e bar that shares its canonical_target with a followed-up comparison
    is labelled with that comparison's label ("Cold exposure",
    GSE201135), as in panels b-e, so the two cold-exposure labels read together

Print size (2026-11-06): the panel is drawn at its final size, 183 x 86 mm,
and placed at 100% in Figure 3. The July version was drawn at 18 x 10 in and
printed at about 0.44x; marker areas and line widths below are the July values
multiplied by that print scale (PRINT_SCALE, PRINT_SCALE**2 for areas), so the
printed points and lines keep their July appearance while the type is set at
6-7 pt.

Outputs:
  figures/volcano_combined_modalities.{pdf,png}
"""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties
import os
from matplotlib.ticker import FixedLocator, MultipleLocator
from matplotlib.transforms import blended_transform_factory
import importlib.util

# Colour roles of Figure 3 (FIG3_PALETTE; palette B when unset, "july" reproduces the July colours).
_spec = importlib.util.spec_from_file_location(
    "fig3_palette", Path(__file__).resolve().parent / "fig3_palette_2026_11_08.py")
fig3_palette = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fig3_palette)
PAL = fig3_palette.active()
DIVERGING = fig3_palette.diverging(PAL)

MM = 1 / 25.4
FIG_W_MM, FIG_H_MM = 183.0, 86.0
# July figure: 18 x 10 in saved with a tight box of 419.7 mm width, placed at 183 mm.
PRINT_SCALE = 183.0 / 419.7

LABEL_PT = 6.0          # intervention labels
LABEL_PITCH_PT = 6.7    # minimum baseline distance between two neighbouring labels
LABEL_LEADING_PT = 5.9  # baseline distance inside a two-line label (tighter than the pitch)
KEY_PT = 6.5            # category and modality keys
RING_S = 3.6 ** 2       # labelled point, redrawn in its modality shape: 3.6 pt across (cloud markers about 2.8 pt)
RING_LW = 0.6           # black edge of a labelled point
TICK_PT = 7.0
AXIS_PT = 7.0
# July scatter edges: matplotlib default (lines.linewidth 1.5 pt, edge in the face colour), scaled to print
SCATTER_LW = 1.5 * PRINT_SCALE

try:
    from aquarel import load_theme
    _theme = (load_theme("scientific").set_grid(draw=False).set_font(size=AXIS_PT)
              .set_ticks(direction="out", width_major=0.6, width_minor=0.4,
                         size_major=2.5, size_minor=1.5, pad_major=1.5, pad_minor=1.5)
              .set_axis_labels(pad=2)
              .set_transforms(trim="both", offset=4))
    _theme.apply()
except Exception:
    _theme = None

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"
# The panel is saved at its page size without a tight box, so pad_inches does not touch the PDF;
# a small pad keeps any tight-box measurement render (preflight) close to the saved page.
plt.rcParams["savefig.pad_inches"] = 0.02
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["mathtext.it"] = "Arial:italic"
plt.rcParams["mathtext.bf"] = "Arial:bold"
plt.rcParams["font.size"] = AXIS_PT
plt.rcParams["axes.labelsize"] = AXIS_PT
plt.rcParams["axes.titlesize"] = AXIS_PT
plt.rcParams["xtick.labelsize"] = TICK_PT
plt.rcParams["ytick.labelsize"] = TICK_PT
plt.rcParams["legend.fontsize"] = KEY_PT
plt.rcParams["legend.title_fontsize"] = KEY_PT
plt.rcParams["axes.linewidth"] = 0.6


CATEGORY_COLORS = dict(PAL["category"])   # used only for annotation text color
CATEGORY_ORDER = ["Genetic", "Drug", "Disease", "Environment", "Other"]

MODALITY_MARKERS = {
    "mouse RNA-seq": "o",
    "human RNA-seq": "s",
    "human DNAm":    "^",
}
MODALITY_ORDER = list(MODALITY_MARKERS.keys())

SCORE_FLOOR = {"mouse RNA-seq": 100, "human RNA-seq": 100, "human DNAm": 80}

# Canonical aging-field interventions to force-annotate. Each entry carries an
# `expect` direction — "neg" for canonically anti-aging, "pos" for canonically
# pro-aging. A row is force-annotated ONLY if its log2FoldChange matches the
# expected direction; wrong-direction cases (e.g., a study where 5-azacytidine
# shows +log2FC) are left to the normal budget path where the plausibility
# gate + per-case q0 corrections handle them.
FORCE_INCLUDE = [
    # Life-extending / senolytic / clock-reducing — expect anti-aging (neg)
    ("rapamycin", "neg"), ("metformin", "neg"), ("ouabain", "neg"),
    ("fenofibrate", "neg"), ("nicotinamide riboside", "neg"),
    ("resveratrol", "neg"), ("acarbose", "neg"), ("Serpina3n", "neg"),
    ("SGI-110", "neg"), ("Zebularine", "neg"), ("Decitabine", "neg"),
    ("5-aza", "neg"),
    # Yamanaka / reprogramming — expect anti-aging (neg)
    ("4F", "neg"), ("OSKM", "neg"), ("iPSC reprogramming", "neg"),
    ("MPTR", "neg"),
    # Pro-aging canonical — expect +
    ("BMH21", "pos"), ("3M-052", "pos"),
    ("doxorubicin", "pos"), ("bleomycin", "pos"),
    # Accelerated-aging syndromes — expect +
    ("HGPS", "pos"), ("Werner", "pos"), ("Cockayne", "pos"),
]

# Category keywords that disqualify an anti-aging annotation as biologically
# implausible (same set s01_visualize_results_human.py uses)
IMPLAUSIBLE_KEYWORDS = (
    "cancer|carcinoma|metastasis|metastatic|tumor|tumour|glioma|glioblastoma|melanoma|lymphoma|"
    "leukemia|leukaemia|sarcoma|adenoma|neoplas|senescen|passage|Alzheimer|Parkinson|COPD|AMD|"
    "macular degeneration|cirrhosis|fibrosis|atherosclero|colitis|inflammat|infection|viral|"
    "influenza|COVID|sepsis|fatal|mutation|mutant| R[0-9]|oncogen|UPD|imprinting|Paternal UPD|"
    "Maternal UPD|Amsterdam|Lynch|Silver-Russell|SRS|vs \\(Age\\)|vs \\(Survival\\)|\\(Age\\)|"
    "\\(Survival\\)|Proxy|clinical criteria|clinical SRS|syndrome|resistance|disorder"
)

# The three comparisons followed up in aged mice, labelled at their own points and marked by a dagger.
TESTED_KEYS = [("GSE122080", "Ouabain vs Saline in old mice"),
               ("GSE134151", "5xFAD Hippocampus PLX5622 vs Vehicle"),
               ("GSE165974", "Cold exposure vs Room temperature in DPP4_pos APCs")]
TESTED_MARK = " †"
SUB_PT = 6.0            # the subscript of the y title
FDR_LINE_PT = 6.0       # the label of the dashed FDR line


def qualify_with_tested(top, judged):
    """Same rule as create_combined_modality_top_hits_v2.py: the label of a b-e bar that shares its canonical_target
    with a followed-up comparison reads that comparison's label. Idempotent."""
    keys = ["gse_id", "comparison_name", "condition_name"]
    tested = pd.DataFrame(TESTED_KEYS, columns=["gse_id", "comparison_name"]).merge(
        judged[["gse_id", "comparison_name", "canonical_target", "label"]], on=["gse_id", "comparison_name"],
        how="left", validate="one_to_one")
    assert tested["canonical_target"].notna().all() and tested["canonical_target"].is_unique
    target = top.merge(judged[keys + ["canonical_target"]], on=keys, how="left", validate="one_to_one")
    display, changed = top["display"].astype(str).copy(), []
    for t in tested.itertuples():
        hit = (target["canonical_target"].eq(t.canonical_target)
               & ~(target["gse_id"].eq(t.gse_id) & target["comparison_name"].eq(t.comparison_name))).to_numpy()
        for i in top.index[hit]:
            display[i] = t.label
            changed.append(i)
    return display, changed


def _canon(df):
    df = df.copy()
    df["condition_category"] = df["condition_category"].where(
        df["condition_category"].isin(CATEGORY_COLORS), "Other"
    )
    return df


def normalize_size(fs: pd.Series) -> pd.Series:
    """Cubic scaling matching s01's normalize_trust_scores behaviour."""
    v = fs.fillna(50).clip(-30, 140)
    normalized = (v + 30) / 170
    return (normalized ** 3) * 120 + 4


def load_human(modality_suffix, label, clock_filter=None):
    master = pd.read_csv(f"result/main_results_hs_{modality_suffix}_outlier_marked.csv", low_memory=False)
    master = master[master["outlier_status"] == "Normal"].copy()
    if clock_filter is not None:
        master = master[master["clock"] == clock_filter].copy()
    master = _canon(master)
    master = master.drop(columns=[c for c in ["final_score", "intervention"] if c in master.columns])
    scored = pd.read_csv(f"result/main_results_hs_{modality_suffix}_with_scores.csv", low_memory=False)
    key = ["gse_id", "comparison_name", "condition_name"]
    scored = scored[key + ["final_score", "intervention"]].drop_duplicates(subset=key)
    out = master.merge(scored, on=key, how="left")
    out["modality"] = label
    return out


def load_mouse():
    master = pd.read_csv("result/main_results_result_mm_rna_outlier_marked_clean.csv", low_memory=False)
    master = master[master["outlier_status"] == "Normal"].copy()
    master = _canon(master)
    try:
        cm = pd.read_csv("data/case_mapping.csv", low_memory=False)
        key = ["gse_id", "comparison_name", "condition_name"]
        keep = [k for k in key if k in cm.columns] + ["case_id"]
        cm = cm[keep].drop_duplicates(subset=[k for k in key if k in cm.columns])
        master = master.merge(cm, on=[k for k in key if k in cm.columns], how="left")
        scores = pd.read_csv("data/scores_summary.csv")[["dataset", "final_score", "intervention"]]
        master = master.merge(scores, left_on="case_id", right_on="dataset", how="left")
        master = master.drop(columns=["dataset"], errors="ignore")
    except FileNotFoundError:
        master["final_score"] = np.nan
        master["intervention"] = ""
    master["modality"] = "mouse RNA-seq"
    return master


def main():
    print("loading combined Normal atlas...")
    keep = ["modality", "gse_id", "comparison_name", "condition_name",
            "condition_category", "log2FoldChange", "fdr", "pvalue",
            "final_score", "intervention"]
    # Single canonical source: plot ALL Normal (non-outlier) comparisons across the
    # three modalities (104,243 rows; the same outlier_status=='Normal' subset used for
    # the manuscript headline), so the plotted points match the title totals exactly.
    df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)
    df = df[df["outlier_status"] == "Normal"].copy()
    df = _canon(df)
    df = df[[c for c in keep if c in df.columns]]
    df = df.dropna(subset=["log2FoldChange", "fdr"])
    # Clip fdr to the open unit interval so boundary rows (fdr==1 non-significant,
    # fdr==0 maximally significant) are retained rather than dropped.
    df["fdr"] = df["fdr"].clip(lower=1e-300, upper=1.0)
    df["_nlf"] = -np.log10(df["fdr"])
    # Display range of the July panel: Cohen's d within +/-30 and -log10 FDR up to 15. The 0.4% of rows beyond
    # are drawn at the boundary; no row is removed and the header totals count every plotted row.
    XCAP, YCAP = 30.0, 15.0
    df["_d_source"] = df["log2FoldChange"]      # unclipped Cohen's d, used to match labelled rows
    df["log2FoldChange"] = df["log2FoldChange"].clip(-XCAP, XCAP)
    df["_nlf"] = df["_nlf"].clip(upper=YCAP)
    # July point area, scaled to the print size (see module docstring)
    df["_size"] = normalize_size(pd.Series(np.nan, index=df.index)) * PRINT_SCALE ** 2
    for _m in ("mouse RNA-seq", "human RNA-seq", "human DNAm"):
        print(f"  {_m} : {(df['modality'] == _m).sum():,}")
    print(f"combined (Normal only): {len(df):,} rows")

    # ---- figure, drawn at its print size ----
    fig = plt.figure(figsize=(FIG_W_MM * MM, FIG_H_MM * MM))
    AX_L, AX_R, AX_B, AX_T = 10.2, 14.5, 8.5, 5.0   # margins in mm
    ax = fig.add_axes([AX_L / FIG_W_MM, AX_B / FIG_H_MM,
                       1 - (AX_L + AX_R) / FIG_W_MM, 1 - (AX_B + AX_T) / FIG_H_MM])

    # Layer 0 — faint density contour from ALL points. -log10 FDR cannot be negative, so the
    # shading is clipped at 0 (a clip path; the contour levels themselves are unchanged).
    try:
        kde_df = df[["log2FoldChange", "_nlf"]].rename(columns={"_nlf": "neg_log10_fdr"})
        before = set(ax.collections)
        sns.kdeplot(
            data=kde_df, x="log2FoldChange", y="neg_log10_fdr",
            fill=True, cmap=fig3_palette.density(PAL), alpha=PAL["density_alpha"], levels=6,
            linewidths=1.0 * PRINT_SCALE, ax=ax, zorder=4,
        )
        kde_layers = [c for c in ax.collections if c not in before]
        assert kde_layers
    except Exception as e:
        kde_layers = []
        print(f"  kde skipped: {e}")

    # Layer 1 — non-significant as light-gray background
    ns = df[df["fdr"] >= 0.05]
    # Every point lies inside the axes (display caps), so the point layers need no clipping.
    ax.scatter(ns["log2FoldChange"], ns["_nlf"],
               s=ns["_size"], c=PAL["nonsignificant"], alpha=0.3, zorder=1, rasterized=True,
               linewidths=SCATTER_LW, clip_on=False)

    # Layer 2 & 3 — significant points, one pass per modality
    sig = df[df["fdr"] < 0.05].copy()
    vmax = np.abs(sig["log2FoldChange"]).max() if not sig.empty else 1
    vmin = -vmax
    scatter_handle = None
    for modality in MODALITY_ORDER:
        m = sig[sig["modality"] == modality]
        if m.empty:
            continue
        marker = MODALITY_MARKERS[modality]
        # black outline
        ax.scatter(m["log2FoldChange"], m["_nlf"],
                   s=m["_size"] * 1.5, c=PAL["point_outline"], alpha=1.0,
                   marker=marker, zorder=5, rasterized=True, linewidths=SCATTER_LW, clip_on=False)
        # colored by effect size
        scatter_handle = ax.scatter(
            m["log2FoldChange"], m["_nlf"],
            s=m["_size"], c=m["log2FoldChange"],
            cmap=DIVERGING, vmin=vmin, vmax=vmax,
            alpha=0.9, marker=marker, zorder=6, rasterized=True, linewidths=SCATTER_LW,
            clip_on=False,
        )

    # ---- axis ranges (July rule: data range widened by 22% on each side for the label columns) ----
    xlim = ax.get_xlim()
    # -log10 FDR starts at 0 and the density shading is clipped there, so the y range starts just
    # below 0 (room for the keys under the FDR line) instead of at the autoscaled -0.68.
    ylim = (-0.35, ax.get_ylim()[1])
    x_range = xlim[1] - xlim[0]
    new_xlim = (xlim[0] - x_range * 0.22, xlim[1] + x_range * 0.22)
    ax.set_xlim(new_xlim)
    ax.set_ylim(ylim)
    # density clip: the full axes width, from -log10 FDR = 0 to the top of the axes (kept inside the axes,
    # so the clipped layer's extent stays inside the panel)
    kde_clip = Rectangle((0.0, 0.0), 1.0, ylim[1], transform=blended_transform_factory(ax.transAxes, ax.transData))
    for c in kde_layers:
        c.set_clip_path(kde_clip)
    # Label columns run from just above the FDR line to just below the row of points capped at
    # -log10 FDR = 15 (the top label's connector passes about 4 pt under that row); the keys sit
    # below the FDR line.
    # The lowest label's centre sits 0.34 y-units (about 4.3 pt) above the FDR line, so its
    # descenders (about 2.8-3.4 pt below the centre of a 6 pt label) stay clear of the dashes.
    y_bottom = -np.log10(0.05) + 0.34
    y_top    = 14.6
    ax_h_pt = (1 - (AX_B + AX_T) / FIG_H_MM) * FIG_H_MM / 25.4 * 72
    pt_per_y = ax_h_pt / (ylim[1] - ylim[0])
    capacity = int(np.floor((y_top - y_bottom) * pt_per_y / LABEL_PITCH_PT)) + 1
    print(f"  label column capacity at {LABEL_PT:g} pt: {capacity} lines per side")

    # ---- labelled set: exactly the interventions of panels b-e, plus the three followed-up comparisons ----
    # Panels b-e show, per category, the five interventions with higher and the five with lower outputs
    # (result/top_hits_combined_modalities_v2.csv). Each is labelled with its b-e display text at the plotted
    # point of that comparison (gse_id and comparison_name; where a series holds several rows of the same
    # comparison, the row with the same Cohen's d and FDR). The three comparisons followed up in aged mice are
    # labelled at their own points and marked by a dagger.
    top = pd.read_csv("result/top_hits_combined_modalities_v2.csv", low_memory=False)
    assert len(top) == 40 and top.groupby(["category", "direction"]).size().eq(5).all()
    judged_path = os.path.join(os.environ.get("CLOCKBASE_JUDGEMENT") or os.path.join(os.environ["CLOCKBASE_VERIFIER_ROOT"], "followup-judgment-001"), "ranked_with_judgement.csv")
    if not os.path.exists(judged_path):
        raise SystemExit(f"Cannot proceed - required data file not found: {judged_path}")
    judged = pd.read_csv(judged_path, low_memory=False)
    top["display"], changed = qualify_with_tested(top, judged)
    assert top.loc[changed, ["gse_id", "display"]].values.tolist() == [["GSE201135", "Cold exposure"]]
    assert np.array_equal(np.sign(top["log2FoldChange"]), top["direction"])
    lab_rows = []

    def _point(gse, comp, d=None, fdr=None):
        m = df[(df["gse_id"] == gse) & (df["comparison_name"] == comp)]
        if len(m) > 1 and d is not None:
            m = m[np.isclose(m["_d_source"], d) & np.isclose(m["fdr"], fdr)]
        assert len(m) == 1, (gse, comp, len(m))
        return m.iloc[0]

    for _, r in top.iterrows():
        p = _point(r["gse_id"], r["comparison_name"], r["log2FoldChange"], r["fdr"])
        assert p["modality"] == r["modality"] and p["condition_name"] == r["condition_name"]
        lab_rows.append(dict(modality=r["modality"], condition_category=r["category"], gse_id=r["gse_id"],
                             comparison_name=r["comparison_name"], condition_name=r["condition_name"],
                             intervention=r["intervention"], log2FoldChange=r["log2FoldChange"], fdr=r["fdr"],
                             final_score=r["final_score"], label=str(r["display"]).replace("Barretts", "Barrett's esophagus").replace("10 Gy x3", "10 Gy ×3"),
                             x=p["log2FoldChange"], y=p["_nlf"]))
    for gse, comp in TESTED_KEYS:
        if ((top["gse_id"] == gse) & (top["comparison_name"] == comp)).any():
            continue    # already labelled as a b-e intervention
        p = _point(gse, comp)
        lab_rows.append(dict(modality=p["modality"], condition_category=p["condition_category"], gse_id=gse,
                             comparison_name=comp, condition_name=p["condition_name"],
                             intervention=p["intervention"], log2FoldChange=p["_d_source"], fdr=p["fdr"],
                             final_score=p["final_score"], label=str(p["intervention"]) + TESTED_MARK,
                             x=p["log2FoldChange"], y=p["_nlf"]))
    lab = pd.DataFrame(lab_rows)
    lab["_side"] = np.where(lab["x"] > 0, "pos", "neg")
    expected = ({(g, c) for g, c in zip(top["gse_id"], top["comparison_name"])} | set(TESTED_KEYS))
    assert set(zip(lab["gse_id"], lab["comparison_name"])) == expected and len(lab) == len(expected)
    assert lab["fdr"].lt(0.05).all()            # every labelled point is in the coloured layer

    # Label geometry at print size: text widths are measured with the renderer, so a label is set
    # on two lines only when its one-line form would sit over plotted points.
    renderer = fig.canvas.get_renderer()
    label_prop = FontProperties(family="Arial", weight="bold", size=LABEL_PT)
    ax_w_pt = (1 - (AX_L + AX_R) / FIG_W_MM) * FIG_W_MM / 25.4 * 72
    pt_per_x = ax_w_pt / (new_xlim[1] - new_xlim[0])
    line_y = LABEL_PITCH_PT / pt_per_y
    lead_y = LABEL_LEADING_PT / pt_per_y
    sig_x = sig["log2FoldChange"].to_numpy()
    sig_y = sig["_nlf"].to_numpy()

    def _width_pt(s):
        w, _, _ = renderer.get_text_width_height_descent(s, label_prop, ismath=False)
        return w * 72 / fig.dpi

    def _split(label):
        words = label.split(" ")
        if words[-1] == TESTED_MARK.strip():   # the dagger stays with the word before it
            words = words[:-2] + [" ".join(words[-2:])]
        if len(words) < 2:
            return [label]
        best = min(range(1, len(words)),
                   key=lambda i: max(_width_pt(" ".join(words[:i])), _width_pt(" ".join(words[i:]))))
        return [" ".join(words[:best]), " ".join(words[best:])]

    def _x_text(side):
        # columns start at the axes edge (4 pt inside the offset y spine), which leaves the long
        # labels room before the column of points capped at |d| = 30
        return new_xlim[1] if side == "pos" else new_xlim[0]

    def _layout(rows, wrap):
        # Rows in column order (see _order); labels stacked from y_bottom to y_top. A two-line label is a
        # block whose lines are LABEL_LEADING_PT apart; the gap between neighbouring labels is the
        # same everywhere in the column and at least LABEL_PITCH_PT, so a block reads as one label.
        parts = [_split(t) if w else [t] for t, w in zip(rows["label"], wrap)]
        inner = [(len(p) - 1) * lead_y for p in parts]
        n = len(parts)
        gap = (y_top - y_bottom - sum(inner)) / (n - 1) if n > 1 else 0.0
        fits = n <= 1 or gap >= line_y - 1e-9
        y0, out = y_bottom, []
        for p, h in zip(parts, inner):
            out.append((p, y0 + h / 2, h))     # (lines, block centre, block height between baselines)
            y0 += h + gap
        return out, fits, gap

    def _covered(p, y_pos, h, side):
        w = max(_width_pt(s) for s in p) / pt_per_x
        x0 = _x_text(side)
        xa, xb = (x0 - w, x0) if side == "pos" else (x0, x0 + w)
        # a point counts as under the label when its marker (1.4 pt half-width) comes within 1.6 pt of
        # the text horizontally, or within half a marker of the label row vertically
        mx, my = 3.0 / pt_per_x, 1.2 / pt_per_y
        half = h / 2 + line_y / 2
        return int(((sig_x > xa - mx) & (sig_x < xb + mx) & (sig_y > y_pos - half - my) & (sig_y < y_pos + half + my)).sum())

    def _leader(x, y, y_pos, side):
        # July connector: a short level stub from the point, a steep run to the label row, then level
        # to the label at the column edge
        s = 1 if side == "pos" else -1
        xh1 = x + s * x_range * 0.01
        xsl = xh1 + s * x_range * 0.025
        return [(x, y), (xh1, y), (xsl, y_pos), (_x_text(side) - s * x_range * 0.005, y_pos)]

    def _n_cross(rows, lay, side):
        # number of crossings between the drawn connectors of one column (segment intersection test)
        P = np.array([_leader(x, y, y_pos, side) for x, y, (_, y_pos, _) in zip(rows["x"], rows["y"], lay)])
        a, b = P[:, :-1].reshape(-1, 2), P[:, 1:].reshape(-1, 2)
        d = b - a
        # o[i, j, e]: side of segment i's line on which end e of segment j lies
        o = np.sign(d[:, None, 0, None] * (np.stack([a, b], 1)[None, :, :, 1] - a[:, None, None, 1])
                    - d[:, None, 1, None] * (np.stack([a, b], 1)[None, :, :, 0] - a[:, None, None, 0]))
        straddle = o[..., 0] * o[..., 1] < 0
        owner = np.repeat(np.arange(len(P)), 3)
        return int((straddle & straddle.T & (owner[:, None] != owner[None, :])).sum() // 2)

    sig_pt = np.c_[sig_x * pt_per_x, sig_y * pt_per_y]
    over_cache = {}

    def _over(x, y, y_pos, side):
        # coloured points that a connector runs over (centre within 1.6 pt of the line), not counting
        # those under the white rim of its own point (2.6 pt)
        key = (x, y, y_pos, side)
        if key not in over_cache:
            L = np.array(_leader(x, y, y_pos, side)) * [pt_per_x, pt_per_y]
            a, d = L[:-1], np.diff(L, axis=0)
            t = np.clip(((sig_pt[:, None] - a) * d).sum(-1) / (d * d).sum(-1), 0, 1)
            dist = np.hypot(*(sig_pt[:, None] - a - t[..., None] * d).transpose(2, 0, 1)).min(1)
            over_cache[key] = int(((dist < 1.6) & (np.hypot(*(sig_pt - L[0]).T) > 2.6)).sum())
        return over_cache[key]

    def _untangle(rows, wrap, side):
        # move one label to another row, or swap two labels, whenever that removes a crossing, or keeps
        # the crossings and runs the connectors over fewer points, or over as many with shorter drops
        def cost(r, w):
            lay, _, _ = _layout(r, w)
            rl = list(zip(r["x"], r["y"], [y_pos for _, y_pos, _ in lay]))
            return (_n_cross(r, lay, side), sum(_over(x, y, y_pos, side) for x, y, y_pos in rl),
                    sum(abs(y_pos - y) for _, y, y_pos in rl))
        n = len(rows)
        moves = [[*range(i)] + [j] + [*range(i + 1, j)] + [i] + [*range(j + 1, n)]
                 for i in range(n) for j in range(i + 1, n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    rest = [k for k in range(n) if k != i]
                    moves.append(rest[:j] + [i] + rest[j:])
        best, improved = cost(rows, wrap), True
        while improved:
            improved = False
            for perm in moves:
                r2, w2 = rows.iloc[perm].reset_index(drop=True), [wrap[k] for k in perm]
                c = cost(r2, w2)
                if c < best:
                    rows, wrap, best, improved = r2, w2, c, True
        return rows, wrap

    def _order(rows, wrap):
        # A connector runs level from its label to just beyond its point and then drops to the point,
        # so it crosses another connector only when its row lies within the drop of a connector that
        # reaches farther out. The point farthest out therefore takes the free row nearest its own
        # height, then the next farthest, and so on; points at the same distance keep the order of
        # their heights. The wrap flag of each label moves with it. _untangle then removes the
        # crossings left between neighbouring points, whose stubs and drops overlap.
        for _ in range(20):
            lay, _, _ = _layout(rows, wrap)
            slots = np.array([y_pos for _, y_pos, _ in lay])
            xs, ys = np.abs(rows["x"].to_numpy()), rows["y"].to_numpy()
            free, take = list(range(len(rows))), {}
            for i in sorted(range(len(rows)), key=lambda i: (-xs[i], ys[i])):
                k = min(free, key=lambda k: (abs(slots[k] - ys[i]), k))
                take[k] = i
                free.remove(k)
            for tie in pd.Series(range(len(rows))).groupby(xs).groups.values():
                ks = sorted(k for k, i in take.items() if i in set(tie))
                for k, i in zip(ks, sorted(tie, key=lambda i: ys[i])):
                    take[k] = i
            perm = [take[k] for k in range(len(rows))]
            if perm == list(range(len(rows))):
                break
            rows = rows.iloc[perm].reset_index(drop=True)
            wrap = [wrap[i] for i in perm]
        return rows, wrap

    def _fit_side(rows, side):
        rows = rows.sort_values(["y", "x"], kind="stable").reset_index(drop=True)
        wrap = [False] * len(rows)
        for _ in range(10):
            rows, wrap = _untangle(*_order(rows, wrap), side)
            lay, fits, gap = _layout(rows, wrap)
            changed = False
            for i, (t, (p, y_pos, h)) in enumerate(zip(rows["label"], lay)):
                if wrap[i] or len(_split(t)) < 2:
                    continue
                p2 = _split(t)
                if _covered(p, y_pos, h, side) > _covered(p2, y_pos, (len(p2) - 1) * lead_y, side):
                    wrap[i] = True
                    changed = True
            if not changed:
                break
        lay, fits, gap = _layout(rows, wrap)
        assert fits, f"{side} label column overflows: {len(rows)} labels, pitch {gap * pt_per_y:.2f} pt"
        return rows, lay, gap

    layouts = {}
    for side in ("pos", "neg"):
        rows, lay, gap = _fit_side(lab[lab["_side"] == side], side)
        layouts[side] = (rows, lay)
        n_lines = sum(len(p) for p, _, _ in lay)
        print(f"  {side}: {len(rows)} labels in {n_lines} lines, label pitch {gap * pt_per_y:.2f} pt, "
              f"{_n_cross(rows, lay, side)} connector crossings, connectors over "
              f"{sum(_over(x, y, y_pos, side) for x, y, (_, y_pos, _) in zip(rows['x'], rows['y'], lay))} points")
        for (_, r), (p, y_pos, h) in zip(rows.iterrows(), lay):
            c = _covered(p, y_pos, h, side)
            print(f"     {' / '.join(p)!r} [{r['condition_category']}, {r['gse_id']}]"
                  + (f": {len(p)} lines" if len(p) > 1 else "")
                  + (f", {c} points under the label" if c else ""))

    # Persist the labelled set for audit: the existing columns plus comparison_name and the printed label.
    # log2FoldChange and fdr are the source values (Cohen's d before the +/-30 display cap).
    annotations_cols = ["_side", "modality", "condition_category", "gse_id", "condition_name", "intervention",
                        "log2FoldChange", "fdr", "final_score", "comparison_name", "label"]
    out_ann = pd.concat([layouts["pos"][0], layouts["neg"][0]], ignore_index=True)[annotations_cols]
    assert out_ann["intervention"].astype(str).str.casefold().eq("ouabain").sum() == 1
    out_ann.to_csv("result/combined_volcano_annotations.csv", index=False)
    print(f"  audit list: result/combined_volcano_annotations.csv")
    print(f"annotations: {len(layouts['pos'][0])} higher + {len(layouts['neg'][0])} lower")
    for side in ("pos", "neg"):
        print(f"    {side}: {layouts[side][0]['condition_category'].value_counts().to_dict()}")

    # Labels are stacked from y_bottom upward and spread evenly, in the column order set by _order and
    # _untangle;
    # each connector runs from the label to its point, and the point is redrawn above the cloud and
    # the connectors with a black ring, so the connector ends at the ring.
    label_groups, drawn = [], []
    for side in ("pos", "neg"):
        rows, lay = layouts[side]
        x_text = _x_text(side)
        for (_, r), (p, y_pos, h) in zip(rows.iterrows(), lay):
            color = CATEGORY_COLORS.get(r["condition_category"], "#333")
            ha = "right" if side == "pos" else "left"
            ax.plot(*zip(*_leader(r["x"], r["y"], y_pos, side)),
                    color="black", ls=":", alpha=0.4, lw=1.0 * PRINT_SCALE, zorder=7)
            # A two-line label is set as a block centred on the connector: lines LABEL_LEADING_PT
            # apart (tighter than the gap to the neighbouring labels), both flush with the column edge.
            z = 11 if r["label"].endswith(TESTED_MARK) else 10
            group = []
            for i, txt in enumerate(p):
                t = ax.text(x_text, y_pos + h / 2 - i * lead_y, txt, ha=ha, va="center",
                            fontsize=LABEL_PT, color=color, weight="bold", zorder=z)
                group.append(t)
            label_groups.append(group)
            drawn.append((r["gse_id"], r["comparison_name"], " ".join(t.get_text() for t in group)))
    assert {(g, c) for g, c, _ in drawn} == expected and len(drawn) == len(expected)
    assert dict(((g, c), s) for g, c, s in drawn) == dict(zip(zip(lab["gse_id"], lab["comparison_name"]), lab["label"]))
    # Each labelled point is redrawn in its effect-size colour and modality shape, 3.6 pt across with a
    # 0.6 pt black edge, over a 0.5 pt white rim that separates the ring from the dark outlines of the cloud.
    for modality in MODALITY_ORDER:
        m = lab[lab["modality"] == modality]
        if not m.empty:
            rim = ax.scatter(m["x"], m["y"], s=RING_S, c="white", marker=MODALITY_MARKERS[modality],
                             edgecolors="white", linewidths=RING_LW + 2 * 0.5, zorder=7.9, clip_on=False)
            rim.set_joinstyle("round")
            ax.scatter(m["x"], m["y"], s=RING_S, c=m["x"], cmap=DIVERGING, vmin=vmin, vmax=vmax,
                       marker=MODALITY_MARKERS[modality], edgecolors="black", linewidths=RING_LW,
                       zorder=8, clip_on=False)

    # ---- labels & title ----
    ax.axhline(-np.log10(0.05), ls="--", color=PAL["fdr_line"], lw=0.5, zorder=0)
    ax.set_xlabel("Cohen's d")
    # Cohen's d is drawn within +/-XCAP: the ticks end at the cap, and the end ticks read "≤−30" and "≥30" where rows
    # lie beyond it (as "≥15" on y). The label columns keep the July axis range.
    beyond_lo, beyond_hi = bool((df["_d_source"] < -XCAP).any()), bool((df["_d_source"] > XCAP).any())
    assert beyond_lo and beyond_hi
    x_ticks = np.linspace(-XCAP, XCAP, 5)
    ax.xaxis.set_major_locator(FixedLocator(x_ticks))
    ax.xaxis.set_minor_locator(FixedLocator(np.arange(-XCAP, XCAP + 1e-9, 5.0)))
    cap_text = lambda v: (("≤" if v < 0 and beyond_lo else "≥" if v > 0 and beyond_hi else "") + f"{v:g}").replace("-", "−")
    ax.set_xticklabels([cap_text(v) if abs(v) == XCAP else f"{v:g}".replace("-", "−") for v in x_ticks])
    assert new_xlim[0] < -XCAP and new_xlim[1] > XCAP
    # -log10 FDR is drawn up to 15 (display cap); the tick at the cap reads "≥15" so the capped row reads as capped.
    ax.set_yticks([0, 5, 10, 15])
    ax.set_yticklabels(["0", "5", "10", "≥15"])
    ax.yaxis.set_minor_locator(MultipleLocator(1))
    ax.set_ylim(ylim)
    # Manuscript-canonical totals: outlier_status='Normal' rows.
    # Mouse: additionally deduplicated by (gse_id, condition_name, case_id).
    # Human RNA-seq and DNAm: Normal rows directly.
    # Title totals derive from the SAME plotted Normal atlas so figure and caption
    # cannot diverge (outlier_status=='Normal' across the three modalities).
    by_mod = df["modality"].value_counts()
    n_mouse = int(by_mod.get("mouse RNA-seq", 0))
    n_hr = int(by_mod.get("human RNA-seq", 0))
    n_hd = int(by_mod.get("human DNAm", 0))
    total = n_mouse + n_hr + n_hd
    sig_count = int((df["fdr"] < 0.05).sum())
    fig.text(0.5, 1 - 0.7 / FIG_H_MM, "ClockBase Agent, combined aging-intervention atlas",
             ha="center", va="top", fontsize=7.0, weight="bold")
    fig.text(0.5, 1 - 3.6 / FIG_H_MM,
             f"Rows supporting inference: {total:,} "
             f"(mouse RNA-seq: {n_mouse:,} · "
             f"human RNA-seq: {n_hr:,} · "
             f"human DNAm: {n_hd:,})   |   "
             f"Significant rows (FDR < 0.05): {sig_count:,} ({100*sig_count/total:.1f}%)",
             ha="center", va="top", fontsize=6.5)

    # ---- colorbar (Effect Size), to the right of the axes, spanning the label columns ----
    if scatter_handle is not None:
        ax_pos = ax.get_position()
        fy = lambda yv: ax_pos.y0 + (yv - ylim[0]) / (ylim[1] - ylim[0]) * ax_pos.height
        cb_y0, cb_y1 = fy(y_bottom), fy(y_top)
        cax = fig.add_axes([ax_pos.x1 + 1.6 / FIG_W_MM, cb_y0, 2.2 / FIG_W_MM, cb_y1 - cb_y0])
        cbar = fig.colorbar(scatter_handle, cax=cax)
        if fig3_palette.name() != "july":
            cbar.solids.set_alpha(1.0)  # the scale keys the stacked points, which read as opaque
        cbar.outline.set_linewidth(0.5)
        cbar.ax.tick_params(labelsize=6.5)
        # the colour scale runs to the display cap, so its ends read as capped too
        assert np.isclose(vmax, XCAP), vmax
        cb_ticks = np.arange(-XCAP, XCAP + 1e-9, 10.0)
        cbar.set_ticks(cb_ticks)
        cbar.set_ticklabels([cap_text(v) if abs(v) == XCAP else f"{v:g}".replace("-", "−") for v in cb_ticks])
        cbar.set_label("Cohen's d", rotation=270, labelpad=7, fontsize=7.0)

    # ---- keys: category (color-code for annotations) + modality (shape), below the FDR line ----
    # The key lists only categories that colour a drawn label.
    shown_categories = set(lab["condition_category"])
    cat_handles = [
        Line2D([0], [0], marker="s", color="w",
               markerfacecolor=CATEGORY_COLORS[c], markeredgecolor=PAL["category_key_edge"],
               markeredgewidth=0.4, markersize=4.6, linestyle="None", label=c)
        for c in CATEGORY_ORDER if c in shown_categories
    ]
    mod_handles = [
        Line2D([0], [0], marker=MODALITY_MARKERS[m], color="black",
               markerfacecolor=PAL["modality_key_face"], markeredgewidth=0.5, markersize=4.2,
               linestyle="None", label=m)
        for m in MODALITY_ORDER
    ]
    key_kw = dict(frameon=False, fontsize=KEY_PT, title_fontsize=KEY_PT, borderaxespad=0.3,
                  borderpad=0.0, handlelength=0.9, handletextpad=0.35, columnspacing=1.0,
                  labelspacing=0.25)
    # The anchor box is the axes widened by the legend's border pad on each side, so the key's outer
    # edges line up with the label columns, which start at the axes edges.
    pad_ax = key_kw["borderaxespad"] * KEY_PT / ax_w_pt
    key_kw["bbox_to_anchor"] = (-pad_ax, 0.0, 1.0 + 2 * pad_ax, 1.0)
    cat_legend = ax.legend(
        handles=cat_handles, title="Intervention category (annotation color)",
        loc="lower left", ncol=len(cat_handles), **key_kw,
    )
    cat_legend._legend_box.align = "left"
    ax.add_artist(cat_legend)
    mod_legend = ax.legend(
        handles=mod_handles, title="Modality (point shape)",
        loc="lower right", ncol=len(mod_handles), **key_kw,
    )
    mod_legend._legend_box.align = "left"

    if _theme is not None:
        try:
            _theme.apply_transforms()
        except Exception:
            pass

    # The y title is set 2 pt left of the tick labels beside it (0, 5, 10). The wider "≥15" sits at the
    # top of the axis, above the title's extent, so it does not push the title off the panel edge.
    # It is set from three pieces, "−log" and " FDR" at 7 pt and the subscript "10" at 6 pt, 0.2 em below the baseline
    # (mathtext would set the subscript at 4.9 pt). The pieces run upwards, centred on the axes.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ticks = [t for t in ax.get_yticklabels() if t.get_visible() and t.get_text()]
    side_ticks = [t.get_window_extent(renderer) for t in ticks if t.get_text() != "≥15"]
    cap_tick = [t.get_window_extent(renderer) for t in ticks if t.get_text() == "≥15"]
    x_right = min(b.x0 for b in side_ticks) - 2.0 * fig.dpi / 72
    px_pt = fig.dpi / 72
    pieces = [("−log", AXIS_PT, 0.0), ("10", SUB_PT, -0.2 * AXIS_PT), (" FDR", AXIS_PT, 0.0)]
    widths = [renderer.get_text_width_height_descent(t, FontProperties(family="Arial", size=sz), ismath=False)[0]
              for t, sz, _ in pieces]
    ax_box = ax.get_window_extent(renderer)
    y_start = (ax_box.y0 + ax_box.y1) / 2 - sum(widths) / 2
    title_parts, y_run = [], y_start
    for (t, sz, shift), w in zip(pieces, widths):
        # rotated by 90 degrees: the baseline is vertical and "below the baseline" points right, towards the axis
        title_parts.append(fig.text(0, y_run / fig.bbox.height, t, fontsize=sz, rotation=90, rotation_mode="anchor",
                                    ha="left", va="baseline", color=ax.yaxis.label.get_color()))
        title_parts[-1]._shift_px = -shift * px_pt
        y_run += w
    ax.set_ylabel("")
    fig.canvas.draw()
    # the right edge of the main pieces (their descent side) at x_right, as for the July label
    main_right = max(p.get_window_extent(renderer).x1 for p in title_parts if p.get_fontsize() == AXIS_PT)
    base_x = x_right - (main_right - 0)  # main pieces are placed at x = 0: their right edge is their descent
    for p in title_parts:
        p.set_x((base_x + p._shift_px) / fig.bbox.width)
    fig.canvas.draw()
    boxes = [p.get_window_extent(renderer) for p in title_parts]
    yl_x0, yl_x1 = min(b.x0 for b in boxes), max(b.x1 for b in boxes)
    yl_y0, yl_y1 = min(b.y0 for b in boxes), max(b.y1 for b in boxes)
    assert yl_x0 >= 0.3 / 25.4 * fig.dpi, "y title runs off the panel"
    for tb in side_ticks + cap_tick:   # no piece of the title touches a tick label
        assert all(b.x1 < tb.x0 - 0.5 * px_pt or b.y1 < tb.y0 or b.y0 > tb.y1 for b in boxes), "y title meets a tick"
    assert all(yl_y1 < b.y0 - 1.0 * px_pt for b in cap_tick), "y title meets the ≥15 tick label"
    assert all(b0.y1 <= b1.y0 + 0.5 for b0, b1 in zip(boxes, boxes[1:])), "y title pieces overlap"
    assert [p.get_fontsize() for p in title_parts] == [AXIS_PT, SUB_PT, AXIS_PT]

    # The dashed line at FDR = 0.05 is labelled at its right end, just below it, right-aligned with the axes; the
    # modality key title on the same side ends well left of it.
    fdr_y = -np.log10(0.05)
    fdr_label = ax.text(new_xlim[1], fdr_y, "FDR = 0.05", ha="right", va="top", fontsize=FDR_LINE_PT, color=PAL["fdr_line"],
                        zorder=10, transform=ax.transData)
    fdr_label.set_y(fdr_y - 0.45 * px_pt / (ax_box.height / (ylim[1] - ylim[0])))
    fig.canvas.draw()
    fb = fdr_label.get_window_extent(renderer)
    line_px = ax.transData.transform((0, fdr_y))[1]
    assert fb.y1 <= line_px - 0.3 * px_pt, "FDR label touches the dashed line"
    for leg in (cat_legend, mod_legend):
        for t in leg.get_texts() + [leg.get_title()]:
            tb = t.get_window_extent(renderer)
            assert fb.x0 > tb.x1 + 2.0 * px_pt or fb.y0 > tb.y1 + 0.5 * px_pt or fb.y1 < tb.y0 - 0.5 * px_pt, \
                f"FDR label meets the key text {t.get_text()!r}"

    # One white backing per label (both lines of a two-line label together), drawn under all
    # label text and over the points and connectors, so neighbouring labels never cover each
    # other's glyphs and a connector ends at the edge of its label.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    pad_px = 0.6 * fig.dpi / 72
    for group in label_groups:
        bbs = [t.get_window_extent(renderer) for t in group]
        (x0, y0), (x1, y1) = inv.transform([(min(b.x0 for b in bbs) - pad_px, min(b.y0 for b in bbs)),
                                            (max(b.x1 for b in bbs) + pad_px, max(b.y1 for b in bbs))])
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor="white", edgecolor="none",
                               zorder=9.5))

    Path("figures").mkdir(exist_ok=True)
    pdf = "figures/volcano_combined_modalities.pdf"
    png = "figures/volcano_combined_modalities.png"
    # Saved at the exact print size (no tight box); rasterized point layers at 600 dpi.
    fig.savefig(pdf, dpi=600)
    fig.savefig(png, dpi=600)
    print(f"\n✓ {pdf}")
    print(f"✓ {png}")


if __name__ == "__main__":
    main()
