#!/usr/bin/env python3
"""
Sample-level box-plot + experimental-design schematic figures for selected
ClockBase Agent top examples.

For each example: LEFT = simple vector schematic of the experimental design,
RIGHT = box plot of real per-sample DNAm biological age (case vs control) with
a Welch t-test annotation.

DATA INTEGRITY: every value is a REAL per-sample clock prediction loaded from the
ClockBase human-DNAm store. Groups are reconstructed from GEO sample titles /
harmonized metadata, matching the published comparison. No simulated data.

Per-sample ages : ClockBase/data/human_dnam/sample_age/<GSE>.csv  (GSM-keyed)
Sample metadata : ClockBase/data/human_dnam/processed_metadata_all.csv
Created 2026-06-13.
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy import stats
warnings.filterwarnings("ignore")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nature_figure_utils import (
    retro_style, trim_axes, figsize_nature, save_for_review,
    RUST, GRAY, BLUE, GOLD, SAGE, INK, ALPHA,
)

ROOT = "/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base"
CB = f"{ROOT}/ClockBase/data/human_dnam"
OUT = f"{ROOT}/figure_rev1/sample_level_examples_2026-06-13"

# ----------------------------------------------------------------------------
# Data access
# ----------------------------------------------------------------------------
_META = None
def meta_for(gse):
    """Return harmonized sample metadata rows for one GSE."""
    global _META
    if _META is None:
        _META = pd.read_csv(f"{CB}/processed_metadata_all.csv", low_memory=False)
    return _META[_META["GEO"] == gse].copy()

def ages_for(gse):
    """Return per-sample clock ages for one GSE (GSM-keyed). Real values only."""
    df = pd.read_csv(f"{CB}/sample_age/{gse}.csv")
    df = df.rename(columns={"sample_id": "GSM"})
    return df

def build_groups(spec):
    """Merge metadata + ages, assign case/control via spec['group_fn'].
    Returns DataFrame with columns GSM, group ('case'/'control'), <clock>."""
    m = meta_for(spec["gse"])
    a = ages_for(spec["gse"])
    df = m.merge(a, on="GSM", how="inner")
    df["group"] = spec["group_fn"](df)
    df = df[df["group"].isin(["case", "control"])].copy()
    clock = spec["clock"]
    # keep only REAL (non-imputed placeholder) clock values: drop if column is a
    # single repeated constant (uniq==1 => global-mean placeholder)
    if df[clock].round(4).nunique() <= 1:
        raise ValueError(f"{spec['gse']}: clock {clock} is an imputed placeholder")
    df = df.dropna(subset=[clock])
    return df, clock

# ----------------------------------------------------------------------------
# Left panel: experimental-design schematic (pure vector)
# ----------------------------------------------------------------------------
def _box(ax, x, y, w, h, text, fc, ec=INK, fs=7.5, tc=INK, lw=0.9, weight="normal"):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h,
                 boxstyle="round,pad=0.012,rounding_size=0.04",
                 facecolor=fc, edgecolor=ec, linewidth=lw, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=tc,
            zorder=4, weight=weight, linespacing=1.25)

def _arrow(ax, x0, y0, x1, y1, lw=1.0, color=INK):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1),
                 arrowstyle="-|>,head_width=2.6,head_length=4.2",
                 mutation_scale=1, lw=lw, color=color,
                 shrinkA=1, shrinkB=1, zorder=2))

def draw_schematic(ax, spec, n_case, n_control):
    """Top-to-bottom design flow: source -> two arms -> assay -> clock -> age."""
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    cc, ctrlc = spec["case_color"], spec["ctrl_color"]
    # source / cohort
    _box(ax, 0.5, 0.92, 0.74, 0.115, spec["source_label"], fc="white", fs=7.5,
         weight="bold")
    # two arms
    yarm = 0.66
    _box(ax, 0.27, yarm, 0.42, 0.135,
         f"{spec['case_label']}\n(n = {n_case})", fc=cc, fs=7.2)
    _box(ax, 0.73, yarm, 0.42, 0.135,
         f"{spec['ctrl_label']}\n(n = {n_control})", fc=ctrlc, fs=7.2)
    _arrow(ax, 0.42, 0.862, 0.30, yarm + 0.072)
    _arrow(ax, 0.58, 0.862, 0.70, yarm + 0.072)
    # assay
    yas = 0.40
    _box(ax, 0.5, yas, 0.78, 0.105, spec["assay_label"], fc="white", fs=7.2)
    _arrow(ax, 0.27, yarm - 0.072, 0.42, yas + 0.057)
    _arrow(ax, 0.73, yarm - 0.072, 0.58, yas + 0.057)
    # clock
    ycl = 0.155
    _box(ax, 0.5, ycl, 0.78, 0.11,
         f"{spec['clock_label']}\n→ biological age", fc="white", fs=7.2)
    _arrow(ax, 0.5, yas - 0.053, 0.5, ycl + 0.057)

# ----------------------------------------------------------------------------
# Right panel: box plot of biological age, case vs control, with Welch t-test
# ----------------------------------------------------------------------------
def _pfmt(p):
    if p < 1e-300: return "p < 1e-300"
    if p < 0.001:  return f"p = {p:.0e}".replace("e-0", "e-")
    return f"p = {p:.3f}"

def draw_boxplot(ax, df, clock, spec):
    rng = np.random.default_rng(0)
    case = df[df.group == "case"][clock].values
    ctrl = df[df.group == "control"][clock].values
    # Control on the LEFT, intervention/case on the RIGHT
    data = [ctrl, case]
    cols = [spec["ctrl_color_solid"], spec["case_color_solid"]]
    pos = [1, 1.34]                      # thin boxes, small gap between them
    bp = ax.boxplot(data, positions=pos, widths=0.24, patch_artist=True,
                    showfliers=False, medianprops=dict(color="black", lw=1.0),
                    whiskerprops=dict(lw=0.8, color=INK),
                    capprops=dict(lw=0.8, color=INK),
                    boxprops=dict(lw=0.8, edgecolor=INK))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(ALPHA["box"])
    for i, (d, c) in enumerate(zip(data, cols)):
        jit = rng.uniform(-0.06, 0.06, len(d))
        ax.scatter(np.full(len(d), pos[i]) + jit, d, s=11, facecolors=c,
                   edgecolors="black", linewidths=0.4, alpha=0.85, zorder=4,
                   clip_on=False)
    ax.set_ylabel(spec["ylabel"], fontsize=8.5)
    ax.set_xticks(pos)
    ax.set_xticklabels([spec["ctrl_short"], spec["case_short"]], fontsize=8.5)
    ax.set_xlim(0.45, 1.89)             # clear gap between y-axis and left box
    ax.tick_params(axis="y", labelsize=8)
    # Welch t-test
    tt = stats.ttest_ind(case, ctrl, equal_var=False)
    lo = min(case.min(), ctrl.min()); hi = max(case.max(), ctrl.max())
    span = hi - lo if hi > lo else 1.0
    ax.set_ylim(lo - 0.07 * span, hi + 0.24 * span)
    # significance bracket
    yb = hi + 0.085 * span; h = 0.038 * span
    ax.plot([pos[0], pos[0], pos[1], pos[1]], [yb, yb + h, yb + h, yb],
            lw=0.9, c=INK)
    ax.text((pos[0] + pos[1]) / 2, yb + h + 0.01 * span,
            f"Welch t-test\n{_pfmt(tt.pvalue)}", ha="center", va="bottom",
            fontsize=6.8, color=INK, linespacing=1.2)
    return tt, len(case), len(ctrl)

# ----------------------------------------------------------------------------
# Compose one figure
# ----------------------------------------------------------------------------
def make_figure(spec):
    retro_style()
    df, clock = build_groups(spec)
    n_case = int((df.group == "case").sum())
    n_control = int((df.group == "control").sum())
    fig, axd = plt.subplot_mosaic(
        [["schema", "box"]], figsize=figsize_nature("double", aspect=0.5),
        layout="constrained", gridspec_kw=dict(width_ratios=[1.18, 0.82]))
    draw_schematic(axd["schema"], spec, n_case, n_control)
    tt, nc, nk = draw_boxplot(axd["box"], df, clock, spec)
    trim_axes(fig, exclude=[axd["schema"]])
    fig.suptitle(spec["title"], fontsize=8.5, weight="bold", y=1.02)
    base = f"{OUT}/{spec['key']}"
    fig.savefig(f"{base}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{base}.pdf", bbox_inches="tight")
    save_for_review(fig, f"{base}_review.png", dpi=110)
    plt.close(fig)
    print(f"[{spec['key']}] case n={nc} mean={df[df.group=='case'][clock].mean():.1f} | "
          f"control n={nk} mean={df[df.group=='control'][clock].mean():.1f} | "
          f"Welch p={tt.pvalue:.2e}  -> {base}.png")
    return df, tt

# ----------------------------------------------------------------------------
# Specs
# ----------------------------------------------------------------------------
def _contains(series, *subs):
    s = series.astype(str)
    m = np.zeros(len(s), dtype=bool)
    for x in subs: m |= s.str.contains(x, case=False, regex=False)
    return m

def _eq(series, *vals):
    s = series.astype(str)
    m = np.zeros(len(s), dtype=bool)
    for v in vals: m |= s.eq(v)
    return m

# Colours: blue = case younger (age-lowering), rust = case older (age-increasing);
# control always grey.
LOW = dict(case_color=BLUE["200"], case_color_solid=BLUE["300"])
HIGH = dict(case_color=RUST["200"], case_color_solid=RUST["400"])
CTRL = dict(ctrl_color=GRAY["200"], ctrl_color_solid=GRAY["300"])
ASSAY = "Genome-wide DNA\nmethylation array"
CLOCKL = "DNAm PhenoAge clock"

def spec(key, gse, title, group_fn, source_label, case_label, ctrl_label,
         case_short, ctrl_short, direction):
    d = dict(key=key, gse=gse, title=title, group_fn=group_fn, clock="PhenoAge",
             ylabel="DNAm PhenoAge (years)", source_label=source_label,
             case_label=case_label, ctrl_label=ctrl_label,
             case_short=case_short, ctrl_short=ctrl_short,
             assay_label=ASSAY, clock_label=CLOCKL, **CTRL)
    d.update(LOW if direction == "low" else HIGH)
    return d

SPECS = {}

# ---- AGE-LOWERING (case younger) -------------------------------------------
SPECS["idh1"] = spec(
    "01_IDH1_glioma", "GSE104723",
    "IDH1 mutation: lower epigenetic age in glioma",
    lambda d: np.where(_contains(d["title"], "IDH1 mutant"), "case",
               np.where(_contains(d["title"], "wild"), "control", None)),
    "Human high-grade\nglioma tumors", "IDH1-mutant", "IDH1 wild-type",
    "IDH1-mut", "WT", "low")

SPECS["chd7"] = spec(
    "02_CHD7_CHARGE", "GSE97362",
    "CHD7 variant (CHARGE): lower epigenetic age",
    lambda d: np.where(_eq(d["harmonized.disease"], "CHD7 variant"), "case",
               np.where(_eq(d["harmonized.disease"], "Control"), "control", None)),
    "Whole blood,\nchromatin-disorder cohort", "CHD7 variant", "Control",
    "CHD7", "Control", "low")

SPECS["imprint"] = spec(
    "03_11p15_imprint", "GSE104451",
    "11p15 methylation loss: lower epigenetic age",
    lambda d: np.where(_eq(d["harmonized.genetic_information"], "11p15 LOM"), "case",
               np.where(_eq(d["harmonized.genetic_information"], "control"), "control", None)),
    "Whole blood,\nimprinting-disorder cohort", "11p15 methylation loss", "Control",
    "11p15 LOM", "Control", "low")

# ---- AGE-INCREASING (case older) -------------------------------------------
SPECS["down"] = spec(
    "04_Down_trisomy21", "GSE140344",
    "Down syndrome (trisomy 21): higher epigenetic age",
    lambda d: np.where(_contains(d["harmonized.disease"], "trisomy 21"), "case",
               np.where(_contains(d["harmonized.disease"], "control"), "control", None)),
    "Whole blood cohort", "Trisomy 21", "Euploid control",
    "Trisomy 21", "Control", "high")

SPECS["sotos"] = spec(
    "05_Sotos_NSD1", "GSE74432",
    "Sotos syndrome (NSD1): higher epigenetic age",
    lambda d: np.where(_contains(d["harmonized.disease"], "Sotos"), "case",
               np.where(_eq(d["harmonized.disease"], "Control"), "control", None)),
    "Whole blood cohort", "Sotos syndrome", "Control",
    "Sotos", "Control", "high")

SPECS["williams"] = spec(
    "06_Williams", "GSE119778",
    "Williams syndrome: higher epigenetic age",
    lambda d: np.where(_contains(d["harmonized.disease"], "Williams"), "case",
               np.where(_contains(d["harmonized.disease"], "Control"), "control", None)),
    "Peripheral blood\ncohort", "Williams syndrome", "Control",
    "Williams", "Control", "high")

SPECS["obesity"] = spec(
    "07_Obesity", "GSE67024",
    "Obesity: higher epigenetic age",
    lambda d: np.where(_eq(d["harmonized.disease"], "obese"), "case",
               np.where(_eq(d["harmonized.disease"], "control"), "control", None)),
    "Subcutaneous\nadipocytes", "Obese", "Lean control",
    "Obese", "Lean", "high")

SPECS["athero"] = spec(
    "08_Atherosclerosis", "GSE107143",
    "Atherosclerosis: higher epigenetic age",
    lambda d: np.where(_contains(d["harmonized.disease"], "Atherosclerosis"), "case",
               np.where(_contains(d["harmonized.disease"], "Healthy"), "control", None)),
    "Whole blood cohort", "Atherosclerosis", "Healthy",
    "Athero.", "Healthy", "high")

SPECS["aorta"] = spec(
    "09_Aortic_dissection", "GSE84274",
    "Aortic dissection: higher epigenetic age",
    lambda d: np.where(_eq(d["harmonized.disease"], "aortic dissection"), "case",
               np.where(_eq(d["harmonized.disease"], "normal"), "control", None)),
    "Ascending aorta\ntissue", "Aortic dissection", "Normal aorta",
    "Dissection", "Normal", "high")

SPECS["afib"] = spec(
    "10_Atrial_fibrillation", "GSE62727",
    "Atrial fibrillation: higher epigenetic age",
    lambda d: np.where(_contains(d["source_name"], "atrial fibrillation"), "case",
               np.where(_contains(d["source_name"], "normal"), "control", None)),
    "Left atrial\nheart tissue", "Atrial fibrillation", "Normal rhythm",
    "AF", "Normal", "high")

if __name__ == "__main__":
    keys = sys.argv[1:] if len(sys.argv) > 1 else list(SPECS)
    for k in keys:
        make_figure(SPECS[k])
