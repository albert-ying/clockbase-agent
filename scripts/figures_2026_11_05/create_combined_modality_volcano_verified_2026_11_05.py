#!/usr/bin/env python3
"""
Combined-modality volcano, rendered in the same style as
volcano_plot_{method}_all_categories_new.pdf in the paper.

Data: all three modalities pooled (mouse RNA-seq, human RNA-seq, human DNAm
primary-clock GrimAgeV2).

Visual encodings (matches s01_visualize_results_human.py Figure-3 block):
  - x            = log2 effect size
  - y            = −log10(FDR)
  - color        = log2FoldChange, RdBu_r divergent cmap (via colorbar)
  - shape        = MODALITY  (● mouse RNA, ■ human RNA, ▲ human DNAm)
  - size         = composite final_score (cubic-scaled, same as s01)
  - black outline layer behind every coloured point
  - non-significant points rendered as faint light-gray background
  - KDE density contour overlay (rocket_r, α = 0.4)
  - selective bracket-connector annotations on the periphery of the plot,
    top N pro-aging on the right, top N anti-aging on the left
  - Effect Size colorbar on the right; modality legend upper-left

Outputs:
  figures/volcano_combined_modalities.{pdf,png}
"""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

try:
    from aquarel import load_theme
    _theme = (load_theme("scientific").set_grid(draw=False).set_font(size=15)
              .set_ticks(direction="out").set_axis_labels(pad=10))
    _theme.apply()
except Exception:
    _theme = None

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"


CATEGORY_COLORS = {   # used only for annotation text color
    "Disease":     "#DC143C",
    "Drug":        "#4B0082",
    "Genetic":     "#228B22",
    "Environment": "#FF8C00",
    "Other":       "#808080",
}
CATEGORY_ORDER = ["Genetic", "Drug", "Disease", "Environment", "Other"]

MODALITY_MARKERS = {
    "mouse RNA-seq": "o",
    "human RNA-seq": "s",
    "human DNAm":    "^",
}
MODALITY_ORDER = list(MODALITY_MARKERS.keys())

# Per-category annotation budget, pooled across modalities. Same structure s01
# uses so the combined plot has the same category balance as the mouse paper.
ANNOTATION_COUNTS = {
    "Disease":     {"positive": 10, "negative": 3},
    "Drug":        {"positive": 10, "negative": 20},
    "Genetic":     {"positive": 8,  "negative": 8},
    "Environment": {"positive": 8,  "negative": 5},
    "Other":       {"positive": 0,  "negative": 0},
}
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
    df["_size"] = normalize_size(df["final_score"])
    for _m in ("mouse RNA-seq", "human RNA-seq", "human DNAm"):
        print(f"  {_m} : {(df['modality'] == _m).sum():,}")
    print(f"combined (Normal only): {len(df):,} rows")

    # ---- figure ----
    fig = plt.figure(figsize=(18, 10))
    ax = plt.gca()

    # Layer 0 — faint density contour from ALL points
    try:
        kde_df = df[["log2FoldChange", "_nlf"]].rename(columns={"_nlf": "neg_log10_fdr"})
        sns.kdeplot(
            data=kde_df, x="log2FoldChange", y="neg_log10_fdr",
            fill=True, cmap="rocket_r", alpha=0.35, levels=6,
            linewidths=1.0, ax=ax, zorder=4,
        )
    except Exception as e:
        print(f"  kde skipped: {e}")

    # Layer 1 — non-significant as light-gray background
    ns = df[df["fdr"] >= 0.05]
    ax.scatter(ns["log2FoldChange"], ns["_nlf"],
               s=5, c="lightgray", alpha=0.3, zorder=1, rasterized=True)

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
                   s=m["_size"] * 1.5, c="black", alpha=1.0,
                   marker=marker, zorder=5, rasterized=True)
        # colored by effect size
        scatter_handle = ax.scatter(
            m["log2FoldChange"], m["_nlf"],
            s=m["_size"], c=m["log2FoldChange"],
            cmap="RdBu_r", vmin=vmin, vmax=vmax,
            alpha=0.9, marker=marker, zorder=6, rasterized=True,
        )

    # ---- bracket-style annotations, per-category budget (pooled across modalities) ----
    ann_pool = df[df["fdr"] < 0.05].dropna(subset=["final_score"]).copy()
    ann_pool = ann_pool[ann_pool.apply(
        lambda r: r["final_score"] >= SCORE_FLOOR.get(r["modality"], 100), axis=1
    )]
    # Plausibility gate — drop wrong-direction implausibility just like s01
    implaus_mask = ann_pool["log2FoldChange"] < 0
    kw_hit = False
    for c in ("condition_name", "comparison_name", "intervention"):
        if c in ann_pool.columns:
            h = ann_pool[c].astype(str).str.contains(IMPLAUSIBLE_KEYWORDS, case=False, na=False, regex=True)
            kw_hit = h if kw_hit is False else (kw_hit | h)
    implaus_mask = implaus_mask & kw_hit
    ann_pool = ann_pool[~implaus_mask]

    ann_pool["_rank"] = ann_pool["final_score"] * np.abs(ann_pool["log2FoldChange"])
    # Normalised dedup key: lowercase + drop non-alphanumerics (handles "WY-14643" vs "WY-14,643")
    import re as _re
    def _intv_norm(s):
        return _re.sub(r"[^a-z0-9]", "", str(s).lower())
    ann_pool["_intv_key"] = ann_pool["intervention"].apply(_intv_norm)

    # 1) FORCE-INCLUDE famous interventions: pick the single highest-ranked hit
    #    per keyword, regardless of category budget — BUT only if direction
    #    matches the expected canonical direction.
    forced = []
    for kw, expect in FORCE_INCLUDE:
        mask = (
            ann_pool["intervention"].astype(str).str.contains(kw, case=False, na=False, regex=False) |
            ann_pool["condition_name"].astype(str).str.contains(kw, case=False, na=False, regex=False)
        )
        hits = ann_pool[mask]
        if expect == "neg":
            hits = hits[hits["log2FoldChange"] < 0]
        elif expect == "pos":
            hits = hits[hits["log2FoldChange"] > 0]
        hits = hits.sort_values("_rank", ascending=False)
        if not hits.empty:
            forced.append(hits.head(1))
    forced_df = pd.concat(forced, ignore_index=True) if forced else pd.DataFrame(columns=ann_pool.columns)
    forced_df = forced_df.drop_duplicates(subset=["_intv_key"])
    forced_keys = set(forced_df["_intv_key"])
    print(f"force-included famous interventions: {len(forced_df)}")

    # 2) Per-category, per-direction picking for the rest
    pos_picks, neg_picks = [forced_df[forced_df["log2FoldChange"] > 0]], [forced_df[forced_df["log2FoldChange"] < 0]]
    for cat, budget in ANNOTATION_COUNTS.items():
        cat_rows = ann_pool[(ann_pool["condition_category"] == cat)
                            & (~ann_pool["_intv_key"].isin(forced_keys))]
        if cat_rows.empty:
            continue
        cat_rows = cat_rows.sort_values("_rank", ascending=False).drop_duplicates(subset=["_intv_key"], keep="first")
        # How many force-picks are already in this category+direction? Subtract from budget.
        forced_pos_in_cat = ((forced_df["condition_category"] == cat) & (forced_df["log2FoldChange"] > 0)).sum()
        forced_neg_in_cat = ((forced_df["condition_category"] == cat) & (forced_df["log2FoldChange"] < 0)).sum()
        pos = cat_rows[cat_rows["log2FoldChange"] > 0].head(max(0, budget["positive"] - forced_pos_in_cat))
        neg = cat_rows[cat_rows["log2FoldChange"] < 0].head(max(0, budget["negative"] - forced_neg_in_cat))
        pos_picks.append(pos)
        neg_picks.append(neg)
    ann_pos = pd.concat(pos_picks, ignore_index=True) if pos_picks else pd.DataFrame()
    ann_neg = pd.concat(neg_picks, ignore_index=True) if neg_picks else pd.DataFrame()
    if not ann_pos.empty: ann_pos = ann_pos.drop_duplicates(subset=["_intv_key"], keep="first")
    if not ann_neg.empty: ann_neg = ann_neg.drop_duplicates(subset=["_intv_key"], keep="first")
    # Persist picked annotations for audit
    _all = pd.concat([ann_pos.assign(_side="pos"), ann_neg.assign(_side="neg")], ignore_index=True)
    all_annotations_cols = ["_side","modality","condition_category","gse_id",
                             "condition_name","intervention","log2FoldChange","fdr","final_score"]
    all_annotations_cols = [c for c in all_annotations_cols if c in _all.columns]
    _all[all_annotations_cols].to_csv("result/combined_volcano_annotations.csv", index=False)
    print(f"  audit list: result/combined_volcano_annotations.csv")
    print(f"annotations picked: {len(ann_pos)} pro-aging + {len(ann_neg)} anti-aging")
    print("  by category:")
    for df_side, side in [(ann_pos, "pos"), (ann_neg, "neg")]:
        if df_side.empty: continue
        print(f"    {side}: {df_side['condition_category'].value_counts().to_dict()}")

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    new_xlim = (xlim[0] - x_range * 0.22, xlim[1] + x_range * 0.22)
    ax.set_xlim(new_xlim)
    # Cap annotation y-range at 10 so labels don't overlap with the upper-left
    # legends (category + modality), and don't drift off into ultra-significant
    # territory where the label spacing would be wasted.
    y_bottom = -np.log10(0.05) + 0.2   # just above the FDR line
    y_top    = 10.0

    # right side (pro-aging) — connector line is hidden behind opaque white
    # label bbox (drawn AFTER the line via higher zorder).
    if len(ann_pos):
        ann_pos_sorted = ann_pos.sort_values("_nlf")
        y_positions = np.linspace(y_bottom, y_top, len(ann_pos_sorted))
        x_text = new_xlim[1] - x_range * 0.02
        for (_, r), y_pos in zip(ann_pos_sorted.iterrows(), y_positions):
            label = str(r.get("intervention") or r.get("condition_name") or "")[:40]
            color = CATEGORY_COLORS.get(r["condition_category"], "#333")
            xh1 = r["log2FoldChange"] + x_range * 0.01
            xsl = xh1 + x_range * 0.025
            xt  = x_text - x_range * 0.005
            ax.plot([r["log2FoldChange"], xh1, xsl, xt],
                    [r["_nlf"], r["_nlf"], y_pos, y_pos],
                    color="black", ls=":", alpha=0.4, lw=1.0, zorder=7)
            ax.text(x_text, y_pos, label, ha="right", va="center",
                    fontsize=9, color=color, weight="bold", zorder=10,
                    path_effects=[pe.withStroke(linewidth=3, foreground="white")],
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                              edgecolor="none", alpha=1.0))

    # left side (anti-aging)
    if len(ann_neg):
        ann_neg_sorted = ann_neg.sort_values("_nlf")
        y_positions = np.linspace(y_bottom, y_top, len(ann_neg_sorted))
        x_text = new_xlim[0] + x_range * 0.02
        for (_, r), y_pos in zip(ann_neg_sorted.iterrows(), y_positions):
            label = str(r.get("intervention") or r.get("condition_name") or "")[:40]
            color = CATEGORY_COLORS.get(r["condition_category"], "#333")
            xh1 = r["log2FoldChange"] - x_range * 0.01
            xsl = xh1 - x_range * 0.025
            xt  = x_text + x_range * 0.005
            ax.plot([r["log2FoldChange"], xh1, xsl, xt],
                    [r["_nlf"], r["_nlf"], y_pos, y_pos],
                    color="black", ls=":", alpha=0.4, lw=1.0, zorder=7)
            ax.text(x_text, y_pos, label, ha="left", va="center",
                    fontsize=9, color=color, weight="bold", zorder=10,
                    path_effects=[pe.withStroke(linewidth=3, foreground="white")],
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                              edgecolor="none", alpha=1.0))

    # ---- labels & title ----
    ax.axhline(-np.log10(0.05), ls="--", color="#666", lw=0.7, zorder=0)
    ax.set_xlabel("Effect Size")
    ax.set_ylabel("−log10(Adjusted P-value)")
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
    ax.set_title(
        f"ClockBase Agent, combined aging-intervention atlas\n"
        f"Total comparisons: {total:,}   "
        f"(mouse RNA-seq: {n_mouse:,} · "
        f"human RNA-seq: {n_hr:,} · "
        f"human DNAm: {n_hd:,})   |   "
        f"Significant (FDR<0.05): {sig_count:,} ({100*sig_count/total:.1f}%)",
        fontsize=12,
    )

    # ---- colorbar (Effect Size) + size-reference legend (Final Score) ----
    if scatter_handle is not None:
        cbar = plt.colorbar(scatter_handle, ax=ax, pad=0.02, shrink=0.55,
                            anchor=(0.0, 0.95), panchor=(0.0, 1.0))
        cbar.set_label("Effect Size", rotation=270, labelpad=20)

    # Size-reference legend: marker sizes at representative final_score values.
    # Positioned upper-right (inside plot, above colorbar zone).
    size_values = [70, 90, 110, 130]
    size_handles = [
        Line2D([0], [0], marker="o", color="black",
               markerfacecolor="#bbb",
               markersize=np.sqrt(normalize_size(pd.Series([v])).iloc[0]),
               linestyle="None", label=f"{v}")
        for v in size_values
    ]
    size_legend = ax.legend(
        handles=size_handles, title="Final score",
        loc="upper right", frameon=False,
        fontsize=9.5, title_fontsize=10,
        ncol=len(size_values), columnspacing=1.5, handletextpad=0.4, borderpad=0.6,
    )
    ax.add_artist(size_legend)

    # ---- legends: category (color-code for annotations) + modality (shape) ----
    cat_handles = [
        Line2D([0], [0], marker="s", color="w",
               markerfacecolor=CATEGORY_COLORS[c], markeredgecolor="#333",
               markersize=11, linestyle="None", label=c)
        for c in CATEGORY_ORDER if c in CATEGORY_COLORS
    ]
    mod_handles = [
        Line2D([0], [0], marker=MODALITY_MARKERS[m], color="black",
               markerfacecolor="#999", markersize=10, linestyle="None", label=m)
        for m in MODALITY_ORDER
    ]
    cat_legend = ax.legend(
        handles=cat_handles, title="Intervention category (annotation color)",
        loc="upper left", frameon=False, fontsize=9.5, title_fontsize=10,
    )
    ax.add_artist(cat_legend)
    ax.legend(
        handles=mod_handles, title="Modality (point shape)",
        loc="upper left", bbox_to_anchor=(0.22, 1.0),
        frameon=False, fontsize=9.5, title_fontsize=10,
    )

    if _theme is not None:
        try:
            _theme.apply_transforms()
        except Exception:
            pass
    plt.tight_layout()

    Path("figures").mkdir(exist_ok=True)
    pdf = "figures/volcano_combined_modalities.pdf"
    png = "figures/volcano_combined_modalities.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    print(f"\n✓ {pdf}")
    print(f"✓ {png}")


if __name__ == "__main__":
    raise SystemExit("Run restore_figure3_2026_11_05.py; the original entry point is disabled to protect July outputs.")

# Verified restoration adapter. July constants, density overlay, black-outlined
# modality marks and peripheral dotted leaders above are retained. Historical
# scores and plausibility-based annotation selection are replaced by released
# statistical ranks; unavailable scores are not assigned surrogate values.
def draw_verified(fig, ax, full, df, ranking):
    from matplotlib.colors import Normalize
    from scipy.stats import gaussian_kde
    from matplotlib.cm import ScalarMappable
    from matplotlib.ticker import FixedLocator
    sig=df[df.fdr<.05];ns=df[df.fdr>=.05]
    norm=Normalize(-45,45)
    cloud=df[['effect_size','neg_log10_fdr']].to_numpy().T
    # Same Scott bandwidth KDE as July's seaborn KDE, evaluated within the
    # readable main domain. This estimates density from observations only.
    gx=np.linspace(-45,45,160);gy=np.linspace(0,12,100)
    xx,yy=np.meshgrid(gx,gy)
    import hashlib
    cache=Path(__file__).resolve().parents[2]/'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-05/figure3/panel_a_kde.npz'
    digest=hashlib.sha256(cloud.tobytes()).hexdigest()
    if cache.exists() and str(np.load(cache)['source_sha256'])==digest:
        zz=np.load(cache)['density']
    else:
        kernel=gaussian_kde(cloud)
        zz=kernel(np.vstack([xx.ravel(),yy.ravel()])).reshape(xx.shape)
        np.savez_compressed(cache,x=gx,y=gy,density=zz,source_sha256=digest,bandwidth='Scott')
    lev=np.linspace(zz.max()*.08,zz.max(),6)
    # Filled density scale in the July rocket family has its own explicit key.
    import seaborn as sns
    density=ax.contourf(xx,yy,zz,levels=lev,cmap='rocket_r',alpha=.35,zorder=0)
    ax.scatter(ns.effect_size,ns.neg_log10_fdr,s=1.8,c='lightgray',alpha=.30,zorder=1,edgecolors='none')
    for m in MODALITY_ORDER:
        s=sig[sig.modality==m]
        ax.scatter(s.effect_size,s.neg_log10_fdr,s=10,c='black',marker=MODALITY_MARKERS[m],edgecolors='none',zorder=3)
        scatter_handle=ax.scatter(s.effect_size,s.neg_log10_fdr,s=8,c=s.effect_size,cmap='RdBu_r',norm=norm,marker=MODALITY_MARKERS[m],edgecolors='none',alpha=.90,zorder=4)
    ax.set_xlim(-45,45);ax.set_ylim(0,12)
    ax.set_xticks([-40,-20,0,20,40]);ax.set_yticks([0,2,4,6,8,10,12])
    ax.set_xlabel("Cohen’s d",labelpad=2);ax.set_ylabel('−log10(BH-FDR)',labelpad=3)
    ax.axhline(-np.log10(.05),ls='--',color='#666',lw=.45,zorder=2)
    for sp in ['top','right']:ax.spines[sp].set_visible(False)
    ax.tick_params(length=2,width=.5,pad=2)
    # Rank-selected label pool, three per category in each observed direction.
    picks=[]
    for cat in ['Genetic','Drug','Disease','Environment']:
        catrows=ranking[ranking.condition_category==cat].sort_values('rank_bound',ascending=False)
        for direction in [-1,1]:
            picked=catrows[np.sign(catrows.hedges_g)==direction].head(3).copy()
            picked['side']=direction;picks.append(picked)
    labels=pd.concat(picks,ignore_index=True)
    # Anchor to matching verified plotting rows instead of ranking's legacy
    # aliases (notably its fdr column, which is historical).
    keys=['modality','gse_id','comparison_name','condition_name']
    prov=pd.read_csv(Path(__file__).resolve().parents[2]/'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-05/data/verified_plotting_row_provenance_2026-11-05.csv',low_memory=False)
    looked=[]
    for _,r in labels.iterrows():
        matched=prov
        for k in keys+['clock']:matched=matched[matched[k].fillna('')==(''if pd.isna(r[k])else r[k])]
        matched=matched[np.isclose(matched.v3_p,r.v3_p,rtol=1e-8,atol=0)]
        if matched.empty:raise ValueError('Ranked label has no verified atlas match: '+r.gse_id+' '+r.label)
        rr=r.to_dict();rr['cohen_d']=float(matched.iloc[0].v3_cohen_d);rr['bh_fdr']=float(matched.iloc[0].v3_fdr);rr['neg_log10_bh_fdr']=-np.log10(rr['bh_fdr']);looked.append(rr)
    labels=pd.DataFrame(looked)
    for direction in [-1,1]:
        side=labels[labels.side==direction].sort_values('neg_log10_bh_fdr')
        ys=np.linspace(2.0,10.6,len(side))
        for (_,r),y in zip(side.iterrows(),ys):
            textx=direction*43.5
            # Deterministic spaced peripheral labels reproduce July's bracket
            # placement and leave the point cloud open at the centre.
            label=str(r.label)
            import textwrap
            label='\n'.join(textwrap.wrap(label,width=30,break_long_words=False))
            endx=textx
            ax.plot([r.cohen_d,r.cohen_d+direction*.7,r.cohen_d+direction*1.8,endx],[r.neg_log10_bh_fdr,r.neg_log10_bh_fdr,y,y],ls=':',lw=.35,color='#444',alpha=.5,zorder=5,clip_on=True)
            ax.text(textx,y,label,ha='right'if direction>0 else'left',va='center',fontsize=6,color=CATEGORY_COLORS[r.condition_category],fontweight='bold',zorder=6,bbox=dict(facecolor='white',edgecolor='none',pad=.3,alpha=.94))
    # Every observation appears in the full-range view; no winsorisation.
    for m in MODALITY_ORDER:
        s=df[df.modality==m];col=np.where(s.fdr<.05,np.where(s.effect_size<0,'#2166ac','#b2182b'),'#d5d5d5')
        full.scatter(s.effect_size,s.neg_log10_fdr,s=.6,c=col,marker=MODALITY_MARKERS[m],edgecolors='none',alpha=.45)
    full.set_xlim(-75,140);full.set_ylim(0,235);full.set_xticks([-50,0,100]);full.set_yticks([0,100,200]);full.tick_params(length=2,width=.5,pad=1)
    full.axhline(-np.log10(.05),ls='--',lw=.3,color='#777')
    full.text(.98,.96,'Full range',transform=full.transAxes,fontsize=6,ha='right',va='top')
    for sp in ['top','right']:full.spines[sp].set_visible(False)
    # Horizontal effect scale outside the volcano, never over the cloud.
    cax=fig.add_axes([.420,.626,.13,.009])
    cb=fig.colorbar(scatter_handle,cax=cax,orientation='horizontal',ticks=[-40,0,40],extend='both')
    cax._ying_is_colorbar=True;cb.set_label("Cohen’s d",fontsize=7,labelpad=1);cb.ax.xaxis.set_label_position('top');cb.ax.tick_params(labelsize=7,length=2,pad=1);cb.outline.set_linewidth(.35);cb.solids.set_rasterized(False)
    dax=fig.add_axes([.420,.579,.13,.009])
    dcb=fig.colorbar(density,cax=dax,orientation='horizontal',ticks=[lev[0],lev[-1]])
    dax._ying_is_colorbar=True;dcb.set_ticklabels([f'{lev[0]:.3g}',f'{lev[-1]:.3g}']);dcb.set_label('Density',fontsize=7,labelpad=1);dcb.ax.xaxis.set_label_position('top');dcb.ax.tick_params(labelsize=7,length=2,pad=1);dcb.outline.set_linewidth(.35);dcb.solids.set_rasterized(False);cb.solids.set_rasterized(False)
    # Category labels retain the July annotation colours without another
    # competing legend inside the cloud.
    for x,cat in zip([.10,.30,.50,.75],['Genetic','Drug','Disease','Environment']):
        ax.text(x,.985,cat,transform=ax.transAxes,fontsize=7,color=CATEGORY_COLORS[cat],ha='center',va='top',fontweight='bold',bbox=dict(facecolor='white',edgecolor='none',pad=.2))
    return labels
