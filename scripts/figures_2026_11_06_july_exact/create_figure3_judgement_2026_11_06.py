#!/usr/bin/env python3
"""Figure 3g shows admitted comparisons in statistical order.

The source is followup-judgment-001/judgement_age-lowering.csv. Admission uses
agreed answers on age at intervention start, health and untreated controls.
The first comparison of each intervention identity and then each GEO series is
drawn in statistical order. Bars give the absolute adjusted confidence bound
nearest zero, as in the statistical panels. Ouabain keeps its dagger and bold
label. Full source comparisons and statistical ranks are recorded in CSV.
"""
import importlib.util
import json
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import blended_transform_factory

# Colour roles of Figure 3 (FIG3_PALETTE; palette B when unset, "july" reproduces the July colours).
_spec = importlib.util.spec_from_file_location(
    "fig3_palette", Path(__file__).resolve().parent / "fig3_palette_2026_11_08.py")
fig3_palette = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fig3_palette)
PAL = fig3_palette.active()

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"
# Print size (7 November 2026): drawn at 99.6 x 38.2 mm beside panel f and placed at 100% in Figure 3, Arial throughout.
plt.rcParams.update({
    "font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "grid.linewidth": 0.4, "hatch.linewidth": 0.4,
})
MM = 1 / 25.4  # inches per millimetre
# Shared with panel f (create_figure3_replication_2026_11_07.py) except the width: same height, axes top and bottom,
# header baseline and row spacing, so the two axes have the same height and baseline in the assembled row.
PANEL = {"width_mm": 99.6, "height_mm": 38.2, "ax_top_mm": 9.6, "ax_bottom_mm": 6.7, "right_mm": 1.0,
         "header_x_mm": 4.0, "header_baseline_mm": 7.8 / 72 * 25.4, "header": 8, "tick": 7, "axis": 7,
         "column": 6.5, "column_gap_mm": 2.0, "head_above_axes_mm": 0.9, "marker": 3.5, "marker_edge": 0.4,
         "bar_edge": 0.3, "edge_mm": 1.7, "rows": 5}
HEADER = ("Judgement ranking", "in old, otherwise healthy animals vs untreated controls")
XLABEL = "|Hedges' g| (CI bound nearest zero)"
INK = PAL["ink"]
KEYS = ["gse_id", "comparison_name", "condition_name"]
SPECIES = {"mouse RNA-seq": "mice"}
TESTED_KEY = ("GSE122080", "Ouabain vs Saline in old mice")  # followed up in aged mice (panel a, dagger)
TESTED_MARK = " †"
RANK_HEAD = "Statistical rank"

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
COLOR_PRO  = PAL["higher"]  # higher clock output (July: #D55E00)
COLOR_ANTI = PAL["lower"]   # lower clock output (July: #0072B2)
EDGE_PRO, EDGE_ANTI = PAL["higher_edge"], PAL["lower_edge"]  # bar outline and hatch strokes (July: #111)


def require(*paths):
    for p in paths:
        if not Path(p).exists():
            raise SystemExit(f"Cannot proceed - required data file not found: {p}")


def box_mm(fig, artist):
    """Window extent in mm from the panel's top-left corner: (left, top, right, bottom)."""
    fig.canvas.draw()
    bb = artist.get_window_extent(fig.canvas.get_renderer())
    px, H = fig.dpi * MM, PANEL["height_mm"]
    return np.array([bb.x0 / px, H - bb.y1 / px, bb.x1 / px, H - bb.y0 / px])


def width_mm(fig, text, **kw):
    probe = fig.text(0, 0, text, **kw)
    left, _, right, _ = box_mm(fig, probe)
    probe.remove()
    return right - left


def ranking_header(fig, bold_text, rest_text):
    """One header line in the panel's top line: the ranking's name in bold, its rule in regular weight."""
    W, H = PANEL["width_mm"], PANEL["height_mm"]
    y = 1 - PANEL["header_baseline_mm"] / H
    size = PANEL["header"]
    bold = fig.text(PANEL["header_x_mm"] / W, y, bold_text, fontsize=size, fontweight="bold", color=INK,
                    ha="left", va="baseline")
    gap = width_mm(fig, "a a", fontsize=size) - width_mm(fig, "aa", fontsize=size)
    x2 = PANEL["header_x_mm"] + box_mm(fig, bold)[2] - box_mm(fig, bold)[0] + gap
    rest = fig.text(x2 / W, y, rest_text, fontsize=size, color=INK, ha="left", va="baseline")
    assert box_mm(fig, rest)[2] <= W - 0.5, "header runs past the panel"
    return [bold, rest]


def column_widths(fig, columns):
    kw = dict(fontsize=PANEL["column"])
    return [max(width_mm(fig, s, **kw) for s in list(heads) + list(values) if s) for heads, values, _ in columns]


def draw_columns(fig, ax, columns, widths, y_values):
    """Columns right of the axes, ending right_mm from the panel's right edge, one value per bar row; each column is
    centred or right-aligned (heads and values alike). Heads are bottom-aligned: their last line sits
    head_above_axes_mm above the axes top."""
    W, H = PANEL["width_mm"], PANEL["height_mm"]
    trans = blended_transform_factory(fig.transFigure, ax.transData)
    lead = PANEL["column"] * 1.2 / 72 * 25.4
    base = PANEL["ax_top_mm"] - PANEL["head_above_axes_mm"]
    x_right = W - PANEL["right_mm"]
    spans = []
    for w in reversed(widths):
        spans.append((x_right - w, x_right))
        x_right -= w + PANEL["column_gap_mm"]
    artists = []
    for (heads, values, align), (x0, x1) in zip(columns, spans[::-1]):
        x = (x0 + x1) / 2 if align == "center" else x1
        for k, line in enumerate(reversed(heads)):
            artists.append(fig.text(x / W, 1 - (base - k * lead) / H, line, fontsize=PANEL["column"], color=INK,
                                    ha=align, va="baseline"))
        for v, y in zip(values, y_values):
            artists.append(fig.text(x / W, y, v, transform=trans, fontsize=PANEL["column"], color=INK,
                                    ha=align, va="center_baseline"))
    return artists, x_right + PANEL["column_gap_mm"]


def place_axis(fig, ax, right_mm):
    """Fix the axes in millimetres: it begins after the widest bar label and ends a column gap before the column."""
    W, H = PANEL["width_mm"], PANEL["height_mm"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label = max(t.get_window_extent(renderer).width for t in ax.get_yticklabels()) / (fig.dpi * MM)
    left = 0.3 + label + (plt.rcParams["ytick.major.size"] + plt.rcParams["ytick.major.pad"]) / 72 * 25.4
    ax.set_position([left / W, PANEL["ax_bottom_mm"] / H, (right_mm - PANEL["column_gap_mm"] - left) / W,
                     (H - PANEL["ax_top_mm"] - PANEL["ax_bottom_mm"]) / H])


def row_geometry(n, bar_mm):
    """Bar height (data units) and y limits that draw n bars at panel b-e thickness on this panel's row spacing."""
    axis_mm = PANEL["height_mm"] - PANEL["ax_top_mm"] - PANEL["ax_bottom_mm"]
    pitch_mm = (axis_mm - 2 * PANEL["edge_mm"] - bar_mm) / (PANEL["rows"] - 1)
    half = axis_mm / (2 * pitch_mm)
    return bar_mm / pitch_mm, ((n - 1) / 2 - half, (n - 1) / 2 + half)


def check_text_layout(fig, ax):
    """No two texts overlap, none leaves the panel, none is under 6 pt, none touches a bar end marker."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    lo, hi = ax.get_xlim()
    texts = list(fig.texts) + [ax.xaxis.label] + list(ax.get_yticklabels())
    texts += [t for v, t in zip(ax.get_xticks(), ax.get_xticklabels()) if lo - 1e-9 <= v <= hi + 1e-9]
    texts = [t for t in texts if t.get_visible() and t.get_text().strip()]
    assert all(t.get_fontsize() >= 6 for t in texts)
    boxes = [t.get_window_extent(renderer) for t in texts]
    eps = 0.1 * fig.dpi * MM
    for t, b in zip(texts, boxes):
        assert b.x0 >= 0 and b.y0 >= 0 and b.x1 <= fig.bbox.x1 and b.y1 <= fig.bbox.y1, t.get_text()
    others = boxes + [ln.get_window_extent(renderer) for ln in ax.lines]
    names = [t.get_text() for t in texts] + ["end marker"] * len(ax.lines)
    for i in range(len(boxes)):
        for j in range(i + 1, len(others)):
            a, b = boxes[i], others[j]
            assert not (a.x0 + eps < b.x1 and b.x0 + eps < a.x1 and a.y0 + eps < b.y1 and b.y0 + eps < a.y1), \
                (names[i], names[j])


def qualifier(row):
    """'(old mice)' or '(aged mice)' from the comparison name, for a label that panels b-e use for another comparison."""
    m = re.search(r"\b(old|aged)\b", row["comparison_name"], flags=re.IGNORECASE)
    assert m, f"no age word in {row['comparison_name']!r} to tell the bar apart from panels b-e"
    return f"{row['display']} ({m.group(1).lower()} {SPECIES[row['modality']]})"


def make_label(row):
    iv = row.get("intervention")
    if isinstance(iv, str) and iv.strip():
        return iv.strip()[:45]
    return str(row.get("condition_name", ""))[:45]


def main():
    verifier = Path(os.environ["CLOCKBASE_VERIFIER_ROOT"])
    source = Path(os.environ.get('CLOCKBASE_JUDGEMENT') or verifier / 'followup-judgment-001') / 'judgement_age-lowering.csv'
    judged_path = Path(os.environ.get('CLOCKBASE_JUDGEMENT') or verifier / 'followup-judgment-001') / 'ranked_with_judgement.csv'
    be_path = Path("result/top_hits_combined_modalities_v2.csv")
    f_path = Path("result/figure3_replication_picks.csv")
    geometry_path = Path("result/figure3_be_geometry.json")
    require(source, judged_path, be_path, f_path, geometry_path)
    bar_mm = json.loads(geometry_path.read_text())["bar_mm"]
    import dedup_interventions_2026_11_06 as dedup
    meta_path = Path(os.environ.get('CLOCKBASE_META') or verifier / 'intervention-meta-2026-10-17') / 'primary/meta_input_effects.csv'
    require(meta_path)
    admitted = pd.read_csv(source).sort_values("statistical_rank", kind="stable")
    assert admitted["judgement_admitted"].all()
    assert admitted["judgement_rank"].tolist() == list(range(1, len(admitted) + 1))
    assert admitted["rank_bound"].abs().is_monotonic_decreasing
    sub = dedup.identity(admitted.drop(columns=["canonical_target"], errors="ignore"), pd.read_csv(meta_path))
    assert (sub["hedges_g"] < 0).all()
    assert sub["age_at_start"].eq("OLD").all()
    assert sub["background"].eq("HEALTHY").all()
    assert sub["control"].eq("UNTREATED").all()
    nearest = np.where(sub["g_ci_low"] > 0, sub["g_ci_low"],
                       np.where(sub["g_ci_high"] < 0, -sub["g_ci_high"], 0))
    assert np.allclose(nearest, sub["rank_bound"].abs(), rtol=0, atol=1e-12)
    sub["statistical_score"] = sub["rank_bound"].abs()
    sub["log2FoldChange"] = sub["hedges_g"]
    # The same identity and one-per-series display rule as the statistical panels.
    sub = sub.drop_duplicates("display_key").drop_duplicates("gse_id").copy()
    sub["display"] = sub["label"].str.replace(r" at ZT\d+", "", regex=True)
    sub.loc[sub["comparison_name"].eq("Overload vs Sham (Senolytic)"), "display"] = "Overload in senolytic-treated mice"
    # A label that panels b-e give to a different comparison is qualified from the comparison name.
    be = pd.read_csv(be_path)
    same_label = [be[be["display"].str.casefold() == d.casefold()] for d in sub["display"]]
    other = [len(b) > 0 and not ((b["gse_id"] == r["gse_id"]) & (b["comparison_name"] == r["comparison_name"])).any()
             for b, (_, r) in zip(same_label, sub.iterrows())]
    sub["display"] = [qualifier(r) if o else r["display"] for o, (_, r) in zip(other, sub.iterrows())]
    assert not set(sub.loc[other, "display"]) & set(be["display"])
    # the followed-up comparison carries the dagger of panel a
    tested = (sub["gse_id"] == TESTED_KEY[0]) & (sub["comparison_name"] == TESTED_KEY[1])
    assert tested.sum() == 1 and sub.loc[tested, "display"].iloc[0] == "Ouabain"
    sub["tested"] = tested
    sub.loc[tested, "display"] = sub.loc[tested, "display"] + TESTED_MARK
    # Statistical rank of the comparison drawn, the Table 3 rank: its place among the lower-output comparisons
    judged = pd.read_csv(judged_path, low_memory=False)
    lower = judged[judged["hedges_g"] < 0].copy()
    order = lower.iloc[np.argsort(-lower["rank_bound"].abs().to_numpy(), kind="stable")]
    assert order["statistical_rank"].tolist() == list(range(1, len(lower) + 1))
    n_lower = len(lower)
    check = sub.merge(lower[KEYS + ["statistical_rank", "rank_bound"]], on=KEYS, how="left", suffixes=("", "_ranked"),
                      validate="one_to_one")
    assert check["statistical_rank_ranked"].notna().all()
    assert (check["statistical_rank_ranked"].astype(int).values == sub["statistical_rank"].astype(int).values).all()
    assert np.allclose(check["rank_bound_ranked"], sub["rank_bound"])
    sub["rank_text"] = sub["statistical_rank"].astype(int).astype(str)
    sub["rank_of"] = n_lower
    sub["rank_source_file"] = str(judged_path)
    sub.to_csv("result/figure3_judgement_picks.csv", index=False)

    # The rank printed beside every bar of f and g, with the comparison it comes from.
    f = pd.read_csv(f_path)
    assert (f["rank_of"] == n_lower).all()
    ranks = pd.concat([
        pd.DataFrame({"panel": "f", "bar": f["display"], "bar_source": f["canonical_target"].map(
            lambda t: f"intervention-meta-2026-10-17/primary/ranked_replicated_age-lowering.csv canonical_target={t}"),
            "statistical_rank": f["statistical_rank"].astype("Int64"), "printed": f["rank_text"], "of": f["rank_of"],
            "rank_gse_id": f["rank_gse_id"], "rank_comparison_name": f["rank_comparison_name"],
            "rank_condition_name": f["rank_condition_name"],
            "rule": "smallest statistical_rank among lower-output ranked comparisons with this canonical_target",
            "rank_source_file": f["rank_source_file"]}),
        pd.DataFrame({"panel": "g", "bar": sub["display"], "bar_source": [
            f"followup-judgment-001/judgement_age-lowering.csv {g} | {c}"
            for g, c in zip(sub["gse_id"], sub["comparison_name"])],
            "statistical_rank": sub["statistical_rank"].astype("Int64"), "printed": sub["rank_text"], "of": sub["rank_of"],
            "rank_gse_id": sub["gse_id"], "rank_comparison_name": sub["comparison_name"],
            "rank_condition_name": sub["condition_name"],
            "rule": "statistical_rank of the comparison drawn", "rank_source_file": sub["rank_source_file"]}),
    ], ignore_index=True)
    ranks.to_csv("result/figure3_statistical_ranks.csv", index=False)

    fig, ax = plt.subplots(figsize=(PANEL["width_mm"] * MM, PANEL["height_mm"] * MM))
    sub = sub.sort_values("statistical_rank", ascending=False)
    labels = [f"{r['display']}" for _, r in sub.iterrows()]
    values = sub["statistical_score"].values
    colors = [COLOR_PRO if v > 0 else COLOR_ANTI for v in sub["log2FoldChange"].values]
    edges = [EDGE_PRO if v > 0 else EDGE_ANTI for v in sub["log2FoldChange"].values]
    hatches = [MODALITY_HATCH[m] for m in sub["modality"].values]

    y_pos = np.arange(len(sub))
    assert len(sub) == PANEL["rows"]
    bar_height, ylim = row_geometry(len(sub), bar_mm)
    bars = ax.barh(y_pos, values, height=bar_height, color=colors, edgecolor=edges, linewidth=PANEL["bar_edge"])
    ax.set_ylim(*ylim)
    for bar, h in zip(bars, hatches):
        bar.set_hatch(h * 2)  # hatch spacing is fixed in inches; doubled to keep the July density at print size

    # Modality markers match the statistical panels.
    x_pad = (values.max() - values.min()) * 0.02 if values.max() > values.min() else 0.5
    for y, (_, r) in zip(y_pos, sub.iterrows()):
        ax.plot(r["statistical_score"] + x_pad, y,
                marker=MODALITY_MARKERS[r["modality"]],
                markersize=PANEL["marker"], color=PAL["marker_edge"], markerfacecolor=PAL["marker_face"],
                markeredgewidth=PANEL["marker_edge"], zorder=5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=PANEL["tick"])
    for tick, (_, r) in zip(ax.get_yticklabels(), sub.iterrows()):
        if r["tested"]:
            tick.set_fontweight("bold")
    ax.set_xlabel(XLABEL, fontsize=PANEL["axis"], labelpad=1.5)
    ax.grid(axis="x", alpha=0.3, linestyle=":")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=PANEL["tick"])
    ax.tick_params(axis="x", labelsize=PANEL["tick"])

    columns = [([RANK_HEAD, f"(of {n_lower} decreases)"], list(sub["rank_text"]), "right")]
    widths = column_widths(fig, columns)
    _, first_column = draw_columns(fig, ax, columns, widths, y_pos)
    place_axis(fig, ax, first_column)
    header = ranking_header(fig, *HEADER)
    heads_top = min(box_mm(fig, t)[1] for t in fig.texts if t.get_text() == RANK_HEAD)
    assert heads_top >= max(box_mm(fig, t)[3] for t in header) + 0.5
    check_text_layout(fig, ax)
    Path("figures").mkdir(exist_ok=True)
    fig.savefig("figures/judgement_barplot_combined_modalities.pdf", dpi=300)
    fig.savefig("figures/judgement_barplot_combined_modalities.png", dpi=300)
    print(ranks[["panel", "bar", "printed", "rank_gse_id", "rank_comparison_name"]].to_string())


if __name__ == "__main__":
    main()
