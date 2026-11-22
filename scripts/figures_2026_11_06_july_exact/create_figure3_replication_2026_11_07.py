#!/usr/bin/env python3
"""Figure 3f: the replication ranking, drawn with a copy of the July horizontal-bar block.

The interventions with lower clock outputs that are robust across independent projects (BH-FDR < 0.05 and a 95%
interval of the pooled Hedges' g excluding zero, results from at least two projects; Results, Replication across
independent experiments), read from the primary meta-analysis
(intervention-meta-2026-10-17/primary/ranked_replicated_age-lowering.csv). Bar length is the bound of the pooled
interval nearest zero, the quantity of panels b-e and g. Three columns right of the bars give the number of projects
and GEO series pooled and the best single rank of the intervention, headed "Best single rank (of N decreases)": a
pooled effect has no rank of its own, so the column gives the best rank of a single comparison of the intervention
among the N ranked comparisons with lower outputs. The followed-up intervention (PLX5622) is bold and carries the
dagger of panel a; the key above panels b-e defines it.

Statistical rank: the best (smallest) statistical_rank among the ranked comparisons with lower outputs
(followup-judgment-001/ranked_with_judgement.csv, hedges_g < 0, ranked 1 to N by the absolute CI bound of Hedges' g
nearest zero) whose canonical_target equals the meta-analysis canonical_target. canonical_target is the column the
meta-analysis pools by (meta_input_effects.csv); for every ranked comparison it equals meta_target of
result/volcano_ranked_comparisons_july.csv where that is set, and the comparison giving the rank is one of the
meta-analysis input effects, from a series the intervention pooled. An intervention with no ranked lower-output
comparison reads "not ranked" (none of its comparisons is among the N). Nothing is typed in: the rows, their order and every printed number come from these
tables.

Outputs: figures/replication_barplot_combined_modalities.{pdf,png}, result/figure3_replication_picks.csv
"""
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator
from matplotlib.transforms import blended_transform_factory
import importlib.util

# Colour roles of Figure 3 (FIG3_PALETTE; palette B when unset, "july" reproduces the July colours).
_spec = importlib.util.spec_from_file_location(
    "fig3_palette", Path(__file__).resolve().parent / "fig3_palette_2026_11_08.py")
fig3_palette = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fig3_palette)
PAL = fig3_palette.active()

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["savefig.facecolor"] = "white"
# Print size (7 November 2026): drawn at 80.4 x 38.2 mm beside panel g and placed at 100% in Figure 3, Arial throughout.
plt.rcParams.update({
    "font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "grid.linewidth": 0.4, "hatch.linewidth": 0.4,
})
MM = 1 / 25.4  # inches per millimetre
# Panels f and g share the height, the axes top and bottom, the header baseline and the row spacing, so that their
# axes have the same height and baseline in the assembled row. The bar thickness is that of panels b-e
# (result/figure3_be_geometry.json); edge_mm is the gap between the outer bar and the end of the axis; both panels space
# their bars for panel g's five rows. f is 80.4 mm wide and g 99.6 mm, 3 mm apart: the g header needs 99 mm at 8 pt.
# header_baseline_mm is the baseline of the 10 pt panel letter set by assemble_figure3_vector.py (7.8 pt).
PANEL = {"width_mm": 80.4, "height_mm": 38.2, "ax_top_mm": 9.6, "ax_bottom_mm": 6.7, "right_mm": 1.0,
         "header_x_mm": 4.0, "header_baseline_mm": 7.8 / 72 * 25.4, "header": 8, "tick": 7, "axis": 7,
         "column": 6.5, "column_gap_mm": 1.6, "head_above_axes_mm": 0.9, "marker": 3.5, "marker_edge": 0.4,
         "bar_edge": 0.3, "edge_mm": 1.7, "rows": 5, "label_edge_mm": 2.0}
HEADER = ("Replication ranking", "pooled across independent projects")
INK = PAL["ink"]

MODALITY_MARKERS = {"mouse RNA-seq": "o", "human RNA-seq": "s", "human DNAm": "^"}
MODALITY_HATCH = {"mouse RNA-seq": "", "human RNA-seq": "////", "human DNAm": "...."}
COLOR_PRO = PAL["higher"]   # July: #D55E00
COLOR_ANTI = PAL["lower"]   # July: #0072B2
EDGE_PRO, EDGE_ANTI = PAL["higher_edge"], PAL["lower_edge"]  # bar outline and hatch strokes (July: #111)
# Display names of the canonical targets (the meta-analysis key is lower case).
DISPLAY = {"sleep deprivation": "Sleep Deprivation", "tlr9": "Tlr9 loss", "plx5622": "PLX5622",
           "erythropoietin": "Erythropoietin"}
TESTED = {"plx5622"}  # followed up in aged mice (Figure 7k-m)
TESTED_MARK = " †"    # the mark of panel a; defined in the key above panels b-e
RANK_HEAD = "Best single rank"
NOT_RANKED = "not ranked"
KEYS = ["gse_id", "comparison_name", "condition_name"]


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
    """Width (mm) of each value column: its widest head line or value."""
    kw = dict(fontsize=PANEL["column"])
    return [max(width_mm(fig, s, **kw) for s in list(heads) + list(values) if s) for heads, values, _ in columns]


def draw_columns(fig, ax, columns, widths, y_values):
    """Columns right of the axes, ending right_mm from the panel's right edge, one value per bar row; each column is
    centred or right-aligned (heads and values alike). Heads are bottom-aligned: their last line sits
    head_above_axes_mm above the axes top."""
    W, H = PANEL["width_mm"], PANEL["height_mm"]
    trans = blended_transform_factory(fig.transFigure, ax.transData)
    lead = PANEL["column"] * 1.2 / 72 * 25.4  # line spacing of the heads
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
    return artists, x_right + PANEL["column_gap_mm"]  # left edge of the first column


def place_axis(fig, ax, right_mm):
    """Fix the axes in millimetres: it begins after the widest bar label and ends a column gap before the columns."""
    W, H = PANEL["width_mm"], PANEL["height_mm"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    label = max(t.get_window_extent(renderer).width for t in ax.get_yticklabels()) / (fig.dpi * MM)
    left = PANEL["label_edge_mm"] + label + (plt.rcParams["ytick.major.size"] + plt.rcParams["ytick.major.pad"]) / 72 * 25.4
    ax.set_position([left / W, PANEL["ax_bottom_mm"] / H, (right_mm - PANEL["column_gap_mm"] - left) / W,
                     (H - PANEL["ax_top_mm"] - PANEL["ax_bottom_mm"]) / H])


def row_geometry(n, bar_mm):
    """Bar height (data units) and y limits that draw n bars at panel b-e thickness on panel g's row spacing."""
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


def statistical_ranks(sub, judged_path, ranked_path, meta_input_path):
    """Statistical rank of each replicated intervention: its best-ranked lower-output comparison (module docstring)."""
    judged = pd.read_csv(judged_path, low_memory=False)
    assert not judged.duplicated(KEYS).any()
    lower = judged[judged["hedges_g"] < 0].copy()
    # statistical_rank orders the lower-output comparisons by the absolute CI bound nearest zero, from 1
    order = lower.iloc[np.argsort(-lower["rank_bound"].abs().to_numpy(), kind="stable")]
    assert order["statistical_rank"].tolist() == list(range(1, len(lower) + 1))
    shown = pd.read_csv(ranked_path, low_memory=False, usecols=KEYS + ["meta_target"])
    linked = lower.merge(shown, on=KEYS, how="left", validate="one_to_one")
    both = linked["canonical_target"].notna() & linked["meta_target"].notna()
    assert (linked.loc[both, "canonical_target"] == linked.loc[both, "meta_target"]).all()
    inputs = pd.read_csv(meta_input_path, low_memory=False, usecols=KEYS + ["canonical_target"])
    rows = []
    for _, r in sub.iterrows():
        hits = linked[linked["canonical_target"] == r["canonical_target"]].sort_values("statistical_rank")
        # the same comparisons are found through meta_target
        assert hits[KEYS].values.tolist() == \
            linked[linked["meta_target"] == r["canonical_target"]].sort_values("statistical_rank")[KEYS].values.tolist()
        out = {"canonical_target": r["canonical_target"], "n_lower_ranked": len(hits),
               "n_ranked_any_direction": int((judged["canonical_target"] == r["canonical_target"]).sum())}
        if len(hits):
            best = hits.iloc[0]
            pooled = inputs[(inputs["canonical_target"] == r["canonical_target"])]
            assert (pooled[KEYS] == best[KEYS].values).all(axis=1).any(), "rank comparison is not a pooled effect"
            assert best["gse_id"] in str(r["series"]).split(";")
            out.update({"statistical_rank": int(best["statistical_rank"]), "rank_text": str(int(best["statistical_rank"])),
                        "rank_gse_id": best["gse_id"], "rank_comparison_name": best["comparison_name"],
                        "rank_condition_name": best["condition_name"], "rank_label": best["label"],
                        "rank_hedges_g": best["hedges_g"], "rank_rank_bound": best["rank_bound"]})
        else:
            out.update({"statistical_rank": np.nan, "rank_text": NOT_RANKED})
        rows.append(out)
    return pd.DataFrame(rows), len(lower)


def main():
    verifier = Path(os.environ["CLOCKBASE_VERIFIER_ROOT"])
    source = Path(os.environ.get('CLOCKBASE_META') or verifier / 'intervention-meta-2026-10-17') / 'primary/ranked_replicated_age-lowering.csv'
    judged_path = Path(os.environ.get('CLOCKBASE_JUDGEMENT') or verifier / 'followup-judgment-001') / 'ranked_with_judgement.csv'
    meta_input_path = Path(os.environ.get('CLOCKBASE_META') or verifier / 'intervention-meta-2026-10-17') / 'primary/meta_input_effects.csv'
    ranked_path = Path("result/volcano_ranked_comparisons_july.csv")
    geometry_path = Path("result/figure3_be_geometry.json")
    require(source, judged_path, meta_input_path, ranked_path, geometry_path)
    bar_mm = json.loads(geometry_path.read_text())["bar_mm"]
    sub = pd.read_csv(source)
    assert (sub["estimate"] < 0).all() and (sub["ci_high"] < 0).all() and (sub["fdr"] < 0.05).all()
    assert (sub["k"] >= 2).all() and set(sub["modalities"]) == {"mouse RNA-seq"}
    assert np.allclose(sub["rank_bound"].abs(), -sub["ci_high"])
    # the printed project and series counts are the counts of the projects and series the table lists
    assert (sub["k"] == sub["studies"].str.split(";").str.len()).all()
    assert (sub["n_series"] == sub["series"].str.split(";").str.len()).all()
    sub["composite"] = sub["rank_bound"].abs()
    sub["log2FoldChange"] = sub["estimate"]
    sub["modality"] = sub["modalities"]
    sub["display"] = sub["canonical_target"].map(DISPLAY)
    assert sub["display"].notna().all(), "a new replicated intervention needs a display name"
    sub["display"] = [d + TESTED_MARK if t in TESTED else d for d, t in zip(sub["display"], sub["canonical_target"])]
    assert sub["display"].str.endswith(TESTED_MARK).sum() == len(TESTED) == 1
    sub = sub.sort_values("composite", ascending=False, kind="stable")
    ranks, n_lower = statistical_ranks(sub, judged_path, ranked_path, meta_input_path)
    sub = sub.merge(ranks, on="canonical_target", how="left", validate="one_to_one")
    sub["rank_of"] = n_lower
    sub["rank_source_file"] = str(judged_path)
    sub.to_csv("result/figure3_replication_picks.csv", index=False)
    fig, ax = plt.subplots(figsize=(PANEL["width_mm"] * MM, PANEL["height_mm"] * MM))
    sub = sub.sort_values("composite", ascending=True)
    labels = [f"{r['display']}" for _, r in sub.iterrows()]
    values = sub["composite"].values
    colors = [COLOR_PRO if v > 0 else COLOR_ANTI for v in sub["log2FoldChange"].values]
    edges = [EDGE_PRO if v > 0 else EDGE_ANTI for v in sub["log2FoldChange"].values]
    hatches = [MODALITY_HATCH[m] for m in sub["modality"].values]

    y_pos = np.arange(len(sub))
    assert len(sub) <= PANEL["rows"]
    bar_height, ylim = row_geometry(len(sub), bar_mm)
    bars = ax.barh(y_pos, values, height=bar_height, color=colors, edgecolor=edges, linewidth=PANEL["bar_edge"])
    ax.set_ylim(*ylim)
    for bar, h in zip(bars, hatches):
        bar.set_hatch(h * 2)  # hatch spacing is fixed in inches; doubled to keep the July density at print size

    # Modality marker at end of each bar
    x_pad = (values.max() - values.min()) * 0.02 if values.max() > values.min() else 0.5
    for y, (_, r) in zip(y_pos, sub.iterrows()):
        ax.plot(r["composite"] + x_pad, y,
                marker=MODALITY_MARKERS[r["modality"]],
                markersize=PANEL["marker"], color=PAL["marker_edge"], markerfacecolor=PAL["marker_face"],
                markeredgewidth=PANEL["marker_edge"], zorder=5)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=PANEL["tick"])
    for tick, (_, r) in zip(ax.get_yticklabels(), sub.iterrows()):
        if r["canonical_target"] in TESTED:
            tick.set_fontweight("bold")
    ax.set_xlabel("Pooled |Hedges' g| (CI bound nearest zero)", fontsize=PANEL["axis"], labelpad=1.5)
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.grid(axis="x", alpha=0.3, linestyle=":")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="y", labelsize=PANEL["tick"])
    ax.tick_params(axis="x", labelsize=PANEL["tick"])

    # Projects, series and statistical rank in columns right of the bars; the x axis holds the bars only.
    columns = [(["Projects"], [f"{int(k)}" for k in sub["k"]], "center"),
               (["Series"], [f"{int(n)}" for n in sub["n_series"]], "center"),
               ([RANK_HEAD, f"(of {n_lower} decreases)"], list(sub["rank_text"]), "right")]
    widths = column_widths(fig, columns)
    _, first_column = draw_columns(fig, ax, columns, widths, y_pos)
    place_axis(fig, ax, first_column)
    header = ranking_header(fig, *HEADER)
    fig.canvas.draw()
    lo, hi = ax.get_xlim()
    # tight to the data: the axis ends at the July autoscale margin past the longest bar's end marker
    assert lo == 0 and values.max() + x_pad < hi < (values.max() + x_pad) * 1.1, (lo, hi)
    # header above the column heads, heads above the axes
    heads_top = min(box_mm(fig, t)[1] for t in fig.texts if t.get_text() in ("Projects", "Series", RANK_HEAD))
    assert heads_top >= max(box_mm(fig, t)[3] for t in header) + 0.5
    # the bar labels start at least label_edge_mm in from the left page edge
    assert min(box_mm(fig, t)[0] for t in ax.get_yticklabels()) >= PANEL["label_edge_mm"] - 0.05
    check_text_layout(fig, ax)
    Path("figures").mkdir(exist_ok=True)
    fig.savefig("figures/replication_barplot_combined_modalities.pdf", dpi=300)
    fig.savefig("figures/replication_barplot_combined_modalities.png", dpi=300)
    print(sub[["display", "k", "n_series", "rank_text", "rank_gse_id", "rank_comparison_name"]].iloc[::-1].to_string())


if __name__ == "__main__":
    main()
