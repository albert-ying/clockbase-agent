#!/usr/bin/env python3
"""
Figure 3b-e — the statistical ranking: per category, the five top-ranked interventions with higher and the five
with lower clock outputs, by the absolute CI bound of Hedges' g nearest zero (one entry per intervention and GEO
series, dedup_interventions_2026_11_06.py). Each bar with lower outputs carries three cells, one per design condition
of the judgement version, from the answers agreed by both judging models
(followup-judgment-001/ranked_with_judgement.csv): old at the start of the intervention, otherwise healthy (the
background answer: HEALTHY when the animals or people carry no disease apart from the compared exposure, so a disease
compared with healthy controls is HEALTHY and a treatment given to a disease model is DISEASE_MODEL), untreated
control. The three heads repeat the words of the panel g header. A filled cell marks a condition met, a crossed white
cell a condition not met and a light grey cell an answer not stated in the source or not agreed by the two models.
The unmet conditions are written left of the cells; the one comparison meeting all three reads "all three met". The
July bar drawing code and modality encodings are retained; the cells are drawn after it.

Labels: the July label (intervention name) of a bar whose canonical_target is that of one of the three comparisons
followed up in aged mice (panel a, dagger) starts with that comparison's label, so the two can be read together
("4°C", GSE201135, becomes "Cold exposure"; the followed-up comparison is GSE165974). The second key line ends
with the dagger, which marks the followed-up interventions in panels a, f and g.

Outputs:
  figures/top_hits_barplot_combined_modalities_v2.{pdf,png}
  result/top_hits_combined_modalities_v2.csv
  result/figure3_condition_grid.csv   one row per lower-output bar: the three states, the reason text and the answers
  result/figure3_be_geometry.json     bar thickness and cell column, read by the panel f and g generators
"""
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerBase
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import to_rgba
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
# Print size (7 November 2026): the grid is drawn at 183 x 106.7 mm and placed at 100% in Figure 3, Arial throughout.
# The 2.7 mm added to the print-size grid of 7 November (104 mm) holds the second key line; the bars keep their 2.78 mm.
plt.rcParams.update({
    "font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.0, "ytick.major.size": 2.0, "xtick.major.pad": 1.5, "ytick.major.pad": 1.5,
    "grid.linewidth": 0.4, "hatch.linewidth": 0.4,
})
MM = 1 / 25.4  # inches per millimetre
PRINT = {"width_mm": 183.0, "height_mm": 106.7, "tick": 7, "title": 8, "axis": 7, "reason": 6, "key": 6.5,
         "column_head": 6.5, "header": 8, "marker": 3.5, "marker_edge": 0.4, "bar_edge": 0.3}

CATEGORIES = ["Genetic", "Drug", "Disease", "Environment"]
TOP_N = 5  # per direction
# Cap on -log10(p) when forming the composite. Extreme-p outliers (e.g. sorted
# immune-cell comparisons with p < 1e-100) are retained but their bar length is
# bounded so they do not swamp the panel. The 95th percentile of -log10(p) across
# qualifying hits is ~20, so a cap of 30 leaves nearly all hits untouched.
CAP_NEGLOG10P = 30
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
INK = PAL["ink"]
REASON_GREY = PAL["reason_text"]
MET_FILL, MET_EDGE = PAL["met_fill"], PAL["met_edge"]          # condition met (July: ClockBase slate)
UNKNOWN_FILL, UNKNOWN_EDGE = PAL["unknown_fill"], PAL["unknown_edge"]

# The ranking header that replaces the July suptitle: the name in bold, the rule in regular weight, one noun phrase.
HEADER = ("Statistical ranking", "of single comparisons with the {n} strongest increases and decreases per category")
HEADER_X_MM = 4.0         # left edge shared with the panel f and g headers
# Design conditions of the judgement version (Methods) and the state of each agreed answer. NaN means the two judging
# models did not agree; it is shown like an answer that is not stated.
FIELDS = ("age_at_start", "background", "control")
COLUMN_HEADS = ("Old at start", "Otherwise healthy", "Untreated control")
STATE_COLUMNS = ("state_old_at_start", "state_otherwise_healthy", "state_untreated_control")
STATE = {
    "age_at_start": {"OLD": "met", "ADULT": "not met", "YOUNG": "not met", "GERMLINE_OR_PRENATAL": "not met",
                     "NOT_STATED": "unknown"},
    "background": {"HEALTHY": "met", "DISEASE_MODEL": "not met", "NOT_STATED": "unknown"},
    "control": {"UNTREATED": "met", "ACTIVE_COMPARATOR": "not met"},
}
REASON = {"GERMLINE_OR_PRENATAL": "germline", "ADULT": "adult", "YOUNG": "young", "DISEASE_MODEL": "disease model",
          "ACTIVE_COMPARATOR": "active control"}
ALL_MET = "all three met"
CELL_KEY = (("met", "condition met"), ("not met", "not met"), ("unknown", "not stated or not agreed"))
DAGGER_KEY = "† followed up in aged mice"   # the mark of panels a, f and g, as a separate group after the cell key
DAGGER_GROUP_GAP_MM = 5.0                    # between the cell key and the dagger entry
# The three comparisons followed up in aged mice (panel a; create_combined_modality_volcano.py TESTED_KEYS).
TESTED_KEYS = [("GSE122080", "Ouabain vs Saline in old mice"),
               ("GSE134151", "5xFAD Hippocampus PLX5622 vs Vehicle"),
               ("GSE165974", "Cold exposure vs Room temperature in DPP4_pos APCs")]
# Cell column geometry (mm): cell side = bar thickness
CELL_GAP_MM = 0.55        # between neighbouring cells
AXIS_TO_CELLS_MM = 1.6    # from the end of the x axis to the first cell
REASON_PAD_MM = 0.7       # from the reason text to the first cell
MARKER_CLEAR_MM = 2.0     # least distance from a bar's end marker to its reason text
HEAD_ABOVE_CELL_MM = 0.5  # from the top cell of a subplot to the rotated column heads
KEY_TOP_MM, KEY_GAP_MM, KEY_MARGIN_MM = 0.4, 0.6, 3.0
HEADER_GAP_MM = 1.2       # from the second key line to the top of the header
LABEL_EDGE_MM = 2.0       # least distance from the left page edge to the bar labels of b and d
COLUMN_GAP_MM = 2.3       # from the cells of b and d to the middle of the page (the c and e labels begin 0.3 mm past it)


def qualify_with_tested(picks, judged):
    """Prefix the label of a bar that shares its canonical_target with a followed-up comparison by that comparison's
    label: "4°C" (GSE201135, cold exposure) reads "Cold exposure". Idempotent; returns the labels and the
    rows changed."""
    keys = ["gse_id", "comparison_name", "condition_name"]
    tested = pd.DataFrame(TESTED_KEYS, columns=["gse_id", "comparison_name"]).merge(
        judged[["gse_id", "comparison_name", "canonical_target", "label"]], on=["gse_id", "comparison_name"],
        how="left", validate="one_to_one")
    assert tested["canonical_target"].notna().all() and tested["canonical_target"].is_unique
    target = picks.merge(judged[keys + ["canonical_target"]], on=keys, how="left", validate="one_to_one")
    target.index = picks.index
    display, changed = picks["display"].copy(), []
    for t in tested.itertuples():
        hit = target["canonical_target"].eq(t.canonical_target) & ~(
            target["gse_id"].eq(t.gse_id) & target["comparison_name"].eq(t.comparison_name))
        for i in picks.index[hit]:
            display[i] = t.label
            changed.append(i)
    return display, changed


def condition_states(row):
    """States of the three conditions for one comparison; an unexpected answer stops the build."""
    return ["unknown" if pd.isna(row[f]) else STATE[f][row[f]] for f in FIELDS]


def reason_text(row):
    """The unmet conditions of one comparison in the order of the cells, or 'all three met'."""
    states = condition_states(row)
    if states == ["met"] * 3:
        return ALL_MET
    return " · ".join(REASON[row[f]] for f, s in zip(FIELDS, states) if s == "not met")


def load_modality_scored():
    """Scored rows for all three modalities, including the extreme-effect-size /
    p-value flagged outliers, read from the same combined table that feeds
    Tables 1-2 (figures/combined_3modality_plot_data.csv) so the bar plot and the
    tables draw from the same data. Rows with a null outlier_status are unscored
    and dropped."""
    df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)
    df = df[df["outlier_status"].notna()].copy()
    keep = ["modality", "gse_id", "condition_category", "condition_name",
            "intervention", "log2FoldChange", "fdr", "final_score", "direction",
            "pvalue", "outlier_status"]
    return df[[c for c in keep if c in df.columns]]


class ModalityKey:
    """Key entry for one modality: the bar hatching of that modality and the marker at the bar end."""
    def __init__(self, modality):
        self.modality = modality


class HandlerModalityKey(HandlerBase):
    """Draw a short blue bar swatch with the modality's hatching, its end marker just past it, as on the bars.
    Swatch and gap are in font-size units; the entry box is widened to hold both."""
    SWATCH, GAP = 1.1, 0.2  # swatch width and swatch-to-marker gap, in font sizes

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        radius = (PRINT["marker"] + 0.5) / 2
        handlebox.width = (self.SWATCH + self.GAP) * fontsize + 2 * radius
        return super().legend_artist(legend, orig_handle, fontsize, handlebox)

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        m = orig_handle.modality
        radius = (PRINT["marker"] + 0.5) / 2
        swatch = Rectangle((-xdescent, -ydescent), self.SWATCH * fontsize, height, facecolor=COLOR_ANTI,
                           edgecolor=EDGE_ANTI, linewidth=PRINT["bar_edge"], hatch=MODALITY_HATCH[m] * 2,
                           transform=trans)
        x = -xdescent + (self.SWATCH + self.GAP) * fontsize + radius
        marker = Line2D([x], [-ydescent + height / 2], marker=MODALITY_MARKERS[m], color=PAL["key_marker_edge"],
                        markerfacecolor=PAL["key_marker_face"], markersize=PRINT["marker"] + 0.5,
                        markeredgewidth=PRINT["marker_edge"], linestyle="None", transform=trans)
        return [swatch, marker]


def cell_artists(x0, y0, w, h, state, transform):
    """One condition cell: met = slate fill, not met = white with a cross, unknown = light grey."""
    if state == "met":
        arts = [Rectangle((x0, y0), w, h, facecolor=MET_FILL, edgecolor=MET_EDGE, linewidth=0.3)]
    elif state == "not met":
        ix, iy = 0.2 * w, 0.2 * h
        arts = [Rectangle((x0, y0), w, h, facecolor=PAL["notmet_fill"], edgecolor=INK, linewidth=0.4),
                Line2D([x0 + ix, x0 + w - ix], [y0 + iy, y0 + h - iy], color=INK, linewidth=0.6,
                       solid_capstyle="butt"),
                Line2D([x0 + ix, x0 + w - ix], [y0 + h - iy, y0 + iy], color=INK, linewidth=0.6,
                       solid_capstyle="butt")]
    elif state == "unknown":
        arts = [Rectangle((x0, y0), w, h, facecolor=UNKNOWN_FILL, edgecolor=UNKNOWN_EDGE, linewidth=0.3)]
    else:
        raise ValueError(state)
    for a in arts:
        a.set_transform(transform)
        a.set_clip_on(False)
    return arts


class CellKey:
    """Key entry for one cell state."""
    def __init__(self, state):
        self.state = state


class HandlerCellKey(HandlerBase):
    """A square cell of the key height, drawn as on the bars."""
    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        handlebox.width = handlebox.height
        return super().legend_artist(legend, orig_handle, fontsize, handlebox)

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        return cell_artists(-xdescent, -ydescent, height, height, orig_handle.state, trans)


def print_size(axes):
    """Type, marker and stroke sizes of the July grid at the 183 mm print width (the July drawing block is unchanged)."""
    for ax in axes.flat:
        ax.tick_params(axis="both", labelsize=PRINT["tick"])
        for bar in ax.patches:
            bar.set_linewidth(PRINT["bar_edge"])
            # outline and hatch strokes in the palette's edge colour of the bar's direction (July: #111 for both)
            bar.set_edgecolor(EDGE_ANTI if np.allclose(bar.get_facecolor(), to_rgba(COLOR_ANTI)) else EDGE_PRO)
            # hatch spacing is fixed in inches; doubling the pattern keeps the July density at about half the size
            if bar.get_hatch():
                bar.set_hatch(bar.get_hatch() * 2)
        for line in ax.lines:
            line.set_markersize(PRINT["marker"])
            line.set_markeredgewidth(PRINT["marker_edge"])


def box_mm(fig, artist):
    """Window extent of an artist in mm from the top-left page corner: (left, top, right, bottom)."""
    fig.canvas.draw()
    bb = artist.get_window_extent(fig.canvas.get_renderer())
    px, H = fig.dpi * MM, PRINT["height_mm"]
    return np.array([bb.x0 / px, H - bb.y1 / px, bb.x1 / px, H - bb.y0 / px])


def text_width_mm(fig, artist):
    left, _, right, _ = box_mm(fig, artist)
    return right - left


def key_and_header(fig):
    """Two key lines centred at the top (bar encodings, then the cell states) and the ranking header below them."""
    W, H = PRINT["width_mm"], PRINT["height_mm"]
    direction_handles = [
        Patch(facecolor=COLOR_ANTI, edgecolor=EDGE_ANTI, linewidth=PRINT["bar_edge"], label="lower clock output"),
        Patch(facecolor=COLOR_PRO,  edgecolor=EDGE_PRO, linewidth=PRINT["bar_edge"], label="higher clock output"),
    ]
    # Modality is shown on the bars by both hatching and end marker, so each key entry shows both.
    modality_handles = [ModalityKey(m) for m in MODALITY_MARKERS]
    style = dict(loc="upper center", frameon=False, fontsize=PRINT["key"], borderaxespad=0, borderpad=0,
                 handlelength=1.6, handleheight=0.8, handletextpad=0.5, columnspacing=1.8)
    bars_key = fig.legend(handles=direction_handles + modality_handles,
                          labels=[h.get_label() for h in direction_handles] + list(MODALITY_MARKERS),
                          handler_map={ModalityKey: HandlerModalityKey()}, ncol=5,
                          bbox_to_anchor=(0.5, 1 - KEY_TOP_MM / H), **style)
    first = box_mm(fig, bars_key)
    cells_key = fig.legend(handles=[CellKey(s) for s, _ in CELL_KEY], labels=[t for _, t in CELL_KEY],
                           handler_map={CellKey: HandlerCellKey()}, ncol=3,
                           bbox_to_anchor=(0.5, 1 - (first[3] + KEY_GAP_MM) / H), **style)
    # the dagger entry is a separate group after the cell key, on the same line; the line is centred as a whole
    cells_box = box_mm(fig, cells_key)
    dagger = fig.text(0, 0, DAGGER_KEY, fontsize=PRINT["key"], color="black", ha="left", va="baseline")
    dagger_w = text_width_mm(fig, dagger)
    line_w = cells_box[2] - cells_box[0] + DAGGER_GROUP_GAP_MM + dagger_w
    shift = (W - line_w) / 2 - cells_box[0]
    cells_key.set_bbox_to_anchor((0.5 + shift / W, 1 - (first[3] + KEY_GAP_MM) / H), transform=fig.transFigure)
    second = box_mm(fig, cells_key)
    key_text = cells_key.get_texts()[0]
    fig.canvas.draw()
    baseline_px = key_text.get_transform().transform(key_text.get_position())[1]
    dagger.set_position(((second[2] + DAGGER_GROUP_GAP_MM) / W, baseline_px / fig.bbox.height))
    d_box = box_mm(fig, dagger)
    assert abs(d_box[3] - box_mm(fig, key_text)[3]) < 0.05, "dagger entry is not level with the cell key"
    second = np.array([second[0], min(second[1], d_box[1]), d_box[2], max(second[3], d_box[3])])
    for b in (first, second):
        assert b[0] >= KEY_MARGIN_MM and b[2] <= W - KEY_MARGIN_MM, f"key line leaves less than 3 mm: {b}"
    assert abs((second[0] + second[2]) / 2 - W / 2) < 0.05, "second key line is not centred"
    assert second[1] - first[3] >= KEY_GAP_MM - 0.05
    bold_text, rest_text = HEADER[0], HEADER[1].format(n=TOP_N)
    probe = fig.text(0, 0, "Hg", fontsize=PRINT["header"], fontweight="bold", va="baseline")
    _, top, _, _ = box_mm(fig, probe)
    ascent = H - top  # mm from the baseline to the top of the line box
    probe.remove()
    baseline = second[3] + HEADER_GAP_MM + ascent
    y = 1 - baseline / H
    bold = fig.text(HEADER_X_MM / W, y, bold_text, fontsize=PRINT["header"], fontweight="bold", color=INK,
                    ha="left", va="baseline")
    space = fig.text(0, 0, "a a", fontsize=PRINT["header"])
    nospace = fig.text(0, 0, "aa", fontsize=PRINT["header"])
    gap = text_width_mm(fig, space) - text_width_mm(fig, nospace)
    space.remove()
    nospace.remove()
    rest = fig.text((HEADER_X_MM + text_width_mm(fig, bold) + gap) / W, y, rest_text, fontsize=PRINT["header"],
                    color=INK, ha="left", va="baseline")
    assert box_mm(fig, rest)[2] <= W - 1.0, "header runs past the page"
    return [bars_key, cells_key], [bold, rest], dagger


def place_axes(fig, axes, header_bottom_mm, reserve_mm):
    """Place the 2 x 2 grid in millimetres: each column's axes begin after that column's widest bar label and end
    reserve_mm before the column's right edge (room for the cells), and the two rows share the height left under
    the key and the header line. Returns the right edge of each column (where the cell column ends) and the axes
    height."""
    W, H = PRINT["width_mm"], PRINT["height_mm"]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    px = fig.dpi * MM  # pixels per millimetre
    label_w = [max(t.get_window_extent(renderer).width for ax in axes[:, c] for t in ax.get_yticklabels()) / px
               for c in range(2)]
    tick = (plt.rcParams["ytick.major.size"] + plt.rcParams["ytick.major.pad"]) / 72 * 25.4
    boxes = [(ax.get_tightbbox(renderer), ax.get_window_extent(renderer)) for ax in axes.flat]
    above = max(tight.y1 - box.y1 for tight, box in boxes) / px   # subplot title
    below = max(box.y0 - tight.y0 for tight, box in boxes) / px   # x ticks and axis title
    top = header_bottom_mm + 1.8 + above
    row_gap, bottom_margin = 1.8, 0.3
    height = (H - top - 2 * below - row_gap - above - bottom_margin) / 2
    tops = [top, top + height + below + row_gap + above]
    col_right = [(c + 1) * W / 2 - (COLUMN_GAP_MM if c == 0 else 1.6) for c in range(2)]
    for (r, c), ax in np.ndenumerate(axes):
        left = c * W / 2 + (LABEL_EDGE_MM if c == 0 else 0.3) + label_w[c] + tick
        ax.set_position([left / W, 1 - (tops[r] + height) / H, (col_right[c] - reserve_mm - left) / W, height / H])
    return col_right, height


def add_cells(fig, ax, rows, right_mm, cell_mm):
    """Cells of the lower-output rows (rows = [(y, row)]) in a column ending at right_mm, the three column heads
    rotated above the top cell, and the reason text right-aligned left of the cells. Returns the cell axes, the
    column heads and the reasons [(y, text artist)]."""
    W, H = PRINT["width_mm"], PRINT["height_mm"]
    pos = ax.get_position()
    col_w = 3 * cell_mm + 2 * CELL_GAP_MM
    left_mm = right_mm - col_w
    cax = fig.add_axes([left_mm / W, pos.y0, col_w / W, pos.height])
    cax.set_xlim(0, col_w)
    cax.set_ylim(*ax.get_ylim())
    cax.set_axis_off()
    pitch_mm = pos.height * H / (ax.get_ylim()[1] - ax.get_ylim()[0])
    h = cell_mm / pitch_mm  # cell height in row units
    for y, r in rows:
        for i, s in enumerate(condition_states(r)):
            for a in cell_artists(i * (cell_mm + CELL_GAP_MM), y - h / 2, cell_mm, h, s, cax.transData):
                cax.add_artist(a)
    y_head = max(y for y, _ in rows) + h / 2 + HEAD_ABOVE_CELL_MM / pitch_mm
    heads = [cax.text(i * (cell_mm + CELL_GAP_MM) + cell_mm / 2, y_head, name, rotation=90, rotation_mode="anchor",
                      ha="left", va="center", fontsize=PRINT["column_head"], color=INK, clip_on=False)
             for i, name in enumerate(COLUMN_HEADS)]
    trans = blended_transform_factory(fig.transFigure, ax.transData)
    reasons = []
    for y, r in rows:
        text = reason_text(r)
        if not text:
            continue  # only answers not stated or not agreed: the grey cells carry it
        passes = text == ALL_MET
        reasons.append((y, ax.text((left_mm - REASON_PAD_MM) / W, y, text, transform=trans, ha="right", va="center",
                                   fontsize=PRINT["reason"], color=PAL["met_text"] if passes else REASON_GREY,
                                   fontweight="bold" if passes else "normal", clip_on=False)))
    return cax, heads, reasons


def marker_clearance(fig, ax, reasons):
    """Distance (mm) from each reason text to the end marker of its own bar."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    gaps = []
    for y, t in reasons:
        marker = [ln for ln in ax.lines if np.isclose(ln.get_ydata()[0], y)]
        assert len(marker) == 1
        gaps.append((t.get_window_extent(renderer).x0 - marker[0].get_window_extent(renderer).x1) / (fig.dpi * MM))
    return gaps


def clear_markers(fig, ax, reasons):
    """Lengthen the x axis of a subplot, in steps of its own tick spacing, until every reason clears its bar's end
    marker by MARKER_CLEAR_MM; the bars keep their values and order, only the scale changes."""
    for _ in range(40):
        gaps = marker_clearance(fig, ax, reasons)
        if min(gaps, default=np.inf) >= MARKER_CLEAR_MM:
            return gaps
        lo, hi = ax.get_xlim()
        ax.set_xlim(lo, hi * 1.025)
    raise AssertionError(f"reason text cannot clear the end markers: {gaps}")


def check_text_layout(fig, obstacles):
    """No two visible texts of the grid overlap, no text runs off the page, no text is under 6 pt, and no text
    touches a cell or a bar end marker (obstacles: artists whose boxes are checked against every text)."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = list(fig.texts)
    for leg in fig.legends:
        texts += list(leg.get_texts())
    for ax in fig.axes:
        texts += list(ax.texts)
        if not ax.axison:
            continue  # the cell columns: no axis is drawn
        texts += [ax.title, ax.xaxis.label]
        for locs, labels, (lo, hi) in ((ax.get_xticks(), ax.get_xticklabels(), ax.get_xlim()),
                                       (ax.get_yticks(), ax.get_yticklabels(), ax.get_ylim())):
            texts += [t for v, t in zip(locs, labels) if lo - 1e-9 <= v <= hi + 1e-9]  # drawn ticks only
    texts = [t for t in texts if t.get_visible() and t.get_text().strip()]
    assert all(t.get_fontsize() >= 6 for t in texts), [t.get_text() for t in texts if t.get_fontsize() < 6]
    boxes = [t.get_window_extent(renderer) for t in texts]
    page = fig.bbox
    for t, b in zip(texts, boxes):
        assert b.x0 >= page.x0 and b.x1 <= page.x1 and b.y0 >= page.y0 and b.y1 <= page.y1, t.get_text()
    eps = 0.1 * fig.dpi * MM  # 0.1 mm
    for i in range(len(texts)):
        for j in range(i + 1, len(texts)):
            a, b = boxes[i], boxes[j]
            if a.x0 + eps < b.x1 and b.x0 + eps < a.x1 and a.y0 + eps < b.y1 and b.y0 + eps < a.y1:
                raise AssertionError(f"texts overlap: {texts[i].get_text()!r} and {texts[j].get_text()!r}")
    for o in obstacles:
        ob = o.get_window_extent(renderer)
        for t, b in zip(texts, boxes):
            if b.x0 + eps < ob.x1 and ob.x0 + eps < b.x1 and b.y0 + eps < ob.y1 and ob.y0 + eps < b.y1:
                raise AssertionError(f"text {t.get_text()!r} touches a cell or marker")
    return len(texts)


def make_label(row):
    iv = row.get("intervention")
    if isinstance(iv, str) and iv.strip():
        return iv.strip()[:45]
    return str(row.get("condition_name", ""))[:45]


def main():
    all_ = load_modality_scored()
    print(f"loaded {len(all_):,} scored rows across modalities")

    # Statistical category ranking supplied with verified confidence intervals.
    qual = pd.read_csv("result/volcano_ranked_comparisons_july.csv", low_memory=False)
    qual["composite"] = qual["rank_bound"].abs()
    picks = {}
    for cat in CATEGORIES:
        # one entry per intervention and per GEO series, both directions (dedup_interventions_2026_11_06.py)
        sub = qual[qual["display_best"] & (qual["display_category"] == cat)].copy()
        sub = sub.sort_values("composite", ascending=False, kind="stable")
        sub["display"] = sub.apply(make_label, axis=1)
        assert ((sub["log2FoldChange"] > 0) == (sub["hedges_g"] > 0)).all()
        sub = pd.concat([sub[sub["hedges_g"] > 0].head(TOP_N), sub[sub["hedges_g"] < 0].head(TOP_N)])
        assert (sub["hedges_g"] > 0).sum() == TOP_N and (sub["hedges_g"] < 0).sum() == TOP_N
        picks[cat] = sub.sort_values("composite", ascending=False, kind="stable")

    # Design answers of the judgement version for the bars with lower outputs.
    judged_path = Path(os.environ.get('CLOCKBASE_JUDGEMENT') or Path(os.environ["CLOCKBASE_VERIFIER_ROOT"]) / 'followup-judgment-001') / 'ranked_with_judgement.csv'
    if not judged_path.exists():
        raise SystemExit(f"Cannot proceed - required data file not found: {judged_path}")
    judged = pd.read_csv(judged_path, low_memory=False)
    keys = ["gse_id", "comparison_name", "condition_name"]
    assert not judged.duplicated(keys).any()
    # A bar sharing its canonical_target with a followed-up comparison is labelled so the two read together.
    qualified = []
    for cat in CATEGORIES:
        picks[cat]["display"], changed = qualify_with_tested(picks[cat], judged)
        qualified += [(cat, picks[cat].loc[i, "gse_id"], picks[cat].loc[i, "display"]) for i in changed]
    assert qualified == [("Environment", "GSE201135", "Cold exposure")], qualified
    for cat in CATEGORIES:
        merged = picks[cat].merge(judged[keys + ["age_at_start", "background", "control", "judgement_admitted"]],
                                  on=keys, how="left", validate="one_to_one")
        assert len(merged) == len(picks[cat]) and merged["judgement_admitted"].notna().all()
        merged.index = picks[cat].index
        lower = merged["hedges_g"] < 0
        states = [condition_states(r) if g else ["", "", ""] for (_, r), g in zip(merged.iterrows(), lower)]
        for i, col in enumerate(STATE_COLUMNS):
            merged[col] = [s[i] for s in states]
        merged["reason"] = [reason_text(r) if g else "" for (_, r), g in zip(merged.iterrows(), lower)]
        picks[cat] = merged
    all_picks = pd.concat([p.assign(category=c) for c, p in picks.items() if not p.empty], ignore_index=True)
    low = all_picks[all_picks["hedges_g"] < 0]
    assert len(low) == len(CATEGORIES) * TOP_N
    all_met = low[low["reason"] == ALL_MET]
    # A comparison is admitted to the judgement version exactly when it meets all three conditions; among the bars
    # with lower outputs only Mechanical Overload (GSE195707) does.
    assert ((low["reason"] == ALL_MET) == low["judgement_admitted"].astype(bool)).all()
    assert len(all_met) == 1 and all_met["gse_id"].iloc[0] == "GSE195707", all_met[keys]
    assert all_met["display"].iloc[0] == "Mechanical Overload"

    # Validate plotted values against the selected run, whose FDR family can change.
    source_keys = ["gse_id", "comparison_name", "condition_name", "clock", "modality"]
    source = qual.set_index(source_keys)
    selected = all_picks.set_index(source_keys)
    assert source.index.is_unique and selected.index.is_unique
    matched = source.reindex(selected.index)
    assert matched["fdr"].notna().all()
    for c in ("composite", "hedges_g", "g_ci_low", "g_ci_high", "log2FoldChange", "fdr"):
        assert np.allclose(matched[c].values, selected[c].values, rtol=0, atol=1e-12), c

    # Figure: 2×2 of categories
    fig, axes = plt.subplots(2, 2, figsize=(PRINT["width_mm"] * MM, PRINT["height_mm"] * MM))
    for ax, cat in zip(axes.flat, CATEGORIES):
        sub = picks[cat]
        if sub.empty:
            ax.set_title(f"{cat}  (no scored cases)", fontsize=PRINT["title"])
            ax.axis("off")
            continue
        sub = sub.sort_values("composite", ascending=True)
        labels = [f"{r['display']}" for _, r in sub.iterrows()]
        values = sub["composite"].values
        colors = [COLOR_PRO if v > 0 else COLOR_ANTI for v in sub["log2FoldChange"].values]
        hatches = [MODALITY_HATCH[m] for m in sub["modality"].values]

        y_pos = np.arange(len(sub))
        bars = ax.barh(y_pos, values, color=colors, edgecolor="#111", linewidth=0.5)
        for bar, h in zip(bars, hatches):
            bar.set_hatch(h)

        # Modality marker at end of each bar
        x_pad = (values.max() - values.min()) * 0.02 if values.max() > values.min() else 0.5
        for y, (_, r) in zip(y_pos, sub.iterrows()):
            ax.plot(r["composite"] + x_pad, y,
                    marker=MODALITY_MARKERS[r["modality"]],
                    markersize=7, color="#222", markerfacecolor="#f7f7f7",
                    markeredgewidth=0.8, zorder=5)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=12)
        ax.set_xlabel("|Hedges' g| (CI bound nearest zero)", fontsize=PRINT["axis"], labelpad=1.5)
        ax.set_title(cat, fontsize=PRINT["title"], fontweight="bold", pad=2.5)
        ax.grid(axis="x", alpha=0.3, linestyle=":")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="y", labelsize=11)
        ax.tick_params(axis="x", labelsize=11)

    # Judgement-version tags: the condition cells of the bars with lower outputs, drawn after the July block.
    print_size(axes)
    for ax in axes.flat:  # the narrower axes would otherwise take 2.5-unit steps (0.0, 2.5, ...)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6, steps=[1, 2, 5, 10], integer=True))
    keys_drawn, header, dagger_key = key_and_header(fig)
    header_bottom = max(box_mm(fig, t)[3] for t in header)
    # cell side = bar thickness: lay out once to find it, then reserve the cell column and lay out again
    _, height = place_axes(fig, axes, header_bottom, reserve_mm=0)
    ylims = {tuple(np.round(ax.get_ylim(), 9)) for ax in axes.flat}
    assert len(ylims) == 1, ylims  # same row pitch in all four subplots, so one cell size
    yl = axes[0, 0].get_ylim()
    cell_mm = 0.8 * height / (yl[1] - yl[0])  # the July bar height is 0.8 row units
    col_right, height2 = place_axes(fig, axes, header_bottom, reserve_mm=3 * cell_mm + 2 * CELL_GAP_MM + AXIS_TO_CELLS_MM)
    assert np.isclose(height, height2)
    obstacles, clearances, head_boxes = [], {}, {}
    for (r, c), ax in np.ndenumerate(axes):
        cat = CATEGORIES[2 * r + c]
        sub = picks[cat].sort_values("composite", ascending=True)
        rows = [(y, row) for y, (_, row) in zip(np.arange(len(sub)), sub.iterrows()) if row["hedges_g"] < 0]
        assert len(rows) == TOP_N
        cax, heads, reasons = add_cells(fig, ax, rows, col_right[c], cell_mm)
        clearances[cat] = clear_markers(fig, ax, reasons)
        obstacles += list(cax.artists) + list(ax.lines)
        head_boxes[cat] = [box_mm(fig, h).round(2).tolist() for h in heads]
        # the heads stay clear of the subplot title and below the header line (or the row above)
        title = box_mm(fig, ax.title)
        for hb in head_boxes[cat]:
            assert hb[0] > title[2] + 0.5 or hb[2] < title[0] - 0.5 or hb[3] < title[1] - 0.5 or hb[1] > title[3] + 0.5
            assert hb[1] >= header_bottom + 0.5
    n_text = check_text_layout(fig, obstacles)
    # the bar labels of b and d start at least LABEL_EDGE_MM in from the left page edge
    label_left = min(box_mm(fig, t)[0] for ax in axes[:, 0] for t in ax.get_yticklabels())
    assert label_left >= LABEL_EDGE_MM - 0.05, label_left
    print(f"{n_text} texts checked: no overlap, none under 6 pt; reason clearance from end markers (mm): "
          + "; ".join(f"{k} {min(v):.1f}" for k, v in clearances.items()))

    Path("figures").mkdir(exist_ok=True)
    Path("result").mkdir(exist_ok=True)
    pdf = "figures/top_hits_barplot_combined_modalities_v2.pdf"
    png = "figures/top_hits_barplot_combined_modalities_v2.png"
    fig.savefig(pdf, dpi=300)
    fig.savefig(png, dpi=300)

    # Save the picked table
    all_picks.to_csv("result/top_hits_combined_modalities_v2.csv", index=False)
    grid = all_picks[all_picks["hedges_g"] < 0].copy()
    grid["panel"] = grid["category"].map(dict(zip(CATEGORIES, "bcde")))
    grid = grid.rename(columns={"age_at_start": "source_age_at_start", "background": "source_background",
                                "control": "source_control"})
    for f in ("source_age_at_start", "source_background", "source_control"):
        grid[f] = grid[f].fillna("NaN (models did not agree)")
    grid["source_file"] = str(judged_path)
    grid[["panel", "category", "display", "gse_id", "comparison_name", "condition_name", "modality", "hedges_g",
          "composite", *STATE_COLUMNS, "reason",
          "source_age_at_start", "source_background", "source_control", "judgement_admitted", "source_file"]
         ].to_csv("result/figure3_condition_grid.csv", index=False)
    geometry = {"page_mm": [PRINT["width_mm"], PRINT["height_mm"]], "bar_mm": cell_mm, "cell_mm": cell_mm,
                "row_pitch_mm": cell_mm / 0.8, "axes_height_mm": height, "cell_gap_mm": CELL_GAP_MM,
                "cell_column_right_mm": col_right, "header_x_mm": HEADER_X_MM, "header_pt": PRINT["header"],
                "header": [t.get_text() for t in header], "column_heads": list(COLUMN_HEADS),
                "column_head_boxes_mm": head_boxes, "reason_marker_clearance_mm": clearances,
                "key_lines": [[t.get_text() for t in leg.get_texts()] for leg in keys_drawn[:1]]
                             + [[t.get_text() for t in keys_drawn[1].get_texts()] + [dagger_key.get_text()]],
                "label_edge_mm": LABEL_EDGE_MM,
                "xlim": {cat: list(ax.get_xlim()) for cat, ax in zip(CATEGORIES, axes.flat)}}
    Path("result/figure3_be_geometry.json").write_text(json.dumps(geometry, indent=2) + "\n")
    print(f"\n✓ {pdf}")
    print(f"✓ {png}")
    print(f"✓ result/top_hits_combined_modalities_v2.csv")
    print(f"✓ result/figure3_condition_grid.csv ({len(grid)} lower-output bars)")
    print(f"bar thickness {cell_mm:.3f} mm")

    # Quick modality-coverage summary
    print("\n--- modality coverage across picked top-hits ---")
    print(all_picks.groupby(["category", "modality"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    main()
