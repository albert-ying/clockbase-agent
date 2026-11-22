#!/usr/bin/env python3
"""Figure 4i,j: AgingAgentBench performance of ClockBase and of language models without the data, as bar plots.

Both panels score every system by the same three-way call (lower and significant, higher and significant, or not
significant). Drawn at print size and placed at 100 % in Figure 4 (i at x = 81 mm, j at x = 151.6 mm):
  i  benchmark_identification_bars.pdf  69.1 x 54 mm  three narrow columns, one per published call (Lower, Higher,
                                                      Not significant, with the call markers of panel h), on the
                                                      seven system rows of j. Each bar is the recall of that call:
                                                      the share of the published calls the system matched, with
                                                      hits/n printed beside it (ClockBase 7 pt bold, the others
                                                      6 pt). A 0.4 pt ClockBase line stands at ClockBase's value
                                                      from its bar down through the five model rows; the model
                                                      counts start just right of it. One axis title spans the
                                                      three columns. The Always 'higher' baseline sits below the
                                                      model rows after a gap of ROW_GAP rows, in i and in j.
  j  benchmark_accuracy_bars.pdf        31.4 x 54 mm  balanced accuracy (the mean of the three columns of i) with
                                                      the 95 % interval and the chance line at one third; no row
                                                      labels, rows asserted level with the rows of i
Sources:
  i  benchmark/published_clock_2026_10_05/groupbench/results_2026_10_10/three_way_by_class.json ('recall': [hits, n]
     per published call; 'published': n per call, printed in the column headers), written by
     groupbench/three_way_by_class_2026_11_07.py
  j  ed23/panel_c_d_models.csv (copied from the checked ED23 package by build_ed23_24_july_2026_11_06.py); each
     balanced accuracy is checked against the mean of the three recalls of the JSON (to 1e-9)
  h  figures/benchmark_schematic.pdf is read only to check that the column headers of i sit on h's header line.
Every bar, count, interval and printed number is read from these files; nothing is typed in except the expected
counts asserted below. The plotted values, the count placement and the row centres of both panels (points from the
page bottom) are written to result/figure4_benchmark_bars.csv.
"""
import ast
import itertools
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.text import Text  # noqa: E402
from matplotlib.transforms import Affine2D  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = ROOT / "paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact"
SOURCE = WORK.parent / "ed23" / "panel_c_d_models.csv"
THREE_WAY = Path(os.environ.get("CLOCKBASE_BENCH_RESULTS") or ROOT / "benchmark/published_clock_2026_10_05/groupbench/results_2026_10_10") / "three_way_by_class.json"
PANEL_H = WORK / "figures" / "benchmark_schematic.pdf"
FIGDIR = WORK / "figures"
RESULT = WORK / "result" / "figure4_benchmark_bars.csv"

# Style tokens shared with the schematic panel h.
INK = "#222222"
CLOCKBASE = "#33414E"
MODEL = "#A6A6A6"
RULE = "#D9D9D9"
RULE_EDGE = "#8C8C8C"
LOWER = "#0072B2"
HIGHER = "#D55E00"
NS = "#8C8C8C"
GRID = "#C8C8C8"
CHANCE = "#595959"

SYSTEMS = ["ClockBase", "GPT-6 Astra", "GPT-6 Sol", "GPT-6 Luna", "GPT-5.5", "GPT-4.1", "Always 'higher'"]
MODELS = SYSTEMS[1:6]
DISPLAY = {s: s for s in SYSTEMS}
DISPLAY["Always 'higher'"] = "Always ‘higher’"  # typographic quotes on the page; ASCII name in the sources
KEYS = {"ClockBase": "S3c", "GPT-6 Astra": "S4:openai:gpt-6-astra", "GPT-6 Sol": "S4:openai:gpt-6-sol",
        "GPT-6 Luna": "S4:openai:gpt-6-luna", "GPT-5.5": "S4:openai:gpt-5.5", "GPT-4.1": "S4:openai:gpt-4.1",
        "Always 'higher'": "S0"}
CLASSES = ("lower", "higher", "none")  # columns of i, left to right (the order of h's key)
CLASS_HEADER = {"lower": "Lower", "higher": "Higher", "none": "Not significant"}  # panel h's words
# Panel h's call markers (create_figure4_benchmark_schematic_2026_11_07.py, Canvas.marker): 3.4 pt, 0.6 pt edge.
CLASS_MARKER = {"lower": dict(marker="v", mfc=LOWER, mec=LOWER),
                "higher": dict(marker="^", mfc=HIGHER, mec=HIGHER),
                "none": dict(marker="o", mfc="white", mec=NS)}
MARKER_PT, MARKER_EW = 3.4, 0.6
EXPECTED_N = {"lower": 16, "higher": 36, "none": 58}
EXPECTED_HITS = {"ClockBase": (10, 30, 51), "GPT-6 Astra": (5, 16, 47), "GPT-6 Sol": (5, 15, 43),
                 "GPT-6 Luna": (5, 12, 42), "GPT-5.5": (6, 17, 38), "GPT-4.1": (6, 28, 23),
                 "Always 'higher'": (0, 36, 0)}
MM = 1 / 25.4

plt.rcParams.update({
    "font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white", "figure.facecolor": "white",
    "axes.linewidth": 0.5, "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2.0, "ytick.major.size": 0.0, "xtick.major.pad": 1.5, "ytick.major.pad": 2.0,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.labelsize": 7, "axes.labelpad": 1.5,
    "xtick.direction": "out",
})

# Geometry in mm. Both panels share the page height, the axes bottom and top and the y limits, so the rows of j are
# level with the rows of i. Page row: h 0-78 mm, i 81-150.1 mm, j 151.6-183 mm.
FIG_H = 54.0
BOTTOM_MM, TOP_STRIP_MM = 6.9, 7.3  # x axis 6.9 mm above the panel bottom; 7.3 mm strip above the axes
PLOT_H = FIG_H - BOTTOM_MM - TOP_STRIP_MM
W_I, W_J, GAP_IJ, X_I = 69.1, 31.4, 1.5, 81.0
LBL_W = 17.8      # row labels of i (the widest, Always 'higher', is 16.6 mm at 7 pt plus the 2 pt tick pad)
COL_W = 11.0      # one column: 0-100 %
COL_GAP = 6.2     # holds the counts that pass 100 and keeps one column's 100 clear of the next column's 0
TAIL_W = W_I - LBL_W - 3 * COL_W - 2 * COL_GAP  # right of the last column: its longest count
COL_X0 = [LBL_W + k * (COL_W + COL_GAP) for k in range(3)]
J_LEFT, J_RIGHT = 2.2, 5.8
HEADER_BASE_MM = 3.64  # baseline of h's 7 pt bold headers, mm from the top (checked against h's PDF)
N_LINE_PITCH_MM = 2.75  # header name to its n line
AXIS_TITLE_BASE_MM = 1.1  # baseline of the axis titles of i and j, mm from the bottom
ROW_GAP = 0.4  # extra row units between GPT-4.1 and the Always 'higher' baseline, in i and j alike
YLIM = (-0.55, len(SYSTEMS) - 0.55 + ROW_GAP)  # the bottom margin keeps the lowest count off the x axis
BAR_H = 0.6  # bar height, row units, as in j
TICK_FS, AXIS_FS, HEADER_FS, CB_COUNT_FS, COUNT_FS, NOTE_FS = 7, 7, 7, 7, 6, 6
COUNT_PAD_PT = 1.4  # bar end to count
COUNT_CLEAR_MM = 1.0  # a count ends at least this far before the next column's axis
EDGE_CLEAR_MM = 0.3  # and before the panel edge in the last column
TICK_GAP_MIN_MM = 3.0  # clear space between one column's '100' and the next column's '0'
# ClockBase line in each column of i, at ClockBase's value from its bar down through the five model rows.
#   "aligned" (used) the line unbroken, with each model count starting just right of ClockBase's value instead of
#             at its own bar end, so the counts of a column form one left-aligned column beside the line and no
#             count crosses it. In the Higher column GPT-4.1 (28/36) stops 0.6 mm short of ClockBase (30/36).
#   "behind"  the line drawn behind the counts at their bar ends (white knockout under each count). Every model
#             bar ends within one count width of ClockBase's value, so the line crosses all fifteen model counts
#             and the knockouts cut it into stubs inside the numbers.
#   "off"     not drawn.
# A fine dash (0, (1.6, 1.0)) was tried in place of the solid line; at print size it reads as a third grid line.
REF_LINE = "aligned"
REF_LW, REF_LINESTYLE = 0.4, "solid"
REF_CLEAR_PT = 0.8  # clear space between the line and the left edge of any count beside it
KNOCKOUT = dict(boxstyle="square,pad=0.05", fc="white", ec="none", lw=0)  # breaks the line under a count


def load():
    for path in (SOURCE, THREE_WAY, PANEL_H):
        if not path.exists():
            raise SystemExit(f"Cannot proceed - required data file not found: {path}")
    cd = pd.read_csv(SOURCE)
    assert list(cd["system"]) == SYSTEMS, f"unexpected row order: {list(cd['system'])}"
    # Only the balanced accuracy and its interval are used; the direction-only counts of ED23 are not plotted.
    cd = cd[["system", "balanced_accuracy", "interval"]].copy()
    iv = cd["interval"].map(ast.literal_eval)
    cd["interval_low"] = [float(a) for a, _ in iv]
    cd["interval_high"] = [float(b) for _, b in iv]
    assert ((cd["interval_low"] <= cd["balanced_accuracy"]) & (cd["balanced_accuracy"] <= cd["interval_high"])).all()

    tw = json.loads(THREE_WAY.read_text())
    assert tw["published"] == EXPECTED_N, f"published calls {tw['published']}"
    assert tw["records"] == sum(EXPECTED_N.values()), f"records {tw['records']}"
    header_n = {c: int(tw["published"][c]) for c in CLASSES}
    rows = []
    for s in SYSTEMS:
        rec = tw["systems"][KEYS[s]]
        recall = rec["recall"]
        for c in CLASSES:
            assert recall[c][1] == header_n[c], f"{s} {c}: n = {recall[c][1]}"
        hits = tuple(int(recall[c][0]) for c in CLASSES)
        assert hits == EXPECTED_HITS[s], f"{s}: recall hits {hits}, expected {EXPECTED_HITS[s]}"
        row = {"source_key": KEYS[s], "printed_name": DISPLAY[s]}
        for c in CLASSES:
            k, n = (int(v) for v in recall[c])
            row.update({f"{c}_hits": k, f"{c}_total": n, f"{c}_matched_pct": k / n * 100, f"{c}_label": f"{k}/{n}"})
        row["recall_mean"] = float(np.mean([recall[c][0] / recall[c][1] for c in CLASSES]))
        row["json_balanced_accuracy"] = float(rec["balanced_accuracy"])
        row["json_interval_low"], row["json_interval_high"] = (float(v) for v in rec["balanced_accuracy_interval"])
        rows.append(row)
    cd = pd.concat([cd, pd.DataFrame(rows)], axis=1)
    # j plots the balanced accuracy of ED23; it must be the mean of the three recalls that i plots.
    for r in cd.itertuples():
        assert abs(r.balanced_accuracy - r.recall_mean) <= 1e-9, (r.system, r.balanced_accuracy, r.recall_mean)
        assert abs(r.balanced_accuracy - r.json_balanced_accuracy) <= 1e-9, r.system
        assert abs(r.interval_low - r.json_interval_low) <= 1e-9 and abs(r.interval_high - r.json_interval_high) <= 1e-9, r.system
    cd["balanced_accuracy_label"] = [f"{v:.2f}" for v in cd["balanced_accuracy"]]
    cd["y"] = np.arange(len(cd))[::-1].astype(float)  # ClockBase on top
    cd.loc[cd["system"] != "Always 'higher'", "y"] += ROW_GAP  # the baseline row set apart from the model rows
    assert cd["y"].iloc[-1] == 0.0 and cd["system"].iloc[-1] == "Always 'higher'"
    cd = cd.set_index("system", drop=False)
    cd.index.name = None
    # The reference line of i stands at ClockBase's value: every model bar must stop short of it.
    for c in CLASSES:
        best = cd.loc[MODELS, f"{c}_hits"].max()
        assert cd.loc["ClockBase", f"{c}_hits"] > best, (c, cd.loc["ClockBase", f"{c}_hits"], best)
    return cd, header_n


def check_h_header_line():
    """The column headers of i sit on the header line of h; read h's PDF so a change there cannot go unnoticed."""
    import pymupdf
    P = 72 / 25.4
    doc = pymupdf.open(PANEL_H)
    base = [s["origin"][1] / P for b in doc[0].get_text("dict")["blocks"] for line in b.get("lines", [])
            for s in line["spans"] if s["text"].strip() in ("Ground truth", "Scoring")]
    doc.close()
    assert base and all(abs(b - HEADER_BASE_MM) < 0.05 for b in base), ("h header line moved", base)


# ------------------------------------------------------------------------------------------------- helpers
def mm_tr(fig):
    """Figure transform in mm from the bottom-left corner."""
    return Affine2D().scale(MM) + fig.dpi_scale_trans


def add_axes_mm(fig, fig_w, left_mm, axes_w_mm):
    ax = fig.add_axes([left_mm / fig_w, BOTTOM_MM / FIG_H, axes_w_mm / fig_w, PLOT_H / FIG_H])
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_ylim(*YLIM)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color=GRID, linestyle=(0, (1, 1.5)), linewidth=0.4)
    ax.grid(axis="y", visible=False)
    return ax


def bbox_mm(artist, fig):
    fig.canvas.draw()
    bb = artist.get_window_extent(fig.canvas.get_renderer())
    k = 25.4 / fig.dpi
    return bb.x0 * k, bb.y0 * k, bb.x1 * k, bb.y1 * k


def text_width_mm(fig, s, fs, weight):
    t = fig.text(0, 0, s, fontsize=fs, fontweight=weight)
    x0, _, x1, _ = bbox_mm(t, fig)
    t.remove()
    return x1 - x0


def axis_title_baseline(ax, text):
    """Axis title with its baseline AXIS_TITLE_BASE_MM above the panel bottom, centred under the axes."""
    ax.set_xlabel(text, fontsize=AXIS_FS)
    ax.xaxis.label.set_verticalalignment("baseline")
    ax.xaxis.set_label_coords(0.5, (AXIS_TITLE_BASE_MM - BOTTOM_MM) / PLOT_H)


def class_header(fig, centre_mm, c, n):
    """Panel h's call marker and word on h's header line, the n of that call on the line below, centred."""
    tr = mm_tr(fig)
    y = FIG_H - HEADER_BASE_MM
    t = fig.text(0, y, CLASS_HEADER[c], fontsize=HEADER_FS, ha="left", va="baseline", transform=tr)
    w = text_width_mm(fig, CLASS_HEADER[c], HEADER_FS, "normal")
    m, gap = MARKER_PT * 25.4 / 72, 0.7
    left = centre_mm - (m + gap + w) / 2
    t.set_x(left + m + gap)
    cap = HEADER_FS * 0.716 * 25.4 / 72  # Arial cap height
    fig.add_artist(Line2D([left + m / 2], [y + cap / 2], transform=tr, linestyle="none", ms=MARKER_PT,
                          mew=MARKER_EW, **CLASS_MARKER[c]))
    fig.text(centre_mm, FIG_H - (HEADER_BASE_MM + N_LINE_PITCH_MM), f"$n$ = {n}", fontsize=HEADER_FS,
             ha="center", va="baseline", transform=tr)


def row_colours(cd):
    colours = [CLOCKBASE if s == "ClockBase" else RULE if s.startswith("Always") else MODEL for s in cd["system"]]
    edges = [RULE_EDGE if s.startswith("Always") else INK for s in cd["system"]]
    return colours, edges


def check_layout(fig, name):
    """No text off the canvas or under 6 pt; no text touching other text; no text on a bar (a count set inside its
    own bar end carries gid 'inside')."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    W, H = fig.bbox.width, fig.bbox.height
    texts = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        bb = t.get_window_extent(r)
        if bb.width == 0:
            continue
        texts.append((t, bb))
    problems = []
    for t, bb in texts:
        if t.get_fontsize() < 6:
            problems.append(f"{name}: {t.get_text()!r} at {t.get_fontsize()} pt")
        if bb.x0 < -0.5 or bb.y0 < -0.5 or bb.x1 > W + 0.5 or bb.y1 > H + 0.5:
            problems.append(f"{name}: text off canvas {t.get_text()!r} {bb.bounds}")
    for (t1, b1), (t2, b2) in itertools.combinations(texts, 2):
        ix = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
        iy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
        if ix > 0.5 and iy > 0.5:
            problems.append(f"{name}: text overlap {t1.get_text()!r} / {t2.get_text()!r}")
    bars = [p.get_window_extent(r) for ax in fig.axes for p in ax.patches
            if isinstance(p, Rectangle) and p.get_width() != 0 and p.get_height() != 0]
    for t, bb in texts:
        if t.get_gid() == "inside":
            continue
        shrunk = bb.padded(-0.8)  # glyph boxes carry side bearings; ignore sub-pixel touches
        if any(shrunk.overlaps(b) for b in bars):
            problems.append(f"{name}: text {t.get_text()!r} overlaps a bar")
    return problems


# ------------------------------------------------------------------------------------------------- panel i
def panel_identification(cd, header_n, ref_line=REF_LINE, ref_linestyle=REF_LINESTYLE):
    assert TAIL_W > 5.0, TAIL_W
    assert ref_line in ("off", "behind", "aligned"), ref_line
    fig = plt.figure(figsize=(W_I * MM, FIG_H * MM))
    colours, edges = row_colours(cd)
    mm_per_pct = COL_W / 100
    pad_mm = COUNT_PAD_PT * 25.4 / 72
    axes, counts, lines = [], {}, {}
    cb = cd.loc["ClockBase"]
    for k, c in enumerate(CLASSES):
        ax = add_axes_mm(fig, W_I, COL_X0[k], COL_W)
        ax.barh(cd["y"], cd[f"{c}_matched_pct"], height=BAR_H, color=colours, edgecolor=edges, linewidth=0.3,
                zorder=2)
        v_cb = float(cb[f"{c}_matched_pct"])
        lines[c] = v_cb
        if ref_line != "off":  # from the ClockBase bar down through the five model rows (not the Always row)
            y_top = cb["y"] - BAR_H / 2
            y_bot = float(cd.loc[MODELS, "y"].min()) - BAR_H / 2
            ax.plot([v_cb, v_cb], [y_top, y_bot], color=CLOCKBASE, lw=REF_LW, linestyle=ref_linestyle,
                    solid_capstyle="butt", dash_capstyle="butt", zorder=3)
        # A count that would run past the next column's axis (or the panel edge) is set inside its bar end.
        limit = COL_W + COL_GAP - COUNT_CLEAR_MM if k < 2 else COL_W + TAIL_W - EDGE_CLEAR_MM
        for row in cd.itertuples():
            v = getattr(row, f"{c}_matched_pct")
            anchor = v_cb if (ref_line == "aligned" and row.system in MODELS) else v
            text = getattr(row, f"{c}_label")
            is_cb = row.system == "ClockBase"
            fs, weight = (CB_COUNT_FS, "bold") if is_cb else (COUNT_FS, "normal")
            outside = anchor * mm_per_pct + pad_mm + text_width_mm(fig, text, fs, weight) <= limit
            if outside:
                t = ax.annotate(text, xy=(anchor, row.y), xytext=(COUNT_PAD_PT, 0), textcoords="offset points",
                                ha="left", va="center", fontsize=fs, fontweight=weight, annotation_clip=False,
                                bbox=KNOCKOUT if ref_line == "behind" else None, zorder=4)
            else:
                t = ax.annotate(text, xy=(v, row.y), xytext=(-COUNT_PAD_PT, 0), textcoords="offset points",
                                ha="right", va="center", fontsize=fs, fontweight=weight, annotation_clip=False,
                                gid="inside", zorder=4, color="white" if is_cb else INK)
            counts[(row.system, c)] = (t, "outside" if outside else "inside", k)
        ax.set_xlim(0, 100)
        ax.set_xticks([0, 50, 100], ["0", "", "100"])  # 50 keeps its gridline and tick, unlabelled
        if k == 0:
            ax.set_yticks(cd["y"], [DISPLAY[s] for s in cd["system"]], fontsize=TICK_FS)
            for tick in ax.get_yticklabels():
                if tick.get_text() == "ClockBase":
                    tick.set_fontweight("bold")
        else:
            ax.set_yticks([])
        class_header(fig, COL_X0[k] + COL_W / 2, c, header_n[c])
        axes.append(ax)
    # One axis title for the three columns: the middle column is centred on the three.
    assert abs((COL_X0[1] + COL_W / 2) - (COL_X0[0] + COL_X0[2] + COL_W) / 2) < 1e-9
    axis_title_baseline(axes[1], "Published calls matched (%)")

    # Checks at print size.
    problems = check_layout(fig, "i")
    for (s, c), (t, where, k) in counts.items():
        x0, _, x1, _ = bbox_mm(t, fig)
        limit = COL_X0[k] + COL_W + COL_GAP - COUNT_CLEAR_MM if k < 2 else W_I - EDGE_CLEAR_MM
        if x1 > limit + 0.05:
            problems.append(f"i: count {t.get_text()} ({s}, {c}) ends at {x1:.2f} mm, limit {limit:.2f}")
    inside = sorted((s, c) for (s, c), (_, where, _) in counts.items() if where == "inside")
    assert inside == [("Always 'higher'", "higher")], inside
    tick_gaps = []
    for k in range(2):
        hundred = [t for t in axes[k].get_xticklabels() if t.get_text() == "100"][0]
        zero = [t for t in axes[k + 1].get_xticklabels() if t.get_text() == "0"][0]
        gap = bbox_mm(zero, fig)[0] - bbox_mm(hundred, fig)[2]
        tick_gaps.append(gap)
        if gap < TICK_GAP_MIN_MM:
            problems.append(f"i: '100' of column {k + 1} is {gap:.2f} mm from the next column's '0'")
    # Distance from the longest model bar to ClockBase's value, per column (where the line stands).
    line_gaps = {c: (lines[c] - float(cd.loc[MODELS, f"{c}_matched_pct"].max())) * mm_per_pct for c in CLASSES}
    if ref_line != "off":
        for k, c in enumerate(CLASSES):
            if line_gaps[c] < 0.5:
                problems.append(f"i: longest model bar in {c} is {line_gaps[c]:.2f} mm from the ClockBase line")
        if ref_line == "aligned":  # no count touches the line
            for (s, c), (t, where, k) in counts.items():
                if s not in MODELS:
                    continue
                x_line = COL_X0[k] + lines[c] * mm_per_pct + REF_LW / 2 * 25.4 / 72
                x0 = bbox_mm(t, fig)[0]
                if x0 - x_line < REF_CLEAR_PT * 25.4 / 72:
                    problems.append(f"i: count {t.get_text()} ({s}, {c}) starts {x0 - x_line:.2f} mm from the line")
    return fig, axes, problems, counts, tick_gaps, line_gaps


# ------------------------------------------------------------------------------------------------- panel j
def panel_accuracy(cd):
    fig = plt.figure(figsize=(W_J * MM, FIG_H * MM))
    ax = add_axes_mm(fig, W_J, J_LEFT, W_J - J_LEFT - J_RIGHT)
    colours, edges = row_colours(cd)
    ax.barh(cd["y"], cd["balanced_accuracy"], height=BAR_H, color=colours, edgecolor=edges, linewidth=0.3, zorder=2)
    cap = 0.2
    for row, colour in zip(cd.itertuples(), colours):
        v, lo, hi, y = row.balanced_accuracy, row.interval_low, row.interval_high, row.y
        if hi > lo:  # Always 'higher' has a zero-width interval: nothing to draw
            inner = "white" if colour == CLOCKBASE else INK  # the part inside a dark bar is drawn white
            ax.plot([lo, v], [y, y], color=inner, lw=0.6, solid_capstyle="butt", zorder=3)
            ax.plot([lo, lo], [y - cap / 2, y + cap / 2], color=inner, lw=0.6, zorder=3)
            ax.plot([v, hi], [y, y], color=INK, lw=0.6, solid_capstyle="butt", zorder=3)
            ax.plot([hi, hi], [y - cap / 2, y + cap / 2], color=INK, lw=0.6, zorder=3)
        is_cb = row.system == "ClockBase"
        ax.annotate(row.balanced_accuracy_label, xy=(max(v, hi), y), xytext=(1.6, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=CB_COUNT_FS if is_cb else NOTE_FS,
                    fontweight="bold" if is_cb else "normal", annotation_clip=False)
    chance = 1 / 3
    ax.axvline(chance, color=CHANCE, linestyle=(0, (2.2, 1.4)), linewidth=0.6, zorder=2.5)  # under the intervals
    ax.annotate("Chance", xy=(chance, 1), xycoords=("data", "axes fraction"), xytext=(0, 1.2),
                textcoords="offset points", ha="center", va="bottom", fontsize=NOTE_FS, color=CHANCE,
                annotation_clip=False)
    ax.set_xlim(0, 1.0)
    ax.set_xticks([0, 0.5, 1.0], ["0", "0.5", "1.0"])
    axis_title_baseline(ax, "Balanced accuracy")
    ax.set_yticks([])
    return fig, ax, check_layout(fig, "j")


# ---------------------------------------------------------------------------------------------------- main
def row_centres_pt(fig, ax, ys):
    fig.canvas.draw()
    return [float(ax.transData.transform((0, y))[1] / fig.dpi * 72) for y in ys]


def save(fig, stem):
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / f"{stem}.pdf")
    fig.savefig(FIGDIR / f"{stem}.png", dpi=600)


def main():
    cd, header_n = load()
    check_h_header_line()
    fig_i, axes_i, prob_i, counts, tick_gaps, line_gaps = panel_identification(cd, header_n)
    fig_j, ax_j, prob_j = panel_accuracy(cd)
    problems = prob_i + prob_j
    if problems:
        raise SystemExit("Layout problems:\n  " + "\n  ".join(problems))
    yj = row_centres_pt(fig_j, ax_j, cd["y"])
    for ax in axes_i:
        yi = row_centres_pt(fig_i, ax, cd["y"])
        assert np.allclose(yi, yj, atol=1e-6), (yi, yj)
        assert np.allclose(ax.get_position().bounds[1::2], ax_j.get_position().bounds[1::2])
        assert ax.get_ylim() == ax_j.get_ylim() == YLIM
    assert abs(X_I + W_I + GAP_IJ + W_J - 183.0) < 1e-9
    # The two axis titles share a baseline.
    bi = axes_i[1].xaxis.label.get_position()[1], axes_i[1].xaxis.label.get_va()
    bj = ax_j.xaxis.label.get_position()[1], ax_j.xaxis.label.get_va()
    assert bi[1] == bj[1] == "baseline"
    cd["row_centre_pt_i"], cd["row_centre_pt_j"] = np.round(yi, 4), np.round(yj, 4)
    for c in CLASSES:
        cd[f"{c}_count_position"] = [counts[(s, c)][1] for s in cd["system"]]
        cd[f"{c}_count_size_pt"] = [counts[(s, c)][0].get_fontsize() for s in cd["system"]]
        cd[f"{c}_header_n"] = header_n[c]
        cd[f"{c}_count_anchor_pct"] = [float(counts[(s, c)][0].xy[0]) for s in cd["system"]]
    save(fig_i, "benchmark_identification_bars")
    save(fig_j, "benchmark_accuracy_bars")
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["system", "printed_name", "source_key"]
    for c in CLASSES:
        cols += [f"{c}_header_n", f"{c}_hits", f"{c}_total", f"{c}_matched_pct", f"{c}_label",
                 f"{c}_count_position", f"{c}_count_anchor_pct", f"{c}_count_size_pt"]
    cols += ["recall_mean", "balanced_accuracy", "interval_low", "interval_high", "balanced_accuracy_label", "y",
             "row_centre_pt_i", "row_centre_pt_j"]
    cd[cols].to_csv(RESULT, index=False)
    print(f"i {W_I} x {FIG_H} mm: labels {LBL_W}, columns {COL_W} with gaps {COL_GAP}, tail {TAIL_W:.2f} mm; "
          f"j {W_J} x {FIG_H} mm")
    print("clear space '100' to next '0' (mm):", [round(g, 2) for g in tick_gaps])
    print(f"ClockBase line: {REF_LINE}; longest model bar to ClockBase's value (mm):",
          {c: round(g, 2) for c, g in line_gaps.items()})
    for s, a, b in zip(cd["system"], yi, yj):
        print(f"row centre (pt from page bottom)  {s:16s} i = {a:7.3f}  j = {b:7.3f}")


if __name__ == "__main__":
    main()
