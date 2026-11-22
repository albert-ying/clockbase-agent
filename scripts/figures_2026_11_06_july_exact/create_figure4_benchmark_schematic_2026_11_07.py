#!/usr/bin/env python3
"""Figure 4h: native-vector schematic of AgingAgentBench, drawn at print size (78 x 54 mm, placed at 100 %).

It sits to the left of the two bar panels (i, published effects identified; j, balanced accuracy) drawn by
create_figure4_benchmark_bars_2026_11_07.py and explains the benchmark to a reader who has not read the Methods,
reading left to right:
  1. Ground truth: papers report aging-clock effects for GEO series in the atlas; each effect is a record with its
     two groups, direction and significance.
  2. Two systems answer each effect. The two groups go to both: ClockBase computes the clock the paper reported from
     the deposited samples of the paper's two groups and tests the difference with the test of the verification;
     five language models receive the GEO description of the series, the group definitions and a neutral
     description of the clock and the measure, but no data.
  3. Scoring: each system's call (lower, higher or not significant; grid rows, titled 'System call') is compared with
     the published call (grid columns, titled 'Published call'), which runs straight from the record to the column
     title, on the results scored for every system. The published calls of those results (decreases, increases, not
     significant; the order of panel i) are counted under the grid columns. The other published effects could not be
     scored with the paper's groups (reasons in Extended Data Fig. 23a); no funnel is drawn.

Everything is drawn in matplotlib from native shapes and Arial text: no AI-generated or embedded raster imagery.
Every printed number is read at build time from the manuscript text module (manuscript_v3_texts_2026_10_01,
T._B and T._BL) and from ed23/numbers.json (c_d.S3c.call_classes, the published calls of the scored results); the
build fails if any value disagrees with the checked values. The numbers are written to
result/figure4_benchmark_schematic_numbers.csv. Glyph contents (matrix tints, dot positions, clock hands, text
lines on cards) are symbolic and carry no values.
"""
import csv
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_DIR = ROOT / "paper/revision_2026-09-07/source"
WORK = ROOT / "paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact"
NUMBERS_JSON = WORK.parent / "ed23" / "numbers.json"
FIGDIR = WORK / "figures"
RESULT = WORK / "result" / "figure4_benchmark_schematic_numbers.csv"
STEM = "benchmark_schematic"

# Style tokens shared with the bar panels i and j (create_figure4_benchmark_bars_2026_11_07.py).
INK = "#222222"
CLOCKBASE = "#33414E"   # ClockBase bars of j, box edge and glyphs
MODEL = "#A6A6A6"       # language-model bars of j, box edge
MODEL_TEXT = "#6E6E6E"  # language-model name and connectors
LOWER = "#0072B2"
HIGHER = "#D55E00"
NS = "#8C8C8C"
PAPER_LINE = "#B4B4B4"
PAGE_BACK = "#E6E6E6"
CB_FILL = "#F1F3F5"
LM_FILL = "#F4F4F4"
AGREE = "#6B6B6B"
CALL_COLOUR = {"higher": HIGHER, "lower": LOWER, "ns": NS}
CALL_ORDER = ("lower", "higher", "ns")  # key, grid rows and grid columns; panel i lists decreases first
# The two groups of a published effect: neutral greys, so that the ClockBase slate means only ClockBase and the
# model grey (MODEL) means only the language models.
GROUP_A, GROUP_B = "#4D4D4D", "#C6C6C6"
MIN_GAP = 0.5  # mm of clear space between an arrowhead and any label (text or call marker)

W_MM, H_MM = 78.0, 54.0
MM = 1 / 25.4
LABEL_FS, HEAD_FS = 6.5, 7.0
LANE_FS = 6.0   # ClockBase step labels: three two-line labels share the 31.8 mm box
TITLE_FS = 6.0  # axis titles of the scoring grid
LANE_PITCH = 2.3  # mm between the baselines of a two-line step label
LW_ARROW = 0.6
HEAD_L, HEAD_W = 2.4, 1.3  # arrowhead length and half-width, pt
LW_GLYPH = 0.5
LINE_PITCH = 2.45  # mm between label baselines
# Header line (Ground truth, ClockBase, Scoring): text top 1.4 mm below the page top, level with the top of the
# legend of panel i (create_figure4_benchmark_bars_2026_11_07.py), which sits beside this panel.
HEAD_Y = 51.0

plt.rcParams.update({
    "font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.facecolor": "white", "figure.facecolor": "white", "text.color": INK,
    "savefig.pad_inches": 0.02,  # only used by tight-bbox measurements; the panel is saved at its exact page size
})

EXPECTED = {"truth_results": 580, "truth_papers": 93, "truth_series": 125, "records": 110, "papers_scored": 32,
            "decreases": 16, "increases": 36, "models": 5, "not_significant": 58}
WORDS = {2: "Two", 3: "Three", 4: "Four", 5: "Five", 6: "Six", 7: "Seven", 8: "Eight", 9: "Nine"}


# ----------------------------------------------------------------------------------------------- numbers
def load_numbers():
    if not NUMBERS_JSON.exists():
        raise SystemExit(f"Cannot proceed - required data file not found: {NUMBERS_JSON}")
    sys.path.insert(0, str(SOURCE_DIR))
    import manuscript_v3_texts_2026_10_01 as T  # slow (about a minute): builds every manuscript number
    B, BL = T._B, T._BL
    nj = json.loads(NUMBERS_JSON.read_text())
    calls = nj["c_d"]["S3c"]["call_classes"]
    nums = {
        "truth_results": B["truth_results"], "truth_papers": B["truth_papers"], "truth_series": B["truth_series"],
        "records": BL["records"], "decreases": BL["decreases"], "increases": BL["increases"],
        "models": len(BL["models"]), "not_significant": calls["none"],
        "papers_scored": B["groups"]["clusters"],
    }
    # The papers of the scored results: the released-test results ClockBase tested are the results every system
    # answered (same count), so their papers are the papers of the 110.
    if B["groups_results"] != nums["records"]:
        raise SystemExit(f"tested results {B['groups_results']} != results scored for every system {nums['records']}")
    bad = {k: (nums[k], v) for k, v in EXPECTED.items() if nums[k] != v}
    if bad:
        raise SystemExit(f"Benchmark numbers differ from the checked values (found, expected): {bad}")
    if nj["ground_truth_total"] != nums["truth_results"]:
        raise SystemExit(f"ed23/numbers.json ground_truth_total {nj['ground_truth_total']} != {nums['truth_results']}")
    # The same 110 results carry the published calls of panels i and j.
    if (calls["lower"], calls["higher"], sum(calls.values()), nj["c_d"]["S3c"]["records"]) != (
            nums["decreases"], nums["increases"], nums["records"], nums["records"]):
        raise SystemExit(f"ed23/numbers.json call classes {calls} disagree with T._BL")
    rows = [
        ("truth_results", nums["truth_results"], "T._B['truth_results']", "published effects of the ground truth"),
        ("truth_papers", nums["truth_papers"], "T._B['truth_papers']", "papers of the ground truth"),
        ("truth_series", nums["truth_series"], "T._B['truth_series']", "GEO series of the ground truth"),
        ("records", nums["records"], "T._BL['records']", "results scored for every system"),
        ("papers_scored", nums["papers_scored"], "T._B['groups']['clusters'] (T._B['groups_results'] == records)",
         "papers of the results scored for every system"),
        ("increases", nums["increases"], "T._BL['increases']; ed23 c_d.S3c.call_classes.higher",
         "published increases among them (grid column 'higher')"),
        ("decreases", nums["decreases"], "T._BL['decreases']; ed23 c_d.S3c.call_classes.lower",
         "published decreases among them (grid column 'lower')"),
        ("not_significant", nums["not_significant"], "ed23/numbers.json c_d.S3c.call_classes.none",
         "published not significant among them (grid column 'not significant')"),
        ("models", nums["models"], "len(T._BL['models'])", "language models: " + ", ".join(BL["models"])),
        ("ground_truth_total_check", nj["ground_truth_total"], "ed23/numbers.json ground_truth_total", "check"),
    ]
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    with RESULT.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["quantity", "value", "source", "meaning"])
        w.writerows(rows)
    return nums


# ----------------------------------------------------------------------------------------------- drawing helpers
class Canvas:
    def __init__(self):
        # 600 dpi: text extents measured at the default 100 dpi are hinted and up to 0.4 mm short at 6 pt
        self.fig = plt.figure(figsize=(W_MM * MM, H_MM * MM), dpi=600)
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, W_MM)
        self.ax.set_ylim(0, H_MM)
        self.ax.set_axis_off()
        self.ax.set_xticks([])  # no hidden tick labels left in the figure
        self.ax.set_yticks([])
        self.texts = []
        self.connectors = []  # polylines (mm) checked against text boxes
        self.heads = []       # arrowheads, checked for clear space around labels
        self.markers = []     # call markers (labels of the key, the record card and the grid)

    def text(self, x, y, s, size=LABEL_FS, weight="normal", color=INK, ha="left", va="center", **kw):
        t = self.ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va, **kw)
        self.texts.append(t)
        return t

    def line(self, pts, color=INK, lw=LW_ARROW, z=3):
        xs, ys = zip(*pts)
        self.ax.add_line(Line2D(xs, ys, color=color, lw=lw, solid_capstyle="butt", solid_joinstyle="miter",
                                zorder=z))
        self.connectors.append(list(pts))

    def arrow(self, pts, color=INK, lw=LW_ARROW, z=3):
        """Polyline through pts (mm) with a small filled head at the last point."""
        if len(pts) > 2:
            self.line(pts[:-1], color=color, lw=lw, z=z)
        a = FancyArrowPatch(pts[-2], pts[-1], arrowstyle=f"-|>,head_length={HEAD_L},head_width={HEAD_W}",
                            mutation_scale=1,
                            color=color, lw=lw, shrinkA=0, shrinkB=0, joinstyle="miter", capstyle="butt", zorder=z)
        self.ax.add_patch(a)
        self.heads.append((tuple(pts[-2]), tuple(pts[-1])))
        self.connectors.append([pts[-2], pts[-1]])

    def rect(self, x, y, w, h, fc="white", ec=INK, lw=LW_GLYPH, z=2, round_mm=0.0):
        if round_mm:
            p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={round_mm}", fc=fc, ec=ec, lw=lw,
                               zorder=z)
        else:
            p = Rectangle((x, y), w, h, fc=fc, ec=ec, lw=lw, zorder=z)
        self.ax.add_patch(p)
        return p

    def hline(self, x0, x1, y, color=PAPER_LINE, lw=0.45, z=6):
        self.ax.add_line(Line2D([x0, x1], [y, y], color=color, lw=lw, zorder=z, solid_capstyle="butt"))

    def marker(self, x, y, kind, size_pt=3.4, z=7):
        m = {"higher": "^", "lower": "v", "ns": "o"}[kind]
        colr = {"higher": HIGHER, "lower": LOWER, "ns": NS}[kind]
        ln, = self.ax.plot([x], [y], marker=m, ms=size_pt, mfc="white" if kind == "ns" else colr, mec=colr,
                           mew=0.6, ls="none", zorder=z)
        self.markers.append(ln)


def paper_stack(c, x, y, w, h, n=3, step=0.75):
    """Stack of published papers: offset pages; the front one has a title bar and text lines."""
    for i in reversed(range(n)):
        c.rect(x + i * step, y + i * step, w, h, fc="white" if i == 0 else PAGE_BACK, ec=INK, lw=LW_GLYPH,
               z=3 + (n - i))
    c.hline(x + 0.6, x + w - 0.6, y + h - 0.9, color=INK, lw=0.9, z=3 + n + 1)
    k, yy = 0, y + h - 1.9
    while yy > y + 0.7:
        c.hline(x + 0.6, x + w - 0.6 - (1.1 if k % 3 == 2 else 0), yy, z=3 + n + 1)
        k, yy = k + 1, yy - 0.72


def record_stack(c, x, y, w, h, rows_y, n=3, step=0.6):
    """Stack of effect records (cards with the same three fields): two groups, direction, significance."""
    for i in reversed(range(n)):
        c.rect(x + i * step, y + i * step, w, h, fc="white" if i == 0 else PAGE_BACK, ec=INK, lw=LW_GLYPH,
               z=3 + (n - i), round_mm=0.5)
    z = 3 + n + 1
    gy, dy, py = rows_y
    cx = x + w / 2
    for dx, colr in ((-0.65, GROUP_A), (0.65, GROUP_B)):
        c.ax.add_patch(Rectangle((cx + dx - 0.5, gy - 0.5), 1.0, 1.0, fc=colr, ec="none", zorder=z))
    c.marker(cx - 0.9, dy, "lower", z=z)
    c.marker(cx + 0.9, dy, "higher", z=z)
    c.text(cx, py, "$P$", ha="center", va="center", zorder=z)
    for yy in ((gy + dy) / 2, (dy + py) / 2):
        c.hline(x + 0.5, x + w - 0.5, yy, color=PAPER_LINE, lw=0.35, z=z)


# Symbolic tint pattern for the data matrices (no values: positions of five tints only).
PATTERN = [[0, 2, 1, 3, 4, 2], [3, 1, 4, 0, 2, 1], [1, 4, 2, 2, 0, 3], [4, 0, 3, 1, 3, 0], [2, 3, 0, 4, 1, 4]]
TINTS = ["#33414E", "#5E6A76", "#8C959E", "#B9BFC5", "#E1E4E7"]
GREY_TINTS = ["#D2D2D2", "#DBDBDB", "#E3E3E3", "#EAEAEA", "#F1F1F1"]


def data_matrix(c, x, y, cell=1.0, gap=0.5, tints=TINTS, group_bars=True):
    """Sample-by-gene matrix: two blocks of three sample columns (the paper's two groups), five gene rows."""
    for r in range(5):
        for k in range(6):
            gx = x + k * cell + (gap if k >= 3 else 0)
            gy = y + (4 - r) * cell
            c.ax.add_patch(Rectangle((gx, gy), cell, cell, fc=tints[PATTERN[r][k]], ec="white", lw=0.25, zorder=4))
    width = 6 * cell + gap
    if group_bars:
        c.ax.add_patch(Rectangle((x, y + 5 * cell + 0.3), 3 * cell, 0.6, fc=GROUP_A, ec="none", zorder=4))
        c.ax.add_patch(Rectangle((x + 3 * cell + gap, y + 5 * cell + 0.3), 3 * cell, 0.6, fc=GROUP_B, ec="none",
                                 zorder=4))
    return width, 5 * cell + (0.9 if group_bars else 0)


def clock_face(c, cx, cy, r, color=CLOCKBASE):
    c.ax.add_patch(Circle((cx, cy), r, fc="white", ec=color, lw=0.8, zorder=4))
    for k in range(12):
        a = math.pi / 2 - k * math.pi / 6
        r0 = r * (0.70 if k % 3 == 0 else 0.80)
        c.ax.add_line(Line2D([cx + r0 * math.cos(a), cx + 0.92 * r * math.cos(a)],
                             [cy + r0 * math.sin(a), cy + 0.92 * r * math.sin(a)], color=color,
                             lw=0.6 if k % 3 == 0 else 0.4, zorder=5, solid_capstyle="butt"))
    for ang, ln, lw in ((math.radians(60), 0.46, 0.9), (math.radians(-30), 0.64, 0.6)):
        c.ax.add_line(Line2D([cx, cx + ln * r * math.cos(ang)], [cy, cy + ln * r * math.sin(ang)], color=color,
                             lw=lw, zorder=5, solid_capstyle="round"))
    c.ax.add_patch(Circle((cx, cy), 0.24, fc=color, ec="none", zorder=6))


def group_test(c, cx, y0, color=CLOCKBASE):
    """Clock outputs of the two groups (symbolic dots) under a test bracket."""
    a_pts = [(0.0, 0.6), (0.45, 1.4), (-0.4, 2.0), (0.3, 2.7)]
    b_pts = [(0.1, 2.1), (-0.45, 2.9), (0.4, 3.5), (-0.1, 4.2)]
    ax0, bx0 = cx - 1.55, cx + 1.55
    for dx, dy in a_pts:
        c.ax.add_patch(Circle((ax0 + dx, y0 + dy), 0.42, fc=GROUP_A, ec="white", lw=0.25, zorder=5))
    for dx, dy in b_pts:
        c.ax.add_patch(Circle((bx0 + dx, y0 + dy), 0.42, fc=GROUP_B, ec="white", lw=0.25, zorder=5))
    top = y0 + 5.2
    for pts in ([(ax0, top - 0.6), (ax0, top), (bx0, top), (bx0, top - 0.6)], [(cx - 2.7, y0), (cx + 2.7, y0)]):
        xs, ys = zip(*pts)
        c.ax.add_line(Line2D(xs, ys, color=color, lw=0.5, zorder=5, solid_capstyle="butt"))


def description_card(c, x, w, rows_y, edge=MODEL):
    """Text card with one field per input, level with its label: the GEO description of the series (text lines),
    the definitions of the two groups (swatches) and the description of the clock (a small clock); no data."""
    gy, dy, cy = rows_y  # GEO description, group definitions, clock description
    pad = 1.15
    c.rect(x, cy - pad, w, gy - cy + 2 * pad, fc="white", ec=edge, lw=LW_GLYPH, z=4, round_mm=0.5)
    for k in range(2):
        c.hline(x + 0.6, x + w - 0.6 - (1.3 if k == 1 else 0), gy + 0.36 - k * 0.72, z=5)
    for k, colr in enumerate((GROUP_A, GROUP_B)):
        yy = dy + 0.45 - k * 0.9
        c.ax.add_patch(Rectangle((x + 0.6, yy - 0.32), 0.64, 0.64, fc=colr, ec="none", zorder=5))
        c.hline(x + 1.6, x + w - 0.8 - k * 0.9, yy, z=5)
    r, ccx = 0.5, x + 0.92
    c.ax.add_patch(Circle((ccx, cy), r, fc="white", ec=MODEL_TEXT, lw=0.4, zorder=5))
    for ang, ln in ((math.radians(60), 0.5), (math.radians(-30), 0.7)):
        c.ax.add_line(Line2D([ccx, ccx + ln * r * math.cos(ang)], [cy, cy + ln * r * math.sin(ang)],
                             color=MODEL_TEXT, lw=0.35, zorder=6, solid_capstyle="round"))
    c.hline(x + 1.6, x + w - 0.8 - 0.5, cy, z=5)


def no_data(c, x, y):
    w, h = data_matrix(c, x, y, tints=GREY_TINTS, group_bars=False)
    c.ax.add_line(Line2D([x - 0.3, x + w + 0.3], [y - 0.3, y + h + 0.3], color=MODEL_TEXT, lw=0.8, zorder=6,
                         solid_capstyle="butt"))
    return w, h


# ----------------------------------------------------------------------------------------------- composition
def build(nums):
    c = Canvas()

    # Column 1: ground truth --------------------------------------------------------------------------------
    c.text(4.6, HEAD_Y, "AgingAgentBench", size=HEAD_FS, weight="bold")
    paper_stack(c, 0.8, 37.4, w=4.8, h=7.6, n=3, step=0.75)
    tx, ty = 8.2, 42.0
    c.text(tx, ty + LINE_PITCH, f"{nums['truth_papers']:,} papers")
    c.text(tx, ty, f"{nums['truth_series']:,} GEO series")
    c.text(tx, ty - LINE_PITCH, "in the atlas")
    # papers report effects: one record per published effect
    rows_y = (21.0, 17.1, 13.2)  # two groups, direction, significance
    rec_x, rec_y, rec_w, rec_h = 0.8, 11.0, 6.4, 12.4
    arrow_top, arrow_end = 36.8, rec_y + rec_h + 1.2 + 0.5
    c.arrow([(3.2, arrow_top), (3.2, arrow_end)])
    c.text(4.6, (arrow_top + arrow_end) / 2, f"{nums['truth_results']:,} effects")
    record_stack(c, rec_x, rec_y, rec_w, rec_h, rows_y, n=3, step=0.6)
    fx = 9.4
    for yy, lab in zip(rows_y, ("Two groups", "Direction", "Significance")):
        c.text(fx, yy, lab)

    # Column 2: the two systems -----------------------------------------------------------------------------
    bx0, bx1 = 26.8, 58.4
    cb_y0, cb_y1 = 32.8, HEAD_Y + 2.0  # the box title shares the header line
    lm_y0, lm_y1 = 6.0, 30.0
    c.rect(bx0, cb_y0, bx1 - bx0, cb_y1 - cb_y0, fc=CB_FILL, ec=CLOCKBASE, lw=0.6, z=1, round_mm=0.8)
    c.rect(bx0, lm_y0, bx1 - bx0, lm_y1 - lm_y0, fc=LM_FILL, ec=MODEL, lw=0.6, z=1, round_mm=0.8)
    c.text(bx0 + 1.2, cb_y1 - 2.0, "ClockBase", size=HEAD_FS, weight="bold", color=CLOCKBASE)
    c.text(bx0 + 1.2, lm_y1 - 2.0, f"{WORDS[nums['models']]} language models", size=HEAD_FS, weight="bold",
           color=MODEL_TEXT)

    # ClockBase: the paper's two groups in the deposited samples -> the clock the paper reported -> the test of the
    # verification. Step centres are set by the label widths (6 pt): labels 1.4 to 1.8 mm apart and 0.7 mm inside
    # the box.
    cb_mid = 43.6
    s1, s2, s3 = 32.25, 42.35, 52.5
    mw, mh = 6.5, 5.9
    data_matrix(c, s1 - mw / 2, cb_mid - mh / 2)
    clock_face(c, s2, cb_mid, r=2.9)
    group_test(c, s3, cb_mid - 2.7)
    c.arrow([(s1 + mw / 2 + 0.5, cb_mid), (s2 - 2.9 - 0.5, cb_mid)])
    c.arrow([(s2 + 2.9 + 0.5, cb_mid), (s3 - 2.7 - 0.5, cb_mid)])
    lab_y = 37.85  # label block centred between the glyphs and the box bottom
    for sx, l1, l2 in ((s1, "Deposited", "samples"), (s2, "Paper’s", "clock"), (s3, "Verification", "test")):
        c.text(sx, lab_y, l1, size=LANE_FS, ha="center")
        c.text(sx, lab_y - LANE_PITCH, l2, size=LANE_FS, ha="center")

    # Language models: GEO description, group definitions and clock description, no data
    card_x, card_w = 29.4, 6.6
    lm_in = rows_y[0]
    in_rows = (lm_in + LINE_PITCH, lm_in, lm_in - LINE_PITCH)
    description_card(c, card_x, card_w, in_rows)
    no_data(c, card_x + (card_w - mw) / 2, 11.6 - 2.5)
    lx = card_x + card_w + 1.2
    for yy, lab in zip(in_rows, ("GEO description", "Group definitions", "Clock description")):
        c.text(lx, yy, lab)
    c.text(lx, 11.6, "No data")

    # Each published effect's two groups go to both systems
    bus_x = 25.4
    c.line([(fx + 12.4, lm_in), (bus_x, lm_in)])
    c.line([(bus_x, lm_in), (bus_x, cb_mid)])
    c.arrow([(bus_x, cb_mid), (s1 - mw / 2 - 0.5, cb_mid)])
    c.arrow([(bus_x, lm_in), (card_x - 0.5, lm_in)])

    # Column 3: scoring -------------------------------------------------------------------------------------
    x3 = 60.8
    c.text(x3, HEAD_Y, "Scoring", size=HEAD_FS, weight="bold")
    call_label = {"lower": "Lower", "higher": "Higher", "ns": "Not significant"}
    for k, kind in enumerate(CALL_ORDER):
        yy = HEAD_Y - 3.5 - k * LINE_PITCH
        c.marker(x3 + 0.9, yy, kind)
        c.text(x3 + 2.2, yy, call_label[kind])
    for k, lab in enumerate((f"{nums['records']:,} results", f"of {nums['papers_scored']:,} papers", "scored for",
                             "every system")):
        c.text(x3, HEAD_Y - 11.9 - k * LINE_PITCH, lab)

    # 3 x 3 grid: rows = system call, columns = published call; agreement on the diagonal. The grid sits far enough
    # left that the column title, centred under it, ends inside the page.
    cell = 3.4
    gx0, gy0 = 65.9, 21.3 - LINE_PITCH  # one text line lower: the scoring note has four lines
    grid_mid = gy0 + 1.5 * cell
    col_x = gx0 + 1.5 * cell
    for i in range(3):
        for j in range(3):
            c.ax.add_patch(Rectangle((gx0 + j * cell, gy0 + (2 - i) * cell), cell, cell,
                                     fc=CALL_COLOUR[CALL_ORDER[i]] if i == j else "white", ec=INK, lw=0.4,
                                     zorder=4))
    count = {"lower": nums["decreases"], "higher": nums["increases"], "ns": nums["not_significant"]}
    for i, kind in enumerate(CALL_ORDER):
        c.marker(gx0 - 1.1, gy0 + (2 - i) * cell + cell / 2, kind, size_pt=3.0)
        c.marker(gx0 + i * cell + cell / 2, gy0 - 1.1, kind, size_pt=3.0)
        c.text(gx0 + i * cell + cell / 2, gy0 - 3.4, f"{count[kind]:,}", ha="center")
    # Axis titles: columns under the counts, rows (rotated) left of the row markers
    col_title_y = gy0 - 5.8
    col_title = c.text(col_x, col_title_y, "Published call", size=TITLE_FS, ha="center")
    row_title_x = gx0 - 2.95
    row_title = c.text(row_title_x, grid_mid, "System call", size=TITLE_FS, ha="center", va="center",
                       rotation=90, rotation_mode="anchor")

    # System calls: both outputs leave their box at its right edge and join one arrow that ends at the left of the
    # row title
    mx = 59.1
    out_x = s3 + 2.7 + 0.4
    c.line([(out_x, cb_mid), (mx, cb_mid), (mx, grid_mid)])
    c.line([(out_x, lm_in), (mx, lm_in), (mx, grid_mid)])
    c.fig.canvas.draw()
    rt = row_title.get_window_extent(c.fig.canvas.get_renderer())
    c.arrow([(mx, grid_mid), (rt.x0 / c.fig.dpi * 25.4 - MIN_GAP - 0.05, grid_mid)])

    # Published call: from the records, under both systems, up to the column title of the grid
    base_y = 0.7  # level with the bottom of the axis titles of panels i and j
    ct = col_title.get_window_extent(c.fig.canvas.get_renderer())
    c.line([(3.4, rec_y - 0.4), (3.4, base_y), (col_x, base_y)])
    c.arrow([(col_x, base_y), (col_x, ct.y0 / c.fig.dpi * 25.4 - MIN_GAP - 0.05)])
    return c


# ----------------------------------------------------------------------------------------------- checks
def _seg_hits_box(p, q, box, pad=0.25):
    """True if segment p-q (axis-aligned or not) passes through box (x0, y0, x1, y1) grown by pad."""
    x0, y0, x1, y1 = box[0] - pad, box[1] - pad, box[2] + pad, box[3] + pad
    for t in [k / 200 for k in range(201)]:
        x, y = p[0] + (q[0] - p[0]) * t, p[1] + (q[1] - p[1]) * t
        if x0 < x < x1 and y0 < y < y1:
            return True
    return False


def check_layout(c):
    """Fail on text outside the page, text overlapping text, a connector through text, or ink in the letter corner."""
    fig = c.fig
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    to_mm = lambda bb: (bb.x0 / fig.dpi * 25.4, bb.y0 / fig.dpi * 25.4, bb.x1 / fig.dpi * 25.4,
                        bb.y1 / fig.dpi * 25.4)
    boxes = [(t.get_text(), to_mm(t.get_window_extent(r))) for t in c.texts]
    problems = []
    for s, (x0, y0, x1, y1) in boxes:
        if x0 < 0.3 or y0 < 0.3 or x1 > W_MM - 0.3 or y1 > H_MM - 0.3:
            problems.append(f"outside page: {s!r} {x0:.1f},{y0:.1f},{x1:.1f},{y1:.1f}")
        if x0 < 4.0 and y1 > H_MM - 4.0:
            problems.append(f"in letter corner: {s!r}")
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if a[0] < b[2] - 0.05 and b[0] < a[2] - 0.05 and a[1] < b[3] - 0.05 and b[1] < a[3] - 0.05:
                problems.append(f"text overlap: {boxes[i][0]!r} / {boxes[j][0]!r}")
    for pts in c.connectors:
        for p, q in zip(pts[:-1], pts[1:]):
            for s, bb in boxes:
                if _seg_hits_box(p, q, bb):
                    problems.append(f"connector {p}->{q} through text {s!r}")
    for art in c.ax.patches + c.ax.lines:
        bb = to_mm(art.get_window_extent(r))
        if bb[0] < 4.0 and bb[3] > H_MM - 4.0:
            problems.append(f"ink in letter corner: {art}")
    # Arrowheads keep MIN_GAP of clear space from every label: text boxes and call markers (marker box grown by
    # the edge stroke and its miter at the triangle corners).
    pt = 25.4 / 72
    labels = list(boxes)
    for m in c.markers:
        x0, y0, x1, y1 = to_mm(m.get_window_extent(r))
        g = m.get_markeredgewidth() * pt
        labels.append((f"marker {m.get_marker()}", (x0 - g, y0 - g, x1 + g, y1 + g)))
    heads = []
    for (px, py), (qx, qy) in c.heads:
        n = math.hypot(qx - px, qy - py)
        ux, uy = (qx - px) / n, (qy - py) / n
        back = (HEAD_L + LW_ARROW) * pt + 0.5 * LW_ARROW * pt / math.sin(math.atan2(HEAD_W, HEAD_L))
        side = (HEAD_W + LW_ARROW / 2) * pt
        corners = [(qx, qy), (qx - ux * back - uy * side, qy - uy * back + ux * side),
                   (qx - ux * back + uy * side, qy - uy * back - ux * side)]
        hb = (min(p[0] for p in corners), min(p[1] for p in corners), max(p[0] for p in corners),
              max(p[1] for p in corners))
        gap, name = min((max(bb[0] - hb[2], hb[0] - bb[2], bb[1] - hb[3], hb[1] - bb[3]), s) for s, bb in labels)
        heads.append(((qx, qy), gap, name))
        if gap < MIN_GAP:
            problems.append(f"arrowhead at ({qx:.2f}, {qy:.2f}) is {gap:.2f} mm from label {name!r}")
    return boxes, heads, problems


def save(c, figdir=FIGDIR, stem=STEM):
    figdir.mkdir(parents=True, exist_ok=True)
    c.fig.savefig(figdir / f"{stem}.pdf")
    c.fig.savefig(figdir / f"{stem}.svg")
    c.fig.savefig(figdir / f"{stem}.png", dpi=600)
    print("wrote", figdir / f"{stem}.pdf")


def main():
    nums = load_numbers()
    c = build(nums)
    boxes, heads, problems = check_layout(c)
    for s, bb in boxes:
        print(f"text {s!r:28s} x {bb[0]:5.1f}-{bb[2]:5.1f}  y {bb[1]:5.1f}-{bb[3]:5.1f}")
    for (qx, qy), gap, name in heads:
        print(f"arrowhead at ({qx:5.2f}, {qy:5.2f}): {gap:4.2f} mm from nearest label {name!r}")
    if problems:
        raise SystemExit("Layout check failed:\n  " + "\n  ".join(problems))
    save(c)


if __name__ == "__main__":
    main()
