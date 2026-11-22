#!/usr/bin/env python3
"""
Combined-modality radar plot showing intervention categories pooled across all
three modalities (mouse RNA-seq + human RNA-seq + human DNAm). Outputs a single
radar with one polygon for "All interventions" and one polygon for
"Significant only", with category counts pooled across modalities.

Drawn at its print size (56 x 50 mm, Figure 4a) and placed at 100% scale.

Outputs:
  figures/create_radar_plot_combined.png
  figures/create_radar_plot_combined.pdf
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.textpath import TextPath

plt.rcParams["font.family"] = "Arial"
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["mathtext.it"] = "Arial:italic"
plt.rcParams["mathtext.bf"] = "Arial:bold"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

# Load combined 3-modality plot data
df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)

# Analytic set: exclude extreme-tail-flagged rows (mirrors analysis_2/3/4_combined)
if "outlier_status" in df.columns:
    df = df[df["outlier_status"] == "Normal"].copy()

# Canonicalise condition_category to the canonical 5 buckets (others become Other)
CANONICAL = {"Disease", "Drug", "Environment", "Genetic", "Other"}
df["condition_category"] = df["condition_category"].where(
    df["condition_category"].isin(CANONICAL), "Other"
)

# Order kept consistent with create_radar_plot.py
categories = ["Disease", "Drug", "Environment", "Genetic", "Other"]

# Pool counts across all modalities (no stratification)
# "Significant" is FDR<0.05 (matches the figure caption); the is_significant
# column additionally applies the score floor and is not used here.
all_counts = df["condition_category"].value_counts()
sig_df = df[df["fdr"] < 0.05]
sig_counts = sig_df["condition_category"].value_counts()

data_all = [all_counts.get(cat, 0) for cat in categories]
data_sig = [sig_counts.get(cat, 0) for cat in categories]
pd.DataFrame({'condition_category': categories, 'all_rows': data_all,
              'significant_rows': data_sig}).to_csv('result/panel_a_category_counts.csv', index=False)

# Style
BG_WHITE = "#ffffff"
GREY70 = "#b3b3b3"
COLOR_ALL = "#B8B8B8"
COLOR_SIG = "#FF6B6B"

VARIABLES = categories
VARIABLES_N = len(VARIABLES)
ANGLES = [n / VARIABLES_N * 2 * np.pi for n in range(VARIABLES_N)]
ANGLES += ANGLES[:1]

# ---------------------------------------------------------------------------
# Print geometry. Every length below is in points (1 pt = 1/72 inch); the
# figure is created at 72 dpi so display units equal points during layout.
# ---------------------------------------------------------------------------
W_MM, H_MM = 56.0, 50.0
W_PT, H_PT = W_MM / 25.4 * 72, H_MM / 25.4 * 72
R_PT = 58.2            # radius of the outer (25k) ring
CX_PT = 78.5           # radar centre
CY_PT = 73.4
LABEL_PAD = 2.4        # gap between the outer ring and category names
FS_CAT, FS_VAL, FS_TICK, FS_LEG = 7.0, 6.0, 6.0, 6.5
LW_DATA, LW_GRID = 0.9, 0.5
MARKER_D = 3.0         # marker diameter, pt
VALUE_GAP = 2.6        # gap between a grey vertex and its value pair
VALUE_LINE = 6.4       # baseline step between the two lines of a value pair

fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4), dpi=72)
ax = fig.add_axes(
    [(CX_PT - R_PT) / W_PT, (CY_PT - R_PT) / H_PT, 2 * R_PT / W_PT, 2 * R_PT / H_PT],
    polar=True,
)
fig.patch.set_facecolor(BG_WHITE)
ax.set_facecolor(BG_WHITE)

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# Display range only: the largest value is 20,523 (Genetic), so the radial axis
# stops at 25k, keeping the July 10% inner offset (ylim -2,500 to 25,000).
RMAX = 25000
ax.set_ylim(-RMAX * 0.1, RMAX)

data_all_plot = data_all + [data_all[0]]
data_sig_plot = data_sig + [data_sig[0]]

ax.plot(
    ANGLES, data_all_plot, c=COLOR_ALL, linewidth=LW_DATA,
    label="All Interventions", linestyle="--", alpha=0.7,
)
ax.scatter(ANGLES, data_all_plot, s=MARKER_D ** 2, c=COLOR_ALL, zorder=10,
           alpha=0.7, linewidths=0)

ax.plot(ANGLES, data_sig_plot, c=COLOR_SIG, linewidth=LW_DATA, label="Significant Only")
ax.scatter(ANGLES, data_sig_plot, s=MARKER_D ** 2, c=COLOR_SIG, zorder=10,
           linewidths=0)

ax.fill(ANGLES, data_all_plot, COLOR_ALL, alpha=0.1)
ax.fill(ANGLES, data_sig_plot, COLOR_SIG, alpha=0.15)

ax.set_xticks(ANGLES[:-1])
ax.set_xticklabels([])
RTICKS = [5000, 10000, 15000, 20000, 25000]
RTICK_LABELS = ["5k", "10k", "15k", "20k", "25k"]
ax.set_yticks(RTICKS)
ax.set_yticklabels([])
ax.yaxis.grid(True, linestyle="--", alpha=0.3, linewidth=LW_GRID)
ax.xaxis.grid(True, linestyle="--", alpha=0.3, linewidth=LW_GRID)
ax.tick_params(axis="both", length=0, pad=0)

ax.spines["start"].set_color("none")
ax.spines["polar"].set_color("none")

renderer = fig.canvas.get_renderer()


def screen_unit(theta):
    """Unit vector on the page for polar angle theta (offset pi/2, clockwise)."""
    phi = np.pi / 2 - theta
    return np.array([np.cos(phi), np.sin(phi)])


def align_for(v, h_cut=0.3, v_cut=0.45):
    ha = "left" if v[0] > h_cut else ("right" if v[0] < -h_cut else "center")
    va = "bottom" if v[1] > v_cut else ("top" if v[1] < -v_cut else "center")
    return ha, va


def ink_box(t):
    """Glyph-ink box of a text artist (layout box less ascender gap and descent)."""
    bb = t.get_window_extent(renderer)
    size = t.get_fontsize()
    return np.array([bb.x0, bb.y0 + 0.2 * size, bb.x1, bb.y1 - 0.2 * size])


def rect_gap(a, b):
    """Signed gap between two boxes (negative when they overlap)."""
    dx = max(b[0] - a[2], a[0] - b[2])
    dy = max(b[1] - a[3], a[1] - b[3])
    if dx < 0 and dy < 0:
        return max(dx, dy)
    return np.hypot(max(dx, 0), max(dy, 0))


def point_rect_dist(p, r):
    """Distance from points p (n x 2) to box r; negative depth when inside."""
    p = np.atleast_2d(p)
    dx = np.maximum.reduce([r[0] - p[:, 0], np.zeros(len(p)), p[:, 0] - r[2]])
    dy = np.maximum.reduce([r[1] - p[:, 1], np.zeros(len(p)), p[:, 1] - r[3]])
    out = np.hypot(dx, dy)
    inside = (dx == 0) & (dy == 0)
    depth = np.minimum.reduce([p[:, 0] - r[0], r[2] - p[:, 0], p[:, 1] - r[1], r[3] - p[:, 1]])
    return np.where(inside, -depth, out)


CENTER = ax.transData.transform((0, -RMAX * 0.1))
POLY = {}
for name, vals in (("all", data_all_plot), ("sig", data_sig_plot)):
    POLY[name] = ax.transData.transform(np.c_[ANGLES, vals])
SEG_POINTS = np.vstack([
    np.linspace(p[i], p[i + 1], 120) for p in POLY.values() for i in range(len(p) - 1)
])
MARKERS = np.vstack([POLY["all"][:-1], POLY["sig"][:-1]])


def clearance(box, boxes=(), inside_ring=True):
    """Smallest distance from a text box to the polygons, markers, other text and ring."""
    gaps = [point_rect_dist(SEG_POINTS, box).min() - LW_DATA / 2]
    gaps.append(point_rect_dist(MARKERS, box).min() - MARKER_D / 2)
    gaps += [rect_gap(box, b) for b in boxes]
    if inside_ring:
        corners = [(box[0], box[1]), (box[0], box[3]), (box[2], box[1]), (box[2], box[3])]
        far = max(np.hypot(*(np.array(c) - CENTER)) for c in corners)
        gaps.append(R_PT - far)
    return min(gaps)


# Category names outside the outer ring, aligned away from the circle.
cat_boxes = []
placed_texts = []      # category names and value labels, for the ring-label search
for angle, cat in zip(ANGLES[:-1], categories):
    u = screen_unit(angle)
    ha, va = align_for(u)
    t = ax.annotate(cat, xy=(angle, RMAX), xytext=tuple(LABEL_PAD * u),
                    textcoords="offset points", ha=ha, va=va,
                    fontsize=FS_CAT, fontweight="bold", annotation_clip=False)
    cat_boxes.append(ink_box(t))
    placed_texts.append(t)

# Value labels: the pair (all rows in grey above significant rows in red) sits
# just beyond each grey vertex, turned away from both polygons.
value_boxes = []
for idx, angle in enumerate(ANGLES[:-1]):
    u = screen_unit(angle)
    best = None
    for k in (0, 1, -1, 2, -2, 3, -3):
        phi = np.arctan2(u[1], u[0]) + k * np.pi / 6
        v = np.array([np.cos(phi), np.sin(phi)])
        ha, va = align_for(v)
        off = (VALUE_GAP + MARKER_D / 2) * v
        if va == "bottom":
            offs = [(off[0], off[1] + VALUE_LINE), (off[0], off[1])]
        elif va == "top":
            offs = [(off[0], off[1]), (off[0], off[1] - VALUE_LINE)]
        else:
            offs = [(off[0], off[1] + VALUE_LINE / 2), (off[0], off[1] - VALUE_LINE / 2)]
        texts = []
        for (dx, dy), val, col in zip(offs, (data_all[idx], data_sig[idx]),
                                      (COLOR_ALL, COLOR_SIG)):
            texts.append(ax.annotate(
                f"{val:,}", xy=(angle, data_all[idx]), xytext=(dx, dy),
                textcoords="offset points", ha=ha, va=va if va != "center" else (
                    "bottom" if dy > 0 else "top"),
                fontsize=FS_VAL, fontweight="bold", color=col, annotation_clip=False,
                bbox=dict(boxstyle="square,pad=0.05", fc=BG_WHITE, ec="none")))
        boxes = [ink_box(t) for t in texts]
        block = np.array([min(b[0] for b in boxes), min(b[1] for b in boxes),
                          max(b[2] for b in boxes), max(b[3] for b in boxes)])
        score = clearance(block, cat_boxes + value_boxes)
        if best is None or (score > best[0] + 0.05 and (best[0] < 1.0)):
            if best is not None:
                for t in best[2]:
                    t.remove()
            best = (score, block, texts, k)
        else:
            for t in texts:
                t.remove()
    value_boxes.append(best[1])
    placed_texts.extend(best[2])
    print(f"value pair {categories[idx]:12s} rotation step {best[3]:+d}  clearance {best[0]:.2f} pt")

# Radial ring labels: one column on a single ray in the gap between Disease
# and Other (or Disease and Drug), away from every value label, category name
# and vertex. Each label sits just outside its own ring: the label ink nearest
# the centre lies TICK_GAP beyond the ring line. Distances use exact glyph ink.
# Among rays that keep every label at least TICK_TEXT_GAP from other text, the
# search keeps the one with the most room from the polygons and markers.
TICK_GAP = 0.4                     # ring line to the nearest label ink, pt
TICK_TEXT_GAP = 1.2 / 25.4 * 72    # minimum space to other text (1.2 mm), pt
TICK_MUTUAL_GAP = 1.0              # minimum space between neighbouring ring labels, pt
TICK_SECTORS = ((296.0, 352.0), (8.0, 64.0))   # Other-Disease and Disease-Drug, deg
CORNER_PT = 4.0 / 25.4 * 72        # top-left square kept free for the panel letter


def glyph_box(t):
    """Exact glyph-ink box of a single-line text artist, in points."""
    wb = t.get_window_extent(renderer)
    _, _, descent = t._get_layout(renderer)
    e = TextPath((0, 0), t.get_text(), prop=t.get_fontproperties()).get_extents()
    base = wb.y0 + descent
    return np.array([wb.x0 + e.x0, base + e.y0, wb.x0 + e.x1, base + e.y1])


OTHER_INK = [glyph_box(t) for t in placed_texts]
TICK_PROP = FontProperties(family="Arial", size=FS_TICK)
TICK_EXT = {lab: np.array(TextPath((0, 0), lab, prop=TICK_PROP).get_extents().extents)
            for lab in RTICK_LABELS}


def ring_pt(r):
    """Ring radius on the page, in points from the centre."""
    return R_PT * (r + RMAX * 0.1) / (RMAX * 1.1)


def tick_boxes(theta):
    """Ink boxes of the ring labels centred on the ray at theta, each just outside its ring."""
    u = screen_unit(theta)
    out = []
    for r, lab in zip(RTICKS, RTICK_LABELS):
        e = TICK_EXT[lab]
        half = np.r_[(e[:2] - e[2:]) / 2, (e[2:] - e[:2]) / 2]
        target = ring_pt(r) + TICK_GAP
        lo, hi = 0.0, R_PT + 20.0
        for _ in range(40):            # nearest ink to the centre at the target radius
            mid = (lo + hi) / 2
            c = CENTER + mid * u
            if point_rect_dist(CENTER, half + np.r_[c, c])[0] < target:
                lo = mid
            else:
                hi = mid
        c = CENTER + hi * u
        out.append(half + np.r_[c, c])
    return out


best_ray = None
for lo_deg, hi_deg in TICK_SECTORS:
    for deg in np.arange(lo_deg, hi_deg + 0.01, 0.25):
        theta = np.deg2rad(deg)
        boxes = tick_boxes(theta)
        data_gap = min(clearance(b, inside_ring=False) for b in boxes)
        text_gap = min(rect_gap(b, c) for b in boxes for c in OTHER_INK)
        mutual = min(rect_gap(boxes[i], boxes[i + 1]) for i in range(len(boxes) - 1))
        on_page = all(b[0] > 0 and b[1] > 0 and b[2] < W_PT and b[3] < H_PT and
                      not (b[0] < CORNER_PT and b[3] > H_PT - CORNER_PT) for b in boxes)
        if text_gap < TICK_TEXT_GAP or mutual < TICK_MUTUAL_GAP or not on_page:
            continue
        key = (round(data_gap, 2), text_gap)
        if best_ray is None or key > best_ray[0]:
            best_ray = (key, deg, boxes, data_gap, text_gap, mutual)
if best_ray is None:
    raise RuntimeError("no ray keeps the ring labels clear of the other text")
_, tick_deg, tick_ink, tick_data, tick_text, tick_mutual = best_ray
print(f"ring labels on the {tick_deg:.2f} deg ray, just outside their rings; clearance "
      f"to polygons {tick_data:.2f} pt, to other text {tick_text:.2f} pt "
      f"({tick_text * 25.4 / 72:.2f} mm), between labels {tick_mutual:.2f} pt")
for r, lab, b in zip(RTICKS, RTICK_LABELS, tick_ink):
    e = TICK_EXT[lab]
    anchor = ax.transData.transform((np.deg2rad(tick_deg), r))
    ax.annotate(lab, xy=(np.deg2rad(tick_deg), r),
                xytext=(b[0] - e[0] - anchor[0], b[1] - e[1] - anchor[1]),
                textcoords="offset points", ha="left", va="baseline",
                fontsize=FS_TICK, color=GREY70, annotation_clip=False)

handles = [
    Line2D([], [], c=COLOR_SIG, lw=LW_DATA, marker="o", markersize=MARKER_D,
           markeredgewidth=0, label=f"Significant ($n$ = {sum(data_sig):,})"),
    Line2D([], [], c=COLOR_ALL, lw=LW_DATA, marker="o", markersize=MARKER_D,
           markeredgewidth=0, label=f"All ($n$ = {sum(data_all):,})",
           linestyle="--", alpha=0.7),
]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
           ncol=2, frameon=False, fontsize=FS_LEG, handlelength=2.0,
           handletextpad=0.4, columnspacing=1.4, borderaxespad=0.15,
           borderpad=0.1)

# Layout report: union of all ink on the page, in mm from each edge.
fig.canvas.draw()
ink = [a.get_window_extent(renderer) for a in fig.findobj(
    lambda a: hasattr(a, "get_text") and a.get_text().strip() and a.get_visible())]
x0 = min(b.x0 for b in ink); y0 = min(b.y0 for b in ink)
x1 = max(b.x1 for b in ink); y1 = max(b.y1 for b in ink)
y0 = min(y0, CY_PT - R_PT)
to_mm = 25.4 / 72
print(f"text margins mm: left {x0 * to_mm:.2f} right {(W_PT - x1) * to_mm:.2f} "
      f"bottom {y0 * to_mm:.2f} top {(H_PT - y1) * to_mm:.2f}")

plt.savefig("figures/create_radar_plot_combined.pdf", facecolor=BG_WHITE)
plt.savefig("figures/create_radar_plot_combined.png", dpi=600, facecolor=BG_WHITE)
print("Saved figures/create_radar_plot_combined.{pdf,png}")

print(f"\nTotal interventions: {sum(data_all):,}")
print(f"Significant interventions: {sum(data_sig):,}")
print("\nBreakdown by category (pooled):")
print("-" * 40)
for cat, all_count, sig_count in zip(categories, data_all, data_sig):
    pct = (sig_count / all_count * 100) if all_count > 0 else 0
    print(f"{cat:12} - All: {all_count:6,}  Sig: {sig_count:6,} ({pct:.1f}%)")
