#!/usr/bin/env python3
"""Figure 3 in its compact form: one score per panel, ranks listed beside every row.

  a    the combined-modality volcano (panel_a.pdf of the July generator), placed as vector at 0.87 without its
       title line, with the keys beside it
  b-e  statistical ranking per category: the five strongest increases and decreases (the drug decreases continue
       down to ouabain, the first drug comparison that meets all three judgement conditions); the gutter gives the
       rank among the 748 increases or 198 decreases; bars give |Hedges' g| (CI bound nearest zero); three condition
       cells and the reason text beside every decrease
  f    replication ranking: pooled Hedges' g with its 95% interval, projects, series, FDR and best single rank
  g    judgement version: admitted lower-output comparisons in statistical order,
       one entry per intervention and GEO series; bars give the absolute adjusted
       confidence bound nearest zero as in c-e.

Inputs: evidence/.../figure3/fig3_rankings.csv written by build_fig3_rankings_table_2026_11_08.py (one row per
intervention and direction) and that module's comparison frames; panel_a.pdf and the July build records in the
same folder. Every printed number comes from these frames. Colours come from fig3_palette_2026_11_08.py.
Output: evidence/.../figure3/figure3_compact.pdf (+ .png, statistical ranks and build record).
Run (after build_july_exact_2026_11_06.py figure3 and package_figures_2026_11_06.py figure3):
  source env_local.sh; .venv-pipeline/bin/python -B build_fig3_rankings_table_2026_11_08.py
  .venv-pipeline/bin/python -B build_figure3_compact_2026_11_08.py
"""
import sys

sys.dont_write_bytecode = True

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pymupdf  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
try:
    import build_fig3_rankings_table_2026_11_08 as B  # noqa: E402  (inserts the repository script folder into sys.path)
except SystemExit as exc:
    raise SystemExit(f"Cannot proceed - {exc}")
import fig3_palette_2026_11_08 as P3  # noqa: E402

be_fig, g_fig = B.be_fig, B.g_fig
NAME = "compact"
TABLE = B.EVIDENCE / "fig3_rankings.csv"
OUTPUT = B.EVIDENCE.parents[1] / "figure-fixes-2026-11-17/figure3"
OUTPUT.mkdir(parents=True, exist_ok=True)
PANEL_A = OUTPUT / "panel_a.pdf"
STRIP = OUTPUT / "figure3_schematic.pdf"   # method strip under panel a
STRIP_H = 39.5
# Display names read against the GEO records (9 November 2026): GSE148730 compares inflamed with steady-state lymph-node
# lymphatic endothelium of wild-type mice; GSE189688 compares ST2-knockout Tregs transplanted with AML cells with a
# transplant without tumour. The atlas labels named the studies' knockouts.
RELABEL = {"Clever-1 Deficiency": "Lymph-node inflammation", "ST2-KO": "ST2-KO Tregs in AML"}
ARIAL_BOLD = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
OUT_PDF = OUTPUT / f"figure3_{NAME}.pdf"
OUT_PNG = OUTPUT / f"figure3_{NAME}.png"
for p in (TABLE, PANEL_A):
    if not p.exists():
        raise SystemExit(f"Cannot proceed - required data file not found: {p}")

PAL = P3.active()
LOW, HIGH = PAL["lower"], PAL["higher"]
EDGE = {"lower": PAL["lower_edge"], "higher": PAL["higher_edge"]}
FILL = {"lower": LOW, "higher": HIGH}
INK = PAL["ink"]
GREY = PAL["reason_text"]
CAT = PAL["category"]
MARK = be_fig.MODALITY_MARKERS
HATCH = be_fig.MODALITY_HATCH
DAGGER = " †"

# ---------------------------------------------------------------- type and geometry (mm, pt)
W = 183.0
S_A = 0.87                  # uniform scale of panel_a.pdf (its 6 pt labels print at 5.2 pt)
A_OFF = 3.6                 # mm: panel a starts below its letter
A_CLIP_TOP_PT = 9.0         # panel_a.pdf above this line holds only the title "ClockBase Agent, ..."
FS = 6.0                    # row labels, numbers, ticks
FS_HEAD = 5.5               # column heads, reason text, key
FS_AX = 6.0                 # axis titles
FS_TITLE = 7.0              # panel titles
FS_LETTER = 10.0            # panel letters (v25)
P = 2.36                    # row pitch
BAR = 1.48                  # bar height = cell side
CELL_GAP = 0.38
BLOCK_GAP = 1.05            # between the increase and decrease blocks of one list
GUT = 4.6                   # right edge of the rank gutter from the panel's left edge
LAB = 3.4                   # left edge of the labels (no rank gutter since v37)
COL_GAP = 4.0               # between the two columns of panels
CW = (W - COL_GAP) / 2
X_LEFT, X_RIGHT = 0.0, CW + COL_GAP
TITLE_H = 3.5               # panel title line (letter and title) above the first row
HEAD_1 = 2.5                # one line of column heads
AXIS_H = 3.5                # axis line, ticks and tick labels below the last row
AXIS_T = 2.7                # axis title
PANEL_GAP = 2.6
SHOW_CELLS = False          # the judgement conditions define panel g; b-e carry no condition cells
LW_AXIS, LW_TICK, TICK = 0.5, 0.5, 0.75
MARKER_PT, MARKER_EW = 3.0, 0.4
MM_PER_PT = 25.4 / 72

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "Arial",
                     "font.sans-serif": ["Arial"], "hatch.linewidth": 0.35, "savefig.facecolor": "white",
                     "axes.unicode_minus": True})


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def minus(s):
    return s.replace("-", "−")


# ---------------------------------------------------------------- data
def state_key(s):
    return {"met": "met", "not met": "not met", "not stated or not agreed": "unknown"}[s]


def load_data():
    t = pd.read_csv(TABLE, low_memory=False)
    for col in ("statistical_display_entry", "followed_up"):
        t[col] = t[col].astype(str).eq("True")
    src, rec = B.load()
    checks = []
    cmp_ = B.comparisons(src, rec, checks)
    jt = B.judgement_tables(src, cmp_, checks)
    return t, cmp_, jt, checks


def followed_label(text, followed):
    return text + DAGGER if followed else text


def be_rows(t):
    """Per category: the 5 strongest increases and decreases (b-e rule); the drug decreases continue to the first
    drug entry whose comparison meets all three conditions (ouabain), so every drug above it is shown."""
    ent = t[t["statistical_display_entry"]]
    out = {}
    for cat in be_fig.CATEGORIES:
        hi = ent[(ent.category == cat) & (ent.direction == "higher")].sort_values("statistical_rank_in_category")
        lo = ent[(ent.category == cat) & (ent.direction == "lower")].sort_values("statistical_rank_in_category")
        n_lo = be_fig.TOP_N
        hi = hi[hi.statistical_rank_in_category <= be_fig.TOP_N].copy()
        lo = lo[lo.statistical_rank_in_category <= n_lo].copy()
        for part in (hi, lo.head(be_fig.TOP_N)):
            assert part["printed_be"].notna().all(), f"{cat}: a top-5 row is not a v25 b-e bar"
        for part in (hi, lo):
            text = []
            for r in part.itertuples():
                if isinstance(r.printed_be, str):
                    s = r.printed_be
                else:  # rows added below the top 5: the ranking label, or the case arm when the label is its stem
                    arm = r.comparison_name.split(" vs ")[0].strip()
                    s = arm if (" " not in arm and arm != r.label and arm.startswith(r.label)) else r.label
                text.append(followed_label(s, bool(r.followed_up)))
            part["text"] = text
            part["bold"] = [x.endswith(DAGGER) for x in text]
        # the v25 reason text equals the table's reason rule on every v25 decrease
        v25 = lo["printed_reason"].notna()
        assert (lo.loc[v25, "printed_reason"] == lo.loc[v25, "reason"]).all()
        out[cat] = (hi, lo)
    assert sum(len(h) + len(l.head(5)) for h, l in out.values()) == 40
    return out


QUALIFY_DURATION = re.compile(r"\b(\d+\s?w)\b")
QUALIFY_PAREN = re.compile(r"\(([^)]+)\)\s*$")


def qualify(label, comparison_name, modality):
    """A label that another panel prints for a different comparison of the same intervention takes a qualifier
    from its own comparison name: an age word (g rule, '(old mice)'), else a duration ('12w'), else a closing
    parenthesis ('(cortex)')."""
    m = re.search(r"\b(old|aged|young)\b", comparison_name, flags=re.IGNORECASE)
    if m:
        return f"{label} ({m.group(1).lower()} {g_fig.SPECIES[modality]})"
    m = QUALIFY_DURATION.search(comparison_name)
    if m and m.group(1).replace(" ", "") not in label:
        return f"{label} {m.group(1).replace(' ', '')}"
    m = QUALIFY_PAREN.search(comparison_name)
    if m:
        return f"{label} ({m.group(1).lower()})"
    raise AssertionError(f"no qualifier for {label!r} in {comparison_name!r}")


def g_rows(jt, t, printed):
    """All 8 admitted decreases in judgement order. Text: the v25 g rule for every row (label without ' at ZT<n>',
    the senolytic arm named, the age qualifier where b-e prints the label for another comparison, the dagger);
    rows that would read alike take their ZT time inside the qualifier."""
    g = jt["lower"].copy()
    assert g["judgement_rank"].tolist() == list(range(1, len(g) + 1))
    disp = g["label"].str.replace(r" at ZT\d+", "", regex=True)
    disp[g["comparison_name"].eq("Overload vs Sham (Senolytic)")] = "Overload in senolytic-treated mice"
    g["display"] = disp
    out = []
    for _, r in g.iterrows():
        s = r["display"]
        other = [(k, c) for k, c in printed.get(s.casefold(), []) if c != (r["gse_id"], r["comparison_name"])]
        if other:
            s = g_fig.qualifier(r)
        if (r["gse_id"], r["comparison_name"]) == g_fig.TESTED_KEY:
            s += g_fig.TESTED_MARK
        out.append(s)
    g["text"] = out
    kept = t[t["printed_g"].notna()].set_index("judgement_gse_id")
    for r in g[g["kept"]].itertuples():   # the five v25 rows keep their v25 text
        assert r.text == kept.loc[r.gse_id, "printed_g"], (r.text, kept.loc[r.gse_id, "printed_g"])
    dup = g["text"].duplicated(keep=False)
    zt = g["comparison_name"].str.extract(r"\b(ZT\d+)\b")[0]
    g.loc[dup, "text"] = [s[:-1] + f", {z})" if s.endswith(")") else f"{s} ({z})" for s, z in zip(g.loc[dup, "text"],
                                                                                              zt[dup])]
    assert g["text"].is_unique
    g["bold"] = g["text"].str.endswith(DAGGER)
    return g


def f_rows(t):
    f = t[t["printed_f"].notna()].sort_values("replication_rank_unit").copy()
    assert len(f) == len(pd.read_csv(B.SOURCES['replication_lower']))
    assert f["replication_rank_unit"].tolist() == list(range(1, len(f) + 1))
    f["text"] = f["printed_f"]
    f["bold"] = f["text"].str.endswith(DAGGER)
    return f


# ---------------------------------------------------------------- drawing
class Canvas:
    def __init__(self, height):
        self.H = height
        self.fig = plt.figure(figsize=(W / 25.4, height / 25.4))
        self.ax = self.fig.add_axes([0, 0, 1, 1])
        self.ax.set_xlim(0, W)
        self.ax.set_ylim(height, 0)
        self.ax.axis("off")
        self.texts = []
        self.renderer = self.fig.canvas.get_renderer()

    def text(self, x, y, s, size=FS, weight="normal", color=INK, ha="left", va="center_baseline", role="", **kw):
        t = self.ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va, **kw)
        t._role = role
        self.texts.append(t)
        return t

    def rich(self, x, y, parts, size=FS_TITLE, va="baseline"):
        """Consecutive runs [(text, weight, color)] on one baseline."""
        out = []
        for s, wgt, col in parts:
            t = self.text(x, y, s, size=size, weight=wgt, color=col, va=va, role="title")
            out.append(t)
            x += width_mm(s, size, wgt)
        return out

    def line(self, xs, ys, lw=LW_AXIS, color=INK, **kw):
        ln = Line2D(xs, ys, linewidth=lw, color=color, solid_capstyle="butt", **kw)
        self.ax.add_line(ln)
        return ln

    def rect(self, x, y, w, h, **kw):
        r = Rectangle((x, y), w, h, **kw)
        self.ax.add_patch(r)
        return r

    def marker(self, x, y, m, size=MARKER_PT, face="#f7f7f7", edge="#222", ew=MARKER_EW, z=6, alpha=1.0):
        self.ax.plot([x], [y], marker=m, ms=size, mfc=face, mec=edge, mew=ew, ls="none", clip_on=False, zorder=z,
                     alpha=alpha)

    def cells(self, x0, yc, states):
        for i, s in enumerate(states):
            for a in be_fig.cell_artists(x0 + i * (BAR + CELL_GAP), yc - BAR / 2, BAR, BAR, s, self.ax.transData):
                if isinstance(a, Rectangle):
                    self.ax.add_patch(a)
                else:
                    self.ax.add_line(a)


_FP = {}
_MEASURE = plt.figure(figsize=(2, 2))
_R = _MEASURE.canvas.get_renderer()


def width_mm(s, size=FS, weight="normal"):
    key = (size, weight)
    if key not in _FP:
        _FP[key] = FontProperties(family="Arial", size=size, weight=weight)
    w, _, _ = _R.get_text_width_height_descent(s, _FP[key], ismath=False)
    return w / _MEASURE.dpi * 25.4


def bar(cv, x0, yc, length, direction, modality, left=False):
    """A bar of the b-e style (fill, edge, modality hatch) with the modality end marker."""
    x = x0 - length if left else x0
    cv.rect(x, yc - BAR / 2, length, BAR, facecolor=FILL[direction], edgecolor=EDGE[direction], linewidth=0.3,
            hatch=HATCH[modality] or None, zorder=3)
    cv.marker(x if left else x0 + length, yc, MARK[modality])


def letter(cv, x, y, s):
    cv.text(x, y, s, size=FS_LETTER, weight="bold", color="black", va="baseline", role="letter")


def panel_title(cv, x0, y0, letter_s, bold, rest, color=INK):
    letter(cv, x0, y0 + 2.75, letter_s)
    parts = [(bold, "bold", color)] + ([(" " + rest, "normal", INK)] if rest else [])
    cv.rich(x0 + 3.4, y0 + 2.6, parts)


def x_axis(cv, x_of, ticks, y, labels=None, title=None, title_x=None):
    """Horizontal axis line from the first to the last tick, ticks pointing down, labels below."""
    xs = [x_of(v) for v in ticks]
    cv.line([min(xs), max(xs)], [y, y], lw=LW_AXIS)
    for v, x, lab in zip(ticks, xs, labels or [None] * len(ticks)):
        cv.line([x, x], [y, y + TICK], lw=LW_TICK)
        s = lab if lab is not None else (f"{v:g}" if v >= 0 else minus(f"{v:g}"))
        cv.text(x, y + TICK + 0.35, s, size=FS, ha="center", va="top", role="tick")
    if title:
        cv.text(title_x if title_x is not None else (min(xs) + max(xs)) / 2, y + TICK + 0.35 + 2.55, title,
                size=FS_AX, ha="center", va="top", role="axis")


def rotated_heads(cv, x0, y_bottom, heads):
    for i, h in enumerate(heads):
        cv.text(x0 + i * (BAR + CELL_GAP) + BAR / 2, y_bottom, h, size=FS_HEAD, rotation=90, ha="left",
                va="center", rotation_mode="anchor", role="head", color=INK)


def draw_be(cv, x0, y0, cat, letter_s, hi, lo, geo, axis_title):
    color = CAT[cat]
    panel_title(cv, x0, y0, letter_s, cat, "", color=color)
    y = y0 + TITLE_H
    x_of = lambda v: geo["bx"] + v * geo["scale"]  # noqa: E731
    rows_y = []
    for k, part in enumerate((hi, lo)):
        off = 0.0
        for i, r in enumerate(part.itertuples()):
            off += getattr(r, "gap_before", 0.0)
            yc = y + off + (i + 0.5) * P
            rows_y.append(yc)
            cv.text(x0 + LAB, yc, r.text, weight="bold" if r.bold else "normal", role="label")
            bar(cv, geo["bx"], yc, r.statistical_score * geo["scale"], r.direction, r.modality)
            if SHOW_CELLS and r.direction == "lower":
                states = [state_key(getattr(r, c)) for c in be_fig.STATE_COLUMNS]
                cv.cells(geo["cells_x0"], yc, states)
                reason = r.reason if isinstance(r.reason, str) else ""
                if reason:
                    met = reason == be_fig.ALL_MET
                    cv.text(geo["cells_x0"] - 0.7, yc, reason, size=FS_HEAD, ha="right",
                            weight="bold" if met else "normal", color=PAL["met_text"] if met else GREY,
                            role="reason")
                    end = x_of(r.statistical_score) + MARKER_PT * MM_PER_PT / 2
                    assert geo["cells_x0"] - 0.7 - width_mm(reason, FS_HEAD, "bold" if met else "normal") > end + 1.0, \
                        f"reason text meets the bar of {r.text}"
        y += off + len(part) * P + (BLOCK_GAP if k == 0 else 0)
    first_lo_top = y0 + TITLE_H + len(hi) * P + BLOCK_GAP
    if SHOW_CELLS:
        rotated_heads(cv, geo["cells_x0"], first_lo_top - 0.45, be_fig.COLUMN_HEADS)
    top, bottom = y0 + TITLE_H, y
    for v in geo["ticks"]:
        if v > 0:
            cv.line([x_of(v)] * 2, [top, bottom], lw=0.3, color="#D0D0D0", linestyle=(0, (1, 1.6)), zorder=1)
    cv.line([x_of(0)] * 2, [top, bottom + 0.25], lw=LW_AXIS)
    x_axis(cv, x_of, geo["ticks"], bottom + 0.25, title=axis_title)
    return rows_y


def be_height(n_hi, n_lo, axis_title):
    return TITLE_H + (n_hi + n_lo) * P + BLOCK_GAP + 0.25 + AXIS_H + (AXIS_T if axis_title else 0)


# ---------------------------------------------------------------- columns of text beside a list
def column_layout(x_start, specs, gap=1.6):
    """specs: [(key, head_lines, values)] -> {key: (x_left, x_right)} sized to the widest head line or value."""
    pos, x = {}, x_start
    for key, heads, values in specs:
        w = max([width_mm(h, FS_HEAD) for h in heads] + [width_mm(v, FS) for v in values])
        pos[key] = (x, x + w)
        x += w + gap
    return pos, x - gap


def draw_heads(cv, pos, specs, y_bottom, align="right"):
    for key, heads, _ in specs:
        xl, xr = pos[key]
        for i, h in enumerate(reversed(heads)):
            cv.text(xr if align == "right" else xl, y_bottom - i * 2.05, h, size=FS_HEAD, ha=align, va="baseline",
                    role="head")


# ---------------------------------------------------------------- build
def main():
    t, cmp_, jt, checks = load_data()
    be = be_rows(t)
    f = f_rows(t)
    # labels printed for comparisons in b-e and g (for the qualifier rule of g and h)
    printed = {}
    for cat, (hi, lo) in be.items():
        for part in (hi, lo):
            for r in part.itertuples():
                printed.setdefault(r.text.replace(DAGGER, "").casefold(), []).append(
                    (r.display_key, (r.gse_id, r.comparison_name)))
    g = g_rows(jt, t, printed)
    for r in g.itertuples():
        printed.setdefault(r.text.replace(DAGGER, "").casefold(), []).append(
            (r.display_key, (r.gse_id, r.comparison_name)))

    n_low = int(t.loc[t.direction == "lower", "statistical_n_comparisons"].dropna().iloc[0])
    n_high = int(t.loc[t.direction == "higher", "statistical_n_comparisons"].dropna().iloc[0])
    assert n_low == len(cmp_[cmp_.direction.eq("lower")])
    assert n_high == len(cmp_[cmp_.direction.eq("higher")])

    # ---------------- c-e: the five strongest decreases of the genetic, drug and environmental categories, one row
    SHOWN = ["Genetic", "Drug", "Environment"]
    # c-e: five strongest increases and decreases of three categories; below a break, the comparisons that
    # nominated a follow-up experiment and rank lower in their category (ouabain and PLX5622 among the drugs)
    ent = t[t["statistical_display_entry"]]
    shown = {}
    for cat in SHOWN:
        hi, lo = be[cat]
        hi, lo = hi.copy(), lo.copy()
        hi["gap_before"], lo["gap_before"] = 0.0, 0.0
        more = ent[(ent.category == cat) & (ent.direction == "lower") & ent.followed_up
                   & (ent.statistical_rank_in_category > be_fig.TOP_N)].sort_values("statistical_rank_in_category").copy()
        if len(more):
            more["text"] = [x + DAGGER for x in more["label"]]
            more["bold"] = True
            more["gap_before"] = [BLOCK_GAP] + [0.0] * (len(more) - 1)
            lo = pd.concat([lo, more])
        shown[cat] = (hi, lo)
    be = shown
    for hi_, lo_ in be.values():
        for part in (hi_, lo_):
            part["text"] = [RELABEL.get(x, x).replace("Barretts", "Barrett's esophagus").replace("10 Gy x3", "10 Gy ×3") for x in part["text"]]
    assert [x for x in be["Drug"][1]["text"] if x.endswith(DAGGER)] == ["Ouabain †", "PLX5622 †"]
    assert "Cold exposure †" in be["Environment"][1]["text"].tolist()
    allrows = pd.concat([pd.concat(v) for v in be.values()])
    be_ticks = list(np.arange(0, np.ceil(float(allrows["statistical_score"].max()) / 4) * 4 + 1e-9, 4))
    lw = {cat: max(width_mm(s_, FS, "bold" if b_ else "normal") for hi, lo in [be[cat]] for part in (hi, lo)
                   for s_, b_ in zip(part["text"], part["bold"])) for cat in SHOWN}
    fixed = LAB + 1.4 + 2.6
    bar_w = (W - 2 * COL_GAP - sum(fixed + v for v in lw.values())) / 3
    geo, xs, x = {}, {}, 0.0
    for cat in SHOWN:
        xs[cat] = x
        geo[cat] = {"bx": x + LAB + lw[cat] + 1.4, "scale": bar_w / be_ticks[-1], "ticks": be_ticks}
        x += fixed + lw[cat] + bar_w + COL_GAP

    # ---------------- vertical layout
    a_h = (pymupdf.open(PANEL_A)[0].rect.height - A_CLIP_TOP_PT) * S_A * MM_PER_PT + A_OFF
    header_y = a_h + 1.2 + STRIP_H + 1.0
    y_rank = header_y + 3.4
    h_row = max(be_height(len(hi), len(lo), True) + float(lo["gap_before"].sum()) for hi, lo in be.values())
    hf = TITLE_H + 2 * 2.05 + 0.7 + len(f) * P + 0.25 + AXIS_H + AXIS_T
    hgm = TITLE_H + 2 * 2.05 + 0.7 + int(g["kept"].sum()) * P + 0.4
    H = y_rank + h_row + PANEL_GAP + max(hf, hgm) + 0.4
    cv = Canvas(H)

    a_w = pymupdf.open(PANEL_A)[0].rect.width * S_A * MM_PER_PT
    letter(cv, 0.2, 3.0, "a")
    draw_key(cv, a_w + 2.2, a_h)

    cv.rich(X_LEFT, header_y + 2.5, [("Statistical ranking", "bold", INK),
                                     (" of single comparisons", "normal", INK)])
    for cat, letter_s in zip(SHOWN, "cde"):
        draw_be(cv, xs[cat], y_rank, cat, letter_s, *be[cat], geo[cat], "|Hedges' g| bound")
    yf = y_rank + h_row + PANEL_GAP
    draw_f(cv, X_LEFT, yf, f, n_low)
    draw_g_judgement(cv, X_RIGHT, yf, g, be_ticks)

    check_overlaps(cv, a_box=(0, 0, a_w, a_h))
    tmp = OUTPUT / f"_{NAME}_ranking_layer.pdf"
    cv.fig.savefig(tmp, dpi=600)
    place_panel_a(tmp, H)
    tmp.unlink()
    render_png()
    records(t, be, f, g, checks, H, n_low, n_high)
    print(f"{OUT_PDF}\n{OUT_PNG}\nheight {H:.1f} mm")


def draw_key(cv, x, a_bottom):
    """Figure key beside panel a, bottom-aligned with it: bar colours, modality, condition cells, dagger."""
    entries = [("swatch", "lower", "lower clock output"), ("swatch", "higher", "higher clock output"),
               ("mod", "mouse RNA-seq", "mouse RNA-seq"), ("mod", "human RNA-seq", "human RNA-seq"),
               ("mod", "human DNAm", "human DNAm"),
               ("text", "†", "followed up in aged mice")]
    pitch = 2.35
    avail = W - x - 4.6
    lines = []
    for kind, what, s in entries:
        words, cur, out = s.split(), "", []
        for wd in words:
            trial = (cur + " " + wd).strip()
            if width_mm(trial, FS_HEAD) > avail and cur:
                out.append(cur)
                cur = wd
            else:
                cur = trial
        out.append(cur)
        lines.append((kind, what, out))
    n = sum(len(o) for _, _, o in lines) + 2  # two group gaps of half a line
    y = a_bottom - 1.4 - n * pitch + pitch / 2
    groups = {"mod": 1, "cell": 2, "text": 3}
    last = 0
    for kind, what, out in lines:
        gnum = groups.get(kind, 0)
        if gnum != last and gnum in (1, 2, 3) and last != 0:
            y += pitch * (0.5 if gnum != 3 else 0.5)
        last = gnum
        if kind == "swatch":
            cv.rect(x, y - BAR / 2, 3.4, BAR, facecolor=FILL[what], edgecolor=EDGE[what], linewidth=0.3)
        elif kind == "mod":
            cv.rect(x, y - BAR / 2, 2.4, BAR, facecolor="#9A9A9A", edgecolor="#333", linewidth=0.3,
                    hatch=HATCH[what] or None)
            cv.marker(x + 2.4, y, MARK[what], size=MARKER_PT)
        elif kind == "cell":
            for a in be_fig.cell_artists(x + 1.0, y - BAR / 2, BAR, BAR, what, cv.ax.transData):
                (cv.ax.add_patch if isinstance(a, Rectangle) else cv.ax.add_line)(a)
        else:
            cv.text(x + 1.7, y, what, size=FS, ha="center", role="key")
        for i, s in enumerate(out):
            cv.text(x + 4.2, y + i * pitch * 0.82, s, size=FS_HEAD, role="key")
        y += pitch + (len(out) - 1) * pitch * 0.82


def draw_f(cv, x0, y0, f, n_low):
    panel_title(cv, x0, y0, "f", "Replication ranking", "pooled across independent projects")
    labels_w = max(width_mm(s, FS, "bold" if b else "normal") for s, b in zip(f["text"], f["bold"]))
    fx0 = x0 + LAB + labels_w + 1.6
    ci_txt = [minus(f"{g:.2f} [{lo:.2f}, {hi:.2f}]") for g, lo, hi in
              zip(f.replication_pooled_g, f.replication_ci_low, f.replication_ci_high)]
    rank_txt = [f"{int(v)}" if not np.isnan(v) else "not ranked" for v in f.replication_best_single_rank]
    fdr_txt = [f"{v:.2g}" for v in f.replication_fdr]
    specs = [("k", ["Projects"], [f"{int(v)}" for v in f.replication_k]),
             ("series", ["Series"], [f"{int(v)}" for v in f.replication_series]),
             ("g", ["Pooled g", "(95% CI)"], ci_txt),
             ("fdr", ["FDR"], fdr_txt)]
    widths = sum(max([width_mm(hh, FS_HEAD) for hh in hs] + [width_mm(v, FS) for v in vals]) for _, hs, vals in specs)
    gaps = 1.6 * (len(specs) - 1)
    forest_w = (x0 + CW - 0.2) - fx0 - widths - gaps - 2.2
    lo_v = float(np.floor(f.replication_ci_low.min()))
    ticks = [lo_v, 0.0]
    sc = forest_w / (0 - lo_v)
    x_of = lambda v: fx0 + (v - lo_v) * sc  # noqa: E731
    pos, _ = column_layout(x_of(0) + 2.2, specs)
    heads_y = y0 + TITLE_H + 2 * 2.05 - 0.2
    draw_heads(cv, pos, specs, heads_y)
    ytop = y0 + TITLE_H + 2 * 2.05 + 0.7
    for i, r in enumerate(f.itertuples()):
        yc = ytop + (i + 0.5) * P
        cv.text(x0 + LAB, yc, r.text, weight="bold" if r.bold else "normal", color=CAT[r.replication_category],
                role="label")
        mod = r.replication_modalities
        assert mod in MARK, mod
        # bar from 0 to the bound nearest zero (the ranking score), interval, pooled g
        cv.rect(x_of(r.replication_ci_high), yc - BAR / 2, x_of(0) - x_of(r.replication_ci_high), BAR,
                facecolor=LOW, edgecolor=EDGE["lower"], linewidth=0.3, hatch=HATCH[mod] or None, zorder=3)
        cv.line([x_of(r.replication_ci_low), x_of(r.replication_ci_high)], [yc, yc], lw=0.6, color=INK, zorder=4)
        cv.marker(x_of(r.replication_ci_high), yc, MARK[mod], z=6)
        cv.marker(x_of(r.replication_pooled_g), yc, "o", size=2.6, face=INK, edge=INK, ew=0.2, z=7)
        vals = {"k": f"{int(r.replication_k)}", "series": f"{int(r.replication_series)}", "g": ci_txt[i],
                "fdr": fdr_txt[i]}
        for key, v in vals.items():
            cv.text(pos[key][1], yc, v, ha="right", role="value",
                    color=GREY if v == "not ranked" else INK)
    bottom = ytop + len(f) * P
    for v in ticks:
        if v < 0:
            cv.line([x_of(v)] * 2, [ytop, bottom], lw=0.3, color="#D0D0D0", linestyle=(0, (1, 1.6)), zorder=1)
    cv.line([x_of(0)] * 2, [ytop, bottom + 0.25], lw=LW_AXIS)
    x_axis(cv, x_of, ticks, bottom + 0.25, title="Pooled Hedges' g (95% CI)")


def draw_g_judgement(cv, x0, y0, g, ticks):
    """Admitted lower-output interventions in statistical order, drawn as in c-e."""
    rows = g[g["kept"]].sort_values("statistical_rank")
    panel_title(cv, x0, y0, "g", "Judgement ranking", "in old, otherwise healthy animals vs untreated controls")
    labels_w = max(width_mm(re.sub(r", ZT\d+\)$", ")", s), FS, "bold" if b else "normal")
                   for s, b in zip(rows["text"], rows["bold"]))
    bx = x0 + LAB + labels_w + 1.4
    avail = (x0 + CW - 0.6 - width_mm(f"{ticks[-1]:g}", FS) / 2) - bx
    sc = avail / ticks[-1]
    x_of = lambda v: bx + v * sc
    ytop = y0 + TITLE_H
    for i, r in enumerate(rows.itertuples()):
        yc = ytop + (i + 0.5) * P
        label = re.sub(r", ZT\d+\)$", ")", r.text)
        cv.text(x0 + LAB, yc, label, weight="bold" if r.bold else "normal", color=CAT[r.condition_category],
                role="label")
        bar(cv, bx, yc, abs(r.rank_bound) * sc, "lower", r.modality)
    bottom = ytop + len(rows) * P
    for v in ticks:
        if v > 0:
            cv.line([x_of(v)] * 2, [ytop, bottom], lw=0.3, color="#D0D0D0", linestyle=(0, (1, 1.6)), zorder=1)
    cv.line([bx, bx], [ytop, bottom + 0.25], lw=LW_AXIS)
    x_axis(cv, x_of, ticks, bottom + 0.25, title="|Hedges' g| bound")
    return bottom - y0


def check_overlaps(cv, a_box):
    """No two text boxes of the ranking layer meet; nothing but the letter and the key enters panel a."""
    cv.fig.canvas.draw()
    r = cv.fig.canvas.get_renderer()
    inv = cv.ax.transData.inverted()
    boxes = []
    for t in cv.texts:
        bb = t.get_window_extent(r)
        (x0, y1), (x1, y0) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        boxes.append((t, x0, y0, x1, y1))
    bad = []
    for i in range(len(boxes)):
        ti, a0, b0, a1, b1 = boxes[i]
        if a0 < 0 or a1 > W + 1e-6:
            bad.append(f"off page: {ti.get_text()!r}")
        if ti._role not in ("letter", "key") and a0 < a_box[2] and b0 < a_box[3]:
            bad.append(f"inside panel a: {ti.get_text()!r}")
        for j in range(i + 1, len(boxes)):
            tj, c0, d0, c1, d1 = boxes[j]
            if a0 < c1 - 0.05 and c0 < a1 - 0.05 and b0 < d1 - 0.05 and d0 < b1 - 0.05:
                if ti._role == "title" and tj._role == "title":
                    continue  # runs of one title line abut by construction
                bad.append(f"{ti.get_text()!r} meets {tj.get_text()!r}")
    if bad:
        raise SystemExit("text overlaps:\n  " + "\n  ".join(bad))


def place_panel_a(layer_pdf, H):
    doc = pymupdf.open(layer_pdf)
    page = doc[0]
    src = pymupdf.open(PANEL_A)
    sr = src[0].rect
    clip = pymupdf.Rect(0, A_CLIP_TOP_PT, sr.width, sr.height)
    pt = 72 / 25.4
    target = pymupdf.Rect(0, A_OFF * pt, sr.width * S_A, A_OFF * pt + (sr.height - A_CLIP_TOP_PT) * S_A)
    page.show_pdf_page(target, src, 0, clip=clip, overlay=True)
    a_bottom = A_OFF * pt + (sr.height - A_CLIP_TOP_PT) * S_A
    spans = [sp for b in page.get_text("dict")["blocks"] for ln in b.get("lines", []) for sp in ln["spans"]
             if sp["bbox"][3] < a_bottom and sp["text"].strip() in RELABEL]
    assert sorted(sp["text"].strip() for sp in spans) == sorted(RELABEL), [sp["text"] for sp in spans]
    for sp in spans:   # the two labels of panel a, reset in the same face, size and colour
        x0, y0, x1, y1 = sp["bbox"]
        new = RELABEL[sp["text"].strip()]
        width = pymupdf.Font(fontfile=ARIAL_BOLD).text_length(new, fontsize=sp["size"])
        page.draw_rect(pymupdf.Rect(x0 - 0.3, y0 - 0.2, max(x1, x0 + width) + 0.8, y1 + 0.2), color=None, fill=(1, 1, 1), overlay=True)
        c = sp["color"]
        page.insert_text((sp["origin"][0], sp["origin"][1]), new, fontname="arialbd", fontfile=ARIAL_BOLD,
                         fontsize=sp["size"], color=((c >> 16 & 255) / 255, (c >> 8 & 255) / 255, (c & 255) / 255))
    page.show_pdf_page(pymupdf.Rect(0, a_bottom + 1.0 * pt, W * pt, a_bottom + (1.0 + STRIP_H) * pt),
                       pymupdf.open(STRIP), 0, overlay=True)
    doc.set_metadata({"title": "Figure 3", "author": "Kejun Albert Ying", "creator": Path(__file__).name})
    doc.save(OUT_PDF, garbage=3, deflate=True)
    assert abs(page.rect.width - W * pt) < 0.5 and abs(page.rect.height - H * pt) < 0.5


def render_png():
    doc = pymupdf.open(OUT_PDF)
    doc[0].get_pixmap(dpi=200, alpha=False).save(OUT_PNG)


def records(t, be, f, g, checks, H, n_low, n_high):
    rows = []
    for cat, letter_s in zip(be, "cde"):
        for part in be[cat]:
            for r in part.itertuples():
                rows.append({"panel": letter_s, "bar": r.text, "direction": r.direction, "gse_id": r.gse_id,
                             "comparison_name": r.comparison_name,
                             "statistical_rank": int(r.statistical_rank_comparison),
                             "of": n_low if r.direction == "lower" else n_high,
                             "statistical_score": r.statistical_score, "reason": r.reason if r.direction == "lower"
                             else ""})
    for r in f.itertuples():
        rows.append({"panel": "f", "bar": r.text, "direction": "lower", "gse_id": "", "comparison_name": "",
                     "statistical_rank": "" if np.isnan(r.replication_best_single_rank)
                     else int(r.replication_best_single_rank), "of": n_low,
                     "statistical_score": r.replication_score, "reason": ""})
    for r in g[g["kept"]].itertuples():
        rows.append({"panel": "g", "bar": r.text, "direction": "lower", "gse_id": r.gse_id,
                     "comparison_name": r.comparison_name, "statistical_rank": int(r.statistical_rank), "of": n_low,
                     "statistical_score": abs(r.rank_bound), "reason": ""})
    ranks = pd.DataFrame(rows)
    ranks.to_csv(OUTPUT / "figure3_compact_statistical_ranks.csv", index=False)
    fr = ranks[ranks.panel == "f"]
    assert set(fr.loc[fr.statistical_rank == "", "bar"]) == set(
        f.loc[f.replication_best_single_rank.isna(), 'text'])
    source_rank = f.loc[f["text"].eq("PLX5622 †"), "replication_best_single_rank"].iloc[0]
    assert int(fr.loc[fr.bar == "PLX5622 †", "statistical_rank"].iloc[0]) == int(source_rank)
    # condition counts across every decrease drawn in b-e
    lows = pd.concat([lo for _, lo in be.values()])
    counts = {c: lows[c].value_counts().to_dict() for c in be_fig.STATE_COLUMNS}
    all_met = lows.loc[lows["reason"].eq(be_fig.ALL_MET), "text"].tolist()
    rec = {"script": str(Path(__file__).resolve()), "run_utc": datetime.now(timezone.utc).isoformat(),
           "page_mm": [W, round(H, 2)], "panel_a": {"path": str(PANEL_A), "sha256": sha256(PANEL_A),
                                                     "scale": S_A, "clip_top_pt": A_CLIP_TOP_PT},
           "table": {"path": str(TABLE), "sha256": sha256(TABLE)},
           "builder": {"path": str(HERE / "build_fig3_rankings_table_2026_11_08.py"),
                       "sha256": sha256(HERE / "build_fig3_rankings_table_2026_11_08.py")},
           "palette": P3.name(), "matplotlib": matplotlib.__version__,
           "builder_checks_passed": checks,
           "be_decreases_drawn": len(lows), "be_condition_counts": counts, "be_all_three_met": all_met,
           "g_rows": int(g["kept"].sum()), "g_admitted_comparisons": len(g)}
    (OUTPUT / "figure3_compact_build_record.json").write_text(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
