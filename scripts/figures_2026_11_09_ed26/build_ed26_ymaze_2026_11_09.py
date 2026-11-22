"""Extended Data Figure 26 with the Y-maze panels (g, h) added under A.M.'s panels a-f.

Source: evidence/ed26-ymaze-2026-11-09/Cecilia Ymaze Summary.xlsx (sheet Testing; received 9 November 2026). Group
labels: 'control' = control diet, 'exp' = PLX5622 (animal numbers match the NOR sheet of the source-data workbook).
Novel-arm preference = 100 x novel / (novel + familiar), recomputed from the arm values and checked against the
workbook's own columns. Two-sided exact Mann-Whitney test between groups. Panels a-f are page 33 of the v19 package,
placed as vector; panel d is redrawn in place (11 November 2026) to show the novel object visits without the outlying
control C64, as in Figure 7n: the control box, its points and the P value are rebuilt from the NOR and Statistics sheets
of the source-data workbook in the panel's own coordinates and style, and everything else on the page is unchanged.
"""
import sys, json, hashlib
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, openpyxl, pymupdf
from scipy import stats

ROOT = Path(__file__).resolve().parents[2] / "paper/revision_2026-09-07"
EV = ROOT / "evidence/ed26-ymaze-2026-11-09"
SRC = EV / "Cecilia Ymaze Summary.xlsx"
V19 = ROOT / "figures/03_Figures_WORKING_2026-10-10-v19.pdf"
for p in (SRC, V19):
    if not p.exists():
        raise SystemExit(f"Cannot proceed - required data file not found: {p}")
CTRL, PLX, CTRL_E, PLX_E, INK = "#B2F3B2", "#FFD0FF", "#3E8E3E", "#B04CA8", "#1A1A1A"
plt.rcParams.update({"pdf.fonttype": 42, "font.family": "Arial", "font.size": 5.5, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2})

ws = openpyxl.load_workbook(SRC, data_only=True)["Testing"]
rows = []
for r in ws.iter_rows(min_row=7, values_only=True):
    for off in (0, 16):
        if r[off] and str(r[off]).startswith("Trial"):
            rows.append({"id": str(r[off + 1]), "novel": r[off + 4], "group": r[off + 5], "dist": r[off + 6:off + 9],
                         "time": r[off + 9:off + 12], "wb_dist": r[off + 13], "wb_time": r[off + 14]})
assert sum(q["group"] == "control" for q in rows) == 9 and sum(q["group"] == "exp" for q in rows) == 10


SD = ROOT / "deliverables/FINAL_REVISION_2026-11-06/SourceData_Fig7_EDFig25_EDFig26.xlsx"
GREEN, GRID, BOXLINE = (0.698039, 0.952941, 0.698039), (0.870588, 0.870588, 0.870588), (0.2, 0.2, 0.2)


def v19_page_without_panel_d_control():
    """A copy of v19 page 33 with panel d's control box, its points and its P value removed (redaction)."""
    base = pymupdf.open(); base.insert_pdf(pymupdf.open(V19), from_page=32, to_page=32)
    v19page = base[0]
    old = [g for g in v19page.get_drawings() if 282 < g["rect"].x0 and g["rect"].x1 < 337 and g["rect"].y0 > 188]
    circles = [g["rect"] for g in old if g.get("fill") and abs(g["fill"][0] - GREEN[0]) < 1e-3
               and g.get("color") == (0.0, 0.0, 0.0) and g["rect"].width < 5]
    boxes = [g["rect"] for g in old if g.get("fill") and abs(g["fill"][0] - GREEN[0]) < 1e-3 and g["rect"].width > 10]
    assert len(circles) == 6 and len(boxes) == 1, (len(circles), len(boxes))
    y0, y150 = 292.1, 192.5  # value 0 and 150 (grid lines of the panel)
    ymap = lambda v: y0 - v * (y0 - y150) / 150
    wb = openpyxl.load_workbook(SD, data_only=True)
    nor = list(wb["NOR"].iter_rows(values_only=True))
    head = next(i for i, r in enumerate(nor) if r and r[0] == "Mouse")
    cols = {name: j for j, name in enumerate(nor[head]) if name}
    mice = [r for r in nor[head + 1:] if r and r[0]]
    keep = [r for r in mice if r[cols["Group"]] == "Control" and str(r[cols["Outlier for visits and exploration"]] or "") != "Yes"]
    dropped = [r[0] for r in mice if r[cols["Group"]] == "Control" and r not in keep]
    assert dropped == ["C64"], dropped
    ctrl = np.array([float(r[cols["Novel object visits (n)"]]) for r in keep])
    stats_rows = list(wb["Statistics"].iter_rows(values_only=True))
    sh = next(i for i, r in enumerate(stats_rows) if r and "Panel" in r)
    sc = {name: j for j, name in enumerate(stats_rows[sh]) if name}
    row = next(r for r in stats_rows[sh + 1:] if r and r[sc["Panel"]] == "ED Fig. 26d")
    assert row[sc["n"]] == "8 vs 5" and "C64 excluded" in str(row[sc["Note"]]), row
    p_d = float(row[sc["P"]])
    # Match each plotted control point to its value by height; the one at 127 visits is C64.
    centre = lambda r: ((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)
    values_of = {i: (y0 - centre(r)[1]) * 150 / (y0 - y150) for i, r in enumerate(circles)}
    outlier = max(values_of, key=values_of.get)
    assert abs(values_of[outlier] - 127) < 1.5, values_of[outlier]
    kept = [circles[i] for i in values_of if i != outlier]
    assert sorted(round(values_of[i]) for i in values_of if i != outlier) == sorted(round(v) for v in ctrl)
    box = boxes[0]
    ptext = [r for r in v19page.search_for("0.22") if 300 < r.x0 < 320 and 190 < r.y0 < 205]
    assert len(ptext) == 1, ptext
    # Remove only objects lying wholly inside these areas: the control box, whiskers and points, the C64 point and
    # the P value; grid lines, axes and the PLX5622 box extend beyond them and stay.
    for r in (pymupdf.Rect(box.x0 - 1.5, 262.0, box.x1 + 1.5, 291.7), pymupdf.Rect(circles[outlier]) + (-1.2, -1.2, 1.2, 1.2),
              ptext[0] + (-0.5, -0.5, 0.5, 0.5)):
        v19page.add_redact_annot(r, fill=False)
    v19page.apply_redactions(images=pymupdf.PDF_REDACT_IMAGE_NONE, graphics=pymupdf.PDF_REDACT_LINE_ART_REMOVE_IF_COVERED,
                             text=pymupdf.PDF_REDACT_TEXT_REMOVE)
    assert "0.22" in v19page.get_text()  # the CA1 value of panel a stays
    assert not [r for r in v19page.search_for("0.22") if 300 < r.x0 < 320]
    return base, {"ctrl": ctrl, "kept": kept, "box": box, "p": p_d, "dropped": dropped, "ptext": ptext[0],
                  "ymap": ymap}


def redraw_panel_d(page, d):
    """Novel object visits without the outlying control (C64), drawn in panel d's own coordinates and style."""
    ctrl, kept, box, p_d, ymap = d["ctrl"], d["kept"], d["box"], d["p"], d["ymap"]
    centre = lambda r: ((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)
    xc = (box.x0 + box.x1) / 2
    q1, med, q3 = np.percentile(ctrl, [25, 50, 75])  # quantile type 7, as ggplot
    iqr = q3 - q1
    upper = ctrl[ctrl <= q3 + 1.5 * iqr].max(); lower = ctrl[ctrl >= q1 - 1.5 * iqr].min()
    if upper > q3:
        page.draw_line((xc, ymap(q3)), (xc, ymap(upper)), color=BOXLINE, width=0.995846)
    if lower < q1:
        page.draw_line((xc, ymap(q1)), (xc, ymap(lower)), color=BOXLINE, width=0.995846)
    page.draw_rect(pymupdf.Rect(box.x0, ymap(q3), box.x1, ymap(q1)), color=BOXLINE, fill=GREEN, width=0.995846)
    page.draw_line((box.x0, ymap(med)), (box.x1, ymap(med)), color=BOXLINE, width=1.991693)
    for r in kept:
        page.draw_circle(centre(r), r.width / 2, color=(0, 0, 0), fill=GREEN, width=0.661417)
    label = f"{p_d:.2f}"
    width = pymupdf.get_text_length(label, fontname="helv", fontsize=6)
    t = d["ptext"]
    page.insert_text(((t.x0 + t.x1) / 2 - width / 2, t.y0 + 0.905 * 6), label, fontname="helv", fontsize=6)
    return {"control_n": len(ctrl), "control_median": float(med), "p": p_d, "excluded": d["dropped"]}


def pref(q, k):
    left, right = q[k][1], q[k][2]
    nov, fam = (left, right) if q["novel"] == "Left" else (right, left)
    return 100 * nov / (nov + fam)


out = {}
fig = plt.figure(figsize=(510.2 / 72, 112 / 72))
for i, (k, wbk, lab, ylabel) in enumerate((("time", "wb_time", "g", "Time in novel arm, %"),
                                             ("dist", "wb_dist", "h", "Distance in novel arm, %"))):
    for q in rows:
        assert abs(pref(q, k) - q[wbk]) < 1e-6
    c = np.array([pref(q, k) for q in rows if q["group"] == "control"])
    e = np.array([pref(q, k) for q in rows if q["group"] == "exp"])
    p = stats.mannwhitneyu(e, c, alternative="two-sided", method="exact").pvalue
    out[k] = {"control_n": len(c), "plx_n": len(e), "control_median": float(np.median(c)), "plx_median": float(np.median(e)),
              "p_between": float(p), "p_control_vs_50": float(stats.wilcoxon(c - 50).pvalue),
              "p_plx_vs_50": float(stats.wilcoxon(e - 50).pvalue)}
    ax = fig.add_axes([0.085 + i * 0.2, 0.2, 0.105, 0.62])
    rng = np.random.default_rng(7 + i)
    for pos, v, fc, ec in ((0, c, CTRL, CTRL_E), (1, e, PLX, PLX_E)):
        bp = ax.boxplot([v], positions=[pos], widths=0.62, patch_artist=True, showfliers=False,
                        medianprops={"color": INK, "lw": 0.9}, whiskerprops={"lw": 0.5, "color": INK},
                        capprops={"lw": 0}, boxprops={"lw": 0.5, "edgecolor": INK, "facecolor": fc})
        ax.scatter(pos + rng.uniform(-0.17, 0.17, len(v)), v, s=7, facecolor=fc, edgecolor=ec, linewidth=0.7, zorder=3,
                   clip_on=False)
    ax.axhline(50, lw=0.5, color="#888888", ls=(0, (2, 2)), zorder=0)
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, 118); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks([0, 1], ["Control", "PLX5622"])
    ax.set_ylabel(ylabel, fontweight="bold", fontsize=5.5)
    ax.plot([0, 0, 1, 1], [106, 109, 109, 106], lw=0.5, color=INK)
    ax.text(0.5, 110.5, f"{p:.2f}", ha="center", va="bottom", fontsize=5.5)
    fig.text(0.085 + i * 0.2 - 0.065, 0.93, lab, fontsize=8, fontweight="bold", va="top")
    if i == 0:
        fig.text(0.085 - 0.035, 0.93, "Y-maze novel-arm test", fontsize=6.5, va="top")
layer = EV / "_ymaze_row.pdf"
fig.savefig(layer); plt.close(fig)
doc = pymupdf.open(); page = doc.new_page(width=510.2, height=322 + 112)
base, panel_d = v19_page_without_panel_d_control()
page.show_pdf_page(pymupdf.Rect(0, 0, 510.2, 322), base, 0)
out_d = redraw_panel_d(page, panel_d)
page.show_pdf_page(pymupdf.Rect(0, 322, 510.2, 434), pymupdf.open(layer), 0)
doc.save(EV / "ed26.pdf", garbage=3, deflate=True); layer.unlink()
pymupdf.open(EV / "ed26.pdf")[0].get_pixmap(dpi=200).save(EV / "ed26.png")
out["total_distance_cm"] = {"control_median": float(np.median([sum(q["dist"]) for q in rows if q["group"] == "control"])),
                            "plx_median": float(np.median([sum(q["dist"]) for q in rows if q["group"] == "exp"]))}
out["source"] = {"path": str(SRC), "sha256": hashlib.sha256(SRC.read_bytes()).hexdigest()}
out["panel_d"] = out_d
out["panel_d_source"] = {"path": str(SD), "sha256": hashlib.sha256(SD.read_bytes()).hexdigest()}
(EV / "numbers.json").write_text(json.dumps(out, indent=1))
print(json.dumps({k: out[k] for k in ("time", "dist")}, indent=1))
