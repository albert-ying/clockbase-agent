"""Figure 3 method strip: one pool of verified comparisons, three rankings (native vector SVG, no raster).

Layout after ImageGen concept A (evidence/figure3-schematic-2026-11-09/concept_A.png; no generated pixels used).
Mouse icon: Servier Medical Art via Bioicons (CC BY 3.0), clip rectangles removed. Ranks printed beside the nominees
are read from fig3_rankings.csv and the judgement table; the build stops if they differ from the values below.
  python3 build_figure3_schematic_2026_11_09.py
"""
import os
import hashlib, re, csv
from pathlib import Path
try:
    import cairosvg
except ModuleNotFoundError:
    import sys
    # Reuse the locally installed pure-Python renderer with the pipeline CFFI backend.
    local = Path('/opt/homebrew/lib/python3.11/site-packages')
    if not (local / 'cairosvg/__init__.py').is_file():
        raise
    sys.path.append(str(local))
    import cairosvg

HERE = Path(__file__).resolve().parent
EV = HERE.parents[1] / "paper/revision_2026-09-07/evidence"
OUT = EV / "figure-fixes-2026-11-17/figure3"
ICON = Path.home() / "Library/Application Support/rpkg_icon/bioicons/cc-by-3.0/Animals/Servier/mouse-gray.svg"
INK, GREY, HAIR = "#181114", "#6B6560", "#B8B2AD"
LOW, LOWT, HIGH, HIGHT, NEU = "#0072B2", "#CFE3F1", "#D55E00", "#F6D5BE", "#D9D5D1"
W, H = 183.0, 39.5
PT = 25.4 / 72

# ---- source check of the three nominees
rows = list(csv.DictReader(open(EV / "figures-rebuild-2026-11-06/figure3/fig3_rankings.csv")))
def one(pred):
    r = [x for x in rows if pred(x)]
    assert len(r) == 1, len(r)
    return r[0]
cold = one(lambda x: x["gse_id"] == "GSE201135" and x["label"].startswith("Cold exposure"))
with open(Path(os.environ.get('CLOCKBASE_JUDGEMENT') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'followup-judgment-001') / 'ranked_with_judgement.csv') as fh:
    ranked = list(csv.DictReader(fh))
cold_source = next(r for r in ranked if r['gse_id'] == cold['gse_id']
                   and r['comparison_name'] == cold['comparison_name'])
assert int(float(cold['statistical_rank_comparison'])) == int(cold_source['statistical_rank'])
plx = one(lambda x: (x["printed_f"] or "").startswith("PLX5622"))
with open(Path(os.environ.get('CLOCKBASE_META') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'intervention-meta-2026-10-17') / 'primary/ranked_replicated_age-lowering.csv') as fh:
    replicated = list(csv.DictReader(fh))
plx_rank = next(i for i, r in enumerate(replicated, start=1) if r['canonical_target'] == 'plx5622')
assert int(float(plx['replication_rank_unit'])) == plx_rank
oua = one(lambda x: (x["printed_g"] or "").startswith("Ouabain"))
with open(Path(os.environ.get('CLOCKBASE_JUDGEMENT') or Path(__import__('os').environ['CLOCKBASE_VERIFIER_ROOT']) / 'followup-judgment-001') / 'judgement_age-lowering.csv') as fh:
    admitted = sorted(csv.DictReader(fh), key=lambda r: int(r['statistical_rank']))
ouabain_rank = next(i for i, r in enumerate(admitted, start=1) if r['gse_id'] == 'GSE122080')
assert int(float(oua["judgement_rank_comparison"])) == ouabain_rank

icon = ICON.read_text()
paths = re.sub(r'\s*clip-path="[^"]*"', "", "".join(re.findall(r"<path[^>]*/>", icon.split("</defs>")[1])))
vb = [float(v) for v in re.search(r'viewBox="([^"]+)"', icon).group(1).split()]

E = []
def text(x, y, s, size=7, weight="normal", fill=INK, anchor="start", style="normal"):
    E.append(f'<text x="{x:.2f}" y="{y:.2f}" font-family="Arial" font-size="{size*PT:.3f}" font-weight="{weight}" '
             f'font-style="{style}" fill="{fill}" text-anchor="{anchor}">{s}</text>')
def line(pts, col=INK, lw=0.2, dash=None):
    d = "M" + " L".join(f"{x:.2f} {y:.2f}" for x, y in pts)
    E.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="{lw}"' + (f' stroke-dasharray="{dash}"' if dash else "") + "/>")
def arrow(pts, col=INK, lw=0.2):
    line(pts, col, lw)
    (x0, y0), (x1, y1) = pts[-2], pts[-1]
    dx, dy = x1 - x0, y1 - y0
    n = (dx * dx + dy * dy) ** 0.5
    ux, uy = dx / n, dy / n
    a, b = 1.1, 0.42
    E.append(f'<path d="M{x1:.2f} {y1:.2f} L{x1-a*ux-b*uy:.2f} {y1-a*uy+b*ux:.2f} L{x1-a*ux+b*uy:.2f} {y1-a*uy-b*ux:.2f}Z" fill="{col}"/>')
def dot(x, y, r, fill, stroke="none"):
    E.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="0.12"/>')
def rect(x, y, w, h, fill, stroke=HAIR, rx=0.5, lw=0.18):
    E.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{lw}"/>')
def mouse(x, y, h, opacity=1.0):
    s = h / vb[3]
    E.append(f'<g transform="translate({x:.2f} {y:.2f}) scale({s:.5f})" opacity="{opacity}">{paths}</g>')
    return vb[2] * s
def group(x, y, col, n=6):
    for i in range(n):
        dot(x + (i % 3) * 1.25, y + (i // 3) * 1.25, 0.45, col)
def ranked(x, yc, hi):
    """three stacked bars; the bar of the nominee's place is filled"""
    for i, wd in enumerate((9.0, 7.0, 5.2)):
        y = yc - 2.55 + i * 1.75
        rect(x, y, wd, 1.25, "#7FB4D6", stroke="none", rx=0.2)

# ---- pool
import random
rnd = random.Random(3)
cx, cy, r = 13.5, 21.5, 8.3
E.append(f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="#F4F2F0" stroke="{HAIR}" stroke-width="0.18"/>')
placed = []
while len(placed) < 46:
    x, y = cx + rnd.uniform(-r, r), cy + rnd.uniform(-r, r)
    if (x - cx) ** 2 + (y - cy) ** 2 < (r - 1.0) ** 2 and all((x - a) ** 2 + (y - b) ** 2 > 1.5 for a, b in placed):
        placed.append((x, y))
for i, (x, y) in enumerate(placed):
    dot(x, y, 0.5, LOW if i % 9 == 0 else HIGH if i % 4 == 0 else "#A9A29C")
text(cx, 5.2, "Verified", 8, "bold", anchor="middle"); text(cx, 8.5, "comparisons", 8, "bold", anchor="middle")
text(0.2, 3.4, "b", 10, "bold", fill="black")

LANES = [(7.5, "Statistical ranking", "Effect size in one comparison", "c–e"),
         (19.5, "Replication ranking", "Agreement across independent studies", "f"),
         (32.0, "Judgement ranking", "AI models mimic expert judgement", "g")]
XT, XP, XL, XN, XM = 29.0, 81.0, 124.0, 137.0, 163.0
for yc, title, what, panel in LANES:
    arrow([(cx + r + 0.6, cy + (yc - cy) * 0.35), (23.5, yc), (26.6, yc)])
    text(XT, yc - 0.5, title, 8, "bold")
    text(XT, yc + 3.0, what, 7, fill=GREY)
# lane 1: treated against control in one study
y = LANES[0][0]
rect(XP, y - 4.2, 31, 8.4, "#F7F6F4")
group(XP + 3.4, y - 0.1, LOW); group(XP + 23.6, y - 0.1, "#A9A29C")
line([(XP + 8.4, y + 0.55), (XP + 22.2, y + 0.55)], INK, 0.2); line([(XP + 8.4, y - 0.2), (XP + 8.4, y + 1.3)]); line([(XP + 22.2, y - 0.2), (XP + 22.2, y + 1.3)])
text(XP + 15.3, y - 0.3, "difference", 6, fill=GREY, anchor="middle")
text(XP + 4.7, y - 1.7, "treated", 6, fill=GREY, anchor="middle"); text(XP + 24.9, y - 1.7, "control", 6, fill=GREY, anchor="middle")
# lane 2: three studies pooled
y = LANES[1][0]
for k in (-1, 0, 1):
    yy = y + k * 3.1
    rect(XP, yy - 1.3, 14.0, 2.6, "#F7F6F4")
    for i in range(3):
        dot(XP + 2.0 + i * 1.3, yy, 0.45, LOW); dot(XP + 9.2 + i * 1.3, yy, 0.45, "#A9A29C")
    line([(XP + 14.0, yy), (XP + 22.0, y)], INK, 0.18)
E.append(f'<path d="M{XP+22.0} {y} l2.6 -2.0 l2.6 2.0 l-2.6 2.0Z" fill="{LOW}" stroke="{INK}" stroke-width="0.15"/>')
text(XP + 24.6, y + 4.6, "pooled", 6, fill=GREY, anchor="middle")
# lane 3: two AI models read the study record and answer the questions an expert would ask
y = LANES[2][0]
rect(XP, y - 3.0, 4.2, 6.0, "#FFFFFF", stroke=GREY, rx=0.3)
for i in range(4):
    line([(XP + 0.7, y - 1.8 + i * 1.2), (XP + 3.5, y - 1.8 + i * 1.2)], HAIR, 0.18)
for k in (-1, 1):
    rect(XP + 6.4, y + k * 1.75 - 1.15, 9.6, 2.3, NEU, stroke=GREY, rx=0.4)
    text(XP + 11.2, y + k * 1.75 + 0.75, "AI model", 6, anchor="middle")
    line([(XP + 4.2, y), (XP + 6.4, y + k * 1.75)], INK, 0.18)
    line([(XP + 16.0, y + k * 1.75), (XP + 17.6, y)], INK, 0.18)
for k, q in zip((-1, 0, 1), ("old at start", "otherwise healthy", "untreated control")):
    yy = y + k * 2.3
    rect(XP + 18.6, yy - 0.75, 1.5, 1.5, "#FFFFFF", stroke=GREY, rx=0.2, lw=0.15)
    line([(XP + 18.9, yy + 0.05), (XP + 19.25, yy + 0.45), (XP + 19.85, yy - 0.45)], LOW, 0.28)
    text(XP + 20.8, yy + 0.75, q, 6, fill=GREY)
for (yc, _, _, panel), hi, nom in zip(LANES, (2, 2, 0), ("Cold exposure", "PLX5622", "Ouabain")):
    if panel == "g":
        # The admitted comparison proceeds to the nominee without a score item.
        arrow([(XP + 39.6, yc), (XN - 0.8, yc)])
    else:
        arrow([({LANES[0][0]: XP + 33.5, LANES[1][0]: XP + 28.4}[yc], yc), (XL - 0.8, yc)])
        ranked(XL, yc, hi)
        text(XL, yc - 3.5, {"c–e": "Effect size", "f": "Pooled effect"}[panel], 6, fill=GREY)
    text(XN, yc - 0.4, nom + " †", 6, "bold")
    text(XN, yc + 2.2, f"panel{'s' if '–' in panel else ''} {panel}", 6, fill=GREY)
    line([(XN + 19.5, yc), (XM - 4.0, yc), (XM - 4.0, LANES[1][0]), (XM - 0.6, LANES[1][0])] if yc != LANES[1][0]
          else [(XN + 19.5, yc), (XM - 0.6, yc)])
wm = mouse(XM + 1.0, LANES[1][0] - 6.4, 9.0)
text(XM + 1.0 + wm / 2, LANES[1][0] + 5.6, "Tested in", 7, "bold", anchor="middle")
text(XM + 1.0 + wm / 2, LANES[1][0] + 8.5, "aged mice", 7, "bold", anchor="middle")

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H}mm" viewBox="0 0 {W} {H}">'
       f'<rect width="{W}" height="{H}" fill="white"/>' + "".join(E) + "</svg>")
(OUT / "figure3_schematic.svg").write_text(svg)
cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(OUT / "figure3_schematic.pdf"))
cairosvg.svg2png(bytestring=svg.encode(), write_to=str(OUT / "figure3_schematic.png"), dpi=400)
with open(OUT / "asset_licences.csv", "w") as fh:
    fh.write("asset,source,licence,sha256\n")
    fh.write(f"mouse-gray.svg,Servier Medical Art via Bioicons,CC BY 3.0,{hashlib.sha256(ICON.read_bytes()).hexdigest()}\n")
print(OUT / "figure3_schematic.pdf")
