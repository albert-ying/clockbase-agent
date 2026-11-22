"""Figure 2a without the mouse DNA methylation bar, and Extended Data Figure 16 without the n in its annotation.

Both pages of the v19 package are raster panels. Figure 2a is redrawn from the three sample counts that the legend and
the Results give (852,381; 964,083; 230,516) and placed over the old panel; the annotation of Extended Data Figure 16
is reset without its count. Everything else on the two pages is unchanged.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pymupdf

ROOT = Path(__file__).resolve().parents[2] / "paper/revision_2026-09-07"
OUT = ROOT / "evidence/figure-numbers-2026-11-09"
V19 = pymupdf.open(ROOT / "figures/03_Figures_WORKING_2026-10-10-v19.pdf")
COUNTS = [("RNA human", 852381), ("RNA mouse", 964083), ("DNAm human", 230516)]
ORANGE = "#F28E1C"

# ---- Figure 2a
box = pymupdf.Rect(86, 126, 240, 254)
plt.rcParams.update({"pdf.fonttype": 42, "font.size": 6})
fig = plt.figure(figsize=(box.width / 72, box.height / 72))
ax = fig.add_axes([0.37, 0.25, 0.42, 0.68])
ys = range(len(COUNTS))[::-1]
ax.barh(list(ys), [c for _, c in COUNTS], color=ORANGE, height=0.62)
for y, (_, c) in zip(ys, COUNTS):
    ax.text(c * 1.15, y, f"{c:,}", va="center", ha="left", color=ORANGE, fontsize=6, clip_on=False)
ax.set_xscale("log"); ax.set_xlim(1e3, 1.2e6); ax.set_xticks([1e4, 1e5, 1e6])
ax.set_yticks(list(ys), [n for n, _ in COUNTS]); ax.set_ylim(-0.6, len(COUNTS) - 0.4)
ax.set_xlabel("Number of Samples (log scale)", fontsize=6.5)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.tick_params(length=2, width=0.5, which="both"); ax.minorticks_off()
layer = OUT / "_fig2a.pdf"; fig.savefig(layer, transparent=False); plt.close(fig)
doc = pymupdf.open(); doc.insert_pdf(V19, from_page=1, to_page=1); page = doc[0]
page.draw_rect(pymupdf.Rect(78, 122, 240, 256), color=None, fill=(1, 1, 1), overlay=True)
page.show_pdf_page(box, pymupdf.open(layer), 0, overlay=True)
doc.save(OUT / "figure2.pdf", garbage=3, deflate=True); layer.unlink()
doc[0].get_pixmap(dpi=140, clip=pymupdf.Rect(70, 95, 420, 270)).save(OUT / "figure2a_check.png")

# ---- Extended Data Figure 16
doc = pymupdf.open(); doc.insert_pdf(V19, from_page=22, to_page=22); page = doc[0]
box = pymupdf.Rect(428, 360, 534, 396)
page.draw_rect(box, color=None, fill=(1, 1, 1), overlay=True)
left = page.insert_textbox(pymupdf.Rect(428, 362, 534, 396),
                           "Pooled three-modality atlas\n(mouse RNA-seq, human RNA-seq,\nhuman DNA methylation)\nx-axis clipped to ±16",
                           fontname="helv", fontsize=5.0, align=pymupdf.TEXT_ALIGN_CENTER, color=(0.2, 0.2, 0.2))
assert left >= 0, left
doc.save(OUT / "ed16.pdf", garbage=3, deflate=True)
doc[0].get_pixmap(dpi=140, clip=pymupdf.Rect(380, 280, 560, 420)).save(OUT / "ed16_check.png")
print("ok")
