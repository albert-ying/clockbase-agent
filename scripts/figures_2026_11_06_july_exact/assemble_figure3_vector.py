#!/usr/bin/env python3
"""
Assemble Figure 3 as a TRUE VECTOR PDF from the individual vector panels.

Background
----------
figure_rev1/Figure3.pdf was a single flattened raster (0 fonts, 1 image) -- even
the axis text was rasterised by whatever external layout tool produced it. Both
panels exist as clean vector PDFs (the volcano carries a raster scatter layer,
which is inherent to a ~100k-point plot, but all of its text/axes/labels are
vector; the bar grid is fully vector). This script re-composites the exact vector
panels so the result is vector (selectable text, infinite zoom on everything
except the volcano point cloud), mirroring assemble_figure4_vector.py.

Panel -> source PDF (both already carry the re-synced numbers, 09:16 build):
  a    combined volcano        Aging-AI-Agent/figures/volcano_combined_modalities.pdf
  b-e  top-15 bar grid         Aging-AI-Agent/figures/top_hits_barplot_combined_modalities_v2.pdf
       (b Genetic, c Drug, d Disease, e Environment -- 2x2; letters added here)

Layout: volcano full width on top, bar grid full width below (matches the
published Figure3.png stack, ratio ~1.30).

Output overwrites the canonical filenames:
  figure/Figure3.pdf, figure/Figure3.png (600 dpi raster for the manuscript .qmd)
  figure_rev1/Figure3.pdf, figure_rev1/Figure3.png
  submission_revision1_2026-05-23/Figures/Figure3.pdf
"""
import os, shutil, sys
import fitz  # PyMuPDF

ROOT = '/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base/clockbase_revision_local_2026-11-03/clockbase-agent/paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact'
FIGSRC = os.path.join(ROOT, "figures")

VOLCANO = os.path.join(FIGSRC, "volcano_combined_modalities.pdf")
BARS = os.path.join(FIGSRC, "top_hits_barplot_combined_modalities_v2.pdf")

JUDGEMENT = os.path.join(FIGSRC, "judgement_barplot_combined_modalities.pdf")
REPLICATION = os.path.join(FIGSRC, "replication_barplot_combined_modalities.pdf")

# Page geometry (7 November 2026): Nature double-column width, panels drawn at print size and placed at 100%.
# b-e is 106.7 mm high (two key lines); f and g are 38.2 mm high, with the ranking header in their top line on the
# baseline of the panel letter, and g is 99.6 mm wide to hold its header at 8 pt, so f is 80.4 mm.
MM = 72 / 25.4
PAGE_W = 183 * MM
LETTER_PT = 10
LETTER_BASELINE_PT = LETTER_PT * 0.78   # below the top of the panel
GAP_MM = 3
BE_H, FG_H, F_W, G_W = 106.7, 38.2, 80.4, 99.6
FG_Y = 89 + BE_H + GAP_MM
# (source, (x, y, width, height) in mm)
LAYOUT = {"a": (VOLCANO, (0, 0, 183, 86)), "b_e": (BARS, (0, 89, 183, BE_H)),
          "f": (REPLICATION, (0, FG_Y, F_W, FG_H)), "g": (JUDGEMENT, (F_W + GAP_MM, FG_Y, G_W, FG_H))}
PAGE_H = (FG_Y + FG_H) * MM
assert abs(F_W + GAP_MM + G_W - 183) < 1e-9 and PAGE_H <= 247 * MM
HEADERS = {"f": "Replication ranking", "g": "Judgement ranking"}


def build(out_pdf):
    for p in (VOLCANO, BARS, REPLICATION, JUDGEMENT):
        if not os.path.exists(p):
            sys.exit(f"Cannot proceed - required panel PDF not found: {p}")
    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=PAGE_H)
    for key, (path, (x, y, w, h)) in LAYOUT.items():
        with fitz.open(path) as s:
            r = s[0].rect
            assert abs(r.width - w * MM) < 1.5 and abs(r.height - h * MM) < 1.5, (key, r.width / MM, r.height / MM, w, h)
            page.show_pdf_page(fitz.Rect(x * MM, y * MM, x * MM + r.width, y * MM + r.height), s, 0)
            if key in HEADERS:
                # the ranking header sits on the baseline of the panel letter, in the panel's top line
                hits = [sp for blk in s[0].get_text("dict")["blocks"] if blk["type"] == 0 for ln in blk["lines"]
                        for sp in ln["spans"] if sp["text"].strip() == HEADERS[key]]
                assert len(hits) == 1 and abs(hits[0]["origin"][1] - LETTER_BASELINE_PT) < 0.05, (key, hits)
    # letters: a, f, g at their panel corners; b-e level with the subplot titles found in the grid PDF
    with fitz.open(BARS) as bars_doc:
        spans = [sp for blk in bars_doc[0].get_text("dict")["blocks"] if blk["type"] == 0
                 for line in blk["lines"] for sp in line["spans"]]
    heads = {}
    for t in ("Genetic", "Drug", "Disease", "Environment"):  # exact subplot titles, not tag text
        hits = [fitz.Rect(sp["bbox"]) for sp in spans if sp["text"].strip() == t]
        assert len(hits) == 1, (t, len(hits))
        heads[t] = hits[0]
    gx, gy = LAYOUT["b_e"][1][0] * MM, LAYOUT["b_e"][1][1] * MM
    half = PAGE_W / 2
    letters = [("a", 0.4, LETTER_BASELINE_PT),
               ("f", 0.4, LAYOUT["f"][1][1] * MM + LETTER_BASELINE_PT),
               ("g", LAYOUT["g"][1][0] * MM + 0.4, LAYOUT["g"][1][1] * MM + LETTER_BASELINE_PT)]
    for letter, title in zip("bcde", ("Genetic", "Drug", "Disease", "Environment")):
        x = gx + 0.4 if heads[title].x0 < half - gx else half + 0.4
        letters.append((letter, x, gy + heads[title].y1))
    for lab, lx, ly in letters:
        page.insert_text((lx, ly), lab, fontsize=LETTER_PT, fontname="hebo", color=(0, 0, 0))
    doc.save(out_pdf, garbage=4, deflate=True)
    doc.close()
    return PAGE_W, PAGE_H

def render_png(pdf_path, png_path, dpi=600):
    doc = fitz.open(pdf_path)
    pix = doc[0].get_pixmap(dpi=dpi, alpha=False)
    pix.save(png_path)
    doc.close()

DESTS = [
    (os.path.join(ROOT, "figure", "Figure3.pdf"), os.path.join(ROOT, "figure", "Figure3.png")),
    (os.path.join(ROOT, "figure_rev1", "Figure3.pdf"), os.path.join(ROOT, "figure_rev1", "Figure3.png")),
    (os.path.join(ROOT, "submission_revision1_2026-05-23", "Figures", "Figure3.pdf"), None),
]
BACKUP_DIR = os.path.join(ROOT, "figure", "outdated")

def deploy():
    master = "/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base/clockbase_revision_local_2026-11-03/clockbase-agent/paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact/assembly_build_3/Figure3_vector.pdf"
    os.makedirs("/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base/clockbase_revision_local_2026-11-03/clockbase-agent/paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact/assembly_build_3", exist_ok=True)
    w, h = build(master)
    print(f"built vector master: {w:.0f}x{h:.0f}pt ratio {h/w:.3f}")
    os.makedirs(BACKUP_DIR, exist_ok=True)
    for pdf_dst, png_dst in DESTS:
        for f in (pdf_dst, png_dst):
            if f and os.path.exists(f):
                bak = os.path.join(BACKUP_DIR, os.path.basename(os.path.dirname(f)) + "_" +
                                   os.path.basename(f).replace(".", "_raster_2026-06-15."))
                if not os.path.exists(bak):
                    shutil.copy2(f, bak)
        shutil.copy2(master, pdf_dst)
        if png_dst:
            render_png(pdf_dst, png_dst, dpi=600)
        print("wrote:", pdf_dst, "+ png" if png_dst else "")
    print("Backups in:", BACKUP_DIR)

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] != "--deploy":
        out = sys.argv[1]
        os.makedirs(os.path.dirname(out), exist_ok=True)
        build(out)
        render_png(out, out.replace(".pdf", ".png"), dpi=150)
        print("Built (preview):", out)
    else:
        deploy()
