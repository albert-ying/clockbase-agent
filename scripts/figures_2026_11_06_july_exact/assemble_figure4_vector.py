#!/usr/bin/env python3
"""
Assemble Figure 4 as a TRUE VECTOR PDF from the individual matplotlib panel PDFs.

Background
----------
The previously delivered figure/Figure4.pdf was a single flattened 600-dpi JPEG
(0 fonts, 1 image) because the external layout tool rasterised the composite.
Every panel, however, exists as a clean vector PDF. This script re-composites the
*exact* vector panels that match the published figure (preserving every number),
so the result is fully vector (selectable text, infinite zoom).

Panel -> source PDF (6 November 2026, night: every panel drawn at its print size on the verified atlas and placed at
100% on a 183 mm page; numbers in the panels come from their generators' verified inputs):
  a  radar (combined)            create_radar_plot_combined.pdf          64,054 rows / 10,243 significant
  b  ridge (combined)            create_ridge_plot_combined.pdf          10,243 significant rows
  c  genetic-type barplot        analysis_2_genetic_types_barplot.pdf    knockout n = 657, baseline 25.0%
  d  genetic-type density        analysis_2_genetic_types_density.pdf
  e  pathway ORA                 analysis_3_pathway_enrichment.pdf       12 KEGG pathways, FDR < 0.05
  f  clinical-stage barplot      analysis_4_drug_clinical_barplot.pdf    FDA approved 59 lower / 102 higher
     FDA-approved drugs inset    analysis_4_fda_antiaging_drugs.pdf      15 drugs with the largest decreases
  g  GenAge waffle (mouse)       geneage_waffle_nature.pdf               62 squares, 39 agree, P = 0.028
  h  AgingAgentBench schematic  benchmark_schematic.pdf                 580 effects, 93 papers, 125 series
  i  published calls matched     benchmark_identification_bars.pdf       one column per published call, rows of j
  j  balanced accuracy           benchmark_accuracy_bars.pdf             mean of i's columns, 95% intervals
The DrugAge waffle is no longer placed; its six matched compounds are reported in the Results and Table S21.

Output overwrites the canonical filenames so downstream picks them up:
  figure/Figure4.pdf, figure/Figure4.png (600 dpi raster for the manuscript .qmd)
  figure_rev1/Figure4.pdf, figure_rev1/Figure4.png
  submission_revision1_2026-05-23/Figures/Figure4.pdf
"""
import os, shutil, sys
import fitz  # PyMuPDF

ROOT = '/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base/clockbase_revision_local_2026-11-03/clockbase-agent/paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact'
FIGSRC = os.path.join(ROOT, "figures")

# Page geometry (6 November 2026, night): Nature double-column width, every panel drawn at its print size and placed at
# 100% on a millimetre grid.
MM = 72 / 25.4
PAGE_W, PAGE_H = 183 * MM, 234 * MM
LETTER_PT = 10

def src(name):
    p = name if os.path.isabs(name) else os.path.join(FIGSRC, name)
    if not os.path.exists(p):
        sys.exit(f"Cannot proceed - required panel PDF not found: {p}")
    return p

# (source pdf, (x, y, width, height) in mm from the top-left, panel letter)
PANELS = [
    ("create_radar_plot_combined.pdf",         (  0,   0,  56, 50), "a"),
    ("create_ridge_plot_combined.pdf",         ( 59,   0, 124, 50), "b"),
    ("analysis_2_genetic_types_barplot.pdf",   (  0,  53,  48, 44), "c"),
    ("analysis_2_genetic_types_density.pdf",   ( 51,  53,  46, 44), "d"),
    ("analysis_3_pathway_enrichment.pdf",      (100,  53,  83, 44), "e"),
    ("analysis_4_drug_clinical_barplot.pdf",   (  0, 100, 104, 46), "f"),
    ("analysis_4_fda_antiaging_drugs.pdf",     (107, 100,  76, 46), None),
    ("geneage_waffle_nature.pdf",              (  0, 149, 183, 28), "g"),
    ("benchmark_schematic.pdf",                (  0, 180,  78, 54), "h"),
    ("benchmark_identification_bars.pdf",      ( 81, 180, 69.1, 54), "i"),
    ("benchmark_accuracy_bars.pdf",            (151.6, 180, 31.4, 54), "j"),
]


def build(out_pdf):
    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=PAGE_H)
    for fname, (x, y, w, h), letter in PANELS:
        s = fitz.open(src(fname))
        r = s[0].rect
        # panels are drawn at their print size: place at 100%, never rescale
        assert abs(r.width - w * MM) < 1.5 and abs(r.height - h * MM) < 1.5, (fname, r.width / MM, r.height / MM, w, h)
        page.show_pdf_page(fitz.Rect(x * MM, y * MM, x * MM + r.width, y * MM + r.height), s, 0)
        s.close()
    for fname, (x, y, w, h), letter in PANELS:
        if letter:
            page.insert_text((x * MM + 0.4, y * MM + LETTER_PT * 0.78), letter, fontsize=LETTER_PT,
                             fontname="hebo", color=(0, 0, 0))
    doc.save(out_pdf, garbage=4, deflate=True)
    doc.close()

def render_png(pdf_path, png_path, dpi=600):
    doc = fitz.open(pdf_path)
    pix = doc[0].get_pixmap(dpi=dpi, alpha=False)
    pix.save(png_path)
    doc.close()

# Canonical destinations to overwrite (PDF + matching 600-dpi PNG where applicable)
DESTS = [
    (os.path.join(ROOT, "figure", "Figure4.pdf"),
     os.path.join(ROOT, "figure", "Figure4.png")),
    (os.path.join(ROOT, "figure_rev1", "Figure4.pdf"),
     os.path.join(ROOT, "figure_rev1", "Figure4.png")),
    (os.path.join(ROOT, "submission_revision1_2026-05-23", "Figures", "Figure4.pdf"),
     None),
]
BACKUP_DIR = os.path.join(ROOT, "figure", "outdated")

def deploy():
    # 1. build the vector master once
    master = "/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base/clockbase_revision_local_2026-11-03/clockbase-agent/paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact/assembly_build_4/Figure4_vector.pdf"
    os.makedirs("/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base/clockbase_revision_local_2026-11-03/clockbase-agent/paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact/assembly_build_4", exist_ok=True)
    build(master)
    # 2. back up existing raster originals, then overwrite each destination
    os.makedirs(BACKUP_DIR, exist_ok=True)
    for pdf_dst, png_dst in DESTS:
        for f in (pdf_dst, png_dst):
            if f and os.path.exists(f):
                bak = os.path.join(BACKUP_DIR,
                                   os.path.basename(os.path.dirname(f)) + "_" +
                                   os.path.basename(f).replace(".", "_raster_2026-06-04."))
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
        render_png(out, out.replace(".pdf", ".png"), dpi=200)
        print("Built (preview):", out)
    else:
        deploy()
