#!/usr/bin/env python3
"""Assemble Extended Data Figure 22 from its existing vector panels without the sensitivity panel.

The panels are the PDFs written by scripts/figures_2026_11_06_july_exact/build_ed21_22_july_2026_11_06.py
(evidence/figures-rebuild-2026-11-06/ed22/panel_*.pdf). No panel is redrawn or scaled. The sensitivity grid
(former panel b) is left out; every specification it showed is a column of Supplementary Table S19.

Old letter -> new letter: a -> a, d -> b, c -> c, e -> d, f -> e, g -> f, h -> g.
Run: .venv-pipeline/bin/python scripts/figures_2026_11_20/assemble_ed22_main_2026_11_20.py
"""
from pathlib import Path
import hashlib
import json

import pymupdf as fitz

REPO = Path(__file__).resolve().parents[2]
BASE = REPO / 'paper/revision_2026-09-07'
SRC = BASE / 'evidence/figures-rebuild-2026-11-06/ed22'
OUT = BASE / 'evidence/ed22-main-2026-11-20'
PAGE_W = 595.68
HALF = PAGE_W / 2
# new letter: (old panel file letter, rect); heights are the native heights of the panel PDFs
LAYOUT = {'a': ('a', (0, 0, HALF, 145)), 'b': ('d', (HALF, 0, PAGE_W, 140)),
          'c': ('c', (0, 145, HALF, 370)), 'd': ('e', (HALF, 140, PAGE_W, 315)),
          'e': ('f', (HALF, 315, PAGE_W, 505)), 'f': ('g', (0, 370, HALF, 745)),
          'g': ('h', (HALF, 505, PAGE_W, 765))}
HEIGHT = 765


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    doc = fitz.open()
    page = doc.new_page(width=PAGE_W, height=HEIGHT)
    sources = {}
    for letter, (old, rect) in LAYOUT.items():
        path = SRC / f'panel_{old}.pdf'
        with fitz.open(path) as src:
            native = src[0].rect
            width, height = rect[2] - rect[0], rect[3] - rect[1]
            if abs(native.width - width) > 0.6 or abs(native.height - height) > 0.6:
                raise SystemExit(f'panel {old}: native {native.width:.1f} x {native.height:.1f}, slot {width:.1f} x {height:.1f}')
            page.show_pdf_page(fitz.Rect(rect), src, 0, keep_proportion=True)
        sources[letter] = {'panel': path.name, 'sha256': sha(path), 'rect_pt': rect}
    for letter, (_, rect) in LAYOUT.items():
        page.insert_text((rect[0] + 5, rect[1] + 16), letter, fontsize=13, fontname='hebo', color=(0, 0, 0))
    out = OUT / 'ed22.pdf'
    doc.save(out, garbage=4, deflate=True)
    doc.close()
    with fitz.open(out) as check:
        pg = check[0]
        pg.get_pixmap(dpi=200, alpha=False).save(OUT / 'ed22_review.png')
        spans = [s for b in pg.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans']
                 if s['text'].strip()]
        facts = {'width_pt': pg.rect.width, 'height_pt': pg.rect.height, 'raster_images': len(pg.get_images()),
                 'minimum_text_pt': min(s['size'] for s in spans), 'panels': sources,
                 'left_out': {'panel': 'panel_b.pdf (sensitivity grid)', 'reported_in': 'Supplementary Table S19'}}
        if facts['raster_images']:
            raise SystemExit('raster image in the assembled page')
        text = pg.get_text()
        if 'H–K' in text or 'Sensitivity' in text:
            raise SystemExit('the sensitivity panel is still on the page')
    (OUT / 'assembly.json').write_text(json.dumps(facts, indent=1))
    print('wrote', out, f"{facts['width_pt']:.1f} x {facts['height_pt']:.1f} pt; minimum text {facts['minimum_text_pt']:.1f} pt")


if __name__ == '__main__':
    main()
