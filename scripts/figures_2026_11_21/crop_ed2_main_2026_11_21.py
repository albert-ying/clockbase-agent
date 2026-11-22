#!/usr/bin/env python3
"""Extended Data Figure 2 with panels a and b.

The page of evidence/figure-fixes-2026-11-17/ed2/ed2.pdf is shown clipped above panels c and d (the comparison of
GrimAge age and sex inputs, reported in Supplementary Note 5). Nothing is redrawn or scaled.
Run: .venv-pipeline/bin/python scripts/figures_2026_11_21/crop_ed2_main_2026_11_21.py
"""
from pathlib import Path
import hashlib
import json

import pymupdf as fitz

REPO = Path(__file__).resolve().parents[2]
BASE = REPO / 'paper/revision_2026-09-07'
SRC = BASE / 'evidence/figure-fixes-2026-11-17/ed2/ed2.pdf'
OUT = BASE / 'evidence/ed2-main-2026-11-21'


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src = fitz.open(SRC)
    page = src[0]
    spans = [s for b in page.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans']
             if s['text'].strip()]
    letters = {s['text'].strip(): s['bbox'] for s in spans if s['text'].strip() in 'abcd' and s['size'] > 8}
    if sorted(letters) != ['a', 'b', 'c', 'd']:
        raise SystemExit(f'panel letters found: {sorted(letters)}')
    top_c = min(letters['c'][1], letters['d'][1])
    above = [s for s in spans if s['bbox'][3] < top_c - 1]
    below = [s for s in spans if s['bbox'][1] >= top_c - 1]
    bottom_b = max(s['bbox'][3] for s in above)
    drawn = [d['rect'] for d in page.get_drawings()]
    crossing = [r for r in drawn if r.y0 < bottom_b + 2 and r.y1 > top_c - 2 and r.height < page.rect.height * 0.9]
    if crossing:
        raise SystemExit(f'{len(crossing)} drawings cross the cut')
    cut = round((bottom_b + top_c) / 2 + 2, 1)
    out = fitz.open()
    new = out.new_page(width=page.rect.width, height=cut)
    clip = fitz.Rect(0, 0, page.rect.width, cut)
    new.show_pdf_page(clip, src, 0, clip=clip)
    target = OUT / 'ed2.pdf'
    out.save(target, garbage=4, deflate=True)
    out.close()
    check = fitz.open(target)[0]
    text = check.get_text()
    kept = [s['text'].strip() for b in check.get_text('dict')['blocks'] if 'lines' in b for line in b['lines']
            for s in line['spans'] if s['text'].strip()]
    for word in ('Horvath age', 'All male', 'Recorded age', 'Coded sex', 'Residual ρ with Horvathv1'):
        if word in text:
            raise SystemExit(f'{word!r} is still on the page')
    if len(kept) != len(above):
        raise SystemExit(f'text spans kept {len(kept)}, expected {len(above)}')
    check.get_pixmap(dpi=200, alpha=False).save(OUT / 'ed2_review.png')
    facts = {'source': str(SRC.relative_to(BASE)), 'source_sha256': hashlib.sha256(SRC.read_bytes()).hexdigest(),
             'source_height_pt': page.rect.height, 'cut_pt': cut, 'bottom_of_panel_b_text_pt': bottom_b,
             'top_of_panels_c_d_pt': top_c, 'text_spans_kept': len(kept), 'text_spans_left_out': len(below),
             'raster_images': len(check.get_images()),
             'output_sha256': hashlib.sha256(target.read_bytes()).hexdigest()}
    (OUT / 'crop.json').write_text(json.dumps(facts, indent=1))
    print('wrote', target, f"{check.rect.width:.1f} x {check.rect.height:.1f} pt; spans kept {len(kept)}, left out {len(below)}")


if __name__ == '__main__':
    main()
