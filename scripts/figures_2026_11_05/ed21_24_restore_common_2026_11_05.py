"""Shared mechanics for the ED21–24 restoration; July source defines style.

Style record: July Figure 3 / create_combined_modality_top_hits_v2.py supplies
direction colours, modality shapes, outlined bars and dotted grids. July ED7 /
src/visualize_validation_results.py supplies category colours and value labels.
July ED19 supplies bold panel titles, substantial marks and compact row spacing.
No retro_style, layout template, provenance footer or in-image figure title.
"""
import os
import hashlib
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.text import Text
from PIL import Image, ImageDraw
import pymupdf as fitz

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'scripts/figures_2026_10_03'))
from audit_system import paths
from nature_figure_utils import strict_figure_qc
from ed21_24_july_top_hits_copy_2026_11_05 import COLOR_PRO, COLOR_ANTI, MODALITY_MARKERS

FINAL = Path(os.environ.get('CLOCKBASE_FINAL') or paths.VERIFIER_ROOT / 'final-2026-10-10')
META = Path(os.environ.get('CLOCKBASE_META') or paths.VERIFIER_ROOT / 'intervention-meta-2026-10-17')
BASE = REPO / 'paper/revision_2026-09-07'
OUTBASE = BASE / 'evidence/figures-rebuild-2026-11-05'
JULY = BASE / 'baseline/03_Figures_2026-07-06.pdf'
WIDTH = 183 / 25.4
HIGHER, LOWER = COLOR_PRO, COLOR_ANTI
INK = '#222222'
GRAY = '#999999'
# Palette read from July ED7's generator, rather than the October soft palette.
CLASSIC = ['#4C72B0', '#DD8452', '#55A868', '#8172B3', '#CCB974', '#64B5CD']
MODS = {'mm_rna': ('Mouse RNA-seq', 'o', '#4C72B0'),
        'hs_rna': ('Human RNA-seq', 's', '#C44E52'),
        'dnam': ('Human DNAm', '^', '#55A868')}


def style():
    plt.rcdefaults()
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8,
                         'axes.labelsize': 8, 'axes.titlesize': 8,
                         'axes.titleweight': 'bold', 'xtick.labelsize': 7,
                         'ytick.labelsize': 7, 'legend.fontsize': 7,
                         'axes.linewidth': .6, 'xtick.major.width': .6,
                         'ytick.major.width': .6, 'xtick.major.size': 2.5,
                         'ytick.major.size': 2.5, 'pdf.fonttype': 42,
                         'ps.fonttype': 42, 'savefig.facecolor': 'white',
                         'figure.facecolor': 'white', 'axes.facecolor': 'white',
                         'text.color': INK, 'axes.labelcolor': INK,
                         'xtick.color': INK, 'ytick.color': INK})


def read_json(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}')
    try:
        return json.loads(path.read_text())
    except Exception as exc:
        raise RuntimeError(f'Cannot read {path}: {exc}') from exc


def axes(fig, box, title=None, letter=None, label_x=None):
    ax = fig.add_axes(box)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_axisbelow(True)
    if title:
        ax.set_title(title, loc='left', pad=8)
    if letter:
        fig.text(box[0] if label_x is None else label_x, box[1]+box[3]+.017,
                 letter, fontsize=11, weight='bold', color='black', va='bottom')
    return ax


def grid(ax, axis='x'):
    ax.grid(axis=axis, color='#bbbbbb', lw=.4, alpha=.5, linestyle=':')


def sci(value):
    if value == 0:
        return '0'
    if value < .01:
        exponent = int(np.floor(np.log10(value)))
        return f'{value / 10. ** exponent:.1f}×10$^{{{exponent}}}$'
    return f'{value:.3f}'


def bar_values(ax, ys, values, *, offset=0, fontsize=7, fmt=None):
    for yy, value in zip(ys, values):
        ax.text(value+offset, yy, fmt(value) if fmt else f'{value:,}',
                va='center', ha='left', fontsize=fontsize)


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def source_records(files):
    records = []
    for file in dict.fromkeys(Path(p).resolve() for p in files):
        if not file.is_file():
            raise FileNotFoundError(f'Cannot proceed - required data file not found: {file}')
        records.append({'path': str(file), 'sha256': sha256(file),
                        'mtime_ns': file.stat().st_mtime_ns, 'bytes': file.stat().st_size})
    return records


def save(fig, name, axd, july_page, *, sources, facts, script):
    """Check defects without imposing the skill's layout or palette rules."""
    out = OUTBASE / name
    out.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    report = strict_figure_qc(fig, axd, check_text_figure_edges=False,
                             forbid_uniform_grid=False, check_text_spine=False,
                             check_data_occupancy=False, check_max_fontsize=False,
                             check_matched_limits=False,
                             max_tall_panel_aspect=100,
                             max_sparse_panel_aspect=100)
    # Only the named defects control the restoration.
    defect_codes = {'font-size', 'overlap', 'text-clipped', 'text-clipping',
                    'missing-colorbar', 'legend-overflow'}
    defects = [i for i in report['issues'] if i['code'] in defect_codes]
    # Ink-mask verification renders at a different dpi. Redraw at the figure
    # dpi before measuring the fixed-page bounds of labels and legends.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    clipped = []
    fonts = []
    inactive_tick_text = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks()+axis.get_minor_ticks():
                if not lo <= tick.get_loc() <= hi:
                    inactive_tick_text.update([id(tick.label1),id(tick.label2)])
    for txt in fig.findobj(Text):
        if id(txt) in inactive_tick_text or not txt.get_visible() or not txt.get_text():
            continue
        bb = txt.get_window_extent(renderer)
        if bb.width == 0 or bb.height == 0:
            continue
        fonts.append(float(txt.get_fontsize()))
        if bb.x0 < -1 or bb.y0 < -1 or bb.x1 > fig.bbox.x1+1 or bb.y1 > fig.bbox.y1+1:
            clipped.append(txt.get_text())
    qc = {'defects': defects, 'canvas_clipping': clipped,
          'minimum_text_size_pt': min(fonts), 'size_mm': [183, fig.get_figheight()*25.4],
          'skill_scope': 'Legibility, overlap, clipping and colour scales only; July appearance retained.'}
    (out / 'checks.json').write_text(json.dumps(qc, indent=2, default=str)+'\n')
    pdf = out / f'{name}.pdf'
    fig.savefig(pdf)
    fig.savefig(out / f'{name}_review.png', dpi=200)
    doc = fitz.open(pdf)
    page = doc[0]
    qc['pdf_images'] = len(page.get_images())
    qc['pdf_vector_paths'] = len(page.get_drawings())
    qc['pdf_size_pt'] = list(page.rect)[2:]
    qc['embedded_fonts'] = [list(f) for f in page.get_fonts()]
    qc['passed'] = not defects and not clipped and qc['pdf_images'] == 0
    (out / 'checks.json').write_text(json.dumps(qc, indent=2, default=str)+'\n')
    (out / 'numbers.json').write_text(json.dumps(facts, indent=2, default=str)+'\n')
    (out / 'provenance.json').write_text(json.dumps({
        'inputs': source_records(sources), 'builder': source_records([script, __file__]),
        'output_sha256': sha256(pdf), 'july_style_reference': str(JULY),
        'july_style_page': july_page}, indent=2)+'\n')
    # Full July page at its original page size beside the new 183-mm-wide figure.
    # Both are rasterised at 150 dpi, so one inch has the same size on either side.
    ref = fitz.open(JULY)[july_page-1].get_pixmap(dpi=150)
    old = Image.frombytes('RGB', [ref.width, ref.height], ref.samples)
    pix = page.get_pixmap(dpi=150)
    new = Image.frombytes('RGB', [pix.width, pix.height], pix.samples)
    canvas = Image.new('RGB', (old.width+new.width+36, max(old.height,new.height)+48), 'white')
    draw = ImageDraw.Draw(canvas)
    draw.text((10,10), f'July page {july_page} at 150 dpi', fill='black')
    draw.text((old.width+28,10), f'{name.upper()} at 150 dpi', fill='black')
    canvas.paste(old,(0,48)); canvas.paste(new,(old.width+36,48))
    canvas.save(out/f'{name}_vs_july.png')
    print(name, 'QC', qc['passed'], 'defects', len(defects), 'clipped', clipped)
    for issue in defects[:20]:
        print(issue['panel'], issue['code'], issue['message'])
    plt.close(fig)
    return out


def table_md(frame):
    return frame.to_markdown(index=False, floatfmt='.6g')


def readme(out, content):
    try:
        (out / 'README.md').write_text(content.strip()+'\n')
    except OSError as exc:
        raise RuntimeError(f'Cannot write restoration report to {out}: {exc}') from exc
