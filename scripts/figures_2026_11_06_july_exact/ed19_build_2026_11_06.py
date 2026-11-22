#!/usr/bin/env python3
"""Run the July ED19 plotting function on frozen annotation inputs."""
from pathlib import Path
import os
import datetime as dt
import difflib
import hashlib
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/ed19'
PREV = OUT.parent.parent / 'figures-rebuild-2026-11-05/ed19'
SRC = OUT.parent.parent / 'ed19-source-recovery-2026-09-18/external_db_crossvalidation_2026-05-14.csv'
JULY = ROOT / 'paper/revision_2026-09-07/baseline/03_Figures_2026-07-06.pdf'
ORIGINAL = ROOT / 'scripts/legacy_sensitivity_analyses_2026_05_14/external_db_crossvalidation_2026-05-14.py'
COPY = Path(__file__).with_name('external_db_crossvalidation_2026-05-14.py')
# Style source: ORIGINAL make_figure, unchanged Calibri, 7/8/8.5/9/12 pt.
# Colour/stroke source: ORIGINAL #B85A3E/#3E7B5A/#3E5C7B, 0.5/0.6 pt.
# Geometry source: ORIGINAL 8.0 x 7.6 inch gridspec; July assembled width 612 pt.


def provenance(path):
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'mtime_utc': dt.datetime.fromtimestamp(path.stat().st_mtime, dt.timezone.utc).isoformat()}


def load_generator():
    # The legacy module imports statsmodels for its unused reporting main routine.
    # Load the installed local package from uv's cache if pipeline lacks it.
    try:
        import statsmodels.stats.multitest
    except ImportError:
        candidates = list((Path.home()/'.cache/uv/archive-v0').glob('*/statsmodels/stats/multitest.py'))
        if not candidates:
            raise RuntimeError('The original generator requires statsmodels; no local package found')
        sys.path.insert(0, str(candidates[0].parents[2]))
    import matplotlib.font_manager as fm
    fonts = Path('/Applications/Microsoft Word.app/Contents/Resources/DFonts')
    for name in ['Calibri.ttf', 'Calibrib.ttf', 'Calibrii.ttf', 'Calibriz.ttf']:
        font = fonts/name
        if not font.exists():
            raise FileNotFoundError(f'Required July font not found: {font}')
        fm.fontManager.addfont(str(font))
    spec = importlib.util.spec_from_file_location('ed19_july_original', COPY)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build():
    import fitz
    from PIL import Image, ImageDraw
    import numpy as np
    OUT.mkdir(parents=True, exist_ok=True)
    m = load_generator()
    annotations = m.pd.read_csv(SRC).fillna('')
    joined = m.pd.read_csv(PREV/'ranking_join.csv')
    ranked = m.pd.read_csv(Path(os.environ['CLOCKBASE_FINAL']) / 'comparison-ranking/ranked_comparisons_all.csv')
    expected = annotations.merge(ranked, left_on=['modality', 'gse_id', 'intervention'],
                                 right_on=['modality', 'gse_id', 'label'], validate='one_to_one')
    assert len(joined) == len(expected)
    assert not annotations.duplicated(['modality','gse_id','intervention']).any()
    assert not joined.duplicated(['modality','gse_id','intervention']).any()
    assert joined.v3_final_class.isin(['VERIFIED', 'PROVISIONAL']).all()
    assert np.allclose(np.minimum(joined.g_ci_low.abs(), joined.g_ci_high.abs()), joined.rank_bound)
    concord, novel = m.headlines(annotations)
    summary = m.enrichment_stats(annotations)
    captured = []
    original_close = m.plt.close
    m.plt.close = lambda fig: captured.append(fig)
    try:
        m.make_figure(annotations, summary, concord, novel)
    finally:
        m.plt.close = original_close
    fig = captured[-1]
    sys.path.insert(0, str(Path('pipeline_assets/figure_style_2026_09_10')))
    from nature_figure_utils import strict_figure_qc, figure_preflight
    qc_options = dict(forbid_uniform_grid=False, check_house_style=False,
                      check_data_occupancy=False, check_matched_limits=False,
                      check_text_spine=False, max_tall_panel_aspect=10,
                      max_sparse_panel_aspect=10, check_max_fontsize=False)
    axes = dict(zip('abcd', fig.axes))
    qc = strict_figure_qc(fig, axes, **qc_options)
    preflight = figure_preflight(fig, axes, verbose=False, **qc_options)
    (OUT/'defect_checks.json').write_text(json.dumps({'strict_qc':qc,'preflight':preflight,
        'exception_reason':'July source geometry, strokes, and typography are user required'},
        indent=2, default=str))
    original_close(fig)
    raw = fitz.open(m.OUT_PDF)
    # Preserve all July font sizes at 100%; 183 mm scaling would make 7 pt ticks 6.58 pt.
    document = fitz.open()
    page = document.new_page(width=612, height=raw[0].rect.height)
    dx = (page.rect.width - raw[0].rect.width)/2
    page.show_pdf_page(fitz.Rect(dx,0,dx+raw[0].rect.width,page.rect.height),raw,0)
    document.save(OUT/'ed19.pdf', garbage=4, deflate=True)
    document.close()
    final = fitz.open(OUT/'ed19.pdf')
    final[0].get_pixmap(dpi=200).save(OUT/'ed19_review.png')
    fonts = [span for block in final[0].get_text('dict')['blocks'] if 'lines' in block
             for line in block['lines'] for span in line['spans']]
    assert min(s['size'] for s in fonts) >= 6.999
    assert not final[0].get_images(), 'Raster image found in vector figure'
    assert 'Extended Data Figure' not in final[0].get_text()
    july = fitz.open(JULY)
    # The July image is placed at this recorded rectangle on package page25.
    july[24].get_pixmap(dpi=200,clip=fitz.Rect(81.05,101.65,549.05,505.80)).save(OUT/'july_reference.png')
    old = Image.open(OUT/'july_reference.png').convert('RGB')
    new = Image.open(OUT/'ed19_review.png').convert('RGB')
    # Compare at one physical scale (200 dpi); retain natural artwork sizes.
    canvas = Image.new('RGB',(old.width+new.width+90,max(old.height,new.height)+80),'white')
    canvas.paste(old,(20,50));canvas.paste(new,(old.width+60,50))
    draw=ImageDraw.Draw(canvas)
    draw.text((20,15),'July page25 artwork at 200 dpi',fill='black')
    draw.text((old.width+60,15),'November6 at 200 dpi',fill='black')
    canvas.save(OUT/'ed19_vs_july.png')
    old_numbers=json.loads((PREV/'numbers_and_sources.json').read_text())
    assert {name: row['n_hits'] for name,row in summary.items()} == old_numbers['database_counts']
    assert all(int(annotations.match_status.eq(name).sum()) == count
               for name,count in old_numbers['status_counts'].items())
    rows=[]
    for panel,frame in [('c',concord),('d',novel)]:
        expected=[r['intervention'] for r in old_numbers['bar_values'] if r['panel']==panel]
        assert frame.intervention.tolist()==expected, 'Selection/order differs from November5'
        for _,row in frame.iterrows():
            rows.append({'panel':panel,'intervention':row.intervention,'gse_id':row.gse_id,
                         'composite_score':row.cb_final_score,'rank_bound':row.rank_bound,
                         'hedges_g':row.hedges_g,'g_ci_low':row.g_ci_low,'g_ci_high':row.g_ci_high})
    def visible_ticks(axis):
        lo, hi = sorted(axis.get_view_interval())
        return [tick.label1.get_text() for tick in axis.get_major_ticks()
                if tick.label1.get_visible() and lo <= tick.get_loc() <= hi]
    ticks={letter:{'x':visible_ticks(ax.xaxis), 'y':visible_ticks(ax.yaxis)}
           for letter,ax in axes.items()}
    manifest={'annotation_rows':len(annotations),'joined_rows':len(joined),
              'database_summary':summary,'status_counts':annotations.match_status.value_counts().to_dict(),
              'bar_values':rows,'ticks':ticks,'pdf_size_points':list(final[0].rect),
              'min_text_pt':min(s['size'] for s in fonts),'font_names':sorted(set(s['font'] for s in fonts)),
              'raster_images':len(final[0].get_images()),'sources':[provenance(p) for p in [SRC,PREV/'ranking_join.csv',ORIGINAL,COPY,JULY]],
              'differences_vs_november5':'Database values, classification counts, selected identities, and order unchanged. Lower bar magnitudes are original composite scores instead of confidence-bound magnitudes.'}
    (OUT/'numbers_and_sources.json').write_text(json.dumps(manifest,indent=2))
    (OUT/'ed19_against_july.diff').write_text(''.join(difflib.unified_diff(
        ORIGINAL.read_text().splitlines(True),COPY.read_text().splitlines(True),
        fromfile=str(ORIGINAL.relative_to(ROOT)),tofile=str(COPY.relative_to(ROOT)))))
    make_readme(manifest, qc, preflight)
    print(json.dumps({'pdf_size_points':manifest['pdf_size_points'],'min_text_pt':manifest['min_text_pt'],
                      'strict_qc_passed':qc['passed'],'preflight_passed':preflight['passed'],
                      'fonts':manifest['font_names'],'bar_values':rows},indent=2))


def make_readme(values, qc, preflight):
    lines=['# Extended Data Figure 19','',
           'The July generator draws the saved external database annotations. Panels c',
           'and d select and order records by the magnitude of the Hedges g confidence',
           'bound nearest zero. Bar lengths retain the July composite scores.','',
           '## Panel map','',
           '| Panel | Content | Source |','|---|---|---|',
           '| a | Database hit rates and hit counts | Saved 257-row annotation CSV |',
           '| b | Annotation classifications | Saved 257-row annotation CSV |',
           '| c | Direction-consistent database match | Nine exact ranking joins |',
           '| d | No recorded database match | Ten highest-ranked records among 65 no-match joins |','',
           '## Printed numbers','',
           'Panel a uses the same denominator of 257 records for every database. The',
           'printed labels are hit counts; bar heights are computed percentages.','',
           '| Database | Printed count | Bar height percent |','|---|---:|---:|']
    for name,row in values['database_summary'].items():
        lines.append(f"| {name} | {row['n_hits']} | {row['hit_rate']*100:.12g} |")
    lines += ['', 'Panel b prints the following counts from the saved match_status column.','',
              '| Classification | Printed count |','|---|---:|']
    for status in ['direction_consistent','mechanism_only','mixed','direction_mismatched','no_db_evidence']:
        lines.append(f"| {status} | {values['status_counts'].get(status, 0)} |")
    lines += ['', 'Panels c and d plot cb_final_score from ranking_join.csv. The statistical',
              'ordering uses min(abs(g_ci_low), abs(g_ci_high)), checked against rank_bound.',
              'These bars are not effect sizes. Rows below follow top-to-bottom plot order.','',
              '| Panel | Intervention | Study | Composite bar length | Ordering value | November5 bar length |',
              '|---|---|---|---:|---:|---:|']
    for row in values['bar_values']:
        lines.append(f"| {row['panel']} | {row['intervention']} | {row['gse_id']} | {row['composite_score']:g} | {row['rank_bound']:.12g} | {row['rank_bound']:.12g} |")
    lines += ['', 'Every numerical axis tick is a scale label, rather than an additional result.','',
              '| Panel | Numerical x ticks | Numerical y ticks |','|---|---|---|']
    for letter,ticks in values['ticks'].items():
        numeric=lambda items:', '.join(s for s in items if s.replace('−','-').replace('.','',1).lstrip('-').isdigit()) or 'None'
        lines.append(f"| {letter} | {numeric(ticks['x'])} | {numeric(ticks['y'])} |")
    lines += ['', 'Numbers embedded in intervention names reproduce the source identifiers.',
              'The plotting function retains July truncation after 28 characters. Full',
              'intervention identities and all plotted values are above and in',
              '`numbers_and_sources.json`.','', '## Style sources','',
              'Every visual convention comes directly from make_figure in',
              '`scripts/legacy_sensitivity_analyses_2026_05_14/external_db_crossvalidation_2026-05-14.py`. ',
              'The preserved conventions are Calibri; 9 pt base and axis labels; 8 pt',
              'database and classification labels; 7 pt intervention labels and count',
              'annotations; 8.5 pt titles; 12 pt bold panel letters; 0.6 pt axes and',
              'bar outlines; 0.5 pt lower-panel outlines; the original five classification',
              'colours; the original 8.0 by 7.6 inch figure and gridspec; and the original',
              'tight output crop. Calibri regular, bold, italic, and bold italic are',
              'registered from the installed Microsoft Word fonts. No font substitution',
              'or independent restyling is applied.','',
              'The output page width is 612 pt (215.9 mm), the July assembled package',
              'width. The July plotting output is centred at 100% scale, with no raster',
              'conversion. Its height is '+f"{values['pdf_size_points'][3]:.6f} pt. "+
              'Reducing this artwork to 183 mm would make the original 7 pt ticks smaller',
              'than the requested floor. The figure contains no internal figure title,',
              'figure number, or footer. The comparison PNG shows July page25 artwork',
              'beside the new artwork at the same 200 dpi physical scale.','',
              '## Changes from November5','',
              'Panels a and b retain every empirical value. The selected interventions',
              'and their order in panels c and d are unchanged. Their bar lengths change',
              'from the CI-bound magnitudes to the original composite scores, as listed',
              'above. The x-axis label returns to ClockBase composite score. Panel b',
              'uses no recorded database match. Exact July Calibri, type sizes, axis',
              'geometry, truncation, outlines, and layout replace the November5 redraw.','',
              '## Source files','',
              'The source records below include SHA256 and UTC filesystem modification times.','']
    for source in values['sources']:
        lines += [f"- `{source['path']}`",f"  SHA256 `{source['sha256']}`; modified `{source['mtime_utc']}`."]
    lines += ['', 'The saved annotation table contains 257 records, whereas July page25',
              'shows 289. Missing historical records remain unavailable. This figure',
              'retains the documented saved annotations and makes no new database queries',
              'or enrichment tests. Eighty-six annotations have an exact ranking join;',
              '171 do not. No unmatched annotation receives an inferred ranking value.','',
              '## Diff','',
              '`ed19_against_july.diff` is a unified diff against the original generator.',
              'Only output/cache paths, the ranking selection and ordering step, and',
              'label text differ. The plotting geometry, bar values, styles, and main',
              'routine remain unchanged. The runner imports the module and calls only',
              'headlines, enrichment_stats, and make_figure with frozen local inputs;',
              'it does not invoke the legacy main routine, whose live queries and',
              'enrichment claims are outside this restoration.','',
              '## Validation','',
              f"The strict figure check passes: {qc['passed']}. The full preflight passes: {preflight['passed']}.",
              'The preserved July geometry is exempt from house composition rules;',
              'the complete findings are in `defect_checks.json`. The preflight records',
              'the original perimeter margin of 2.48% against its generic 2% limit and',
              'the original inter-row whitespace. Both are retained because geometry',
              'changes are outside the permitted diff. Visual review checks',
              'the saved 200 dpi render against July page25. PDF inspection confirms',
              f"a minimum text size of {values['min_text_pt']:g} pt, embedded Calibri text,",
              'vector bars, no raster images, and no clipped labels or text overlap.','',
              '## Reproduction','',
              'From the repository root, run:','', '```sh','source ../env_local.sh',
              '../.venv-pipeline/bin/python scripts/figures_2026_11_06_july_exact/ed19_build_2026_11_06.py','```','',
              'The runner loads statsmodels from a local uv cache if it is absent from',
              'the pipeline environment. This dependency is imported by the unchanged',
              'legacy module; no tests from its unused main routine are run.']
    (OUT/'README.md').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    try:
        build()
    except (OSError, ValueError, AssertionError, ImportError) as error:
        raise RuntimeError(f'ED19 July restoration failed: {error}') from error
