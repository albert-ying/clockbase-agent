#!/usr/bin/env python3
"""Redraw the unchanged November 5 ED23/24 data with July plotting conventions."""
from pathlib import Path
import ast
import difflib
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image

REPO = Path(__file__).resolve().parents[2]
EVIDENCE = REPO / 'paper/revision_2026-09-07/evidence'
OLD = EVIDENCE / 'figures-rebuild-2026-11-05'
NEW = EVIDENCE / 'figure-fixes-2026-11-17'
JULY = REPO / 'paper/revision_2026-09-07/baseline/03_Figures_2026-07-06.pdf'
LEGACY = REPO / 'paper/legacy_reference/figure_rev1'
TOP = REPO.parents[1] / 'Aging-AI-Agent/create_combined_modality_top_hits_v2.py'
# July ED19 supplies this rcParams block, font, ticks and spine geometry.
# July ED7 supplies black-edged category bars, alpha .8 and blue/orange palette.
# July Figure 3b-e supplies the directional key and outlined modality marks.
STYLE = {
    'ed7': TOP.parent / 'src/visualize_validation_results.py',
    'ed14': REPO / 'scripts/legacy_sensitivity_analyses_2026_05_14/extreme_tail_dedicated_analysis.py',
    'ed16': LEGACY / 'make_ed16_pooled_volcano.py',
    'ed17': LEGACY / 'sensitivity_analyses/effect_size_distribution_2026-05-14.py',
    'ed19': REPO / 'scripts/legacy_sensitivity_analyses_2026_05_14/external_db_crossvalidation_2026-05-14.py',
    'figure3': TOP,
    'assembly': LEGACY / 'assemble_figure4_vector.py',
    'sample_boxes': LEGACY / 'sample_level_examples_2026-06-13/make_example_figures.py',
}
W = 183 / 25.4 * 72
BLUE, ORANGE = '#0072B2', '#D55E00'
CLASSIC = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3', '#937860']
PREFLIGHT = {}
PRINTED = {}


def checked(path):
    if not path.is_file():
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}')
    return path


def sha(path):
    return hashlib.sha256(checked(path).read_bytes()).hexdigest()


def style():
    for p in Path('/Applications/Microsoft Word.app/Contents/Resources/DFonts').glob('Calibr*.ttf'):
        font_manager.fontManager.addfont(str(p))
    plt.rcParams.update({'font.family': 'Calibri', 'font.size': 9,
        'axes.linewidth': .6, 'axes.spines.top': False, 'axes.spines.right': False,
        'xtick.major.width': .6, 'ytick.major.width': .6,
        'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'axes.labelsize': 8.5,
        'axes.titlesize': 9, 'axes.titleweight': 'bold',
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white'})


def panel(name, letter, width, height, bounds):
    # Half panels retain the July Figure3b-e 7:5 figure proportion; the full
    # rows retain the ED7 14:10 page grammar, with height set by category count.
    fig = plt.figure(figsize=(width / 72, height / 72), dpi=200)
    ax = fig.add_axes(bounds)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, pad=2)
    return fig, ax


def grid(ax, axis='x'):
    ax.grid(axis=axis, color='#cccccc', linewidth=.5, alpha=.5, linestyle=':')
    ax.tick_params(axis='y' if axis == 'x' else 'x', length=0)


def sci(value):
    # Preserve the November 5 printed precision as well as the source value.
    if value == 0:
        return '0'
    if value < .01:
        exponent = int(np.floor(np.log10(value)))
        superscript = str(exponent).translate(str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹'))
        return f'{value / 10. ** exponent:.1f}×10{superscript}'
    return f'{value:.3f}'


def save_panel(fig, ax, name, letter):
    out = NEW / name / 'panels'
    out.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    # July references override the house's palette, grid and title defaults.
    # Keep the full report for checking only the structural findings.
    from nature_figure_utils import figure_preflight
    report = figure_preflight(fig, {letter: ax}, verbose=False)
    PREFLIGHT[f'{name}_{letter}'] = report['text']
    txt = []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            txt.append({'text': t.get_text(), 'font_pt': t.get_fontsize()})
    PRINTED[f'{name}_{letter}'] = txt
    fig.savefig(out / f'panel_{letter}.pdf')
    plt.close(fig)


def assemble(name, rects, height, july_page, letters=None):
    """rects maps a drawn panel file to its place; letters maps it to the printed letter (default: the same)."""
    letters = letters or {}
    out = NEW / name
    doc = fitz.open()
    page = doc.new_page(width=W, height=height)
    for letter, rect in rects.items():
        with fitz.open(out / 'panels' / f'panel_{letter}.pdf') as src:
            page.show_pdf_page(fitz.Rect(rect), src, 0, keep_proportion=True)
        # Original July Figure4 vector assembly: fixed rects, bold Helvetica.
        page.insert_text((rect[0] + 2, rect[1] + 14), letters.get(letter, letter),
                         fontsize=12, fontname='hebo', color=(0, 0, 0))
    pdf = out / f'{name}.pdf'
    doc.save(pdf, garbage=4, deflate=True)
    page.get_pixmap(dpi=200, alpha=False).save(out / f'{name}_review.png')
    # Compare the figure area in July at precisely the same physical scale.
    with fitz.open(JULY) as j:
        jp = j[july_page - 1]
        # Page text heading is outside the plotting region; preserve the page
        # figure at its printed scale and omit the reference heading.
        clip = fitz.Rect(64, 95, 553, 445 if july_page == 13 else 415)
        pix = jp.get_pixmap(dpi=200, alpha=False, clip=clip)
        ref = Image.frombytes('RGB', [pix.width, pix.height], pix.samples)
    now = Image.open(out / f'{name}_review.png').convert('RGB')
    gap = 48
    comparison = Image.new('RGB', (now.width + gap + ref.width, max(now.height, ref.height)), 'white')
    comparison.paste(now, (0, 0)); comparison.paste(ref, (now.width + gap, 0))
    comparison.save(out / f'{name}_vs_july.png', dpi=(200, 200))
    spans = [s for b in page.get_text('dict')['blocks'] if b['type'] == 0
             for line in b['lines'] for s in line['spans']]
    checks = {'page_width_pt': W, 'page_height_pt': height,
              'embedded_images': len(page.get_images()), 'vector_paths': len(page.get_drawings()),
              'minimum_printed_font_pt': min(s['size'] for s in spans),
              'all_spans_inside_page': all(page.rect.contains(fitz.Rect(s['bbox'])) for s in spans),
              'source_csvs_byte_identical': all(sha(p) == sha(out / p.name) for p in (OLD / name).glob('*.csv')),
              'numbers_json_byte_identical': sha(OLD / name / 'numbers.json') == sha(out / 'numbers.json'),
              'review_dpi': 200, 'numeric_changes': []}
    assert checks['embedded_images'] == 0 and checks['minimum_printed_font_pt'] >= 6
    assert checks['all_spans_inside_page'] and checks['source_csvs_byte_identical']
    (out / 'checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    doc.close()


def prepare(name):
    out = NEW / name
    out.mkdir(parents=True, exist_ok=True)
    for src in (OLD / name).glob('*.csv'):
        shutil.copy2(src, out / src.name)
    for fn in ['numbers.json', 'provenance.json']:
        shutil.copy2(checked(OLD / name / fn), out / fn)
    return out


def ed23():
    name = 'ed23'; out = prepare(name)
    numbers = json.loads((out / 'numbers.json').read_text())
    half = W / 2; h = half * 5 / 7
    # Panels c and d (ClockBase against language models without the data) are drawn here and placed in
    # Figure 4h,i (assemble_figure4_vector.py); this page prints the drawn e and f as c and d.
    rects = {'a': (0, 0, half, h), 'b': (half, 0, W, h),
             'e': (0, h + 9, W, h + 213),
             'f': (0, h + 222, W, h + 426)}
    printed_letters = {'e': 'c', 'f': 'd'}
    a = pd.read_csv(out / 'panel_a_funnel.csv')
    fig, ax = panel(name, 'a', half, h, [.44, .23, .51, .62])
    y = np.arange(len(a))[::-1]
    ax.barh(y, a.n, color=[CLASSIC[2],CLASSIC[0],CLASSIC[4],CLASSIC[1],CLASSIC[3],CLASSIC[5],'#999999','#aaaaaa','#cccccc'],
            edgecolor='black', linewidth=.5, alpha=.8, height=.68)
    for yi, n in zip(y, a.n):
        ax.text(n+3, yi, str(n), va='center', fontsize=6.5)
    ax.set_yticks(y, a.reason); ax.set_xlim(0, 200); ax.set_xticks([0,100,200])
    ax.set_xlabel('Published effects'); ax.set_title(f'Ground truth ({numbers["ground_truth_total"]} effects)', pad=7)
    grid(ax); save_panel(fig, ax, name, 'a')

    fig, ax = panel(name, 'b', half, h, [.26, .45, .69, .36])
    outcomes = [('significant, same direction','Sig., same',BLUE),
                ('not significant, same direction','NS, same','#9fc2df'),
                ('not significant, opposite direction','NS, opposite','#f2c49b'),
                ('significant, opposite direction','Sig., opposite',ORANGE),
                ('not testable','Not testable','#cccccc')]
    for yi, key in zip([1,0], ['S3','S3c']):
        left = 0
        for field, label, color in outcomes:
            n = numbers['b'][key].get(field,0)
            if n:
                ax.barh(yi,n,left=left,color=color,height=.55,edgecolor='black',linewidth=.5,alpha=.8)
                if n >= 3:
                    ax.text(left+n/2,yi,str(n),ha='center',va='center',fontsize=6.5)
                left += n
    ax.set_yticks([1,0],['Welch test','Verification test']); ax.set_xlim(0,55); ax.set_xticks([0,25,50])
    ax.set_xlabel('Published significant effects')
    ax.set_title('Clock outputs in published groups', pad=9)
    # Dedicated lower band, following the July ED7 small in-panel key.
    ax.legend(handles=[Patch(fc=c,ec='black',lw=.5,label=l) for _,l,c in outcomes],
              loc='upper left',bbox_to_anchor=(-.32,-.66),ncol=2,frameon=False,
              fontsize=7,handlelength=1,columnspacing=.8,labelspacing=.25)
    grid(ax); save_panel(fig, ax, name, 'b')

    cd = pd.read_csv(out / 'panel_c_d_models.csv')
    y = np.arange(len(cd))[::-1]
    fig, ax = panel(name, 'c', half, h, [.36,.24,.59,.57])
    for yi, row in zip(y, cd.itertuples()):
        lo, hi = row.lower_n/row.lower_total*100, row.higher_n/row.higher_total*100
        ax.plot([lo,hi],[yi,yi],color='#aaaaaa',lw=.8)
        ax.plot(lo,yi,'o',ms=4.5,mfc=BLUE,mec='#111',mew=.4,clip_on=False)
        ax.plot(hi,yi,'o',ms=4.5,mfc=ORANGE,mec='#111',mew=.4,clip_on=False)
    ax.set_yticks(y,cd.system);ax.set_xlim(-3,103);ax.set_xticks([0,50,100])
    ax.set_xlabel('Published effects identified (%)');ax.set_title('Published effect identification',pad=23)
    ax.plot([],[],'o',color=BLUE,label='Lower (16)',ms=4.5)
    ax.plot([],[],'o',color=ORANGE,label='Higher (36)',ms=4.5)
    ax.legend(loc='lower left',bbox_to_anchor=(-.04,1.0),ncol=2,frameon=False,fontsize=7,
              handletextpad=.25,columnspacing=.7)
    grid(ax);save_panel(fig,ax,name,'c')

    fig,ax=panel(name,'d',half,h,[.36,.24,.59,.57])
    for yi,row in zip(y,cd.itertuples()):
        lo,hi=ast.literal_eval(row.interval);v=row.balanced_accuracy
        color=BLUE if yi==len(cd)-1 else '#777777'
        ax.errorbar(v,yi,xerr=[[v-lo],[hi-v]],fmt='o',ms=4.5,color=color,
                    markeredgecolor='#111',markeredgewidth=.4,elinewidth=.8,capsize=0)
    ax.set_yticks(y,cd.system);ax.set_xlim(.3,.9);ax.set_xticks([.3,.6,.9])
    ax.axvline(1/3,color='gray',ls='--',lw=.6);ax.set_xlabel('Balanced accuracy')
    ax.set_title('Three-way call',pad=23);grid(ax);save_panel(fig,ax,name,'d')

    e=pd.read_csv(out/'panel_e_agreement.csv');y=np.arange(len(e))[::-1]
    fig,ax=panel(name,'e',W,204,[.32,.20,.61,.67])
    for yi,row in zip(y,e.itertuples()):
        lo,hi=ast.literal_eval(row.interval) if pd.notna(row.interval) else (row.agreement,row.agreement)
        color=CLASSIC[2] if yi==len(e)-1 else BLUE if yi==len(e)-2 else ORANGE if yi>=5 else '#777777'
        ax.errorbar(row.agreement*100,yi,xerr=[[(row.agreement-lo)*100],[(hi-row.agreement)*100]],
                    fmt='o',ms=4.5,color=color,markeredgecolor='#111',markeredgewidth=.4,elinewidth=.8,capsize=0)
        ax.text(1.025,yi,str(row.n),transform=ax.get_yaxis_transform(),va='center',fontsize=6.5)
    ax.text(1.025,len(e)-.1,'n',transform=ax.get_yaxis_transform(),fontsize=7,va='bottom')
    ax.set_yticks(y,e.system.replace({'Historical agent':'ClockBase agent','Released rows':'Verified rows'}));ax.set_xlim(0,100);ax.set_xticks([0,25,50,75,100]);ax.set_ylim(-.6,len(e)-.3)
    ax.set_xlabel('Agreement with published call (%)');ax.set_title(f'Benchmark series ({numbers["e_series"]})',pad=8)
    grid(ax);save_panel(fig,ax,name,'e')

    f=pd.read_csv(out/'panel_f_ranked_counts.csv');y=np.arange(len(f))[::-1]
    fig,ax=panel(name,'f',W,204,[.32,.26,.61,.62])
    for j,(field,label,color) in enumerate([
        ('higher_found','Higher found','#f2c49b'),('higher_recurs','Higher recurs',ORANGE),
        ('lower_found','Lower found','#9fc2df'),('lower_recurs','Lower recurs',BLUE)]):
        pos=y+.3-j*.2
        ax.barh(pos,f[field],height=.19,color=color,edgecolor='black',lw=.5,alpha=.8,label=label)
    # One paired count label per direction leaves the bar marks in the July
    # grouped style while avoiding four tightly stacked numerical labels.
    for yi,row in zip(y,f.itertuples()):
        for direction,offset in [('higher',.2),('lower',-.2)]:
            found=getattr(row,direction+'_found');recurs=getattr(row,direction+'_recurs')
            ax.text(max(found,recurs)+.25,yi+offset,f'{found} / {recurs}',va='center',fontsize=6)
    model={'gpt-6-astra':'GPT-6 Astra','gpt-6-sol':'GPT-6 Sol','gpt-6-luna':'GPT-6 Luna','gpt-5.5':'GPT-5.5','claude-opus-5-5':'Opus 5.5'}
    ax.set_yticks(y,[r.framework+' + '+model[r.model] for r in f.itertuples()]);ax.set_xlim(0,32);ax.set_xticks([0,10,20,30])
    ax.set_xlabel('Ranked comparisons (of 30)');ax.set_title('Tables 1 and 2 reanalyzed',pad=8)
    ax.legend(loc='upper center',bbox_to_anchor=(.5,-.31),ncol=4,frameon=False,fontsize=7,
              handlelength=.9,columnspacing=.8,handletextpad=.3)
    grid(ax);save_panel(fig,ax,name,'f')
    assemble(name,rects,h+426,13,printed_letters)
    return rects


def ed24():
    name='ed24';out=prepare(name)
    stats=pd.read_csv(out/'example_statistics.csv');values=pd.read_csv(out/'plotted_sample_values.csv')
    # July ED16 uses a 3+2 panel grammar. Each independent panel has its
    # printed 1:1.1-ish proportion; the final row is centred and stays this size.
    width=W/3; height=width*1.23
    rects={k:(i*width,0,(i+1)*width,height) for i,k in enumerate('abc')}
    rects.update({k:((i+.5)*width,height+10,(i+1.5)*width,2*height+10) for i,k in enumerate('de')})
    titles=['VEXAS syndrome','Critical limb\nischemia','Ataxia-telangiectasia','Gonadotropin\nstimulation','Varicose veins']
    cases=['VEXAS','CLI','A-T','Stimulated','Varicose']
    unadjusted=values[values.panel.isin(['a','b','d'])].clock_value
    for i,(letter,row) in enumerate(zip('abcde',stats.itertuples())):
        fig,ax=panel(name,letter,width,height,[.27,.24,.68,.48])
        sub=values[values.panel==letter]
        data=[sub[sub.group==group].clock_value.to_numpy() for group in ['control','case']]
        assert [len(v) for v in data]==[row.n_control,row.n_case]
        assert abs(data[1].mean()-data[0].mean()-row.mean_diff)<1e-6
        # Exact July sample-example median/whisker/cap style; Nov5 cosmetic
        # jitter is retained verbatim rather than drawing new random values.
        bp=ax.boxplot(data,positions=[0,1],widths=.24/.34,patch_artist=True,showfliers=False,
                      medianprops=dict(color='black',lw=1),whiskerprops=dict(lw=.8,color='#181114'),
                      capprops=dict(lw=.8,color='#181114'),boxprops=dict(lw=.8,edgecolor='#181114'))
        for box,color in zip(bp['boxes'],['#cccccc',ORANGE]):box.set_facecolor(color);box.set_alpha(.70)
        marker='^' if letter=='e' else 's'
        for group,color in [('control','#666666'),('case',ORANGE)]:
            ss=sub[sub.group==group]
            ax.scatter(ss.display_x,ss.clock_value,s=11,marker=marker,facecolors=color,edgecolors='black',
                       linewidths=.4,alpha=.85,zorder=4,clip_on=False)
        scale=unadjusted if letter in 'abd' else sub.clock_value
        lo,hi=scale.min(),scale.max();span=hi-lo
        ax.set_ylim(lo-.09*span,hi+.42*span);ax.set_xlim(-.55,1.55)
        control='Natural\ncycle' if letter=='d' else 'Control'
        ax.set_xticks([0,1],[f'{control}\nn = {row.n_control}',f'{cases[i]}\nn = {row.n_case}'])
        ax.set_ylabel('Horvath clock (years)' if letter=='e' else 'Clock output\nadjusted for age' if letter=='c' else 'Clock output')
        ax.yaxis.set_major_locator(MaxNLocator(4));grid(ax,'y')
        # Split title, tissue and accession are the original panel's content.
        ax.set_title(titles[i]+'\n'+row.tissue,pad=27,fontsize=8.5)
        ax.text(0,1.08,row.gse,transform=ax.transAxes,fontsize=7,ha='left',va='bottom')
        local_hi=max(v.max() for v in data);bracket=local_hi+.105*span
        ax.plot([0,0,1,1],[bracket-.025*span,bracket,bracket,bracket-.025*span],color='black',lw=.75)
        # Stat block is anchored to the bracket, as in the July sample plots.
        fdr=sci(row.fdr)
        ax.text(.5,bracket+.025*span,f'g = {row.hedges_g:.2f}\nFDR = {fdr}',ha='center',va='bottom',fontsize=6.5,linespacing=1.2)
        save_panel(fig,ax,name,letter)
    assemble(name,rects,2*height+10,22)
    return rects


def provenance(name,rects):
    out=NEW/name
    inputs=list((OLD/name).glob('*.csv'))+[OLD/name/'numbers.json',OLD/name/'provenance.json']
    entries=[]
    for path in inputs+list(STYLE.values())+[JULY,Path(__file__)]:
        checked(path)
        entries.append({'path':str(path),'sha256':sha(path),
                        'modified_utc':datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat()})
    (out/'build_provenance.json').write_text(json.dumps({'files':entries,'panel_rects_pt':rects,
                    'data_changes':[]},indent=2)+'\n')
    (out/'printed_text.json').write_text(json.dumps({k:v for k,v in PRINTED.items() if k.startswith(name)},indent=2)+'\n')
    (out/'preflight.txt').write_text('\n\n'.join(k+'\n'+v for k,v in PREFLIGHT.items() if k.startswith(name)))
    oldscript=REPO/f'scripts/figures_2026_11_05/restore_{name}_2026_11_05.py'
    diff=''.join(difflib.unified_diff(oldscript.read_text().splitlines(True),Path(__file__).read_text().splitlines(True),
                                    fromfile=str(oldscript.relative_to(REPO)),tofile=str(Path(__file__).relative_to(REPO))))
    (out/f'{name}_plotting.diff').write_text(diff)


def validate_sources():
    out=NEW/'ed23';n=json.loads((out/'numbers.json').read_text())
    a=pd.read_csv(out/'panel_a_funnel.csv')
    assert dict(zip(a.reason,a.n))==n['a'] and int(a.n.sum())==n['ground_truth_total']
    c=pd.read_csv(out/'panel_c_d_models.csv')
    for (_,row),source in zip(c.iterrows(),n['c_d'].values()):
        assert [row.lower_n,row.lower_total]==source['published_lower_identified']
        assert [row.higher_n,row.higher_total]==source['published_higher_identified']
        assert np.isclose(row.balanced_accuracy,source['call_balanced_accuracy'],rtol=0,atol=1e-15)
        assert ast.literal_eval(row.interval)==source['call_balanced_accuracy_interval']
    e=pd.read_csv(out/'panel_e_agreement.csv')
    for (_,row),source in zip(e.iterrows(),n['e']):
        assert row.system==source['label'] and row.n==source['call_n']
        assert np.isclose(row.agreement,source['call_agreement'],rtol=0,atol=1e-15)
        assert ast.literal_eval(row.interval)==source['call_agreement_interval']
    f=pd.read_csv(out/'panel_f_ranked_counts.csv')
    for (_,row),source in zip(f.iterrows(),n['f']):
        assert row.framework==source['framework'] and row.model==source['model']
        for direction in ['higher','lower']:
            for metric in ['found','recurs']:
                assert row[f'{direction}_{metric}']==source[direction][metric]
    (out/'validation.json').write_text(json.dumps({'all_source_values_match':True,
       'coverage_rows':len(a),'system_rows':len(c),'agreement_rows':len(e),
       'ranked_configurations':len(f),'numeric_changes':[]},indent=2)+'\n')
    out=NEW/'ed24';n=json.loads((out/'numbers.json').read_text())
    stats=pd.read_csv(out/'example_statistics.csv');values=pd.read_csv(out/'plotted_sample_values.csv')
    for (_,row),source in zip(stats.iterrows(),n['examples']):
        assert row.panel==source['panel'] and row.gse==source['gse']
        for field in ['n_case','n_control','hedges_g','fdr','mean_diff']:
            assert np.isclose(row[field],source[field],rtol=0,atol=1e-14)
        sub=values[values.panel==row.panel]
        assert len(sub[sub.group=='control'])==row.n_control and len(sub[sub.group=='case'])==row.n_case
        delta=sub[sub.group=='case'].clock_value.mean()-sub[sub.group=='control'].clock_value.mean()
        assert abs(delta-row.mean_diff)<1e-6
    (out/'validation.json').write_text(json.dumps({'all_source_values_match':True,
       'examples':len(stats),'plotted_observations':len(values),
       'all_mean_differences_recovered':True,'numeric_changes':[]},indent=2)+'\n')


def readmes():
    shared='''## Style sources

- `scripts/legacy_sensitivity_analyses_2026_05_14/external_db_crossvalidation_2026-05-14.py`, lines 430–438 supplies Calibri, base 9 pt, 0.6 pt axes and tick strokes, and absent top/right spines. Its compact titles and horizontal rows supply the ED23 comparison panels.
- `Aging-AI-Agent/src/visualize_validation_results.py` is the July ED7 validation generator. Its bars use black 0.5 pt outlines, alpha 0.8, `#4C72B0`/`#DD8452` and the named correctness/error colours. Its 14 × 10 in composition and July Figure3b–e's 7 × 5 in subplot proportions determine ED23 half-panel proportions. Final text sizes are set at print size to keep tick labels at 7.5 pt.
- `Aging-AI-Agent/create_combined_modality_top_hits_v2.py` supplies `#D55E00` for higher and `#0072B2` for lower, `#111` 0.5 pt bar edges, substantial outlined markers, square human RNA-seq and triangle human DNAm markers, and bold centred panel titles. These exact directional colours replace only the previous rendering's category colours.
- `paper/legacy_reference/figure_rev1/make_ed16_pooled_volcano.py` supplies the five-panel 3+2 layout for ED24; its printed proportions are retained with an empty position removed by centring the last two panels.
- `paper/legacy_reference/figure_rev1/sample_level_examples_2026-06-13/make_example_figures.py` supplies the July sample-example generator: quartile boxes at alpha 0.70, black 1 pt medians, 0.8 pt whiskers and caps, and black-edged 11 pt² alpha-0.85 sample marks. The July 0.24 box width and 0.34 arm separation are normalized to a 1-unit arm separation while preserving their ratio; sample x coordinates retain November 5 jitter.
- `paper/legacy_reference/figure_rev1/sensitivity_analyses/effect_size_distribution_2026-05-14.py` supplies compact 8 pt ticks, 9 pt axis-label proportions and top/right spine removal. ED14's generator was inspected as a volcano reference; no volcano is added to these pages.
- `paper/legacy_reference/figure_rev1/assemble_figure4_vector.py` supplies independent PDF placement in fixed PyMuPDF rectangles and bold Helvetica panel letters. The output is 183 mm wide; letters are 12 pt. `build_provenance.json` records every style-source hash, timestamp and panel rectangle.

## Verification

All November 5 source CSVs, `numbers.json` and their original provenance are copied byte-for-byte. `checks.json` checks those identities, page bounds, minimum text size and absence of embedded raster images. Every axis tick and annotation is listed in `printed_text.json`. The review render is 200 dpi. In the comparison PNG, the current figure is on the left and the July figure is on the right. ED23 uses package page 13 and ED24 uses page 22. The reference figure area is shown at the same physical scale; the PDF page heading is outside the cropped reference figure area.

The figure preflight report is in `preflight.txt`. Its perimeter and sparse-category proportion findings reflect the named July panel geometry; the explicit source reference governs those choices. The final vector PDFs, review renders and same-scale July comparisons were visually reviewed for overlaps and clipping. No unresolved text collision or clipping was found.

## Changes from November 5

No data, observations, intervals, sample counts, identifiers or numerical annotations change. The plots are independent vector panels assembled at fixed rectangles. Font, stroke, marker, palette, title and layout conventions now come from the July scripts. No footer, figure number or figure title appears in the image. `ed23_plotting.diff` or `ed24_plotting.diff` records the complete plotting-script difference against the November 5 builder. No manuscript or other figure is edited.
'''
    a=pd.read_csv(NEW/'ed23/panel_a_funnel.csv');c=pd.read_csv(NEW/'ed23/panel_c_d_models.csv');e=pd.read_csv(NEW/'ed23/panel_e_agreement.csv');f=pd.read_csv(NEW/'ed23/panel_f_ranked_counts.csv')
    n=json.loads((NEW/'ed23/numbers.json').read_text())
    b=pd.DataFrame([{'test':key,**v} for key,v in n['b'].items() if key in ['S3','S3c']]).fillna(0)
    text='''# Benchmark performance

ED23 retains the November 5 panels and their system order. The drawn panels c and d (panel_c.pdf, panel_d.pdf) are placed in Figure 4h,i; the page prints a, b and the drawn e and f as c and d. The style comparison uses July ED7 on package page 13.

## Panel sources

| Panel | Data source | Content |
| --- | --- | --- |
| a | `panel_a_funnel.csv`, `reason` and `n` | Nine ground-truth coverage categories |
| b | `numbers.json`, `b.S3` and `b.S3c` | Welch-test and released-test outcomes |
| c | `panel_c_d_models.csv`, lower/higher identified counts and totals | Seven systems and both directions |
| d | `panel_c_d_models.csv`, `balanced_accuracy` and `interval` | The same systems with source intervals |
| e | `panel_e_agreement.csv`, `agreement`, `interval` and `n` | Thirteen systems on 61 benchmark series |
| f | `panel_f_ranked_counts.csv`, all four found/recurs columns | Six rerun configurations |

## Printed numbers

Panel a's title prints 580 from `numbers.json.ground_truth_total`; its bar counts are the `n` column below.

'''+a.to_markdown(index=False)+'''

Panel b prints counts of at least 3 inside their corresponding segments; all five outcomes remain drawn. Exact counts, including the narrow unlabelled 1/2 segments, follow.

'''+b.to_markdown(index=False)+'''

Panel c prints the direction denominators 16 and 36 from `numbers.json.c_direction_denominators`. Its values are lower_n/lower_total × 100 and higher_n/higher_total × 100. Panel d draws the supplied estimates and intervals without recalculation; the dashed 1/3 line is three-class chance. The following source table holds every plotted estimate.

'''+c.to_markdown(index=False)+'''

Panel e's title prints 61 from `numbers.json.e_series`. Its right-hand n column prints `panel_e_agreement.csv.n` in the displayed system order; agreement and intervals are multiplied by 100.

'''+e.to_markdown(index=False)+'''

Panel f prints each direction's counts as found / recurs next to its paired bars, using every value in the following table. The denominator 30 in the x-axis label is the same ranked-comparison denominator used on November 5. Tables 1 and 2 in the panel title identify the manuscript tables.

'''+f.to_markdown(index=False)+'\n\n'+shared
    (NEW/'ed23/README.md').write_text(text)
    st=pd.read_csv(NEW/'ed24/example_statistics.csv')
    printed=st[['panel','gse','n_case','n_control','hedges_g','fdr','mean_diff','method']].copy()
    printed['printed_g']=printed.hedges_g.map(lambda x:f'{x:.2f}')
    printed['printed_FDR']=printed.fdr.map(sci)
    text='''# Human intervention examples

ED24 retains the five November 5 human examples, observations and source statistics. The style comparison uses the July five-panel ED16 layout on package page 22.

## Panel sources

| Panel | Comparison | Sample-value source |
| --- | --- | --- |
| a | VEXAS syndrome in whole blood | `plotted_sample_values.csv`, panel a |
| b | Critical limb ischaemia in gastrocnemius | `plotted_sample_values.csv`, panel b |
| c | Ataxia-telangiectasia in cerebellar vermis | `plotted_sample_values.csv`, panel c |
| d | Gonadotropin stimulation in granulosa cells | `plotted_sample_values.csv`, panel d |
| e | Varicose veins in vein wall | `plotted_sample_values.csv`, panel e |

The observation y values and cosmetic x positions are read directly from `clock_value` and `display_x`, including November 5's jitter. Panels a and c use individual means; c additionally adjusts for chronological age. Other panels use the same sample values. All five case-minus-control means reproduce `example_statistics.csv.mean_diff` within 10⁻⁶.

## Printed numbers

Each accession, n, g and FDR comes from the same row of `example_statistics.csv`; its full values and printed formatting follow. `n_control` and `n_case` are printed under the corresponding arms. The source effects and FDR are preserved; no new tests are computed. Numeric parts of GEO accessions are identifiers. Tick positions are scale values, listed in `printed_text.json`.

'''+printed.to_markdown(index=False)+'\n\n'+shared
    (NEW/'ed24/README.md').write_text(text)


def main():
    try:
        sys.path.insert(0,'pipeline_assets/figure_style_2026_09_10')
        style()
        rect23=ed23();rect24=ed24()
        validate_sources()
        provenance('ed23',rect23);provenance('ed24',rect24)
        readmes()
        print('Built ED23 and ED24; all numerical sources unchanged.')
    except Exception as exc:
        print(f'ED23/24 build failed: {type(exc).__name__}: {exc}',file=sys.stderr)
        raise


if __name__=='__main__':
    main()
