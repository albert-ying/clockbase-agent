#!/usr/bin/env python3
"""Restore July Figure 3 using verified effects and statistical rankings.

Reference: July page 3 and copied July generators. RdBu_r volcano, category
annotation colours, black outlined modality shapes, orange/blue hatched bars.
Typography: July sans-serif appearance at 7 pt ticks/labels and 6 pt annotations.
Composition: July full-width volcano above Genetic/Drug then Disease/Environment.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib
import importlib.util
import json
import os
import sys
import textwrap
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import Normalize
from scipy.stats import gaussian_kde
import fitz
from PIL import Image,ImageOps,ImageDraw
ROOT=Path(__file__).resolve().parents[2]
HERE=Path(__file__).parent
OUT=ROOT/'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-05/figure3'
DATA=OUT.parent/'data/combined_3modality_plot_data_verified_2026-11-05.csv'
FINAL=Path(os.environ['CLOCKBASE_FINAL'])
RANK=FINAL/'comparison-ranking'
REFERENCE=ROOT/'paper/revision_2026-09-07/baseline/03_Figures_2026-07-06.pdf'

def load_module(name):
    spec=importlib.util.spec_from_file_location(name,HERE/(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def source(path):
    return {'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
            'mtime_utc':datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat()}

def save_comparison():
    july=fitz.open(REFERENCE)
    # The package page contains a heading and large page margins; the figure
    # content itself is compared at the same 183-mm physical width.
    pix=july[2].get_pixmap(dpi=200)
    original=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
    crop=original.crop((int(original.width*.13),int(original.height*.125),int(original.width*.84),int(original.height*.835)))
    current=Image.open(OUT/'figure3_review.png').convert('RGB')
    crop=crop.resize((current.width,round(crop.height*current.width/crop.width)),Image.Resampling.LANCZOS)
    canvas=Image.new('RGB',(current.width*2+40,max(crop.height,current.height)+50),'white')
    draw=ImageDraw.Draw(canvas);draw.text((10,10),'July reference',fill='black');draw.text((current.width+40,10),'Verified restoration',fill='black')
    canvas.paste(crop,(0,40));canvas.paste(current,(current.width+40,40));canvas.save(OUT/'figure3_vs_july.png',dpi=(200,200))
    original.save(OUT/'july_page3.png',dpi=(200,200))

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    volcano=load_module('create_combined_modality_volcano_verified_2026_11_05')
    bars=load_module('create_combined_modality_top_hits_v2_verified_2026_11_05')
    plt.rcdefaults()
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,'axes.titlesize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.55,'savefig.facecolor':'white','hatch.linewidth':.25})
    df=pd.read_csv(DATA,low_memory=False)
    rank=pd.read_csv(RANK/'ranked_comparisons_all.csv')
    top=pd.read_csv(RANK/'figure3_category_top15.csv')
    summary=json.loads((RANK/'summary.json').read_text())
    fig=plt.figure(figsize=(183/25.4,9.7))
    axes={}
    axes['a']=fig.add_axes([.080,.690,.865,.257])
    axes['a_full']=fig.add_axes([.080,.585,.29,.065])
    axes['f']=fig.add_axes([.690,.580,.255,.065])
    # Labels occupy the same left-margin role as in July, with GEO on a
    # second line where necessary. Separate positions avoid tight-layout
    # shrinking the data axes under long biological names.
    for name,rect in {'b':[.270,.320,.220,.210],'c':[.750,.320,.240,.210],
                      'd':[.270,.035,.220,.227],'e':[.750,.035,.240,.227]}.items():
        axes[name]=fig.add_axes(rect)
    annotations=volcano.draw_verified(fig,axes['a'],axes['a_full'],df,rank)
    bar_rows=bars.draw_verified(fig,{k:axes[k] for k in 'bcde'},top)
    funnel=summary['funnel']
    steps=[('FDR < 0.05',funnel['usable_significant_n3_rankable_category']['comparisons']),
           ('Agreed contrast',funnel['agreed_matching_contrast_type']['comparisons']),
           ('No design concern',funnel['no_agreed_adjudication_design_concern']['comparisons']),
           ('Design admitted',funnel['context_admitted']['comparisons']),
           ('Ranked',summary['outputs']['all']['ranked'])]
    ax=axes['f'];yp=np.arange(5)
    ax.barh(yp,[v for _,v in steps],color=bars.COLOR_PRO,height=.72,edgecolor='#111',lw=.35)
    ax.invert_yaxis();ax.set_xlim(0,8500);ax.set_ylim(4.65,-.65)
    ax.set_yticks(yp);ax.set_yticklabels([s for s,_ in steps],fontsize=7)
    ax.tick_params(axis='y',length=0,pad=3);ax.set_xticks([])
    for y,(_,n) in zip(yp,steps):ax.text(n+160,y,f'{n:,}',fontsize=6,ha='left',va='center')
    for spine in ax.spines.values():spine.set_visible(False)
    fig.text(.560,.650,'f',fontsize=11,fontweight='bold',ha='left')
    fig.text(.790,.650,'Ranking eligibility',fontsize=8,fontweight='bold',ha='center')
    mod_counts=df.modality.value_counts()
    fig.text(.08,.984,f'Rows supporting inference: mouse RNA-seq {mod_counts["mouse RNA-seq"]:,}  ·  human RNA-seq {mod_counts["human RNA-seq"]:,}  ·  human DNAm {mod_counts["human DNAm"]:,}',fontsize=7,va='top')
    fig.text(.012,.956,'a',fontsize=11,fontweight='bold')
    # Same direction and shape encodings as July. Fixed top placement is
    # outside the axes and cannot cover a point or a label.
    handles=[Patch(facecolor=bars.COLOR_ANTI,edgecolor='#111',lw=.35,label='Lower clock output'),Patch(facecolor=bars.COLOR_PRO,edgecolor='#111',lw=.35,label='Higher clock output')]
    handles += [Line2D([0],[0],marker=volcano.MODALITY_MARKERS[m],color='#222',markerfacecolor='#eee',markersize=3.8,ls='None',label=m) for m in volcano.MODALITY_ORDER]
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.53,.974),frameon=False,ncol=5,columnspacing=.9,handlelength=1.2,handletextpad=.4)
    fig.canvas.draw()
    # Only requested defect checks apply. The named July composition takes
    # precedence over the skill's alternative palette and layout rules.
    sys.path.insert(0,'pipeline_assets/figure_style_2026_09_10')
    from nature_figure_utils import strict_figure_qc
    qc=strict_figure_qc(fig,axes,min_fontsize=7,forbid_uniform_grid=False,check_data_occupancy=False,check_house_style=False,check_max_fontsize=False,check_matched_limits=False,check_matrix_row_labels=False,check_text_figure_edges=False,max_tall_panel_aspect=10,max_sparse_panel_aspect=10)
    qc['requested_defect_issues']=[i for i in qc.get('issues',[]) if i['code'] in ['overlap','font-size','text-clipping','missing-colorbar','text-spine']]
    qc['requested_defect_checks_passed']=not qc['requested_defect_issues']
    qc['scope_note']='Per brief, clip_on and baseline point-edge findings are not style instructions. The main volcano is a view window and the full-range view retains all measurements.'
    (OUT/'nature_figure_defect_checks.json').write_text(json.dumps(qc,indent=2,default=str)+'\n')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    canvas_issues=[]
    for artist in fig.findobj(matplotlib.text.Text):
        if artist.get_visible() and artist.get_text():
            box=artist.get_window_extent(renderer)
            if box.x0 < fig.bbox.x0-1 or box.y0 < fig.bbox.y0-1 or box.x1 > fig.bbox.x1+1 or box.y1 > fig.bbox.y1+1:
                canvas_issues.append(artist.get_text())
    (OUT/'canvas_edge_checks.json').write_text(json.dumps({'outside_canvas':canvas_issues},indent=2)+'\n')
    assert not canvas_issues, f'Figure 3 text extends beyond canvas: {canvas_issues}'
    assert qc['requested_defect_checks_passed'], f'Figure 3 has requested defect failures: {qc["requested_defect_issues"]}'
    fig.savefig(OUT/'figure3.pdf')
    fig.savefig(OUT/'figure3_review.png',dpi=200)
    # PDF readback verifies vector text, physical size and embedded image count.
    pdf=fitz.open(OUT/'figure3.pdf');p=pdf[0]
    spans=[s for b in p.get_text('dict')['blocks'] if 'lines'in b for l in b['lines'] for s in l['spans']]
    checks={'width_mm':p.rect.width/72*25.4,'height_mm':p.rect.height/72*25.4,'embedded_raster_images':len(p.get_images()),'fonts':len(p.get_fonts()),'minimum_text_pt':min(s['size'] for s in spans),'text_spans':len(spans),'nature_qc':qc,'requested_visual_checks_passed':qc['requested_defect_checks_passed'] and not canvas_issues,'all_bar_rows':len(bar_rows),'main_view_rows':int(df.effect_size.between(-45,45).mul(df.neg_log10_fdr.between(0,12)).sum()),'full_range_rows':len(df)}
    (OUT/'checks.json').write_text(json.dumps(checks,indent=2,default=str)+'\n')
    annotations=annotations.drop(columns=['effect_size','pvalue','fdr','direction','final_score'],errors='ignore')
    annotations.to_csv(OUT/'panel_a_annotations.csv',index=False)
    bar_rows=bar_rows.drop(columns=['effect_size','pvalue','fdr','direction','final_score'],errors='ignore')
    bar_rows.to_csv(OUT/'panels_b_to_e_values.csv',index=False)
    pd.DataFrame(steps,columns=['stage','comparisons']).to_csv(OUT/'panel_f_values.csv',index=False)
    provenance={'created_utc':datetime.now(timezone.utc).isoformat(),'sources':[source(p) for p in [DATA,RANK/'figure3_category_top15.csv',RANK/'ranked_comparisons_all.csv',RANK/'summary.json',REFERENCE,Path(__file__),HERE/'create_combined_modality_volcano_verified_2026_11_05.py',HERE/'create_combined_modality_top_hits_v2_verified_2026_11_05.py']],'modality_rows':mod_counts.to_dict(),'significant_rows':int(df.is_significant.sum()),'annotation_count':len(annotations),'uniform_marker_area_pt2':8,'main_view':{'x':[-45,45],'y':[0,12]},'full_range':{'x':[-75,140],'y':[0,235]}}
    (OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    save_comparison()
    write_readme(df,top,steps,annotations,checks)
    print(json.dumps({k:v for k,v in checks.items() if k!='nature_qc'},indent=2));print('QC',qc.get('passed'),len(qc.get('issues',[])))

def write_readme(df,top,steps,annotations,checks):
    rows=['# Figure 3 restoration','','The figure restores the July volcano and category bars using verified atlas statistics and the released comparison ranking. All panels are vector artwork at 183 mm width and 246.38 mm height.','','## Panel map','','| July panel | Restored panel | Change |','|---|---|---|','| a | a | Modality shapes, RdBu_r effect colours, density contours, category-coloured leader labels and FDR threshold retained. Verified Cohen’s d replaces the historical effect. Labels use statistical rank. Equal marker area replaces unavailable composite scores. A small full-range view preserves every row while the main view retains July’s readable scale. |','| b | b | Genetic top 15 ordered by the released CI-bound rank; bars show Hedges’ g with capless adjusted exact confidence limits. |','| c | c | Drug top 15 with the same changes. |','| d | d | Disease top 15 with the same changes. |','| e | e | Environment top 15 with the same changes. |','| none | f | Horizontal eligibility bars added below a beside its full-range view. |','','## Data sources','','- Panel a uses `../data/combined_3modality_plot_data_verified_2026-11-05.csv` and `ranked_comparisons_all.csv` in CLOCKBASE_FINAL/comparison-ranking.','- Panels b–e use CLOCKBASE_FINAL/comparison-ranking/figure3_category_top15.csv, without relevance filtering or a new score.','- Panel f uses CLOCKBASE_FINAL/comparison-ranking/summary.json.','- Absolute paths, source modification times and SHA-256 hashes are recorded in provenance.json.','','## Printed numbers','','The header reports 26,999 mouse RNA-seq rows, 24,152 human RNA-seq rows and 12,903 human DNAm rows, derived by grouping the verified plotting table. Together these are 64,054 rows from 53,865 distinct comparisons. There are 10,243 rows with BH-FDR < 0.05. The FDR threshold line is at −log10(0.05).','','The full-range view contains all 64,054 rows. The main view covers Cohen’s d from −45 to 45 and −log10(BH-FDR) from 0 to 12; '+f'{checks["main_view_rows"]:,} rows lie within both limits. Values outside these limits remain visible in the full-range view. Effect and FDR coordinates retain their source values. Colour limits are −45 and 45, with extensions marking values outside the colour range. The effect scale prints −40, 0 and 40; density bounds are the rounded first and last KDE contour levels. panel_a_kde.npz records the evaluated grid, Scott bandwidth and source-array SHA-256 digest.','','Each category contains 15 comparisons, ordered by descending rank_bound (the adjusted exact CI bound nearest zero). Full labels, GEO accessions, effects, CI limits, P values, FDR values, sample counts and multiplicity factors are in panels_b_to_e_values.csv. The 60 selected comparisons all belong to mouse RNA-seq; marker shapes retain the July convention. Labels are line-wrapped rather than truncated.','','| Eligibility stage | Comparisons |','|---|---:|']
    rows.extend(f'| {label} | {n:,} |'for label,n in steps)
    rows+=['','These stages use distinct comparisons through design admission; the final 946 are deduplicated study-and-label rankings after adjusted-CI exclusion of zero. The source additionally records 1,133 rows whose adjusted CI excludes zero, before deduplication.','','## Draft panel legends','','**a.** Combined volcano of 64,054 VERIFIED or PROVISIONAL atlas rows from 53,865 comparisons. The horizontal axis shows case-minus-control Cohen’s d and the vertical axis shows −log10 of BH-FDR adjusted within modality. Circles, squares and triangles denote mouse RNA-seq, human RNA-seq and human DNAm, respectively. Significant rows are coloured by Cohen’s d using the blue-to-red scale; other rows are pale grey. The dashed line marks FDR = 0.05. Contours show the pooled two-dimensional kernel density. Category-coloured labels identify the three highest-ranked raising and three highest-ranked lowering comparisons per category. The small view displays the complete effect and FDR ranges. All marker areas are equal.','','**b–e.** The top 15 Genetic (b), Drug (c), Disease (d) and Environment (e) comparisons ranked by the adjusted exact confidence bound nearest zero. Bars show Hedges’ g and thin horizontal lines show exact 95% confidence intervals adjusted for the number of eligible tests sharing a study and intervention label. Orange bars indicate higher clock output and blue bars indicate lower clock output. Marker shape and hatch identify modality using the convention in a. Labels show the intervention and GEO accession. Sample counts and test methods are provided in the source table and row provenance.','','**f.** Eligibility of comparisons for statistical ranking. Counts show comparisons with FDR < 0.05, at least three samples per arm and a ranked category; agreement of two models on contrast type; absence of an agreed design concern; admission by the design assessment; and the final study-and-label ranking after adjusted confidence intervals and deduplication.','','## Validation','','The build checks the exact row count and schema in the shared export, reads the released CI limits without recomputation, records every selected row, reads back the vector PDF, and runs nature-figure defect checks without applying a replacement style or layout. The PDF has '+str(checks['embedded_raster_images'])+' embedded raster images and '+str(checks['fonts'])+' embedded fonts. Tick and axis text use at least 7 pt; in-panel annotations use at least 6 pt. The volcano carries a Cohen’s d colour scale and the density contours carry a separate density scale. Render and same-width July comparison are supplied.','','## Open issues','','The raw nature-figure report also flags point clipping and proximity to the zero-FDR baseline. These findings concern the documented view window and bounded FDR values; the requested text, scale and clipping checks are reported separately in nature_figure_defect_checks.json and canvas_edge_checks.json.','','The historical composite score has no verified counterpart. The plot therefore uses equal marker areas and statistical ranks. The main volcano uses a readable window because the verified atlas extends to Cohen’s d = 127.54 and −log10(FDR) = 221.59; the adjacent full-range view retains those measurements. Source values, not these rounded extrema, are plotted. The ranking export also contains historical effect_size and fdr aliases; these are excluded from the figure sidecars and are never used for volcano coordinates. panel_a_annotations.csv records verified cohen_d and bh_fdr after matching the atlas by modality, GEO accession, comparison, condition, clock and verified P value.']
    (OUT/'README.md').write_text('\n'.join(rows)+'\n')

if __name__=='__main__':
    try:main()
    except (OSError,ValueError,AssertionError) as exc:
        raise SystemExit(f'Figure 3 restoration failed: {exc}') from exc
