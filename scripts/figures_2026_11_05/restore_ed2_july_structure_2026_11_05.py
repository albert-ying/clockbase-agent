"""Restore ED2's July square correlation matrix and labelled clock PCA.

Style source: src/historical_geo_pipeline/analysis/1_dnam_analysis.R,
lines 195–216 and 273–296; July package page 8. A dated byte-for-byte
copy of that R source accompanies this bounded Python port because the
installed Rscript executable cannot run on this Mac (bad CPU type).
Tokens: July Helvetica/Arial, black axes, blue-positive/red-negative RdBu,
outlined squares with area |rho|, pals::stepped fills extracted from July PDF.
Additional panels retain October's data and plot types, in July typography.
Run after sourcing ../env_local.sh with ../.venv-pipeline/bin/python.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymupdf
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle
from PIL import Image, ImageDraw
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'paper/revision_2026-09-07/evidence/figure-fixes-2026-11-17/ed2'
BASE = ROOT / 'paper/revision_2026-09-07'
sys.path.insert(0, str(ROOT))
sys.path.insert(0, 'pipeline_assets/figure_style_2026_09_10')
from audit_system.blood_clock_structure_2026_09_30 import CLOCKS, age_residuals, pooled_within_study
from nature_figure_utils import smart_annotate, figure_preflight

plt.rcParams.update({'font.family': 'Arial', 'font.size': 7,
    'axes.labelsize': 8, 'axes.titlesize': 8, 'xtick.labelsize': 7,
    'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.linewidth': .65,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.facecolor': 'white'})


def receipt(path):
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'mtime_utc': datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()}


def clean_axes(ax):
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    ax.tick_params(width=.65, length=3)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    run = Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'r6-6-blood'
    cohort_path = run / 'structure-002/cohort.csv'
    grim_path = run / 'grimage-age-input-002/grimage_age_input.csv'
    summary_path = run / 'structure-002/summary.json'
    saved_path = run / 'structure-002/correlation_repaired_residual_spearman_chronological_age_input.csv'
    old_cor = BASE / 'evidence/blood-source-audit-2026-09-18/correlation_historical_gsm_global_iqr.csv'
    july = BASE / 'baseline/03_Figures_2026-07-06.pdf'
    old_pca = ROOT / 'paper/legacy_reference/plot/human_blood_clock_pca.pdf'
    r_source = ROOT / 'src/historical_geo_pipeline/analysis/1_dnam_analysis.R'
    for p in [cohort_path, grim_path, summary_path, saved_path, old_cor, july, old_pca, r_source]:
        if not p.exists():
            raise FileNotFoundError(f'Cannot proceed - required data file not found: {p}')
    copy_path = Path(__file__).with_name('ed2_july_dnam_analysis_source_copy_2026_11_05.R')
    if not copy_path.exists():
        shutil.copyfile(r_source, copy_path)
    elif copy_path.read_bytes() != r_source.read_bytes():
        raise RuntimeError('Existing July source copy differs; refusing overwrite')
    cohort = pd.read_csv(cohort_path, dtype={'gsm_id': str})
    grim = pd.read_csv(grim_path, dtype={'gsm_id': str})
    frame = cohort.loc[cohort.source.eq('stream_repaired')].merge(
        grim.drop(columns=['age', 'Horvathv1', 'predicted_sex']),
        on=['gse_id', 'gsm_id'], validate='one_to_one')
    sizes = frame.groupby('gse_id').size()
    frame = frame[frame.gse_id.isin(sizes[sizes >= 20].index)].reset_index(drop=True)
    assert len(frame) == 18098 and frame.gse_id.nunique() == 48
    variants = {'historical_horvath_age_input': frame[CLOCKS],
        'sex_corrected_horvath_age_input': frame[CLOCKS].assign(
            GrimAgeV1=frame.GrimAgeV1_sex_corrected, GrimAgeV2=frame.GrimAgeV2_sex_corrected),
        'chronological_age_input': frame[CLOCKS].assign(
            GrimAgeV1=frame.GrimAgeV1_chronological, GrimAgeV2=frame.GrimAgeV2_chronological)}
    residuals, matrices = {}, {}
    saved_summary = json.loads(summary_path.read_text())['structures']['grimage_age_input_contrast']
    contrast = []
    for name, values in variants.items():
        residuals[name] = age_residuals(values.assign(gse_id=frame.gse_id, age=frame.age))
        matrices[name] = pooled_within_study(residuals[name], frame.gse_id, 'spearman')
        for clock in ['GrimAgeV1', 'GrimAgeV2']:
            rho = float(matrices[name].loc['Horvathv1', clock])
            assert abs(rho - saved_summary[name]['horvath_grimage_residual_r'][clock]) < 1e-12
            contrast.append({'age_input': name, 'clock': clock, 'rho': rho})
    corr = matrices['chronological_age_input']
    saved = pd.read_csv(saved_path, index_col=0)
    delta = float(np.nanmax(abs(corr - saved)))
    assert delta < 1e-12
    corr.to_csv(OUT / 'correlation_recorded_age.csv')
    contrast = pd.DataFrame(contrast)
    contrast.to_csv(OUT / 'grimage_input_correlations.csv', index=False)
    # R corrplot default hclust uses complete linkage on 1-r. Preserve that
    # July clustering rule; October's average linkage is not its visual order.
    tree = linkage(squareform(1-corr.to_numpy(), checks=False), method='complete')
    order = [CLOCKS[i] for i in leaves_list(tree)]
    # July prcomp: predictors are observations; samples are variables.
    # Center and scale each sample across the 17 clocks, then SVD.
    res = residuals['chronological_age_input']
    complete = res.dropna()
    x = complete.to_numpy().T
    centered = x-x.mean(axis=0)
    standardized = centered / centered.std(axis=0, ddof=1)
    u, singular, _ = np.linalg.svd(standardized, full_matrices=False)
    pc = u[:, :2]*singular[:2]
    explained = 100*singular[:2]**2/(singular**2).sum()
    pca = pd.DataFrame({'clock': CLOCKS, 'PC1': pc[:, 0], 'PC2': pc[:, 1]})
    pca.to_csv(OUT / 'clock_pca_coordinates.csv', index=False)
    with pymupdf.open(old_pca) as doc:
        fills = [d['fill'] for d in doc[0].get_drawings() if d.get('fill')]
    assert len(fills) == len(CLOCKS)
    colors = dict(zip(CLOCKS, fills))
    fig = plt.figure(figsize=(183/25.4, 9.65))
    # July's full matrix stays above its labelled PCA; additions sit below.
    a = fig.add_axes([.235, .505, .53, .395])
    b = fig.add_axes([.17, .273, .72, .195])
    c = fig.add_axes([.255, .067, .26, .128])
    d1 = fig.add_axes([.605, .067, .17, .128])
    d2 = fig.add_axes([.81, .067, .17, .128])
    for ax in [b, c, d1, d2]:
        clean_axes(ax)
    # Full corrplot square glyphs, outlined cells, blue positive and red negative.
    values = corr.loc[order, order].to_numpy()
    cm = plt.get_cmap('RdBu')
    norm = Normalize(-1, 1)
    for i in range(len(order)):
        for j in range(len(order)):
            a.add_patch(Rectangle((j-.5, i-.5), 1, 1, facecolor='white',
                edgecolor='#A9A9A9', linewidth=.35))
            side = np.sqrt(abs(values[i, j]))*.94
            a.add_patch(Rectangle((j-side/2, i-side/2), side, side,
                facecolor=cm(norm(values[i, j])), edgecolor='#333333', linewidth=.35))
    a.set(xlim=(-.5,16.5), ylim=(16.5,-.5), aspect='equal')
    a.set_xticks(range(17), ["PedBE" if clock == "PEDBE" else clock for clock in order], rotation=90, ha='left', va='center', rotation_mode='anchor')
    a.set_yticks(range(17), ["PedBE" if clock == "PEDBE" else clock for clock in order])
    a.tick_params(axis='x', top=True, labeltop=True, bottom=False, labelbottom=False, length=0, pad=3)
    a.tick_params(axis='y', length=0, pad=3)
    for spine in a.spines.values():
        spine.set_visible(False)
    cbax = fig.add_axes([.79, .515, .018, .375])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cm), cax=cbax,
        ticks=[-1,-.5,0,.5,1])
    cb.set_label('Spearman ρ', fontsize=8, labelpad=5)
    cb.outline.set_linewidth(.5)
    cb.solids.set_rasterized(False)
    b.scatter(pc[:,0],pc[:,1], s=17, c=[colors[k] for k in CLOCKS],
        edgecolors='#333333', linewidths=.5, zorder=3, clip_on=False)
    b.set(xlim=(-265,280),ylim=(-100,150),
        xlabel=f'PC1 ({explained[0]:.1f}%)', ylabel=f'PC2 ({explained[1]:.1f}%)')
    b.set_xticks([-200,-100,0,100,200]); b.set_yticks([-100,0,100])
    # The residual PCA puts nine clocks nearly on top of each other. Explicit
    # label positions preserve July's labelled points and short repel leaders;
    # no value or point location is changed to make labels fit.
    positions = {'PhenoAge':(-45,127),'Lin':(-12,88),
        'YingAdaptAge':(219,86),'YingDamAge':(-200,54),
        'Hannum':(50,49),'Horvathv1':(-64,25),'YingCausAge':(-65,62),
        'Horvathv2':(63,12),'GrimAgeV2':(-65,-4),'GrimAgeV1':(60,-16),
        'PEDBE':(66,-36),'Zhang_10':(-58,-35),
        'DunedinPoAm38':(133,-55),'DunedinPACE':(-108,-55),
        'DNAmTL':(57,-78),'Knight':(-12,-92),'LeeRobust':(-105,-78)}
    for clock, xy in zip(CLOCKS, pc):
        b.annotate("PedBE" if clock == "PEDBE" else clock,xy=xy,xytext=positions[clock],ha='center',va='center',
            fontsize=6.6,arrowprops={'arrowstyle':'-','color':'#555555','lw':.4,
            'shrinkA':2,'shrinkB':3},annotation_clip=False)
    # Added October panel c: two GrimAge versions under three input recipes.
    hist_color = '#D55E00'
    chrono_color = '#0072B2'
    keys = list(variants)
    y = np.array([2,1,0])
    for clock, marker, offset in [('GrimAgeV1','o',.10),('GrimAgeV2','s',-.10)]:
        vals = [matrices[k].loc['Horvathv1',clock] for k in keys]
        c.plot(vals,y+offset,color='#999999',lw=.65,zorder=1)
        c.scatter(vals,y+offset,s=21,marker=marker,c=[hist_color,hist_color,chrono_color],
            edgecolors='black',linewidths=.4,zorder=3,clip_on=False)
    c.set_yticks(y,['Horvath age\nAll male','Horvath age\nCoded sex','Recorded age\nCoded sex'])
    c.set(xlim=(0,.75),ylim=(-.45,2.45),xlabel='Residual ρ with Horvathv1')
    c.set_xticks([0,.2,.4,.6])
    handles=[plt.Line2D([],[],ls='',marker=m,mfc='white',mec='black',ms=4,label=l)
        for m,l in [('o','GrimAgeV1'),('s','GrimAgeV2')]]
    # Dedicated above-panel space holds the legend clear of all data.
    c.legend(handles=handles,loc='lower left',bbox_to_anchor=(-.02,1.02),
        frameon=False,ncol=2,columnspacing=.8,handletextpad=.3,borderaxespad=0)
    # Added October panel d: identical source selection and 99.9% display bound.
    histres = residuals['historical_horvath_age_input']
    scatter_values = pd.DataFrame({'Horvathv1':histres.Horvathv1,
        'GrimAgeV1':histres.GrimAgeV1, 'GrimAgeV1_chronological':res.GrimAgeV1})
    limit = float(np.nanpercentile(abs(scatter_values.to_numpy()),99.9))
    shown = scatter_values.abs().le(limit).all(axis=1)
    span = 1.08*limit
    for ax, name, col, title in [(d1,'historical_horvath_age_input',hist_color,'Horvath age'),
                                (d2,'chronological_age_input',chrono_color,'Recorded age')]:
        ax.scatter(histres.loc[shown,'Horvathv1'],residuals[name].loc[shown,'GrimAgeV1'],
            s=.7,c=col,alpha=.20,edgecolors='none',clip_on=False)
        ax.set(xlim=(-span,span),ylim=(-span,span))
        ax.set_xlabel('Horvath residual\n(years)',fontsize=7)
        ax.set_xticks([-20,0,20]);ax.set_yticks([-20,0,20])
        ax.set_title(title,loc='left',pad=4)
        ax.text(.05,.96,f'ρ = {matrices[name].loc["Horvathv1","GrimAgeV1"]:.2f}',
            transform=ax.transAxes,ha='left',va='top',fontsize=7)
    d1.set_ylabel('GrimAge residual (years)',fontsize=7,labelpad=2)
    d2.tick_params(labelleft=False)
    for label, x0, y0 in [('a',.035,.982),('b',.035,.485),('c',.035,.215),('d',.55,.215)]:
        fig.text(x0,y0,label,fontweight='bold',fontsize=10,ha='left',va='top')
    fig.canvas.draw()
    report = figure_preflight(fig,{'a':a,'b':b,'c':c,'d1':d1,'d2':d2},verbose=False,
        check_house_style=False,check_max_fontsize=False,forbid_uniform_grid=False,
        check_data_occupancy=False,check_matched_limits=False)
    (OUT/'skill_preflight.txt').write_text(report['text'])
    fig.savefig(OUT/'ed2.pdf')
    fig.savefig(OUT/'ed2_review.png',dpi=200)
    plt.close(fig)
    with pymupdf.open(july) as doc:
        page = doc[7]
        # Compare the plotted content at the same 183 mm width, with package
        # heading removed; never shrink the new figure to July page whitespace.
        ref = page.get_pixmap(matrix=pymupdf.Matrix(2.2,2.2),alpha=False)
        im = Image.frombytes('RGB',(ref.width,ref.height),ref.samples)
        mask = np.asarray(im).min(axis=2)<245
        mask[:int(95/792*im.height)] = False
        yy,xx=np.where(mask)
        im=im.crop((int(xx.min()),int(yy.min()),int(xx.max()+1),int(yy.max()+1)))
    new=Image.open(OUT/'ed2_review.png').convert('RGB')
    w=new.width
    im=im.resize((w,round(im.height*w/im.width)),Image.Resampling.LANCZOS)
    canvas=Image.new('RGB',(w*2+36,max(im.height,new.height)+40),'white')
    canvas.paste(im,(0,40));canvas.paste(new,(w+36,40))
    draw=ImageDraw.Draw(canvas);draw.text((12,10),'July reference',fill='black')
    draw.text((w+48,10),'Restored figure',fill='black')
    canvas.save(OUT/'ed2_vs_july.png')
    absent=[c for c in pd.read_csv(old_cor,index_col=0).columns if c not in CLOCKS]
    stats={'samples':len(frame),'studies':int(frame.gse_id.nunique()),
        'pca_samples':len(complete),'pca_predictors':17,'pca_variance_percent':explained.tolist(),
        'scatter_samples_drawn':int(shown.sum()),'scatter_samples_omitted':int((~shown).sum()),
        'scatter_99_9_percentile_years':limit,'scatter_axis_limit_years':span,
        'correlation_max_abs_difference_from_saved':delta,'heatmap_order':order,
        'absent_july_predictors':absent,'pca_colors':colors,
        'sources':[receipt(p) for p in [cohort_path,grim_path,summary_path,saved_path,old_cor,july,old_pca,r_source]],
        'script':receipt(Path(__file__))}
    (OUT/'numbers_and_sources.json').write_text(json.dumps(stats,indent=2)+'\n')
    with pymupdf.open(OUT/'ed2.pdf') as doc:
        stats['pdf_width_mm']=float(doc[0].rect.width/72*25.4)
        stats['pdf_height_mm']=float(doc[0].rect.height/72*25.4)
        stats['pdf_embedded_images']=len(doc[0].get_images())
    (OUT/'checks.json').write_text(json.dumps({k:v for k,v in stats.items() if k!='sources'},indent=2)+'\n')
    print(json.dumps({k:v for k,v in stats.items() if k not in ['sources','pca_colors','script']},indent=2))


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        raise SystemExit(f'ED2 restoration failed: {exc}') from exc
