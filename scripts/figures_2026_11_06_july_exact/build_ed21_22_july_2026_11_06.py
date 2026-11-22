#!/usr/bin/env python3
"""Rebuild ED21/22 from Nov 5 CSVs as independent July-style vector panels.

Style sources: July ED7 visualize_validation_results.py (Arial, #4C72B0,
#DD8452, outlined .5pt bars); ED19 external_db_crossvalidation (Calibri 9,
.6pt spines); Figure3 top_hits_v2 (orange/blue, end markers, dotted grid);
combined_modality_volcano (RdBu_r, black marker outline, category labels).
Panel geometry follows native July ED19 / 7x5-inch top-hits panel proportions,
adapted for the requested 7pt final tick floor. Fixed rect assembly follows
assemble_figure4_vector.py; no panel PDF is scaled during assembly.
"""
from pathlib import Path
import hashlib
import json
import sys
import re
import os
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.text import Text
import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image, ImageDraw

REPO=Path(__file__).resolve().parents[2]
BASE=REPO/'paper/revision_2026-09-07'
OLD=BASE/'evidence/figures-rebuild-2026-11-05'
OUT=BASE/'evidence/figure-fixes-2026-11-17'
JULY=BASE/'baseline/03_Figures_2026-07-06.pdf'
AGENT=REPO.parents[1]/'Aging-AI-Agent'
SKILL=Path('pipeline_assets/figure_style_2026_09_10')
sys.path.insert(0,str(SKILL))
from nature_figure_utils import strict_figure_qc, figure_preflight

INK='#222222'; HIGH='#D55E00'; LOW='#0072B2'
CLASSIC=['#4C72B0','#DD8452','#55A868','#8172B3','#CCB974','#64B5CD']
MODS={'mm_rna':('Mouse RNA-seq','o','#4C72B0'),
      'hs_rna':('Human RNA-seq','s','#C44E52'),
      'dnam':('Human DNAm','^','#55A868')}
CATEGORY={'Disease':'#DC143C','Drug':'#4B0082','Genetic':'#228B22','Environment':'#FF8C00','Other':'#808080'}
PAGE_W=595.68
SOURCE_FILES=[]
CHECKS={}
DEFECTS = [
    ('STORED_P_UNUSABLE', 'Stored P missing, 0 or 1'),
    ('ARMS_INVERTED_CONFIRMED', 'Arms swapped in two models'),
    ('ARMS_DEFINED_BY_CLOCK_OUTPUT', 'Arms set by clock outputs'),
    ('ARMS_DO_NOT_MATCH_GEO_SAMPLE_RECORDS', 'Arms contradict sample records'),
    ('CLOCK_OUTPUTS_OUTSIDE_PLAUSIBLE_RANGE', 'Implausible clock outputs'),
    ('SINGLE_CELLS_AS_INDEPENDENT_SAMPLES', 'Cells counted as individuals'),
    ('EFFECT_NOT_REPRODUCED', 'Effect not reproduced'),
    ('EMPTY_OR_ARTIFACT_INPUT', 'Empty or artifactual input'),
    ('TEST_NOT_A_GROUP_CONTRAST', 'Not a group contrast'),
    ('DETACHED_INPUT_PROFILES', 'Detached DNAm profiles'),
    ('NON_STANDARD_ASSAY_IN_ARMS', 'Non-bisulfite assay'),
    ('FDR_BELOW_RAW_P', 'FDR below raw P'),
]


def current_defects():
    final = Path(os.environ['CLOCKBASE_FINAL']).resolve()
    assert (final / 'reclassification_provenance.json').is_file(), final
    path = final / 'error-census/census.json'
    try:
        census = json.loads(path.read_text())
        provenance = json.loads((OLD / 'ed21/provenance.json').read_text())
        record = next(r for r in provenance['inputs'] if Path(r['path']).resolve() == path)
        assert record['sha256'] == hashlib.sha256(path.read_bytes()).hexdigest(), 'ED21 source census is stale'
    except (OSError, ValueError, StopIteration) as exc:
        raise RuntimeError(f'Cannot verify ED21 census {path}: {exc}') from exc
    SOURCE_FILES.append(path)
    rows = []
    for modality in ('dnam', 'hs_rna', 'mm_rna'):
        for code, _ in DEFECTS:
            value = census['modalities'][modality]['defects'][code]
            rows.append({'modality': modality, 'defect': code, **value,
                         'percent': 100 * value['n'] / value['of']})
    return pd.DataFrame(rows)

def read(path):
    path=Path(path)
    if not path.is_file(): raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}')
    SOURCE_FILES.append(path)
    try: return pd.read_csv(path)
    except Exception as exc: raise RuntimeError(f'Cannot read {path}: {exc}') from exc


def copy_file(source, target):
    try:
        subprocess.run(['/bin/cp', '-f', str(source), str(target)], check=True)
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError(f'Cannot copy {source} to {target}: {exc}') from exc

def setup(font):
    plt.rcdefaults()
    for name in ['Calibri.ttf','Calibrib.ttf','Calibrii.ttf','Calibriz.ttf','Arial.ttf','Arial Bold.ttf']:
        p=Path('/Applications/Microsoft Word.app/Contents/Resources/DFonts')/name
        if p.exists(): font_manager.fontManager.addfont(str(p))
    plt.rcParams.update({'font.family':font,'font.size':9,'axes.labelsize':8,
       'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
       'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,
       'xtick.major.width':.6,'ytick.major.width':.6,'xtick.major.size':3,
       'ytick.major.size':3,'pdf.fonttype':42,'ps.fonttype':42,
       'savefig.facecolor':'white','figure.facecolor':'white','axes.facecolor':'white'})

def panel(w,h,title,box=(.23,.23,.72,.60),bold=False):
    fig=plt.figure(figsize=(w/72,h/72)); ax=fig.add_axes(box)
    ax.set_title(title,fontsize=9,fontweight='bold' if bold else 'normal',pad=7)
    ax.set_axisbelow(True)
    return fig,ax

def grid(ax): ax.grid(axis='x',alpha=.3,linestyle=':',linewidth=.5)

def save(fig,out,letter):
    fig.canvas.draw()
    # July's fixed bar labels, outside legends and fixed bracket annotations are
    # intentional placements copied from the named generators. House composition
    # styling is advisory here because the user explicitly requires July style.
    rep=figure_preflight(fig,{'panel':fig.axes[0]},strict=False)
    sizes=[t.get_fontsize() for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    geometry_codes={'aspect','outer_margin_pct','sparse-panel-aspect'}
    unresolved=[item for item in rep['blockers'] if item['code'] not in geometry_codes]
    tick_sizes=[t.get_fontsize() for ax in fig.axes for t in ax.get_xticklabels()+ax.get_yticklabels() if t.get_visible() and t.get_text().strip()]
    CHECKS[f'{out.name}_{letter}']={'preflight':rep,'minimum_matplotlib_text_pt':min(sizes),
        'minimum_tick_size_pt':min(tick_sizes),'july_geometry_override_codes':sorted(geometry_codes),
        'passed_after_user_geometry_override':not unresolved}
    if min(sizes)<6: raise ValueError(f'Sub-floor annotation in {out.name} {letter}')
    if min(tick_sizes)<7: raise ValueError(f'Sub-floor tick in {out.name} {letter}')
    if unresolved: raise ValueError(f'Unresolved figure defects in {out.name} {letter}: {unresolved}')
    fig.savefig(out/f'panel_{letter}.pdf')
    plt.close(fig)

def stacks(ax,df,spec,labels,print_segments=True):
    groups=df.group.drop_duplicates().tolist(); yy=np.arange(len(groups))[::-1]
    for y,key in zip(yy,groups):
        sub=df[df.group==key].set_index('class'); left=0
        for code,title,color in spec:
            r=sub.loc[code]; width=float(r.percent)
            if r.n:
                ax.barh(y,width,left=left,height=.67,color=color,edgecolor='black',linewidth=.5,alpha=.8)
                if print_segments and width>=9:
                    ax.text(left+width/2,y,f'{int(r.n):,}',ha='center',va='center',fontsize=7,color='white' if code=='VERIFIED' else INK)
            left+=width
        ax.text(102,y,f'{int(sub.total.iloc[0]):,}',va='center',fontsize=7)
    ax.set_yticks(yy,labels); ax.set_ylim(-.55,len(groups)-.45)
    ax.set_xlim(0,116);ax.set_xticks([0,25,50,75,100]);ax.set_xlabel('Atlas rows (%)' if print_segments else 'Candidates (%)')
    ax.spines['bottom'].set_bounds(0,100);ax.tick_params(axis='y',length=0)

def build21():
    SOURCE_FILES.clear()
    setup('Arial'); out=OUT/'ed21';out.mkdir(parents=True,exist_ok=True)
    a=read(OLD/'ed21/panel_a_classification.csv');b=current_defects()
    c=read(OLD/'ed21/panel_c_published_fates.csv');d=read(OLD/'ed21/panel_d_sample_values.csv')
    e=read(OLD/'ed21/panel_e_histogram.csv');f=read(OLD/'figure3/panel_f_values.csv')
    classes=[('VERIFIED','Verified',CLASSIC[0]),('PROVISIONAL','Provisional','#9fc2df'),('DESCRIPTIVE','Descriptive',CLASSIC[2]),
       ('DUPLICATE','Duplicate',CLASSIC[3]),('PENDING','Withheld',CLASSIC[4]),('NOT_REPAIRABLE','Not repairable',CLASSIC[1])]
    fates=[('CONFIRMED_CLAIMED_DIRECTION','Confirmed',CLASSIC[0]),('NOT_SIGNIFICANT_AFTER_REPAIR','Not significant','#b9d3e7'),
       ('DESCRIPTIVE_ONLY_ARM_BELOW_3','Descriptive',CLASSIC[2]),('PENDING','Pending',CLASSIC[4]),
       ('CONTRADICTS_CLAIMED_DIRECTION','Reversed',HIGH),('NO_EXACT_ROW','No atlas row','white')]
    fig,ax=panel(PAGE_W,125,'Final classification',(.15,.26,.81,.51),True)
    stacks(ax,a,classes,['Human DNAm','Human RNA-seq','Mouse RNA-seq'])
    ax.legend(handles=[Patch(fc=co,ec='black',lw=.5,label=t) for _,t,co in classes],ncol=6,frameon=False,
       loc='upper center',bbox_to_anchor=(.5,1.48),columnspacing=1,handlelength=1)
    save(fig,out,'a')
    defect_height=265
    fig,ax=panel(355,defect_height,'Defects found by the verification',(.48,.21,.48,.65),True)
    codes=[code for code,_ in DEFECTS]; yy=np.arange(len(codes))[::-1]
    labels=[label for _,label in DEFECTS]
    for off,key in zip((.21,0,-.21),['mm_rna','hs_rna','dnam']):
        title,marker,color=MODS[key];sub=b[b.modality==key].set_index('defect').loc[codes]
        ax.scatter(sub.percent,yy+off,s=26,marker=marker,c=color,ec='black',lw=.4,clip_on=False,label=title)
    upper=max(10,float(np.ceil(b.percent.max()/2)*2))
    ax.set_yticks(yy,labels);ax.set_ylim(-.6,len(codes)-.4);ax.set_xlim(-.2,upper);ax.set_xticks(np.arange(0,upper+1,2));ax.set_xlabel('Rows with defect (%)');grid(ax)
    ax.tick_params(axis='y',length=0)
    ax.legend(ncol=3,frameon=False,loc='upper right',bbox_to_anchor=(1.01,-.23),columnspacing=.7,handletextpad=.3,handlelength=.8)
    save(fig,out,'b')
    fig,ax=panel(PAGE_W-355,225,'Candidates of the submitted version',(.39,.40,.55,.46),True)
    fate_labels={'Figure3_bars':'Suppl. Fig. 1 bars','S14':'Table S14','Table1_pro_aging':'Table S20, higher','Table2_anti_aging':'Table S20, lower'}
    stacks(ax,c,fates,[fate_labels[key] for key in c.group.drop_duplicates()],False);ax.set_xticks([0,50,100])
    ax.legend(handles=[Patch(fc=co,ec='black',lw=.5,label=t) for _,t,co in fates],ncol=2,frameon=False,
       loc='upper center',bbox_to_anchor=(.35,-.44),handlelength=.9,columnspacing=.7,labelspacing=.4)
    save(fig,out,'c')
    fig,ax=panel(PAGE_W/2,195,'GSE73950',(.20,.24,.74,.60),True)
    for name,color in [('Endometriosis',HIGH),('Control',LOW)]:
        sub=d[d.arm==name];ax.scatter(sub.r_reference,sub.historical_Horvathv1,s=24,marker='^',c=color,ec='black',lw=.4,alpha=.8,clip_on=False,label=f'{name} (n = {len(sub)})')
    ax.set(xlim=(-.15,1.03),ylim=(-12,67),xticks=[0,.5,1],yticks=[0,20,40,60],xlabel='Reference methylome correlation (r)',ylabel='Horvath age before repair (years)');grid(ax)
    ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(.0,1.00),handletextpad=.3,labelspacing=.3)
    ctrl=d[d.arm=='Control'];rmed=ctrl.r_reference.median();amed=ctrl.historical_Horvathv1.median()
    ax.annotate(f'{len(ctrl)} controls overlap\nr = {rmed:.3f}; {amed:.1f} years',xy=(rmed,amed),xytext=(.22,3),fontsize=6.5,color=LOW,arrowprops={'arrowstyle':'-','color':LOW,'lw':.6})
    save(fig,out,'d')
    fig,ax=panel(PAGE_W/2,195,'DNAm input integrity',(.20,.24,.74,.60),True)
    width=e.bin_right-e.bin_left;ax.bar(e.bin_left,e.n,width=width,align='edge',color=CLASSIC[0],edgecolor='white',lw=.35,log=True)
    ax.bar(e.bin_left[e['below_0.10_n']>0],e.loc[e['below_0.10_n']>0,'below_0.10_n'],width=width[e['below_0.10_n']>0],align='edge',color=HIGH,ec='white',lw=.35,log=True)
    for t in [.1,.5]:ax.axvline(t,color=INK,lw=.7,ls='--')
    ax.set(xlim=(-.5,1.03),xticks=[-.5,0,.5,1],ylabel='Samples',xlabel='Reference methylome correlation (r)')
    ax.set_yticks([1,10,100,1000,10000],['1','10','100','1,000','10,000']);ax.minorticks_off()
    facts=json.loads((OLD/'ed21/numbers.json').read_text())['e_screen']
    ax.text(.04,.96,f'n = {facts["samples"]:,}\n{facts["studies"]:,} studies',transform=ax.transAxes,va='top',fontsize=7)
    ax.text(.04,.63,f'{facts["below_0.10"]:,} below 0.10',transform=ax.transAxes,fontsize=7,color=HIGH)
    for t in [.1,.5]:ax.text(t,.97,f'{t:.2f}',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=6.5)
    ax.grid(axis='y',lw=.4,ls=':',alpha=.4);save(fig,out,'e')
    fig,ax=panel(PAGE_W,165,'Ranking eligibility',(.34,.25,.61,.59),True)
    yy=np.arange(5)[::-1];vals=f.comparisons
    ax.barh(yy,vals,color=CLASSIC[0],ec='black',lw=.5,alpha=.8)
    for y,v in zip(yy,vals):ax.text(v+95,y,f'{v:,}',va='center',fontsize=7)
    ax.set_yticks(yy,['FDR < 0.05','Contrast type agreed by two models','No design concern','Design audit admitted','Ranked'])
    ax.set_xlim(0,7900);ax.set_xticks([0,2000,4000,6000]);ax.set_xlabel('Comparisons');grid(ax);ax.tick_params(axis='y',length=0)
    save(fig,out,'f')
    rects={'a':(0,0,PAGE_W,125),'b':(0,125,355,390),'c':(355,125,PAGE_W,350),
       'd':(0,390,PAGE_W/2,585),'e':(PAGE_W/2,390,PAGE_W,585),'f':(0,585,PAGE_W,750)}
    finish(out,rects,750,13,{'run':Path(os.environ['CLOCKBASE_FINAL']).name,
        'b_defects':b.to_dict('records'),'b_defect_labels':dict(DEFECTS),
        'f_ranking_eligibility':f.to_dict('records')})

DISPLAY={'3;4 dihydrocollidine;5 diethoxycarbonyl 1':'DDC',"theiler's murine encephalomyelitis virus":"Theiler's virus",'respiratory syncytial virus':'RSV','hepatitis virus':'Mouse hepatitis virus','escherichia coli':'E. coli','lipopolysaccharide':'LPS','dextran sodium sulfate':'DSS','anti cd40 antibody':'Anti-CD40','heligmosomoides polygyrus':'H. polygyrus','conduritol b epoxide':'Conduritol B epoxide'}
def names(rows,counts=False,projects=False):
    multi=rows.groupby(['species','canonical_target']).condition_category.transform('nunique')>1
    labels=[]
    for repeated,r in zip(multi,rows.itertuples()):
        target=str(r.canonical_target)
        if target in DISPLAY: name=DISPLAY[target]
        elif r.intervention_class=='GENE':
            name=(target.upper() if r.species=='human' else target[:1].upper()+target[1:])+{'DECREASE':' loss','INCREASE':' gain'}.get(r.manipulation,'')
        elif re.match(r'^[a-z]{2,5}[ -]?\d',target):name=target.upper()
        else:name=target[:1].upper()+target[1:]
        if getattr(r,'large_arm',False):name+='†'
        if repeated:name+=', '+r.condition_category.lower().replace('environment','env.')
        if r.species=='human':name+=' (human)'
        if counts:name+=f' [{r.k}/{r.n_series}]'
        if projects:name+=f' ({r.k})'
        labels.append(name)
    return labels

def forest(ax,rows,color,counts=True):
    yy=np.arange(len(rows))[::-1]
    ax.barh(yy,rows.estimate,height=.62,color=color,ec='#111',lw=.5)
    ax.errorbar(rows.estimate,yy,xerr=[rows.estimate-rows.ci_low,rows.ci_high-rows.estimate],fmt='o',ms=4,mfc='#f7f7f7',mec='#222',mew=.8,ecolor='#222',lw=.75,capsize=2,zorder=3)
    ax.set_yticks(yy,names(rows,counts=counts));ax.set_ylim(-.6,len(rows)-.25);ax.tick_params(axis='y',length=0)
    ax.axvline(0,color=INK,lw=.7,ls='--');grid(ax);ax.set_xlabel("Pooled Hedges' g (95% CI)")
    ax.set_xlim(np.floor(min(rows.ci_low.min(),0))-.1,np.ceil(max(rows.ci_high.max(),0))+.1)

def project_panel(out,letter,title,df,pooled,height):
    fig,ax=panel(PAGE_W/2,height,title,(.47,.23,.48,.61))
    yy=np.arange(len(df))[::-1]+1.3
    ax.errorbar(df.g,yy,xerr=[df.g-df.ci_low_plotted,df.ci_high_plotted-df.g],fmt='s',ms=4,mec=INK,mew=.45,color=LOW,ecolor=LOW,lw=1,capsize=2,capthick=.6,clip_on=False)
    ax.add_patch(Polygon([[pooled['ci_low'],0],[pooled['estimate'],.29],[pooled['ci_high'],0],[pooled['estimate'],-.29]],fc=LOW,ec=INK,lw=.6))
    ax.set_yticks(list(yy)+[0],df.display_label.tolist()+[f'Pooled ({pooled["k"]} projects)'])
    ax.set_ylim(-.65,max(yy)+.7);ax.set_xlim(-4.7,2.5);ax.set_xticks([-4,-2,0,2]);ax.tick_params(axis='y',length=0)
    ax.axvline(0,color=INK,lw=.7,ls='--');grid(ax);ax.set_xlabel("Hedges' g (95% CI)");save(fig,out,letter)

def build22():
    SOURCE_FILES.clear()
    setup('Calibri');out=OUT/'ed22';out.mkdir(parents=True,exist_ok=True)
    data={l:read(OLD/f'ed22/{name}.csv') for l,name in [('a','panel_a_lower'),('b','panel_b_sensitivity'),('c','panel_c_higher'),('d','panel_d_counting_unit'),('e','panel_e_sleep_deprivation'),('f','panel_f_rapamycin'),('g','panel_g_lifespan_matched')]}
    facts=json.loads((OLD/'ed22/numbers.json').read_text())
    a,b,c,d,e,f,g=[data[l] for l in 'abcdefg']
    fig,ax=panel(PAGE_W/2,145,'Lower clock outputs',(.45,.27,.49,.49));forest(ax,a,LOW);ax.set_xticks([-3,-2,-1,0]);save(fig,out,'a')
    fig,ax=panel(PAGE_W/2,145,'Sensitivity',(.22,.32,.73,.43))
    cols=b.specification.drop_duplicates().tolist(); interventions=b.intervention.drop_duplicates().tolist()
    for r in b.itertuples():
        x=cols.index(r.specification);y=len(interventions)-1-interventions.index(r.intervention)
        if r.state=='under two units':ax.plot([x-.16,x+.16],[y,y],color='#999',lw=1)
        else:ax.scatter(x,y,s=28,fc=LOW if r.state=='robust' else 'white',ec=LOW,lw=.7,clip_on=False)
    ax.set_xticks(range(5),['Primary','Series','All\ncontexts','No\ncontact','H–K']);ax.tick_params(axis='x',length=0)
    ax.set_yticks(range(len(interventions)));ax.tick_params(axis='y',length=0)
    # Link sensitivity rows to exactly the source labels, including the category.
    ax.set_yticklabels(list(reversed(names(a))))
    ax.set_ylim(-.6,len(interventions)-.4);ax.set_xlim(-.5,4.5);ax.spines[:].set_visible(False)
    ax.legend(handles=[Line2D([],[],marker='o',ls='',mfc=LOW,mec=LOW,ms=4,label='Robust'),Line2D([],[],marker='o',ls='',mfc='white',mec=LOW,ms=4,label='Not robust')],ncol=2,frameon=False,loc='lower center',bbox_to_anchor=(.5,-.64),handlelength=.9,columnspacing=.7,handletextpad=.3)
    save(fig,out,'b')
    fig,ax=panel(PAGE_W/2,225,'Higher clock outputs',(.51,.17,.43,.73));forest(ax,c,HIGH);ax.set_xticks([0,5,10]);save(fig,out,'c')
    fig,ax=panel(PAGE_W/2,140,'Counting unit',(.27,.25,.67,.55))
    yy=np.arange(4)[::-1]
    for off,col,color,title in [(.18,'series',CLASSIC[0],'Series'),(-.18,'projects',CLASSIC[1],'Projects')]:
        vals=d[col];ax.barh(yy+off,vals-1,left=1,height=.34,color=color,ec='black',lw=.5)
        for y,v in zip(yy+off,vals):ax.text(v*1.15,y,str(v),fontsize=6.5,va='center')
    ax.set_xscale('log');ax.set_xlim(1,1800);ax.set_xticks([1,10,100,1000],['1','10','100','1,000']);ax.minorticks_off()
    ax.set_yticks(yy,d.group);ax.set_ylim(-.6,3.6);ax.tick_params(axis='y',length=0);ax.set_xlabel('Interventions');grid(ax)
    ax.legend(handles=[Patch(fc=CLASSIC[0],ec='black',lw=.5,label='Series'),Patch(fc=CLASSIC[1],ec='black',lw=.5,label='Projects')],ncol=2,frameon=False,loc='upper right',bbox_to_anchor=(1,1.40),handlelength=.9)
    save(fig,out,'d')
    project_panel(out,'e','Sleep deprivation',e,facts['sleep_pooled'],175)
    project_panel(out,'f','Rapamycin',f,facts['rapamycin_pooled'],190)
    fig,ax=panel(PAGE_W/2,375,'Lifespan matches',(.48,.20,.46,.72))
    yy=np.arange(len(g))[::-1]
    for y,r in zip(yy,g.itertuples()):
        co=HIGH if r.expected>0 else LOW
        ax.errorbar(r.estimate,y,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt='o',ms=4.7,mfc=co if r.concordant else 'white',mec=co,mew=.8,ecolor=co,lw=1,capsize=2,capthick=.7,clip_on=False)
    ax.set_yticks(yy,names(g,projects=True));ax.set_ylim(-.7,len(g)-.3);ax.tick_params(axis='y',length=0)
    for tick,r in zip(ax.get_yticklabels(),g.itertuples()):
        if r.source=='GenAge':tick.set_fontstyle('italic')
    ax.set_xlim(-6.1,6.1);ax.set_xticks([-6,-3,0,3,6]);ax.axvline(0,color=INK,ls='--',lw=.75);grid(ax);ax.set_xlabel("Pooled Hedges' g (95% CI)")
    ax.legend(handles=[Line2D([],[],marker='o',ls='',mfc=HIGH,mec=HIGH,ms=4,label='Higher'),Line2D([],[],marker='o',ls='',mfc=LOW,mec=LOW,ms=4,label='Lower'),Line2D([],[],marker='o',ls='',mfc='white',mec=INK,ms=4,label='Opposite')],title='Expected direction',title_fontsize=7,ncol=3,frameon=False,loc='upper center',bbox_to_anchor=(.35,-.16),handlelength=.8,columnspacing=.6,handletextpad=.3)
    save(fig,out,'g')
    h=read(OLD/'figure4/panel_i_replication_points.csv');named=read(OLD/'figure4/panel_i_named_points.csv')
    replication = json.loads((OLD/'figure4/numbers.json').read_text())['replication']
    assert len(h) == replication['two_or_more_projects']
    assert int((h.robust & h.estimate.gt(0)).sum()) == replication['robust_higher']
    assert int((h.robust & h.estimate.lt(0)).sum()) == replication['robust_lower']
    fig,ax=panel(PAGE_W/2,260,'Replication across independent projects',(.18,.22,.76,.62))
    nlf=-np.log10(h.fdr.clip(lower=1e-300));sig=h.fdr<.05
    ax.scatter(h.loc[~sig,'estimate'],nlf[~sig],s=7,c='lightgray',alpha=.3,clip_on=False)
    vmax=h.loc[sig,'estimate'].abs().max();norm=Normalize(-vmax,vmax)
    ax.scatter(h.loc[sig,'estimate'],nlf[sig],s=22,c='black',clip_on=False,zorder=5)
    mapped=ax.scatter(h.loc[sig,'estimate'],nlf[sig],s=15,c=h.loc[sig,'estimate'],cmap='RdBu_r',norm=norm,alpha=.9,clip_on=False,zorder=6)
    ax.axhline(-np.log10(.05),ls='--',color='#666',lw=.7,zorder=0);ax.axvline(0,color='#aaa',lw=.5,zorder=0)
    ax.set_xlabel("Pooled Hedges' g");ax.set_ylabel('−log10(BH-FDR)');ax.set_xlim(h.estimate.min()-.5,h.estimate.max()+.5);ax.set_ylim(-.2,nlf.max()+1)
    cb=fig.colorbar(mapped,ax=ax,pad=.03,fraction=.04,shrink=.6);cb.set_label("Hedges' g",fontsize=7);cb.ax.tick_params(labelsize=7)
    cb.solids.set_rasterized(False)
    # Fixed dotted bracket labels copied from July combined volcano, with source
    # category colours. The selected points are fixed and do not overlap here.
    locations={'lipopolysaccharide':(1.1,16.3),'bleomycin':(4.3,5.1),'sleep deprivation':(-3.7,6.0),'plx5622':(-3.7,4.4)}
    for r in named.itertuples():
        target=r.canonical_target.lower();text={'lipopolysaccharide':'LPS','plx5622':'PLX5622'}.get(target,target)
        if target not in locations:continue
        x,y=locations[target];point=(r.estimate,-np.log10(max(r.fdr,1e-300)))
        arrow={'arrowstyle':'-','color':'black','ls':':','alpha':.4,'lw':1}
        if target=='sleep deprivation':arrow['relpos']=(1,0)
        ax.annotate(text,xy=point,xytext=(x,y),fontsize=6.5,fontweight='bold',color=CATEGORY[r.condition_category],arrowprops=arrow,bbox={'boxstyle':'round,pad=0.3','fc':'white','ec':'none','alpha':1},zorder=10)
    ax.text(.98,.98,f"{len(h)} interventions\n{replication['robust_higher']} robust higher\n{replication['robust_lower']} robust lower",transform=ax.transAxes,ha='right',va='top',fontsize=7)
    save(fig,out,'h')
    rects={'a':(0,0,PAGE_W/2,145),'b':(PAGE_W/2,0,PAGE_W,145),
       'c':(0,145,PAGE_W/2,370),'d':(PAGE_W/2,145,PAGE_W,285),
       'e':(0,370,PAGE_W/2,545),'f':(PAGE_W/2,285,PAGE_W,475),
       'g':(PAGE_W/2,475,PAGE_W,850),'h':(0,545,PAGE_W/2,805)}
    finish(out,rects,850,25,{'h_volcano':{'interventions':len(h),'robust_higher':replication['robust_higher'],'robust_lower':replication['robust_lower'],'threshold':.05,'named':named.to_dict('records')}})

def finish(out,rects,height,july_page,additions):
    for p in (OLD/out.name).glob('panel_*.csv'):
        if out.name=='ed21' and p.name=='panel_b_defects.csv':continue
        copy_file(p,out/p.name)
    copy_file(OLD/out.name/'provenance.json',out/'source_provenance_2026_11_05.json')
    oldfacts=json.loads((OLD/out.name/'numbers.json').read_text());oldfacts.update(additions)
    if out.name=='ed21':pd.DataFrame(additions['b_defects']).to_csv(out/'panel_b_defects.csv',index=False)
    if out.name=='ed21':copy_file(OLD/'figure3/panel_f_values.csv',out/'panel_f_ranking_eligibility.csv')
    else:
        copy_file(OLD/'figure4/panel_i_replication_points.csv',out/'panel_h_replication_points.csv')
        copy_file(OLD/'figure4/panel_i_named_points.csv',out/'panel_h_named_points.csv')
    (out/'numbers.json').write_text(json.dumps(oldfacts,indent=2))
    doc=fitz.open();page=doc.new_page(width=PAGE_W,height=height)
    for letter,rect in rects.items():
        with fitz.open(out/f'panel_{letter}.pdf') as src:page.show_pdf_page(fitz.Rect(rect),src,0,keep_proportion=True)
    for letter,rect in rects.items():page.insert_text((rect[0]+5,rect[1]+16),letter,fontsize=13,fontname='hebo',color=(0,0,0))
    doc.save(out/f'{out.name}.pdf',garbage=4,deflate=True);doc.close()
    with fitz.open(out/f'{out.name}.pdf') as doc:
        pg=doc[0];pg.get_pixmap(dpi=200,alpha=False).save(out/f'{out.name}_review.png')
        spans=[s for b in pg.get_text('dict')['blocks'] if 'lines' in b for line in b['lines'] for s in line['spans'] if s['text'].strip()]
        pdfchecks={'width_pt':pg.rect.width,'height_pt':pg.rect.height,'raster_images':len(pg.get_images()),'vector_paths':len(pg.get_drawings()),'minimum_pdf_text_pt':min(s['size'] for s in spans),'native_panel_rects_pt':rects}
        if pdfchecks['minimum_pdf_text_pt']<6.49:raise ValueError('PDF text fell below native annotation size')
        if pdfchecks['raster_images']:raise ValueError('Raster plot embedded in vector PDF')
    with fitz.open(JULY) as july:
        pg=july[july_page-1]; pix=pg.get_pixmap(dpi=200,alpha=False); reference=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
        # Crop the July figure artwork below its page title, preserve physical scale.
        scale=200/72; crop={13:(65,95,550,420),25:(70,90,560,515)}[july_page]
        reference=reference.crop(tuple(round(v*scale) for v in crop))
    new=Image.open(out/f'{out.name}_review.png').convert('RGB'); gap=40
    comp=Image.new('RGB',(new.width+reference.width+gap,max(new.height,reference.height)+50),'white')
    draw=ImageDraw.Draw(comp);draw.text((10,10),f'{out.name.upper()} rebuilt 2026-11-06',fill='black');draw.text((new.width+gap,10),f'July Extended Data Figure {july_page-6}',fill='black')
    comp.paste(new,(0,50));comp.paste(reference,(new.width+gap,50));comp.save(out/f'{out.name}_vs_july.png')
    unchanged=[]
    for p in (OLD/out.name).glob('panel_*.csv'):
        if out.name=='ed21' and p.name=='panel_b_defects.csv':continue
        assert p.read_bytes()==(out/p.name).read_bytes();unchanged.append(p.name)
    provenance=[]
    for p in dict.fromkeys(SOURCE_FILES+[OLD/out.name/'numbers.json']):
        provenance.append({'path':str(p),'mtime_ns':p.stat().st_mtime_ns,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2))
    (out/'checks.json').write_text(json.dumps({'pdf':pdfchecks,'unchanged_source_csvs':unchanged,'panel_checks':{k:v for k,v in CHECKS.items() if k.startswith(out.name)}},indent=2,default=str))
    write_readme(out,rects,additions,pdfchecks)

def write_readme(out,rects,additions,checks):
    ed=out.name;is21=ed=='ed21'
    panelmap=([('a','Final classification','panel_a_classification.csv'),('b','Defects found by the verification','panel_b_defects.csv'),('c','Published results','panel_c_published_fates.csv'),('d','GSE73950','panel_d_sample_values.csv'),('e','DNAm input integrity','panel_e_histogram.csv'),('f','Ranking eligibility','panel_f_ranking_eligibility.csv')] if is21 else
      [('a','Lower clock outputs','panel_a_lower.csv'),('b','Sensitivity','panel_b_sensitivity.csv'),('c','Higher clock outputs','panel_c_higher.csv'),('d','Counting unit','panel_d_counting_unit.csv'),('e','Sleep deprivation','panel_e_sleep_deprivation.csv'),('f','Rapamycin','panel_f_rapamycin.csv'),('g','Lifespan matches','panel_g_lifespan_matched.csv'),('h','Replication across independent projects','panel_h_replication_points.csv')])
    lines=[f'# {"Atlas verification" if is21 else "Intervention replication"}','',
      ('Panel b reads every listed defect category from the current final error census; percentages are 100 × n / of. The other supplied panel CSVs are preserved byte for byte. ' if is21 else 'The supplied panel CSVs are preserved byte for byte. ') + 'Each panel is drawn in its own vector PDF and placed in a fixed PyMuPDF rectangle at native size. The page width is 595.68 pt, matching the July vector assembly.','',
      '## Panel map','', '| Panel | Content | Data source |','| --- | --- | --- |']
    lines += [f'| {l} | {title} | [{name}]({name}) |' for l,title,name in panelmap]
    lines += ['','## Style sources','',
      '- July ED7 uses `Aging-AI-Agent/src/visualize_validation_results.py`: Arial, #4C72B0 and #DD8452 bars at alpha 0.8, black 0.5-pt outlines and numerical bar labels. The remaining classification colours follow its seaborn muted palette.',
      '- July ED19 uses `scripts/legacy_sensitivity_analyses_2026_05_14/external_db_crossvalidation_2026-05-14.py`: Calibri 9 pt, 0.6-pt axes and ticks, top/right spines removed, compact horizontal bars and 7–9-pt labels. ED22 copies this typography; ED21 follows ED7 Arial.',
      '- July Figure 3b–e uses `Aging-AI-Agent/create_combined_modality_top_hits_v2.py`: #D55E00 higher and #0072B2 lower, #111 0.5-pt bar outlines, pale filled end markers with #222 edges, dotted x grids, bold categorical titles. Native source panels are 7 × 5 inches; panel ratios follow that convention where content permits.',
      '- July ED14 uses `scripts/legacy_sensitivity_analyses_2026_05_14/extreme_tail_dedicated_analysis.py`: 8.5 × 6.2-inch volcano, dashed FDR threshold and light zero guide. July ED16 uses `paper/legacy_reference/figure_rev1/make_ed16_pooled_volcano.py`: multiple compact category volcanoes. July ED17 uses `scripts/legacy_sensitivity_analyses_2026_05_14/effect_size_distribution_2026-05-14.py`: 12 × 4.5-inch three-panel histogram, white-edged bins, modality colours and dashed reference lines. ED21e follows the ED17 histogram treatment.',
      '- ED22h copies `Aging-AI-Agent/create_combined_modality_volcano.py`: lightgray nonsignificant marks, RdBu_r effect colours, black outlines, bold category-coloured annotation text on white with dotted leaders. The colourbar preserves the continuous effect scale. All marks stay vector; the July generator rasterized its dense layers.',
      '- `paper/legacy_reference/figure_rev1/assemble_figure4_vector.py` supplies 595.68-pt page width, fixed rectangles, `show_pdf_page` and bold Helvetica 13-pt letters. Panel PDFs are placed without downscaling. Final tick labels are 7 pt or greater; annotations are 6.5 pt or greater. Layout differences from July accommodate the additional panels and print-size text requirement.','',
      '## Numerical provenance','',
      'Every original panel value and annotation is retained in [numbers.json](numbers.json). Its original empirical source paths and hashes are retained in [source_provenance_2026_11_05.json](source_provenance_2026_11_05.json). [provenance.json](provenance.json) records the immediate CSV source paths, timestamps and SHA-256 hashes. Axis ticks are coordinates in the labelled units. No value is simulated.','']
    if is21:
        lines += ['Panel a prints segment counts where the segment is at least 9% wide and the modality totals. Panel c prints the published-item totals. Both are computed only from the corresponding CSV n and total columns. The full count tables follow.','']
        for p in ['panel_a_classification.csv','panel_c_published_fates.csv']:
            df=pd.read_csv(out/p);lines += [df.pivot(index='group',columns='class',values='n').to_markdown(), '']
        facts=json.loads((out/'numbers.json').read_text());dd=facts['d_index_study'];ee=facts['e_screen']
        lines += [f'Panel d prints 26 endometriosis and 17 control samples. Its overlapping-control annotation rounds the source median correlation {dd["control_median_r"]:.12g} to −0.002 and historical age {dd["control_median_age"]:.12g} to −0.7 years. Panel e prints {ee["samples"]:,} samples, {ee["studies"]:,} studies and {ee["below_0.10"]:,} correlations below 0.10. Reference thresholds 0.10 and 0.50 and all histogram bin counts remain unchanged.','',
          'Panel f adds the ranking funnel from `figures-rebuild-2026-11-05/figure3/panel_f_values.csv`. The five printed counts are:','',pd.read_csv(out/'panel_f_ranking_eligibility.csv').to_markdown(index=False),'']
    else:
        lines += ['Panels a/c print project/series counts in brackets; panel g prints project counts in parentheses. They come directly from k and n_series in each corresponding CSV. Panels e/f print the supplied project display labels and pooled project counts. Their pooled estimates and confidence intervals are preserved in numbers.json. Panel b preserves all source sensitivity decisions, including the < 2 units key and H–K interval sensitivity. A dagger retains the original arm-size threshold of 120. Panel d prints:','',pd.read_csv(out/'panel_d_counting_unit.csv').to_markdown(index=False),'',
          f"Panel h shows {len(pd.read_csv(out/'panel_h_replication_points.csv'))} pooled independent-project rows. The robustness criterion remains k ≥ 2, BH-FDR < 0.05 and rank_bound > 0. Counts and plotted values come from the current Figure 4 source tables; LPS, sleep deprivation, PLX5622 and bleomycin coordinates come from panel_i_named_points.csv.",'',
          'The project and series counts printed beside intervention labels are recorded below.','']
        for p in ['panel_a_lower.csv','panel_c_higher.csv','panel_g_lifespan_matched.csv']:
            df=pd.read_csv(out/p);lines += [df[['canonical_target','condition_category','k','n_series']].to_markdown(index=False),'']
    lines += ['## Differences from the previous build','',
      ('Panel f uses the current ranking eligibility counts listed above. Panel b reports the current census, including implausible clock outputs and single-cell units.' if is21 else 'Panel h uses the current pooled intervention effects and robustness counts. These repeat the current Figure 4 data.'),
      'The plotting implementation now draws independent panels, copies July script conventions and uses vector PDF composition instead of a freely styled single matplotlib canvas. The original scripts remain untouched. [Source diff](build_vs_2026_11_05.diff) records the implementation changes.','',
      '## Validation','',
      ('Panel b was recomputed from the current census; the other original panel CSVs are byte-identical to the supplied sources. ' if is21 else 'All original panel CSV files are byte-identical to the supplied sources. ') + f'The assembled PDF has {checks["vector_paths"]} vector paths, no raster plot images and minimum PDF text size {checks["minimum_pdf_text_pt"]:.2f} pt. Native panel rectangles and automated findings are in [checks.json](checks.json). Generic house style and aspect findings retain the existing July geometry overrides. The review PNG is rendered at 200 dpi. The side-by-side image compares artwork at the same physical scale; the July page title is outside the reference crop.','',
      'The final 200-dpi render and same-scale July comparison were inspected. No clipped legend, overlapping text or truncated confidence interval remains.','']
    (out/'README.md').write_text('\n'.join(lines))
    import difflib
    before=(REPO/f'scripts/figures_2026_11_05/restore_{ed}_2026_11_05.py').read_text().splitlines(True)
    after=Path(__file__).read_text().splitlines(True)
    (out/'build_vs_2026_11_05.diff').write_text(''.join(difflib.unified_diff(before,after,fromfile=f'restore_{ed}_2026_11_05.py',tofile=Path(__file__).name)))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--figure',choices=['ed21','ed22','both'],default='both')
    args=parser.parse_args()
    try:
        if args.figure in ('ed21','both'):build21()
        if args.figure in ('ed22','both'):build22()
        print('Built',args.figure,'in',OUT)
    except Exception as exc:
        print(f'ED21/22 build failed: {exc}',file=sys.stderr)
        raise
