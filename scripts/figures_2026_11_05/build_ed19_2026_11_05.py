#!/usr/bin/env python3
"""Restore July ED19 bar panels from frozen annotations and exact ranking joins."""
from pathlib import Path
import os,json,hashlib
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fitz
from PIL import Image,ImageOps,ImageDraw
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-05/ed19'
SRC=ROOT/'paper/revision_2026-09-07/evidence/ed19-source-recovery-2026-09-18/external_db_crossvalidation_2026-05-14.csv'
RANK=Path(os.environ['CLOCKBASE_FINAL'])/'comparison-ranking/ranked_comparisons_all.csv'
JULY=ROOT/'paper/revision_2026-09-07/baseline/03_Figures_2026-07-06.pdf'
STATUS=['direction_consistent','mechanism_only','mixed','direction_mismatched','no_db_evidence']
COLORS=['#3E7B5A','#9C9C9C','#D9A441','#A53636','#3E5C7B']
def build():
    x=pd.read_csv(SRC).fillna('');r=pd.read_csv(RANK)
    assert not x.duplicated(['intervention','modality']).any()
    keys=['modality','gse_id','label']
    assert not r.duplicated(keys).any(), 'Ambiguous ranking label key'
    m=x.merge(r,left_on=['modality','gse_id','intervention'],right_on=keys,validate='one_to_one')
    assert m.v3_final_class.isin(['VERIFIED','PROVISIONAL']).all()
    m.to_csv(OUT/'ranking_join.csv',index=False)
    unmatched=x.merge(r[keys],left_on=['modality','gse_id','intervention'],right_on=keys,how='left',indicator=True)
    unmatched[unmatched._merge.eq('left_only')].drop(columns=['_merge','label']).to_csv(OUT/'unranked_annotations.csv',index=False)
    c=m[m.match_status.eq('direction_consistent')].sort_values('rank_bound',ascending=False).head(10)
    d=m[m.match_status.eq('no_db_evidence')].sort_values('rank_bound',ascending=False).head(10)
    dbs=['LongevityMap','AgingAtlas','OpenTargets','DrugAge','GenAge']
    hit=[x.lmap_hit.ne('').sum(),(x.aa_gene_hit.ne('')|x.aa_drug_hit.ne('')).sum(),x.ot_hit_id.ne('').sum(),x.drugage_hit.ne('').sum(),x.genage_hit.ne('').sum()]
    counts=x.match_status.value_counts().reindex(STATUS,fill_value=0)
    # July generator make_figure supplies the bar types, palette, order, strokes and font.
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,'xtick.major.width':.6,'ytick.major.width':.6,'pdf.fonttype':42,'ps.fonttype':42})
    fig=plt.figure(figsize=(183/25.4,7.1))
    a=fig.add_axes([.12,.64,.31,.29]);b=fig.add_axes([.70,.64,.28,.29]);cc=fig.add_axes([.23,.10,.25,.34]);dd=fig.add_axes([.74,.10,.24,.34])
    rates=[v/len(x)*100 for v in hit]
    bars=a.bar(range(5),rates,color='#B85A3E',edgecolor='black',linewidth=.6)
    for bar,n in zip(bars,hit):a.text(bar.get_x()+bar.get_width()/2,bar.get_height()+.55,f'n={n}',ha='center',fontsize=7)
    a.set_xticks(range(5),dbs,rotation=30,ha='right',fontsize=7);a.set_ylabel('Hit rate (%)');a.set_ylim(0,max(rates)*1.25+5)
    bars=b.barh(range(5),counts,color=COLORS,edgecolor='black',linewidth=.6)
    for bar,n in zip(bars,counts):b.text(n+2,bar.get_y()+bar.get_height()/2,str(n),va='center',fontsize=7)
    b.set_yticks(range(5),['direction consistent','mechanism only','mixed','direction mismatched','no recorded match'],fontsize=7);b.invert_yaxis();b.set_xlabel('Number of interventions');b.set_xlim(0,220)
    printed=[]
    for ax,z,title,col,panel in [(cc,c,'Direction-consistent database match','#3E7B5A','c'),(dd,d,'No recorded database match','#3E5C7B','d')]:
        z=z.iloc[::-1]
        ax.barh(range(len(z)),z.rank_bound,color=col,edgecolor='black',linewidth=.5)
        # Shortened label spellings do not alter the documented intervention identity.
        labels=[s.replace('Intestinal Epithelial RIPK3 Knockout','Intestinal epithelial RIPK3 KO').replace('Liver Insulin Receptor Knockout','Liver insulin receptor KO').replace('Atg7 motor-neuron cKO in SOD1G93A background','Atg7 motor-neuron cKO').replace('Foxo1 KO in uterine epithelium','Uterine Foxo1 KO').replace('Nicotinamide Riboside','Nicotinamide riboside').replace('Inactive (remission) SJIA monocytes','Inactive SJIA monocytes') for s in z.intervention]
        ax.set_yticks(range(len(z)),labels,fontsize=7)
        ax.set_xlabel("CI bound nearest zero (|g|)",fontsize=7)
        ax.set_title(title,fontsize=7.3,pad=10)
        ax.set_xlim(0,3.05 if panel=='c' else 5.5);ax.set_xticks([0,1,2,3] if panel=='c' else [0,2,4])
        for _,row in z.iloc[::-1].iterrows(): printed.append({'panel':panel,'intervention':row.intervention,'gse_id':row.gse_id,'comparison_name':row.comparison_name,'clock':row.clock,'hedges_g':row.hedges_g,'g_ci_low':row.g_ci_low,'g_ci_high':row.g_ci_high,'rank_bound':row.rank_bound,'n_case':row.v3_n_case,'n_control':row.v3_n_control,'v3_final_class':row.v3_final_class})
    for ax,label in zip([a,b,cc,dd],'abcd'):fig.text(ax.get_position().x0-(.04 if label in 'ab' else .18),ax.get_position().y1+.025,label,fontsize=12,fontweight='bold')
    fig.savefig(OUT/'ed19.pdf');fig.savefig(OUT/'ed19_review.png',dpi=200)
    import sys
    sys.path.insert(0,str(Path('pipeline_assets/figure_style_2026_09_10')))
    from nature_figure_utils import strict_figure_qc
    qc=strict_figure_qc(fig,dict(zip('abcd',[a,b,cc,dd])),check_text_figure_edges=True,forbid_uniform_grid=False,check_house_style=False,check_data_occupancy=False,check_matched_limits=False,check_text_spine=False,max_tall_panel_aspect=10,max_sparse_panel_aspect=10,check_max_fontsize=False)
    (OUT/'defect_checks.json').write_text(json.dumps(qc,indent=2,default=str))
    plt.close(fig)
    pdf=fitz.open(OUT/'ed19.pdf');fonts=[s['size'] for block in pdf[0].get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for s in line['spans']]
    manifest={'annotation_rows':len(x),'database_counts':dict(zip(dbs,map(int,hit))),'database_rates_percent':dict(zip(dbs,rates)),'status_counts':{k:int(v) for k,v in counts.items()},'joined':len(m),'unjoined':len(x)-len(m),'joined_by_status':m.match_status.value_counts().to_dict(),'bar_values':printed,'pdf_size_points':list(pdf[0].rect),'min_pdf_text_pt':min(fonts),'sources':{str(p):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'mtime':p.stat().st_mtime} for p in [SRC,RANK,JULY]}}
    (OUT/'numbers_and_sources.json').write_text(json.dumps(manifest,indent=2))
    make_comparison()
    print(json.dumps(manifest,indent=2))
def make_comparison():
    from PIL import ImageChops
    july=fitz.open(JULY)
    # Page 25's placed image contains the complete July artwork, excluding
    # the package header. Crop the surrounding white pixels in both drawings.
    image_rect=fitz.Rect(81.05,101.65,549.05,505.80)
    pix=july[24].get_pixmap(matrix=fitz.Matrix(3,3),clip=image_rect)
    old=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
    new=Image.open(OUT/'ed19_review.png').convert('RGB')
    def artwork(im):
        grey=ImageChops.difference(im,Image.new('RGB',im.size,'white')).convert('L')
        box=grey.point(lambda v:255 if v>25 else 0).getbbox()
        return im.crop(box)
    old,new=artwork(old),artwork(new)
    width=max(old.width,new.width)
    old=old.resize((width,round(old.height*width/old.width)),Image.Resampling.LANCZOS)
    new=new.resize((width,round(new.height*width/new.width)),Image.Resampling.LANCZOS)
    pad,gap,head=35,55,55
    canvas=Image.new('RGB',(2*width+2*pad+gap,max(old.height,new.height)+head+pad),'white')
    canvas.paste(old,(pad,head));canvas.paste(new,(pad+width+gap,head))
    draw=ImageDraw.Draw(canvas)
    draw.text((pad,15),'July reference',fill='black')
    draw.text((pad+width+gap,15),'Restored figure',fill='black')
    canvas.save(OUT/'ed19_vs_july.png')

if __name__=='__main__':
    try:build()
    except (OSError,ValueError,AssertionError) as exc:raise RuntimeError(f'ED19 build failed: {exc}') from exc
