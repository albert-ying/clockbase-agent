#!/usr/bin/env python3
"""Restore Figure 4 from July generator methods and verified measurements.

Source tokens: July page 4 and copied create_* / analysis_* generators.
Category palette: Disease #A23B72, Drug #2E86AB, Environment #C73E1D,
Genetic #F18F01, Other #6B7A8F. Waffle concordance #3B8EA5/#E64B35.
Arial type; white ground; original panel order and plot types preserved.
Only Nature defect checks apply; no retro_style or stock mosaic is used.
"""
from pathlib import Path
import os,sys,json,re,hashlib,importlib.util,subprocess
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle,Patch
from matplotlib.colors import Normalize
from scipy import stats
from scipy.stats import gaussian_kde,binomtest,fisher_exact

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-05/figure4'
TABLE=OUT.parent/'data/combined_3modality_plot_data_verified_2026-11-05.csv'
META=Path(os.environ.get('CLOCKBASE_META') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'intervention-meta-2026-10-17')
SCRIPT=Path(__file__).parent
COL={'Disease':'#A23B72','Drug':'#2E86AB','Environment':'#C73E1D','Genetic':'#F18F01','Other':'#6B7A8F'}
AGREE,DISAGREE='#3B8EA5','#E64B35'
plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.titlesize':8,'axes.linewidth':.55,'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white','figure.facecolor':'white'})
NUM={};INPUTS=[TABLE,META/'primary/intervention_meta.csv',META/'concordance/lifespan_concordance.csv']

def bh(p):
 p=np.array(p,float);idx=np.argsort(p);out=np.empty(len(p));out[idx]=np.minimum.accumulate((p[idx]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];return np.minimum(out,1)

def save_csv(d,name):
 d.to_csv(OUT/name,index=False)

def base(ax):
 ax.spines[['top','right']].set_visible(False);ax.grid(False);ax.tick_params(width=.55,length=2,pad=2)

def letter(ax,s,x=-.14,y=1.07):ax.text(x,y,s,transform=ax.transAxes,fontweight='bold',fontsize=11,va='top',ha='left')

def title(ax,s):ax.set_title(s,loc='left',pad=5,fontweight='bold')

def fmtp(p):return f'{p:.1e}' if p<.001 else f'{p:.3f}'

def load_module(name,path):
 spec=importlib.util.spec_from_file_location(name,path);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod

def compute_genetic(df):
 # Literal keyword order and matching expressions from the July generator.
 g=df[(df.condition_category=='Genetic')&(df.fdr<.05)].copy();g['intervention_type']='Other'
 keywords={'Knockout':['knockout','knock-out','ko','null','deletion','deleted','-/-'],'Knockdown':['knockdown','knock-down','kd','shrna','sirna','rnai'],'Overexpression':['overexpression','overexpress','transgenic','tg','oe'],'Mutation':['mutant','mutation','mut','variant']}
 for name,words in keywords.items():g.loc[g.condition_name.astype(str).str.contains('|'.join(words),case=False,na=False),'intervention_type']=name
 rows=[]
 for typ in keywords:
  p=g[g.intervention_type==typ];rows.append({'Type':typ,'n':len(p),'lower':int((p.effect_size<0).sum()),'higher':int((p.effect_size>0).sum()),'median_d':p.effect_size.median()})
 r=pd.DataFrame(rows);r=r[r.n>0].copy();baseline=r.lower.sum()/r.n.sum();r['lower_fraction']=r.lower/r.n;r['fold_enrichment']=r.lower_fraction/baseline;r['p_binomial']=[binomtest(int(x.lower),int(x.n),baseline).pvalue for x in r.itertuples()];r=r.sort_values('fold_enrichment')
 chi=stats.chi2_contingency(r[['lower','higher']].values);NUM['genetic']={'all_significant':len(g),'classified':int(r.n.sum()),'baseline_lower_fraction':baseline,'chi2':float(chi.statistic),'chi2_p':float(chi.pvalue),'dof':int(chi.dof),'types':r.to_dict('records')}
 save_csv(g,'panel_cd_genetic_rows.csv');save_csv(r,'panel_c_genetic_summary.csv');return g,r,baseline

def compute_ora(df):
 # Exact July gene-name stripping/validation and custom observed-gene background.
 ref=set()
 for fn in ['gene_reference_mouse.csv','gene_reference_human.csv']:
  p=ROOT/'result'/fn;INPUTS.append(p);ref.update(pd.read_csv(p).symbol.astype(str).str.upper())
 suffixes=[' KNOCKOUT',' KNOCKDOWN',' OVEREXPRESSION',' KO',' KD',' OE',' DELETION',' MUTANT',' -/-',' +/-',' TRANSGENIC',' DEFICIENT',' DEFICIENCY',' NULL',' HETEROZYGOUS']
 def extract(s):
  cn=str(s).upper()
  for a in suffixes:cn=cn.replace(a,'')
  for w in cn.split():
   w=re.sub('[^A-Z0-9]','',w)
   if len(w)>=2 and w[0].isalpha() and w in ref:return w
  return None
 g=df[df.condition_category=='Genetic'].copy();g['gene']=g.condition_name.map(extract);g=g.dropna(subset=['gene']);bg=set(g.gene);query=set(g[g.fdr<.05].gene)
 save_csv(g[['gse_id','comparison_name','condition_name','effect_size','fdr','gene']],'panel_e_gene_assignments.csv')
 rows=[];library_counts={}
 for lib in ['KEGG_2019_Mouse','KEGG_2019_Human']:
  p=OUT/f'{lib}.gmt';INPUTS.append(p);sets={}
  for line in p.read_text().splitlines():
   a=line.strip().split('\t');sets[a[0]]=set(x.upper() for x in a[2:])
  rr=[]
  for term,genes in sets.items():
   members=genes&bg;overlap=query&members
   # July Enrichr custom background method; exact one-sided Fisher test.
   if not overlap:continue
   a=len(overlap);b=len(query)-a;c=len(members)-a;d=len(bg)-a-b-c
   odds,pv=fisher_exact([[a,b],[c,d]],alternative='greater')
   rr.append({'Pathway':term,'Database':lib,'N_Overlap':a,'Pathway_Size':len(members),'Library_Pathway_Size':len(genes),'P_value':pv,'Odds_Ratio':odds,'Genes':';'.join(sorted(overlap))})
  q=bh([r['P_value'] for r in rr])
  for r,v in zip(rr,q):r['FDR']=v
  rows.extend(rr);library_counts[lib]={'terms':len(sets),'tested_nonzero_overlap':len(rr),'FDR_significant':int((q<.05).sum())}
 allr=pd.DataFrame(rows);save_csv(allr,'panel_e_all_ora_results.csv')
 filt=allr[(allr.N_Overlap>=3)&(allr.Pathway_Size<200)].sort_values('FDR').drop_duplicates('Pathway')
 filt=filt[~filt.Pathway.str.contains('cancer|carcinoma|tumor|tumour|neoplasm|oncogene|metastasis|malignancy|leukemia|lymphoma|melanoma|glioma|sarcoma|myeloma',case=False,regex=True)]
 plot=filt.head(10).copy();plot['neg_log10_FDR']=-np.log10(plot.FDR.clip(lower=1e-300));save_csv(plot,'panel_e_plotted_pathways.csv')
 NUM['pathways']={'background_unique_genes':len(bg),'significant_unique_genes':len(query),'libraries':library_counts,'eligible_non_cancer_pathways':len(filt),'FDR_significant_eligible_non_cancer':int((filt.FDR<.05).sum()),'plotted':len(plot),'FDR_significant_plotted':int((plot.FDR<.05).sum()),'background_method':'observed verified Genetic genes; pathway members intersected with background; one-sided Fisher; BH within library over nonzero-overlap terms'}
 return plot

def compute_fda(df):
 # Import copied July matcher; adding cached pure Python statsmodels only.
 sys.path.append('/Users/A.Y/.cache/uv/archive-v0/tU7tOHgnLOjnWV6m')
 module=load_module('july_fda',SCRIPT/'analysis_4_fda_approved_drugs_figure4_july_source_2026_11_05.py')
 # Import resets seaborn global defaults; restore the measured physical type.
 plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.titlesize':8,'axes.linewidth':.55,'pdf.fonttype':42,'ps.fonttype':42})
 annotation=OUT/'opentargets_annotation.json';INPUTS.extend([ROOT/'data/drug_molecule/part-00000-35323264-2241-4256-b643-3a92cd0230e5-c000.snappy.parquet',annotation]);ot=pd.DataFrame(json.loads(annotation.read_text()));lookup=module.build_opentargets_lookup(ot)
 g=df[(df.condition_category=='Drug')&(df.fdr<.05)].copy();cachefile=OUT/'panel_f_match_cache.json'
 cache=json.loads(cachefile.read_text()) if cachefile.exists() else {}
 for idx,name in enumerate(g.condition_name.dropna().unique()):
  if name not in cache:
   matched,status,info=module.match_drug_to_opentargets(name,ot,lookup)
   cache[name]={'matched':bool(matched),'status':status,'info':info}
  if idx%25==0:cachefile.write_text(json.dumps(cache,default=str));print('drug matching',idx,flush=True)
 cachefile.write_text(json.dumps(cache,default=str));g['matched']=g.condition_name.map(lambda x:cache.get(x,{}).get('matched',False));g['clinical_status']=g.condition_name.map(lambda x:cache.get(x,{}).get('status'));g['matched_name']=g.condition_name.map(lambda x:(cache.get(x,{}).get('info') or {}).get('name'));g['match_type']=g.condition_name.map(lambda x:(cache.get(x,{}).get('info') or {}).get('match_type'))
 save_csv(g,'panel_f_drug_matches.csv');matched=g[g.matched].copy();rr=[]
 for status,p in matched.groupby('clinical_status'):
  lo=int((p.effect_size<0).sum());hi=int((p.effect_size>0).sum());rr.append({'status':status,'n':len(p),'lower':lo,'higher':hi,'p_binomial':binomtest(lo,lo+hi,.5).pvalue})
 r=pd.DataFrame(rr).sort_values('n');save_csv(r,'panel_f_phase_counts.csv');chi=stats.chi2_contingency(r[['lower','higher']].values)
 NUM['drugs']={'significant_drug_rows':len(g),'matched_rows':len(matched),'unique_queries':g.condition_name.nunique(),'matched_unique_queries':matched.condition_name.nunique(),'chi2':float(chi.statistic),'chi2_p':float(chi.pvalue),'chi2_df':int(chi.dof),'stages':r.to_dict('records')};return r

def approved_inset():
 # Exact July high-confidence selection and strongest row per standard drug name.
 import ast
 raw=pd.read_csv(OUT/'panel_f_drug_matches.csv');part=raw[(raw.clinical_status=='FDA_Approved')&(raw.effect_size<0)&raw.match_type.isin(['exact','synonym','tradename'])].sort_values('effect_size').drop_duplicates('matched_name').head(15).copy()
 tree=ast.parse((SCRIPT/'analysis_4_fda_approved_drugs_figure4_july_source_2026_11_05.py').read_text());mechanisms={}
 for node in ast.walk(tree):
  if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='drug_mechanisms' for t in node.targets):mechanisms=ast.literal_eval(node.value)
 part['mechanism']=part.matched_name.str.upper().map(mechanisms).fillna('');save_csv(part,'panel_f_approved_lower_drug_inset.csv')
 fig,ax=plt.subplots(figsize=(4.3,3.4));fig.subplots_adjust(left=.35,right=.97,bottom=.16,top=.86);y=np.arange(len(part));ax.barh(y,part.effect_size,color='#27AE60',alpha=.3,edgecolor='white',lw=.5);ax.set_yticks(y,part.matched_name.str.upper(),fontsize=7);ax.set_xlabel("Cohen's d",fontsize=7);ax.axvline(0,color='black',ls='--',lw=.6);fig.text(.35,.95,'Approved drugs with lower clock output',fontsize=8,fontweight='bold');ax.text(0,1.025,f'Top {len(part)} high-confidence matched drugs',transform=ax.transAxes,fontsize=7)
 for yy,r in zip(y,part.itertuples()):
  if r.mechanism:ax.text(-.1,yy,r.mechanism,ha='right',va='center',fontsize=6)
 base(ax);fig.savefig(OUT/'figure4_f_approved_drug_inset.pdf');fig.savefig(OUT/'figure4_f_approved_drug_inset_review.png',dpi=200);plt.close(fig);NUM['approved_inset']={'n':len(part),'selection':'FDR<0.05; Approved; d<0; exact/synonym/tradename; strongest d per matched OpenTargets name; first15','rows':part[['matched_name','effect_size','gse_id','comparison_name','match_type','mechanism']].to_dict('records')}

def waffle(ax,source,letter_s,life,ncol):
 p=life[life.source==source].sort_values(['rank_bound','canonical_target'],ascending=[False,True]).copy();p['display_name']=p.canonical_target.str.upper() if source=='GenAge' else p.canonical_target.str.capitalize()
 n=len(p);k=int(p.concordant.sum());pv=binomtest(k,n,.5).pvalue
 rows=int(np.ceil(n/ncol));size=.9
 for j,r in enumerate(p.itertuples()):
  col=j%ncol;row=rows-1-j//ncol;ax.add_patch(Rectangle((col,row),size,size,facecolor=AGREE if r.concordant else DISAGREE,edgecolor='white',linewidth=.3))
 # July leader-line labels, selected most certain matched genes for legibility.
 labelled=p.head(12 if source=='GenAge' else n)
 for j,r in enumerate(labelled.itertuples()):
  col=j%ncol;row=rows-1-j//ncol;color=AGREE if r.concordant else DISAGREE
  yy=rows+.48 if j%2==0 else -.25;anchor=row+.95 if j%2==0 else row
  ax.plot([col+.45,col+.45],[anchor,yy],color=color,lw=.6);ax.text(col+.45,yy,({'Nicotinamide riboside':'NR'}.get(r.display_name,r.display_name)),rotation=45,ha='left' if j%2==0 else 'right',va='bottom' if j%2==0 else 'top',fontsize=6,color=color,fontstyle='italic' if source=='GenAge' else 'normal')
 ax.set_xlim(-.7,ncol+.6);ax.set_ylim(-1.0,rows+2.4);ax.axis('off');letter(ax,letter_s,x=-.05,y=1.02)
 ax.text(.5,1.2,f'{source} lifespan concordance',ha='center',transform=ax.transAxes,fontsize=8,fontweight='bold');ax.text(.28 if source=='GenAge' else .5,1.0,f'{k}/{n} ({k/n:.0%})   P = {pv:.3f}',ha='center',transform=ax.transAxes,fontsize=7)
 
 if source=='GenAge':ax.legend(handles=[Patch(color=AGREE,label='Expected direction'),Patch(color=DISAGREE,label='Opposite direction')],loc='lower right',bbox_to_anchor=(1,.75),ncol=2,frameon=False,handlelength=.7,columnspacing=1,fontsize=7)
 save_csv(p,f'panel_{letter_s}_{source.lower()}_waffle.csv');NUM[source]={'matched':n,'concordant':k,'share':k/n,'two_sided_binomial_p':pv,'labelled_names':labelled.display_name.tolist()}

def main():
 OUT.mkdir(parents=True,exist_ok=True);df=pd.read_csv(TABLE,low_memory=False);df['effect_size']=pd.to_numeric(df.effect_size);sig=df[df.fdr<.05].copy();print('verified',len(df),'significant',len(sig),flush=True)
 g,gen,baseline=compute_genetic(df);pathways=compute_ora(df);phases=compute_fda(df);approved_inset()
 pairs=[]
 for j,t1 in enumerate(gen.Type):
  for t2 in gen.Type.iloc[j+1:]:
   a=g[g.intervention_type==t1].effect_size;b=g[g.intervention_type==t2].effect_size;u,pv=stats.mannwhitneyu(a,b,alternative='two-sided');pairs.append({'type1':t1,'type2':t2,'U':u,'P_value':pv})
 pairs=pd.DataFrame(pairs);pairs['FDR']=bh(pairs.P_value);save_csv(pairs,'panel_cd_genetic_pairwise_tests.csv')
 H,pv=stats.kruskal(*[p.effect_size for _,p in sig.groupby('condition_category')]);pos=sig.loc[sig.effect_size>0,'effect_size'];neg=-sig.loc[sig.effect_size<0,'effect_size'];T,Tp=stats.ttest_ind(pos,neg);NUM['ridge_global_tests']={'kruskal_H':H,'kruskal_p':pv,'equal_variance_t_abs_effects':T,'t_p':Tp,'higher_n':len(pos),'lower_n':len(neg)}
 life=pd.read_csv(META/'concordance/lifespan_concordance.csv');primary=pd.read_csv(META/'primary/intervention_meta.csv');rep=primary[primary.k>=2].copy();rep['robust']=(rep.fdr<.05)&(rep.rank_bound>0)
 NUM['all_rows']=len(df);NUM['significant_rows']=len(sig);NUM['concordance_total']={'matched':len(life),'concordant':int(life.concordant.sum()),'two_sided_p':binomtest(int(life.concordant.sum()),len(life),.5).pvalue,'constant_higher_matches':int((life.expected>0).sum()),'constant_higher_share':float((life.expected>0).mean())}
 NUM['replication']={'two_or_more_projects':len(rep),'robust':int(rep.robust.sum()),'robust_higher':int((rep.robust&(rep.estimate>0)).sum()),'robust_lower':int((rep.robust&(rep.estimate<0)).sum())}
 # Manual physical boxes track the July a-h composition, with i appended below.
 fig=plt.figure(figsize=(183/25.4,247/25.4));axes={}
 ax=fig.add_axes([.07,.785,.38,.16],projection='polar');axes['a']=ax
 cats=['Disease','Drug','Environment','Genetic','Other'];theta=np.arange(5)*2*np.pi/5;ang=np.r_[theta,theta[0]];ac=df.condition_category.value_counts().reindex(cats,fill_value=0);sc=sig.condition_category.value_counts().reindex(cats,fill_value=0);rmax=int(np.ceil(ac.max()/5000)*5000)
 ax.set_theta_offset(np.pi/2);ax.set_theta_direction(-1);ax.set_ylim(-rmax*.1,rmax*1.14)
 for vals,color,style,label in [(ac,'#B8B8B8','--',f'All ({len(df):,})'),(sc,'#FF6B6B','-',f'FDR < 0.05 ({len(sig):,})')]:
  v=np.r_[vals.values,vals.values[0]];ax.plot(ang,v,color=color,lw=1.2,ls=style,marker='o',ms=3,label=label);ax.fill(ang,v,color=color,alpha=.12)
  for t,c in zip(theta,vals):ax.text(t,c+rmax*(.18 if color=='#FF6B6B' or t==theta[4] else (.22 if t==theta[2] else .06)),f'{c:,}',fontsize=6,color=color,ha='center')
 ax.set_xticks(theta);ax.set_xticklabels(cats,fontsize=7,fontweight='bold');ax.tick_params(pad=0);ax.set_yticks(np.arange(5000,rmax+1,5000));ax.set_yticklabels(['' if v==5000 else f'{v//1000}k' for v in np.arange(5000,rmax+1,5000)],fontsize=7,color='#999999');ax.spines['polar'].set_visible(False);ax.grid(linestyle='--',alpha=.2);ax.legend(loc='upper center',bbox_to_anchor=(.5,-.11),frameon=False,ncol=2,fontsize=7,handlelength=1.1,columnspacing=.8);letter(ax,'a',-.15,1.12);title(ax,'Intervention categories')
 save_csv(pd.DataFrame({'category':cats,'all':ac.values,'FDR_significant':sc.values}),'panel_a_category_counts.csv');NUM['radar']={'all':ac.to_dict(),'significant':sc.to_dict(),'radial_max':rmax}
 # b ridges retain category order based on the verified median and July colors.
 order=sig.groupby('condition_category').effect_size.median().sort_values().index.tolist();b_stats=[];x=np.linspace(sig.effect_size.min(),sig.effect_size.max(),1800)
 for j,cat in enumerate(order):
  ax=fig.add_axes([.60,.797+(4-j)*.033,.335,.041]);axes['b'+str(j)]=ax;p=sig[sig.condition_category==cat].effect_size.dropna();kde=gaussian_kde(p);kde.set_bandwidth(kde.factor*.5);dens=kde(x);ax.fill_between(x,dens,color=COL[cat],alpha=.7);ax.plot(x,dens,color=COL[cat],lw=.9);med=p.median();pv=stats.wilcoxon(p[p!=0]).pvalue;ax.axvline(0,ls='--',lw=.5,color='black');ax.axvline(med,lw=.8,color=COL[cat]);ax.set_xlim(-25,25);ax.set_ylim(0,max(dens)*1.2);ax.set_yticks([]);ax.text(-.035,.5,cat,ha='right',va='center',transform=ax.transAxes,fontweight='bold',fontsize=7);ax.text(-.035,.02,f'n={len(p):,}',ha='right',transform=ax.transAxes,fontsize=6);ax.text(.98,.45,f'P={fmtp(pv)}',ha='right',transform=ax.transAxes,fontsize=6,color='white',bbox={'boxstyle':'round,pad=.18','facecolor':COL[cat],'edgecolor':'none','alpha':.95});ax.spines[['left','top','right']].set_visible(False);ax.grid(False)
  if j<4:ax.spines['bottom'].set_visible(False);ax.set_xticks([])
  else:ax.set_xlabel("Cohen's d",labelpad=0);ax.tick_params(length=2,pad=1)
  b_stats.append({'category':cat,'n':len(p),'median_d':med,'wilcoxon_p':pv})
  if j==0:letter(ax,'b',-.29,1.5);title(ax,'Significant effect distributions')
 save_csv(pd.DataFrame(b_stats),'panel_b_ridge_statistics.csv');NUM['ridge']=b_stats;NUM['ridge_display']={'xlim':[-25,25],'outside_display_n':int((sig.effect_size.abs()>25).sum()),'KDE_uses_all_significant_rows':True}
 # c genetic fold enrichment using July binomial method.
 ax=fig.add_axes([.14,.602,.32,.132]);axes['c']=ax;types=gen.Type.tolist();colors={}
 for r in gen.itertuples():
  intensity=min(abs(r.fold_enrichment-1)/.5,1);colors[r.Type]=plt.cm.Blues(.5+.5*intensity) if r.fold_enrichment>1 else plt.cm.YlOrRd(.5+.5*intensity)
 yy=np.arange(len(gen));ax.barh(yy,gen.fold_enrichment,color=[colors[t] for t in types],alpha=.6,edgecolor='black',lw=.5);ax.axvline(1,color='black',ls='--',lw=.8);ax.set_yticks(yy,types,fontsize=7);ax.set_xlabel('Enrichment of lower-output signatures',labelpad=2);ax.set_xlim(0,max(1.6,gen.fold_enrichment.max()*1.6));base(ax);letter(ax,'c',-.37,1.18)
 for y,r in zip(yy,gen.itertuples()):ax.text(r.fold_enrichment+.025,y,f'{r.fold_enrichment:.2f}  n={r.n}\nP={fmtp(r.p_binomial)}',va='center',fontsize=6)
 ax.text(.0,-.33,f'Baseline lower-output fraction = {baseline:.1%}',transform=ax.transAxes,fontsize=7)
 # d genetic density colors linked to c as July source.
 ax=fig.add_axes([.58,.602,.355,.132]);axes['d']=ax
 for typ in types:
  vals=g[g.intervention_type==typ].effect_size.dropna();xx=np.linspace(vals.min(),vals.max(),500);den=gaussian_kde(vals)(xx);ax.fill_between(xx,den,color=colors[typ],alpha=.2);ax.plot(xx,den,color=colors[typ],lw=1,label=f'{typ} ({len(vals)})')
 ax.axvline(0,color='black',ls='--',lw=.8);ax.set_xlabel("Cohen's d");ax.set_ylabel('Density');ax.legend(frameon=False,fontsize=7,loc='lower left',bbox_to_anchor=(-.02,1.0),ncol=2,columnspacing=.5,handlelength=.7,handletextpad=.3);base(ax);letter(ax,'d',-.15,1.18)
 # Shared colour semantics are categorical by genetic type; legend carries scale.
 # e pathway bars preserve July rocket rank colors with explicit rank scale.
 ax=fig.add_axes([.265,.378,.22,.157]);axes['e']=ax
 if len(pathways):
  import seaborn as sns
  cmap=sns.color_palette('rocket',as_cmap=True);norm=Normalize(1,max(len(pathways),2));n=len(pathways);ax.barh(np.arange(n),pathways.neg_log10_FDR,color=[cmap(norm(v)) for v in range(1,n+1)],edgecolor='black',lw=.3,alpha=.9);ax.set_yticks(np.arange(n),[{'Thyroid hormone signaling pathway':'Thyroid hormone signalling','Th17 cell differentiation':'Th17 cell differentiation','Longevity regulating pathway':'Longevity regulation','Signaling pathways regulating pluripotency of stem cells':'Stem cell pluripotency','AGE-RAGE signaling pathway in diabetic complications':'AGE–RAGE diabetic signalling'}.get(str(v),str(v)) for v in pathways.Pathway],fontsize=7);ax.invert_yaxis();ax.axvline(-np.log10(.05),color='#E64B35',ls='--',lw=.8);ax.set_xlabel('-log10 BH-FDR',labelpad=2);base(ax)
  cax=fig.add_axes([.494,.378,.007,.157]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cax,orientation='vertical',ticks=[1,n]);cb.ax.tick_params(labelsize=7,length=1,pad=1);cb.set_label('Rank',size=7,labelpad=0);cb.solids.set_rasterized(False)
 else:ax.text(.5,.5,'No pathway meets\nthe overlap criteria',ha='center',transform=ax.transAxes);ax.axis('off')
 letter(ax,'e',-1.0,1.1);title(ax,'Pathway enrichment')
 # f divergence counts, original phase colors and direction-opacity split.
 ax=fig.add_axes([.64,.378,.295,.157]);axes['f']=ax;phcol={'FDA_Approved':'#27AE60','Phase_3':'#2874A6','Phase_2':'#3498DB','Phase_1':'#5DADE2','Preclinical':'#95A5A6'};yy=np.arange(len(phases));ax.barh(yy,-phases.lower,color=[phcol[t] for t in phases.status],alpha=.8,edgecolor='black',lw=.5);ax.barh(yy,phases.higher,color=[phcol[t] for t in phases.status],alpha=.4,edgecolor='black',lw=.5);ax.axvline(0,color='black',lw=.8);ax.set_yticks(yy,[t.replace('FDA_Approved','Approved').replace('Phase_','Phase ') for t in phases.status],fontsize=7);ax.set_xlabel('Number of atlas rows',labelpad=2);base(ax);letter(ax,'f',-.32,1.1);title(ax,'Clinical drug development')
 for y,r in zip(yy,phases.itertuples()):
  ax.text(-r.lower/2 if r.lower>=35 else -r.lower-12,y,str(r.lower),ha='center' if r.lower>=35 else 'right',va='center',fontsize=6,color='white' if r.lower>=35 else 'black');ax.text(r.higher/2 if r.higher>=25 else r.higher+12,y,str(r.higher),ha='center' if r.higher>=25 else 'left',va='center',fontsize=6);ax.text(.99,y,f'P={fmtp(r.p_binomial)}',ha='right',va='center',fontsize=6,transform=ax.get_yaxis_transform())
 ax.set_xlim(-210,490)
 ax.text(.0,.98,'Lower output',transform=ax.transAxes,fontsize=6);ax.text(1,.98,'Higher output',transform=ax.transAxes,fontsize=6,ha='right')
 # g/h match only verified mouse interventions; no invented unmatched denominator.
 ax=fig.add_axes([.08,.251,.86,.069]);axes['g']=ax;waffle(ax,'GenAge','g',life,31)
 ax=fig.add_axes([.08,.067,.25,.099]);axes['h']=ax;waffle(ax,'DrugAge','h',life,6)
 # i added replication volcano, using same statistical input as October panel.
 ax=fig.add_axes([.45,.048,.485,.125]);axes['i']=ax;Y=-np.log10(rep.fdr.clip(lower=1e-300));cols=np.where(rep.robust,np.where(rep.estimate<0,AGREE,'#E64B35'),'#B8B8B8');orderp=np.argsort(rep.robust.values);ax.scatter(rep.estimate.iloc[orderp],Y.iloc[orderp],c=cols[orderp],s=(3+2*np.sqrt(rep.k)).iloc[orderp],alpha=.75,edgecolors='none',clip_on=False);ax.axhline(-np.log10(.05),ls='--',lw=.6,color='black');ax.axvline(0,ls='--',lw=.6,color='black');ax.set_xlabel("Pooled Hedges' g across independent projects",labelpad=1);ax.set_ylabel('-log10 BH-FDR',labelpad=1);base(ax);letter(ax,'i',-.13,1.28)
 wanted={'lipopolysaccharide':'LPS','sleep deprivation':'Sleep deprivation','plx5622':'PLX5622','bleomycin':'Bleomycin'};named=rep[(rep.species=='mouse')&rep.canonical_target.isin(wanted)].sort_values('k',ascending=False).drop_duplicates('canonical_target')
 # Labels occupy distinct fixed offsets from observed points.
 offsets={'lipopolysaccharide':(-20,8),'sleep deprivation':(25,20),'plx5622':(-22,10),'bleomycin':(26,7)}
 for r in named.itertuples():ax.annotate(wanted[r.canonical_target],(r.estimate,-np.log10(r.fdr)),xytext=offsets[r.canonical_target],textcoords='offset points',fontsize=6,ha='left',arrowprops={'arrowstyle':'-','lw':.5,'color':'#555555'})
 ax.legend(handles=[Line2D([],[],marker='o',ls='',color='#E64B35',ms=3,label=f'Robust higher ({NUM["replication"]["robust_higher"]})'),Line2D([],[],marker='o',ls='',color=AGREE,ms=3,label=f'Robust lower ({NUM["replication"]["robust_lower"]})'),Line2D([],[],marker='o',ls='',color='#B8B8B8',ms=3,label='Not robust')],loc='upper right',bbox_to_anchor=(1,1.14),ncol=1,frameon=False,fontsize=7,handletextpad=.2,columnspacing=.8)
 ax.text(0,1.26,'Replication across projects',transform=ax.transAxes,fontsize=8,fontweight='bold');ax.text(0,1.12,f'n = {len(rep)} interventions',transform=ax.transAxes,fontsize=7);save_csv(rep,'panel_i_replication_points.csv');save_csv(named,'panel_i_named_points.csv')
 # Record all printed values and immutable input hashes alongside render.
 (OUT/'numbers.json').write_text(json.dumps(NUM,indent=2,default=str)+'\n');INPUTS.extend(sorted(SCRIPT.glob('*figure4*2026_11_05.py')));provenance={str(p):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'mtime_utc':datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat()} for p in INPUTS};(OUT/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 # Only defect checks are used; retain July composition/palette.
 fig.canvas.draw();small=[]
 for text in fig.findobj(matplotlib.text.Text):
  if text.get_visible() and text.get_text() and text.get_fontsize()<6:small.append(text.get_text())
 if small:raise RuntimeError('Text below annotation floor: '+repr(small))
 sys.path.insert(0,'pipeline_assets/figure_style_2026_09_10')
 from nature_figure_utils import strict_figure_qc
 qc=strict_figure_qc(fig,axes,check_text_figure_edges=False,forbid_uniform_grid=False,check_text_spine=False,check_data_occupancy=False,check_max_fontsize=False,check_matched_limits=False,check_smoother_domain=False,check_house_style=False,max_tall_panel_aspect=100,max_sparse_panel_aspect=100)
 (OUT/'defect_check.json').write_text(json.dumps(qc,indent=2,default=str)+'\n')
 fig.canvas.draw();rdr=fig.canvas.get_renderer();edge_issues=[]
 for txt in fig.findobj(matplotlib.text.Text):
  if not txt.get_visible() or not txt.get_text():continue
  bb=txt.get_window_extent(rdr);box=fig.bbox
  if bb.x0 < box.x0-1 or bb.x1 > box.x1+1 or bb.y0 < box.y0-1 or bb.y1 > box.y1+1:edge_issues.append(txt.get_text())
 (OUT/'canvas_edge_check.json').write_text(json.dumps({'outside_canvas':edge_issues},indent=2)+'\n')
 fig.savefig(OUT/'figure4.pdf');fig.savefig(OUT/'figure4_review.png',dpi=200)
 import fitz
 checkdoc=fitz.open(OUT/'figure4.pdf');spans=[span for block in checkdoc[0].get_text('dict')['blocks'] if 'lines' in block for line in block['lines'] for span in line['spans']];outside=[span['text'] for span in spans if span['bbox'][0]<-1 or span['bbox'][1]<-1 or span['bbox'][2]>checkdoc[0].rect.width+1 or span['bbox'][3]>checkdoc[0].rect.height+1];(OUT/'pdf_checks.json').write_text(json.dumps({'page_size_points':list(checkdoc[0].rect),'raster_image_count':len(checkdoc[0].get_images()),'smallest_extracted_font_pt':min(span['size'] for span in spans),'outside_page_text':outside,'image_title_absent':not any('Figure 4' in span['text'] for span in spans)},indent=2)+'\n')
 doc=fitz.open(ROOT/'paper/revision_2026-09-07/baseline/03_Figures_2026-07-06.pdf');page=doc[3]
 # Compare the July figure content, removing package margins/header only.
 page.get_pixmap(matrix=fitz.Matrix(2,2),clip=fitz.Rect(72,101,552,663)).save(OUT/'july_figure4_content.png')
 from PIL import Image,ImageOps,ImageDraw
 old=Image.open(OUT/'july_figure4_content.png').convert('RGB');new=Image.open(OUT/'figure4_review.png').convert('RGB');w=new.width;old=old.resize((w,round(old.height*w/old.width)));height=max(old.height,new.height);both=Image.new('RGB',(w*2+40,height+40),'white');both.paste(old,(0,40));both.paste(new,(w+40,40));draw=ImageDraw.Draw(both);draw.text((12,8),'July reference',fill='black');draw.text((w+52,8),'Verified restoration',fill='black');both.save(OUT/'figure4_vs_july.png')
 print(json.dumps(NUM,indent=2,default=str),flush=True)

if __name__=='__main__':
 try:main()
 except Exception as e:
  print(f'Figure 4 restoration failed: {type(e).__name__}: {e}',file=sys.stderr);raise
