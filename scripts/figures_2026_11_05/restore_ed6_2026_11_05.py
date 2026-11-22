#!/usr/bin/env python3
"""Restore July ED6 geometry from saved predictions and registered OLS models.
The July plotting script is absent; the archived September builder is copied
alongside this file. Its data sources, not its visual theme, are reused.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pymupdf as fitz
from PIL import Image,ImageOps,ImageDraw
from scipy.stats import t as student_t
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'paper/revision_2026-09-07'
SRC=BASE/'evidence/gse162005-covariate-audit-2026-09-18'
OUT=BASE/'evidence/figure-fixes-2026-11-17/ed6'
OUT.mkdir(parents=True,exist_ok=True)

def main():
 p=pd.read_csv(SRC/'prepared-run-001/saved_predictions.csv').set_index('sample_id')
 d=pd.read_csv(SRC/'numerical-run-001/sample_design_and_predictions.csv').set_index('sample_id').loc[p.index]
 assert len(p)==15 and p.index.is_unique
 cols=[c for c in p if c.endswith('_EN_tAge')]
 x=p[cols].to_numpy();z=(x-x.mean(axis=0))/x.std(axis=0,ddof=0)
 u,s,v=np.linalg.svd(z,full_matrices=False)
 scores=u[:,:2]*s[:2]
 # Deterministic orientation matching July: first library positive on both axes.
 for j in range(2):
  if scores[0,j]<0:scores[:,j]*=-1
 pd.DataFrame(scores,index=p.index,columns=['PC1','PC2']).join(d[['batch','pool','treatment']]).to_csv(OUT/'ed6_pca_scores.csv')
 models=pd.read_csv(SRC/'numerical-run-001/conditional_models.csv')
 for r in models.itertuples():
  parts=[np.ones((len(d),1)),d[['treatment']].to_numpy()]
  for term in ['batch','pool']:
   if f'C({term})' in r.formula:parts.append(pd.get_dummies(d[term],drop_first=True,dtype=float).to_numpy())
  design=np.column_stack(parts);response=d[r.clock].to_numpy();pinv=np.linalg.pinv(design);coef=pinv@response
  df=len(d)-np.linalg.matrix_rank(design);sigma=((response-design@coef)**2).sum()/df;se=np.sqrt(sigma*(pinv@pinv.T)[1,1]);ci=coef[1]+np.array([-1,1])*student_t.ppf(.975,df)*se
  assert abs(coef[1]-r.condition_coefficient)<1e-12
  assert np.max(abs(ci-[r.ci_lower,r.ci_upper]))<1e-9
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.7,'pdf.fonttype':42,'ps.fonttype':42})
 fig=plt.figure(figsize=(183/25.4,3.2))
 ax=fig.add_axes([.085,.26,.365,.62]);b=fig.add_axes([.625,.26,.355,.62])
 batchcolors=['#1f77b4','#ff7f0e','#2ca02c']
 # July uses muted semi-transparent default colours with black outlines.
 for batch in [1,2,3]:
  for treatment,mark in [(0,'o'),(1,'^')]:
   mask=(d.batch==batch)&(d.treatment==treatment)
   ax.scatter(scores[mask,0],scores[mask,1],s=25,c=batchcolors[batch-1],marker=mark,edgecolor='#333333',linewidth=.6,alpha=.8)
 handles=[Line2D([],[],marker=m,linestyle='',color=batchcolors[k-1],markeredgecolor='#333333',alpha=.8,label=f'Batch {k} | {cond}',markersize=4) for k in [1,2,3] for m,cond in [('o','Indolent'),('^','Aggressive')]]
 ax.legend(handles=handles,ncol=1,loc='lower left',bbox_to_anchor=(.07,0),frameon=False,handletextpad=.3,fontsize=7, borderaxespad=.4, labelspacing=.25)
 ax.axhline(0,color='#bdbdbd',lw=.6,zorder=0);ax.axvline(0,color='#bdbdbd',lw=.6,zorder=0)
 ax.set(xlabel='PC1',ylabel='PC2',xlim=(-10.2,7),ylim=(-8.8,6.8))
 ax.set_title('Clock-space PCA',loc='left',pad=8)
 clocks=['Hazard_Mouse_All_All_scaleddiff_EN_tAge','Hazard_Rodents_All_All_scaleddiff_EN_tAge']
 formulas=list(dict.fromkeys(models.formula))
 labels=['Condition only','With batch','With pool','With batch + pool']
 colors=['#888888','#3774bb','#dd8b27','#bb382d']
 offsets=[-.30,-.10,.10,.30]
 for k,c in enumerate(clocks):
  for j,f in enumerate(formulas):
   r=models[(models.clock==c)&(models.formula==f)].iloc[0];y=k+offsets[j]
   b.errorbar(r.condition_coefficient,y,xerr=[[r.condition_coefficient-r.ci_lower],[r.ci_upper-r.condition_coefficient]],fmt='o',color=colors[j],markersize=3.8,elinewidth=1,capsize=2,label=labels[j] if k==0 else None)
 b.axvline(0,color='#888888',ls='--',lw=.7,zorder=0)
 b.set(yticks=[0,1],yticklabels=['Mouse\nmortality','Rodent\nmortality'],ylim=(1.7,-.62),xlim=(-.4,.25),xlabel='Lesion-class coefficient (native clock units)')
 b.set_xticks([-.4,-.2,0,.2]);b.set_title('Effect estimates across covariate models',loc='left',pad=8)
 b.legend(loc='lower center',bbox_to_anchor=(.5,-.39),ncol=2,frameon=False,columnspacing=1,handletextpad=.4,fontsize=7)
 for letter,axis in [('a',ax),('b',b)]:axis.text(-.20,1.08,letter,transform=axis.transAxes,fontweight='bold',fontsize=11)
 fig.savefig(OUT/'ed6.pdf');fig.savefig(OUT/'ed6_review.png',dpi=200)
 fig.canvas.draw()
 import sys
 sys.path.insert(0,str(Path('pipeline_assets/figure_style_2026_09_10')))
 from nature_figure_utils import detect_overlaps
 overlaps=detect_overlaps(fig,{'a':ax,'b':b},padding=1)
 (OUT/'ed6_overlap_check.json').write_text(json.dumps(overlaps,indent=2,default=str))
 texts=[t for t in fig.findobj(matplotlib.text.Text) if t.get_text()]
 # Axis/legend floors; in-panel legend annotation retained at 6 pt to match July.
 assert all(t.get_fontsize()>=6 for t in texts)
 page=fitz.open(OUT/'ed6.pdf')[0];assert abs(page.rect.width-183/25.4*72)<.01;assert not page.get_images()
 july=fitz.open(BASE/'baseline/03_Figures_2026-07-06.pdf')[11]
 # Use the July figure region, not the enclosing package's empty paper margins.
 old=Image.open(__import__('io').BytesIO(july.get_pixmap(matrix=fitz.Matrix(200/72,200/72),clip=fitz.Rect(70,100,550,300)).tobytes('png'))).convert('RGB')
 new=Image.open(OUT/'ed6_review.png').convert('RGB')
 # Equal 183-mm figure widths for comparison.
 old=old.resize((new.width,round(old.height*new.width/old.width)))
 canvas=Image.new('RGB',(new.width*2,max(old.height,new.height)+40),'white');canvas.paste(old,(0,40));canvas.paste(new,(new.width,40));draw=ImageDraw.Draw(canvas);draw.text((10,12),'July page 12',(0,0,0));draw.text((new.width+10,12),'Restored ED6',(0,0,0));canvas.save(OUT/'ed6_vs_july.png')
 meta={'n_libraries':len(d),'n_aggressive':int(d.treatment.sum()),'n_indolent':int((1-d.treatment).sum()),'n_predictors_pca':len(cols),'pca_predictors':cols,'pca_variance_fraction':(s[:2]**2/(s@s)).tolist(),'batch_counts':d.batch.value_counts().sort_index().to_dict(),'model_verification':'8 coefficients independently refit to 1e-12; 16 interval endpoints to 1e-9','pdf_images':len(page.get_images()),'size_mm':[183,3.2*25.4]}
 (OUT/'ed6_checks.json').write_text(json.dumps(meta,indent=2))
 models.to_csv(OUT/'ed6_model_values.csv',index=False)
 print(json.dumps(meta,indent=2))
if __name__=='__main__':
 try:main()
 except Exception as exc:raise RuntimeError(f'ED6 restoration failed: {exc}') from exc
