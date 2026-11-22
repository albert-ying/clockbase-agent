"""Restore replication panels and consolidate Figure 4b/d into ED22."""
from ed21_24_restore_common_2026_11_05 import *
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
from matplotlib.ticker import MaxNLocator
import build_ed22_replication_2026_10_14 as original


def names(rows, *, project_counts=False, flags=False):
    rows=original.with_category_suffix(rows)
    labels=[]
    for r in rows.itertuples():
        text=original.label(r).replace(' (mm)','').replace(' (hs)',' (human)')
        if not flags: text=text.replace('†','')
        text=text.replace(', environment',', env.').replace(', other',', other')
        if project_counts: text+=f' ({r.k})'
        labels.append(text)
    return labels


def forest(ax, rows, color):
    yy=np.arange(len(rows))[::-1]
    ax.barh(yy, rows.estimate, height=.52,color=color,alpha=.8,
            edgecolor=INK,lw=.45,zorder=2)
    ax.errorbar(rows.estimate, yy,xerr=[rows.estimate-rows.ci_low,rows.ci_high-rows.estimate],
                fmt='o',ms=4,color=INK,ecolor=INK,lw=.75,capsize=2,capthick=.6,zorder=3)
    ax.set_yticks(yy,names(rows,flags=True));ax.set_ylim(-.6,len(rows)-.25)
    ax.tick_params(axis='y',length=0);grid(ax)
    ax.axvline(0,color=INK,lw=.7,ls='--',zorder=1)
    ax.set_xlabel("Pooled Hedges' g (95% CI)")
    lo=min(rows.ci_low.min(),0);hi=max(rows.ci_high.max(),0)
    ax.set_xlim(np.floor(lo)-.1,np.ceil(hi)+.1)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    return yy


def projects(ax, units, primary, project_records, key, color):
    part=units[(units[original.KEYS]==pd.Series(dict(zip(original.KEYS,key)))).all(axis=1)].sort_values('g').copy()
    row=primary[(primary[original.KEYS]==pd.Series(dict(zip(original.KEYS,key)))).all(axis=1)].iloc[0]
    yy=np.arange(len(part))[::-1]+1.3;half=1.959964*np.sqrt(part['var'])
    ax.errorbar(part.g,yy,xerr=half,fmt='s',ms=4,mec=INK,mew=.45,color=color,
                ecolor=color,lw=1.0,capsize=2,capthick=.6,zorder=3)
    label=[]
    for r in part.itertuples():
        series=r.series.split(';');who=original.lab(project_records['records'].get(series[0]))
        more=f' +{len(series)-1}' if len(series)>1 else ''
        label.append(f'{who + ", " if who else ""}{series[0]}{more}')
    ax.add_patch(Polygon([[row.ci_low,0],[row.estimate,.29],[row.ci_high,0],[row.estimate,-.29]],
                         fc=color,ec=INK,lw=.6,zorder=3))
    ax.set_yticks(list(yy)+[0],label+[f'Pooled ({row.k} projects)'])
    ax.set_ylim(-.65,max(yy)+.7);ax.tick_params(axis='y',length=0)
    ax.axvline(0,color=INK,lw=.7,ls='--');grid(ax)
    ax.set_xlabel("Hedges' g (95% CI)")
    part['ci_low_plotted']=part.g-half;part['ci_high_plotted']=part.g+half
    part['display_label']=label
    return part,row,(min((part.g-half).min(),row.ci_low),max((part.g+half).max(),row.ci_high))


def main():
    style();tables,units,project_records=original.load(META)
    primary=tables['primary'];rob=original.robust(primary)
    lower=rob[rob.estimate<0].head(original.TOP).copy();higher=rob[rob.estimate>0].head(original.TOP).copy()
    fig=plt.figure(figsize=(WIDTH,9.7))
    a=axes(fig,[.265,.813,.23,.105],'Lower clock outputs','a',.015)
    b=axes(fig,[.525,.813,.135,.105],None,'b',.515)
    c=axes(fig,[.285,.525,.245,.215],'Higher clock outputs','c',.015)
    d=axes(fig,[.805,.835,.17,.095],'Counting unit','d',.69)
    e=axes(fig,[.285,.315,.245,.145],'Sleep deprivation','e',.015)
    f=axes(fig,[.285,.07,.245,.18],'Rapamycin','f',.015)
    g=axes(fig,[.775,.09,.205,.64],'Lifespan matches','g',.69)
    yy=forest(a,lower,LOWER);forest(c,higher,HIGHER)
    # Counts remain linked to each row without competing with the CI drawing.
    for ax,rows in [(a,lower),(c,higher)]:
        for y,r in zip(np.arange(len(rows))[::-1],rows.itertuples()):
            ax.text(1.025,y,f'{r.k} ({r.n_series})',transform=ax.get_yaxis_transform(),
                    fontsize=6.5,ha='left',va='center')
    # Only panel c has a separate count column; the aligned sensitivity takes that area in a.
    for t in list(a.texts): t.remove()
    lower_labels=[n+f' [{r.k}/{r.n_series}]' for n,r in zip(names(lower,flags=True),lower.itertuples())]
    a.set_yticklabels(lower_labels)
    for t in list(c.texts):t.remove()
    c.set_yticklabels([n+f' [{r.k}/{r.n_series}]' for n,r in zip(names(higher,flags=True),higher.itertuples())])
    a.set_xticks([-3,-2,-1,0])
    columns=list(original.RUNS)+['hartung-knapp'];sens=[]
    for j,name in enumerate(columns):
        table=tables['primary' if name=='hartung-knapp' else name]
        for y,r in zip(yy,lower.itertuples()):
            hit=table[(table[original.KEYS]==pd.Series({k:getattr(r,k) for k in original.KEYS})).all(axis=1)]
            if hit.empty or hit.iloc[0].k<2:
                b.plot([j-.16,j+.16],[y,y],color=GRAY,lw=1)
                state='under two units';h=None
            else:
                h=hit.iloc[0]
                ok=h.ci_high_hk<0 if name=='hartung-knapp' else h.fdr<.05 and h.rank_bound>0
                state='robust' if ok else 'not robust'
                b.scatter(j,y,s=26,marker='o',fc=LOWER if ok else 'white',ec=LOWER,lw=.7,zorder=3)
            sens.append({'intervention':original.label(r),'specification':name,'state':state,
                         'k':None if h is None else int(h.k),
                         'fdr':None if h is None else float(h.fdr),
                         'ci_high_hk':None if h is None else float(h.ci_high_hk)})
    b.set_xticks(range(5),['Primary','Series','All contexts','No contact','H–K'],rotation=90)
    b.xaxis.tick_top();b.tick_params(axis='x',length=0,pad=3)
    b.set_yticks([]);b.set_xlim(-.5,4.5);b.set_ylim(a.get_ylim());b.spines[:].set_visible(False)
    b.legend(handles=[Line2D([],[],marker='o',ls='',mfc=LOWER,mec=LOWER,ms=4,label='Robust'),
                      Line2D([],[],marker='o',ls='',mfc='white',mec=LOWER,ms=4,label='Not robust'),
                      Line2D([],[],ls='-',color=GRAY,lw=1,label='< 2 units')],
             loc='upper left',bbox_to_anchor=(-.14,-.06),frameon=False,handlelength=.9,
             handletextpad=.3,labelspacing=.25)
    series=tables['by-series']
    groups=[('Replicated',int((series.k>=2).sum()),int((primary.k>=2).sum())),
            ('Robust',len(original.robust(series)),len(rob)),
            ('Higher',int((original.robust(series).estimate>0).sum()),int((rob.estimate>0).sum())),
            ('Lower',int((original.robust(series).estimate<0).sum()),int((rob.estimate<0).sum()))]
    yd=np.arange(len(groups))[::-1]
    for off,idx,color,title in [(.18,1,CLASSIC[0],'Series'),(-.18,2,CLASSIC[1],'Projects')]:
        vv=[group[idx] for group in groups]
        d.barh(yd+off,np.array(vv)-1,left=1,height=.34,color=color,ec=INK,lw=.45,label=title)
        for y,value in zip(yd+off,vv):d.text(value*1.16,y,str(value),fontsize=6.5,va='center')
    d.set_xscale('log');d.set_xlim(1,1900);d.set_xticks([1,10,100,1000],['1','10','100','1,000'])
    d.minorticks_off();d.set_yticks(yd,[row[0] for row in groups]);d.tick_params(axis='y',length=0)
    d.set_ylim(-.6,len(groups)-.4);d.set_xlabel('Interventions');grid(d)
    d.legend(loc='lower center',bbox_to_anchor=(.4,1.29),ncol=2,frameon=False,
             handlelength=.85,columnspacing=.6,handletextpad=.3)
    most=rob[rob.estimate<0].sort_values(['k','rank_bound'],ascending=False).iloc[0]
    ep,er,ex=projects(e,units,primary,project_records,tuple(most[original.KEYS]),LOWER)
    fp,fr,fx=projects(f,units,primary,project_records,original.REFERENCE,LOWER)
    for ax in [e,f]:ax.set_xlim(np.floor(min(ex[0],fx[0]))-.1,np.ceil(max(ex[1],fx[1]))+.1);ax.set_xticks([-4,-2,0,2])
    lifespan_file=META/'concordance/lifespan_concordance.csv'
    matched=original.with_category_suffix(pd.read_csv(lifespan_file))
    matched=matched[matched.k>=2].sort_values(['expected','estimate'],ascending=[False,True]).copy()
    yg=np.arange(len(matched))[::-1]
    for y,r in zip(yg,matched.itertuples()):
        color=HIGHER if r.expected>0 else LOWER
        g.errorbar([r.estimate],[y],xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],
                   fmt='o',ms=4.7,mfc=color if r.concordant else 'white',mec=color,mew=.8,
                   ecolor=color,lw=1,capsize=2,capthick=.7,zorder=3)
    g.set_yticks(yg,names(matched,project_counts=True));g.set_ylim(-.7,len(matched)-.3)
    for tick,r in zip(g.get_yticklabels(),matched.itertuples()):
        if r.source=='GenAge':tick.set_fontstyle('italic')
    g.set_xlim(-6.1,6.1);g.set_xticks([-6,-3,0,3,6]);g.tick_params(axis='y',length=0)
    g.axvline(0,color=INK,ls='--',lw=.75);grid(g);g.set_xlabel("Pooled Hedges' g (95% CI)")
    fig.legend(handles=[Line2D([],[],marker='o',ls='',mfc=HIGHER,mec=HIGHER,ms=4,label='Higher'),
                        Line2D([],[],marker='o',ls='',mfc=LOWER,mec=LOWER,ms=4,label='Lower'),
                        Line2D([],[],marker='o',ls='',mfc='white',mec=INK,ms=4,label='Opposite')],
               title='Expected direction',title_fontsize=7,ncol=3,loc='lower center',
               bbox_to_anchor=(.79,.008),frameon=False,handletextpad=.2,columnspacing=.6,labelspacing=.3)
    facts={'counting_unit':[{'group':n,'series':s,'projects':p} for n,s,p in groups],
           'lower':lower.to_dict('records'),'higher':higher.to_dict('records'),
           'sensitivity':sens,'project_sleep_deprivation':ep.to_dict('records'),
           'sleep_pooled':er.to_dict(),'project_rapamycin':fp.to_dict('records'),
           'rapamycin_pooled':fr.to_dict(),'lifespan_matched':matched.to_dict('records')}
    sources=[META/k/'intervention_meta.csv' for k in original.RUNS]
    sources += [META/'primary/unit_effects.csv',META/'primary/meta_input_effects.csv',META/'projects.json',lifespan_file]
    out=save(fig,'ed22',dict(zip('abcdefg',[a,b,c,d,e,f,g])),4,sources=sources,facts=facts,script=__file__)
    for name,df in [('panel_a_lower',lower),('panel_b_sensitivity',pd.DataFrame(sens)),('panel_c_higher',higher),
                    ('panel_d_counting_unit',pd.DataFrame(facts['counting_unit'])),('panel_e_sleep_deprivation',ep),
                    ('panel_f_rapamycin',fp),('panel_g_lifespan_matched',matched)]:
        df.to_csv(out/f'{name}.csv',index=False)
    details=matched[['canonical_target','manipulation','estimate','ci_low','ci_high','k','expected','concordant']]
    readme(out,f'''# Intervention replication

ED22 retains October panels a–f and adds the matched-intervention forest from
Figure 4d as panel g. The existing counting-unit panel d consolidates Figure
4b and adds the robust total. July Figure 3 supplies outlined bars and the
orange/blue direction key. July Figure 4 supplies density and category-colour
conventions. The comparison uses July page 4 at the same physical scale.

## Panel map

| Original panel | Restored panel | Content |
| --- | --- | --- |
| ED22a | a | All four robust lower-output interventions, with bars and exact 95% CIs |
| ED22b | b | Five-specification sensitivity of the same four rows |
| ED22c | c | The same twelve highest-ranked robust higher-output interventions |
| ED22d and Figure 4b | d | Series-versus-project counts, with robust total added |
| ED22e | e | Sleep deprivation by independent project and pooled estimate |
| ED22f | f | Rapamycin by independent project and pooled estimate |
| Figure 4d | g | All {len(matched)} lifespan-matched rows with at least two projects |

## Data and numbers

All panels read `intervention-meta-2026-10-17`. Absolute paths, modification
times and input hashes are in `provenance.json`; every plotted estimate, CI,
sample-unit count and decision is in `numbers.json` and the panel CSV files.
Axis ticks are in the stated units. Point coordinates are never simulated.

Panel a/c read `primary/intervention_meta.csv`, keep k ≥ 2, BH-FDR < 0.05
and rank_bound > 0, split by estimate sign, and retain the first twelve rows
per direction in the source ordering. Four rows meet the lower-output
criteria. The primary source ordering ranks by the CI bound nearest zero.
The bracketed counts in a and c are projects/series. A dagger retains the October large-arm flag
(an input arm exceeds 120 observations), including the possibility of
single-cell measurements rather than independently sampled individuals.

Panel b reads the primary, by-series, all-contexts and no-contact-link
`intervention_meta.csv` files. Filled dots meet k ≥ 2, FDR < 0.05 and
rank_bound > 0. The H–K column uses primary `ci_high_hk < 0`; it does not
claim a separate FDR correction. Open dots fail the relevant criterion and
horizontal strokes indicate fewer than two units.

Panel d gives these exact counts:

{table_md(pd.DataFrame(facts['counting_unit']))}

The robust totals are 88 by series and 77 by projects, partitioned into 84/73
higher and 4/4 lower. Replicated means at least two counting units. Bars
start at 1 on a logarithmic axis, as in both October source panels.

Panels e/f read `primary/unit_effects.csv` (`g`, `var`, `series`) and
`projects.json` (GEO contact surnames), with pooled estimates from
`primary/intervention_meta.csv`. Project 95% CIs use
g ± 1.959964 × sqrt(var), exactly as in October. Panel e is
{original.label(most)} with {int(er.k)} projects, not a substituted example.
Panel f is rapamycin with {int(fr.k)} projects. Pooled estimates and CIs are
{er.estimate:.9g} [{er.ci_low:.9g}, {er.ci_high:.9g}] and
{fr.estimate:.9g} [{fr.ci_low:.9g}, {fr.ci_high:.9g}], respectively.
The two plots share their horizontal limits.

Panel g reads `concordance/lifespan_concordance.csv`, keeps k ≥ 2 and sorts
by expected direction then estimate. Unlike the former clipped forest, all
26 CIs fit within −6 to 6. Gene labels are italic; project counts appear in
parentheses. `expected` = +1 predicts higher output, −1 predicts lower output.

{table_md(details)}

## Draft legend

**a,** Robust replicated interventions associated with lower clock outputs,
ranked by the 95% confidence bound nearest zero. Bars show pooled Hedges' g
and points and whiskers show estimates and 95% CIs. The four rows pass
BH-FDR < 0.05 and have at least two independent projects; blue indicates
lower output. Brackets give projects/series. **b,** Sensitivity of these four
rows to the counting unit, context inclusion, project-linkage rule and
Hartung–Knapp interval. Filled dots indicate robustness under the specified
criterion, open dots its absence, and strokes fewer than two units.
**c,** The twelve highest-ranked robust replicated interventions with higher
outputs, drawn as in a in orange. Brackets give projects/series.
A dagger marks an intervention with an input arm exceeding 120 observations.
**d,** Counts of replicated and robust intervention–category combinations
when the independent unit is a GEO series (blue) or linked project (orange).
The horizontal scale is logarithmic. **e,f,** Project-level estimates for
sleep deprivation ({int(er.k)} projects) and rapamycin ({int(fr.k)} projects).
Squares and whiskers show project Hedges' g and normal-approximation 95% CIs;
diamonds show random-effects pooled estimates and 95% CIs. **g,** Pooled
effects for {len(matched)} GenAge- or DrugAge-matched intervention–category
combinations represented by at least two projects. Colours indicate the
clock-output direction predicted from the lifespan database. Filled points
agree with that direction and open points have the opposite sign; whiskers
show 95% CIs. Parentheses give the project count. Gene names are italicised.
All effects are case minus control. Normal random-effects tests and BH-FDR
follow the supplied intervention meta-analysis; H–K denotes the separate
interval sensitivity, not an additional FDR test.

## Checks

`checks.json` records the nature-figure defect checks, canvas bounds, text
floor, exact PDF size, fonts and vector-path inspection. The page is
183 × 246.38 mm. The source panel rows, sensitivity states and counts are
preserved; the counting panel is consolidated rather than duplicated.
The 200-dpi render and the July comparison are inspected for label collisions.
No original script, manuscript or combined figure package is changed.
''')


if __name__=='__main__':
    main()
