"""Restore ED21 with July validation colours and Figure 3 modality symbols."""
from ed21_24_restore_common_2026_11_05 import *
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import build_ed21_2026_10_03 as original


def stacks(ax, counts, labels, spec, *, print_segments=True):
    yy = np.arange(len(labels))[::-1]
    rows = []
    for y, (key, row), label in zip(yy, counts.items(), labels):
        total = sum(row.values()); left = 0
        for name, title, color in spec:
            count = row.get(name, 0); width = 100*count/total
            if count:
                ax.barh(y, width, left=left, height=.67, color=color,
                        edgecolor=INK, linewidth=.4)
                if print_segments and width >= 9:
                    ax.text(left+width/2, y, f'{count:,}', ha='center', va='center',
                            fontsize=6.5, color='white' if name in ('VERIFIED', 'CONFIRMED_CLAIMED_DIRECTION') else INK)
            rows.append({'group': key, 'class': name, 'n': count, 'percent': width, 'total': total})
            left += width
        ax.text(102, y, f'{total:,}', va='center', fontsize=7)
    ax.set_yticks(yy, labels)
    ax.set_xlim(0, 116); ax.set_xticks([0,25,50,75,100]); ax.set_ylim(-.55,len(labels)-.45)
    ax.spines['bottom'].set_bounds(0,100)
    ax.tick_params(axis='y', length=0)
    return rows


def main():
    style(); original.FINAL = FINAL
    counts = original.class_counts(); fates=original.fate_counts()
    census_file = FINAL/'error-census/census.json'
    census=read_json(census_file)['modalities']
    classes=[('VERIFIED','Verified',CLASSIC[0]),('PROVISIONAL','Provisional','#9fc2df'),
             ('DESCRIPTIVE','Descriptive',CLASSIC[2]),('DUPLICATE','Duplicate',CLASSIC[3]),
             ('PENDING','Withheld',CLASSIC[4]),('NOT_REPAIRABLE','Not repairable',CLASSIC[1])]
    fate_spec=[('CONFIRMED_CLAIMED_DIRECTION','Confirmed',CLASSIC[0]),
               ('NOT_SIGNIFICANT_AFTER_REPAIR','Not significant','#b9d3e7'),
               ('DESCRIPTIVE_ONLY_ARM_BELOW_3','Descriptive',CLASSIC[2]),
               ('PENDING','Pending',CLASSIC[4]),
               ('CONTRADICTS_CLAIMED_DIRECTION','Reversed',HIGHER),
               ('NO_EXACT_ROW','No atlas row','white')]
    fig=plt.figure(figsize=(WIDTH,6.6))
    a=axes(fig,[.17,.805,.745,.095],'Final classification', 'a',.015)
    b=axes(fig,[.285,.425,.255,.285],'Historical defects', 'b',.015)
    c=axes(fig,[.705,.485,.225,.215],'Published results', 'c',.575)
    d=axes(fig,[.095,.06,.36,.23],'GSE73950', 'd',.015)
    e=axes(fig,[.605,.06,.36,.23],'DNAm input integrity', 'e',.525)
    rows_a=stacks(a,counts,[x[1] for x in original.MODALITIES],classes)
    a.set_xlabel('Atlas rows (%)')
    fig.legend(handles=[Patch(fc=color,ec=INK,lw=.4,label=title) for _,title,color in classes],
             ncol=3,frameon=False,loc='upper center',bbox_to_anchor=(.54,.997),
             columnspacing=1.2,handlelength=1.1,labelspacing=.4)
    rows_b=[]; yy=np.arange(len(original.DEFECTS))[::-1]
    defect_labels=['Stored P missing, 0 or 1','Arms swapped in two models','Effect not reproduced',
                   'Empty or artefactual input','Not a group contrast','Detached DNAm profiles',
                   'Non-bisulfite assay','FDR below raw P']
    for off,(key,label) in zip((.21,0,-.21),original.MODALITIES):
        _,marker,color=MODS[key]
        values=[100*census[key]['defects'][code]['fraction'] for code,_ in original.DEFECTS]
        b.scatter(values,yy+off,s=22,marker=marker,c=color,edgecolors=INK,lw=.35,zorder=3,clip_on=False,label=label)
        for code,value in zip([x[0] for x in original.DEFECTS],values):
            rows_b.append({'modality':key,'defect':code,**census[key]['defects'][code],'percent':value})
    b.set_yticks(yy,defect_labels);b.set_ylim(-.6,len(yy)-.4);b.set_xlim(-.2,10)
    b.set_xticks([0,2,4,6,8,10]);b.set_xlabel('Rows with defect (%)');grid(b)
    b.legend(frameon=False,loc='upper left',bbox_to_anchor=(-.88,-.19),ncol=3,
             columnspacing=.65,handletextpad=.25,handlelength=1)
    rows_c=stacks(c,fates,[x[1] for x in original.OUTPUTS],fate_spec,print_segments=False)
    c.set_xlabel('Published items (%)');c.set_xticks([0,50,100])
    c.legend(handles=[Patch(fc=color,ec=INK,lw=.4,label=title) for _,title,color in fate_spec],
             ncol=2,frameon=False,loc='upper center',bbox_to_anchor=(.40,-.26),
             columnspacing=.8,handlelength=.9,labelspacing=.35)
    screen=pd.read_csv(original.SCREEN)
    study=screen[screen.gse_id=='GSE73950'].set_index('gsm')
    metadata_file=paths.GEO/'GEO_hs_dnam/metadata/GSE73950.csv'
    meta=pd.read_csv(metadata_file,low_memory=False).set_index('geo_accession')
    arm=meta.reindex(study.index).title.str[0].map({'E':'Endometriosis','H':'Control'})
    historical=pd.read_csv(original.INDEX_STUDY).set_index('gsm_id')['Horvathv1'].reindex(study.index)
    case_n={}
    for name,color in [('Endometriosis',HIGHER),('Control',LOWER)]:
        mask=arm==name;case_n[name]=int(mask.sum())
        d.scatter(study.r_reference[mask],historical[mask],s=23,c=color,edgecolors=INK,lw=.35,
                  marker='^',alpha=.8,zorder=3,clip_on=False,label=f'{name} (n = {mask.sum()})')
    d.set_xlim(-.15,1.03);d.set_ylim(-12,67);d.set_xticks([0,.5,1]);d.set_yticks([0,20,40,60])
    d.set_xlabel('Reference methylome correlation (r)');d.set_ylabel('Historical Horvath age (years)');grid(d)
    d.legend(frameon=False,loc='upper left',bbox_to_anchor=(-.015,1.00),handletextpad=.25,labelspacing=.35)
    ctrl=arm=='Control'; rmed=float(study.r_reference[ctrl].median()); amed=float(historical[ctrl].median())
    d.annotate(f'{ctrl.sum()} controls overlap\nr = {rmed:.3f}; {amed:.1f} years',xy=(rmed,amed),xytext=(.19,6),
               fontsize=6.5,color=LOWER,arrowprops={'arrowstyle':'-','color':LOWER,'lw':.6})
    r=screen.r_reference.dropna().to_numpy(); bins=np.linspace(-.45,1,59)
    e.hist(r,bins=bins,color=CLASSIC[0],edgecolor='white',lw=.35,log=True)
    e.hist(r[r<.10],bins=bins,color=HIGHER,edgecolor='white',lw=.35,log=True)
    for threshold in [.1,.5]: e.axvline(threshold,color=INK,lw=.7,ls='--')
    e.set_xlim(-.5,1.03);e.set_xticks([-.5,0,.5,1]);e.set_ylabel('Samples');e.set_xlabel('Reference methylome correlation (r)')
    e.text(.04,.96,f'n = {len(r):,}\n{screen.dropna(subset=["r_reference"]).gse_id.nunique():,} studies',
           transform=e.transAxes,va='top',fontsize=7)
    e.text(.04,.69,f'{int((r<.1).sum()):,} below 0.10',transform=e.transAxes,fontsize=7,color=HIGHER)
    e.text(.1,.97,'0.10',transform=e.get_xaxis_transform(),ha='center',va='top',fontsize=6.5)
    e.text(.5,.97,'0.50',transform=e.get_xaxis_transform(),ha='center',va='top',fontsize=6.5)
    e.grid(axis='y',lw=.4,ls=':',alpha=.4)
    facts={'a_classification':rows_a,'b_defects':rows_b,'c_published_fates':rows_c,
           'd_index_study':{'gse':'GSE73950','n':case_n,'control_median_r':rmed,'control_median_age':amed},
           'e_screen':{'samples':len(r),'studies':int(screen.dropna(subset=['r_reference']).gse_id.nunique()),
                       'below_0.10':int((r<.1).sum()),'thresholds':[.1,.5], 'bins':bins.tolist()}}
    sources=[FINAL/f'repairs-{key}/manifest.json' for key,_ in original.MODALITIES]
    sources += [census_file, FINAL/'published-status/summary.json',original.SCREEN,original.INDEX_STUDY,metadata_file]
    out=save(fig,'ed21',dict(zip('abcde',[a,b,c,d,e])),13,sources=sources,facts=facts,script=__file__)
    pd.DataFrame(rows_a).to_csv(out/'panel_a_classification.csv',index=False)
    pd.DataFrame(rows_b).to_csv(out/'panel_b_defects.csv',index=False)
    pd.DataFrame(rows_c).to_csv(out/'panel_c_published_fates.csv',index=False)
    pd.DataFrame({'gsm':study.index,'arm':arm.values,'r_reference':study.r_reference.values,
                  'historical_Horvathv1':historical.values}).to_csv(out/'panel_d_sample_values.csv',index=False)
    pd.DataFrame({'bin_left':bins[:-1],'bin_right':bins[1:],'n':np.histogram(r,bins)[0],
                  'below_0.10_n':np.histogram(r[r<.1],bins)[0]}).to_csv(out/'panel_e_histogram.csv',index=False)
    readme(out,f'''# Atlas verification

ED21 retains all five October panels. July ED7 supplies the stacked-bar colours
and outlined marks; July Figure 3 supplies modality symbols, and July ED17
supplies the histogram treatment. This is a new figure, so no July data panel
is replaced. The comparison uses July page 13 at the same physical scale.

## Panel map

| October panel | Restored panel | Change |
| --- | --- | --- |
| a | a | Retain classification proportions; add counts inside wide segments |
| b | b | Retain defect percentages; circle for mouse RNA-seq, square for human RNA-seq and triangle for DNAm |
| c | c | Retain published-item fates with category colours |
| d | d | Retain every index-study point and overlapping-control annotation |
| e | e | Retain histogram bins and integrity thresholds with blue and orange fills |

## Data and numbers

Absolute source paths, modification times and SHA-256 hashes are in
`provenance.json`. `numbers.json` records every numerical annotation and every
plotted count, percentage and threshold. Axis ticks indicate the stated units.
The adjacent CSV files hold panel-level plotted values.

Panel a reads `repairs-*/manifest.json`, field `final_classes`. Pending
subclasses are combined as withheld. Counts and totals are:

{table_md(pd.DataFrame(rows_a).pivot(index='group',columns='class',values='n').reset_index())}

Panel b reads `error-census/census.json`, field
`modalities[modality].defects[defect].fraction`, multiplied by 100. These are
historical defects, which can coexist on one row; the percentages are not
mutually exclusive. Panel c reads `published-status/summary.json`; pending
unresolved outcomes and missing exact matches retain the October definitions.

{table_md(pd.DataFrame(rows_c).pivot(index='group',columns='class',values='n').reset_index())}

Panel d reads `input-integrity-dnam-001/samples.csv`, historical predictions
from `replay-dnam-001/inputs/GSE73950.csv`, and sample arms from
`GEO_hs_dnam/metadata/GSE73950.csv`. Endometriosis n = {case_n['Endometriosis']};
control n = {case_n['Control']}. The control median r is {rmed:.9g}; median
historical age is {amed:.9g} years. Printed values round to three and one
decimal places, respectively. These historical predictions illustrate an
input defect and are not used for an inferential disease claim.

Panel e reads `input-integrity-dnam-001/samples.csv`. The finite correlations
comprise {len(r):,} samples from {screen.dropna(subset=['r_reference']).gse_id.nunique():,}
studies; {int((r<.1).sum()):,} are below 0.10. Dashed lines are the pre-existing
0.10 detached-input and 0.50 atypical-input thresholds. Bins follow the October
builder exactly: 59 edges from −0.45 to 1.00.

## Draft legend

**a,** Final classification of atlas rows by modality. Bar segments show the
percentage in each class; labels within sufficiently wide segments give row
counts, and labels at right give all rows in the modality. Blue denotes
verified, pale blue provisional, green descriptive, purple duplicate, yellow
withheld and orange non-repairable records. **b,** Historical defects as a
percentage of each modality's atlas rows. Circles indicate mouse RNA-seq,
squares human RNA-seq and triangles human DNAm. Defects are not mutually
exclusive. **c,** Outcomes of re-evaluating published Table S14, Figure 3 and
Tables 1 and 2 items against their claimed direction. Numbers at right give
the numbers of published items. **d,** Historical Horvath predictions against
reference-methylome correlation for GSE73950, with {case_n['Endometriosis']}
endometriosis and {case_n['Control']} control samples. Orange indicates
endometriosis and blue controls; the control predictions overlap. **e,** Distribution of
reference-methylome correlations across {len(r):,} screened samples from
{screen.dropna(subset=['r_reference']).gse_id.nunique():,} studies. Orange marks
samples below 0.10. The vertical scale is logarithmic; dashed lines indicate
correlations of 0.10 and 0.50. These panels are descriptive counts and input
diagnostics; no new hypothesis tests are performed.

## Checks

The build runs the nature-figure defect checks for text size, overlap,
clipping and colour scales, with layout rules disabled to preserve July
styling. `checks.json` records findings. PDF dimensions are 183 × 167.64 mm.
PDF vector paths, embedded fonts and absence of raster plot images are
checked. The 200-dpi render and the July comparison require visual inspection.
No plotted value is simulated. The original builders and July generators
remain unchanged.
''')


if __name__=='__main__':
    main()
