"""Restore five human examples with July box-and-point styling."""
from ed21_24_restore_common_2026_11_05 import *
from matplotlib.ticker import MaxNLocator
import build_ed24_human_examples_2026_10_18 as original
from audit_system.row_repairs_2026_09_28 import individuals_by_comparison


def main():
    style()
    ranking_file=FINAL/'comparison-ranking/ranked_comparisons_all.csv'
    ranking=pd.read_csv(ranking_file)
    audit=paths.VERIFIER_ROOT/'individuals-audit-001'
    individuals=individuals_by_comparison(audit)
    fig=plt.figure(figsize=(WIDTH,5.95))
    positions=[[.095,.59,.22,.285],[.415,.59,.22,.285],[.735,.59,.22,.285],
               [.13,.115,.34,.285],[.625,.115,.33,.285]]
    axd={letter:axes(fig,box,letter=letter,label_x=box[0]-.075)
         for letter,box in zip('abcde',positions)}
    rng=np.random.default_rng(20261105)
    facts=[];values_rows=[];sources=[ranking_file]
    plotted=[]
    for letter,spec in zip('abcde',original.EXAMPLES):
        modality,folder,gse,comparison,title,tissue=spec
        case,control,row,ranked,method=original.tested_values(FINAL,individuals,modality,folder,gse,comparison,ranking)
        if row['final_class'] not in ['VERIFIED','PROVISIONAL']:
            raise ValueError(f'{gse}: ineligible class {row["final_class"]}')
        ax=axd[letter]
        marker=MODS[modality][1]
        # July sample-example generator's quartile box, median, caps and points.
        bp=ax.boxplot([control,case],positions=[0,1],widths=.48,patch_artist=True,
                      showfliers=False,medianprops={'color':INK,'lw':1.1},
                      whiskerprops={'color':INK,'lw':.7},capprops={'color':INK,'lw':.7},
                      boxprops={'edgecolor':INK,'lw':.7})
        for patch,color in zip(bp['boxes'],['#cccccc',HIGHER]):
            patch.set_facecolor(color);patch.set_alpha(.42)
        for x,data,color,group in [(0,control,'#666666','control'),(1,case,HIGHER,'case')]:
            offset=rng.uniform(-.13,.13,len(data))
            ax.scatter(x+offset,data,s=20,marker=marker,c=color,ec=INK,lw=.35,
                       alpha=.85,zorder=3,clip_on=False)
            values_rows.extend({'panel':letter,'gse':gse,'group':group,'unit_index':i,
                                'clock_value':float(value),'display_x':float(x+jitter),
                                'method':method or 'WELCH'}
                               for i,(value,jitter) in enumerate(zip(data,offset)))
        ax.set_xlim(-.55,1.55)
        low=min(case.min(),control.min());high=max(case.max(),control.max());span=high-low
        ax.set_ylim(low-.09*span,high+.42*span)
        xlabel0=original.CONTROL_LABEL.get(gse,'Control')
        xlabel1=original.CASE_LABEL[gse]
        ax.set_xticks([0,1],[f'{xlabel0}\nn = {len(control)}',f'{xlabel1}\nn = {len(case)}'])
        ylabel=original.CLOCK_LABEL.get(ranked.clock,'Clock output')
        if method=='AGE_ADJUSTED':ylabel+='\nadjusted for age'
        ax.set_ylabel(ylabel)
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_title(f'{title}\n{tissue}',loc='left',pad=19)
        ax.text(0,1.015,gse,transform=ax.transAxes,ha='left',va='bottom',fontsize=7)
        # The bracket and statistics use the verified row's test rather than a new t test.
        bracket=high+.105*span
        ax.plot([0,0,1,1],[bracket-.025*span,bracket,bracket,bracket-.025*span],color=INK,lw=.75)
        ax.text(.5,bracket+.025*span,f'g = {ranked.hedges_g:.2f}\nFDR = {sci(row["recomputed"]["fdr_family"])}',
                va='bottom',ha='center',fontsize=6.5,linespacing=1.25)
        grid(ax,axis='y')
        unit='individuals' if method in ['UNIT_MEANS','AGE_ADJUSTED'] else 'samples'
        facts.append({'panel':letter,'gse':gse,'comparison':comparison,'tissue':tissue,
                      'clock':ranked.clock,'modality':ranked.modality,'method':method or 'WELCH',
                      'n_case':len(case),'n_control':len(control),'unit':unit,
                      'hedges_g':float(ranked.hedges_g),'fdr':float(row['recomputed']['fdr_family']),
                      'mean_diff':float(row['recomputed']['mean_diff']),
                      'plotted_mean_diff':float(case.mean()-control.mean()),'final_class':row['final_class'],
                      'design':row.get('design')})
        sources += [FINAL/f'public-export/studies/{modality}/{gse}.json',
                    paths.resolve(row['predictions']['path'])]
        sources += original.design.load_series(paths.GEO,folder,gse,True)['files']
        plotted.append((ax,case,control,method,ranked.clock))
    # Match the unadjusted RNA clock scale across panels, as in October.
    unadjusted=[t for t in plotted if t[4]==plotted[0][4] and t[3]!='AGE_ADJUSTED']
    low=min(min(t[1].min(),t[2].min()) for t in unadjusted)
    high=max(max(t[1].max(),t[2].max()) for t in unadjusted)
    span=high-low
    for ax,_,_,_,_ in unadjusted:ax.set_ylim(low-.09*span,high+.42*span)
    # Include the specific individuals-audit files used to recover units.
    sources += [audit/'comparisons.jsonl']
    out=save(fig,'ed24',axd,3,sources=sources,
             facts={'examples':facts,'jitter_seed':20261105,'jitter_range':[-.13,.13]},script=__file__)
    pd.DataFrame(values_rows).to_csv(out/'plotted_sample_values.csv',index=False)
    df=pd.DataFrame(facts)
    df.drop(columns='design').to_csv(out/'example_statistics.csv',index=False)
    readme(out,f'''# Human intervention examples

ED24 retains the five October human examples and their verified observations.
The box-and-point drawing follows the copied July sample-example generator,
with the Figure 3 orange direction colour and square/triangle modality key.
The five plots occupy two filled rows. July Figure 3 on page 3 is the
same-scale colour and density comparison; no July scientific panel is replaced.

## Panel map

| October panel | Restored panel | Content |
| --- | --- | --- |
| a | a | VEXAS syndrome in whole blood |
| b | b | Critical limb ischaemia in gastrocnemius |
| c | c | Ataxia-telangiectasia in cerebellar vermis |
| d | d | Gonadotropin stimulation in granulosa cells |
| e | e | Varicose veins in vein wall |

Horizontal box-and-point strips are redrawn vertically using the July sample
example generator's boxes, medians, caps and overlaid observations. Every
observation, arm, clock, effect and FDR remains the October panel's value.

## Data and numbers

`provenance.json` lists absolute file paths, modification times and SHA-256
hashes. `plotted_sample_values.csv` contains every displayed clock value and
the cosmetic point position. `example_statistics.csv` and `numbers.json`
hold every printed n, Hedges' g and FDR before rounding. Axis ticks state
the clock-output scale and are not measurements.

Each panel reads the named per-study JSON under
`final-2026-10-10/public-export/studies/`, uses
`comparison-ranking/ranked_comparisons_all.csv` to select the same clock as
the October builder, and loads the record's prediction file. Metadata under
`GEO_*/metadata/` and the individuals audit recover the tested units. The
October `tested_values` function is reused unchanged. Repeated samples are
collapsed to individual means where the verified design requires it; the
age-adjusted panel uses the recorded age slope and the comparison mean age.

{table_md(df[['panel','gse','n_case','n_control','unit','hedges_g','fdr','mean_diff','method','final_class']])}

For every panel the plotted case mean minus control mean must reproduce the
recorded mean difference within 10⁻⁶. All five assertions pass during loading.
Effects and FDR labels come directly from the verified ranking and row;
no replacement tests are computed. The point jitter uses seed 20261105 and
uniform offsets within ±0.13 along the categorical axis. It changes only
point placement, not the measured clock values. The shared unadjusted RNA
clock limits retain cross-panel comparability; the adjusted RNA output and
DNAm age have their own scales.

## Draft legend

**a–e,** Clock outputs in five human comparisons: VEXAS syndrome in whole
blood (a), critical limb ischaemia in gastrocnemius (b), ataxia-telangiectasia
in cerebellar vermis (c), gonadotropin stimulation in granulosa cells (d),
and varicose veins in vein wall (e). Points show the tested observations,
with squares for RNA-seq and triangles for DNAm. Grey denotes controls and
orange the case or stimulated arm. Boxes span the 25th to 75th percentiles,
centre lines show medians and whiskers extend to observations within 1.5
interquartile ranges. All observations are overlaid, including those outside
the whiskers. The plotted sample sizes are {facts[0]['n_case']} cases and
{facts[0]['n_control']} controls in a, {facts[1]['n_case']} and
{facts[1]['n_control']} in b, {facts[2]['n_case']} and {facts[2]['n_control']}
in c, {facts[3]['n_case']} stimulated and {facts[3]['n_control']} natural-cycle
samples in d, and {facts[4]['n_case']} and {facts[4]['n_control']} in e.
Panels a and c use individual means; c additionally adjusts clock output
to the comparison mean chronological age using the fitted age coefficient.
Panels b, d and e use samples. Labels give case-minus-control Hedges' g and
BH-FDR within modality from the verified analysis. The corresponding source
tests are Welch two-sample tests on individual means for a and on samples
for b and d, the arm coefficient of a clock-on-arm-and-age linear model
for c, and a paired t-test on the ten matched pairs for e. The Horvath clock
output in e is in years.

## Checks

The build enforces VERIFIED or PROVISIONAL eligibility and exact recovery of
each recorded mean difference. `checks.json` records the text-floor,
overlap, canvas-clipping and vector checks. PDF width is 183 mm and height
151.13 mm. Render review uses 200 dpi and a same-scale July comparison.
Original generators, source records, manuscript files and the combined
figure package remain unchanged.
''')


if __name__=='__main__':
    main()
