#!/usr/bin/env python3
"""Package vector panels, same-scale comparisons and source ledgers."""
import os
from pathlib import Path
import ast
import collections
import difflib
import hashlib
import importlib.util
import json
import re
import shutil
import sys

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('builder', HERE/'build_july_exact_2026_11_06.py')
b = importlib.util.module_from_spec(spec); spec.loader.exec_module(b)
ROOT, OUT, WORK, ORIGINAL = b.ROOT, b.OUT, b.WORK, b.ORIGINAL
PREVIOUS = b.PREVIOUS
import fitz
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from scipy import stats

PANELS = {
    3: [('a','volcano_combined_modalities'), ('b_e','top_hits_barplot_combined_modalities_v2')],
    4: [('a','create_radar_plot_combined'), ('b','create_ridge_plot_combined'),
        ('c','analysis_2_genetic_types_barplot'), ('d','analysis_2_genetic_types_density'),
        ('e','analysis_3_pathway_enrichment'), ('f_clinical','analysis_4_drug_clinical_barplot'),
        ('f_approved','analysis_4_fda_antiaging_drugs'), ('g','geneage_waffle_nature'),
        ('h','drugage_waffle_nature')],
}
SOURCES = {
    'volcano_combined_modalities': 'Verified combined table; Cohen d and supplied BH-FDR; statistical annotation ranking. Axis ticks follow unchanged July autoscaling. Score-reference ticks 70/90/110/130 are inherited style constants, not verified observations.',
    'top_hits_barplot_combined_modalities_v2': 'comparison-ranking/figure3_category_top15.csv, 15 per category; rank_bound=absolute Hedges g confidence bound nearest zero. Numbers in intervention labels are source names; axis ticks use July autoscaling.',
    'create_radar_plot_combined': 'Verified combined table; category row counts and FDR<0.05 counts; 35k radial cap and 5k ticks are July constants.',
    'create_ridge_plot_combined': 'Verified combined table, FDR<0.05; category row counts; scipy Wilcoxon against zero; ticks -20/-10/0/10/20 are July constants; view limit +/-25.',
    'analysis_2_genetic_types_barplot': 'Verified mouse RNA rows with Genetic and FDR<0.05; original keyword classification; analysis_2_genetic_types_fold_change.csv gives counts, fold enrichment and exact two-sided binomial P. Ticks use July autoscaling.',
    'analysis_2_genetic_types_density': 'Same significant mouse classified Genetic rows as c; legend counts from those rows; KDE and tick positions use unchanged July code.',
    'analysis_3_pathway_enrichment': 'Verified mouse Genetic rows; original mouse symbol reference and extract_genes; cached KEGG_2019_Mouse; result/panel_e_mouse_all_ora_results.csv supplies Fisher P/BH-FDR; July threshold 0.05 and axis autoscaling.',
    'analysis_4_drug_clinical_barplot': 'Previous panel_f_drug_matches.csv filtered to mouse RNA, FDR<0.05, matched OpenTargets; result/analysis_4_drug_clinical_status.csv gives direction counts; two-sided binomial P versus 0.5; ticks derived from max direction counts in July code.',
    'analysis_4_fda_antiaging_drugs': 'Same matched significant mouse Drug rows; FDA_Approved and d<0; exact/synonym/tradename; strongest d per standardized name; first 15. Names and mechanism text are the original July mappings; axis ticks from July autoscaling.',
    'geneage_waffle_nature': 'Verified lifespan_concordance.csv, mouse GenAge; one target-manipulation unit per row; sum(concordant), size and one-sided binomtest versus 0.5. Unmatched=0 reflects matched-only input, not measured coverage.',
    'drugage_waffle_nature': 'Verified lifespan_concordance.csv, mouse DrugAge; sum(concordant), size and one-sided binomtest versus 0.5. Unmatched=0 reflects matched-only input, not measured coverage.',
}


def raster(pdf, dpi=200):
    with fitz.open(pdf) as d:
        p = d[0].get_pixmap(dpi=dpi, alpha=False)
        return Image.frombytes('RGB', (p.width,p.height), p.samples)


def pair(july, new, dest):
    # Both panels are rasterised at exactly 200dpi, then padded, never resized.
    old_im, new_im = raster(july), raster(new)
    width = max(old_im.width,new_im.width); height=max(old_im.height,new_im.height)
    canvas=Image.new('RGB',(width*2+24,height+36),'white')
    draw=ImageDraw.Draw(canvas);draw.text((10,10),'July',fill='black');draw.text((width+34,10),'Verified',fill='black')
    canvas.paste(old_im,(0,36));canvas.paste(new_im,(width+24,36))
    canvas.save(dest, dpi=(200,200))


def crop_pdf(source, rect, dest):
    with fitz.open(source) as src:
        doc=fitz.open();p=doc.new_page(width=rect.width,height=rect.height)
        p.show_pdf_page(p.rect,src,0,clip=rect,keep_proportion=True)
        doc.save(dest,garbage=4,deflate=True);doc.close()


def text_lines(pdf):
    with fitz.open(pdf) as d:
        return [line.strip() for line in d[0].get_text().splitlines() if line.strip()]


def diff_audit():
    rows=[]
    for name in b.GENERATORS + ['assemble_figure3_vector.py','assemble_figure4_vector.py']:
        src=ORIGINAL/name if name in b.GENERATORS else ROOT/'paper/legacy_reference/figure_rev1'/name
        old=src.read_text();new=(HERE/name).read_text()
        changes=list(difflib.unified_diff(old.splitlines(True),new.splitlines(True),fromfile=str(src),tofile=str(HERE/name)))
        assert ''.join(changes)==(HERE/Path(name).with_suffix('.diff')).read_text()
        added=sum(x.startswith('+') and not x.startswith('+++') for x in changes)
        removed=sum(x.startswith('-') and not x.startswith('---') for x in changes)
        # Apart from explicitly allowed membership blocks, erase only allowed
        # text/path substitutions and demand identical Python ASTs.
        ta,tb=ast.parse(old),ast.parse(new)
        if name=='create_combined_modality_volcano.py':
            def remove_selection(tree):
                fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
                start=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='ann_pool' for t in n.targets))
                end=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='_all' for t in n.targets))
                del fn.body[start:end]
            remove_selection(ta);remove_selection(tb)
        elif name=='create_combined_modality_top_hits_v2.py':
            def remove_selection(tree):
                fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
                start=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in ['thresholds','qual'] for t in n.targets))
                end=next(i for i,n in enumerate(fn.body) if isinstance(n,ast.Assign) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Attribute) and n.value.func.attr=='subplots')
                del fn.body[start:end]
            remove_selection(ta);remove_selection(tb)
        class Normalise(ast.NodeTransformer):
            def visit_Assign(self,node):
                if name.startswith('assemble_figure') and any(isinstance(t,ast.Name) and t.id in ['ROOT','FIGSRC'] for t in node.targets):
                    node.value=ast.Constant(value='<path>')
                return self.generic_visit(node)
            def visit_Constant(self,node):
                if isinstance(node.value,str):
                    # Individual string changes are additionally checked against
                    # the human-readable unified diff and documented whitelist.
                    return ast.copy_location(ast.Constant(value='<string>'),node)
                return node
        equal=ast.dump(Normalise().visit(ta),include_attributes=False)==ast.dump(Normalise().visit(tb),include_attributes=False)
        assert equal, f'Non-text/path/selection AST changed: {name}'
        rows.append({'script':name,'added_lines':added,'removed_lines':removed,
            'unchanged_nontext_plotting_ast':equal,'original_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
            'copy_sha256':hashlib.sha256((HERE/name).read_bytes()).hexdigest()})
    pd.DataFrame(rows).to_csv(OUT/'script_diff_summary.csv',index=False)
    return pd.DataFrame(rows)


def package():
    b.assemble()
    references=[]
    for _, panels in PANELS.items():
        for _,name in panels:
            path=ORIGINAL/'figures'/(name+'.pdf')
            with fitz.open(path) as dd:
                references.append({'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                    'mtime_utc':b.datetime.fromtimestamp(path.stat().st_mtime,b.timezone.utc).isoformat(),
                    'page_size_pt':list(dd[0].rect)[2:],'metadata':dd.metadata})
    b.write_json(OUT/'july_reference/reference_provenance.json',references)
    allchecks=[]
    for num,panels in PANELS.items():
        folder=OUT/f'figure{num}'
        pdf=folder/f'figure{num}.pdf'
        raster(pdf).save(folder/f'figure{num}_review.png',dpi=(200,200))
        pair(OUT/f'july_reference/figure{num}_july.pdf',pdf,folder/f'figure{num}_vs_july.png')
        ledger=[];changes=[]
        for letter,name in panels:
            source=WORK/'figures'/(name+'.pdf');dest=folder/f'panel_{letter}.pdf'
            shutil.copy2(source,dest)
            pair(ORIGINAL/'figures'/(name+'.pdf'),dest,folder/f'panel_{letter}_vs_july.png')
            raster(dest).save(folder/f'panel_{letter}_review.png',dpi=(200,200))
            old_lines=text_lines(ORIGINAL/'figures'/(name+'.pdf'));new_lines=text_lines(dest)
            removed=list((collections.Counter(old_lines)-collections.Counter(new_lines)).elements())
            added=list((collections.Counter(new_lines)-collections.Counter(old_lines)).elements())
            changes.append({'panel':letter,'source':name,'july_only_text':removed,'new_only_text':added})
            (folder/f'panel_{letter}_text.txt').write_text('\n'.join(new_lines)+'\n')
            for line in new_lines:
                if re.search(r'\d',line):
                    ledger.append({'panel':letter,'printed_text':line,'source_or_rule':SOURCES[name]})
            with fitz.open(dest) as d:
                page=d[0]
                allchecks.append({'figure':num,'panel':letter,'page_size_pt':list(page.rect)[2:],
                    'font_count':len(page.get_fonts()),'image_count':len(page.get_images()),
                    'producer':d.metadata.get('producer'),'missing_glyph': '\ufffd' in page.get_text() or '\u25a1' in page.get_text()})
        if num==3:
            # Isolate b-e as vector clips of the original 2x2 grid. The full
            # panel_b_e.pdf retains the shared legend and statistical heading.
            for label,x,y in [('b',0,0),('c',1,0),('d',0,1),('e',1,1)]:
                paths=[]
                for old,source in [(True,ORIGINAL/'figures/top_hits_barplot_combined_modalities_v2.pdf'),
                                   (False,folder/'panel_b_e.pdf')]:
                    with fitz.open(source) as dd:w,h=dd[0].rect.width,dd[0].rect.height
                    top=62.;mid=(top+h)/2
                    rect=fitz.Rect(x*w/2,top if y==0 else mid,(x+1)*w/2,mid if y==0 else h)
                    dest=(OUT/'july_reference'/f'figure3_panel_{label}.pdf') if old else folder/f'panel_{label}.pdf'
                    crop_pdf(source,rect,dest);paths.append(dest)
                pair(*paths,folder/f'panel_{label}_vs_july.png')
        pd.DataFrame(ledger).to_csv(folder/'printed_numbers.csv',index=False)
        b.write_json(folder/'text_changes_vs_july.json',changes)
    b.write_json(OUT/'pdf_checks.json',allchecks)
    assert all(c['font_count']>0 and not c['missing_glyph'] for c in allchecks)


def metrics_and_readmes(diffs):
    table=WORK/'figures/combined_3modality_plot_data.csv'
    df=pd.read_csv(table,low_memory=False);sig=df[df.fdr<.05];mouse=df[df.modality=='mouse RNA-seq']
    rank=pd.read_csv(WORK/'result/figure3_category_top15_july.csv')
    gen=pd.read_csv(WORK/'result/analysis_2_genetic_types_fold_change.csv')
    clinical=pd.read_csv(WORK/'result/analysis_4_drug_clinical_status.csv')
    life=pd.read_csv(WORK/'data/lifespan_concordance_source.csv')
    paths=pd.read_csv(WORK/'result/analysis_3_pathway_enrichment_ora.csv')
    cancer='cancer|carcinoma|tumor|tumour|neoplasm|oncogene|metastasis|malignancy|leukemia|lymphoma|melanoma|glioma|sarcoma|myeloma'
    plotted=paths[(paths.FDR<.05)&~paths.Pathway.str.contains(cancer,case=False,na=False)].head(30)
    plotted.to_csv(WORK/'result/panel_e_plotted_pathways.csv',index=False)
    baseline=int(gen.N_Anti_Aging.sum())/int(gen.N_Interventions.sum())
    catcounts=pd.DataFrame({'all_rows':df.condition_category.value_counts(),
        'significant_rows':sig.condition_category.value_counts()}).sort_index()
    catcounts.to_csv(WORK/'result/panel_a_category_counts.csv')
    ridges=[]
    for cat,g in sig.groupby('condition_category'):
        d=g.effect_size;stat,p=stats.wilcoxon(d[d!=0]);ridges.append({'category':cat,'n':len(d),'median_d':float(d.median()),'wilcoxon_p':p})
    ridge=pd.DataFrame(ridges);ridge.to_csv(WORK/'result/panel_b_ridge_statistics.csv',index=False)
    life_numbers={s:{'units':len(g),'genes_or_drugs':g.canonical_target.nunique(),'agree':int(g.concordant.sum()),
        'disagree':int((~g.concordant).sum()),'p_greater':stats.binomtest(int(g.concordant.sum()),len(g),.5,alternative='greater').pvalue} for s,g in life.groupby('source')}
    figuresize={}
    for n in [3,4]:
        with fitz.open(OUT/f'figure{n}/figure{n}.pdf') as d:figuresize[n]=list(d[0].rect)[2:]
    numbers={'atlas_rows':len(df),'unique_comparisons':len(df.drop_duplicates(['modality','gse_id','comparison_name','condition_name'])),
        'significant_rows':len(sig),'significant_percent':100*len(sig)/len(df),'modality_counts':df.modality.value_counts().to_dict(),
        'category_counts':catcounts.reset_index().to_dict('records'),'d_range':[df.effect_size.min(),df.effect_size.max()],
        'max_neg_log10_fdr':float(-np.log10(df.fdr.min())),
        'ridge_outside_fixed_xlim':int((sig.effect_size.abs()>25).sum()),
        'radar_above_35000':int((catcounts.all_rows>35000).sum()),
        'genetic_baseline':baseline,'genetic':gen.to_dict('records'),
        'clinical':clinical.to_dict('records'),'lifespan':life_numbers,
        'pathways':plotted.to_dict('records'),'page_size_pt':figuresize}
    b.write_json(OUT/'numbers.json',numbers)
    for n in [3,4]:
        folder=OUT/f'figure{n}'
        rows=[]
        for letter,name in PANELS[n]:
            generator=next((x for x in b.GENERATORS if x.replace('.py','')==name),None)
            if name=='volcano_combined_modalities':generator='create_combined_modality_volcano.py'
            elif name=='top_hits_barplot_combined_modalities_v2':generator='create_combined_modality_top_hits_v2.py'
            if name.startswith('analysis_2'):generator='analysis_2_genetic_intervention_types.py'
            elif name.startswith('analysis_3'):generator='analysis_3_pathway_enrichment.py'
            elif name.startswith('analysis_4'):generator='analysis_4_fda_approved_drugs.py'
            elif name=='geneage_waffle_nature':generator='create_geneage_waffle_plot_nature.py'
            elif name=='drugage_waffle_nature':generator='create_drugage_waffle_plot_nature.py'
            rows.append(f'| {letter.replace("b_e","b–e")} | `{name}.pdf` | [{generator}]({HERE/generator}) |')
        scripts=diffs[diffs.script.isin(['create_combined_modality_volcano.py','create_combined_modality_top_hits_v2.py','assemble_figure3_vector.py'])] if n==3 else diffs[~diffs.script.isin(['create_combined_modality_volcano.py','create_combined_modality_top_hits_v2.py','assemble_figure3_vector.py'])]
        diff_lines='\n'.join(f'| [{r.script}]({HERE/r.script}) | +{r.added_lines} / −{r.removed_lines} | [{Path(r.script).with_suffix(".diff")}]({HERE/Path(r.script).with_suffix(".diff")}) |' for r in scripts.itertuples())
        intro=f'''# Figure {n} reproduction

The vector PDF uses the July generators and assembly geometry with verified inputs. `figure{n}_review.png` is rendered at 200 dpi. `figure{n}_vs_july.png` places the July assembly on the left and the verified assembly on the right at the same physical scale. Comparison images are padded without resizing.

## Panel map

| Panel | July panel PDF | Copied generator |
|---|---|---|
'''+ '\n'.join(rows)+'''

The per-panel PDFs preserve the generated vectors. Figure 3 b–e also have individual vector clips; `panel_b_e.pdf` retains their shared heading and legend.

## Data sources

'''
        if n==4:
            intro=intro.replace(' Figure 3 b–e also have individual vector clips; `panel_b_e.pdf` retains their shared heading and legend.','')
        source_text='\n\n'.join(f'**{letter.replace("b_e","b–e")}**. {SOURCES[name]}' for letter,name in PANELS[n])
        common=f'''

The original generators are under `{ORIGINAL}`. The assembly originals are under `{ROOT / 'paper/legacy_reference/figure_rev1'}`. July reference PDFs were assembled from `{ORIGINAL / 'figures'}` with the original functions. The reference is the supplied panel PDFs, not a rerender of the previous attempt.

The empirical plotting source is `{PREVIOUS / 'data/combined_3modality_plot_data_verified_2026-11-05.csv'}`. Its row-provenance companion supplies clocks for exact joins. Source paths, UTC modification times and SHA-256 hashes are recorded in [input_provenance.json](../july_exact/data/input_provenance.json). [adapter_rules.json](../july_exact/data/adapter_rules.json) records every compatibility field. No measurements are simulated.

## Printed numbers

[printed_numbers.csv](printed_numbers.csv) lists every PDF text line containing a digit and its source or generation rule, including tick labels and digits in intervention names. [text_changes_vs_july.json](text_changes_vs_july.json) records all removed and added text against July. [numbers.json](../numbers.json) stores the unrounded results. The original ticks and July formatting determine rounding in the image.

'''
        if n==3:
            details=f'''The title now reports {len(df):,} atlas rows, replacing July's 101,872. These rows represent {numbers['unique_comparisons']:,} distinct modality/GEO/comparison/condition keys and retain separate clocks. Significant rows change from 11,664 (11.4%) to {len(sig):,} ({100*len(sig)/len(df):.1f}%). Modality counts are {len(df[df.modality.eq('mouse RNA-seq')]):,} / {len(df[df.modality.eq('human RNA-seq')]):,} / {len(df[df.modality.eq('human DNAm')]):,} for mouse RNA / human RNA / human DNAm. The inherited title calls these rows comparisons.

The 60 bars use the supplied 15-per-category statistical ranking. Their lengths are absolute Hedges' g CI bounds nearest zero. [Selected rows](../july_exact/result/top_hits_combined_modalities_v2.csv) include both CI endpoints, rank_bound, verified row index, Cohen's d and BH-FDR. No composite score is calculated. Volcano labels use the ranked-comparison file and July's per-category direction budgets, yielding 36 higher and 36 lower labels.

| Category | First ranked comparison | CI bound nearest zero |
|---|---|---:|
'''+ '\n'.join(f'| {r.condition_category} | {r.condition_name} | {r.rank_bound:.9g} |' for r in rank.groupby('condition_category',sort=False).head(1).itertuples())+'''

The ranking file's historical `effect_size` alias is not used as a measurement. Some aliases differ in magnitude or sign from the corrected values. Joins use the corrected Hedges' g with the sample-size correction and verified p-values to resolve repeated identifiers, then use the supplied verified table for plotted Cohen's d, direction and FDR.
'''
            limits=f'''## Reproduction limits

The original volcano chooses axis limits from its input and expands them by the original formula. The verified d range is {numbers['d_range'][0]:.6g} to {numbers['d_range'][1]:.6g}; maximum −log10(FDR) is {numbers['max_neg_log10_fdr']:.6g}. July's annotation strip remains at y=1.50103–10. Labels therefore crowd near the bottom of the larger significance range. Changing the strip or freezing its autoscaled limits would change plotting code outside the permitted edits.

`final_score` is missing throughout the verified input. July's unchanged `normalize_size` uses its default 50, so every significant point has size `(80/170)^3*120+4`. The inherited Final score legend still shows 70/90/110/130 reference symbols. It is an obsolete score scale and does not describe verified measurements. It is retained under the instruction to preserve July's plotting code; the figure is unsuitable for publication with that legend unless its removal is authorized.

All {len(df):,} supplied inference-supporting rows receive the input compatibility literal `outlier_status=Normal`, so the July inclusion filter retains them. This is not a historical outlier classification. No old scores or exclusions are imported.

Figure 3 retains PAGE_W=1001 pt. Its height is computed by the unchanged July formula from each panel's tight PDF bounding box, so longer labels and altered ticks can change the assembled height. The new page is {figuresize[3][0]:.6f} × {figuresize[3][1]:.6f} pt. Nominal figsize, fonts, colours, marker shapes, label-placement logic, layout functions and limits logic are unchanged.
'''
        else:
            old_counts={'Disease':(16325,2585),'Drug':(22336,1907),'Environment':(9747,1108),'Genetic':(35323,3300),'Other':(18141,2764)}
            details='''| Radar category | July all / significant | Verified all / significant |
|---|---:|---:|
'''+ '\n'.join(f'| {cat} | {old_counts[cat][0]:,} / {old_counts[cat][1]:,} | {r.all_rows:,} / {r.significant_rows:,} |' for cat,r in catcounts.iterrows())+f'''

The radar reports {len(df):,} rows and {len(sig):,} significant rows ({100*len(sig)/len(df):.1f}%). Ridge category counts are the significant counts above; [ridge statistics](../july_exact/result/panel_b_ridge_statistics.csv) contain every median and Wilcoxon P. The significant set is pooled across the three modalities only in a–b.

Mouse genetic counts change from KO 775 / KD 85 / OE 140 / mutation 70 to KO 657 / KD 64 / OE 121 / mutation 55. The classified significant mouse baseline changes from 27.1% to {baseline:.1%}. Enrichment and binomial P use the original July calculation.

| Type | n | Lower / higher | Fold enrichment | Two-sided P |
|---|---:|---:|---:|---:|
'''+ '\n'.join(f'| {r.Type} | {r.N_Interventions} | {r.N_Anti_Aging} / {r.N_Pro_Aging} | {r.Fold_Change:.8g} | {r.P_value:.8g} |' for r in gen.itertuples())+'''

The mouse KEGG analysis uses 1,497 background genes, 403 significant genes and 237 tested nonzero-overlap terms. Seventeen terms pass size<100, overlap≥15 and FDR<0.05; twelve remain after July's cancer-name exclusion, replacing seven July bars. [Pathway values](../july_exact/result/panel_e_plotted_pathways.csv) give every bar's −log10(FDR), FDR, raw P, overlap and pathway size. Fisher contingency counts for every tested term are in [all results](../july_exact/result/panel_e_mouse_all_ora_results.csv). The cached library and local Fisher/BH calculation replace an online Enrichr call; the plotted function is byte-identical to July. This reproduces the stated custom-background method, not an archived Enrichr server response.

| Clinical stage | July lower / higher | Verified lower / higher | Verified n | Two-sided P |
|---|---:|---:|---:|---:|
'''
            old_phases={'FDA_Approved':(80,136),'Phase_3':(9,13),'Phase_2':(7,24),'Phase_1':(2,11),'Preclinical':(4,11)}
            for r in clinical.itertuples():
                old=old_phases[r.Clinical_Status];p=stats.binomtest(r.N_Anti_Aging,r.N_Drugs,.5).pvalue
                details+=f'| {r.Clinical_Status} | {old[0]} / {old[1]} | {r.N_Anti_Aging} / {r.N_Pro_Aging} | {r.N_Drugs} | {p:.9g} |\n'
            details+='''
These are atlas signatures, not unique drugs. The 161 FDA signatures include 75 distinct standardized names. The approved-drug inset uses the original strongest-effect and high-confidence matching selection, with 15 names. [Inset rows](../july_exact/result/panel_f_mouse_approved_lower_top15.csv) give every plotted effect, name, match type and source comparison.

| Lifespan reference | July matched agreement / disagreement | Verified agreement / disagreement | One-sided P |
|---|---:|---:|---:|
'''
            for name,old in [('GenAge','16 / 4'),('DrugAge','4 / 0')]:
                r=life_numbers[name];details+=f'| {name} | {old} | {r["agree"]} / {r["disagree"]} | {r["p_greater"]:.12g} |\n'
            details+='''
GenAge now has 62 target–manipulation units across 49 genes, and DrugAge has six units. The GenAge inputs append manipulation to names only when a gene has multiple rows, preserving distinct expected directions through July's name-based grouping. The July greater-than-50% binomial test is unchanged. Different manipulations of one gene need not be independent; the inherited test does not account for that dependence.
'''
            limits='''## Reproduction limits

The radar cap is unchanged at 35,000; no new category count exceeds it. The ridge remains at −25 to +25; 132 significant input rows lie outside the display and remain in the KDE and tests. These rows are not removed from the data.

The lifespan source contains matched units only. The old waffle totals 555/204, which included 535/200 unmatched entries, cannot be reconstructed from it. The new totals 62/6 and unmatched=0 describe the supplied matched-only input. The grey legend's zero is not a coverage estimate. No unmatched denominator is invented. Compatibility `trust_score=+inf` passes the unchanged >80 input filter and is explicitly a sentinel, not a measured trust score.

GenAge has more matched labels than July and fewer waffle rows. The original 100-column geometry and label positions are preserved, which crowds labels and brings them into the title area. No label thinning, font shrinking, axis expansion or new layout is applied. The fixed pathway panel also inherits July's long-name pressure on the plotting area. The clinical panel's longer corrected labels inherit its fixed placement.

Figure 4 keeps PAGE_W=595.68 and PAGE_H=706.68 pt, the nine original source rectangles and letters a–h. Panels c–h use mouse-only inputs. The FDA density PDF generated as an original function side effect stays in the isolated working folder and is not assembled.
'''
        ending=f'''
## Script diffs

| Script | Lines added / removed | Unified diff |
|---|---:|---|
{diff_lines}

The non-text plotting AST matches the July originals after removing the two authorized membership/ranking blocks. Numeric plotting constants are unchanged. Diff summaries and original/copy hashes are in [script_diff_summary.csv](../script_diff_summary.csv). Radar and pathway scripts have empty diffs because they are byte-identical. The assembly diffs change paths only.

## Verification

The renderer uses Matplotlib 3.10.1, matching the July PDF producer metadata, with aquarel 0.0.7 and the supplied Python environment. Other July package versions were not recorded. Per-source [figure records](../july_exact/logs/) retain nominal sizes, axes positions and limits, ticks, text, package versions and nature-figure preflight findings. Preflight style changes are not applied because the July reference is required.

The generic preflight does not pass these preserved figures. It flags inherited large type and layout conventions as well as visible data-dependent overlaps; the polar radar check also raises a missing `left` spine error. The recorded findings are review evidence, not a claim of publication readiness. Same-data pixel identity to July is not established: the current historical mouse master yields KO n=768 on replay, whereas the archived July panel prints n=775. Reference comparisons therefore use the archived PDFs. Their hashes and producer metadata are in [reference_provenance.json](../july_reference/reference_provenance.json).

All delivered panels have embedded fonts and selectable text. The volcano keeps July's rasterized scatter layer; its axes and labels remain vector. All other assembled panels remain vector. [PDF checks](../pdf_checks.json) record page sizes, fonts and image counts. No footers, extra panels or figure numbers were added.

## Reproduction command

```sh
source {ROOT.parent / 'env_local.sh'}
python {HERE / 'build_july_exact_2026_11_06.py'} prepare
python {HERE / 'build_july_exact_2026_11_06.py'} run
python {HERE / 'package_figures_2026_11_06.py'}
```
'''
        (folder/'README.md').write_text(intro+source_text+common+details+'\n'+limits+ending)


# Figure 3 colours: the palette chosen by FIG3_PALETTE (B when unset), as in the generators.
_pspec = importlib.util.spec_from_file_location('fig3_palette', HERE / 'fig3_palette_2026_11_08.py')
fig3_palette = importlib.util.module_from_spec(_pspec)
_pspec.loader.exec_module(fig3_palette)
PALETTE_NAME = fig3_palette.name()
_PAL = fig3_palette.active()
BAR_COLOR_RULE = (f"{_PAL['names']['lower']} {_PAL['lower']} for lower, "
                  f"{_PAL['names']['higher']} {_PAL['higher']} for higher.")
BAR_COLOR_SENTENCE = (f"{_PAL['names']['lower'][:1].upper()}{_PAL['names']['lower'][1:]} (`{_PAL['lower']}`) marks lower outputs and "
                      f"{_PAL['names']['higher']} (`{_PAL['higher']}`) higher outputs")
PALETTE_RULE = f'Colour roles of fig3_palette_2026_11_08.py, palette {PALETTE_NAME} (FIG3_PALETTE).'


def _map_text(pal):
    """The volcano colour map of the palette, for the README: a matplotlib name, or its nodes with the bar fills named."""
    spec = pal['diverging']
    if isinstance(spec, str):
        return f'the {spec} map'
    role = {pal['lower']: ' (lower fill)', pal['higher']: ' (higher fill)'}
    return 'a diverging map through ' + ', '.join(f"`{c}`{role.get(c, '')}" for _, c in spec)


_DENSITY = _PAL['density']
VOLCANO_COLOR_SENTENCE = (
    f"Rows with FDR < 0.05 are coloured by Cohen's d on {_map_text(_PAL)}, the others light grey, over "
    + (f"{_DENSITY} density shading" if isinstance(_DENSITY, str) else
       f"density shading of all rows from `{_DENSITY[0][1]}` to `{_DENSITY[-1][1]}`")
    + f" at alpha {_PAL['density_alpha']:g}; label colours give the intervention category ("
    + ', '.join(f"{k} `{v}`" for k, v in _PAL['category'].items() if k != 'Other') + ').')

# Figure 3g uses admitted comparisons in statistical order.
G_XLABEL = "|Hedges' g| (CI bound nearest zero)"
G_MAP_SELECTION = ('Admitted judgement rows in statistical order, the first of each intervention and GEO series; '
                   'the column gives the statistical rank of the comparison drawn')
G_RULE = ('Admitted comparisons from followup-judgment-001 in statistical order; the first comparison of each '
          'display_key and then each GEO series is drawn. Bar length is the absolute adjusted confidence bound '
          'nearest zero. Ouabain is bold and carries the dagger of panel a.')
ORDINAL = {1: 'first', 2: 'second', 3: 'third', 4: 'fourth', 5: 'fifth'}


def G_SOURCE_PATH():
    import os
    return Path(os.environ.get('CLOCKBASE_JUDGEMENT') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'followup-judgment-001') / 'judgement_age-lowering.csv'


def check_figure3_g(gpicks, judgements, keys):
    """Verify statistical order, display identity and adjusted confidence bounds."""
    import dedup_interventions_2026_11_06 as dedup
    import os
    assert judgements.judgement_admitted.all()
    meta = pd.read_csv(Path(os.environ.get('CLOCKBASE_META') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'intervention-meta-2026-10-17') / 'primary/meta_input_effects.csv')
    source = judgements.sort_values('statistical_rank', kind='stable')
    assert source.rank_bound.abs().is_monotonic_decreasing
    source = dedup.identity(source.drop(columns=['canonical_target'], errors='ignore'), meta)
    expected = source.drop_duplicates('display_key').drop_duplicates('gse_id')
    assert gpicks[keys].values.tolist() == expected[keys].values.tolist()
    for col in ['hedges_g', 'g_ci_low', 'g_ci_high', 'rank_bound', 'v3_p', 'v3_fdr', 'judgement_rank']:
        assert np.allclose(gpicks[col], expected[col], rtol=1e-12, atol=0), col
    assert np.allclose(gpicks.statistical_score, expected.rank_bound.abs(), rtol=0, atol=1e-12)
    return gpicks.label.tolist().index('Ouabain') + 1


def figure3_g_readme(gpicks):
    """Source comparisons and bounds of the displayed judgement entries."""
    lines = '\n'.join(
        f'| {int(r.judgement_rank)} | {r.display} | {r.gse_id} | {r.statistical_score:.6f} | '
        f'{r.hedges_g:.6f} | {int(r.v3_n_case)}/{int(r.v3_n_control)} | {r.v3_fdr:.9g} | '
        f'{r.rank_text} |' for r in gpicks.itertuples())
    return f'''## Panel g numbers

Admitted comparisons remain in statistical order. The panel draws the first
comparison of each intervention identity and then each GEO series. Bars give
the absolute adjusted confidence bound nearest zero, as in the statistical
panels. Values come from `{G_SOURCE_PATH()}`. Sample sizes are case/control
and FDR is `v3_fdr`. Full precision remains in the CSV and numbers record.

| Rank | Comparison | GEO | CI bound nearest zero | Hedges' g | n | FDR | Statistical rank |
|---:|---|---|---:|---:|---:|---:|---:|
{lines}'''


def package_figure3_followup():
    """Package only Figure 3."""
    import io
    import os
    import subprocess
    from datetime import datetime, timezone

    def subprocess_show(rel, commit):
        """A tracked CSV as committed in an earlier round."""
        text = subprocess.run(['git', '-C', str(ROOT), 'show', f'{commit}:{rel}'], capture_output=True, text=True,
                              check=True).stdout
        return pd.read_csv(io.StringIO(text), low_memory=False)

    folder = OUT / 'figure3'
    rank_path = Path(os.environ['CLOCKBASE_FINAL']) / 'comparison-ranking/ranked_comparisons_all.csv'
    judgement_path = Path(os.environ.get('CLOCKBASE_JUDGEMENT') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'followup-judgment-001') / 'judgement_age-lowering.csv'
    ranked = pd.read_csv(rank_path)
    judgements = pd.read_csv(judgement_path)
    replication_path = Path(os.environ.get('CLOCKBASE_META') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'intervention-meta-2026-10-17') / 'primary/ranked_replicated_age-lowering.csv'
    answers_path = Path(os.environ.get('CLOCKBASE_JUDGEMENT') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'followup-judgment-001') / 'ranked_with_judgement.csv'
    picks = pd.read_csv(WORK / 'result/top_hits_combined_modalities_v2.csv')
    fpicks = pd.read_csv(WORK / 'result/figure3_replication_picks.csv')
    gpicks = pd.read_csv(WORK / 'result/figure3_judgement_picks.csv')
    annotations = pd.read_csv(WORK / 'result/combined_volcano_annotations.csv')
    answers = pd.read_csv(answers_path, low_memory=False)
    keys = ['modality', 'gse_id', 'comparison_name', 'condition_name', 'clock']
    # b-e: per category, the five top-ranked with higher and the five with lower outputs, one entry per
    # intervention and GEO series (dedup_interventions_2026_11_06.py), recomputed here from the ranking.
    shown = pd.read_csv(WORK / 'result/volcano_ranked_comparisons_july.csv', low_memory=False)
    grid = pd.read_csv(WORK / 'result/figure3_condition_grid.csv')
    ranks = pd.read_csv(WORK / 'result/figure3_statistical_ranks.csv')
    geometry = json.loads((WORK / 'result/figure3_be_geometry.json').read_text())
    # b-e cells, recomputed here from the agreed answers: met, not met, or unknown (not stated, or NaN where the two
    # models did not agree); the reason text lists the unmet conditions, or reads 'all three met'.
    state = {'age_at_start': {'OLD': 'met', 'ADULT': 'not met', 'YOUNG': 'not met', 'GERMLINE_OR_PRENATAL': 'not met',
                              'NOT_STATED': 'unknown'},
             'background': {'HEALTHY': 'met', 'DISEASE_MODEL': 'not met', 'NOT_STATED': 'unknown'},
             'control': {'UNTREATED': 'met', 'ACTIVE_COMPARATOR': 'not met'}}
    words = {'GERMLINE_OR_PRENATAL': 'germline', 'ADULT': 'adult', 'YOUNG': 'young', 'DISEASE_MODEL': 'disease model',
             'ACTIVE_COMPARATOR': 'active control'}
    fields = ('age_at_start', 'background', 'control')
    state_cols = ('state_old_at_start', 'state_otherwise_healthy', 'state_untreated_control')
    # Validate selected rows against this run's ranking rather than an earlier release.
    assert len(grid) == int((picks.hedges_g < 0).sum())
    assert grid[['gse_id', 'comparison_name']].values.tolist() == \
        picks.loc[picks.hedges_g < 0, ['gse_id', 'comparison_name']].values.tolist()
    for cat in ['Genetic', 'Drug', 'Disease', 'Environment']:
        expected = shown[shown.display_best & shown.display_category.eq(cat)].copy()
        expected = expected.iloc[np.argsort(-expected.rank_bound.abs().to_numpy(), kind='stable')]
        expected = pd.concat([expected[expected.hedges_g > 0].head(5), expected[expected.hedges_g < 0].head(5)])
        expected = expected.iloc[np.argsort(-expected.rank_bound.abs().to_numpy(), kind='stable')]
        selected = picks[picks.category.eq(cat)]
        assert len(selected) == len(expected)
        assert selected[keys].fillna('').values.tolist() == expected[keys].fillna('').values.tolist()
        assert np.allclose(selected.hedges_g, expected.hedges_g)
        assert np.allclose(selected.composite, expected.rank_bound.abs())
        assert np.allclose(selected.g_ci_low, expected.g_ci_low)
        assert np.allclose(selected.g_ci_high, expected.g_ci_high)
        # every cell and reason restates the two models' agreed answers for that comparison
        for _, r in selected.iterrows():
            a = answers[(answers.gse_id == r.gse_id) & (answers.comparison_name == r.comparison_name)
                        & (answers.condition_name == r.condition_name)]
            assert len(a) == 1
            a = a.iloc[0]
            if r.hedges_g > 0:
                assert pd.isna(r.reason) and all(pd.isna(r[c]) for c in state_cols)
                continue
            states = ['unknown' if pd.isna(a[f]) else state[f][a[f]] for f in fields]
            reason = 'all three met' if states == ['met'] * 3 else \
                ' · '.join(words[a[f]] for f, s in zip(fields, states) if s == 'not met')
            g = grid[(grid.gse_id == r.gse_id) & (grid.comparison_name == r.comparison_name)].iloc[0]
            assert [r[c] for c in state_cols] == states == [g[c] for c in state_cols], (r.display, states)
            assert (r.reason if isinstance(r.reason, str) else '') == reason == (g.reason if isinstance(g.reason, str) else '')
            for f in fields:
                source = 'NaN (models did not agree)' if pd.isna(a[f]) else a[f]
                assert g['source_' + f] == source
            assert (states == ['met'] * 3) == bool(a.judgement_admitted) == bool(g.judgement_admitted)
    all_met = grid[grid.reason == 'all three met']
    assert all_met.judgement_admitted.all() and not grid.loc[grid.reason != 'all three met', 'judgement_admitted'].any()
    # f: every robust replicated intervention with lower outputs, in the order of the pooled bound
    replicated = pd.read_csv(replication_path)
    replicated = replicated.iloc[np.argsort(-replicated.rank_bound.abs().to_numpy(), kind='stable')]
    assert fpicks.canonical_target.tolist() == replicated.canonical_target.tolist() and len(fpicks) == len(replicated)
    for col in ['estimate', 'ci_low', 'ci_high', 'k', 'n_series', 'fdr']:
        assert np.allclose(fpicks[col], replicated[col])
    assert np.allclose(fpicks.composite, -replicated.ci_high)
    plx5622_replication_rank = replicated.canonical_target.tolist().index('plx5622') + 1
    # g: the judgement version in statistical order, one entry per intervention and GEO series
    ouabain_rank = check_figure3_g(gpicks, judgements, keys)
    assert len(annotations[annotations.intervention.str.casefold().eq('ouabain')]) == 1
    # the followed-up interventions carry the dagger of panel a in f and g
    assert [d for d in fpicks.display if d.endswith(' †')] == ['PLX5622 †']
    assert [d for d in gpicks.display if d.endswith(' †')] == ['Ouabain †']
    assert gpicks.loc[gpicks.display.eq('Ouabain †'), 'gse_id'].tolist() == ['GSE122080']
    # Statistical ranks of f and g, recomputed from the judged ranking: the lower-output comparisons ordered by the
    # absolute CI bound nearest zero. f takes the best rank among the comparisons whose canonical_target is the
    # meta-analysis canonical_target; g takes the rank of the comparison drawn.
    lower = answers[answers.hedges_g < 0]
    lower = lower.iloc[np.argsort(-lower.rank_bound.abs().to_numpy(), kind='stable')].reset_index(drop=True)
    lower['recomputed_rank'] = np.arange(1, len(lower) + 1)
    assert (lower.recomputed_rank == lower.statistical_rank).all()
    expected_f = {}
    for t in fpicks.canonical_target:
        hit = lower[lower.canonical_target == t]
        expected_f[t] = str(int(hit.recomputed_rank.min())) if len(hit) else 'not ranked'
    assert fpicks.rank_text.tolist() == [expected_f[t] for t in fpicks.canonical_target]
    plx = fpicks[fpicks.canonical_target == 'plx5622'].iloc[0]
    assert plx.rank_text == expected_f[plx.canonical_target]
    plx_source = lower[lower.canonical_target.eq(plx.canonical_target)].iloc[0]
    assert plx.rank_gse_id == plx_source.gse_id and plx.rank_comparison_name == plx_source.comparison_name
    expected_g_rank = [str(int(lower.loc[(lower.gse_id == r.gse_id) & (lower.comparison_name == r.comparison_name)
                                          & (lower.condition_name == r.condition_name), 'recomputed_rank'].iloc[0]))
                       for r in gpicks.itertuples()]
    assert gpicks.rank_text.astype(str).tolist() == expected_g_rank
    assert ranks.printed.astype(str).tolist() == fpicks.rank_text.tolist() + expected_g_rank
    assert (ranks['of'] == len(lower)).all()
    # g labels that panels b-e give to another comparison carry the age word of their own comparison name
    for r in gpicks.itertuples():
        base = re.sub(r' at ZT\d+$', '', r.label)
        other = picks[(picks.display.str.casefold() == base.casefold())
                      & ~((picks.gse_id == r.gse_id) & (picks.comparison_name == r.comparison_name))]
        if len(other):
            word = re.search(r'\b(old|aged)\b', r.comparison_name, flags=re.I).group(1).lower()
            assert r.display == f'{base} ({word} mice)', r.display
        else:
            assert r.display in (base, base + ' †', 'Overload in senolytic-treated mice'), r.display
    # a: the labels are the 40 b-e bars at their display text plus the three followed-up comparisons (dagger)
    labelled = set(zip(annotations.gse_id, annotations.comparison_name, annotations.label))
    bars_labelled = set(zip(picks.gse_id, picks.comparison_name, picks.display.astype(str)))
    tested_keys = next(ast.literal_eval(n.value) for n in ast.parse(
        (HERE / 'create_combined_modality_volcano.py').read_text()).body
        if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'TESTED_KEYS' for t in n.targets))
    assert bars_labelled <= labelled and len(labelled) == len(annotations)
    assert set(zip(annotations.gse_id, annotations.comparison_name)) == \
        set(zip(picks.gse_id, picks.comparison_name)) | set(tested_keys)
    assert all(lab.endswith(' †') for g_, c_, lab in labelled - bars_labelled)
    picks.to_csv(folder / 'selected_rankings.csv', index=False)
    grid.to_csv(folder / 'condition_grid.csv', index=False)
    ranks.to_csv(folder / 'statistical_ranks.csv', index=False)
    fpicks.to_csv(folder / 'panel_f_rows.csv', index=False)
    gpicks.to_csv(folder / 'panel_g_rows.csv', index=False)
    inputs = []
    for path in [rank_path, judgement_path, G_SOURCE_PATH(), replication_path, answers_path,
                 WORK/'result/volcano_ranked_comparisons_july.csv',
                 WORK/'figures/combined_3modality_plot_data.csv',
                 b.PREVIOUS/'data/verified_plotting_row_provenance_2026-11-05.csv']:
        inputs.append({'path':str(path.resolve()), 'mtime_utc':datetime.fromtimestamp(
            path.stat().st_mtime, timezone.utc).isoformat(),
            'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    b.write_json(folder/'input_provenance.json', inputs)
    rules = {
        'b_c_e':'Per category (display_best rows of dedup_interventions_2026_11_06.py): the first 5 with hedges_g>0 and the first 5 with hedges_g<0 by abs(rank_bound), stable sort descending.',
        'd':'As b, c and e.',
        'tags':'Bars with hedges_g<0 carry three cells from the agreed judgement answers: old at start (OLD met; ADULT, YOUNG, GERMLINE_OR_PRENATAL not met; NOT_STATED or NaN unknown), otherwise healthy (background HEALTHY met, which includes a disease compared with healthy controls; DISEASE_MODEL not met; NOT_STATED or NaN unknown), untreated control (UNTREATED met; ACTIVE_COMPARATOR not met; NaN unknown). NaN means the two models did not agree. The reason text lists the unmet conditions (germline, adult, young, disease model, active control) or reads all three met.',
        'f':f'All {len(fpicks)} robust replicated interventions with lower outputs (primary meta-analysis), abs(rank_bound) = -ci_high descending. Columns: k projects, n_series series, best single rank (of {len(lower)} decreases) = smallest statistical_rank among the {len(lower)} lower-output comparisons of ranked_with_judgement.csv whose canonical_target equals the meta-analysis canonical_target (equal to meta_target of volcano_ranked_comparisons_july.csv wherever both are set), or not ranked when none of them is among the {len(lower)}. PLX5622, followed up in aged mice, is bold and marked with the dagger of panel a.',
        'g':G_RULE,
        'bar_length':'b-e: absolute confidence bound nearest zero from the supplied corrected Hedges g CI; f: the bound of the pooled interval nearest zero; g: the absolute adjusted confidence bound nearest zero.',
        'bar_color':BAR_COLOR_RULE,
        'palette':PALETTE_RULE,
        'g_label_changes':{r.label:r.display for r in gpicks.itertuples() if r.label != r.display},
        'a':'Labels are the 40 b-e bars at their display text plus the three comparisons followed up in aged mice (dagger), as in commit 038eaa416 except 4°C (GSE201135), which shares canonical_target cold exposure with the followed-up GSE165974 and reads Cold exposure in a and e. Constant July unscored-point size, including nonsignificant points; no Final score legend. Display caps marked on the axes: x ticks at -30, -15, 0, 15, 30 with the ends reading ≤−30 and ≥30, colour-bar ends ≤−30 and ≥30, y cap ≥15; the dashed line is labelled FDR = 0.05; the y-title subscript is 6 pt.',
        'volcano_label_columns':f'July column formulas, fixed y strip, fonts{", colours" if PALETTE_NAME == "july" else ""} and axis limits retained. Labels are ordered within each column so that the point farthest out takes the free row nearest its height, then pairwise moves remove remaining connector crossings (_order, _untangle).',
        'layout':'183 mm page, panels at 100%: a 183 x 86 mm; b-e 183 x 106.7 mm (two key lines, the second ending with the dagger entry; ranking header; 2 x 2 grid); f 80.4 x 38.2 mm and g 99.6 x 38.2 mm side by side, sharing axes height and baseline, headers on the panel-letter baseline. Bar labels at the left page edge (b, d, f) start at least 2 mm in.',
    }
    b.write_json(folder/'selection_rules.json',rules)
    script_map = [('create_combined_modality_volcano.py', ORIGINAL/'create_combined_modality_volcano.py'),
                  ('create_combined_modality_top_hits_v2.py', ORIGINAL/'create_combined_modality_top_hits_v2.py'),
                  ('create_figure3_replication_2026_11_07.py', ORIGINAL/'create_combined_modality_top_hits_v2.py'),
                  ('create_figure3_judgement_2026_11_06.py', ORIGINAL/'create_combined_modality_top_hits_v2.py'),
                  ('assemble_figure3_vector.py', ROOT/'paper/legacy_reference/figure_rev1/assemble_figure3_vector.py')]
    diffs=[]
    for name, src in script_map:
        dest=HERE/name
        delta=list(difflib.unified_diff(src.read_text().splitlines(True),dest.read_text().splitlines(True),
                   fromfile=str(src),tofile=str(dest)))
        dest.with_suffix('.diff').write_text(''.join(delta))
        diffs.append({'script':name,'added_lines':sum(x.startswith('+') and not x.startswith('+++') for x in delta),
                      'removed_lines':sum(x.startswith('-') and not x.startswith('---') for x in delta),
                      'original_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),
                      'copy_sha256':hashlib.sha256(dest.read_bytes()).hexdigest()})
    pd.DataFrame(diffs).to_csv(folder/'script_diff_summary.csv',index=False)
    shared_path=OUT/'script_diff_summary.csv'
    shared=pd.read_csv(shared_path)
    updates=pd.DataFrame(diffs)
    updates['unchanged_nontext_plotting_ast']=False
    # Numeric figure edits are explicitly authorized and recorded in the diffs.
    shared=shared[~shared.script.isin(updates.script)]
    pd.concat([shared,updates],ignore_index=True).to_csv(shared_path,index=False)

    # Verify the July grid's bar drawing code character for character, excluding titles and axis text.
    def drawing_block(path):
        text=path.read_text()
        start=text.index('        sub = sub.sort_values("composite", ascending=True)')
        marker='    # Judgement-version tags' if '    # Judgement-version tags' in text else '    # Shared legend on the figure'
        end=text.index(marker,start)
        block=text[start:end]
        # sizes were set for print on 6 November 2026 (night); the drawing structure must still match July
        block=re.sub(r'(fontsize|markersize|linewidth|markeredgewidth|labelsize|height)=[0-9.]+',r'\1=N',block)
        return '\n'.join(l for l in block.splitlines() if 'ax.set_title(' not in l and 'ax.set_xlabel(' not in l)
    assert drawing_block(ORIGINAL/'create_combined_modality_top_hits_v2.py') == drawing_block(HERE/'create_combined_modality_top_hits_v2.py')
    pdf=folder/'figure3.pdf'
    raster(pdf).save(folder/'figure3_review.png',dpi=(200,200))
    pair(OUT/'july_reference/figure3_july.pdf',pdf,folder/'figure3_vs_july.png')
    panel_list=[('a','volcano_combined_modalities'),('b_e','top_hits_barplot_combined_modalities_v2'),
                ('f','replication_barplot_combined_modalities'),('g','judgement_barplot_combined_modalities')]
    source_rules={name:rules[letter] if letter in ('a','f','g') else rules['b_c_e']+' '+rules['tags']
                  for letter,name in panel_list}
    ledger=[];changes=[];checks=[]
    for letter,name in panel_list:
        dest=folder/f'panel_{letter}.pdf'
        shutil.copy2(WORK/'figures'/(name+'.pdf'),dest)
        raster(dest).save(folder/f'panel_{letter}_review.png',dpi=(200,200))
        new_lines=text_lines(dest)
        (folder/f'panel_{letter}_text.txt').write_text('\n'.join(new_lines)+'\n')
        old_lines=text_lines(ORIGINAL/'figures'/(name+'.pdf')) if letter not in ('f','g') else []
        if letter not in ('f','g'): pair(ORIGINAL/'figures'/(name+'.pdf'),dest,folder/f'panel_{letter}_vs_july.png')
        changes.append({'panel':letter,'source':name,
            'july_only_text':list((collections.Counter(old_lines)-collections.Counter(new_lines)).elements()),
            'new_only_text':list((collections.Counter(new_lines)-collections.Counter(old_lines)).elements())})
        for line in new_lines:
            if re.search(r'\d',line):ledger.append({'panel':letter,'printed_text':line,'source_or_rule':source_rules[name]})
        with fitz.open(dest) as d:
            checks.append({'panel':letter,'page_size_pt':list(d[0].rect)[2:],
                'font_count':len(d[0].get_fonts()),'image_count':len(d[0].get_images()),
                'missing_glyph':'\ufffd' in d[0].get_text() or '\u25a1' in d[0].get_text()})
            if letter=='a': assert 'Final score' not in d[0].get_text()
    source=folder/'panel_b_e.pdf'
    with fitz.open(source) as d:
        w,h=d[0].rect.width,d[0].rect.height
        spans=[sp for blk in d[0].get_text('dict')['blocks'] if blk['type']==0 for ln in blk['lines'] for sp in ln['spans']]
    # rows start just above the subplot titles (print-size grid of 6 November 2026, night)
    title_top={t:min(sp['bbox'][1] for sp in spans if sp['text'].strip()==t) for t in ('Genetic','Disease')}
    top,mid=title_top['Genetic']-2,title_top['Disease']-2
    for label,x,y in [('b',0,0),('c',1,0),('d',0,1),('e',1,1)]:
        crop_pdf(source,fitz.Rect(x*w/2,top if y==0 else mid,(x+1)*w/2,mid if y==0 else h),folder/f'panel_{label}.pdf')
        pair(OUT/'july_reference'/f'figure3_panel_{label}.pdf',folder/f'panel_{label}.pdf',folder/f'panel_{label}_vs_july.png')
    pd.DataFrame(ledger).to_csv(folder/'printed_numbers.csv',index=False)
    b.write_json(folder/'text_changes_vs_july.json',changes)
    b.write_json(folder/'pdf_checks.json',checks)
    assert all(c['font_count']>0 and not c['missing_glyph'] for c in checks)
    with fitz.open(pdf) as d:size=list(d[0].rect)[2:]
    # Validate the new f label and every rendered bar against its source number.
    for name,selected in [('top_hits_barplot_combined_modalities_v2',picks),
                          ('replication_barplot_combined_modalities',fpicks),
                          ('judgement_barplot_combined_modalities',gpicks)]:
        record=json.loads((WORK/'logs'/(name+'_figure_record.json')).read_text())
        if name.startswith('judgement'):
            assert next(t for t in record['axes'][0]['ytick_labels'] if t['text']=='Ouabain †')['weight']=='bold'
            assert np.allclose([x['width'] for x in record['axes'][0]['bars']],gpicks.statistical_score.iloc[::-1])
        elif name.startswith('replication'):
            assert next(t for t in record['axes'][0]['ytick_labels'] if t['text']=='PLX5622 †')['weight']=='bold'
            assert np.allclose([x['width'] for x in record['axes'][0]['bars']],fpicks.composite.iloc[::-1])
        else:
            hatch={'mouse RNA-seq':'','human RNA-seq':'////////','human DNAm':'........'}
            for axis,cat in zip(record['axes'],['Genetic','Drug','Disease','Environment']):
                expected=selected[selected.category.eq(cat)].iloc[::-1]
                assert np.allclose([x['width'] for x in axis['bars']],expected.composite)
                assert [x['hatch'] or '' for x in axis['bars']]==[hatch[m] for m in expected.modality]
                assert [t['text'] for t in axis['ytick_labels']]==expected.display.astype(str).tolist()
                assert axis['xlabel']=="|Hedges' g| (CI bound nearest zero)"
    # Headers, column heads, reasons and ranks as printed in the panel PDFs.
    def pdf_spans(path):
        with fitz.open(path) as d:
            return [sp for blk in d[0].get_text('dict')['blocks'] if blk['type']==0 for ln in blk['lines'] for sp in ln['spans']]
    def column_values(spans, head):
        h=[sp for sp in spans if sp['text'].strip()==head];assert len(h)==1
        x0,x1,y=h[0]['bbox'][0],h[0]['bbox'][2],h[0]['bbox'][3]
        rows=[sp for sp in spans if sp is not h[0] and sp['bbox'][1]>=y-1 and x0-1<=(sp['bbox'][0]+sp['bbox'][2])/2<=x1+1]
        return [sp['text'].strip() for sp in sorted(rows,key=lambda sp:sp['bbox'][1])]
    headers={'b_e':('Statistical ranking','of single comparisons with the 5 strongest increases and decreases per category'),
             'f':('Replication ranking','pooled across independent projects'),
             'g':('Judgement ranking','in old, otherwise healthy animals vs untreated controls')}
    texts={}
    for letter,name in panel_list[1:]:
        spans=pdf_spans(folder/f'panel_{letter}.pdf');texts[letter]=[sp['text'].strip() for sp in spans]
        bold,rest=headers[letter]
        hb=[sp for sp in spans if sp['text'].strip()==bold];hr=[sp for sp in spans if sp['text'].strip()==rest]
        assert len(hb)==len(hr)==1 and 'Bold' in hb[0]['font'] and 'Bold' not in hr[0]['font']
        assert abs(hb[0]['origin'][1]-hr[0]['origin'][1])<0.01 and hb[0]['size']==hr[0]['size']==8
    assert geometry['header']==list(headers['b_e'])
    assert [texts['b_e'].count(h) for h in ('Old at start','Otherwise healthy','Untreated control')]==[4,4,4]
    assert 'No disease model' not in texts['b_e']
    assert geometry['column_heads']==['Old at start','Otherwise healthy','Untreated control']
    assert texts['b_e'].count("|Hedges' g| (CI bound nearest zero)")==4
    for key in ('condition met','not met','not stated or not agreed','lower clock output','higher clock output',
                '† followed up in aged mice'):
        assert texts['b_e'].count(key)==1, key
    assert 'Cold exposure' in texts['b_e'] and '4°C' not in texts['b_e']
    # bar labels at the left page edge start at least 2 mm in (b and d of the grid, f)
    mm=72/25.4
    left_labels=[sp for sp in pdf_spans(folder/'panel_b_e.pdf') if sp['text'].strip() in set(picks.display)
                 and sp['bbox'][0]<w/2]
    assert len(left_labels)==20 and min(sp['bbox'][0] for sp in left_labels)>=2.0*mm-0.1
    assert min(sp['bbox'][0] for sp in pdf_spans(folder/'panel_f.pdf') if sp['text'].strip() in set(fpicks.display))>=2.0*mm-0.1
    # panel a: the display caps are marked on both axes and the colour bar, the FDR line is labelled and the y-title
    # subscript is 6 pt
    a_spans=pdf_spans(folder/'panel_a.pdf');a_text=[sp['text'].strip() for sp in a_spans]
    assert a_text.count('≤−30')==2 and a_text.count('≥30')==2 and a_text.count('≥15')==1, a_text
    assert '−40' not in a_text and '40' not in a_text and a_text.count('FDR = 0.05')==1
    sub_spans=[sp for sp in a_spans if sp['text'].strip()=='10' and abs(sp['size']-6.0)<0.01]
    assert len(sub_spans)==1 and min(sp['size'] for sp in a_spans)>=6.0-1e-6, min(sp['size'] for sp in a_spans)
    assert {'−log','FDR'}<=set(a_text)
    reasons=sorted(r for r in grid.reason.dropna() if r)
    printed=sorted(t for t in texts['b_e'] if t in set(reasons))
    assert printed==reasons, (printed,reasons)
    assert all(min(v)>=2.0 for v in geometry['reason_marker_clearance_mm'].values())
    assert column_values(pdf_spans(folder/'panel_f.pdf'),'Best single rank')==[f'(of {len(lower)} decreases)']+fpicks.rank_text.tolist()
    assert column_values(pdf_spans(folder/'panel_f.pdf'),'Projects')==[str(int(k)) for k in fpicks.k]
    assert column_values(pdf_spans(folder/'panel_f.pdf'),'Series')==[str(int(n)) for n in fpicks.n_series]
    assert column_values(pdf_spans(folder/'panel_g.pdf'),'Statistical rank')==[f'(of {len(lower)} decreases)']+gpicks.rank_text.astype(str).tolist()
    assert 'Statistical rank' not in texts['f'] and f'(of {len(lower)})' not in texts['f'] + texts['g']
    assert "Pooled |Hedges' g| (CI bound nearest zero)" in texts['f'] and G_XLABEL in texts['g']
    assert G_XLABEL in texts['g']
    for r in gpicks.itertuples(): assert r.display in texts['g']
    page_text=[sp['text'].strip() for sp in pdf_spans(pdf)]
    for bold,rest in headers.values(): assert bold in page_text and rest in page_text
    for old in ('Replicated across','independent projects','Old, otherwise healthy animals','vs untreated controls',
                'Top 5 interventions per category and direction by Hedges\' g (CI bound nearest zero)'):
        assert old not in page_text, old
    with fitz.open(folder/'panel_f.pdf') as df_, fitz.open(folder/'panel_g.pdf') as dg_:
        assert abs(df_[0].rect.height-dg_[0].rect.height)<0.01
    frec=json.loads((WORK/'logs/replication_barplot_combined_modalities_figure_record.json').read_text())['axes'][0]
    grec=json.loads((WORK/'logs/judgement_barplot_combined_modalities_figure_record.json').read_text())['axes'][0]
    fh,gh=frec['rect_fraction'],grec['rect_fraction']
    # same axes height and baseline, same row pitch
    assert np.isclose(fh[1],gh[1]) and np.isclose(fh[3],gh[3]) and np.isclose(np.diff(frec['ylim']),np.diff(grec['ylim'])).all()
    plotted = pd.read_csv(WORK / 'figures/combined_3modality_plot_data.csv', low_memory=False)
    modality_counts = plotted.modality.value_counts().to_dict()
    numbers={'atlas_rows':len(plotted),'significant_rows':int(plotted.fdr.lt(.05).sum()),
             'category_bars':len(picks),'replication_bars':len(fpicks),
             'judgement_bars':len(gpicks),'plx5622_replication_rank':plx5622_replication_rank,'ouabain_judgement_rank':ouabain_rank,'page_size_pt':size,
             'headers':{k:' '.join(v) for k,v in headers.items()},
             'lower_output_ranked_comparisons':int(len(lower)),
             'statistical_ranks':ranks.replace({np.nan:None}).to_dict('records'),
             'condition_states':{c:grid[c].value_counts().to_dict() for c in state_cols},
             'condition_grid':grid.replace({np.nan:None}).to_dict('records'),
             'reason_marker_clearance_mm':geometry['reason_marker_clearance_mm'],
             'f_rows':fpicks.replace({np.nan:None}).to_dict('records'),
             'g_rows':gpicks.replace({np.nan:None}).to_dict('records')}
    b.write_json(folder/'numbers.json',numbers)
    diff_lines='\n'.join(f'| [{r["script"]}]({HERE/r["script"]}) | +{r["added_lines"]} / −{r["removed_lines"]} | [Diff]({HERE/Path(r["script"]).with_suffix(".diff")}) |' for r in diffs)
    f_lines='\n'.join(f'| {i+1} | {r.display} | {int(r.k)} | {int(r.n_series)} | {r.estimate:.6f} | [{r.ci_low:.6f}, {r.ci_high:.6f}] | {r.composite:.6f} | {r.fdr:.6g} | {r.rank_text} | {"" if pd.isna(r.rank_gse_id) else f"{r.rank_gse_id}, {r.rank_comparison_name}"} |' for i,r in enumerate(fpicks.itertuples()))
    g_section=figure3_g_readme(gpicks)
    unranked=' and '.join(fpicks.loc[fpicks.rank_text.eq('not ranked'),'display'])
    shown_state={'met':'met','not met':'not met','unknown':'not stated or not agreed'}
    state_counts='; '.join(f"{h}: " + ', '.join(f'{n} {shown_state[s]}' for s,n in grid[c].value_counts().items())
                           for h,c in zip(('old at start','otherwise healthy','untreated control'),state_cols))
    first_lines='\n'.join(f'| {cat} | {group.iloc[0].intervention} | {group.iloc[0].composite:.8g} |' for cat,group in picks.groupby('category',sort=False))
    readme=f'''# Figure 3 reproduction

Figure 3 shows three versions of the ranking: the statistical ranking (b–e), the
replication ranking (f) and the judgement ranking (g). It uses copies of the
July generators; f and g copy the July bar-drawing block. The assembly keeps the
July full-width volcano and bar-grid formulas and places f and g side by side
below. The page is {size[0]:.6f} × {size[1]:.6f} pt. Review images use 200 dpi.

## Panel map

| Panel | Content | Selection | Panel PDF |
|---|---|---|---|
| a | Combined-modality volcano | The 40 b–e bars at their display text plus the three comparisons followed up in aged mice (†) | [a](panel_a.pdf) |
| b–e | {' '.join(headers['b_e'])} | Per category the top 5 with Hedges' g > 0 and the top 5 with g < 0, one per intervention and GEO series; each lower bar carries three cells for the judgement conditions and its unmet conditions as text | [b](panel_b.pdf) [c](panel_c.pdf) [d](panel_d.pdf) [e](panel_e.pdf) |
| f | {' '.join(headers['f'])} | The {len(fpicks)} robust replicated interventions with lower outputs; columns give projects, series and best single rank | [f](panel_f.pdf) |
| g | {' '.join(headers['g'])} | {G_MAP_SELECTION} | [g](panel_g.pdf) |

Bars in b–e give the absolute CI bound of Hedges' g nearest zero, bars in f the
bound of the pooled interval nearest zero and bars in g the adjusted confidence bound nearest zero.
{BAR_COLOR_SENTENCE}; marker shapes and
hatching give the modality. PLX5622 (f, {ORDINAL.get(plx5622_replication_rank, str(plx5622_replication_rank))}) and ouabain (g, {ORDINAL.get(ouabain_rank, str(ouabain_rank))}), followed
up in aged mice, are bold and carry the dagger that marks the three followed-up
comparisons in a; the second key line ends with "† followed up in aged mice". The
bar 4°C in e (GSE201135) shares the canonical target cold exposure with the
followed-up comparison GSE165974 and reads "Cold exposure" in a and e.

The three cells beside each lower bar in b–e are old at start, otherwise healthy
and untreated control, the three conditions of the g header: filled (`{_PAL['met_fill']}`) when the
agreed answer meets the condition, crossed when it does not, light grey when the
answer is not stated or the two models did not agree. Otherwise healthy is the
background answer HEALTHY, which the judging models also give to a disease
compared with healthy controls (NASH vs Normal), and not met is DISEASE_MODEL, a
treatment given to a disease model. Across the 20 lower bars, {state_counts}.
Only Mechanical Overload (GSE195707) meets all three, and it is the only one of
these bars admitted to the judgement ranking.

The rank columns count the {len(lower)} ranked comparisons with lower outputs
in the judged ranking (`statistical_rank`, ordered by the absolute CI bound
nearest zero), headed "(of {len(lower)} decreases)". In f, headed "Best single
rank", it is the best rank among the comparisons whose `canonical_target` is the
meta-analysis `canonical_target`; that column equals `meta_target` of the July
ranking table wherever both are set. {unranked} have no comparison among the
{len(lower)} and read not ranked. In g, headed "Statistical rank", it is the rank
of the comparison drawn, as in Table 3. Labels in g that b–e use for a different
comparison take the age word of their own comparison name.

In a, Cohen's d is drawn within ±30 and −log₁₀ FDR up to 15; the x ticks end at
the cap and read ≤−30 and ≥30, as do the colour-bar ends, matching ≥15 on y. The
dashed line is labelled FDR = 0.05, and the subscript of the y title is set at
6 pt. {VOLCANO_COLOR_SENTENCE}

## Data sources

The category selections come from `{rank_path}` through the display table of
`dedup_interventions_2026_11_06.py`. The b–e cells restate the agreed answers in
`{answers_path}`. Panel f comes from `{replication_path}` and panel g from
`{judgement_path}`, ordered by statistical rank. These are supplied model judgements and meta-analysis
results; nothing was judged or pooled during plotting.

The volcano uses the supplied verified atlas with {len(plotted):,} rows and {numbers['significant_rows']:,} rows
with FDR < 0.05. Modality counts are {modality_counts.get('mouse RNA-seq', 0):,} mouse RNA-seq, {modality_counts.get('human RNA-seq', 0):,} human RNA-seq,
and {modality_counts.get('human DNAm', 0):,} human DNAm. Rows retain separate clocks. Verified Cohen's d and
BH-FDR supply the volcano positions. Historical `effect_size`, `pvalue`, `fdr`,
and `final_score` aliases in ranking and judgement exports do not replace
corrected Hedges' g or the `v3_` statistics reported below. The pre-existing
atlas adapter uses `outlier_status=Normal` as an inclusion literal, not a
rerun of the historical outlier classifier. No measurements are simulated.

[Input provenance](input_provenance.json) records paths, UTC modification times,
and SHA-256 hashes. [Selection rules](selection_rules.json) records filters and
ordering. [Selected rankings](selected_rankings.csv) contains the 40 category
rows with their condition states and reason text, [condition grid](condition_grid.csv)
the 20 lower-output bars with their source answers, [statistical ranks](statistical_ranks.csv)
the rank beside each f and g bar with its source comparison, [panel f rows](panel_f_rows.csv)
the four replicated interventions and [panel g rows](panel_g_rows.csv) the five
judgement rows.
The original plotting-source provenance remains in
[the July input record](../july_exact/data/input_provenance.json).

The first ranked comparison in each category is shown below.

| Category | First comparison | Bar length |
|---|---|---:|
{first_lines}

## Panel f numbers

Copied from the replication source; k is projects, series is GEO series.

| Rank | Intervention | k | Series | Pooled g | 95% CI | Bar length | FDR | Best single rank | From |
|---:|---|---:|---:|---:|---|---:|---:|---:|---|
{f_lines}

{g_section}

## Script diffs

The unified diffs compare each generator against its July original. Panels f
and g are compared against the July top-hits generator from which their bar
code was copied. [Diff summary](script_diff_summary.csv) includes hashes and counts.

| Script | Added / removed | Unified diff |
|---|---:|---|
{diff_lines}

The b–e bar drawing block matches July character for character after excluding
axis-label and title text; the condition cells are drawn after it in a fixed
column at the right of each subplot, and an x axis is lengthened only where a
reason text would come within 2 mm of its bar's end marker. The f and g scripts
copy that drawing block, set the ranking header on the panel-letter baseline and
the value columns right of the bars, and bold the followed-up intervention and mark
it with the dagger. Assembly places f and g
side by side below the grid and sets the b–e letters level with the subplot
titles found in the grid PDF.

The volcano adds the three comparisons followed up in aged mice without
displacing the category-budget labels, {"with the same July font, colour and bracket style" if PALETTE_NAME == "july" else "in the July font and bracket style and the category colour of the other labels"}. All points now use the July unscored default size,
`(80/170)^3*120+4` = {((80/170)**3*120+4):.9f} pt²; the inherited outline
multiplier is retained. The unsupported Final score legend is removed.
Label columns, axis limits{", fonts and colours are July's" if PALETTE_NAME == "july" else " and fonts are July's"}; within each column the
labels are ordered so that fewer connectors cross.

## Verification

Category membership and order were independently recomputed from the supplied
ranking CSV. Every rendered bar width was checked against its source confidence bound.
Every b–e cell state and reason text was recomputed from the agreed judgement
answers, and the bars equal those of commit 038eaa416. The statistical ranks of
f and g were recomputed from the judged ranking and read back from the panel
PDFs, as were the headers, column heads, key entries, x-axis titles, the
capped tick labels of a and the left label margin. Panel f membership, order,
pooled estimates, intervals, project and series counts and the bold PLX5622 †
label, and panel g membership, order, confidence bounds nearest zero, P values,
FDR values and the bold Ouabain † label were checked against their sources and the rendered figure
records.
[PDF checks](pdf_checks.json) record embedded fonts, image counts, page sizes,
and missing-glyph checks. [Printed numbers](printed_numbers.csv) and
[text changes](text_changes_vs_july.json) list PDF-extractable text; July's
path-effect volcano annotations render as paths and are checked through the
annotation CSV and review image instead of PDF text extraction.

Generic nature-figure preflight findings are recorded without applying style
changes beyond the requested scope. This Figure 3-only command does not
write Figure 4 or Extended Data outputs.

## Reproduction command

Run the Figure 3-only path to preserve the other figure outputs.

```sh
source {ROOT.parent/'env_local.sh'}
python {HERE/'build_july_exact_2026_11_06.py'} figure3
python {HERE/'package_figures_2026_11_06.py'} figure3
```
'''
    (folder/'README.md').write_text(readme)
    old_helpers=folder/'outdated/before_direction_selection_2026-11-06/scripts'
    for helper in ['build_july_exact_2026_11_06.py','package_figures_2026_11_06.py']:
        delta=''.join(difflib.unified_diff((old_helpers/helper).read_text().splitlines(True),
            (HERE/helper).read_text().splitlines(True),fromfile=str(old_helpers/helper),tofile=str(HERE/helper)))
        (folder/(helper.removesuffix('.py')+'_followup.diff')).write_text(delta)

    # Verify the outer PDFs for a-e have not been scaled down to accommodate f.
    old_panels=folder/'outdated/before_direction_selection_2026-11-06/artifacts'
    scale_checks=[]
    for name in ['panel_a.pdf','panel_b_e.pdf']:
        with fitz.open(old_panels/name) as old_doc, fitz.open(folder/name) as new_doc:
            old_rect,new_rect=old_doc[0].rect,new_doc[0].rect
            ratio=old_rect.width/new_rect.width
            assert ratio>=1-1e-4, f'Existing panel shrunk during assembly: {name}'  # tight-bbox jitter below 0.1 pt
            scale_checks.append({'panel':name,'old_source_size_pt':list(old_rect)[2:],
                'new_source_size_pt':list(new_rect)[2:],'new_to_old_assembly_scale':ratio})
    # Confirm all pre-existing protected outputs retain their bytes.
    baseline=json.loads((folder/'outdated/before_direction_selection_2026-11-06/protected_hashes.json').read_text())
    concurrent_changes=[]
    for name,sha in baseline.items():
        if hashlib.sha256((OUT/name).read_bytes()).hexdigest()!=sha:
            concurrent_changes.append(name)
    b.write_json(folder/'verification.json',{'category_membership_and_order':'passed',
        'judgement_membership_and_numbers':'passed','rendered_bar_widths':'passed',
        'bold_ouabain_in_f':'passed','july_bar_drawing_block':'passed',
        'bars_match_selected_run':'passed','condition_grid_vs_source':'passed',
        'admitted_condition_states_match_source':'passed','statistical_ranks_vs_source':'passed',
        'plx5622_statistical_rank_matches_source':'passed',
        'volcano_label_set_unchanged_except_cold_exposure_4C':'passed','dagger_in_a_f_g_and_key':'passed',
        'otherwise_healthy_column_head':'passed','panel_a_caps_fdr_label_subscript_6pt':'passed',
        'left_label_margin_2mm':'passed',
        'headers_column_heads_axis_titles_in_pdf':'passed','reason_marker_clearance_min_mm':
            min(min(v) for v in geometry['reason_marker_clearance_mm'].values()),
        'protected_outputs_unchanged':len(baseline)-len(concurrent_changes),
        'other_outputs_changed_concurrently':concurrent_changes,
        'a_e_scale_checks':scale_checks,
        'scope':'This command writes only Figure 3 outputs and associated source diffs.',
        'page_size_pt':size})
    print(pd.DataFrame(diffs)[['script','added_lines','removed_lines']].to_string(index=False))
    print('Packaged Figure 3 only:',folder)


def main():
    try:
        if len(sys.argv)>1 and sys.argv[1]=="figure3":
            package_figure3_followup()
            return
        diffs=diff_audit();package();metrics_and_readmes(diffs)
        print(diffs[['script','added_lines','removed_lines']].to_string(index=False))
        print('Packaged',OUT)
    except Exception as exc:
        print(f'Packaging failed: {exc}',file=sys.stderr)
        raise


if __name__=='__main__': main()
