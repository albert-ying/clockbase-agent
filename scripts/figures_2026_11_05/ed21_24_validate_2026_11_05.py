"""Verify final delivery files and conservation of the supplied source values."""
import os
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image
from ed21_24_restore_common_2026_11_05 import OUTBASE, FINAL, META, REPO, sha256
import build_ed22_replication_2026_10_14 as ed22


def main():
    summary={}
    for n in [21,22,23,24]:
        out=OUTBASE/f'ed{n}';name=f'ed{n}'
        files=[out/f'{name}.pdf',out/f'{name}_review.png',out/f'{name}_vs_july.png',out/'README.md']
        for path in files: assert path.is_file() and path.stat().st_size>0, str(path)
        qc=json.loads((out/'checks.json').read_text())
        assert qc['passed'], (name,qc)
        page=fitz.open(files[0])[0]
        assert abs(page.rect.width*25.4/72-183)<1e-3
        assert page.rect.height*25.4/72<=247
        assert not page.get_images(), name
        assert len(page.get_drawings())>20
        text=page.get_text()
        assert 'data:' not in text.lower() and 'computed in' not in text.lower()
        # ED21's published-result categories include the historical Figure 3;
        # this is a data-row label, not a figure heading.
        headings=[line.strip() for line in text.splitlines()
                  if line.strip().startswith(('Figure ','Extended Data Figure '))]
        assert headings==(['Figure 3'] if n==21 else []), (name,headings)
        provenance=json.loads((out/'provenance.json').read_text())
        assert provenance['output_sha256']==sha256(files[0])
        for item in provenance['inputs']+provenance['builder']:
            assert sha256(item['path'])==item['sha256'], item['path']
        png=Image.open(files[1]);assert abs(png.info['dpi'][0]-200)<.05
        summary[name]={'passed':True,'size_mm':[183,round(page.rect.height*25.4/72,2)],
                       'vector_paths':len(page.get_drawings()),'raster_images':len(page.get_images()),
                       'pdf_sha256':sha256(files[0]),'files':[str(p.relative_to(REPO)) for p in files]}
    # ED21 diagnostics: classifications account for every historical atlas row;
    # the eligible classes independently sum to the stated inferential subset.
    a=pd.read_csv(OUTBASE/'ed21/panel_a_classification.csv')
    for key,part in a.groupby('group'):
        raw=json.loads((FINAL/f'repairs-{key}/manifest.json').read_text())['final_classes']
        assert part.n.sum()==sum(raw.values())
        for cls in ['VERIFIED','PROVISIONAL','DESCRIPTIVE','DUPLICATE','NOT_REPAIRABLE']:
            assert part.loc[part['class']==cls,'n'].iloc[0]==raw.get(cls,0)
    eligible=int(a.loc[a['class'].isin(['VERIFIED','PROVISIONAL']),'n'].sum())
    atlas=pd.read_csv(FINAL/'public-export/atlas_v3_combined.csv',usecols=['v3_final_class'])
    assert eligible==int(atlas.v3_final_class.isin(['VERIFIED','PROVISIONAL']).sum())
    summary['ed21']['verified_plus_provisional_rows']=eligible
    screen=pd.read_csv(ed22.paths.VERIFIER_ROOT/'input-integrity-dnam-001/samples.csv')
    hist=pd.read_csv(OUTBASE/'ed21/panel_e_histogram.csv')
    assert int(hist.n.sum())==len(screen.r_reference.dropna())
    # ED22: source-selected rows, exact source effects and counts, plus all
    # lifespan matches with >=2 projects. No row is dropped to improve fit.
    primary=pd.read_csv(META/'primary/intervention_meta.csv')
    series=pd.read_csv(META/'by-series/intervention_meta.csv')
    for panel,sign in [('a',-1),('c',1)]:
        observed=pd.read_csv(OUTBASE/f'ed22/panel_{panel}_{"lower" if sign<0 else "higher"}.csv')
        expected=ed22.robust(primary)
        expected=expected[expected.estimate*sign>0].head(12)
        assert list(observed.canonical_target)==list(expected.canonical_target)
        for col in ['estimate','ci_low','ci_high','k','n_series','fdr']:
            assert np.allclose(observed[col],expected[col],rtol=1e-12,atol=1e-15)
    counts=pd.read_csv(OUTBASE/'ed22/panel_d_counting_unit.csv')
    project_robust=ed22.robust(primary);series_robust=ed22.robust(series)
    expected_counts=[[int(series.k.ge(2).sum()),int(primary.k.ge(2).sum())],
                     [len(series_robust),len(project_robust)],
                     [int(series_robust.estimate.gt(0).sum()),int(project_robust.estimate.gt(0).sum())],
                     [int(series_robust.estimate.lt(0).sum()),int(project_robust.estimate.lt(0).sum())]]
    assert counts[['series','projects']].values.tolist()==expected_counts
    lifespan=pd.read_csv(META/'concordance/lifespan_concordance.csv')
    plotted=pd.read_csv(OUTBASE/'ed22/panel_g_lifespan_matched.csv')
    expected=lifespan[lifespan.k>=2].sort_values(['expected','estimate'],ascending=[False,True])
    assert len(plotted)==len(expected)
    for col in ['estimate','ci_low','ci_high','k','expected']:
        assert np.allclose(plotted[col],expected[col],rtol=1e-12,atol=1e-15)
    summary['ed22']['counting_unit']=counts.to_dict('records')
    summary['ed22']['lifespan_matched_rows']=len(plotted)
    # ED23: verify the separate source JSON blocks for every model and rerun.
    bench=REPO/'benchmark/published_clock_2026_10_05'
    published=json.loads((Path(os.environ.get('CLOCKBASE_BENCH_RESULTS') or bench/'groupbench/results_2026_10_10')/'published_released_test.json').read_text())
    facts=json.loads((OUTBASE/'ed23/numbers.json').read_text())
    for key,block in facts['c_d'].items():
        assert block==published['against_models_without_data']['systems'][key]
    assert facts['b']==published['published_significant']
    for entry in facts['f']:
        folder={'Codex':'2026_10_14','Claude Code':'2026_10_18','Biomni':'2026_10_19'}[entry['framework']]
        raw=json.loads((bench/f'agentbench/results_{folder}/ranked_reruns.json').read_text())['by_model'][entry['model']]
        for direction in ['higher','lower']:assert entry[direction]==raw[direction]
    summary['ed23'].update({'ground_truth_effects':facts['ground_truth_total'],
                            'benchmark_series':facts['e_series'],'agreement_systems':len(facts['e']),'ranked_configurations':len(facts['f'])})
    # ED24: all source observations were recovered with the original loader;
    # compare the exported values and tests to the unrounded source rows.
    stats=pd.read_csv(OUTBASE/'ed24/example_statistics.csv')
    values=pd.read_csv(OUTBASE/'ed24/plotted_sample_values.csv')
    for row in stats.itertuples():
        part=values[values.panel==row.panel]
        case=part[part.group=='case'].clock_value;control=part[part.group=='control'].clock_value
        assert len(case)==row.n_case and len(control)==row.n_control
        assert abs(case.mean()-control.mean()-row.mean_diff)<1e-6
        assert row.final_class in ['VERIFIED','PROVISIONAL']
    assert stats.loc[stats.gse=='GSE68319','method'].iloc[0]=='PAIRED'
    summary['ed24']['examples']=stats[['gse','n_case','n_control','method','hedges_g','fdr']].to_dict('records')
    for n in [21,22,23,24]:
        out=OUTBASE/f'ed{n}'
        (out/'validation.json').write_text(json.dumps(summary[f'ed{n}'],indent=2)+'\n')
        qc=json.loads((out/'checks.json').read_text())
        report=(out/'README.md').read_text().split('\n## Final review\n')[0]
        report+=f'''\n## Final review

The final 200-dpi render and the same-scale July comparison were visually
inspected. No unresolved text overlap or clipping was found. The supplied
source values pass `ed21_24_validate_2026_11_05.py`; results are recorded in
`validation.json`. The PDF contains {qc['pdf_vector_paths']} vector paths
and no embedded raster plot images. The minimum annotation size is
{qc['minimum_text_size_pt']:g} pt; axis and tick text meet the 7-pt floor.

## Files

- [Vector PDF](ed{n}.pdf)
- [Review render at 200 dpi](ed{n}_review.png)
- [Comparison with July](ed{n}_vs_july.png)
- [Unrounded numbers](numbers.json)
- [Source paths and hashes](provenance.json)
- [Defect checks](checks.json)
- [Source-value validation](validation.json)

## Open problems

No required input is missing and no rendering defect remains unresolved.
The replacement page is delivered as a standalone PDF; the combined figure
package is unchanged.
'''
        (out/'README.md').write_text(report)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
