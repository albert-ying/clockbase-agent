#!/usr/bin/env python3
"""Export the verified atlas in the July plotting schema without historical statistics."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-05/data'
FINAL = Path(os.environ['CLOCKBASE_FINAL'])
KEY = ['modality', 'gse_id', 'comparison_name', 'condition_name']

def receipt(path):
    return {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'mtime_utc': datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat()}

def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src = FINAL/'public-export/atlas_v3_combined.csv'
    cls = FINAL/'public-export/classification/comparisons.jsonl'
    schema = ROOT/'figures/combined_3modality_plot_data.csv'
    atlas = pd.read_csv(src, low_memory=False)
    atlas = atlas[atlas.v3_final_class.isin(['VERIFIED','PROVISIONAL'])].copy()
    categories = pd.read_json(cls, lines=True)[KEY+['category']]
    atlas = atlas.merge(categories,on=KEY,how='left',validate='many_to_one',sort=False)
    columns = pd.read_csv(schema,nrows=0).columns.tolist()
    out = pd.DataFrame(index=atlas.index,columns=columns)
    mapping = {k:k for k in KEY}
    mapping.update(effect_size='v3_cohen_d',log2FoldChange='v3_cohen_d',
                   fdr='v3_fdr',pvalue='v3_p',intervention='intervention',case_id='case_id')
    for target,source in mapping.items(): out[target]=atlas[source]
    out['condition_category'] = atlas['category'].where(atlas.category.isin(['Genetic','Drug','Disease','Environment']),'Other')
    out['direction'] = np.sign(atlas.v3_cohen_d)
    out['neg_log10_fdr'] = -np.log10(atlas.v3_fdr)
    out['is_significant'] = atlas.v3_fdr < .05
    assert len(out)==len(atlas)
    assert list(out.columns)==columns
    assert np.isfinite(out[['effect_size','fdr','pvalue','neg_log10_fdr']].astype(float)).all().all()
    assert out.final_score.isna().all()
    name = 'combined_3modality_plot_data_verified_2026-11-05.csv'
    out.to_csv(OUT/name,index=False)
    provenance = atlas[['modality','gse_id','comparison_name','condition_name','clock','atlas_row_ordinal','v3_final_class','v3_n_case','v3_n_control','v3_design_method','v3_cohen_d','v3_p','v3_fdr','condition_category','category']].copy()
    provenance.insert(0,'plot_row_index',np.arange(len(provenance)))
    provenance.to_csv(OUT/'verified_plotting_row_provenance_2026-11-05.csv',index=False)
    counts=atlas.groupby(['modality','v3_final_class']).size().rename('rows').reset_index()
    counts.to_csv(OUT/'counts_by_modality_and_class.csv',index=False)
    summary={'rows':len(out),'comparisons':len(out[KEY].drop_duplicates()),'significant_rows':int(out.is_significant.sum()),'category_missing_rows':int(atlas.category.isna().sum()),'counts':counts.to_dict('records'),'sources':[receipt(p) for p in [src,cls,schema]],'created_utc':datetime.now(timezone.utc).isoformat(),'mapping':mapping,'empty_columns':[c for c in out if out[c].isna().all()]}
    (OUT/'provenance.json').write_text(json.dumps(summary,indent=2)+'\n')
    lines=['# Verified plotting table','',f'The table contains {len(out):,} inference-supporting atlas rows from {summary["comparisons"]:,} comparisons. Only VERIFIED and PROVISIONAL rows are included. Row order follows the atlas. Multiple clocks remain separate rows.','', '## Column mapping','', '| July column | Verified source or rule |','|---|---|']
    for c in columns:
        rule = {'condition_category':'classification/comparisons.jsonl category, joined by modality, GEO accession, comparison and condition; Genetic, Drug, Disease and Environment retained, all others and missing values become Other','direction':'sign(v3_cohen_d); positive is case greater than control','neg_log10_fdr':'−log10(v3_fdr), without a floor or clipping','is_significant':'v3_fdr < 0.05','log2FoldChange':'v3_cohen_d compatibility alias; this column is Cohen’s d, not a log fold change'}.get(c,mapping.get(c,'empty; no verified counterpart used'))
        lines.append(f'| {c} | {rule} |')
    lines+=['','The final_score column is empty. Historical score, outlier and review fields are not carried forward. Intervention and case identifiers are atlas annotations, not statistical measurements. The separate row-provenance table preserves final class, clock, sample counts, test method and atlas row ordinal without changing the July schema.','','## Row counts','','| Modality | VERIFIED | PROVISIONAL | Total |','|---|---:|---:|---:|']
    for m,g in counts.groupby('modality'):
        d=g.set_index('v3_final_class').rows
        lines.append(f'| {m} | {d.get("VERIFIED",0):,} | {d.get("PROVISIONAL",0):,} | {g.rows.sum():,} |')
    lines+=['',f'{summary["significant_rows"]:,} rows have BH-FDR below 0.05. BH-FDR is the supplied v3_fdr value, adjusted within modality; no FDR is recomputed here. {summary["category_missing_rows"]:,} rows lack the classification category and receive Other.','','## Sources','','Source timestamps and SHA-256 hashes are in provenance.json. Exact paths follow.','']
    lines += [f'- `{p}`' for p in [src,cls,schema]]
    lines+=['','## Validation','','Assertions check the selected inference-supporting row count, the exact July column order, finite effect sizes and probabilities, and an empty final_score column. Each export row has a matching row in verified_plotting_row_provenance_2026-11-05.csv.']
    (OUT/'README.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    try: main()
    except (OSError,ValueError,AssertionError) as exc:
        raise SystemExit(f'Verified plotting export failed: {exc}') from exc
