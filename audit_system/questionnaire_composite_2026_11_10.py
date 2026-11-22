"""Record the saved questionnaire score and its composite with verified v3_p.

Scope: scored Genetic, Drug and Environment comparisons across both directions.
questionnaire_score = final_score. The questionnaire composite is a record,
never an admission rule, filter or displayed ranking. No model calls.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from audit_system import paths
# Same formula constants as the archived composite implementation.
CAP, P_FLOOR = 30, 1e-300
DIRECTIONS = ('age-lowering', 'age-raising')
KEY = ['modality', 'gse_id', 'comparison_name', 'condition_name']
SCOPE = ('Genetic', 'Drug', 'Environment')
from audit_system.membership_recovery_2026_09_27 import sha256

FORMULA = 'questionnaire_composite = min(-log10(max(v3_p, 1e-300)), 30) * (final_score - m) / (M - m)'


def questionnaire_table(table):
    required = KEY + ['condition_category', 'direction_v3', 'final_score', 'v3_p']
    if set(required) - set(table):
        raise ValueError(f'missing columns: {sorted(set(required)-set(table))}')
    fields = list(dict.fromkeys(KEY + ['condition_category', 'direction_v3', 'v3_p'] +
                               [c for c in ('label',) if c in table]))
    out = table[fields].copy()
    if table.duplicated(KEY).any():
        raise ValueError('duplicate comparison keys')
    if set(table.direction_v3) - set(DIRECTIONS):
        raise ValueError('unknown direction')
    out['statistical_rank_overall'] = np.arange(1, len(out)+1)
    out['statistical_rank'] = out.groupby('direction_v3').cumcount()+1
    raw_score = pd.to_numeric(table.final_score, errors='raise')
    out['questionnaire_score'] = raw_score.where(table.condition_category.isin(SCOPE))
    scored = out.condition_category.isin(SCOPE) & out.questionnaire_score.notna()
    scores = out.loc[scored, 'questionnaire_score']
    if scores.empty or not np.isfinite(scores).all():
        raise ValueError('no finite scored questionnaire scope')
    m, M = float(scores.min()), float(scores.max())
    if not np.isfinite([m, M]).all() or M == m:
        raise ValueError('questionnaire range undefined')
    p = pd.to_numeric(out.loc[scored, 'v3_p'], errors='raise')
    if not (np.isfinite(p) & p.between(0, 1)).all():
        raise ValueError('invalid verified P value in scored scope')
    out['questionnaire_composite'] = np.nan
    out.loc[scored, 'questionnaire_composite'] = np.minimum(-np.log10(np.maximum(p, P_FLOOR)), CAP)*(scores-m)/(M-m)
    out['questionnaire_composite_rank'] = pd.Series(pd.NA, index=out.index, dtype='Int64')
    for direction in DIRECTIONS:
        ordered = out[scored & out.direction_v3.eq(direction)].sort_values(
            ['questionnaire_composite', 'statistical_rank'], ascending=[False, True], kind='stable')
        out.loc[ordered.index, 'questionnaire_composite_rank'] = np.arange(1, len(ordered)+1)
    return out, m, M


def record(row):
    fields = list(row.index)
    return {k: None if pd.isna(v) else v.item() if hasattr(v, 'item') else v
            for k, v in row[fields].items()}


def run(verifier, out_dir, final=None, overwrite=False):
    verifier, out_dir = Path(verifier), Path(out_dir)
    final = Path(final) if final else verifier / 'final-2026-10-10'
    inputs = {'ranked': final / 'comparison-ranking/ranked_comparisons_all.csv',
              'module': Path(__file__)}
    for path in inputs.values():
        if not path.exists():
            raise SystemExit(f'Cannot proceed - required data file not found: {path}')
    outputs = {'ranked': out_dir / 'ranked_with_questionnaire_composite.csv'}
    if not overwrite and any(p.exists() for p in [*outputs.values(), out_dir / 'summary.json']):
        raise SystemExit(f'refusing overwrite in {out_dir}; pass --overwrite')
    table, m, M = questionnaire_table(pd.read_csv(inputs['ranked'], low_memory=False))
    counts = {}
    for direction in DIRECTIONS:
        part = table[table.direction_v3.eq(direction)]
        scoped = part[part.condition_category.isin(SCOPE)]
        counts[direction] = {'ranked': len(part), 'scoped': len(scoped),
                            'scored': int(scoped.questionnaire_composite.notna().sum())}
    selectors = {'ouabain': table.gse_id.eq('GSE122080') & table.condition_name.str.contains('ouabain', case=False, na=False),
                 'PLX5622_hippocampus': table.gse_id.eq('GSE134151') & table.comparison_name.str.contains('Hippocampus', case=False, na=False),
                 'cold_exposure': table.gse_id.eq('GSE165974')}
    followed = {}
    for name, selected in selectors.items():
        if not selected.any():
            raise SystemExit(f'Cannot find required comparison: {name}')
        followed[name] = [record(row) for _, row in table[selected].iterrows()]
    summary = {'schema': 'clockbase-questionnaire-composite-v1', 'final': final.name, 'formula': FORMULA,
               'scope': list(SCOPE), 'm': m, 'M': M, 'counts': counts,
               'ranked_comparisons': len(table), 'questionnaire_score': 'final_score',
               'use': 'record only; not used for ranking, filtering or admission',
               'questionnaire_composite_rank': 'within direction among scored, descending composite, ties by statistical_rank',
               'followed_up': followed,
               'inputs': {k: {'path': str(p.resolve()), 'sha256': sha256(p), 'mtime_ns': p.stat().st_mtime_ns}
                          for k, p in inputs.items()},
               'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
        table.to_csv(outputs['ranked'], index=False)
        summary['outputs'] = {k: {'path': str(p.resolve()), 'sha256': sha256(p)} for k, p in outputs.items()}
        (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    except OSError as error:
        raise SystemExit(f'Cannot write questionnaire outputs: {error}') from error
    print(json.dumps({'m': m, 'M': M, 'counts': counts, 'followed_up': followed}, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--verifier', type=Path, default=paths.VERIFIER_ROOT.resolve())
    parser.add_argument('--final', type=Path)
    parser.add_argument('--out', type=Path)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    run(args.verifier, args.out or args.verifier / 'questionnaire-composite-001', args.final, args.overwrite)


if __name__ == '__main__':
    main()
