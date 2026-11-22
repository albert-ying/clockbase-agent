"""Composite version of the comparison ranking: the historical composite formula on the verified comparisons.

The comparison ranking (comparison_ranking_2026_09_30) orders verified comparisons by the bound of the adjusted 95%
confidence interval of Hedges' g nearest zero. The submitted tables were ordered by a composite of significance and
the agent score (generate_tables_1_2.py; create_combined_modality_top_hits_v2.py):
    min(-log10(clip(P, 1e-300)), 30) x scaled agent score.
This module gives that composite for the verified comparisons, with the verified P value (v3_p) for P and the blinded
relevance score (mean of the two models' rubric scores, 0-10; intervention_classification_2026_09_29) divided by 10
for the scaled agent score:

    composite = min(-log10(max(v3_p, 1e-300)), 30) x relevance_score / 10

Scope: comparisons of the Genetic, Drug and Environment categories. The rubric scores interventions, so
non-interventions score 0-2 by design; Disease comparisons get no composite (empty) and keep the statistical ranking
only. Within each direction (lower and higher clock outputs, direction_v3) the scored comparisons are ranked by
descending composite, ties broken by the statistical rank. The statistical rank is the position of a comparison among
the comparisons of its direction in ranked_comparisons_all.csv (as in followup_judgment_2026_10_18); the eligible set,
the sign convention and the statistical ranking are not changed.

Input:  FINAL/comparison-ranking/ranked_comparisons_all.csv
Output: DIR/ranked_with_composite.csv  every ranked comparison in the order of the input file (keys modality, gse_id,
                                       comparison_name, condition_name; category, direction, statistical rank, v3_p,
                                       relevance_score, composite, composite_rank)
        DIR/summary.json               counts; the top 15 comparisons with lower clock outputs; the composite ranks of
                                       ouabain (GSE122080) and PLX5622 (GSE134151); Spearman correlation of the
                                       composite and statistical ranks; overlap with Table 2; input sha256

  python -m audit_system.composite_ranking_2026_11_08 --final final-2026-10-10 --out DIR
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import sha256

SCHEMA = 'clockbase-composite-ranking-v1'
KEY = ['modality', 'gse_id', 'comparison_name', 'condition_name']
P_FLOOR, CAP = 1e-300, 30
SCOPE = ('Genetic', 'Drug', 'Environment')
DIRECTIONS = {'age-lowering': 'lower clock outputs', 'age-raising': 'higher clock outputs'}
FOLLOWED_UP = {'ouabain': ('GSE122080', 'ouabain'), 'PLX5622': ('GSE134151', None)}  # None: every ranked row
TABLE_N, TOP_N = 30, 15


def composite(p, relevance):
    """min(-log10(max(P, 1e-300)), 30) x relevance / 10."""
    return np.minimum(-np.log10(np.maximum(p, P_FLOOR)), CAP) * relevance / 10


def ranked_rows(final):
    """The ranked comparisons with their direction and statistical rank (within direction and overall)."""
    source = Path(final) / 'comparison-ranking' / 'ranked_comparisons_all.csv'
    try:
        table = pd.read_csv(source, low_memory=False)
    except OSError:
        raise SystemExit(f'Cannot proceed - required data file not found: {source}')
    if table.duplicated(KEY).any():
        raise SystemExit('ranked_comparisons_all.csv has more than one row per comparison')
    unknown = set(table.direction_v3) - set(DIRECTIONS)
    if unknown or table.v3_p.isna().any() or table.relevance_score.isna().any():
        raise SystemExit(f'unexpected ranked rows: directions {unknown}, missing P or relevance score')
    table['statistical_rank_overall'] = np.arange(1, len(table) + 1)
    table['statistical_rank'] = table.groupby('direction_v3').cumcount() + 1
    return table, source


def composite_table(table):
    """Composite and composite rank (within direction, ties by statistical rank) for the comparisons in scope."""
    out = table.copy()
    scored = out.condition_category.isin(SCOPE)
    out['neg_log10_p_capped'] = np.where(scored, np.minimum(-np.log10(np.maximum(out.v3_p, P_FLOOR)), CAP), np.nan)
    out['composite'] = np.where(scored, composite(out.v3_p, out.relevance_score), np.nan)
    out['composite_rank'] = pd.Series(pd.NA, index=out.index, dtype='Int64')
    for direction in DIRECTIONS:
        part = out[scored & (out.direction_v3 == direction)].sort_values(['composite', 'statistical_rank'],
                                                                            ascending=[False, True])
        out.loc[part.index, 'composite_rank'] = np.arange(1, len(part) + 1)
    return out


def spearman(frame):
    result = stats.spearmanr(frame.composite_rank.astype(float), frame.statistical_rank)
    return {'n': len(frame), 'rho': float(result.statistic), 'p': float(result.pvalue)}


def summarize(table, source, final):
    scored = table[table.composite.notna()]
    lower = scored[scored.direction_v3 == 'age-lowering']
    higher = scored[scored.direction_v3 == 'age-raising']
    lower_all = table[table.direction_v3 == 'age-lowering']
    table2 = lower_all[lower_all.statistical_rank <= TABLE_N]
    table2_scored = table2[table2.composite.notna()]
    ties = int(sum(part.composite.duplicated(keep=False).sum() for _, part in scored.groupby('direction_v3')))

    def row(r):
        return {'gse_id': r.gse_id, 'label': r.label, 'category': r.condition_category, 'modality': r.modality,
                'comparison_name': r.comparison_name, 'direction': DIRECTIONS[r.direction_v3],
                'statistical_rank': int(r.statistical_rank),
                'composite_rank': int(r.composite_rank) if pd.notna(r.composite_rank) else None,
                'composite': float(r.composite) if pd.notna(r.composite) else None,
                'v3_p': float(r.v3_p), 'relevance_score': float(r.relevance_score)}

    followed = {}
    for name, (gse, term) in FOLLOWED_UP.items():
        rows = table[table.gse_id == gse]
        if term:
            rows = rows[rows.label.str.contains(term, case=False, na=False)
                        | rows.condition_name.str.contains(term, case=False, na=False)]
        if rows.empty:
            raise SystemExit(f'no ranked row for {name} ({gse})')
        followed[name] = [row(r) for r in rows.itertuples()]
    # Single order across both directions (the S17 Rank against the composite over all scored comparisons).
    single = scored.sort_values(['composite', 'statistical_rank_overall'], ascending=[False, True])
    single_rho = stats.spearmanr(np.arange(1, len(single) + 1), single.statistical_rank_overall)
    return {
        'schema': SCHEMA, 'final': Path(final).name,
        'formula': 'composite = min(-log10(max(v3_p, 1e-300)), 30) * relevance_score / 10',
        'scope': list(SCOPE), 'not_scored': sorted(set(table.condition_category) - set(SCOPE)),
        'ranking': 'within direction by descending composite; ties broken by the statistical rank',
        'ranked_comparisons': len(table),
        'n_scored': len(scored), 'n_scored_lower': len(lower), 'n_scored_higher': len(higher),
        'scored_by_category': {d: part.condition_category.value_counts().to_dict()
                               for d, part in ((DIRECTIONS[k], scored[scored.direction_v3 == k]) for k in DIRECTIONS)},
        'composite_values_tied_within_direction': ties,
        'lower_top15': [row(r) for r in lower.sort_values('composite_rank').head(TOP_N).itertuples()],
        'followed_up': followed,
        'spearman_composite_vs_statistical_rank': {
            'scored_lower': spearman(lower),
            'scored_higher': spearman(higher),
            'all_scored': {**spearman(scored), 'ranks': 'composite_rank and statistical_rank, both within direction'},
            'all_scored_single_order': {'n': len(single), 'rho': float(single_rho.statistic),
                                        'p': float(single_rho.pvalue),
                                        'ranks': 'composite over all scored comparisons (ties by statistical rank) '
                                                 'against the overall statistical rank (S17 Rank)'}},
        'table2': {'comparisons': len(table2), 'scored': len(table2_scored),
                   f'scored_in_composite_top{TABLE_N}_lower': int((table2_scored.composite_rank <= TABLE_N).sum()),
                   'definition': f'top {TABLE_N} comparisons with lower clock outputs by statistical rank'},
        'input': paths.portable(source), 'input_sha256': sha256(source),
        'module_sha256': sha256(__file__), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}


def run(final, out_dir):
    out_dir = Path(out_dir)
    if out_dir.exists():
        raise SystemExit(f'refusing to overwrite {out_dir}')
    table, source = ranked_rows(final)
    table = composite_table(table)
    summary = summarize(table, source, final)
    columns = (KEY + ['label', 'condition_category', 'direction_v3', 'statistical_rank', 'statistical_rank_overall',
                      'v3_p', 'neg_log10_p_capped', 'relevance_score', 'composite', 'composite_rank'])
    out = table[columns].rename(columns={'condition_category': 'category', 'direction_v3': 'direction'})
    out['direction'] = out.direction.map(DIRECTIONS)
    out_dir.mkdir(parents=True)
    out.to_csv(out_dir / 'ranked_with_composite.csv', index=False)
    summary['output_sha256'] = sha256(out_dir / 'ranked_with_composite.csv')
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str) + '\n')
    print(json.dumps(summary, indent=1, default=str))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--final', required=True, help='verifier final run (name under the verifier root, or a path)')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    final = Path(args.final) if Path(args.final).is_dir() else paths.VERIFIER_ROOT / args.final
    run(final, args.out)


if __name__ == '__main__':
    main()
