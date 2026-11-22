"""Judgement ranking by the composite score: the judgement-admitted comparisons ordered by the composite of S17.

The judgement version (followup_judgment_2026_10_18) admits the comparisons whose design two models judged as old at
the start, otherwise healthy and against untreated controls (lower clock outputs), or otherwise healthy and against
untreated controls (higher clock outputs). Its judgement_rank orders the admitted comparisons by the bound of the
adjusted 95% confidence interval of Hedges' g nearest zero. This module orders the same admitted comparisons by the
composite score of composite_ranking_2026_11_08 (Supplementary Table S17):

    composite = min(-log10(max(v3_p, 1e-300)), 30) x relevance_score / 10

where v3_p is the verified P value and relevance_score the blinded rubric score (mean of the two models, 0-10). The
composite is defined for Genetic, Drug and Environment comparisons only; admitted Disease comparisons keep an empty
composite and no composite rank. Within each direction the admitted comparisons with a composite are ranked
1..n by descending composite, ties broken by the statistical rank (judgement_composite_rank). The composite is
recomputed here with composite_ranking_2026_11_08.composite() from the judgement file's v3_p and relevance_score and
must equal the stored composite of ranked_with_composite.csv; the statistical rank, v3_p and relevance score of the
two files must agree as well. Any disagreement stops the run.

Display (Figure 3g, comparisons with lower clock outputs only): one entry per intervention and GEO series, the rule
of dedup_interventions_2026_11_06.py applied to the composite order. The intervention identity (display_key) is that
script's identity(): the canonical target assigned for the replication analysis
(intervention-meta-2026-10-17/primary/meta_input_effects.csv), else the ranking table's canonical target, else the
label, normalized and merged through its synonym table. Walking down the composite order, the first comparison of
each display_key is kept, then the first of each GEO series. display_label is the bar text of
create_figure3_judgement_2026_11_06.py before its panel b-e qualifier and dagger: the label with " at ZT<n>" removed,
and "Overload in senolytic-treated mice" for the comparison "Overload vs Sham (Senolytic)". displayed_bound marks the
rows the same rule keeps under the earlier order (descending |rank_bound|), for comparison.

Inputs:  VERIFIER/followup-judgment-001/judgement_age-lowering.csv, judgement_age-raising.csv
         VERIFIER/composite-ranking-001/ranked_with_composite.csv
         VERIFIER/intervention-meta-2026-10-17/primary/meta_input_effects.csv  (display identity)
Outputs: OUT/judgement_composite_age-lowering.csv, OUT/judgement_composite_age-raising.csv
             every judgement-admitted row with the judgement file's columns (judgement_rank renamed
             judgement_rank_bound) plus neg_log10_p_capped, composite, composite_rank (within direction among all
             scored comparisons, from the composite table), judgement_composite_rank, display_key, display_label,
             displayed_bound and displayed; rows in judgement_composite_rank order, unscored rows last by
             judgement_rank_bound
         OUT/judgement_composite_summary.json   formula, ranking and display rules, counts, displayed rows, rank
                                                changes, input and output sha256

VERIFIER is CLOCKBASE_VERIFIER_ROOT when set, else <repository>/../data/clockbase-verifier-v3-2026-09-27 (the local
revision layout); OUT defaults to VERIFIER/composite-ranking-001.

  python -m audit_system.judgement_composite_2026_11_08 [--verifier DIR] [--out DIR] [--overwrite]
"""
import argparse
import importlib.util
import json
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import paths
from audit_system import composite_ranking_2026_11_08 as cr
from audit_system.membership_recovery_2026_09_27 import sha256

SCHEMA = 'clockbase-judgement-composite-v1'
KEY = cr.KEY
DIRECTIONS = cr.DIRECTIONS  # 'age-lowering' / 'age-raising' -> composite table direction text
DISPLAYED_DIRECTION = 'age-lowering'  # Figure 3g draws the comparisons with lower clock outputs
DEDUP_SCRIPT = paths.REPO / 'scripts/figures_2026_11_06_july_exact/dedup_interventions_2026_11_06.py'
FIGURE_SCRIPT = paths.REPO / 'scripts/figures_2026_11_06_july_exact/create_figure3_judgement_2026_11_06.py'
RENAMED = {'Overload vs Sham (Senolytic)': 'Overload in senolytic-treated mice'}  # as in FIGURE_SCRIPT
LOCAL_VERIFIER = paths.REPO.parent / 'data' / 'clockbase-verifier-v3-2026-09-27'
RTOL = 1e-12


def default_verifier():
    return paths.VERIFIER_ROOT if os.environ.get('CLOCKBASE_VERIFIER_ROOT') else LOCAL_VERIFIER


def read_csv(path):
    try:
        return pd.read_csv(path, low_memory=False)
    except FileNotFoundError:
        raise SystemExit(f'Cannot proceed - required data file not found: {path}')
    except (OSError, pd.errors.ParserError) as error:
        raise SystemExit(f'Cannot read {path}: {error}')


def where(path, verifier):
    """'verifier:relative/path' for a path under the verifier root, 'repo:...' under the repository, else the path."""
    path = Path(path).resolve()
    for label, root in (('verifier', Path(verifier).resolve()), ('repo', paths.REPO)):
        try:
            return f'{label}:{path.relative_to(root).as_posix()}'
        except ValueError:
            continue
    return str(path)


def dedup_identity(verifier):
    """identity() of dedup_interventions_2026_11_06.py, loaded from that script so the display rule has one source."""
    if not DEDUP_SCRIPT.exists():
        raise SystemExit(f'Cannot proceed - required script not found: {DEDUP_SCRIPT}')
    os.environ.setdefault('CLOCKBASE_VERIFIER_ROOT', str(verifier))  # the script resolves its own paths on import
    spec = importlib.util.spec_from_file_location('dedup_interventions_2026_11_06', DEDUP_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except Exception as error:  # noqa: BLE001 - report which script failed and why
        raise SystemExit(f'Cannot load the display rule from {DEDUP_SCRIPT}: {error!r}')
    return module.identity


def join_composite(judged, table, direction, source):
    """Admitted rows of one direction with the composite columns; checks the join and the recomputed composite."""
    if judged.duplicated(KEY).any():
        raise SystemExit(f'{source} has more than one row per comparison')
    if 'judgement_admitted' not in judged or 'judgement_rank' not in judged:
        raise SystemExit(f'{source} lacks judgement_admitted or judgement_rank')
    admitted = judged[judged.judgement_admitted.astype(bool)].copy()
    if not admitted.direction_v3.eq(direction).all():
        raise SystemExit(f'{source} has rows outside direction {direction}')
    checked = ['category', 'direction', 'statistical_rank', 'v3_p', 'relevance_score']  # compared, then dropped
    kept = ['neg_log10_p_capped', 'composite', 'composite_rank']
    part = table[KEY + checked + kept].rename(columns={c: f'ct_{c}' for c in checked})
    merged = admitted.merge(part, on=KEY, how='left', validate='one_to_one', indicator=True)
    missing = merged[merged._merge != 'both']
    if len(missing):
        raise SystemExit(f'{len(missing)} admitted rows of {source} are not in ranked_with_composite.csv: '
                         f'{missing[KEY].to_dict("records")}')
    checks = {
        'direction': merged.ct_direction.eq(DIRECTIONS[direction]),
        'category': merged.ct_category.eq(merged.condition_category),
        'statistical_rank': merged.ct_statistical_rank.astype(int).eq(merged.statistical_rank.astype(int)),
        'v3_p': pd.Series(np.isclose(merged.ct_v3_p, merged.v3_p, rtol=RTOL, atol=0), index=merged.index),
        'relevance_score': pd.Series(np.isclose(merged.ct_relevance_score, merged.relevance_score, rtol=RTOL, atol=0),
                                     index=merged.index),
    }
    recomputed = pd.Series(np.where(merged.condition_category.isin(cr.SCOPE),
                                    cr.composite(merged.v3_p, merged.relevance_score), np.nan), index=merged.index)
    both_empty = recomputed.isna() & merged.composite.isna()
    checks['composite'] = both_empty | pd.Series(np.isclose(recomputed, merged.composite, rtol=RTOL, atol=0),
                                                 index=merged.index)
    checks['composite_rank'] = merged.composite.isna().eq(merged.composite_rank.isna())
    for name, ok in checks.items():
        if not ok.all():
            bad = merged.loc[~ok, KEY + [c for c in merged.columns if name in c]]
            raise SystemExit(f'{source}: {name} disagrees with ranked_with_composite.csv for {(~ok).sum()} rows:\n'
                             f'{bad.to_string()}')
    merged = merged.drop(columns=['_merge'] + [f'ct_{c}' for c in checked])
    merged['composite_rank'] = merged.composite_rank.astype('Int64')  # read as float: empty for Disease rows
    return merged.rename(columns={'judgement_rank': 'judgement_rank_bound'}), int(both_empty.sum())


def rank_by_composite(frame):
    """judgement_composite_rank: 1..n by descending composite, ties by statistical rank; unscored rows get none."""
    out = frame.copy()
    out['judgement_composite_rank'] = pd.Series(pd.NA, index=out.index, dtype='Int64')
    scored = out[out.composite.notna()].sort_values(['composite', 'statistical_rank'], ascending=[False, True],
                                                    kind='stable')
    out.loc[scored.index, 'judgement_composite_rank'] = np.arange(1, len(scored) + 1)
    unscored = out[out.composite.isna()].sort_values('judgement_rank_bound', kind='stable')
    return out.loc[list(scored.index) + list(unscored.index)].reset_index(drop=True)


def first_per_intervention_and_series(ordered):
    """Index of the first row of each display_key, then of each GEO series, walking down the given order."""
    return ordered.drop_duplicates('display_key').drop_duplicates('gse_id').index


def add_display(frame, identity, meta, direction):
    out = identity(frame, meta)
    if len(out) != len(frame) or not (out[KEY].values == frame[KEY].values).all():
        raise SystemExit('the display identity changed the rows or their order')
    out = out.drop(columns=['meta_target'])
    out['display_label'] = [RENAMED.get(c, re.sub(r' at ZT\d+', '', str(label)))
                            for c, label in zip(out.comparison_name, out.label)]
    out['displayed_bound'] = False
    out['displayed'] = False
    if direction == DISPLAYED_DIRECTION:
        by_bound = out.assign(_abs=out.rank_bound.abs()).sort_values('_abs', ascending=False, kind='stable')
        out.loc[first_per_intervention_and_series(by_bound), 'displayed_bound'] = True
        by_composite = out[out.judgement_composite_rank.notna()].sort_values('judgement_composite_rank')
        out.loc[first_per_intervention_and_series(by_composite), 'displayed'] = True
    return out


def row(r):
    def num(v, cast):
        return cast(v) if pd.notna(v) else None
    return {'judgement_composite_rank': num(r.judgement_composite_rank, int),
            'judgement_rank_bound': num(r.judgement_rank_bound, int), 'label': r.label,
            'display_label': r.display_label, 'display_key': r.display_key, 'category': r.condition_category,
            'gse_id': r.gse_id, 'comparison_name': r.comparison_name, 'condition_name': r.condition_name,
            'modality': r.modality, 'composite': num(r.composite, float), 'v3_p': float(r.v3_p),
            'relevance_score': float(r.relevance_score), 'statistical_rank': int(r.statistical_rank),
            'composite_rank': num(r.composite_rank, int), 'abs_rank_bound': float(abs(r.rank_bound)),
            'displayed': bool(r.displayed), 'displayed_bound': bool(r.displayed_bound)}


def summarize(frames, unscored_checked, inputs, verifier, outputs):
    directions = {}
    for direction, frame in frames.items():
        scored = frame[frame.composite.notna()]
        moved = scored[scored.judgement_composite_rank.astype(int) != scored.judgement_rank_bound.astype(int)]
        entry = {
            'admitted': len(frame), 'scored': len(scored), 'unscored': len(frame) - len(scored),
            'unscored_by_category': frame[frame.composite.isna()].condition_category.value_counts().to_dict(),
            'admitted_by_category': frame.condition_category.value_counts().to_dict(),
            'unscored_rows_checked_empty_in_composite_table': unscored_checked[direction],
            'rank_changed_vs_bound': len(moved),
        }
        if direction == DISPLAYED_DIRECTION:
            entry['rows'] = [row(r) for r in frame.itertuples()]
            entry['displayed'] = [row(r) for r in frame[frame.displayed].itertuples()]
            entry['displayed_same_rows_as_bound_rule'] = bool(frame.displayed.eq(frame.displayed_bound).all())
        else:
            entry['top10'] = [row(r) for r in frame[frame.composite.notna()].head(10).itertuples()]
        directions[direction] = entry
    return {
        'schema': SCHEMA,
        'formula': 'composite = min(-log10(max(v3_p, 1e-300)), 30) * relevance_score / 10',
        'composite_scope': list(cr.SCOPE),
        'ranking': 'judgement-admitted comparisons within direction by descending composite; ties broken by the '
                   'statistical rank; comparisons without a composite (Disease) get no judgement_composite_rank',
        'previous_ranking': 'judgement_rank_bound: descending |rank_bound| (CI bound of Hedges\' g nearest zero)',
        'display': {'figure': 'Figure 3g', 'direction': DIRECTIONS[DISPLAYED_DIRECTION],
                    'rule': 'walking down judgement_composite_rank, the first comparison of each display_key '
                            '(identity() of dedup_interventions_2026_11_06.py), then the first of each GEO series',
                    'label': 'display_label = label without " at ZT<n>"; "Overload vs Sham (Senolytic)" -> '
                             '"Overload in senolytic-treated mice"; the figure adds its panel b-e qualifier and '
                             'the dagger of the followed-up comparison'},
        'directions': directions,
        'inputs': {name: {'path': where(p, verifier), 'sha256': sha256(p)} for name, p in inputs.items()},
        'outputs': {name: {'path': where(p, verifier), 'sha256': sha256(p)} for name, p in outputs.items()},
        'module_sha256': sha256(__file__), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}


def run(verifier, out_dir, overwrite=False, judgement=None, composite=None, meta=None):
    verifier, out_dir = Path(verifier), Path(out_dir)
    judgement = Path(judgement or os.environ.get('CLOCKBASE_JUDGEMENT') or verifier / 'followup-judgment-001')
    composite = Path(composite or os.environ.get('CLOCKBASE_COMPOSITE') or verifier / 'composite-ranking-001')
    meta = Path(meta or os.environ.get('CLOCKBASE_META') or verifier / 'intervention-meta-2026-10-17')
    inputs = {
        'judgement_age-lowering': judgement / 'judgement_age-lowering.csv',
        'judgement_age-raising': judgement / 'judgement_age-raising.csv',
        'ranked_with_composite': composite / 'ranked_with_composite.csv',
        'meta_input_effects': meta / 'primary' / 'meta_input_effects.csv',
        'composite_module': Path(cr.__file__), 'dedup_script': DEDUP_SCRIPT, 'figure_script': FIGURE_SCRIPT}
    for p in inputs.values():
        if not Path(p).exists():
            raise SystemExit(f'Cannot proceed - required data file not found: {p}')
    outputs = {d: out_dir / f'judgement_composite_{d}.csv' for d in DIRECTIONS}
    outputs['summary'] = out_dir / 'judgement_composite_summary.json'
    existing = [str(p) for p in outputs.values() if p.exists()]
    if existing and not overwrite:
        raise SystemExit(f'refusing to overwrite {existing} (pass --overwrite)')
    table = read_csv(inputs['ranked_with_composite'])
    if table.duplicated(KEY).any():
        raise SystemExit('ranked_with_composite.csv has more than one row per comparison')
    meta = read_csv(inputs['meta_input_effects'])
    identity = dedup_identity(verifier)
    frames, unscored = {}, {}
    for direction in DIRECTIONS:
        source = inputs[f'judgement_{direction}']
        merged, unscored[direction] = join_composite(read_csv(source), table, direction, source.name)
        frames[direction] = add_display(rank_by_composite(merged), identity, meta, direction)
    out_dir.mkdir(parents=True, exist_ok=True)
    for direction, frame in frames.items():
        frame.to_csv(outputs[direction], index=False)
    summary_path = outputs.pop('summary')
    summary = summarize(frames, unscored, inputs, verifier, outputs)
    summary_path.write_text(json.dumps(summary, indent=2, default=str) + '\n')
    shown = frames[DISPLAYED_DIRECTION]
    print(shown[['judgement_composite_rank', 'judgement_rank_bound', 'label', 'gse_id', 'composite',
                 'statistical_rank', 'composite_rank', 'displayed_bound', 'displayed']].to_string(index=False))
    for direction, entry in summary['directions'].items():
        print(direction, {k: v for k, v in entry.items() if k not in ('rows', 'displayed', 'top10')})


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--verifier', type=Path, default=None,
                        help=f'verifier root (default: CLOCKBASE_VERIFIER_ROOT, else {LOCAL_VERIFIER})')
    parser.add_argument('--out', type=Path, default=None, help='output directory (default: VERIFIER/composite-ranking-001)')
    parser.add_argument('--judgement', type=Path, help='judgement directory (default: CLOCKBASE_JUDGEMENT or verifier default)')
    parser.add_argument('--composite', type=Path, help='composite directory (default: CLOCKBASE_COMPOSITE or verifier default)')
    parser.add_argument('--meta', type=Path, help='meta-analysis directory (default: CLOCKBASE_META or verifier default)')
    parser.add_argument('--overwrite', action='store_true', help='replace earlier outputs of this module')
    args = parser.parse_args(argv)
    verifier = args.verifier or default_verifier()
    run(verifier, args.out or verifier / 'composite-ranking-001', args.overwrite,
        judgement=args.judgement, composite=args.composite, meta=args.meta)


if __name__ == '__main__':
    main()
