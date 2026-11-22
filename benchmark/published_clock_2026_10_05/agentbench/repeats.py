"""Score repeated runs of one configuration (PROTOCOL.md of groupbench, addendum item 27).

Each series was analysed three times by the online agent of the ClockBase website (gpt-6-sol), in separate model
sessions. Every replicate is scored with the functions of items 16 (published results, `score.py`) and 20 (ranked
comparisons, `ranked.py`); the stability endpoints of item 27 (f) compare the replicates.

  python -m benchmark.published_clock_2026_10_05.agentbench.repeats --published RUNS --ranked RUNS --out DIR

RUNS are website run roots (one directory per series with run.json and rep-<k>/out/ result files).
"""
import argparse
import json
import shutil
from collections import Counter
from pathlib import Path

from audit_system import paths
from audit_system.row_repairs_2026_09_28 import individuals_by_comparison

from ..summarize_2026_10_05 import call, judged
from . import ranked as R
from . import score as S
from .prepare import published_records

BENCH = paths.DATA / 'benchmark'
CHECKS = BENCH / 'groupbench/design-checks/final-2026-10-10'
FINAL = 'final-2026-10-10'
REPLICATES = (1, 2, 3)


def website_runs(root):
    """{gse: (run.json record, run directory)} of a website run root."""
    out = {}
    for run_json in sorted(Path(root).glob('*/run.json')):
        record = json.loads(run_json.read_text())
        out[record['gse']] = (record, run_json.parent)
    return out


def stage(run_dir, replicate, target):
    """Copy one replicate's result files to target/results/ (the layout score.read_run and ranked.run_arms read)."""
    source = run_dir / f'rep-{replicate}' / 'out'
    (target / 'results').mkdir(parents=True, exist_ok=True)
    for name in ('aging_analysis.csv', 'comparison_samples.csv'):
        if (source / name).is_file():
            shutil.copyfile(source / name, target / 'results' / name)
    return target


def replicate_status(run_record, replicate):
    for entry in run_record.get('replicate_results') or []:
        if entry.get('replicate') == replicate:
            return entry.get('status')
    return None


def status(record, system):
    """Item 27 (f): not covered, untested, or the call (higher, lower, none) of one replicate."""
    j = judged(record, system)
    if j is None:
        return 'NOT_COVERED'
    c = call(j['direction'], j['significant'])
    return c if c is not None else 'UNTESTED'


def published(root, out):
    runs = website_runs(root)
    series = sorted(runs)
    records = [r for r in published_records(BENCH) if r['gse'] in set(series)]
    checks = {(r['gse'], r['row_ordinal']): r
              for r in map(json.loads, (CHECKS / 'dnam/rows.jsonl').read_text().splitlines())}
    reader, cache = S.PredictionReader(), {}
    replicate_runs = {}
    for k in REPLICATES:
        replicate_runs[k] = {gse: S.read_run(stage(runs[gse][1], k, out / 'staged' / 'published' / f'rep-{k}' / gse))
                             for gse in series}
    for record in records:
        record['systems']['S2c'] = S.released_result(record, checks)
        for k in REPLICATES:
            run = replicate_runs[k][record['gse']]
            raw_key, rec_key = f'W{k}:raw', f'W{k}:recomputed'
            found = S.match(record, run)
            if found is None:
                state = 'NO_RESULTS' if not run['rows'] else 'NO_MATCHING_COMPARISON'
                record['systems'][raw_key] = record['systems'][rec_key] = {'covered': False, 'status': state}
                continue
            row, swapped, overlap = found
            arms = run['arms'][(row['comparison'], row['condition'])]
            record['systems'][raw_key] = S.raw_result(row, swapped, overlap)
            record['systems'][rec_key] = S.recomputed_result(record, row, arms, swapped, overlap, reader, cache)
    statuses = [[status(r, f'W{k}:recomputed') for k in REPLICATES] for r in records]
    covered_all = [s for s in statuses if all(x not in ('NOT_COVERED',) for x in s)]
    tested_all = [s for s in covered_all if 'UNTESTED' not in s]
    same_tested = [s for s in tested_all if len(set(s)) == 1]
    per_replicate = {f'W{k}:{kind}': S.describe(records, f'W{k}:{kind}') for k in REPLICATES for kind in ('raw', 'recomputed')}
    for system in ('S1', 'S2c', 'S3'):
        per_replicate[system] = S.describe(records, system)
    patterns = Counter('/'.join(s) for s in statuses)
    return {
        'series': len(series), 'records': len(records),
        'replicates_completed': {k: sum(replicate_status(runs[g][0], k) == 'COMPLETED' for g in series) for k in REPLICATES},
        'primary': {'covered_by_all_three': len(covered_all), 'tested_by_all_three': len(tested_all),
                    'same_call_in_all_three': len(same_tested),
                    'fraction': len(same_tested) / len(tested_all) if tested_all else None},
        'secondary': {'same_status_in_all_three': sum(len(set(s)) == 1 for s in statuses), 'of': len(statuses),
                      'fraction': sum(len(set(s)) == 1 for s in statuses) / len(statuses) if statuses else None},
        'status_patterns': dict(patterns.most_common()),
        'per_replicate': per_replicate,
        'cost_usd': round(sum(runs[g][0].get('cost_usd') or 0 for g in series), 4),
    }, records


def ranked_scores(root, out):
    runs = website_runs(root)
    targets = R.ranked(FINAL)
    read_individuals = individuals_by_comparison(paths.VERIFIER_ROOT / 'individuals-audit-001')
    cache, series_cache = {}, {}
    check = []
    for target in targets:
        values = R.clock_values(target['gse'], cache)
        reference = R.statistic(target, target['case'], target['control'], values,
                                R.target_series(target, series_cache, read_individuals))
        released = target['released']
        check.append(reference['p'] is not None and abs(reference['p'] - released['p']) <= 1e-9 * max(1, released['p'])
                     and abs(reference['mean_diff'] - released['mean_diff']) <= 1e-9 * max(1, abs(released['mean_diff'])))
    records = []
    for k in REPLICATES:
        arms_by_series = {gse: R.run_arms(stage(runs[gse][1], k, out / 'staged' / 'ranked' / f'rep-{k}' / gse))
                          for gse in runs}
        for target in targets:
            expected = 1 if target['table'] == 'higher' else -1
            record = {'replicate': k, 'table': target['table'], 'rank': target['rank'], 'gse': target['gse'],
                      'comparison': target['comparison'], 'found': False, 'same_direction': False, 'recurs': False,
                      'opposite': False}
            arms = arms_by_series.get(target['gse'], {'arms': {}})
            match = R.best_match(target, arms['arms'])
            if match is not None:
                (overlap, _, _), key, swapped, case, control = match
                result = R.statistic(target, case, control, R.clock_values(target['gse'], cache),
                                     R.target_series(target, series_cache, read_individuals))
                diff = result['mean_diff']
                record.update({'found': True, 'jaccard_min': overlap, 'swapped': swapped, 'run_comparison': key[0],
                               'method': result['method'], 'mean_diff': diff, 'p': result['p'],
                               'n_case': result['n_case'], 'n_control': result['n_control']})
                record['same_direction'] = diff is not None and diff * expected > 0
                record['opposite'] = diff is not None and diff * expected < 0
                record['recurs'] = record['same_direction'] and result['p'] is not None and result['p'] < 0.05
            records.append(record)
    by_target = {}
    for record in records:
        by_target.setdefault((record['table'], record['rank']), []).append(record)
    recur_counts = Counter(sum(r['recurs'] for r in group) for group in by_target.values())
    per_replicate = {}
    for k in REPLICATES:
        part = [r for r in records if r['replicate'] == k]
        per_replicate[k] = {table: {key: sum(r[key] for r in part if r['table'] == table)
                                    for key in ('found', 'same_direction', 'recurs', 'opposite')}
                            for table in ('higher', 'lower')}
    return {
        'comparisons': len(targets), 'series': len(runs),
        'recomputation_check_equal': sum(check),
        'replicates_completed': {k: sum(replicate_status(runs[g][0], k) == 'COMPLETED' for g in runs) for k in REPLICATES},
        'per_replicate': per_replicate,
        'recurring_in_all_three': recur_counts.get(3, 0),
        'recurring_in_at_least_two': recur_counts.get(3, 0) + recur_counts.get(2, 0),
        'recurring_in_at_least_one': sum(v for n, v in recur_counts.items() if n >= 1),
        'recurrence_distribution': {n: recur_counts.get(n, 0) for n in (0, 1, 2, 3)},
        'found_with_opposite_sign_in_any': sum(any(r['opposite'] for r in g) for g in by_target.values()),
        'cost_usd': round(sum(runs[g][0].get('cost_usd') or 0 for g in runs), 4),
    }, records


def website_labels(root):
    labels = Counter()
    for gse, (record, run_dir) in website_runs(root).items():
        analysis = run_dir / 'analysis.json'
        labels[json.loads(analysis.read_text())['stability']['label'] if analysis.is_file() else 'NO_ANALYSIS'] += 1
    return dict(labels)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--published', type=Path, required=True)
    parser.add_argument('--ranked', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    pub, pub_records = published(args.published, args.out)
    rank, rank_records = ranked_scores(args.ranked, args.out)
    summary = {'item': 27, 'published': pub, 'ranked': rank,
               'website_labels': {'published': website_labels(args.published), 'ranked': website_labels(args.ranked)}}
    (args.out / 'repeats.json').write_text(json.dumps(summary, indent=1, default=str))
    with (args.out / 'published_records.jsonl').open('w') as handle:
        for record in pub_records:
            systems = {s: record['systems'].get(s) for s in ('S1', 'S2c') + tuple(f'W{k}:{t}' for k in REPLICATES
                                                                                for t in ('raw', 'recomputed'))}
            handle.write(json.dumps({k: record.get(k) for k in ('gse', 'clock_clockbase', 'direction', 'significant')}
                                    | {'systems': systems}, default=str) + '\n')
    with (args.out / 'ranked_records.jsonl').open('w') as handle:
        for record in rank_records:
            handle.write(json.dumps(record, default=str) + '\n')
    print(json.dumps({k: v for k, v in summary.items()}, indent=1, default=str)[:6000])


if __name__ == '__main__':
    main()
