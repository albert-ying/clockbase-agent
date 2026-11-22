"""Agent reruns of the ranked comparisons (PROTOCOL.md of groupbench, addendum item 20).

`prepare` writes one run directory per series of Tables 1 and 2 and model, as prepare.py does for item 16, with the
task text of the released mouse template `all_samples_local` and the series' study design from its public GEO record.
`score` matches each ranked comparison to the runs' comparisons by their arms and recomputes the verifier's statistic
on the run's arms.

  python -m benchmark.published_clock_2026_10_05.agentbench.ranked prepare --final NAME --root RUNS --python PYTHON \
      --geo-cache DIR
  python -m benchmark.published_clock_2026_10_05.agentbench.ranked score --final NAME --root RUNS --out DIR
"""
import argparse
import csv
import json
import re
import shutil
import time
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

from audit_system import design_estimand_2026_10_10 as design
from audit_system import paths
from audit_system.row_repairs_2026_09_28 import MIN_ARM, individuals_by_comparison, recompute
from src.parallel_mm_rna import TASK_TEMPLATES

from .prepare import EXTRA, MODELS, NOTE
from .score import audit_run, jaccard, number, run_dir

CLOCK = 'Hazard_Mouse_All_All_scaleddiff_EN_tAge'
TEMPLATE = 'all_samples_local'
JACCARD = 0.5
GEO_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=self&form=text&view=brief'


def ranked(final):
    """The 60 comparisons of Tables 1 and 2 with their verified arms."""
    root = paths.VERIFIER_ROOT / final
    rows = []
    for table, label in (('table1_pro_aging.csv', 'higher'), ('table2_anti_aging.csv', 'lower')):
        frame = pd.read_csv(root / 'comparison-ranking' / table)
        for rank, row in enumerate(frame.itertuples(), start=1):
            record = json.loads((root / 'public-export/studies/mm_rna' / f'{row.gse_id}.json').read_text())
            hits = [r for r in record['rows'] if r['comparison_name'] == row.comparison_name
                    and r['condition_name'] == row.condition_name and r['clock'] == row.clock]
            if len(hits) > 1:  # repeated labels: the record whose released statistic is the table row's
                hits = [r for r in hits if r.get('duplicate_of') is None and (r.get('recomputed') or {}).get('p') is not None
                        and abs(r['recomputed']['p'] - row.v3_p) <= 1e-12 * max(1, row.v3_p)
                        and r['recomputed']['n_case'] == row.v3_n_case and r['recomputed']['n_control'] == row.v3_n_control]
            if len(hits) != 1:
                raise SystemExit(f'{row.gse_id} {row.comparison_name}: {len(hits)} records')
            rows.append({'table': label, 'rank': rank, 'gse': row.gse_id, 'comparison': row.comparison_name,
                         'condition': row.condition_name, 'clock': row.clock, 'hedges_g': float(row.hedges_g),
                         'p': float(row.v3_p), 'case': hits[0]['case_gsms'], 'control': hits[0]['control_gsms'],
                         'released': hits[0]['recomputed']})
    return rows


def study_design(gse, cache):
    """Title, summary and overall design of a series from its public GEO record (cached as text)."""
    path = Path(cache) / f'{gse}.txt'
    if not path.exists():
        for attempt in range(4):
            try:
                with urllib.request.urlopen(GEO_URL.format(gse), timeout=40) as response:
                    text = response.read().decode('utf-8', 'replace')
                if '!Series_title' in text:
                    path.write_text(text)
                    break
            except Exception:  # retried, then reported by the caller
                pass
            time.sleep(2 + 2 * attempt)
    text = path.read_text() if path.exists() else ''

    def field(name):
        return ' '.join(v.strip() for v in re.findall(rf'^!Series_{name} = (.*)$', text, flags=re.M))
    if not field('title'):
        raise ValueError(f'{gse}: GEO record not read')
    return {'title': field('title'), 'summary': field('summary'), 'overall_design': field('overall_design')}


def task_text(gse, design_block, run, python):
    text = TASK_TEMPLATES[TEMPLATE].format(study_id=gse, study_design=design_block, data_dir=f'{run}/data',
                                           RESULTS_DIR=f'{run}/results', OUTPUT_DIR=f'{run}/results')
    text = text.replace(f'{run}/results/{gse}/', f'{run}/results/')
    if '/dfs/' in text or str(paths.REPO) in text:
        raise ValueError(f'{gse}: server path left in the task text')
    return NOTE.format(run=run, python=python) + text.strip('\n') + '\n' + EXTRA.format(run=run)


def write_inputs(gse, run):
    overview = paths.GEO / 'gse_rna_mm.csv'
    metadata = paths.GEO / 'GEO_mm_rna/processed_metadata_preharmonize' / f'{gse}.csv'
    clocks = paths.VERIFIER_ROOT / 'replay-mm_rna-001/inputs' / f'{gse}.csv'
    for path in (overview, metadata, clocks):
        if not path.is_file():
            raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}')
    (run / 'data/dataset_level').mkdir(parents=True)
    (run / 'data/GEO_mm_rna/processed_metadata_preharmonize').mkdir(parents=True)
    (run / 'data/clock_output').mkdir(parents=True)
    (run / 'results').mkdir()
    with overview.open(newline='') as source, (run / 'data/dataset_level/gse_rna_mm.csv').open('w', newline='') as sink:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(sink, fieldnames=reader.fieldnames)
        writer.writeheader()
        rows = [row for row in reader if row.get('accession') == gse]
        writer.writerows(rows)
    if not rows:
        raise ValueError(f'{gse}: no row in the dataset overview')
    shutil.copyfile(metadata, run / 'data/GEO_mm_rna/processed_metadata_preharmonize' / f'{gse}.csv')
    shutil.copyfile(clocks, run / 'data/clock_output/mouse_rna_age_short.csv')
    return {'overview': str(overview), 'metadata': str(metadata), 'clocks': str(clocks)}


def prepare(args):
    rows = ranked(args.final)
    series = sorted({r['gse'] for r in rows}, key=lambda g: int(g[3:]))
    Path(args.geo_cache).mkdir(parents=True, exist_ok=True)
    args.root.mkdir(parents=True, exist_ok=True)
    plan = {'final': args.final, 'comparisons': len(rows), 'series': series, 'models': args.models,
            'framework': args.framework,
            'template': TEMPLATE, 'python': args.python, 'runs': []}
    for gse in series:
        try:
            block = study_design(gse, args.geo_cache)
        except ValueError as error:
            plan['runs'].append({'gse': gse, 'status': 'NOT_PREPARED', 'error': str(error)})
            continue
        for model in args.models:
            run = run_dir(args.root, args.framework, model, gse).resolve()
            if run.exists():
                continue
            try:
                sources = write_inputs(gse, run)
                (run / 'TASK.md').write_text(task_text(gse, block, run, args.python))
            except (OSError, ValueError) as error:
                shutil.rmtree(run, ignore_errors=True)
                plan['runs'].append({'model': model, 'gse': gse, 'status': 'NOT_PREPARED', 'error': str(error)})
                continue
            plan['runs'].append({'model': model, 'gse': gse, 'status': 'PREPARED', 'run': str(run), 'sources': sources})
    (args.root / 'plan.json').write_text(json.dumps(plan, indent=1) + '\n')
    print(len(series), 'series;', sum(r['status'] == 'PREPARED' for r in plan['runs']), 'runs prepared;',
          [r for r in plan['runs'] if r['status'] != 'PREPARED'][:5])


def run_arms(run):
    """The arms of every comparison of one run, and its reported rows on the ranked clock."""
    out = {'arms': {}, 'reported': {}, 'problems': []}
    samples, table = run / 'results/comparison_samples.csv', run / 'results/aging_analysis.csv'
    if not samples.is_file():
        out['problems'].append('NO_SAMPLES_FILE')
        return out
    try:
        with samples.open(newline='', encoding='utf-8-sig') as handle:
            for row in csv.DictReader(handle):
                key = ((row.get('comparison_name') or '').strip(), (row.get('condition_name') or '').strip())
                arm, gsm = (row.get('arm') or '').strip().lower(), (row.get('gsm_id') or '').strip()
                if arm in ('case', 'control') and gsm:
                    out['arms'].setdefault(key, {'case': [], 'control': []})[arm].append(gsm)
        if table.is_file():
            with table.open(newline='', encoding='utf-8-sig') as handle:
                for row in csv.DictReader(handle):
                    if (row.get('clock') or '').strip() == CLOCK:
                        key = ((row.get('comparison_name') or '').strip(), (row.get('condition_name') or '').strip())
                        out['reported'][key] = {'p': number(row.get('p_value')), 'effect': number(row.get('effect_size')),
                                                'direction': number(row.get('direction'))}
        else:
            out['problems'].append('NO_RESULTS_TABLE')
    except (OSError, ValueError, csv.Error) as error:
        out['problems'].append(f'UNREADABLE: {error}')
    return out


def best_match(target, arms):
    best = None
    for key, arm in arms.items():
        for swapped in (False, True):
            case, control = (arm['control'], arm['case']) if swapped else (arm['case'], arm['control'])
            a, b = jaccard(target['case'], case), jaccard(target['control'], control)
            rank = (min(a, b), a + b, not swapped)
            if min(a, b) >= JACCARD and (best is None or rank > best[0]):
                best = (rank, key, swapped, case, control)
    return best


def clock_values(gse, cache):
    if gse not in cache:
        table = pd.read_csv(paths.VERIFIER_ROOT / 'replay-mm_rna-001/inputs' / f'{gse}.csv', low_memory=False)
        values = pd.to_numeric(table[CLOCK], errors='coerce')
        values.index = table['gsm_id'].astype(str)
        cache[gse] = values[~values.index.duplicated(keep=False)]
    return cache[gse]


def target_series(target, series_cache, read_individuals):
    """The design inputs of the target's series; where GEO has no individual field, the two readers' individuals of
    the ranked comparison (the verifier's rule for ranking candidates)."""
    gse = target['gse']
    if gse not in series_cache:
        series_cache[gse] = design.load_series(paths.GEO, 'GEO_mm_rna', gse, False)
    series = series_cache[gse]
    key = (gse, target['comparison'], target['condition'])
    if not series['subjects'] and key in read_individuals:
        series = {**series, 'subjects': read_individuals[key], 'subject_field': None,
                  'subject_source': 'two-reader individuals audit'}
    return series


def statistic(target, case, control, values, series):
    dropped = set(design.overflow(values, list(case) + list(control)))
    case, control = [g for g in case if g not in dropped], [g for g in control if g not in dropped]
    estimate = design.estimate(values, case, control, series, 'mm_rna', CLOCK, (target['comparison'], target['condition']), ())
    if estimate is not None and estimate['recomputed'] is not None:
        result, method = estimate['recomputed'], estimate['design']['method']
    else:
        result, method = recompute(values, case, control), 'WELCH'
    tested = (result.get('mean_diff') is not None and min(result['n_case'], result['n_control']) >= MIN_ARM
              and result.get('p') is not None)
    return {'method': method, 'mean_diff': result.get('mean_diff'), 'p': result.get('p') if tested else None,
            'n_case': result['n_case'], 'n_control': result['n_control']}


def score(args):
    targets = ranked(args.final)
    plan = json.loads((args.root / 'plan.json').read_text())
    python = plan['python']
    framework = plan.get('framework', 'codex')
    cache, series_cache, records, runs = {}, {}, [], []
    read_individuals = individuals_by_comparison(paths.VERIFIER_ROOT / 'individuals-audit-001')
    checks = []
    for target in targets:
        values = clock_values(target['gse'], cache)
        reference = statistic(target, target['case'], target['control'], values,
                              target_series(target, series_cache, read_individuals))
        released = target['released']
        checks.append({'gse': target['gse'], 'comparison': target['comparison'],
                       'equal': reference['p'] is not None and abs(reference['p'] - released['p']) <= 1e-9 * max(1, released['p'])
                       and abs(reference['mean_diff'] - released['mean_diff']) <= 1e-9 * max(1, abs(released['mean_diff']))})
    for model in plan['models']:
        for gse in plan['series']:
            run = run_dir(args.root, framework, model, gse)
            info = audit_run(framework, run, python)
            arms = run_arms(run)
            # Web access excludes a run as in the benchmark scorer; before item 24 no ranked run had any.
            excluded = bool(info['outside_paths'] or info['web_searches']) or info['status'] is None
            runs.append({'model': model, 'gse': gse, **{k: info[k] for k in ('status', 'commands', 'seconds',
                         'web_searches', 'network_attempts')}, 'outside_paths': len(info['outside_paths']),
                         'comparisons': len(arms['arms']), 'problems': ';'.join(arms['problems']), 'excluded': excluded,
                         'input_tokens': (info['usage'] or {}).get('input_tokens'),
                         'cached_input_tokens': (info['usage'] or {}).get('cached_input_tokens'),
                         'output_tokens': (info['usage'] or {}).get('output_tokens')})
            for target in (t for t in targets if t['gse'] == gse):
                expected = 1 if target['table'] == 'higher' else -1
                record = {'model': model, 'table': target['table'], 'rank': target['rank'], 'gse': gse,
                          'comparison': target['comparison'], 'found': False, 'same_direction': False, 'recurs': False}
                match = None if excluded else best_match(target, arms['arms'])
                if match is not None:
                    (overlap, _, _), key, swapped, case, control = match
                    result = statistic(target, case, control, clock_values(gse, cache),
                                       target_series(target, series_cache, read_individuals))
                    diff = result['mean_diff']
                    record.update({'found': True, 'jaccard_min': overlap, 'swapped': swapped, 'run_comparison': key[0],
                                   'method': result['method'], 'mean_diff': diff, 'p': result['p'],
                                   'n_case': result['n_case'], 'n_control': result['n_control']})
                    record['same_direction'] = diff is not None and diff * expected > 0
                    record['recurs'] = record['same_direction'] and result['p'] is not None and result['p'] < 0.05
                    reported = arms['reported'].get(key)
                    if reported:
                        sign = reported['direction'] if reported['direction'] in (1.0, -1.0) else reported['effect']
                        record['reported_p'] = reported['p']
                        record['reported_same_direction'] = (sign is not None and sign != 0
                                                             and ((sign > 0) != swapped) == (expected > 0))
                records.append(record)
    frame = pd.DataFrame(records)
    summary = {'final': args.final, 'comparisons': len(targets), 'series': len(plan['series']), 'models': plan['models'],
               'recomputation_check': {'rows': len(checks), 'equal': sum(c['equal'] for c in checks),
                                       'unequal': [c for c in checks if not c['equal']]},
               'runs': {'planned': len(runs), 'status_0': sum(r['status'] == '0' for r in runs),
                        'timed_out': sum(r['status'] == '124' for r in runs), 'excluded': sum(r['excluded'] for r in runs),
                        'no_comparison': sum(r['comparisons'] == 0 for r in runs),
                        'web_searches': sum(r['web_searches'] for r in runs),
                        'input_tokens': int(sum(r['input_tokens'] or 0 for r in runs)),
                        'cached_input_tokens': int(sum(r['cached_input_tokens'] or 0 for r in runs)),
                        'output_tokens': int(sum(r['output_tokens'] or 0 for r in runs)),
                        'median_seconds': float(pd.Series([r['seconds'] for r in runs if r['seconds']]).median())},
               'by_model': {}, 'across_models': {}}
    for model, part in frame.groupby('model'):
        summary['by_model'][model] = {table: {key: int(part[part.table == table][key].sum())
                                              for key in ('found', 'same_direction', 'recurs')}
                                      for table in ('higher', 'lower')}
    per = frame.groupby(['table', 'gse', 'comparison']).agg(found=('found', 'sum'), recurs=('recurs', 'sum'))
    for table in ('higher', 'lower'):
        part = per.loc[table]
        summary['across_models'][table] = {'found_by_any': int((part.found >= 1).sum()),
                                           'recurs_in_any': int((part.recurs >= 1).sum()),
                                           'recurs_in_two_or_more': int((part.recurs >= 2).sum()),
                                           'recurs_in_all': int((part.recurs == len(plan['models'])).sum())}
    args.out.mkdir(parents=True, exist_ok=False)
    frame.to_csv(args.out / 'ranked_reruns.csv', index=False)
    pd.DataFrame(runs).to_csv(args.out / 'runs.csv', index=False)
    (args.out / 'ranked_reruns.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps({k: summary[k] for k in ('recomputation_check', 'runs', 'by_model', 'across_models')}, indent=1))


def combine(args):
    """Ranked comparisons across the configurations of several scored runs (PROTOCOL.md items 20, 22 and 24).

    A configuration is a scored directory and a model, since a model can run under several frameworks (gpt-6-astra in
    Codex and in Biomni)."""
    tables = []
    for directory in args.scored:
        table = pd.read_csv(Path(directory) / 'ranked_reruns.csv')
        table['configuration'] = Path(directory).name + ':' + table.model.astype(str)
        tables.append(table)
    rows = pd.concat(tables, ignore_index=True)
    rows['recurs'] = rows.recurs.fillna(False).astype(bool)
    rows['found'] = rows.found.fillna(False).astype(bool)
    rows['reversed'] = rows.found & ~rows.same_direction.fillna(True).astype(bool)
    key = ['table', 'rank', 'gse', 'comparison']
    per = rows.groupby(key).agg(recurs=('recurs', 'sum'), reversed=('reversed', 'sum')).reset_index()
    configurations = sorted(rows.configuration.unique())
    out = {'scored': [str(d) for d in args.scored], 'models': sorted(rows.model.unique()),
           'configurations': configurations, 'comparisons': len(per)}
    for table in ('higher', 'lower'):
        part = per[per.table == table]
        out[table] = {'recurs_in_any': int((part.recurs >= 1).sum()), 'recurs_in_two_or_more': int((part.recurs >= 2).sum()),
                      'recurs_in_all': int((part.recurs == len(configurations)).sum()), 'reversed': int((part.reversed > 0).sum())}
    out['recurs_in_two_or_more'] = out['higher']['recurs_in_two_or_more'] + out['lower']['recurs_in_two_or_more']
    out['recurs_in_all'] = out['higher']['recurs_in_all'] + out['lower']['recurs_in_all']
    Path(args.out).write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    c = sub.add_parser('combine')
    c.add_argument('--scored', type=Path, nargs='+', required=True, help='scored directories (ranked_reruns.csv each)')
    c.add_argument('--out', type=Path, required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--final', required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--python', required=True)
    p.add_argument('--geo-cache', required=True)
    p.add_argument('--models', nargs='*', default=list(MODELS))
    p.add_argument('--framework', choices=['codex', 'claude', 'biomni'], default='codex',
                   help='agent framework; run directories are <framework>-<model>/<GSE> (PROTOCOL.md items 22, 24)')
    s = sub.add_parser('score')
    s.add_argument('--final', required=True)
    s.add_argument('--root', type=Path, required=True)
    s.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    {'prepare': prepare, 'score': score, 'combine': combine}[args.command](args)


if __name__ == '__main__':
    main()
