"""Summarize engine runs across systems: shared record sets, significance at a stated or default threshold, S4 by
publication recognition, and T2 expected directions.

The engine scores each system on the records it covers, so systems differ in their records. This summary adds the
comparisons on records that every compared system covers. Significance uses the paper's stated threshold and P < 0.05
where the paper states none (labeled in the output). S4 results are split by whether the model said it recognized the
publication. T2 records carry an expected direction (higher in the condition) and no significance claim; a covered
record agrees when the system's direction is higher. Intervals resample papers (T1) or series (T2). S0 is a constant
rule (every contrast higher in the case and significant) scored on the same records, the base rate any system must beat.

  python -m benchmark.published_clock_2026_10_05.summarize_2026_10_05 --t1 RUN --t2 RUN --out DIR
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from .engine.metrics import paper_key
from .engine.stats import threshold

DEFAULT_ALPHA = 0.05
REPLICATES, SEED = 1000, 20261005


def read_records(run):
    return [json.loads(line) for line in (Path(run) / 'records.jsonl').read_text().splitlines()]


BASELINE = 'S0'


def systems_of(records):
    found = {s for r in records for s in r.get('systems', {})}
    return [BASELINE] + [s for s in ('S1', 'S2', 'S3') if s in found] + sorted(s for s in found if s.startswith('S4:'))


def judged(record, system):
    """Direction and significance of one covered system result, or None when the system did not cover the record."""
    if system == BASELINE:
        return {'direction': 'higher', 'significant': True, 'alpha': None, 'alpha_assumed': False,
                'recognizes_publication': None}
    result = record.get('systems', {}).get(system, {})
    if result.get('covered') is not True or result.get('direction') not in ('higher', 'lower', 'none'):
        return None
    significant = result.get('significant') if isinstance(result.get('significant'), bool) else None
    stated = threshold(record)
    if significant is None and isinstance(result.get('p'), (int, float)):
        significant = result['p'] < (stated or DEFAULT_ALPHA)
    return {'direction': result['direction'], 'significant': significant,
            'alpha': None if system.startswith('S4:') else (stated or DEFAULT_ALPHA),
            'alpha_assumed': not system.startswith('S4:') and stated is None,
            'recognizes_publication': result.get('recognizes_publication')}


def rate(pairs):
    return float(np.mean([a == b for a, b in pairs])) if pairs else None


def call(direction, significant):
    """Three-class call: the direction when significant, else none."""
    return direction if significant is True and direction in ('higher', 'lower') else 'none' if significant is False else None


def balanced(pairs):
    """Mean recall over the published classes present (higher, lower, none)."""
    classes = sorted({a for a, _ in pairs})
    recalls = [np.mean([b == a for x, b in pairs if x == a]) for a in classes]
    return float(np.mean(recalls)) if recalls else None


def t1_metrics(items, system):
    judgments = [(r, judged(r, system)) for r in items]
    covered = [(r, j) for r, j in judgments if j]
    calls = [(call(r.get('direction'), r.get('significant')), call(j['direction'], j['significant'])) for r, j in covered]
    calls = [(a, b) for a, b in calls if a is not None and b is not None]
    direction = [(r['direction'], j['direction']) for r, j in covered
                 if r.get('significant') is True and r.get('direction') in ('higher', 'lower')]
    significance = [(r['significant'], j['significant']) for r, j in covered
                    if isinstance(r.get('significant'), bool) and isinstance(j['significant'], bool)]
    return {'records': len(items), 'covered': len(covered), 'direction_agreement': rate(direction),
            'direction_n': len(direction), 'call_agreement': rate(calls), 'call_balanced_accuracy': balanced(calls),
            'call_n': len(calls), 'call_classes': {c: sum(a == c for a, _ in calls) for c in ('higher', 'lower', 'none')},
            'significance_agreement': rate(significance),
            'significance_n': len(significance),
            'significance_alpha_assumed': sum(j['alpha_assumed'] for r, j in covered
                                              if isinstance(r.get('significant'), bool) and isinstance(j['significant'], bool))}


def t2_metrics(items, system):
    covered = [(r, j) for r, j in ((r, judged(r, system)) for r in items) if j]
    expected = [(r['direction'], j['direction']) for r, j in covered]
    significant_expected = [(True, j['significant'] is True and j['direction'] == r['direction']) for r, j in covered]
    return {'records': len(items), 'covered': len(covered), 'direction_agreement': rate(expected),
            'direction_n': len(expected), 'significant_in_expected_direction': rate(significant_expected)}


def with_interval(items, system, compute, unit, keys=('direction_agreement',)):
    point = compute(items, system)
    groups = defaultdict(list)
    for record in items:
        groups[unit(record)].append(record)
    names = sorted(groups)
    rng, draws = np.random.default_rng(SEED), {k: [] for k in keys}
    if len(names) >= 2:
        for _ in range(REPLICATES):
            sample = [r for i in rng.integers(0, len(names), len(names)) for r in groups[names[i]]]
            value = compute(sample, system)
            for key in keys:
                if value[key] is not None:
                    draws[key].append(value[key])
    for key in keys:
        point[f'{key}_interval'] = ([float(np.quantile(draws[key], 0.025)), float(np.quantile(draws[key], 0.975))]
                                    if draws[key] else None)
    point['units'] = len(names)
    return point


def shared(records, systems):
    return [r for r in records if all(judged(r, s) for s in systems)]


def summarize(t1, t2):
    output = {'default_alpha': DEFAULT_ALPHA, 'bootstrap': {'replicates': REPLICATES, 'seed': SEED}}
    if t1 is not None:
        systems = systems_of(t1)
        s4 = [s for s in systems if s.startswith('S4:')]
        mapped = [r for r in t1 if r.get('clock_clockbase')]
        keys = ('direction_agreement', 'call_balanced_accuracy')
        tables = {'all_records': {s: with_interval(t1, s, t1_metrics, paper_key, keys) for s in systems},
                  'clockbase_clock_records': {s: with_interval(mapped, s, t1_metrics, paper_key, keys) for s in systems}}
        groups = [('shared_S1_S2_S3', [BASELINE, 'S1', 'S2', 'S3'])]
        if s4:
            groups += [('shared_S3_S4', [BASELINE, 'S3'] + s4), ('shared_all', systems)]
        for name, members in groups:
            if all(m in systems for m in members):
                subset = shared(t1, members)
                tables[name] = {'systems': members, 'records': len(subset),
                                'papers': len({paper_key(r) for r in subset}),
                                'metrics': {s: with_interval(subset, s, t1_metrics, paper_key, keys) for s in members}}
        recognition = {}
        for system in s4:
            for flag in (True, False):
                subset = [r for r in t1 if (judged(r, system) or {}).get('recognizes_publication') is flag]
                recognition[f'{system}|recognized={flag}'] = with_interval(subset, system, t1_metrics, paper_key, keys)
        tables['s4_by_recognition'] = recognition
        # S3 against each model on the records that model covered without recognizing the publication.
        unrecognized = {}
        for system in s4:
            subset = [r for r in t1 if judged(r, 'S3') and (judged(r, system) or {}).get('recognizes_publication') is False]
            unrecognized[system] = {'records': len(subset), 'papers': len({paper_key(r) for r in subset}),
                                    'metrics': {s: with_interval(subset, s, t1_metrics, paper_key, keys)
                                                for s in (BASELINE, 'S3', system)}}
        tables['s3_vs_s4_unrecognized'] = unrecognized
        output['T1'] = tables
    if t2 is not None:
        systems = [s for s in systems_of(t2) if s not in (BASELINE,) and not s.startswith('S4:')]
        series = lambda r: r['gse']
        keys = ('direction_agreement', 'significant_in_expected_direction')
        output['T2'] = {'all_records': {s: with_interval(t2, s, t2_metrics, series, keys) for s in systems},
                        'shared_S1_S2_S3': {'records': len(shared(t2, systems)),
                                            'metrics': {s: with_interval(shared(t2, systems), s, t2_metrics, series, keys)
                                                        for s in systems}}}
    return output


def fmt(value, interval=None, n=None):
    if value is None:
        return 'undefined' + (f' ({n})' if n is not None else '')
    text = f'{value:.3f}'
    if interval:
        text += f' [{interval[0]:.3f}, {interval[1]:.3f}]'
    return text + (f' ({n})' if n is not None else '')


def balanced_cell(m):
    classes = m['call_classes']
    return (fmt(m['call_balanced_accuracy'], m.get('call_balanced_accuracy_interval'))
            + f" ({classes['higher']}/{classes['lower']}/{classes['none']})")


def markdown(summary):
    lines = ['# Published-clock benchmark: shared-record comparison', '',
             'Direction agreement: published-significant contrasts (T1) or expected direction (T2). Call balanced '
             'accuracy: mean recall over the published calls present (higher, lower, or none = not significant), each '
             'system calling its direction when significant and none otherwise; S0 (constant higher, significant) scores '
             'at most 1/k with k classes present. Significance '
             f'agreement applies the paper\'s stated threshold, else P < {DEFAULT_ALPHA} (count of assumed thresholds '
             'shown). Brackets: 95% bootstrap intervals over papers (T1) or series (T2); (n) = denominator.', '']
    t1 = summary.get('T1')
    if t1:
        for name, title in (('all_records', 'T1, all records'), ('clockbase_clock_records', 'T1, records with a ClockBase clock')):
            lines += [f'## {title}', '', '| System | Covered / records | Direction agreement | Call balanced accuracy '
                      '(higher/lower/none n) | Significance agreement | Assumed thresholds |',
                      '| --- | --- | --- | --- | --- | --- |']
            for system, m in t1[name].items():
                lines.append(f"| {system} | {m['covered']} / {m['records']} | "
                             f"{fmt(m['direction_agreement'], m['direction_agreement_interval'], m['direction_n'])} | "
                             f"{balanced_cell(m)} | "
                             f"{fmt(m['significance_agreement'], None, m['significance_n'])} | {m['significance_alpha_assumed']} |")
            lines.append('')
        for name in ('shared_S1_S2_S3', 'shared_S3_S4', 'shared_all'):
            if name in t1:
                block = t1[name]
                lines += [f"## T1, records covered by every system of {', '.join(s for s in block['systems'] if s != BASELINE)} "
                          f"({block['records']} records, {block['papers']} papers)", '',
                          '| System | Direction agreement | Call balanced accuracy (higher/lower/none n) | Significance agreement |',
                          '| --- | --- | --- | --- |']
                for system, m in block['metrics'].items():
                    lines.append(f"| {system} | {fmt(m['direction_agreement'], m['direction_agreement_interval'], m['direction_n'])} | "
                                 f"{balanced_cell(m)} | {fmt(m['significance_agreement'], None, m['significance_n'])} |")
                lines.append('')
        if t1.get('s4_by_recognition'):
            lines += ['## T1, LLM-only predictions by whether the model said it recognized the publication', '',
                      '| Model and recognition | Covered | Direction agreement | Call balanced accuracy | Significance agreement |',
                      '| --- | --- | --- | --- | --- |']
            for key, m in t1['s4_by_recognition'].items():
                lines.append(f"| {key} | {m['covered']} | {fmt(m['direction_agreement'], m['direction_agreement_interval'], m['direction_n'])} | "
                             f"{balanced_cell(m)} | {fmt(m['significance_agreement'], None, m['significance_n'])} |")
            lines.append('')
    if t1 and t1.get('s3_vs_s4_unrecognized'):
        lines += ['## T1, S3 against each model on records the model did not recognize (and S3 covers)', '',
                  '| Model | Records (papers) | System | Direction agreement | Call balanced accuracy | Significance agreement |',
                  '| --- | --- | --- | --- | --- | --- |']
        for model, block in t1['s3_vs_s4_unrecognized'].items():
            for system, m in block['metrics'].items():
                lines.append(f"| {model} | {block['records']} ({block['papers']}) | {system} | "
                             f"{fmt(m['direction_agreement'], m['direction_agreement_interval'], m['direction_n'])} | "
                             f"{balanced_cell(m)} | {fmt(m['significance_agreement'], None, m['significance_n'])} |")
        lines.append('')
    t2 = summary.get('T2')
    if t2:
        lines += ['## T2 (ComputAgeBench expected directions, one record per ClockBase clock)', '',
                  '| System | Covered / records | Expected direction | Significant in expected direction |',
                  '| --- | --- | --- | --- |']
        for system, m in t2['all_records'].items():
            lines.append(f"| {system} | {m['covered']} / {m['records']} | "
                         f"{fmt(m['direction_agreement'], m['direction_agreement_interval'], m['direction_n'])} | "
                         f"{fmt(m['significant_in_expected_direction'], m['significant_in_expected_direction_interval'])} |")
        lines.append('')
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--t1', type=Path, nargs='+', help='one or more engine runs (their records are combined)')
    parser.add_argument('--ground-truth', type=Path, help='score only the T1 records whose ids are in this file')
    parser.add_argument('--t2', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    t1 = [r for run in args.t1 for r in read_records(run)] if args.t1 else None
    if t1 is not None and args.ground_truth:
        keep = {json.loads(line)['id'] for line in args.ground_truth.read_text().splitlines()}
        t1 = [r for r in t1 if r['id'] in keep]
    if t1 is not None and len({r['id'] for r in t1}) != len(t1):
        raise SystemExit('the T1 runs share record ids')
    summary = summarize(t1, read_records(args.t2) if args.t2 else None)
    summary['inputs'] = {'t1': [str(run) for run in args.t1] if args.t1 else None, 't2': str(args.t2) if args.t2 else None,
                         'ground_truth': str(args.ground_truth) if args.ground_truth else None,
                         't1_records': len(t1) if t1 is not None else None}
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    (args.out / 'report.md').write_text(markdown(summary))
    print(markdown(summary))


if __name__ == '__main__':
    main()
