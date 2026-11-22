"""Published results recomputed for the paper's groups with ClockBase's released test (PROTOCOL.md, addendum item 19).

S3 of the published-clock benchmark tests the paper's groups with a Welch test. The released verifier tests a
comparison with the design-respecting estimand (individual means, a paired test, age as a covariate), drops overflow
clock values and, for comparisons of interest in a series without a GEO individual field, reads the individuals from
the sample records with two models. This module applies those rules to the paper's groups:

  S3    Welch test (as scored before)
  S3c   released test, individuals from a GEO sample field
  S3r   released test, individuals from a GEO sample field or, where the series has none, from the two readers when both
        report the same shared individuals

  read   builds the reading package of every comparison without an individual field and asks the two readers
  score  writes published_released_test.json and published_released_test.csv (one row per published result); the JSON
         also compares the systems with the no-data models (S4) on the results every one of them answers: three-way call,
         and the direction of published increases and of published decreases separately

  python -m benchmark.published_clock_2026_10_05.groupbench.published_released_test read --bench DATA --out DIR
  python -m benchmark.published_clock_2026_10_05.groupbench.published_released_test score --bench DATA --out DIR
"""
import argparse
import csv
import hashlib
import json
import os
import re
import threading
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from audit_system import design_estimand_2026_10_10 as design
from audit_system import paths
from audit_system.individuals_audit_2026_10_10 import MAX_SAMPLES, MAX_TOKENS, MODELS, package, parse, reconcile, series_text
from audit_system.label_adjudication_2026_09_27 import Ledger, call as ledger_call
from audit_system.row_repairs_2026_09_28 import MIN_ARM, recompute

from ..engine.metrics import paper_key
from ..engine.predictions import PredictionReader, atlas_rows
from ..engine.stats import threshold
from ..summarize_2026_10_05 import call as three_way, judged, read_records, t1_metrics, with_interval
from .design_eval import auc, finding_metrics, signed_z

FOLDER = 'GEO_hs_dnam'
SUMMARY_NAME = os.environ.get('CLOCKBASE_BENCH_SUMMARY') or os.environ.get('CLOCKBASE_FINAL_NAME', 'final-2026-10-05')
ONE_SIDED = re.compile(r'one[- ]sided|one[- ]tailed', re.I)


def published(bench):
    """Published human DNAm results with a call that S3 scored (the paper's groups, ClockBase clock outputs)."""
    inputs = json.loads((bench / 'summaries' / SUMMARY_NAME / 'summary.json').read_text())['inputs']
    final = {json.loads(line)['id'] for line in Path(inputs['ground_truth']).read_text().splitlines()}
    out = []
    for run in inputs['t1']:
        for record in read_records(run):
            s3 = record.get('systems', {}).get('S3', {})
            if (record['id'] in final and record['modality'] == 'human DNAm'
                    and three_way(record.get('direction'), record.get('significant')) is not None
                    and s3.get('covered') is True and isinstance(s3.get('p'), (int, float))
                    and s3.get('direction') in ('higher', 'lower', 'none')):
                out.append(record)
    return out


def unit_key(record):
    groups = record['groups']
    text = json.dumps([record['gse'], sorted(groups['case_gsms']), sorted(groups['control_gsms'])])
    return hashlib.sha256(text.encode()).hexdigest()[:16]


def read(bench, out, max_requests, workers):
    """The two readers on every comparison whose series has no GEO individual field (the verifier's rules)."""
    out.mkdir(parents=True, exist_ok=True)
    units, frames = {}, {}
    for record in published(bench):
        units.setdefault(unit_key(record), record)
    ledger, budget = Ledger(out / 'ledger.jsonl'), {'max': max_requests, 'sent': 0, 'lock': threading.Lock()}

    def one(item):
        key, record = item
        gse, groups = record['gse'], record['groups']
        case, control = [str(g) for g in groups['case_gsms']], [str(g) for g in groups['control_gsms']]
        base = {'unit': key, 'gse': gse, 'case_gsms': case, 'control_gsms': control}
        raw, field = frames[gse]
        if raw is None:
            return {**base, 'status': 'METADATA_UNAVAILABLE'}
        if field:
            return {**base, 'status': 'GEO_INDIVIDUAL_FIELD', 'subject_field': field}
        if len(case) + len(control) > MAX_SAMPLES:
            return {**base, 'status': 'NOT_READ_ARMS_ABOVE_LIMIT'}
        unit = {'gse': gse, 'comparison_name': f"{record.get('case_label')} vs {record.get('control_label')}",
                'condition_name': str(record.get('case_label')), 'case': case, 'control': control}
        prompt = package(unit, raw, series_text(paths.GEO / FOLDER, gse))
        answers, receipts = [], []
        for model in MODELS:
            reply, cached = ledger_call(model, prompt, ledger, budget, parse=parse, max_tokens=MAX_TOKENS)
            receipts.append({'model': model, 'status': reply.get('status'), 'request_sha256': reply.get('request_sha256'), 'cached': cached})
            answers.append(reply.get('parsed') if reply.get('status') == 'OK' else None)
        return {**base, **reconcile(case, control, answers), 'prompt_sha256': hashlib.sha256(prompt.encode()).hexdigest(),
                'evidence': [None if a is None else a.get('individual_evidence') for a in answers], 'receipts': receipts}

    for record in units.values():
        if record['gse'] not in frames:
            raw = design.read_frame(paths.GEO / FOLDER / 'processed_metadata_preharmonize' / f"{record['gse']}.csv")
            frames[record['gse']] = (raw, design.subject_table(raw)[0] if raw is not None else None)
    with ThreadPoolExecutor(workers) as pool:
        results = list(pool.map(one, sorted(units.items())))
    (out / 'individuals.jsonl').write_text(''.join(json.dumps(r, default=str) + '\n' for r in results))
    summary = {'comparisons': len(results), 'status': dict(Counter(r['status'] for r in results)), 'requests_sent': budget['sent'],
               'models': list(MODELS), 'max_samples': MAX_SAMPLES}
    (out / 'individuals_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


def statistic(record, values, series):
    """The released test of one published contrast on the paper's groups."""
    groups = record['groups']
    case, control = [str(g) for g in groups['case_gsms']], [str(g) for g in groups['control_gsms']]
    dropped = set(design.overflow(values, case + control))
    case, control = [g for g in case if g not in dropped], [g for g in control if g not in dropped]
    names = (str(record.get('case_label') or ''), str(record.get('control_label') or ''))
    estimate = design.estimate(values, case, control, series, 'dnam', record['clock_clockbase'], names, ())
    if estimate is not None and estimate['recomputed'] is not None:
        result, method = estimate['recomputed'], estimate['design']['method']
    else:
        result, method = recompute(values, case, control), 'WELCH'
    diff, p = result.get('mean_diff'), result.get('p')
    tested = p is not None and diff is not None and min(result['n_case'], result['n_control']) >= MIN_ARM
    return {'method': method, 'tested': tested, 'p': p if tested else None, 'mean_diff': diff,
            'direction': None if diff is None else 'higher' if diff > 0 else 'lower' if diff < 0 else 'none',
            'n_case': result.get('n_case'), 'n_control': result.get('n_control'), 'overflow_samples_dropped': len(dropped)}


def score(bench, out):
    records = published(bench)
    individuals = {}
    path = out / 'individuals.jsonl'
    if path.is_file():
        individuals = {r['unit']: r for r in map(json.loads, path.read_text().splitlines())}
    reader, series_cache, rows, items = PredictionReader(), {}, [], {'S3': [], 'S3c': [], 'S3r': []}
    for record in records:
        gse, truth, alpha = record['gse'], three_way(record['direction'], record['significant']), threshold(record) or 0.05
        prediction = reader.load(record['modality'], gse, record['clock_clockbase'], atlas_rows(record['modality'], gse))
        values = prediction['samples']['clock']
        if gse not in series_cache:
            series_cache[gse] = design.load_series(paths.GEO, FOLDER, gse, True)
        series = series_cache[gse]
        read_unit = individuals.get(unit_key(record), {})
        with_readers = series
        if not series.get('subject_field') and read_unit.get('status') == 'SHARED_AGREED' and read_unit.get('individuals'):
            with_readers = {**series, 'subjects': read_unit['individuals'], 'subject_source': 'two-reader individuals audit'}
        s3 = record['systems']['S3']
        results = {'S3': {'method': 'WELCH', 'tested': True, 'p': s3['p'], 'direction': s3['direction'], 'mean_diff': s3.get('mean_diff'),
                          'n_case': s3.get('n_case'), 'n_control': s3.get('n_control')},
                   'S3c': statistic(record, values, series), 'S3r': statistic(record, values, with_readers)}
        row = {'id': record['id'], 'gse': gse, 'clock': record['clock_clockbase'], 'paper': paper_key(record), 'published_call': truth,
               'published_p': record.get('p_value'), 'published_test': record.get('test'),
               'published_one_sided': bool(ONE_SIDED.search(str(record.get('test') or ''))),
               'published_n_case': record.get('n_case'), 'published_n_control': record.get('n_control'),
               'individuals_read': read_unit.get('status')}
        for system in ('S3c', 'S3r'):  # as systems of the record, for the functions of the published-clock summary
            record['systems'][system] = {'covered': bool(results[system]['tested']), 'direction': results[system]['direction'],
                                         'p': results[system]['p']}
        for system, result in results.items():
            call = None
            if result['tested'] and result['direction'] is not None:
                sign = {'higher': 1.0, 'lower': -1.0, 'none': 0.0}[result['direction']]
                call = result['direction'] if result['p'] < alpha and sign else 'none'
                items[system].append({'cluster': paper_key(record), 'truth': truth, 'z': signed_z(sign, result['p']), 'call': call,
                                      'id': record['id']})
            row.update({f'{system}_method': result['method'], f'{system}_call': call, f'{system}_p': result['p'],
                        f'{system}_direction': result['direction'], f'{system}_n_case': result['n_case'],
                        f'{system}_n_control': result['n_control']})
        rows.append(row)
    significant = [r for r in rows if r['published_call'] != 'none']

    def outcome(row, system):
        call = row[f'{system}_call']
        if call is None:
            return 'not testable'
        if call == row['published_call']:
            return 'significant, same direction'
        if call != 'none':
            return 'significant, opposite direction'
        return 'not significant, same direction' if row[f'{system}_direction'] == row['published_call'] else 'not significant, opposite direction'

    result = {'records': len(rows), 'papers': len({r['paper'] for r in rows}), 'series': len({r['gse'] for r in rows}),
              'systems': {s: finding_metrics(i) for s, i in items.items()},
              'methods': {s: dict(Counter(r[f'{s}_method'] for r in rows)) for s in ('S3c', 'S3r')},
              'not_testable': {s: sum(r[f'{s}_call'] is None for r in rows) for s in items},
              'published_significant': {s: dict(Counter(outcome(r, s) for r in significant)) for s in items},
              'individuals_read': dict(Counter(r['status'] for r in individuals.values())),
              'one_sided_published_tests': {
                  'results': sum(r['published_one_sided'] for r in rows),
                  'not_recovered_by_S3r_two_sided_but_by_one_sided': sum(
                      r['published_one_sided'] and r['published_call'] != 'none' and r['S3r_call'] == 'none'
                      and r['S3r_direction'] == r['published_call'] and r['S3r_p'] is not None and r['S3r_p'] / 2 < 0.05 for r in rows)},
              'not_recovered_by_S3r': [
                  {k: r[k] for k in ('id', 'gse', 'clock', 'published_call', 'published_p', 'published_test', 'published_one_sided',
                                     'published_n_case', 'published_n_control', 'individuals_read', 'S3r_method', 'S3r_direction',
                                     'S3r_p', 'S3r_n_case', 'S3r_n_control')} | {'outcome': outcome(r, 'S3r')}
                  for r in significant if r['S3r_call'] != r['published_call']]}
    # Against the no-data models, on the results that S3, the released test and every model answer.
    models = sorted({name for record in records for name in record['systems'] if name.startswith('S4:')})
    compared = ['S0', 'S3', 'S3c', 'S3r'] + models
    answered = [r for r in records if all(judged(r, name) for name in compared)]
    signed = [r for r in answered if r.get('significant') is True and r.get('direction') in ('higher', 'lower')]
    block = {'records': len(answered), 'papers': len({paper_key(r) for r in answered}),
             'published_significant': dict(Counter(r['direction'] for r in signed)), 'systems': {}}
    for name in compared:
        entry = with_interval(answered, name, t1_metrics, paper_key, ('direction_agreement', 'call_agreement', 'call_balanced_accuracy'))
        for sign in ('higher', 'lower'):
            members = [r for r in signed if r['direction'] == sign]
            entry[f'published_{sign}_identified'] = [sum(judged(r, name)['direction'] == sign for r in members), len(members)]
        entry['direction_balanced_accuracy'] = sum(k / n for k, n in (entry['published_higher_identified'], entry['published_lower_identified'])) / 2
        if name.startswith('S4:'):  # the model's stated confidence, signed by its direction
            scores = [(r['direction'], {'higher': 1.0, 'lower': -1.0}.get(r['systems'][name]['direction'], 0.0)
                       * float(r['systems'][name].get('confidence') or 0.0)) for r in signed]
            entry['direction_auroc_signed_confidence'] = auc([z for d, z in scores if d == 'higher'], [z for d, z in scores if d == 'lower'])
        block['systems'][name] = entry
    result['against_models_without_data'] = block
    out.mkdir(parents=True, exist_ok=True)
    (out / 'published_released_test.json').write_text(json.dumps(result, indent=1, default=float) + '\n')
    with (out / 'published_released_test.csv').open('w', newline='') as sink:
        writer = csv.DictWriter(sink, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for system, m in result['systems'].items():
        print(f"{system:4s} direction AUROC {m['direction_auroc']:.3f} {[round(x, 2) for x in m['direction_auroc_interval']]}  "
              f"same direction both significant {m['same_direction_both_significant']}  replicated {m['replicated']}  "
              f"published null and ClockBase significant {m['significant_only_in_clockbase']}  detection AUROC {m['detection_auroc']:.3f}  "
              f"not testable {result['not_testable'][system]}")
    print(json.dumps({k: result[k] for k in ('methods', 'published_significant', 'one_sided_published_tests')}, indent=1))
    print('against models without data:', block['records'], 'results;', block['published_significant'])
    for name, m in block['systems'].items():
        print(f"  {name:24s} call {m['call_agreement']:.3f} (n {m['call_n']})  balanced {m['call_balanced_accuracy']:.3f}  "
              f"increases {m['published_higher_identified']}  decreases {m['published_lower_identified']}  "
              f"AUROC(confidence) {m.get('direction_auroc_signed_confidence')}")


def funnel(bench, out):
    """Why each published result of the ground truth is or is not scored with the paper's groups (S3, S3c)."""
    inputs = json.loads((bench / 'summaries' / SUMMARY_NAME / 'summary.json').read_text())['inputs']
    final = {json.loads(line)['id'] for line in Path(inputs['ground_truth']).read_text().splitlines()}
    scored = {r['id'] for r in published(bench)}
    reasons, rna_clocks, rna_papers, seen = {}, {}, set(), set()
    for run in inputs['t1']:
        for record in read_records(run):
            if record['id'] not in final or record['id'] in seen:
                continue
            seen.add(record['id'])
            s3 = record.get('systems', {}).get('S3', {})
            if record['modality'] != 'human DNAm':
                reason = f"{record['modality']}: clock not computed by ClockBase" if record.get('clock_clockbase') is None \
                    else f"{record['modality']}: not scored"
                rna_clocks[record.get('clock_reported')] = rna_clocks.get(record.get('clock_reported'), 0) + 1
                rna_papers.add(paper_key(record))
            elif three_way(record.get('direction'), record.get('significant')) is None:
                reason = 'DNAm: no published direction or significance'
            elif s3.get('covered') is not True:
                reason = f"DNAm: {s3.get('reason') or s3.get('status') or 'not covered'}"
            elif record['id'] not in scored:
                reason = 'DNAm: covered, no statistic'
            else:
                reason = 'DNAm: scored with the paper\'s groups'
            reasons[reason] = reasons.get(reason, 0) + 1
    result = json.loads((out / 'published_released_test.json').read_text()) if (out / 'published_released_test.json').exists() else {}
    summary = {'ground_truth': len(final), 'read': len(seen), 'reasons': dict(sorted(reasons.items(), key=lambda kv: -kv[1])),
               'scored_s3': len(scored),
               'testable_s3c': result['records'] - result['not_testable']['S3c'] if result else None,
               'rna_papers': len(rna_papers), 'rna_clocks_reported': dict(sorted(rna_clocks.items(), key=lambda kv: -kv[1]))}
    (out / 'ground_truth_funnel.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1)[:1500])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('command', choices=('read', 'score', 'funnel'))
    parser.add_argument('--bench', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-requests', type=int, default=600)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()
    if args.command == 'read':
        read(args.bench, args.out, args.max_requests, args.workers)
    elif args.command == 'funnel':
        funnel(args.bench, args.out)
    else:
        score(args.bench, args.out)


if __name__ == '__main__':
    main()
