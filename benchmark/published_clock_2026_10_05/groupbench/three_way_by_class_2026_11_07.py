"""Recall of each published call under the three-way call, for ClockBase and the models without the data.

The comparison with language models in published_released_test.py reports, per system, the balanced accuracy of the
three-way call (higher, lower, not significant) and the number of published increases and decreases to which the
system assigned the published direction. For ClockBase that direction is the sign of its estimate whatever its
significance, whereas a model may answer that there is no difference. This script scores every system the same way:
for each published call it counts the results the system called identically (lower and significant, higher and
significant, or not significant), on the results that every system answered. The mean of the three recalls is the
balanced accuracy, and the script fails unless it reproduces the stored values of published_released_test.json.

  python -m benchmark.published_clock_2026_10_05.groupbench.three_way_by_class_2026_11_07

Inputs: the benchmark records of the summary's t1 runs (CLOCKBASE_BENCHMARK_DATA; server paths remapped to it), the
released-test calls of results_2026_10_10/published_released_test.csv (S3c, S3r) and the stored summary
results_2026_10_10/published_released_test.json. Output: results_2026_10_10/three_way_by_class.json.
"""
import csv
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

from ..summarize_2026_10_05 import call, judged, read_records

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results_2026_10_10'
SERVER = '/Entropy/kying/clockbase-benchmark-2026-10-05'
CLASSES = ('lower', 'higher', 'none')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    bench = Path(os.environ['CLOCKBASE_BENCHMARK_DATA'])
    local = lambda p: Path(str(p).replace(SERVER, str(bench)))
    summary = bench / 'summaries/final-2026-10-05/summary.json'
    inputs = json.loads(summary.read_text())['inputs']
    released_csv = RESULTS / 'published_released_test.csv'
    released = {r['id']: r for r in csv.DictReader(released_csv.open())}
    records = {}
    for run in inputs['t1']:
        for record in read_records(local(run)):
            if record['id'] in released and record['id'] not in records:
                records[record['id']] = record
    if len(records) != len(released):
        raise SystemExit(f'{len(released) - len(records)} released-test records not found in the benchmark runs')
    for record in records.values():  # the released test as scored in published_released_test.score()
        row = released[record['id']]
        for system in ('S3c', 'S3r'):
            if row[f'{system}_direction'] in ('higher', 'lower', 'none') and row[f'{system}_p'] != '':
                record['systems'][system] = {'covered': True, 'direction': row[f'{system}_direction'],
                                             'p': float(row[f'{system}_p']), 'significant': row[f'{system}_call'] != 'none'}
    models = sorted({name for r in records.values() for name in r['systems'] if name.startswith('S4:')})
    compared = ['S0', 'S3', 'S3c', 'S3r'] + models
    answered = [r for r in records.values() if all(judged(r, name) for name in compared)]
    published = {r['id']: call(r.get('direction'), r.get('significant')) for r in answered}
    stored = json.loads((RESULTS / 'published_released_test.json').read_text())['against_models_without_data']
    if len(answered) != stored['records']:
        raise SystemExit(f"{len(answered)} answered results, stored {stored['records']}")
    out = {'records': len(answered), 'published': dict(Counter(published.values())), 'systems': {}}
    for name in compared:
        recall = {}
        for c in CLASSES:
            members = [r for r in answered if published[r['id']] == c]
            hits = sum(call(judged(r, name)['direction'], judged(r, name)['significant']) == c for r in members)
            recall[c] = [hits, len(members)]
        balanced = sum(k / n for k, n in recall.values()) / len(CLASSES)
        if abs(balanced - stored['systems'][name]['call_balanced_accuracy']) > 1e-9:
            raise SystemExit(f'{name}: balanced accuracy {balanced} differs from the stored value')
        out['systems'][name] = {'recall': recall, 'balanced_accuracy': balanced,
                                'balanced_accuracy_interval': stored['systems'][name].get('call_balanced_accuracy_interval')}
    out['inputs'] = {str(p): sha(p) for p in (summary, released_csv, RESULTS / 'published_released_test.json')}
    (RESULTS / 'three_way_by_class.json').write_text(json.dumps(out, indent=1) + '\n')
    for name, m in out['systems'].items():
        print(f"{name:24s} {m['recall']}  balanced accuracy {m['balanced_accuracy']:.3f}")


if __name__ == '__main__':
    main()
