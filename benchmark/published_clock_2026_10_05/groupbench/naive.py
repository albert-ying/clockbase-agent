"""Naive-designer reference for Benchmark A (PROTOCOL.md v2): every single-field split of a series' GEO metadata.

For each series with benchmark units and each GEO field with 2 to 8 values (every value held by at least 2 samples,
fields unique per sample skipped), each pair of values forms one comparison whose declared contrast is that field. The
design rules (not the reported-direction rule, which needs a reported result) are applied with the frozen inputs.

  python -m benchmark.published_clock_2026_10_05.groupbench.naive --root G --prepared DIR --frozen DIR --out DIR
"""
import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

from .classify_v2 import classify_unit
from .run_v2 import A1_FACTORS, A2_FACTORS, unit_samples
from .reconcile_inputs import series_samples


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    units = [json.loads(line) for line in (args.prepared / 'units.jsonl').read_text().splitlines()]
    frozen = {name: json.loads((args.frozen / f'{name}.json').read_text()) for name in ('declared', 'coding', 'refined', 'samples')}
    series = sorted({(u['set'], u['gse']) for u in units})
    cache, results = {}, []
    for set_name, gse in series:
        samples = series_samples(args.root, set_name, gse)
        cache[(set_name, gse)] = samples
        fields = {}
        for gsm, record in samples.items():
            for field, value in record['geo'].items():
                fields.setdefault(field, {}).setdefault(value, []).append(gsm)
        for field, groups in fields.items():
            if not 2 <= len(groups) <= 8 or any(len(m) < 2 for m in groups.values()):
                continue
            for v1, v2 in itertools.combinations(sorted(groups), 2):
                declared = {'declared_factor': 'geo_field', 'status': 'agreed',
                            'case_rule': [[{'source': 'geo', 'field': field, 'values': [v1]}]],
                            'reference_rule': [[{'source': 'geo', 'field': field, 'values': [v2]}]]}
                unit = {'unit_id': f'NAIVE:{set_name}:{gse}:{field}:{v1}|{v2}', 'set': set_name, 'gse': gse,
                        'arm_a': sorted(groups[v1]), 'arm_b': sorted(groups[v2]),
                        'variants': [{'index': 0, 'declared': declared, 'case_is_arm_a': True, 'clocks': [],
                                      'series_case_level': groups[v1], 'series_reference_level': groups[v2]}]}
                arm_samples, _ = unit_samples(unit, frozen, args.root, cache)
                factors = A1_FACTORS if set_name == 'A1' else A2_FACTORS
                result = classify_unit(unit, arm_samples, set_name, factors, lambda d: None)
                results.append({'unit_id': unit['unit_id'], 'set': set_name, 'gse': gse, 'class': result['class'],
                                'error_types': result.get('error_types', [])})
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'naive_designer.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in results))
    summary = {s: {'splits': sum(r['set'] == s for r in results),
                   'classes': dict(Counter(r['class'] for r in results if r['set'] == s)),
                   'error_types': dict(Counter(t for r in results if r['set'] == s for t in r['error_types']))}
               for s in ('A1', 'A2')}
    (args.out / 'naive_designer_summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
