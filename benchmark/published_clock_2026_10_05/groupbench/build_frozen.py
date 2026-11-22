"""Freeze the Benchmark A inputs (PROTOCOL.md v2): agreed reads plus third-read decisions, with hashes.

Outputs (frozen-v2/): declared.json ('unit#variant' -> rules and status), coding.json (intervention components and
strain backgrounds), refined.json ('dataset|gsm' -> GEO subtype), samples.json (per series: unlabelled categories,
subjects, confirmed label conflicts mapped to the annotation's field names), FROZEN.json (hashes and provenance).

  python -m benchmark.published_clock_2026_10_05.groupbench.build_frozen --reconciled DIR --third THIRD.json --out DIR
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

A1_FIELDS = {'intervention': 'Intervention', 'intervention type': 'Intervention.type', 'genotype': 'Intervention',
             'compound': 'Intervention', 'diet': 'Intervention', 'treatment': 'Intervention', 'strain': 'Strain',
             'background': 'Strain', 'sex': 'Sex', 'tissue': 'Tissue.details', 'tissue.details': 'Tissue.details',
             'cell type': 'Tissue.details', 'cell_type': 'Tissue.details', 'age': 'Age.days', 'age.days': 'Age.days'}
A2_FIELDS = {'condition': 'condition', 'class': 'class', 'disease': 'condition', 'diagnosis': 'condition', 'age': 'age',
             'tissue': 'tissue', 'cell type': 'cell_type', 'cell_type': 'cell_type', 'celltype': 'cell_type', 'sex': 'sex',
             'gender': 'sex'}


def field_name(raw, set_name):
    key = raw.strip().lower().replace('-', ' ')
    table = A1_FIELDS if set_name == 'A1' else A2_FIELDS
    if key in table:
        return table[key]
    for word, target in table.items():
        if word in key:
            return target
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--reconciled', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--third', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    agreed = {name: json.loads((args.reconciled / f'{name}_agreed.json').read_text())
              for name in ('declared', 'coding', 'refined', 'samples')}
    third = json.loads(args.third.read_text())
    units = {u['unit_id']: u for u in map(json.loads, (args.prepared / 'units.jsonl').read_text().splitlines())}
    log = {'unmapped_conflict_fields': [], 'declared_from_third_read': 0, 'ambiguous': 0}

    declared = {}
    for key, value in agreed['declared'].items():
        declared[key] = {**value, 'status': value.get('status', 'agreed')}
    for item in (x for r in third if r['kind'] == 'declared' for x in r['declared']):
        key = f"{item['unit_id']}#{item['variant']}"
        declared[key] = {**item, 'status': 'ambiguous' if item['ambiguous'] else 'third_read'}
        log['declared_from_third_read'] += 1
    missing = [f"{u}#{i}" for u, unit in units.items() for i in range(len(unit['variants'])) if f"{u}#{i}" not in declared]
    log['declared_missing'] = missing
    log['ambiguous'] = sum(v['status'] == 'ambiguous' for v in declared.values())

    strains = agreed['coding']['strains']
    coding = {'interventions': agreed['coding']['interventions'],
              'backgrounds': {s: x['background'] for s, x in strains['A'].items()} if strains['agreed'] else None}

    refined = dict(agreed['refined'])
    samples = {gse: {'unlabelled': value['unlabelled'], 'subjects': value['subjects'] or {}, 'conflicts': []}
               for gse, value in agreed['samples'].items()}
    set_of = {u['gse']: u['set'] for u in units.values()}
    for result in (r for r in third if r['kind'] == 'samples'):
        for series in result['series']:
            gse = series['gse']
            entry = samples.setdefault(gse, {'unlabelled': {}, 'subjects': {}, 'conflicts': []})
            for x in series['unlabelled']:
                entry['unlabelled'][x['gsm']] = {'gsm': x['gsm'], 'category': x['category'], 'label_like': x['label_like'],
                                                 'reason': x['reason']}
            if series['subjects'] or 'none' in series['subjects_resolution'].lower() or 'no subject' in series['subjects_resolution'].lower():
                entry['subjects'] = {x['gsm']: x['subject_id'] for x in series['subjects'] if x['subject_id']}
            for c in series['confirmed_conflicts']:
                field = field_name(c['field'], set_of.get(gse, 'A2'))
                if field is None:
                    log['unmapped_conflict_fields'].append({'gse': gse, **c})
                    continue
                entry['conflicts'].append({'gsm': c['gsm'], 'field': field, 'expert_value': c['expert_value'],
                                           'geo_value': c['geo_value']})
            for r in series['refined']:
                if r['refined_condition']:
                    refined[f"{r['dataset']}|{r['gsm']}"] = r['refined_condition']
    args.out.mkdir(parents=True, exist_ok=False)
    tables = {'declared': declared, 'coding': coding, 'refined': refined, 'samples': samples}
    hashes = {}
    for name, value in tables.items():
        text = json.dumps(value, indent=1, sort_keys=True) + '\n'
        (args.out / f'{name}.json').write_text(text)
        hashes[name] = hashlib.sha256(text.encode()).hexdigest()
    frozen = {'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'sha256': hashes, 'log': log,
              'counts': {'declared': len(declared), 'interventions': len(coding['interventions']),
                         'refined': len(refined), 'series': len(samples),
                         'conflicts': sum(len(v['conflicts']) for v in samples.values())}}
    (args.out / 'FROZEN.json').write_text(json.dumps(frozen, indent=1) + '\n')
    print(json.dumps(frozen, indent=1))


if __name__ == '__main__':
    main()
