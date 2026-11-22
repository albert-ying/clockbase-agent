"""Reconcile the two independent reads of the frozen inputs (PROTOCOL.md v2) and list disagreements for a third read.

- Declared contrasts agree when both readers name the same factor family and their case rules, and their reference
  rules, select the same samples among all samples of the series (expert labels and GEO metadata).
- A1 coding agrees per intervention string on control status and the set of (class, name) components, and per strain
  on the partition of strains into backgrounds.
- A2 refinement agrees when both readers refine the same samples to the same subtypes.
- Sample checks agree per sample on the unlabelled category and on the partition of samples into subjects; every label
  conflict reported by either reader goes to the third read.

  python -m benchmark.published_clock_2026_10_05.groupbench.reconcile_inputs --prepared DIR --reads READS.json --out DIR
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

from audit_system import paths
from .classify_v2 import holds

GEO_META = {'A1': paths.GEO / 'GEO_mm_rna/processed_metadata_preharmonize',
            'A2': paths.GEO / 'GEO_hs_dnam/processed_metadata_preharmonize'}
FAMILY = {'intervention_compound': 'intervention', 'intervention_genotype': 'intervention',
          'intervention_diet': 'intervention', 'intervention_other': 'intervention'}


def series_samples(root, set_name, gse):
    """Every GEO sample of the series with its GEO fields and expert labels (for evaluating declared rules)."""
    geo = pd.read_csv(GEO_META[set_name] / f'{gse}.csv', dtype=str).fillna('')
    out = {r['GSM']: {'geo': {k: v for k, v in r.items() if k != 'GSM'}, 'expert': {}} for r in geo.to_dict('records')}
    if set_name == 'A1':
        annotation = pd.read_excel(root / 'a1/source/Data_annotation_absolute_rodents.xlsx').set_index('Sample')
        mapping = pd.read_csv(root / f'a1/mapping/{gse}.csv', dtype=str).fillna('')
        for sample, gsm in zip(mapping['sample'], mapping['gsm']):
            if gsm and gsm in out:
                row = annotation.loc[sample]
                out[gsm]['expert'] = {k: str(row[k]) for k in ('Strain', 'Sex', 'Tissue', 'Tissue.details',
                                                               'Intervention', 'Intervention.type', 'Age.days')}
    else:
        labels = pd.read_csv(root / 'a2/labels.csv', dtype=str).fillna('')
        meta = pd.read_csv(root / 'a2/computage_bench_meta_3cf5268.tsv', sep='\t', dtype=str).fillna('')
        cell = {(r['DatasetID'], r['SampleID']): r['CellType'] for r in meta.to_dict('records')}
        for r in labels.to_dict('records'):
            if r['gsm'] not in out:
                continue
            value = {k: r[k] for k in ('condition', 'class', 'age', 'tissue', 'sex')}
            value['cell_type'] = cell.get((r['dataset'], r['gsm']), '')
            previous = out[r['gsm']]['expert']
            if previous:  # listed under two datasets: fields that disagree become missing
                value = {k: (v if previous.get(k) == v else None) for k, v in value.items()}
            out[r['gsm']]['expert'] = value
    return out


def selected(rule, samples):
    return frozenset(g for g, s in samples.items() if holds(rule, s))


def reconcile_declared(root, reads, units_by_id):
    by_key = defaultdict(dict)
    for read in reads:
        for item in read['declared']:
            by_key[(item['unit_id'], item['variant'])][read['reader']] = item
    frozen, disagreements, cache = {}, [], {}
    for key, pair in sorted(by_key.items()):
        unit = units_by_id.get(key[0])
        if unit is None or len(pair) < 2:
            disagreements.append({'unit_id': key[0], 'variant': key[1], 'reason': 'missing read', 'reads': pair})
            continue
        a, b = pair['A'], pair['B']
        if a['ambiguous'] and b['ambiguous']:
            frozen[key] = {**a, 'status': 'ambiguous'}
            continue
        series = (unit['set'], unit['gse'])
        if series not in cache:
            cache[series] = series_samples(root, *series)
        samples = cache[series]
        same_family = FAMILY.get(a['declared_factor'], a['declared_factor']) == FAMILY.get(b['declared_factor'], b['declared_factor'])
        sets = {r: (selected(pair[r]['case_rule'], samples), selected(pair[r]['reference_rule'], samples)) for r in 'AB'}
        agree = same_family and not a['ambiguous'] and not b['ambiguous'] and sets['A'] == sets['B'] and all(sets['A'])
        if agree:
            frozen[key] = {**a, 'status': 'agreed', 'case_n': len(sets['A'][0]), 'reference_n': len(sets['A'][1])}
        else:
            disagreements.append({'unit_id': key[0], 'variant': key[1], 'reads': pair,
                                  'selected_counts': {r: [len(sets[r][0]), len(sets[r][1])] for r in 'AB'},
                                  'case_sets_equal': sets['A'][0] == sets['B'][0],
                                  'reference_sets_equal': sets['A'][1] == sets['B'][1], 'same_family': same_family})
    return frozen, disagreements


def components(entry):
    return (entry['is_control'], frozenset((c['class'], c['name'].strip().lower()) for c in entry['components']))


def reconcile_coding(reads):
    a, b = ({x['intervention']: x for x in read['interventions']} for read in reads)
    frozen, disagreements = {}, []
    for name in sorted(set(a) | set(b)):
        if name in a and name in b and components(a[name]) == components(b[name]):
            frozen[name] = a[name]
        else:
            disagreements.append({'intervention': name, 'reads': {'A': a.get(name), 'B': b.get(name)}})
    sa, sb = ({x['strain']: x for x in read['strains']} for read in reads)
    def partition(table):
        groups = defaultdict(set)
        for strain, x in table.items():
            groups[x['background'].strip().lower()].add(strain)
        return {frozenset(g) for g in groups.values()}
    strains = {'agreed': partition(sa) == partition(sb), 'A': sa, 'B': sb}
    return frozen, disagreements, strains


def reconcile_refine(reads):
    a, b = ({(x['dataset'], x['gsm']): x['refined_condition'] for x in read['refined']} for read in reads)
    if a == b:
        return {f'{k[0]}|{k[1]}': v for k, v in a.items()}, []
    keys = set(a) | set(b)
    return ({f'{k[0]}|{k[1]}': a[k] for k in keys if a.get(k) == b.get(k)},
            [{'dataset': k[0], 'gsm': k[1], 'A': a.get(k), 'B': b.get(k)} for k in sorted(keys) if a.get(k) != b.get(k)])


def subject_partition(entries):
    groups = defaultdict(set)
    for x in entries:
        if x['subject_id']:
            groups[x['subject_id']].add(x['gsm'])
    return {frozenset(g) for g in groups.values()}


def reconcile_samples(reads):
    by_series = defaultdict(dict)
    for read in reads:
        for series in read['series']:
            by_series[series['gse']][read['reader']] = series
    frozen, disagreements = {}, []
    for gse, pair in sorted(by_series.items()):
        if len(pair) < 2:
            disagreements.append({'gse': gse, 'reason': 'missing read'})
            continue
        a, b = pair['A'], pair['B']
        ua, ub = ({x['gsm']: x for x in s['unlabelled']} for s in (a, b))
        categories, mismatched = {}, []
        for gsm in sorted(set(ua) | set(ub)):
            if gsm in ua and gsm in ub and ua[gsm]['category'] == ub[gsm]['category']:
                categories[gsm] = ua[gsm]
            else:
                mismatched.append({'gsm': gsm, 'A': ua.get(gsm), 'B': ub.get(gsm)})
        subjects_agree = subject_partition(a['subjects']) == subject_partition(b['subjects'])
        conflicts = {(c['gsm'], c['field']): c for s in (a, b) for c in s['label_conflicts']}
        frozen[gse] = {'unlabelled': categories, 'subjects': {x['gsm']: x['subject_id'] for x in a['subjects']}
                       if subjects_agree else None, 'conflicts_to_check': list(conflicts.values())}
        if mismatched or not subjects_agree or conflicts:
            disagreements.append({'gse': gse, 'unlabelled': mismatched, 'subjects_agree': subjects_agree,
                                  'subjects': {'A': a['subjects'], 'B': b['subjects']} if not subjects_agree else None,
                                  'label_conflicts': list(conflicts.values())})
    return frozen, disagreements


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True, help='groupbench data directory')
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--reads', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    reads = json.loads(args.reads.read_text())
    units = {u['unit_id']: u for u in map(json.loads, (args.prepared / 'units.jsonl').read_text().splitlines())}
    declared, declared_dis = reconcile_declared(args.root, reads['declared'], units)
    coding, coding_dis, strains = reconcile_coding(sorted(reads['coding'], key=lambda r: r['reader']))
    refined, refine_dis = reconcile_refine(sorted(reads['refine'], key=lambda r: r['reader']))
    samples, samples_dis = reconcile_samples(reads['samples'])
    args.out.mkdir(parents=True, exist_ok=False)
    write = lambda name, value: (args.out / name).write_text(json.dumps(value, indent=1, default=list) + '\n')
    write('declared_agreed.json', {f'{k[0]}#{k[1]}': v for k, v in declared.items()})
    write('coding_agreed.json', {'interventions': coding, 'strains': strains})
    write('refined_agreed.json', refined)
    write('samples_agreed.json', samples)
    write('disagreements.json', {'declared': declared_dis, 'coding': coding_dis, 'strains_agreed': strains['agreed'],
                                 'refine': refine_dis, 'samples': samples_dis})
    summary = {'declared_agreed': len(declared), 'declared_disagreements': len(declared_dis),
               'coding_agreed': len(coding), 'coding_disagreements': len(coding_dis), 'strain_partition_agreed': strains['agreed'],
               'refined_agreed': len(refined), 'refine_disagreements': len(refine_dis),
               'sample_series': len(samples), 'sample_series_needing_third_read': len(samples_dis)}
    write('summary.json', summary)
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
