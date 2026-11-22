"""Run Benchmark A of PROTOCOL.md: label agent comparisons with expert annotations and classify them.

A1 (mouse RNA-seq): the Tyshkovskiy et al. 2026 per-sample annotation, mapped to GSMs (a1/mapping/<GSE>.csv). A
series is used only when its mapping is one-to-one, every GSM belongs to the series and the mapped counts per expert
group equal Supplementary Table 1B. Design factors: intervention, tissue (detail), sex, strain, age.
A2 (human DNAm): ComputAgeBench labels (a2/labels.csv). Design factors: condition and tissue; sex and age are cohort
characteristics the agent cannot choose and are reported, not classified.

Writes classified.jsonl (every comparison of the covered series, with its exclusion reason or class), summary.json
and adjudication_queue.jsonl (all error candidates and a seeded 10% sample of Correct comparisons).

  python -m benchmark.published_clock_2026_10_05.groupbench.run_a --root GROUPBENCH --out DIR
"""
import argparse
import json
import random
from collections import Counter
from pathlib import Path

import pandas as pd

from .classify import ERROR_CANDIDATES, classify, completeness

SEED = 20261006
A1_FACTORS = ['intervention', 'tissue', 'sex', 'strain', 'age']
A2_FACTORS = ['condition', 'tissue']


def st1b_counts(path):
    table = pd.read_excel(path, sheet_name='(B) Aggregated rodent data', header=0)
    table.columns = [str(c).strip() for c in table.columns]
    return table


def a1_labels(root, log):
    """GSM -> expert label for every A1 series whose mapping validates; log the series that do not."""
    annotation = pd.read_excel(root / 'a1/source/Data_annotation_absolute_rodents.xlsx')
    st1b = st1b_counts(root / 'a1/source/Supplementary_Table_1_MOESM3.xlsx')
    labels, valid = {}, []
    for path in sorted((root / 'a1/mapping').glob('GSE*.csv')):
        gse = path.stem
        mapping = pd.read_csv(path, dtype=str).fillna('')
        mapped = mapping[mapping.gsm != '']
        rows = annotation[annotation.Source == gse].set_index('Sample')
        problems = []
        if mapped.gsm.duplicated().any() or mapped['sample'].duplicated().any():
            problems.append('mapping not one-to-one')
        if len(mapped) != len(rows):
            problems.append(f'{len(rows) - len(mapped)} of {len(rows)} annotated samples unmapped')
        joined = rows.loc[[s for s in mapped['sample'] if s in rows.index]]
        expected = st1b[st1b.Source == gse].groupby(['Strain', 'Sex', 'Tissue', 'Intervention'])['Number of samples'].sum()
        observed = joined.groupby(['Strain', 'Sex', 'Tissue', 'Intervention']).size()
        if not observed.sort_index().equals(expected.sort_index().astype(observed.dtype)):
            problems.append('mapped counts differ from Supplementary Table 1B')
        if problems:
            log.append({'set': 'A1', 'gse': gse, 'excluded': problems})
            continue
        valid.append(gse)
        for sample, gsm in zip(mapped['sample'], mapped.gsm):
            row = rows.loc[sample]
            labels[gsm] = {'gse': gse, 'intervention': str(row['Intervention']), 'type': str(row['Intervention.type']),
                           'tissue': str(row['Tissue.details']), 'sex': str(row['Sex']), 'strain': str(row['Strain']),
                           'age': float(row['Age.days'])}
    return labels, valid


def a2_labels(root):
    table = pd.read_csv(root / 'a2/labels.csv', dtype=str).fillna('')
    return {r.gsm: {'gse': r.dataset, 'condition': r.condition, 'class': r['class'], 'tissue': r.tissue,
                    'sex': r.sex, 'age': r.age} for _, r in table.iterrows()}, sorted(set(table.dataset))


def run(root, out, a2_control_classes):
    log = []
    comparisons = [json.loads(line) for line in (root / 'agent_comparisons_2026-10-06.jsonl').read_text().splitlines()]
    sets = {}
    labels1, series1 = a1_labels(root, log)
    sets['A1'] = ('mouse RNA-seq', labels1, set(series1), A1_FACTORS, 'type', {'Control'}, 'age')
    if (root / 'a2/labels.csv').exists():
        labels2, series2 = a2_labels(root)
        sets['A2'] = ('human DNAm', labels2, set(series2), A2_FACTORS, 'class', set(a2_control_classes), None)
    classified = []
    for name, (modality, labels, series, factors, control_factor, control_values, age_factor) in sets.items():
        by_series = {}
        for gsm, value in labels.items():
            by_series.setdefault(value['gse'], {})[gsm] = value
        for item in comparisons:
            if item['modality'] != modality or item['gse'] not in series:
                continue
            entry = {'set': name, **{k: item[k] for k in ('comparison_id', 'gse', 'comparison_name', 'condition_name',
                                                           'membership_source', 'membership_status')},
                     'final_classes': sorted({c['final_class'] for c in item['clocks'] if c['final_class']}),
                     'n_case': len(item['case_gsms']), 'n_control': len(item['control_gsms'])}
            if not item['agent_arms']:
                entry['class'] = 'EXCLUDED_NOT_AGENT_ARMS'
            else:
                result = classify(item['case_gsms'], item['control_gsms'], labels, factors, control_factor,
                                  control_values, age_factor)
                entry.update(result)
                if result['class'] == 'CORRECT':
                    series_labels = by_series.get(item['gse'], {})
                    entry['completeness'] = {arm: completeness(item[f'{arm}_gsms'], labels, series_labels, factors)
                                             for arm in ('case', 'control')}
                    if any(v is not None and v < 0.8 for v in entry['completeness'].values()):
                        entry['flags'].append('incomplete')
                if name == 'A2':
                    entry['cohort'] = {f: {arm: dict(Counter(str(labels[g].get(f)) for g in item[f'{arm}_gsms']
                                                              if g in labels)) for arm in ('case', 'control')}
                                       for f in ('sex',)}
            classified.append(entry)
    rng = random.Random(SEED)
    correct = [e for e in classified if e.get('class') == 'CORRECT']
    sampled = set(e['comparison_id'] for e in rng.sample(correct, max(1, round(0.1 * len(correct))))) if correct else set()
    queue = [dict(e, queue_reason='error_candidate' if e['class'] in ERROR_CANDIDATES else 'correct_sample')
             for e in classified if e.get('class') in ERROR_CANDIDATES or e['comparison_id'] in sampled]
    summary = {name: {'series_used': sorted(sets[name][2]), 'classes': dict(Counter(e['class'] for e in classified if e['set'] == name)),
                      'classes_by_source': {f"{e_src}|{cls}": n for (e_src, cls), n in sorted(Counter(
                          (e['membership_source'], e['class']) for e in classified if e['set'] == name).items(), key=str)}}
               for name in sets}
    summary['excluded_series'] = log
    summary['queue'] = dict(Counter(e['queue_reason'] for e in queue))
    out.mkdir(parents=True, exist_ok=False)
    (out / 'classified.jsonl').write_text(''.join(json.dumps(e) + '\n' for e in classified))
    (out / 'adjudication_queue.jsonl').write_text(''.join(json.dumps(e) + '\n' for e in queue))
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--a2-control-class', nargs='+', default=['HC'])
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.out, args.a2_control_class), indent=1))


if __name__ == '__main__':
    main()
