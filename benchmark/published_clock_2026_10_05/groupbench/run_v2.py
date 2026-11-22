"""Run Benchmark A of PROTOCOL.md v2 with the frozen inputs.

E1: each unit's own arms and the agent's reported direction. E2: each released row (VERIFIED or PROVISIONAL) on its
post-repair arms, with the direction of the verifier's recomputation. Writes rule results, the E2 results and the
blinded census queue (all assessable units, shuffled; no rule outcome and no verifier output in the packages).

  python -m benchmark.published_clock_2026_10_05.groupbench.run_v2 --root GROUPBENCH --prepared DIR --frozen DIR --out DIR
"""
import argparse
import json
import random
import statistics
from collections import Counter
from pathlib import Path

import pandas as pd

from ..engine.cli import atlas_rows
from ..engine.predictions import PredictionReader
from .classify_v2 import classify_unit, holds
from .reconcile_inputs import series_samples

SEED = 20261006
MODALITY = {'A1': 'mouse RNA-seq', 'A2': 'human DNAm'}
A1_FACTORS = ['compound', 'genotype', 'diet', 'exercise', 'hypoxia', 'other', 'tissue', 'sex', 'background', 'age']
A2_FACTORS = ['condition', 'tissue_cell', 'sex', 'age']
DECLARED_TO_FACTOR = {'A1': {'intervention_compound': 'compound', 'intervention_genotype': 'genotype',
                             'intervention_diet': 'diet', 'intervention_other': 'other', 'tissue': 'tissue',
                             'cell_type': 'tissue', 'sex': 'sex', 'background_strain': 'background', 'age_group': 'age'},
                      'A2': {'condition': 'condition', 'tissue': 'tissue_cell', 'cell_type': 'tissue_cell', 'sex': 'sex',
                             'age_group': 'age'}}
ADJUSTED_WORDS = ('adjust', 'covariat', 'linear model', 'ancova', 'regression', 'mixed')
MISSING_SEX = {'Unknown', 'U', '', 'nan', 'None'}


def coded_a1(expert, coding):
    entry = coding['interventions'].get(expert.get('Intervention')) or {'components': []}
    names = {}
    for component in entry['components']:
        names.setdefault(component['class'], []).append(component['name'].strip().lower())
    out = {k: '+'.join(sorted(names.get(k, []))) or 'none' for k in ('compound', 'genotype', 'diet', 'exercise', 'hypoxia', 'other')}
    out['tissue'] = expert.get('Tissue.details')
    out['sex'] = None if str(expert.get('Sex')) in MISSING_SEX else expert.get('Sex')
    out['background'] = coding['backgrounds'].get(expert.get('Strain'), expert.get('Strain'))
    out['age'] = float(expert['Age.days']) if expert.get('Age.days') not in (None, '', 'nan') else None
    return out


def coded_a2(expert, refined):
    age = expert.get('age')
    return {'condition': refined or expert.get('condition'),
            'tissue_cell': f"{expert.get('tissue')}|{expert.get('cell_type', '')}",
            'sex': None if str(expert.get('sex')) in MISSING_SEX else expert.get('sex'),
            'age': float(age) if age not in (None, '', 'nan') else None}


def unit_samples(unit, frozen, root, cache):
    """Per-sample records for a unit's arms, with the frozen categories, subjects and label-audit removals."""
    set_name, gse = unit['set'], unit['gse']
    if (set_name, gse) not in cache:
        cache[(set_name, gse)] = series_samples(root, set_name, gse)
    series = cache[(set_name, gse)]
    checks = frozen['samples'].get(gse, {'unlabelled': {}, 'subjects': {}, 'conflicts': []})
    removed = {(c['gsm'], c['field']) for c in checks.get('conflicts', [])}
    samples = {}
    for arm in ('arm_a', 'arm_b'):
        for gsm in unit[arm]:
            base = series.get(gsm, {'geo': {}, 'expert': {}})
            expert = {k: (None if (gsm, k) in removed else v) for k, v in base['expert'].items()}
            record = {'geo': base['geo'], 'expert': expert, 'subject': (checks.get('subjects') or {}).get(gsm) or None,
                      'arm': arm}
            if expert:
                record['category'] = 'labelled'
            else:
                record['category'] = (checks['unlabelled'].get(gsm) or {}).get('category', 'unassignable')
            samples[gsm] = record
    for gsm, record in samples.items():
        if record['category'] == 'same_group':  # the curator excluded it from a group the arm already holds
            peers = [samples[g]['expert'] for g in unit[record['arm']] if samples[g]['category'] == 'labelled']
            if peers:
                modal = Counter(json.dumps(p, sort_keys=True) for p in peers).most_common(1)[0][0]
                record['expert'] = json.loads(modal)
                record['flags'] = ['curator_excluded']
    for gsm, record in samples.items():
        expert = record['expert'] or {}
        if not expert:
            record['coded'] = {}
        elif set_name == 'A1':
            record['coded'] = coded_a1(expert, frozen['coding'])
        else:
            refined = next((v for k, v in frozen['refined'].items() if k.endswith(f'|{gsm}')), None)
            record['coded'] = coded_a2(expert, refined)
    return samples, series


def clock_values(reader, modality, gse, clock, rows):
    loaded = reader.load(modality, gse, clock, rows)
    if loaded.get('status') != 'OK':
        return {}
    return {str(k): float(v) for k, v in loaded['samples']['clock'].dropna().items()}


def prepare_unit(unit, frozen, root, cache, reader, rows_cache):
    samples, series = unit_samples(unit, frozen, root, cache)
    modality = MODALITY[unit['set']]
    if unit['gse'] not in rows_cache:
        rows_cache[unit['gse']] = atlas_rows(modality, unit['gse'])
    rows = rows_cache[unit['gse']]
    variants = []
    for index, variant in enumerate(unit['variants']):
        declared = frozen['declared'].get(f"{unit['unit_id']}#{index}")
        test = ' '.join(str((c.get('historical') or {}).get('statistical_test') or '') for c in variant['clocks']).lower()
        clocks = [{**c, 'values': clock_values(reader, modality, unit['gse'], c['clock'], rows)} for c in variant['clocks']]
        entry = {'index': index, 'declared': declared, 'case_is_arm_a': variant['case_is_arm_a'], 'clocks': clocks,
                 'adjusted': any(w in test for w in ADJUSTED_WORDS)}
        if declared and declared.get('status') != 'ambiguous':
            entry['series_case_level'] = [g for g, s in series.items() if holds(declared['case_rule'], s)]
            entry['series_reference_level'] = [g for g, s in series.items() if holds(declared['reference_rule'], s)]
        variants.append(entry)
    return {**unit, 'variants': variants}, samples, rows


def released_units(unit, prepared, rows):
    """E2: one unit per released row, on its post-repair arms with the direction of the verifier's recomputation."""
    by_ordinal = {r['row_ordinal']: r for r in rows}
    out = []
    for variant in prepared['variants']:
        for clock in variant['clocks']:
            row = by_ordinal.get(clock['row_ordinal'])
            if not row or row.get('final_class') not in ('VERIFIED', 'PROVISIONAL'):
                continue
            recomputed = row.get('recomputed') or {}
            mean_diff = recomputed.get('mean_diff')
            if not row.get('case_gsms') or not row.get('control_gsms') or mean_diff is None:
                continue
            released_clock = {'clock': clock['clock'], 'row_ordinal': clock['row_ordinal'], 'values': clock['values'],
                              'historical': {'pvalue': recomputed.get('p'), 'effect_size': mean_diff,
                                             'direction': 1.0 if mean_diff > 0 else -1.0 if mean_diff < 0 else None,
                                             'statistical_test': 'verifier recomputation'}}
            out.append({'unit_id': f"{prepared['unit_id']}|E2|{clock['row_ordinal']}", 'set': prepared['set'],
                        'gse': prepared['gse'], 'arm_a': sorted(row['case_gsms']), 'arm_b': sorted(row['control_gsms']),
                        'variants': [{**variant, 'case_is_arm_a': True, 'clocks': [released_clock]}],
                        'final_class': row['final_class'], 'actions': [a.get('action') for a in row.get('actions') or []]})
    return out


def census_package(prepared, samples):
    """What reviewers see: names, arms, per-sample labels, reported results with arm means; nothing else."""
    variants = []
    for variant in prepared['variants']:
        original = prepared['raw_variants'][variant['index']]
        agent_case = prepared['arm_a'] if variant['case_is_arm_a'] else prepared['arm_b']
        reported = []
        for clock in variant['clocks']:
            values, historical = clock['values'], clock.get('historical') or {}
            mean = lambda arm: statistics.mean([values[g] for g in arm if g in values]) if any(g in values for g in arm) else None
            reported.append({'clock': clock['clock'], 'test': historical.get('statistical_test'),
                             'reported_direction': historical.get('direction'), 'reported_p': historical.get('pvalue'),
                             'mean_arm_1': mean(prepared['arm_a']), 'mean_arm_2': mean(prepared['arm_b'])})
        variants.append({'variant': variant['index'], 'comparison_name': original['comparison_name'],
                         'condition_name': original['condition_name'],
                         'agent_case_arm': 'arm 1' if variant['case_is_arm_a'] else 'arm 2', 'reported': reported,
                         'agent_case_size': len(agent_case)})
    rows = [{'gsm': g, 'arm': 'arm 1' if s['arm'] == 'arm_a' else 'arm 2', 'expert_label': s['expert'] or None,
             'unlabelled_category': None if s['category'] == 'labelled' else s['category'], 'subject': s['subject'],
             'geo': s['geo']} for g, s in sorted(samples.items(), key=lambda x: (x[1]['arm'], x[0]))]
    return {'set': prepared['set'], 'gse': prepared['gse'], 'variants': variants, 'samples': rows}


def run(root, prepared_dir, frozen_dir, out):
    units = [json.loads(line) for line in (prepared_dir / 'units.jsonl').read_text().splitlines()]
    frozen = {name: json.loads((frozen_dir / f'{name}.json').read_text()) for name in ('declared', 'coding', 'refined', 'samples')}
    reader, cache, rows_cache = PredictionReader(), {}, {}
    results, e2, census = [], [], []
    for unit in units:
        prepared, samples, rows = prepare_unit(unit, frozen, root, cache, reader, rows_cache)
        prepared['raw_variants'] = unit['variants']
        factors = A1_FACTORS if unit['set'] == 'A1' else A2_FACTORS
        mapping = DECLARED_TO_FACTOR[unit['set']]
        factor_of = lambda d: mapping.get(d.get('declared_factor'))
        result = classify_unit(prepared, samples, unit['set'], factors, factor_of)
        result.update(set=unit['set'], gse=unit['gse'], n_arm_a=len(unit['arm_a']), n_arm_b=len(unit['arm_b']),
                      membership_sources=sorted({v['membership_source'] for v in unit['variants']}),
                      orientation_conflict=unit['orientation_conflict'])
        results.append(result)
        for released in released_units(unit, prepared, rows):
            e2_result = classify_unit(released, samples_for(released, frozen, root, cache), unit['set'], factors, factor_of)
            e2_result.update(set=unit['set'], gse=unit['gse'], parent=unit['unit_id'],
                             final_class=released['final_class'], actions=released['actions'])
            e2.append(e2_result)
        if result['class'] in ('ERROR', 'CORRECT'):
            census.append({'unit_id': unit['unit_id'], 'package': census_package(prepared, samples)})
    random.Random(SEED).shuffle(census)
    for position, item in enumerate(census):
        item['queue_position'] = position
    out.mkdir(parents=True, exist_ok=False)
    (out / 'e1_rules.jsonl').write_text(''.join(json.dumps(r, default=str) + '\n' for r in results))
    (out / 'e2_rules.jsonl').write_text(''.join(json.dumps(r, default=str) + '\n' for r in e2))
    (out / 'census_queue.jsonl').write_text(''.join(json.dumps(c, default=str) + '\n' for c in census))
    summary = {'E1': {s: dict(Counter(r['class'] for r in results if r['set'] == s)) for s in ('A1', 'A2')},
               'E1_error_types': {s: dict(Counter(t for r in results if r['set'] == s for t in r['error_types'])) for s in ('A1', 'A2')},
               'E2_rows': {s: dict(Counter(r['class'] for r in e2 if r['set'] == s)) for s in ('A1', 'A2')},
               'census_units': len(census)}
    (out / 'summary_rules.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


def samples_for(released, frozen, root, cache):
    samples, _ = unit_samples({'set': released['set'], 'gse': released['gse'], 'arm_a': released['arm_a'],
                               'arm_b': released['arm_b']}, frozen, root, cache)
    return samples


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.prepared, args.frozen, args.out), indent=1))


if __name__ == '__main__':
    main()
