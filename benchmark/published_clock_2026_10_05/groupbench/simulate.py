"""SIMULATED planted errors for validating the Benchmark A rules (PROTOCOL.md v2), before the agent is scored.

Every output of this module is a simulation, labelled SIMULATED. Copies of real units that the rules classify Correct
receive one planted error each:
  swap_arms          the agent's case arm is set to the arm holding the declared reference level and the reported
                     direction is recomputed from the agent's arms, as the agent's code would (arm order alone is a sign
                     convention); expected REPORTED_RESULT_INVERSION
  merge_levels       samples of another declared level of the series (at least two, a third of the arm) join the case
                     arm; expected a design error
  add_nuisance       samples of another tissue or age group of the series join the case arm; expected a design error
  reverse_labels     the declared case and reference levels are exchanged (names stating the opposite), the reported
                     direction kept; expected REPORTED_RESULT_INVERSION
  drop_half          half of the case arm is removed; expected the incomplete-arm flag
  duplicate_subject  two samples of one arm are given one subject identifier; expected REPEATED_SUBJECT
Seed 20261006; the planted units and outcomes are written with the code version.

  python -m benchmark.published_clock_2026_10_05.groupbench.simulate --root G --prepared DIR --frozen DIR --out DIR
"""
import argparse
import copy
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

from .classify_v2 import ERRORS, classify_unit, holds, orient
from .run_v2 import A1_FACTORS, A2_FACTORS, DECLARED_TO_FACTOR, PredictionReader, prepare_unit, unit_samples

SEED = 20261006
DESIGN_ERRORS = ERRORS - {'REPORTED_RESULT_INVERSION', 'SMALL_ARM_INFERENCE', 'PAIRED_AS_INDEPENDENT'}


def case_arm_key(variant):
    return 'arm_a' if variant['case_is_arm_a'] else 'arm_b'


def recompute_direction(unit):
    """After an arm exchange, the agent's code reports case minus control on the new arms."""
    for variant in unit['variants']:
        case = unit[case_arm_key(variant)]
        control = unit['arm_b' if variant['case_is_arm_a'] else 'arm_a']
        for clock in variant['clocks']:
            values = clock['values']
            a = [values[g] for g in case if g in values]
            b = [values[g] for g in control if g in values]
            if a and b:
                diff = sum(a) / len(a) - sum(b) / len(b)
                clock['historical'] = {**(clock.get('historical') or {}), 'direction': 1.0 if diff > 0 else -1.0,
                                       'effect_size': diff}


def plant(kind, unit, samples, series, rng):
    """Return (planted unit, planted samples) or None when the series cannot host this error."""
    u, s = copy.deepcopy(unit), copy.deepcopy(samples)
    variant = u['variants'][0]
    case_key = case_arm_key(variant)
    if kind == 'swap_arms':
        # The agent's case arm becomes the arm holding the declared reference level, and the reported direction is
        # recomputed from the agent's arms, as the agent's code would; the reported result then contradicts the
        # declared groups whatever the stored arm order was.
        for v in u['variants']:
            if not v['declared'] or v['declared'].get('status') == 'ambiguous':
                return None
            orientation = orient(u, v, s)
            if orientation['case_arm'] is None:
                return None
            v['case_is_arm_a'] = orientation['reference_arm'] == 'arm_a'
        recompute_direction(u)
        return u, s
    if kind == 'reverse_labels':
        for v in u['variants']:
            if v['declared'] and v['declared'].get('status') != 'ambiguous':
                v['declared'] = {**v['declared'], 'case_rule': v['declared']['reference_rule'],
                                 'reference_rule': v['declared']['case_rule']}
                v['series_case_level'], v['series_reference_level'] = v.get('series_reference_level'), v.get('series_case_level')
        return u, s
    if kind == 'drop_half':
        members = sorted(u[case_key])
        if len(members) < 2:
            return None
        u[case_key] = members[: len(members) // 2]
        return u, s
    if kind == 'duplicate_subject':
        members = sorted(u[case_key])
        if len(members) < 2:
            return None
        for g in members[:2]:
            s[g]['subject'] = 'SIMULATED-subject'
        return u, s
    in_arms = set(u['arm_a']) | set(u['arm_b'])
    declared = variant['declared']
    if kind == 'merge_levels':
        pool = [g for g, rec in series.items() if g not in in_arms and rec['expert']
                and not holds(declared['case_rule'], rec)]
    else:  # add_nuisance: another tissue or age group than the case arm holds
        case_tissue = {samples[g]['coded'].get('tissue') or samples[g]['coded'].get('tissue_cell') for g in u[case_key]}
        case_ages = [samples[g]['coded'].get('age') for g in u[case_key] if samples[g]['coded'].get('age') is not None]
        pool = []
        for g, rec in series.items():
            if g in in_arms or not rec['expert']:
                continue
            expert = rec['expert']
            tissue = expert.get('Tissue.details') or expert.get('tissue')
            age = expert.get('Age.days') or expert.get('age')
            other_tissue = tissue not in case_tissue
            other_age = bool(case_ages) and age not in (None, '', 'nan') and \
                abs(float(age) - sorted(case_ages)[len(case_ages) // 2]) > 0.5 * max(sorted(case_ages)[len(case_ages) // 2], 1)
            if other_tissue or other_age:
                pool.append(g)
    need = max(2, round(len(u[case_key]) / 3))
    if len(pool) < need:
        return None
    added = sorted(rng.sample(sorted(pool), need))
    u[case_key] = sorted(u[case_key] + added)
    return u, s


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    units = [json.loads(line) for line in (args.prepared / 'units.jsonl').read_text().splitlines()]
    frozen = {name: json.loads((args.frozen / f'{name}.json').read_text()) for name in ('declared', 'coding', 'refined', 'samples')}
    reader, cache, rows_cache, rng = PredictionReader(), {}, {}, random.Random(SEED)
    outcomes, planted_records = defaultdict(Counter), []
    kinds = ('swap_arms', 'merge_levels', 'add_nuisance', 'reverse_labels', 'drop_half', 'duplicate_subject')
    for unit in units:
        prepared, samples, _ = prepare_unit(unit, frozen, args.root, cache, reader, rows_cache)
        factors = A1_FACTORS if unit['set'] == 'A1' else A2_FACTORS
        mapping = DECLARED_TO_FACTOR[unit['set']]
        factor_of = lambda d: mapping.get(d.get('declared_factor'))
        if classify_unit(prepared, samples, unit['set'], factors, factor_of)['class'] != 'CORRECT':
            continue
        series = cache[(unit['set'], unit['gse'])]
        for kind in kinds:
            planted = plant(kind, prepared, samples, series, rng)
            if planted is None:
                outcomes[kind]['not_applicable'] += 1
                continue
            u, s = planted
            if kind in ('merge_levels', 'add_nuisance', 'drop_half'):
                fresh, _ = unit_samples(u, frozen, args.root, cache)
                s = {**fresh, **{g: v for g, v in s.items() if g in fresh}}
            result = classify_unit(u, s, unit['set'], factors, factor_of)
            types = set(result['error_types'])
            flags = {f for v in result['variants'] for f in v.get('flags', [])}
            detected = {'swap_arms': 'REPORTED_RESULT_INVERSION' in types,
                        'reverse_labels': 'REPORTED_RESULT_INVERSION' in types,
                        'merge_levels': bool(types & DESIGN_ERRORS), 'add_nuisance': bool(types & DESIGN_ERRORS),
                        'drop_half': any(f.startswith('incomplete') for f in flags),
                        'duplicate_subject': 'REPEATED_SUBJECT' in types}[kind]
            outcomes[kind]['detected' if detected else 'missed'] += 1
            planted_records.append({'SIMULATED': True, 'kind': kind, 'source_unit': unit['unit_id'], 'set': unit['set'],
                                    'detected': detected, 'error_types': sorted(types), 'flags': sorted(flags)})
    table = {kind: {**counts, 'detection_rate': (counts['detected'] / (counts['detected'] + counts['missed'])
                                                 if counts['detected'] + counts['missed'] else None)}
             for kind, counts in outcomes.items()}
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out / 'SIMULATED_planted_errors.jsonl').write_text(''.join(json.dumps(r) + '\n' for r in planted_records))
    (args.out / 'SIMULATED_detection.json').write_text(json.dumps({'SIMULATED': True, 'seed': SEED, 'table': table}, indent=1) + '\n')
    print(json.dumps(table, indent=1))


if __name__ == '__main__':
    main()
