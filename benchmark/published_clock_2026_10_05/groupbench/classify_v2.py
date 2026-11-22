"""Rules of PROTOCOL.md v2 for one unit (a distinct pair of arms) and its name variants.

Inputs per unit: the two arms; per sample, coded nuisance factors and the fields the declared rules read (expert
labels and GEO values); per variant, the declared contrast (rules selecting the declared case and reference levels),
which arm the agent made the case, the recorded test and, per clock, the stored direction, effect and P and the
per-sample clock values. The functions return rule outcomes; adjudication decides afterwards.
"""
import math
import random
import statistics
from collections import Counter

PURITY = 0.9
TV = 0.2
MIN_MOVE = 2
PERMUTATIONS = 2000
SEED = 20261006
UNASSIGNABLE = 0.2
INCOMPLETE = 0.8
AGE_RELATIVE = 0.25
A2_AGE_YEARS, A2_AGE_SMD = 5.0, 0.5
RATE_CLOCKS = {'DunedinPACE'}
# Cohort factors can be imbalanced by recruitment or chance, so a partial imbalance needs a permutation P < 0.05.
# Every other factor is set by the investigators' design; a partial imbalance of at least MIN_MOVE samples counts.
COHORT_FACTORS = {'sex'}
ERRORS = {'DEGENERATE_TEST', 'WRONG_FACTOR', 'LABEL_MISMATCH', 'IMPURE_ARM', 'CONFOUNDED', 'CONFOUNDED_AGE',
          'REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT', 'SMALL_ARM_INFERENCE', 'REPORTED_RESULT_INVERSION'}
DESIGN = {'DEGENERATE_TEST', 'WRONG_FACTOR', 'LABEL_MISMATCH', 'IMPURE_ARM', 'CONFOUNDED', 'CONFOUNDED_AGE',
          'REPEATED_SUBJECT'}


def holds(rule, sample):
    """A declared-level rule: OR over groups, AND over predicates; a predicate holds on an exact field value."""
    for group in rule or []:
        if all(str(sample.get(p['source'], {}).get(p['field'])) in set(map(str, p['values'])) for p in group):
            return True
    return False


def total_variation(a, b):
    ca, cb = Counter(a), Counter(b)
    return 0.5 * sum(abs(ca[v] / len(a) - cb[v] / len(b)) for v in set(ca) | set(cb))


def permutation_p(a, b, statistic, seed=SEED):
    observed, pooled, n = statistic(a, b), list(a) + list(b), len(a)
    rng, hits = random.Random(seed), 0
    for _ in range(PERMUTATIONS):
        rng.shuffle(pooled)
        hits += statistic(pooled[:n], pooled[n:]) >= observed - 1e-12
    return (hits + 1) / (PERMUTATIONS + 1)


def median_gap(a, b):
    ma, mb = statistics.median(a), statistics.median(b)
    return abs(ma - mb) / max(abs(ma), abs(mb), 1e-9)


def smd(a, b):
    if len(a) < 2 or len(b) < 2:
        return 0.0
    pooled = math.sqrt((statistics.variance(a) + statistics.variance(b)) / 2)
    return abs(statistics.mean(a) - statistics.mean(b)) / pooled if pooled else 0.0


def nuisance(arm_a, arm_b, samples, factors, set_name, exempt):
    """Nuisance-factor findings: complete separation or a significant imbalance that needs at least two samples moved."""
    findings = []
    for factor in factors:
        if factor in exempt:
            continue
        a = [samples[g]['coded'].get(factor) for g in arm_a]
        b = [samples[g]['coded'].get(factor) for g in arm_b]
        a, b = [v for v in a if v is not None], [v for v in b if v is not None]
        if len(a) < 1 or len(b) < 1:
            continue
        if factor == 'age':
            a, b = [float(v) for v in a], [float(v) for v in b]
            if set_name == 'A2':
                gap = abs(statistics.median(a) - statistics.median(b))
                if gap > A2_AGE_YEARS or smd(a, b) > A2_AGE_SMD:
                    findings.append({'factor': 'age', 'kind': 'age_imbalance', 'median_gap_years': round(gap, 2),
                                     'smd': round(smd(a, b), 3)})
                continue
            complete = max(a) < min(b) or max(b) < min(a)
            if median_gap(a, b) > AGE_RELATIVE:  # rodent age groups are set by design
                findings.append({'factor': 'age', 'kind': 'complete' if complete else 'imbalance',
                                 'median_gap_relative': round(median_gap(a, b), 3)})
            continue
        tv = total_variation(a, b)
        if tv <= TV:
            continue
        complete = not (set(a) & set(b))
        moved = tv * min(len(a), len(b))
        significant = factor not in COHORT_FACTORS or permutation_p(a, b, total_variation) < 0.05
        if complete or (moved >= MIN_MOVE - 1e-9 and significant):
            findings.append({'factor': factor, 'kind': 'complete' if complete else 'imbalance', 'tv': round(tv, 3)})
        else:
            findings.append({'factor': factor, 'kind': 'minor_imbalance', 'tv': round(tv, 3)})
    return findings


def orient(unit, variant, samples):
    """Which arm carries the declared case and reference levels, with purity; None when neither orientation holds."""
    purity = {}
    for arm in ('arm_a', 'arm_b'):
        members = unit[arm]
        purity[arm] = {level: sum(holds(variant['declared'][f'{level}_rule'], samples[g]) for g in members) / len(members)
                       for level in ('case', 'reference')}
    straight = min(purity['arm_a']['case'], purity['arm_b']['reference'])
    crossed = min(purity['arm_b']['case'], purity['arm_a']['reference'])
    if max(straight, crossed) >= PURITY:
        case_arm = 'arm_a' if straight >= crossed else 'arm_b'
        return {'case_arm': case_arm, 'reference_arm': 'arm_b' if case_arm == 'arm_a' else 'arm_a', 'purity': purity}
    return {'case_arm': None, 'purity': purity}


def clock_sign(values, case, reference):
    a = [values[g] for g in case if g in values and values[g] == values[g]]
    b = [values[g] for g in reference if g in values and values[g] == values[g]]
    if not a or not b:
        return None
    diff = statistics.mean(a) - statistics.mean(b)
    return 1 if diff > 0 else -1 if diff < 0 else 0


def classify_unit(unit, samples, set_name, factors, factor_of_declared):
    """All rule findings for one unit. samples: {gsm: {'expert': {...}, 'geo': {...}, 'coded': {...},
    'category': labelled|same_group|repeat_subject|other_group|unassignable, 'subject': id or None}}."""
    out = {'unit_id': unit['unit_id'], 'findings': [], 'variants': []}
    arms = {arm: unit[arm] for arm in ('arm_a', 'arm_b')}
    if set(arms['arm_a']) == set(arms['arm_b']) or not arms['arm_a'] or not arms['arm_b']:
        out['findings'].append({'type': 'DEGENERATE_TEST', 'detail': 'identical or empty arms'})
    for arm, members in arms.items():
        share = sum(samples[g]['category'] == 'unassignable' for g in members) / len(members) if members else 1
        if share > UNASSIGNABLE:
            out['class'] = 'NOT_ASSESSABLE'
            out['reason'] = f'{arm}: {share:.0%} of samples unassignable'
            return out
    for arm, members in arms.items():
        subjects = [samples[g]['subject'] for g in members if samples[g].get('subject')]
        repeated = [s for s, n in Counter(subjects).items() if n > 1]
        if repeated:
            out['findings'].append({'type': 'REPEATED_SUBJECT', 'arm': arm, 'subjects': repeated})
    both = {samples[g]['subject'] for g in arms['arm_a'] if samples[g].get('subject')} & {
        samples[g]['subject'] for g in arms['arm_b'] if samples[g].get('subject')}
    all_subjects = {samples[g]['subject'] for g in arms['arm_a'] + arms['arm_b'] if samples[g].get('subject')}
    paired = bool(all_subjects) and len(both) >= 0.5 * len(all_subjects)
    for variant in unit['variants']:
        declared = variant['declared']
        record = {'variant': variant['index'], 'findings': []}
        if declared is None or declared.get('status') == 'ambiguous' or declared['declared_factor'] == 'unclear':
            record['class'] = 'AMBIGUOUS_LABEL'
            out['variants'].append(record)
            continue
        fields = {(p['source'], p['field']) for level in ('case', 'reference')
                  for group in declared[f'{level}_rule'] or [] for p in group}
        recorded = {arm: sum(any(samples[g].get(src, {}).get(f) not in (None, '') for src, f in fields)
                             for g in members) / len(members) for arm, members in arms.items()}
        if not fields or min(recorded.values()) < 1 - UNASSIGNABLE:
            record['class'] = 'NOT_ASSESSABLE_UNENCODED'
            out['variants'].append(record)
            continue
        orientation = orient(unit, variant, samples)
        record['purity'] = orientation['purity']
        declared_factor = factor_of_declared(declared)
        if orientation['case_arm'] is None:
            separating = [f for f in nuisance(arms['arm_a'], arms['arm_b'], samples, factors, set_name, {declared_factor})
                          if f['kind'] == 'complete']
            if separating:
                record['findings'].append({'type': 'WRONG_FACTOR', 'factors': separating})
            elif max(orientation['purity']['arm_a'].values()) >= PURITY or max(orientation['purity']['arm_b'].values()) >= PURITY:
                record['findings'].append({'type': 'IMPURE_ARM', 'purity': orientation['purity']})
            else:
                record['findings'].append({'type': 'LABEL_MISMATCH', 'purity': orientation['purity']})
        else:
            case, reference = arms[orientation['case_arm']], arms[orientation['reference_arm']]
            for label, members, level in (('case', case, 'series_case_level'), ('reference', reference, 'series_reference_level')):
                pool = set(variant.get(level) or [])
                if pool and len(set(members) & pool) / len(pool) < INCOMPLETE:
                    record.setdefault('flags', []).append(f'incomplete {label} arm')
            exempt = {declared_factor}
            if variant.get('adjusted'):
                exempt |= {'sex', 'age'}
            complete_groups = (set(case) == set(variant.get('series_case_level', [])) and
                               set(reference) == set(variant.get('series_reference_level', [])))
            for finding in nuisance(case, reference, samples, factors, set_name, exempt):
                if finding['kind'] == 'minor_imbalance':
                    record.setdefault('flags', []).append(f"minor imbalance in {finding['factor']}")
                elif complete_groups:
                    record.setdefault('flags', []).append(f"study-design imbalance in {finding['factor']}")
                elif finding['kind'] == 'age_imbalance':
                    raw_clocks = [c['clock'] for c in variant['clocks'] if c['clock'] not in RATE_CLOCKS]
                    if raw_clocks:
                        record['findings'].append({'type': 'CONFOUNDED_AGE', 'clocks': raw_clocks, **finding})
                else:
                    record['findings'].append({'type': 'CONFOUNDED', **finding})
            agent_case = unit['arm_a'] if variant['case_is_arm_a'] else unit['arm_b']
            for clock in variant['clocks']:
                historical = clock.get('historical') or {}
                if historical.get('pvalue') is not None and min(len(case), len(reference)) < 3:
                    record['findings'].append({'type': 'SMALL_ARM_INFERENCE', 'clock': clock['clock']})
                if historical.get('pvalue') == 1 and historical.get('effect_size') == 0:
                    record['findings'].append({'type': 'DEGENERATE_TEST', 'clock': clock['clock'], 'detail': 'P = 1, effect 0'})
                values = clock.get('values') or {}
                sign = clock_sign(values, case, reference)
                direction = historical.get('direction')
                if sign is None or direction is None:
                    continue
                if sign != 0 and math.copysign(1, direction) != sign:
                    record['findings'].append({'type': 'REPORTED_RESULT_INVERSION', 'clock': clock['clock'],
                                               'declared_sign': sign, 'reported_direction': direction,
                                               'agent_case_is_declared_case': set(agent_case) == set(case)})
                elif historical.get('effect_size') is not None and sign != 0 and historical['effect_size'] != 0 and \
                        math.copysign(1, historical['effect_size']) != sign:
                    record.setdefault('record_errors', []).append({'type': 'EFFECT_SIGN_INCONSISTENCY', 'clock': clock['clock']})
            test = ' '.join(str((c.get('historical') or {}).get('statistical_test') or '') for c in variant['clocks']).lower()
            if paired and 'paired' not in test and 'mixed' not in test:
                record['findings'].append({'type': 'PAIRED_AS_INDEPENDENT', 'subjects_in_both_arms': len(both)})
        record['class'] = 'ERROR' if record['findings'] else 'CORRECT'
        out['variants'].append(record)
    errors = [f for f in out['findings']] + [f for v in out['variants'] for f in v.get('findings', [])]
    classes = {v['class'] for v in out['variants']}
    if errors:
        out['class'] = 'ERROR'
    elif classes and classes <= {'NOT_ASSESSABLE_UNENCODED'}:
        out['class'] = 'NOT_ASSESSABLE_UNENCODED'
    elif classes and classes <= {'AMBIGUOUS_LABEL', 'NOT_ASSESSABLE_UNENCODED'}:
        out['class'] = 'AMBIGUOUS_LABEL'
    else:
        out['class'] = 'CORRECT'
    out['error_types'] = sorted({f['type'] for f in errors})
    return out
