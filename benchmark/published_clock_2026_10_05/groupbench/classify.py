"""Classify agent comparisons against expert sample labels (Benchmark A of PROTOCOL.md).

Each comparison's arms are labelled sample by sample with the expert annotation. A factor differs between arms when
the arms' distributions over its values differ by more than TV_THRESHOLD in total variation distance. The rules, in
order: not assessable (label coverage below COVERAGE in an arm), no contrast, confounded (more than one factor
differs), mixed arm (the contrast factor takes several values within an arm, except pooled doses of one compound and
pooled non-overlapping ages), possible inversion (the case arm is the expert control or the younger group), correct.
Error candidates go to adjudication; this module decides nothing beyond the rules.
"""
from collections import Counter

COVERAGE = 0.8
TV_THRESHOLD = 0.2
INCOMPLETE = 0.8
ERROR_CANDIDATES = {'NO_CONTRAST', 'CONFOUNDED', 'MIXED_ARM', 'POSSIBLE_INVERSION'}


def total_variation(a, b):
    """Total variation distance between the value distributions of two lists."""
    ca, cb = Counter(a), Counter(b)
    return 0.5 * sum(abs(ca[v] / len(a) - cb[v] / len(b)) for v in set(ca) | set(cb))


def compound(value):
    """The compound of an intervention label ('Rapamycin; 14ppm' -> 'Rapamycin')."""
    return str(value).split(';')[0].strip().lower()


def classify(case, control, labels, factors, control_factor, control_values, age_factor=None):
    """Classify one comparison.

    case, control: GSM lists. labels: {gsm: {factor: value}}. factors: design factors to compare (ordered).
    control_factor: the factor whose values mark the expert control (for example intervention type), with
    control_values the values meaning control. age_factor: the factor holding numeric age (pooled-age rule).
    """
    result = {'n_case': len(case), 'n_control': len(control), 'flags': []}
    labelled = {arm: [g for g in gsms if g in labels] for arm, gsms in (('case', case), ('control', control))}
    result['coverage'] = {arm: (len(labelled[arm]) / len(gsms) if gsms else 0.0)
                          for arm, gsms in (('case', case), ('control', control))}
    if min(result['coverage'].values()) < COVERAGE:
        return {**result, 'class': 'NOT_ASSESSABLE'}
    values = {f: {arm: [labels[g].get(f) for g in labelled[arm]] for arm in labelled} for f in factors}
    result['tv'] = {f: round(total_variation(values[f]['case'], values[f]['control']), 3) for f in factors}
    differing = [f for f in factors if result['tv'][f] > TV_THRESHOLD]
    result['differing'] = differing
    result['arm_values'] = {f: {arm: dict(Counter(map(str, v))) for arm, v in values[f].items()} for f in factors}
    if not differing:
        return {**result, 'class': 'NO_CONTRAST'}
    if len(differing) > 1:
        return {**result, 'class': 'CONFOUNDED'}
    factor = differing[0]
    result['contrast_factor'] = factor
    for arm in ('case', 'control'):
        distinct = set(map(str, values[factor][arm]))
        if len(distinct) <= 1:
            continue
        if factor == age_factor:
            continue  # judged below: pooled ages are acceptable when the arms do not overlap
        if factor != control_factor and len({compound(v) for v in distinct}) == 1:
            result['flags'].append(f'pooled_doses_{arm}')
            continue
        return {**result, 'class': 'MIXED_ARM', 'mixed_arm': arm}
    if factor == age_factor:
        case_ages = [float(v) for v in values[factor]['case']]
        control_ages = [float(v) for v in values[factor]['control']]
        if not (max(case_ages) < min(control_ages) or max(control_ages) < min(case_ages)):
            return {**result, 'class': 'MIXED_ARM', 'mixed_arm': 'overlapping ages'}
        if len(set(case_ages)) > 1 or len(set(control_ages)) > 1:
            result['flags'].append('pooled_ages')
        if min(case_ages) < min(control_ages):
            return {**result, 'class': 'POSSIBLE_INVERSION', 'inversion': 'case arm is the younger group'}
    else:
        is_control = {arm: [labels[g].get(control_factor) in control_values for g in labelled[arm]]
                      for arm in labelled}
        if all(is_control['case']) and not any(is_control['control']):
            return {**result, 'class': 'POSSIBLE_INVERSION', 'inversion': 'case arm is the expert control'}
        if not any(is_control['control']):
            result['flags'].append('reference_is_not_expert_control')
    for f in factors:
        if f != factor and len({str(v) for arm in values[f].values() for v in arm}) > 1:
            result['flags'].append(f'pooled_{f}')
    return {**result, 'class': 'CORRECT'}


def completeness(arm, labels, all_labels, factors):
    """Share of the expert stratum's samples (in the series) that the arm includes, when the arm is one stratum."""
    strata = {tuple(labels[g].get(f) for f in factors) for g in arm if g in labels}
    if len(strata) != 1:
        return None
    stratum = strata.pop()
    members = [g for g, value in all_labels.items() if tuple(value.get(f) for f in factors) == stratum]
    return len(set(arm) & set(members)) / len(members) if members else None
