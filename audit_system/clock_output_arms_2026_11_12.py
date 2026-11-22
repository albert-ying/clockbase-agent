"""Apply checked comparison decisions at the verifier's final step.

Detection helpers document candidate selection only. Labels, saved code and
arm separation never decide a final class. The checked decision list supplies
the class and reason for every withheld comparison.
"""
import json
import re
from pathlib import Path

import pandas as pd
import numpy as np

REASON = 'ARMS_DEFINED_BY_CLOCK_OUTPUT'
MEANING = 'arms defined by a clock output'
INFERENTIAL = ('VERIFIED', 'PROVISIONAL')
CLOCK = re.compile(r'\bclock\b|scaleddiff|Horvath|PhenoAge|GrimAge|Dunedin|YingDam|YingAdapt|\bKnight\b|age acceleration', re.I)
SPLIT = re.compile(r'\bhigh(?:est)?\b|\blow(?:est)?\b|accelerated|decelerated|quartile|tertile|outlier|older|younger|cluster', re.I)
ACCEL = re.compile(r'accelerated.*(?:decelerated|normal).*aging|accelerated.*aging.*cluster', re.I)
SELECTION = re.compile(r'quantile\(|percentile\(|median\(|qcut\(|pd.cut\(|\.apply\(|np.where\(|sort_values\(|nlargest\(|nsmallest\(|fit_predict\(|KMeans\(|(?:\[clock\]|clock_values|\[clock_name\]).*[<>]|_category.*=', re.I)


def labels(record):
    return ' '.join(str(record.get(k) or '') for k in ('comparison_name', 'condition_name'))


def label_candidate(record):
    text = labels(record)
    return bool(CLOCK.search(text) or ACCEL.search(text))


def code_evidence(trace):
    """Keep source blocks with selection expressions and clock input context."""
    blocks = trace.get('blocks') or []
    return [{'block': i, 'code': b} for i, b in enumerate(blocks)
            if SELECTION.search(b) and (CLOCK.search(b) or re.search(r'\[clock\]|clock_values|clock_name|scaled_data|X_scaled', b))]


def arm_evidence(record, predictions):
    inputs = record.get('inputs') or {}
    case, control = inputs.get('case_gsms') or [], inputs.get('control_gsms') or []
    result = {'case_gsms': case, 'control_gsms': control,
              'shared_gsms': sorted(set(case) & set(control)), 'clocks': {}}
    if predictions is None or not case or not control:
        return result
    id_col = next(c for c in predictions if c.lower() in ('gsm_id', 'gsm', 'sample_id'))
    table = predictions.set_index(id_col)
    table = table.loc[~table.index.duplicated(keep=False)]
    for clock in table.columns:
        if not CLOCK.search(clock) and clock != record['clock']:
            continue
        values = pd.to_numeric(table[clock], errors='coerce')
        a, b = values.reindex(case).dropna(), values.reindex(control).dropna()
        if not len(a) or not len(b):
            continue
        # Pairwise ordering is descriptive evidence, not a significance test.
        ordered = np.sort(b.to_numpy())
        left = np.searchsorted(ordered, a.to_numpy(), side='left')
        right = np.searchsorted(ordered, a.to_numpy(), side='right')
        below = int((len(b) - right).sum())
        tied = int((right - left).sum())
        result['clocks'][clock] = {
            'n_case': len(a), 'n_control': len(b),
            'case_min': float(a.min()), 'case_max': float(a.max()), 'case_mean': float(a.mean()),
            'control_min': float(b.min()), 'control_max': float(b.max()), 'control_mean': float(b.mean()),
            'range_overlap': bool(max(a.min(), b.min()) <= min(a.max(), b.max())),
            'strict_separation': bool(a.min() > b.max() or b.min() > a.max()),
            'case_below_control_fraction': below / (len(a) * len(b)),
            'ties_fraction': tied / (len(a) * len(b))}
    return result


def assess(record, trace, predictions):
    """Return a candidate decision; source review remains explicit in the table.

    A labelled clock split plus recovered arms is sufficient when the saved code
    establishes output selection. Without code, require strict separation on the
    tested clock. Cross-clock splits do not require separation on the tested
    clock when code identifies the selecting clock.
    """
    text = labels(record)
    candidate = label_candidate(record)
    if not candidate:
        return None
    evidence = code_evidence(trace or {})
    arms = arm_evidence(record, predictions)
    decision = {'row_ordinal': record['row_ordinal'], 'gse_id': record['gse_id'],
                'comparison_name': record['comparison_name'], 'condition_name': record['condition_name'],
                'clock': record['clock'], 'old_class': record['final_class'],
                'trace_available': bool(trace and trace.get('blocks')),
                'code_evidence': evidence, 'arm_evidence': arms}
    if record['final_class'] not in INFERENTIAL + ('DUPLICATE',):
        decision.update(decision='OUT_OF_SCOPE', reason='already withheld or descriptive')
        return decision
    explicit = bool((CLOCK.search(text) and SPLIT.search(text)) or ACCEL.search(text))
    # Merely naming the tested clock in a biological comparison is insufficient.
    # CLOCK gene perturbations likewise do not select samples by a clock output.
    explicit = explicit and not re.search(r'siCLOCK|sh-CLOCK', text, re.I)
    selected = bool(evidence and explicit)
    tested = arms['clocks'].get(record['clock']) or arms['clocks'].get(record['clock'] + '_mean') or {}
    separated = tested.get('strict_separation', False)
    recovered = bool(arms['case_gsms'] and arms['control_gsms'] and not arms['shared_gsms'])
    keep = recovered and explicit and (selected or separated)
    decision.update(decision='WITHHOLD' if keep else 'RETAIN',
                    reason=MEANING if keep else 'label does not establish output-defined recovered arms',
                    evidence_used='labels and replayed selection code and recovered arms' if selected else
                                  'labels and recovered arm separation on tested clock' if keep else
                                  'labels and available code and arm evidence',
                    tested_clock_range_overlap=tested.get('range_overlap'),
                    tested_clock_strict_separation=separated)
    return decision


REASON_CODES = {MEANING: REASON,
                'arms do not match the GEO sample records': 'ARMS_DO_NOT_MATCH_GEO_SAMPLE_RECORDS'}


def load_decisions(path):
    """Validate comparison identities and preserve the supplied class and reason."""
    source = json.loads(Path(path).read_text())
    withheld = {}
    for group in source.values():
        if not isinstance(group, dict) or 'comparisons' not in group:
            continue
        if group['final_class'] not in ('NOT_REPAIRABLE', 'PENDING_DESIGN'):
            raise ValueError(f'Unsupported decision class: {group["final_class"]}')
        if group['reason'] not in REASON_CODES:
            raise ValueError(f'Unsupported decision reason: {group["reason"]}')
        for key in group['comparisons']:
            if len(key.split('|', 2)) != len(('modality', 'gse_id', 'comparison_name')):
                raise ValueError(f'Invalid comparison identity: {key}')
            if key in withheld:
                raise ValueError(f'Repeated withheld comparison: {key}')
            withheld[key] = {'comparison_key': key, 'final_class': group['final_class'],
                             'reason': group['reason'], 'reason_code': REASON_CODES[group['reason']]}
    retained = source['retain_checked']
    if len(set(retained)) != len(retained) or set(retained) & set(withheld):
        raise ValueError('Retained decisions overlap or repeat')
    return withheld, set(retained)


def comparison_key(modality, record):
    return '|'.join(str(x) for x in (modality, record['gse_id'], record['comparison_name']))


def withhold(record, evidence):
    """Preserve estimates as historical evidence and remove family membership."""
    record['previous_final_class'] = record['final_class']
    record['final_class'] = evidence['final_class']
    record['reasons'] = sorted(set(record['reasons']) | {evidence['reason_code']})
    record['pending_reason'] = [evidence['reason']]
    record['comparison_reclassification'] = dict(evidence)
    action = 'WITHHOLD_CLOCK_OUTPUT_ARMS' if evidence['reason_code'] == REASON else 'WITHHOLD_GEO_ARM_MISMATCH'
    record['actions'].append({'action': action, 'reason': evidence['reason']})
    recomputed = record.get('recomputed', {})
    record['withheld_family'] = {k: recomputed.pop(k) for k in ('fdr_family', 'fdr_family_size') if k in recomputed}


def apply(records, decisions, modality):
    """Apply explicit decisions to every inferential clock and its duplicate rows."""
    kept = {r['row_ordinal']: decisions[comparison_key(modality, r)] for r in records
            if r.get('final_class_before_verification_rules', r['final_class']) in INFERENTIAL
            and comparison_key(modality, r) in decisions}
    for record in records:
        if record['final_class'] == 'DUPLICATE' and record.get('duplicate_of') in kept:
            evidence = {**kept[record['duplicate_of']], 'representative': record['duplicate_of']}
            kept[record['row_ordinal']] = evidence
            record['representative_class'] = evidence['final_class']
    for record in records:
        if record['row_ordinal'] in kept:
            withhold(record, kept[record['row_ordinal']])
    return kept
