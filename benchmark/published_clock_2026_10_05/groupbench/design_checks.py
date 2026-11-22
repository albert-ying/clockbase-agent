"""Prototype of three design checks the released verifier does not make (not part of the code release).

For each released row (VERIFIED or PROVISIONAL) of the final verifier run, on its post-repair arms:

- REPEATED_SUBJECT: one subject contributes more than one sample to an arm (at least two subjects in the comparison).
- PAIRED_AS_INDEPENDENT: at least half of the subjects have samples in both arms; the released statistic is an
  unpaired Welch test.
- AGE_IMBALANCE (human series): the arms differ in recorded age (median gap above 5 years, or a standardized mean
  difference above 0.5), the clock is not a pace clock, and age or time is not the declared contrast (neither of the
  verifier's two contrast classifications says AGE or DEVELOPMENT_OR_TIME, and the names carry no age-group word).

Subjects come from a GEO sample field whose name is a subject identifier (individual, subject, patient, donor,
participant, animal, mouse, person, volunteer, with an optional id or number suffix). A field qualifies when it is
recorded for at least 80% of the series' samples, has at least two values, at least one value shared by two samples,
and at least 80% of its values look like identifiers. Ages are the harmonized age; for human series without one, a GEO
field named as the age at sampling whose values parse as years.

Variant v2 (--variant v2) also used the harmonized subject field when no GEO field qualified and applied the age check
to mouse series (relative median gap above 25%). Two blinded readers confirmed 0 of 2 sampled flags from the harmonized
subject field and 9 of 15 from mouse ages (wrong harmonized ages), so the default (v3) leaves both out.

Each flagged row is recomputed with a statistic that respects the design: subject means per arm; a paired t-test over
the subjects in both arms when the design is paired, otherwise a Welch test on subject means; under age imbalance an
ordinary least-squares model of the clock on arm and age. The released call (direction, two-sided P < 0.05) is
compared with the recomputed call.

  python -m benchmark.published_clock_2026_10_05.groupbench.design_checks --out DIR [--modality dnam hs_rna mm_rna]
      [--variant v3|v2]
"""
import argparse
import hashlib
import json
import math
import re
import statistics
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

COVERAGE = 0.8
ID_SHARE = 0.8
PAIRED_SHARE = 0.5
MIN_PAIRS = 3
AGE_YEARS, AGE_SMD, AGE_RELATIVE = 5.0, 0.5, 0.25
ALPHA = 0.05
RELEASED = ('VERIFIED', 'PROVISIONAL')
HUMAN = {'dnam', 'hs_rna'}
SUBJECT_WORDS = 'individual|subject|patient|donor|participant|animal|mouse|person|volunteer'
SUBJECT_NAME = re.compile(rf'^(source |original |biological )?({SUBJECT_WORDS})s?'
                          r'([ _.-]?(id|ids|identifier|number|no|nr|num|code|name|#))?$')
MISSING = {'', 'na', 'n/a', 'nan', 'none', 'null', 'unknown', 'not available', 'not applicable', 'not provided', '-', '--'}
AGE_WORDS = re.compile(
    r'\b(aged|aging|ageing|old|older|young|younger|elderly|adult|adults|juvenile|neonat\w*|newborn\w*|infant\w*|'
    r'child\w*|pediatric|paediatric|adolescen\w*|fetal|foetal|fetus|embryo\w*|postnatal|prenatal|centenarian\w*)\b', re.I)
AGE_CONTRASTS = {'AGE', 'DEVELOPMENT_OR_TIME'}
AGE_FIELD = re.compile(
    r'^((donor|patient|patients|subject|participant|individual|visit|baseline|sample)[ _.]?)?age'
    r'([ _.]?(of donor|at (sampling|sample collection|collection|draw|blood draw|visit|biopsy|surgery|enrol\w*|baseline\w*)|'
    r'sampling|draw/autopsy))?[ _.]?(in[ _.])?[(\[]?(years?|yrs?|y)?[)\]]?( at enrolment)?$')
YEARS = re.compile(r'\s*(\d+(\.\d+)?)\s*(y|yr|yrs|year|years|years old|yo)?\s*')
PACE = re.compile(r'pace of aging', re.I)


NOT_SAMPLE_AGE = ('maternal', 'diagnosis', 'gestational')


def formula_has_age(formula):
    """Whether a model formula carries the sample's age as a term on its right-hand side.

    The response is a clock whose name ends in "Age", and terms such as agent, passage or stage contain the letters,
    so the test reads the right-hand side term by term.
    """
    text = str(formula or '').lower()
    if '~' not in text:
        return False
    for term in re.findall(r'[a-z_][a-z0-9_.]*', text.split('~', 1)[1]):
        if any(word in term for word in NOT_SAMPLE_AGE):
            continue
        if 'age' in re.split(r'[_.]', term) or term in ('agebin', 'agegroup'):
            return True
    return False


def clean(value):
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    return None if text.lower() in MISSING else text


def identifier(value):
    """A value that reads as an identifier: it has a digit, or is a one- or two-letter code after an optional word."""
    return bool(re.search(r'\d', value)) or bool(re.fullmatch(r'(\w+[ _-])?[A-Za-z]{1,2}', value))


def qualifies(values, n_samples):
    recorded = [v for v in values.values() if v is not None]
    distinct = set(recorded)
    return (n_samples > 0 and len(recorded) >= COVERAGE * n_samples and len(distinct) >= 2
            and len(distinct) < len(recorded)
            and sum(identifier(v) for v in distinct) >= ID_SHARE * len(distinct))


def subject_table(raw, processed, harmonized=False):
    """{'field', 'source', 'subjects': {gsm: id}} from the finest qualifying subject field, or None."""
    candidates = []
    if raw is not None:
        for column in raw.columns:
            if SUBJECT_NAME.match(str(column).strip().lower()):
                values = {str(g): clean(v) for g, v in raw[column].items()}
                if qualifies(values, len(raw)):
                    candidates.append((len(set(v for v in values.values() if v)), 'geo', str(column), values))
    if harmonized and not candidates and processed is not None and 'harmonized.subject' in processed:
        values = {str(g): clean(v) for g, v in processed['harmonized.subject'].items()}
        if qualifies(values, len(processed)):
            candidates.append((len(set(v for v in values.values() if v)), 'harmonized', 'harmonized.subject', values))
    if not candidates:
        return None
    _, source, field, values = max(candidates, key=lambda c: (c[0], c[2]))
    return {'field': field, 'source': source, 'subjects': {g: v for g, v in values.items() if v is not None}}


def age_declared(row):
    """Why age or time is the declared contrast of a row: 'classification', 'name', both, or None."""
    answers = ((row.get('ranking_evidence') or {}).get('classification') or {}).get('answers') or []
    reasons = []
    if any(a.get('contrast_type') in AGE_CONTRASTS for a in answers):
        reasons.append('classification')
    if AGE_WORDS.search(' '.join(t for t in (row.get('comparison_name'), row.get('condition_name')) if t)):
        reasons.append('name')
    return '+'.join(reasons) or None


def age_table(raw, processed, human):
    """({gsm: age}, source field): the harmonized age, else for human series a GEO age-at-sampling field in years."""
    if processed is not None and 'harmonized.age' in processed:
        numeric = pd.to_numeric(processed['harmonized.age'], errors='coerce').dropna()
        if len(numeric):
            return {str(g): float(v) for g, v in numeric.items()}, 'harmonized.age'
    if human and raw is not None:
        for column in raw.columns:
            if not AGE_FIELD.match(str(column).strip().lower()):
                continue
            parsed = {}
            for gsm, value in raw[column].items():
                match = YEARS.fullmatch(str(value).lower()) if clean(value) else None
                if match:
                    parsed[str(gsm)] = float(match.group(1))
            if len(parsed) >= COVERAGE * len(raw) and all(0 <= v <= 120 for v in parsed.values()):
                return parsed, str(column)
    return {}, None


def subject_findings(case, control, subjects):
    """Counts of subjects per arm and the repeated-subject and pairing findings."""
    by_arm = {'case': [subjects[g] for g in case if g in subjects], 'control': [subjects[g] for g in control if g in subjects]}
    known = len(by_arm['case']) + len(by_arm['control'])
    if known < COVERAGE * (len(case) + len(control)):
        return {'assessable': False}
    everyone = set(by_arm['case']) | set(by_arm['control'])
    both = set(by_arm['case']) & set(by_arm['control'])
    out = {'assessable': True, 'subjects': len(everyone), 'subjects_case': len(set(by_arm['case'])),
           'subjects_control': len(set(by_arm['control'])), 'subjects_both': len(both), 'findings': []}
    repeated = {arm: sorted(s for s, n in Counter(ids).items() if n > 1) for arm, ids in by_arm.items()}
    if len(everyone) >= 2 and (repeated['case'] or repeated['control']):
        out['findings'].append('REPEATED_SUBJECT')
        out['repeated'] = {arm: len(ids) for arm, ids in repeated.items()}
    if len(both) >= 2 and len(both) >= PAIRED_SHARE * len(everyone):
        out['findings'].append('PAIRED_AS_INDEPENDENT')
    return out


def smd(a, b):
    if len(a) < 2 or len(b) < 2:
        return 0.0
    pooled = math.sqrt((statistics.variance(a) + statistics.variance(b)) / 2)
    return abs(statistics.mean(a) - statistics.mean(b)) / pooled if pooled else 0.0


def age_finding(case, control, ages, human):
    """The age gap between arms and which criterion it meets; None when age is not recorded for 80% of each arm."""
    a = [ages[g] for g in case if g in ages]
    b = [ages[g] for g in control if g in ages]
    if not a or not b or len(a) < COVERAGE * len(case) or len(b) < COVERAGE * len(control):
        return None
    gap = statistics.median(a) - statistics.median(b)
    out = {'median_case': statistics.median(a), 'median_control': statistics.median(b), 'median_gap': gap,
           'smd': round(smd(a, b), 3), 'separated': max(a) < min(b) or max(b) < min(a), 'criterion': None}
    if human:
        if abs(gap) > AGE_YEARS:
            out['criterion'] = 'median_gap'
        elif out['smd'] > AGE_SMD:
            out['criterion'] = 'smd_only'
    else:
        scale = max(abs(statistics.median(a)), abs(statistics.median(b)), 1e-9)
        if abs(gap) / scale > AGE_RELATIVE:
            out['criterion'] = 'median_gap'
    return out


def call(estimate, p):
    if estimate is None or p is None or not math.isfinite(p):
        return None
    if p >= ALPHA or estimate == 0:
        return 'ns'
    return 'higher' if estimate > 0 else 'lower'


def welch(a, b):
    if len(a) < 2 or len(b) < 2:
        return None
    result = stats.ttest_ind(a, b, equal_var=False)
    return {'estimate': float(np.mean(a) - np.mean(b)), 'p': float(result.pvalue), 'n': [len(a), len(b)]}


def units(arm, values, subjects, ages=None):
    """Subject means of the clock (and age) within one arm; a sample without a subject is its own unit."""
    grouped = {}
    for gsm in arm:
        if gsm in values:
            grouped.setdefault(subjects.get(gsm, f'sample:{gsm}'), []).append(gsm)
    out = {}
    for subject, members in grouped.items():
        recorded = [ages[g] for g in members if ages and g in ages]
        out[subject] = (float(np.mean([values[g] for g in members])), float(np.mean(recorded)) if recorded else None)
    return out


def adjusted(case, control):
    """Arm coefficient of clock ~ arm + age; case and control are lists of (clock, age)."""
    rows = [(v, 1.0, a) for v, a in case if a is not None] + [(v, 0.0, a) for v, a in control if a is not None]
    case_ages = [a for _, a in case if a is not None]
    control_ages = [a for _, a in control if a is not None]
    if len(case_ages) < 2 or len(control_ages) < 2 or len(rows) < 5:
        return {'status': 'TOO_FEW_SAMPLES'}
    if max(case_ages) < min(control_ages) or max(control_ages) < min(case_ages):
        return {'status': 'AGE_SEPARATES_ARMS'}
    y = np.array([r[0] for r in rows])
    x = np.array([[1.0, r[1], r[2]] for r in rows])
    if np.linalg.matrix_rank(x) < 3:
        return {'status': 'AGE_NOT_ESTIMABLE'}
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    df = len(rows) - 3
    sigma2 = float(np.sum((y - x @ beta) ** 2) / df)
    se = math.sqrt(sigma2 * np.linalg.inv(x.T @ x)[1, 1])
    if se == 0:
        return {'status': 'AGE_NOT_ESTIMABLE'}
    return {'status': 'OK', 'estimate': float(beta[1]), 'p': float(2 * stats.t.sf(abs(beta[1] / se), df)),
            'n': [len(case_ages), len(control_ages)], 'age_slope': float(beta[2])}


def recompute(case, control, values, subjects, ages, findings):
    """The design-respecting statistic for a flagged row."""
    a, b = units(case, values, subjects, ages), units(control, values, subjects, ages)
    if 'PAIRED_AS_INDEPENDENT' in findings:
        shared = sorted(set(a) & set(b))
        shared = [s for s in shared if not s.startswith('sample:')]
        if len(shared) < MIN_PAIRS:
            return {'method': 'paired t-test on subject means', 'status': 'TOO_FEW_PAIRS', 'n': [len(shared)]}
        x, y = [a[s][0] for s in shared], [b[s][0] for s in shared]
        if np.allclose(np.array(x) - np.array(y), (np.array(x) - np.array(y))[0]):
            return {'method': 'paired t-test on subject means', 'status': 'NO_VARIANCE', 'n': [len(shared)]}
        result = stats.ttest_rel(x, y)
        return {'method': 'paired t-test on subject means', 'status': 'OK', 'estimate': float(np.mean(x) - np.mean(y)),
                'p': float(result.pvalue), 'n': [len(shared)]}
    if 'AGE_IMBALANCE' in findings:
        result = adjusted(list(a.values()), list(b.values()))
        return {'method': 'least squares on arm and age' + (', subject means' if 'REPEATED_SUBJECT' in findings else ''),
                **result}
    result = welch([v for v, _ in a.values()], [v for v, _ in b.values()])
    if result is None:
        return {'method': 'Welch t-test on subject means', 'status': 'ONE_SUBJECT_ARM', 'n': [len(a), len(b)]}
    return {'method': 'Welch t-test on subject means', 'status': 'OK', **result}


def change(released, recomputed):
    if released is None or recomputed is None:
        return None
    if released == recomputed:
        return 'UNCHANGED'
    if released == 'ns':
        return 'SIGNIFICANCE_GAINED'
    if recomputed == 'ns':
        return 'SIGNIFICANCE_LOST'
    return 'DIRECTION_REVERSED'


def unit_key(gse, case, control):
    arms = sorted(['|'.join(sorted(case)), '|'.join(sorted(control))])
    return f"{gse}:{hashlib.sha1('#'.join(arms).encode()).hexdigest()[:12]}"


SUBJECT_FINDINGS = {'REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT'}
NOT_ESTIMABLE = {'AGE_NOT_ESTIMABLE', 'TOO_FEW_SAMPLES'}


def candidate_statistic(findings, age_usable, case, control, values, subjects, ages, released):
    """The candidate release estimand of PROTOCOL.md addendum item 10 for one row."""
    kept = {'method': 'Welch t-test (released)', 'status': 'OK', 'estimate': released['estimate'], 'p': released['p']}
    if set(findings) & SUBJECT_FINDINGS:
        use = [f for f in findings if f in SUBJECT_FINDINGS] + (['AGE_IMBALANCE'] if age_usable else [])
        result = recompute(case, control, values, subjects, ages, use)
        if result.get('status') in NOT_ESTIMABLE and age_usable:  # the age term cannot be fitted on subject means
            result = recompute(case, control, values, subjects, ages, [f for f in use if f != 'AGE_IMBALANCE'])
        return result
    if age_usable:
        result = recompute(case, control, values, {}, ages, ['AGE_IMBALANCE'])
        return kept if result.get('status') in NOT_ESTIMABLE else result
    return kept


def check_row(row, modality, table, ages, clock_values, mouse_age=False):
    """Findings for one released row; clock_values is called only when the row is flagged or its ages are usable."""
    case, control = [str(g) for g in row['case_gsms']], [str(g) for g in row['control_gsms']]
    released = row.get('recomputed') or {}
    test = (row.get('membership') or {}).get('test')
    formula = str(((row.get('replay') or {}).get('test') or {}).get('formula') or '')
    target = (row.get('evidence') or {}).get('clock_target') or ''
    record = {'modality': modality, 'gse': row['gse'], 'row_ordinal': row['row_ordinal'], 'clock': row['clock'],
              'comparison_name': row.get('comparison_name'), 'condition_name': row.get('condition_name'),
              'final_class': row['final_class'], 'unit': unit_key(row['gse'], case, control),
              'n_case': len(case), 'n_control': len(control), 'agent_test': test,
              'agent_paired': test in ('ttest_rel', 'sm.MixedLM'), 'agent_age_adjusted': formula_has_age(formula),
              'released': {'estimate': released.get('mean_diff'), 'p': released.get('p'),
                           'fdr': released.get('fdr_family'), 'call': call(released.get('mean_diff'), released.get('p'))},
              'findings': []}
    if table:
        subject = subject_findings(case, control, table['subjects'])
        record['subject'] = {'field': table['field'], 'source': table['source'], **subject}
        record['findings'] += subject.get('findings', [])
    age = age_finding(case, control, ages, modality in HUMAN) if ages else None
    if age:
        record['age'] = age
        record['age_declared'] = age_declared(row)
        record['pace_clock'] = bool(PACE.search(target))
        record['age_checked'] = modality in HUMAN or mouse_age
        if record['age_checked'] and age['criterion'] and not record['age_declared'] and not record['pace_clock']:
            record['findings'].append('AGE_IMBALANCE')
    age_usable = bool(age) and modality in HUMAN and not record['age_declared'] and not record['pace_clock']
    record['candidate'] = {'method': 'Welch t-test (released)', 'status': 'OK', 'estimate': record['released']['estimate'],
                           'p': record['released']['p'], 'call': record['released']['call']}
    if not record['findings'] and not age_usable:
        return record
    values = clock_values(row['clock'])
    if not values:
        if record['findings']:
            record['recomputed'] = {'status': 'NO_PREDICTIONS'}
        return record
    candidate = candidate_statistic(record['findings'], age_usable, case, control, values,
                                    table['subjects'] if table else {}, ages or {}, record['released'])
    candidate['call'] = call(candidate.get('estimate'), candidate.get('p')) if candidate.get('status') == 'OK' else None
    record['candidate'] = candidate
    if not record['findings']:
        return record
    check = welch([values[g] for g in case if g in values], [values[g] for g in control if g in values])
    record['released_reproduced'] = bool(check and released.get('p') is not None
                                         and math.isclose(check['p'], released['p'], rel_tol=1e-6, abs_tol=1e-12))
    recomputed = recompute(case, control, values, table['subjects'] if table else {}, ages or {}, record['findings'])
    recomputed['call'] = call(recomputed.get('estimate'), recomputed.get('p')) if recomputed.get('status') == 'OK' else None
    record['recomputed'] = recomputed
    record['change'] = change(record['released']['call'], recomputed['call'])
    return record


def read_frame(path):
    if not path.exists():
        return None
    try:
        frame = pd.read_csv(path, dtype=str)
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot read sample metadata {path}: {error}') from error
    if 'GSM' not in frame:
        return None
    frame = frame.set_index('GSM')
    frame.index = frame.index.astype(str)
    return frame[~frame.index.duplicated(keep=False)]


def series_inputs(modality, gse, harmonized=False):
    """The subject table, the recorded ages and the age field of one series."""
    from audit_system import paths
    from audit_system import membership_recovery_2026_09_27 as membership
    folder = paths.GEO / membership.MODALITIES[modality]['folder']
    processed = read_frame(folder / 'processed_metadata' / f'{gse}.csv')
    raw = read_frame(folder / 'processed_metadata_preharmonize' / f'{gse}.csv')
    ages, age_field = age_table(raw, processed, modality in HUMAN)
    return subject_table(raw, processed, harmonized), ages, age_field


def run(out, modalities, variant='v3'):
    from audit_system import paths
    from ..engine.predictions import PredictionReader
    reader = PredictionReader()
    out.mkdir(parents=True, exist_ok=True)
    from ..engine.predictions import final_export
    root = final_export() / 'studies'
    with (out / 'rows.jsonl').open('w') as rows_file, (out / 'series.jsonl').open('w') as series_file:
        for modality in modalities:
            for path in sorted((root / modality).glob('*.json')):
                gse = path.stem
                try:
                    rows = json.loads(path.read_text())['rows']
                except (OSError, ValueError) as error:
                    raise ValueError(f'Cannot read study export {path}: {error}') from error
                released = [r for r in rows if r.get('final_class') in RELEASED and r.get('case_gsms') and r.get('control_gsms')]
                if not released:
                    continue
                table, ages, age_field = series_inputs(modality, gse, harmonized=variant == 'v2')
                cache = {}

                def clock_values(clock):
                    if clock not in cache:
                        loaded = reader.load(modality, gse, clock, rows)
                        cache[clock] = ({str(k): float(v) for k, v in loaded['samples']['clock'].dropna().items()}
                                        if loaded.get('status') == 'OK' else {})
                    return cache[clock]

                for row in released:
                    record = check_row({**row, 'gse': gse}, modality, table, ages, clock_values, mouse_age=variant == 'v2')
                    rows_file.write(json.dumps(record) + '\n')
                series_file.write(json.dumps({'modality': modality, 'gse': gse, 'released_rows': len(released),
                                              'subject_field': table['field'] if table else None,
                                              'subject_source': table['source'] if table else None,
                                              'subjects': len(set(table['subjects'].values())) if table else None,
                                              'samples_with_age': len(ages), 'age_field': age_field}) + '\n')


def tally(records):
    """Counts of rows, distinct arm pairs and series for a set of records."""
    return {'rows': len(records), 'units': len({r['unit'] for r in records}), 'series': len({r['gse'] for r in records})}


def summarize(records):
    summary = {}
    for modality in sorted({r['modality'] for r in records}):
        rows = [r for r in records if r['modality'] == modality]
        flagged = [r for r in rows if r['findings']]
        block = {'released': tally(rows),
                 'subjects_recorded': tally([r for r in rows if r.get('subject', {}).get('assessable')]),
                 'age_recorded': tally([r for r in rows if r.get('age')]),
                 'flagged': tally(flagged), 'findings': {}, 'age_exempt': {}}
        for finding in ('REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT', 'AGE_IMBALANCE'):
            members = [r for r in flagged if finding in r['findings']]
            entry = tally(members)
            entry['only_this_finding'] = sum(r['findings'] == [finding] for r in members)
            entry['change'] = dict(Counter(r.get('change') or (r.get('recomputed') or {}).get('status') or 'NOT_COMPUTED'
                                           for r in members))
            entry['released_significant'] = sum(r['released']['call'] in ('higher', 'lower') for r in members)
            entry['released_fdr_significant'] = sum((r['released']['fdr'] or 1) < ALPHA for r in members)
            if finding == 'PAIRED_AS_INDEPENDENT':
                entry['agent_test_was_paired'] = sum(r['agent_paired'] for r in members)
            if finding == 'AGE_IMBALANCE':
                entry['agent_model_had_age'] = sum(r['agent_age_adjusted'] for r in members)
                entry['criterion'] = dict(Counter(r['age']['criterion'] for r in members))
            block['findings'][finding] = entry
        imbalanced = [r for r in rows if r.get('age') and r['age']['criterion'] and r.get('age_checked', True)]
        block['age_exempt'] = {'age_or_time_declared': tally([r for r in imbalanced if r['age_declared']]),
                               'declared_by': dict(Counter(r['age_declared'] for r in imbalanced if r['age_declared'])),
                               'pace_clock': tally([r for r in imbalanced if r['pace_clock'] and not r['age_declared']])}
        # Rows whose released Welch statistic replaced a paired or age-adjusted test of the agent.
        paired = [r for r in rows if r['agent_paired']]
        adjusted_rows = [r for r in rows if r['agent_age_adjusted']]
        block['agent_design'] = {
            'agent_test_paired': {**tally(paired), 'flagged_paired': sum('PAIRED_AS_INDEPENDENT' in r['findings'] for r in paired),
                                  'subjects_recorded': sum(bool(r.get('subject', {}).get('assessable')) for r in paired)},
            'agent_model_with_age': {**tally(adjusted_rows), 'flagged_age': sum('AGE_IMBALANCE' in r['findings'] for r in adjusted_rows),
                                     'age_recorded': sum(bool(r.get('age')) for r in adjusted_rows)}}
        block['change_any_finding'] = dict(Counter(r.get('change') or (r.get('recomputed') or {}).get('status') or 'NOT_COMPUTED'
                                                   for r in flagged))
        outcome = lambda r: r.get('change') or (r.get('recomputed') or {}).get('status') or 'NOT_COMPUTED'
        significant = [r for r in rows if r['released']['call'] in ('higher', 'lower')]
        fdr = [r for r in significant if (r['released']['fdr'] or 1) < ALPHA]
        block['released_significant'] = {'rows': len(significant), 'flagged': dict(Counter(outcome(r) for r in significant if r['findings']))}
        block['released_fdr_significant'] = {'rows': len(fdr), 'flagged': dict(Counter(outcome(r) for r in fdr if r['findings']))}
        block['released_not_reproduced'] = sum(r.get('released_reproduced') is False for r in flagged)
        summary[modality] = block
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--modality', nargs='+', default=['dnam', 'hs_rna', 'mm_rna'])
    parser.add_argument('--variant', choices=['v3', 'v2'], default='v3')
    parser.add_argument('--summarize-only', action='store_true')
    args = parser.parse_args()
    if not args.summarize_only:
        run(args.out, args.modality, args.variant)
    records = [json.loads(line) for line in (args.out / 'rows.jsonl').read_text().splitlines()]
    summary = summarize(records)
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
