"""Design-respecting estimand of a two-arm comparison, from the sample records of its series.

The declared estimand of the verifier is a Welch two-sample t-test on the arms. It assumes that every sample is an
independent unit and that nothing but the labelled condition separates the arms. Three designs break that, and the
sample records often say so:

  repeated individual   one individual contributes several samples to an arm
  paired arms           at least half of the individuals have samples in both arms
  recorded ages         the arms of a human comparison carry ages, and the clock tracks age

For such a row the statistic is recomputed on the same arms:

  repeated individual   individual means replace the samples (Welch test on individual means)
  paired arms           paired t-test on the individual means of the individuals in both arms
  recorded ages         least squares of the clock on arm and age; the arm coefficient is the estimate

Individuals come from a GEO sample field named as an individual identifier (individual, subject, patient, donor,
participant, animal, mouse, person, volunteer, with an optional id or number suffix), or from the two-reader
individuals audit of a comparison when the series has no such field. A field qualifies when it is recorded for at
least 80% of the series' samples, has at least two values, at least one value shared by two samples, and at least 80%
of its values read as identifiers. Ages are the harmonized age; for a human series without one, a GEO field named as
the age at sampling whose values parse as years. Age enters only for human series, clocks in age units, ages recorded
for at least 80% of each arm, and comparisons that do not contrast age or time (no age-group word in the names and,
when the contrast classification is given, neither model answering AGE or DEVELOPMENT_OR_TIME). When age separates
the arms completely the age term cannot be fitted: the row keeps its statistic and carries AGE_SEPARATES_ARMS.

Clock values with an absolute value above 1,000 are overflow artefacts of the clock software (values of the order of
1e39 occur); `overflow` lists them so that the caller can drop the samples.

Validation: benchmark/published_clock_2026_10_05/groupbench (PROTOCOL.md addendum, items 10 to 13).
"""
import math
import re
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

OVERFLOW = 1000.0
MIN_UNITS = 3
COVERAGE = 0.8
ID_SHARE = 0.8
PAIRED_SHARE = 0.5
HUMAN = {'dnam', 'hs_rna'}
RATE_CLOCKS = {'DunedinPACE'}
AGE_CONTRASTS = {'AGE', 'DEVELOPMENT_OR_TIME'}
SUBJECT_WORDS = 'individual|subject|patient|donor|participant|animal|mouse|person|volunteer'
SUBJECT_NAME = re.compile(rf'^(source |original |biological )?({SUBJECT_WORDS})s?'
                          r'([ _.-]?(id|ids|identifier|number|no|nr|num|code|name|#))?$')
MISSING = {'', 'na', 'n/a', 'nan', 'none', 'null', 'unknown', 'not available', 'not applicable', 'not provided', '-', '--'}
AGE_WORDS = re.compile(
    r'\b(aged|aging|ageing|old|older|young|younger|elderly|adult|adults|juvenile|neonat\w*|newborn\w*|infant\w*|'
    r'child\w*|pediatric|paediatric|adolescen\w*|fetal|foetal|fetus|embryo\w*|postnatal|prenatal|centenarian\w*)\b', re.I)
AGE_FIELD = re.compile(
    r'^((donor|patient|patients|subject|participant|individual|visit|baseline|sample)[ _.]?)?age'
    r'([ _.]?(of donor|at (sampling|sample collection|collection|draw|blood draw|visit|biopsy|surgery|enrol\w*|baseline\w*)|'
    r'sampling|draw/autopsy))?[ _.]?(in[ _.])?[(\[]?(years?|yrs?|y)?[)\]]?( at enrolment)?$')
YEARS = re.compile(r'\s*(\d+(\.\d+)?)\s*(y|yr|yrs|year|years|years old|yo)?\s*')
ESTIMANDS = {
    'UNIT_MEANS': 'Welch two-sample t-test on individual means; Cohen d with pooled SD; case minus control',
    'PAIRED': ('Paired t-test on the individual means of individuals in both arms; mean within-individual difference; '
               'Cohen d with the pooled SD of the two conditions; case minus control'),
    'AGE_ADJUSTED': ('Least squares of the clock on arm and age; arm coefficient; Cohen d with the residual SD; '
                     'case minus control'),
    'AGE_ADJUSTED_UNIT_MEANS': ('Least squares of individual means of the clock on arm and age; arm coefficient; '
                                'Cohen d with the residual SD; case minus control'),
}


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


def subject_table(raw):
    """(field, {gsm: individual}) from the finest qualifying GEO individual field, or (None, {})."""
    candidates = []
    if raw is not None:
        for column in raw.columns:
            if SUBJECT_NAME.match(str(column).strip().lower()):
                values = {str(g): clean(v) for g, v in raw[column].items()}
                if qualifies(values, len(raw)):
                    candidates.append((len(set(v for v in values.values() if v)), str(column), values))
    if not candidates:
        return None, {}
    _, field, values = max(candidates, key=lambda c: (c[0], c[1]))
    return field, {g: v for g, v in values.items() if v is not None}


def age_table(raw, processed, human):
    """(field, {gsm: age}): the harmonized age, else for a human series a GEO age-at-sampling field in years."""
    if processed is not None and 'harmonized.age' in processed:
        numeric = pd.to_numeric(processed['harmonized.age'], errors='coerce').dropna()
        if len(numeric):
            return 'harmonized.age', {str(g): float(v) for g, v in numeric.items()}
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
                return str(column), parsed
    return None, {}


def read_frame(path):
    """A sample table indexed by GSM (duplicated GSMs dropped), or None when the file or the column is absent."""
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


def load_series(geo, folder, gse, human):
    """The design inputs of one series: individual field and map, age field and map, and the files read."""
    processed_path = geo / folder / 'processed_metadata' / f'{gse}.csv'
    raw_path = geo / folder / 'processed_metadata_preharmonize' / f'{gse}.csv'
    processed, raw = read_frame(processed_path), read_frame(raw_path)
    subject_field, subjects = subject_table(raw)
    age_field, ages = age_table(raw, processed, human)
    return {'subject_field': subject_field, 'subjects': subjects, 'subject_source': 'GEO sample field' if subject_field else None,
            'age_field': age_field, 'ages': ages,
            'files': [str(p) for p in (raw_path, processed_path) if p.exists()]}


def overflow(values, gsms):
    """GSMs whose clock value is an overflow artefact."""
    out = []
    for gsm in gsms:
        if gsm in values.index:
            value = values[gsm]
            if isinstance(value, (int, float, np.floating)) and np.isfinite(value) and abs(value) > OVERFLOW:
                out.append(gsm)
    return out


def age_declared(names, contrast_types=()):
    """Whether the comparison contrasts age or time: by its names, or by a contrast classification."""
    return bool(AGE_WORDS.search(' '.join(str(n) for n in names if isinstance(n, str)))) \
        or bool(set(contrast_types) & AGE_CONTRASTS)


def subject_findings(case, control, subjects):
    """Repeated-individual and pairing findings for arms whose individuals are recorded for at least 80% of samples."""
    by_arm = {'case': [subjects[g] for g in case if g in subjects], 'control': [subjects[g] for g in control if g in subjects]}
    if not case or not control or len(by_arm['case']) + len(by_arm['control']) < COVERAGE * (len(case) + len(control)):
        return []
    everyone = set(by_arm['case']) | set(by_arm['control'])
    both = set(by_arm['case']) & set(by_arm['control'])
    findings = []
    if len(everyone) >= 2 and any(n > 1 for ids in by_arm.values() for n in Counter(ids).values()):
        findings.append('REPEATED_INDIVIDUAL')
    if len(both) >= 2 and len(both) >= PAIRED_SHARE * len(everyone):
        findings.append('PAIRED_ARMS')
    return findings


def unit_means(values, arm, subjects, ages):
    """{unit: (mean clock value, mean recorded age or None, samples)}; a sample without an individual is its own unit."""
    grouped = {}
    for gsm in arm:
        if gsm in values.index and pd.notna(values[gsm]):
            grouped.setdefault(subjects.get(gsm, f'sample:{gsm}'), []).append(gsm)
    out = {}
    for unit, members in grouped.items():
        recorded = [ages[g] for g in members if g in ages]
        out[unit] = (float(np.mean([values[g] for g in members])), float(np.mean(recorded)) if recorded else None, len(members))
    return out


CONSTANT_OUTPUT_REASON = 'outputs constant within both arms'
CONSTANT_OUTPUT_CODE = 'OUTPUTS_CONSTANT_WITHIN_BOTH_ARMS'
CONSTANT_OUTPUT_RTOL = 1e-12


def constant_output(a, b):
    """Whether both arms are constant, including pooled SD at the round-off floor.

    Empty arms do not establish constancy. A singleton has zero within-arm
    variation; exact single-value arms are also covered when no SD is estimable.
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if not len(a) or not len(b):
        return False
    if np.all(a == a[0]) and np.all(b == b[0]):
        return True
    df = len(a) + len(b) - 2
    if df <= 0:
        return False
    ss = sum(float(np.sum((v - v.mean()) ** 2)) for v in (a, b))
    pooled = math.sqrt(ss / df)
    return pooled <= CONSTANT_OUTPUT_RTOL * max(1.0, abs(float(a.mean())), abs(float(b.mean())))


def constant_summary(a, b):
    """Descriptive mean difference, with minimum replication taking precedence.

    Small arms retain the ordinary minimum-unit reason in the caller. Only
    arms meeting the declared minimum receive the constant-output reason.
    """
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    # An exact constant's mean is its value, without summation round-off.
    mean = lambda v: float(v[0]) if np.all(v == v[0]) else float(v.mean())
    ma, mb = mean(a), mean(b)
    out = {'n_case': len(a), 'n_control': len(b), 'mean_case': ma, 'mean_control': mb,
           'mean_diff': ma - mb}
    if len(a) < MIN_UNITS or len(b) < MIN_UNITS:
        return out
    out['undefined_reason'] = CONSTANT_OUTPUT_REASON
    return out


def two_sample(a, b):
    """The declared estimand on two arrays of independent units (the formulas of row_repairs.recompute)."""
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if constant_output(a, b):
        return constant_summary(a, b)
    out = {'n_case': int(len(a)), 'n_control': int(len(b)),
           'mean_case': float(a.mean()) if len(a) else None, 'mean_control': float(b.mean()) if len(b) else None}
    if len(a) and len(b):
        out['mean_diff'] = float(a.mean() - b.mean())
    if len(a) >= 2 and len(b) >= 2:
        pooled = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
        out['cohen_d'] = float(out['mean_diff'] / pooled) if pooled > 0 else None
    if len(a) >= MIN_UNITS and len(b) >= MIN_UNITS:
        test = stats.ttest_ind(a, b, equal_var=False)
        out['p'] = float(test.pvalue) if np.isfinite(test.pvalue) else None
        out['t'] = float(test.statistic) if np.isfinite(test.statistic) else None
    return out


def paired(case_units, control_units):
    """Paired statistics over the individuals in both arms (samples without an individual cannot be paired)."""
    shared = sorted(u for u in set(case_units) & set(control_units) if not u.startswith('sample:'))
    x = np.array([case_units[u][0] for u in shared], dtype=float)
    y = np.array([control_units[u][0] for u in shared], dtype=float)
    if constant_output(x, y):
        return constant_summary(x, y), len(set(case_units) | set(control_units)) - len(shared)
    out = {'n_case': int(len(shared)), 'n_control': int(len(shared)),
           'mean_case': float(x.mean()) if len(x) else None, 'mean_control': float(y.mean()) if len(y) else None}
    if len(shared):
        out['mean_diff'] = float((x - y).mean())
    if len(shared) >= 2:
        pooled = math.sqrt((x.var(ddof=1) + y.var(ddof=1)) / 2)
        out['cohen_d'] = float(out['mean_diff'] / pooled) if pooled > 0 else None
    if len(shared) >= MIN_UNITS and not np.allclose(x - y, (x - y)[0]):
        test = stats.ttest_rel(x, y)
        out['p'] = float(test.pvalue) if np.isfinite(test.pvalue) else None
        out['t'] = float(test.statistic) if np.isfinite(test.statistic) else None
    return out, len(set(case_units) | set(control_units)) - len(shared)


def age_adjusted(case_units, control_units):
    """('OK', statistics, slope) of clock ~ arm + age over units with a recorded age, or (reason, None, None)."""
    case = [(v, a) for v, a, _ in case_units.values() if a is not None]
    control = [(v, a) for v, a, _ in control_units.values() if a is not None]
    if len(case) < MIN_UNITS or len(control) < MIN_UNITS:
        return 'TOO_FEW_UNITS_WITH_AGE', None, None
    case_ages, control_ages = [a for _, a in case], [a for _, a in control]
    case_values, control_values = [v for v, _ in case], [v for v, _ in control]
    if constant_output(case_values, control_values):
        return 'OK', constant_summary(case_values, control_values), None
    if max(case_ages) < min(control_ages) or max(control_ages) < min(case_ages):
        return 'AGE_SEPARATES_ARMS', None, None
    y = np.array([v for v, _ in case + control], dtype=float)
    x = np.array([[1.0, 1.0, a] for _, a in case] + [[1.0, 0.0, a] for _, a in control], dtype=float)
    if np.linalg.matrix_rank(x) < 3:
        return 'AGE_NOT_ESTIMABLE', None, None
    # Remove the response baseline before fitting. Otherwise a tiny arm effect
    # on a nearly constant response loses precision against the intercept.
    centered_y = y - y[0]
    beta, *_ = np.linalg.lstsq(x, centered_y, rcond=None)
    df = len(y) - 3
    sigma2 = float(np.sum((centered_y - x @ beta) ** 2) / df)
    se = math.sqrt(sigma2 * np.linalg.inv(x.T @ x)[1, 1])
    if se == 0 or sigma2 <= 0:
        return 'AGE_NOT_ESTIMABLE', None, None
    t = float(beta[1] / se)
    out = {'n_case': len(case), 'n_control': len(control),
           'mean_case': float(np.mean([v for v, _ in case])), 'mean_control': float(np.mean([v for v, _ in control])),
           'mean_diff': float(beta[1]), 'cohen_d': float(beta[1] / math.sqrt(sigma2)),
           'p': float(2 * stats.t.sf(abs(t), df)), 't': t}
    return 'OK', out, float(beta[2])


def estimate(values, case, control, series, modality, clock, names=(), contrast_types=()):
    """The design-respecting statistic of one row, or None when the declared Welch estimand stands unchanged.

    Returns {'recomputed', 'estimand', 'design', 'soft'}; 'recomputed' is None when only a soft flag applies
    (age separates the arms of a comparison on independent samples).
    """
    subjects, ages = series.get('subjects') or {}, series.get('ages') or {}
    present = lambda arm: [g for g in arm if g in values.index and pd.notna(values[g])]
    case, control = present(case), present(control)
    findings = subject_findings(case, control, subjects) if subjects else []
    if series.get('force_individuals') and subjects:
        if not all(g in subjects for g in case + control):
            raise ValueError('Single-cell estimand requires an individual for every sample')
        findings = sorted(set(findings) | {'RECORDED_SINGLE_CELL_INDIVIDUALS'})
    recorded = lambda arm: sum(g in ages for g in arm)
    usable = (modality in HUMAN and clock not in RATE_CLOCKS and bool(ages) and bool(case) and bool(control)
              and recorded(case) >= COVERAGE * len(case) and recorded(control) >= COVERAGE * len(control)
              and not age_declared(names, contrast_types))
    if not findings and not usable:
        return None
    collapse = subjects if findings else {}
    case_units, control_units = unit_means(values, case, collapse, ages), unit_means(values, control, collapse, ages)
    design = {'findings': findings, 'samples_case': len(case), 'samples_control': len(control),
              'units_case': len(case_units), 'units_control': len(control_units)}
    if findings:
        design.update(subject_field=series.get('subject_field'), subject_source=series.get('subject_source'))
    soft = []
    if 'PAIRED_ARMS' in findings:
        result, dropped = paired(case_units, control_units)
        design.update(method='PAIRED', pairs=result['n_case'], unpaired_individuals_dropped=dropped)
        return {'recomputed': result, 'estimand': ESTIMANDS['PAIRED'], 'design': design, 'soft': soft}
    if usable:
        status, result, slope = age_adjusted(case_units, control_units)
        if status == 'OK':
            method = 'AGE_ADJUSTED_UNIT_MEANS' if findings else 'AGE_ADJUSTED'
            design.update(method=method, age_field=series.get('age_field'), age_slope=slope,
                          units_without_age=len(case_units) + len(control_units) - result['n_case'] - result['n_control'])
            return {'recomputed': result, 'estimand': ESTIMANDS[method], 'design': design, 'soft': soft}
        design['age_term'] = status
        if status == 'AGE_SEPARATES_ARMS':
            soft.append('AGE_SEPARATES_ARMS')
            design['age_field'] = series.get('age_field')
    if findings:
        result = two_sample([v for v, _, _ in case_units.values()], [v for v, _, _ in control_units.values()])
        design['method'] = 'UNIT_MEANS'
        return {'recomputed': result, 'estimand': ESTIMANDS['UNIT_MEANS'], 'design': design, 'soft': soft}
    if soft:
        return {'recomputed': None, 'estimand': None, 'design': design, 'soft': soft}
    return None
