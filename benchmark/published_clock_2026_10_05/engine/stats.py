"""Published estimands followed by the verifier's Welch test and pooled Cohen's d."""
import re

import numpy as np
import pandas as pd

from audit_system.row_repairs_2026_09_28 import ESTIMAND, recompute
from .io import clean


def threshold(record):
    value = record.get('threshold')
    if isinstance(value, (int, float)) and not isinstance(value, bool) and 0 < value < 1:
        return float(value)
    if isinstance(value, str):
        found = re.fullmatch(r'\s*(?:[pP]\s*)?(?:<|<=|≤|=)?\s*(0?\.\d+|\d+(?:\.\d+)?[eE]-\d+)\s*', value)
        if found and 0 < float(found[1]) < 1:
            return float(found[1])
    return None


def direction(value):
    if value is None or not np.isfinite(value):
        return None
    return 'higher' if value > 0 else 'lower' if value < 0 else 'none'


def evaluate(samples, groups, measure, alpha=None):
    if groups.get('status') not in ('AGREED', 'MANUAL'):
        return {'status': 'UNRESOLVED_GROUPS', 'covered': False}
    if measure not in ('age', 'acceleration_residual', 'acceleration_difference', 'rate'):
        return {'status': 'UNSUPPORTED_MEASURE', 'covered': False}
    values = pd.to_numeric(samples['clock'], errors='coerce').replace([np.inf, -np.inf], np.nan)
    fit = None
    if measure.startswith('acceleration_'):
        age = pd.to_numeric(samples['age'], errors='coerce').replace([np.inf, -np.inf], np.nan)
        if measure == 'acceleration_difference':
            values = values - age
        else:
            usable = values.notna() & age.notna()
            if usable.sum() < 2 or age[usable].nunique() < 2:
                return {'status': 'INSUFFICIENT_RECORDED_AGE', 'covered': False}
            design = np.column_stack([np.ones(usable.sum()), age[usable]])
            intercept, slope = np.linalg.lstsq(design, values[usable], rcond=None)[0]
            values = (values - (intercept + slope * age)).where(usable)
            fit = {'method': 'OLS with intercept', 'population': 'all series samples with finite clock and recorded age',
                   'n': int(usable.sum()), 'intercept': float(intercept), 'slope': float(slope),
                   'gsms': sorted(map(str, samples.index[usable]))}
    case, control = groups['case_gsms'], groups['control_gsms']
    result = clean(recompute(values, case, control))
    covered = result['n_case'] > 0 and result['n_control'] > 0
    inferential = result['n_case'] >= 3 and result['n_control'] >= 3
    p = result.get('p')
    return {'status': 'INFERENTIAL' if inferential else 'DESCRIPTIVE' if covered else 'EMPTY_ARM',
            'covered': covered, 'measure': measure, 'estimand': ESTIMAND, **result,
            'direction': direction(result.get('mean_diff')),
            'significant': bool(p < alpha) if p is not None and alpha is not None else None,
            'threshold': alpha, 'fit': fit,
            'excluded_case_gsms': sorted(g for g in case if pd.isna(values.get(g, np.nan))),
            'excluded_control_gsms': sorted(g for g in control if pd.isna(values.get(g, np.nan)))}
