"""Source-bound statistical reruns with immutable, resumable comparison receipts.

Numerical verification is conditional on the declared design. Neither agreement
between formulas nor a source quote establishes scientific validity. The runner
keeps the source review, original result, proposed calculation and release
decision separate. It never executes model-generated code or overwrites inputs.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from numbers import Real
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import scipy
from scipy import linalg, optimize, stats
from scipy.special import stdtr

from .agent_workflow_2026_09_10 import sha, write_new
from .source_packet import read_source_csv
from src.pipeline_runtime_2026_09_10 import StageStore, PipelineIntegrityError

SCHEMA = 'clockbase-source-bound-statistical-rerun-v1'
WELCH = {'welch', 'welch t test', 'welch s t test', 'welch two sample t test',
         'welch s two sample t test', 'welch independent samples t test',
         'welch independent t test', 'two sided welch t test'}
STUDENT = {'student t test', 'student s t test', 'student two sample t test'}
PAIRED = {'paired t test', 'two sided paired t test'}
OLS = {'ols adjusted contrast', 'ols hc3 adjusted contrast'}
POOLED = 'case-minus-control divided by pooled within-group sample SD'
PAIRED_DZ = 'mean within-unit case-minus-control divided by SD of within-unit differences'


def record_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _critical_t(df, confidence_level):
    """Invert the lower t tail with bounded, checked convergence.

    Some supported SciPy builds lose precision in t.ppf. Retain the existing
    two-sided probability arithmetic, but solve against stdtr independently
    of that inverse. Failure remains a numerical stop, never a method change.
    """
    if any(isinstance(v, (bool, np.bool_)) or not isinstance(v, Real)
           or not math.isfinite(v) for v in (df, confidence_level)):
        raise ValueError('Critical-t inputs must be finite real numbers')
    if df <= 0 or not 0 < confidence_level < 1:
        raise ValueError('Critical-t requires positive df and confidence in (0, 1)')
    probability = float((1 + confidence_level) / 2)
    if not 0.5 < probability < 1:
        raise ValueError('Critical-t probability is not representable inside (0.5, 1)')
    tail = 1 - probability

    def residual(value):
        cdf = float(stdtr(df, -value))
        if not math.isfinite(cdf) or not 0 <= cdf <= 0.5:
            raise ValueError('Critical-t CDF is nonfinite or outside its lower tail')
        return cdf - tail

    upper = 1.
    for _ in range(1024):
        if residual(upper) <= 0:
            break
        upper *= 2
        if not math.isfinite(upper):
            raise ValueError('Critical-t has no finite bracket')
    else:
        raise ValueError('Critical-t bracketing did not converge')
    try:
        critical, convergence = optimize.brentq(
            residual, 0., upper, xtol=math.nextafter(0., 1.),
            rtol=8 * sys.float_info.epsilon, maxiter=256, full_output=True)
    except (RuntimeError, ValueError) as exc:
        raise ValueError('Critical-t inversion failed: ' + str(exc)) from exc
    if not convergence.converged or not math.isfinite(critical) or critical <= 0:
        raise ValueError('Critical-t inversion did not meet its accuracy check')
    error, limit = abs(residual(critical)), 64 * sys.float_info.epsilon * tail
    if error > limit:
        # Brent's relative stopping interval can span several floats. The
        # 8-epsilon bound spans at most 16 usual spacings; 32 also covers a
        # power-of-two transition. Keep the same tail-accuracy requirement.
        below = above = critical
        for _ in range(32):
            below, above = math.nextafter(below, 0.), math.nextafter(above, math.inf)
            for candidate in (below, above):
                if math.isfinite(candidate) and candidate > 0:
                    candidate_error = abs(residual(candidate))
                    if candidate_error < error:
                        critical, error = candidate, candidate_error
            if error <= limit:
                break
    if error > limit:
        raise ValueError('Critical-t inversion did not meet its accuracy check')
    return float(critical)


def recompute_contrast(frame, case_ids, control_ids, clock, test, *, design=None,
                       confidence_level=0.95):
    """Recompute a prespecified two-sided test; no automatic method selection.

    Paired tests require an explicit one-to-one biological-unit map. OLS accepts
    a numeric design fixed before calculation; the source runner derives its
    columns from exact metadata fields. Missing observations are never dropped.
    """
    normalized = re.sub(r'[^a-z0-9]+', ' ', test.lower()).strip()
    if normalized not in WELCH | STUDENT | PAIRED | OLS:
        raise NotImplementedError('Declared statistical test has no verified adapter: ' + test)
    if normalized in PAIRED | OLS and design is None:
        raise NotImplementedError('Declared statistical test requires an explicit source-bound design: ' + test)
    if not 0 < confidence_level < 1:
        raise ValueError('Confidence level must lie strictly between zero and one')
    case_ids, control_ids = list(case_ids), list(control_ids)
    members = case_ids + control_ids
    if not frame.index.is_unique or len(members) != len(set(members)) or set(members) - set(frame.index):
        raise ValueError('Declared comparison duplicates or invents sample membership')
    if min(len(case_ids), len(control_ids)) < 2:
        raise ValueError('Declared comparison has fewer than two samples per group')
    if clock not in frame or list(frame.columns).count(clock) != 1:
        raise ValueError('Declared clock is absent or ambiguous in predictions')
    a = frame.loc[case_ids, clock].to_numpy(float)
    b = frame.loc[control_ids, clock].to_numpy(float)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Declared comparison contains nonfinite predictions')
    difference = float(a.mean() - b.mean())
    va, vb = float(a.var(ddof=1)), float(b.var(ddof=1))
    pooled = ((len(a)-1)*va + (len(b)-1)*vb) / (len(a)+len(b)-2)
    result = {'clock': clock, 'case_ids': case_ids, 'control_ids': control_ids,
              'n_case': len(a), 'n_control': len(b), 'mean_difference': difference,
              'case_mean': float(a.mean()), 'control_mean': float(b.mean()),
              'alternative': 'two-sided', 'confidence_level': confidence_level,
              'missingness_policy': 'stop; no complete-case deletion or imputation',
              'assumptions_status': 'CANNOT VERIFY',
              'assumptions_note': 'Calculations do not establish independence, normality, covariate adequacy or clock applicability.',
              'scientific_status': 'CANNOT VERIFY'}
    if normalized in PAIRED:
        pairs = design.get('pairs', [])
        if len(pairs) != len(a) or len(a) != len(b):
            raise ValueError('Explicit pairs must cover every sample exactly once')
        if (set(p.get('case_id') for p in pairs) != set(case_ids) or
                set(p.get('control_id') for p in pairs) != set(control_ids) or
                any(not isinstance(p.get('unit_id'), str) or not p['unit_id'].strip() for p in pairs) or
                len(set(p['unit_id'] for p in pairs)) != len(pairs)):
            raise ValueError('Explicit pairs lack unique units or exact arm membership')
        paired_a = frame.loc[[p['case_id'] for p in pairs], clock].to_numpy(float)
        paired_b = frame.loc[[p['control_id'] for p in pairs], clock].to_numpy(float)
        changes = paired_a - paired_b
        difference = float(changes.mean())
        change_sd = float(changes.std(ddof=1))
        if change_sd <= 0:
            raise ValueError('Degenerate within-pair variance; t inference unavailable')
        se2, df = change_sd**2 / len(pairs), len(pairs)-1
        library = stats.ttest_rel(paired_a, paired_b, nan_policy='raise', alternative='two-sided')
        p = float(library.pvalue)
        result.update(statistical_test='Paired t-test', pairs=pairs, n_biological_units=len(pairs),
                      mean_difference=difference, within_unit_differences=changes.tolist(),
                      effect_size=difference/change_sd, effect_size_definition=PAIRED_DZ,
                      pooled_standardized_difference=difference/math.sqrt(pooled) if pooled > 0 else None,
                      pooled_standardized_difference_definition=POOLED)
    elif normalized in OLS:
        covariates = design.get('covariates')
        if not isinstance(covariates, list) or not covariates:
            raise ValueError('Adjusted test requires explicit covariates; no automatic covariate search')
        names = ['intercept', 'case']
        columns = [np.ones(len(members)), np.array([1.]*len(a) + [0.]*len(b))]
        for covariate in covariates:
            name, values = covariate.get('name'), covariate.get('values', {})
            if not isinstance(name, str) or not name.strip() or name in names or set(values) != set(members):
                raise ValueError('Covariate names or sample identities are invalid')
            names.append(name)
            columns.append(np.array([float(values[s]) for s in members]))
        x, y = np.column_stack(columns), np.concatenate([a, b])
        if not np.isfinite(x).all():
            raise ValueError('Covariate contains nonfinite or missing values')
        n, k = x.shape
        if n <= k or np.linalg.matrix_rank(x) != k or np.linalg.cond(x) > 1e10:
            raise ValueError('Confounded, rank-deficient, ill-conditioned or saturated covariate design')
        q, r = linalg.qr(x, mode='economic')
        coefficients = linalg.solve_triangular(r, q.T @ y)
        independent_coefficients = np.linalg.lstsq(x, y, rcond=None)[0]
        if not np.allclose(coefficients, independent_coefficients, rtol=1e-9, atol=1e-11):
            raise ValueError('Independent least-squares solutions disagree')
        residuals = y - x @ coefficients
        inv_r = linalg.solve_triangular(r, np.eye(k))
        bread = inv_r @ inv_r.T
        df = n-k
        leverage = np.square(q).sum(axis=1)
        if normalized == 'ols hc3 adjusted contrast':
            if np.any(1-leverage <= 1e-10):
                raise ValueError('Unit leverage prevents HC3 covariance estimation')
            covariance = bread @ (x.T @ (x * ((residuals/(1-leverage))**2)[:, None])) @ bread
            variance = 'HC3 with residual-df t reference'
        else:
            covariance = float(residuals @ residuals)/df * bread
            variance = 'classical homoskedastic OLS'
        difference, se2 = float(coefficients[1]), float(covariance[1, 1])
        if se2 <= 0:
            raise ValueError('Degenerate adjusted contrast variance')
        p = float(2 * stats.t.sf(abs(difference/math.sqrt(se2)), df))
        result.update(statistical_test='OLS adjusted contrast', variance_estimator=variance,
                      adjusted_mean_difference=difference, effect_size=None,
                      effect_size_definition='unstandardized case coefficient; standardized effect not defined',
                      design_columns=names, design_matrix=x.tolist(), coefficients=coefficients.tolist(),
                      residuals=residuals.tolist(), leverage=leverage.tolist(),
                      covariates=covariates, residual_degrees_of_freedom=df)
    else:
        equal_var = normalized in STUDENT
        se2 = pooled*(1/len(a)+1/len(b)) if equal_var else va/len(a)+vb/len(b)
        if se2 <= 0 or pooled <= 0:
            raise ValueError('Degenerate within-group variance; standardized contrast unavailable')
        df = len(a)+len(b)-2 if equal_var else se2**2/((va/len(a))**2/(len(a)-1)+(vb/len(b))**2/(len(b)-1))
        p = float(stats.ttest_ind(a, b, equal_var=equal_var, nan_policy='raise').pvalue)
        result.update(statistical_test='Student t-test' if equal_var else 'Welch t-test',
                      effect_size=difference/math.sqrt(pooled), effect_size_definition=POOLED,
                      case_variance=va, control_variance=vb)
    se = math.sqrt(se2)
    statistic = difference/se
    independent_p = float(2 * stdtr(df, -abs(statistic)))
    critical = _critical_t(df, confidence_level)
    ci = [difference-critical*se, difference+critical*se]
    if not all(math.isfinite(v) for v in (difference, se, df, p, independent_p, *ci)):
        raise ValueError('Nonfinite recomputed contrast')
    if not math.isclose(p, independent_p, rel_tol=1e-10, abs_tol=1e-14):
        raise ValueError('Independent t-reference formula disagrees with library result')
    result.update(p_value=p, independent_p_value=independent_p, standard_error=se,
                  statistic=statistic, degrees_of_freedom=float(df), confidence_interval=ci,
                  confidence_interval_estimand=('adjusted case coefficient' if normalized in OLS else 'case-minus-control mean'),
                  numerical_verification='PASS', verification_scope='conditional numerical consistency; independent source verification is separate')
    return result


def source_record(record):
    path = Path(record['path'])
    if not path.is_file():
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    if sha(path) != record['sha256']:
        raise PipelineIntegrityError('Frozen source changed: ' + str(path))
    return path


def select_groups(metadata, proposal):
    """Retain every source row matching the prespecified exact filters."""
    def arm(filters):
        if not isinstance(filters, dict) or not filters or any(not isinstance(v, str) for v in filters.values()):
            raise ValueError('Nonempty exact source filters are required')
        selected = [row for row in metadata if all(row.get(k) == v for k, v in filters.items())]
        return [row.get('GSM', row.get('geo_accession')) for row in selected]
    case, control = arm(proposal['case_filter']), arm(proposal['control_filter'])
    if not case or not control or None in case+control or len(set(case+control)) != len(case+control):
        raise ValueError('Source groups are empty, overlap or lack unique sample identity')
    return case, control


def _label_pairs(metadata, field, rule, case, control):
    """Parse a frozen label rule; biological interpretation is assessed separately.

    Every record in the declared literal scope is parsed. Other strata/arms are
    recorded explicitly, and the reconstructed groups must equal the supplied
    complete selectors. No raw field is changed and no incomplete pair is dropped.
    """
    keys = {'schema', 'delimiter', 'token_count', 'fixed_tokens', 'stratum_tokens',
            'unit_tokens', 'arm_token', 'case_token', 'control_token', 'scope_filters'}
    if not isinstance(rule, dict) or set(rule) != keys or rule['schema'] != 'clockbase-delimited-pair-v1':
        raise ValueError('Unknown or incomplete source-label pair rule')
    delimiter, count = rule['delimiter'], rule['token_count']
    if not isinstance(delimiter, str) or not 1 <= len(delimiter) <= 8 or type(count) is not int or not 2 <= count <= 32:
        raise ValueError('Pair-label delimiter or token count is invalid')

    def positions(values, name):
        if not isinstance(values, dict):
            raise ValueError(name + ' must declare exact token positions and values')
        parsed = {}
        for key, value in values.items():
            if (not isinstance(key, str) or not re.fullmatch(r'0|[1-9][0-9]*', key)
                    or int(key) >= count or not isinstance(value, str) or not value):
                raise ValueError(name + ' contains an invalid position or literal value')
            parsed[int(key)] = value
        return parsed

    fixed = positions(rule['fixed_tokens'], 'fixed_tokens')
    strata = positions(rule['stratum_tokens'], 'stratum_tokens')
    units, arm = rule['unit_tokens'], rule['arm_token']
    if (not strata or not isinstance(units, list) or not units
            or any(type(i) is not int or not 0 <= i < count for i in units)
            or len(units) != len(set(units)) or type(arm) is not int or not 0 <= arm < count):
        raise ValueError('Pair-label stratum, unit or arm positions are invalid')
    if (set(strata) & (set(units) | set(fixed) | {arm}) or arm in set(units) | set(fixed)
            or set(fixed) | set(strata) | set(units) | {arm} != set(range(count))):
        raise ValueError('Pair-label token roles overlap or leave an undeclared token')
    case_token, control_token = rule['case_token'], rule['control_token']
    if (not isinstance(case_token, str) or not case_token or not isinstance(control_token, str)
            or not control_token or case_token == control_token):
        raise ValueError('Distinct exact case and control labels are required')
    filters = rule['scope_filters']
    if not isinstance(filters, dict) or any(not isinstance(k, str) or not k or not isinstance(v, str) or not v
                                           for k, v in filters.items()):
        raise ValueError('Pair-label scope filters must be exact nonempty source fields and values')
    if not isinstance(field, str) or not field:
        raise ValueError('Source-label pairing requires a literal pair_field')
    if (not case or not control or len(set(case + control)) != len(case + control)
            or any(not isinstance(sid, str) or not sid for sid in case + control)):
        raise ValueError('Pair-label groups lack unique, disjoint sample identities')

    maps = {'case': {}, 'control': {}}
    witnesses, excluded, seen = [], [], set()
    scoped = 0
    for row in metadata:
        sid = row.get('GSM', row.get('geo_accession'))
        if not isinstance(sid, str) or not sid or sid in seen:
            raise ValueError('Source-label metadata sample identities are absent or duplicated')
        seen.add(sid)
        if any(row.get(k) != v for k, v in filters.items()):
            excluded.append({'sample_id': sid, 'reason': 'OUTSIDE_LITERAL_SCOPE',
                             'scope_field_values': {k: row.get(k) for k in filters}})
            continue
        scoped += 1
        raw = row.get(field)
        if not isinstance(raw, str) or not raw:
            raise ValueError('Missing pair-label source field for ' + sid)
        tokens = raw.split(delimiter)
        if len(tokens) != count or any(not token for token in tokens) or any(tokens[i] != v for i, v in fixed.items()):
            raise ValueError('Pair label does not match the complete declared token structure: ' + sid)
        witness = {'sample_id': sid, 'field': field, 'value': raw, 'tokens': tokens}
        if any(tokens[i] != v for i, v in strata.items()):
            excluded.append({**witness, 'reason': 'OTHER_DECLARED_STRATUM'})
            continue
        if tokens[arm] not in {case_token, control_token}:
            excluded.append({**witness, 'reason': 'OTHER_ARM_LABEL'})
            continue
        side = 'case' if tokens[arm] == case_token else 'control'
        unit = json.dumps({'stratum': [tokens[i] for i in sorted(strata)],
                           'unit': [tokens[i] for i in units]}, sort_keys=True, separators=(',', ':'))
        if unit in maps[side]:
            raise ValueError('Repeated source-derived unit within an arm; no aggregation is allowed')
        maps[side][unit] = sid
        witnesses.append({**witness, 'arm': side, 'derived_unit_id': unit})
    if (set(maps['case'].values()) != set(case) or set(maps['control'].values()) != set(control)):
        raise ValueError('Derived complete source groups differ from the selected arms')
    if set(maps['case']) != set(maps['control']):
        raise ValueError('Unmatched derived source units; no complete-pair deletion')
    pairs = [{'unit_id': unit, 'case_id': maps['case'][unit], 'control_id': maps['control'][unit]}
             for unit in sorted(maps['case'])]
    return {'pairs': pairs, 'pair_identity': {
        'schema': 'clockbase-derived-pair-map-v1', 'source_field': field,
        'rule': json.loads(json.dumps(rule, allow_nan=False)),
        'raw_field_witnesses': sorted(witnesses, key=lambda row: row['sample_id']),
        'excluded_samples': sorted(excluded, key=lambda row: row['sample_id']),
        'counts': {'source_records': len(seen), 'scoped_records': scoped,
                   'paired_label_units': len(pairs), 'included_samples': len(witnesses),
                   'excluded_samples': len(excluded)},
        'scientific_status': 'CANNOT VERIFY',
        'scope': 'Exact source-label derivation and membership only; physical sample derivation, independence, statistical assumptions and clock applicability are separate judgments.'}}


def bind_design(metadata, proposal, case, control):
    """Derive declared units and covariates from metadata, without outcomes."""
    by_id = {row.get('GSM', row.get('geo_accession')): row for row in metadata}
    if len(by_id) != len(metadata):
        raise ValueError('Source metadata sample identities are not unique')
    design = {}
    if proposal.get('pair_identity') is not None:
        design.update(_label_pairs(metadata, proposal.get('pair_field'), proposal['pair_identity'], case, control))
    elif proposal.get('pair_field'):
        field = proposal['pair_field']
        def units(ids):
            result = {}
            for sid in ids:
                value = by_id[sid].get(field, '').strip()
                if not value or value in result:
                    raise ValueError('Missing or repeated biological unit within an arm')
                result[value] = sid
            return result
        a, b = units(case), units(control)
        if set(a) != set(b):
            raise ValueError('Unmatched source biological units; no complete-pair deletion')
        design['pairs'] = [{'unit_id': unit, 'case_id': a[unit], 'control_id': b[unit]} for unit in sorted(a)]
    if proposal.get('covariates'):
        columns = []
        for spec in proposal['covariates']:
            if not spec.get('justification') or not spec.get('role') in {'baseline_confounder', 'design_stratum'}:
                raise ValueError('Each covariate requires a prespecified baseline/design role and justification')
            values = {sid: by_id[sid].get(spec['field'], '') for sid in case+control}
            if any(not v.strip() for v in values.values()):
                raise ValueError('Missing source covariate; no imputation or row deletion')
            if spec['type'] == 'numeric':
                prefix = spec.get('prefix', '')
                if any(not v.startswith(prefix) for v in values.values()):
                    raise ValueError('Source covariate prefix does not match')
                columns.append({'name': spec['field'], 'values': {sid: float(v[len(prefix):]) for sid, v in values.items()}})
            elif spec['type'] == 'categorical':
                levels = spec['levels']
                if len(levels) < 2 or len(levels) != len(set(levels)) or set(values.values()) != set(levels):
                    raise ValueError('Categorical levels must exactly match the prespecified source levels')
                for level in levels[1:]:
                    columns.append({'name': spec['field']+'='+level, 'values': {sid: float(v == level) for sid, v in values.items()}})
            else:
                raise ValueError('Unsupported source covariate type')
        design['covariates'] = columns
    return design


def clock_frame(path, clock, *, sample_id_column=None):
    """Read an exact unique clock column; never silently mangle duplicate names."""
    try:
        with Path(path).open(newline='') as handle:
            rows = list(csv.reader(handle))
    except (OSError, csv.Error, UnicodeError) as error:
        raise ValueError('Cannot read prediction CSV: ' + str(error)) from error
    if not rows or rows[0].count(clock) != 1:
        raise ValueError('Requested clock column is absent or ambiguous: ' + clock)
    header = rows[0]
    id_names = [name for name in ('sample_id', 'gsm_id', 'GSM', 'geo_accession') if name in header]
    id_columns = [i for i, name in enumerate(header) if name in id_names]
    if sample_id_column is not None:
        index = sample_id_column.get('index')
        if (type(index) is not int or not 0 <= index < len(header) or
                header[index] != sample_id_column.get('header') or
                sample_id_column.get('format') != 'GSM accession'):
            raise ValueError('Explicit source ID-column binding is invalid')
        id_columns = list(dict.fromkeys([index]+id_columns))
    if not id_columns:
        raise ValueError('Prediction CSV lacks an explicit sample identity column')
    clock_index = header.index(clock)
    data = []
    for row in rows[1:]:
        if len(row) != len(header) or len({row[i] for i in id_columns}) != 1:
            raise ValueError('Ragged prediction CSV or conflicting sample identity columns')
        if sample_id_column is not None and not re.fullmatch(r'GSM[0-9]+', row[id_columns[0]]):
            raise ValueError('Explicit source sample ID is not a literal GSM accession')
        data.append((row[id_columns[0]], row[clock_index]))
    frame = pd.DataFrame(data, columns=['sample_id', clock]).set_index('sample_id')
    if not frame.index.is_unique or any(not sid for sid in frame.index):
        raise ValueError('Prediction IDs are blank or duplicated')
    return frame


def execute_proposal(plan, proposal, reference):
    """Detect, propose and calculate; independent verification remains separate."""
    sources = plan['sources'][proposal['study_id']]
    paths = {key: source_record(value) for key, value in sources.items()}
    metadata, _ = read_source_csv(paths['metadata'])
    case, control = select_groups(metadata, proposal)
    reviewed = reference
    for key in proposal['reference_target']:
        reviewed = reviewed[key]
    if set(reviewed['case_ids']) != set(case) or set(reviewed['control_ids']) != set(control):
        raise PipelineIntegrityError('Source-defined groups disagree with the frozen independent reference')
    if reviewed.get('clock', proposal['clock']) != proposal['clock']:
        raise PipelineIntegrityError('Requested clock disagrees with the frozen reference')
    method_reference = reference
    for key in proposal['reference_method']:
        method_reference = method_reference[key]
    methods = {
        'paired_t_by_donor': ('inferential_conditional', PAIRED),
        'welch_independent_conditional_on_independence': ('inferential_conditional', WELCH),
        'student_independent_conditional_on_equal_variance': ('inferential_conditional', STUDENT),
        'ols_adjusted_conditional_on_design': ('inferential_conditional', {'ols adjusted contrast'}),
        'ols_hc3_adjusted_conditional_on_design': ('inferential_conditional', {'ols hc3 adjusted contrast'}),
        'descriptive_technical_replicates': ('descriptive_technical', set()),
        'missing_input_stop': ('missing_input_stop', set()),
    }
    if method_reference not in methods:
        raise ValueError('Independent source reference does not support an executable method')
    mode, tests = methods[method_reference]
    test = re.sub(r'[^a-z0-9]+', ' ', proposal.get('statistical_test', '').lower()).strip()
    if proposal['mode'] != mode or (tests and test not in tests):
        raise PipelineIntegrityError('Proposed mode/test contradicts the independent source reference')
    originals, _ = read_source_csv(paths['original'])
    ordinal = proposal['original_row_ordinal']
    if type(ordinal) is not int or ordinal < 1 or ordinal > len(originals):
        raise ValueError('Invalid original logical row ordinal')
    original = originals[ordinal-1]
    if original.get('gse_id') != proposal['study_id'] or original.get('clock') != proposal['clock']:
        raise PipelineIntegrityError('Original row identity differs from the proposal')
    if record_hash(original) != proposal['original_record_sha256']:
        raise PipelineIntegrityError('Original comparison/condition record differs from the frozen proposal')
    identity = sources['original']['sha256'] + ':' + str(ordinal)
    record = {'proposal_id': proposal['proposal_id'], 'study_id': proposal['study_id'],
              'clock': proposal['clock'], 'original_row_id': identity, 'original': original,
              'original_row_ordinal': ordinal, 'original_record_number': ordinal+1,
              'sources': sources, 'case_ids': case, 'control_ids': control,
              'proposal': proposal, 'design_review_status': 'REVIEWED',
              'scientific_status': 'CANNOT VERIFY', 'scientific_limitations': proposal['scientific_limitations'],
              'historical_membership_status': 'CANNOT VERIFY',
              'measurement_applicability': 'CANNOT VERIFY', 'independent_verification': 'NOT_YET_VERIFIED',
              'release_status': 'WITHHELD', 'accepted_scientific_repair': False,
              'adjusted_p_value': None, 'multiplicity_status': 'CANNOT VERIFY',
              'multiplicity_reason': plan['multiplicity']['reason'],
              'clock_execution': False, 'original_overwritten': False,
              'detections': [], 'candidate': None}
    for field in ('p_value', 'effect_size'):
        try:
            value = float(original[field])
            valid = math.isfinite(value) and (field != 'p_value' or 0 <= value <= 1)
        except (KeyError, ValueError, TypeError):
            valid = False
        if not valid:
            record['detections'].append({'status': 'FAIL', 'field': field,
                'scope': 'literal original record', 'reason': 'Not a finite number in the permitted domain',
                'proposed_action': 'withhold invalid field; do not invent a replacement'})
    try:
        frame = clock_frame(paths['predictions'], proposal['clock'],
                            sample_id_column=sources['predictions'].get('sample_id_column'))
        missing = sorted(set(case+control) - set(frame.index))
        if missing != sorted(reviewed.get('missing_prediction_ids', [])):
            raise PipelineIntegrityError('Missing prediction IDs disagree with the frozen independent reference')
        if missing:
            record.update(calculation_status='CANNOT VERIFY', execution_state='UNRESOLVED',
                          missing_prediction_ids=missing,
                          reason='Required source-defined samples lack predictions; no samples omitted')
            return record
        if proposal['mode'] == 'descriptive_technical':
            a, b = frame.loc[case, proposal['clock']].to_numpy(float), frame.loc[control, proposal['clock']].to_numpy(float)
            if not np.isfinite(np.concatenate([a, b])).all():
                raise ValueError('Nonfinite descriptive input')
            record['candidate'] = {'case_mean': float(a.mean()), 'control_mean': float(b.mean()),
                'mean_difference': float(a.mean()-b.mean()), 'n_case_assays': len(a), 'n_control_assays': len(b),
                'p_value': None, 'effect_size': None, 'confidence_interval': None,
                'inference_status': 'WITHHELD', 'scope': proposal['estimand']}
            record['detections'].append({'status': 'FAIL', 'scope': 'donor-population inference from source-defined assays',
                'reason': proposal['technical_replication_reason'], 'proposed_action': 'descriptive_only'})
        elif proposal['mode'] == 'inferential_conditional':
            design = bind_design(metadata, proposal, case, control)
            if test in PAIRED:
                expected = [{'unit_id': p['donor'], 'case_id': p['case_gsm'], 'control_id': p['control_gsm']}
                            for p in reviewed.get('pair_map', [])]
                if sorted(design.get('pairs', []), key=record_hash) != sorted(expected, key=record_hash) or not expected:
                    raise PipelineIntegrityError('Source-derived pair map disagrees with the frozen independent reference')
            if test in OLS and (not reviewed.get('covariates') or proposal.get('covariates') != reviewed['covariates']):
                raise PipelineIntegrityError('Covariate design differs from the frozen independent reference')
            candidate = recompute_contrast(frame, case, control, proposal['clock'], proposal['statistical_test'],
                design=design, confidence_level=plan['confidence_level'])
            record['candidate'] = candidate
            # The saved effect definition is a separately declared comparison;
            # paired dz must never be compared to a pooled-group-SD effect.
            fields = {'p_value': 'p_value'}
            if proposal.get('saved_effect_definition') == POOLED:
                fields['effect_size'] = ('pooled_standardized_difference' if 'pairs' in candidate else 'effect_size')
            record['conditional_original_agreement'] = {}
            for field, candidate_field in fields.items():
                try:
                    prior, current = float(original[field]), float(candidate[candidate_field])
                    agrees = math.isfinite(prior) and math.isclose(prior, current,
                        rel_tol=plan['relative_tolerance'], abs_tol=plan['absolute_tolerance'])
                except (ValueError, TypeError, KeyError):
                    agrees = False
                record['conditional_original_agreement'][field] = 'PASS' if agrees else 'FAIL'
            record['conditional_agreement_scope'] = 'source-defined companion data only; historical input lineage is unresolved'
        elif proposal['mode'] == 'missing_input_stop':
            raise ValueError('Frozen missing-input stop now has complete data; obtain a new source review before inference')
        else:
            raise ValueError('Unsupported proposal mode')
        record.update(calculation_status='PASS', execution_state='RERUN')
    except PipelineIntegrityError:
        raise
    except (ValueError, TypeError, KeyError, NotImplementedError) as error:
        record.update(calculation_status='CANNOT VERIFY', execution_state='UNRESOLVED',
                      candidate=None, reason=type(error).__name__+': '+str(error))
    return record


def run(plan_path, output):
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = json.loads(plan_path.read_text())
    if plan.get('schema') != SCHEMA or plan.get('evaluation_role') != 'exposed_development':
        raise ValueError('This runner requires an explicit exposed-development plan; held-out execution uses the frozen evaluation runner')
    if plan.get('multiplicity', {}).get('action') != 'withhold_adjusted_inference':
        raise ValueError('Complete testing-family validation is required in the separate family gate')
    proposals = plan['proposals']
    ids = [p['proposal_id'] for p in proposals]
    if len(ids) != len(set(ids)) or any(not re.fullmatch(r'[A-Za-z0-9_-]+', p) for p in ids):
        raise ValueError('Proposal IDs must be unique safe stage names')
    records = [plan['independent_design_reference']]
    records += [source for study in plan['sources'].values() for source in study.values()]
    for record in records:
        source_record(record)
    reference = json.loads(source_record(plan['independent_design_reference']).read_text())
    code_paths = [Path(__file__), Path(__file__).with_name('source_packet.py'),
                  Path(__file__).with_name('agent_workflow_2026_09_10.py'),
                  Path(__file__).resolve().parents[1]/'src/pipeline_runtime_2026_09_10.py']
    config = {'schema': SCHEMA, 'plan_sha256': sha(plan_path), 'plan': plan,
              'implementation': {str(p): sha(p) for p in code_paths},
              'versions': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__}}
    store = StageStore(output, config)
    results = []
    for proposal in proposals:
        def calculate(directory, proposal=proposal):
            result = execute_proposal(plan, proposal, reference)
            for source in plan['sources'][proposal['study_id']].values():
                source_record(source)
            write_new(directory/'result.json', result)
            return result
        results.append(store.run(proposal['proposal_id'], calculate))
    for record in records:
        source_record(record)
    def summarize(directory):
        write_new(directory/'calculations.json', results)
        table = []
        for result in results:
            c = result['candidate'] or {}
            table.append({'proposal_id': result['proposal_id'], 'study_id': result['study_id'],
                'clock': result['clock'], 'execution_state': result['execution_state'],
                'calculation_status': result['calculation_status'], 'original_p_value': result['original'].get('p_value'),
                'rerun_p_value': c.get('p_value'), 'mean_difference': c.get('mean_difference'),
                'confidence_interval': json.dumps(c.get('confidence_interval')),
                'scientific_status': result['scientific_status'], 'release_status': result['release_status']})
        if table:
            with (directory/'validation-table.csv').open('x', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(table[0])); writer.writeheader(); writer.writerows(table)
        summary = {'schema': SCHEMA, 'recorded_at_utc': datetime.now(timezone.utc).isoformat(),
            'studies': len(set(r['study_id'] for r in results)), 'planned_targets': len(results),
            'original_rows': len(set(r['original_row_id'] for r in results)),
            'rerun_targets': sum(r['execution_state'] == 'RERUN' for r in results),
            'unresolved_targets': sum(r['execution_state'] == 'UNRESOLVED' for r in results),
            'inferential_conditional_targets': sum(bool(r['candidate'] and r['candidate'].get('p_value') is not None) for r in results),
            'descriptive_only_targets': sum(bool(r['candidate'] and r['candidate'].get('inference_status') == 'WITHHELD') for r in results),
            'accepted_scientific_repairs': 0, 'independently_verified_targets': 0,
            'independent_verification': 'NOT_YET_VERIFIED', 'scientific_status': 'CANNOT VERIFY',
            'originals_preserved': True, 'clock_executions': 0, 'model_calls': 0, 'api_cost_usd': 0,
            'scope': 'source-bound exposed-development statistical reruns; source clock outputs are reused'}
        write_new(directory/'summary.json', summary)
        return summary
    return store.run('summary', summarize)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.plan, args.output), indent=2))
