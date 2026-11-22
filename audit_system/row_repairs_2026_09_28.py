"""Repair every atlas row the verifier withheld, and record the full repair path per row.

Rows with a recovered membership are recomputed on the saved per-sample predictions with one
declared estimand (Welch two-sample t-test; Cohen's d with pooled SD; case minus control), so every
repaired value is comparable across studies and categories. Repairs applied in order:

  NON_STANDARD_ASSAY          drop the flagged GSMs from both arms, then recompute
  INPUT_PROFILE_DETACHED      use the stream-repaired input when the study was re-derived from GEO; otherwise
                              drop the detached GSMs, and route the row to input repair if an arm empties
  LABEL_ORIENTATION_INVERTED  swap arms when both adjudicating models confirmed the inversion
  label match not confirmed   a row whose label match two models checked without agreeing keeps LABEL_UNCONFIRMED (soft)
  EFFECT_NOT_REPRODUCED       replaced by the recomputed effect under the declared estimand
  ARM_TOO_SMALL (<3 per arm)  DESCRIPTIVE: mean difference and d reported, no P, excluded from inference
  overflow clock values       drop GSMs whose clock value is an overflow artefact (absolute value above 1,000)
  design (--design-estimand)  when the sample records show a repeated individual, paired arms or recorded ages,
                              the statistic is recomputed with the design-respecting estimand of
                              design_estimand_2026_10_10 (individual means, paired t-test, age as a covariate); the
                              Welch result is kept beside it; fewer than 3 independent units per arm is DESCRIPTIVE
  FDR_BELOW_RAW_P / any FDR   BH recomputed within the declared family (all inferential repaired rows)
  repeated tests              rows with the same study, clock and arms as another row become DUPLICATE and
                              point to one representative; each test enters the BH family once

Rows still lacking membership or valid inputs are routed, with their reason, to design reconstruction
(PENDING_DESIGN) or input repair (PENDING_INPUT). Nothing is dropped. Original values stay beside the
repaired ones, and every row carries: original values, reasons, actions, inputs (paths/hashes),
recomputed values, final class and verdict.

  python -m audit_system.row_repairs_2026_09_28 --modality dnam --membership DIR --verdicts DIR \
      --replay DIR --out DIR [--adjudication DIR ...] [--design-estimand [--classification DIR] [--individuals DIR]]

Without --design-estimand the module behaves as before 2026-10-10 (every row under the Welch estimand, no overflow
drop); passes 1 and 2 of the verifier run that way, so the rows they select for the model stages do not change.
"""
from audit_system.adjudication_files import decision_files
import argparse
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system import design_estimand_2026_10_10 as design_estimand
from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import CORRECTIONS, MODALITIES, sha256
from audit_system.row_verdicts_2026_09_27 import NON_STANDARD_ASSAY, series_sample_text

SCHEMA = 'clockbase-row-repair-v1'
ESTIMAND = 'Welch two-sample t-test; Cohen d with pooled SD; case minus control'
MIN_ARM = 3
DESIGN_REASONS = {'MEMBERSHIP_UNRECOVERED', 'TEST_NOT_RECONSTRUCTIBLE', 'STORED_P_UNUSABLE', 'MEMBERSHIP_AMBIGUOUS',
                  'ONE_SAMPLE_TEST_NOT_CONTRAST', 'CONTINUOUS_TERM_NOT_GROUP_CONTRAST', 'ARMS_UNRESOLVED',
                  'PREDICTIONS_UNAVAILABLE'}
INPUT_REASONS = {'INPUT_EMPTY', 'RNA_INPUT_ARTIFACT', 'CLOCK_COVERAGE_LOW', 'COVERAGE_UNAVAILABLE', 'CLOCK_UNREGISTERED'}
SOFT = {'TISSUE_OUT_OF_DOMAIN', 'CLOCK_COVERAGE_MID', 'LABEL_UNCONFIRMED', 'EFFECT_SIGN_ONLY', 'CLOCK_AGE_INPUT_FROM_HORVATH',
        'PROFILE_ATYPICAL'}
# Historical input defects that a stream repair from the GEO source replaces.
REPAIRABLE_INPUT = {'INPUT_EMPTY', 'INPUT_PROFILE_DETACHED'}
DECISIONS = {'LABEL_SUPPORTED', 'LABEL_INVERTED_CONFIRMED', 'LABEL_UNCONFIRMED'}


def recompute(values, case, control):
    a, b = values.reindex(case).dropna(), values.reindex(control).dropna()
    if design_estimand.constant_output(a, b):
        return design_estimand.constant_summary(a, b)
    out = {'n_case': int(len(a)), 'n_control': int(len(b)),
           'mean_case': float(a.mean()) if len(a) else None, 'mean_control': float(b.mean()) if len(b) else None}
    if len(a) and len(b):
        out['mean_diff'] = float(a.mean() - b.mean())
    if len(a) >= 2 and len(b) >= 2:
        pooled = np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) / (len(a) + len(b) - 2))
        out['cohen_d'] = float(out['mean_diff'] / pooled) if pooled > 0 else None
    if len(a) >= MIN_ARM and len(b) >= MIN_ARM:
        test = stats.ttest_ind(a, b, equal_var=False)
        out['p'] = float(test.pvalue) if np.isfinite(test.pvalue) else None
        out['t'] = float(test.statistic) if np.isfinite(test.statistic) else None
    return out


def study_values(replay_root, gse, clock, cache):
    key = (gse, clock)
    if key not in cache:
        path = Path(replay_root) / 'inputs' / f'{gse}.csv'
        values = None
        if path.exists():
            table = pd.read_csv(path, low_memory=False)
            column = 'gsm_id' if 'gsm_id' in table.columns else next(c for c in table.columns if c.lower() in ('gsm', 'sample_id'))
            for name in (clock, clock + '_mean'):
                if name in table.columns:
                    values = pd.to_numeric(table[name], errors='coerce')
                    values.index = table[column].astype(str)
                    values = values[~values.index.duplicated(keep=False)]
                    break
        if len(cache) > 500:
            cache.clear()
        cache[key] = (values, {'path': str(path), 'sha256': sha256(path)} if path.exists() else None)
    return cache[key]


def repaired_dnam_values(root, gse, clock, cache):
    """Clock values from a stream repair (study re-derived from the GEO source), QC-retained samples only."""
    key = ('repaired', str(root), gse, clock)
    if key not in cache:
        study = Path(root) / 'studies' / gse
        status_path = study / 'status.json'
        result = (None, None, None)
        if status_path.exists():
            status = json.loads(status_path.read_text())
            if status.get('state') == 'REPAIRED':
                table = pd.read_csv(study / 'clock_predictions.csv', index_col=0)
                table = table.rename(columns={'DNAmGrimAge': 'GrimAgeV1', 'DNAmGrimAge.1': 'GrimAgeV2'})
                qc = pd.read_csv(study / 'sample_status.csv', index_col=0)
                retained = qc.index[qc.status == 'RETAINED'].astype(str)
                verdict = next((v for v in status.get('clock_verdicts', []) if v['clock'] == clock), None)
                if clock in table.columns:
                    values = pd.to_numeric(table[clock], errors='coerce')
                    values.index = values.index.astype(str)
                    result = (values.reindex(retained), {'path': str(study / 'clock_predictions.csv'),
                              'sha256': sha256(study / 'clock_predictions.csv'), 'status_sha256': sha256(status_path)},
                              verdict)
        cache[key] = result
    return cache[key]


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(p)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def apply_design(record, values, case, control, series, modality, contrast_types):
    """Replace the Welch result by the design-respecting statistic when the sample records call for one."""
    estimate = design_estimand.estimate(values, case, control, series, modality, record['clock'],
                                        (record['comparison_name'], record['condition_name']), contrast_types)
    if estimate is None:
        return []
    record['design'] = {**estimate['design'], 'metadata': [paths.portable(f) for f in series.get('files', [])]}
    if estimate['recomputed'] is not None:
        record['recomputed_declared'] = record['recomputed']
        record['recomputed'] = estimate['recomputed']
        record['estimand'] = estimate['estimand']
        record['actions'].append({'action': 'RECOMPUTE_DESIGN_ESTIMAND', 'method': estimate['design']['method']})
    return estimate['soft']


def repair_row(row, member, verdict, adjudicated, values, value_source, assay_text, adjudicated_issues=(),
               repaired=None, series=None, modality=None, contrast_types=(), verification_rules=False):
    reasons = set(verdict['reasons'])
    extra_soft = []
    record = {'row_ordinal': verdict['row_ordinal'], 'gse_id': row['gse_id'], 'comparison_name': row['comparison_name'],
              'condition_name': row['condition_name'], 'clock': row['clock'],
              'original': {k: row.get(k) for k in ('statistical_test', 'pvalue', 'effect_size', 'direction', 'fdr')},
              'original_verdict': verdict['verdict'], 'reasons': sorted(reasons), 'actions': [],
              'membership_source': (member or {}).get('source'), 'estimand': ESTIMAND}
    if (member or {}).get('status') != 'UNIQUE':
        record['final_class'] = 'PENDING_DESIGN'
        record['pending_reason'] = sorted(reasons & DESIGN_REASONS) or ['MEMBERSHIP_UNRECOVERED']
        return record
    if member.get('source') == 'design_reconstruction':
        record['actions'].append({'action': 'DESIGN_RECONSTRUCTED_DUAL_MODEL', 'rules': member.get('design_rules')})
        adjudicated = adjudicated or 'LABEL_SUPPORTED'  # arms were assigned to the labels by both models
    if reasons & REPAIRABLE_INPUT and repaired and repaired[0] is not None:
        values, value_source, clock_verdict = repaired
        record['actions'].append({'action': 'USE_REPAIRED_INPUT', 'clock_verdict': clock_verdict})
        # Historical coverage described the historical matrix; the repaired clock verdict now governs coverage.
        reasons = reasons - REPAIRABLE_INPUT - {'CLOCK_COVERAGE_LOW', 'CLOCK_COVERAGE_MID', 'COVERAGE_UNAVAILABLE'}
        if clock_verdict and clock_verdict.get('verdict') == 'WITHHELD':
            record['final_class'] = 'PENDING_INPUT'
            record['pending_reason'] = ['REPAIRED_CLOCK_WITHHELD'] + list(clock_verdict.get('reasons', []))
            return record
        if clock_verdict and clock_verdict.get('verdict') == 'PROVISIONAL':
            extra_soft.append('REPAIRED_CLOCK_PROVISIONAL')
    if reasons & INPUT_REASONS:
        record['final_class'] = 'PENDING_INPUT'
        record['pending_reason'] = sorted(reasons & INPUT_REASONS)
        return record
    if values is None:
        record['final_class'] = 'PENDING_INPUT'
        record['pending_reason'] = ['PREDICTIONS_UNAVAILABLE']
        return record
    case, control = list(member['case_gsms']), list(member['control_gsms'])
    # case_gsms/control_gsms are the arms actually used, after any drop or swap below; the stored arms are kept.
    record['inputs'] = {'predictions': value_source, 'case_gsms': case, 'control_gsms': control,
                        'stored_case_gsms': case, 'stored_control_gsms': control}
    if 'INPUT_PROFILE_DETACHED' in reasons:
        detached = set((verdict.get('evidence') or {}).get('detached_gsms') or [])
        case, control = [g for g in case if g not in detached], [g for g in control if g not in detached]
        record['actions'].append({'action': 'DROP_DETACHED_INPUT_SAMPLES', 'removed_gsms': sorted(detached)})
        if not case or not control:
            record['final_class'] = 'PENDING_INPUT'
            record['pending_reason'] = ['INPUT_PROFILE_DETACHED']
            return record
    if 'NON_STANDARD_ASSAY' in reasons:
        flagged = [g for g in case + control if NON_STANDARD_ASSAY.search(assay_text.get(g, ''))]
        case, control = [g for g in case if g not in flagged], [g for g in control if g not in flagged]
        record['actions'].append({'action': 'DROP_NON_STANDARD_ASSAY_SAMPLES', 'removed_gsms': flagged})
        if not case or not control:
            # The comparison contrasts assays (e.g. BS vs oxBS), not biological groups: nothing left to repair.
            record['final_class'] = 'NOT_REPAIRABLE'
            record['pending_reason'] = ['ARM_ENTIRELY_NON_STANDARD_ASSAY']
            return record
    if adjudicated == 'LABEL_INVERTED_CONFIRMED':
        # Both models placed the labelled condition in the other arm, whatever the heuristic said.
        case, control = control, case
        record['actions'].append({'action': 'SWAP_ARMS_CONFIRMED_INVERSION', 'evidence': 'dual-model adjudication'})
    elif 'LABEL_ORIENTATION_INVERTED' in reasons:
        if adjudicated == 'LABEL_SUPPORTED':
            record['actions'].append({'action': 'KEEP_ORIENTATION_ADJUDICATED_SUPPORTED'})
        else:
            record['final_class'] = 'PENDING_ADJUDICATION'
            record['pending_reason'] = ['LABEL_ORIENTATION_INVERTED']
            return record
    if series is not None:
        artefacts = design_estimand.overflow(values, case + control)
        if artefacts:
            case, control = [g for g in case if g not in artefacts], [g for g in control if g not in artefacts]
            record['actions'].append({'action': 'DROP_OVERFLOW_CLOCK_VALUES', 'removed_gsms': sorted(artefacts)})
            if not case or not control:
                record['final_class'] = 'PENDING_INPUT'
                record['pending_reason'] = ['CLOCK_VALUES_OVERFLOW']
                return record
    record['inputs'].update(case_gsms=case, control_gsms=control)
    result = recompute(values, case, control)
    record['recomputed'] = result
    record['actions'].append({'action': 'RECOMPUTE_DECLARED_ESTIMAND'})
    if 'EFFECT_NOT_REPRODUCED' in reasons:
        record['actions'].append({'action': 'REPLACE_UNREPRODUCIBLE_EFFECT'})
    if series is not None and not result.get('undefined_reason'):
        extra_soft += apply_design(record, values, case, control, series, modality, contrast_types)
        result = record['recomputed']
    soft = sorted(reasons & SOFT) + extra_soft
    if adjudicated in ('LABEL_SUPPORTED', 'LABEL_INVERTED_CONFIRMED') and 'LABEL_UNCONFIRMED' in soft:
        soft.remove('LABEL_UNCONFIRMED')
    elif adjudicated == 'LABEL_UNCONFIRMED' and 'LABEL_UNCONFIRMED' not in soft:
        soft.append('LABEL_UNCONFIRMED')  # the label match was checked and the two models did not confirm it
    # Design concerns both models raised are kept as soft flags: the row stays usable, with the concern visible.
    # For a confirmed inversion, CONTROL_INAPPROPRIATE describes the stored control arm, which is the labelled
    # condition and becomes the case after the swap; it restates the inversion and is not carried.
    restated = {'CASE_NOT_CONDITION', 'OTHER', 'NON_STANDARD_ASSAY'}
    if adjudicated == 'LABEL_INVERTED_CONFIRMED':
        restated.add('CONTROL_INAPPROPRIATE')
    soft += sorted('ADJUDICATED_' + issue for issue in adjudicated_issues if issue not in restated)
    record['soft_flags'] = soft
    if result.get('p') is None:
        record['final_class'] = 'DESCRIPTIVE'
        if result.get('undefined_reason') == design_estimand.CONSTANT_OUTPUT_REASON:
            record['pending_reason'] = [design_estimand.CONSTANT_OUTPUT_REASON]
            record['reasons'] = sorted(set(record['reasons']) | {design_estimand.CONSTANT_OUTPUT_CODE})
            record['actions'].append({'action': 'DESCRIPTIVE_ONLY_CONSTANT_OUTPUTS',
                                      'reason': design_estimand.CONSTANT_OUTPUT_REASON})
            if verification_rules:
                from audit_system.final_verification_rules_2026_11_13 import apply
                apply(record, values, series, modality, contrast_types)
            return record
        # Under a design estimand the units are individuals or pairs; any other missing P keeps the general action.
        few_units = 'recomputed_declared' in record and min(result['n_case'], result['n_control']) < MIN_ARM
        action = 'DESCRIPTIVE_ONLY_UNITS_BELOW_3' if few_units else 'DESCRIPTIVE_ONLY_ARM_BELOW_3'
        record['actions'].append({'action': action, 'n_case': result['n_case'], 'n_control': result['n_control']})
    else:
        record['final_class'] = 'PROVISIONAL' if soft else 'VERIFIED'
    if verification_rules:
        from audit_system.final_verification_rules_2026_11_13 import apply
        apply(record, values, series, modality, contrast_types)
    return record


def mark_duplicates(records):
    """One test per study, clock and pair of arms (after swaps and drops).

    The agent saved some tests under several labels, and design reconstruction can re-arm rows whose historical
    analysis differed (for example covariate-adjusted) with the labelled arms. Repeats point to one representative
    (VERIFIED before PROVISIONAL before DESCRIPTIVE, then the lowest row) and leave the BH family and the counts.
    """
    rank = {'VERIFIED': 0, 'PROVISIONAL': 1, 'DESCRIPTIVE': 2}
    groups = defaultdict(list)
    for r in records:
        if r['final_class'] in rank:
            inputs = r['inputs']
            groups[(r['gse_id'], r['clock'], tuple(sorted(inputs['case_gsms'])),
                    tuple(sorted(inputs['control_gsms'])))].append(r)
    marked = 0
    for members in groups.values():
        members.sort(key=lambda r: (rank[r['final_class']], r['row_ordinal']))
        for r in members[1:]:
            r['duplicate_of'] = members[0]['row_ordinal']
            r['representative_class'] = members[0]['final_class']
            r['final_class'] = 'DUPLICATE'
            r['actions'].append({'action': 'MARK_DUPLICATE_TEST', 'duplicate_of': members[0]['row_ordinal']})
            marked += 1
    return marked


def contrast_types_by_comparison(classification_dir):
    """{(study, comparison, condition): contrast types answered by either classifying model}."""
    out = {}
    path = Path(classification_dir) / 'comparisons.jsonl'
    try:
        lines = open(path)
    except OSError as error:
        raise ValueError(f'Cannot read contrast classification {path}: {error}') from error
    with lines:
        for line in lines:
            entry = json.loads(line)
            types = {a.get('contrast_type') for a in entry.get('answers') or [] if isinstance(a, dict)}
            if entry.get('contrast_type'):
                types.add(entry['contrast_type'])
            out[(entry.get('gse_id'), entry.get('comparison_name'), entry.get('condition_name'))] = sorted(t for t in types if t)
    return out


def individuals_by_comparison(individuals_dir):
    """{(study, comparison, condition): {gsm: individual}} from the two-reader individuals audit (agreed only)."""
    out = {}
    path = Path(individuals_dir) / 'comparisons.jsonl'
    try:
        lines = open(path)
    except OSError as error:
        raise ValueError(f'Cannot read individuals audit {path}: {error}') from error
    with lines:
        for line in lines:
            entry = json.loads(line)
            if entry.get('status') == 'SHARED_AGREED' and entry.get('individuals'):
                out[(entry['gse_id'], entry['comparison_name'], entry['condition_name'])] = entry['individuals']
    return out


def deposited_individuals_by_comparison(modality):
    """Read source-supported grouping decisions without claiming a model review."""
    path = paths.REPO / 'audit_system/policy/individual_grouping_decisions_2026-11-13.json'
    if not path.exists():
        return {}, path
    try:
        policy = json.loads(path.read_text())
        records = {}
        for entry in policy['decisions']:
            if entry['modality'] != modality:
                continue
            if entry['decision'] != 'USE_RECORDED_INDIVIDUALS':
                raise ValueError(f'Unsupported individual grouping decision in {path}')
            key = (entry['gse_id'], entry['comparison_name'], entry['condition_name'])
            if key in records or not entry['individuals']:
                raise ValueError(f'Duplicate or empty individual grouping in {path}: {key}')
            records[key] = entry['individuals']
        return records, path
    except (OSError, ValueError, KeyError) as error:
        raise ValueError(f'Cannot read deposited individual grouping {path}: {error}') from error


def assign_fdr(records):
    """One BH family per modality, with duplicate tests represented once."""
    for record in records:
        for field in ('fdr_family', 'fdr_family_size'):
            record.get('recomputed', {}).pop(field, None)
    inferential = [r for r in records if r['final_class'] in ('VERIFIED', 'PROVISIONAL')]
    fdr = bh([r['recomputed']['p'] for r in inferential]) if inferential else []
    for r, q in zip(inferential, fdr):
        r['recomputed']['fdr_family'] = float(q)
        r['recomputed']['fdr_family_size'] = len(inferential)
    by_ordinal = {r['row_ordinal']: r for r in records}
    for r in records:
        if r['final_class'] == 'DUPLICATE' and 'fdr_family' in by_ordinal[r['duplicate_of']].get('recomputed', {}):
            r['recomputed']['fdr_family'] = by_ordinal[r['duplicate_of']]['recomputed']['fdr_family']  # the same test
    return inferential


def flat_records(records, atlas):
    """Flatten repaired records using the verifier export schema."""
    flat = []
    for r, row in zip(records, atlas.to_dict('records')):
        rec = r.get('recomputed', {})
        flat.append({**row, 'v3_final_class': r['final_class'], 'v3_original_verdict': r['original_verdict'],
                     'v3_reasons': ';'.join(r['reasons']), 'v3_actions': ';'.join(a['action'] for a in r['actions']),
                     'v3_soft_flags': ';'.join(r.get('soft_flags', [])), 'v3_pending': ';'.join(r.get('pending_reason', [])),
                     'v3_n_case': rec.get('n_case'), 'v3_n_control': rec.get('n_control'),
                     'v3_mean_diff': rec.get('mean_diff'), 'v3_cohen_d': rec.get('cohen_d'), 'v3_p': rec.get('p'),
                     'v3_fdr': rec.get('fdr_family'), 'v3_design_method': (r.get('design') or {}).get('method'),
                     'v3_membership_source': r.get('membership_source'),
                     'v3_duplicate_of': r.get('duplicate_of'), 'v3_representative_class': r.get('representative_class')})
    return flat


def run(modality, membership_dir, verdicts_dir, replay_root, out_dir, adjudication_dirs=(), design_dirs=(),
        repaired_dnam_root=None, design=False, classification_dir=None, individuals_dir=None, clock_output_decisions=None,
        verification_rules=True):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    table = CORRECTIONS / f"2026-09-08-{MODALITIES[modality]['table']}-atlas_WORKING_NOT_VALIDATED.csv"
    atlas = pd.read_csv(table, low_memory=False)
    members = {r['row_ordinal']: r for r in map(json.loads, open(Path(membership_dir) / 'rows.jsonl'))}
    designed = 0
    for folder in design_dirs:
        path = Path(folder) / 'rows.jsonl'
        if path.exists():
            for line in open(path):
                d = json.loads(line)
                if (members.get(d['row_ordinal']) or {}).get('status') != 'UNIQUE':
                    members[d['row_ordinal']] = d
                    designed += 1
    verdicts = {r['row_ordinal']: r for r in map(json.loads, open(Path(verdicts_dir) / 'row_verdicts.jsonl'))}
    adjudication, adjudicated_issues = {}, {}
    for folder in adjudication_dirs:
        for path in decision_files(folder):
            for line in open(path):
                d = json.loads(line)
                # Later runs override earlier ones only with a decision; an incomplete or uncached run never erases one.
                if d.get('decision') in DECISIONS or d['row_ordinal'] not in adjudication:
                    adjudication[d['row_ordinal']] = d.get('decision')
                    adjudicated_issues[d['row_ordinal']] = d.get('issues_agreed') or []
    contrasts = contrast_types_by_comparison(classification_dir) if design and classification_dir else {}
    read_individuals = individuals_by_comparison(individuals_dir) if design and individuals_dir else {}
    deposited_individuals, individual_policy = deposited_individuals_by_comparison(modality) if design else ({}, None)
    series_cache = {}
    cache, assay_cache, records = {}, {}, []
    for ordinal, row in enumerate(atlas.to_dict('records')):
        member, verdict = members.get(ordinal), verdicts[ordinal]
        values, source = (None, None)
        if (member or {}).get('status') == 'UNIQUE':
            values, source = study_values(replay_root, row['gse_id'], row['clock'], cache)
        repaired = None
        if repaired_dnam_root and REPAIRABLE_INPUT & set(verdict['reasons']):
            # Stream-repair roots in order (empty matrices, then detached profiles); the first repaired study wins.
            for root in ([repaired_dnam_root] if isinstance(repaired_dnam_root, (str, Path)) else repaired_dnam_root):
                repaired = repaired_dnam_values(root, row['gse_id'], row['clock'], cache)
                if repaired[0] is not None:
                    break
        text = {}
        if modality == 'dnam' and 'NON_STANDARD_ASSAY' in verdict['reasons']:
            if row['gse_id'] not in assay_cache:
                assay_cache[row['gse_id']] = series_sample_text(row['gse_id'])
            text = assay_cache[row['gse_id']]
        series = None
        if design:
            gse = row['gse_id']
            if gse not in series_cache:
                series_cache[gse] = design_estimand.load_series(paths.GEO, MODALITIES[modality]['folder'], gse,
                                                                modality in design_estimand.HUMAN)
                if verification_rules:
                    from audit_system.final_verification_rules_2026_11_13 import metadata
                    series_cache[gse]['verification_metadata'] = metadata(paths.GEO, modality, gse)
            series = series_cache[gse]
            key = (gse, row['comparison_name'], row['condition_name'])
            if not series['subjects'] and key in deposited_individuals:
                series = {**series, 'subjects': deposited_individuals[key], 'subject_field': None,
                          'subject_source': paths.portable(individual_policy),
                          'files': [*series['files'], str(individual_policy)]}
            elif not series['subjects'] and key in read_individuals:
                # No GEO individual field: the individuals both readers of this comparison identified.
                series = {**series, 'subjects': read_individuals[key], 'subject_field': None,
                          'subject_source': 'two-reader individuals audit'}
        records.append(repair_row(row, member, verdict, adjudication.get(ordinal), values, source, text,
                                  adjudicated_issues.get(ordinal, ()), repaired, series, modality,
                                  contrasts.get((row['gse_id'], row['comparison_name'], row['condition_name']), ()),
                                  verification_rules=verification_rules))
    duplicates = mark_duplicates(records)
    if clock_output_decisions:
        from audit_system.clock_output_arms_2026_11_12 import apply, load_decisions
        decisions, _ = load_decisions(clock_output_decisions)
        apply(records, decisions, modality)
    inferential = assign_fdr(records)
    with open(out_dir / 'repairs.jsonl', 'w') as sink:
        for r in records:
            sink.write(json.dumps(r, default=str, sort_keys=True) + '\n')
    flat = flat_records(records, atlas)
    pd.DataFrame(flat).to_csv(out_dir / 'atlas_repaired.csv', index=False)
    classes = pd.Series([r['final_class'] for r in records]).value_counts().to_dict()
    actions = pd.Series([a['action'] for r in records for a in r['actions']]).value_counts().to_dict()
    methods = pd.Series([(r.get('design') or {}).get('method') for r in records]).dropna().value_counts().to_dict()
    manifest = {'schema': SCHEMA, 'modality': modality, 'estimand': ESTIMAND, 'rows': len(records),
                'final_classes': classes, 'actions': actions, 'fdr_family_size': len(inferential),
                'design_estimand': {'applied': bool(design), 'methods': methods,
                                    'module_sha256': sha256(design_estimand.__file__) if design else None,
                                    'classification': str(classification_dir) if design and classification_dir else None,
                                    'individuals': str(individuals_dir) if design and individuals_dir else None},
                'duplicate_tests_marked': duplicates,
                'table_sha256': sha256(table), 'membership_sha256': sha256(Path(membership_dir) / 'rows.jsonl'),
                'verdicts_sha256': sha256(Path(verdicts_dir) / 'row_verdicts.jsonl'), 'module_sha256': sha256(__file__),
                'adjudication_rows': len(adjudication), 'design_reconstructed_rows': designed,
                'repaired_dnam_root': [str(r) for r in ([repaired_dnam_root] if isinstance(repaired_dnam_root, (str, Path))
                                                        else repaired_dnam_root)] if repaired_dnam_root else None,
                'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                'repairs_sha256': sha256(out_dir / 'repairs.jsonl'), 'atlas_sha256': sha256(out_dir / 'atlas_repaired.csv')}
    manifest['final_verification_rules'] = {'enabled': verification_rules}
    if verification_rules:
        from audit_system.final_verification_rules_2026_11_13 import provenance
        manifest['final_verification_rules']['sources'] = provenance()
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True))
    print(json.dumps({'final_classes': classes, 'actions': actions}, indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--modality', choices=sorted(MODALITIES), required=True)
    for name in ('membership', 'verdicts', 'replay', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--adjudication', type=Path, action='append', default=[])
    parser.add_argument('--design', type=Path, action='append', default=[])
    parser.add_argument('--repaired-dnam', type=Path, action='append', help='stream-repair root; repeat for several')
    parser.add_argument('--design-estimand', action='store_true',
                        help='drop overflow clock values and use the design-respecting estimand where the sample records call for it')
    parser.add_argument('--classification', type=Path, help='contrast classification run (age and time contrasts keep the Welch test)')
    parser.add_argument('--individuals', type=Path, help='two-reader individuals audit run')
    parser.add_argument('--clock-output-decisions', type=Path, help='checked comparison decisions with explicit final classes and reasons')
    parser.add_argument('--no-final-verification-rules', action='store_true', help='reproduce the preceding final without output, cell-unit or title-pair checks')
    args = parser.parse_args(argv)
    run(args.modality, args.membership, args.verdicts, args.replay, args.out, args.adjudication, args.design,
        args.repaired_dnam, args.design_estimand, args.classification, args.individuals, args.clock_output_decisions,
        verification_rules=not args.no_final_verification_rules)


if __name__ == '__main__':
    main()
