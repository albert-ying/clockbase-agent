"""Score the prototype design checks against the census outcomes of the group-choice benchmark (in-sample).

For every released benchmark row (E2 of PROTOCOL.md v2), the census-confirmed error types the checks target
(REPEATED_SUBJECT, PAIRED_AS_INDEPENDENT, CONFOUNDED_AGE) are compared with the findings of design_checks on the same
row. Units are the benchmark units (distinct arm pairs). The benchmark is where these error types were found, so this
measures what the checks catch there; performance on other series is measured separately.

With --validation, the findings both model readers of design_validate gave for a unit are added to the findings from
the recorded fields (the two-layer variant), and both variants are reported.

With --atlas, two atlas-wide tables are written instead: for released rows whose own agent test was paired or carried
age as a covariate and which the checks flag, the agent's call against the released call and the design-respecting
call (agent_design.json); and the check findings of the ranked comparisons (ranked.json).

With --corrected, a released unit also counts as an error when the rules find pairing, repeated subjects or age
imbalance on its released statistic although the agent's own test had handled it (the released Welch test does not);
outputs carry the suffix _corrected.

With --age, the agent's handling of age is tabulated for released human rows (the covariates it declared, the test
it ran, whether its model formula carried age) and every row with recorded ages on a clock that tracks age is
recomputed as clock ~ arm + age and compared with the released call (age_instruction.json).

With --external, the candidate release estimand (PROTOCOL.md addendum item 10) is scored against published results:
the published-clock benchmark T1 with systems S1 (historical row), S2 (released row) and S2c (candidate on the
released row), and Benchmark B with the agent's, the released and the candidate call (external.json).

With --auroc, the significant-finding metrics of PROTOCOL.md addendum item 12 are computed by modality (auroc.json):
direction AUROC, same direction among results significant on both sides, detection AUROC, replication and results
significant only in ClockBase, against published results and against the expert-defined groups.

  python -m benchmark.published_clock_2026_10_05.groupbench.design_eval --root GROUPBENCH --checks DIR --out DIR
      [--validation DIR] [--corrected] [--atlas] [--age] [--external] [--auroc]
"""
import argparse
import json
import os
import re

import numpy as np
from collections import Counter
from pathlib import Path

from .final_report import load_jsonl, load_released_arms, released_rows
from .stats_a import estimate

TARGET = {'REPEATED_SUBJECT': 'REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT': 'PAIRED_AS_INDEPENDENT',
          'CONFOUNDED_AGE': 'AGE_IMBALANCE'}
MODALITY = {'A1': 'mm_rna', 'A2': 'dnam'}
RELEASED = ('VERIFIED', 'PROVISIONAL')


def corrected(rows, e2):
    """Released rows with the design types the rules find on the released statistic added to their error types."""
    rules = {r['unit_id']: set(r.get('error_types', [])) for r in e2}
    out = []
    for row in rows:
        types = set(row['types']) | (rules.get(row['unit_id'], set()) & set(TARGET))
        out.append({**row, 'types': sorted(types), 'error': int(bool(types))})
    return out


def score(set_name, rows, checks, readers=None):
    """Row and unit agreement between census error types and check findings (plus reader findings when given)."""
    by_type = {finding: Counter() for finding in TARGET.values()}
    units = {}
    for row in rows:
        check = checks[(row['gse'], int(row['unit_id'].rsplit('|', 1)[1]))]
        found = set(check['findings'])
        if readers:
            extra = readers.get(check['unit'], set())
            found |= extra - ({'AGE_IMBALANCE'} if check['clock'] == 'DunedinPACE' else set())
        truth = {TARGET[t] for t in row['types'] if t in TARGET}
        other = {t for t in row['types'] if t not in TARGET}
        for finding in TARGET.values():
            by_type[finding][('census' if finding in truth else 'no census') + ', ' +
                             ('check' if finding in found else 'no check')] += 1
        unit = units.setdefault(row['parent'], {'gse': row['gse'], 'error': 0, 'truth': set(), 'found': set(), 'other': set()})
        unit['error'] |= row['error']
        unit['truth'] |= truth
        unit['found'] |= found
        unit['other'] |= other
    errors = [u for u in units.values() if u['error']]
    correct = [u for u in units.values() if not u['error']]
    # Withheld: every flagged unit leaves the release. Recomputed: a flagged finding is repaired, so a unit stays an
    # error when it has an error type outside the three or a targeted type the checks did not find.
    kept = [u for u in units.values() if not u['found']]
    after_recompute = [{'gse': u['gse'], 'error': int(bool(u['error'] and (u['other'] or u['truth'] - u['found'])))}
                       for u in units.values()]
    return {'set': set_name, 'released_rows': len(rows), 'released_units': len(units), 'rows_by_type': {k: dict(v) for k, v in by_type.items()},
            'error_units': len(errors), 'error_units_flagged': sum(bool(u['found']) for u in errors),
            'error_units_fully_addressed': sum(not u['other'] and not (u['truth'] - u['found']) for u in errors),
            'correct_units': len(correct), 'correct_units_flagged': sum(bool(u['found']) for u in correct),
            'as_released': estimate([{'gse': u['gse'], 'error': u['error']} for u in units.values()]),
            'if_flagged_withheld': estimate([{'gse': u['gse'], 'error': u['error']} for u in kept]) if kept else None,
            'if_flagged_recomputed': estimate(after_recompute),
            'missed': sorted((name, sorted(u['truth'] - u['found']), sorted(u['other'])) for name, u in units.items()
                             if u['error'] and (u['other'] or u['truth'] - u['found']))}


def agent_call(row):
    """The agent's stored call in the released orientation (mirrored when the verifier swapped the arms)."""
    historical = row.get('historical') or {}
    direction, p = historical.get('direction'), historical.get('pvalue')
    if direction is None or p is None or p != p:
        return None
    if any(str(a.get('action', '')).startswith('SWAP_ARMS') for a in row.get('actions') or []):
        direction = -direction
    return 'ns' if p >= 0.05 or direction == 0 else ('higher' if direction > 0 else 'lower')


def atlas(checks, out):
    from audit_system import paths
    from ..engine.predictions import final_export
    export = (paths.VERIFIER_ROOT / os.environ['CLOCKBASE_BENCH_SCORED'] / 'public-export'
              if os.environ.get('CLOCKBASE_BENCH_SCORED') else final_export()) / 'studies'
    design = {'agent_model_with_age': ('agent_age_adjusted', 'AGE_IMBALANCE'),
              'agent_test_paired': ('agent_paired', 'PAIRED_AS_INDEPENDENT')}
    tables = {name: Counter() for name in design}
    ranked = {'rows': 0, 'subjects_recorded': 0, 'age_recorded': 0, 'flagged': 0, 'findings': Counter(), 'outcome': Counter(),
              'flagged_rows': []}
    for modality in ('dnam', 'hs_rna', 'mm_rna'):
        by_series = {}
        for record in load_jsonl(checks / modality / 'rows.jsonl'):
            by_series.setdefault(record['gse'], []).append(record)
        for gse, records in by_series.items():
            exported = {r['row_ordinal']: r for r in json.loads((export / modality / f'{gse}.json').read_text())['rows']}
            for record in records:
                row = exported[record['row_ordinal']]
                recomputed = record.get('recomputed') or {}
                for name, (flag, finding) in design.items():
                    if record[flag] and finding in record['findings'] and recomputed.get('status') == 'OK':
                        agent = agent_call(row)
                        if agent is None:
                            tables[name]['agent call not stored'] += 1
                            continue
                        tables[name]['rows'] += 1
                        tables[name]['released call equals agent call'] += record['released']['call'] == agent
                        tables[name]['design-respecting call equals agent call'] += recomputed['call'] == agent
                        if agent == 'ns' and record['released']['call'] != 'ns':
                            tables[name]['agent not significant, released significant'] += 1
                        if agent != 'ns' and record['released']['call'] == 'ns':
                            tables[name]['agent significant, released not significant'] += 1
                ranking = ((row.get('ranking_evidence') or {}).get('ranking') or {})
                if ranking.get('ranked'):
                    ranked['rows'] += 1
                    ranked['subjects_recorded'] += bool(record.get('subject', {}).get('assessable'))
                    ranked['age_recorded'] += bool(record.get('age'))
                    if record['findings']:
                        outcome = record.get('change') or recomputed.get('status') or 'NOT_COMPUTED'
                        ranked['flagged'] += 1
                        ranked['outcome'][outcome] += 1
                        for finding in record['findings']:
                            ranked['findings'][finding] += 1
                        ranked['flagged_rows'].append({'modality': modality, 'gse': gse, 'comparison_name': record['comparison_name'],
                                                       'rank_within_direction': ranking.get('rank_within_direction'),
                                                       'direction': ranking.get('direction_v3'), 'findings': record['findings'],
                                                       'outcome': outcome, 'released_p': record['released']['p'],
                                                       'recomputed_p': recomputed.get('p')})
    out.mkdir(parents=True, exist_ok=True)
    (out / 'agent_design.json').write_text(json.dumps({k: dict(v) for k, v in tables.items()}, indent=1) + '\n')
    ranked['findings'], ranked['outcome'] = dict(ranked['findings']), dict(ranked['outcome'])
    ranked['flagged_rows'].sort(key=lambda r: (r['rank_within_direction'] or 0))
    (out / 'ranked.json').write_text(json.dumps(ranked, indent=1) + '\n')
    print(json.dumps({k: dict(v) for k, v in tables.items()}, indent=1))
    print(json.dumps({k: v for k, v in ranked.items() if k != 'flagged_rows'}, indent=1))


AGE_WORD = re.compile(r'(^|[^a-z])age([^a-z]|$)')
MODEL_TEST = re.compile(r'linear|regress|ancova|ols|mixed|covar|adjust|glm')
PAIRED_TEST = re.compile(r'paired|repeated')


def age_instruction(checks, out):
    """How the agent handled age in released human rows, and the effect of age as a covariate wherever it is recorded."""
    import pandas as pd
    from audit_system import paths
    from ..engine.predictions import PredictionReader
    from . import design_checks as dc
    from ..engine.predictions import final_export
    export = final_export()
    atlas = pd.read_csv(export / 'atlas_v3_combined.csv', low_memory=False,
                        usecols=['modality', 'atlas_row_ordinal', 'historical_covariate', 'historical_covariates',
                                 'historical_statistical_test'])
    names = {'human DNAm': 'dnam', 'human RNA-seq': 'hs_rna', 'mouse RNA-seq': 'mm_rna'}
    declared = {}
    for row in atlas.itertuples(index=False):
        text = ' '.join(str(v).lower() for v in (row.historical_covariate, row.historical_covariates) if isinstance(v, str))
        test = str(row.historical_statistical_test).lower()
        declared[(names[row.modality], int(row.atlas_row_ordinal))] = (
            bool(AGE_WORD.search(text)), 'model' if MODEL_TEST.search(test) else 'paired' if PAIRED_TEST.search(test) else 'two-sample')
    reader, result = PredictionReader(), {}
    for modality in ('dnam', 'hs_rna'):
        records = load_jsonl(checks / modality / 'rows.jsonl')
        tally, changes = Counter(), {'age_imbalanced': Counter(), 'age_balanced': Counter()}
        agreement = Counter()
        by_series = {}
        for record in records:
            says_age, family = declared[(modality, record['row_ordinal'])]
            in_model = record['agent_age_adjusted']
            tally['released rows'] += 1
            tally['agent declared age as a covariate'] += says_age
            tally['agent declared age and ran a model'] += says_age and family == 'model'
            tally['agent declared age and ran a two-sample or paired test'] += says_age and family != 'model'
            tally['agent model formula carries age'] += in_model
            if record.get('age'):
                tally['age recorded for both arms'] += 1
                tally['age recorded, agent model carries age'] += in_model
                tally['age recorded, agent declared age without a model'] += says_age and not in_model
                tally['age recorded, age neither declared nor modelled'] += not says_age and not in_model
                if not record.get('pace_clock') and not record.get('age_declared'):
                    by_series.setdefault(record['gse'], []).append(record)
            if 'AGE_IMBALANCE' in record['findings']:
                tally['age-imbalanced arms'] += 1
                tally['age-imbalanced, agent model carries age'] += in_model
                tally['age-imbalanced, agent declared age without a model'] += says_age and not in_model
                tally['age-imbalanced, age neither declared nor modelled'] += not says_age and not in_model
        for gse, members in by_series.items():
            exported = {r['row_ordinal']: r for r in json.loads((export / 'studies' / modality / f'{gse}.json').read_text())['rows']}
            _, ages, _ = dc.series_inputs(modality, gse)
            cache = {}
            for record in members:
                row = exported[record['row_ordinal']]
                if record['clock'] not in cache:
                    loaded = reader.load(modality, gse, record['clock'], list(exported.values()))
                    cache[record['clock']] = ({str(k): float(v) for k, v in loaded['samples']['clock'].dropna().items()}
                                              if loaded.get('status') == 'OK' else {})
                values = cache[record['clock']]
                case, control = [str(g) for g in row['case_gsms']], [str(g) for g in row['control_gsms']]
                fit = dc.adjusted([(values[g], ages.get(g)) for g in case if g in values],
                                  [(values[g], ages.get(g)) for g in control if g in values]) if values else {'status': 'NO_PREDICTIONS'}
                group = 'age_imbalanced' if 'AGE_IMBALANCE' in record['findings'] else 'age_balanced'
                if fit['status'] != 'OK':
                    changes[group][fit['status']] += 1
                    continue
                call = dc.call(fit['estimate'], fit['p'])
                changes[group][dc.change(record['released']['call'], call)] += 1
                if record['agent_age_adjusted']:
                    agent = agent_call(row)
                    if agent is not None:
                        agreement['rows'] += 1
                        agreement['released call equals agent call'] += record['released']['call'] == agent
                        agreement['age-adjusted call equals agent call'] += call == agent
        result[modality] = {'agent': dict(tally), 'age_as_covariate_wherever_recorded': {k: dict(v) for k, v in changes.items()},
                            'rows_where_agent_model_carried_age': dict(agreement)}
    out.mkdir(parents=True, exist_ok=True)
    (out / 'age_instruction.json').write_text(json.dumps(result, indent=1) + '\n')
    print(json.dumps(result, indent=1))


MODALITY_KEY = {'human DNAm': 'dnam', 'human RNA-seq': 'hs_rna', 'mouse RNA-seq': 'mm_rna'}
SUMMARY_NAME = os.environ.get('CLOCKBASE_BENCH_SUMMARY') or os.environ.get('CLOCKBASE_FINAL_NAME', 'final-2026-10-05')
GUARD_YEARS = 1000.0


def signed_z(estimate, p):
    """Sign of the estimate times the normal quantile of 1 - P/2."""
    from scipy import stats
    if estimate is None or p is None or p != p or estimate != estimate:
        return None
    if estimate == 0:
        return 0.0
    return float(np.copysign(stats.norm.isf(max(float(p), 1e-300) / 2), estimate))


def auc(positive, negative):
    """P(score of a positive > score of a negative), ties one half; None when a class is empty."""
    from scipy import stats
    if not positive or not negative:
        return None
    return float(stats.mannwhitneyu(positive, negative, alternative='two-sided', method='asymptotic').statistic
                 / (len(positive) * len(negative)))


def finding_metrics(items):
    """Item 12 for items {'cluster', 'truth' (higher|lower|none), 'z', 'call' (higher|lower|none)}."""
    def point(members):
        significant = [i for i in members if i['truth'] in ('higher', 'lower')]
        null = [i for i in members if i['truth'] == 'none']
        both = [i for i in significant if i['call'] in ('higher', 'lower')]
        return {'direction_auroc': auc([i['z'] for i in significant if i['truth'] == 'higher'],
                                       [i['z'] for i in significant if i['truth'] == 'lower']),
                'detection_auroc': auc([abs(i['z']) for i in significant], [abs(i['z']) for i in null]),
                'same_direction_both_significant': [sum(i['call'] == i['truth'] for i in both), len(both)],
                'replicated': [sum(i['call'] == i['truth'] for i in significant), len(significant)],
                'sign_matches_significant': [sum((i['z'] > 0) == (i['truth'] == 'higher') and i['z'] != 0 for i in significant),
                                             len(significant)],
                'significant_only_in_clockbase': [sum(i['call'] != 'none' for i in null), len(null)],
                'reference': {c: sum(i['truth'] == c for i in members) for c in ('higher', 'lower', 'none')}}
    out = point(items)
    clusters = {}
    for item in items:
        clusters.setdefault(item['cluster'], []).append(item)
    names, rng = sorted(clusters), np.random.default_rng(20261009)
    draws = {'direction_auroc': [], 'detection_auroc': []}
    if len(names) >= 2:
        for _ in range(2000):
            sample = point([i for k in rng.integers(0, len(names), len(names)) for i in clusters[names[k]]])
            for key in draws:
                if sample[key] is not None:
                    draws[key].append(sample[key])
    for key, values in draws.items():
        out[f'{key}_interval'] = [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))] if values else None
    out['clusters'] = len(names)
    return out


def local_summary_inputs(inputs, bench):
    """Resolve archived benchmark paths under the explicitly selected local bundle."""
    def local(value):
        if isinstance(value, list):
            return [local(item) for item in value]
        if not isinstance(value, str):
            return value
        source = Path(value)
        if source.exists():
            return value
        if bench.name in source.parts:
            candidate = bench.joinpath(*source.parts[source.parts.index(bench.name) + 1:])
            if candidate.exists():
                return str(candidate)
            raise FileNotFoundError(f'Cannot proceed - required data file not found: {candidate}')
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {source}')
    return {key: local(value) for key, value in inputs.items()}


def auroc(root, checks, out):
    """Significant-finding metrics by modality against published results and expert-defined groups."""
    from ..engine.metrics import paper_key
    from ..engine.stats import threshold
    from ..summarize_2026_10_05 import call as three_way, judged, read_records
    bench = root.parent
    candidates = {}
    for modality in MODALITY_KEY.values():
        for record in load_jsonl(checks / modality / 'rows.jsonl'):
            candidates[(modality, record['gse'], record['row_ordinal'])] = record
    inputs = json.loads((bench / 'summaries' / SUMMARY_NAME / 'summary.json').read_text())['inputs']
    inputs = local_summary_inputs(inputs, bench)
    final_ids = {r['id'] for r in load_jsonl(inputs['ground_truth'])}
    records = [r for run in inputs['t1'] for r in read_records(run) if r['id'] in final_ids]
    # Overflow guard (item 13b): published-group results whose estimate is non-physical are recomputed without the
    # overflow clock values and kept as a second system, S3g.
    from ..engine.predictions import PredictionReader, atlas_rows
    from ..engine.stats import evaluate
    reader, guarded = PredictionReader(), []
    for record in records:
        s3 = record.get('systems', {}).get('S3', {})
        record['systems']['S3g'] = s3
        if isinstance(s3.get('mean_diff'), (int, float)) and abs(s3['mean_diff']) > GUARD_YEARS:
            prediction = reader.load(record['modality'], record['gse'], record['clock_clockbase'], atlas_rows(record['modality'], record['gse']))
            samples = prediction['samples'].copy()
            dropped = int((samples['clock'].abs() > GUARD_YEARS).sum())
            samples.loc[samples['clock'].abs() > GUARD_YEARS, 'clock'] = np.nan
            record['systems']['S3g'] = evaluate(samples, record['groups'], record['measure'], threshold(record))
            guarded.append({'id': record['id'], 'gse': record['gse'], 'clock': record['clock_clockbase'], 'samples_dropped': dropped,
                            'before': [s3['mean_diff'], s3.get('p')],
                            'after': [record['systems']['S3g'].get('mean_diff'), record['systems']['S3g'].get('p')]})
    published = {}
    effects = {'S3': [], 'S3g': []}
    effect_dropped = Counter()
    for record in records:
        effect = record.get('effect') or {}
        if effect.get('value') is not None and record['modality'] == 'human DNAm':
            unit = str(effect.get('unit', '')).lower()
            if not isinstance(effect['value'], (int, float)) or isinstance(effect['value'], bool):
                effect_dropped['effect not numeric'] += 1
            elif not unit.startswith('year'):
                effect_dropped['effect not in years'] += 1
            elif 'uncertain' in str(effect.get('type', '')).lower():
                effect_dropped['sign convention uncertain'] += 1
            elif record.get('direction') in ('higher', 'lower') and effect['value'] != 0 \
                    and (effect['value'] > 0) != (record['direction'] == 'higher'):
                effect_dropped['sign contradicts the stated direction'] += 1
            else:
                for system in effects:
                    result = record['systems'].get(system, {})
                    if result.get('covered') is True and isinstance(result.get('mean_diff'), (int, float)):
                        effects[system].append({'cluster': paper_key(record), 'published': float(effect['value']),
                                                'clockbase': float(result['mean_diff']), 'significant': record.get('significant')})
                if record['systems'].get('S3', {}).get('covered') is not True:
                    effect_dropped['no ClockBase estimate for the published groups'] += 1
    for record in records:
        truth = three_way(record.get('direction'), record.get('significant'))
        if truth is None:
            continue
        systems = record.get('systems', {})
        alpha = threshold(record) or 0.05
        entries = {}
        for system in ('S1', 'S2', 'S3', 'S3g'):
            result = systems.get(system, {})
            if result.get('covered') is True and isinstance(result.get('p'), (int, float)) and result.get('direction') in ('higher', 'lower', 'none'):
                sign = {'higher': 1.0, 'lower': -1.0, 'none': 0.0}[result['direction']]
                entries[system] = (signed_z(sign, result['p']), result['direction'] if result['p'] < alpha and sign else 'none')
        s2 = systems.get('S2', {})
        row = candidates.get((MODALITY_KEY[record['modality']], record['gse'], s2.get('row_ordinal'))) if s2.get('covered') is True else None
        if row is not None and row['candidate']['status'] == 'OK' and row['candidate'].get('p') is not None:
            estimate = -row['candidate']['estimate'] if s2.get('swapped') else row['candidate']['estimate']
            direction = 'higher' if estimate > 0 else 'lower' if estimate < 0 else 'none'
            entries['S2c'] = (signed_z(estimate, row['candidate']['p']),
                              direction if row['candidate']['p'] < alpha and direction != 'none' else 'none')
        for system, (z, call_) in entries.items():
            if z is not None:
                published.setdefault((record['modality'], system), []).append(
                    {'cluster': paper_key(record), 'truth': truth, 'z': z, 'call': call_, 'id': record['id']})
    result = {'published': {}, 'expert_groups': {}}
    for (modality, system), items in sorted(published.items()):
        result['published'].setdefault(modality, {})[system] = finding_metrics(items)

    def effect_metrics(pairs):
        from scipy.stats import pearsonr, spearmanr
        def point(members):
            x, y = np.array([m['published'] for m in members]), np.array([m['clockbase'] for m in members])
            if len(x) < 3 or len(set(x)) < 2 or len(set(y)) < 2:
                return None, None
            return float(spearmanr(x, y).statistic), float(pearsonr(x, y).statistic)
        x, y = np.array([m['published'] for m in pairs]), np.array([m['clockbase'] for m in pairs])
        rho, r = point(pairs)
        clusters = {}
        for m in pairs:
            clusters.setdefault(m['cluster'], []).append(m)
        names, rng, draws = sorted(clusters), np.random.default_rng(20261009), {'spearman': [], 'pearson': []}
        for _ in range(2000 if len(names) >= 2 else 0):
            a, b = point([m for k in rng.integers(0, len(names), len(names)) for m in clusters[names[k]]])
            if a is not None and a == a and b == b:
                draws['spearman'].append(a)
                draws['pearson'].append(b)
        interval = lambda v: [float(np.quantile(v, 0.025)), float(np.quantile(v, 0.975))] if v else None
        return {'pairs': len(pairs), 'papers': len(names), 'spearman': rho, 'spearman_interval': interval(draws['spearman']),
                'pearson': r, 'pearson_interval': interval(draws['pearson']),
                'same_sign': [int(np.sum(np.sign(x) == np.sign(y))), len(x)],
                'median_absolute_difference_years': float(np.median(np.abs(x - y))),
                'median_absolute_published_effect_years': float(np.median(np.abs(x)))}
    result['effect_agreement_human_DNAm'] = {system: effect_metrics(pairs) for system, pairs in effects.items()}
    result['effect_agreement_human_DNAm']['left_out'] = dict(effect_dropped)
    result['overflow_guard'] = {'threshold_years': GUARD_YEARS, 'recomputed': guarded}
    shared_ids = None
    for system in ('S1', 'S2', 'S2c'):
        ids = {i['id'] for i in published.get(('human DNAm', system), [])}
        shared_ids = ids if shared_ids is None else shared_ids & ids
    result['published']['human DNAm']['shared_S1_S2_S2c'] = {
        system: finding_metrics([i for i in published[('human DNAm', system)] if i['id'] in shared_ids]) for system in ('S1', 'S2', 'S2c')}

    names = {'A1': 'mouse RNA-seq', 'A2': 'human DNAm'}
    stage_rows = load_jsonl(checks / 'a_result_stages.jsonl')
    for set_name, modality in names.items():
        block = result['expert_groups'].setdefault(modality, {})
        for stage in ('agent', 'released', 'candidate'):
            items = []
            for row in stage_rows:
                stats_ = row.get('stats') or {}
                if row['set'] != set_name or row['status'] != 'ok' or not stats_.get('expert_design') or not stats_.get(stage):
                    continue
                if stage != 'agent' and not row.get('released'):
                    continue
                z = signed_z(stats_[stage][0], stats_[stage][1])
                if z is None or row.get(f'{stage}_call') is None:
                    continue
                items.append({'cluster': row['gse'], 'truth': row['expert_design_call'], 'z': z, 'call': row[f'{stage}_call']})
            block[stage] = finding_metrics(items)
            if set_name == 'A2' and stage != 'agent':
                confirmed = [i for i, row in zip(items, [r for r in stage_rows if r['set'] == set_name and r['status'] == 'ok'
                                                         and (r.get('stats') or {}).get('expert_design') and (r.get('stats') or {}).get(stage)
                                                         and r.get('released') and r.get(f'{stage}_call') is not None
                                                         and signed_z(r['stats'][stage][0], r['stats'][stage][1]) is not None])
                             if row.get('label_unconfirmed') is False]
                block[f'{stage}, labels confirmed'] = finding_metrics(confirmed)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'auroc.json').write_text(json.dumps(result, indent=1) + '\n')

    def line(name, m):
        f = lambda v, i: 'n/a' if v is None else f"{v:.3f}" + (f" [{i[0]:.2f}, {i[1]:.2f}]" if i else '')
        r = m['reference']
        return (f"   {name:30s} ref {r['higher']}/{r['lower']}/{r['none']}  direction AUROC {f(m['direction_auroc'], m['direction_auroc_interval'])}  "
                f"same direction {m['same_direction_both_significant'][0]}/{m['same_direction_both_significant'][1]}  "
                f"detection AUROC {f(m['detection_auroc'], m['detection_auroc_interval'])}  replicated {m['replicated'][0]}/{m['replicated'][1]}  "
                f"ClockBase-only {m['significant_only_in_clockbase'][0]}/{m['significant_only_in_clockbase'][1]}")
    print('overflow guard recomputed', guarded)
    print('effect agreement', json.dumps(result['effect_agreement_human_DNAm'], indent=1))
    for reference in ('published', 'expert_groups'):
        block = result[reference]
        for modality, systems in block.items():
            print(reference, modality)
            for name, m in systems.items():
                if name == 'shared_S1_S2_S2c':
                    for system, mm in m.items():
                        print(line('shared ' + system, mm))
                else:
                    print(line(name, m))


def external(root, checks, out):
    """The candidate estimand against published results (T1 and Benchmark B), with the existing metrics."""
    from ..engine.metrics import paper_key
    from ..engine.stats import threshold
    from ..summarize_2026_10_05 import call as three_way, judged, read_records, shared, t1_metrics, with_interval
    from . import b_stats
    bench = root.parent
    candidates = {}
    for modality in MODALITY_KEY.values():
        for record in load_jsonl(checks / modality / 'rows.jsonl'):
            candidates[(modality, record['gse'], record['row_ordinal'])] = record

    # T1: S2c is the candidate statistic of the released row that S2 uses, oriented like S2.
    inputs = json.loads((bench / 'summaries' / SUMMARY_NAME / 'summary.json').read_text())['inputs']
    inputs = local_summary_inputs(inputs, bench)
    final_ids = {r['id'] for r in load_jsonl(inputs['ground_truth'])}  # the audited ground truth, as in the summary
    records = [r for run in inputs['t1'] for r in read_records(run) if r['id'] in final_ids]
    applied = set()
    for record in records:
        s2 = record.get('systems', {}).get('S2', {})
        if s2.get('covered') is not True:
            record.setdefault('systems', {})['S2c'] = {'covered': False, 'status': s2.get('status')}
            continue
        row = candidates.get((MODALITY_KEY[record['modality']], record['gse'], s2.get('row_ordinal')))
        if row is None:  # a descriptive row: no released statistic to replace
            record['systems']['S2c'] = dict(s2)
            continue
        candidate = row['candidate']
        if candidate['status'] != 'OK' or candidate.get('p') is None:
            record['systems']['S2c'] = {'covered': False, 'status': candidate['status'], 'method': candidate['method']}
            continue
        estimate = -candidate['estimate'] if s2.get('swapped') else candidate['estimate']
        alpha = threshold(record)
        record['systems']['S2c'] = {'covered': True, 'row_ordinal': s2['row_ordinal'], 'p': candidate['p'],
                                    'direction': 'higher' if estimate > 0 else 'lower' if estimate < 0 else 'none',
                                    'significant': bool(candidate['p'] < alpha) if alpha is not None else None,
                                    'method': candidate['method']}
        if not candidate['method'].startswith('Welch t-test (released)'):
            applied.add(record['id'])
    keys = ('direction_agreement', 'call_balanced_accuracy')

    def block(items, systems):
        return {s: with_interval(items, s, t1_metrics, paper_key, keys) for s in systems}

    def published(record):
        return three_way(record.get('direction'), record.get('significant'))

    def system_call(record, system):
        j = judged(record, system)
        return three_way(j['direction'], j['significant']) if j else None

    common = shared(records, ['S2', 'S2c'])

    def gain(items):
        """Call agreement of S2c minus S2 on the same records."""
        a, b = t1_metrics(items, 'S2c')['call_agreement'], t1_metrics(items, 'S2')['call_agreement']
        return None if a is None or b is None else a - b

    groups = {}
    for record in common:
        groups.setdefault(paper_key(record), []).append(record)
    names_, rng, draws = sorted(groups), np.random.default_rng(20261009), []
    for _ in range(2000):
        value = gain([r for i in rng.integers(0, len(names_), len(names_)) for r in groups[names_[i]]])
        if value is not None:
            draws.append(value)
    difference = {'call_agreement_S2c_minus_S2': gain(common), 'papers': len(names_),
                  'interval': [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))],
                  'bootstrap': {'unit': 'paper', 'replicates': 2000, 'seed': 20261009}}
    moved = Counter()
    for record in common:
        truth, before, after = published(record), system_call(record, 'S2'), system_call(record, 'S2c')
        if None in (truth, before, after) or before == after:
            continue
        moved['call changed'] += 1
        moved['now agrees with the paper' if after == truth else 'no longer agrees with the paper' if before == truth
              else 'disagrees before and after'] += 1
    t1 = {'records': len(records), 'own_coverage': block(records, ['S1', 'S2', 'S2c']),
          'shared_S1_S2_S2c': {'records': len(shared(records, ['S1', 'S2', 'S2c'])),
                               'metrics': block(shared(records, ['S1', 'S2', 'S2c']), ['S1', 'S2', 'S2c'])},
          'shared_S2_S2c': {'records': len(common), 'metrics': block(common, ['S2', 'S2c'])},
          'candidate_estimand_applied': {'records': len([r for r in common if r['id'] in applied]),
                                         'metrics': block([r for r in common if r['id'] in applied], ['S1', 'S2', 'S2c'])},
          'S2_to_S2c': dict(moved), 'paired_difference': difference,
          'S2c_not_covered': dict(Counter(r['systems']['S2c'].get('status') for r in records
                                          if r['systems'].get('S2', {}).get('covered') is True and r['systems']['S2c'].get('covered') is not True))}

    # Benchmark B: the same 3x3 statistics for the agent's, the released and the candidate call.
    linked = load_jsonl(root / 'b-v2/linked.jsonl')
    items = {'agent': [], 'released': [], 'candidate': []}
    transitions = Counter()
    for entry in linked:
        link = entry.get('primary')
        if not link or link['tier'] != 'contained' or link['result']['status'] != 'OK':
            continue
        truth = b_stats.published_call(entry)
        agent = b_stats.agent_call(link['result'], 'historical')
        released = b_stats.agent_call(link['result'], 'verified')
        row = candidates.get((MODALITY_KEY[entry['modality']], entry['gse'], link['result']['row_ordinal']))
        if truth is None or agent is None or released is None or row is None or row['candidate']['status'] != 'OK':
            transitions['not scored: ' + ('row not released' if row is None else 'no candidate call' if truth and agent and released else 'no call')] += 1
            continue
        estimate = -row['candidate']['estimate'] if link['swapped'] else row['candidate']['estimate']
        direction = 'higher' if estimate > 0 else 'lower' if estimate < 0 else 'none'
        candidate = direction if row['candidate']['p'] < link['result']['alpha'] and direction != 'none' else 'none'
        paper = b_stats.paper_of(entry['id'])
        for name, value in (('agent', agent), ('released', released), ('candidate', candidate)):
            items[name].append({'paper': paper, 'pair': (truth, value), 'id': entry['id']})
        if candidate != released:
            transitions['call changed from released'] += 1
            transitions['now agrees with the paper' if candidate == truth else 'no longer agrees with the paper' if released == truth
                        else 'disagrees before and after'] += 1
    b = {'transitions_released_to_candidate': dict(transitions)}
    for name, members in items.items():
        pairs = [i['pair'] for i in members]
        entry = b_stats.metrics(pairs)
        entry['papers'] = len({i['paper'] for i in members})
        entry['call_agreement_interval'] = b_stats.cluster_bootstrap(members, lambda sample: sum(x == y for x, y in sample) / len(sample) if sample else None)
        entry['kappa_interval'] = b_stats.cluster_bootstrap(members, b_stats.kappa)
        entry['discordant'] = [i['id'] for i in members if i['pair'][0] != i['pair'][1]]
        b[name] = entry
    result = {'T1': t1, 'B': b}
    out.mkdir(parents=True, exist_ok=True)
    (out / 'external.json').write_text(json.dumps(result, indent=1) + '\n')

    def line(m):
        d = f"{m['direction_agreement']:.3f} ({m['direction_n']})" if m['direction_agreement'] is not None else 'n/a'
        b_ = (f"{m['call_balanced_accuracy']:.3f} [{m['call_balanced_accuracy_interval'][0]:.2f}, {m['call_balanced_accuracy_interval'][1]:.2f}]"
              if m['call_balanced_accuracy'] is not None and m.get('call_balanced_accuracy_interval') else 'n/a')
        return f"covered {m['covered']}/{m['records']} direction {d} call agreement {m['call_agreement']} balanced {b_} significance {m['significance_agreement']} ({m['significance_n']})"

    for name in ('own_coverage', 'shared_S1_S2_S2c', 'shared_S2_S2c', 'candidate_estimand_applied'):
        section = t1[name]
        metrics = section if name == 'own_coverage' else section['metrics']
        print('T1', name, '' if name == 'own_coverage' else f"({section['records']} records)")
        for system, m in metrics.items():
            print('   ', system, line(m))
    print('T1 S2 -> S2c', t1['S2_to_S2c'], '| not covered', t1['S2c_not_covered'], '| difference', t1['paired_difference'])
    for name in ('agent', 'released', 'candidate'):
        e = b[name]
        print('B', name, 'n', e['n'], 'agreement', round(e['call_agreement'], 3), e['call_agreement_interval'], 'kappa', round(e['kappa'], 3),
              e['kappa_interval'], 'direction both significant', e['direction_agreement_both_significant'], e['n_both_significant'])
    print('B transitions', b['transitions_released_to_candidate'])


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--checks', type=Path, required=True, help='directory with <modality>/rows.jsonl')
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--validation', type=Path, help='directory with scored.jsonl of design_validate')
    parser.add_argument('--atlas', action='store_true')
    parser.add_argument('--corrected', action='store_true')
    parser.add_argument('--age', action='store_true')
    parser.add_argument('--external', action='store_true')
    parser.add_argument('--auroc', action='store_true')
    args = parser.parse_args()
    if args.atlas:
        atlas(args.checks, args.out)
        return
    if args.age:
        age_instruction(args.checks, args.out)
        return
    if args.external:
        external(args.root, args.checks, args.out)
        return
    if args.auroc:
        auroc(args.root, args.checks, args.out)
        return
    e2 = load_jsonl(args.root / 'rules-v2/e2_rules.jsonl')
    units = {u['unit_id']: u for u in load_jsonl(args.root / 'prepared-v2/units.jsonl')}
    outcomes = json.loads((args.root / 'final-v2/outcomes.json').read_text())
    arms = load_released_arms(e2)
    readers, read_units = None, set()
    if args.validation:
        scored = load_jsonl(args.validation / 'scored.jsonl')
        readers = {u['unit']: {f for f in TARGET.values() if u[f'reference_{f}']} for u in scored}
        read_units = set(readers)
    result = {}
    for set_name, modality in MODALITY.items():
        checks = {(r['gse'], r['row_ordinal']): r for r in load_jsonl(args.checks / modality / 'rows.jsonl')}
        rows = released_rows(set_name, e2, outcomes, units, arms)
        if args.corrected:
            rows = corrected(rows, e2)
        # Records written by export_as_checks name the class of each row in the scored run: rows that run no longer
        # releases are counted and left out, so the rates describe what the scored run releases.
        scored_class = lambda r: checks[(r['gse'], int(r['unit_id'].rsplit('|', 1)[1]))].get('scored_final_class', 'VERIFIED')
        withheld = Counter(scored_class(r) for r in rows if scored_class(r) not in RELEASED)
        rows = [r for r in rows if scored_class(r) in RELEASED]
        result[set_name] = {**score(set_name, rows, checks), 'rows_not_released_in_scored_run': dict(withheld)}
        if readers is not None:
            block = score(f'{set_name} with readers', rows, checks, readers)
            keys = {checks[(r['gse'], int(r['unit_id'].rsplit('|', 1)[1]))]['unit'] for r in rows}
            block['released_unit_keys'] = len(keys)
            block['released_unit_keys_read'] = len(keys & read_units)
            result[f'{set_name}_with_readers'] = block
    args.out.mkdir(parents=True, exist_ok=True)
    name = ('in_sample_with_readers' if readers is not None else 'in_sample') + ('_corrected' if args.corrected else '') + '.json'
    (args.out / name).write_text(json.dumps(result, indent=1) + '\n')
    for block in result.values():
        print(block['set'], 'units', block['released_units'], 'errors', block['error_units'],
              'flagged', block['error_units_flagged'], 'fully addressed', block['error_units_fully_addressed'],
              'correct flagged', block['correct_units_flagged'], 'of', block['correct_units'])
        for finding, counts in block['rows_by_type'].items():
            print('  rows', finding, counts)
        for label in ('as_released', 'if_flagged_withheld', 'if_flagged_recomputed'):
            e = block[label]
            if e:
                print(f"  {label}: {e['errors']}/{e['units']} primary {e['primary']['interval']} ({e['primary']['method']})")
        for name, targeted, other in block['missed']:
            print('  left:', name, 'not found', targeted, 'other', other)


if __name__ == '__main__':
    main()
