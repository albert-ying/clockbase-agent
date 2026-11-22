"""Final analysis of the group-choice benchmark (PROTOCOL.md v2): census outcomes, E1 and E2, verifier, Benchmark B.

Run with a Python that has scipy and statsmodels-free dependencies only (stats_a uses numpy and scipy).

  python -m benchmark.published_clock_2026_10_05.groupbench.final_report --root GROUPBENCH --out DIR
"""
import argparse
import json
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from scipy import stats as sstats

from audit_system import paths

from .combine_census import combine, kappa
from .stats_a import estimate, exclusion_bounds

CONTRAST_TYPE = {'age_group': 'age contrast', 'intervention_compound': 'intervention', 'intervention_diet': 'intervention',
                 'intervention_other': 'intervention', 'intervention_genotype': 'genotype', 'condition': 'disease condition',
                 'tissue': 'tissue or cell type', 'cell_type': 'tissue or cell type', 'sex': 'other',
                 'background_strain': 'other', 'geo_field': 'other (GEO field)'}
STATISTICAL = {'PAIRED_AS_INDEPENDENT', 'SMALL_ARM_INFERENCE'}
ORIENTATION = {'REPORTED_RESULT_INVERSION', 'LABEL_MISMATCH'}
EXPORT = Path(os.environ.get('CLOCKBASE_BENCH_EXPORT') or paths.VERIFIER_ROOT / 'final-2026-10-05/public-export/studies')
MOD_DIR = {'A1': 'mm_rna', 'A2': 'dnam'}


GEO_DNAM = paths.GEO / 'GEO_hs_dnam/processed_metadata_preharmonize'


def median_age_gap(unit):
    """Absolute difference of the arms' median GEO ages (years) for a human DNAm unit; None when GEO has no age."""
    import pandas as pd
    path = GEO_DNAM / f"{unit['gse']}.csv"
    if not path.exists():
        return None
    meta = pd.read_csv(path, dtype=str).set_index('GSM')
    for column in [c for c in meta.columns if 'age' in c.lower()]:
        values = pd.to_numeric(meta[column], errors='coerce')
        a = values.reindex(unit['arm_a']).dropna()
        b = values.reindex(unit['arm_b']).dropna()
        if len(a) >= 0.8 * len(unit['arm_a']) and len(b) >= 0.8 * len(unit['arm_b']) and len(a) and len(b):
            return float(abs(a.median() - b.median()))
    return None


def load_released_arms(e2):
    """Released (post-repair) arms of every E2 row, from the public export."""
    cache, out = {}, {}
    for row in e2:
        key = (row['set'], row['gse'])
        if key not in cache:
            path = EXPORT / MOD_DIR[row['set']] / f"{row['gse']}.json"
            cache[key] = {r['row_ordinal']: r for r in json.loads(path.read_text())['rows']} if path.exists() else {}
        ordinal = int(row['unit_id'].rsplit('|', 1)[1])
        exported = cache[key].get(ordinal)
        if exported:
            out[row['unit_id']] = (exported.get('case_gsms') or [], exported.get('control_gsms') or [])
    return out


def load_jsonl(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines()]


def census_by_unit(reviews):
    out = defaultdict(list)
    for review in reviews:
        for unit in review['units']:
            out[unit['unit_id']].append({'brief': review['brief'], 'variants': unit['variants']})
    return out


def contrast_type(unit_id, declared):
    factors = {d['declared_factor'] for k, d in declared.items() if k.startswith(unit_id + '#')}
    types = sorted({CONTRAST_TYPE.get(f, 'other') for f in factors})
    return types[0] if len(types) == 1 else 'mixed'


def a_result(rule, unit_record):
    """Welch test on the declared levels' complete series groups vs the agent's reported call (rule-correct units)."""
    return None  # computed in run_v2 outputs when values exist (see a_result_rows)


def released_rows(set_name, e2, outcomes, units, released_arms):
    """E2: released rows, keeping only error types the census confirmed for the parent unit (inversions are
    recomputed on the released arms and kept as found)."""
    e2_rows = []
    for row in (r for r in e2 if r['set'] == set_name):
        parent = outcomes.get(row['parent'], {})
        if parent.get('outcome') not in ('CONFIRMED', 'CONFIRMED_MISSED_BY_RULES', 'UNRESOLVED', 'OVERTURNED', 'CORRECT'):
            continue
        unit = units[row['parent']]
        released = released_arms.get(row['unit_id'])
        same_arms = released is not None and {frozenset(released[0]), frozenset(released[1])} == {
            frozenset(unit['arm_a']), frozenset(unit['arm_b'])}
        confirmed = set(parent.get('confirmed_types', [])) if parent.get('outcome') in ('CONFIRMED', 'CONFIRMED_MISSED_BY_RULES') else set()
        # Orientation is judged on the released arms (a repair can fix it); other confirmed errors persist when the
        # released arms are the agent's arms, because the recomputation keeps the arms and the unpaired raw test.
        kept = set(row.get('error_types', [])) & (confirmed | {'REPORTED_RESULT_INVERSION', 'SMALL_ARM_INFERENCE'})
        if same_arms:
            kept |= confirmed - ORIENTATION
        e2_rows.append({'gse': row['gse'], 'parent': row['parent'], 'error': int(bool(kept)), 'types': sorted(kept),
                        'final_class': row['final_class'], 'same_arms': same_arms, 'unit_id': row['unit_id']})
    return e2_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--census', type=Path, required=True, help='census reviews JSON')
    parser.add_argument('--b-adjudication', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    e1 = {r['unit_id']: r for r in load_jsonl(root / 'rules-v2/e1_rules.jsonl')}
    e2 = load_jsonl(root / 'rules-v2/e2_rules.jsonl')
    excluded = load_jsonl(root / 'prepared-v2/excluded.jsonl')
    declared = json.loads((root / 'frozen-v2/declared.json').read_text())
    units = {u['unit_id']: u for u in load_jsonl(root / 'prepared-v2/units.jsonl')}
    reviews = census_by_unit(json.loads(args.census.read_text()))
    released_arms = load_released_arms(e2)

    outcomes = {}
    for unit_id, rule in e1.items():
        if rule['class'] not in ('ERROR', 'CORRECT'):
            outcomes[unit_id] = {'outcome': rule['class'], 'rule_class': rule['class']}
            continue
        combined = combine(rule['class'], reviews.get(unit_id, []))
        confirmed_types = sorted(set(rule['error_types']) & {t for t, n in combined['reviewer_error_types'].items() if n >= 2}) \
            if combined['outcome'] == 'CONFIRMED' else sorted(t for t, n in combined['reviewer_error_types'].items() if n >= 2)
        outcomes[unit_id] = {**combined, 'rule_class': rule['class'], 'rule_error_types': rule['error_types'],
                             'confirmed_types': confirmed_types, 'n_reviews': len(reviews.get(unit_id, []))}

    report = {'kappa': kappa([reviews[u] for u in e1 if u in reviews]), 'sets': {}}
    for set_name in ('A1', 'A2'):
        ids = [u for u in e1 if e1[u]['set'] == set_name]
        flow = Counter(outcomes[u]['outcome'] for u in ids)
        assessable = [u for u in ids if outcomes[u]['outcome'] in ('CONFIRMED', 'CONFIRMED_MISSED_BY_RULES', 'UNRESOLVED',
                                                                  'OVERTURNED', 'CORRECT')]
        error = lambda u, upper=False: outcomes[u]['outcome'] in ('CONFIRMED', 'CONFIRMED_MISSED_BY_RULES') or (
            upper and outcomes[u]['outcome'] == 'UNRESOLVED')
        rows = [{'gse': e1[u]['gse'], 'error': int(error(u))} for u in assessable]
        rows_upper = [{'gse': e1[u]['gse'], 'error': int(error(u, True))} for u in assessable]
        design_rows = [{'gse': e1[u]['gse'], 'error': int(error(u) and bool(set(outcomes[u]['confirmed_types']) - STATISTICAL))}
                       for u in assessable]
        # Result-sign record errors: the stored effect size has the opposite sign of the declared groups while the
        # stored direction agrees (deterministic, from the rules). The adjudication packages described the reported
        # direction as relative to the stored case arm, which is false for these units; reviewer calls of label
        # mismatch or inversion on them are reported both ways.
        sign = {u for u in ids if any(r for v in e1[u].get('variants', []) for r in v.get('record_errors', []))}
        artifact = {u for u in assessable if u in sign and outcomes[u]['outcome'] == 'CONFIRMED_MISSED_BY_RULES'
                    and set(outcomes[u]['confirmed_types']) <= {'LABEL_MISMATCH', 'REPORTED_RESULT_INVERSION'}}
        rows_reclassified = [{'gse': e1[u]['gse'], 'error': int(error(u) and u not in artifact)} for u in assessable]
        rows_any = [{'gse': e1[u]['gse'], 'error': int(error(u) or u in sign)} for u in assessable]
        # Sensitivity (stated with the results): age confounding counted only when the arms' median GEO ages differ
        # by more than 5 years, instead of 5 years or a standardized mean difference above 0.5.
        age_only = set()
        if set_name == 'A2':
            for u in assessable:
                types = set(outcomes[u].get('confirmed_types') or [])
                if error(u) and types and types <= {'CONFOUNDED_AGE'}:
                    gap = median_age_gap(units[u])
                    if gap is not None and gap <= 5:
                        age_only.add(u)
        rows_age_strict = [{'gse': e1[u]['gse'], 'error': int(error(u) and u not in age_only)} for u in assessable]
        block = {'units': len(ids), 'flow': dict(flow), 'assessable': len(assessable),
                 'E1': estimate(rows), 'E1_upper_with_unresolved': estimate(rows_upper),
                 'E1_package_artifact_reclassified': estimate(rows_reclassified),
                 'package_artifact_units': sorted(artifact),
                 'result_sign_record_errors': {'units': len(sign & set(assessable)), 'of': len(assessable),
                                               'estimate': estimate([{'gse': e1[u]['gse'], 'error': int(u in sign)} for u in assessable])},
                 'E1_any_error_including_result_sign': estimate(rows_any),
                 'E1_age_confounding_median_gap_over_5y_only': estimate(rows_age_strict),
                 'age_only_units_within_5y': sorted(age_only),
                 'E1_design_errors_only': estimate(design_rows),
                 'error_types': dict(Counter(t for u in assessable if error(u) for t in outcomes[u]['confirmed_types'])),
                 'rule_based': estimate([{'gse': e1[u]['gse'], 'error': int(e1[u]['class'] == 'ERROR')} for u in assessable]),
                 'by_stratum': {}, 'by_contrast_type': {}}
        for stratum in ('replay', 'recovery', 'replay+recovery'):
            sel = [u for u in assessable if stratum in e1[u]['membership_sources']]
            block['by_stratum'][stratum] = {'units': len(sel), 'errors': sum(error(u) for u in sel)}
        for u in assessable:
            t = contrast_type(u, declared)
            entry = block['by_contrast_type'].setdefault(t, {'units': 0, 'errors': 0})
            entry['units'] += 1
            entry['errors'] += int(error(u))
        n_excluded = sum(1 for e in excluded if e['set'] == set_name)
        block['excluded_comparisons'] = n_excluded
        block['exclusion_bounds'] = exclusion_bounds(sum(r['error'] for r in rows), len(rows), n_excluded)
        e2_rows = released_rows(set_name, e2, outcomes, units, released_arms)
        e2_units = {}
        for r in e2_rows:
            e2_units.setdefault(r['parent'], {'gse': r['gse'], 'error': 0})
            e2_units[r['parent']]['error'] |= r['error']
        block['E2_rows'] = estimate(e2_rows) if e2_rows else None
        block['E2_units'] = estimate(list(e2_units.values())) if e2_units else None
        block['E2_error_types'] = dict(Counter(t for r in e2_rows for t in r['types']))
        # Verifier: what happened to the rows of each confirmed E1 error.
        caught = Counter()
        for u in assessable:
            if not error(u):
                continue
            released = [r for r in e2_rows if r['parent'] == u]
            if not released:
                caught['withheld (no released row)'] += 1
            elif any(r['error'] for r in released):
                caught['released with the error'] += 1
            else:
                caught['released after repair'] += 1
        false_alarm = Counter()
        for u in assessable:
            if error(u):
                continue
            released = [r for r in e2_rows if r['parent'] == u]
            false_alarm['withheld' if not released else 'released'] += 1
        block['verifier_on_confirmed_errors'] = dict(caught)
        block['verifier_on_correct_units'] = dict(false_alarm)
        report['sets'][set_name] = block

    if args.b_adjudication and args.b_adjudication.exists():
        adjudications = json.loads(args.b_adjudication.read_text())
        by_item = defaultdict(list)
        for review in adjudications:
            for item in review['items']:
                by_item[item['id']].append(item)
        stats_b = json.loads((root / 'b-v2/b_stats_contained.json').read_text())
        attributed = {}
        for item_id, items in by_item.items():
            votes = Counter(i['category'] for i in items)
            agent_votes = sum(i['agent_error'] for i in items)
            top, count = votes.most_common(1)[0]
            attributed[item_id] = {'categories': dict(votes), 'majority': top if count >= 2 else None,
                                   'agent_error': agent_votes >= 2 or (count < 2)}
        n = stats_b['historical']['n']
        agent_errors = sum(v['agent_error'] for v in attributed.values())
        papers = {d['id']: d['id'].split('-GSE')[0] for d in stats_b['historical']['discordant']}
        report['B'] = {'pairs': n, 'papers': stats_b['historical']['papers'], 'call_agreement': stats_b['historical']['call_agreement'],
                       'call_agreement_interval': stats_b['historical']['call_agreement_interval'],
                       'kappa': stats_b['historical']['kappa'], 'kappa_interval': stats_b['historical']['kappa_interval'],
                       'direction_agreement_both_significant': stats_b['historical']['direction_agreement_both_significant'],
                       'n_both_significant': stats_b['historical']['n_both_significant'],
                       'discordant': len(stats_b['historical']['discordant']), 'attributions': attributed,
                       'agent_error_rate': agent_errors / n if n else None, 'agent_errors': agent_errors,
                       'verified': {k: stats_b['verified'][k] for k in ('n', 'call_agreement', 'kappa',
                                                                        'direction_agreement_both_significant', 'n_both_significant')},
                       'coverage': stats_b['coverage']}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'outcomes.json').write_text(json.dumps(outcomes, indent=1) + '\n')
    (args.out / 'final_report.json').write_text(json.dumps(report, indent=1, default=str) + '\n')
    print(json.dumps(report, indent=1, default=str)[:6000])


if __name__ == '__main__':
    main()
