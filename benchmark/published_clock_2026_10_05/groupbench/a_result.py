"""A-result endpoint (PROTOCOL.md v2, amendment 3): the agent's reported result against the right groups.

For every unit and name variant whose declared contrast is settled, the declared levels' complete groups in the series
(all samples the declared rules select) are compared with a two-sided Welch test on the same atlas clock values, and
the result is set against the agent's reported call for that clock (stored direction; significant when the stored P is
below 0.05). The final report keeps the units the census found correct. Units whose arms are the complete declared
groups give the agent's own test back; the endpoint therefore measures incompleteness of arms and test choice.

With --stages CHECKS_DIR, the same reference is set against three stages of the same row (the agent's reported
call, the released call, the candidate call of PROTOCOL.md addendum item 10), each oriented to the declared case, and a
second reference is added: the declared groups analysed with the candidate estimand (subject means, pairing, age as a
covariate where recorded) instead of the Welch test. Output: FILE with one record per unit, variant and clock.

  python -m benchmark.published_clock_2026_10_05.groupbench.a_result --root G --prepared DIR --frozen DIR --out FILE
      [--stages CHECKS_DIR]
"""
import argparse
import json
from pathlib import Path

from scipy import stats

from .run_v2 import PredictionReader, prepare_unit


def call(direction, p):
    if direction is None or p is None:
        return None
    if p >= 0.05:
        return 'none'
    return 'higher' if direction > 0 else 'lower' if direction < 0 else 'none'


MODALITY = {'A1': 'mm_rna', 'A2': 'dnam'}


def stages(args):
    """Agent, released and candidate calls against the declared groups (Welch and design-respecting references)."""
    from . import design_checks as dc
    units = [json.loads(line) for line in (args.prepared / 'units.jsonl').read_text().splitlines()]
    frozen = {name: json.loads((args.frozen / f'{name}.json').read_text()) for name in ('declared', 'coding', 'refined', 'samples')}
    checks = {m: {(r['gse'], r['row_ordinal']): r for r in map(json.loads, (args.stages / m / 'rows.jsonl').read_text().splitlines())}
              for m in set(MODALITY.values())}
    reader, cache, rows_cache, inputs, out = PredictionReader(), {}, {}, {}, []
    for unit in units:
        prepared, _, rows = prepare_unit(unit, frozen, args.root, cache, reader, rows_cache)
        modality = MODALITY[unit['set']]
        by_ordinal = {r['row_ordinal']: r for r in rows}
        if unit['gse'] not in inputs:
            inputs[unit['gse']] = dc.series_inputs(modality, unit['gse'])
        table, ages, _ = inputs[unit['gse']]
        for variant in prepared['variants']:
            declared = variant['declared']
            if not declared or declared.get('status') == 'ambiguous':
                continue
            case, reference = variant.get('series_case_level') or [], variant.get('series_reference_level') or []
            for clock in variant['clocks']:
                values = clock['values']
                a = [values[g] for g in case if g in values]
                b = [values[g] for g in reference if g in values]
                historical = clock.get('historical') or {}
                agent = call(historical.get('direction'), historical.get('pvalue'))
                record = {'unit_id': unit['unit_id'], 'set': unit['set'], 'gse': unit['gse'], 'variant': variant['index'],
                          'clock': clock['clock'], 'row_ordinal': clock['row_ordinal'], 'agent_call': agent}
                if len(a) < 2 or len(b) < 2:
                    out.append({**record, 'status': 'not_computable'})
                    continue
                test = stats.ttest_ind(a, b, equal_var=False)
                diff = sum(a) / len(a) - sum(b) / len(b)
                record.update(status='ok', expert_welch_call=call(diff, float(test.pvalue)))
                # The declared groups under the candidate estimand: the design-respecting reference.
                subjects = table['subjects'] if table else {}
                found = dc.subject_findings(case, reference, subjects).get('findings', []) if subjects else []
                age = dc.age_finding(case, reference, ages, modality in dc.HUMAN) if ages else None
                usable = bool(age) and modality in dc.HUMAN and clock['clock'] != 'DunedinPACE' \
                    and declared.get('declared_factor') != 'age_group'
                fit = dc.candidate_statistic(found, usable, case, reference, values, subjects, ages,
                                             {'estimate': diff, 'p': float(test.pvalue)})
                record['expert_design_method'] = fit['method']
                record['expert_design_call'] = call(fit.get('estimate'), fit.get('p')) if fit.get('status') == 'OK' else None
                record['stats'] = {'expert_welch': [diff, float(test.pvalue)],
                                   'expert_design': [fit['estimate'], fit['p']] if fit.get('status') == 'OK' else None,
                                   'agent': [historical.get('direction'), historical.get('pvalue')]}
                row = by_ordinal.get(clock['row_ordinal'])
                check = checks[modality].get((unit['gse'], clock['row_ordinal']))
                if row and check and row.get('final_class') in ('VERIFIED', 'PROVISIONAL') and row.get('case_gsms'):
                    released_case = set(map(str, row['case_gsms']))
                    orient = 1 if len(released_case & set(case)) >= len(released_case & set(reference)) else -1
                    recomputed = row.get('recomputed') or {}
                    record['released_call'] = call(orient * recomputed['mean_diff'], recomputed['p']) \
                        if recomputed.get('mean_diff') is not None and recomputed.get('p') is not None else None
                    if record['released_call'] is not None:
                        record['stats']['released'] = [orient * recomputed['mean_diff'], recomputed['p']]
                    candidate = check['candidate']
                    record['candidate_method'] = candidate['method']
                    record['candidate_call'] = call(orient * candidate['estimate'], candidate['p']) \
                        if candidate['status'] == 'OK' and candidate.get('p') is not None else None
                    if record['candidate_call'] is not None:
                        record['stats']['candidate'] = [orient * candidate['estimate'], candidate['p']]
                    record['released'] = True
                    record['final_class'] = row['final_class']
                    record['label_unconfirmed'] = 'LABEL_UNCONFIRMED' in (row.get('soft_flags') or [])
                    record['case_not_condition'] = 'CASE_NOT_CONDITION' in ((row.get('adjudication') or {}).get('issues_agreed') or [])
                else:
                    record['released'] = False
                out.append(record)
    args.out.write_text(''.join(json.dumps(r) + '\n' for r in out))
    from .stats_a import estimate
    strata = {'all released rows': lambda r: True,
              'VERIFIED': lambda r: r.get('final_class') == 'VERIFIED',
              'PROVISIONAL': lambda r: r.get('final_class') == 'PROVISIONAL',
              'labels confirmed': lambda r: r.get('label_unconfirmed') is False,
              'LABEL_UNCONFIRMED': lambda r: r.get('label_unconfirmed') is True,
              'no agreed CASE_NOT_CONDITION': lambda r: r.get('case_not_condition') is False}
    summary = {}
    for set_name in ('A1', 'A2'):
        ok = [r for r in out if r['set'] == set_name and r['status'] == 'ok']
        block = summary.setdefault(set_name, {})
        for reference in ('expert_welch_call', 'expert_design_call'):
            for stage in ('agent_call', 'released_call', 'candidate_call'):
                for name, keep in strata.items():
                    if stage == 'agent_call' and name != 'all released rows':
                        continue
                    members = [r for r in ok if (stage == 'agent_call' or (r.get('released') and keep(r)))
                               and r.get(reference) and r.get(stage)]
                    if not members:
                        continue
                    result = estimate([{'gse': r['gse'], 'error': int(r[reference] != r[stage])} for r in members])
                    result['opposite_significant'] = sum({r[reference], r[stage]} == {'higher', 'lower'} for r in members)
                    label = 'all rows of the agent' if stage == 'agent_call' else name
                    block.setdefault(reference, {}).setdefault(stage, {})[label] = result
                    if name in ('all released rows', 'labels confirmed'):
                        low, high = result['effective']['interval']
                        print(f"{set_name} {reference:19s} {stage:15s} {label:22s} {result['units'] - result['errors']}/{result['units']} agree "
                              f"({100 * (1 - result['rate']):.1f}%; disagreement CI {100 * low:.1f}-{100 * high:.1f}%), opposite {result['opposite_significant']}")
    args.out.with_suffix('.summary.json').write_text(json.dumps(summary, indent=1) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--frozen', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--stages', type=Path, help='design-checks directory with <modality>/rows.jsonl')
    args = parser.parse_args()
    if args.stages:
        stages(args)
        return
    units = [json.loads(line) for line in (args.prepared / 'units.jsonl').read_text().splitlines()]
    frozen = {name: json.loads((args.frozen / f'{name}.json').read_text()) for name in ('declared', 'coding', 'refined', 'samples')}
    reader, cache, rows_cache, out = PredictionReader(), {}, {}, []
    for unit in units:
        prepared, _, _ = prepare_unit(unit, frozen, args.root, cache, reader, rows_cache)
        for variant in prepared['variants']:
            declared = variant['declared']
            if not declared or declared.get('status') == 'ambiguous':
                continue
            case, reference = variant.get('series_case_level') or [], variant.get('series_reference_level') or []
            agent_case = set(unit['arm_a'] if variant['case_is_arm_a'] else unit['arm_b'])
            agent_reference = set(unit['arm_b'] if variant['case_is_arm_a'] else unit['arm_a'])
            for clock in variant['clocks']:
                values = clock['values']
                a = [values[g] for g in case if g in values]
                b = [values[g] for g in reference if g in values]
                historical = clock.get('historical') or {}
                agent = call(historical.get('direction'), historical.get('pvalue'))
                if len(a) < 2 or len(b) < 2 or agent is None:
                    out.append({'unit_id': unit['unit_id'], 'variant': variant['index'], 'clock': clock['clock'],
                                'status': 'not_computable', 'n_case': len(a), 'n_reference': len(b)})
                    continue
                test = stats.ttest_ind(a, b, equal_var=False)
                diff = sum(a) / len(a) - sum(b) / len(b)
                expert = call(diff, float(test.pvalue))
                out.append({'unit_id': unit['unit_id'], 'set': unit['set'], 'gse': unit['gse'], 'variant': variant['index'],
                            'clock': clock['clock'], 'status': 'ok', 'n_case': len(a), 'n_reference': len(b),
                            'expert_groups_p': float(test.pvalue), 'expert_groups_diff': diff, 'expert_call': expert,
                            'agent_call': agent, 'agree': expert == agent,
                            'agent_arms_are_complete_groups': agent_case == set(case) and agent_reference == set(reference)})
    args.out.write_text(''.join(json.dumps(r) + '\n' for r in out))
    ok = [r for r in out if r['status'] == 'ok']
    print(len(out), 'rows;', len(ok), 'computable;', sum(r['agree'] for r in ok), 'agree')


if __name__ == '__main__':
    main()
