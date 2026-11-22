"""Benchmark summary figure (2026-10-07): agreement with published results, comparison with language models alone,
error rates by stage, where the human DNAm errors sit, the released atlas under the design checks, and the validation
of the checks. Every number is read from the benchmark result files or computed here from them; nothing is simulated.
With --table, the values are written as one summary table by stage (Table_benchmark_2026-10-09.md); with
--scorecard, a two-panel figure of the stages (raw agent, released, candidate estimand) is drawn. Neither redraws
the seven-panel figure.

With CLOCKBASE_BENCH_SCORED naming a verifier run, the tables and the stage figure describe that run: the result
files are those computed on the records export_as_checks wrote for it, where the "released" stage is the earlier
release (CLOCKBASE_FINAL_NAME, a Welch test for every row) and the "candidate" stage is the scored run.

  MPLBACKEND=Agg python3 -m benchmark.published_clock_2026_10_05.groupbench.figure_benchmark_2026_10_07 \
      --bench /Entropy/kying/clockbase-benchmark-2026-10-05 --out benchmark/published_clock_2026_10_05/groupbench/results_2026_10_07
"""
import argparse
import csv
import json
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from ..engine.metrics import paper_key
from .design_eval import TARGET, corrected, local_summary_inputs
from .final_report import load_jsonl, load_released_arms, released_rows
from .stats_a import clopper_pearson, estimate

DESIGN = set(TARGET)
FINAL_NAME = os.environ.get('CLOCKBASE_FINAL_NAME', 'final-2026-10-05')
SUMMARY_NAME = os.environ.get('CLOCKBASE_BENCH_SUMMARY', FINAL_NAME)
SCORED = os.environ.get('CLOCKBASE_BENCH_SCORED')
NAME = 'ClockEffectBench'  # the benchmark's name in the tables of a scored run and in the manuscript
# Scored agent reruns (agentbench/score.py; PROTOCOL.md addendum items 16, 22 and 24), listed in the tidy table of a
# scored run: Codex with four models, Claude Code and Biomni, each scored on the same published results.
RERUN_SETS = [(framework, Path(__file__).resolve().parents[1] / 'agentbench' / folder) for framework, folder in (
    ('Codex', os.environ.get('CLOCKBASE_BENCH_RERUNS', 'results_2026_10_11')), ('Claude Code', 'results_2026_10_18'),
    ('Biomni', 'results_2026_10_19'))]
RERUNS = RERUN_SETS[0][1]
# The paper's groups tested with a Welch test (S3) and with the released test (S3c, S3r; PROTOCOL.md addendum item 19).
WELCH_GROUPS = "Paper's groups, Welch test"
RELEASED_GROUPS = "Paper's groups, released test"
READER_GROUPS = "Paper's groups, released test with individuals read by two models"
MODEL_NAME = {'gpt-6-astra': 'GPT-6 Astra', 'gpt-6-sol': 'GPT-6 Sol', 'gpt-6-luna': 'GPT-6 Luna', 'gpt-5.5': 'GPT-5.5',
              'gpt-4.1': 'GPT-4.1', 'claude-opus-5-5': 'Claude Opus 5.5'}
CHECKS = SCORED or 'v3'
STAGE = ({'agent': 'Raw agent', 'released': 'Verifier, Welch test only', 'candidate': 'Released'} if SCORED else
         {'agent': 'Raw agent', 'released': 'Released', 'candidate': 'Candidate'})
STAMP = SCORED.replace('final-', '') if SCORED else '2026-10-09'
SYSTEMS = {'S3': 'ClockBase', 'S4:openai:gpt-6-astra': 'GPT-6 Astra', 'S4:openai:gpt-4.1': 'GPT-4.1',
           'S4:openai:gpt-5.5': 'GPT-5.5', 'S4:openai:gpt-6-sol': 'GPT-6 Sol', 'S4:openai:gpt-6-luna': 'GPT-6 Luna'}
KINDS = {'CONFOUNDED_AGE': 'Age imbalance', 'PAIRED_AS_INDEPENDENT': 'Paired, tested unpaired',
         'REPEATED_SUBJECT': 'Repeated individual', 'CONFOUNDED': 'Another factor differs',
         'LABEL_MISMATCH': 'Wrong groups', 'IMPURE_ARM': 'Mixed arm',
         'REPORTED_RESULT_INVERSION': 'Direction inverted', 'SMALL_ARM_INFERENCE': 'Arm under 3 samples'}
MODALITY = {'dnam': 'Human DNAm', 'hs_rna': 'Human RNA-seq', 'mm_rna': 'Mouse RNA-seq'}
CHANGED = ('SIGNIFICANCE_LOST', 'SIGNIFICANCE_GAINED', 'DIRECTION_REVERSED')


def read(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise SystemExit(f'Cannot proceed - required data file not found: {path} ({error})')


def rate_row(label, block, hypothetical=False):
    """A forest row from an estimate() block: the observed share and the design-effect-adjusted interval."""
    low, high = block['effective']['interval']
    return {'label': label, 'k': block['errors'], 'n': block['units'], 'low': low, 'high': high, 'open': hypothetical}


def binomial_row(label, k, n):
    low, high = clopper_pearson(k, n)
    return {'label': label, 'k': k, 'n': n, 'low': low, 'high': high, 'open': False}


def collect(bench, results):
    group = bench / 'groupbench'
    checks = group / 'design-checks' / CHECKS
    data = {'sources': {}}

    path = group / 'b-v2/b_stats_contained.json'
    b = read(path)['historical']
    order = ['higher', 'none', 'lower']
    data['published'] = {'matrix': [[b['table_published_by_agent'][r][c] for c in order] for r in order], 'n': b['n'],
                         'papers': b['papers'], 'agreement': b['call_agreement'], 'kappa': b['kappa']}
    data['sources']['published'] = str(path)

    path = bench / 'summaries' / SUMMARY_NAME / 'summary.json'
    shared = read(path)['T1']['shared_S3_S4']
    data['systems'] = {'records': shared['records'], 'papers': shared['papers'],
                       'constant_rule': shared['metrics']['S0']['call_balanced_accuracy'],
                       'constant_rule_direction': shared['metrics']['S0']['direction_agreement'],
                       'rows': [{'label': SYSTEMS[k], 'value': shared['metrics'][k]['call_balanced_accuracy'],
                                 'low': shared['metrics'][k]['call_balanced_accuracy_interval'][0],
                                 'high': shared['metrics'][k]['call_balanced_accuracy_interval'][1],
                                 'call_n': shared['metrics'][k]['call_n'],
                                 'direction': shared['metrics'][k]['direction_agreement'],
                                 'direction_n': shared['metrics'][k]['direction_n'],
                                 'direction_interval': shared['metrics'][k]['direction_agreement_interval']} for k in SYSTEMS]}
    constant = shared['metrics']['S0']['call_balanced_accuracy']
    data['systems']['rows'].append({'label': 'Constant rule', 'value': constant, 'low': constant, 'high': constant,
                                    'open': True})
    data['sources']['systems'] = str(path)

    path = group / 'final-v2/final_report.json'
    final = read(path)['sets']
    plain = read(checks / 'in_sample_corrected.json')
    readers = read(checks / 'in_sample_with_readers_corrected.json')
    data['stages'] = [
        {'header': 'Mouse RNA-seq'},
        rate_row('Raw agent', final['A1']['E1']),
        rate_row('Released', plain['A1']['as_released']),
        rate_row('+ checks, readers', readers['A1_with_readers']['if_flagged_recomputed'], True),
        {'header': 'Human DNAm'},
        rate_row('Raw agent', final['A2']['E1']),
        rate_row('Released', plain['A2']['as_released']),
        rate_row('+ field checks', plain['A2']['if_flagged_recomputed'], True),
        rate_row('+ checks, readers', readers['A2_with_readers']['if_flagged_recomputed'], True)]
    data['sources']['stages'] = [str(path), str(checks / 'in_sample_corrected.json'),
                                 str(checks / 'in_sample_with_readers_corrected.json')]
    # The same rates by stage for the tables; for a scored run the corrected rows are what that run releases.
    data['in_sample'] = {name: {'agent': rate_row(STAGE['agent'], final[name]['E1']),
                                'released': rate_row(STAGE['released'], plain[name]['as_released']),
                                'candidate': rate_row(STAGE['candidate'], plain[name]['if_flagged_recomputed'], not SCORED),
                                'readers': rate_row(f"{STAGE['candidate']} + two readers",
                                                    readers[f'{name}_with_readers']['if_flagged_recomputed'], True)}
                         for name in ('A1', 'A2')}
    data['not_released_in_scored_run'] = {name: plain[name].get('rows_not_released_in_scored_run', {}) for name in ('A1', 'A2')}

    # Released human DNAm units by series and kind of error (corrected error set).
    e2 = load_jsonl(group / 'rules-v2/e2_rules.jsonl')
    units = {u['unit_id']: u for u in load_jsonl(group / 'prepared-v2/units.jsonl')}
    outcomes = read(group / 'final-v2/outcomes.json')
    rows = corrected(released_rows('A2', e2, outcomes, units, load_released_arms(e2)), e2)
    per_unit = {}
    for row in rows:
        unit = per_unit.setdefault(row['parent'], {'gse': row['gse'], 'types': set()})
        unit['types'] |= set(row['types'])
    series = {}
    for unit in per_unit.values():
        entry = series.setdefault(unit['gse'], {'gse': unit['gse'], 'design_only': 0, 'other': 0, 'correct': 0})
        entry['correct' if not unit['types'] else 'design_only' if unit['types'] <= DESIGN else 'other'] += 1
    ranked = sorted(series.values(), key=lambda s: (-(s['design_only'] + s['other'] + s['correct']),
                                                    -(s['design_only'] + s['other']), s['gse']))
    errors = sum(s['design_only'] + s['other'] for s in ranked)
    data['series'] = {'rows': ranked, 'units': len(per_unit), 'errors': errors,
                      'top4_errors': sum(s['design_only'] + s['other'] for s in ranked[:4]),
                      'top4_units': sum(s['design_only'] + s['other'] + s['correct'] for s in ranked[:4]),
                      'estimate': estimate([{'gse': u['gse'], 'error': int(bool(u['types']))} for u in per_unit.values()])['rate']}
    kinds = Counter(t for u in per_unit.values() for t in u['types'])
    data['kinds'] = [{'label': KINDS[t], 'units': n, 'design': t in DESIGN} for t, n in kinds.most_common()]
    data['sources']['series'] = [str(group / 'rules-v2/e2_rules.jsonl'), str(group / 'final-v2/outcomes.json')]

    data['atlas'] = []
    if SCORED:  # the atlas under the prototype checks and the validation of those checks belong to the prototype
        return data
    for modality, label in MODALITY.items():
        path = results / f'atlas_summary_{modality}.json'
        block = read(path)[modality]
        change = block['change_any_finding']
        released = block['released']['rows']
        data['atlas'].append({'label': label, 'released': released, 'flagged': block['flagged']['rows'],
                              'unchanged': change.get('UNCHANGED', 0), 'changed': sum(change.get(k, 0) for k in CHANGED),
                              'not_testable': sum(v for k, v in change.items() if k not in CHANGED + ('UNCHANGED',))})
        data['sources'].setdefault('atlas', []).append(str(path))

    path = results / 'validation.json'
    outside = read(path)['outside']
    confirmed = {f: [sum(outside[m][f]['flagged_confirmed'][k] for m in outside) for k in ('count', 'of')]
                 for f in ('REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT', 'AGE_IMBALANCE')}
    data['validation'] = (
        [{'header': 'Flagged: confirmed'},
         binomial_row('Repeated individual', *confirmed['REPEATED_SUBJECT']),
         binomial_row('Paired arms', *confirmed['PAIRED_AS_INDEPENDENT']),
         binomial_row('Age imbalance', *confirmed['AGE_IMBALANCE']),
         {'header': 'Unflagged: missed'}]
        + [{**binomial_row(MODALITY[m], outside[m]['unflagged_any_missed']['count'], outside[m]['unflagged_any_missed']['of']),
            'missed': True} for m in MODALITY])
    data['sources']['validation'] = str(path)
    return data


def share(k, n, low=None, high=None):
    text = f'{k:,} / {n:,} ({100 * k / n:.1f}%)' if n else f'{k} / {n}'
    return text if low is None else f'{text}; 95% CI {100 * low:.0f}-{100 * high:.0f}%'


def stages(data, bench, results):
    """The staged comparison: raw agent, released atlas, candidate estimand (and with two readers where scored)."""
    ext = read(results / 'external.json')
    t1, b = ext['T1']['shared_S1_S2_S2c']['metrics'], ext['B']
    sample = data['in_sample']
    conclusion = read(results / 'conclusion_stages.json')

    def agree(set_name, stage, stratum, reference='expert_design_call'):
        e = conclusion[set_name][reference][stage][stratum]
        return {'k': e['units'] - e['errors'], 'n': e['units'], 'value': 100 * (1 - e['rate']),
                'low': 100 * (1 - e['effective']['interval'][1]), 'high': 100 * (1 - e['effective']['interval'][0]),
                'unadjusted': [100 * (1 - e['unadjusted_interval'][1]), 100 * (1 - e['unadjusted_interval'][0])],
                'opposite': e['opposite_significant'], 'series': e['series']}

    rows_c = []
    for label, set_name, stratum in (('Mouse RNA-seq', 'A1', 'all released rows'),
                                     ('Human DNAm, labels confirmed', 'A2', 'labels confirmed'),
                                     ('Human DNAm, all released', 'A2', 'all released rows')):
        released, candidate = agree(set_name, 'released_call', stratum), agree(set_name, 'candidate_call', stratum)
        row = {'label': label, 'n': f"{candidate['n']} results", 'released': released['value'], 'candidate': candidate['value'],
               'detail': {'released': released, 'candidate': candidate}}
        if stratum == 'all released rows':
            agent = agree(set_name, 'agent_call', 'all rows of the agent')
            row['agent'] = agent['value']
            row['detail']['agent'] = agent
        rows_c.append(row)
    return {
        'conclusion': rows_c, 'conclusion_source': conclusion,
        'agreement': [
            {'label': 'Call, published results', 'n': f"{ext['T1']['shared_S1_S2_S2c']['records']} results",
             'agent': 100 * t1['S1']['call_agreement'], 'released': 100 * t1['S2']['call_agreement'],
             'candidate': 100 * t1['S2c']['call_agreement']},
            {'label': 'Direction, significant results', 'n': f"{t1['S2c']['direction_n']} results",
             'agent': 100 * t1['S1']['direction_agreement'], 'released': 100 * t1['S2']['direction_agreement'],
             'candidate': 100 * t1['S2c']['direction_agreement']},
            {'label': 'Call, same-sample contrasts', 'n': f"{b['candidate']['n']} contrasts",
             'agent': 100 * b['agent']['call_agreement'], 'released': 100 * b['released']['call_agreement'],
             'candidate': 100 * b['candidate']['call_agreement']}],
        'errors': [{'label': label, 'n': f"{sample[name]['candidate']['n']} comparisons",
                    **{stage: 100 * row['k'] / row['n'] for stage, row in sample[name].items()}}
                   for label, name in (('Mouse RNA-seq', 'A1'), ('Human DNAm', 'A2'))],
        'external': ext}


def table(data, bench, results):
    """One summary table of the benchmark by stage, from the figure data and the result files beside it."""
    group = bench / 'groupbench'
    final = read(group / 'final-v2/final_report.json')
    shared = read(bench / 'summaries' / SUMMARY_NAME / 'summary.json')['T1']['shared_S3_S4']
    systems = {k: shared['metrics'][k] for k in SYSTEMS}
    models = [v for k, v in systems.items() if k != 'S3']
    staged = stages(data, bench, results)
    ext = staged['external']
    t1 = ext['T1']['shared_S1_S2_S2c']
    b = ext['B']
    gain = ext['T1']['paired_difference']
    moved, moved_b = ext['T1']['S2_to_S2c'], b['transitions_released_to_candidate']
    mouse, human = data['in_sample']['A1'], data['in_sample']['A2']
    caught = final['sets']['A1']['verifier_on_confirmed_errors']
    candidate = Counter()
    for modality in MODALITY:
        for record in load_jsonl(group / 'design-checks' / CHECKS / modality / 'rows.jsonl'):
            c = record['candidate']
            significant = record['released']['call'] in ('higher', 'lower')
            candidate['rows'] += 1
            candidate['significant'] += significant
            candidate['applies'] += not c['method'].startswith('Welch t-test (released)')
            if c['status'] != 'OK':
                candidate['no call'] += 1
            elif c['call'] != record['released']['call']:
                candidate['significant lost' if significant else 'gained'] += 1
    source = staged['conclusion_source']
    welch = {'A1': source['A1']['expert_welch_call']['released_call']['all released rows'],
             'A2': source['A2']['expert_welch_call']['released_call']['all released rows'],
             'A2c': source['A2']['expert_welch_call']['released_call']['labels confirmed']}
    label_counts = Counter()
    from audit_system import paths
    atlas_csv = paths.VERIFIER_ROOT / Path(SCORED or FINAL_NAME) / 'public-export/atlas_v3_combined.csv'
    try:
        with atlas_csv.open(newline='') as handle:
            for row in csv.DictReader(handle):
                if row['v3_final_class'] in ('VERIFIED', 'PROVISIONAL'):
                    flagged = 'LABEL_UNCONFIRMED' in row['v3_soft_flags']
                    label_counts['rows'] += 1
                    label_counts['flagged'] += flagged
                    label_counts['dnam'] += flagged and row['modality'] == 'human DNAm'
    except OSError as error:
        raise SystemExit(f'Cannot proceed - required data file not found: {atlas_csv} ({error})')

    def cell(m, key, n_key):
        return f"{100 * m[key]:.0f}% ({m[n_key]})"

    def error(row):
        return f"{row['k']} / {row['n']} ({100 * row['k'] / row['n']:.0f}%; {100 * row['low']:.0f}-{100 * row['high']:.0f}%)"

    def kappa(e):
        return f"{e['kappa']:.2f} ({e['kappa_interval'][0]:.2f}-{e['kappa_interval'][1]:.2f})"

    names = [STAGE[k] for k in ('agent', 'released', 'candidate')]
    header = ' | '.join(names)
    if SCORED:
        intro = [
            f'Generated by `../figure_benchmark_2026_10_07.py --table` with CLOCKBASE_BENCH_SCORED={SCORED}. Stages: **raw agent**',
            f'(its stored result), **verifier, Welch test only** (the verifier\'s repairs with a Welch test for every row; the',
            f'release of {FINAL_NAME}), **released** (the verifier run {SCORED}: the same repairs with the design-respecting',
            'estimand of PROTOCOL.md addendum item 10 - individual means, pairing, age as a covariate where recorded - the',
            'overflow rule, and the two-reader individuals audit of ranking candidates). The estimand was fixed before it',
            'was scored against published results.']
    else:
        intro = [
            'Generated by `../figure_benchmark_2026_10_07.py --table`. Stages: **raw agent** (its stored result), **released**',
            '(the verifier\'s repairs and its Welch recomputation), **candidate** (the released arms with the design-respecting',
            'estimand of PROTOCOL.md addendum item 10: subject means, pairing, age as a covariate where recorded). The',
            'candidate was fixed before it was scored against published results and is not in the released verifier.']
    latest = STAGE['candidate']
    lines = [f"# {NAME + ': summary table' if SCORED else 'Benchmark summary table'} by stage ({STAMP})", ''] + intro + [
        '',
        '## Conclusion matches the analysis of expert-defined groups (higher is better)',
        '',
        'The call of each comparison (direction, P < 0.05) against the call obtained from the complete groups an expert',
        'annotation defines, analysed with the design-respecting estimand. Cells: agreeing / scored (percent; 95% CI',
        'for the percent, design-effect-adjusted; unadjusted in brackets).',
        '',
        f'| Expert-annotated benchmark | {header} |',
        '|---|---|---|---|',
    ] + [f"| {row['label']} | " + ' | '.join(
            (lambda d: f"{d['k']} / {d['n']} ({d['value']:.1f}%; {d['low']:.0f}-{d['high']:.0f}% [{d['unadjusted'][0]:.0f}-{d['unadjusted'][1]:.0f}%])")(row['detail'][k])
            if k in row['detail'] else 'not defined' for k in ('agent', 'released', 'candidate')) + ' |'
         for row in staged['conclusion']] + [
        '',
        f"Opposite significant conclusions: {staged['conclusion'][0]['detail']['agent']['opposite']} of the raw agent's mouse results and "
        f"{staged['conclusion'][2]['detail']['agent']['opposite']} of its human DNAm results; none in the {'rows of the two verifier stages' if SCORED else 'released or candidate rows'}. "
        "\"Labels confirmed\" excludes released rows the verifier itself marked LABEL_UNCONFIRMED "
        f"({label_counts['flagged']:,} of {label_counts['rows']:,} released rows; {label_counts['dnam']:,} of the human DNAm rows). "
        f"Against the Welch reference of protocol amendment 3, the {'rows under the Welch test' if SCORED else 'released rows'} agree in "
        f"{welch['A1']['units'] - welch['A1']['errors']} / {welch['A1']['units']} (mouse) and "
        f"{welch['A2c']['units'] - welch['A2c']['errors']} / {welch['A2c']['units']} (human DNAm, labels confirmed; "
        f"{welch['A2']['units'] - welch['A2']['errors']} / {welch['A2']['units']} for all released rows).",
        '',
        '## A. Agreement with published results (higher is better)',
        '',
        f'| Benchmark | {header} |',
        '|---|---|---|---|',
        f"| Three-way call, {t1['records']} published results on released rows (calls scored) | {cell(t1['metrics']['S1'], 'call_agreement', 'call_n')} | "
        f"{cell(t1['metrics']['S2'], 'call_agreement', 'call_n')} | **{cell(t1['metrics']['S2c'], 'call_agreement', 'call_n')}** |",
        f"| Direction, published significant results | {cell(t1['metrics']['S1'], 'direction_agreement', 'direction_n')} | "
        f"{cell(t1['metrics']['S2'], 'direction_agreement', 'direction_n')} | **{cell(t1['metrics']['S2c'], 'direction_agreement', 'direction_n')}** |",
        f"| Balanced accuracy of the call | {t1['metrics']['S1']['call_balanced_accuracy']:.2f} | {t1['metrics']['S2']['call_balanced_accuracy']:.2f} | "
        f"**{t1['metrics']['S2c']['call_balanced_accuracy']:.2f}** |",
        f"| Three-way call, {b['candidate']['n']} contrasts on the same samples | {100 * b['agent']['call_agreement']:.1f}% | "
        f"{100 * b['released']['call_agreement']:.1f}% | **{100 * b['candidate']['call_agreement']:.1f}%** |",
        f"| Cohen's kappa, same contrasts | {kappa(b['agent'])} | {kappa(b['released'])} | **{kappa(b['candidate'])}** |",
        '',
        (f"Released against the Welch test only on all {ext['T1']['shared_S2_S2c']['records']} results with a row in both: call agreement " if SCORED else
         f"Candidate against released on all {ext['T1']['shared_S2_S2c']['records']} released-row results: call agreement ")
        + f"{100 * ext['T1']['shared_S2_S2c']['metrics']['S2']['call_agreement']:.1f}% to "
        f"{100 * ext['T1']['shared_S2_S2c']['metrics']['S2c']['call_agreement']:.1f}%, a gain of {100 * gain['call_agreement_S2c_minus_S2']:.1f} points "
        f"(95% CI {100 * gain['interval'][0]:.1f} to {100 * gain['interval'][1]:.1f}; bootstrap over {gain['papers']} papers). Of the "
        f"{moved['call changed']} calls that change, {moved['now agrees with the paper']} move to the paper's call and "
        f"{moved.get('no longer agrees with the paper', 0)} away from it; in the same-sample contrasts, "
        f"{moved_b['now agrees with the paper']} of {moved_b['call changed from released']} move to it and "
        f"{moved_b.get('no longer agrees with the paper', 0)} away.",
        '',
        f"ClockBase against language models given only the paper's description ({shared['records']} results, {shared['papers']} papers): "
        f"balanced accuracy {systems['S3']['call_balanced_accuracy']:.2f} "
        f"({systems['S3']['call_balanced_accuracy_interval'][0]:.2f}-{systems['S3']['call_balanced_accuracy_interval'][1]:.2f}) against "
        f"{min(m['call_balanced_accuracy'] for m in models):.2f}-{max(m['call_balanced_accuracy'] for m in models):.2f} for five models and "
        f"{shared['metrics']['S0']['call_balanced_accuracy']:.2f} for a constant rule; direction {systems['S3']['direction_agreement']:.2f} against "
        f"{min(m['direction_agreement'] for m in models):.2f}-{max(m['direction_agreement'] for m in models):.2f}.",
        '',
        '## B. Comparisons with a design error (lower is better)',
        '',
        f"| Expert-annotated benchmark | {header} | {latest} + two readers{' on every comparison' if SCORED else ''} |",
        '|---|---|---|---|---|',
        f"| Mouse RNA-seq ({final['sets']['A1']['E1']['series']} series) | " + ' | '.join(error(mouse[k]) for k in ('agent', 'released', 'candidate'))
        + f" | **{error(mouse['readers'])}** |",
        f"| Human DNAm, disease cohorts ({len(data['series']['rows'])} series) | " + ' | '.join(error(human[k]) for k in ('agent', 'released', 'candidate'))
        + f" | **{error(human['readers'])}** |",
        '',
        f"The verifier withheld or repaired {caught['withheld (no released row)'] + caught['released after repair']} of "
        f"{sum(caught.values())} raw-agent errors in mouse. "
        + ("The columns after the Welch test are in-sample: the estimand was designed from this " if SCORED else "The candidate columns are in-sample: the checks were designed from this ")
        + f"benchmark. In human DNAm, {data['series']['top4_errors']} of the {data['series']['errors']} "
        + ('errors under the Welch test' if SCORED else 'released errors') + " sit in the four largest series."
        + (f" A comparison counts as corrected when the released statistic respects the design the census names; rows the run no longer "
           f"releases are left out ({sum(sum(v.values()) for v in data['not_released_in_scored_run'].values())} rows). The last column adds the "
           "two readers to every comparison; the run reads only ranking candidates." if SCORED else ''),
        '',
        f"## C. What the {'design-respecting estimand' if SCORED else 'candidate'} changes in the {'atlas' if SCORED else 'released atlas'}",
        '',
        '| | Rows |',
        '|---|---|',
        f"| {'Rows released under the Welch test' if SCORED else 'Released rows'} | {candidate['rows']:,} |",
        f"| Rows where the {'released statistic' if SCORED else 'candidate estimand'} differs from the Welch test | {share(candidate['applies'], candidate['rows'])} |",
        f"| Significant rows that lose or reverse their call | {share(candidate['significant lost'], candidate['significant'])} |",
        f"| Rows that become significant | {share(candidate['gained'], candidate['rows'] - candidate['significant'])} |",
        (f"| Rows no longer released (fewer than three individuals in an arm, overflow clock values) | {share(candidate['no call'], candidate['rows'])} |" if SCORED else
         f"| Rows with no candidate call (one subject in an arm, too few pairs, age separating the arms) | {share(candidate['no call'], candidate['rows'])} |"),
        '']
    if SCORED:
        (results / f'Table_benchmark_{STAMP}.md').write_text('\n'.join(lines))
        print('\n'.join(lines))
        return
    validation = [r for r in data['validation'] if 'header' not in r]
    flags, missed = validation[:3], validation[3:]
    headline = read(results / 'headline.json')

    def ranked(label):
        block = headline[label]
        r = block['REPEATED_SUBJECT']
        parts = []
        if r['both_readers'] or r['one_reader']:
            parts.append(f"repeated individuals in {r['both_readers']} (both readers) + {r['one_reader']} (one reader)")
        if block['PAIRED_AS_INDEPENDENT']['both_readers']:
            parts.append(f"paired arms in {block['PAIRED_AS_INDEPENDENT']['both_readers']}")
        if block['AGE_IMBALANCE']['both_readers']:
            parts.append(f"age imbalance in {block['AGE_IMBALANCE']['both_readers']}")
        return f"{block['units_read']} read; " + ('; '.join(parts) if parts else 'no finding')

    lines += [
        '## D. Independent checks',
        '',
        '| | Result |',
        '|---|---|',
        f"| Flags confirmed by two blinded readers, series outside the benchmark | {share(sum(r['k'] for r in flags), sum(r['n'] for r in flags))} |",
    ] + [f"| Unflagged comparisons where both readers find shared individuals, {r['label']} | {share(r['k'], r['n'], r['low'], r['high'])} |"
         for r in missed] + [
        f"| Ranked comparisons read by two readers, Table 1 | {ranked('Table 1')} |",
        f"| Table 2 | {ranked('Table 2')} |",
        f"| Figure 3 | {ranked('Figure 3')} |",
        f"| Table S14 | {ranked('Table S14')} |",
        '']
    (results / f'Table_benchmark_{STAMP}.md').write_text('\n'.join(lines))
    print('\n'.join(lines))


def modality_table(bench, results):
    """The benchmark by modality: direction, significance and effect size against published and expert references."""
    metrics = read(results / 'auroc.json')
    linked = [json.loads(line) for line in (bench / 'groupbench/b-v2/linked.jsonl').read_text().splitlines()]
    mouse = [e for e in linked if e['modality'] == 'mouse RNA-seq' and (e.get('primary') or {}).get('tier') == 'contained'
             and e['primary']['result']['status'] == 'OK']
    mouse_same = sum(e['published']['direction'] == e['primary']['result']['verified']['direction'] for e in mouse)
    mouse_call = sum(e['published'].get('call') == e['primary']['result']['verified']['call'] for e in mouse)

    def auroc(m, key):
        value, interval = m[key], m[f'{key}_interval']
        return 'n/a' if value is None else f"{value:.2f}" + (f" ({interval[0]:.2f}-{interval[1]:.2f})" if interval else '')

    def frac(pair):
        k, n = pair
        return f"{k} / {n}" + (f" ({100 * k / n:.0f}%)" if n else '')

    def null_agree(m):
        k, n = m['significant_only_in_clockbase']
        return [n - k, n]

    def significance(m):
        both = m['same_direction_both_significant'][1]
        agree, null = null_agree(m)
        return [both + agree, m['replicated'][1] + null]

    def block(columns, systems):
        rows = [('Direction', 'Direction AUROC (95% CI)', lambda m: auroc(m, 'direction_auroc')),
                ('', 'Effect sign matches a significant reference effect', lambda m: frac(m['sign_matches_significant'])),
                ('', 'Same direction when both are significant', lambda m: frac(m['same_direction_both_significant'])),
                ('Significance', 'Reference significant: ClockBase significant, same direction', lambda m: frac(m['replicated'])),
                ('', 'Reference not significant: ClockBase not significant', lambda m: frac(null_agree(m))),
                ('', 'Significance agrees overall', lambda m: frac(significance(m))),
                ('', 'Detection AUROC (95% CI)', lambda m: auroc(m, 'detection_auroc'))]
        out = ['| Question | Measure | ' + ' | '.join(columns) + ' |', '|---|---|' + '---|' * len(columns)]
        for question, measure, cell in rows:
            out.append(f"| {('**' + question + '**') if question else ''} | {measure} | " + ' | '.join(cell(m) for m in systems) + ' |')
        return out

    dnam_p, dnam_e = metrics['published']['human DNAm'], metrics['expert_groups']['human DNAm']
    # The published ground truth the results were scored on, counted for the table and for the manuscript.
    truth_path = Path(local_summary_inputs(read(bench / 'summaries' / SUMMARY_NAME / 'summary.json')['inputs'], bench)['ground_truth'])
    truth = load_jsonl(truth_path)
    counts = {'ground_truth': truth_path.name, 'results': len(truth), 'papers': len({paper_key(r) for r in truth}),
              'series': len({r['gse'] for r in truth}), 'by_modality': dict(Counter(r['modality'] for r in truth))}
    (results / 'ground_truth_counts.json').write_text(json.dumps(counts, indent=1) + '\n')
    mouse_e = metrics['expert_groups']['mouse RNA-seq']
    effect, guard = metrics['effect_agreement_human_DNAm'], metrics['overflow_guard']
    released_test = read(results / 'published_released_test.json') if (results / 'published_released_test.json').is_file() else None
    blank = ' |' * (4 if released_test else 3)

    def effect_cell(e, key):
        return f"{e[key]:.2f} ({e[key + '_interval'][0]:.2f}-{e[key + '_interval'][1]:.2f})"

    names = [STAGE[k] for k in ('agent', 'released', 'candidate')]
    lines = [
        f"# {NAME + ': table' if SCORED else 'Benchmark table'} by modality ({STAMP})",
        '',
        'What the benchmark asks: when ClockBase compares an intervention or condition with its control, does it find the',
        'effect the reference reports - the same direction, the same significance, a similar size? Generated by',
        '`../figure_benchmark_2026_10_07.py --modality-table` from `auroc.json` (PROTOCOL.md, addendum items 12 and 13).',
        (f'Stages: raw agent (stored result), verifier with a Welch test for every row (the release of {FINAL_NAME}), released'
         if SCORED else 'Stages: raw agent (stored result), released (current verifier), candidate (design-respecting estimand; not in the'),
        (f'(the verifier run {SCORED}, with the design-respecting estimand). Intervals: bootstrap over papers or series.'
         if SCORED else 'released verifier). Intervals: bootstrap over papers or series.'),
        '',
        '## DNA methylation clocks against published papers',
        '',
        f"Human series; Horvath, PhenoAge, GrimAge, DunedinPACE and related clocks. The ground truth holds {counts['results']} published "
        f"results from {counts['papers']} papers on {counts['series']} series ({counts['by_modality'].get('human DNAm', 0)} for DNA methylation clocks). "
        f"Published findings scored: "
        f"{dnam_p['S3']['reference']['higher']} higher, {dnam_p['S3']['reference']['lower']} lower, {dnam_p['S3']['reference']['none']} not significant "
        f"({dnam_p['S3']['clusters']} papers) for the paper's groups; up to {dnam_p['S2c']['clusters']} papers for atlas rows.",
        ''] + block([WELCH_GROUPS] + ([RELEASED_GROUPS] if released_test else []) + ['Atlas row, ' + n[0].lower() + n[1:] for n in names],
                    [dnam_p['S3']] + ([released_test['systems']['S3c']] if released_test else []) + [dnam_p['S1'], dnam_p['S2'], dnam_p['S2c']]) + [
        f"| **Effect size** | Spearman correlation of effects in years ({effect['S3g']['pairs']} pairs, {effect['S3g']['papers']} papers) | "
        f"{effect_cell(effect['S3g'], 'spearman')}; without the overflow guard {effect_cell(effect['S3'], 'spearman')} |{blank}",
        f"| | Pearson correlation | {effect_cell(effect['S3g'], 'pearson')} |{blank}",
        f"| | Same sign | {frac(effect['S3g']['same_sign'])}; without the guard {frac(effect['S3']['same_sign'])} |{blank}",
        f"| | Median absolute difference | {effect['S3g']['median_absolute_difference_years']:.1f} years "
        f"(median published effect {effect['S3g']['median_absolute_published_effect_years']:.1f} years) |{blank}",
        '',
        f"Effect sizes can be compared for few papers: most report no effect in years. The overflow guard treats clock values above "
        f"{guard['threshold_years']:.0f} years in absolute value as missing; {len(guard['recomputed'])} published-group results in "
        f"{len({g['gse'] for g in guard['recomputed']})} series carried such values (of the order of 1e38 to 1e39 years). "
        + ('The released verifier drops such values before every test (action DROP_OVERFLOW_CLOCK_VALUES).' if SCORED else 'The released atlas has no such guard.'),
        '',
        '## DNA methylation clocks against expert-defined groups',
        '',
        f"{dnam_e['agent']['clusters']} series with expert sample annotation; the reference is the complete expert groups analysed with the "
        'design-respecting estimand. Same samples and clock values as ClockBase: this tests the choice of groups, not replication.',
        ''] + block(names, [dnam_e['agent'], dnam_e['released'], dnam_e['candidate']]) + [
        '',
        '## RNA-seq clock (mouse) against expert-defined groups',
        '',
        f"Transcriptomic mortality clock; {mouse_e['agent']['clusters']} series annotated for Tyshkovskiy et al. (2026). Same samples and clock values as "
        'ClockBase.',
        ''] + block(names, [mouse_e['agent'], mouse_e['released'], mouse_e['candidate']]) + [
        '',
        f"Against the published results of Tyshkovskiy et al. (2026, Fig. 5a) for the same samples: {mouse_same} of {len(mouse)} contrasts in the "
        f"same direction and {mouse_call} of {len(mouse)} with the same call (all significant); too few for an AUROC, and the figure gives no effect "
        'size to compare. Other published RNA-seq clock results use different clocks. Human RNA-seq has no reference and is not benchmarked.',
        '']
    (results / f'Table_benchmark_by_modality_{STAMP}.md').write_text('\n'.join(lines))
    print('\n'.join(lines))


def tidy(data, bench, results):
    """The benchmark as one tidy table (a row per reference, modality, measure and stage) for the supplementary workbook."""
    metrics, ext = read(results / 'auroc.json'), read(results / 'external.json')
    conclusion = read(results / 'conclusion_stages.json')
    rows = []

    def add(reference, modality, measure, stage, k=None, n=None, estimate=None, interval=None, clusters=None, note=''):
        if estimate is None and k is not None and n:
            estimate = k / n
        low, high = interval if interval else (None, None)
        rows.append({'reference': reference, 'modality': modality, 'measure': measure, 'stage': stage, 'numerator': k,
                     'denominator': n, 'estimate': None if estimate is None else round(estimate, 4),
                     'ci_low': None if low is None else round(low, 4), 'ci_high': None if high is None else round(high, 4),
                     'clusters': clusters, 'note': note})

    def finding_rows(reference, modality, stage, m, unit):
        note = f'bootstrap over {unit}'
        add(reference, modality, 'Direction AUROC', stage, estimate=m['direction_auroc'], interval=m['direction_auroc_interval'],
            clusters=m['clusters'], note=note)
        add(reference, modality, 'Effect sign matches a significant reference effect', stage, *m['sign_matches_significant'])
        add(reference, modality, 'Same direction when both are significant', stage, *m['same_direction_both_significant'])
        add(reference, modality, 'Reference significant: ClockBase significant, same direction', stage, *m['replicated'])
        k, n = m['significant_only_in_clockbase']
        add(reference, modality, 'Reference not significant: ClockBase not significant', stage, n - k, n)
        add(reference, modality, 'Detection AUROC', stage, estimate=m['detection_auroc'], interval=m['detection_auroc_interval'],
            clusters=m['clusters'], note=note)

    published = 'Published results'
    dnam = metrics['published']['human DNAm']
    released_test = read(results / 'published_released_test.json') if (results / 'published_released_test.json').is_file() else None
    if released_test:
        finding_rows(published, 'Human DNAm', RELEASED_GROUPS, released_test['systems']['S3c'], 'papers')
        finding_rows(published, 'Human DNAm', READER_GROUPS, released_test['systems']['S3r'], 'papers')
        for stage, key in ((WELCH_GROUPS, 'S3'), (RELEASED_GROUPS, 'S3c'), (READER_GROUPS, 'S3r')):
            for outcome, count in sorted(released_test['published_significant'][key].items()):
                add(published, 'Human DNAm', f'Published significant results: {outcome}', stage, count,
                    sum(released_test['published_significant'][key].values()))
    for key, stage in (('S3', WELCH_GROUPS), ('S1', 'Atlas row, ' + STAGE['agent'].lower()),
                       ('S2', 'Atlas row, ' + STAGE['released'][0].lower() + STAGE['released'][1:]),
                       ('S2c', 'Atlas row, ' + STAGE['candidate'].lower())):
        finding_rows(published, 'Human DNAm', stage, dnam[key], 'papers')
    shared = ext['T1']['shared_S2_S2c']
    for key, stage in (('S2', STAGE['released']), ('S2c', STAGE['candidate'])):
        m = shared['metrics'][key]
        add(published, 'Human DNAm', 'Three-way call agrees (higher, lower, not significant)', 'Atlas row, ' + stage[0].lower() + stage[1:],
            round(m['call_agreement'] * m['call_n']), m['call_n'], note=f"{shared['records']} results with a row at both verifier stages")
    gain = ext['T1']['paired_difference']
    add(published, 'Human DNAm', 'Three-way call agreement, difference between the verifier stages', STAGE['candidate'] + ' minus ' + STAGE['released'][0].lower() + STAGE['released'][1:],
        estimate=gain['call_agreement_S2c_minus_S2'], interval=gain['interval'], clusters=gain['papers'], note='bootstrap over papers')
    effect = metrics['effect_agreement_human_DNAm']
    for key, note in (('S3g', 'overflow clock outputs excluded'), ('S3', 'overflow clock outputs kept')):
        e = effect[key]
        stage = WELCH_GROUPS
        add(published, 'Human DNAm', 'Effect in years: Spearman correlation', stage, n=e['pairs'], estimate=e['spearman'],
            interval=e['spearman_interval'], clusters=e['papers'], note=note + '; bootstrap over papers')
        add(published, 'Human DNAm', 'Effect in years: same sign', stage, *e['same_sign'], note=note)
        add(published, 'Human DNAm', 'Effect in years: median absolute difference (years)', stage, n=e['pairs'],
            estimate=e['median_absolute_difference_years'], note=note + f"; median published effect {e['median_absolute_published_effect_years']} years")

    constant = 'Constant rule: higher and significant'
    if released_test:
        against = released_test['against_models_without_data']
        note = (f"{against['records']} results from {against['papers']} papers answered by every system "
                f"({against['published_significant']['higher']} published increases, {against['published_significant']['lower']} decreases); "
                'bootstrap over papers')

        def baseline_name(name):
            if name.startswith('S4:'):
                return f"Language model without data: {MODEL_NAME[name.split(':')[-1]]}"
            return {'S0': constant, 'S3': WELCH_GROUPS, 'S3c': RELEASED_GROUPS, 'S3r': READER_GROUPS}[name]

        for name, m in against['systems'].items():
            stage = baseline_name(name)
            add(published, 'Human DNAm', 'Systems: three-way call agrees', stage, round(m['call_agreement'] * m['call_n']), m['call_n'],
                interval=m.get('call_agreement_interval'), clusters=against['papers'], note=note)
            add(published, 'Human DNAm', 'Systems: three-way call, balanced accuracy', stage, n=m['call_n'],
                estimate=m['call_balanced_accuracy'], interval=m.get('call_balanced_accuracy_interval'), clusters=against['papers'], note=note)
            add(published, 'Human DNAm', 'Systems: published increases identified', stage, *m['published_higher_identified'], note=note)
            add(published, 'Human DNAm', 'Systems: published decreases identified', stage, *m['published_lower_identified'], note=note)
            if m.get('direction_auroc_signed_confidence') is not None:
                add(published, 'Human DNAm', "Systems: direction AUROC from the model's stated confidence", stage,
                    estimate=m['direction_auroc_signed_confidence'], note=note)
    else:
        systems = data['systems']
        shared_note = f"{systems['records']} results from {systems['papers']} papers scored for every system; bootstrap over papers"
        for row in systems['rows']:
            if row.get('open'):  # the constant rule, added below with its direction
                continue
            stage = (WELCH_GROUPS if row['label'] == 'ClockBase'
                     else f"Language model without data: {row['label']}")
            add(published, 'Human DNAm', 'Three-way call, balanced accuracy', stage, n=row['call_n'], estimate=row['value'],
                interval=[row['low'], row['high']], clusters=systems['papers'], note=shared_note)
            add(published, 'Human DNAm', 'Direction agrees, published significant results', stage,
                round(row['direction'] * row['direction_n']), row['direction_n'], interval=row['direction_interval'],
                clusters=systems['papers'], note=shared_note)
        add(published, 'Human DNAm', 'Three-way call, balanced accuracy', constant, estimate=systems['constant_rule'], note=shared_note)
        add(published, 'Human DNAm', 'Direction agrees, published significant results', constant,
            round(systems['constant_rule_direction'] * systems['rows'][0]['direction_n']), systems['rows'][0]['direction_n'], note=shared_note)


    for framework, folder in RERUN_SETS:
        if not (SCORED and (folder / 'agent_reruns.json').is_file()):
            continue
        reruns = read(folder / 'agent_reruns.json')
        first = folder == RERUNS  # the baselines (S0-S4) are listed once, from the Codex scoring

        def system_name(system, framework=framework):
            if system == 'S0':
                return constant
            if system == 'S1':
                return 'Historical agent (CodingAgent with Claude 3.7 Sonnet), stored result'
            if system == 'S2c':
                return 'Released row of the historical agent'
            if system == 'S3':
                return WELCH_GROUPS
            kind, model = system.split(':')[0], system.split(':')[-2 if system.startswith('R:') else -1]
            if kind == 'S4':
                return f'Language model without data: {MODEL_NAME[model]}'
            return f"{framework} agent with {MODEL_NAME[model]}, " + ('reported result' if system.endswith(':raw') else 'verifier recomputation')

        def rerun_rows(block, note, system_name=system_name):
            for system, m in block.items():
                stage, f = system_name(system), m.get('findings') or {}
                add(published, 'Human DNAm', 'Agent reruns: published results covered', stage, m['covered'], m['records'], note=note)
                if m['call_n']:
                    add(published, 'Human DNAm', 'Agent reruns: three-way call agrees', stage, round(m['call_agreement'] * m['call_n']),
                        m['call_n'], interval=m.get('call_agreement_interval'), clusters=m.get('units'), note=note)
                if m['direction_n']:
                    add(published, 'Human DNAm', 'Agent reruns: direction agrees, published significant results', stage,
                        round(m['direction_agreement'] * m['direction_n']), m['direction_n'],
                        interval=m.get('direction_agreement_interval'), clusters=m.get('units'), note=note)
                if f.get('direction_auroc') is not None:
                    add(published, 'Human DNAm', 'Agent reruns: direction AUROC', stage, estimate=f['direction_auroc'],
                        interval=f.get('direction_auroc_interval'), clusters=f.get('clusters'), note=note)

        base = (f"{reruns['records']} published results from {reruns['papers']} papers in {reruns['series']} series on the five clocks "
                'of the task; bootstrap over papers')
        own = reruns['own_coverage'] if first else {k: v for k, v in reruns['own_coverage'].items() if k.startswith('R:')}
        rerun_rows(own, base + '; each system on the results it covers')
        shared_block = reruns['shared_agent_systems']
        shared = (shared_block['metrics'] if first else
                  {k: v for k, v in shared_block['metrics'].items() if k in ('S0', 'S1', 'S2c') or k.startswith('R:')})
        rerun_rows(shared, base + f"; the {shared_block['records']} results covered by the historical agent, "
                   + ('its released row and every rerun' if first else f'its released row and the {framework} rerun'))

    expert = 'Expert-defined groups'
    for set_name, modality in (('A1', 'Mouse RNA-seq'), ('A2', 'Human DNAm')):
        block = conclusion[set_name]['expert_design_call']
        for stage_key, stage, strata in (('agent_call', STAGE['agent'], ('all rows of the agent',)),
                                         ('released_call', STAGE['released'], ('all released rows', 'labels confirmed')),
                                         ('candidate_call', STAGE['candidate'], ('all released rows', 'labels confirmed'))):
            for stratum in strata:
                e = block[stage_key][stratum]
                add(expert, modality, 'Conclusion agrees (direction and significance)', stage, e['units'] - e['errors'], e['units'],
                    interval=[1 - e['effective']['interval'][1], 1 - e['effective']['interval'][0]], clusters=e['series'],
                    note=f"{stratum}; Clopper-Pearson on the design-effect-adjusted sample size; {e['opposite_significant']} significant in the opposite direction")
        groups = metrics['expert_groups'][modality[0].lower() + modality[1:] if modality.startswith('Mouse') else 'human DNAm']
        for key, stage in (('agent', STAGE['agent']), ('released', STAGE['released']), ('candidate', STAGE['candidate'])):
            finding_rows(expert, modality, stage, groups[key], 'series')
        sample = data['in_sample'][set_name]
        for key in ('agent', 'released', 'candidate', 'readers'):
            row = sample[key]
            add(expert, modality, 'Comparisons with a design flaw of any kind (strict count)', row['label'], row['k'], row['n'],
                interval=[row['low'], row['high']], note='in-sample: the estimand was designed from this set'
                + ('; two readers applied to every comparison (the run reads ranking candidates)' if key == 'readers' else ''))
    path = results / f'Table_benchmark_{STAMP}.csv'
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'{len(rows)} rows -> {path}')


def scorecard(data, bench, results, utils):
    """Two panels: agreement with published results and design errors, by stage."""
    import matplotlib.pyplot as plt
    U = utils
    U.retro_style()
    BLUE, GRAY, INK = U.BLUE, U.GRAY, U.INK
    staged = stages(data, bench, results)
    # Rings are drawn larger than the filled marks and the released mark is a square, so coinciding stages stay visible.
    styles = {'agent': dict(label=STAGE['agent'], marker='o', ms=7.6, mfc='none', mec=GRAY['400'], mew=1.0),
              'released': dict(label=STAGE['released'], marker='s', ms=5.4, mfc=GRAY['300'], mec=INK, mew=0.4),
              'candidate': dict(label='Released' if SCORED else 'Candidate estimand', marker='o', ms=4.6, mfc=BLUE['400'], mec=INK, mew=0.4),
              'readers': dict(label=f"{STAGE['candidate']} + two readers", marker='o', ms=7.6, mfc='none', mec=BLUE['400'], mew=1.0)}
    fig, axd = plt.subplot_mosaic([['a'], ['b'], ['c']], figsize=U.figsize_nature('1.5col', aspect=1.02), layout='constrained',
                                  gridspec_kw=dict(height_ratios=[3.3, 3.3, 2.5]))
    fig.get_layout_engine().set(h_pad=0.10)
    for key, rows, xlabel, title, ticks in (
            ('a', staged['conclusion'], 'Conclusion matches the expert-group analysis (%)', 'Expert-defined groups', [0, 25, 50, 75, 100]),
            ('b', staged['agreement'], 'Agreement with the published result (%)', 'Published results', [0, 25, 50, 75, 100]),
            ('c', staged['errors'], 'Comparisons with any design error (%)', 'Design errors', [0, 25, 50, 75, 100])):
        ax = axd[key]
        ys = np.arange(len(rows))
        for y, row in zip(ys, rows):
            present = [k for k in styles if k in row]
            xs = [row[k] for k in present]
            ax.plot([min(xs), max(xs)], [y, y], color=GRAY['300'], lw=0.9, zorder=1, solid_capstyle='butt')
            for k in present:
                style = {a: b for a, b in styles[k].items() if a != 'label'}
                ax.plot([row[k]], [y], ls='none', zorder=3, clip_on=False, **style)
        for k, style in styles.items():
            if any(k in row for row in rows):
                ax.plot([], [], ls='none', **style)
        ax.set_yticks(ys)
        ax.set_yticklabels([f"{row['label']}\n({row['n']})" for row in rows])
        for tick in ax.get_yticklabels():
            tick.set_linespacing(1.05)
        ax.invert_yaxis()
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y', length=0, pad=3.5)
        ax.set_xlim(ticks[0], ticks[-1])
        ax.set_xticks(ticks)
        ax.spines['bottom'].set_bounds(ticks[0], ticks[-1])
        U.pad_categorical_axis(ax, 'y', pad=0.6)
        order = [k for k in styles if any(k in row for row in rows)]
        values = [' | '.join(f"{row[k]:.0f}" if k in row else '\u2013' for k in order) for row in rows]
        U.add_row_value_column(ax, ys, values, header='% by stage', x=1.03, fontsize=U.TYPE['tick'])
    for k in ('readers',):  # the legend of panel a names every stage, including the one drawn only in panel b
        axd['a'].plot([], [], ls='none', **styles[k])
    U.nature_legend(axd['a'], title='Stage', loc='center left')
    U.add_panel_labels(fig, axd)
    U.save_figure(fig, str(results / f'Figure_benchmark_stages_{STAMP}.pdf'), axes_dict=axd,
                  footer=f'Data and intervals: Table_benchmark_{STAMP}.md',
                  review_png=str(results / f'Figure_benchmark_stages_{STAMP}_review.png'))
    fig.savefig(results / f'Figure_benchmark_stages_{STAMP}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def forest(ax, rows, colors, xlabel, value_header, utils, percent=True):
    """Point-range rows with bold group headers on the category axis and a k / n column on the right."""
    scale = 100 if percent else 1
    ys, labels, headers, values = [], [], [], []
    for y, row in enumerate(rows):
        ys.append(y)
        if 'header' in row:
            labels.append(row['header'])
            headers.append(y)
            values.append('')
            continue
        labels.append(row['label'])
        x = scale * (row['k'] / row['n'] if 'k' in row else row['value'])
        color = colors(row)
        ax.errorbar([x], [y], xerr=[[x - scale * row['low']], [scale * row['high'] - x]], fmt='none', ecolor=utils.INK,
                    zorder=2, **utils.ERRORBAR['forest'])
        ax.plot([x], [y], marker='o', ms=4.6, mfc='white' if row.get('open') else color, mec=color if row.get('open') else utils.INK,
                mew=1.0 if row.get('open') else 0.4, ls='none', zorder=3, clip_on=False)
        values.append(f"{row['k']} / {row['n']}" if 'k' in row else f"{row['value']:.2f}")
    ax.set_yticks(ys)
    ax.set_yticklabels(labels)
    for tick, y in zip(ax.get_yticklabels(), ys):
        if y in headers:
            tick.set_fontweight('bold')
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=3.5)
    return ys, values


def finish_forest(ax, ys, values, value_header, utils, xticks):
    ax.set_xlim(xticks[0], xticks[-1])
    ax.set_xticks(xticks)
    ax.spines['bottom'].set_bounds(xticks[0], xticks[-1])
    utils.pad_categorical_axis(ax, 'y', pad=0.5)
    if value_header:
        utils.add_row_value_column(ax, ys, values, header=value_header, x=1.03, fontsize=utils.TYPE['tick'], header_y=1.0)


def draw(data, out, utils, draft=False):
    import matplotlib.pyplot as plt
    U = utils
    U.retro_style()
    RUST, BLUE, GRAY, INK = U.RUST, U.BLUE, U.GRAY, U.INK
    # One nested grid per row, so a wide tick label is charged to its own row only; e is a full-width band.
    fig, axd = plt.subplot_mosaic([[[['a', 'b']]], [[['c', 'd']]], ['e'], [[['f', 'g']]]],
                                  figsize=U.figsize_nature('double', aspect=1.0), layout='constrained',
                                  gridspec_kw=dict(height_ratios=[1.12, 1.16, 0.74, 1.10]))
    for key, ratios in (('a', [0.72, 1.5]), ('c', [1.42, 1.0]), ('f', [1.0, 1.02])):
        axd[key].get_subplotspec().get_gridspec().set_width_ratios(ratios)
    fig.get_layout_engine().set(h_pad=0.09, w_pad=0.06)

    # a: published call against the ClockBase call for the same samples
    ax = axd['a']
    p = data['published']
    U.confusion_matrix_panel(ax, np.array(p['matrix'], dtype=float), ['Higher', 'None', 'Lower'], normalize=False,
                             aspect='auto', show_zeros=True)
    ax.set_xlabel(f"ClockBase call\nagreement {100 * p['agreement']:.0f}%, " + rf"$\kappa$ = {p['kappa']:.2f}")
    ax.set_ylabel('Published call')
    ax.set_title(f"Published results ($n$ = {p['n']})")

    # b: three-way call accuracy against language models given the paper's description only
    ax = axd['b']
    rows = sorted(data['systems']['rows'], key=lambda r: -r['value'])
    color_b = lambda r: BLUE['400'] if r['label'] == 'ClockBase' else GRAY['400'] if r.get('open') else GRAY['300']
    ys_b, values_b = forest(ax, rows, color_b, f"Balanced accuracy of the call ($n$ = {data['systems']['records']} published results)",
                            None, U, percent=False)
    ax.set_title('ClockBase and language models alone')

    # c: share of comparisons with an error, by stage
    ax = axd['c']
    ys_c, values_c = forest(ax, data['stages'], lambda r: RUST['400'], 'Comparisons with an error (%)', None, U)
    ax.plot([], [], marker='o', ms=4.6, mfc=RUST['400'], mec=INK, mew=0.4, ls='none', label='Measured')
    ax.plot([], [], marker='o', ms=4.6, mfc='white', mec=RUST['400'], mew=1.0, ls='none', label='With checks (in-sample)')
    ax.set_title('Error rate by stage')

    # d: kinds of error in the released human DNAm comparisons
    ax = axd['d']
    kinds = data['kinds']
    yk = np.arange(len(kinds))[::-1]
    bars = ax.barh(yk, [k['units'] for k in kinds], height=0.62,
                   color=[RUST['200'] if k['design'] else RUST['400'] for k in kinds],
                   edgecolor=RUST['400'], linewidth=0.5)
    ax.set_yticks(yk)
    ax.set_yticklabels([k['label'] for k in kinds])
    ax.set_xlabel('Released human DNAm comparisons')
    ax.set_title('Kinds of error')

    # e: the same comparisons by series
    ax = axd['e']
    series = data['series']['rows']
    x = np.arange(len(series))
    design = np.array([s['design_only'] for s in series])
    other = np.array([s['other'] for s in series])
    correct = np.array([s['correct'] for s in series])
    ax.bar(x, design, width=0.68, color=RUST['200'], edgecolor=RUST['400'], linewidth=0.5,
           label='Pairing, repeats or age')
    ax.bar(x, other, width=0.68, bottom=design, color=RUST['400'], edgecolor=RUST['400'], linewidth=0.5,
           label='Other error')
    ax.bar(x, correct, width=0.68, bottom=design + other, color=BLUE['100'], edgecolor=BLUE['400'], linewidth=0.5,
           label='No error')
    ax.set_xticks(x)
    ax.set_xticklabels([s['gse'] for s in series])
    ax.set_ylabel('Released comparisons')
    ax.set_title(f"Human DNAm benchmark by series ($n$ = {data['series']['units']} comparisons)")

    # f: released atlas rows flagged by the checks on recorded fields, by outcome of the recomputation
    ax = axd['f']
    atlas = data['atlas']
    ya = np.arange(len(atlas))[::-1]
    left = np.zeros(len(atlas))
    for key, label, fill, edge in (('unchanged', 'Unchanged', BLUE['100'], BLUE['400']),
                                   ('changed', 'Changed', RUST['400'], RUST['400']),
                                   ('not_testable', 'Not testable', GRAY['200'], GRAY['400'])):
        share = np.array([100 * a[key] / a['released'] for a in atlas])
        ax.barh(ya, share, left=left, height=0.56, color=fill, edgecolor=edge, linewidth=0.5, label=label)
        left += share
    ax.set_yticks(ya)
    ax.set_yticklabels([a['label'] for a in atlas])
    ax.set_xlabel('Released rows flagged (%)')
    ax.set_title('Released atlas, field checks')

    # g: the checks on series outside the benchmark
    ax = axd['g']
    ys_g, values_g = forest(ax, data['validation'], lambda r: RUST['400'] if r.get('missed') else BLUE['400'],
                            'Comparisons (%)', None, U)
    ax.set_title('Checks on other series')

    U.trim_axes(fig, exclude=[axd['a'], axd['b'], axd['c'], axd['g']])
    finish_forest(axd['b'], ys_b, values_b, 'Accuracy', U, [0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
    finish_forest(axd['c'], ys_c, values_c, 'Errors / total', U, [0, 25, 50, 75, 100])
    finish_forest(axd['g'], ys_g, values_g, 'Count / total', U, [0, 25, 50, 75, 100])
    axd['f'].set_xticks([0, 5, 10, 15, 20])
    U.set_domain(axd['f'], xlim=(0, 24))  # room for the legend beside the two shorter bars
    axd['e'].set_yticks([0, 2, 4, 6, 8, 10])  # counts
    U.set_domain(axd['e'], ylim=axd['e'].get_ylim())
    U.nature_legend(axd['c'], title='Estimate', loc='upper right')
    U.nature_legend(axd['e'], title='Released comparison', loc='upper right')
    U.nature_legend(axd['f'], title='Recomputed call', loc='lower right')
    plt.setp(axd['e'].get_xticklabels(), rotation=40, ha='right', rotation_mode='anchor')
    U.add_bar_value_labels(axd['d'], bars, fmt='{:.0f}')
    U.add_row_value_column(axd['f'], ya, [f"{a['flagged']:,} / {a['released']:,}" for a in atlas],
                           header='Flagged / released', x=1.02, fontsize=U.TYPE['tick'])
    s = data['series']
    axd['e'].text(0.30, 0.95, f"Four series hold {s['top4_errors']} of the {s['errors']} errors",
                  transform=axd['e'].transAxes, ha='left', va='top', fontsize=U.TYPE['tick'], color=INK)
    U.add_panel_labels(fig, axd)
    if draft:  # a look at the layout while the gates still fail
        report = U.figure_preflight(fig, axd, verbose=False)
        print(report['text'])
        U.save_for_review(fig, str(out / 'Figure_benchmark_2026-10-07_review.png'), dpi=110)
        plt.close(fig)
        return
    U.save_figure(fig, str(out / 'Figure_benchmark_2026-10-07.pdf'), axes_dict=axd,
                  footer='Data: ClockBase verifier run final-2026-10-05 and benchmark results of 2026-10-06 and 2026-10-07; '
                         'computed in figure_benchmark_2026_10_07.py',
                  review_png=str(out / 'Figure_benchmark_2026-10-07_review.png'))
    fig.savefig(out / 'Figure_benchmark_2026-10-07.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--bench', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--draft', action='store_true', help='print the gate report and write only the review PNG')
    parser.add_argument('--table', action='store_true', help='write the summary table and stop')
    parser.add_argument('--scorecard', action='store_true', help='draw the stage figure and stop')
    parser.add_argument('--modality-table', action='store_true', help='write the benchmark table by modality and stop')
    parser.add_argument('--tidy', action='store_true', help='write the benchmark as one tidy CSV and stop')
    args = parser.parse_args()
    if args.modality_table:
        modality_table(args.bench, args.out)
        return
    if SCORED and not (args.table or args.scorecard or args.tidy):
        raise SystemExit('With CLOCKBASE_BENCH_SCORED, use --table, --modality-table or --scorecard: the seven-panel figure '
                         'describes the design checks before they entered the verifier.')
    data = collect(args.bench, args.out)
    (args.out / 'figure_benchmark_data.json').write_text(json.dumps(data, indent=1) + '\n')
    if args.table:
        table(data, args.bench, args.out)
        return
    if args.tidy:
        tidy(data, args.bench, args.out)
        return
    if args.scorecard:
        sys.path.insert(0, str(args.out))
        import nature_figure_utils
        scorecard(data, args.bench, args.out, nature_figure_utils)
        return
    sys.path.insert(0, str(args.out))
    import nature_figure_utils
    draw(data, args.out, nature_figure_utils, args.draft)


if __name__ == '__main__':
    main()
