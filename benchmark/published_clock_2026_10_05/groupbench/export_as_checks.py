"""Score a verifier run that already applies the design-respecting estimand with the benchmark's existing functions.

The benchmark functions (design_eval, a_result --stages, the table and figure builders) read one record per released
row with its findings and a `candidate` statistic. For a run whose row repairs use the design estimand, the released
statistic is that candidate. This module writes those records from the run's public export, so that every benchmark
function can be pointed at the run without change:

  released    the row's statistic in the reference run of the benchmark (CLOCKBASE_FINAL_NAME, default
              final-2026-10-05); kept so that both releases stay side by side
  candidate   the row's statistic in the scored run; rows the scored run no longer releases have no candidate call
  findings    what the scored run corrected, in the benchmark's names: repeated individuals and paired arms it
              analysed by individual (recorded in GEO or found by the two-reader audit) and age it entered as a
              covariate. A soft flag without a corrected statistic (age separating the arms) is listed under
              `flags` and does not count as a finding

  python -m benchmark.published_clock_2026_10_05.groupbench.export_as_checks --scored final-2026-10-10 --out DIR
"""
import argparse
import json
from collections import Counter
from pathlib import Path

from audit_system import paths

from ..engine.predictions import FINAL_NAME
from .design_checks import call, formula_has_age, unit_key

RELEASED = ('VERIFIED', 'PROVISIONAL')
FINDING = {'REPEATED_INDIVIDUAL': 'REPEATED_SUBJECT', 'PAIRED_ARMS': 'PAIRED_AS_INDEPENDENT'}
AGE_METHODS = {'AGE_ADJUSTED', 'AGE_ADJUSTED_UNIT_MEANS'}
SOFT = {'AGE_SEPARATES_ARMS'}


def findings_of(design):
    """The benchmark's names for what a production design block corrected (empty when the Welch result stands)."""
    design = design or {}
    if not design.get('method'):
        return []
    out = [FINDING[f] for f in design.get('findings', []) if f in FINDING]
    if design['method'] in AGE_METHODS:
        out.append('AGE_IMBALANCE')
    return out


def statistic(row):
    recomputed = row.get('recomputed') or {}
    return {'estimate': recomputed.get('mean_diff'), 'p': recomputed.get('p'), 'fdr': recomputed.get('fdr_family'),
            'call': call(recomputed.get('mean_diff'), recomputed.get('p'))}


def records(reference_rows, scored_rows, modality, gse):
    """One record per row released in the reference run."""
    scored = {r['row_ordinal']: r for r in scored_rows}
    for row in reference_rows:
        if row.get('final_class') not in RELEASED or not row.get('case_gsms') or not row.get('control_gsms'):
            continue
        new = scored.get(row['row_ordinal']) or {}
        test = (row.get('membership') or {}).get('test')
        formula = ((row.get('replay') or {}).get('test') or {}).get('formula')
        record = {'modality': modality, 'gse': gse, 'row_ordinal': row['row_ordinal'], 'clock': row['clock'],
                  'comparison_name': row.get('comparison_name'), 'condition_name': row.get('condition_name'),
                  'final_class': row['final_class'],
                  'unit': unit_key(gse, [str(g) for g in row['case_gsms']], [str(g) for g in row['control_gsms']]),
                  'n_case': len(row['case_gsms']), 'n_control': len(row['control_gsms']), 'agent_test': test,
                  'agent_paired': test in ('ttest_rel', 'sm.MixedLM'), 'agent_age_adjusted': formula_has_age(formula),
                  'released': statistic(row), 'findings': findings_of(new.get('design')),
                  'flags': sorted(set(new.get('soft_flags') or []) & SOFT), 'scored_final_class': new.get('final_class')}
        design = new.get('design') or {}
        if new.get('final_class') in RELEASED and (new.get('recomputed') or {}).get('p') is not None:
            stats = statistic(new)
            record['candidate'] = {'method': design.get('method') or 'Welch t-test (released)', 'status': 'OK',
                                   'estimate': stats['estimate'], 'p': stats['p'], 'call': stats['call'],
                                   'subject_source': design.get('subject_source')}
        else:
            record['candidate'] = {'method': design.get('method') or 'Welch t-test (released)',
                                   'status': 'NOT_RELEASED_' + str(new.get('final_class')), 'estimate': None, 'p': None, 'call': None}
        yield record


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--scored', required=True, help='name of the verifier run to score')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    reference = paths.VERIFIER_ROOT / FINAL_NAME / 'public-export/studies'
    scored = paths.VERIFIER_ROOT / args.scored / 'public-export/studies'
    summary = {'reference': FINAL_NAME, 'scored': args.scored}
    for modality in ('dnam', 'hs_rna', 'mm_rna'):
        target = args.out / modality
        target.mkdir(parents=True, exist_ok=True)
        counts = Counter()
        with (target / 'rows.jsonl').open('w') as sink:
            for path in sorted((reference / modality).glob('*.json')):
                other = scored / modality / path.name
                try:
                    reference_rows = json.loads(path.read_text())['rows']
                    scored_rows = json.loads(other.read_text())['rows'] if other.exists() else []
                except (OSError, ValueError) as error:
                    raise ValueError(f'Cannot read study export {path} or {other}: {error}') from error
                for record in records(reference_rows, scored_rows, modality, path.stem):
                    sink.write(json.dumps(record) + '\n')
                    counts['rows'] += 1
                    counts[record['candidate']['status']] += 1
                    counts['method ' + record['candidate']['method']] += 1
                    counts['call differs'] += record['candidate']['status'] == 'OK' and record['candidate']['call'] != record['released']['call']
        summary[modality] = dict(counts)
    (args.out / 'export_as_checks.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


if __name__ == '__main__':
    main()
