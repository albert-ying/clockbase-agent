"""Package a stratified random 40 units for blind review by two domain experts (PROTOCOL.md v2, adjudication).

Strata: A1 rule-error, A1 rule-correct, A2 rule-error, A2 rule-correct; 10 per stratum (all when fewer), seed 20261006.
The workbook shows what the model reviewers saw (names, arms, per-sample expert labels and GEO metadata, reported
results with arm means) and nothing else: no rule outcome, model judgment or verifier output. Reviewers fill the
Verdicts sheet independently.

  python -m benchmark.published_clock_2026_10_05.groupbench.human_review --root G --out FILE.xlsx
"""
import argparse
import json
import random
from pathlib import Path

from openpyxl import Workbook

SEED = 20261006
INSTRUCTIONS = [
    'Blind review of comparisons made by an AI agent (ClockBase). Each unit is a pair of sample groups (arm 1, arm 2)',
    'from one GEO series, compared on an aging clock. The Units sheet gives the comparison names (condition_name names',
    'the case level), which arm the agent treated as the case, and per clock the reported direction (+1 higher, -1',
    'lower in the agent case), reported P and the mean clock value of each arm. The Samples sheet lists every sample with',
    'its expert label (Tyshkovskiy et al. 2026 rodent annotation for A1; ComputAgeBench for A2) and GEO metadata.',
    'For each unit and name variant, record in Verdicts: correct, error or not assessable; for errors, the type',
    '(wrong factor, label mismatch, impure arm, confounded, age-confounded, repeated subject, paired analysed as',
    'independent, P value from an arm under 3 samples, reported result inverted relative to the named groups,',
    'degenerate test, other) and the evidence (GSMs and field values). Work independently; do not consult other',
    'reviewers or any other output of the benchmark.',
]


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    rules = [json.loads(line) for line in (args.root / 'rules-v2/e1_rules.jsonl').read_text().splitlines()]
    packages = {}
    for path in (args.root / 'census-v2/packages').glob('unit_*.json'):
        package = json.loads(path.read_text())
        packages[package['unit_id']] = package
    rng = random.Random(SEED)
    chosen = []
    for set_name in ('A1', 'A2'):
        for rule_class in ('ERROR', 'CORRECT'):
            pool = sorted(r['unit_id'] for r in rules if r['set'] == set_name and r['class'] == rule_class and r['unit_id'] in packages)
            chosen += rng.sample(pool, min(10, len(pool)))
    rng.shuffle(chosen)
    wb = Workbook()
    sheet = wb.active
    sheet.title = 'Instructions'
    for line in INSTRUCTIONS:
        sheet.append([line])
    units = wb.create_sheet('Units')
    units.append(['review_id', 'set', 'series', 'variant', 'comparison_name', 'condition_name', 'agent_case_arm', 'clock',
                  'test', 'reported_direction', 'reported_p', 'mean_arm_1', 'mean_arm_2'])
    samples = wb.create_sheet('Samples')
    samples.append(['review_id', 'arm', 'gsm', 'expert_label', 'unlabelled_category', 'subject', 'geo_metadata'])
    verdicts = wb.create_sheet('Verdicts')
    verdicts.append(['review_id', 'variant', 'reviewer', 'verdict', 'error_type', 'evidence'])
    key = {}
    for number, unit_id in enumerate(chosen, 1):
        review_id = f'R{number:02d}'
        key[review_id] = unit_id
        package = packages[unit_id]
        for variant in package['variants']:
            for reported in variant['reported']:
                units.append([review_id, package['set'], package['gse'], variant['variant'], variant['comparison_name'],
                              variant['condition_name'], variant['agent_case_arm'], reported['clock'], reported['test'],
                              reported['reported_direction'], reported['reported_p'], reported['mean_arm_1'], reported['mean_arm_2']])
            verdicts.append([review_id, variant['variant'], '', '', '', ''])
        for s in package['samples']:
            samples.append([review_id, s['arm'], s['gsm'], json.dumps(s['expert_label']) if s['expert_label'] else '',
                            s['unlabelled_category'] or '', s['subject'] or '', json.dumps(s['geo'])])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    wb.save(args.out)
    args.out.with_suffix('.key.json').write_text(json.dumps({'seed': SEED, 'key': key}, indent=1) + '\n')
    print(len(chosen), 'units packaged ->', args.out)


if __name__ == '__main__':
    main()
