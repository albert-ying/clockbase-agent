"""Validate the prototype design checks on series outside the group-choice benchmark, with two blinded model readers.

sample   Draws released units (distinct arm pairs, one per series) outside the benchmark series, by modality and
         stratum: SUBJECT (flagged for repeated subjects or pairing), AGE (flagged for age imbalance only) and NONE
         (no finding), plus every released benchmark unit (stratum BENCHMARK). Units whose arms hold more than
         MAX_SAMPLES samples are not drawn. Each package gives the readers the series text, the comparison names and
         the GEO fields of the samples in the two arms; no check output, no clock value and no harmonized field.
read     Two models answer each package: which samples come from the same individual, each sample's age, and whether
         age or time is the declared contrast.
score    Each reader's answer is turned into findings with the rule functions of design_checks; a finding stands as
         the reference when both readers give it, is absent when neither does, and is unresolved otherwise.
headline Packages for every comparison of the ranked tables (Table 1, Table 2, Figure 3, Table S14) instead of a
         sample; read and score then run on that directory.

  python -m benchmark.published_clock_2026_10_05.groupbench.design_validate sample --root G --checks DIR --out DIR
  python -m benchmark.published_clock_2026_10_05.groupbench.design_validate read --out DIR --max-requests N
  python -m benchmark.published_clock_2026_10_05.groupbench.design_validate score --root G --checks DIR --out DIR
  python -m benchmark.published_clock_2026_10_05.groupbench.design_validate headline --checks DIR --out DIR
"""
import argparse
import csv
import gzip
import hashlib
import json
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from audit_system.individuals_audit_2026_10_10 import PROMPT, package, series_text  # the verifier's reader prompt

from . import design_checks as dc
from .stats_a import clopper_pearson

SEED = 20261007
MAX_SAMPLES = 120
QUOTA = {'SUBJECT': 20, 'AGE': 15, 'NONE': 60}
SUBJECT_FINDINGS = {'REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT'}
FINDINGS = ('REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT', 'AGE_IMBALANCE')
YEAR = {'years': 1.0, 'months': 1 / 12, 'weeks': 1 / 52.18, 'days': 1 / 365.25}
def load_jsonl(path):
    return [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]


def benchmark_series(root):
    """Series of the group-choice benchmark (annotated for A1 or A2), by modality."""
    mouse = {p.stem for p in (root / 'a1/mapping').glob('*.csv')}
    with (root / 'a2/labels.csv').open() as handle:
        human = {row['dataset'] for row in csv.DictReader(handle)}
    return {'mm_rna': mouse, 'dnam': human, 'hs_rna': set()}


def stratum(findings):
    if findings & SUBJECT_FINDINGS:
        return 'SUBJECT'
    return 'AGE' if 'AGE_IMBALANCE' in findings else 'NONE'


def collect_units(checks, modality):
    """Released units of one modality with the union of their rows' findings and their arms from the export."""
    from audit_system import paths
    units = {}
    for row in load_jsonl(checks / modality / 'rows.jsonl'):
        unit = units.setdefault(row['unit'], {'unit': row['unit'], 'modality': modality, 'gse': row['gse'],
                                              'comparison_name': row['comparison_name'], 'condition_name': row['condition_name'],
                                              'row_ordinals': [], 'findings': set(), 'clocks': [], 'nonpace': False,
                                              'n': row['n_case'] + row['n_control']})
        unit['row_ordinals'].append(row['row_ordinal'])
        unit['findings'] |= set(row['findings'])
        unit['clocks'].append(row['clock'])
        unit['nonpace'] |= row['clock'] != 'DunedinPACE'
    from ..engine.predictions import final_export
    return units, final_export() / 'studies' / modality


def arms_of(export, unit):
    rows = json.loads((export / f"{unit['gse']}.json").read_text())['rows']
    row = next(r for r in rows if r['row_ordinal'] == unit['row_ordinals'][0])
    return [str(g) for g in row['case_gsms']], [str(g) for g in row['control_gsms']]


def sample(root, checks, out):
    from audit_system import paths
    from audit_system import membership_recovery_2026_09_27 as membership
    rng = random.Random(SEED)
    excluded = benchmark_series(root)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'prompts').mkdir(exist_ok=True)
    packages, frame_counts = [], {}
    for modality in ('dnam', 'hs_rna', 'mm_rna'):
        units, export = collect_units(checks, modality)
        folder = paths.GEO / membership.MODALITIES[modality]['folder']
        outside = [u for u in units.values() if u['gse'] not in excluded[modality]]
        benchmark = [u for u in units.values() if u['gse'] in excluded[modality]]
        chosen = []
        for name in ('SUBJECT', 'AGE', 'NONE'):
            pool = sorted((u for u in outside if stratum(u['findings']) == name), key=lambda u: u['unit'])
            eligible = [u for u in pool if u['n'] <= MAX_SAMPLES]
            frame_counts[f'{modality}:{name}'] = {'units': len(pool), 'units_within_size_limit': len(eligible),
                                                  'series': len({u['gse'] for u in eligible})}
            rng.shuffle(eligible)
            seen = set()
            for unit in eligible:
                if unit['gse'] in seen:
                    continue
                seen.add(unit['gse'])
                chosen.append({**unit, 'stratum': name})
                if len(seen) == QUOTA[name]:
                    break
        chosen += [{**u, 'stratum': 'BENCHMARK'} for u in sorted(benchmark, key=lambda u: u['unit']) if u['n'] <= MAX_SAMPLES]
        frame_counts[f'{modality}:BENCHMARK'] = {'units': len(benchmark),
                                                 'units_within_size_limit': sum(u['n'] <= MAX_SAMPLES for u in benchmark)}
        for unit in chosen:
            raw = dc.read_frame(folder / 'processed_metadata_preharmonize' / f"{unit['gse']}.csv")
            if raw is None:
                continue
            unit['case'], unit['control'] = arms_of(export, unit)
            prompt = package(unit, raw, series_text(folder, unit['gse']))
            digest = hashlib.sha256(prompt.encode()).hexdigest()[:16]
            (out / 'prompts' / f'{digest}.txt').write_text(prompt)
            packages.append({**{k: v for k, v in unit.items() if k != 'findings'}, 'prompt': digest})
    (out / 'packages.jsonl').write_text(''.join(json.dumps(p) + '\n' for p in packages))
    (out / 'frame.json').write_text(json.dumps({'seed': SEED, 'max_samples': MAX_SAMPLES, 'quota': QUOTA,
                                                'frame': frame_counts,
                                                'drawn': dict(Counter(f"{p['modality']}:{p['stratum']}" for p in packages))},
                                               indent=1) + '\n')
    print(json.dumps(json.loads((out / 'frame.json').read_text()), indent=1))


HEADLINE = {'table1_pro_aging': 'Table 1', 'table2_anti_aging': 'Table 2', 'figure3_category_top15': 'Figure 3',
            's14_age_lowering': 'Table S14'}
MODALITY_NAME = {'human DNAm': 'dnam', 'human RNA-seq': 'hs_rna', 'mouse RNA-seq': 'mm_rna'}


def headline(checks, out):
    """Packages for the comparisons of the ranked tables; a unit lists every table it appears in."""
    from audit_system import paths
    from audit_system import membership_recovery_2026_09_27 as membership
    from ..engine.predictions import final_export
    ranking = final_export() / 'ranking'
    rows = {}
    for modality in ('dnam', 'hs_rna', 'mm_rna'):
        for row in load_jsonl(checks / modality / 'rows.jsonl'):
            rows.setdefault((modality, row['gse'], row['comparison_name'], row['condition_name'], row['clock']), []).append(row)
    wanted, unmatched = {}, []
    for name, label in HEADLINE.items():
        with (ranking / f'{name}.csv').open() as handle:
            for entry in csv.DictReader(handle):
                key = (MODALITY_NAME[entry['modality']], entry['gse_id'], entry['comparison_name'], entry['condition_name'], entry['clock'])
                p = float(entry['v3_p'])
                # The ranked arm sizes count samples with a clock value; the arms can list one more sample.
                match = [r for r in rows.get(key, []) if abs(r['released']['p'] - p) <= 1e-12 + 1e-6 * abs(p)]
                if not match:
                    unmatched.append({'table': label, 'gse': entry['gse_id'], 'comparison_name': entry['comparison_name']})
                    continue
                wanted.setdefault((key[0], match[0]['unit']), {'tables': [], 'row': match[0]})['tables'].append(label)
    out.mkdir(parents=True, exist_ok=True)
    (out / 'prompts').mkdir(exist_ok=True)
    packages, too_large = [], []
    units_by_modality = {m: collect_units(checks, m) for m in ('dnam', 'hs_rna', 'mm_rna')}
    for (modality, key), entry in sorted(wanted.items()):
        units, export = units_by_modality[modality]
        unit = {**units[key], 'stratum': '+'.join(sorted(set(entry['tables']))), 'ranked_clock': entry['row']['clock']}
        unit['nonpace'] = entry['row']['clock'] != 'DunedinPACE'
        if unit['n'] > MAX_SAMPLES:
            too_large.append({'unit': key, 'gse': unit['gse'], 'comparison_name': unit['comparison_name'], 'n': unit['n'],
                              'tables': unit['stratum']})
            continue
        folder = paths.GEO / membership.MODALITIES[modality]['folder']
        raw = dc.read_frame(folder / 'processed_metadata_preharmonize' / f"{unit['gse']}.csv")
        if raw is None:
            unmatched.append({'table': unit['stratum'], 'gse': unit['gse'], 'comparison_name': unit['comparison_name'],
                              'reason': 'no GEO sample fields'})
            continue
        unit['case'], unit['control'] = arms_of(export, unit)
        prompt = package(unit, raw, series_text(folder, unit['gse']))
        digest = hashlib.sha256(prompt.encode()).hexdigest()[:16]
        (out / 'prompts' / f'{digest}.txt').write_text(prompt)
        packages.append({**{k: v for k, v in unit.items() if k != 'findings'}, 'prompt': digest})
    (out / 'packages.jsonl').write_text(''.join(json.dumps(p) + '\n' for p in packages))
    frame = {'max_samples': MAX_SAMPLES, 'units': len(wanted), 'packages': len(packages), 'not_matched': unmatched,
             'above_size_limit': too_large, 'by_table': dict(Counter(t for p in packages for t in p['stratum'].split('+')))}
    (out / 'frame.json').write_text(json.dumps(frame, indent=1) + '\n')
    print(json.dumps({k: v for k, v in frame.items() if k != 'above_size_limit'}, indent=1), 'above size limit', len(too_large))


def read(out, max_requests, workers, limit=None):
    from ..engine.groups import Adjudicator
    adjudicator = Adjudicator(out / 'ledger.jsonl', max_requests=max_requests, cache_only=max_requests == 0)
    packages = load_jsonl(out / 'packages.jsonl')[:limit]

    def one(item):
        responses = adjudicator.ask((out / 'prompts' / f"{item['prompt']}.txt").read_text())
        return {'unit': item['unit'], 'answers': [{'model': r.get('model_requested'), 'status': r.get('status'),
                                                   'parsed': r.get('parsed')} for r in responses]}

    with ThreadPoolExecutor(max_workers=workers) as pool:
        answers = list(pool.map(one, packages))
    (out / 'reads.jsonl').write_text(''.join(json.dumps(a) + '\n' for a in answers))
    print(dict(Counter(a['status'] for r in answers for a in r['answers'])), 'models', adjudicator.models)


def reader_findings(item, answer):
    """Findings implied by one reader's answer, with the rule functions of design_checks; None when unusable."""
    if not isinstance(answer, dict) or not isinstance(answer.get('individuals'), dict):
        return None
    case, control = item['case'], item['control']
    subjects = {g: str(v) for g, v in answer['individuals'].items() if v is not None and dc.clean(v) is not None}
    result = dc.subject_findings(case, control, subjects)
    shared = answer.get('shared_individuals')
    found = set(result.get('findings', [])) if result.get('assessable') else set()
    cell_line = answer.get('material') == 'CELL_LINE'  # an established cell line is not an individual
    if shared == 'NO' or cell_line:
        found = set()
    out = {'shared': shared, 'material': answer.get('material'),
           'subject_determined': shared == 'NO' or cell_line or bool(result.get('assessable')),
           'REPEATED_SUBJECT': 'REPEATED_SUBJECT' in found, 'PAIRED_AS_INDEPENDENT': 'PAIRED_AS_INDEPENDENT' in found}
    scale = YEAR.get(answer.get('age_unit')) if item['modality'] in dc.HUMAN else 1.0
    ages = {}
    if scale and isinstance(answer.get('ages'), dict):
        for gsm, value in answer['ages'].items():
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                ages[gsm] = float(value) * scale
    age = dc.age_finding(case, control, ages, item['modality'] in dc.HUMAN) if ages else None
    out['age_recorded'] = age is not None
    out['age_gap'] = age['median_gap'] if age else None
    out['age_declared'] = bool(answer.get('age_or_time_is_the_contrast'))
    out['AGE_IMBALANCE'] = bool(age and age['criterion'] and not out['age_declared'] and item['nonpace'])
    return out


def reference(readers, finding):
    """True or False when both readers agree; None when one is unusable or they differ."""
    if any(r is None for r in readers):
        return None
    if finding != 'AGE_IMBALANCE' and not all(r['subject_determined'] for r in readers):
        return None
    values = {r[finding] for r in readers}
    return values.pop() if len(values) == 1 else None


def rate(k, n):
    low, high = clopper_pearson(k, n)
    return {'count': k, 'of': n, 'rate': k / n if n else None, 'interval': [low, high]}


def reader_recomputation(item, answers, values):
    """The statistic on subject means for each reader that reports shared individuals in the arms."""
    out = []
    for answer in answers:
        parsed = answer.get('parsed') if answer.get('status') == 'OK' else None
        found = reader_findings(item, parsed)
        if not found or not (found['REPEATED_SUBJECT'] or found['PAIRED_AS_INDEPENDENT']):
            out.append(None)
            continue
        subjects = {g: str(v) for g, v in parsed['individuals'].items() if v is not None and dc.clean(v) is not None}
        findings = [f for f in ('REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT') if found[f]]
        result = dc.recompute(item['case'], item['control'], values, subjects, {}, findings)
        result['call'] = dc.call(result.get('estimate'), result.get('p')) if result.get('status') == 'OK' else None
        result['subjects'] = [len({subjects[g] for g in arm if g in subjects}) for arm in (item['case'], item['control'])]
        out.append(result)
    return out


def score(root, checks, out, suffix=''):
    packages = {p['unit']: p for p in load_jsonl(out / 'packages.jsonl')}
    reads = {r['unit']: r for r in load_jsonl(out / 'reads.jsonl')}
    found = {}
    for modality in ('dnam', 'hs_rna', 'mm_rna'):
        for row in load_jsonl(checks / modality / 'rows.jsonl'):
            found.setdefault(row['unit'], set()).update(row['findings'])
    units = []
    for unit, item in packages.items():
        answers = reads.get(unit, {}).get('answers', [])
        readers = [reader_findings(item, a.get('parsed')) if a.get('status') == 'OK' else None for a in answers]
        record = {'unit': unit, 'modality': item['modality'], 'gse': item['gse'], 'stratum': item['stratum'],
                  'comparison_name': item['comparison_name'], 'check': sorted(found.get(unit, set())), 'readers': readers}
        for finding in FINDINGS:
            record[f'reference_{finding}'] = reference(readers, finding) if len(readers) == 2 else None
        units.append(record)
    (out / f'scored{suffix}.jsonl').write_text(''.join(json.dumps(u) + '\n' for u in units))
    if any(label in u['stratum'] for u in units for label in HEADLINE.values()):
        from ..engine.predictions import PredictionReader, atlas_rows
        reader, released = PredictionReader(), {}
        for modality in ('dnam', 'hs_rna', 'mm_rna'):
            for row in load_jsonl(checks / modality / 'rows.jsonl'):
                released[(row['unit'], row['clock'])] = row['released']
        for u in units:
            item = packages[u['unit']]
            if not any(r and (r['REPEATED_SUBJECT'] or r['PAIRED_AS_INDEPENDENT']) for r in u['readers']):
                continue
            loaded = reader.load(item['modality'], item['gse'], item['ranked_clock'], atlas_rows(item['modality'], item['gse']))
            values = ({str(k): float(v) for k, v in loaded['samples']['clock'].dropna().items()}
                      if loaded.get('status') == 'OK' else {})
            u['released'] = released.get((u['unit'], item['ranked_clock']))
            u['arms'] = [len(item['case']), len(item['control'])]
            u['reader_recomputation'] = reader_recomputation(item, reads[u['unit']]['answers'], values) if values else None
        summary = {}
        for label in HEADLINE.values():
            members = [u for u in units if label in u['stratum'].split('+')]
            block = {'units_read': len(members),
                     'both_readers_determined_individuals': sum(u['reference_REPEATED_SUBJECT'] is not None for u in members),
                     'both_readers_recorded_age': sum(len(u['readers']) == 2 and all(r and r['age_recorded'] for r in u['readers'])
                                                      for u in members)}
            for finding in FINDINGS:
                block[finding] = {'both_readers': sum(u[f'reference_{finding}'] is True for u in members),
                                  'one_reader': sum(u[f'reference_{finding}'] is None and len(u['readers']) == 2
                                                    and any(r and r[finding] for r in u['readers']) for u in members),
                                  'recorded_field_check': sum(finding in u['check'] for u in members)}
            block['units_with_a_finding'] = [
                {'gse': u['gse'], 'comparison_name': u['comparison_name'], 'check': u['check'], 'arms': u.get('arms'),
                 'released': u.get('released'), 'reader_recomputation': u.get('reader_recomputation'),
                 'both_readers': [f for f in FINDINGS if u[f'reference_{f}'] is True],
                 'one_reader': [f for f in FINDINGS if u[f'reference_{f}'] is None and len(u['readers']) == 2
                                and any(r and r[f] for r in u['readers'])]}
                for u in members if u['check'] or any(u[f'reference_{f}'] is True or (
                    u[f'reference_{f}'] is None and len(u['readers']) == 2 and any(r and r[f] for r in u['readers']))
                    for f in FINDINGS)]
            summary[label] = block
        (out / f'headline{suffix}.json').write_text(json.dumps(summary, indent=1) + '\n')
        print(json.dumps(summary, indent=1))
        return
    summary = {'reader_agreement': {}, 'outside': {}, 'benchmark': {}}
    usable = [u for u in units if len(u['readers']) == 2 and all(u['readers'])]
    for finding in FINDINGS:
        same = sum(u['readers'][0][finding] == u['readers'][1][finding] for u in usable)
        summary['reader_agreement'][finding] = rate(same, len(usable))
    for scope, members in (('outside', [u for u in units if u['stratum'] != 'BENCHMARK']),
                           ('benchmark', [u for u in units if u['stratum'] == 'BENCHMARK'])):
        for modality in sorted({u['modality'] for u in members}):
            block = summary[scope].setdefault(modality, {})
            for finding in FINDINGS:
                resolved = [u for u in members if u['modality'] == modality and u[f'reference_{finding}'] is not None]
                flagged = [u for u in resolved if finding in u['check']]
                strata = ('SUBJECT',) if finding != 'AGE_IMBALANCE' else ('SUBJECT', 'AGE')
                unflagged = [u for u in resolved if finding not in u['check']
                             and (scope == 'benchmark' or u['stratum'] == 'NONE')]
                block[finding] = {
                    'units_resolved': len(resolved),
                    'units_unresolved': sum(u['modality'] == modality and u[f'reference_{finding}'] is None for u in members),
                    'flagged_confirmed': rate(sum(u[f'reference_{finding}'] for u in flagged if scope == 'benchmark' or u['stratum'] in strata),
                                              sum(1 for u in flagged if scope == 'benchmark' or u['stratum'] in strata)),
                    'unflagged_with_finding': rate(sum(u[f'reference_{finding}'] for u in unflagged), len(unflagged))}
            none = [u for u in members if u['modality'] == modality and (scope == 'benchmark' or u['stratum'] == 'NONE')]
            resolved = [u for u in none if all(u[f'reference_{f}'] is not None for f in FINDINGS)]
            block['unflagged_any_missed'] = rate(sum(any(u[f'reference_{f}'] and f not in u['check'] for f in FINDINGS)
                                                     for u in resolved), len(resolved))
    (out / f'validation{suffix}.json').write_text(json.dumps(summary, indent=1) + '\n')
    print(json.dumps(summary, indent=1))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('command', choices=['sample', 'read', 'score', 'headline'])
    parser.add_argument('--root', type=Path)
    parser.add_argument('--checks', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--max-requests', type=int, default=0, help='0 replays cached answers only')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--limit', type=int, help='read only the first N packages')
    parser.add_argument('--suffix', default='', help='suffix of the score output files')
    args = parser.parse_args()
    if args.command == 'sample':
        sample(args.root, args.checks, args.out)
    elif args.command == 'headline':
        headline(args.checks, args.out)
    elif args.command == 'read':
        read(args.out, args.max_requests, args.workers, args.limit)
    else:
        score(args.root, args.checks, args.out, args.suffix)


if __name__ == '__main__':
    main()
