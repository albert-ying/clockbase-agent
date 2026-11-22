"""Planted-defect self-test of the verifier. SIMULATED: defects are planted into copies of real atlas rows.

Data-generating process. From a finished verifier run (REFERENCE), take rows that passed every check with no
repair (VERIFIED, membership recovered or replay-linked, label orientation supported by the heuristic without
model adjudication), one row per study, a seeded random sample per modality. Each sampled row gets exactly one
treatment, in rotation:

  clean                copied unchanged                                  expect VERIFIED, statistics as REFERENCE
  p_altered            stored P x 1.37 (FDR raised to at least P)        expect MEMBERSHIP_UNRECOVERED, PENDING_DESIGN
  p_zero               stored P set to 0                                 expect STORED_P_UNUSABLE, PENDING_DESIGN
  fdr_below_p          stored FDR set to P / 10                          expect FDR_BELOW_RAW_P, FDR recomputed (>= P)
  effect_scaled        stored effect x 1.5                               expect EFFECT_NOT_REPRODUCED, effect replaced
                                                                         by the REFERENCE value
  effect_sign_flipped  stored effect, direction and log2FC negated       expect LABEL_ORIENTATION_INVERTED,
                       (the stored sign now contradicts the labels)      PENDING_ADJUDICATION
  wrong_values         DNAm: clock name replaced by another clock;       expect PENDING_DESIGN or PENDING_INPUT
                       RNA (one clock): study accession replaced by an   (membership or predictions not recovered)
                       unsampled eligible study (values from elsewhere)
  duplicate            the row saved twice                               expect one VERIFIED, one DUPLICATE

The sandbox atlas tables are written to ROOT/corrections and every deterministic stage (membership recovery,
replay links, merge, verdicts, repairs) runs on them through CLOCKBASE_CORRECTIONS, reusing the REFERENCE run's
replay traces and DNAm integrity screen. RNA recovery can read saved replay predictions when a clock_result
file is absent from the deposit. Repairs use the final chain's design estimand, classification, individuals
and checked clock-output decisions. No model is called. BH is recomputed over the simulated inferential
family, not the full reference family; clean equality concerns the unadjusted statistics listed below.

A second simulation tests the DNAm input-integrity screen: in copies of real DNAm matrices, a seeded third of the
samples have their beta values permuted across probes and a third are inverted (1 - beta); the screen should put
both below r = 0.10 and leave untouched samples at their REFERENCE r.

  python -m audit_system.verifier_selftest_2026_10_01 --root SANDBOX --reference VERIFIER_ROOT/FINAL [--per-modality 80]
"""
import argparse
import gzip
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import paths

REPO = Path(__file__).resolve().parents[1]
SOURCE = paths.CORRECTIONS  # the three atlas correction tables (CLOCKBASE_CORRECTIONS)
TABLES = {'dnam': 'DNAm', 'hs_rna': 'human-RNA', 'mm_rna': 'mouse-RNA'}
DEFECTS = ('clean', 'p_altered', 'p_zero', 'fdr_below_p', 'effect_scaled', 'effect_sign_flipped', 'wrong_values',
           'duplicate')
SEED = 20261002
P_FACTOR, EFFECT_FACTOR = 1.37, 1.5
SIMPLE_SOURCES = {'recovery', 'replay', 'replay+recovery'}
SCHEMA = 'clockbase-verifier-selftest-v2'


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def jsonl(path):
    with open(path) as handle:
        return [json.loads(line) for line in handle]


def table_path(folder, modality):
    return Path(folder) / f'2026-09-08-{TABLES[modality]}-atlas_WORKING_NOT_VALIDATED.csv'


def eligible(reference, modality):
    """REFERENCE rows that passed every check without a repair: the clean material defects are planted into."""
    verdicts = {r['row_ordinal']: r for r in jsonl(Path(reference) / f'verdicts-{modality}/row_verdicts.jsonl')}
    out = []
    for r in jsonl(Path(reference) / f'repairs-{modality}/repairs.jsonl'):
        v = verdicts[r['row_ordinal']]
        rec = r.get('recomputed') or {}
        if (r['final_class'] == 'VERIFIED' and r.get('membership_source') in SIMPLE_SOURCES
                and [a['action'] for a in r['actions']] == ['RECOMPUTE_DECLARED_ESTIMAND']
                and v['verdict'] == 'VERIFIED' and v['evidence'].get('label_orientation') == 'SUPPORTED'
                and rec.get('p') is not None):
            out.append(r)
    return out


def plant(row, defect, other_clocks, other_study):
    """Return the planted copies of one atlas row (a dict of the table's columns)."""
    out = dict(row)
    p, fdr = float(row['pvalue']), pd.to_numeric(row.get('fdr'), errors='coerce')
    if defect == 'p_altered':
        out['pvalue'] = p * P_FACTOR if p * P_FACTOR < 1 else p / P_FACTOR
        if np.isfinite(fdr) and fdr < out['pvalue']:
            out['fdr'] = out['pvalue']
    elif defect == 'p_zero':
        out['pvalue'] = 0.0
    elif defect == 'fdr_below_p':
        out['fdr'] = p / 10
    elif defect == 'effect_scaled':
        out['effect_size'] = float(row['effect_size']) * EFFECT_FACTOR
    elif defect == 'effect_sign_flipped':
        for column in ('effect_size', 'direction', 'log2FoldChange'):
            value = pd.to_numeric(row.get(column), errors='coerce')
            if np.isfinite(value):
                out[column] = -value
    elif defect == 'wrong_values':
        if other_clocks:
            out['clock'] = other_clocks[0]
        else:
            out['gse_id'] = other_study
    return [out, dict(out)] if defect == 'duplicate' else [out]


def build(root, reference, per_modality, rng):
    """Write the sandbox atlas tables and the plant key (sandbox ordinal -> base row and defect)."""
    corrections = Path(root) / 'corrections'
    corrections.mkdir(parents=True)
    key = []
    for modality in TABLES:
        atlas = pd.read_csv(table_path(SOURCE, modality), low_memory=False)
        clocks = sorted(atlas.clock.dropna().unique())
        pool = eligible(reference, modality)
        by_study = {}
        for r in pool:
            by_study.setdefault(r['gse_id'], []).append(r)
        studies = sorted(by_study)
        rng.shuffle(studies)
        chosen = [by_study[g][rng.integers(len(by_study[g]))] for g in studies[:per_modality]]
        rows = []
        for i, base in enumerate(chosen):
            defect = DEFECTS[i % len(DEFECTS)]
            source_row = atlas.iloc[base['row_ordinal']].to_dict()
            if source_row['gse_id'] != base['gse_id'] or source_row['clock'] != base['clock']:
                raise SystemExit(f'{modality} row {base["row_ordinal"]}: REFERENCE and source table disagree')
            others = [c for c in clocks if c != base['clock']]
            for copy, planted in enumerate(plant(source_row, defect, others, studies[per_modality + i])):
                key.append({'modality': modality, 'sandbox_ordinal': len(rows), 'base_ordinal': base['row_ordinal'],
                            'gse_id': base['gse_id'], 'defect': defect, 'copy': copy,
                            'reference': {k: base['recomputed'].get(k) for k in ('n_case', 'n_control', 'mean_diff', 'cohen_d', 'p')},
                            'planted': {k: planted.get(k) for k in ('gse_id', 'clock', 'pvalue', 'effect_size', 'fdr')}})
                rows.append(planted)
        pd.DataFrame(rows, columns=atlas.columns).to_csv(table_path(corrections, modality), index=False)
        print(f'{modality}: {len(pool)} eligible rows in {len(by_study)} studies; planted {len(rows)} rows from {len(chosen)}')
    return corrections, key


def repair_options(verifier_root, modality):
    """The deterministic repair settings of verifier_v3_final_chain_2026_09_30.sh."""
    verifier_root = Path(verifier_root)
    decisions = os.environ.get('CLOCKBASE_CLOCK_OUTPUT_DECISIONS') or str(
        REPO / 'audit_system/policy/clock_output_arm_decisions_2026-11-12.json')
    if not Path(decisions).is_file():
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {decisions}')
    options = ['--design-estimand', '--clock-output-decisions', decisions]
    for name, flag in (('classification-001', '--classification'), ('individuals-audit-001', '--individuals')):
        if (verifier_root / name / 'comparisons.jsonl').is_file():
            options += [flag, str(verifier_root / name)]
    if modality == 'dnam':
        options += ['--repaired-dnam', str(paths.STREAM_EMPTY)]
        if (paths.STREAM_DETACHED / 'studies').is_dir():
            options += ['--repaired-dnam', str(paths.STREAM_DETACHED)]
    return options


def run_stages(root, verifier_root, corrections, workers):
    env = {**os.environ, 'CLOCKBASE_CORRECTIONS': str(corrections)}
    py = [sys.executable, '-B', '-m']
    log = open(Path(root) / 'stages.log', 'a')
    for m in TABLES:
        replay = f'{verifier_root}/replay-{m}-001'
        verdicts = py + ['audit_system.row_verdicts_2026_09_27', '--modality', m, '--membership', f'{root}/merged-{m}',
                         '--out', f'{root}/verdicts-{m}']
        if m == 'dnam':
            verdicts += ['--integrity', f'{verifier_root}/input-integrity-dnam-001']
        for command in (
                py + ['audit_system.membership_recovery_2026_09_27', 'run', '--modality', m, '--out', f'{root}/membership-{m}',
                      '--workers', str(workers), '--prediction-replay', replay],
                py + ['audit_system.replay_link_2026_09_27', '--modality', m, '--replay', replay, '--out', f'{root}/link-{m}'],
                py + ['audit_system.membership_merge_2026_09_27', '--modality', m, '--recovery', f'{root}/membership-{m}',
                      '--link', f'{root}/link-{m}', '--replay', replay, '--out', f'{root}/merged-{m}'],
                verdicts,
                py + ['audit_system.row_repairs_2026_09_28', '--modality', m, '--membership', f'{root}/merged-{m}',
                      '--verdicts', f'{root}/verdicts-{m}', '--replay', replay, '--out', f'{root}/repairs-{m}',
                      *repair_options(verifier_root, m)]):
            log.write(f'$ {" ".join(command)}\n')
            log.flush()
            result = subprocess.run(command, cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise SystemExit(f'stage failed: {" ".join(command[3:5])}; see {root}/stages.log')
    log.close()


def same(a, b, tolerance=1e-9):
    return a is not None and b is not None and abs(a - b) <= tolerance * max(1.0, abs(a), abs(b))


def score(entry, verdict, repair, pair=None):
    """Outcome of one planted row: EXPECTED, FLAGGED_OTHER (stopped or flagged another way) or MISSED."""
    defect, reasons, cls = entry['defect'], set(verdict['reasons']), repair['final_class']
    rec, ref = repair.get('recomputed') or {}, entry['reference']
    stats_match = all(same(rec.get(k), ref[k]) for k in ('mean_diff', 'cohen_d', 'p')) and \
        (rec.get('n_case'), rec.get('n_control')) == (ref['n_case'], ref['n_control'])
    if defect == 'clean':
        return 'EXPECTED' if cls == 'VERIFIED' and stats_match else 'FALSE_ALARM'
    expected = {
        'p_altered': 'MEMBERSHIP_UNRECOVERED' in reasons and cls == 'PENDING_DESIGN',
        'p_zero': 'STORED_P_UNUSABLE' in reasons and cls == 'PENDING_DESIGN',
        'fdr_below_p': 'FDR_BELOW_RAW_P' in reasons and cls == 'VERIFIED' and stats_match
                       and rec.get('fdr_family', -1) >= rec.get('p', 2) * (1 - 1e-12),
        'effect_scaled': 'EFFECT_NOT_REPRODUCED' in reasons and cls == 'VERIFIED' and stats_match
                         and any(a['action'] == 'REPLACE_UNREPRODUCIBLE_EFFECT' for a in repair['actions']),
        'effect_sign_flipped': 'LABEL_ORIENTATION_INVERTED' in reasons and cls == 'PENDING_ADJUDICATION',
        'wrong_values': cls in ('PENDING_DESIGN', 'PENDING_INPUT'),
        'duplicate': pair is not None and sorted([cls, pair]) == ['DUPLICATE', 'VERIFIED'],
    }[defect]
    if expected:
        return 'EXPECTED'
    if defect == 'duplicate':
        return 'MISSED'
    if cls != 'VERIFIED' or verdict['verdict'] != 'VERIFIED':
        return 'FLAGGED_OTHER'  # stopped or flagged by another check
    # No flag: the row is published with recomputed statistics; they either undo the planted value or carry it.
    return 'CORRECTED_UNFLAGGED' if stats_match else 'MISSED'


def evaluate(root, key):
    verdicts, repairs = {}, {}
    for m in TABLES:
        verdicts[m] = {r['row_ordinal']: r for r in jsonl(Path(root) / f'verdicts-{m}/row_verdicts.jsonl')}
        repairs[m] = {r['row_ordinal']: r for r in jsonl(Path(root) / f'repairs-{m}/repairs.jsonl')}
    rows = []
    for entry in key:
        m, o = entry['modality'], entry['sandbox_ordinal']
        pair = None
        if entry['defect'] == 'duplicate':
            if entry['copy'] == 1:
                continue
            pair = repairs[m][o + 1]['final_class']
        outcome = score(entry, verdicts[m][o], repairs[m][o], pair)
        rows.append({**{k: entry[k] for k in ('modality', 'sandbox_ordinal', 'base_ordinal', 'gse_id', 'defect')},
                     'outcome': outcome, 'final_class': repairs[m][o]['final_class'],
                     'reasons': ';'.join(verdicts[m][o]['reasons'])})
    frame = pd.DataFrame(rows)
    table = frame.groupby('defect').outcome.value_counts().unstack(fill_value=0).reindex(DEFECTS)
    return frame, table


def detached_simulation(root, verifier_root, rng, studies=3, max_bytes=40_000_000):
    """Permute or invert a seeded third of the samples each in copies of real DNAm matrices and rescreen."""
    from audit_system.input_integrity_2026_09_30 import DNAM_INPUTS, REFERENCE, screen_dnam
    screened = pd.read_csv(Path(verifier_root) / 'input-integrity-dnam-001/samples.csv', usecols=['gse_id', 'gsm', 'r_reference'])
    clean = screened.groupby('gse_id').r_reference.agg(['min', 'size'])
    clean = clean[(clean['min'] >= 0.9) & (clean['size'] >= 6)].index
    candidates = sorted(g for g in clean if (DNAM_INPUTS / f'{g}.csv.gz').exists()
                        and (DNAM_INPUTS / f'{g}.csv.gz').stat().st_size < max_bytes)
    chosen = [candidates[i] for i in rng.choice(len(candidates), size=min(studies, len(candidates)), replace=False)]
    reference = pd.read_csv(REFERENCE, index_col=0, converters={0: str})['mean']
    out_dir = Path(root) / 'detached-simulation'
    out_dir.mkdir()
    results = []
    for gse in chosen:
        matrix = pd.read_csv(DNAM_INPUTS / f'{gse}.csv.gz', index_col=0, low_memory=False)
        order = rng.permutation(matrix.shape[1])
        third = len(order) // 3
        treatment = {matrix.columns[i]: 'permuted' for i in order[:third]}
        treatment.update({matrix.columns[i]: 'inverted' for i in order[third:2 * third]})
        values = matrix.apply(pd.to_numeric, errors='coerce')
        for column, kind in treatment.items():
            column_values = values[column].to_numpy()
            values[column] = rng.permutation(column_values) if kind == 'permuted' else 1 - column_values
        path = out_dir / f'{gse}_SIMULATED.csv.gz'
        with gzip.open(path, 'wt') as handle:
            values.to_csv(handle)
        frame = screen_dnam(path, reference)
        before = screened[screened.gse_id == gse].set_index('gsm').r_reference
        for gsm, r in zip(frame.gsm, frame.r_reference):
            kind = treatment.get(gsm, 'untouched')
            results.append({'gse_id': gse, 'gsm': gsm, 'treatment': kind, 'r_reference': r, 'r_reference_before': before.get(gsm),
                            'outcome': ('EXPECTED' if (r < 0.10 if kind != 'untouched' else same(r, before.get(gsm), 1e-6))
                                        else 'MISSED' if kind != 'untouched' else 'CHANGED')})
        # The simulated matrix is regenerable from the seed; retain screen results.
        subprocess.run(['/bin/rm', '-f', str(path)], check=True)
    frame = pd.DataFrame(results)
    frame.to_csv(out_dir / 'screen_SIMULATED.csv', index=False)
    return frame, chosen


def run(root, reference, per_modality=80, workers=8):
    root, reference = Path(root), Path(reference)
    verifier_root = reference.parent
    root.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(SEED)
    started = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    print('SIMULATED: planted defects in copies of real atlas rows (see module docstring)')
    corrections, key = build(root, reference, per_modality, rng)
    with open(root / 'plant_key_SIMULATED.jsonl', 'w') as sink:
        for entry in key:
            sink.write(json.dumps(entry, default=str) + '\n')
    run_stages(root, verifier_root, corrections, workers)
    rows, table = evaluate(root, key)
    rows.to_csv(root / 'outcomes_SIMULATED.csv', index=False)
    detached, detached_studies = detached_simulation(root, verifier_root, rng)
    planted = rows[rows.defect != 'clean']
    clean = rows[rows.defect == 'clean']
    summary = {
        'schema': SCHEMA, 'simulation': True,
        'label': 'SIMULATED: defects planted into copies of real atlas rows; not empirical error rates',
        'reference': str(reference), 'seed': SEED, 'per_modality': per_modality,
        'parameters': {'p_factor': P_FACTOR, 'effect_factor': EFFECT_FACTOR, 'fdr_below_p': 'P / 10',
                       'detached_threshold_r': 0.10},
        'source_tables_sha256': {m: sha256(table_path(SOURCE, m)) for m in TABLES},
        'sandbox_tables_sha256': {m: sha256(table_path(corrections, m)) for m in TABLES},
        'module_sha256': sha256(__file__),
        'verifier_modules_sha256': {name: sha256(REPO / 'audit_system' / name) for name in (
            'membership_recovery_2026_09_27.py', 'row_verdicts_2026_09_27.py', 'row_repairs_2026_09_28.py',
            'design_estimand_2026_10_10.py', 'clock_output_arms_2026_11_12.py', 'final_verification_rules_2026_11_13.py')},
        'repair_options': {m: repair_options(verifier_root, m) for m in TABLES},
        'clock_output_decisions_sha256': sha256(repair_options(verifier_root, 'dnam')[2]),
        'reference_repairs_sha256': {m: sha256(reference / f'repairs-{m}/repairs.jsonl') for m in TABLES},
        'outcomes_by_defect': {d: {k: int(v) for k, v in table.loc[d].items()} for d in DEFECTS if d in table.index},
        'planted_rows': len(planted),
        'caught_any_route': int((planted.outcome != 'MISSED').sum()),
        'caught_as_expected': int((planted.outcome == 'EXPECTED').sum()),
        'clean_rows': len(clean), 'clean_verified_unchanged': int((clean.outcome == 'EXPECTED').sum()),
        'detached_simulation': {'studies': detached_studies,
                                'by_treatment': {k: {o: int(n) for o, n in v.items()} for k, v in
                                                 detached.groupby('treatment').outcome.value_counts().unstack(fill_value=0)
                                                 .to_dict('index').items()}},
        'started': started, 'finished': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
    }
    (root / 'summary_SIMULATED.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps({k: summary[k] for k in ('outcomes_by_defect', 'planted_rows', 'caught_any_route', 'caught_as_expected',
                                              'clean_rows', 'clean_verified_unchanged', 'detached_simulation')}, indent=1))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True, help='new sandbox directory')
    parser.add_argument('--reference', type=Path, required=True, help='finished final chain directory')
    parser.add_argument('--per-modality', type=int, default=80)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args(argv)
    run(args.root, args.reference, args.per_modality, args.workers)


if __name__ == '__main__':
    main()
