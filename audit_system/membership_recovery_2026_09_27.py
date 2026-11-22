"""Recover historical case/control membership of atlas rows by reproducing their stored statistics.

The saved atlas rows record gse_id, comparison and condition labels, the clock, the test name,
P and Cohen's d, but no sample identifiers. For two-sample t-test rows, every candidate split of
the study's samples is enumerated from its GEO metadata fields and the GPT-4 per-sample
annotations (every pair of levels of one field, one level against the rest, and pairs of levels
within one level of a second field). Welch and Student t-tests are recomputed on the saved
per-sample clock predictions. A split whose P matches the stored P to relative 1e-6 recovers the
membership. The stored effect sign fixes which arm is the case arm; the effect definition that
reproduces the stored effect is recorded.

  python -m audit_system.membership_recovery_2026_09_27 run --modality dnam --out DIR [--limit N]

Outcomes per row: UNIQUE (one membership reproduces P), AMBIGUOUS (several distinct memberships
do), NONE (no candidate does), NO_PREDICTIONS, NOT_TWO_SAMPLE_T (other test families).
"""
import argparse
import hashlib
import itertools
import json
import re
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system import paths

GEO = paths.GEO
# CLOCKBASE_CORRECTIONS points every verifier stage at another set of atlas tables (the planted-defect self-test).
CORRECTIONS = paths.CORRECTIONS
MODALITIES = {
    'dnam': {'table': 'DNAm', 'folder': 'GEO_hs_dnam', 'gpt4': 'hs_dnam_meta.csv'},
    'hs_rna': {'table': 'human-RNA', 'folder': 'GEO_hs_rna', 'gpt4': 'hs_rna_meta.csv'},
    'mm_rna': {'table': 'mouse-RNA', 'folder': 'GEO_mm_rna', 'gpt4': 'mm_rna_meta.csv'},
}
P_RELATIVE_TOLERANCE = 1e-6
EFFECT_TOLERANCE = 1e-4
MAX_LEVELS, MAX_STRATA_LEVELS = 40, 12
SCHEMA = 'clockbase-membership-recovery-v1'


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def is_two_sample_t(test):
    test = str(test).lower()
    return bool(re.search(r't-test|t test|welch|student', test)) and not re.search(r'paired|one[- ]sample', test)


def candidate_splits(meta):
    """Yield (description, index_a, index_b) for every candidate two-group split."""
    fields = [c for c in meta.columns if 1 < meta[c].nunique(dropna=True) <= MAX_LEVELS]
    for field in fields:
        levels = meta[field].dropna().unique()
        for a, b in itertools.combinations(levels, 2):
            yield (field, a, b, None), meta.index[meta[field] == a], meta.index[meta[field] == b]
        labelled = meta.index[meta[field].notna()]
        for a in levels:
            case = meta.index[meta[field] == a]
            yield (field, a, '<rest>', None), case, labelled.difference(case)
        for stratum_field in fields:
            if stratum_field == field or meta[stratum_field].nunique() > MAX_STRATA_LEVELS:
                continue
            for stratum in meta[stratum_field].dropna().unique():
                sub = meta[meta[stratum_field] == stratum]
                for a, b in itertools.combinations(sub[field].dropna().unique(), 2):
                    yield (field, a, b, (stratum_field, stratum)), sub.index[sub[field] == a], sub.index[sub[field] == b]


def effect_sizes(a, b):
    diff = a.mean() - b.mean()
    pooled = np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) / (len(a) + len(b) - 2))
    average = np.sqrt((a.var() + b.var()) / 2)
    d = diff / pooled if pooled else np.nan
    return {'mean_diff': float(diff), 'cohen_d_pooled': float(d),
            'cohen_d_avg_sd': float(diff / average) if average else np.nan,
            'hedges_g': float(d * (1 - 3 / (4 * (len(a) + len(b)) - 9))),
            'glass_delta': float(diff / b.std()) if b.std() else np.nan}


def unique_samples(values):
    """Drop every sample identifier that occurs more than once: its prediction is not attributable."""
    duplicated = values.index[values.index.duplicated(keep=False)]
    return values[~values.index.isin(duplicated)], sorted(set(map(str, duplicated)))


def recover_row(values, meta, stored_p, stored_effect):
    """Match one stored (P, effect) against every candidate split of one prediction vector."""
    values, duplicated = unique_samples(values)
    meta = meta[~meta.index.duplicated(keep=False)]
    matches, seen = {}, set()
    for key, case, control in candidate_splits(meta):
        a, b = values.reindex(case).dropna(), values.reindex(control).dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        identity = frozenset([tuple(sorted(a.index)), tuple(sorted(b.index))])
        if identity in seen:
            continue
        seen.add(identity)
        for equal_var in (False, True):
            p = stats.ttest_ind(a, b, equal_var=equal_var).pvalue
            if not (np.isfinite(p) and abs(p - stored_p) <= P_RELATIVE_TOLERANCE * stored_p):
                continue
            if np.isfinite(stored_effect) and np.sign(a.mean() - b.mean()) != np.sign(stored_effect):
                a, b, key = b, a, (key[0], key[2], key[1], key[3])
            effects = effect_sizes(a, b)
            matched = [k for k, v in effects.items()
                       if np.isfinite(v) and abs(v - stored_effect) <= EFFECT_TOLERANCE * max(1.0, abs(stored_effect))]
            matches[identity] = {
                'field': str(key[0]), 'case_level': str(key[1]), 'control_level': str(key[2]),
                'stratum': None if key[3] is None else [str(key[3][0]), str(key[3][1])],
                'test': 'student' if equal_var else 'welch', 'recomputed_p': float(p),
                'n_case': len(a), 'n_control': len(b), 'effect_definitions_matched': matched,
                'recomputed_cohen_d': effects['cohen_d_pooled'],
                'case_gsms': sorted(map(str, a.index)), 'control_gsms': sorted(map(str, b.index))}
            break
    if not matches:
        return {'status': 'NONE', 'candidates_tested': len(seen), 'duplicated_sample_ids': duplicated[:20]}
    result = next(iter(matches.values()))
    return {'status': 'UNIQUE' if len(matches) == 1 else 'AMBIGUOUS', 'matching_memberships': len(matches),
            'candidates_tested': len(seen), 'duplicated_sample_ids': duplicated[:20], **result}


# ------------------------------------------------------------------ data access (worker globals)

_STATE = {}


def _init(modality, prediction_replay=None):
    config = MODALITIES[modality]
    gpt4 = pd.read_csv(GEO / 'result' / config['gpt4'], low_memory=False, dtype=str)
    gpt4['gse'] = gpt4.GSE_ID.str.extract(r'(GSE\d+)', expand=False)
    gpt4 = gpt4.drop(columns=['GSE_ID']).rename(columns={'GSM_ID': 'gsm'})
    _STATE.update(modality=modality, config=config, prediction_replay=prediction_replay,
                  prediction_inputs={}, gpt4={g: f.drop(columns='gse').set_index('gsm')
                                                            for g, f in gpt4.groupby('gse')})
    if modality == 'dnam':
        raw = pd.read_csv(GEO / 'result/hs_dnam_age.csv.gz').rename(columns={'GSE_ID': 'gse', 'GSM_ID': 'gsm'})
        dup = [c for c in raw.columns if c.startswith('DNAmGrimAge')]
        if len(dup) >= 2:  # the historical table stores GrimAge V1 then V2 under one name
            raw = raw.rename(columns={dup[0]: 'GrimAgeV1', dup[1]: 'GrimAgeV2'})
        qced = pd.read_csv(GEO / 'result/qced_human_dnam_age.csv').rename(columns={'gse_id': 'gse', 'gsm_id': 'gsm'})
        _STATE['dnam'] = {name: {g: f.set_index('gsm') for g, f in frame.groupby('gse')}
                          for name, frame in (('raw', raw), ('qced', qced))}


def _prediction_sources(gse, clock):
    if _STATE['modality'] == 'dnam':
        for name, by_study in _STATE['dnam'].items():
            frame = by_study.get(gse)
            if frame is not None and clock in frame.columns:
                yield name, pd.to_numeric(frame[clock], errors='coerce')
        return
    path = GEO / _STATE['config']['folder'] / 'clock_result' / f'{gse}.csv'
    if path.exists():
        frame = pd.read_csv(path).set_index('sample_id')
        for column in (clock, clock + '_mean'):
            if column in frame.columns:
                yield 'clock_result', pd.to_numeric(frame[column], errors='coerce')
                return
    # A traced deposit can omit clock_result files not opened by its original
    # self-test sample. The saved replay inputs contain per-sample predictions,
    # not repaired atlas statistics. Recovery still enumerates metadata splits
    # and must reproduce the planted row's P; no reference arms are copied.
    replay = _STATE.get('prediction_replay')
    if replay is not None:
        path = Path(replay) / 'inputs' / f'{gse}.csv'
        if not path.is_file():
            raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}')
        frame = pd.read_csv(path, low_memory=False)
        key = next(c for c in frame if c.lower() in ('gsm_id', 'gsm', 'sample_id'))
        frame = frame.set_index(key)
        for column in (clock, clock + '_mean'):
            if column in frame:
                _STATE['prediction_inputs'][(gse, clock)] = {
                    'path': str(path), 'sha256': sha256(path), 'column': column}
                yield 'replay_inputs', pd.to_numeric(frame[column], errors='coerce')
                return


def _metadata(gse):
    path = GEO / _STATE['config']['folder'] / 'processed_metadata' / f'{gse}.csv'
    meta = pd.read_csv(path, dtype=str).set_index('GSM') if path.exists() else pd.DataFrame()
    extra = _STATE['gpt4'].get(gse, pd.DataFrame()).copy()
    extra.columns = ['gpt4.' + c for c in extra.columns]
    if meta.empty:
        return extra
    return meta.join(extra, how='outer') if not extra.empty else meta


def _study(task):
    gse, rows = task
    try:
        return _study_rows(gse, rows)
    except Exception as error:  # one malformed study must not stop the corpus run
        records = [{'row_ordinal': row['row_ordinal'], 'gse_id': gse, 'comparison_name': row['comparison_name'],
                    'condition_name': row['condition_name'], 'clock': row['clock'],
                    'stored_test': row['statistical_test'], 'stored_p': row['pvalue'],
                    'stored_effect': row['effect_size'], 'status': 'STUDY_ERROR',
                    'error': f'{type(error).__name__}: {error}'[:300]} for row in rows]
        return gse, records, 0.0


def _study_rows(gse, rows):
    warnings.filterwarnings('ignore')
    started, out = time.time(), []
    meta = None
    for row in rows:
        record = {'row_ordinal': row['row_ordinal'], 'gse_id': gse, 'comparison_name': row['comparison_name'],
                  'condition_name': row['condition_name'], 'clock': row['clock'],
                  'stored_test': row['statistical_test'], 'stored_p': row['pvalue'],
                  'stored_effect': row['effect_size']}
        if not is_two_sample_t(row['statistical_test']):
            record['status'] = 'NOT_TWO_SAMPLE_T'
        elif not (np.isfinite(row['pvalue']) and 0 < row['pvalue'] < 1):
            record['status'] = 'STORED_P_UNUSABLE'
        else:
            meta = _metadata(gse) if meta is None else meta
            record['status'] = 'NO_PREDICTIONS'
            for source, values in _prediction_sources(gse, row['clock']):
                result = recover_row(values, meta, row['pvalue'], row['effect_size'])
                if result['status'] != 'NONE' or record['status'] == 'NO_PREDICTIONS':
                    record.update(result, prediction_source=source)
                    if source == 'replay_inputs':
                        record['prediction_input'] = _STATE['prediction_inputs'][(gse, row['clock'])]
                if result['status'] != 'NONE':
                    break
        out.append(record)
    return gse, out, time.time() - started


def run(modality, out_dir, workers=16, limit=None, resume=False, prediction_replay=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=resume)
    done = set()
    if resume and (out_dir / 'rows.jsonl').exists():
        with open(out_dir / 'rows.jsonl') as handle:
            for line in handle:
                try:
                    done.add(json.loads(line)['gse_id'])
                except json.JSONDecodeError:
                    break  # a truncated final line from an interrupted run
    table = CORRECTIONS / f"2026-09-08-{MODALITIES[modality]['table']}-atlas_WORKING_NOT_VALIDATED.csv"
    atlas = pd.read_csv(table, low_memory=False)
    atlas['row_ordinal'] = np.arange(len(atlas))
    atlas['pvalue'] = pd.to_numeric(atlas.pvalue, errors='coerce')
    atlas['effect_size'] = pd.to_numeric(atlas.effect_size, errors='coerce')
    columns = ['row_ordinal', 'comparison_name', 'condition_name', 'clock', 'statistical_test', 'pvalue', 'effect_size']
    tasks = [(gse, frame[columns].to_dict('records')) for gse, frame in atlas.groupby('gse_id', sort=True)]
    if limit:
        tasks = tasks[:limit]
    if done:
        # A study is written in one block, so only the last recorded study can be partial: redo it.
        last = None
        with open(out_dir / 'rows.jsonl') as handle:
            kept = []
            for line in handle:
                try:
                    kept.append(json.loads(line))
                except json.JSONDecodeError:
                    break
        last = kept[-1]['gse_id'] if kept else None
        expected = {gse: len(r) for gse, r in tasks}
        counts = {}
        for record in kept:
            counts[record['gse_id']] = counts.get(record['gse_id'], 0) + 1
        complete = {g for g, n in counts.items() if n == expected.get(g)}
        with open(out_dir / 'rows.jsonl', 'w') as handle:
            for record in kept:
                if record['gse_id'] in complete:
                    handle.write(json.dumps(record, default=str) + '\n')
        tasks = [(g, r) for g, r in tasks if g not in complete]
    manifest = {'schema': SCHEMA, 'modality': modality, 'table': str(table), 'table_sha256': sha256(table),
                'module_sha256': sha256(__file__), 'studies': len(tasks), 'rows': sum(len(r) for _, r in tasks),
                'tolerances': {'p_relative': P_RELATIVE_TOLERANCE, 'effect': EFFECT_TOLERANCE},
                'candidate_limits': {'max_levels': MAX_LEVELS, 'max_strata_levels': MAX_STRATA_LEVELS},
                'started': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    if done:
        manifest['resumed_skipping_studies'] = len(expected) - len(tasks)
    manifest['prediction_replay'] = str(prediction_replay) if prediction_replay is not None else None
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    counts = {}
    _init(modality, prediction_replay)  # fail here before starting pool initializers
    with open(out_dir / 'rows.jsonl', 'a' if done else 'w') as sink, Pool(
            workers, initializer=_init, initargs=(modality, prediction_replay)) as pool:
        for i, (gse, records, seconds) in enumerate(pool.imap_unordered(_study, tasks, chunksize=4)):
            for record in records:
                sink.write(json.dumps(record, default=str) + '\n')
                counts[record['status']] = counts.get(record['status'], 0) + 1
            if i % 500 == 0:
                print(f'{i}/{len(tasks)} studies {counts}', flush=True)
    counts = {}
    with open(out_dir / 'rows.jsonl') as handle:
        for line in handle:
            status = json.loads(line)['status']
            counts[status] = counts.get(status, 0) + 1
    manifest.update(finished=time.strftime('%Y-%m-%dT%H:%M:%S%z'), status_counts=counts,
                    rows_sha256=sha256(out_dir / 'rows.jsonl'))
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps(counts))
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    commands = parser.add_subparsers(dest='command', required=True)
    command = commands.add_parser('run')
    command.add_argument('--modality', choices=sorted(MODALITIES), required=True)
    command.add_argument('--out', type=Path, required=True)
    command.add_argument('--workers', type=int, default=16)
    command.add_argument('--limit', type=int)
    command.add_argument('--resume', action='store_true', help='keep complete studies already in OUT/rows.jsonl')
    command.add_argument('--prediction-replay', type=Path,
                         help='saved replay inputs to use when an RNA clock_result source is absent')
    args = parser.parse_args(argv)
    run(args.modality, args.out, args.workers, args.limit, args.resume, args.prediction_replay)


if __name__ == '__main__':
    main()
