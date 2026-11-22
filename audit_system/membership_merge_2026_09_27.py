"""Merge replay links and numerical recovery into one membership record per atlas row.

Replay links (the agent's own executed test) take precedence; numerical recovery fills rows the
replay could not link. Output uses the membership-recovery schema so row verdicts consume it
unchanged, with added provenance: source, misfiled, kind, covariates.

  python -m audit_system.membership_merge_2026_09_27 --modality dnam --recovery DIR --link DIR \
      --replay DIR --out DIR
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system.membership_recovery_2026_09_27 import GEO, MODALITIES, effect_sizes, sha256

SCHEMA = 'clockbase-membership-merged-v1'


def separating_levels(meta, case, control):
    """(field, case_level, control_level) for the metadata field that splits the arms cleanly, if any."""
    if meta is None:
        return None, None, None
    a, b = meta.reindex(case), meta.reindex(control)
    best = None
    for column in meta.columns:
        if column == 'GSM':
            continue
        va, vb = set(a[column].dropna().astype(str)), set(b[column].dropna().astype(str))
        if len(va) == 1 and len(vb) == 1 and va != vb:
            return column, next(iter(va)), next(iter(vb))
        if len(va) == 1 and va.isdisjoint(vb) and vb and best is None:
            best = (column, next(iter(va)), '|'.join(sorted(vb))[:120])
    return best or (None, None, None)


def load_meta(modality, gse):
    folder = MODALITIES[modality]['folder']
    for sub in ('processed_metadata_preharmonize', 'processed_metadata'):
        path = GEO / folder / sub / f'{gse}.csv'
        if path.exists():
            frame = pd.read_csv(path, dtype=str)
            key = 'GSM' if 'GSM' in frame.columns else frame.columns[0]
            return frame.drop_duplicates(key).set_index(key)
    return None


def study_values(replay_root, gse, clock):
    path = Path(replay_root) / 'inputs' / f'{gse}.csv'
    if not path.exists():
        return None
    table = pd.read_csv(path, low_memory=False)
    column = 'gsm_id' if 'gsm_id' in table.columns else next(c for c in table.columns if c.lower() in ('gsm', 'sample_id'))
    for name in (clock, clock + '_mean'):
        if name in table.columns:
            values = pd.to_numeric(table[name], errors='coerce')
            values.index = table[column].astype(str)
            return values[~values.index.duplicated(keep=False)]
    return None


def from_link(record, modality, replay_root, cache):
    link = record['link']
    gse, clock = record['gse_id'], record['clock']
    base = {'source': 'replay', 'kind': link['kind'], 'test': link['fn'], 'misfiled': record.get('misfiled', False),
            'link_traces': record.get('link_traces'), 'covariates': link.get('covariates', [])}
    if link['kind'] == 'one_sample':
        return {**base, 'status': 'ONE_SAMPLE_TEST', 'members': link['members']}
    if link.get('continuous_term'):
        return {**base, 'status': 'CONTINUOUS_TERM', 'members': link['members'], 'term': link.get('term')}
    case, control = link.get('case_gsms'), link.get('control_gsms')
    if not case or not control:
        return {**base, 'status': 'ARMS_UNRESOLVED'}
    key = (gse, clock)
    if key not in cache:
        cache[key] = study_values(replay_root, gse, clock)
    values = cache[key]
    matched = []
    if values is not None:
        a, b = values.reindex(case).dropna(), values.reindex(control).dropna()
        if len(a) >= 2 and len(b) >= 2:
            effects = effect_sizes(a, b)
            matched = [k for k, v in effects.items() if np.isfinite(v)
                       and abs(v - record['stored_effect']) <= 1e-4 * max(1.0, abs(record['stored_effect']))]
            if link['kind'] == 'model' and np.sign(effects['mean_diff']) == np.sign(record['stored_effect']):
                matched.append('model_direction_sign')
            stored = record['stored_effect']
            close = lambda v: np.isfinite(v) and abs(v - stored) <= 1e-4 * max(1.0, abs(stored))
            if link.get('paired_case'):
                diff = values.reindex(link['paired_case']).to_numpy(float) - values.reindex(link['paired_control']).to_numpy(float)
                diff = diff[np.isfinite(diff)]
                if len(diff) >= 2 and diff.std(ddof=1) > 0:
                    if close(diff.mean() / diff.std(ddof=1)):
                        matched.append('paired_d_z')
                    if close(diff.mean() / np.sqrt((a.var() + b.var()) / 2)):
                        matched.append('paired_d_av')
            param = link.get('term_param')
            if param is not None:
                pooled_sd = pd.concat([a, b]).std()
                if pooled_sd and close(param / pooled_sd):
                    matched.append('model_coefficient_standardized')
    param = link.get('term_param')
    if param is not None and abs(param - record['stored_effect']) <= 1e-4 * max(1.0, abs(record['stored_effect'])):
        matched.append('model_coefficient')
    if gse not in cache:
        cache[gse] = load_meta(modality, gse)
    field, case_level, control_level = separating_levels(cache[gse], case, control)
    return {**base, 'status': 'UNIQUE', 'case_gsms': case, 'control_gsms': control, 'n_case': len(case),
            'n_control': len(control), 'effect_definitions_matched': matched, 'field': field,
            'case_level': case_level, 'control_level': control_level, 'stratum': None, 'term': link.get('term')}


def run(modality, recovery_dir, link_dir, replay_root, out_dir, design_dirs=()):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    recovery = {r['row_ordinal']: r for r in map(json.loads, open(Path(recovery_dir) / 'rows.jsonl'))}
    designed = {}
    for folder in design_dirs:
        path = Path(folder) / 'rows.jsonl'
        if path.exists():
            for line in open(path):
                d = json.loads(line)
                designed[d['row_ordinal']] = d
    cache, counts = {}, {}
    with open(out_dir / 'rows.jsonl', 'w') as sink:
        for line in open(Path(link_dir) / 'row_links.jsonl'):
            link = json.loads(line)
            ordinal = link['row_ordinal']
            base = {k: link[k] for k in ('row_ordinal', 'gse_id', 'comparison_name', 'condition_name', 'clock',
                                          'stored_test', 'stored_p', 'stored_effect')}
            rec = recovery.get(ordinal, {})
            if link['link_status'] == 'LINKED':
                merged = from_link(link, modality, replay_root, cache)
                if rec.get('status') == 'UNIQUE' and merged.get('status') == 'UNIQUE':
                    merged['source'] = 'replay+recovery'
                    merged['recovery_agrees'] = sorted(rec['case_gsms']) == sorted(merged['case_gsms'])
            elif rec.get('status') == 'UNIQUE':
                merged = {k: rec[k] for k in ('case_gsms', 'control_gsms', 'n_case', 'n_control', 'field', 'case_level',
                                              'control_level', 'stratum', 'test', 'effect_definitions_matched')}
                merged.update(status='UNIQUE', source='recovery', kind='two_sample', misfiled=False)
            else:
                merged = {'status': {'STORED_P_UNUSABLE': 'STORED_P_UNUSABLE', 'AMBIGUOUS': 'AMBIGUOUS'}.get(
                    link['link_status'], rec.get('status') if rec.get('status') not in (None, 'NONE') else 'NONE'),
                          'source': None, 'replay_link_status': link['link_status']}
            if merged.get('status') != 'UNIQUE' and ordinal in designed:
                d = designed[ordinal]
                merged = {'status': 'UNIQUE', 'source': 'design_reconstruction', 'kind': 'two_sample',
                          'case_gsms': d['case_gsms'], 'control_gsms': d['control_gsms'], 'n_case': d['n_case'],
                          'n_control': d['n_control'], 'design_rules': d.get('design_rules'), 'field': None,
                          'case_level': None, 'control_level': None, 'stratum': None,
                          'historical_status': merged.get('status')}
            record = {**base, **merged}
            counts[record['status']] = counts.get(record['status'], 0) + 1
            sink.write(json.dumps(record, default=str) + '\n')
            if len(cache) > 400:
                cache.clear()
    manifest = {'schema': SCHEMA, 'modality': modality, 'recovery': str(recovery_dir), 'link': str(link_dir),
                'recovery_sha256': sha256(Path(recovery_dir) / 'rows.jsonl'),
                'link_sha256': sha256(Path(link_dir) / 'row_links.jsonl'), 'module_sha256': sha256(__file__),
                'status_counts': counts, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                'rows_sha256': sha256(out_dir / 'rows.jsonl')}
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps(counts))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--modality', choices=sorted(MODALITIES), required=True)
    for name in ('recovery', 'link', 'replay', 'out'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--design', type=Path, action='append', default=[])
    args = parser.parse_args(argv)
    run(args.modality, args.recovery, args.link, args.replay, args.out, args.design)


if __name__ == '__main__':
    main()
