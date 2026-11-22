"""Read the verifier's prediction sources and recorded sample demographics."""
import os
import threading
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import paths
from audit_system import membership_recovery_2026_09_27 as membership
from audit_system.row_repairs_2026_09_28 import repaired_dnam_values, study_values
from .io import read_json, receipt

# The verifier run the benchmark scores (CLOCKBASE_FINAL_NAME) and the folder of benchmark summaries built on it.
FINAL_NAME = os.environ.get('CLOCKBASE_FINAL_NAME', 'final-2026-10-05')
SUMMARY_NAME = os.environ.get('CLOCKBASE_BENCH_SUMMARY', FINAL_NAME)


def final_export():
    return paths.VERIFIER_ROOT / FINAL_NAME / 'public-export'


MODALITIES = {'human DNAm': 'dnam', 'human RNA-seq': 'hs_rna', 'mouse RNA-seq': 'mm_rna',
              'dnam': 'dnam', 'hs_rna': 'hs_rna', 'mm_rna': 'mm_rna'}
SOURCE_LOCK = threading.Lock()


def modality_key(modality):
    if modality not in MODALITIES:
        raise ValueError(f'Unsupported modality: {modality}')
    return MODALITIES[modality]


def export_path(modality, gse, root=None):
    root = Path(root) if root is not None else final_export()
    return root / 'studies' / modality_key(modality) / f'{gse}.json'


def atlas_rows(modality, gse, root=None):
    path = export_path(modality, gse, root)
    return read_json(path).get('rows', []) if path.exists() else []


class PredictionReader:
    def __init__(self):
        self.cache, self.tables, self.receipts = {}, {}, {}

    def source_receipt(self, path):
        key = str(path)
        if key not in self.receipts:
            self.receipts[key] = receipt(path)
        return self.receipts[key]

    def metadata(self, modality, gse):
        folder = membership.MODALITIES[modality_key(modality)]['folder']
        path = paths.GEO / folder / 'processed_metadata' / f'{gse}.csv'
        if not path.exists():
            return pd.DataFrame(), None
        try:
            frame = pd.read_csv(path, dtype=str).set_index('GSM')
        except (OSError, ValueError, KeyError) as error:
            raise ValueError(f'Cannot read sample metadata {path}: {error}') from error
        frame.index = frame.index.astype(str)
        return frame[~frame.index.duplicated(keep=False)], self.source_receipt(path)

    def historical_sources(self, modality, gse, clock):
        """Delegate column choice to membership recovery without leaving its globals changed."""
        key = modality_key(modality)
        state = {'modality': key, 'config': membership.MODALITIES[key]}
        source_paths = {}
        if key == 'dnam':
            state['dnam'] = {}
            for name, filename in [('raw', 'hs_dnam_age.csv.gz'), ('qced', 'qced_human_dnam_age.csv')]:
                path = paths.GEO / 'result' / filename
                if not path.exists():
                    continue
                if str(path) not in self.tables:
                    try:
                        frame = pd.read_csv(path, low_memory=False).rename(
                            columns={'GSE_ID': 'gse', 'GSM_ID': 'gsm', 'gse_id': 'gse', 'gsm_id': 'gsm'})
                    except (OSError, ValueError) as error:
                        raise ValueError(f'Cannot read predictions {path}: {error}') from error
                    dup = [c for c in frame.columns if c.startswith('DNAmGrimAge')]
                    if name == 'raw' and len(dup) >= 2:
                        frame = frame.rename(columns={dup[0]: 'GrimAgeV1', dup[1]: 'GrimAgeV2'})
                    self.tables[str(path)] = {g: f.set_index('gsm') for g, f in frame.groupby('gse')}
                state['dnam'][name] = self.tables[str(path)]
                source_paths[name] = path
        else:
            source_paths['clock_result'] = paths.GEO / state['config']['folder'] / 'clock_result' / f'{gse}.csv'
        with SOURCE_LOCK:
            previous, previous_geo = membership._STATE.copy(), membership.GEO
            try:
                membership._STATE.clear()
                membership._STATE.update(state)
                membership.GEO = paths.GEO
                results = list(membership._prediction_sources(gse, clock))
            finally:
                membership._STATE.clear()
                membership._STATE.update(previous)
                membership.GEO = previous_geo
        return [(name, values, source_paths[name]) for name, values in results]

    def load(self, modality, gse, clock, rows=()):
        """Prefer the exported verifier source; otherwise repaired, QCed, raw, or RNA study outputs."""
        meta, meta_source = self.metadata(modality, gse)
        if not clock:
            return {'status': 'NO_CLOCK_MAPPING', 'samples': pd.DataFrame(), 'metadata': meta,
                    'metadata_source': meta_source}
        selected = [r for r in rows if r.get('clock') == clock]
        recorded = {}
        conflicts = []
        for row in selected:
            item = row.get('predictions') or {}
            if item.get('path'):
                previous = recorded.get(item['path'])
                if previous and previous.get('sha256') != item.get('sha256'):
                    conflicts.extend([previous, item])
                recorded[item['path']] = item
        values, source, verdict, name = None, None, None, None
        if len(recorded) > 1 or conflicts:
            return {'status': 'AMBIGUOUS_PREDICTION_SOURCE', 'sources': conflicts or list(recorded.values()),
                    'samples': pd.DataFrame(), 'metadata': meta, 'metadata_source': meta_source}
        if recorded:
            stored = next(iter(recorded.values()))
            path = paths.resolve(stored['path'])
            if not path.exists():
                return {'status': 'NO_PREDICTIONS', 'error': f'Cannot proceed - required data file not found: {path}',
                        'samples': pd.DataFrame(), 'metadata': meta, 'metadata_source': meta_source}
            source = dict(self.source_receipt(path))
            if stored.get('sha256') and stored['sha256'] != source['sha256']:
                return {'status': 'SOURCE_HASH_MISMATCH', 'source': source, 'expected_sha256': stored['sha256'],
                        'samples': pd.DataFrame(), 'metadata': meta, 'metadata_source': meta_source}
            if path.name == 'clock_predictions.csv':
                status_source = self.source_receipt(path.parent / 'status.json')
                if stored.get('status_sha256') and stored['status_sha256'] != status_source['sha256']:
                    return {'status': 'SOURCE_HASH_MISMATCH', 'source': source, 'status_source': status_source,
                            'expected_sha256': stored['status_sha256'], 'samples': pd.DataFrame(),
                            'metadata': meta, 'metadata_source': meta_source}
                values, _, verdict = repaired_dnam_values(path.parents[2], gse, clock, self.cache)
                for filename in ('status.json', 'sample_status.csv'):
                    source[filename] = self.source_receipt(path.parent / filename)
                name = 'stream_repaired'
            elif path.parent.name == 'inputs':
                values, _ = study_values(path.parent.parent, gse, clock, self.cache)
                name = 'verifier_replay'
            else:
                available = self.historical_sources(modality, gse, clock)
                values = next((v for _, v, p in available if p.resolve() == path.resolve()), None)
                name = 'recorded_historical'
        else:
            if modality_key(modality) == 'dnam':
                for root in (paths.STREAM_EMPTY, paths.STREAM_DETACHED):
                    values, _, verdict = repaired_dnam_values(root, gse, clock, self.cache)
                    if values is not None:
                        path = root / 'studies' / gse / 'clock_predictions.csv'
                        source = dict(self.source_receipt(path))
                        for filename in ('status.json', 'sample_status.csv'):
                            source[filename] = self.source_receipt(path.parent / filename)
                        name = 'stream_repaired'
                        break
            if values is None:
                available = self.historical_sources(modality, gse, clock)
                preferred = next((r.get('evidence', {}).get('prediction_source') for r in selected
                                  if r.get('evidence', {}).get('prediction_source')), None)
                available.sort(key=lambda x: (x[0] != (preferred or 'qced'), x[0]))
                if available:
                    name, values, path = available[0]
                    source = self.source_receipt(path)
        if values is None:
            return {'status': 'NO_PREDICTIONS', 'source': source, 'metadata_source': meta_source,
                    'samples': pd.DataFrame(), 'metadata': meta}
        values.index = values.index.astype(str)
        values, duplicated = membership.unique_samples(values)
        samples = pd.DataFrame({'clock': values.replace([np.inf, -np.inf], np.nan)})
        age_field = next((c for c in ('harmonized.age', 'age', 'chronological_age') if c in meta), None)
        sex_field = next((c for c in ('harmonized.sex', 'sex', 'gender') if c in meta), None)
        samples['age'] = pd.to_numeric(meta[age_field], errors='coerce') if age_field else np.nan
        samples['sex'] = meta[sex_field] if sex_field else None
        withheld = verdict and verdict.get('verdict') == 'WITHHELD'
        return {'status': 'CLOCK_WITHHELD' if withheld else 'OK', 'samples': samples, 'metadata': meta,
                'source': source, 'source_kind': name, 'metadata_source': meta_source,
                'age_field': age_field, 'sex_field': sex_field, 'clock_verdict': verdict,
                'duplicated_sample_ids': duplicated}
