"""Input-integrity screen: is every sample's clock input attached to the right features?

A bisulfite methylation profile from any human tissue correlates strongly with the population-mean
methylome (the GSE40279 means that the historical recipe also uses for imputation), because CpG-island
and non-island probes sit at opposite ends of the beta range. A sample whose values have been detached
from their probe identifiers (for example, a file merged by row position instead of probe ID) keeps a
plausible beta distribution and full probe coverage, so every coverage check passes, yet every clock
computed from it is meaningless. GSE73950 is the index case: its 17 control samples correlate with the
reference at r = -0.002 while its 26 cases correlate at 0.88-0.95, and the controls' historical
PhenoAge is -72 years.

For each historical DNAm input matrix (GEO_hs_dnam/processed_data/{GSE}.csv.gz) the screen streams
probe chunks and accumulates, per sample, the Pearson correlation with the reference means and with
the study's per-probe median, plus the share of values outside [0, 1]. Nothing is thresholded here;
the verdict layer applies the published threshold. No model calls.

  python -m audit_system.input_integrity_2026_09_30 dnam --out DIR [--workers 16] [--gse GSE ...]
  python -m audit_system.input_integrity_2026_09_30 summary --out DIR
"""
import argparse
import csv
import gzip
import json
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import GEO, sha256

SCHEMA = 'clockbase-input-integrity-v1'
DNAM_INPUTS = GEO / 'GEO_hs_dnam' / 'processed_data'
REFERENCE = paths.METHYLOME_REFERENCE
CELLS_PER_CHUNK = 8_000_000  # bounds worker memory: rows per chunk = cells / samples


class Correlation:
    """Per-sample Pearson correlation with a reference, accumulated over feature chunks (pairwise complete)."""

    def __init__(self, samples):
        self.n, self.sx, self.sy, self.sxx, self.syy, self.sxy = (np.zeros(samples) for _ in range(6))

    def add(self, x, y):
        y = np.broadcast_to(y[:, None], x.shape) if y.ndim == 1 else y
        mask = np.isfinite(x) & np.isfinite(y)
        xz, yz = np.where(mask, x, 0.0), np.where(mask, y, 0.0)
        self.n += mask.sum(0)
        self.sx += xz.sum(0)
        self.sy += yz.sum(0)
        self.sxx += (xz * xz).sum(0)
        self.syy += (yz * yz).sum(0)
        self.sxy += (xz * yz).sum(0)

    def r(self):
        numerator = self.n * self.sxy - self.sx * self.sy
        denominator = np.sqrt(np.clip(self.n * self.sxx - self.sx ** 2, 0, None) * np.clip(self.n * self.syy - self.sy ** 2, 0, None))
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(denominator > 0, numerator / denominator, np.nan)


def numeric(chunk):
    try:
        return chunk.to_numpy(dtype=float)
    except (TypeError, ValueError):
        return chunk.apply(pd.to_numeric, errors='coerce').to_numpy(dtype=float)


def screen_dnam(path, reference):
    """Per-sample integrity metrics for one historical DNAm input matrix."""
    with gzip.open(path, 'rt') as handle:
        samples = next(csv.reader([handle.readline()]))[1:]
    rows = max(1000, CELLS_PER_CHUNK // max(1, len(samples)))
    to_reference, to_median = Correlation(len(samples)), Correlation(len(samples))
    finite, outside, features, truncated = np.zeros(len(samples)), np.zeros(len(samples)), 0, False
    try:
        for chunk in pd.read_csv(path, index_col=0, chunksize=rows, low_memory=False):
            x = numeric(chunk)
            features += len(chunk)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                median = np.nanmedian(x, axis=1)
            to_reference.add(x, reference.reindex(chunk.index.astype(str)).to_numpy(dtype=float))
            to_median.add(x, median)
            is_finite = np.isfinite(x)
            finite += is_finite.sum(0)
            outside += (is_finite & ((x < 0) | (x > 1))).sum(0)
    except EOFError:
        truncated = True  # the gzip stream ends early (GSE51032 on disk); metrics use the chunks read
    return pd.DataFrame({'gsm': samples, 'features': features, 'input_truncated': truncated, 'finite': finite.astype(int),
                         'reference_pairs': to_reference.n.astype(int), 'r_reference': to_reference.r(),
                         'median_pairs': to_median.n.astype(int), 'r_study_median': to_median.r(),
                         'share_outside_unit': np.where(finite > 0, outside / np.maximum(finite, 1), np.nan)})


def _dnam_job(gse, out_dir):
    target = Path(out_dir) / 'studies' / f'{gse}.csv'
    started = time.time()
    try:
        reference = pd.read_csv(REFERENCE, index_col=0, converters={0: str})['mean']
        frame = screen_dnam(DNAM_INPUTS / f'{gse}.csv.gz', reference)
        frame.insert(0, 'gse_id', gse)
        partial = target.with_suffix('.partial')
        frame.to_csv(partial, index=False)
        partial.rename(target)
        return {'gse_id': gse, 'status': 'OK', 'samples': len(frame), 'seconds': round(time.time() - started, 1)}
    except Exception as error:  # recorded per study; the run continues
        return {'gse_id': gse, 'status': 'ERROR', 'error': f'{type(error).__name__}: {error}'[:500],
                'trace': traceback.format_exc()[-2000:], 'seconds': round(time.time() - started, 1)}


def run_dnam(out_dir, workers=16, gses=None):
    out_dir = Path(out_dir)
    (out_dir / 'studies').mkdir(parents=True, exist_ok=True)
    available = sorted(p.name[:-len('.csv.gz')] for p in DNAM_INPUTS.glob('*.csv.gz'))
    todo = [g for g in (gses or available) if not (out_dir / 'studies' / f'{g}.csv').exists()]
    run = {'schema': SCHEMA, 'modality': 'dnam', 'inputs': str(DNAM_INPUTS), 'studies_available': len(available),
           'studies_todo': len(todo), 'reference': str(REFERENCE), 'reference_sha256': sha256(REFERENCE),
           'module_sha256': sha256(__file__), 'started': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / f'run-{time.strftime("%Y%m%dT%H%M%S")}.json').write_text(json.dumps(run, indent=2))
    # Largest matrices first so the tail of the run is short.
    todo.sort(key=lambda g: -(DNAM_INPUTS / f'{g}.csv.gz').stat().st_size)
    with ProcessPoolExecutor(workers) as pool, open(out_dir / 'events.jsonl', 'a') as events:
        futures = [pool.submit(_dnam_job, g, out_dir) for g in todo]
        for future in as_completed(futures):
            events.write(json.dumps(future.result()) + '\n')
            events.flush()
    summarize(out_dir)


def summarize(out_dir):
    out_dir = Path(out_dir)
    parts = [pd.read_csv(p) for p in sorted((out_dir / 'studies').glob('*.csv'))]
    samples = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    samples.to_csv(out_dir / 'samples.csv', index=False)
    events = [json.loads(line) for line in open(out_dir / 'events.jsonl')] if (out_dir / 'events.jsonl').exists() else []
    errors = {}
    for event in events:  # the latest attempt per study decides
        errors[event['gse_id']] = event if event['status'] == 'ERROR' else None
    bins = [-1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    summary = {'schema': SCHEMA, 'studies': int(samples.gse_id.nunique()) if len(samples) else 0,
               'samples': int(len(samples)),
               'studies_with_errors': sorted(g for g, e in errors.items() if e),
               'r_reference_histogram': (pd.cut(samples.r_reference, bins).value_counts().sort_index()
                                         .rename(str).to_dict() if len(samples) else {}),
               'r_reference_missing': int(samples.r_reference.isna().sum()) if len(samples) else 0,
               'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'module_sha256': sha256(__file__)}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=1, default=str))


def detached_list(screen_dir, out_path, threshold=0.10):
    """Studies with any sample below the detached threshold, in the format of the GEO source inventory."""
    samples = pd.read_csv(Path(screen_dir) / 'samples.csv')
    studies = sorted(samples[samples.r_reference < threshold].gse_id.unique())
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as sink:
        for gse in studies:
            path = DNAM_INPUTS / f'{gse}.csv.gz'
            rows = samples[samples.gse_id == gse]
            sink.write(json.dumps({'study_id': gse, 'sample_columns': int(len(rows)), 'probe_rows': int(rows.features.iloc[0]),
                                   'matrix_path': str(path), 'matrix_sha256': sha256(path),
                                   'detached_samples': int((rows.r_reference < threshold).sum())}) + '\n')
    print(len(studies), 'studies ->', out_path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    dnam = sub.add_parser('dnam')
    dnam.add_argument('--out', type=Path, required=True)
    dnam.add_argument('--workers', type=int, default=16)
    dnam.add_argument('--gse', nargs='*')
    summary = sub.add_parser('summary')
    summary.add_argument('--out', type=Path, required=True)
    listing = sub.add_parser('detached-list')
    listing.add_argument('--screen', type=Path, required=True)
    listing.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == 'dnam':
        run_dnam(args.out, args.workers, args.gse)
    elif args.command == 'detached-list':
        detached_list(args.screen, args.out)
    else:
        summarize(args.out)


if __name__ == '__main__':
    main()
