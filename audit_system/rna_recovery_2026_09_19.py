"""Conditional reconstruction of a saved YuGene differential table.

This adapter cannot certify raw normalization, historical invocation or clock
applicability. Original files are read only; existing output paths are refused.
"""
from __future__ import annotations

import csv
import gzip
import io
import json
from pathlib import Path
import platform
import re
import time

import numpy as np

from audit_system.hashing import inspect_file

CONTRACT = 'saved-yugene-cohort-median-recovery-2026-09-19-v2'
NUMBER = re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?')
GSM = re.compile(r'GSM[0-9]+')


def file_record(path):
    record = inspect_file(path)
    return {'path': str(record.path), 'sha256': record.sha256,
            'size_bytes': record.byte_count, 'mtime_ns': Path(path).stat().st_mtime_ns}


def _open(path):
    return gzip.open(path, 'rt', newline='') if str(path).endswith('.gz') else Path(path).open(newline='')


def _fields(line):
    if '\x00' in line:
        raise ValueError('NUL byte in source CSV')
    return next(csv.reader([line], strict=True))


def _header(line):
    if not line.endswith('\n'):
        raise ValueError('Incomplete header')
    fields = _fields(line)
    if not fields or fields[0] not in ('', 'sample_id') or len(fields) < 2:
        raise ValueError('Expected explicit sample_id or retained R row-name header')
    genes = fields[1:]
    if len(set(genes)) != len(genes) or any(not g or g.strip() != g for g in genes):
        raise ValueError('Ambiguous feature identities')
    return genes


def _values(tokens, *, yugene=False):
    if any(t != 'NA' and not NUMBER.fullmatch(t) for t in tokens):
        raise ValueError('Only finite numeric tokens or the documented R NA token are accepted')
    values = np.fromiter((np.nan if t == 'NA' else float(t) for t in tokens), dtype=float)
    finite = values[~np.isnan(values)]
    if not np.isfinite(finite).all():
        raise ValueError('Infinite source value')
    if yugene and ((finite < 0).any() or (finite > 1).any()):
        raise ValueError('Observed YuGene values are outside [0,1]')
    return values


def read_saved_table(path, *, expected_samples=None, expected_genes=None, keep_values=True, yugene=False):
    """Strict row-width, identity and lexical scan before numerical use."""
    before = file_record(path)
    samples, arrays = [], []
    with _open(path) as handle:
        genes = _header(handle.readline())
        if expected_genes is not None and genes != list(expected_genes):
            raise ValueError('Feature identity/order mismatch: ' + str(path))
        for number, line in enumerate(handle, 2):
            if not line.endswith('\n'):
                raise ValueError(f'Unterminated source record at line {number}')
            row = _fields(line)
            if len(row) != len(genes) + 1:
                raise ValueError(f'Ragged source record at line {number}')
            if not GSM.fullmatch(row[0]) or row[0] in samples:
                raise ValueError(f'Invalid or duplicate sample identity at line {number}')
            samples.append(row[0])
            values = _values(row[1:], yugene=yugene)
            if keep_values:
                arrays.append(values)
    if not samples:
        raise ValueError('No source samples')
    if expected_samples is not None and samples != list(expected_samples):
        raise ValueError('Sample identity/order mismatch: ' + str(path))
    if file_record(path)['sha256'] != before['sha256']:
        raise ValueError('Source changed while reading: ' + str(path))
    return {'source': before, 'samples': samples, 'genes': genes,
            'values': np.vstack(arrays) if keep_values else None}


def r_median_differences(values):
    """R median(..., na.rm=FALSE) across every source sample for each gene."""
    if values.ndim != 2 or not len(values):
        raise ValueError('Expected a nonempty sample-by-gene matrix')
    # np.median, unlike nanmedian or pandas' default, propagates missingness.
    medians = np.median(values, axis=0)
    return values - medians, medians


def compare_saved_prefix(path, samples, genes, candidate, *, atol):
    """Compare terminated fields; retain the final EOF token without validating it."""
    before = file_record(path)
    result = {'source': before, 'complete_rows': 0, 'partial_rows': 0,
              'checked_cells': 0, 'checked_numeric_cells': 0, 'checked_na_cells': 0,
              'max_abs_difference': 0.0, 'unverified_final_token': None}
    with _open(path) as handle:
        if _header(handle.readline()) != genes:
            raise ValueError('Legacy prefix feature identity/order mismatch')
        row_index = 0
        for line in handle:
            complete = line.endswith('\n')
            row = _fields(line)
            if row_index >= len(samples) or not row or row[0] != samples[row_index]:
                raise ValueError('Legacy prefix sample identity/order mismatch')
            if complete:
                if len(row) != len(genes) + 1:
                    raise ValueError('Ragged terminated legacy row')
                tokens = row[1:]
                result['complete_rows'] += 1
            else:
                if len(row) < 2 or len(row) > len(genes) + 1 or result['partial_rows']:
                    raise ValueError('Ambiguous final legacy record')
                # The last token may parse as a number but is not delimiter-terminated.
                tokens = row[1:-1]
                result['partial_rows'] += 1
                result['partial_sample_id'] = row[0]
                result['partial_terminated_cells'] = len(tokens)
                result['unverified_final_token'] = row[-1]
                result['unverified_final_token_feature'] = genes[len(tokens)]
            observed = _values(tokens)
            expected = candidate[row_index, :len(tokens)]
            if not np.array_equal(np.isnan(observed), np.isnan(expected)):
                raise ValueError('Legacy prefix NA mask differs from reconstruction')
            finite = ~np.isnan(observed)
            error = np.abs(observed[finite] - expected[finite])
            if (error > atol).any():
                raise ValueError(f'Legacy prefix values differ beyond atol={atol} at {row[0]}')
            result['checked_cells'] += len(tokens)
            result['checked_numeric_cells'] += int(finite.sum())
            result['checked_na_cells'] += int((~finite).sum())
            if len(error):
                result['max_abs_difference'] = max(result['max_abs_difference'], float(error.max()))
            row_index += 1
    if result['checked_cells'] == 0:
        raise ValueError('No intact legacy cells to verify')
    if file_record(path)['sha256'] != before['sha256']:
        raise ValueError('Legacy prefix changed while reading')
    result['status'] = 'PASS'
    return result


def validate_plan_scope(plan):
    """V2 can inspect arithmetic, but cannot turn a source-domain failure into PASS."""
    if plan.get('schema') == 'clockbase-rna-truncation-recovery-plan-v1':
        return False
    if plan.get('schema') != 'clockbase-rna-truncation-arithmetic-plan-v2':
        raise ValueError('Unsupported recovery plan')
    if (plan.get('arithmetic_only') is not True
            or plan.get('source_value_policy') != 'preserve'
            or plan.get('repair_must_remain_unverified') is not True):
        raise ValueError('Arithmetic-only plan must preserve values and withhold repair certification')
    prior = plan['v1_plan']
    if file_record(prior['path'])['sha256'] != prior['sha256']:
        raise ValueError('Original v1 plan changed')
    original = json.loads(Path(prior['path']).read_text())
    for key in ('inputs', 'definition_sources', 'studies'):
        if plan[key] != original[key]:
            raise ValueError('Arithmetic plan changed original ' + key)
    for key in ('comparison_atol', 'comparison_rtol', 'missingness', 'definition', 'comparison_scope'):
        if plan['reconstruction'][key] != original['reconstruction'][key]:
            raise ValueError('Arithmetic plan changed original ' + key)
    return True


def recover_study(plan_path, study, output):
    """Execute one hash-pinned study, preserving failures and original sources."""
    plan_path, output = Path(plan_path), Path(output)
    plan = json.loads(plan_path.read_text())
    arithmetic_only = validate_plan_scope(plan)
    if study not in plan['studies']:
        raise ValueError('Study is absent from recovery plan')
    definition = plan['reconstruction']
    if definition['comparison_atol'] != 1e-12 or definition['comparison_rtol'] != 0:
        raise ValueError('This adapter implements only the frozen tolerance')
    sources = {s['role']: s for s in plan['inputs'] if s['study_id'] == study}
    for source in list(sources.values()) + plan['definition_sources']:
        path = Path(source['path'])
        if path.exists() != source['available']:
            raise ValueError('Pinned source availability changed: ' + str(path))
        if source['available'] and file_record(path)['sha256'] != source['sha256']:
            raise ValueError('Pinned source content changed: ' + str(path))
    output.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    receipt = {'contract': CONTRACT, 'study_id': study, 'plan': file_record(plan_path),
               'started_unix': time.time(),
               'arithmetic_only': arithmetic_only,
               'implementation': file_record(__file__), 'inputs': sources,
               'versions': {'python': platform.python_version(), 'numpy': np.__version__},
               'detection': 'FAIL', 'scientific_status': 'CANNOT VERIFY',
               'independent_verification': 'CANNOT VERIFY', 'clock_execution': False,
               'scope': 'Conditional reconstruction from saved normalization; no historical runtime, raw preprocessing or clock applicability certification.'}
    try:
        table = read_saved_table(sources['yugene_data']['path'], yugene=not arithmetic_only)
        samples, genes, values = table['samples'], table['genes'], table['values']
        outside = (values < 0) | (values > 1)
        observed = values[~np.isnan(values)]
        if not len(observed):
            raise ValueError('No observed finite source values')
        receipt['source_domain'] = {'status': 'FAIL' if outside.any() else 'PASS',
                                    'rule': 'Every observed value in [0,1]',
                                    'out_of_range_cells': int(outside.sum()),
                                    'out_of_range_samples': int(outside.any(axis=1).sum()),
                                    'observed_min': float(observed.min()),
                                    'observed_max': float(observed.max()),
                                    'values_clamped_or_changed': False,
                                    'negative_token_origin': 'CANNOT VERIFY'}
        for role in ('scaled_data', 'scaled_diff_data'):
            read_saved_table(sources[role]['path'], expected_samples=samples,
                             expected_genes=genes, keep_values=False)
        with _open(sources['raw_counts']['path']) as handle:
            raw_line = handle.readline()
            if not raw_line.endswith('\n'):
                raise ValueError('Raw count header is not terminated')
            raw = _fields(raw_line)[1:]
        if len(set(raw)) != len(raw) or set(raw) != set(samples):
            raise ValueError('Raw count header and saved sample identities differ')
        candidate, medians = r_median_differences(values)
        prefix = compare_saved_prefix(sources['yugene_diff_data']['path'], samples, genes,
                                      candidate, atol=definition['comparison_atol'])
        if prefix['partial_rows'] != 1:
            raise ValueError('Recovery requires one observed truncated final record')
        target = output / (study + '_yugene_diff_candidate.csv.gz')
        with target.open('xb') as binary:
            with gzip.GzipFile(fileobj=binary, filename='', mode='wb', mtime=0, compresslevel=6) as compressed:
                with io.TextIOWrapper(compressed, newline='') as text:
                    writer = csv.writer(text, lineterminator='\n')
                    writer.writerow(['sample_id'] + genes)
                    for sample, row in zip(samples, candidate):
                        writer.writerow([sample] + ['NA' if np.isnan(v) else format(v, '.17g') for v in row])
        for source in list(sources.values()) + plan['definition_sources']:
            if Path(source['path']).exists() != source['available'] or (
                source['available'] and file_record(source['path'])['sha256'] != source['sha256']):
                raise ValueError('Pinned source changed during reconstruction: ' + source['path'])
        receipt.update(repair_execution='CANNOT VERIFY' if arithmetic_only else 'PASS',
                       conditional_arithmetic='PASS', prefix=prefix, candidate=file_record(target),
                       sample_count=len(samples), feature_count=len(genes),
                       candidate_cell_count=int(candidate.size),
                       candidate_na_cells=int(np.isnan(candidate).sum()),
                       genes_with_na_median=int(np.isnan(medians).sum()),
                       raw_header_same_order=raw == samples,
                       candidate_sample_ids=samples)
    except Exception as error:
        receipt.update(repair_execution='CANNOT VERIFY' if arithmetic_only else 'FAIL',
                       conditional_arithmetic='FAIL', error=f'{type(error).__name__}: {error}')
    receipt['elapsed_seconds'] = time.monotonic() - start
    receipt['finished_unix'] = time.time()
    with (output / 'recovery-receipt.json').open('x') as handle:
        json.dump(receipt, handle, indent=2, allow_nan=False)
        handle.write('\n')
    return receipt


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', required=True, type=Path)
    parser.add_argument('--study', required=True)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = recover_study(args.plan, args.study, args.output)
    print(json.dumps({'study_id': args.study, 'repair_execution': result['repair_execution'],
                      'conditional_arithmetic': result.get('conditional_arithmetic'),
                      'error': result.get('error')}))
    raise SystemExit(0 if result.get('conditional_arithmetic') == 'PASS' else 1)
