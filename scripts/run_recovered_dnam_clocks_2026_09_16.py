"""Run existing DNAm producers after verified source recovery, without model calls.

All samples and all measured rows are retained. Numerical verification does not
establish clock applicability, biological units, treatment effects or historical
producing lineage. Run from the repository with the pinned pipeline Python.
"""
import argparse
import csv
import fcntl
import gzip
from importlib.metadata import version
import json
from pathlib import Path
import re
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import biolearn
import numpy as np
import pandas as pd

from audit_system.agent_workflow_2026_09_10 import sha, write_new
from audit_system.report_2026_09_10 import recover_dnam_input
from src.pipeline_inputs_2026_09_10 import preprocess_dnam, predict_dnam
from src.pipeline_runtime_2026_09_10 import StageStore


DIRECT_PLAN_SCHEMA = 'clockbase-direct-series-clock-repair-plan-v1'
DIRECT_PLAN_SCHEMA_V2 = 'clockbase-direct-series-clock-repair-plan-v2'


def summarize(predictions, recovery, historical, directory, source_mode='recovery'):
    fresh = pd.read_csv(predictions['path'], index_col=0)
    if (fresh.index.has_duplicates or len(fresh) != recovery['source_samples']
            or not all(re.fullmatch(r'GSM[0-9]+', str(s)) for s in fresh.index)
            or not np.isfinite(fresh.to_numpy(float)).all()):
        raise ValueError('Fresh clock sample identity, dimension or values invalid')
    old = None
    if historical:
        # Duplicate unrelated historical columns remain untouched. Only the three
        # explicitly named, unique clock columns and GSM IDs enter this comparison.
        with historical.open(newline='') as handle:
            header = next(csv.reader(handle))
        if any(header.count(name) != 1 for name in predictions['clocks']):
            raise ValueError('Historical selected clock columns are missing/ambiguous')
        old = pd.read_csv(historical, index_col=0)[predictions['clocks']]
        if old.index.has_duplicates or set(old.index) != set(fresh.index):
            raise ValueError('Historical/fresh GSM membership differs')
        old = old.loc[fresh.index]
        if not np.isfinite(old.to_numpy(float)).all():
            raise ValueError('Selected historical clock values are nonfinite')
    rows = []
    for check in predictions['checks']:
        name = check['clock']
        row = {'study': recovery['study_id'], 'clock': name, 'samples': len(fresh),
               'required_cpgs': check['required_features'],
               'observed_cpgs_min': check['observed_features_min'],
               'imputed_cpgs_max': check['required_features'] - check['observed_features_min'],
               'fresh_min': float(fresh[name].min()), 'fresh_max': float(fresh[name].max()),
               'fresh_unique_values': int(fresh[name].nunique()),
               'coefficient_verification': check['numerical_status'],
               'max_independent_difference': check['max_abs_independent_difference'],
               'measurement_unit_applicability': 'CANNOT VERIFY',
               'biological_sample_unit_status': 'CANNOT VERIFY',
               'scientific_status': 'CANNOT VERIFY'}
        if old is not None:
            row.update(historical_min=float(old[name].min()), historical_max=float(old[name].max()),
                       historical_unique_values=int(old[name].nunique()),
                       max_absolute_change=float((fresh[name]-old[name]).abs().max()))
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(directory/'clock-validation.csv', index=False, mode='x')
    fields = ['clock','samples','observed_cpgs_min','required_cpgs',
              'fresh_unique_values','coefficient_verification','scientific_status']
    if old is not None: fields.insert(4, 'historical_unique_values')
    lines = ['| ' + ' | '.join(fields) + ' |', '| ' + ' | '.join(['---']*len(fields)) + ' |']
    lines.extend('| ' + ' | '.join(str(row[key]) for key in fields) + ' |' for row in rows)
    (directory/'clock-validation.md').write_text('\n'.join(lines)+'\n')
    scope = ('Fresh canonical direct-Series beta-to-clock arithmetic using the existing producer; '
             'numerical candidate only, with no historical vector comparison'
             if source_mode == 'direct_series' else
             'Fresh beta-to-clock arithmetic using existing producer; no group contrast, '
             'historical overwrite or accepted-result error rate')
    result = {'study_id': recovery['study_id'], 'sample_count': len(fresh), 'rows': rows,
              'generation_calls': 0, 'scientific_status': 'CANNOT VERIFY',
              'historical_producing_lineage': 'CANNOT VERIFY',
              'measurement_unit_applicability': 'CANNOT VERIFY',
              'biological_sample_unit_status': 'CANNOT VERIFY',
              'source_mode': source_mode, 'candidate_status': 'NUMERICAL_CANDIDATE',
              'scope': scope}
    write_new(directory/'summary.json', result)
    return result


def _bound_file(record, label):
    path = Path(record['path']).resolve()
    if not path.is_file():
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
    if sha(path) != record['sha256']:
        raise ValueError('Direct-Series bound file changed: ' + label)
    if 'size_bytes' in record and path.stat().st_size != record['size_bytes']:
        raise ValueError('Direct-Series bound file size changed: ' + label)
    return path


def _primary_series_sample_ids(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.reader(iter(handle.readline, ''), delimiter='\t', strict=True)
        for row in reader:
            if row == ['!series_matrix_table_begin']:
                header = next(reader, [])
                break
        else:
            raise ValueError('AUDIT_FAIL: GEO matrix begin marker missing')
    if len(header) < 2 or header[0] != 'ID_REF':
        raise ValueError('AUDIT_FAIL: invalid GEO primary matrix header')
    return header[1:]


def _metadata_sample_ids(path):
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or 'GSM' not in reader.fieldnames:
            raise ValueError('AUDIT_FAIL: metadata lacks the GSM identity column')
        identities = [row.get('GSM', '') for row in reader]
    if not identities:
        raise ValueError('AUDIT_FAIL: metadata has no sample records')
    if (len(identities) != len(set(identities))
            or any(re.fullmatch(r'GSM[0-9]+', value) is None for value in identities)):
        raise ValueError('AUDIT_FAIL: metadata GSM identities are missing, noncanonical or duplicated')
    return identities


def _repository_path(value):
    path = Path(value)
    return (path if path.is_absolute() else ROOT/path).resolve()


def _validate_independent_readback(plan, studies, record):
    receipt_path = _bound_file(record, 'independent source cross-readback')
    receipt = json.loads(receipt_path.read_text())
    required_checks = {
        'source_sha256', 'sample_count', 'feature_count', 'sample_ids',
        'full_header_sha256', 'feature_order_sha256', 'numeric_stream_sha256',
        'finite_count', 'missing_count', 'out_of_beta_count',
        'source_contract_status'}
    if (record.get('status_scope') != 'saved-source-cross-readback-only'
            or record.get('scientific_status') != 'CANNOT VERIFY'
            or record.get('source_origin_plan_sha256')
            != plan['source_origin_plan']['sha256']
            or record.get('root_summary_sha256')
            != plan['source_origin_summary']['sha256']
            or record.get('dnam_study_count') != len(studies)
            or set(record.get('required_checks', [])) != required_checks):
        raise ValueError('Independent source cross-readback plan binding is invalid')
    if (receipt.get('schema') != record.get('schema')
            or receipt.get('status') != 'PASS'
            or receipt.get('status_scope') != 'saved-source-cross-readback-only'
            or receipt.get('source_only') is not True
            or receipt.get('scientific_status') != 'CANNOT VERIFY'
            or any(receipt.get(key) != 0
                   for key in ('clock_executions', 'model_calls', 'model_api_calls'))
            or receipt.get('failures') != []
            or receipt.get('source_status_partition')
            != record.get('source_status_partition_retained')):
        raise ValueError('Independent source cross-readback scope/status is invalid')
    if (receipt.get('plan_sha256') != plan['source_origin_plan']['sha256']
            or _repository_path(receipt.get('plan_path', ''))
            != Path(plan['source_origin_plan']['path']).resolve()
            or receipt.get('root_summary_sha256') != plan['source_origin_summary']['sha256']
            or _repository_path(receipt.get('root_summary_path', ''))
            != Path(plan['source_origin_summary']['path']).resolve()):
        raise ValueError('Independent source cross-readback has stale root bindings')
    source_run = _bound_file(record['source_run'], 'independent source readback run')
    if (receipt.get('run_readback_sha256') != record['source_run']['sha256']
            or _repository_path(receipt.get('run_readback_path', '')) != source_run):
        raise ValueError('Independent source cross-readback has a stale source-run binding')

    expected = {entry['study_id']: entry for entry in studies}
    comparisons = receipt.get('dnam_comparisons')
    if (not isinstance(comparisons, list)
            or receipt.get('dnam_study_count') != len(expected)
            or len(comparisons) != len(expected)
            or {row.get('study_id') for row in comparisons} != set(expected)):
        raise ValueError('Independent source cross-readback study denominator differs from plan')
    for row in comparisons:
        entry = expected[row['study_id']]
        expected_status = entry['beta_input_contract']['status']
        checks = row.get('checks', {})
        if (row.get('all_checks_pass') is not True
                or any(checks.get(key) is not True for key in required_checks)
                or row.get('ours_source_status') != expected_status
                or row.get('root_contract_status') != expected_status):
            raise ValueError('Independent source cross-readback differs for ' + row['study_id'])

    source_receipt = json.loads(source_run.read_text())
    if (source_receipt.get('source_only') is not True
            or source_receipt.get('clock_executions') != 0
            or source_receipt.get('model_api_calls') != 0):
        raise ValueError('Independent source readback run has invalid execution controls')
    source_rows = [row for row in source_receipt.get('records', [])
                   if row.get('role') == 'series_matrix'
                   and row.get('study_id') in expected]
    if (len(source_rows) != len(expected)
            or {row.get('study_id') for row in source_rows} != set(expected)):
        raise ValueError('Independent source readback run lacks the planned DNAm studies')
    for row in source_rows:
        entry = expected[row['study_id']]
        counts = entry['numeric_counts']
        expected_status = entry['beta_input_contract']['status']
        bindings = {
            'sha256_before': entry['series_matrix']['sha256'],
            'sha256_after': entry['series_matrix']['sha256'],
            'sample_ids': entry['ordered_source_sample_ids'],
            'sample_id_count': entry['source_samples'],
            'feature_count': entry['source_features'],
            'header_identity_sha256': entry['source_header_sha256'],
            'feature_order_sha256': entry['source_feature_order_sha256'],
            'value_stream_sha256': entry['numeric_stream_sha256'],
            'finite_value_count': counts['finite'],
            'missing_value_count': counts.get('missing', 0),
            'out_of_beta_count': counts['outside_beta_range'],
            'source_status': expected_status,
        }
        if (any(row.get(key) != value for key, value in bindings.items())
                or row.get('source_mutated_during_read') is not False
                or row.get('arithmetic_status') != 'PASS'
                or row.get('scientific_status') != 'CANNOT VERIFY'):
            raise ValueError('Independent source readback values differ for ' + row['study_id'])
    return receipt_path, source_run


def load_direct_series_plan(plan_path, study_id, *, for_execution=False):
    plan_path = Path(plan_path).resolve()
    plan = json.loads(plan_path.read_text())
    if plan.get('schema') not in {DIRECT_PLAN_SCHEMA, DIRECT_PLAN_SCHEMA_V2}:
        raise ValueError('Unsupported direct-Series clock plan schema')
    if plan.get('schema') == DIRECT_PLAN_SCHEMA_V2:
        prior_path = _bound_file(plan.get('prior_plan', {}), 'held direct-Series plan v1')
        prior = json.loads(prior_path.read_text())
        if (prior.get('schema') != DIRECT_PLAN_SCHEMA
                or prior.get('status') != 'HOLD_PENDING_INDEPENDENT_SOURCE_READBACK'
                or plan['prior_plan'].get('status') != prior['status']):
            raise ValueError('Direct-Series plan v2 does not preserve its held v1 predecessor')
    studies = plan.get('studies')
    if not isinstance(studies, list) or not studies:
        raise ValueError('Direct-Series clock plan has no studies')
    ids = [entry.get('study_id') for entry in studies]
    if len(ids) != len(set(ids)) or any(re.fullmatch(r'GSE[0-9]+', str(v)) is None for v in ids):
        raise ValueError('Direct-Series study identities are missing, noncanonical or duplicated')
    counts = {
        'selected_studies': len(studies),
        'eligible_studies': sum(x.get('execution_status') == 'ELIGIBLE' for x in studies),
        'stopped_studies': sum(x.get('execution_status') == 'STOPPED' for x in studies),
        'eligible_source_samples': sum(x.get('source_samples', 0) for x in studies
                                       if x.get('execution_status') == 'ELIGIBLE'),
        'stopped_source_samples': sum(x.get('source_samples', 0) for x in studies
                                      if x.get('execution_status') == 'STOPPED'),
    }
    if plan.get('denominator') != counts:
        raise ValueError('Direct-Series plan denominator differs from its study records')
    matches = [entry for entry in studies if entry['study_id'] == study_id]
    if len(matches) != 1:
        raise ValueError('Requested study is absent or duplicated in the direct-Series plan')
    entry = matches[0]

    for name in ('source_origin_plan', 'source_origin_summary',
                 'source_origin_validation_table'):
        _bound_file(plan[name], name)
    series = _bound_file(entry['series_matrix'], study_id + ' series matrix')
    metadata = _bound_file(entry['metadata'], study_id + ' metadata')
    readback_path = _bound_file(entry['source_readback'], study_id + ' source readback')
    readback = json.loads(readback_path.read_text())
    ordered = entry.get('ordered_source_sample_ids')
    if (not isinstance(ordered, list) or len(ordered) != entry.get('source_samples')
            or len(ordered) != len(set(ordered))
            or any(re.fullmatch(r'GSM[0-9]+', str(value)) is None for value in ordered)):
        raise ValueError('Direct-Series ordered sample identities are invalid')
    if _primary_series_sample_ids(series) != ordered:
        raise ValueError('Direct-Series primary header differs from the planned identity order')
    metadata_order = _metadata_sample_ids(metadata)
    if set(metadata_order) != set(ordered):
        raise ValueError('Direct-Series metadata membership differs from the primary header')
    for key in ('study_id', 'source_samples', 'source_features', 'source_header_sha256',
                'source_feature_order_sha256', 'numeric_stream_sha256',
                'source_parser_status', 'beta_input_contract'):
        if readback.get(key) != entry.get(key):
            raise ValueError('Direct-Series source readback differs from plan: ' + key)
    if readback.get('source_sample_ids') != ordered:
        raise ValueError('Direct-Series source readback identity order differs from plan')
    if readback.get('source_sha256') != entry['series_matrix']['sha256']:
        raise ValueError('Direct-Series source readback file hash differs from plan')

    reference = plan['reference']
    reference_source = _bound_file(
        {'path': reference['source'], 'sha256': reference['source_sha256']},
        'GSE40279 source')
    cache_path = _bound_file(
        {'path': reference['cache_path'], 'sha256': reference['cache_sha256']},
        'GSE40279 versioned cache')
    reference_receipt = _bound_file(
        {'path': reference['receipt_path'], 'sha256': reference['receipt_sha256']},
        'GSE40279 reference receipt')
    receipt = json.loads(reference_receipt.read_text())
    receipt_source = receipt.get('source_record', {})
    if (reference.get('receipt_status') != 'PASS' or receipt.get('status') != 'PASS'
            or receipt.get('reader_contract') != reference.get('reader_contract')
            or receipt_source.get('source_sha256') != reference['source_sha256']
            or receipt_source.get('means_sha256') != reference['cache_sha256']
            or receipt.get('versioned_cache_sha256') != reference['cache_sha256']
            or Path(receipt.get('versioned_cache_path', '')).resolve() != cache_path):
        raise ValueError('GSE40279 reference receipt is not PASS')
    execution = plan['clock_execution']
    if execution.get('clocks') != ['Horvathv1', 'Hannum', 'PhenoAge']:
        raise ValueError('Direct-Series selected clocks differ from the approved plan')
    if execution.get('coefficient_arithmetic_tolerance') != {
            'rtol': 1e-10, 'atol': 1e-10}:
        raise ValueError('Direct-Series coefficient tolerance differs from the producer contract')
    if (execution.get('historical_prediction_join') is not False
            or execution.get('historical_prediction_comparison') is not False):
        raise ValueError('Direct-Series plan must prohibit historical prediction joins')
    limits = plan['resource_limits']
    if (limits.get('study_concurrency') != 1 or limits.get('threads_per_process') != 1
            or limits.get('automatic_retries') != 0 or limits.get('network_calls') != 0
            or limits.get('paid_model_api_calls') != 0):
        raise ValueError('Direct-Series resource limits differ from the approved plan')

    independent = plan['independent_source_readback']
    independent_path = None
    independent_source_run = None
    if independent.get('status') == 'PASS':
        independent_path, independent_source_run = _validate_independent_readback(
            plan, studies, independent)
    elif for_execution and entry.get('execution_status') == 'ELIGIBLE':
        raise ValueError('Direct-Series clocks remain on hold pending independent source readback')
    if (for_execution and entry.get('execution_status') == 'ELIGIBLE'
            and plan.get('status') != 'READY_FOR_EXECUTION'):
        raise ValueError('Direct-Series clock plan is not released for execution')
    if entry.get('execution_status') not in {'ELIGIBLE', 'STOPPED'}:
        raise ValueError('Direct-Series study execution status is invalid')
    if entry.get('execution_status') == 'ELIGIBLE':
        if (entry.get('source_parser_status') != 'PASS'
                or entry.get('beta_input_contract', {}).get('status') != 'PASS'):
            raise ValueError('Eligible direct-Series study lacks a passing source/beta contract')
    else:
        if entry.get('beta_input_contract', {}).get('status') != 'FAIL':
            raise ValueError('Stopped direct-Series study lacks its retained beta-contract failure')

    return {'plan': plan, 'plan_path': plan_path, 'plan_sha256': sha(plan_path),
            'entry': entry, 'series': series, 'metadata': metadata,
            'metadata_source_order': metadata_order,
            'metadata_reindexed_to_primary_header': metadata_order != ordered,
            'reference_source': reference_source, 'reference_cache': cache_path,
            'reference_receipt': reference_receipt,
            'independent_readback': independent_path,
            'independent_source_run': independent_source_run}


def snapshot_direct_sources(binding, directory):
    copied = {}
    for label, source, destination in (
            ('series_matrix', binding['series'], directory/'series_matrix.txt.gz'),
            ('metadata', binding['metadata'], directory/'metadata.csv')):
        if destination.exists():
            raise FileExistsError('Direct-Series snapshot destination exists: ' + str(destination))
        shutil.copyfile(source, destination)
        expected = binding['entry'][label]['sha256']
        if sha(destination) != expected or sha(source) != expected:
            raise ValueError('Direct-Series source changed during snapshot: ' + label)
        copied[label] = {'path': str(destination.resolve()), 'sha256': expected,
                         'source_path': str(source)}
    metadata_view = directory/'metadata-primary-header-order.csv'
    with Path(copied['metadata']['path']).open(newline='', encoding='utf-8') as source:
        reader = csv.DictReader(source)
        rows = {row['GSM']: row for row in reader}
        fields = reader.fieldnames
    with metadata_view.open('x', newline='', encoding='utf-8') as target:
        writer = csv.DictWriter(target, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows[gsm] for gsm in binding['entry']['ordered_source_sample_ids'])
    manifest = {'study_id': binding['entry']['study_id'], 'status': 'PASS',
                'ordered_source_sample_ids': binding['entry']['ordered_source_sample_ids'],
                'metadata_source_order': binding['metadata_source_order'],
                'metadata_reindexed_to_primary_header': binding['metadata_reindexed_to_primary_header'],
                'metadata_primary_header_order_view': {
                    'path': str(metadata_view.resolve()), 'sha256': sha(metadata_view),
                    'join_key': 'GSM'},
                'files': copied, 'originals_modified': False}
    write_new(directory/'snapshot.json', manifest)
    return manifest


def _csv_sample_header(path):
    with Path(path).open(newline='', encoding='utf-8') as handle:
        header = next(csv.reader(handle), [])
    if len(header) < 2:
        raise ValueError('Clock intermediate has no sample columns')
    return header[1:]


def run_direct_series(args):
    if not args.study or re.fullmatch(r'GSE[0-9]+', args.study) is None:
        raise ValueError('--study with a canonical GSE accession is required in direct-Series mode')
    if args.historical:
        raise ValueError('--historical is prohibited in direct-Series mode')
    binding = load_direct_series_plan(args.direct_series_plan, args.study,
                                      for_execution=not args.prepare)
    entry = binding['entry']
    prepared = {'prepared': True, 'study_id': args.study,
                'execution_status': entry['execution_status'],
                'plan_status': binding['plan']['status'], 'model_calls': 0,
                'source_samples': entry['source_samples']}
    if args.prepare:
        return prepared
    if entry['execution_status'] == 'STOPPED':
        return {'study_id': args.study, 'status': 'STOPPED', 'execution_started': False,
                'model_calls': 0, 'source_samples': entry['source_samples'],
                'reason': entry['stop_reason'], 'candidate_status': 'NOT_EXECUTED',
                'measurement_unit_applicability': 'CANNOT VERIFY',
                'biological_sample_unit_status': 'CANNOT VERIFY',
                'scientific_status': 'CANNOT VERIFY'}
    if args.reference.resolve() != binding['reference_source']:
        raise ValueError('--reference differs from the direct-Series plan')
    if args.reference_cache.resolve() != binding['reference_cache'].parent:
        raise ValueError('--reference-cache differs from the direct-Series versioned cache directory')

    code = [Path(__file__).resolve(), ROOT/'src/pipeline_inputs_2026_09_10.py',
            ROOT/'src/pipeline_runtime_2026_09_10.py',
            ROOT/'audit_system/agent_workflow_2026_09_10.py']
    package = Path(biolearn.__file__).parent
    code.extend(package/name for name in ['model.py','model_gallery.py','imputation.py',
                                          'data/Horvath1.csv','data/Hannum.csv','data/PhenoAge.csv'])
    identities = {str(path): sha(path) for path in code}
    sources = {str(path): {'sha256': sha(path), 'size_bytes': path.stat().st_size}
               for path in [binding['plan_path'], binding['series'], binding['metadata'],
                            binding['reference_source'], binding['reference_cache'],
                            binding['reference_receipt']]}
    if binding['independent_readback']:
        for path in (binding['independent_readback'], binding['independent_source_run']):
            sources[str(path)] = {'sha256': sha(path), 'size_bytes': path.stat().st_size}
    run_config = {'schema': 'direct-series-clocks-v1', 'study_id': args.study,
        'plan_sha256': binding['plan_sha256'], 'entry': entry, 'sources': sources,
        'implementation_sha256': identities,
        'versions': {name: version(name) for name in ('biolearn','numpy','pandas','scipy')},
        'python': sys.version, 'model_calls': 0,
        'reference_cache': str(args.reference_cache.resolve()),
        'clocks': ['Horvathv1','Hannum','PhenoAge'],
        'sample_selection': 'all ordered canonical direct-Series source samples',
        'historical_prediction_comparison': False,
        'metadata_reindexed_to_primary_header': binding['metadata_reindexed_to_primary_header'],
        'measurement_unit_applicability': 'CANNOT VERIFY',
        'biological_sample_unit_status': 'CANNOT VERIFY',
        'scientific_status': 'CANNOT VERIFY'}
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/'.run.lock').open('a') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('Another process owns this clock run')
        store = StageStore(args.output, run_config)
        snapshot = store.run('00_source_snapshot', lambda directory:
            snapshot_direct_sources(binding, directory))
        recovery = {'status': 'PASS', 'study_id': args.study,
                    'source_samples': entry['source_samples'],
                    'source_features': entry['source_features'],
                    'source_sample_ids': entry['ordered_source_sample_ids'],
                    'recovered_input': snapshot['files']['series_matrix']['path'],
                    'source_mode': 'direct_series'}
        processed = store.run('01_preprocessing', lambda directory:
            preprocess_dnam(Path(recovery['recovered_input']), args.reference,
                            args.reference_cache, directory))
        if _csv_sample_header(processed['path']) != entry['ordered_source_sample_ids']:
            raise ValueError('Preprocessed sample order differs from the direct-Series plan')
        clocks = store.run('02_clock_predictions', lambda directory:
            predict_dnam(processed, directory))
        result = store.run('03_validation_table', lambda directory:
            summarize(clocks, recovery, None, directory, source_mode='direct_series'))
        for path, digest in identities.items():
            if sha(path) != digest: raise ValueError('Clock implementation changed during execution: '+path)
        for path, record in sources.items():
            if sha(path) != record['sha256']: raise ValueError('Clock source changed during execution: '+path)
        return result


def run(args):
    if args.direct_series_plan:
        return run_direct_series(args)
    if args.study:
        raise ValueError('--study is only valid with --direct-series-plan')
    config_path = args.recovery.resolve()/'config.json'
    config = json.loads(config_path.read_text())
    p = config['proposal']
    inputs = p['source_files']
    recovery = recover_dnam_input(p['study_id'], inputs['failed_matrix']['path'],
        inputs['supplement']['path'], inputs['metadata']['path'], inputs['series']['path'],
        args.recovery, sample_source=inputs.get('sample_records', {}).get('path'))
    if recovery['status'] != 'PASS': raise ValueError('Source reconstruction is not verified')
    recovery['study_id'] = p['study_id']
    code = [Path(__file__).resolve(), ROOT/'src/pipeline_inputs_2026_09_10.py',
            ROOT/'src/pipeline_runtime_2026_09_10.py', ROOT/'audit_system/report_2026_09_10.py',
            ROOT/'audit_system/agent_workflow_2026_09_10.py', ROOT/'audit_system/source_packet.py']
    package = Path(biolearn.__file__).parent
    code.extend(package/name for name in ['model.py','model_gallery.py','imputation.py',
                                          'data/Horvath1.csv','data/Hannum.csv','data/PhenoAge.csv'])
    identities = {str(path): sha(path) for path in code}
    source_paths = [config_path, args.reference.resolve()]
    if args.historical: source_paths.append(args.historical.resolve())
    sources = {str(path): {'sha256': sha(path), 'size_bytes': path.stat().st_size,
                         'mtime_ns': path.stat().st_mtime_ns} for path in source_paths}
    run_config = {'schema':'source-recovered-clocks-v1', 'recovery':recovery,
        'sources':sources, 'implementation_sha256':identities,
        'versions': {name:version(name) for name in ('biolearn','numpy','pandas','scipy')},
        'python':sys.version, 'model_calls':0, 'reference_cache':str(args.reference_cache.resolve()),
        'clocks':['Horvathv1','Hannum','PhenoAge'], 'sample_selection':'all recovered source samples',
        'preprocessing':'Existing GSE40279 population mean imputation; independently checked',
        'idat_processing':False, 'held_out':False}
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output/'.run.lock').open('a') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('Another process owns this clock run')
        store = StageStore(args.output, run_config)
        if args.prepare: return {'prepared':True, 'study_id':p['study_id'], 'model_calls':0}
        processed = store.run('01_preprocessing', lambda directory:
            preprocess_dnam(Path(recovery['recovered_input']), args.reference,
                            args.reference_cache, directory))
        clocks = store.run('02_clock_predictions', lambda directory: predict_dnam(processed, directory))
        result = store.run('03_validation_table', lambda directory:
            summarize(clocks, recovery, args.historical, directory))
        for path, digest in identities.items():
            if sha(path) != digest: raise ValueError('Clock implementation changed during execution: '+path)
        for path, record in sources.items():
            if sha(path) != record['sha256']: raise ValueError('Clock source changed during execution: '+path)
        return result


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--recovery', type=Path)
    source.add_argument('--direct-series-plan', type=Path)
    parser.add_argument('--study')
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--reference-cache', type=Path, default=ROOT/'.clockbase-cache/dnam-reference')
    parser.add_argument('--historical', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--prepare', action='store_true')
    return parser.parse_args(argv)


if __name__ == '__main__':
    args = parse_args()
    try:
        print(json.dumps(run(args), indent=2, allow_nan=False), flush=True)
    except Exception as error:
        print(f'Clock recovery validation failed: {type(error).__name__}: {error}', file=sys.stderr, flush=True)
        raise
