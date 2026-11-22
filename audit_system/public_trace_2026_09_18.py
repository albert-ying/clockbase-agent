"""Export source, molecular and numerical evidence from a completed pipeline.

Preparation pins archived receipts before export. Export verifies their payloads
and independently reads sample predictions; it never executes archived code,
calls a model, reruns a clock or upgrades scientific acceptance.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import gzip
import json
import math
from pathlib import Path
import re
import shutil
import sys
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import scipy
from scipy import stats

from .agent_workflow_2026_09_10 import AUDIT_DOMAINS, sha, source_packet, write_new
from .source_packet import source_grouping_diagnostics
from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError


STAGES = ('01_metadata', '02_molecular_input', '03_blinded_design',
          '04_preprocessing', '05_clock_predictions', '06_coding_agent',
          '07_audit_and_repair', '08_report')
FAILED_STAGE = '03_blinded_design'
STATUS = {'PASS', 'FAIL', 'CANNOT VERIFY'}
EFFECT = 'cohen_d_pooled_sample_sd_case_minus_control'


def _export_identity():
    root = Path(__file__).resolve().parents[1]
    files = ('audit_system/public_trace_2026_09_18.py',
             'audit_system/agent_workflow_2026_09_10.py',
             'audit_system/source_packet.py', 'src/pipeline_runtime_2026_09_10.py',
             'src/pipeline_2026_09_10.py')
    return {'implementation_sha256': {name: sha(root/name) for name in files},
            'python': sys.version, 'numpy': np.__version__,
            'pandas': pd.__version__, 'scipy': scipy.__version__}


def _digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def _inside(root, path):
    root, path = Path(root).resolve(), Path(path).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        raise PipelineIntegrityError('Trace input is missing or outside its recorded run')
    return path


def _snapshot(run):
    run = Path(run).resolve()
    paths = [run/'config.json']
    for name in STAGES:
        receipt = run/name/'complete.json'
        if not receipt.is_file():
            raise ValueError('Full pipeline trace requires a completed stage: ' + name)
        record = json.loads(receipt.read_text())
        sidecar = receipt.with_suffix('.sha256')
        if sidecar.exists() and sidecar.read_text().strip() != sha(receipt):
            raise PipelineIntegrityError('Stage receipt hash changed: ' + name)
        paths.append(receipt)
        if sidecar.exists(): paths.append(sidecar)
        for relative, expected in record['outputs'].items():
            artifact = _inside(run/name, run/name/relative)
            if sha(artifact) != expected:
                raise PipelineIntegrityError('Stage payload hash changed: ' + name + '/' + relative)
            paths.append(artifact)
    paths.extend(p for p in (run/'api').rglob('*') if p.is_file())
    return {str(p.relative_to(run)): sha(p) for p in sorted(set(paths))}


def prepare(run, plan_path):
    run, plan_path = Path(run).resolve(), Path(plan_path).resolve()
    snapshot = _snapshot(run)
    config = json.loads((run/'config.json').read_text())
    if not re.fullmatch(r'GSE\d+', config.get('study', '')):
        raise ValueError('A public trace requires a GEO series identifier')
    plan = {'schema': 'clockbase-pipeline-trace-plan-v1', 'study_id': config['study'],
            'run': str(run), 'snapshot_sha256': snapshot, 'split': 'development',
            'export_identity': _export_identity(), 'max_download_bytes': 50_000_000,
            'scope': 'Archived pipeline evidence and conditional numerical readback; no clock rerun'}
    write_new(plan_path, plan)
    return {'status': 'PASS', 'files_pinned': len(snapshot), 'plan_sha256': sha(plan_path)}


def _positive_int(value, name):
    if type(value) is not int or value <= 0:
        raise ValueError(name + ' must be a positive integer')
    return value


def _finite_nonnegative(value, name, *, allow_none=False):
    if value is None and allow_none:
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(name + ' must be a finite nonnegative number')
    if not math.isfinite(float(value)) or value < 0:
        raise ValueError(name + ' must be a finite nonnegative number')


def _archive_prompt_cap(config):
    cap = config.get('max_prompt_bytes', 48000)
    return _positive_int(cap, 'archived max_prompt_bytes')


def _failed_field_reference_rejection(run, config, packet, transport, failure, slots):
    """Replay a zero-retry rejection with trusted current code, never archived code.

    An exception during materialization has no saved candidate. An exception from
    strict parsing has the exact candidate and complete materialization record.
    Neither shape supplies a historical response-validation record or retry.
    """
    from .agent_workflow_2026_09_10 import parse_design, review_messages
    from src.pipeline_2026_09_10 import (SOURCE_FIELD_REFERENCE_LIMITS,
        SOURCE_FIELD_REFERENCE_POLICY, source_contract_candidate,
        source_contract_instructions, source_contract_visible_packet)

    def unique(pairs):
        value = {}
        for key, child in pairs:
            if key in value:
                raise ValueError('Duplicate archived JSON key: ' + key)
            value[key] = child
        return value

    def invalid_constant(value):
        raise ValueError('Nonfinite archived JSON value: ' + value)

    def read(path):
        try:
            return json.loads(path.read_text(), object_pairs_hook=unique,
                              parse_constant=invalid_constant)
        except OSError as error:
            raise PipelineIntegrityError('Cannot read retained field-reference artifact') from error

    if (config.get('source_contract_policy') != SOURCE_FIELD_REFERENCE_POLICY
            or type(config.get('audit_validation_retries')) is not int
            or config['audit_validation_retries'] != 0
            or config.get('audit_profile') not in {'source_domains_v3', 'source_domains_v4'}
            or len(slots) != 1):
        raise ValueError('Field-reference failure requires one response and zero validation retries')
    attempt = run / FAILED_STAGE / 'attempt-000'
    adjustment_dir = attempt / 'contract-adjustments'
    slot = slots[0]
    candidate_path = adjustment_dir / '000-candidate.json'
    expected_files = {'config.json'}
    completed = {}
    for name in STAGES[:2]:
        completed[name] = read(run / name / 'complete.json')
        expected_files.update(name + '/' + item for item in
                              ['complete.json', 'complete.sha256', *completed[name]['outputs']])
    expected_files.update(FAILED_STAGE + '/attempt-000/' + name for name in
        ['failure.json', 'source-packet.json', 'source-transport.json',
         'contract-adjustments/000-original-response.txt', 'contract-adjustments/000-receipt.json'])
    if candidate_path.exists():
        expected_files.add(str(candidate_path.relative_to(run)))
    expected_files.update('api/design-000/' + name for name in
                          ['request.json', 'reservation.json', 'response.json', 'receipt.json'])
    if (run / '.run.lock').exists():
        expected_files.add('.run.lock')
    expected_directories = {str(parent) for name in expected_files
                            for parent in Path(name).parents if str(parent) != '.'}
    observed_files, observed_directories = set(), set()
    for path in run.rglob('*'):
        if path.is_symlink():
            raise PipelineIntegrityError('Field-reference failure cannot contain symlinks')
        relative = str(path.relative_to(run))
        if path.is_file():
            observed_files.add(relative)
        elif path.is_dir():
            observed_directories.add(relative)
        else:
            raise PipelineIntegrityError('Field-reference failure requires regular files')
    if observed_files != expected_files or observed_directories != expected_directories:
        raise ValueError('Field-reference failure contains missing, extra or mixed artifacts')
    # Read strict JSON again so duplicates cannot hide behind the legacy reader.
    for path, value in [(run/'config.json', config), (attempt/'source-packet.json', packet),
                        (attempt/'source-transport.json', transport), (attempt/'failure.json', failure)]:
        if _digest(read(path)) != _digest(value):
            raise PipelineIntegrityError('Field-reference control artifact changed during readback')
    dependencies = ('src/pipeline_2026_09_10.py', 'audit_system/agent_workflow_2026_09_10.py',
                    'audit_system/source_packet.py', 'src/pipeline_runtime_2026_09_10.py')
    current = _export_identity()['implementation_sha256']
    if (not isinstance(config.get('implementation'), dict)
            or any(config['implementation'].get(name) != current[name] for name in dependencies)):
        raise PipelineIntegrityError('Field-reference replay implementation differs from the recorded producer')
    metadata_stage = run / STAGES[0]
    metadata_result = completed[STAGES[0]].get('result', {})
    def source_role(name, required=False):
        value = metadata_result.get(name)
        if value is None and not required:
            return None
        if not isinstance(value, str) or not value:
            raise ValueError('Field-reference metadata source role is missing: ' + name)
        path = _inside(metadata_stage, Path(value) if Path(value).is_absolute() else metadata_stage / value)
        if str(path.relative_to(metadata_stage)) not in completed[STAGES[0]]['outputs']:
            raise PipelineIntegrityError('Field-reference source role is not a completed-stage output')
        return path
    if not isinstance(config.get('target'), str) or not config['target'].strip():
        raise ValueError('Field-reference source target is missing')
    reconstructed_packet = source_packet(source_role('metadata', required=True), config['study'], config['target'],
        context=source_role('study_context'), sample_source=source_role('source_sample_records'),
        comparison_scope=config.get('comparison_scope', 'exact_target'))
    if _digest(reconstructed_packet) != _digest(packet):
        raise PipelineIntegrityError('Field-reference source packet differs from completed metadata/context')
    messages = review_messages(packet, mode='independent_audit', primary=None,
        review_profile=config['audit_profile'], source_binding_profile=config.get('source_binding_profile', 'legacy'))
    wire = [{'role': 'system', 'content': messages[0].content + source_contract_instructions(SOURCE_FIELD_REFERENCE_POLICY)},
            {'role': 'user', 'content': messages[1].content}]
    if _digest(transport) != _digest(json.loads(wire[1]['content'])):
        raise PipelineIntegrityError('Field-reference source transport differs from the current source request')
    visible = source_contract_visible_packet(packet, transport)
    model = config.get('audit_model') or config.get('model')
    if not isinstance(model, str) or not model:
        raise ValueError('Field-reference configured source model is missing')
    expected_request = {'model': model, 'messages': wire,
        'max_completion_tokens': _positive_int(config.get('max_output_tokens'), 'max_output_tokens'),
        'reasoning_effort': 'low', 'store': False, 'response_format': {'type': 'json_object'}}
    request, response, receipt, reservation = [read(slot/name) for name in
        ['request.json', 'response.json', 'receipt.json', 'reservation.json']]
    if _digest(request) != _digest(expected_request):
        raise PipelineIntegrityError('Field-reference request differs from the configured source request')
    if response['choices'][0].get('finish_reason') != 'stop':
        raise ValueError('Field-reference rejection requires a completed generation')
    raw = response['choices'][0]['message']['content']
    original_path = adjustment_dir / '000-original-response.txt'
    if original_path.read_bytes() != raw.encode('utf-8'):
        raise PipelineIntegrityError('Field-reference original bytes differ from the API response')
    if (set(failure) != {'status', 'error', 'elapsed_seconds'}
            or set(receipt) != {'usage', 'estimated_cost_usd', 'billed_cost_usd', 'response_sha256',
                                'requested_model', 'returned_model', 'elapsed_seconds'}
            or set(reservation) != {'reserved_usd', 'started_unix', 'model', 'automatic_retries', 'concurrency'}
            or _digest(reservation.get('automatic_retries')) != _digest(0)
            or _digest(reservation.get('concurrency')) != _digest(1)
            or receipt.get('billed_cost_usd') is not None
            or _digest(receipt.get('usage')) != _digest(response.get('usage'))):
        raise ValueError('Field-reference failure, API usage or reservation differs from the supported shape')
    usage = receipt.get('usage')
    if (not isinstance(usage, dict)
            or any(type(usage.get(key)) is not int or usage[key] < 0 for key in ('prompt_tokens', 'completion_tokens'))):
        raise ValueError('Field-reference response lacks recorded token usage')
    rate_prefix = 'audit_' if config.get('audit_model') else ''
    rates = [config.get(rate_prefix + direction + '_usd_per_million') for direction in ('input', 'output')]
    for rate in rates:
        _finite_nonnegative(rate, 'configured token rate')
        if rate == 0:
            raise ValueError('Configured token rate must be positive')
    request_bytes = len(json.dumps(request, ensure_ascii=False).encode())
    estimated = (usage['prompt_tokens'] * rates[0] + usage['completion_tokens'] * rates[1]) / 1e6
    reserved = ((request_bytes + 100 * len(wire)) * rates[0] + config['max_output_tokens'] * rates[1]) / 1e6
    if receipt['estimated_cost_usd'] != estimated or reservation['reserved_usd'] != reserved:
        raise PipelineIntegrityError('Field-reference recorded estimate or reservation differs from frozen rates')
    adjustment = read(adjustment_dir / '000-receipt.json')
    expected = {'schema': 'clockbase-source-contract-adjustment-v1', 'policy': SOURCE_FIELD_REFERENCE_POLICY,
        'attempt': 0, 'original_response_sha256': sha(original_path), 'changes': [],
        'scientific_status': 'CANNOT VERIFY', 'accepted_result': False,
        'scope': 'Exact selected-field materialization and conservative advisory statuses; strict evidence and scientific checks still apply',
        'source_packet_sha256': sha(attempt/'source-packet.json'),
        'source_transport_sha256': sha(attempt/'source-transport.json'),
        'implementation_sha256': current['src/pipeline_2026_09_10.py'],
        'normalization_limits': SOURCE_FIELD_REFERENCE_LIMITS}
    rejection_step = 'materialization'
    try:
        value, changes = source_contract_candidate(raw, visible, SOURCE_FIELD_REFERENCE_POLICY)
    except (ValueError, TypeError, KeyError) as error:
        rejection_error = f'{type(error).__name__}: {error}'
        if candidate_path.exists():
            raise PipelineIntegrityError('Rejected materialization must not retain a candidate') from error
    else:
        rejection_step = 'strict_parse'
        expected_bytes = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode('utf-8')
        if not candidate_path.is_file() or candidate_path.read_bytes() != expected_bytes:
            raise PipelineIntegrityError('Field-reference candidate differs from exact current materialization')
        expected.update(changes=changes, candidate_sha256=sha(candidate_path))
        try:
            parse_design(json.dumps(value, ensure_ascii=False, allow_nan=False), packet,
                         primary=None, review_profile=config['audit_profile'])
        except (ValueError, TypeError, KeyError) as error:
            rejection_error = f'{type(error).__name__}: {error}'
        else:
            raise PipelineIntegrityError('Currently valid field-reference candidate was labeled rejected')
    expected.update(contract_status='REJECTED', error=rejection_error)
    if _digest(adjustment) != _digest(expected) or failure['error'] != rejection_error:
        raise PipelineIntegrityError('Recorded field-reference rejection differs from trusted replay')
    return {'rejection_step': rejection_step, 'error': rejection_error,
            'adjustment': adjustment, 'adjustment_sha256': sha(adjustment_dir/'000-receipt.json'),
            'current_implementation_sha256': {name: current[name] for name in dependencies}}


def _failed_stage_snapshot(run):
    """Pin the deliberately narrow pre-analysis failure profile.

    This profile accepts complete metadata and molecular-input stages followed by
    a failed blinded-design stage.  It never treats a later stage as merely
    unavailable: any later receipt, output, or completed stage is rejected.
    """
    run = Path(run).resolve()
    config_path = run / 'config.json'
    if not config_path.is_file():
        raise ValueError('Failed trace requires config.json')
    allowed_top_level = {'config.json', '01_metadata', '02_molecular_input', '03_blinded_design', 'api', '.run.lock'}
    if any(path.name not in allowed_top_level for path in run.iterdir()):
        raise ValueError('Failed trace contains unsupported top-level artifacts')
    config = json.loads(config_path.read_text())
    archived_cap = _archive_prompt_cap(config)
    study = config.get('study')
    if not isinstance(study, str) or not re.fullmatch(r'GSE\d+', study):
        raise ValueError('A failed public trace requires a GEO series identifier')
    snapshot = {str(config_path.relative_to(run)): sha(config_path)}
    completed = {}
    for name in STAGES[:2]:
        stage_dir = run / name
        receipt_path = stage_dir / 'complete.json'
        sidecar = stage_dir / 'complete.sha256'
        if not receipt_path.is_file() or not sidecar.is_file():
            raise ValueError('Failed trace requires completed ' + name)
        if sidecar.read_text().strip() != sha(receipt_path):
            raise PipelineIntegrityError('Completed stage receipt hash changed: ' + name)
        receipt = json.loads(receipt_path.read_text())
        completed[name] = receipt
        outputs = receipt.get('outputs')
        if not isinstance(outputs, dict) or not outputs:
            raise ValueError('Completed stage lacks a pinned output manifest: ' + name)
        snapshot[str(receipt_path.relative_to(run))] = sha(receipt_path)
        snapshot[str(sidecar.relative_to(run))] = sha(sidecar)
        for relative, expected in outputs.items():
            artifact = _inside(stage_dir, stage_dir / relative)
            if not isinstance(expected, str) or sha(artifact) != expected:
                raise PipelineIntegrityError('Completed stage payload hash changed: ' + name + '/' + relative)
            snapshot[str(artifact.relative_to(run))] = expected
    design_dir = run / FAILED_STAGE
    if (not design_dir.is_dir() or (design_dir / 'complete.json').exists()
            or (design_dir / 'complete.sha256').exists()):
        raise ValueError('Failed trace requires an incomplete blinded-design stage')
    attempts = sorted(path for path in design_dir.glob('attempt-*') if path.is_dir())
    if len(attempts) != 1 or attempts[0].name != 'attempt-000':
        raise ValueError('Failed trace supports exactly one attempt-000 blinded-design failure')
    attempt = attempts[0]
    failure_path = attempt / 'failure.json'
    packet_path, transport_path = attempt / 'source-packet.json', attempt / 'source-transport.json'
    if not failure_path.is_file() or not packet_path.is_file() or not transport_path.is_file():
        raise ValueError('Failed blinded-design attempt lacks its failure and source artifacts')
    failure = json.loads(failure_path.read_text())
    if (failure.get('status') != 'CANNOT VERIFY' or not isinstance(failure.get('error'), str)
            or not failure['error'].strip() or type(failure.get('elapsed_seconds')) not in (int, float)
            or not math.isfinite(failure['elapsed_seconds']) or failure['elapsed_seconds'] < 0):
        raise ValueError('Failed blinded-design attempt lacks a valid recorded failure')
    packet = json.loads(packet_path.read_text())
    transport = json.loads(transport_path.read_text())
    if packet.get('study_id') != study or transport.get('study_id') != study:
        raise ValueError('Failed source packet study identity changed')
    if packet.get('source_sha256') != transport.get('source_sha256'):
        raise ValueError('Failed source packet and transport source identities differ')
    metadata_result = completed['01_metadata'].get('result', {})
    metadata_names = [name for name in completed['01_metadata']['outputs']
                      if Path(name).name == Path(str(metadata_result.get('metadata', ''))).name]
    if len(metadata_names) != 1 or packet.get('source_sha256') != completed['01_metadata']['outputs'][metadata_names[0]]:
        raise PipelineIntegrityError('Failed source packet metadata identity differs from completed metadata stage')
    declared_ids = metadata_result.get('sample_ids')
    packet_ids = [row.get('sample_id') for row in packet.get('samples', [])]
    if (not isinstance(declared_ids, list) or packet_ids != declared_ids
            or len(set(packet_ids)) != len(packet_ids)):
        raise ValueError('Failed source packet sample roster differs from completed metadata stage')
    for path in sorted(attempt.rglob('*')):
        if path.is_file():
            snapshot[str(path.relative_to(run))] = sha(path)
    api_dir = run / 'api'
    if not api_dir.is_dir():
        raise ValueError('Failed trace requires a pinned API artifact directory')
    slots = sorted(path for path in api_dir.glob('design-*') if path.is_dir())
    if not slots or len(slots) > 2 or [p.name for p in slots] != [f'design-{i:03d}' for i in range(len(slots))]:
        raise ValueError('Failed trace supports one initial call and at most one validation retry')
    if any(not path.is_dir() or path not in slots for path in api_dir.iterdir()):
        raise ValueError('Failed trace API directory contains unsupported artifacts')
    for slot in slots:
        expected_files = {'request.json', 'reservation.json', 'response.json', 'receipt.json'}
        if {p.name for p in slot.iterdir() if p.is_file()} != expected_files:
            raise ValueError('API failure slot has unsupported or missing artifacts: ' + slot.name)
        request = json.loads((slot / 'request.json').read_text())
        response = json.loads((slot / 'response.json').read_text())
        receipt = json.loads((slot / 'receipt.json').read_text())
        reservation = json.loads((slot / 'reservation.json').read_text())
        if not isinstance(request.get('messages'), list) or not isinstance(request.get('model'), str):
            raise ValueError('API failure request is malformed')
        if receipt.get('response_sha256') != sha(slot / 'response.json'):
            raise PipelineIntegrityError('API response hash mismatch: ' + slot.name)
        if (receipt.get('requested_model') != request['model'] or receipt.get('returned_model') != response.get('model')
                or reservation.get('model') != request['model'] or response.get('model') != request['model']):
            raise ValueError('API failure model identity changed: ' + slot.name)
        if (len(response.get('choices', [])) != 1 or not isinstance(response['choices'][0], dict)
                or not isinstance(response['choices'][0].get('message', {}).get('content'), str)):
            raise ValueError('API failure response lacks one retained message')
        request_bytes = len(json.dumps(request, ensure_ascii=False).encode())
        if request_bytes > archived_cap:
            raise ValueError('Sent API request exceeds the archived prompt byte limit')
        _finite_nonnegative(receipt.get('estimated_cost_usd'), 'estimated_cost_usd')
        _finite_nonnegative(receipt.get('billed_cost_usd'), 'billed_cost_usd', allow_none=True)
        _finite_nonnegative(receipt.get('elapsed_seconds'), 'elapsed_seconds')
        _finite_nonnegative(reservation.get('reserved_usd'), 'reserved_usd')
        _finite_nonnegative(reservation.get('started_unix'), 'started_unix')
        for path in sorted(slot.iterdir()):
            if path.is_file(): snapshot[str(path.relative_to(run))] = sha(path)
    if (config.get('source_contract_policy') == 'literal_field_refs_v1'
            and config.get('audit_validation_retries') == 0):
        _failed_field_reference_rejection(run, config, packet, transport, failure, slots)
        return snapshot, config, packet, failure, slots, []
    validation_dir = attempt / 'response-validation'
    validations = sorted(validation_dir.glob('*-validation.json')) if validation_dir.is_dir() else []
    response_texts = sorted(validation_dir.glob('*-response.txt')) if validation_dir.is_dir() else []
    if len(validations) != len(slots) or len(response_texts) != len(slots):
        raise ValueError('Every retained API response needs one response and validation artifact')
    for index, path in enumerate(validations):
        if path.name != f'{index:03d}-validation.json' or response_texts[index].name != f'{index:03d}-response.txt':
            raise ValueError('Validation artifacts are not contiguous')
        feedback = json.loads(path.read_text())
        if feedback.get('validation_status') != 'FAIL' or not isinstance(feedback.get('instruction'), str):
            raise ValueError('Failed trace requires a failed response-validation record')
        issues = feedback.get('citation_issues', [])
        omitted = feedback.get('citation_issues_omitted', 0)
        if (not isinstance(issues, list) or type(feedback.get('citation_issue_count')) is not int
                or type(omitted) is not int or omitted < 0
                or feedback['citation_issue_count'] != len(issues) + omitted):
            raise ValueError('Response-validation citation issue counts are malformed')
        _failed_coverage_details(feedback, packet)
        response = json.loads((slots[index] / 'response.json').read_text())
        if response_texts[index].read_text() != response['choices'][0]['message']['content']:
            raise PipelineIntegrityError('Saved response text differs from API response: ' + str(index))
        snapshot[str(path.relative_to(run))] = sha(path)
        snapshot[str(response_texts[index].relative_to(run))] = sha(response_texts[index])
    for index in range(1, len(slots)):
        previous_request = json.loads((slots[index - 1] / 'request.json').read_text())
        previous_response = json.loads((slots[index - 1] / 'response.json').read_text())
        previous_feedback = json.loads(validations[index - 1].read_text())
        request = json.loads((slots[index] / 'request.json').read_text())
        expected = _failed_retry_request(previous_request, previous_response, previous_feedback)
        if request != expected:
            raise PipelineIntegrityError('Saved validation retry request differs from its retained feedback: ' + slots[index].name)
    for name in STAGES[3:]:
        directory = run / name
        if directory.exists():
            raise ValueError('Failed trace cannot include downstream stage artifacts: ' + name)
    return snapshot, config, packet, failure, slots, validations


def prepare_failed(run, plan_path, prompt_bytes=None):
    """Prepare an immutable plan for a pre-analysis blinded-design failure."""
    run, plan_path = Path(run).resolve(), Path(plan_path).resolve()
    snapshot, config, packet, failure, slots, validations = _failed_stage_snapshot(run)
    archived_cap = _archive_prompt_cap(config)
    cap = archived_cap if prompt_bytes is None else _positive_int(prompt_bytes, 'prompt_bytes')
    if cap != archived_cap:
        raise ValueError('prompt_bytes must equal the archived prompt byte limit')
    plan = {'schema': 'clockbase-failed-pipeline-trace-plan-v1', 'study_id': config['study'],
            'run': str(run), 'split': 'development', 'failure_stage': FAILED_STAGE,
            'failure_attempt': 'attempt-000', 'prompt_byte_limit': cap,
            'snapshot_sha256': snapshot, 'export_identity': _export_identity(),
            'max_download_bytes': 50_000_000,
            'scope': 'Recorded metadata and molecular-input stages followed by a failed blinded-design response; no downstream execution'}
    write_new(plan_path, plan)
    return {'status': 'PASS', 'files_pinned': len(snapshot), 'api_responses': len(slots),
            'validation_records': len(validations), 'plan_sha256': sha(plan_path)}


def _failed_public_sources(packet):
    documents = packet.get('study_documents', [])
    source_url = _public_url(packet['source_url'])
    public_documents = []
    for doc in documents:
        public_documents.append({'source_id': doc['source_id'], 'source_url': _public_url(doc['source_url']),
                                 'source_sha256': doc['source_sha256'], 'fields': doc['fields']})
    return {'schema': 'clockbase-public-source-records-v1', 'study_id': packet['study_id'],
            'source_metadata_sha256': packet['source_sha256'], 'source_url': source_url,
            'samples': packet['samples'], 'study_documents': public_documents}


def _failed_retry_request(request, response, feedback):
    retry = dict(request)
    retry['messages'] = list(request['messages']) + [
        {'role': 'assistant', 'content': response['choices'][0]['message']['content']},
        {'role': 'user', 'content': json.dumps(feedback, ensure_ascii=False)}]
    return retry


def _failed_coverage_details(feedback, packet):
    """Expose bounded v3 diagnostics without changing historical v1/v2 exports."""
    if feedback.get('schema') != 'clockbase-design-validation-feedback-v3':
        return {}
    issues = feedback.get('coverage_issues')
    count, omitted = feedback.get('coverage_issue_count'), feedback.get('coverage_issues_omitted')
    if (not isinstance(issues, list) or type(count) is not int or count < 0
            or type(omitted) is not int or len(issues) != min(count, 20)
            or omitted != max(0, count - 20)):
        raise ValueError('Response-validation coverage issue counts are malformed')
    source_ids = {row['sample_id'] for row in packet['samples']}
    source_fields = {field for row in packet['samples'] for field, value in row.items()
                     if field != 'sample_id' and isinstance(value, str)}
    index = r'(?:0|[1-9][0-9]*)'
    comparison = r'(?:\.comparisons\[' + index + r'\])?'
    kinds = {'COVARIATE_FIELD_QUOTATION_MISSING': ('missing_source_fields',
                 r'\$\.covariate_plan\[' + index + r'\]\.source_fields'),
             'CONDITION_SAMPLE_QUOTATION_MISSING': ('missing_sample_ids',
                 r'\$' + comparison + r'\.findings\[' + index + r'\]\.condition'),
             'ASSIGNED_SAMPLE_QUOTATION_MISSING': ('missing_sample_ids',
                 r'\$' + comparison + r'\.(?:case_ids|control_ids)')}
    for item in issues:
        field, pattern = (kinds.get(item.get('kind'), (None, None))
                          if isinstance(item, dict) and isinstance(item.get('kind'), str) else (None, None))
        if (field is None or set(item) != {'path', 'kind', field, 'missing_count', 'omitted_count'}
                or not isinstance(item['path'], str) or re.fullmatch(pattern, item['path']) is None
                or not isinstance(item[field], list)
                or any(not isinstance(value, str) or not value for value in item[field])
                or any(value not in (source_ids if field == 'missing_sample_ids' else source_fields) for value in item[field])
                or item[field] != sorted(set(item[field]))
                or type(item['missing_count']) is not int or item['missing_count'] < 1
                or type(item['omitted_count']) is not int
                or len(item[field]) != min(item['missing_count'], 20)
                or item['omitted_count'] != max(0, item['missing_count'] - 20)):
            raise ValueError('Response-validation coverage issue details are malformed')
    return {'validation_feedback_schema': feedback['schema'], 'coverage_issue_count': count,
            'coverage_issues': issues, 'coverage_issues_omitted': omitted}


def _failed_history(run, plan, config, packet, failure, slots, validations, snapshot):
    if not validations:
        attempt = run / FAILED_STAGE / 'attempt-000'
        transport = json.loads((attempt / 'source-transport.json').read_text())
        replay = _failed_field_reference_rejection(run, config, packet, transport, failure, slots)
        slot = slots[0]
        request = json.loads((slot / 'request.json').read_text())
        receipt = json.loads((slot / 'receipt.json').read_text())
        adjustment = replay['adjustment']
        recorded = {'schema': adjustment['schema'], 'policy': adjustment['policy'],
            'contract_status': adjustment['contract_status'], 'error': adjustment['error'],
            'receipt_sha256': replay['adjustment_sha256'],
            'original_response_sha256': adjustment['original_response_sha256'],
            'candidate_sha256': adjustment.get('candidate_sha256'),
            'source_packet_sha256': adjustment['source_packet_sha256'],
            'source_transport_sha256': adjustment['source_transport_sha256'],
            'implementation_sha256': adjustment['implementation_sha256'],
            'recorded_change_count': len(adjustment['changes'])}
        diagnostic = {'status': 'REPRODUCED_REJECTION', 'rejection_step': replay['rejection_step'],
            'error': replay['error'], 'current_implementation_sha256': replay['current_implementation_sha256'],
            'scope': 'Export-time replay of the trusted materializer and strict parser; not a recorded response-validation artifact, transmitted feedback, or scientific assessment'}
        return {'schema': 'clockbase-failure-history-v1', 'study_id': packet['study_id'],
            'failure_profile': 'zero_retry_literal_field_reference_rejection_v1',
            'failure_stage': FAILED_STAGE, 'failure_attempt': 'attempt-000',
            'failure': {'status': failure['status'], 'error': failure['error'],
                        'elapsed_seconds': failure['elapsed_seconds']},
            'prompt_byte_limit': plan['prompt_byte_limit'],
            'attempts': [{'index': 0, 'api_slot': slot.name,
                'request_bytes': len(json.dumps(request, ensure_ascii=False).encode()),
                'request_sha256': sha(slot/'request.json'), 'prompt_byte_limit': plan['prompt_byte_limit'],
                'prompt_limit_check': 'WITHIN_LIMIT', 'response_sha256': receipt['response_sha256'],
                'requested_model': receipt['requested_model'], 'returned_model': receipt['returned_model'],
                'estimated_cost_usd': receipt['estimated_cost_usd'], 'billed_cost_usd': receipt.get('billed_cost_usd'),
                'elapsed_seconds': receipt['elapsed_seconds'], 'recorded_rejection': recorded,
                'export_time_diagnostic': diagnostic, 'retry_transmitted_after_this_response': False}],
            'automatic_validation_retries': 0, 'response_validation_records': 0,
            'next_retry': None, 'source_sha256': dict(sorted(snapshot.items())),
            'source_links': [{'source_id': doc['source_id'], 'source_url': _public_url(doc['source_url'])}
                             for doc in packet.get('study_documents', [])],
            'scope': 'Recorded contract-adjustment rejection with separate export-time diagnostic; raw response content and scientific model assessments are not accepted or published'}
    attempts, retry_request = [], None
    for index, slot in enumerate(slots):
        request = json.loads((slot / 'request.json').read_text())
        response = json.loads((slot / 'response.json').read_text())
        receipt = json.loads((slot / 'receipt.json').read_text())
        feedback = json.loads(validations[index].read_text())
        serialized_bytes = len(json.dumps(request, ensure_ascii=False).encode())
        issues = feedback.get('citation_issues', [])
        attempts.append({'index': index, 'api_slot': slot.name, 'request_bytes': serialized_bytes,
                         'prompt_byte_limit': plan['prompt_byte_limit'],
                         'prompt_limit_check': 'WITHIN_LIMIT' if serialized_bytes <= plan['prompt_byte_limit'] else 'EXCEEDS_LIMIT',
                         'response_sha256': receipt['response_sha256'],
                         'requested_model': receipt['requested_model'], 'returned_model': receipt['returned_model'],
                         'estimated_cost_usd': receipt['estimated_cost_usd'],
                         'billed_cost_usd': receipt.get('billed_cost_usd'),
                         'elapsed_seconds': receipt['elapsed_seconds'],
                         'validation_status': feedback['validation_status'],
                         'validation_error': feedback.get('error'),
                         'citation_issue_count': feedback['citation_issue_count'],
                         'citation_issues': issues,
                         'citation_issues_omitted': feedback.get('citation_issues_omitted', 0),
                         'retry_transmitted_after_this_response': index + 1 < len(slots)})
        attempts[-1].update(_failed_coverage_details(feedback, packet))
        if index + 1 == len(slots):
            retry_request = _failed_retry_request(request, response, feedback)
    next_retry = None
    if retry_request is not None and len(slots) == 1:
        retry_bytes = len(json.dumps(retry_request, ensure_ascii=False).encode())
        next_retry = {'request_bytes': retry_bytes,
                      'prompt_byte_limit': plan['prompt_byte_limit'],
                      'transmitted': False,
                      'reason': ('reconstructed_retry_exceeds_prompt_byte_limit'
                                 if retry_bytes > plan['prompt_byte_limit']
                                 else 'validation_feedback_saved_without_a_second_api_slot')}
    public_hashes = {relative: digest for relative, digest in sorted(snapshot.items())}
    return {'schema': 'clockbase-failure-history-v1', 'study_id': packet['study_id'],
            'failure_stage': FAILED_STAGE, 'failure_attempt': 'attempt-000',
            'failure': {'status': failure['status'], 'error': failure['error'],
                        'elapsed_seconds': failure['elapsed_seconds']},
            'prompt_byte_limit': plan['prompt_byte_limit'], 'attempts': attempts,
            'next_retry': next_retry, 'source_sha256': public_hashes,
            'source_links': [{'source_id': doc['source_id'], 'source_url': _public_url(doc['source_url'])}
                             for doc in packet.get('study_documents', [])],
            'scope': 'Recorded response-contract failure; response content and scientific model assessments are not accepted'}


def export_failed(plan_path, output):
    """Export the bounded pre-analysis failure as a source-only v1 trace."""
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = json.loads(plan_path.read_text())
    if plan.get('schema') != 'clockbase-failed-pipeline-trace-plan-v1' or plan.get('split') != 'development':
        raise ValueError('Unsupported failed pipeline trace plan')
    run = Path(plan['run']).resolve()
    snapshot, config, packet, failure, slots, validations = _failed_stage_snapshot(run)
    if (snapshot != plan.get('snapshot_sha256') or _export_identity() != plan.get('export_identity')
            or config.get('study') != plan.get('study_id')
            or plan.get('failure_stage') != FAILED_STAGE or plan.get('failure_attempt') != 'attempt-000'):
        raise PipelineIntegrityError('Pinned failed run or exporter changed')
    cap = _positive_int(plan.get('prompt_byte_limit'), 'prompt_byte_limit')
    if cap != _archive_prompt_cap(config):
        raise PipelineIntegrityError('Failed plan prompt byte limit differs from archived configuration')
    if output.exists():
        receipt_path = output / 'export-receipt.json'
        if not receipt_path.is_file():
            raise PipelineIntegrityError('Existing failed public export lacks its receipt')
        receipt = json.loads(receipt_path.read_text())
        if receipt.get('plan_sha256') != sha(plan_path):
            raise PipelineIntegrityError('Public failed export belongs to a different plan')
        observed = {str(path.relative_to(output)): sha(path) for path in output.rglob('*')
                    if path.is_file() and path.name != 'export-receipt.json'}
        if observed != receipt.get('output_sha256'):
            raise PipelineIntegrityError('Public failed export changed')
        return {'status': 'PASS', 'study_id': plan['study_id'], 'trace_sha256': sha(output/'audit_trace.json'), 'resume': True}
    if config.get('study') != packet.get('study_id') or packet.get('study_id') != plan['study_id']:
        raise ValueError('Failed source identity differs from pinned plan')
    source = _failed_public_sources(packet)
    history = _failed_history(run, plan, config, packet, failure, slots, validations, snapshot)
    domains = []
    for name in AUDIT_DOMAINS:
        domains.append({'id': name, 'title': name.replace('_', ' ').capitalize(),
            'model_assessment': {'status': 'CANNOT VERIFY',
                'summary': 'The retained blinded-design response failed its contract; no model assessment is accepted',
                'evidence': [], 'unresolved': ['No validated source-domain assessment is available from the failed response'],
                'proposed_repair': None, 'verification_needed': ['A separately budgeted validated source review']},
            'verification': {'status': 'CANNOT VERIFY', 'scope': 'No validated assessment to check', 'checks': []}})
    trace = {'schema': 'clockbase-audit-trace-v1', 'study_id': plan['study_id'],
        'run_id': 'sha256:' + _digest(snapshot), 'recorded_at_utc': datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
        'scope': 'Source-only development record of a failed blinded-design response; no downstream execution or scientific result',
        'evaluation_split': 'development',
        'models': [{'provider': 'openai', 'model': history['attempts'][0]['requested_model'],
                    'role': 'source-design reviewer'}],
        'overall': {'status': 'CANNOT VERIFY', 'reason': 'The blinded-design response failed validation; no source-domain judgment, comparison, numerical result or biological inference is accepted'},
        'domains': domains, 'comparisons': [], 'covariates': [], 'artifacts': [],
        'events': [{'stage': 'metadata and molecular input', 'status': 'PASS',
                    'summary': 'The two completed upstream stage receipts and listed payload hashes were verified'},
                   {'stage': 'blinded design response contract', 'status': 'FAIL',
                    'summary': 'The retained response-validation record(s) failed; model response content is not treated as an assessment'},
                   {'stage': 'preprocessing and clocks', 'status': 'CANNOT VERIFY',
                    'summary': 'These downstream stages were not reached in the pinned failure profile'},
                   {'stage': 'repair application', 'status': 'CANNOT VERIFY',
                    'summary': 'No repair or scientific result was applied'}],
        'resource_usage': {'generation_attempts': len(slots), 'responses_received': len(slots),
            'estimated_cost_usd': float(sum((Decimal(str(item['estimated_cost_usd'])) for item in history['attempts']), Decimal(0))),
            'billed_cost_usd': None},
        'limitations': ['This is a bounded development failure record, not a completed pipeline analysis.',
            'Published source prose may contain findings; the archived source packet is not certified outcome-blinded.',
            'No validated model assessment, comparison, molecular analysis, clock execution, causal claim or error rate is supplied.',
            'The response, validation, API and upstream artifacts are retained by hash; raw model/API content is not published in the failure history.',
            'Estimated costs are copied from archived receipts and are not invoices.']}
    if not validations:
        trace['events'][1]['summary'] = ('The recorded zero-retry contract adjustment was REJECTED; '
            'trusted export-time replay reproduces the failure without creating historical validation feedback')
        trace['limitations'][3] = ('The original response, rejected adjustment, API and upstream artifacts are retained by hash. '
            'Raw model/API content is omitted; the export-time diagnostic is separate from the recorded failure. '
            'No response-validation artifact or untransmitted retry is synthesized.')
    artifacts = {'source-records': source, 'failure-history': history,
                 'stage-manifest': {'schema': 'clockbase-failed-stage-manifest-v1',
                                    'study_id': plan['study_id'], 'files': snapshot,
                                    'completed_stages': list(STAGES[:2]), 'failed_stage': FAILED_STAGE,
                                    'api_response_count': len(slots), 'prompt_byte_limit': cap}}
    for value in [trace, *artifacts.values()]: _safe(value)
    output.mkdir(parents=True, exist_ok=False); (output / 'public').mkdir()
    for label, value in artifacts.items():
        path = output / 'public' / (label + '.json')
        write_new(path, value)
        trace['artifacts'].append({'id': label, 'kind': 'public failure evidence',
            'label': label.replace('-', ' ').capitalize(), 'sha256': sha(path), 'download_path': 'public/' + path.name})
    trace['artifacts'].append({'id': 'source-metadata', 'kind': 'source hash',
                               'label': 'Original source metadata digest', 'sha256': packet['source_sha256'],
                               'url': _public_url(packet['source_url'])})
    write_new(output / 'audit_trace.json', trace)
    write_new(output / 'export-receipt.json', {'schema': 'clockbase-failed-trace-export-v1',
        'plan_sha256': sha(plan_path), 'source_snapshot_sha256': _digest(snapshot),
        'new_model_calls': 0, 'clocks_rerun': False,
        'output_sha256': {str(path.relative_to(output)): sha(path) for path in output.rglob('*')
                          if path.is_file() and path.name != 'export-receipt.json'}})
    return {'status': 'PASS', 'study_id': plan['study_id'], 'trace_sha256': sha(output/'audit_trace.json'),
            'scientific_status': 'CANNOT VERIFY', 'new_model_calls': 0,
            'api_responses': len(slots), 'next_retry_transmitted': len(slots) > 1}


# Explicit names for callers that keep preparation and export workflows separate.
prepare_failed_run = prepare_failed
export_failed_run = export_failed


def _number(value):
    if value is None or value == '': return None
    if isinstance(value, bool): raise ValueError('A boolean is not a numerical estimate')
    number = float(value)
    if not math.isfinite(number): raise ValueError('Nonfinite estimate cannot enter a public trace')
    return number


def _ids(values, source_ids):
    if (not isinstance(values, list)
            or any(not isinstance(v, str) or not re.fullmatch(r'GSM\d+', v) for v in values)
            or len(set(values)) != len(values)
            or set(values) - source_ids):
        raise ValueError('Invalid source sample identities in numerical record')
    return list(values)


def numerical_record(row):
    values = {key: _number(row.get(key)) for key in ('p_value', 'mean_difference', 'effect_size')}
    if values['p_value'] is not None and not 0 <= values['p_value'] <= 1:
        raise ValueError('P value outside [0,1]')
    present = any(value is not None for value in values.values())
    declared = row.get('status', '')
    state = ('COMPUTED' if declared == 'RESULT' and all(value is not None for value in values.values())
             else 'INCONSISTENT' if present or declared == 'RESULT'
             else 'NOT_REPORTED' if declared == 'NOT_REPORTED' else 'ABSTAINED')
    text = {key: row.get(key) or '' for key in ('statistical_test', 'reason')}
    definition = row.get('effect_size_definition') or None
    if any(not isinstance(value, str) for value in text.values()) or (definition is not None and not isinstance(definition, str)):
        raise ValueError('Numerical labels must be text')
    return {'state': state, **text, **values, 'effect_size_definition': definition}


def verify_numbers(row, frame):
    """Conditional sample-level arithmetic, independently of pipeline helpers."""
    record = numerical_record(row)
    scope = 'Recomputed during export from saved sample predictions; independence and clock validity unverified'
    def result(name, status, reported=None, recomputed=None, reason=''):
        return {'name': name, 'status': status, 'scope': scope, 'reported': reported,
                'recomputed': recomputed,
                'absolute_difference': abs(reported-recomputed) if reported is not None and recomputed is not None else None,
                'reason': reason}
    if record['state'] != 'COMPUTED':
        return [result('declared calculation', 'CANNOT VERIFY',
                       reason='No complete RESULT declaration; no numerical result inferred from an abstention')]
    cases, controls = row['case_ids'], row['control_ids']
    if min(len(cases), len(controls)) < 2 or set(cases) & set(controls):
        return [result('sample membership', 'FAIL', reason='Unpaired calculation needs disjoint groups with at least two sample records')]
    if row['clock'] not in frame or set(cases+controls) - set(frame.index):
        return [result('sample membership', 'FAIL', reason='Requested samples or clock are absent from the saved predictions')]
    a, b = frame.loc[cases, row['clock']].to_numpy(float), frame.loc[controls, row['clock']].to_numpy(float)
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return [result('finite predictions', 'FAIL', reason='Saved group predictions are not finite')]
    difference = float(a.mean()-b.mean())
    output = [result('mean difference', 'PASS' if math.isclose(record['mean_difference'], difference, rel_tol=1e-8, abs_tol=1e-10) else 'FAIL',
                     record['mean_difference'], difference)]
    test = record['statistical_test'].casefold()
    # Do not reinterpret adjusted, paired, repeated or multigroup declarations.
    unadjusted = not re.search(r'\b(paired|adjusted|ancova|mixed|repeated|regression|one[- ]sided|one[- ]tailed)\b', test)
    if unadjusted and 'welch' in test:
        p = float(stats.ttest_ind(a, b, equal_var=False, alternative='two-sided').pvalue)
    elif unadjusted and ('student' in test or test in {'t-test', 'two-sample t-test'}):
        p = float(stats.ttest_ind(a, b, equal_var=True, alternative='two-sided').pvalue)
    else: p = None
    if p is None or not math.isfinite(p):
        output.append(result('P value', 'CANNOT VERIFY', record['p_value'],
                             reason='Declared method is unsupported or yields an undefined conditional test'))
    else:
        output.append(result('P value', 'PASS' if math.isclose(record['p_value'], p, rel_tol=1e-7, abs_tol=1e-12) else 'FAIL',
                             record['p_value'], p))
    if record['effect_size_definition'] == EFFECT:
        variance = ((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2)
        effect = difference/math.sqrt(variance) if variance > 0 else None
        output.append(result('standardized effect', 'CANNOT VERIFY' if effect is None else
            'PASS' if math.isclose(record['effect_size'], effect, rel_tol=1e-8, abs_tol=1e-10) else 'FAIL',
            record['effect_size'], effect, 'Explicit pooled-sample-SD definition only'))
    else:
        output.append(result('effect definition', 'CANNOT VERIFY', record['effect_size'],
                             reason='Effect definition was not supplied or is unsupported; no definition inferred from its value'))
    return output


def _tabular_format(path):
    """Retain the declared delimiter and detected compression in public names."""
    path = Path(path)
    with path.open('rb') as raw:
        compressed = raw.read(2) == b'\x1f\x8b'
    compressed = compressed or path.suffix == '.gz'
    tabular_suffix = path.with_suffix('').suffix if path.suffix == '.gz' else path.suffix
    extension = '.tsv' if tabular_suffix == '.tsv' else '.csv'
    return compressed, extension


def _matrix_qc(path, source_ids, modality):
    """Inspect values without retaining the full matrix in memory."""
    compressed, extension = _tabular_format(path)
    opener = gzip.open if compressed else open
    with opener(path, 'rt', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t' if extension == '.tsv' else ',')
        header = next(reader, None)
        if not header or len(header)<2:
            raise ValueError('Molecular matrix is empty or lacks sample columns')
        columns = header[1:]
        _ids(columns, source_ids)
        if set(columns) != source_ids: raise ValueError('Molecular matrix and GEO sample roster differ')
        genes, count, missing, out_of_range = set(), 0, 0, 0
        for row in reader:
            if len(row) != len(header) or not row[0] or row[0] in genes:
                raise ValueError('Molecular matrix feature identifiers or row lengths are invalid')
            genes.add(row[0]); count += 1
            for value in row[1:]:
                if value.strip().lower() in {'', 'na', 'nan'}:
                    missing += 1
                    continue
                number = _number(value)
                if modality == 'human_dnam':
                    out_of_range += int(not 0 <= number <= 1)
                elif modality in {'human_rna', 'mouse_rna'}:
                    out_of_range += int(number < 0 or not number.is_integer())
                else:
                    raise ValueError('Unrecognized molecular measurement profile')
    if not count: raise ValueError('Molecular matrix contains no feature measurements')
    return count, len(columns), {'name': 'molecular values',
        'status': 'FAIL' if out_of_range else 'CANNOT VERIFY' if missing else 'PASS',
        'scope': 'Rechecked raw input during export; missing values were not replaced',
        'reason': f'{count*len(columns)} cells inspected; {missing} missing; {out_of_range} outside declared '
                  + ('beta range [0,1]' if modality == 'human_dnam' else 'nonnegative integer count domain')}


def _public_url(url):
    parts = urlparse(url)
    if (parts.scheme != 'https' or parts.hostname not in {'www.ncbi.nlm.nih.gov', 'pmc.ncbi.nlm.nih.gov',
            'pubmed.ncbi.nlm.nih.gov', 'www.ebi.ac.uk'} or parts.username or parts.password):
        raise ValueError('Unrecognized public evidence URL')
    if parts.hostname == 'www.ebi.ac.uk' and not re.fullmatch(r'/europepmc/webservices/rest/PMC\d+/fullTextXML', parts.path):
        raise ValueError('Unrecognized Europe PMC source path')
    return url


def _safe(value):
    if isinstance(value, dict):
        for key, item in value.items(): _safe(key); _safe(item)
    elif isinstance(value, list):
        for item in value: _safe(item)
    elif isinstance(value, str) and re.search(
            r'/Entropy/|/Lasagna/|/home/|/Users/|Bearer\s+|sk-[A-Za-z0-9_-]{20,}'
            r'|(?:^|[\s\x22\x27(=])/[A-Za-z0-9._-]+/[A-Za-z0-9._/-]+'
            r'|[A-Za-z0-9.!#$%&\x27*+/=?^_`{|}~-]+@[A-Za-z0-9-]+(?:\.[A-Za-z0-9-]+)+', value):
        raise ValueError('Potentially private text refused; source quotes were not rewritten')


def grouping_readback(packet, comparisons):
    """Public diagnostics preserve request checks beside unapplied witnesses."""
    records, checks = [], []
    for comparison in comparisons:
        scope = {**packet, 'target': comparison['target']}
        record = source_grouping_diagnostics(scope, comparison['case_ids'], comparison['control_ids'])
        record.update(comparison_id=comparison['comparison_id'], comparison_name=comparison['comparison_name'])
        records.append(record)
        binding = record['target_binding']
        reasons = '; '.join(binding['failures'] + binding['unresolved'])
        checks.append({'status': binding['status'],
            'description': comparison['comparison_name'] + ': literal target binding to complete directional source groups. '
                           + (reasons or 'Exact membership matches; target meaning and control roles remain unverified.')})
        for witness in record['source_membership_witnesses']:
            checks.append({'status': witness['source_binding']['status'],
                'description': comparison['comparison_name'] + ': unapplied source-label witness: '
                    + witness['proposed_target'] + '. Conditional on the supplied groups; no scientific approval or repair.'})
    return {'schema': 'clockbase-public-grouping-readback-v1', 'study_id': packet['study_id'],
            'scientific_status': 'CANNOT VERIFY', 'reference_judgment': False,
            'repairs_applied': [], 'comparisons': records}, checks


def export(plan_path, output):
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = json.loads(plan_path.read_text())
    if plan.get('schema') == 'clockbase-failed-pipeline-trace-plan-v1':
        return export_failed(plan_path, output)
    if plan.get('schema') != 'clockbase-pipeline-trace-plan-v1' or plan.get('split') != 'development':
        raise ValueError('Unsupported pipeline trace plan')
    run = Path(plan['run']).resolve()
    snapshot = _snapshot(run)
    if snapshot != plan['snapshot_sha256'] or _export_identity() != plan['export_identity']:
        raise PipelineIntegrityError('Pinned run or exporter changed')
    if output.exists():
        receipt = json.loads((output/'export-receipt.json').read_text())
        if receipt['plan_sha256'] != sha(plan_path): raise PipelineIntegrityError('Public export belongs to a different plan')
        observed = {str(p.relative_to(output)): sha(p) for p in output.rglob('*')
                    if p.is_file() and p.name != 'export-receipt.json'}
        if observed != receipt['output_sha256']: raise PipelineIntegrityError('Public export changed')
        return {'status': 'PASS', 'trace_sha256': sha(output/'audit_trace.json'), 'resume': True}
    stage = {name: json.loads((run/name/'complete.json').read_text()) for name in STAGES}
    data = {name: value['result'] for name, value in stage.items()}
    config = json.loads((run/'config.json').read_text())
    study = config['study']
    if study != plan['study_id'] or not re.fullmatch(r'GSE\d+', study): raise ValueError('Study identity changed')
    metadata, design, raw, clocks, primary, audit = (data[name] for name in
        ('01_metadata', '03_blinded_design', '02_molecular_input', '05_clock_predictions', '06_coding_agent', '07_audit_and_repair'))
    def registered(path):
        path = _inside(run, path)
        if str(path.relative_to(run)) not in snapshot: raise ValueError('Input is absent from the pinned stage manifest')
        return path
    packet = source_packet(registered(metadata['metadata']), study, config['target'],
                           context=registered(metadata['study_context']) if metadata.get('study_context') else None)
    ids = [row['sample_id'] for row in packet['samples']]
    source_ids = set(ids)
    _ids(ids, source_ids)
    if set(metadata['sample_ids']) != source_ids: raise ValueError('GEO metadata receipt and sample packet differ')
    raw_path, predictions_path = registered(raw['path']), registered(clocks['path'])
    features, molecular_samples, molecular_check = _matrix_qc(raw_path, source_ids, config['modality'])
    with predictions_path.open(newline='') as handle:
        prediction_header = next(csv.reader(handle), None)
    if (not prediction_header or prediction_header[0]!='sample_id'
            or len(prediction_header)!=len(set(prediction_header))
            or not isinstance(clocks['clocks'],list)
            or any(not isinstance(name,str) or not name for name in clocks['clocks'])
            or len(clocks['clocks'])!=len(set(clocks['clocks']))):
        raise ValueError('Prediction columns or declared clocks are empty or duplicated')
    frame = pd.read_csv(predictions_path).set_index('sample_id', verify_integrity=True)
    if set(frame.index) != source_ids or set(frame.columns) != set(clocks['clocks']) or not np.isfinite(frame.to_numpy(float)).all():
        raise ValueError('Prediction identities or finite values differ from the declared source and clocks')
    if clocks['sample_count'] != len(ids): raise ValueError('Clock sample count differs from source roster')
    documents = {doc['source_id']: doc for doc in packet.get('study_documents', [])}
    def evidence(items):
        from .agent_workflow_2026_09_10 import _check_evidence
        checked, _ = _check_evidence(items, packet)
        return [{**item, 'source_url': _public_url('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+item['sample_id']
                    if 'sample_id' in item else documents[item['source_id']]['source_url'])} for item in checked]
    source = {'schema': 'clockbase-public-source-records-v1', 'study_id': study,
        'source_metadata_sha256': packet['source_sha256'], 'source_url': _public_url(packet['source_url']),
        'samples': packet['samples'], 'study_documents': [{key: doc[key] for key in
            ('source_id', 'source_url', 'source_sha256', 'fields')} for doc in documents.values()]}
    for doc in source['study_documents']: _public_url(doc['source_url'])
    domains = []
    for name in AUDIT_DOMAINS:
        assessment = dict(design.get('audit_domains', {}).get(name) or {
            'status': 'CANNOT VERIFY', 'summary': 'This archived review did not record a separate assessment for this domain',
            'evidence': [], 'unresolved': ['A domain-specific source assessment is unavailable in this run'],
            'proposed_repair': None, 'verification_needed': ['Independent source and method assessment']})
        assessment['evidence'] = evidence(assessment['evidence'])
        domains.append({'id': name, 'title': name.replace('_', ' ').capitalize(), 'model_assessment': assessment,
            'verification': {'status': 'CANNOT VERIFY', 'scope': 'Export verifies recorded artifacts and conditional calculations only',
                'checks': [{'status': 'PASS', 'description': 'All registered stage artifact hashes verified'}] if name == 'data_provenance' else []}})
    comparisons, seen = [], set()
    children = design.get('comparisons', [design])
    for child in children:
        case, control = _ids(child['case_ids'], source_ids), _ids(child['control_ids'], source_ids)
        cid = child.get('comparison_id') or 'cmp-'+_digest([study,case,control,child.get('target',config['target'])])[:20]
        if cid in seen: raise ValueError('Duplicate archived comparison identity')
        seen.add(cid)
        comparisons.append({'comparison_id': cid, 'comparison_name': child.get('comparison_name', config['target']),
            'target': child.get('target', config['target']), 'status': 'CANNOT VERIFY', 'case_ids': case,
            'control_ids': control, 'excluded_ids': _ids(child['excluded_ids'],source_ids),
            'biological_unit': child['biological_unit'], 'pairing': child['pairing'],
            'unresolved': list(design.get('unresolved',[])) if child is design else
                          list(design.get('unresolved',[]))+list(child.get('unresolved',[]))})
    grouping, grouping_checks = grouping_readback(packet, comparisons)
    grouping_domain = next(domain for domain in domains if domain['id'] == 'grouping')
    grouping_domain['verification'] = {
        'status': 'FAIL' if any(check['status'] == 'FAIL' for check in grouping_checks) else 'CANNOT VERIFY',
        'scope': 'Exact source labels and conditional membership only; no target interpretation or scientific reference',
        'checks': grouping_checks}
    results, result_ids = [], set()
    for original in primary['rows']:
        case, control = _ids(original['case_ids'],source_ids), _ids(original['control_ids'],source_ids)
        if set(case) & set(control): raise ValueError('Overlapping original case/control identities')
        if original['clock'] not in frame: raise ValueError('Unrecognized original clock')
        cid = original.get('comparison_id') or 'cmp-'+_digest([study,case,control,original.get('target',config['target'])])[:20]
        rid = cid+'/'+original['clock']
        if rid in result_ids: raise ValueError('Duplicate original numerical result')
        result_ids.add(rid)
        def matching(rows):
            matches = []
            for row in rows:
                if row['clock'] != original['clock']:
                    continue
                if 'comparison_id' not in row and len(children) > 1:
                    raise ValueError('Legacy repair lacks an unambiguous comparison identity')
                if row.get('comparison_id', cid) == cid:
                    matches.append(row)
            return matches
        proposed, selected = matching(audit.get('candidates',[])), matching(audit.get('final_outputs',[]))
        if len(proposed)>1 or len(selected)>1: raise ValueError('Ambiguous candidate or selected result identity')
        checks = verify_numbers(original,frame)
        proposed_groups = None
        if proposed:
            proposed_groups = {key: _ids(proposed[0][key],source_ids) for key in ('case_ids','control_ids')}
            if set(proposed_groups['case_ids']) & set(proposed_groups['control_ids']):
                raise ValueError('Overlapping proposed case/control identities')
            checks.extend({**r,'name':'proposal '+r['name']} for r in verify_numbers(proposed[0],frame))
        final = selected[0] if selected else {}
        selection = {'output': final.get('selected_output','unrecorded'),
            'status': final.get('scientific_status','CANNOT VERIFY'),
            'accepted_result': final.get('accepted_result',False),
            'repair_applied': final.get('repair_applied'),
            'reason': final.get('selection_reason','No separate selection rationale was saved by this archived version')}
        if (selection['output'] not in {'original','candidate','none','unrecorded'}
                or selection['status'] not in STATUS or type(selection['accepted_result']) is not bool
                or (selection['repair_applied'] is not None and type(selection['repair_applied']) is not bool)):
            raise ValueError('Invalid archived selection declaration')
        selected_row = original if selection['output']=='original' else proposed[0] if selection['output']=='candidate' and proposed else None
        if (selection['output']=='candidate' and not proposed or selection['repair_applied'] and selection['output']!='candidate'
                or selection['accepted_result'] and (selection['status']!='PASS' or selected_row is None or numerical_record(selected_row)['state']!='COMPUTED')):
            raise ValueError('Inconsistent archived selection declaration')
        if final and selected_row is not None:
            if (any(final[key] != selected_row[key] for key in ('case_ids','control_ids'))
                    or numerical_record(final) != numerical_record(selected_row)):
                raise ValueError('Selected output differs from its declared original or proposal')
        results.append({'result_id':rid,'comparison_id':cid,'clock':original['clock'],
            'case_ids':case,'control_ids':control,'original':numerical_record(original),
            'proposal':numerical_record(proposed[0]) if proposed else None,
            'proposal_groups': proposed_groups, 'selection': selection,
            'checks':checks})
    covariates = [{**row,'evidence':evidence(row['evidence'])} for row in design.get('covariate_plan',[])]
    receipts = [json.loads(path.read_text()) for path in sorted((run/'api').glob('*/receipt.json'))]
    for path in (run/'api').glob('*/receipt.json'):
        record = json.loads(path.read_text())
        if sha(path.parent/'response.json') != record['response_sha256']: raise PipelineIntegrityError('API response hash mismatch')
    processing = data['04_preprocessing']
    keys = ('log_transform','scaling','differential_center','RLE_applied','YuGene_applied','filter_order','normalization',
            'raw_methylation_preprocessing_recomputed','independent_imputation_verification','idat_processing')
    multiplicity = audit.get('multiplicity',{})
    trace = {'schema':'clockbase-audit-trace-v2','study_id':study,'run_id':'sha256:'+_digest(snapshot),
        'recorded_at_utc':datetime.now(timezone.utc).isoformat().replace('+00:00','Z'),
        'scope':'Saved full-pipeline evidence with new conditional numerical readback; no clock rerun during export',
        'evaluation_split':'development',
        'models':[{'provider':'openai','model':config['model'],'role':'archived original CodingAgent'},
                  {'provider':'openai','model':config.get('audit_model') or config['model'],'role':'archived source-design reviewer'}],
        'overall':{'status':'CANNOT VERIFY','reason':'Independent scientific reference, biological units and clock applicability are not established by artifact or arithmetic checks'},
        'domains':domains,'comparisons':comparisons,'covariates':covariates,'artifacts':[],
        'events':[{'stage':'archived pipeline','status':'PASS','summary':'Eight recorded stages and all pinned payload hashes verified'},
                  {'stage':'export numerical readback','status':'FAIL' if any(c['status']=='FAIL' for r in results for c in r['checks']) else 'PASS',
                   'summary':'Saved sample predictions checked against explicit original declarations; unperformed and undefined checks retained'}],
        'execution':{'stage_count':len(STAGES),'completed_stage_count':len(stage),'clocks_rerun_during_export':False,
            'stages':[{'stage':name,'status':'PASS','elapsed_seconds':stage[name]['elapsed_seconds'],
                       'artifact_count':len(stage[name]['outputs'])} for name in STAGES]},
        'data':{'sample_ids':ids,'sample_count':len(ids),'molecular_feature_count':features,
            'molecular_sample_count':molecular_samples,'clock_names':list(frame.columns),
            'preprocessing':[{'name':key,'value':str(processing[key])} for key in keys if key in processing],
            'checks':[{'name':'molecular matrix','status':'PASS','scope':'Rechecked during export',
                       'reason':'Nonempty unique features and exact source sample roster'}, molecular_check,
                      {'name':'prediction identity','status':'PASS','scope':'Rechecked during export',
                       'reason':'Unique source samples, declared clocks and finite saved predictions'}]},
        'analysis':{'scope':'Original and proposed statistics; saved-output readback does not rerun preprocessing or clocks',
            'multiplicity':{'method':multiplicity.get('method'),
                'attempted_tests':multiplicity.get('attempted_endpoints'),
                'denominator_status':multiplicity.get('denominator_status','CANNOT VERIFY')},'results':results},
        'resource_usage':{'generation_attempts':len(list((run/'api').glob('*/reservation.json'))),
            'responses_received':len(receipts),'estimated_cost_usd':float(sum((Decimal(str(r['estimated_cost_usd'])) for r in receipts),Decimal(0))),
            'billed_cost_usd':None},
        'limitations':['This export reads completed development runs; it is not a new full pipeline execution.',
            'Scientific acceptance and organismal aging effects are not certified by numerical agreement.',
            'Unspecified effect definitions are retained; standardized effects are not age in years.',
            'Source data and derived predictions are public-study records, not proof of independent biological units.',
            'Estimated costs describe the archived model calls; export made zero model calls and invoices are unavailable.',
            'The original scripts and dependency versions are identified by archived configuration hashes; they were not executed by export.']}
    artifacts = {'source-records':source,'numerical-readback':trace['analysis'], 'grouping-readback':grouping,
        'stage-manifest':{'schema':'clockbase-public-stage-manifest-v1','files':snapshot,
                          'source_implementation_sha256':config.get('implementation',{}),
                          'export_identity':plan['export_identity'],
                          'source_config_sha256':sha(run/'config.json')}}
    for value in [trace,*artifacts.values()]: _safe(value)
    output.mkdir(parents=True,exist_ok=False);(output/'public').mkdir()
    for label,value in artifacts.items():
        path=output/'public'/(label+'.json');write_new(path,value)
        trace['artifacts'].append({'id':label,'kind':'verified export','label':label.replace('-',' ').capitalize(),
            'sha256':sha(path),'download_path':'public/'+path.name})
    for label,path in (('sample-predictions',predictions_path),('molecular-matrix',raw_path)):
        record={'id':label,'kind':'saved data','label':label.replace('-',' ').capitalize(),'sha256':sha(path)}
        if path.stat().st_size <= plan['max_download_bytes']:
            compressed, extension = _tabular_format(path)
            target=output/'public'/(label+extension+('.gz' if compressed else ''));shutil.copyfile(path,target)
            if sha(target)!=record['sha256']:raise PipelineIntegrityError('Export copy changed')
            record['download_path']='public/'+target.name
        trace['artifacts'].append(record)
    if _snapshot(run)!=snapshot: raise PipelineIntegrityError('Source run changed during export')
    write_new(output/'audit_trace.json',trace)
    write_new(output/'export-receipt.json',{'schema':'clockbase-pipeline-trace-export-v1',
        'plan_sha256':sha(plan_path),'source_snapshot_sha256':_digest(snapshot),
        'new_model_calls':0,'clocks_rerun':False,
        'output_sha256':{str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file()}})
    return {'status':'PASS','study_id':study,'trace_sha256':sha(output/'audit_trace.json'),
            'original_results':len(results),'scientific_status':'CANNOT VERIFY','new_model_calls':0}


def _corpus_projection_helpers(observed):
    def checked(path, expected=None):
        path = Path(path).resolve()
        actual = sha(path)
        if expected is not None and actual != expected:
            raise PipelineIntegrityError('Corpus export source hash changed')
        observed[path] = actual
        return path
    def read(path):
        return json.loads(checked(path).read_text())
    def clean(value):
        if isinstance(value, dict):
            return {str(key): clean(item) for key, item in value.items() if key not in {
                'path', 'source_path', 'command', 'stdout', 'stderr', 'safe_text', 'previous_response', 'raw_response'}}
        if isinstance(value, list):
            return [clean(item) for item in value]
        if isinstance(value, str):
            try:
                _safe(value)
            except ValueError:
                return '[Text withheld from public export; SHA256 ' + hashlib.sha256(value.encode()).hexdigest() + ']'
        return value
    def project(value, keys):
        return clean({key: value[key] for key in keys if key in value})
    def stage(study, name):
        directory = study / name
        path = directory / 'complete.json'
        if not path.exists():
            return None
        expected = checked(directory/'complete.sha256').read_text().strip()
        value = read(checked(path, expected))
        for relative, digest in value['outputs'].items():
            target = _inside(directory, directory/relative)
            checked(target, digest)
        return value
    def artifact_from_stage(study, name, basename):
        complete = stage(study, name)
        if complete is None:
            return None
        hits = [relative for relative in complete['outputs'] if Path(relative).name == basename]
        if len(hits) != 1:
            raise ValueError('Corpus stage does not identify one ' + basename)
        return read(study/name/hits[0])
    return checked, read, clean, project, stage, artifact_from_stage


def _corpus_record_hash(value):
    """Digest convention used by the saved corpus contracts (includes spaces)."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _measurement_projection(study, terminal, task, source_index, frozen_design, aggregate, helpers):
    """Export method records only when the bound reviewer requests received them."""
    checked, read, clean, project, stage, artifact_from_stage = helpers
    completed = stage(study, '02c_measurement_evidence')
    if completed is None:
        if terminal.get('measurement_evidence_sha256'):
            raise ValueError('Terminal measurement registry lacks its completed stage')
        return None, {}
    value = completed['result']; registry_hash = _corpus_record_hash(value)
    if (value.get('schema') != 'clockbase-pipeline-measurement-evidence-v1'
            or value.get('study_id') != task['study_id'] or value.get('modality') != task['modality']
            or value.get('scientific_status') != 'CANNOT VERIFY'
            or terminal.get('measurement_evidence_sha256') != registry_hash
            or terminal.get('measurement_evidence_records') != len(value.get('evidence', {}))):
        raise ValueError('Measurement registry differs from its source terminal')
    records = value['evidence']
    if not records or any(not re.fullmatch(r'M[0-9]{5}', key) or key in source_index.get('evidence', {}) for key in records):
        raise ValueError('Measurement registry identifiers are missing or collide with source evidence')
    context = task.get('pipeline_measurement_context') or {}
    pipeline = Path(context.get('root', '')).resolve()
    source_hashes = {}
    for role, declaration in value['source_files'].items():
        checked(declaration['path'], declaration['sha256'])
        source_hashes[role] = declaration['sha256']
    stage_hashes = {}
    required_stages = {'02_molecular_input', '04_preprocessing', '05_clock_predictions'}
    if any(row.get('kind') == 'source_linked_dnam_input_recovery' for row in records.values()):
        required_stages.add('01_metadata')
    if set(value['stage_receipts']) != required_stages:
        raise ValueError('Measurement registry must bind input, preprocessing and clock stages')
    for name, declaration in value['stage_receipts'].items():
        path = checked(declaration['path'], declaration['sha256'])
        if path != pipeline/name/'complete.json' or stage(pipeline, name) is None:
            raise ValueError('Measurement receipt belongs to a different pipeline stage')
        stage_hashes[name] = declaration['sha256']
    source_records = clean(records); request_bindings = []; limit = None
    for scope in aggregate['scoped_scientific_checks']:
        name = scope['review_stage']
        if not re.fullmatch(r'04_review_[0-9]{4}', name):
            raise ValueError('Unsupported measurement review stage')
        reviewed = stage(study, name)
        if reviewed is None or reviewed['result'].get('checks') != scope['checks']:
            raise ValueError('Measurement review scope differs from its completed review')
        accepted = []
        for relative in reviewed['outputs']:
            if re.fullmatch(r'contract-[0-9]{2}\.json', Path(relative).name):
                contract = read(study/name/relative)
                if contract.get('status') in {'PASS', 'PARTIAL_WITH_CONSERVATIVE_REJECTIONS'}:
                    accepted.append(contract)
        if len(accepted) != 1 or type(accepted[0].get('attempt')) is not int or not 0 <= accepted[0]['attempt'] <= 1:
            raise ValueError('Measurement review does not identify one accepted request')
        slot = name[3:] + f"_{accepted[0]['attempt']:02d}-000"
        request_rel = 'api/'+slot+'/request.json'; binding_rel = 'api/'+slot+'.aggregate-binding.json'
        outputs = terminal['output_sha256']
        if request_rel not in outputs or binding_rel not in outputs:
            raise ValueError('Reviewer request or aggregate binding is absent from the terminal ledger')
        request_path = checked(study/request_rel, outputs[request_rel])
        binding_path = checked(study/binding_rel, outputs[binding_rel])
        request, binding = read(request_path), read(binding_path)
        if (binding.get('schema') != 'clockbase-aggregate-call-binding-v1'
                or binding.get('request_file_sha256') != sha(request_path) or binding.get('call_slot') != slot):
            raise ValueError('Reviewer API request binding has changed')
        messages = request.get('messages', [])
        if len(messages) != 2 or messages[0].get('role') != 'system' or messages[1].get('role') != 'user':
            raise ValueError('Unsupported reviewer API message envelope')
        payload = json.loads(messages[1]['content']); evidence = payload['source']['evidence']
        if binding.get('aggregate_key') != _corpus_record_hash({'study': str(study.resolve()), 'actor': slot[:-4],
                'system': messages[0]['content'], 'payload': payload}):
            raise ValueError('Reviewer request no longer matches its aggregate identity')
        supplied = {key: item for key, item in evidence.items() if re.fullmatch(r'M[0-9]{5}', key)}
        if (payload.get('study_id') != task['study_id'] or payload['source'].get('study_id') != task['study_id']
                or payload.get('measurement_evidence_sha256') != registry_hash or supplied != records
                or payload.get('frozen_source_design') != frozen_design
                or set(scope['row_ids']) != {key for key in evidence if re.fullmatch(r'R[0-9a-f]{24}', key)}
                or evidence.get('L00000') != source_index.get('evidence', {}).get('L00000')):
            raise ValueError('Reviewer payload did not receive the exact measurement registry, rows and frozen design')
        current_limit = evidence.get('L00001')
        if (not isinstance(current_limit, dict) or current_limit.get('field') != 'post_design_measurement_scope'
                or not isinstance(current_limit.get('value'), str) or not current_limit['value'].strip()
                or (limit is not None and current_limit != limit)):
            raise ValueError('Reviewer post-design evidence scope is missing or inconsistent')
        limit = current_limit
        request_bindings.append({'review_stage': name, 'request_sha256': sha(request_path),
            'aggregate_binding_sha256': sha(binding_path), 'measurement_evidence_sha256': registry_hash,
            'measurement_evidence_ids': sorted(records), 'post_design_scope_id': 'L00001'})
    if limit is not None:
        source_records['L00001'] = clean(limit)
    if 'L00000' in source_index.get('evidence', {}):
        source_records['L00000'] = clean(source_index['evidence']['L00000'])
    public = {'schema': 'clockbase-public-measurement-evidence-v1', 'registry_sha256': registry_hash,
        'stage_sha256': sha(study/'02c_measurement_evidence/complete.json'), 'evidence': clean(records),
        'source_sha256': source_hashes, 'stage_sha256_by_name': stage_hashes,
        'review_request_bindings': request_bindings, 'scientific_status': 'CANNOT VERIFY',
        'scope': 'Measurement records supplied after the unchanged source-design assessment; completeness and applicability are separate judgments.',
        'limitations': clean(value.get('limitations', []))}
    _safe(public)
    return public, source_records


def _intermediate_corpus_counts(value):
    """Count interrupted-stage evidence without adding it to terminal coverage."""
    numerical = value['independent_numerics']
    return {'calculated_targets': sum(c['status'] == 'PASS' for c in value['calculations']),
            'unresolved_targets': sum(c['status'] != 'PASS' for c in value['calculations']),
            'independent_numerical_counts': numerical['counts'] if numerical else {s: 0 for s in STATUS},
            'review_row_assignments': 0, 'reviewed_unique_rows': 0, 'repair_stages': 0,
            'completed_stages': len(value['completed_stages']),
            'failed_attempts': len(value['failed_attempts'])}


def _project_corpus_study(study, result, task, root, helpers, *, report_sha256=None):
    from .independent_statistics_2026_09_20 import model_verification_view
    from .legacy_analysis_review_2026_09_19 import SCIENTIFIC_CHECKS
    checked, read, clean, project, stage, artifact_from_stage = helpers
    row_fields = ('gse_id', 'comparison_name', 'condition_name', 'condition_category', 'clock',
                  'statistical_test', 'covariate', 'p_value', 'pvalue', 'effect_size', 'direction', 'trust_score')
    calculation_fields = ('statistical_test', 'variance_estimator', 'case_ids', 'control_ids', 'n_case', 'n_control',
        'n_case_assays', 'n_control_assays', 'n_biological_units', 'case_mean', 'control_mean', 'mean_difference',
        'adjusted_mean_difference', 'p_value', 'effect_size', 'effect_size_definition', 'confidence_interval',
        'confidence_interval_estimand', 'confidence_level', 'standard_error', 'degrees_of_freedom',
        'pairs', 'covariates', 'scope', 'pooled_standardized_difference', 'pooled_standardized_difference_definition')
    key = result['study_key']
    if key != task['modality'] + ':' + task['study_id'] or not re.fullmatch(r'GSE[0-9]+', task['study_id']):
        raise ValueError('Unsupported corpus study identity')
    terminal_path = study/'terminal.json'
    if not terminal_path.exists():
        if result['status'] != 'CONTROLLER_FAILURE':
            raise ValueError('Completed report lacks its terminal study receipt')
        terminal = result
        terminal_hash = None
    else:
        terminal = read(checked(terminal_path, checked(study/'terminal.sha256').read_text().strip()))
        if terminal != result:
            raise ValueError('Report result differs from the terminal study receipt')
        terminal_hash = sha(terminal_path)
        for relative, expected in terminal['output_sha256'].items():
            checked(_inside(study, study/relative), expected)
    source_config = read(study/'config.json') if (study/'config.json').exists() else {'artifacts': [], 'task': task}
    sources = []
    for role, source in source_config['task']['source_files'].items():
        if source.get('available'):
            checked(source['path'], source['sha256'])
        sources.append({'id': task['study_id']+'_'+role, 'role': role, 'available': source.get('available', False),
                        'sha256': source.get('sha256')})
    index = stage(study, '01_source_index')
    design = stage(study, '02_blinded_design')
    calculation_stage = stage(study, '03_conditional_reruns')
    numeric_stage = stage(study, '03b_independent_r')
    aggregation = stage(study, '05_aggregation')
    repair_stage = stage(study, '06_versioned_repairs')
    # The stage name is part of the exporter contract; no guessed candidate
    # directory or latest attempt is substituted for a completion receipt.
    if repair_stage is None and terminal.get('repair_summary', {}).get('status') == 'PASS':
        raise ValueError('Terminal repair summary has no completed 06_versioned_repairs stage')
    source_index = index['result'] if index else {}
    frozen_design = design['result'] if design else None
    calculations = calculation_stage['result'] if calculation_stage else []
    aggregate = aggregation['result'] if aggregation else {'rows': [], 'scoped_scientific_checks': [], 'design_disagreements': []}
    if len(aggregate['rows']) != terminal.get('reviewed_rows', 0):
        raise ValueError('Terminal reviewed-row count differs from preserved review decisions')
    for scoped in aggregate['scoped_scientific_checks']:
        if set(scoped['checks']) != set(SCIENTIFIC_CHECKS):
            raise ValueError('Public corpus review must retain all sixteen scoped checks')
    original_rows, original_parse_issues = [], []
    for artifact in source_config['artifacts']:
        path = (root/artifact['path']).resolve()
        if not path.is_relative_to(root):
            raise ValueError('Original analysis artifact escapes its source root')
        checked(path, artifact['sha256'])
        try:
            with path.open(newline='') as handle:
                reader = csv.reader(handle, strict=True); header = next(reader, [])
                for number, values in enumerate(reader, 2):
                    rid = 'R' + hashlib.sha256(json.dumps([artifact['path'], artifact['sha256'], number],
                        sort_keys=True, allow_nan=False).encode()).hexdigest()[:24]
                    fields = dict(zip(header, values)) if len(header) == len(set(header)) == len(values) else None
                    original_rows.append({'row_id': rid, 'artifact_id': artifact['artifact_id'],
                        'artifact_sha256': artifact['sha256'], 'logical_record': number,
                        'fields': project(fields, row_fields) if fields else None,
                        'record_status': 'PARSED' if fields else 'MALFORMED_OR_AMBIGUOUS'})
            if not header:
                original_parse_issues.append({'artifact_id': artifact['artifact_id'], 'status': 'FAIL',
                    'reason': 'Original CSV has no header; no numerical or scientific claim is inferred.'})
        except (OSError, csv.Error, UnicodeError) as error:
            original_parse_issues.append({'artifact_id': artifact['artifact_id'], 'status': 'CANNOT VERIFY',
                'reason': clean(type(error).__name__ + ': ' + str(error)),
                'scope': 'Previously parsed logical rows are retained; unparsed records are not invented.'})
    if set(r['row_id'] for r in aggregate['rows']) - set(r['row_id'] for r in original_rows):
        raise ValueError('Public review decision is not bound to an original source row')
    referenced = set()
    def refs(value):
        if isinstance(value, dict):
            for name, item in value.items():
                if name in {'evidence_ids', 'case_all', 'case_any', 'control_all', 'control_any'} and isinstance(item, list):
                    referenced.update(x for x in item if isinstance(x, str) and x.startswith(('S', 'D', 'L')))
                else:
                    refs(item)
        elif isinstance(value, list):
            for item in value: refs(item)
    refs(frozen_design); refs(aggregate)
    source_evidence = {name: project(source_index.get('evidence', {})[name],
        ('artifact', 'field', 'value', 'source_sha256', 'value_sha256', 'sample_ids'))
        for name in sorted(referenced) if name in source_index.get('evidence', {})}
    measurement, method_records = _measurement_projection(study, terminal, source_config['task'], source_index,
        frozen_design, aggregate, helpers)
    source_evidence.update(method_records)
    trace_stage = stage(study, '02b_original_trace')
    trace = project(trace_stage['result'], ('availability', 'source_sha256', 'entry_canonical_sha256',
        'request_binding', 'eligible_for_post_design_context', 'identity_scope', 'identity_limit',
        'evidence_count', 'evidence_characters', 'scope', 'withheld_reason')) if trace_stage else None
    public_calculations = [project(c, ('comparison_id', 'clock', 'case_ids', 'control_ids', 'target', 'status',
        'scientific_status', 'source_design_sha256', 'reason', 'meaning')) |
        {'calculation': project(c['calculation'], calculation_fields) if c.get('calculation') else None}
        for c in calculations]
    numeric = model_verification_view(numeric_stage['result']) if numeric_stage else None
    if numeric_stage:
        raw_numeric = numeric_stage['result'].get('results', [])
        if len(raw_numeric) != len(calculations) or len(numeric['results']) != len(calculations):
            raise ValueError('Public independent calculation coverage differs from its source')
        for index, (calculation, original, public) in enumerate(zip(
                calculations, raw_numeric, numeric['results'], strict=True)):
            if (type(original.get('input_index')) is not int or original['input_index'] != index
                    or any(field not in original or field not in calculation
                           or original[field] != calculation[field] for field in ('comparison_id', 'clock'))):
                raise ValueError('Public independent calculation identity differs from its source')
            # The reviewer view hashes labels containing a slash. This one exact
            # stopped-target label is public and means that no clock was resolved.
            # Restore it only after the sealed R/calculation identities agree;
            # other unsupported/private labels remain withheld by that view.
            if original['clock'] == 'N/A':
                if calculation.get('status') != 'CANNOT VERIFY' or original.get('numerical_status') != 'CANNOT VERIFY':
                    raise ValueError('Unresolved clock label cannot identify a verified calculation')
                public['clock'] = original['clock']
                public['clock_label_status'] = 'SOURCE_BOUND_LITERAL_STOP_LABEL'
    # N identifiers restart inside each review chunk. Preserve that scope;
    # V identifiers refer to the full study's independent result ordering.
    evidence_bindings = []
    by_row = {row['row_id']: row for row in original_rows}
    for scoped in aggregate['scoped_scientific_checks']:
        clocks = {by_row[r]['fields'].get('clock') for r in scoped['row_ids'] if by_row[r]['fields']}
        relevant = [i for i, c in enumerate(calculations) if c['clock'] in clocks]
        evidence_bindings.append({'review_stage': scoped['review_stage'], 'row_ids': scoped['row_ids'],
            'conditional_evidence': {f'N{i:05d}': {'calculation_index': index} for i, index in enumerate(relevant)},
            'independent_evidence': {f'V{i:05d}': {'verification_index': i} for i in range(len(numeric['results']) if numeric else 0)},
            'trace_evidence_scope': 'Trace text stays private; source and entry hashes plus identity limits are shown separately.'})
    repairs = None
    if repair_stage:
        repair_inputs = artifact_from_stage(study, '06_versioned_repairs', 'inputs.json')
        row_decisions = artifact_from_stage(study, '06_versioned_repairs', 'row-decisions.json')
        manifest = artifact_from_stage(study, '06_versioned_repairs', 'withholding-manifest.json')
        readback = artifact_from_stage(study, '06_versioned_repairs', 'independent-readback.json')
        repair_receipt = repair_stage['result']
        if repair_receipt.get('accepted_scientific_repairs') != 0 or repair_receipt.get('numeric_replacements') != 0:
            raise ValueError('Unsupported scientific or numerical repair release')
        repairs = {'summary': project(repair_receipt, ('status', 'original_rows', 'original_files', 'originals_preserved',
            'row_action_counts', 'withheld_rows', 'withheld_files', 'attempted_withholding_rows',
            'conditional_numeric_candidates', 'unresolved_candidates', 'accepted_scientific_repairs', 'numeric_replacements')),
            'rows': clean(row_decisions), 'independent_readback': project(readback, ('status', 'check_counts', 'scope')),
            'copies': []}
        # Candidate linkage is preserved separately from its already exported
        # conditional calculation, avoiding duplicate prediction arrays.
        candidates = artifact_from_stage(study, '06_versioned_repairs', 'correction-candidates.json')
        repairs['candidates'] = [project(c, ('row_id', 'status', 'reason', 'source_comparison_id', 'clock',
            'calculation_index', 'historical_sample_lineage', 'original_orientation', 'clock_applicability',
            'scientific_status', 'accepted_scientific_repair', 'original_replaced', 'independent_receipt_sha256')) for c in candidates]
        for file in manifest['files']:
            if file.get('withheld_copy'):
                repairs['copies'].append({'artifact_id': file['artifact_id'], 'source_sha256': file['before_sha256'],
                    'copy_sha256': file['after_sha256'], 'operations': [
                        project(op, ('row_id', 'logical_record', 'changes', 'defects')) |
                        {'before': project(op.get('before', {}), row_fields)} for op in file['operations']],
                    'scope': 'Versioned local copy; public packet shows exact row/field differences, not private source paths'})
        bound_rows = {r['row_id']: r for r in repair_inputs['original_rows']}
        for row in original_rows:
            if row['row_id'] in bound_rows:
                row['identity_binding'] = clean(bound_rows[row['row_id']].get('identity_binding'))
    packet = {'schema': 'clockbase-public-corpus-study-v1', 'study_key': key, 'study_id': task['study_id'],
        'modality': task['modality'], 'operational_status': terminal['status'], 'scientific_status': 'CANNOT VERIFY',
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+task['study_id'],
        'terminal_sha256': terminal_hash, 'source_files': sources, 'source_evidence': source_evidence,
        'frozen_source_design': clean(frozen_design), 'trace_identity': trace, 'original_rows': original_rows,
        'original_parse_issues': original_parse_issues,
        'calculations': public_calculations, 'independent_numerics': numeric,
        'review': clean(aggregate), 'review_evidence_bindings': evidence_bindings, 'repairs': repairs,
        'repair_stage_limit': clean(terminal.get('repair_summary')) if repair_stage is None else None,
        'contract_warning_count': terminal.get('contract_warning_count'),
        'failure': clean(terminal.get('model_error') or terminal.get('error')),
        'accepted_scientific_repairs': 0, 'scientific_error_rate': None}
    if measurement is not None:
        packet['measurement_evidence'] = measurement
    _safe(packet)
    name = key.replace(':', '-')
    summary = {'study_key': key, 'study_id': task['study_id'], 'modality': task['modality'],
            'operational_status': terminal['status'], 'reviewed_rows': len(aggregate['rows']),
            'original_rows': len(original_rows), 'calculated_targets': sum(c.get('status') == 'PASS' for c in calculations),
            'unresolved_targets': sum(c.get('status') != 'PASS' for c in calculations),
            'independent_numerical_counts': numeric['counts'] if numeric else None,
            'withheld_rows': repairs['summary']['withheld_rows'] if repairs else 0,
            'repair_stage_status': repairs['summary']['status'] if repairs else 'NOT_EXECUTED',
            'contract_warning_count': terminal.get('contract_warning_count'),
            'scientific_status': 'CANNOT VERIFY', 'artifact_id': name, 'terminal_sha256': terminal_hash}
    if terminal_hash is None:
        if not re.fullmatch(r'[0-9a-f]{64}', report_sha256 or ''):
            raise ValueError('Interrupted corpus evidence requires its exact report hash')
        # This first version supports the observed pre-review interruptions.
        # Preserve later-stage cases privately and refuse release until their
        # partial review/repair mapping can be represented without ambiguity.
        allowed = {'01_source_index', '02_blinded_design', '02b_original_trace',
                   '03_conditional_reruns', '03b_independent_r'}
        receipts = {}
        for path in sorted(study.glob('*/complete.json')):
            stage_name = path.parent.name
            if stage_name not in allowed:
                raise ValueError('Unsupported interrupted corpus stage: ' + stage_name)
            receipt = stage(study, stage_name)
            receipts[stage_name] = {'receipt_sha256': sha(checked(path)),
                'seal_sha256': sha(checked(path.with_name('complete.sha256'))),
                'output_sha256': dict(receipt['outputs'])}
        failed_attempts = []
        for path in sorted(study.glob('*/attempt-*/failure.json')):
            checked(_inside(study, path))
            failure = read(path)
            failed_attempts.append({'stage': path.parent.parent.name, 'attempt': path.parent.name,
                'receipt_sha256': sha(path), 'failure': project(failure, ('status', 'error', 'elapsed_seconds'))})
        intermediate = {'schema': 'clockbase-interrupted-corpus-evidence-v1',
            'study_key': key, 'source_report_sha256': report_sha256,
            'counted_in_terminal_summary': False, 'scientific_status': 'CANNOT VERIFY',
            'accepted_scientific_repairs': 0, 'scientific_error_rate': None,
            'scope': 'Completed-stage outputs saved before controller failure. Excluded from terminal review, calculation and repair totals.',
            'completed_stages': receipts, 'failed_attempts': failed_attempts,
            'calculations': packet['calculations'], 'independent_numerics': packet['independent_numerics'],
            'completed_review_stages': [], 'aggregation': None, 'repairs': None}
        _safe(intermediate)
        packet['intermediate_evidence'] = intermediate
        packet['calculations'] = []
        packet['independent_numerics'] = None
        packet['review'] = {'rows': [], 'scoped_scientific_checks': [], 'design_disagreements': []}
        packet['review_evidence_bindings'] = []
        packet['repairs'] = None
        summary.update(reviewed_rows=0, calculated_targets=0, unresolved_targets=0,
                       independent_numerical_counts=None, withheld_rows=0, repair_stage_status='NOT_EXECUTED')
        summary['intermediate_evidence_summary'] = {
            'sha256': _digest(intermediate), 'source_report_sha256': report_sha256,
            'counted_in_terminal_summary': False, **_intermediate_corpus_counts(intermediate)}
    return packet, summary


def export_corpus(run, output, *, report_path=None, source_root=None, independent_review_path=None):
    """Export a dated operational corpus snapshot with lazy study packets.

    Only an immutable completed report and its terminal/stage receipts are read.
    Original clock outputs are reused; this exporter performs no new audit,
    statistical calculation, model call, or scientific acceptance.
    """
    from .independent_statistics_2026_09_20 import model_verification_view
    from .legacy_analysis_review_2026_09_19 import SCIENTIFIC_CHECKS
    run, output = Path(run).resolve(), Path(output).resolve()
    root = Path(source_root or Path(__file__).resolve().parents[1]).resolve()
    if output.exists() or run.is_relative_to(output):
        raise ValueError('Corpus export requires a new directory outside its source run')
    observed = {}
    helpers = _corpus_projection_helpers(observed)
    checked, read, clean, project, stage, artifact_from_stage = helpers
    config = read(checked(run/'config.json', checked(run/'config.sha256').read_text().strip()))
    queue_receipt = read(run/'queue-receipt.json')
    checked(run/'queue.jsonl.gz', queue_receipt['sha256'])
    with gzip.open(run/'queue.jsonl.gz', 'rt') as handle:
        queue = [json.loads(line) for line in handle]
    if len(queue) != queue_receipt['studies'] or len({r['study_key'] for r in queue}) != len(queue):
        raise ValueError('Corpus queue denominator or identities differ from the frozen receipt')
    queued = {r['study_key']: r for r in queue}
    reports = sorted((run/'reports').glob('report-*.json'))
    if report_path is None and not reports:
        raise ValueError('Corpus export requires a completed dated report')
    report_path = _inside(run, report_path or reports[-1])
    report = read(report_path)
    independent_review = None
    if independent_review_path is not None:
        independent = read(independent_review_path)
        if (independent.get('run', {}).get('report_sha256') != sha(report_path)
                or independent['run'].get('config_sha256') != sha(run/'config.json')
                or independent['run'].get('queue_sha256') != queue_receipt['sha256']):
            raise ValueError('Independent corpus review does not bind this exact run and report')
        independent_review = project(independent, ('status', 'status_scope', 'decision', 'decision_scope',
            'scientific_findings', 'counts', 'limitations', 'recorded_at_utc'))
        independent_review['sha256'] = sha(independent_review_path)
        independent_review['source_report_sha256'] = sha(report_path)
    if report.get('accepted_scientific_repairs') != 0 or report.get('scientific_error_rate') is not None:
        raise ValueError('This corpus export has no scientific-acceptance or error-rate reference gate')
    results = report['results']
    if (len(results) != report['selected_studies'] or len({r['study_key'] for r in results}) != len(results)
            or any(r['study_key'] not in queued for r in results)
            or report['whole_corpus_studies'] != len(queue)):
        raise ValueError('Report study membership is inconsistent with the frozen queue')
    packets, study_summaries, downloads = [], [], {}
    report_digest = sha(report_path)
    for result in sorted(results, key=lambda row: row['study_key']):
        key = result['study_key']
        packet, summary = _project_corpus_study(run/'studies'/key.replace(':', '-'), result, queued[key], root, helpers,
                                             report_sha256=report_digest)
        downloads['studies/'+summary['artifact_id']+'.json'] = packet
        packets.append(packet); study_summaries.append(summary)
    reviewed_ids = {r['row_id'] for p in packets for r in p['review']['rows']}
    summary = {'queued_study_tasks': len(queue), 'attempted_study_tasks': len(results),
        'pending_study_tasks': len(queue)-len(results), 'status_counts': dict(Counter(r['status'] for r in results)),
        'reviewed_unique_rows': len(reviewed_ids),
        'rows_without_completed_review': config['denominators']['parsed_logical_rows']-len(reviewed_ids),
        'conditional_calculations': sum(s['calculated_targets'] for s in study_summaries),
        'unresolved_calculations': sum(s['unresolved_targets'] for s in study_summaries),
        'independent_numerical_counts': {status: sum((s['independent_numerical_counts'] or {}).get(status, 0) for s in study_summaries) for status in STATUS},
        'verified_withheld_rows': sum(s['withheld_rows'] for s in study_summaries if s['repair_stage_status'] == 'PASS'),
        'accepted_scientific_repairs': report['accepted_scientific_repairs'], 'scientific_error_rate': report['scientific_error_rate'],
        'new_clock_executions': 0,
        'scope': 'Operational review coverage under this dated report; agent judgments are not independent scientific reference labels'}
    if (len(reviewed_ids) != report['unique_reviewed_rows'] or summary['conditional_calculations'] != report['calculations']
            or ('unresolved_calculations' in report and summary['unresolved_calculations'] != report['unresolved_calculations'])):
        raise ValueError('Corpus report summary disagrees with source study stages')
    assessment = {'schema': 'clockbase-continuous-corpus-assessment-v1', 'assessment_id': output.name,
        'assessment_type': 'continuous_corpus_audit', 'title': 'Corpus audit: source design, calculations and verified safeguards',
        'recorded_at_utc': report['created_at_utc'], 'scope': report['scope'],
        'overall_status': 'CANNOT VERIFY', 'scientific_status': 'CANNOT VERIFY',
        'status_scope': 'Completed review steps and conditional numerical checks do not establish scientific acceptance',
        'inventory': config['denominators'], 'summary': summary, 'studies': study_summaries,
        'execution': {'provider': config['config']['provider'], 'model': config['config']['api']['model'],
            'concurrency': config['config']['concurrency'],
            'validation_scope': config['config'].get('validation_scope', config['source_scope'])},
        'budget': project(report['budget'], ('max_cost_usd', 'max_calls', 'call_reservations', 'estimated_or_reserved_usd', 'by_status', 'billed_cost_usd')),
        'coverage': clean(report.get('coverage')), 'source_report_sha256': sha(report_path),
        'independent_review': independent_review,
        'interpretation': ['The frozen inventory and this execution have different denominators.',
            'PASS, FAIL and CANNOT VERIFY are retained for each named check. Partial contracts retain controller restrictions.',
            'Original outputs, proposed comparisons, independent R calculations and versioned withholding remain separate.'],
        'limitations': ['This is exposed operational development, not held-out evidence of a scientific error rate.',
            'Historical membership, comparison orientation, measurement applicability and testing-family completeness can remain unresolved.',
            'Provider usage estimates and reserved costs are not an invoice.',
            'Private text is explicitly withheld with a content hash; raw prompts, responses and unpublished source material are not exported.'],
        'artifacts': []}
    intermediate_summaries = [s['intermediate_evidence_summary'] for s in study_summaries if 'intermediate_evidence_summary' in s]
    if intermediate_summaries:
        assessment['intermediate_evidence_summary'] = {
            'study_tasks': len(intermediate_summaries), 'counted_in_terminal_summary': False,
            'scientific_status': 'CANNOT VERIFY',
            **{field: sum(s[field] for s in intermediate_summaries) for field in (
                'calculated_targets', 'unresolved_targets', 'review_row_assignments', 'reviewed_unique_rows',
                'repair_stages', 'completed_stages', 'failed_attempts')},
            'independent_numerical_counts': {status: sum(s['independent_numerical_counts'][status] for s in intermediate_summaries) for status in STATUS}}
    downloads['coverage-public.json'] = {'inventory': assessment['inventory'], 'summary': summary, 'coverage': assessment['coverage']}
    if independent_review is not None:
        downloads['independent-review-public.json'] = independent_review
    downloads['source-ledger-public.json'] = {'source_report_sha256': sha(report_path), 'config_sha256': sha(run/'config.json'),
        'queue_sha256': queue_receipt['sha256'], 'source_implementation_sha256': config['implementation_sha256'],
        'terminal_sha256': {s['study_key']: s['terminal_sha256'] for s in study_summaries},
        'export_implementation_sha256': sha(__file__), 'new_model_calls': 0, 'new_clock_executions': 0}
    output.mkdir(parents=True, exist_ok=False)
    for name, value in downloads.items():
        _safe(value)
        path = output/name; path.parent.mkdir(parents=True, exist_ok=True); write_new(path, value)
        label = Path(name).stem
        assessment['artifacts'].append({'id': label, 'label': label.replace('-', ' '),
            'kind': 'corpus study trace' if name.startswith('studies/') else 'corpus audit evidence',
            'sha256': sha(path), 'download_path': name, 'media_type': 'application/json'})
    _safe(assessment)
    write_new(output/'assessment.json', assessment)
    for path, expected in observed.items():
        if sha(path) != expected:
            raise PipelineIntegrityError('Corpus source changed during public export')
    write_new(output/'artifact-manifest.json', {'schema': 'clockbase-data-quality-artifact-manifest-v1',
        'assessment_id': output.name, 'files': {str(p.relative_to(output)): sha(p) for p in sorted(output.rglob('*')) if p.is_file()}})
    return {'status': 'PASS', 'assessment_id': output.name, 'assessment_sha256': sha(output/'assessment.json'),
            'studies_exported': len(packets), 'scientific_status': 'CANNOT VERIFY', 'new_model_calls': 0}


def _future_terminal_result(terminal, result):
    """The future-hook envelope adds accounting, without changing terminal evidence."""
    envelope = {'execution_status', 'budget', 'pipeline_root_sha256', 'original_paper_2025_eligibility'}
    if (not set(terminal) <= set(result) or set(result) - set(terminal) != envelope or any(result.get(key) != value for key, value in terminal.items())
            or result.get('scientific_acceptance') != 'CANNOT VERIFY' or result.get('accepted_scientific_repairs') != 0
            or result.get('execution_status') != 'COMPLETED' or result.get('original_paper_2025_eligibility') != 'QUARANTINED'):
        raise ValueError('Future comprehensive result differs from its bound terminal or scientific limits')


def _geo_source_text(path):
    """Archived HTTP bodies may be gzip payloads without a .gz filename."""
    payload = Path(path).read_bytes()
    if payload.startswith(b'\x1f\x8b'):
        payload = gzip.decompress(payload)
    return payload.decode('utf-8')


def _future_public_date(date, quote, series_text, ingestion_date=None):
    if (quote not in series_text or not re.fullmatch(r'Public on [A-Z][a-z]{2} [0-9]{2} [0-9]{4}', quote)
            or datetime.strptime(quote, 'Public on %b %d %Y').strftime('%Y-%m-%d') != date
            or (ingestion_date is not None and ingestion_date != date)):
        raise ValueError('First-public date lacks its exact retained source quote')


def _future_producer_binding(packet, summary, result, study_id, prediction_sha256, producer_sha256):
    if (packet['study_id'] != study_id or summary['reviewed_rows'] != result['input_rows']
            or summary['original_rows'] != result['input_rows']):
        raise ValueError('Future comprehensive review is not complete for its declared producer rows')
    if packet['independent_numerics']['predictions_sha256'] != prediction_sha256:
        raise ValueError('Future R verification used different predictions')
    for row in packet['original_rows']:
        binding = row.get('identity_binding') or {}
        if (binding.get('identity_status') != 'PASS' or binding.get('bound_study_key') != packet['study_key']
                or binding.get('producer_stage_sha256') != producer_sha256):
            raise ValueError('Future original row lacks a verified producer-stage identity binding')


def _future_prediction_source(pipeline, clock_stage, checked):
    """Bind RNA or DNAm predictions to one exact sealed stage-output member."""
    root = Path(pipeline).resolve()/'05_clock_predictions'
    path = _inside(root, clock_stage['result']['path'])
    relative = str(path.relative_to(root))
    outputs = clock_stage['outputs']
    aliases = [name for name in outputs if (root/name).resolve() == path]
    if aliases != [relative]:
        raise ValueError('Future predictions are not an exact unique clock-stage output')
    expected = outputs[relative]
    declared = clock_stage['result'].get('prediction_sha256')
    if declared is not None and declared != expected:
        raise ValueError('Future prediction hash differs from sealed clock-stage membership')
    return checked(path, expected)


def export_future_update(version_dir, output, *, plan_path, readback_path, lineage_path, ingestion_readback_path=None):
    """Join an executed future-ingestion pipeline and its comprehensive audit.

    This function exports completed source records only. It never calls a model,
    recomputes clocks, accepts a scientific result, or admits a paper extension.
    """
    version_dir, output = Path(version_dir).resolve(), Path(output).resolve()
    plan_path, readback_path, lineage_path = map(lambda p: Path(p).resolve(), (plan_path, readback_path, lineage_path))
    if output.exists() or version_dir.is_relative_to(output) or output.is_relative_to(version_dir):
        raise ValueError('Future trace requires a new output directory outside its source version')
    observed = {}; helpers = _corpus_projection_helpers(observed)
    checked, read, clean, project, stage, artifact_from_stage = helpers
    version = read(checked(version_dir/'version.json', checked(version_dir/'version.json.sha256').read_text().strip()))
    verification = read(readback_path); lineage = read(lineage_path)
    if (version['version_id'] != version_dir.name or verification.get('version_id') != version['version_id']
            or verification.get('version_sha256') != sha(version_dir/'version.json')
            or verification.get('study_id') != version['study_id'] or verification.get('integrity_status') != 'PASS'
            or verification.get('completed_stage_count') != len(STAGES)
            or verification.get('comprehensive_audit_status') != 'COMPLETED'
            or verification.get('original_paper_2025_eligibility') != 'QUARANTINED'):
        raise ValueError('Future trace readback does not bind a complete quarantined version')
    for relative, expected in verification['artifact_sha256'].items():
        checked(_inside(version_dir, version_dir/relative), expected)
    for name, expected in version['implementation'].items():
        checked(_inside(version_dir/'implementation', version_dir/'implementation'/name), expected)
    for role, source in version['sources'].items():
        checked(source['path'], source['sha256'])
        if version['source_sha256'][role] != source['sha256']:
            raise ValueError('Future source role hash differs from its version')
    for name, expected in lineage['evidence_sha256'].items():
        checked(_inside(lineage_path.parent, lineage_path.parent/name), expected)
    current = [row for row in lineage['development_attempts'] if row['version_id'] == version['version_id']]
    if len(current) != 1 or current[0]['api_calls'] != verification['api_call_count']:
        raise ValueError('Future attempt lineage does not identify the completed run')
    lineage_receipts = []
    for name in lineage['evidence_sha256']:
        if re.fullmatch(r'full-run-readback-[0-9]+\.json', name):
            prior = read(lineage_path.parent/name)
            lineage_receipts.append(project(prior, ('version_id', 'version_sha256', 'integrity_status', 'pipeline_status',
                'comprehensive_audit_status', 'completed_stage_count', 'api_call_count', 'estimated_cost_usd',
                'scientific_status', 'original_paper_2025_eligibility', 'all_version_files_unchanged_on_resume',
                'new_calls_on_resume', 'artifact_sha256')) | {'readback_sha256': sha(lineage_path.parent/name)})
    pipeline = version_dir/'pipeline'; study = version_dir/'comprehensive-audit/study'
    terminal = read(study/'terminal.json')
    result = read(version_dir/'comprehensive-audit-result.json')
    _future_terminal_result(terminal, result)
    cfg = read(study/'config.json')
    if (result['pipeline_root_sha256'] != sha(checked(pipeline/'config.json'))
            or cfg['task']['producer_version_id'] != version['version_id']):
        raise ValueError('Comprehensive hook is bound to a different pipeline version')
    packet, summary = _project_corpus_study(study, terminal, cfg['task'], pipeline, helpers)
    stage(pipeline, '06_coding_agent')
    clock_stage = stage(pipeline, '05_clock_predictions')
    matrix_stage = stage(pipeline, '02_molecular_input')
    predictions_path = _future_prediction_source(pipeline, clock_stage, checked)
    _future_producer_binding(packet, summary, result, version['study_id'], sha(predictions_path), sha(pipeline/'06_coding_agent/complete.json'))
    ingestion = None
    if ingestion_readback_path is not None:
        ingestion_path = Path(ingestion_readback_path).resolve(); ing = read(ingestion_path)
        if (ing.get('status') != 'PASS' or ing.get('study_id') != version['study_id']
                or ing.get('counts_sha256') != sha(checked(matrix_stage['result']['path']))
                or ing.get('all_count_cells_exactly_equal') is not True
                or ing.get('original_paper_2025_eligibility') != 'QUARANTINED'):
            raise ValueError('Molecular acquisition readback is not bound to this source matrix')
        checked(ing['counts_path'], ing['counts_sha256'])
        ingestion = project(ing, ('status', 'scope', 'source_features', 'source_samples', 'source_sample_ids',
            'count_cells_compared', 'all_count_cells_exactly_equal', 'counts_sha256', 'full_hdf5_hash_recomputed',
            'hdf5_size_bytes', 'first_public_date', 'original_paper_2025_eligibility'))
        ingestion['readback_sha256'] = sha(ingestion_path)
    date = version['source_design']['series_first_public_date']
    quote = version['source_design']['source_status_quote']
    series = checked(version['sources']['series']['path'], version['source_sha256']['series'])
    _future_public_date(date, quote, _geo_source_text(series), (ingestion or {}).get('first_public_date'))
    prepare(pipeline, plan_path)
    export(plan_path, output)
    trace_path = output/'audit_trace.json'; trace = json.loads(trace_path.read_text())
    # Preserve the base pipeline export before adding the comprehensive hook.
    shutil.copyfile(trace_path, output/'public/pipeline-audit-trace.json')
    (output/'export-receipt.json').rename(output/'pipeline-export-receipt.json')
    clock_result = clock_stage['result']
    fresh = clock_result.get('clock_predictions_recomputed') is True
    ingestion_view = {'schema': 'clockbase-future-ingestion-public-v1', 'version_id': version['version_id'],
        'source_version_sha256': sha(version_dir/'version.json'), 'study_id': version['study_id'],
        'source_first_public_date': date, 'source_first_public_quote': quote, 'series_sha256': sha(series),
        'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+version['study_id'],
        'original_paper_2025_eligibility': 'QUARANTINED', 'scientific_status': 'CANNOT VERIFY',
        'pipeline_stages_completed': len(STAGES), 'sample_count': trace['data']['sample_count'],
        'clock_count': len(trace['data']['clock_names']), 'clock_names': trace['data']['clock_names'],
        'clock_predictions_recomputed_in_source_run': fresh,
        'fresh_prediction_cells': trace['data']['sample_count']*len(trace['data']['clock_names']) if fresh else 0,
        'clock_prediction_checks': [project(row, ('clock', 'model_sha256', 'numerical_status', 'observed_features_min',
            'required_features', 'max_abs_independent_difference', 'measurement_applicability')) for row in clock_result.get('checks', [])],
        'molecular_acquisition_readback': ingestion, 'comprehensive_summary': summary,
        'accepted_scientific_repairs': 0, 'source_pipeline_report_sha256': sha(version_dir/'pipeline-report-result.json'),
        'comprehensive_result_sha256': sha(version_dir/'comprehensive-audit-result.json'),
        'readback_sha256': sha(readback_path), 'lineage_sha256': sha(lineage_path),
        'completed_run_api_calls': verification['api_call_count'], 'completed_run_estimated_cost_usd': verification['estimated_cost_usd'],
        'all_selected_attempt_api_calls': lineage['total_model_calls_in_selected_attempts'],
        'all_selected_attempt_estimated_cost_usd': lineage['total_estimated_cost_usd_in_selected_attempts'],
        'billed_cost_usd': verification['billed_cost_usd'], 'prior_attempts': clean(lineage['development_attempts']),
        'limits': clean(lineage['incomplete']), 'export_model_calls': 0, 'clocks_rerun_during_export': False}
    trace['future_ingestion'] = ingestion_view
    trace['comprehensive_audit'] = summary | {'artifact_id': 'comprehensive-audit'}
    trace['scope'] = 'Executed future-ingestion pipeline plus source-first comprehensive scientific audit; paper-extension eligibility remains quarantined.'
    trace['models'].append({'provider': 'openai', 'model': version['comprehensive_audit']['model'], 'role': 'comprehensive source-design and sixteen-check review'})
    trace['resource_usage'].update(generation_attempts=verification['api_call_count'],responses_received=verification['api_call_count'],
        estimated_cost_usd=verification['estimated_cost_usd'],billed_cost_usd=verification['billed_cost_usd'])
    public_readback = project(verification, ('integrity_status', 'scientific_status', 'completed_stage_count',
        'verified_source_files', 'verified_implementation_files', 'api_call_count', 'estimated_cost_usd',
        'all_version_files_unchanged_on_resume', 'new_calls_on_resume', 'credential_free_resume', 'version_id', 'version_sha256'))
    public_readback['api_calls'] = [project(row, ('phase', 'receipt_sha256', 'request_sha256', 'reservation_sha256',
        'returned_model', 'prompt_tokens', 'completion_tokens', 'estimated_cost_usd', 'billing_status')) for row in verification['api_calls']]
    files = {'comprehensive-audit.json': packet, 'future-ingestion.json': ingestion_view,
        'update-readback.json': public_readback, 'attempt-lineage.json': {'attempts': clean(lineage['development_attempts']), 'receipts': lineage_receipts}}
    for name, value in files.items():
        _safe(value); write_new(output/'public'/name, value)
    for name in (*files, 'pipeline-audit-trace.json'):
        path = output/'public'/name
        trace['artifacts'].append({'id': Path(name).stem, 'kind': 'future ingestion evidence', 'label': Path(name).stem.replace('-', ' ').capitalize(),
            'sha256': sha(path), 'download_path': 'public/'+name})
    _safe(trace);trace_path.write_text(json.dumps(trace,indent=2,allow_nan=False)+'\n')
    for path, expected in observed.items():
        if sha(path) != expected: raise PipelineIntegrityError('Future source changed during export')
    write_new(output/'export-receipt.json', {'schema': 'clockbase-future-update-export-v1',
        'plan_sha256': sha(plan_path), 'readback_sha256': sha(readback_path), 'lineage_sha256': sha(lineage_path),
        'new_model_calls': 0, 'clocks_rerun': False, 'scientific_status': 'CANNOT VERIFY',
        'output_sha256': {str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file() and p.name != 'export-receipt.json'}})
    return {'status': 'PASS', 'study_id': version['study_id'], 'trace_sha256': sha(trace_path),
        'fresh_source_prediction_cells': ingestion_view['fresh_prediction_cells'], 'comprehensive_calculations': summary['calculated_targets'],
        'scientific_status': 'CANNOT VERIFY', 'paper_2025_eligibility': 'QUARANTINED', 'new_model_calls': 0}



def _planned_clock_snapshot(version_dir, binding):
    """Read a separately declared planned checkpoint; never infer intent from absence."""
    version_dir = Path(version_dir).resolve()
    if (binding.get('schema') != 'clockbase-planned-clock-checkpoint-binding-v1'
            or binding.get('integrity_status') != 'PASS'
            or binding.get('scientific_status') != 'CANNOT VERIFY'
            or type(binding.get('accepted_scientific_repairs')) is not int
            or binding['accepted_scientific_repairs'] != 0
            or binding.get('original_paper_2025_eligibility') != 'QUARANTINED'
            or binding.get('completed_stage_count') != 5 or binding.get('pipeline_status') != 'COMPLETED'
            or binding.get('comprehensive_audit_status') != 'NOT_RUN'):
        raise ValueError('Planned clocks require their separate completed checkpoint binding')
    snapshot = {str(p.relative_to(version_dir)): sha(_inside(version_dir, p))
                for p in sorted(version_dir.rglob('*')) if p.is_file()}
    if not snapshot or snapshot != binding.get('artifact_sha256'):
        raise PipelineIntegrityError('Planned clock version membership or bytes changed')
    observed = {}; checked, read, clean, project, stage, artifact = _corpus_projection_helpers(observed)
    def sealed(name):
        path = version_dir/name
        return read(checked(path, checked(str(path)+'.sha256').read_text().strip()))
    def declared(name):
        value = binding.get(name)
        if not isinstance(value, dict) or set(value) != {'path', 'sha256'}:
            raise ValueError('Planned checkpoint requires a path/hash declaration: '+name)
        return checked(value['path'], value['sha256'])
    version = sealed('version.json'); pipeline = version_dir/'pipeline'; config = read(pipeline/'config.json')
    identity = version.get('identity')
    if (not isinstance(identity, dict) or version.get('version_id') != _digest(identity)
            or any(version.get(key) != value for key, value in identity.items())
            or version['version_id'] != version_dir.name or version['version_id'] != binding.get('version_id')
            or sha(version_dir/'version.json') != binding.get('version_sha256')
            or not re.fullmatch(r'GSE\d+', version.get('study_id', ''))
            or version['study_id'] != binding.get('study_id') or config.get('study') != version['study_id']
            or config.get('modality') != version.get('modality') or config.get('target') != version.get('target')):
        raise ValueError('Planned checkpoint version/configuration identity differs')
    plan_path, scope_path, first_path, freeze_path = (declared(name) for name in
        ('update_plan', 'planned_scope', 'first_source_request', 'protocol_freeze'))
    plan, scope, freeze = read(plan_path), read(scope_path), read(freeze_path)
    entries = plan.get('studies')
    if not isinstance(entries, list) or len(entries) != 1:
        raise ValueError('Planned checkpoint requires its single-study frozen plan')
    entry = entries[0]
    if (plan.get('schema') != 'clockbase-continuous-update-plan-v1'
            or scope.get('allowed_until_sequence') != ['input', 'design', 'clocks']
            or scope.get('prohibited_until') != ['agent', 'audit', 'report']
            or scope.get('input_plan', {}).get('sha256') != sha(plan_path)
            or Path(scope.get('registry', '')).resolve() != version_dir.parents[3]
            or scope.get('study_id') != version['study_id']
            or scope.get('analysis_revision') != version.get('explicit_analysis_revision')
            or scope.get('supersedes_version') != version.get('supersedes_version')
            or any(entry.get(key) != version.get(key) for key in ('study_id', 'modality', 'target', 'supersedes_version'))
            or entry.get('analysis_revision') != version.get('explicit_analysis_revision')
            or any(version['pipeline_options'].get(key) != value for key, value in plan.get('pipeline', {}).items())
            or plan.get('comprehensive_audit') != version.get('comprehensive_audit')
            or plan.get('limits', {}).get('concurrency') != 1 or plan['limits'].get('allow_network') is not False):
        raise ValueError('Preexisting planned stop scope differs from the frozen version')
    sources = version.get('sources', {})
    if (not {'series', 'samples'} <= set(sources)
            or {role: value['sha256'] for role, value in sources.items()} != version.get('source_sha256')):
        raise ValueError('Planned checkpoint source declarations are incomplete')
    for role, source in sources.items(): checked(source['path'], source['sha256'])
    for role, source in entry.get('sources', {}).items():
        if (role not in sources or source['sha256'] != sources[role]['sha256']
                or Path(source['path']).resolve() != Path(sources[role]['path']).resolve()):
            raise ValueError('Planned source role differs from the frozen input plan')
    if not version.get('implementation'):
        raise ValueError('Planned checkpoint lacks its producing implementation snapshot')
    for name, digest in version['implementation'].items():
        checked(_inside(version_dir/'implementation', version_dir/'implementation'/name), digest)
        if name in config.get('implementation', {}) and config['implementation'][name] != digest:
            raise ValueError('Planned pipeline implementation differs from its source version')
    options = version['pipeline_options']
    if (any(config.get(key) != value for key, value in options.items())
            or config.get('source_contract_policy') != 'literal_field_refs_v1'
            or config.get('audit_profile') != 'source_domains_v4'
            or config.get('source_binding_profile') != 'literal_source_v1'
            or config.get('comparison_scope') != 'source_defined_interventions'):
        raise ValueError('Planned checkpoint source-review options differ')
    if (freeze.get('schema') not in {'clockbase-planned-clock-protocol-freeze-v1',
                                    'clockbase-dnam-source-clock-development-freeze-v1'}
            or freeze.get('status') != 'FROZEN' or freeze.get('study_id') != version['study_id']
            or freeze.get('version_id') != version['version_id'] or freeze.get('plan_sha256') != sha(plan_path)
            or freeze.get('first_source_request_sha256') != sha(first_path)
            or freeze.get('allowed_until_sequence') != ['input', 'design', 'clocks']
            or freeze.get('prohibited_until') != ['agent', 'audit', 'report']
            or freeze.get('provider') != 'openai' or freeze.get('model') != (config.get('audit_model') or config['model'])
            or freeze.get('serial_concurrency') != 1 or freeze.get('max_new_reserved_calls') != options['max_calls']
            or Decimal(str(freeze.get('max_new_reserved_usd'))) != Decimal(str(options['max_cost_usd']))
            or freeze.get('scientific_status') != 'CANNOT VERIFY' or freeze.get('accepted_scientific_repairs') != 0
            or any(type(freeze.get(key)) is not int or freeze[key] != 0 for key in
                   ('model_calls_at_freeze', 'new_reservations_at_freeze', 'clock_runs_at_freeze'))):
        raise ValueError('Planned prelaunch freeze differs from its source-only checkpoint')
    frozen_at = datetime.fromisoformat(freeze['utc'].replace('Z', '+00:00'))
    if frozen_at.tzinfo is None:
        raise ValueError('Planned prelaunch timestamp must include its timezone')
    if ({p.name for p in pipeline.iterdir()} - (set(STAGES[:5]) | {'api', 'config.json', '.run.lock'})
            or any(version_dir.glob('pipeline-agent*')) or any(version_dir.glob('pipeline-audit*'))
            or any(version_dir.glob('pipeline-report*')) or any(version_dir.glob('comprehensive-audit*'))
            or any(pipeline.rglob('failure.json'))):
        raise ValueError('Planned clocks contain a failure or unpermitted downstream stage')
    stages = {name: stage(pipeline, name) for name in STAGES[:5]}
    if any(value is None for value in stages.values()):
        raise ValueError('Planned clocks require all five completed producer stages')
    for name, value in stages.items():
        files = {str(p.relative_to(pipeline/name)) for p in (pipeline/name).rglob('*') if p.is_file()}
        if files != set(value['outputs']) | {'complete.json', 'complete.sha256'}:
            raise ValueError('Planned stage has unbound or missing artifacts')
    wrapper = sealed('pipeline-clocks.json')
    request = read(version_dir/'pipeline-clocks-request.json'); result_path = version_dir/'pipeline-clocks-result.json'
    result = read(result_path)
    expected_request = {'phase':'pipeline', 'version_path':str(version_dir/'version.json'),
        'pipeline_root':str(pipeline), 'cache':str(version_dir.parents[3]/'cache'), 'until':'clocks', 'result_path':str(result_path)}
    expected_outputs = {str(p): sha(p) for p in pipeline.rglob('*') if p.is_file() and p.name != '.run.lock'}
    expected_outputs[str(result_path)] = sha(result_path)
    if (request != expected_request or wrapper.get('status') != 'COMPLETED' or wrapper.get('until') != 'clocks'
            or wrapper.get('output_sha256') != expected_outputs or result != stages[STAGES[4]]['result']):
        raise ValueError('Planned clocks lack their exact completed controller invocation')
    for path in version_dir.glob('pipeline-*.json'):
        if re.fullmatch(r'pipeline-(metadata|input|design|preprocess|clocks)\.json', path.name):
            prior = sealed(path.name)
            if prior.get('status') != 'COMPLETED' or prior.get('until') != path.stem.removeprefix('pipeline-'):
                raise ValueError('A failed controller phase cannot become a planned checkpoint')
            for name, digest in prior.get('output_sha256', {}).items(): checked(_inside(version_dir, name), digest)
    def member(name, path):
        target = _inside(pipeline/name, path)
        if stages[name]['outputs'].get(str(target.relative_to(pipeline/name))) != sha(target):
            raise ValueError('Planned source/result is not a completed-stage member')
        return checked(target)
    metadata, clocks, design = (stages[name]['result'] for name in (STAGES[0], STAGES[4], STAGES[2]))
    packet = source_packet(member(STAGES[0], metadata['metadata']), version['study_id'], config['target'],
        context=member(STAGES[0], metadata['study_context']) if metadata.get('study_context') else None,
        sample_source=member(STAGES[0], metadata['source_sample_records']), comparison_scope=config['comparison_scope'])
    for role, key in [('series', 'series_source'), ('samples', 'source_sample_records')]:
        if sha(member(STAGES[0], metadata[key])) != sources[role]['sha256']:
            raise ValueError('Planned physical GEO source differs from its metadata stage')
    source_packet_path = member(STAGES[2], pipeline/STAGES[2]/'attempt-000/source-packet.json')
    if read(source_packet_path) != packet:
        raise ValueError('Planned source request packet differs from the physical source stage')
    transport = read(member(STAGES[2], pipeline/STAGES[2]/'attempt-000/source-transport.json'))
    from .agent_workflow_2026_09_10 import review_messages, parse_design
    from src.pipeline_2026_09_10 import (SOURCE_FIELD_REFERENCE_LIMITS, source_contract_candidate,
        source_contract_feedback, source_contract_instructions, source_contract_visible_packet)
    messages = review_messages(packet, review_profile='source_domains_v4', source_binding_profile='literal_source_v1')
    expected_transport = json.loads(messages[1].content)
    if transport != expected_transport:
        raise ValueError('Planned model-visible packet differs from source-only transport')
    expected_messages = [{'role':'system','content':messages[0].content+source_contract_instructions('literal_field_refs_v1')},
                         {'role':'user','content':messages[1].content}]
    reference_packet = source_contract_visible_packet(packet,transport)
    design_path = member(STAGES[2], pipeline/STAGES[2]/'attempt-000/design.json')
    if read(design_path) != design:
        raise ValueError('Planned design differs from the completed source-stage result')
    ids = [row['sample_id'] for row in packet['samples']]
    if ids != metadata['sample_ids'] or ids != version['source_design']['sample_ids'] or len(ids) != len(set(ids)):
        raise ValueError('Planned source roster differs')
    predictions = member(STAGES[4], clocks['path'])
    with predictions.open(newline='') as handle: header = next(csv.reader(handle), None)
    if (not header or header != ['sample_id']+clocks['clocks'] or len(header) != len(set(header))
            or not clocks['clocks'] or any(not isinstance(name,str) or not name for name in clocks['clocks'])):
        raise ValueError('Planned prediction columns are empty or duplicated')
    frame = pd.read_csv(predictions).set_index('sample_id', verify_integrity=True)
    if (list(frame.index) != ids or not np.isfinite(frame.to_numpy(float)).all()
            or clocks.get('sample_count') != len(ids) or clocks.get('clock_predictions_recomputed') is not True):
        raise ValueError('Planned predictions differ from saved source samples or finite clocks')
    api = pipeline/'api'; slots = sorted(p for p in api.iterdir() if p.is_dir())
    if (not slots or {p.name for p in api.iterdir()} != {p.name for p in slots}
            or [p.name for p in slots] != ['design-%03d'%i for i in range(len(slots))]
            or len(slots) > options['max_calls'] or len(slots) > options['audit_validation_retries']+1):
        raise ValueError('Planned checkpoint has non-source or excessive model slots')
    calls=[]; history=[]; previous_request=previous_response=None
    for index, slot in enumerate(slots):
        if {p.name for p in slot.iterdir()} != {'request.json','reservation.json','response.json','receipt.json'}:
            raise ValueError('Planned source slot has missing or uncertain call evidence')
        request, reservation, receipt = (read(slot/name) for name in ('request.json','reservation.json','receipt.json'))
        response = read(checked(slot/'response.json', receipt['response_sha256']))
        if (isinstance(reservation.get('started_unix'), bool)
                or not math.isfinite(float(reservation.get('started_unix', float('nan'))))
                or float(reservation['started_unix']) < frozen_at.timestamp()):
            raise ValueError('Planned scope was not frozen before this source call')
        model = config.get('audit_model') or config['model']
        settings = {'model':model,'max_completion_tokens':options['max_output_tokens'],
            'reasoning_effort':'low','store':False,'response_format':{'type':'json_object'}}
        if ({key:value for key,value in request.items() if key!='messages'} != settings
                or model != reservation.get('model') or model != receipt.get('requested_model')
                or model != receipt.get('returned_model') or response.get('model') != model
                or receipt.get('usage') != response.get('usage') or reservation.get('automatic_retries') != 0
                or reservation.get('concurrency') != 1 or len(response.get('choices',[])) != 1
                or response['choices'][0].get('finish_reason') != 'stop'
                or not isinstance(response['choices'][0].get('message',{}).get('content'),str)):
            raise ValueError('Planned source API model, settings, completion or usage differs')
        if index == 0:
            if sha(slot/'request.json') != sha(first_path) or request['messages'] != expected_messages:
                raise ValueError('Planned first source request differs from its prelaunch capture')
        else:
            feedback = read(member(STAGES[2], pipeline/STAGES[2]/('attempt-000/response-validation/%03d-validation.json'%(index-1))))
            if feedback.get('validation_status') != 'FAIL' or request != _failed_retry_request(previous_request,previous_response,feedback):
                raise ValueError('Planned source correction differs from its retained completed-response feedback')
        raw = json.dumps(request,ensure_ascii=False).encode()
        usage = receipt['usage']
        if any(type(usage.get(key)) is not int or usage[key]<0 for key in ('prompt_tokens','completion_tokens')):
            raise ValueError('Planned source usage is invalid')
        cost=(Decimal(usage['prompt_tokens'])*Decimal(str(options['audit_input_usd_per_million']))
              +Decimal(usage['completion_tokens'])*Decimal(str(options['audit_output_usd_per_million'])))/Decimal(1000000)
        reserved=(Decimal(len(raw)+100*len(request['messages']))*Decimal(str(options['audit_input_usd_per_million']))
                  +Decimal(options['max_output_tokens'])*Decimal(str(options['audit_output_usd_per_million'])))/Decimal(1000000)
        if (len(raw)>options['max_prompt_bytes'] or abs(Decimal(str(receipt['estimated_cost_usd']))-cost)>Decimal('1e-12')
                or abs(Decimal(str(reservation['reserved_usd']))-reserved)>Decimal('1e-12')):
            raise ValueError('Planned source cost/reservation differs from frozen rates')
        calls.append({'slot':slot.name,'phase':'source_design',**{name+'_sha256':sha(slot/(name+'.json'))
            for name in ('request','reservation','response','receipt')},
            **project(receipt,('requested_model','returned_model','usage','estimated_cost_usd','billed_cost_usd','elapsed_seconds'))})
        raw_text = response['choices'][0]['message']['content']
        original = member(STAGES[2], pipeline/STAGES[2]/('attempt-000/contract-adjustments/%03d-original-response.txt'%index))
        if original.read_text() != raw_text:
            raise ValueError('Planned source contract original differs from its API response')
        attempt = pipeline/STAGES[2]/'attempt-000'
        adjustment_path = member(STAGES[2],attempt/('contract-adjustments/%03d-receipt.json'%index))
        adjustment = read(adjustment_path)
        expected_adjustment = {'schema':'clockbase-source-contract-adjustment-v1','policy':'literal_field_refs_v1',
            'attempt':index,'original_response_sha256':sha(original),'changes':[],'scientific_status':'CANNOT VERIFY',
            'accepted_result':False,'source_packet_sha256':sha(source_packet_path),
            'source_transport_sha256':sha(attempt/'source-transport.json'),
            'implementation_sha256':config['implementation']['src/pipeline_2026_09_10.py'],
            'normalization_limits':SOURCE_FIELD_REFERENCE_LIMITS,
            'scope':'Exact selected-field materialization and conservative advisory statuses; strict evidence and scientific checks still apply'}
        rejected = None
        candidate_path = attempt/('contract-adjustments/%03d-candidate.json'%index)
        try:
            candidate,changes = source_contract_candidate(raw_text,reference_packet,'literal_field_refs_v1')
        except (ValueError,TypeError,KeyError) as error:
            if candidate_path.exists():raise ValueError('Rejected materialization has an invented saved candidate') from error
            rejected = error
        else:
            candidate_path = member(STAGES[2],candidate_path)
            expected_adjustment.update(changes=changes,candidate_sha256=sha(candidate_path))
            if _corpus_record_hash(read(candidate_path)) != _corpus_record_hash(candidate):
                raise ValueError('Planned materialized source candidate differs')
            try:
                parsed = parse_design(json.dumps(candidate,ensure_ascii=False,allow_nan=False),packet,review_profile='source_domains_v4')
            except (ValueError,TypeError,KeyError) as error:rejected = error
        if rejected is not None:
            if index == len(slots)-1:raise ValueError('Planned final source response fails the strict source contract') from rejected
            expected_adjustment.update(contract_status='REJECTED',error=f'{type(rejected).__name__}: {rejected}')
        else:
            if index != len(slots)-1 or _corpus_record_hash(parsed) != _corpus_record_hash(design):
                raise ValueError('Planned accepted source design or stop after acceptance differs')
            expected_adjustment['contract_status']='PARTIALLY_ADMISSIBLE' if changes else 'UNCHANGED_VALID'
            if _corpus_record_hash(read(member(STAGES[2],attempt/'source-contract-admissibility.json'))) != _corpus_record_hash(expected_adjustment):
                raise ValueError('Planned final source admissibility differs')
        if _corpus_record_hash(adjustment) != _corpus_record_hash(expected_adjustment):
            raise ValueError('Planned source adjustment provenance or outcome differs')
        validation_sha, error, validation_status = None, None, 'NOT_RECORDED'
        if options['audit_validation_retries']:
            text_path = member(STAGES[2], pipeline/STAGES[2]/('attempt-000/response-validation/%03d-response.txt'%index))
            feedback_path = member(STAGES[2], pipeline/STAGES[2]/('attempt-000/response-validation/%03d-validation.json'%index))
            feedback = read(feedback_path)
            validation_status = 'PASS' if index == len(slots)-1 else 'FAIL'
            if (text_path.read_text() != raw_text or feedback.get('validation_status') != validation_status
                    or (validation_status == 'PASS' and feedback.get('scientific_status') != 'CANNOT VERIFY')):
                raise ValueError('Planned retained source response-validation sequence differs')
            if validation_status == 'FAIL':
                expected_feedback = source_contract_feedback(raw_text,reference_packet,rejected,'literal_field_refs_v1',schema=feedback.get('schema'))
                if _corpus_record_hash(feedback) != _corpus_record_hash(expected_feedback):
                    raise ValueError('Planned rejected source feedback differs from its strict contract replay')
                error = feedback.get('error')
                if not isinstance(error,str) or not error:
                    raise ValueError('Planned rejected source response lacks its recorded error')
                _safe(error)
            elif feedback != {'validation_status':'PASS','scientific_status':'CANNOT VERIFY',
                              'scope':'Response schema and source checks; no scientific repair accepted'}:
                raise ValueError('Planned final validation receipt differs')
            validation_sha = sha(feedback_path)
        elif (pipeline/STAGES[2]/'attempt-000/response-validation').exists():
            raise ValueError('Planned zero-retry source call has conflicting validation history')
        history.append({'slot':slot.name,'recorded_validation_status':validation_status,
            'response_text_sha256':sha(original),'validation_sha256':validation_sha,
            'validation_error':error,'retry_transmitted':index<len(slots)-1,
            'scientific_status':'CANNOT VERIFY','accepted_scientific_repairs':0})
        previous_request,previous_response=request,response
    expected_validation = {str(pipeline/STAGES[2]/('attempt-000/response-validation/%03d-%s'%(index,suffix)))
        for index in range(len(slots)) for suffix in ('response.txt','validation.json')} if options['audit_validation_retries'] else set()
    if {str(p) for p in (pipeline/STAGES[2]/'attempt-000/response-validation').rglob('*') if p.is_file()} != expected_validation:
        raise ValueError('Planned validation history contains unrepresented artifacts')
    if read(member(STAGES[2],pipeline/STAGES[2]/'attempt-000/model-design.json')) != json.loads(previous_response['choices'][0]['message']['content']):
        raise ValueError('Planned advisory model design differs from its final API response')
    cost=float(sum((Decimal(str(call['estimated_cost_usd'])) for call in calls),Decimal(0)))
    if (len(calls)!=binding.get('api_call_count') or abs(cost-binding.get('estimated_cost_usd',-1))>1e-12
            or cost>options['max_cost_usd']):
        raise ValueError('Planned source total accounting differs')
    readback_paths, readback_declarations = {}, {}
    readbacks = binding.get('independent_readbacks', {})
    if not isinstance(readbacks,dict) or set(readbacks)-{'independent-source-readback','independent-clock-readback'}:
        raise ValueError('Planned independent readback roles are unknown')
    for name, declaration in readbacks.items():
        if (not isinstance(declaration,dict) or set(declaration)!={'path','sha256','source_receipt'}
                or not isinstance(declaration['source_receipt'],dict)
                or set(declaration['source_receipt'])!={'path','sha256'}):
            raise ValueError('Planned independent readback lacks its original receipt binding')
        original = declaration['source_receipt']; checked(original['path'],original['sha256'])
        path = checked(declaration['path'],declaration['sha256']); record = read(path); _safe(record)
        kind = 'source_design_sha256' if name=='independent-source-readback' else 'source_predictions_sha256'
        expected = sha(design_path if name=='independent-source-readback' else predictions)
        if (not isinstance(record,dict) or not isinstance(record.get('schema'),str) or not record['schema']
                or record.get('study_id')!=version['study_id'] or record.get('version_id')!=version['version_id']
                or record.get('source_version_sha256')!=sha(version_dir/'version.json')
                or record.get('source_receipt_sha256')!=original['sha256'] or record.get(kind)!=expected
                or record.get('scientific_status',record.get('overall_scientific_status'))!='CANNOT VERIFY'
                or ('overall_scientific_status' in record and record['overall_scientific_status']!='CANNOT VERIFY')
                or type(record.get('accepted_scientific_repairs')) is not int or record['accepted_scientific_repairs']!=0):
            raise ValueError('Planned independent readback differs from the same source/prediction record or scientific scope')
        readback_paths[name]=path
        readback_declarations[name]={'sha256':sha(path),'source_receipt_sha256':original['sha256'],'schema':record['schema']}
    date, quote = (version['source_design'][key] for key in ('series_first_public_date','source_status_quote'))
    _future_public_date(date,quote,_geo_source_text(Path(sources['series']['path'])))
    value={'schema':'clockbase-planned-clock-checkpoint-public-v1','outcome':'planned_clocks_checkpoint',
        'study_id':version['study_id'],'version_id':version['version_id'],'run_key':'planned-clocks-'+version['version_id'][:20],
        'planned_until':'clocks','completed_stages':list(STAGES[:5]),'unexecuted_pipeline_stages':list(STAGES[5:]),
        'comprehensive_status':'NOT_RUN','source_version_sha256':sha(version_dir/'version.json'),
        'source_snapshot_sha256':_digest(snapshot),'planned_scope_sha256':sha(scope_path),'update_plan_sha256':sha(plan_path),
        'protocol_freeze_sha256':sha(freeze_path),
        'controller_completion_sha256':sha(version_dir/'pipeline-clocks.json'),
        'controller_request_sha256':sha(version_dir/'pipeline-clocks-request.json'),
        'controller_result_sha256':sha(result_path),'first_source_request_sha256':sha(first_path),
        'source_design_sha256':sha(design_path),'source_config_sha256':sha(pipeline/'config.json'),
        'source_first_public_date':date,'source_first_public_quote':quote,'series_sha256':sources['series']['sha256'],
        'source_url':packet['source_url'],'sample_ids':ids,'sample_count':len(ids),'clock_names':clocks['clocks'],
        'fresh_prediction_cells':int(frame.size),'predictions_sha256':sha(predictions),'api_calls':calls,
        'source_response_history':history,
        'independent_readbacks':readback_declarations,
        'source_budget':{'model':config.get('audit_model') or config['model'],'max_calls':options['max_calls'],
            'max_cost_usd':options['max_cost_usd'],'input_usd_per_million':options['audit_input_usd_per_million'],
            'output_usd_per_million':options['audit_output_usd_per_million'],
            'validation_retries':options['audit_validation_retries']},
        'estimated_cost_usd':cost,'billed_cost_usd':None,'scientific_status':'CANNOT VERIFY',
        'accepted_scientific_repairs':0,'original_paper_2025_eligibility':'QUARANTINED','scientific_error_rate':None,
        'export_model_calls':0,'clocks_rerun_during_export':False,
        'scope':'Planned completion through saved clock predictions. Pipeline coding, statistical analysis, audit/repair and report stages and comprehensive audit were not executed. Separately downloaded independent readbacks retain their own scope.'}
    return snapshot,version,config,stages,value,packet,predictions,design_path,observed,readback_paths


def export_planned_clock_checkpoint(version_dir, output, *, readback_path, plan_path):
    """Additive planned checkpoint export; old failure profiles remain unchanged."""
    version_dir,output=Path(version_dir).resolve(),Path(output).resolve()
    if output.exists() or output.is_relative_to(version_dir) or version_dir.is_relative_to(output):
        raise ValueError('Planned checkpoint requires a new output outside its source')
    readback_path,plan_path=Path(readback_path).resolve(),Path(plan_path).resolve()
    binding=json.loads(readback_path.read_text()); binding_sha=sha(readback_path)
    snapshot,version,config,stages,value,packet,predictions,design_path,observed,readback_paths=_planned_clock_snapshot(version_dir,binding)
    export_plan={'schema':'clockbase-planned-clock-export-plan-v1','binding_sha256':binding_sha,
        'source_snapshot_sha256':_digest(snapshot),'export_identity':_export_identity()}
    if plan_path.exists():
        if json.loads(plan_path.read_text()) != export_plan:
            raise PipelineIntegrityError('Planned checkpoint export plan differs')
    else: write_new(plan_path,export_plan)
    value['binding_sha256']=binding_sha
    design=json.loads(design_path.read_text())
    if set(design.get('audit_domains',{}))!=set(AUDIT_DOMAINS):
        raise ValueError('Planned source design lacks its complete advisory domain assessment')
    domains=[]
    for name in AUDIT_DOMAINS:
        assessment={key:design['audit_domains'][name][key] for key in
            ('status','summary','evidence','unresolved','proposed_repair','verification_needed')}
        assessment['evidence']=[{key:item[key] for key in ('sample_id','source_id','field','quote','source_url') if key in item}
                                for item in assessment['evidence']]
        domains.append({'id':name,'title':name.replace('_',' ').capitalize(),'model_assessment':assessment,
            'verification':{'status':'CANNOT VERIFY','scope':'Completed source artifact bindings only; advisory model assessment is not scientific acceptance','checks':[]}})
    trace={'schema':'clockbase-audit-trace-v1','study_id':version['study_id'],'run_id':'sha256:'+_digest(snapshot),
        'recorded_at_utc':datetime.now(timezone.utc).isoformat().replace('+00:00','Z'),'evaluation_split':'development',
        'scope':value['scope'],'models':[{'provider':'openai','model':config.get('audit_model') or config['model'],
            'role':'source-design reviewer; no coding or comprehensive audit executed'}],
        'overall':{'status':'CANNOT VERIFY','reason':'The predeclared endpoint was clock prediction; later scientific analysis was not executed'},
        'domains':domains, 'comparisons':[],'covariates':[], 'artifacts':[],
        'events':[{'stage':'planned clock checkpoint','status':'CANNOT VERIFY',
            'summary':'Five producer stages completed at the planned endpoint. Coding, audit and report were not executed.'}],
        'resource_usage':{'generation_attempts':len(value['api_calls']),'responses_received':len(value['api_calls']),
            'estimated_cost_usd':value['estimated_cost_usd'],'billed_cost_usd':None},
        'limitations':['Saved predictions and a completed source contract do not establish biological independence, clock applicability or a treatment effect.',
            'The source-design download and domain assessments preserve the archived advisory reasoning; no statistical comparison, repair or comprehensive audit is represented.',
            'Independent source and clock readback downloads, when present, describe separate verification and its limits; they are not completed pipeline audit or report stages.',
            'This development checkpoint has no scientific error-rate estimate and accepts no repairs.'],
        'planned_checkpoint':value}
    source=_failed_public_sources(packet)
    manifest={'files':snapshot,'study_id':version['study_id'],'source_config_sha256':sha(version_dir/'pipeline/config.json'),
        'completed_stages':list(STAGES[:5]),'planned_until':'clocks'}
    public={'planned-checkpoint':value,'source-records':source,'stage-manifest':manifest}
    for item in public.values():_safe(item)
    _safe(json.loads(design_path.read_text())); _safe(trace)
    output.mkdir(parents=True);(output/'public').mkdir()
    for name,item in public.items():write_new(output/'public'/(name+'.json'),item)
    shutil.copyfile(predictions,output/'public/sample-predictions.csv')
    shutil.copyfile(design_path,output/'public/source-design.json')
    for name,path in readback_paths.items():shutil.copyfile(path,output/'public'/(name+'.json'))
    for name,suffix in [('planned-checkpoint','.json'),('source-records','.json'),('stage-manifest','.json'),
                        ('sample-predictions','.csv'),('source-design','.json')]+[(name,'.json') for name in readback_paths]:
        path=output/'public'/(name+suffix)
        trace['artifacts'].append({'id':name,'kind':'planned checkpoint evidence','label':name.replace('-',' ').capitalize(),
            'sha256':sha(path),'download_path':'public/'+path.name})
    _safe(trace);write_new(output/'audit_trace.json',trace)
    if (sha(readback_path)!=binding_sha or _planned_clock_snapshot(version_dir,binding)[0]!=snapshot
            or any(sha(path)!=digest for path,digest in observed.items())):
        raise PipelineIntegrityError('Planned checkpoint sources changed during export')
    write_new(output/'export-receipt.json',{'schema':'clockbase-planned-clock-export-v1','binding_sha256':binding_sha,
        'source_snapshot_sha256':_digest(snapshot),'new_model_calls':0,'clocks_rerun':False,
        'output_sha256':{str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file()}})
    return {'status':'PASS','trace_sha256':sha(output/'audit_trace.json'),'planned_until':'clocks',
        'completed_stages':5,'scientific_status':'CANNOT VERIFY','accepted_scientific_repairs':0}


def _partial_future_snapshot(version_dir, binding):
    """Validate two explicit partial profiles without relaxing completed exporters."""
    version_dir = Path(version_dir).resolve()
    if (binding.get('schema') != 'clockbase-partial-update-binding-v1'
            or binding.get('integrity_status') != 'PASS'
            or binding.get('scientific_status') != 'CANNOT VERIFY'
            or binding.get('original_paper_2025_eligibility') != 'QUARANTINED'):
        raise ValueError('Partial update requires a quarantined source binding')
    expected = binding.get('artifact_sha256')
    if not isinstance(expected, dict) or not expected:
        raise ValueError('Partial update lacks a complete version manifest')
    snapshot = {str(p.relative_to(version_dir)): sha(_inside(version_dir, p))
                for p in sorted(version_dir.rglob('*')) if p.is_file()}
    if snapshot != expected:
        raise PipelineIntegrityError('Partial update version membership or bytes changed')
    observed = {}; checked, read, clean, project, stage, artifact = _corpus_projection_helpers(observed)
    def sealed(relative, sidecar):
        return read(checked(version_dir/relative, checked(version_dir/sidecar).read_text().strip()))
    version = sealed('version.json', 'version.json.sha256')
    pipeline = version_dir/'pipeline'
    config = read(pipeline/'config.json')
    if (version.get('version_id') != version_dir.name or version['version_id'] != binding.get('version_id')
            or not re.fullmatch(r'[0-9a-f]{64}', version['version_id'])
            or sha(version_dir/'version.json') != binding.get('version_sha256')
            or version.get('study_id') != binding.get('study_id') or config.get('study') != version['study_id']
            or not re.fullmatch(r'GSE\d+', version['study_id'])):
        raise ValueError('Partial update study or version identity differs')
    for role, source in version['sources'].items():
        checked(source['path'], source['sha256'])
        if version['source_sha256'].get(role) != source['sha256']:
            raise ValueError('Partial update source-role identity differs')
    for name, digest in version['implementation'].items():
        checked(_inside(version_dir/'implementation', version_dir/'implementation'/name), digest)
    completed = [name for name in STAGES if (pipeline/name/'complete.json').exists()]
    if completed not in [list(STAGES[:5]), list(STAGES)]:
        raise ValueError('Unsupported partial update stage profile')
    stages = {name: stage(pipeline, name) for name in completed}
    if len(stages) != binding.get('completed_stage_count'):
        raise ValueError('Partial update completed-stage count differs from readback')
    wrapper = sealed('pipeline-report.json', 'pipeline-report.json.sha256')
    contracts = []
    if len(completed) == 5:
        outcome, phase, failed_stage = 'coding_stop', 'pipeline', '06_coding_agent'
        if (wrapper.get('status') != 'FAILED' or binding.get('pipeline_status') != 'FAILED'
                or binding.get('comprehensive_audit_status') != 'MISSING'
                or any((version_dir/name).exists() for name in ('pipeline-report-result.json',
                    'comprehensive-audit', 'comprehensive-audit.json', 'comprehensive-audit-result.json'))
                or any((pipeline/name).exists() for name in STAGES[6:])
                or (pipeline/failed_stage/'complete.sha256').exists()):
            raise ValueError('Coding stop contains a report or downstream evidence')
        failure_path = pipeline/failed_stage/'attempt-000/failure.json'
        if sorted(pipeline.rglob('failure.json')) != [failure_path]:
            raise ValueError('Coding stop requires its single retained failure')
        comprehensive_status = 'MISSING'
    else:
        outcome, phase, failed_stage = 'producer_complete_comprehensive_unresolved', 'comprehensive', '02_blinded_design'
        if wrapper.get('status') != 'COMPLETED' or binding.get('pipeline_status') != 'COMPLETED':
            raise ValueError('Unresolved comprehensive audit requires a completed producer report')
        _snapshot(pipeline)
        report = read(version_dir/'pipeline-report-result.json')
        if (report.get('execution_status') != 'PASS' or report.get('accepted_results') != 0
                or report.get('scientific_status') != 'CANNOT VERIFY'):
            raise ValueError('Producer report exceeds the partial scientific scope')
        audit = version_dir/'comprehensive-audit/study'
        hook = sealed('comprehensive-audit.json', 'comprehensive-audit.json.sha256')
        terminal = sealed('comprehensive-audit/study/terminal.json', 'comprehensive-audit/study/terminal.sha256')
        result = read(version_dir/'comprehensive-audit-result.json')
        envelope = {'execution_status', 'budget', 'pipeline_root_sha256', 'original_paper_2025_eligibility'}
        if (set(result)-set(terminal) != envelope or not set(terminal) <= set(result)
                or any(result[key] != value for key, value in terminal.items())
                or terminal.get('status') != 'UNRESOLVED' or result.get('execution_status') != 'UNRESOLVED'
                or hook.get('status') != 'UNRESOLVED' or hook.get('result') != result
                or binding.get('comprehensive_audit_status') != 'UNRESOLVED'
                or terminal.get('scientific_acceptance') != 'CANNOT VERIFY'
                or terminal.get('accepted_scientific_repairs') != 0 or terminal.get('reviewed_rows') != 0
                or terminal.get('calculations') != 0 or terminal.get('study_id') != version['study_id']
                or result.get('original_paper_2025_eligibility') != 'QUARANTINED'
                or result.get('pipeline_root_sha256') != sha(pipeline/'config.json')
                or read(audit/'config.json')['task']['producer_version_id'] != version['version_id']):
            raise ValueError('Unresolved comprehensive terminal is not bound to this producer')
        for relative, digest in terminal['output_sha256'].items():
            checked(_inside(audit, audit/relative), digest)
        if (stage(audit, '01_source_index') is None
                or {p.name for p in audit.iterdir() if p.is_dir()} != {'01_source_index', '02_blinded_design', 'api'}
                or any((audit/failed_stage/name).exists() for name in ('complete.json', 'complete.sha256'))):
            raise ValueError('Unresolved comprehensive audit contains later completed stages')
        failure_path = audit/failed_stage/'attempt-000/failure.json'
        if sorted(audit.rglob('failure.json')) != [failure_path]:
            raise ValueError('Unresolved audit requires its single retained failure')
        for path in sorted(failure_path.parent.glob('contract-*.json')):
            contract = read(path)
            if contract.get('status') != 'REJECTED_RESPONSE_CONTRACT' or contract.get('attempt') != len(contracts):
                raise ValueError('Unresolved audit contract attempts are incomplete or accepted')
            contracts.append(project(contract, ('status', 'attempt', 'error', 'scientific_status')) | {'sha256': sha(path)})
        if len(contracts) != 2:
            raise ValueError('Unresolved audit requires its two rejected contracts')
        comprehensive_status = 'UNRESOLVED'
    failure = read(failure_path)
    if failure.get('status') != 'CANNOT VERIFY' or not isinstance(failure.get('error'), str) or not failure['error']:
        raise ValueError('Partial update lacks its saved failure')
    _finite_nonnegative(failure.get('elapsed_seconds'), 'partial failure elapsed seconds')
    if {p.name for p in failure_path.parent.parent.iterdir()} != {'attempt-000'}:
        raise ValueError('Partial update has unsupported failure attempts')
    for name in ('pipeline-report', 'comprehensive-audit'):
        path = version_dir/(name+'.json')
        if path.exists():
            for target, digest in read(path).get('output_sha256', {}).items():
                checked(_inside(version_dir, target), digest)
    calls = []
    slots = sorted(p.parent for p in version_dir.rglob('reservation.json'))
    if set(slots) != {p.parent for p in version_dir.rglob('response.json')}:
        raise ValueError('Partial update has unaccounted API responses')
    for slot in slots:
        if slot.parent.name != 'api':
            raise ValueError('Partial update has an unsupported API slot')
        request, reservation, receipt = (read(slot/name) for name in ('request.json', 'reservation.json', 'receipt.json'))
        response = read(checked(slot/'response.json', receipt['response_sha256']))
        if (request['model'] != reservation['model'] or request['model'] != receipt['requested_model']
                or response['model'] != receipt['returned_model'] or response['model'] != request['model']
                or receipt['usage'] != response['usage']):
            raise ValueError('Partial update API identity or usage differs')
        _finite_nonnegative(receipt['estimated_cost_usd'], 'API estimate')
        calls.append({'phase': 'pipeline' if slot.is_relative_to(pipeline) else 'comprehensive',
            'slot': slot.name, 'receipt_sha256': sha(slot/'receipt.json'), 'response_sha256': sha(slot/'response.json'),
            'request_sha256': sha(slot/'request.json'), 'reservation_sha256': sha(slot/'reservation.json'),
            **project(receipt, ('requested_model', 'returned_model', 'estimated_cost_usd', 'billed_cost_usd', 'elapsed_seconds', 'usage'))})
    cost = float(sum((Decimal(str(c['estimated_cost_usd'])) for c in calls), Decimal(0)))
    if len(calls) != binding.get('api_call_count') or abs(cost-binding.get('estimated_cost_usd', -1)) > 1e-12:
        raise ValueError('Partial update API accounting differs from bound readback')
    return snapshot, version, config, stages, {'outcome': outcome, 'comprehensive_status': comprehensive_status,
        'failure': {'phase': phase, 'stage': failed_stage, **project(failure, ('status', 'error', 'elapsed_seconds')),
                    'sha256': sha(failure_path)}, 'contract_attempts': contracts, 'api_calls': calls,
        'estimated_cost_usd': cost}, observed


def export_partial_future_update(version_dir, output, *, readback_path, plan_path):
    """Export a coding stop or an unresolved comprehensive audit from pinned files."""
    version_dir, output = Path(version_dir).resolve(), Path(output).resolve()
    if output.exists() or output.is_relative_to(version_dir) or version_dir.is_relative_to(output):
        raise ValueError('Partial update requires a new output directory outside its source')
    binding = json.loads(Path(readback_path).read_text())
    snapshot, version, config, stages, partial, observed = _partial_future_snapshot(version_dir, binding)
    pipeline = version_dir/'pipeline'
    metadata, clocks = (stages[name]['result'] for name in ('01_metadata', '05_clock_predictions'))
    def member(stage_name, value):
        path = _inside(pipeline/stage_name, value)
        relative = str(path.relative_to(pipeline/stage_name))
        if stages[stage_name]['outputs'].get(relative) != sha(path):
            raise ValueError('Partial source is not an exact stage-output member')
        return path
    packet = source_packet(member('01_metadata', metadata['metadata']), version['study_id'], config['target'],
        context=member('01_metadata', metadata['study_context']) if metadata.get('study_context') else None)
    source = _failed_public_sources(packet)
    ids = [row['sample_id'] for row in packet['samples']]
    if ids != metadata['sample_ids'] or ids != version['source_design']['sample_ids'] or len(ids) != len(set(ids)):
        raise ValueError('Partial source roster differs from metadata and version')
    predictions = member('05_clock_predictions', clocks['path'])
    frame = pd.read_csv(predictions).set_index('sample_id', verify_integrity=True)
    if (set(frame.index) != set(ids) or list(frame.columns) != clocks['clocks']
            or not np.isfinite(frame.to_numpy(float)).all() or clocks['sample_count'] != len(ids)
            or clocks.get('clock_predictions_recomputed') is not True):
        raise ValueError('Partial saved predictions do not bind the source samples and fresh clocks')
    series = Path(version['sources']['series']['path'])
    date, quote = (version['source_design'][key] for key in ('series_first_public_date', 'source_status_quote'))
    _future_public_date(date, quote, _geo_source_text(series))
    partial.update(schema='clockbase-partial-future-update-public-v1', study_id=version['study_id'],
        version_id=version['version_id'], source_version_sha256=sha(version_dir/'version.json'),
        source_first_public_date=date, source_first_public_quote=quote, series_sha256=sha(series),
        source_url=source['source_url'], completed_stages=list(stages),
        missing_pipeline_stages=[name for name in STAGES if name not in stages],
        sample_ids=ids, sample_count=len(ids), clock_names=clocks['clocks'],
        fresh_prediction_cells=int(frame.size), predictions_sha256=sha(predictions),
        source_completion_readback_sha256=binding['source_completion_readback_sha256'],
        binding_sha256=sha(readback_path), source_snapshot_sha256=_digest(snapshot),
        scientific_status='CANNOT VERIFY', original_paper_2025_eligibility='QUARANTINED',
        accepted_scientific_repairs=0, scientific_error_rate=None, billed_cost_usd=None,
        export_model_calls=0, clocks_rerun_during_export=False,
        numerical_scope='Export arithmetic on completed producer rows is separate from the saved independent R readback; no inference is accepted.')
    if len(stages) == 8:
        prepare(pipeline, plan_path)
        export(plan_path, output)
        trace = json.loads((output/'audit_trace.json').read_text())
        (output/'export-receipt.json').rename(output/'pipeline-export-receipt.json')
    else:
        trace = {'schema': 'clockbase-audit-trace-v1', 'study_id': version['study_id'],
            'run_id': 'sha256:'+_digest(snapshot), 'recorded_at_utc': datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
            'scope': 'Five completed producer stages, saved clock predictions and a coding contract failure; no producer report or comprehensive audit',
            'evaluation_split': 'development', 'models': [
                {'provider': 'openai', 'model': config['model'], 'role': 'archived original CodingAgent'},
                {'provider': 'openai', 'model': config.get('audit_model') or config['model'], 'role': 'source-design reviewer'}],
            'overall': {'status': 'CANNOT VERIFY', 'reason': 'Coding stopped before a report; scientific acceptance remains unresolved'},
            'domains': [{'id': name, 'title': name.replace('_', ' ').capitalize(),
                'model_assessment': {'status': 'CANNOT VERIFY', 'summary': 'The partial run does not establish this scientific domain',
                    'evidence': [], 'unresolved': ['No completed comprehensive audit'], 'proposed_repair': None,
                    'verification_needed': ['Independent scientific source and method review']},
                'verification': {'status': 'CANNOT VERIFY', 'scope': 'Completed stage bytes and saved failure only', 'checks': []}}
                for name in AUDIT_DOMAINS], 'comparisons': [], 'covariates': [], 'artifacts': [], 'events': [],
            'limitations': ['A coding response contract failed after clock prediction. Later report and audit outputs are absent.',
                'Saved predictions and source-stage completion do not establish biological units, clock applicability or an error rate.']}
        output.mkdir(parents=True); (output/'public').mkdir()
        for label, value in {'source-records': source, 'stage-manifest': {'files': snapshot, 'source_config_sha256': sha(pipeline/'config.json')}}.items():
            _safe(value); path = output/'public'/(label+'.json'); write_new(path, value)
            trace['artifacts'].append({'id': label, 'kind': 'partial update evidence', 'label': label.replace('-', ' ').capitalize(),
                'sha256': sha(path), 'download_path': 'public/'+path.name})
        target = output/'public/sample-predictions.csv'; shutil.copyfile(predictions, target)
        trace['artifacts'].append({'id': 'sample-predictions', 'kind': 'saved data', 'label': 'Saved sample predictions',
            'sha256': sha(target), 'download_path': 'public/'+target.name})
        write_new(plan_path, {'schema': 'clockbase-partial-update-export-plan-v1', 'binding_sha256': sha(readback_path),
            'source_snapshot_sha256': _digest(snapshot), 'export_identity': _export_identity()})
    trace['partial_update'] = partial
    trace['scope'] = 'Partial future-data update: '+partial['outcome'].replace('_', ' ')+'; scientific status CANNOT VERIFY'
    trace['resource_usage'] = {'generation_attempts': len(partial['api_calls']), 'responses_received': len(partial['api_calls']),
        'estimated_cost_usd': partial['estimated_cost_usd'], 'billed_cost_usd': None}
    trace['events'].append({'stage': partial['failure']['phase']+' '+partial['failure']['stage'],
        'status': 'FAIL', 'summary': partial['failure']['error']})
    _safe(partial); write_new(output/'public/partial-update.json', partial)
    trace['artifacts'].append({'id': 'partial-update', 'kind': 'partial update evidence', 'label': 'Partial run, failure and API receipts',
        'sha256': sha(output/'public/partial-update.json'), 'download_path': 'public/partial-update.json'})
    _safe(trace)
    (output/'audit_trace.json').write_text(json.dumps(trace, indent=2, allow_nan=False)+'\n')
    if _partial_future_snapshot(version_dir, binding)[0] != snapshot or any(sha(p) != digest for p, digest in observed.items()):
        raise PipelineIntegrityError('Partial source changed during export')
    write_new(output/'export-receipt.json', {'schema': 'clockbase-partial-update-export-v1',
        'binding_sha256': sha(readback_path), 'source_snapshot_sha256': _digest(snapshot),
        'new_model_calls': 0, 'clocks_rerun': False,
        'output_sha256': {str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()}})
    return {'status': 'PASS', 'study_id': version['study_id'], 'outcome': partial['outcome'],
        'trace_sha256': sha(output/'audit_trace.json'), 'scientific_status': 'CANNOT VERIFY', 'new_model_calls': 0}


def export_future_followup(followup_dir, parent_trace_dir, output, *, independent_review_path=None):
    """Export a completed, sealed audit-only follow-up in a new public version."""
    followup_dir, parent_trace_dir, output = (Path(p).resolve() for p in (followup_dir, parent_trace_dir, output))
    if output.exists() or any(output.is_relative_to(p) or p.is_relative_to(output) for p in (followup_dir, parent_trace_dir)):
        raise ValueError('Follow-up export requires a new directory outside its sources')
    observed = {}; helpers = _corpus_projection_helpers(observed)
    checked, read, clean, project, stage, artifact_from_stage = helpers
    def sealed(path):
        path = Path(path)
        return read(checked(path, checked(path.with_suffix(path.suffix+'.sha256')).read_text().strip()))
    frozen = sealed(followup_dir/'followup.json'); result = sealed(followup_dir/'result.json')
    if frozen.get('schema') != 'clockbase-pipeline-audit-followup-v1' or frozen.get('new_producer_or_clock_execution') is not False:
        raise ValueError('Unsupported audit-only follow-up manifest')
    request = frozen['request']; pipeline = Path(request['pipeline_root']).resolve(); version_dir = pipeline.parent
    if output.is_relative_to(version_dir) or version_dir.is_relative_to(output):
        raise ValueError('Follow-up export may not modify its producer version')
    if pipeline.name != 'pipeline' or Path(request['output']).resolve() != followup_dir/'audit':
        raise ValueError('Follow-up request is outside its declared producer or audit root')
    version = sealed(version_dir/'version.json')
    if (version['version_id'] != _digest(version['identity']) or version['version_id'] != version_dir.name
            or request['version_id'] != version['version_id'] or request['study_id'] != version['study_id']
            or request['modality'] != version['modality'] or request['source_bindings'] != version['sources']):
        raise ValueError('Follow-up does not bind the immutable producer version')
    for record in frozen['input_files'].values(): checked(record['path'], record['sha256'])
    for name, expected in frozen['implementation_sha256'].items():
        checked(_inside(followup_dir/'implementation', followup_dir/'implementation'/name), expected)
    for source in version['sources'].values(): checked(source['path'], source['sha256'])
    for name, expected in version['implementation'].items():
        checked(_inside(version_dir/'implementation', version_dir/'implementation'/name), expected)
    replay = frozen['source_design_replay']; prior_study = Path(replay['root']).resolve()
    if prior_study != version_dir/'comprehensive-audit/study':
        raise ValueError('Follow-up supersedes an audit outside the original producer version')
    for record in replay['files'].values():
        checked(_inside(prior_study, record['path']), record['sha256'])
    prior_terminal_path = checked(prior_study/'terminal.json', checked(prior_study/'terminal.sha256').read_text().strip())
    if not any(Path(item['path']).resolve() == prior_terminal_path and item['sha256'] == sha(prior_terminal_path)
            for item in frozen['input_files'].values()):
        raise ValueError('Superseded audit terminal is not frozen by the follow-up')
    prior_terminal = read(prior_terminal_path)
    for relative, expected in prior_terminal['output_sha256'].items(): checked(_inside(prior_study, prior_study/relative), expected)
    parent_receipt = read(parent_trace_dir/'export-receipt.json')
    for relative, expected in parent_receipt['output_sha256'].items():
        checked(_inside(parent_trace_dir, parent_trace_dir/relative), expected)
    trace = read(parent_trace_dir/'audit_trace.json'); parent_hash = sha(parent_trace_dir/'audit_trace.json')
    if (trace['study_id'] != version['study_id'] or trace['future_ingestion']['version_id'] != version['version_id']
            or trace['future_ingestion']['source_version_sha256'] != sha(version_dir/'version.json')
            or trace['comprehensive_audit']['terminal_sha256'] != sha(prior_terminal_path)):
        raise ValueError('Parent public trace does not represent the superseded producer audit')
    study = followup_dir/'audit/study'; terminal = read(study/'terminal.json')
    _future_terminal_result(terminal, result)
    cfg = read(study/'config.json'); task = cfg['task']
    if (result['pipeline_root_sha256'] != sha(checked(pipeline/'config.json'))
            or task.get('producer_version_id') != version['version_id'] or task.get('source_design_replay') != replay
            or task.get('audit_followup') != {'revision': frozen['revision'], 'manifest': {'path': str(followup_dir/'followup.json'), 'sha256': sha(followup_dir/'followup.json')}}):
        raise ValueError('Follow-up study does not bind its sealed replay manifest')
    for name in ('01_source_index', '02_blinded_design'):
        before, after = stage(prior_study, name), stage(study, name)
        if before is None or after is None or before['result'] != after['result']:
            raise ValueError('Audit-only follow-up changed the blinded source index or design')
        if name == '02_blinded_design' and (before['result'].get('frozen_before_saved_outcomes') is not True
                or before['result'].get('source_index_sha256') != _corpus_record_hash(stage(prior_study, '01_source_index')['result'])):
            raise ValueError('Replayed design is not bound to the unchanged source-first index')
    replay_provenance = artifact_from_stage(study, '02_blinded_design', 'replay-provenance.json')
    if replay_provenance != {'binding': replay, 'new_model_calls': 0, 'new_independent_design_assessment': False}:
        raise ValueError('Follow-up does not preserve its source-design replay provenance')
    packet, summary = _project_corpus_study(study, terminal, task, pipeline, helpers)
    if not packet.get('measurement_evidence') or not packet['measurement_evidence']['review_request_bindings']:
        raise ValueError('Completed follow-up lacks measurement evidence supplied to its reviewer')
    clock = stage(pipeline, '05_clock_predictions'); stage(pipeline, '06_coding_agent')
    predictions = _future_prediction_source(pipeline, clock, checked)
    _future_producer_binding(packet, summary, result, version['study_id'], sha(predictions), sha(pipeline/'06_coding_agent/complete.json'))
    original_packet = read(parent_trace_dir/next(a['download_path'] for a in trace['artifacts'] if a['id'] == 'comprehensive-audit'))
    if packet['original_rows'] != original_packet['original_rows']:
        raise ValueError('Audit-only follow-up changed the producer rows or sample lineage')
    audit_calls = []
    for relative in sorted(terminal['output_sha256']):
        if not re.fullmatch(r'api/[^/]+/receipt\.json', relative): continue
        receipt = read(study/relative)
        usage = receipt.get('usage', {})
        if any(type(usage.get(key)) is not int or usage[key] < 0 for key in ('prompt_tokens', 'completion_tokens')):
            raise ValueError('Follow-up usage must retain explicit nonnegative token counts')
        expected_cost = (usage['prompt_tokens']*request['budget']['input_usd_per_million']
            + usage['completion_tokens']*request['budget']['output_usd_per_million'])/1e6
        if abs(receipt['estimated_cost_usd']-expected_cost) > 1e-12:
            raise ValueError('Follow-up usage estimate differs from its frozen tariffs')
        audit_calls.append(project(receipt, ('requested_model', 'returned_model', 'usage', 'estimated_cost_usd', 'billed_cost_usd', 'elapsed_seconds'))
            | {'receipt_sha256': sha(study/relative)})
    budget = result['budget']; estimated = sum(row['estimated_cost_usd'] for row in audit_calls)
    request_count = sum(bool(re.fullmatch(r'api/[^/]+/request\.json', relative)) for relative in terminal['output_sha256'])
    if (len(audit_calls) != request_count or len(audit_calls) != budget['call_reservations'] or abs(estimated-budget['estimated_or_reserved_usd']) > 1e-10
            or any(row['requested_model'] != request['budget']['model'] or row.get('billed_cost_usd') is not None for row in audit_calls)
            or estimated > request['budget']['max_cost_usd'] or len(audit_calls) > request['budget']['max_calls']):
        raise ValueError('Completed follow-up API receipts do not reconcile with its separate budget')
    independent = None
    if independent_review_path is not None:
        independent_review_path = checked(independent_review_path); review = read(independent_review_path)
        expected = {'followup_manifest_sha256': sha(followup_dir/'followup.json'), 'result_sha256': sha(followup_dir/'result.json'),
                    'terminal_sha256': sha(study/'terminal.json')}
        if (any(review.get(key) != value for key, value in expected.items()) or review.get('status') != 'PASS'
                or review.get('scientific_status') != 'CANNOT VERIFY'):
            raise ValueError('Independent follow-up readback does not bind this execution')
        independent = project(review, ('status', 'scope', 'scientific_status', 'checks', 'limitations', 'recorded_at_utc', *expected))
        independent['source_receipt_sha256'] = sha(independent_review_path)
    followup = {'schema': 'clockbase-public-audit-followup-v1', 'revision': frozen['revision'],
        'study_id': version['study_id'], 'producer_version_id': version['version_id'],
        'followup_manifest_sha256': sha(followup_dir/'followup.json'), 'result_sha256': sha(followup_dir/'result.json'),
        'terminal_sha256': sha(study/'terminal.json'), 'created_at_utc': frozen['created_at_utc'],
        'supersedes': {'trace_sha256': parent_hash, 'terminal_sha256': sha(prior_terminal_path), 'run_key': parent_trace_dir.parent.name},
        'source_index_reused': True, 'grouping_design_reused': True, 'new_independent_design_assessment': False,
        'source_design_sha256': _corpus_record_hash(stage(study, '02_blinded_design')['result']), 'clocks_rerun': False, 'producer_rerun': False,
        'measurement_evidence_records': len(packet['measurement_evidence']['evidence']),
        'measurement_registry_sha256': packet['measurement_evidence']['registry_sha256'],
        'additional_audit_budget': project(budget, ('max_cost_usd', 'max_calls', 'call_reservations', 'estimated_or_reserved_usd', 'billed_cost_usd', 'by_status')),
        'api_receipts': audit_calls, 'independent_review': independent, 'scientific_status': 'CANNOT VERIFY', 'accepted_scientific_repairs': 0,
        'difference': 'Measurement records added; grouping design reused; clocks not rerun; additional audit cost recorded separately.'}
    trace['audit_followup'] = followup
    trace.pop('future_independent_review', None)  # Its unchanged file remains a prior-review download.
    trace['comprehensive_audit'] = summary | {'artifact_id': 'comprehensive-audit'}
    trace['future_ingestion']['comprehensive_summary'] = summary
    trace['scope'] = 'Executed source pipeline with a separately budgeted audit-only measurement-evidence follow-up; scientific acceptance and paper-extension eligibility remain unresolved.'
    finished = datetime.fromisoformat(result['finished_at_utc'])
    if finished.tzinfo is None or finished.utcoffset() != timezone.utc.utcoffset(finished):
        raise ValueError('Follow-up finish time must identify UTC explicitly')
    trace['recorded_at_utc'] = finished.isoformat().replace('+00:00', 'Z')
    for artifact in trace['artifacts']:
        if artifact['id'] in {'comprehensive-audit', 'future-ingestion', 'independent-future-readback'}:
            artifact['id'] = 'prior-'+artifact['id']; artifact['label'] = 'Prior '+artifact['label']
    additions = {'comprehensive-audit': ('comprehensive-audit-followup.json', packet),
        'future-ingestion': ('future-ingestion-followup.json', trace['future_ingestion']),
        'audit-followup': ('audit-followup.json', followup), 'measurement-evidence': ('measurement-evidence.json', packet['measurement_evidence']),
        'superseded-audit-trace': ('superseded-audit-trace.json', read(parent_trace_dir/'audit_trace.json'))}
    for _, value in additions.values(): _safe(value)
    shutil.copytree(parent_trace_dir, output)
    receipt_target = output/'receipts'/('parent-'+sha(output/'export-receipt.json')+'.json')
    receipt_target.parent.mkdir(exist_ok=True); (output/'export-receipt.json').rename(receipt_target)
    for artifact_id, (name, value) in additions.items():
        path = output/'public'/name; write_new(path, value)
        trace['artifacts'].append({'id': artifact_id, 'kind': 'audit-only follow-up evidence', 'label': artifact_id.replace('-', ' ').capitalize(),
            'sha256': sha(path), 'download_path': 'public/'+name})
    _safe(trace); (output/'audit_trace.json').write_text(json.dumps(trace, indent=2, allow_nan=False)+'\n')
    for path, expected in observed.items():
        if sha(path) != expected: raise PipelineIntegrityError('Follow-up source changed during export')
    write_new(output/'export-receipt.json', {'schema': 'clockbase-future-followup-export-v1', 'source_trace_sha256': parent_hash,
        'followup_manifest_sha256': followup['followup_manifest_sha256'], 'result_sha256': followup['result_sha256'],
        'export_identity': _export_identity(), 'new_model_calls': 0, 'clocks_rerun': False, 'scientific_status': 'CANNOT VERIFY',
        'output_sha256': {str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()}})
    return {'status': 'PASS', 'trace_sha256': sha(output/'audit_trace.json'), 'measurement_evidence_records': followup['measurement_evidence_records'],
        'additional_audit_calls': len(audit_calls), 'additional_audit_estimated_usd': estimated, 'scientific_status': 'CANNOT VERIFY', 'new_model_calls': 0}


def export_audit_reassessment(followup_dir, output, *, plan_path, independent_review_path=None):
    """Expose a new source-first audit after a failed hook, preserving both runs.

    The producer has completed all eight stages. This narrower contract accepts
    the sealed handoff-repair driver without pretending it replayed an old design
    or that the failed original comprehensive audit completed.
    """
    followup_dir, output, plan_path = (Path(p).resolve() for p in (followup_dir, output, plan_path))
    observed = {}; helpers = _corpus_projection_helpers(observed)
    checked, read, clean, project, stage, artifact_from_stage = helpers
    def sealed(path):
        path = Path(path)
        return read(checked(path, checked(path.with_suffix(path.suffix+'.sha256')).read_text().strip()))
    fresh = (followup_dir/'followup.json').is_file()
    manifest_path = followup_dir/('followup.json' if fresh else 'freeze.json')
    if fresh and (followup_dir/'freeze.json').exists():
        raise ValueError('Ambiguous reassessment manifests')
    raw_frozen = sealed(manifest_path); result = sealed(followup_dir/'result.json')
    if fresh:
        if (raw_frozen.get('schema') != 'clockbase-pipeline-audit-followup-v2'
                or raw_frozen.get('mode') != 'fresh_source_v1'
                or raw_frozen.get('new_producer_or_clock_execution') is not False
                or 'source_design_replay' in raw_frozen or independent_review_path is None):
            raise ValueError('Fresh reassessment requires its sealed source-first manifest and independent readback')
        frozen = {'request': raw_frozen['request'],
            'input_sha256': {record['path']:record['sha256'] for record in raw_frozen['input_files'].values()},
            'original_audit_output':raw_frozen['fresh_source']['prior_audit']}
        for name, expected in raw_frozen['implementation_sha256'].items():
            checked(_inside(followup_dir/'implementation', followup_dir/'implementation'/name), expected)
    else:
        frozen = raw_frozen
        execution = sealed(followup_dir/'execution-receipt.json')
        if (frozen.get('schema') != 'clockbase-audit-handoff-followup-v1'
                or execution.get('schema') != 'clockbase-audit-handoff-execution-v1'
                or execution.get('freeze_sha256') != sha(manifest_path)
                or execution.get('result_sha256') != sha(followup_dir/'result.json')
                or execution.get('producer_executions') != 0 or execution.get('clock_executions') != 0
                or execution.get('scientific_status') != 'CANNOT VERIFY'):
            raise ValueError('Reassessment execution does not bind its frozen audit-only inputs')
    request = frozen['request']; pipeline = Path(request['pipeline_root']).resolve(); version_dir = pipeline.parent
    if (output.exists() or any(output.is_relative_to(p) or p.is_relative_to(output) for p in (followup_dir, version_dir))
            or pipeline.name != 'pipeline' or Path(request['output']).resolve() != followup_dir/'audit'
            or request.get('source_design_replay') is not None):
        raise ValueError('Reassessment requires new output and a new source-design assessment')
    for path, expected in frozen['input_sha256'].items(): checked(path, expected)
    if not fresh:
        for name, expected in frozen['implementation_sha256'].items():
            checked(_inside(frozen['runtime'], Path(frozen['runtime'])/name), expected)
        checked(followup_dir/'run_frozen_audit_2026_09_21.py', frozen['driver_sha256'])
    version = sealed(version_dir/'version.json')
    if (version['version_id'] != _digest(version['identity']) or version['version_id'] != version_dir.name
            or request['version_id'] != version['version_id'] or request['study_id'] != version['study_id']
            or request['modality'] != version['modality'] or request['source_bindings'] != version['sources']):
        raise ValueError('Reassessment changed its immutable producer identity')
    for source in version['sources'].values(): checked(source['path'], source['sha256'])
    for name, expected in version['implementation'].items():
        checked(_inside(version_dir/'implementation', version_dir/'implementation'/name), expected)
    prior_study = Path(frozen['original_audit_output']).resolve()/'study'
    if prior_study != version_dir/'comprehensive-audit/study':
        raise ValueError('Prior failed audit belongs to a different producer')
    prior_terminal_path = prior_study/'terminal.json'
    checked(prior_terminal_path, checked(prior_study/'terminal.sha256').read_text().strip())
    if frozen['input_sha256'].get(str(prior_terminal_path)) != sha(prior_terminal_path):
        raise ValueError('Prior terminal is absent from the frozen inputs')
    prior_terminal = read(prior_terminal_path)
    prior_result_path = version_dir/'comprehensive-audit-result.json'
    prior_result = read(checked(prior_result_path, frozen['input_sha256'][str(prior_result_path)]))
    if (prior_terminal.get('status') != 'UNRESOLVED' or prior_result.get('execution_status') == 'COMPLETED'
            or prior_result.get('study_id') != version['study_id']
            or any(prior_result.get(k) != v for k, v in prior_terminal.items())):
        raise ValueError('Reassessment requires the preserved failed original audit')
    for relative, expected in prior_terminal['output_sha256'].items():
        checked(_inside(prior_study, prior_study/relative), expected)
    study = followup_dir/'audit/study'; terminal = read(study/'terminal.json')
    _future_terminal_result(terminal, result)
    cfg = read(study/'config.json'); task = cfg['task']
    if (result['pipeline_root_sha256'] != sha(checked(pipeline/'config.json'))
            or task.get('producer_version_id') != version['version_id'] or task.get('source_design_replay') is not None):
        raise ValueError('Reassessment is not bound to the unchanged producer')
    design = stage(study, '02_blinded_design'); index = stage(study, '01_source_index')
    if (not design or not index or design['result'].get('frozen_before_saved_outcomes') is not True
            or design['result'].get('source_index_sha256') != _corpus_record_hash(index['result'])
            or any(Path(p).name == 'replay-provenance.json' for p in design['outputs'])):
        raise ValueError('Reassessment lacks a new source-first design')
    independent = None
    if fresh:
        expected_task = {'revision':raw_frozen['revision'], 'mode':'fresh_source_v1',
            'manifest':{'path':str(manifest_path),'sha256':sha(manifest_path)}}
        if task.get('audit_followup') != expected_task:
            raise ValueError('Fresh reassessment task does not bind its declared follow-up')
        prior_index = stage(prior_study, '01_source_index')
        if (prior_index is None or prior_index['result'] != index['result']
                or _corpus_record_hash(index['result']) != raw_frozen['fresh_source']['source_index_sha256']):
            raise ValueError('Fresh reassessment changed its frozen source index')
        first = read(study/'api/source_design_00-000/request.json')
        messages = first.get('messages')
        if (not isinstance(messages,list) or len(messages) != 2
                or messages[0].get('role') != 'system' or messages[1].get('role') != 'user'
                or first.get('model') != request['budget']['model']
                or _corpus_record_hash({'system':messages[0].get('content'),
                    'payload':json.loads(messages[1].get('content',''))}) != raw_frozen['fresh_source']['initial_request_sha256']):
            raise ValueError('Fresh source request differs from its prelaunch blinded packet')
        independent_review_path = checked(independent_review_path)
        review = read(independent_review_path)
        expected_review = {'followup_manifest_sha256':sha(manifest_path),
            'result_sha256':sha(followup_dir/'result.json'),'terminal_sha256':sha(study/'terminal.json')}
        if (review.get('status') != 'PASS' or review.get('scientific_status') != 'CANNOT VERIFY'
                or any(review.get(k) != v for k,v in expected_review.items())
                or not isinstance(review.get('source_claim_findings'),list)):
            raise ValueError('Independent source readback does not bind this reassessment')
        independent = project(review, ('status','scope','scientific_status','checks','limitations',
            'recorded_at_utc','source_claim_findings',*expected_review))
        independent['source_receipt_sha256'] = sha(independent_review_path)
        _safe(independent)
    packet, summary = _project_corpus_study(study, terminal, task, pipeline, helpers)
    measurement = packet.get('measurement_evidence')
    if not measurement or not measurement['review_request_bindings']:
        raise ValueError('Reassessment did not supply measurement records to its reviewer')
    clock = stage(pipeline, '05_clock_predictions'); stage(pipeline, '06_coding_agent')
    predictions = _future_prediction_source(pipeline, clock, checked)
    _future_producer_binding(packet, summary, result, version['study_id'], sha(predictions), sha(pipeline/'06_coding_agent/complete.json'))
    calls = []
    for relative in sorted(terminal['output_sha256']):
        if not re.fullmatch(r'api/[^/]+/receipt\.json', relative): continue
        receipt = read(study/relative); usage = receipt['usage']
        if any(type(usage.get(k)) is not int or usage[k] < 0 for k in ('prompt_tokens', 'completion_tokens')):
            raise ValueError('Reassessment usage counts are invalid')
        estimated = sum(usage[k]*request['budget'][r] for k,r in (
            ('prompt_tokens','input_usd_per_million'),('completion_tokens','output_usd_per_million')))/1e6
        if (abs(receipt['estimated_cost_usd']-estimated) > 1e-12
                or receipt['requested_model'] != request['budget']['model'] or receipt.get('billed_cost_usd') is not None):
            raise ValueError('Reassessment cost differs from frozen tariffs or model')
        calls.append(project(receipt, ('requested_model','returned_model','usage','estimated_cost_usd','billed_cost_usd','elapsed_seconds'))
                     | {'receipt_sha256': sha(study/relative), 'role': Path(relative).parent.name})
    budget = result['budget']
    if (len(calls) != budget['call_reservations'] or len(calls) > request['budget']['max_calls']
            or sum(bool(re.fullmatch(r'api/[^/]+/request\.json', p)) for p in terminal['output_sha256']) != len(calls)
            or not any(c['role'].startswith('source_design_') for c in calls)
            or abs(sum(c['estimated_cost_usd'] for c in calls)-budget['estimated_or_reserved_usd']) > 1e-10
            or budget['estimated_or_reserved_usd'] > request['budget']['max_cost_usd']):
        raise ValueError('Reassessment calls do not reconcile with its separate budget')
    date, quote = version['source_design']['series_first_public_date'], version['source_design']['source_status_quote']
    _future_public_date(date, quote, _geo_source_text(checked(version['sources']['series']['path'])))
    reassessment = {'schema':'clockbase-public-audit-reassessment-v1','study_id':version['study_id'],
        'producer_version_id':version['version_id'],'version_sha256':sha(version_dir/'version.json'),
        'freeze_sha256':sha(manifest_path),'result_sha256':sha(followup_dir/'result.json'),
        'terminal_sha256':sha(study/'terminal.json'),
        'source_design_sha256':_corpus_record_hash(design['result']),
        'measurement_registry_sha256':measurement['registry_sha256'],'measurement_evidence_records':len(measurement['evidence']),
        'new_independent_design_assessment':True,'grouping_design_reused':False,'producer_rerun':False,'clocks_rerun':False,
        'source_first_public_date':date,'source_first_public_quote':quote,
        'prior_audit':project(prior_result, ('status','execution_status','model_error','budget','reviewed_rows','input_rows'))
            | {'terminal_sha256':sha(prior_terminal_path),'result_sha256':sha(version_dir/'comprehensive-audit-result.json')},
        'additional_audit_budget':project(budget, ('max_cost_usd','max_calls','call_reservations','estimated_or_reserved_usd','billed_cost_usd','by_status')),
        'api_receipts':calls,'scientific_status':'CANNOT VERIFY','accepted_scientific_repairs':0,
        'original_paper_2025_eligibility':'QUARANTINED',
        'scope':'New source-only design and measurement review after a preserved failed audit hook. Producer outputs are unchanged; descriptive calculations do not establish a biological effect.'}
    if fresh:
        reassessment.update(schema='clockbase-public-audit-reassessment-v2',
            followup_manifest_sha256=sha(manifest_path), independent_review=independent,
            independent_readback_sha256=sha(independent_review_path), source_index_reused=True,
            scope='New source-first audit after a preserved source-design contract failure. Independent source review findings remain visible, including errors in the audit. Completed execution and conditional calculations do not establish scientific acceptance.')
    else:
        reassessment['execution_receipt_sha256'] = sha(followup_dir/'execution-receipt.json')
    _safe(reassessment)
    prepare(pipeline, plan_path); export(plan_path, output)
    trace_path = output/'audit_trace.json'; trace = json.loads(trace_path.read_text())
    shutil.copyfile(trace_path, output/'public/pipeline-audit-trace.json')
    (output/'export-receipt.json').rename(output/'pipeline-export-receipt.json')
    trace['audit_reassessment'] = reassessment
    trace['comprehensive_audit'] = summary | {'artifact_id':'comprehensive-audit'}
    trace['scope'] = reassessment['scope']
    finished = datetime.fromisoformat(result['finished_at_utc'])
    if finished.tzinfo is None or finished.utcoffset() != timezone.utc.utcoffset(finished):
        raise ValueError('Reassessment finish time must identify UTC explicitly')
    trace['recorded_at_utc'] = finished.isoformat().replace('+00:00','Z')
    trace['models'].append({'provider':'openai','model':request['budget']['model'],'role':'separate source-design and comprehensive reassessment'})
    additions = {'comprehensive-audit':packet,'audit-reassessment':reassessment,'measurement-evidence':measurement,
                 'prior-failed-audit':reassessment['prior_audit']}
    if fresh:
        additions['independent-source-readback'] = independent
    for artifact_id, value in additions.items():
        _safe(value); path = output/'public'/(artifact_id+'.json'); write_new(path,value)
        trace['artifacts'].append({'id':artifact_id,'kind':'separate audit reassessment','label':artifact_id.replace('-',' ').capitalize(),
                                  'sha256':sha(path),'download_path':'public/'+path.name})
    trace['artifacts'].append({'id':'pipeline-audit-trace','kind':'preserved producer trace','label':'Original producer trace',
                              'sha256':sha(output/'public/pipeline-audit-trace.json'),'download_path':'public/pipeline-audit-trace.json'})
    _safe(trace); trace_path.write_text(json.dumps(trace,indent=2,allow_nan=False)+'\n')
    for path, expected in observed.items():
        if sha(path) != expected: raise PipelineIntegrityError('Reassessment input changed during export')
    write_new(output/'export-receipt.json',{'schema':'clockbase-audit-reassessment-export-v1','freeze_sha256':sha(manifest_path),
        'plan_sha256':sha(plan_path),'export_identity':_export_identity(),'new_model_calls':0,'clocks_rerun':False,
        'scientific_status':'CANNOT VERIFY','output_sha256':{str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file()}})
    return {'status':'PASS','study_id':version['study_id'],'trace_sha256':sha(trace_path),
            'additional_audit_calls':len(calls),'new_model_calls':0,'scientific_status':'CANNOT VERIFY'}


_SOURCE_REFERENCE_VIEW = (
    ('data_provenance', 'Assay context and study-family links', ('assay_context', 'family_overlap')),
    ('sample_identity', 'Sample record identity', ('sample_identity',)),
    ('biological_units', 'Biological units', ('biological_units',)),
    ('grouping', 'Comparison definitions and pairing', ('comparison_definition', 'pairing')),
    ('covariates', 'Covariates and confounding', ('confounding',)),
    ('measurement_applicability', 'Measurement prerequisites', ('measurement_prerequisites',)),
    ('statistical_method', 'Conditional statistical requirements', ('test_requirements',)),
)


def _source_reference_snapshot(release_path, release_sha256, registry_path, registry_sha256):
    """Read a complete pinned version chain and only its allowed source payloads."""
    from . import evaluation_2026_09_10 as readiness
    for checksum in (release_sha256, registry_sha256):
        if re.fullmatch(r'[0-9a-f]{64}', checksum or '') is None:
            raise ValueError('Source-reference release and registry require explicit SHA256 bindings')
    observed = {}
    checked, _, clean, _, _, _ = _corpus_projection_helpers(observed)
    checked(release_path, release_sha256)
    context = readiness.source_readiness_context(release_path, release_sha256)
    for binding in context['release']['bindings'].values():
        checked(binding['path'], binding['sha256'])
    bindings = readiness._readiness_json(checked(registry_path, registry_sha256))
    if not isinstance(bindings, list):
        raise ValueError('Source-reference registry must be a complete ordered list')
    versions, payloads = [], {}
    for binding in bindings:
        frozen = readiness._readiness_frozen(context, binding)
        versions.append((binding, frozen))
        path = checked(binding['path'], binding['sha256'])
        stage = path.parent.parent
        receipt = readiness._readiness_json(checked(
            stage/'complete.json', checked(stage/'complete.sha256').read_text().strip()))
        checked(stage.parent/'config.json')
        for relative, checksum in receipt['outputs'].items():
            checked(_inside(stage, stage/relative), checksum)
        packet = readiness._readiness_json(path.parent/'packet.json', frozen['artifacts']['packet.json'])
        try:
            submission = readiness._readiness_json(
                path.parent/'submission.json', frozen['artifacts']['submission.json'])
        except (ValueError, UnicodeError) as error:
            if frozen['validation']['validation_status'] != 'REJECTED':
                raise ValueError('Accepted source-reference submission is unreadable') from error
            submission = None
        sid = frozen['study_id']
        payloads.setdefault(sid, []).append({
            'binding': binding, 'frozen': frozen, 'packet': packet, 'submission': submission})
    table = readiness.source_readiness_table(context['enrollment'], versions)
    for sid in payloads:
        # context and frozen checks reject administrative holds before raw access.
        for binding in context['allowed_by_id'][sid]['sources'].values():
            source = checked(binding['path'], binding['sha256'])
            if source.stat().st_size != binding['bytes']:
                raise PipelineIntegrityError('Source-reference raw-source size changed')
    return context, table, payloads, {str(p): h for p, h in sorted(observed.items())}, clean


def prepare_source_reference(release_path, release_sha256, registry_path, registry_sha256, plan_path):
    """Pin a complete source-reference snapshot without computing scientific outcomes."""
    context, table, payloads, snapshot, _ = _source_reference_snapshot(
        release_path, release_sha256, registry_path, registry_sha256)
    plan = {
        'schema': 'clockbase-source-reference-trace-plan-v1',
        'release': {'path': str(Path(release_path).resolve()), 'sha256': release_sha256},
        'registry': {'path': str(Path(registry_path).resolve()), 'sha256': registry_sha256},
        'recorded_at_utc': datetime.now(timezone.utc).isoformat().replace('+00:00', 'Z'),
        'snapshot_sha256': snapshot, 'export_identity': _export_identity(),
        'enrolled_studies': table['enrolled_studies'], 'submitted_studies': len(payloads),
        'source_validated_studies': sum(r['source_validated'] for r in table['studies']),
        'reference_runtime': context['release']['reference_runtime'],
        'max_download_bytes': 50_000_000,
        'scope': 'Source-attribution development readiness; no scientific acceptance or accuracy estimate',
    }
    write_new(Path(plan_path), plan)
    return {'status': 'PASS', 'plan_sha256': sha(plan_path), 'files_pinned': len(snapshot),
            'submitted_studies': len(payloads), 'new_model_calls': 0}


def _source_reference_evidence(judgment, evidence_ids):
    """Show five intact public witnesses; the full filtered judgment is downloadable."""
    ids = list(dict.fromkeys(evidence_ids))
    evidence, withheld = [], 0
    for key in ids:
        witness = judgment['evidence'][key]
        try:
            _safe(witness['quote'])
        except ValueError:
            withheld += 1
            continue
        if not witness['quote'].strip() or len(witness['quote']) > 12000:
            withheld += 1
            continue
        if len(evidence) < 5:
            entity = witness['entity']
            evidence.append({
                'sample_id' if entity.startswith('GSM') else 'source_id': entity,
                'field': witness['field'] + ' (source lines ' + str(witness['line_start'])
                         + '–' + str(witness['line_end']) + ')',
                'quote': witness['quote'],
                'source_url': _public_url('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc='+entity),
            })
    note = (f'Displaying {len(evidence)} of {len(ids)} linked witnesses; {withheld} omitted from the inline '
            'view for privacy, empty content or length. Complete public-safe witnesses and source hashes '
            'are in the version-history download.')
    return evidence, note


def _source_reference_inline(value):
    """Keep the inline string contract without abbreviating scientific source fields."""
    if isinstance(value, dict):
        return {key: _source_reference_inline(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_source_reference_inline(item) for item in value]
    if isinstance(value, str) and (not value.strip() or len(value) > 12000):
        return ('[Inline text omitted: ' + str(len(value)) + ' characters; SHA256 '
                + hashlib.sha256(value.encode()).hexdigest()
                + '. Full public-safe source fields are retained in the version-history download.]')
    return value


def _source_reference_trace(plan, plan_sha256, row, history, artifacts, clean):
    current = history[-1]
    frozen, judgment = current['frozen'], current['submission']
    valid = row['source_validated']
    source_domains = frozen['validation']['domains'] if valid else {}
    unresolved = {r['unresolved_id']: r for r in judgment['unresolved']} if valid else {}
    domains = []
    for public_id, title, names in _SOURCE_REFERENCE_VIEW:
        assertions = [source_domains[name] for name in names] if valid else []
        statuses = {r['status'] for r in assertions}
        status = ('FAIL' if 'FAIL' in statuses else
                  'CANNOT VERIFY' if not statuses or 'CANNOT VERIFY' in statuses else 'PASS')
        ids = [key for assertion in assertions for key in assertion['evidence_ids']]
        evidence, note = _source_reference_evidence(judgment, ids) if valid else ([], '')
        dependencies = list(dict.fromkeys(
            key for assertion in assertions for key in assertion['unresolved_ids']))
        summary = (' '.join(name + ' [' + source_domains[name]['status'] + ']: '
                            + source_domains[name]['rationale'] for name in names) + ' ' + note
                   if valid else 'NOT_ASSESSED: the current source-attribution submission was rejected. '
                   + row['not_assessed_reason'])
        domains.append({
            'id': public_id, 'title': title,
            'model_assessment': {
                'status': status, 'summary': clean(summary), 'evidence': clean(evidence),
                'unresolved': clean([unresolved[k]['question'] + ' ' + unresolved[k]['reason']
                                    for k in dependencies]) if valid else ['Source-attribution validation is incomplete.'],
                'proposed_repair': None,
                'verification_needed': clean([unresolved[k]['needed_evidence'] for k in dependencies])
                    + ['Independent scientific verification of the source interpretation.'],
            },
            'verification': {
                'status': 'PASS' if valid else 'FAIL',
                'scope': 'Archived literal/structural validation and source-file hash readback only; '
                         'this does not certify the scientific judgment.',
                'checks': [
                    {'status': 'PASS', 'description': 'Pinned version chain, stage payloads and allowed source hashes agree.'},
                    {'status': 'PASS' if valid else 'FAIL',
                     'description': 'Current submission: ' + frozen['validation']['validation_status']
                         + '. Quoted-source support does not establish independence, applicability or valid inference.'},
                ],
            },
        })
    comparisons, covariates = [], []
    if valid:
        for comparison in judgment['comparisons']:
            comparisons.append({
                'comparison_id': comparison['comparison_id'], 'comparison_name': comparison['target'],
                'target': comparison['target'], 'status': 'CANNOT VERIFY',
                'case_ids': comparison['case_ids'], 'control_ids': comparison['control_ids'],
                'excluded_ids': comparison['excluded_sample_ids'],
                'biological_unit': comparison['domains']['biological_units']['rationale'],
                'pairing': comparison['domains']['pairing']['rationale'],
                'unresolved': [
                    'Source-declared contrast type: ' + comparison['contrast_type']
                    + '; conditional test family: ' + comparison['conditional_test_family']
                    + '. No statistical outcome is verified here.',
                    *[r['requirement'] for r in comparison['analysis_requirements']],
                ],
            })
        for covariate in judgment['covariates']:
            evidence, note = _source_reference_evidence(judgment, covariate['evidence_ids'])
            covariates.append({
                'comparison_name': 'Source inventory; ' + covariate['scope'],
                'name': covariate['name'], 'decision': 'unknown',
                'source_fields': sorted({judgment['evidence'][k]['field'] for k in covariate['evidence_ids']}),
                'rationale': f"{len(covariate['observed_sample_ids'])} observed and "
                    f"{len(covariate['missing_sample_ids'])} missing sample records; "
                    'a source inventory does not authorize statistical adjustment. ' + note,
                'evidence': evidence,
                'missingness_action': 'Retain missingness and source dependencies; no values imputed or adjustment applied.',
            })
    trace = {
        'schema': 'clockbase-audit-trace-v1', 'study_id': row['study_id'],
        'run_id': 'sha256:' + _digest({'plan_sha256': plan_sha256, 'study_id': row['study_id'],
                                      'current_version_sha256': current['binding']['sha256']}),
        'recorded_at_utc': plan['recorded_at_utc'], 'evaluation_split': 'development', 'models': [],
        'scope': 'Source-only design reconstruction with prior pipeline results withheld. '
                 'Native collaboration review; verified model identity and billable cost are unavailable.',
        'overall': {'status': 'CANNOT VERIFY',
                    'reason': 'Source-attribution readiness is separate from scientific correctness. '
                              'No accepted repair, validated treatment effect or audit error rate is established.'},
        'domains': domains, 'comparisons': clean(comparisons), 'covariates': clean(covariates),
        'artifacts': artifacts,
        'events': [{'stage': 'Source-reference version ' + h['frozen']['version_id'],
                    'status': 'PASS' if h['frozen']['validation']['validation_status'] == 'ACCEPTED_SOURCE_ATTRIBUTION' else 'FAIL',
                    'summary': h['frozen']['validation']['validation_status'] + '. '
                        + ('Current version.' if h is current else 'Preserved superseded version.')}
                   for h in history],
        'resource_usage': {'generation_attempts': None, 'responses_received': None,
                           'estimated_cost_usd': None, 'billed_cost_usd': None},
        'limitations': [
            'These are exposed development studies. Source-reference and method outputs are separate; this is not held-out model accuracy.',
            'Allowed source narratives may contain published findings; withholding prior pipeline results is not total outcome blinding.',
            'The seven display sections project nine named source domains. Combined sections use FAIL, then CANNOT VERIFY, then PASS; the download preserves each original verdict.',
            'All 80 enrolled studies remain in the cohort table. Two administrative holds and unsubmitted or rejected cases are NOT_ASSESSED, not scientific failures or abstentions.',
            'All registered rejected and corrected versions are retained. Prompt-format clarifications and corrections do not establish a model-performance gain.',
            'Native collaboration model/cost are unknown. Direct model API calls by the source controller and exporter are zero; this is not an OpenAI-key trial.',
            'Server paths and potentially private text are omitted or represented by explicit SHA256 withholding descriptors. Other literal quotes are unchanged.',
            'Overlength or empty inline prose is represented by an explicit character count and hash. Full public-safe source fields remain in the version-history download; literal evidence quotes are never shortened.',
            'No outcome values, P values, treatment-effect estimates or clocks are computed by this export. Biological independence and measurement applicability require separate evidence.',
        ],
    }
    trace = _source_reference_inline(trace)
    _safe(trace)
    return trace


def export_source_reference(plan_path, output):
    """Export every submitted study, full source-version history and all enrolled denominators."""
    from . import evaluation_2026_09_10 as readiness
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    plan = readiness._readiness_json(plan_path)
    if (plan.get('schema') != 'clockbase-source-reference-trace-plan-v1'
            or plan.get('export_identity') != _export_identity()
            or plan.get('max_download_bytes') != 50_000_000):
        raise ValueError('Source-reference export plan or implementation changed')
    if output.exists() or any(Path(p).is_relative_to(output) for p in plan['snapshot_sha256']):
        raise ValueError('Source-reference export requires a new directory outside its inputs')
    context, table, payloads, snapshot, clean = _source_reference_snapshot(
        plan['release']['path'], plan['release']['sha256'],
        plan['registry']['path'], plan['registry']['sha256'])
    if (snapshot != plan['snapshot_sha256'] or len(payloads) != plan['submitted_studies']
            or table['enrolled_studies'] != plan['enrolled_studies']
            or sum(r['source_validated'] for r in table['studies']) != plan['source_validated_studies']
            or plan['reference_runtime'] != context['release']['reference_runtime']):
        raise PipelineIntegrityError('Source-reference export inputs differ from the frozen plan')
    public_table = clean(table)
    public_table.update(scientific_error_rate=None, accepted_scientific_repairs=0,
                        registry_sha256=plan['registry']['sha256'], release_sha256=plan['release']['sha256'])
    output.mkdir(parents=True)
    write_new(output/'readiness.json', public_table)
    fields = ['study_id', 'source_access_allowed', 'administrative_state', 'source_record_count',
              'source_validated', 'not_assessed_reason', *table['domains']]
    with (output/'readiness.csv').open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for row in public_table['studies']:
            writer.writerow({**{k: row[k] for k in fields if k not in row['domains']}, **row['domains']})
    with (output/'domain-counts.csv').open('x', newline='') as handle:
        fields = ['domain', 'PASS', 'FAIL', 'CANNOT VERIFY', 'NOT_ASSESSED',
                  'assessed_and_source_validated', 'allowed_studies', 'enrolled_studies', 'administratively_held']
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for domain, counts in table['domains'].items():
            writer.writerow({'domain': domain, **{k: counts[k] for k in fields[1:]}})
    traces = []
    for row in table['studies']:
        sid = row['study_id']
        if sid not in payloads:
            continue
        history = payloads[sid]
        directory = output/sid; directory.mkdir()
        download = directory/'downloads'; download.mkdir()
        public_history = {
            'schema': 'clockbase-public-source-reference-history-v1', 'study_id': sid,
            'scope': plan['scope'], 'reference_runtime': clean(plan['reference_runtime']),
            'release_sha256': plan['release']['sha256'], 'registry_sha256': plan['registry']['sha256'],
            'current_version_sha256': history[-1]['binding']['sha256'],
            'scientific_status': 'CANNOT VERIFY', 'accepted_scientific_repairs': 0, 'scientific_error_rate': None,
            'privacy_policy': 'Paths omitted; unsafe text replaced by explicit withholding hashes; other quotes unchanged.',
            'versions': [{
                'version_id': h['frozen']['version_id'], 'frozen_sha256': h['binding']['sha256'],
                'supersedes': clean(h['frozen']['config']['supersedes']),
                'submission_sha256': h['frozen']['artifacts']['submission.json'],
                'packet_sha256': h['frozen']['artifacts']['packet.json'],
                'validation': clean(h['frozen']['validation']),
                'source_packet': clean(h['packet']), 'judgment': clean(h['submission']),
                'unreadable_submission_preserved_by_hash': h['submission'] is None,
            } for h in history],
        }
        _safe(public_history)
        write_new(download/'source-reference-history.json', public_history)
        for name in ('readiness.json', 'readiness.csv', 'domain-counts.csv'):
            with (download/name).open('xb') as handle:
                handle.write((output/name).read_bytes())
        artifacts = []
        for artifact_id, name, label in (
                ('source-reference-history', 'source-reference-history.json', 'Complete source judgments, evidence, grouping, units, pairs and conditional methods'),
                ('readiness-json', 'readiness.json', 'All enrolled studies and preserved version history'),
                ('readiness-csv', 'readiness.csv', 'Study-level source-readiness table'),
                ('domain-counts', 'domain-counts.csv', 'Source-domain counts with allowed and enrolled denominators')):
            path = download/name
            if path.stat().st_size > plan['max_download_bytes']:
                raise ValueError('Source-reference download exceeds the frozen size limit')
            artifacts.append({'id': artifact_id, 'kind': path.suffix[1:], 'label': label,
                              'sha256': sha(path), 'download_path': 'downloads/'+name})
        if len({a['id'] for a in artifacts}) != len(artifacts):
            raise ValueError('Source-reference download IDs must be unique')
        trace = _source_reference_trace(plan, sha(plan_path), row, history, artifacts, clean)
        write_new(directory/'audit_trace.json', trace)
        traces.append({'study_id': sid, 'trace_sha256': sha(directory/'audit_trace.json'),
                       'source_validated': row['source_validated'],
                       'current_version_sha256': history[-1]['binding']['sha256']})
    for path, checksum in snapshot.items():
        if sha(path) != checksum:
            raise PipelineIntegrityError('Source-reference input changed during export')
    receipt = {
        'schema': 'clockbase-source-reference-public-export-v1', 'status': 'PASS',
        'plan_sha256': sha(plan_path), 'source_snapshot_sha256': _digest(snapshot),
        'enrolled_studies': table['enrolled_studies'], 'allowed_studies': table['allowed_studies'],
        'administrative_holds': table['administrative_holds'], 'traces': traces,
        'scientific_status': 'CANNOT VERIFY', 'accepted_scientific_repairs': 0, 'scientific_error_rate': None,
        'new_model_calls': 0, 'clocks_rerun': False,
        'output_sha256': {str(p.relative_to(output)): sha(p) for p in sorted(output.rglob('*')) if p.is_file()},
    }
    write_new(output/'verification.json', receipt)
    return {'status': 'PASS', 'traces': len(traces), 'enrolled_studies': table['enrolled_studies'],
            'source_validated_studies': plan['source_validated_studies'],
            'scientific_status': 'CANNOT VERIFY', 'new_model_calls': 0,
            'verification_sha256': sha(output/'verification.json')}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,help='Prepare a plan for this completed eight-stage pipeline run')
    parser.add_argument('--failed-run','--prepare-failed','--prepare-failed-run',dest='failed_run',type=Path,
                        help='Prepare a bounded plan for metadata+input complete and blinded-design failure only')
    parser.add_argument('--max-prompt-bytes',dest='prompt_bytes',type=int,
                        help='Failed-run prompt cap used to reconstruct a retry request; default is archived config or 48000')
    parser.add_argument('--plan',type=Path)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--corpus-run',type=Path,help='Export a completed dated corpus-audit report')
    parser.add_argument('--corpus-report',type=Path)
    parser.add_argument('--source-root',type=Path)
    parser.add_argument('--corpus-independent-review',type=Path)
    parser.add_argument('--future-version',type=Path,help='Export an executed future-ingestion version with its comprehensive audit')
    parser.add_argument('--future-readback',type=Path)
    parser.add_argument('--future-lineage',type=Path)
    parser.add_argument('--future-ingestion-readback',type=Path)
    parser.add_argument('--future-followup',type=Path,help='Export a sealed audit-only follow-up after a completed future-ingestion run')
    parser.add_argument('--parent-future-trace',type=Path,help='Preserved public trace for the original source audit')
    parser.add_argument('--followup-independent-review',type=Path)
    parser.add_argument('--audit-reassessment',type=Path,help='Export a new source-first audit after a preserved failed comprehensive hook')
    parser.add_argument('--prepare-source-reference', type=Path, help='Prepare a source-reference trace plan from a frozen release')
    parser.add_argument('--source-reference-release-sha256')
    parser.add_argument('--source-reference-registry', type=Path)
    parser.add_argument('--source-reference-registry-sha256')
    parser.add_argument('--source-reference-plan', type=Path, help='Export a pinned source-reference trace plan')
    args=parser.parse_args()
    reference_args = (args.prepare_source_reference, args.source_reference_release_sha256,
                      args.source_reference_registry, args.source_reference_registry_sha256, args.source_reference_plan)
    if any(value is not None for value in reference_args):
        if any(value is not None for value in (
                args.run, args.failed_run, args.prompt_bytes, args.corpus_run, args.corpus_report, args.source_root,
                args.corpus_independent_review, args.future_version, args.future_readback, args.future_lineage,
                args.future_ingestion_readback, args.future_followup, args.parent_future_trace,
                args.followup_independent_review, args.audit_reassessment)):
            parser.error('Source-reference preparation/export cannot be combined with another export mode')
        try:
            if args.source_reference_plan is not None:
                if args.output is None or args.plan is not None or any(v is not None for v in reference_args[:4]):
                    parser.error('--source-reference-plan requires only --output')
                result = export_source_reference(args.source_reference_plan, args.output)
            else:
                if None in reference_args[:4] or args.plan is None or args.output is not None:
                    parser.error('Source-reference preparation requires release/registry paths, both hashes and --plan')
                result = prepare_source_reference(*reference_args[:4], args.plan)
            print(json.dumps(result, allow_nan=False))
        except (OSError, ValueError, KeyError, TypeError, PipelineIntegrityError) as error:
            raise SystemExit(f'Public source-reference export stopped: {type(error).__name__}: {error}')
        return
    if args.audit_reassessment is not None:
        if (args.output is None or args.plan is None or any(value is not None for value in
                (args.future_followup,args.parent_future_trace,args.future_version,args.run,args.failed_run,args.corpus_run))):
            parser.error('--audit-reassessment requires --plan and --output without another export mode')
        try:
            print(json.dumps(export_audit_reassessment(args.audit_reassessment,args.output,plan_path=args.plan,
                independent_review_path=args.followup_independent_review)))
        except (OSError,ValueError,KeyError,TypeError,PipelineIntegrityError) as error:
            raise SystemExit(f'Public reassessment export stopped: {type(error).__name__}: {error}')
        return
    if args.future_followup is not None:
        if (args.output is None or args.parent_future_trace is None or any(value is not None for value in
                (args.future_version,args.run,args.failed_run,args.corpus_run,args.plan))):
            parser.error('--future-followup requires --parent-future-trace and --output, with no other export mode')
        try:
            print(json.dumps(export_future_followup(args.future_followup,args.parent_future_trace,args.output,
                independent_review_path=args.followup_independent_review)))
        except (OSError,ValueError,KeyError,TypeError,PipelineIntegrityError) as error:
            raise SystemExit(f'Public audit-only follow-up export stopped: {type(error).__name__}: {error}')
        return
    if args.parent_future_trace is not None or args.followup_independent_review is not None:
        parser.error('Follow-up trace arguments require --future-followup')
    if args.future_version is not None:
        if (None in (args.output,args.plan,args.future_readback,args.future_lineage)
                or any(value is not None for value in (args.run,args.failed_run,args.corpus_run,args.corpus_report,args.corpus_independent_review))):
            parser.error('--future-version requires --plan, --output, --future-readback and --future-lineage, with no other export mode')
        try:
            print(json.dumps(export_future_update(args.future_version,args.output,plan_path=args.plan,
                readback_path=args.future_readback,lineage_path=args.future_lineage,ingestion_readback_path=args.future_ingestion_readback)))
        except (OSError,ValueError,KeyError,TypeError,PipelineIntegrityError) as error:
            raise SystemExit(f'Public future-ingestion export stopped: {type(error).__name__}: {error}')
        return
    if any(value is not None for value in (args.future_readback,args.future_lineage,args.future_ingestion_readback)):
        parser.error('Future readbacks require --future-version')
    if args.corpus_run is not None:
        if args.output is None or args.run is not None or args.failed_run is not None or args.plan is not None:
            parser.error('--corpus-run requires --output and no other export mode')
        try:
            print(json.dumps(export_corpus(args.corpus_run,args.output,report_path=args.corpus_report,source_root=args.source_root,
                independent_review_path=args.corpus_independent_review)))
        except (OSError,ValueError,KeyError,TypeError,PipelineIntegrityError) as error:
            raise SystemExit(f'Public corpus export stopped: {type(error).__name__}: {error}')
        return
    if args.plan is None:
        parser.error('--plan is required for pipeline trace preparation/export')
    if args.run is not None and args.failed_run is not None:
        parser.error('Choose either --run or --failed-run')
    if args.prompt_bytes is not None and args.failed_run is None:
        parser.error('--max-prompt-bytes is only valid with --failed-run')
    if args.prompt_bytes is not None and args.prompt_bytes <= 0:
        parser.error('--max-prompt-bytes must be a positive integer')
    if (args.run is not None or args.failed_run is not None) == (args.output is not None):
        parser.error('Choose one prepare mode (--run/--failed-run) or --output to export')
    try:
        if args.run is not None:
            result = prepare(args.run, args.plan)
        elif args.failed_run is not None:
            result = prepare_failed(args.failed_run, args.plan, args.prompt_bytes)
        else:
            result = export(args.plan, args.output)
        print(json.dumps(result,allow_nan=False))
    except (OSError,ValueError,KeyError,TypeError,PipelineIntegrityError) as error:
        raise SystemExit(f'Public pipeline trace stopped: {type(error).__name__}: {error}')


if __name__=='__main__': main()
