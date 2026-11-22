"""Study-level error and coverage from independently adjudicated input records.

Rows are studies, not clocks or comparisons. Related studies share family_id.
No inference is made from agent agreement or successful software execution.
"""
import argparse
import csv
import fcntl
import gzip
import html
from collections import Counter, defaultdict
import hashlib
import json
import math
from itertools import combinations
from pathlib import Path
import random
import re
import time


COMPARISON_METHODS = ('existing_pipeline', 'deterministic', 'independent_audit', 'self_review')
ERROR_SIGNAL_DEFINITION = 'observed-producer-action-v1'
REFERENCE_VERSION = 'source-first-reference-2026-09-17-v2'
REFERENCE_DESIGN_DOMAINS = ('sample_identity', 'biological_units', 'comparison_definition', 'pairing', 'confounding')
REFERENCE_OUTCOME_DOMAINS = ('measurement_applicability', 'statistical_validity', 'numerical_accuracy',
                             'inferential_claim_scope', 'multiplicity_and_ranking')


SOURCE_READINESS_DOMAINS = REFERENCE_DESIGN_DOMAINS + (
    'assay_context', 'measurement_prerequisites', 'test_requirements', 'family_overlap')
SOURCE_READINESS_HOLDS = {
    'GSE295890': 'PROVIDER_REFUSAL_NO_NEW_SOURCE_TRANSMISSION',
    'GSE298220': 'PENDING_AUTHORIZATION_NO_NEW_SOURCE_TRANSMISSION'}
SOURCE_READINESS_FIELDS = {
    'series': frozenset(('title', 'overall_design', 'type', 'geo_accession',
                         'sample_id', 'platform_id', 'relation', 'pubmed_id')),
    'samples': frozenset(('title', 'geo_accession', 'source_name_ch1', 'organism_ch1',
        'molecule_ch1', 'characteristics_ch1', 'description', 'growth_protocol_ch1',
        'treatment_protocol_ch1', 'extract_protocol_ch1', 'data_processing',
        'library_selection', 'library_source', 'library_strategy', 'platform_id',
        'instrument_model', 'taxid_ch1', 'relation', 'series_id'))}
_SOURCE_READINESS_STAGE = '00_source_readiness'


def _readiness_json(path, expected=None):
    """Read a caller-declared administrative artifact, never a discovered file."""
    try:
        raw = Path(path).read_bytes()
        if expected is not None and hashlib.sha256(raw).hexdigest() != expected:
            raise ValueError('Source-readiness input hash changed: ' + str(path))
        def unique(pairs):
            result = {}
            for key, value in pairs:
                if key in result: raise ValueError('Duplicate JSON key: ' + key)
                result[key] = value
            return result
        return json.loads(raw, object_pairs_hook=unique,
                          parse_constant=lambda value: (_ for _ in ()).throw(ValueError('Non-finite JSON value')))
    except FileNotFoundError as error:
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path)) from error
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError('Cannot read source-readiness JSON: ' + str(path)) from error


def _readiness_enrollment(enrollment):
    """Check the controller roster without inspecting any source payload."""
    rows = enrollment.get('studies', [])
    if (enrollment.get('schema') != 'clockbase-rna80-source-readiness-enrollment-v1'
            or enrollment.get('split') != 'development_exposed' or len(rows) != 80
            or [r.get('manifest_ordinal') for r in rows] != list(range(1, 81))
            or len({r.get('study_id') for r in rows}) != 80):
        raise ValueError('The complete ordered RNA80 development enrollment is required')
    for row in rows:
        sid = row.get('study_id')
        held = sid in SOURCE_READINESS_HOLDS
        if (not isinstance(sid, str) or re.fullmatch(r'GSE[0-9]+', sid) is None
                or type(row.get('source_record_count_from_fixed_manifest')) is not int
                or row['source_record_count_from_fixed_manifest'] < 1
                or row.get('new_source_access_allowed') is not (not held)
                or row.get('administrative_state') != (SOURCE_READINESS_HOLDS[sid] if held else 'ENROLLED_NOT_ASSESSED')
                or row.get('reference_judgment') is not None):
            raise ValueError('RNA80 administrative identity, hold or source-record count changed')
    allowed = [r for r in rows if r['new_source_access_allowed']]
    held = [r for r in rows if not r['new_source_access_allowed']]
    counts = {'enrolled_studies': 80, 'enrolled_source_records': sum(r['source_record_count_from_fixed_manifest'] for r in rows),
              'source_access_allowed_studies': len(allowed),
              'source_access_allowed_records': sum(r['source_record_count_from_fixed_manifest'] for r in allowed),
              'administrative_holds': len(held),
              'administrative_hold_source_records': sum(r['source_record_count_from_fixed_manifest'] for r in held)}
    if len(allowed) != 78 or {r['study_id'] for r in held} != set(SOURCE_READINESS_HOLDS) or any(
            enrollment.get(key) != value for key, value in counts.items()):
        raise ValueError('RNA80 enrollment denominators or typed holds changed')
    return rows


def source_readiness_context(release_path, release_sha256):
    """Require a parent-frozen source-only release; preparation drafts do not release sources."""
    from importlib.metadata import version
    from . import agent_workflow_2026_09_10 as workflow
    from src import pipeline_runtime_2026_09_10 as runtime
    from jsonschema import Draft202012Validator
    if not isinstance(release_sha256, str) or re.fullmatch(r'[0-9a-f]{64}', release_sha256) is None:
        raise ValueError('An explicit frozen release SHA256 is required')
    release = _readiness_json(release_path, release_sha256)
    if (release.get('schema') != 'clockbase-rna80-source-readiness-release-v1'
            or release.get('phase') != 'SOURCE_ONLY_FROZEN' or release.get('split') != 'development_exposed'
            or re.fullmatch(r'[0-9a-f]{40}', release.get('protocol_commit', '')) is None
            or not isinstance(release.get('reference_runtime'), dict)
            or set(release['reference_runtime']) != {'runtime', 'model'}
            or not isinstance(release['reference_runtime']['runtime'], str) or not release['reference_runtime']['runtime'].strip()
            or (release['reference_runtime']['model'] is not None and
                (not isinstance(release['reference_runtime']['model'], str) or not release['reference_runtime']['model'].strip()))
            or release.get('jsonschema_version') != version('jsonschema')):
        raise ValueError('A committed, runtime-bound SOURCE_ONLY_FROZEN release is required')
    code = {'audit_system/evaluation_2026_09_10.py': workflow.sha(__file__),
            'audit_system/agent_workflow_2026_09_10.py': workflow.sha(workflow.__file__),
            'src/pipeline_runtime_2026_09_10.py': workflow.sha(runtime.__file__)}
    if release.get('implementation_sha256') != code:
        raise ValueError('Source-readiness implementation differs from its frozen release')
    bindings = release.get('bindings', {})
    if set(bindings) != {'protocol', 'enrollment', 'allowed_sources', 'schedule', 'judgment_schema'}:
        raise ValueError('Release must bind protocol, enrollment, allowed sources, schedule and schema')
    for binding in bindings.values():
        if (not isinstance(binding, dict) or set(binding) != {'path', 'sha256'}
                or not isinstance(binding['path'], str) or not Path(binding['path']).is_absolute()
                or re.fullmatch(r'[0-9a-f]{64}', binding.get('sha256', '')) is None):
            raise ValueError('Malformed administrative release binding')
    try:
        if workflow.sha(bindings['protocol']['path']) != bindings['protocol']['sha256']:
            raise ValueError('Frozen source-readiness protocol changed')
    except OSError as error:
        raise ValueError('Cannot read frozen source-readiness protocol') from error
    loaded = {key: _readiness_json(value['path'], value['sha256']) for key, value in bindings.items() if key != 'protocol'}
    rows = _readiness_enrollment(loaded['enrollment'])
    allowed = loaded['allowed_sources']
    if (allowed.get('schema') != 'clockbase-rna80-source-readiness-allowed-sources-v1'
            or allowed.get('split') != 'development_exposed'
            or allowed.get('source_manifest_sha256') != loaded['enrollment'].get('source_manifest_sha256')
            or [r.get('study_id') for r in allowed.get('studies', [])] != [r['study_id'] for r in rows if r['new_source_access_allowed']]):
        raise ValueError('Allowed source manifest changed the fixed allowed roster')
    by_id, paths = {}, set()
    for row in allowed['studies']:
        original = rows[row['manifest_ordinal'] - 1]
        samples = row.get('sample_ids', [])
        if (original['study_id'] != row['study_id'] or row['study_id'] in SOURCE_READINESS_HOLDS
                or len(samples) != original['source_record_count_from_fixed_manifest']
                or len(set(samples)) != len(samples) or any(re.fullmatch(r'GSM[0-9]+', s) is None for s in samples)
                or set(row.get('sources', {})) != {'series', 'samples'}):
            raise ValueError('Allowed source membership or administrative hold mismatch')
        for binding in row['sources'].values():
            path = binding.get('path')
            if (set(binding) != {'path', 'sha256', 'bytes'} or not isinstance(path, str)
                    or not Path(path).is_absolute() or '..' in Path(path).parts or path in paths
                    or type(binding.get('bytes')) is not int or binding['bytes'] < 1
                    or re.fullmatch(r'[0-9a-f]{64}', binding.get('sha256', '')) is None):
                raise ValueError('Malformed or duplicate allowlisted raw-source identity')
            paths.add(path)
        by_id[row['study_id']] = row
    schedule = loaded['schedule']
    expected = list(by_id)
    batches = schedule.get('batches', [])
    if (schedule.get('schema') != 'clockbase-rna80-source-readiness-schedule-v1' or schedule.get('batch_size') != 5
            or len(batches) != 16):
        raise ValueError('Fixed RNA80 five-study batch schedule required')
    for index, batch in enumerate(batches):
        ids = expected[index * 5:(index + 1) * 5]
        if (batch.get('batch_id') != 'source-reference-batch-%03d' % (index + 1)
                or batch.get('study_ids') != ids
                or batch.get('manifest_ordinals') != [by_id[s]['manifest_ordinal'] for s in ids]
                or batch.get('source_record_count') != sum(len(by_id[s]['sample_ids']) for s in ids)):
            raise ValueError('Batch schedule changes the fixed order, membership or denominator')
    schema = loaded['judgment_schema']
    Draft202012Validator.check_schema(schema)
    if (schema.get('$schema') != 'https://json-schema.org/draft/2020-12/schema'
            or schema.get('properties', {}).get('study_id', {}).get('enum') != expected):
        raise ValueError('Judgment schema does not match the fixed allowed roster')
    return {**loaded, 'release': release, 'release_sha256': release_sha256,
            'release_path': str(Path(release_path).absolute()), 'allowed_by_id': by_id,
            'batches': {b['batch_id']: b for b in batches}}


def _readiness_batch(context, batch_id):
    batch = context['batches'].get(batch_id)
    if batch is None or any(s in SOURCE_READINESS_HOLDS or s not in context['allowed_by_id'] for s in batch['study_ids']):
        raise ValueError('Unlisted or administratively held source batch; no raw source was opened')
    return batch


def _readiness_raw(context, batch_id, study_id, role, binding):
    """Reject access first, including before stat, symlink checks or checksumming."""
    batch = _readiness_batch(context, batch_id)
    if (study_id in SOURCE_READINESS_HOLDS or study_id not in batch['study_ids']
            or role not in {'series', 'samples'}
            or binding != context['allowed_by_id'][study_id]['sources'][role]):
        raise ValueError('Unlisted or held raw source/path/hash; no raw source was opened')
    import os
    import stat
    path = Path(binding['path'])
    try:
        if any(p.is_symlink() for p in (path, *path.parents)):
            raise ValueError('Symlinked raw sources are not released by the literal allowlist')
        with os.fdopen(os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK), 'rb') as handle:
            if not stat.S_ISREG(os.fstat(handle.fileno()).st_mode):
                raise ValueError('Released raw source is not a regular file')
            raw = handle.read(binding['bytes'] + 1)
    except FileNotFoundError as error:
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path)) from error
    except OSError as error:
        raise ValueError('Cannot read released raw source: ' + str(path)) from error
    if len(raw) != binding['bytes'] or hashlib.sha256(raw).hexdigest() != binding['sha256']:
        raise ValueError('Released raw source size/hash changed: ' + str(path))
    return raw


def _readiness_soft(raw, role, study_id, sample_ids, binding):
    """Parse complete SOFT fields; repeated occurrences and line breaks are literal.

    This adopts the old extractor's field/entity boundary approach, without its
    historical case-specific main or validator. Blank continuation lines and
    CRLF within values are retained. No source summaries or contacts are exported.
    """
    try:
        lines = raw.decode('utf-8').splitlines(keepends=True)
    except UnicodeError as error:
        raise ValueError('Released SOFT source is not UTF-8') from error
    prefix, kind = ('!Series_', 'SERIES') if role == 'series' else ('!Sample_', 'SAMPLE')
    expected = [study_id] if role == 'series' else sample_ids
    entities, records, counts = [], [], Counter()
    entity, current, previous_eol, table = None, None, '', False
    for number, physical in enumerate(lines, 1):
        eol = '\r\n' if physical.endswith('\r\n') else ('\n' if physical.endswith('\n') else '')
        line = physical[:-len(eol)] if eol else physical
        if '\r' in line or any(c in line for c in ('\v', '\f', '\x85', '\u2028', '\u2029')):
            raise ValueError('Malformed SOFT line terminator')
        if line.startswith('^'):
            match = re.fullmatch(r'\^' + kind + r' = (GSE[0-9]+|GSM[0-9]+)', line)
            if not match or table or match[1] not in expected or match[1] in entities:
                raise ValueError('Malformed, duplicate or unlisted SOFT entity boundary')
            entity = match[1]; entities.append(entity); current = None
        elif line.startswith('!'):
            current = None
            marker = re.fullmatch(r'!sample_table_(begin|end)', line)
            if marker:
                if role != 'samples' or entity is None or (marker[1] == 'begin') == table:
                    raise ValueError('Malformed SOFT table boundary')
                table = marker[1] == 'begin'
            else:
                match = re.fullmatch(r'(' + prefix + r'[A-Za-z0-9_]+(?:/[A-Za-z0-9_]+)*) = (.*)', line)
                if match is None or entity is None or table:
                    raise ValueError('Malformed SOFT field boundary (%s, %s, line %d)' % (role, entity, number))
                field, value = match.groups()
                short = field[len(prefix):]
                occurrence = counts[(entity, short)]; counts[(entity, short)] += 1
                current = {'source_role': role, 'source_path': binding['path'], 'source_sha256': binding['sha256'],
                           'entity': entity, 'field': field, 'line_start': number, 'line_end': number, 'quote': value,
                           'projection_field': short + ('.' + str(occurrence) if occurrence else '')}
                if short in SOURCE_READINESS_FIELDS[role]: records.append(current)
        elif line.startswith('#') or table:
            current = None
        elif current is not None:
            current['quote'] += previous_eol + line
            current['line_end'] = number
        elif line:
            raise ValueError('Unattributed SOFT text outside a field/table')
        previous_eol = eol
    if table or set(entities) != set(expected):
        raise ValueError('SOFT entities differ from the fixed unique roster')
    for sid in expected:
        accessions = [r['quote'] for r in records if r['entity'] == sid and r['field'] == prefix + 'geo_accession']
        if accessions != [sid]: raise ValueError('SOFT accession does not match its entity')
    if role == 'series':
        members = [r['quote'] for r in records if r['field'] == '!Series_sample_id']
        if len(members) != len(set(members)) or set(members) != set(sample_ids):
            raise ValueError('Series GSM membership differs from the fixed unique sample roster')
    return records


def _readiness_packet(context, batch_id, study_id):
    batch = _readiness_batch(context, batch_id)
    if study_id not in batch['study_ids'] or study_id in SOURCE_READINESS_HOLDS:
        raise ValueError('Study not released in this fixed batch; no raw source was opened')
    row = context['allowed_by_id'][study_id]
    records = []
    for role, binding in row['sources'].items():
        records.extend(_readiness_soft(_readiness_raw(context, batch_id, study_id, role, binding),
                                      role, study_id, row['sample_ids'], binding))
    samples = {s: {'sample_id': s} for s in row['sample_ids']}
    series, witnesses, projections = {}, {}, {}
    for record in records:
        # The shared evidence validator reserves sample_id for an identity, even
        # in document fields. Keep Series fields in their explicit namespace.
        field = record['projection_field'] if record['source_role'] == 'samples' else 'Series_' + record['projection_field']
        witness = {key: value for key, value in record.items() if key != 'projection_field'}
        eid = study_id + ':' + record['source_role'] + ':L' + str(record['line_start'])
        witnesses[eid] = witness
        if record['source_role'] == 'samples':
            samples[record['entity']][field] = record['quote']
            projections[eid] = {'sample_id': record['entity'], 'field': field, 'quote': record['quote']}
        else:
            series[field] = record['quote']
            projections[eid] = {'source_id': study_id + ':series', 'field': field, 'quote': record['quote']}
    bindings = context['release']['bindings']
    packet = {'schema': 'clockbase-rna80-source-readiness-packet-v1', 'phase': 'SOURCE_ONLY',
              'split': 'development_exposed', 'release_sha256': context['release_sha256'],
              'protocol_sha256': bindings['protocol']['sha256'], 'enrollment_sha256': bindings['enrollment']['sha256'],
              'judgment_schema_sha256': bindings['judgment_schema']['sha256'], 'batch_id': batch_id,
              'study_id': study_id, 'sample_ids': row['sample_ids'], 'sources': row['sources'],
              'source_packet': {'study_id': study_id, 'samples': list(samples.values()),
                                'study_documents': [{'source_id': study_id + ':series', 'fields': series}]},
              'literal_witnesses': witnesses,
              'exposure_note': 'Prior pipeline outcomes are withheld; allowed public source narratives may contain published findings.'}
    return packet, projections


def _readiness_store(output, config, operation):
    """Use the existing append-only stage receipts; serialize writers to one version."""
    from src.pipeline_runtime_2026_09_10 import StageStore
    output = Path(output).absolute()
    for value in config.get('input_paths', []):
        if Path(value).absolute().is_relative_to(output):
            raise ValueError('Source-readiness output must not contain its input artifacts')
    output.mkdir(parents=True, exist_ok=True)
    with (output / '.source-readiness.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return StageStore(output, config).run(_SOURCE_READINESS_STAGE, operation)


def prepare_source_readiness_batch(release_path, release_sha256, batch_id, output):
    """Project only a scheduled batch after release; forms contain no default verdicts."""
    from .agent_workflow_2026_09_10 import sha, write_new
    context = source_readiness_context(release_path, release_sha256)
    batch = _readiness_batch(context, batch_id)
    config = {'schema': 'clockbase-source-readiness-projection-store-v1', 'release_sha256': release_sha256,
              'frozen_release': context['release'],
              'batch_id': batch_id, 'sources': [context['allowed_by_id'][s] for s in batch['study_ids']],
              'input_paths': [str(release_path), *[b['path'] for b in context['release']['bindings'].values()],
                              *[b['path'] for s in batch['study_ids'] for b in context['allowed_by_id'][s]['sources'].values()]]}
    # Even a resumed projection checks its original sources, without exposing other batches.
    packets = {sid: _readiness_packet(context, batch_id, sid)[0] for sid in batch['study_ids']}
    def prepare(directory):
        result = {'schema': 'clockbase-source-readiness-projection-result-v1', 'batch_id': batch_id,
                  'release_sha256': release_sha256, 'packets': [], 'reference_labels_completed': False,
                  'new_model_calls': 0, 'source_files_read': 2 * len(batch['study_ids'])}
        for sid, packet in packets.items():
            folder = directory / sid; folder.mkdir()
            path = folder / 'packet.json'; write_new(path, packet)
            form = {'schema': 'clockbase-rna80-source-readiness-judgment-v1', 'completed': False,
                    'split': 'development_exposed', 'protocol_sha256': packet['protocol_sha256'],
                    'enrollment_sha256': packet['enrollment_sha256'], 'packet_sha256': sha(path),
                    'study_id': sid, 'sample_ids': packet['sample_ids'],
                    'reviewer': {'name': '', 'role': 'independent_source_adjudication', 'independent_from_methods': None,
                                 'runtime_model_as_reported': None, 'known_prior_exposure': '', 'source_narrative_exposure_note': ''},
                    'domains': {name: None for name in REFERENCE_DESIGN_DOMAINS},
                    'readiness_domains': {name: None for name in SOURCE_READINESS_DOMAINS if name not in REFERENCE_DESIGN_DOMAINS},
                    'facts': [], 'evidence': {}, 'unresolved': [], 'comparison_inventory_status': None,
                    'comparisons': [], 'unassigned_sample_ids': [], 'unit_map': [], 'unit_unresolved_sample_ids': [],
                    'pair_map': [], 'covariates': [], 'family_links': [], 'family_closure': 'UNRESOLVED_SOURCE_ONLY',
                    'additional_sources_needed': []}
            write_new(folder / 'uncompleted-form.json', form)
            result['packets'].append({'study_id': sid, 'path': str(path), 'sha256': sha(path),
                                      'form': str(folder / 'uncompleted-form.json')})
        write_new(directory / 'projection.json', result)
        return result
    return _readiness_store(output, config, prepare)


def _readiness_validate_submission(context, packet, packet_sha256, judgment):
    """Mechanical source attribution and structural integrity, not scientific adjudication."""
    from jsonschema import Draft202012Validator
    from .agent_workflow_2026_09_10 import _check_evidence, comparison_identity
    errors = sorted(Draft202012Validator(context['judgment_schema']).iter_errors(judgment), key=lambda e: str(list(e.path)))
    if errors:
        error = errors[0]
        raise ValueError('Submission schema: /' + '/'.join(map(str, error.path)) + ': ' + error.message)
    if (judgment['study_id'] != packet['study_id'] or judgment['sample_ids'] != packet['sample_ids']
            or judgment['packet_sha256'] != packet_sha256
            or any(judgment[key] != packet[key] for key in ('protocol_sha256', 'enrollment_sha256'))
            or judgment['reviewer']['runtime_model_as_reported'] != context['release']['reference_runtime']['model']
            or any(not judgment['reviewer'][key].strip() for key in ('name', 'known_prior_exposure', 'source_narrative_exposure_note'))):
        raise ValueError('Submission packet, roster, release or declared reviewer identity mismatch')
    roster = set(packet['sample_ids'])
    def unique(items, key):
        result = {r[key]: r for r in items}
        if len(result) != len(items): raise ValueError('Duplicate ' + key)
        return result
    facts = unique(judgment['facts'], 'fact_id')
    unresolved = unique(judgment['unresolved'], 'unresolved_id')
    comparisons = unique(judgment['comparisons'], 'comparison_id')
    units = unique(judgment['unit_map'], 'unit_id')
    evidence = judgment['evidence']
    # Check all declared identities and cross-references before opening any raw source.
    def walk(value):
        if isinstance(value, dict):
            for key, child in value.items():
                targets = {'fact_ids': facts, 'unresolved_ids': unresolved, 'comparison_ids': comparisons,
                           'evidence_ids': evidence, 'related_evidence_ids': evidence}
                if key in targets and set(child) - targets[key].keys(): raise ValueError('Unknown ' + key)
                if key.endswith('sample_ids') or key in {'case_ids', 'control_ids'}:
                    if set(child) - roster: raise ValueError('Sample identity outside packet roster: ' + key)
                if key in {'case_sample_id', 'control_sample_id'} and child not in roster:
                    raise ValueError('Pair sample identity outside packet roster')
                walk(child)
        elif isinstance(value, list):
            for child in value: walk(child)
    walk(judgment)
    for witness in evidence.values():
        binding = packet['sources'].get(witness['source_role'])
        if (binding is None or witness['source_path'] != binding['path'] or witness['source_sha256'] != binding['sha256']
                or witness['line_end'] < witness['line_start']
                or witness['entity'] not in ([packet['study_id']] if witness['source_role'] == 'series' else roster)):
            raise ValueError('Unlisted witness path/hash/entity; no raw source was opened')
    rebuilt, projections = _readiness_packet(context, packet['batch_id'], packet['study_id'])
    if rebuilt != packet: raise ValueError('Source packet differs from its literal raw-source reconstruction')
    canonical = {json.dumps(v, sort_keys=True): key for key, v in rebuilt['literal_witnesses'].items()}
    checked = {}
    for key, witness in evidence.items():
        match = canonical.get(json.dumps(witness, sort_keys=True))
        if match is None: raise ValueError('Witness must equal one complete literal field at its exact line interval')
        checked[key] = projections[match]
        _check_evidence([checked[key]], packet['source_packet'])
    def supported(ids, *, unit_key=False):
        entries = [checked[key] for key in ids]
        if unit_key:
            entries = [e for e in entries if e['field'].split('.')[0] not in {'title', 'geo_accession'}]
        return _check_evidence(entries, packet['source_packet'])[1]
    def links(row, domain=None, comparison_id=None):
        if not row.get('fact_ids') and not row.get('unresolved_ids'):
            raise ValueError('A conditional requirement needs source facts or an unresolved dependency')
        for key in row.get('fact_ids', []):
            if comparison_id is not None and facts[key]['comparison_ids'] and comparison_id not in facts[key]['comparison_ids']:
                raise ValueError('Source fact is scoped to a different comparison')
            if domain is not None and facts[key]['domain'] != domain:
                raise ValueError('Fact has a different domain from the assertion')
            if 'evidence_ids' in row and not set(facts[key]['evidence_ids']).issubset(row['evidence_ids']):
                raise ValueError('Assertion omits its source-fact witnesses')
        for key in row.get('unresolved_ids', []):
            if comparison_id is not None and unresolved[key]['comparison_ids'] and comparison_id not in unresolved[key]['comparison_ids']:
                raise ValueError('Unresolved dependency is scoped to a different comparison')
            if domain is not None and unresolved[key]['domain'] != domain:
                raise ValueError('Unresolved dependency has a different domain from the assertion')
    for row in unresolved.values():
        for field in row['searched_fields']:
            role = 'series' if field.startswith('!Series_') else 'samples'
            prefix = '!Series_' if role == 'series' else '!Sample_'
            if role not in row['searched_source_roles'] or not field.startswith(prefix) or field[len(prefix):] not in SOURCE_READINESS_FIELDS[role]:
                raise ValueError('Unresolved search must name released SOFT fields and source roles')
    for fact in facts.values():
        if set(fact['sample_ids']) - supported(fact['evidence_ids']):
            raise ValueError('A sample-specific fact lacks a direct source witness')
    for row in [*facts.values(), *unresolved.values()]:
        for cid in row['comparison_ids']:
            contrast = comparisons[cid]
            if set(row['sample_ids']) - set(contrast['case_ids'] + contrast['control_ids'] + contrast['excluded_sample_ids']):
                raise ValueError('A comparison-scoped fact/dependency names samples outside its arms/exclusions')
    def domains(rows, comparison_id=None):
        for domain, row in rows.items():
            links(row, domain, comparison_id)
            if not row['rationale'].strip(): raise ValueError('Empty scoped domain rationale')
            if row['status'] == 'FAIL' and type(row['material']) is not bool:
                raise ValueError('FAIL requires explicit scoped materiality')
            if row['status'] in {'PASS', 'FAIL'} and any(
                    not set(facts[f]['evidence_ids']).issubset(row['evidence_ids']) for f in row['fact_ids']):
                raise ValueError('Resolved domain omits its supporting fact witnesses')
    domains(judgment['domains']); domains(judgment['readiness_domains'])
    assigned = set()
    for cid, row in comparisons.items():
        members = set(row['case_ids'] + row['control_ids'])
        if (set(row['case_ids']) & set(row['control_ids']) or members & set(row['excluded_sample_ids'])
                or cid != comparison_identity(packet['study_id'], row['case_ids'], row['control_ids'], row['target'])
                or members - supported(row['evidence_ids'])):
            raise ValueError('Comparison identity, disjoint arms/exclusions or direct source coverage failed')
        domains(row['domains'], cid)
        for requirement in row['analysis_requirements']: links(requirement, comparison_id=cid)
        assigned.update(members)
    disposition = judgment['comparison_inventory_status']
    if ((disposition == 'SOURCE_DEFINED_CONTRASTS_LISTED') != bool(comparisons)
            or set(judgment['unassigned_sample_ids']) != roster - assigned):
        raise ValueError('Comparison disposition or complete unassigned sample accounting failed')
    def missing_coverage(ids, domain):
        covered = set().union(*(set(r['sample_ids']) for r in unresolved.values() if r['domain'] == domain))
        if set(ids) - covered: raise ValueError('Unassigned samples require explicit ' + domain + ' dependencies')
    missing_coverage(judgment['unassigned_sample_ids'], 'comparison_definition')
    mapped = set()
    for unit in units.values():
        members = set(unit['sample_ids'])
        links(unit, 'biological_units')
        if (mapped & members or members - supported(unit['evidence_ids'], unit_key=True)
                or any(facts[f]['qualification'] != 'LITERAL_SOURCE_DECLARATION' for f in unit['fact_ids'])
                or not members.issubset(set().union(*(set(facts[f]['sample_ids']) for f in unit['fact_ids'])))):
            raise ValueError('Unit map lacks explicit sample-linked unit evidence or overlaps units')
        mapped.update(members)
    if set(judgment['unit_unresolved_sample_ids']) != roster - mapped:
        raise ValueError('Known and unresolved unit samples must partition the roster')
    missing_coverage(judgment['unit_unresolved_sample_ids'], 'biological_units')
    paired = set()
    for pair in judgment['pair_map']:
        cid, uid = pair['comparison_id'], pair['unit_id']
        links(pair, 'pairing', cid)
        key = (cid, pair['case_sample_id'], pair['control_sample_id'])
        members = {pair['case_sample_id'], pair['control_sample_id']}
        if (cid not in comparisons or uid not in units or key in paired
                or pair['case_sample_id'] not in comparisons[cid]['case_ids']
                or pair['control_sample_id'] not in comparisons[cid]['control_ids']
                or not members.issubset(units[uid]['sample_ids'])
                or members - supported(pair['evidence_ids'], unit_key=True)
                or any(facts[f]['qualification'] != 'LITERAL_SOURCE_DECLARATION' for f in pair['fact_ids'])
                or not members.issubset(set().union(*(set(facts[f]['sample_ids']) for f in pair['fact_ids'])))):
            raise ValueError('Pair map lacks explicit unit/arm/source bindings; suffix matching is insufficient')
        paired.add(key)
    for row in comparisons.values():
        if row['conditional_test_family'] == 'PAIRED':
            pairs = [p for p in judgment['pair_map'] if p['comparison_id'] == row['comparison_id']]
            for field, arm in (('case_sample_id', 'case_ids'), ('control_sample_id', 'control_ids')):
                ids = [p[field] for p in pairs]
                if len(ids) != len(set(ids)) or set(ids) != set(row[arm]):
                    raise ValueError('A paired conditional plan requires a complete one-to-one declared pair map')
    names = set()
    for row in judgment['covariates']:
        key = (row['name'], row['scope']); observed, missing = set(row['observed_sample_ids']), set(row['missing_sample_ids'])
        links(row, 'confounding')
        if key in names or observed & missing or observed | missing != roster:
            raise ValueError('Covariate observed/missing IDs must partition the roster without duplicates')
        if observed and (not row['fact_ids'] or observed - supported(row['evidence_ids'])
                         or not observed.issubset(set().union(*(set(facts[f]['sample_ids']) for f in row['fact_ids'])))):
            raise ValueError('Observed covariate values need sample-linked source facts')
        if missing and (not row['unresolved_ids'] or not missing.issubset(set().union(
                *(set(unresolved[u]['sample_ids']) for u in row['unresolved_ids'])))):
            raise ValueError('Missing covariates need sample-specific unresolved dependencies')
        names.add(key)
    family_keys = set()
    for row in judgment['family_links']:
        key = (row['link_type'], row['target_identifier'])
        quotes = [evidence[e]['quote'] for e in row['evidence_ids']]
        if key in family_keys or not any(row['target_identifier'] in quote for quote in quotes):
            raise ValueError('Family identifier must occur literally in its cited source; duplicate link')
        family_keys.add(key)
    return {'validation_status': 'ACCEPTED_SOURCE_ATTRIBUTION', 'study_id': packet['study_id'],
            'domains': {**judgment['domains'], **judgment['readiness_domains']},
            'literal_witnesses_checked': len(evidence), 'comparisons': len(comparisons),
            'scope': 'Schema, literal attribution and structural bindings only. Scientific interpretation, contradiction and unit-key semantics are not independently certified.',
            'family_closure': 'UNRESOLVED_SOURCE_ONLY', 'new_model_calls': 0, 'outcome_values_read': False}


def _readiness_frozen(context, binding):
    """Read one pinned completed readiness version, without reopening raw sources."""
    from .agent_workflow_2026_09_10 import sha
    path = Path(binding['path'])
    if (path.name != 'frozen-source-readiness.json' or not re.fullmatch(r'attempt-[0-9]{3}', path.parent.name)
            or path.parent.parent.name != _SOURCE_READINESS_STAGE):
        raise ValueError('A completed source-readiness artifact path is required')
    frozen = _readiness_json(path, binding['sha256'])
    stage = path.parent.parent
    complete = _readiness_json(stage / 'complete.json', (stage / 'complete.sha256').read_text().strip())
    config = _readiness_json(stage.parent / 'config.json')
    if (config.get('schema') != 'clockbase-source-readiness-submission-store-v1'
            or config.get('release_sha256') != context['release_sha256']
            or frozen.get('config') != config or frozen.get('study_id') in SOURCE_READINESS_HOLDS
            or frozen.get('study_id') not in context['allowed_by_id']
            or complete.get('result', {}).get('sha256') != binding['sha256']
            or complete['result'].get('frozen_submission') != str(path)):
        raise ValueError('Readiness version differs from its frozen stage/release identity')
    for name, checksum in complete['outputs'].items():
        rel = Path(name)
        if rel.is_absolute() or '..' in rel.parts or rel.parts[0] != path.parent.name:
            raise ValueError('Unlisted readiness stage output path')
        if sha(stage / rel) != checksum: raise ValueError('Frozen readiness stage output changed')
    if complete['outputs'].get(str(path.relative_to(stage))) != binding['sha256']:
        raise ValueError('Readiness result is not bound in its stage outputs')
    if set(frozen.get('artifacts', {})) != {'submission.json', 'packet.json', 'literal-validation.json'} or any(
            complete['outputs'].get(path.parent.name + '/' + name) != checksum for name, checksum in frozen['artifacts'].items()):
        raise ValueError('Readiness input/validation artifacts are not bound in its completed stage')
    if _readiness_json(path.parent / 'literal-validation.json') != frozen.get('validation'):
        raise ValueError('Frozen readiness validation receipt differs from its embedded result')
    return frozen


def freeze_source_readiness_submission(release_path, release_sha256, batch_id, study_id,
                                      packet_binding, submission_path, version_id, output, *, supersedes=None):
    """Keep accepted and rejected submissions in separately named input-bound versions."""
    from .agent_workflow_2026_09_10 import sha, write_new
    context = source_readiness_context(release_path, release_sha256)
    batch = _readiness_batch(context, batch_id)
    if study_id in SOURCE_READINESS_HOLDS or study_id not in batch['study_ids']:
        raise ValueError('Study is not released in this batch; no packet/submission/raw file was opened')
    if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,79}', version_id) is None:
        raise ValueError('A stable named submission version is required')
    prior = _readiness_frozen(context, supersedes) if supersedes is not None else None
    if prior and (prior['study_id'] != study_id or prior['version_id'] == version_id):
        raise ValueError('Superseded submission must be the same study with a different version name')
    packet = _readiness_json(packet_binding['path'], packet_binding['sha256'])
    row = context['allowed_by_id'][study_id]
    if (packet.get('schema') != 'clockbase-rna80-source-readiness-packet-v1' or packet.get('study_id') != study_id
            or packet.get('batch_id') != batch_id or packet.get('release_sha256') != release_sha256
            or packet.get('sources') != row['sources'] or packet.get('sample_ids') != row['sample_ids']):
        raise ValueError('Unlisted packet study/batch/raw-source binding; no raw source was opened')
    try:
        submission_raw = Path(submission_path).read_bytes()
    except FileNotFoundError as error:
        raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(submission_path)) from error
    except OSError as error:
        raise ValueError('Cannot read submitted source-readiness judgment') from error
    config = {'schema': 'clockbase-source-readiness-submission-store-v1', 'release_sha256': release_sha256,
              'frozen_release': context['release'],
              'study_id': study_id, 'batch_id': batch_id, 'version_id': version_id,
              'packet': packet_binding, 'submission_sha256': hashlib.sha256(submission_raw).hexdigest(),
              'sources': row['sources'], 'supersedes': supersedes,
              'input_paths': [str(release_path), packet_binding['path'], str(Path(submission_path).absolute()),
                              *[b['path'] for b in context['release']['bindings'].values()],
                              *[b['path'] for b in row['sources'].values()]]}
    def freeze(directory):
        with (directory / 'submission.json').open('xb') as handle: handle.write(submission_raw)
        write_new(directory / 'packet.json', packet)
        try:
            judgment = _readiness_json(directory / 'submission.json')
            validation = _readiness_validate_submission(context, packet, packet_binding['sha256'], judgment)
        except (ValueError, OSError, KeyError, TypeError) as error:
            validation = {'validation_status': 'REJECTED', 'study_id': study_id, 'readiness_status': 'NOT_ASSESSED',
                          'reason': type(error).__name__ + ': ' + str(error), 'new_model_calls': 0,
                          'scope': 'No scientific domain verdict can be counted from an invalid submission.'}
        write_new(directory / 'literal-validation.json', validation)
        value = {'schema': 'clockbase-frozen-source-readiness-v1', 'phase': 'SOURCE_ATTRIBUTION_FROZEN',
                 'study_id': study_id, 'version_id': version_id, 'split': 'development_exposed',
                 'config': config, 'validation': validation,
                 'artifacts': {n: sha(directory / n) for n in ('submission.json', 'packet.json', 'literal-validation.json')}}
        path = directory / 'frozen-source-readiness.json'; write_new(path, value)
        return {'study_id': study_id, 'version_id': version_id, 'validation_status': validation['validation_status'],
                'frozen_submission': str(path), 'sha256': sha(path)}
    # Recheck already-frozen inputs before permitting hash-checked resume.
    if (Path(output) / _SOURCE_READINESS_STAGE / 'complete.json').exists():
        _readiness_packet(context, batch_id, study_id)
    return _readiness_store(output, config, freeze)


def source_readiness_table(enrollment, versions=()):
    """All-enrolled administrative/readiness counts; never an accuracy estimator."""
    rows = _readiness_enrollment(enrollment)
    current, history, seen = {}, defaultdict(list), set()
    enrolled_ids = {r['study_id'] for r in rows}
    for binding, frozen in versions:
        sid = frozen['study_id']; checksum = binding['sha256']
        if sid not in enrolled_ids or sid in SOURCE_READINESS_HOLDS or checksum in seen or frozen['validation']['validation_status'] not in {'REJECTED', 'ACCEPTED_SOURCE_ATTRIBUTION'}:
            raise ValueError('Duplicate, held or unknown readiness version')
        prior = frozen['config']['supersedes']
        if (prior is None and sid in current) or (prior is not None and (sid not in current or prior != current[sid][0])):
            raise ValueError('Readiness history must retain the complete ordered supersession chain')
        if frozen['version_id'] in {v['version_id'] for v in history[sid]}:
            raise ValueError('Duplicate named readiness version')
        seen.add(checksum); current[sid] = (binding, frozen)
        history[sid].append({'version_id': frozen['version_id'], 'binding': binding,
                             'validation_status': frozen['validation']['validation_status']})
    studies = []
    for row in rows:
        sid = row['study_id']; held = not row['new_source_access_allowed']
        frozen = current[sid][1] if sid in current else None
        validation = frozen['validation'] if frozen else None
        accepted = bool(validation and validation['validation_status'] == 'ACCEPTED_SOURCE_ATTRIBUTION')
        domains = {d: validation['domains'][d]['status'] if accepted else 'NOT_ASSESSED' for d in SOURCE_READINESS_DOMAINS}
        if any(v not in {'PASS', 'FAIL', 'CANNOT VERIFY', 'NOT_ASSESSED'} for v in domains.values()):
            raise ValueError('Unknown source-readiness domain verdict')
        studies.append({'study_id': sid, 'manifest_ordinal': row['manifest_ordinal'],
                        'administrative_state': row['administrative_state'], 'source_access_allowed': not held,
                        'source_record_count': row['source_record_count_from_fixed_manifest'], 'domains': domains,
                        'source_validated': accepted, 'versions': history[sid],
                        'not_assessed_reason': (row.get('reason', row['administrative_state']) if held else
                            (validation['reason'] if validation and not accepted else ('NOT_SUBMITTED' if not frozen else None))),
                        'comparison_count': validation['comparisons'] if accepted else None,
                        'literal_witness_count': validation['literal_witnesses_checked'] if accepted else None})
    allowed = [r for r in studies if r['source_access_allowed']]
    counts = {}
    for domain in SOURCE_READINESS_DOMAINS:
        counter = Counter(r['domains'][domain] for r in allowed)
        numerator = counter['PASS'] + counter['FAIL']
        assessed = sum(r['source_validated'] for r in allowed)
        counts[domain] = {s: counter[s] for s in ('PASS', 'FAIL', 'CANNOT VERIFY', 'NOT_ASSESSED')}
        counts[domain].update({'assessed_and_source_validated': assessed, 'allowed_studies': 78, 'enrolled_studies': 80,
            'administratively_held': 2, 'submitted_studies': sum(bool(r['versions']) for r in allowed),
            'source_validated_coverage_allowed': {'numerator': assessed, 'denominator': 78},
            'source_validated_coverage_enrolled': {'numerator': assessed, 'denominator': 80},
            'resolved_coverage_allowed': {'numerator': numerator, 'denominator': 78},
            'resolved_coverage_enrolled': {'numerator': numerator, 'denominator': 80}})
        if sum(counter.values()) != 78: raise ValueError('Incomplete study-domain denominator')
    return {'schema': 'clockbase-rna80-source-readiness-table-v1', 'split': 'development_exposed',
            'enrolled_studies': 80, 'allowed_studies': 78, 'administrative_holds': 2,
            'enrolled_source_records': enrollment['enrolled_source_records'],
            'source_access_allowed_records': enrollment['source_access_allowed_records'],
            'administrative_hold_source_records': enrollment['administrative_hold_source_records'],
            'enrolled_domain_slots': 720, 'allowed_domain_slots': 702, 'held_domain_slots': 18,
            'domains': counts, 'studies': studies, 'new_model_calls': 0, 'outcome_values_read': False,
            'scope': 'Source-reference readiness only; quoted evidence and supplied interpretations are not scientific ground truth. No accuracy or repair estimate is computed.'}


def report_source_readiness(release_path, release_sha256, frozen_bindings, output):
    """Read pinned readiness versions; preserve rejected and superseded versions in counts."""
    from .agent_workflow_2026_09_10 import write_new
    context = source_readiness_context(release_path, release_sha256)
    versions = [(b, _readiness_frozen(context, b)) for b in frozen_bindings]
    table = source_readiness_table(context['enrollment'], versions)
    config = {'schema': 'clockbase-source-readiness-table-store-v1', 'release_sha256': release_sha256,
              'frozen_release': context['release'],
              'frozen_versions': frozen_bindings, 'input_paths': [str(release_path), *[b['path'] for b in frozen_bindings]]}
    def report(directory):
        write_new(directory / 'readiness.json', table)
        return {'table': str(directory / 'readiness.json'), 'enrolled_studies': 80,
                'source_validated_studies': sum(r['source_validated'] for r in table['studies'])}
    return _readiness_store(output, config, report)


def reference_archive_context(archive):
    """Dispatch only the explicit prospective schema to its release validator."""
    config = json.loads((Path(archive) / 'config.json').read_text())
    if config.get('schema') == 'clockbase-prospective-comparison-v1':
        from .heldout_2026_09_17 import validate_archive_context
        context = validate_archive_context(Path(archive))
        if context.get('split') not in {'development', 'held_out'} or not context.get('enrollment'):
            raise ValueError('Prospective archive validator returned no split/enrollment')
        return context
    if config.get('development_only') is not True:
        raise ValueError('This reference workflow requires a development or validated prospective archive')
    return {'split': 'development'}


def reference_packet_context(packet):
    if packet.get('schema') == 'clockbase-reference-design-packet-v2':
        from .heldout_2026_09_17 import validate_reference_packet_context
        validate_reference_packet_context(packet)
        if packet.get('split') not in {'development', 'held_out'}:
            raise ValueError('Unknown prospective source-reference split')
    elif packet.get('schema') != 'clockbase-reference-design-packet-v1' or packet.get('split') != 'development':
        raise ValueError('Source-only reference packet required')
    if packet.get('phase') != 'SOURCE_ONLY': raise ValueError('Source-only reference packet required')
    return packet['split']


def validate_attempted_universe(records, enrollment):
    """Bind denominators to enrollment, not to the surviving reference forms."""
    if (not isinstance(enrollment, dict) or enrollment.get('schema') != 'clockbase-comparison-enrollment-v1'
            or enrollment.get('split') not in {'development', 'held_out'}):
        raise ValueError('A declared comparison enrollment is required')
    methods, studies = enrollment.get('methods'), enrollment.get('studies')
    if (not isinstance(methods, list) or not methods or len(set(methods)) != len(methods)
            or any(not isinstance(method, str) or not method.strip() for method in methods)
            or not isinstance(studies, list) or not studies):
        raise ValueError('Enrollment requires unique methods and a nonempty study roster')
    families = {}
    for row in studies:
        if (set(row) != {'study_id', 'family_id'} or any(not isinstance(value, str) or not value.strip()
                for value in row.values()) or row['study_id'] in families):
            raise ValueError('Enrollment requires unique study identities with their source families')
        families[row['study_id']] = row['family_id']
    expected = {(study, method) for study in families for method in methods}
    observed = [(row.get('study_id'), row.get('method')) for row in records]
    if len(observed) != len(set(observed)) or set(observed) != expected:
        raise ValueError('Records differ from the frozen attempted study-by-method roster; missing studies or methods cannot be dropped')
    if any(row.get('family_id') != families[row['study_id']] for row in records):
        raise ValueError('Record family differs from the enrolled study family')
    return {'studies': len(studies), 'families': len(set(families.values())),
            'methods': len(methods), 'expected_records': len(expected)}


def reference_stage(archive, name, filename):
    """Read one pinned archive artifact without executing a method."""
    from .agent_workflow_2026_09_10 import sha
    from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
    archive = Path(archive).resolve(); folder = (archive / name).resolve()
    if not folder.is_relative_to(archive): raise PipelineIntegrityError('Reference stage escapes its archive')
    if sha(folder / 'complete.json') != (folder / 'complete.sha256').read_text().strip():
        raise PipelineIntegrityError('Reference source stage receipt changed')
    receipt = json.loads((folder / 'complete.json').read_text())
    for relative, digest in receipt['outputs'].items():
        path = (folder / relative).resolve()
        if not path.is_relative_to(folder) or sha(path) != digest:
            raise PipelineIntegrityError('Reference source artifact changed: ' + relative)
    matches = [folder / relative for relative in receipt['outputs'] if Path(relative).name == filename]
    if len(matches) != 1: raise ValueError('Missing or ambiguous reference source artifact: ' + filename)
    return json.loads(matches[0].read_text()), receipt['result']


def reference_cases(archive, _seen=None):
    from .agent_workflow_2026_09_10 import sha
    archive = Path(archive).resolve()
    seen = set() if _seen is None else set(_seen)
    if archive in seen: raise ValueError('Cyclic reference archive ancestry')
    seen.add(archive)
    config = json.loads((archive / 'config.json').read_text())
    context = reference_archive_context(archive)
    if (archive / '00_inputs' / 'complete.json').is_file():
        cases, receipt = reference_stage(archive, '00_inputs', 'prepared-cases.json')
        if not cases or cases != receipt['cases'] or len({case['study_id'] for case in cases}) != len(cases):
            raise ValueError('Reference cases differ from the pinned input stage')
    elif isinstance(config.get('archive'), str) and isinstance(config.get('archive_sha256'), dict):
        parent = Path(config['archive']).resolve()
        identities = {str(path): sha(path) for path in parent.rglob('*') if path.is_file()}
        if not identities or identities != config['archive_sha256']:
            raise ValueError('Replay reference source archive changed')
        cases = reference_cases(parent, seen)
    else:
        raise ValueError('Reference archive has neither prepared inputs nor a pinned replay source')
    enrollment = enrollment_sha = None
    if (archive / '00_enrollment' / 'complete.json').is_file():
        enrollment, enrolled = reference_stage(archive, '00_enrollment', 'enrollment.json')
        enrollment_sha = enrolled['sha256']
        if (enrollment != enrolled['enrollment'] or enrollment.get('split') != context['split']
                or enrollment.get('schema') != 'clockbase-comparison-enrollment-v1'
                or enrollment.get('methods') != list(COMPARISON_METHODS)):
            raise ValueError('Attempted enrollment differs from its receipt')
        if context.get('enrollment') is not None and enrollment != context['enrollment']:
            raise ValueError('Prepared archive differs from the validated prospective enrollment')
        expected = {row['study_id']: row['family_id'] for row in enrollment['studies']}
        if {case['study_id']: case['family_id'] for case in cases} != expected:
            raise ValueError('Prepared or unavailable cases differ from the frozen attempted enrollment')
    elif context.get('enrollment') is not None:
        raise ValueError('Prospective reference archive lacks its sealed enrollment')
    for case in cases:
        unavailable = case.get('preparation_status')
        if unavailable in {'UNAVAILABLE_SHARED_PREREQUISITE', 'UNAVAILABLE_PROSPECTIVE_PREPARATION'}:
            attempt = case.get('attempt', {})
            prospective = unavailable == 'UNAVAILABLE_PROSPECTIVE_PREPARATION'
            expected_schema = 'clockbase-prospective-unavailable-attempt-v1' if prospective else 'clockbase-unavailable-attempt-v1'
            if (enrollment is None or prospective != (context.get('enrollment') is not None)
                    or attempt.get('schema') != expected_schema
                    or attempt.get('enrollment_sha256') != enrollment_sha
                    or any(key in case for key in ('packet', 'primary', 'clocks'))):
                raise ValueError('Unavailable case lacks its enrolled attempt or invents prepared inputs')
            if prospective and (not isinstance(attempt.get('phase'), str) or not attempt['phase'].strip()):
                raise ValueError('Prospective preparation failure requires a pinned failure phase')
            for path, digest in case['source_sha256'].items():
                if sha(path) != digest: raise ValueError('Unavailable-attempt source changed: ' + path)
            if (case['source_sha256'].get(attempt['upstream_config']) != attempt['upstream_config_sha256']
                    or case['source_sha256'].get(attempt['failure_receipt']) != attempt['failure_sha256']
                    or not prospective and json.loads(Path(attempt['upstream_config']).read_text()).get('study') != case['study_id']
                    or json.loads(Path(attempt['failure_receipt']).read_text()) != attempt['recorded_failure']):
                raise ValueError('Unavailable source attempt identity differs from its evidence')
            if prospective:
                failure, returned = reference_stage(archive, attempt['phase'], 'terminal-attempt.json')
                if (Path(attempt['upstream_config']).resolve() != archive / 'config.json'
                        or not attempt['phase'].startswith(case['study_id'] + '/')
                        or failure != attempt['recorded_failure'] or failure.get('status') != 'CANNOT VERIFY'
                        or failure.get('phase') != attempt['phase'] or failure.get('automatic_retry') is not False
                        or Path(returned['terminal_receipt']).resolve() != Path(attempt['failure_receipt']).resolve()
                        or returned.get('terminal_sha256') != attempt['failure_sha256']
                        or {key: value for key, value in returned.items() if key not in {'terminal_receipt', 'terminal_sha256'}} != failure):
                    raise ValueError('Prospective preparation failure differs from its sealed terminal stage')
            continue
        if not case.get('family_id') or case['study_id'] != case['packet']['study_id']:
            raise ValueError('Reference study/family identity is missing or inconsistent')
        for path, digest in case['source_sha256'].items():
            if sha(path) != digest: raise ValueError('Reference empirical source changed: ' + path)
    return cases


def reference_source_view(packet):
    """Deliver source metadata/documents, excluding method outputs and summaries."""
    fields = {'study_id', 'target', 'comparison_scope', 'samples', 'study_documents', 'source_url', 'source_sha256',
              'context_sha256', 'source_context_withheld', 'withheld', 'sample_source', 'metadata_source'}
    view = {key: value for key, value in packet.items() if key in fields}
    allowed = {'sample_id', 'title', 'source_name_ch1', 'organism_ch1', 'molecule_ch1', 'platform_id',
               'treatment_protocol_ch1', 'growth_protocol_ch1', 'extract_protocol_ch1',
               'library_strategy', 'library_source', 'library_selection'}
    view['samples'] = [{key: value for key, value in sample.items() if key in allowed
                       or key.startswith('characteristics_ch1')
                       or re.fullmatch(r'(?:(?:treatment|growth|extract)_protocol_ch1|library_(?:strategy|source|selection))\.\d+', key)}
                      for sample in packet['samples']]
    ids = [sample['sample_id'] for sample in view['samples']]
    if not ids or len(set(ids)) != len(ids): raise ValueError('Reference source sample identities are missing or duplicated')
    return view


def export_reference_designs(archive, output):
    """Prepare uncompleted source-only forms; never manufacture reference labels."""
    from .agent_workflow_2026_09_10 import sha, write_new
    from src.pipeline_runtime_2026_09_10 import StageStore
    archive, output = Path(archive).resolve(), Path(output).resolve()
    if output == archive or output.is_relative_to(archive) or archive.is_relative_to(output):
        raise ValueError('Reference output must be separate from its source archive')
    if reference_archive_context(archive).get('enrollment') is not None:
        raise ValueError('Prospective references require the controller source-first packet export')
    cases = reference_cases(archive)
    config = {'version': REFERENCE_VERSION, 'implementation_sha256': sha(__file__), 'archive': str(archive),
              'archive_sha256': {str(path.relative_to(archive)): sha(path) for path in archive.rglob('*') if path.is_file()},
              'development_only': True, 'generation_calls': 0, 'reference_labels_completed': False}
    store = StageStore(output, config)

    def prepare(directory):
        packets, unavailable = [], []
        for case in cases:
            if case.get('preparation_status') == 'UNAVAILABLE_SHARED_PREREQUISITE':
                unavailable.append({'study_id': case['study_id'], 'family_id': case['family_id'],
                    'reason': 'No prepared source/numerical packet; freeze an unavailable reference separately',
                    'attempt': case['attempt']})
                continue
            folder = directory / case['study_id']; folder.mkdir()
            packet = {'schema': 'clockbase-reference-design-packet-v1', 'phase': 'SOURCE_ONLY',
                      'study_id': case['study_id'], 'family_id': case['family_id'], 'split': 'development',
                      'source_packet': reference_source_view(case['packet'])}
            source_path = folder / 'source-design.json'; write_new(source_path, packet)
            form = {'schema': 'clockbase-reference-design-judgment-v1', 'completed': False,
                    'packet_sha256': sha(source_path), 'study_id': case['study_id'], 'family_id': case['family_id'],
                    'reviewer': {'name': '', 'role': 'independent_source_adjudication', 'independent_from_methods': None},
                    'domains': {name: {'status': 'CANNOT VERIFY', 'material': None, 'rationale': '', 'evidence': []}
                                for name in REFERENCE_DESIGN_DOMAINS}, 'comparisons': [], 'additional_sources_needed': []}
            write_new(folder / 'design-form.json', form)
            content = '<!doctype html><meta charset="utf-8"><title>Source design reference</title><style>body{font:16px system-ui;max-width:1100px;margin:40px auto}pre{white-space:pre-wrap}td,th{padding:8px;text-align:left;vertical-align:top;border-bottom:1px solid #ddd}</style>'
            content += '<h1>' + html.escape(case['study_id']) + ': source-only reference review</h1><p>' + html.escape(packet['source_packet']['target']) + '</p>'
            content += '<p>Complete the design form from source evidence before viewing numerical outputs. Prior numerical outcomes, method identities and decisions are omitted. This is a development case, not a held-out validation observation.</p>'
            for sample in packet['source_packet']['samples']:
                content += '<h2>' + html.escape(sample['sample_id']) + '</h2><table>' + ''.join(
                    '<tr><th>' + html.escape(key) + '</th><td>' + html.escape(str(value)) + '</td></tr>'
                    for key, value in sample.items() if key != 'sample_id') + '</table>'
            for document in packet['source_packet'].get('study_documents', []):
                content += '<h2>Source document</h2><pre>' + html.escape(json.dumps(document, indent=2)) + '</pre>'
            (folder / 'source-design.html').write_text(content)
            packets.append({'study_id': case['study_id'], 'family_id': case['family_id'], 'source_packet': str(source_path),
                            'packet_sha256': sha(source_path), 'form': str(folder / 'design-form.json')})
        result = {'phase': 'SOURCE_ONLY', 'packets': packets, 'reference_labels_completed': False,
                  'generation_calls': 0, 'incremental_api_cost_usd': 0, 'held_out_release': 'CLOSED'}
        if unavailable:
            result['unavailable_cases'] = unavailable
            write_new(directory / 'unavailable-cases.json', unavailable)
        write_new(directory / 'reference-packets.json', result)
        return result
    return store.run('00_source_packets', prepare)


def freeze_reference_design(packet_path, judgment_path, output):
    """Freeze a submitted assessor judgment and validate its source quotations.

    The assessor declaration records provenance, not authenticated identity or
    proof that the scientific interpretation is correct.
    """
    from .agent_workflow_2026_09_10 import sha, write_new, _check_evidence, comparison_identity
    from src.pipeline_runtime_2026_09_10 import StageStore
    packet_path, judgment_path, output = Path(packet_path).resolve(), Path(judgment_path).resolve(), Path(output).resolve()
    if (output.is_relative_to(packet_path.parent) or packet_path.is_relative_to(output)
            or judgment_path.is_relative_to(output)):
        raise ValueError('Reference freeze output must be separate from its input packet and judgment')
    packet, judgment = json.loads(packet_path.read_text()), json.loads(judgment_path.read_text())
    reference_packet_context(packet)
    if judgment.get('schema') != 'clockbase-reference-design-judgment-v1' or judgment.get('completed') is not True:
        raise ValueError('An uncompleted reference form is not a judgment')
    if judgment.get('packet_sha256') != sha(packet_path) or any(judgment.get(key) != packet[key] for key in ('study_id', 'family_id')):
        raise ValueError('Reference judgment names a different source packet or study')
    reviewer = judgment.get('reviewer', {})
    if (not isinstance(reviewer.get('name'), str) or not reviewer['name'].strip()
            or reviewer.get('role') != 'independent_source_adjudication'
            or reviewer.get('independent_from_methods') is not True):
        raise ValueError('Named independent source adjudicator required; method agreement is not a reference')
    if set(judgment.get('domains', {})) != set(REFERENCE_DESIGN_DOMAINS):
        raise ValueError('All reference design domains must be assessed')
    source = packet['source_packet']
    def check_domains(domains):
        if not isinstance(domains, dict) or set(domains) != set(REFERENCE_DESIGN_DOMAINS):
            raise ValueError('All comparison source-design domains must be assessed')
        for name, domain in domains.items():
            if (domain.get('status') not in {'PASS', 'FAIL', 'CANNOT VERIFY'}
                    or not isinstance(domain.get('rationale'), str) or not domain['rationale'].strip()):
                raise ValueError('Missing reference status or rationale: ' + name)
            evidence, _ = _check_evidence(domain.get('evidence', []), source)
            if domain['status'] != 'CANNOT VERIFY' and not evidence:
                raise ValueError('Resolved source judgment requires checked evidence: ' + name)
            if ('material' in domain and domain['material'] is not None
                    and (type(domain['material']) is not bool or
                         domain['status'] != 'FAIL' and domain['material'])):
                raise ValueError('Source materiality must be an explicit failure judgment')
    check_domains(judgment['domains'])
    comparisons = judgment.get('comparisons')
    if not isinstance(comparisons, list): raise ValueError('Reference comparisons must be a list')
    if judgment['domains']['comparison_definition']['status'] == 'PASS' and not comparisons:
        raise ValueError('Resolved comparison definition requires explicit source groups')
    source_ids = {sample['sample_id'] for sample in source['samples']}
    seen, assigned = set(), set()
    for row in comparisons:
        if not isinstance(row.get('target'), str) or not row['target'].strip():
            raise ValueError('Reference comparison requires its estimand/target')
        groups = [row.get('case_ids'), row.get('control_ids')]
        if any(not isinstance(group, list) or not group or any(not isinstance(sid, str) for sid in group) for group in groups):
            raise ValueError('Reference comparison groups must be nonempty source-ID lists')
        members = groups[0] + groups[1]
        if len(set(members)) != len(members) or set(members) - source_ids:
            raise ValueError('Reference comparison duplicates or invents sample identities')
        cid = comparison_identity(packet['study_id'], *groups, row['target'])
        if cid in seen: raise ValueError('Duplicate reference comparison')
        evidence, supported = _check_evidence(row.get('evidence', []), source)
        if not evidence or set(members) - supported:
            raise ValueError('Every reference group sample requires source-linked evidence')
        if 'domains' in row:
            check_domains(row['domains'])
        seen.add(cid); assigned.update(members)
    config = {'version': REFERENCE_VERSION, 'implementation_sha256': sha(__file__),
              'packet_path': str(packet_path), 'packet_sha256': sha(packet_path),
              'judgment_path': str(judgment_path), 'judgment_sha256': sha(judgment_path)}
    store = StageStore(output, config)

    def freeze(directory):
        value = {'schema': ('clockbase-frozen-reference-design-v2' if packet['schema'].endswith('-v2')
                            else 'clockbase-frozen-reference-design-v1'), 'phase': 'DESIGN_FROZEN',
                 'study_id': packet['study_id'], 'family_id': packet['family_id'], 'split': packet['split'],
                 'packet': packet, 'judgment': judgment, 'source_identities': config,
                 'unassigned_ids': sorted(source_ids - assigned),
                 'scientific_reference_complete': False, 'assessor_identity_authenticated': False}
        destination = directory / 'frozen-design.json'; write_new(destination, value)
        (directory / 'frozen-design.sha256').write_text(sha(destination) + '\n')
        return {'phase': 'DESIGN_FROZEN', 'frozen_design': str(destination), 'sha256': sha(destination),
                'scientific_reference_complete': False}
    return store.run('00_design_freeze', freeze)


def export_reference_outcomes(archive, frozen_path, output):
    """Release anonymized outcomes only after a source judgment is frozen."""
    from .agent_workflow_2026_09_10 import sha, write_new
    from src.pipeline_runtime_2026_09_10 import StageStore
    archive, frozen_path, output = Path(archive).resolve(), Path(frozen_path).resolve(), Path(output).resolve()
    if output == archive or output.is_relative_to(archive) or archive.is_relative_to(output):
        raise ValueError('Outcome reference output must be separate from its source archive')
    if output.is_relative_to(frozen_path.parent.parent) or frozen_path.is_relative_to(output):
        raise ValueError('Numerical packet output must be separate from the frozen design stage')
    frozen = read_reference_design(frozen_path)
    context = reference_archive_context(archive)
    if frozen['split'] != context['split']: raise ValueError('Reference design/archive split differs')
    _require_registered_prospective_design(archive, frozen_path, frozen, context)
    cases = reference_cases(archive)
    matches = [case for case in cases if case['study_id'] == frozen['study_id']]
    if len(matches) != 1 or matches[0]['family_id'] != frozen['family_id']:
        raise ValueError('Frozen reference does not identify one study in this archive')
    case = matches[0]
    if case.get('preparation_status') in {'UNAVAILABLE_SHARED_PREREQUISITE', 'UNAVAILABLE_PROSPECTIVE_PREPARATION'}:
        raise ValueError('Unavailable prerequisites cannot supply a numerical reference packet')
    if reference_source_view(case['packet']) != frozen['packet']['source_packet']:
        raise ValueError('Outcome archive uses a different source design packet')
    config = {'version': REFERENCE_VERSION, 'implementation_sha256': sha(__file__), 'archive': str(archive),
              'frozen_design_path': str(frozen_path), 'frozen_design_sha256': sha(frozen_path),
              'archive_sha256': {str(path.relative_to(archive)): sha(path) for path in archive.rglob('*') if path.is_file()}}
    store = StageStore(output, config)

    def prepare(directory):
        outcomes, mapping, _ = reference_output_inventory(archive, case)
        reviewer = directory / 'reviewer'; reviewer.mkdir()
        packet = {'schema': ('clockbase-reference-outcome-packet-v3' if context.get('enrollment') is not None
                             else 'clockbase-reference-outcome-packet-v2'), 'phase': 'NUMERICAL_REVIEW',
                  'study_id': case['study_id'], 'family_id': case['family_id'], 'split': context['split'],
                  'frozen_design_sha256': sha(frozen_path), 'outcomes': [outcomes[key] for key in sorted(outcomes)],
                  'prediction_source': {'path': case['clocks']['path'], 'sha256': sha(case['clocks']['path'])},
                  'clock_declarations': [{key: row[key] for key in ('clock', 'model_sha256', 'required_features') if key in row}
                                         for row in case['clocks'].get('checks', [])],
                  'measurement_sources': [row['measurement_sources'] for row in case['clocks'].get('checks', []) if row.get('measurement_sources')],
                  'measurement_evidence_gaps': [row['clock'] + ': molecular input/model source files are not included in this archived clock declaration'
                                               for row in case['clocks'].get('checks', []) if not row.get('measurement_sources')],
                  'withheld': ['method identities', 'audit findings', 'acceptance decisions', 'reference answers'],
                  'reference_labels_completed': False}
        write_new(reviewer / 'outcomes.json', packet)
        form = {'schema': 'clockbase-reference-outcome-judgment-v2', 'completed': False,
                'packet_sha256': sha(reviewer / 'outcomes.json'), 'frozen_design_sha256': sha(frozen_path),
                'reviewer': frozen['judgment']['reviewer'], 'supporting_artifacts': [],
                'instructions': 'Assess the declared output and its inferential scope. An explicit abstention or absence of a paper-wide FDR claim is not by itself an error. Explain materiality for each FAIL. CANNOT VERIFY is retained.',
                'assessments': [{'output_id': key, 'domains': {name: {
                    'status': 'CANNOT VERIFY', 'material': None, 'rationale': '', 'evidence': []}
                    for name in REFERENCE_OUTCOME_DOMAINS}} for key in sorted(outcomes)]}
        write_new(reviewer / 'outcome-form.json', form)
        write_new(directory / 'controller-output-mapping.json', mapping)
        result = {'phase': 'NUMERICAL_REVIEW', 'study_id': case['study_id'], 'unique_outputs': len(outcomes),
                  'mapped_output_occurrences': len(mapping), 'reviewer_packet': str(reviewer / 'outcomes.json'),
                  'reference_labels_completed': False, 'generation_calls': 0, 'held_out_release': 'CLOSED'}
        if context.get('enrollment') is not None:
            result.update(split=context['split'], held_out_release='OPEN' if context['split'] == 'held_out' else 'CLOSED')
        write_new(directory / 'reference-outcomes.json', result)
        return result
    return store.run('00_numerical_packets', prepare)


def read_reference_design(frozen_path):
    from .agent_workflow_2026_09_10 import sha
    frozen_path = Path(frozen_path).resolve()
    frozen = json.loads(frozen_path.read_text())
    checked_frozen, freeze_receipt = reference_stage(
        frozen_path.parent.parent.parent, frozen_path.parent.parent.name, 'frozen-design.json')
    if (frozen.get('schema') not in {'clockbase-frozen-reference-design-v1', 'clockbase-frozen-reference-design-v2'}
            or frozen.get('phase') != 'DESIGN_FROZEN'
            or sha(frozen_path) != frozen_path.with_suffix('.sha256').read_text().strip()
            or frozen != checked_frozen or freeze_receipt['sha256'] != sha(frozen_path)):
        raise ValueError('Missing or changed frozen source-design judgment')
    for name in ('packet', 'judgment'):
        source = frozen['source_identities']
        if sha(source[name + '_path']) != source[name + '_sha256']:
            raise ValueError('Source packet or submitted reference judgment changed after freeze')
        if json.loads(Path(source[name + '_path']).read_text()) != frozen[name]:
            raise ValueError('Frozen design differs from the submitted reference inputs')
    split = reference_packet_context(frozen['packet'])
    expected = ('clockbase-frozen-reference-design-v2' if frozen['packet']['schema'].endswith('-v2')
                else 'clockbase-frozen-reference-design-v1')
    if frozen['split'] != split or frozen['schema'] != expected:
        raise ValueError('Frozen reference design schema/split differs from its source packet')
    return frozen


def _require_registered_prospective_design(archive, frozen_path, frozen, context):
    """A source judgment cannot be selected or replaced after method execution."""
    from .agent_workflow_2026_09_10 import sha
    if context.get('enrollment') is None: return
    state, returned = reference_stage(archive, frozen['study_id'] + '/reference_state', 'reference-state.json')
    references = [row for row in state.get('references', []) if row.get('phase') == 'SOURCE_DESIGN']
    if (frozen.get('schema') != 'clockbase-frozen-reference-design-v2' or state != returned
            or state.get('study_id') != frozen['study_id'] or len(references) != 1
            or Path(frozen['packet'].get('evaluation_context', {}).get('archive', '')).resolve() != Path(archive).resolve()
            or Path(references[0]['path']).resolve() != Path(frozen_path).resolve()
            or references[0]['sha256'] != sha(frozen_path)):
        raise ValueError('Prospective numerical review requires the pre-execution sealed source judgment')


def reference_output_inventory(archive, case):
    from .agent_workflow_2026_09_10 import sha
    outcomes, mapping, methods = {}, [], {}
    fields = {'clock', 'case_ids', 'control_ids', 'statistical_test', 'p_value', 'mean_difference', 'effect_size',
              'effect_size_definition', 'p_value_by_attempted_family', 'q_value_run_family', 'status',
              'covariate', 'covariates', 'claimed_covariates', 'adjustment', 'estimand'}

    def retain(row, origin, multiplicity=None):
        visible = {key: row[key] for key in fields if key in row}
        visible['target'] = row.get('target', case['packet']['target'])
        visible['prediction_source_sha256'] = sha(case['clocks']['path'])
        if multiplicity is not None:
            visible['multiplicity'] = {key: multiplicity[key] for key in ('method', 'family_denominator', 'attempted_endpoints') if key in multiplicity}
        token = hashlib.sha256(json.dumps(visible, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        output_id = 'output-' + token[:24]
        item = {'output_id': output_id, **visible}
        if output_id in outcomes and outcomes[output_id] != item: raise ValueError('Reference output identity collision')
        outcomes[output_id] = item
        mapping.append({'output_id': output_id, **origin})
    for index, row in enumerate(case['primary']['rows']): retain(row, {'method': 'producer', 'collection': 'original', 'index': index})
    for method in COMPARISON_METHODS:
        result, receipt = reference_stage(archive, case['study_id'] + '/' + method + '/decision', 'method-result.json')
        if result != receipt: raise ValueError('Method output differs from its source-stage result')
        methods[method] = result
        for collection in ('final_outputs', 'candidates'):
            for index, row in enumerate(result.get(collection, [])):
                retain(row, {'method': method, 'collection': collection, 'index': index}, result.get('multiplicity'))
    return outcomes, mapping, methods


def reference_set_status(values):
    """A missing member cannot be silently treated as a correct result."""
    values = list(values)
    return ('FAIL' if 'FAIL' in values else 'PASS' if values and all(v == 'PASS' for v in values)
            else 'CANNOT VERIFY')


def reference_output_design(frozen, output):
    """Carry frozen biological-design limits into a computed output's label.

    Comparison-specific judgments take precedence when provided before numerical
    release. Legacy judgments use their global domains; an unmatched contrast
    cannot acquire design validity from numerical arithmetic alone.
    """
    from .agent_workflow_2026_09_10 import comparison_identity
    computed = any(output.get(key) not in (None, '') for key in ('p_value', 'mean_difference', 'effect_size'))
    if not computed:
        return {'status': 'PASS', 'scope': 'No computed inferential output; no biological-design approval implied'}
    cid = comparison_identity(frozen['study_id'], output.get('case_ids', []),
                              output.get('control_ids', []), output['target'])
    matches = [row for row in frozen['judgment']['comparisons']
               if comparison_identity(frozen['study_id'], row['case_ids'], row['control_ids'], row['target']) == cid]
    if len(matches) != 1:
        return {'status': 'CANNOT VERIFY', 'scope': 'No unique source-frozen comparison matches this output'}
    domains = matches[0].get('domains', frozen['judgment']['domains'])
    effective = [('FAIL' if row.get('material') is True else 'PASS' if row.get('material') is False else 'CANNOT VERIFY')
                 if row['status'] == 'FAIL' else row['status'] for row in domains.values()]
    return {'status': reference_set_status(effective),
            'domain_statuses': {name: row['status'] for name, row in domains.items()},
            'material_failure_domains': [name for name, row in domains.items()
                                        if row['status'] == 'FAIL' and row.get('material') is True],
            'scope': 'Comparison-specific source domains' if 'domains' in matches[0] else 'Global source domains from legacy judgment'}


def _reference_pointer(document, pointer):
    if not isinstance(pointer, str) or (pointer and not pointer.startswith('/')):
        raise ValueError('Reference evidence requires a JSON pointer')
    current = document
    for token in pointer.split('/')[1:]:
        if re.search(r'~(?![01])', token): raise ValueError('Invalid reference JSON pointer escape')
        token = token.replace('~1', '/').replace('~0', '~')
        if isinstance(current, list):
            if not re.fullmatch(r'0|[1-9][0-9]*', token): raise ValueError('Invalid reference array index')
            current = current[int(token)]
        else:
            current = current[token]
    return current


def validate_reference_outcomes(numerical_archive, judgment_path):
    """Check links and assessor declarations; these checks do not establish scientific truth."""
    from .agent_workflow_2026_09_10 import sha, _check_evidence, comparison_identity
    numerical_archive, judgment_path = Path(numerical_archive).resolve(), Path(judgment_path).resolve()
    config = json.loads((numerical_archive / 'config.json').read_text())
    archive = Path(config['archive']).resolve()
    context = reference_archive_context(archive)
    if {str(p.relative_to(archive)): sha(p) for p in archive.rglob('*') if p.is_file()} != config['archive_sha256']:
        raise ValueError('Reference method archive changed after numerical release')
    frozen_path = Path(config['frozen_design_path'])
    if sha(frozen_path) != config['frozen_design_sha256']: raise ValueError('Reference design changed after numerical release')
    frozen = read_reference_design(frozen_path)
    _require_registered_prospective_design(archive, frozen_path, frozen, context)
    packet, release = reference_stage(numerical_archive, '00_numerical_packets', 'outcomes.json')
    mapping, _ = reference_stage(numerical_archive, '00_numerical_packets', 'controller-output-mapping.json')
    expected_schema = ('clockbase-reference-outcome-packet-v3' if context.get('enrollment') is not None
                       else 'clockbase-reference-outcome-packet-v2')
    if (packet.get('schema') != expected_schema or packet.get('phase') != 'NUMERICAL_REVIEW'
            or packet.get('split') != context['split'] or frozen['split'] != context['split']
            or packet['frozen_design_sha256'] != sha(frozen_path)):
        raise ValueError('Version-bound numerical reference packet required')
    cases = [case for case in reference_cases(archive) if case['study_id'] == packet['study_id']]
    if len(cases) != 1 or cases[0]['family_id'] != packet['family_id'] or packet['study_id'] != frozen['study_id']:
        raise ValueError('Numerical reference study/family mismatch')
    case = cases[0]
    if reference_source_view(case['packet']) != frozen['packet']['source_packet']:
        raise ValueError('Reference source design differs from the method source')
    outcomes, expected_mapping, methods = reference_output_inventory(archive, case)
    if mapping != expected_mapping or packet['outcomes'] != [outcomes[key] for key in sorted(outcomes)]:
        raise ValueError('Anonymous reference outputs do not match every archived output')
    judgment = json.loads(judgment_path.read_text())
    if judgment.get('schema') != 'clockbase-reference-outcome-judgment-v2' or judgment.get('completed') is not True:
        raise ValueError('An uncompleted numerical form is not a reference judgment')
    if (judgment.get('packet_sha256') != sha(release['reviewer_packet'])
            or judgment.get('frozen_design_sha256') != sha(frozen_path)):
        raise ValueError('Numerical judgment identifies a different packet or design')
    if judgment.get('reviewer') != frozen['judgment']['reviewer']:
        raise ValueError('Numerical reviewer differs from the frozen independent design reviewer')
    assessments = judgment.get('assessments')
    if (not isinstance(assessments, list) or len(assessments) != len(outcomes)
            or {item.get('output_id') for item in assessments} != set(outcomes)):
        raise ValueError('Every anonymous output needs exactly one assessment, including rejected outputs')
    supporting = judgment.get('supporting_artifacts')
    if not isinstance(supporting, list): raise ValueError('Supporting reference artifacts must be a list')
    artifacts, dependencies = {}, {}
    for item in supporting:
        aid = item.get('artifact_id')
        if not isinstance(aid, str) or not aid.strip() or aid in artifacts:
            raise ValueError('Supporting artifact identifiers must be unique')
        path = Path(item['path']).resolve()
        if item.get('role') not in {'independent_calculation', 'measurement_evidence', 'source_document'}:
            raise ValueError('Unsupported reference artifact role')
        if sha(path) != item.get('sha256'): raise ValueError('Supporting reference artifact changed')
        if not isinstance(item.get('inputs'), dict) or not item['inputs']:
            raise ValueError('Supporting reference artifact needs source file identities')
        if not item.get('output_ids') or set(item['output_ids']) - set(outcomes):
            raise ValueError('Supporting reference artifact must name its assessed outputs')
        for source, digest in item['inputs'].items():
            source = str(Path(source).resolve())
            if sha(source) != digest: raise ValueError('Supporting reference calculation input changed')
            dependencies[source] = digest
        if item['role'] == 'independent_calculation':
            if (item.get('independent_from_methods') is not True
                    or item['inputs'].get(case['clocks']['path']) != sha(case['clocks']['path'])
                    or path.is_relative_to(archive)):
                raise ValueError('Independent calculation must use the named predictions and be separate from method outputs')
        dependencies[str(path)] = sha(path)
        artifacts[aid] = (item, json.loads(path.read_text()))
    labels, design_constraints = {}, {}
    for assessment in assessments:
        output_id, domains = assessment['output_id'], assessment.get('domains', {})
        if set(domains) != set(REFERENCE_OUTCOME_DOMAINS): raise ValueError('All numerical reference domains must be assessed')
        output, statuses = outcomes[output_id], []
        computed = any(output.get(key) not in (None, '') for key in ('p_value', 'mean_difference', 'effect_size'))
        for name, domain in domains.items():
            status = domain.get('status')
            if (status not in {'PASS', 'FAIL', 'CANNOT VERIFY'} or not isinstance(domain.get('rationale'), str)
                    or not domain['rationale'].strip()): raise ValueError('Missing numerical reference status or rationale')
            if status == 'FAIL' and type(domain.get('material')) is not bool:
                raise ValueError('Reference failures require an explicit materiality judgment and rationale')
            if status != 'FAIL' and domain.get('material') not in (None, False):
                raise ValueError('Material failure conflicts with reference domain status')
            evidence = domain.get('evidence')
            if not isinstance(evidence, list) or (status != 'CANNOT VERIFY' and not evidence):
                raise ValueError('Resolved reference domain requires checked evidence')
            roles = set()
            for witness in evidence:
                if 'artifact_id' in witness:
                    if set(witness) != {'artifact_id', 'pointer', 'value'}: raise ValueError('Reference artifact witness schema mismatch')
                    item, document = artifacts[witness['artifact_id']]
                    if output_id not in item['output_ids']: raise ValueError('Reference evidence concerns another output')
                    observed = _reference_pointer(document, witness['pointer'])
                    if json.dumps(observed, sort_keys=True, allow_nan=False) != json.dumps(witness['value'], sort_keys=True, allow_nan=False):
                        raise ValueError('Reference artifact witness differs from its source')
                    roles.add(item['role'])
                else:
                    _check_evidence([witness], frozen['packet']['source_packet'])
            if computed and status != 'CANNOT VERIFY' and name == 'numerical_accuracy' and 'independent_calculation' not in roles:
                raise ValueError('Resolved numerical accuracy requires a separate source-linked calculation')
            if computed and status == 'PASS' and name == 'measurement_applicability' and 'measurement_evidence' not in roles:
                raise ValueError('Measurement PASS requires source-linked molecular/model evidence')
            if computed and status == 'PASS' and name == 'statistical_validity':
                cid = comparison_identity(case['study_id'], output.get('case_ids', []), output.get('control_ids', []), output['target'])
                eligible = {(comparison_identity(case['study_id'], row['case_ids'], row['control_ids'], row['target']), row['target'])
                            for row in frozen['judgment']['comparisons']}
                if (cid, output['target']) not in eligible: raise ValueError('Statistical PASS must identify a source-frozen comparison and target')
            statuses.append('PASS' if status == 'FAIL' and domain['material'] is False else status)
        design_constraints[output_id] = reference_output_design(frozen, output)
        labels[output_id] = reference_set_status(statuses + [design_constraints[output_id]['status']])
    identities = {'numerical_archive': str(numerical_archive),
                  'numerical_archive_sha256': {str(p.relative_to(numerical_archive)): sha(p) for p in numerical_archive.rglob('*') if p.is_file()},
                  'judgment_path': str(judgment_path), 'judgment_sha256': sha(judgment_path),
                  'supporting_source_sha256': dependencies}
    return {'case': case, 'methods': methods, 'mapping': mapping, 'labels': labels, 'judgment': judgment,
            'design_constraints': design_constraints,
            'frozen_design': frozen, 'archive': str(archive), 'source_identities': identities,
            'split': context['split'], 'prospective': context.get('enrollment') is not None}


def reference_method_facts(archive, case, methods, mapping):
    """Project archived method decisions without assuming reference availability.

    Private output-pair fields support both legacy and incomplete-reference
    projections. Existing frozen-reference records retain their original shape.
    """
    from .agent_workflow_2026_09_10 import sha
    archive = Path(archive)
    replay = 'archive' in json.loads((archive / 'config.json').read_text())
    if set(methods) != set(COMPARISON_METHODS):
        raise ValueError('Reference method facts require every enrolled method slot')
    lookup = {(row['method'], row['collection'], row['index']): row['output_id'] for row in mapping}
    if len(lookup) != len(mapping): raise ValueError('Duplicate reference output occurrence')
    original_rows = case.get('primary', {}).get('rows', [])
    originals = [lookup['producer', 'original', i] for i in range(len(original_rows))]
    facts = []
    for method, result in methods.items():
        if (result.get('method') != method or result.get('study_id', case['study_id']) != case['study_id']
                or result.get('family_id', case['family_id']) != case['family_id']):
            raise ValueError('Method result identity differs from the reference case')
        accepted = [i for i, row in enumerate(result.get('final_outputs', [])) if row.get('accepted_result') is True]
        if result.get('accepted_results') != len(accepted): raise ValueError('Accepted output count differs from archived method result')
        if bool(accepted) != (result.get('study_decision_status') == 'PASS'):
            raise ValueError('Archived acceptance differs from its method decision')
        corrections = [row for row in result.get('checks', []) if row.get('repair_kind') == 'CORRECTION_PROPOSAL']
        if type(result.get('repair_proposed')) is not bool or result['repair_proposed'] != bool(corrections):
            raise ValueError('Proposed repair count differs from correction checks')
        pairs = []
        for check in corrections:
            def matches(rows, comparison_id):
                return [i for i, row in enumerate(rows) if row.get('clock') == check['clock']
                        and (comparison_id is None or row.get('comparison_id') == comparison_id)]
            candidates = matches(result.get('candidates', []), check.get('comparison_id'))
            source = matches(original_rows, check.get('primary_analysis_state', {}).get('comparison_id'))
            pairs.append(None if len(candidates) != 1 or len(source) != 1 else
                         (lookup['producer', 'original', source[0]], lookup[method, 'candidates', candidates[0]]))
        durations = []
        for stage in ('review', 'decision'):
            folder = archive / case['study_id'] / method / stage
            if (folder / 'complete.json').is_file():
                if sha(folder / 'complete.json') != (folder / 'complete.sha256').read_text().strip():
                    raise ValueError('Method timing receipt changed')
                durations.append(json.loads((folder / 'complete.json').read_text())['elapsed_seconds'])
        slots = list((archive / 'api').glob(case['study_id'] + '-' + method + '-*'))
        resources = [json.loads((slot / 'receipt.json').read_text()) for slot in slots if (slot / 'receipt.json').is_file()]
        output_ids = {'original': originals,
                      'accepted': [lookup[method, 'final_outputs', i] for i in accepted],
                      'final': [lookup[method, 'final_outputs', i] for i in range(len(result.get('final_outputs', [])))],
                      'candidates': [lookup[method, 'candidates', i] for i in range(len(result.get('candidates', [])))]}
        facts.append({'study_id': case['study_id'], 'family_id': case['family_id'], 'method': method,
            'accepted': bool(accepted), 'decision_status': result['study_decision_status'],
            'flagged_material_error': result['flagged_material_error'], 'repair_proposed': bool(corrections),
            'runtime_seconds': None if replay else sum(durations), 'cost_usd': None if replay or slots else 0,
            'estimated_cost_usd': None if replay or len(resources) != len(slots) else sum(row['estimated_cost_usd'] for row in resources),
            'resource_scope': 'Incremental archived review/decision; excludes shared producer execution. Replays lack comparable total method time/cost.',
            'recorded_current_archive_stage_seconds': sum(durations), 'accepted_output_count': len(accepted),
            'correction_candidate_count': len(corrections),
            'reference_output_ids': {key: output_ids[key] for key in ('original', 'accepted')},
            '_all_output_ids': output_ids, '_correction_pairs': pairs})
    return facts


def _reference_fact_labels(fact, labels):
    repairs = []
    for pair in fact['_correction_pairs']:
        if pair is None:
            repairs.append('CANNOT VERIFY'); continue
        original, candidate = (labels[oid] for oid in pair)
        repairs.append('FAIL' if candidate == 'FAIL' else 'PASS' if candidate == 'PASS' and original == 'FAIL'
                       else 'CANNOT VERIFY')
    return {'reference_status': reference_set_status(labels[oid] for oid in fact['_all_output_ids']['original']),
            'final_reference_status': reference_set_status(labels[oid] for oid in fact['_all_output_ids']['accepted']),
            'repair_status': reference_set_status(repairs)}


def reference_study_records(checked):
    """Derive legacy records while preserving their archived serialized fields."""
    facts = reference_method_facts(checked['archive'], checked['case'], checked['methods'], checked['mapping'])
    records = [{**{key: value for key, value in fact.items() if not key.startswith('_')},
             **_reference_fact_labels(fact, checked['labels']),
             'reference_source': checked['frozen_design']['source_identities']['judgment_path'],
             'reference_adjudicator': checked['judgment']['reviewer']['name'],
             'reference_role': 'independent_source_adjudication',
             'final_reference_source': checked['source_identities']['judgment_path']} for fact in facts]
    if checked.get('prospective'):
        records = [{**_reference_record_resources(Path(checked['archive']), row, prospective=True),
                    'reference_output_ids': fact['_all_output_ids']} for row, fact in zip(records, facts)]
    return records


def _reference_record_resources(archive, row, *, prospective=False):
    """Keep incremental reviews and shared producer costs in separate accounts."""
    resources = comparison_resource_record(archive, row)
    # Replays retain unknown totals, while known usage and billing gaps survive.
    resources.update({key: row[key] for key in ('runtime_seconds', 'cost_usd', 'estimated_cost_usd', 'resource_scope')})
    if not prospective: return resources
    resources['resource_scope'] = 'Incremental archived method review/decision; fresh shared producer execution is reported separately'
    resources['api_evidence_scope'] = 'Incremental method review only; the prospective producer has its own shared API account'
    directory = archive / 'baseline-api' / row['study_id']
    slots = {path.parent for path in directory.rglob('reservation.json')}
    receipt_slots = {path.parent for path in directory.rglob('receipt.json')}
    if not receipt_slots <= slots: raise ValueError('Producer API receipt lacks its request reservation')
    receipts = [json.loads((path / 'receipt.json').read_text()) for path in sorted(receipt_slots)]
    for receipt in receipts:
        if (type(receipt.get('estimated_cost_usd')) not in (int, float)
                or not math.isfinite(receipt['estimated_cost_usd']) or receipt['estimated_cost_usd'] < 0
                or not isinstance(receipt.get('returned_model'), str)
                or any(type(receipt.get('usage', {}).get(key)) is not int or receipt['usage'][key] < 0
                       for key in ('prompt_tokens', 'completion_tokens'))):
            raise ValueError('Invalid shared producer API usage receipt')
    resources['shared_upstream_resources'] = {
        'study_id': row['study_id'], 'scope': 'Shared producer API account; count once per study, separately from incremental method costs',
        'runtime_seconds': None, 'generation_attempts': len(slots), 'responses_with_usage': len(receipts),
        'unknown_billing_attempts': len(slots - receipt_slots), 'billed_cost_usd': None if slots else 0,
        'estimated_cost_usd': None if slots != receipt_slots else sum(item['estimated_cost_usd'] for item in receipts),
        'estimated_known_cost_usd': sum(item['estimated_cost_usd'] for item in receipts),
        'prompt_tokens': sum(item['usage']['prompt_tokens'] for item in receipts),
        'completion_tokens': sum(item['usage']['completion_tokens'] for item in receipts),
        'returned_models': sorted({item['returned_model'] for item in receipts})}
    return resources


def freeze_reference_outcomes(numerical_archive, judgment_path, output):
    from .agent_workflow_2026_09_10 import sha, write_new
    from src.pipeline_runtime_2026_09_10 import StageStore
    output, numerical_archive, judgment_path = Path(output).resolve(), Path(numerical_archive).resolve(), Path(judgment_path).resolve()
    if output.is_relative_to(numerical_archive) or numerical_archive.is_relative_to(output) or judgment_path.is_relative_to(output):
        raise ValueError('Outcome freeze output must be separate from its input archive and judgment')
    checked = validate_reference_outcomes(numerical_archive, judgment_path)
    records = reference_study_records(checked)
    config = {'version': REFERENCE_VERSION, 'implementation_sha256': sha(__file__), **checked['source_identities']}
    store = StageStore(output, config)
    def freeze(directory):
        value = {'schema': ('clockbase-frozen-reference-outcomes-v2' if checked['prospective']
                            else 'clockbase-frozen-reference-outcomes-v1'), 'phase': 'OUTCOMES_FROZEN',
                 'source_identities': checked['source_identities'], 'judgment': checked['judgment'],
                 'labels': checked['labels'], 'records': records, 'split': checked['split'],
                 'source_design_constraints': checked['design_constraints'],
                 'reference_labels_completed': True, 'assessor_identity_authenticated': False,
                 'scope': 'Validated evidence links and assessor declarations; scientific interpretations remain the assessor responsibility'}
        path = directory / 'frozen-outcomes.json'; write_new(path, value)
        path.with_suffix('.sha256').write_text(sha(path) + '\n')
        return {'phase': 'OUTCOMES_FROZEN', 'frozen_outcomes': str(path), 'sha256': sha(path),
                'reference_labels_completed': True, 'reference_resolved_outputs': sum(s != 'CANNOT VERIFY' for s in checked['labels'].values())}
    return store.run('00_outcome_freeze', freeze)


def read_reference_records(frozen_path):
    from .agent_workflow_2026_09_10 import sha
    frozen_path = Path(frozen_path).resolve()
    frozen, receipt = reference_stage(frozen_path.parent.parent.parent, frozen_path.parent.parent.name, 'frozen-outcomes.json')
    if frozen.get('schema') == 'clockbase-reference-incomplete-v1':
        if (sha(frozen_path) != receipt['sha256']
                or sha(frozen_path) != frozen_path.with_suffix('.sha256').read_text().strip()
                or json.loads(frozen_path.read_text()) != frozen
                or frozen.get('reference_labels_completed') is not False):
            raise ValueError('Incomplete reference freeze changed')
        checked = _incomplete_reference_state(**frozen['inputs'])
        if frozen['state'] != checked or frozen.get('split') != checked['split'] or frozen.get('phase') != 'REFERENCE_INCOMPLETE':
            raise ValueError('Incomplete reference source, prior labels or method decisions changed')
        return checked['records']
    if frozen.get('schema') == 'clockbase-reference-unavailable-v1':
        identities = frozen['source_identities']; archive = Path(identities['archive'])
        if (frozen.get('phase') != 'REFERENCE_UNAVAILABLE' or frozen.get('split') != 'development'
                or frozen.get('reference_labels_completed') is not False
                or sha(frozen_path) != receipt['sha256']
                or sha(frozen_path) != frozen_path.with_suffix('.sha256').read_text().strip()
                or {str(p.relative_to(archive)): sha(p) for p in archive.rglob('*') if p.is_file()}
                   != identities['archive_sha256']):
            raise ValueError('Unavailable reference archive or receipt changed')
        records = unavailable_reference_records(archive, frozen['study_id'])
        if frozen['records'] != records:
            raise ValueError('Unavailable reference altered the archived method decisions')
        return records
    if (frozen.get('schema') not in {'clockbase-frozen-reference-outcomes-v1', 'clockbase-frozen-reference-outcomes-v2'}
            or frozen.get('phase') != 'OUTCOMES_FROZEN'
            or frozen.get('schema') == 'clockbase-frozen-reference-outcomes-v1' and frozen.get('split') != 'development'
            or sha(frozen_path) != receipt['sha256']
            or sha(frozen_path) != frozen_path.with_suffix('.sha256').read_text().strip()
            or json.loads(frozen_path.read_text()) != frozen): raise ValueError('Missing or changed frozen outcome judgment')
    identities = frozen['source_identities']
    checked = validate_reference_outcomes(identities['numerical_archive'], identities['judgment_path'])
    records = reference_study_records(checked)
    expected_schema = 'clockbase-frozen-reference-outcomes-v2' if checked['prospective'] else 'clockbase-frozen-reference-outcomes-v1'
    if (frozen['schema'] != expected_schema or frozen['split'] != checked['split']
            or checked['source_identities'] != identities or checked['judgment'] != frozen['judgment']
            or checked['labels'] != frozen['labels'] or records != frozen['records']
            or 'source_design_constraints' in frozen and frozen['source_design_constraints'] != checked['design_constraints']):
        raise ValueError('Reference evidence, judgment or derived study labels changed after freeze')
    return records


def unavailable_reference_records(archive, study_id):
    """Administrative missing-reference records, never scientific judgments."""
    archive = Path(archive).resolve()
    if reference_archive_context(archive).get('enrollment') is not None:
        raise ValueError('Legacy unavailable references cannot label prospective preparation failures')
    matches = [case for case in reference_cases(archive) if case['study_id'] == study_id]
    if len(matches) != 1 or matches[0].get('preparation_status') != 'UNAVAILABLE_SHARED_PREREQUISITE':
        raise ValueError('This unavailable-reference adapter requires a frozen shared-prerequisite failure; completed analyses need independent adjudication')
    case, records = matches[0], []
    for method in COMPARISON_METHODS:
        result, returned = reference_stage(archive, study_id + '/' + method + '/decision', 'method-result.json')
        if (result != returned or result.get('study_id') != study_id or result.get('family_id') != case['family_id']
                or result.get('method') != method or result.get('preparation_attempt') != case['attempt']
                or result.get('execution_status') != 'NOT_RUN_SHARED_PREREQUISITE'
                or result.get('study_decision_status') != 'CANNOT VERIFY' or result.get('accepted_results') != 0
                or any(result.get(key) != [] for key in ('final_outputs', 'candidates', 'checks'))
                or any(result.get(key) is not False for key in ('repair_proposed', 'flagged_material_error', 'raw_model_alarm'))
                or list((archive / 'api').glob(study_id + '-' + method + '-*'))):
            raise ValueError('Unavailable reference cannot erase executed method results, flags, repairs or API attempts')
        records.append({'study_id': study_id, 'family_id': case['family_id'], 'method': method,
            'reference_role': 'reference_unavailable', 'reference_adjudicator': None,
            'reference_availability': 'UNAVAILABLE_SHARED_PREREQUISITE',
            'reference_source': case['attempt']['failure_receipt'],
            'final_reference_source': case['attempt']['failure_receipt'],
            'reference_status': 'CANNOT VERIFY', 'final_reference_status': 'CANNOT VERIFY',
            'accepted': False, 'decision_status': result['study_decision_status'],
            'flagged_material_error': result['flagged_material_error'], 'repair_proposed': result['repair_proposed'],
            'repair_status': 'CANNOT VERIFY', 'runtime_seconds': None, 'cost_usd': 0, 'estimated_cost_usd': 0,
            'resource_scope': 'No incremental method execution or API call; shared upstream runtime/cost are unknown here',
            'upstream_runtime_seconds': None, 'upstream_cost_usd': None,
            'accepted_output_count': 0, 'correction_candidate_count': 0,
            'reference_output_ids': {'original': [], 'accepted': []}})
    return records


def freeze_reference_unavailable(archive, study_id, output):
    """Freeze unavailable assessment separately from completed adjudication."""
    from .agent_workflow_2026_09_10 import sha, write_new
    from src.pipeline_runtime_2026_09_10 import StageStore
    archive, output = Path(archive).resolve(), Path(output).resolve()
    if output == archive or output.is_relative_to(archive) or archive.is_relative_to(output):
        raise ValueError('Unavailable reference must be separate from its method archive')
    records = unavailable_reference_records(archive, study_id)
    identities = {'archive': str(archive), 'archive_sha256': {
        str(p.relative_to(archive)): sha(p) for p in archive.rglob('*') if p.is_file()}}
    store = StageStore(output, {'version': 'clockbase-reference-unavailable-v1',
        'implementation_sha256': sha(__file__), 'study_id': study_id, **identities})
    def freeze(directory):
        value = {'schema': 'clockbase-reference-unavailable-v1', 'phase': 'REFERENCE_UNAVAILABLE',
            'study_id': study_id, 'split': 'development', 'source_identities': identities, 'records': records,
            'reference_labels_completed': False, 'reference_adjudicator': None,
            'scope': 'Missing assessment caused by an unavailable shared prerequisite; no scientific judgment'}
        path = directory / 'frozen-outcomes.json'; write_new(path, value)
        path.with_suffix('.sha256').write_text(sha(path) + '\n')
        return {'phase': 'REFERENCE_UNAVAILABLE', 'frozen_outcomes': str(path), 'sha256': sha(path),
                'reference_labels_completed': False, 'reference_resolved_outputs': 0}
    return store.run('00_unavailable_freeze', freeze)


def _incomplete_reference_state(archive, study_id, missingness_path, frozen_design=None, frozen_outcomes=None):
    """Retain method facts and validated prior labels; never infer missing labels."""
    from .agent_workflow_2026_09_10 import sha
    archive, missingness_path = Path(archive).resolve(), Path(missingness_path).resolve()
    context = reference_archive_context(archive)
    matches = [case for case in reference_cases(archive) if case['study_id'] == study_id]
    if len(matches) != 1: raise ValueError('Incomplete reference must identify one enrolled study')
    case = matches[0]
    if case.get('preparation_status') == 'UNAVAILABLE_SHARED_PREREQUISITE':
        raise ValueError('Use the strict shared-prerequisite unavailable adapter for an unexecuted case')
    unavailable = case.get('preparation_status') == 'UNAVAILABLE_PROSPECTIVE_PREPARATION'
    prospective = context.get('enrollment') is not None
    archive_hashes = {str(p.relative_to(archive)): sha(p) for p in archive.rglob('*') if p.is_file()}
    dependencies = dict(case['source_sha256'])

    def pin(path, expected=None):
        path = str(Path(path).resolve())
        digest = sha(path)
        if expected is not None and digest != expected: raise ValueError('Incomplete reference evidence changed: ' + path)
        if path in dependencies and dependencies[path] != digest:
            raise ValueError('Conflicting incomplete-reference dependency identity: ' + path)
        dependencies[path] = digest
        return path

    def read_missingness(path):
        value = json.loads(Path(path).read_text())
        if (set(value) != {'schema', 'study_id', 'family_id', 'archive_config_sha256', 'phase', 'reason', 'evidence'}
                or value.get('schema') != 'clockbase-reference-missingness-v1'
                or value.get('study_id') != study_id or value.get('family_id') != case['family_id']
                or value.get('archive_config_sha256') != sha(archive / 'config.json')
                or value.get('phase') not in {'SOURCE_DESIGN', 'NUMERICAL_REVIEW'}
                or not isinstance(value.get('reason'), str) or not value['reason'].strip()
                or not isinstance(value.get('evidence'), list) or not value['evidence']):
            raise ValueError('Source-linked administrative reference missingness receipt required')
        for item in value['evidence']:
            if not isinstance(item, dict) or set(item) != {'path', 'sha256'}:
                raise ValueError('Missingness evidence requires exact file identities')
            pin(item['path'], item['sha256'])
        return value

    missingness_path = pin(missingness_path)
    missing = read_missingness(missingness_path)

    prior = {'SOURCE_DESIGN': frozen_design, 'NUMERICAL_REVIEW': frozen_outcomes}
    state_folder = archive / study_id / 'reference_state'
    if (state_folder / 'complete.json').is_file():
        state, returned = reference_stage(archive, study_id + '/reference_state', 'reference-state.json')
        if (state != returned or state.get('study_id') != study_id
                or not isinstance(state.get('references'), list) or not isinstance(state.get('missingness'), list)):
            raise ValueError('Prior reference-state receipt differs from its study')
        seen = set()
        for item in state['references']:
            if not isinstance(item, dict) or set(item) != {'phase', 'path', 'sha256'}:
                raise ValueError('Prior reference requires its phase and exact file identity')
            phase = item['phase']
            if phase not in prior or phase in seen: raise ValueError('Ambiguous prior reference phase')
            seen.add(phase); path = pin(item['path'], item['sha256'])
            if prior[phase] is not None and str(Path(prior[phase]).resolve()) != path:
                raise ValueError('Explicit prior reference would replace a sealed reference')
            prior[phase] = path
        recorded_missing = []
        for item in state['missingness']:
            if (not isinstance(item, dict) or set(item) != {'phase', 'path', 'sha256'}
                    or item['phase'] not in prior or item['phase'] in {phase for phase, _ in recorded_missing}):
                raise ValueError('Unknown or duplicate sealed missing-reference phase')
            path = pin(item['path'], item['sha256'])
            if read_missingness(path)['phase'] != item['phase']:
                raise ValueError('Sealed reference missingness differs from its declared phase')
            recorded_missing.append((item['phase'], path))
        if prospective and 'SOURCE_DESIGN' not in seen and 'SOURCE_DESIGN' not in {phase for phase, _ in recorded_missing}:
            raise ValueError('Prospective source-reference state requires a frozen design or pinned missingness')
        if prospective and frozen_design is not None and 'SOURCE_DESIGN' not in seen:
            raise ValueError('A late source-design judgment cannot replace the pre-numerical source state')
        if missing['phase'] == 'SOURCE_DESIGN' and ('SOURCE_DESIGN', missingness_path) not in recorded_missing:
            raise ValueError('Source-design absence differs from its sealed pre-numerical state')
    elif prospective:
        raise ValueError('Prospective incomplete references require a sealed prior reference-state stage')

    if unavailable:
        outcomes, mapping, methods = {}, [], {}
        for method in COMPARISON_METHODS:
            result, returned = reference_stage(archive, study_id + '/' + method + '/decision', 'method-result.json')
            if (result != returned or result.get('study_id') != study_id or result.get('family_id') != case['family_id']
                    or result.get('method') != method or result.get('preparation_attempt') != case['attempt']
                    or result.get('execution_status') != 'NOT_RUN_SHARED_PREREQUISITE'
                    or result.get('study_decision_status') != 'CANNOT VERIFY' or result.get('accepted_results') != 0
                    or any(result.get(key) != [] for key in ('final_outputs', 'candidates', 'checks'))
                    or any(result.get(key) is not False for key in ('repair_proposed', 'flagged_material_error', 'raw_model_alarm'))
                    or list((archive / 'api').glob(study_id + '-' + method + '-*'))):
                raise ValueError('Prospective preparation failure cannot erase method execution, outputs or flags')
            methods[method] = result
    else:
        outcomes, mapping, methods = reference_output_inventory(archive, case)
    labels = dict.fromkeys(outcomes, 'CANNOT VERIFY')
    design_constraints, assessor, known_sources = {}, None, []
    frozen = None
    if prior['NUMERICAL_REVIEW'] is not None:
        if unavailable: raise ValueError('Unavailable prepared inputs cannot have numerical reference outputs')
        outcome_path = pin(prior['NUMERICAL_REVIEW'])
        prior_outcome = json.loads(Path(outcome_path).read_text())
        if prior_outcome.get('schema') not in {'clockbase-frozen-reference-outcomes-v1', 'clockbase-frozen-reference-outcomes-v2'}:
            raise ValueError('Prior outcome evidence must be an independently frozen judgment')
        # Reject wrapper cycles before recomputing the independent judgment.
        read_reference_records(outcome_path)
        numerical = Path(prior_outcome['source_identities']['numerical_archive'])
        numerical_config = json.loads((numerical / 'config.json').read_text())
        if Path(numerical_config['archive']).resolve() != archive or prior_outcome['split'] != context['split']:
            raise ValueError('Prior outcome reference belongs to a different archive/split')
        if set(prior_outcome['labels']) != set(outcomes): raise ValueError('Prior reference omits or invents output identities')
        declared_design = str(Path(numerical_config['frozen_design_path']).resolve())
        if prospective and prior['SOURCE_DESIGN'] is None:
            raise ValueError('Numerical reference requires the source design sealed before execution')
        if prior['SOURCE_DESIGN'] is not None and str(Path(prior['SOURCE_DESIGN']).resolve()) != declared_design:
            raise ValueError('Prior outcome and source-design freezes differ')
        prior['SOURCE_DESIGN'] = declared_design
        labels.update(prior_outcome['labels'])
        assessor = prior_outcome['judgment']['reviewer']['name']
        known_sources.append(outcome_path)
    if prior['SOURCE_DESIGN'] is not None:
        design_path = pin(prior['SOURCE_DESIGN'])
        frozen = read_reference_design(design_path)
        _require_registered_prospective_design(archive, design_path, frozen, context)
        if (frozen['study_id'] != study_id or frozen['family_id'] != case['family_id']
                or frozen['split'] != context['split']
                or prospective and frozen['schema'] != 'clockbase-frozen-reference-design-v2'
                or not unavailable and frozen['packet']['source_packet'] != reference_source_view(case['packet'])):
            raise ValueError('Prior source reference differs from the enrolled method source')
        assessor = assessor or frozen['judgment']['reviewer']['name']
        known_sources.append(design_path)
        for oid, outcome in outcomes.items():
            design_constraints[oid] = reference_output_design(frozen, outcome)
            if design_constraints[oid]['status'] == 'FAIL': labels[oid] = 'FAIL'
    facts = reference_method_facts(archive, case, methods, mapping)
    records = []
    for fact in facts:
        row = {key: value for key, value in fact.items() if not key.startswith('_')}
        resources = _reference_record_resources(archive, row, prospective=prospective)
        if unavailable:
            resources.update(preparation_status=case['preparation_status'], preparation_attempt=case['attempt'])
        records.append({**resources, **_reference_fact_labels(fact, labels),
            'reference_role': 'reference_incomplete', 'reference_availability': 'MISSING_OR_PARTIAL_REFERENCE',
            'reference_adjudicator': assessor, 'reference_source': missingness_path,
            'final_reference_source': missingness_path, 'reference_missing_phase': missing['phase'],
            'known_reference_sources': known_sources,
            'reference_output_ids': fact['_all_output_ids'],
            'reference_scope': 'Archived decisions retained; labels use only validated prior judgments and material design constraints; remaining labels are unresolved'})
    if archive_hashes != {str(p.relative_to(archive)): sha(p) for p in archive.rglob('*') if p.is_file()}:
        raise ValueError('Method archive changed while retaining incomplete references')
    for path, digest in dependencies.items():
        if sha(path) != digest: raise ValueError('Reference dependency changed while freezing: ' + path)
    return {'study_id': study_id, 'family_id': case['family_id'], 'split': context['split'],
            'archive': str(archive), 'archive_sha256': archive_hashes, 'source_sha256': dependencies,
            'missingness': missing, 'prior_references': prior, 'known_reference_sources': known_sources,
            'known_source_design': ({'domains': frozen['judgment']['domains'],
                                     'comparisons': frozen['judgment']['comparisons']} if frozen is not None else None),
            'labels': labels, 'source_design_constraints': design_constraints, 'records': records}


def freeze_reference_incomplete(archive, study_id, missingness_path, output, *, frozen_design=None, frozen_outcomes=None):
    """Freeze missing/partial references without resetting successful methods."""
    from .agent_workflow_2026_09_10 import sha, write_new
    from src.pipeline_runtime_2026_09_10 import StageStore
    archive, output = Path(archive).resolve(), Path(output).resolve()
    inputs = {'archive': str(archive), 'study_id': study_id,
              'missingness_path': str(Path(missingness_path).resolve()),
              'frozen_design': str(Path(frozen_design).resolve()) if frozen_design is not None else None,
              'frozen_outcomes': str(Path(frozen_outcomes).resolve()) if frozen_outcomes is not None else None}
    if (output == archive or output.is_relative_to(archive) or archive.is_relative_to(output)
            or any(Path(value).is_relative_to(output) for key, value in inputs.items()
                   if key not in {'study_id', 'archive'} and value is not None)):
        raise ValueError('Incomplete reference freeze must be separate from its method archive and input evidence')
    state = _incomplete_reference_state(**inputs)
    store = StageStore(output, {'version': 'clockbase-reference-incomplete-v1',
        'implementation_sha256': sha(__file__), 'inputs': inputs, 'state': state})
    def freeze(directory):
        value = {'schema': 'clockbase-reference-incomplete-v1', 'phase': 'REFERENCE_INCOMPLETE',
            'split': state['split'], 'inputs': inputs, 'state': state,
            'reference_labels_completed': False, 'scope': 'Administrative missingness with retained independent evidence; no new scientific judgment'}
        path = directory / 'frozen-outcomes.json'; write_new(path, value)
        path.with_suffix('.sha256').write_text(sha(path) + '\n')
        return {'phase': 'REFERENCE_INCOMPLETE', 'frozen_outcomes': str(path), 'sha256': sha(path),
                'reference_labels_completed': False, 'reference_resolved_outputs': sum(s != 'CANNOT VERIFY' for s in state['labels'].values())}
    return store.run('00_incomplete_freeze', freeze)


def unresolved_design(study, reason):
    return {'study_id':study,'status':'CANNOT VERIFY','case_ids':[],'control_ids':[],
            'excluded_ids':[],'biological_unit':'CANNOT VERIFY','pairing':'unknown','findings':[],
            'unresolved':[reason],'evidence':[],'proposed_action':'obtain_source_evidence'}


def deterministic_design(packet, primary):
    """Use the same explicit source contract as the audit, without a model."""
    from .source_packet import verify_source_design
    from .agent_workflow_2026_09_10 import DESIGN_SET_SCHEMA
    if primary.get('comparisons') is not None:
        comparisons = []
        for item in primary['comparisons']:
            child = deterministic_design({**packet, 'target': item['target']},
                {'rows': [row for row in primary['rows'] if row['comparison_id'] == item['comparison_id']]})
            comparisons.append({**child, **item})
        assigned = {sid for item in comparisons for sid in item['case_ids'] + item['control_ids']}
        return {'schema': DESIGN_SET_SCHEMA, 'study_id': packet['study_id'],
                'status': 'DESIGN_PROPOSED' if comparisons else 'CANNOT VERIFY',
                'comparisons': comparisons,
                'unassigned_ids': sorted({row['sample_id'] for row in packet['samples']} - assigned),
                'findings': [], 'evidence': [], 'unresolved': [] if comparisons else ['No declared producer comparisons']}
    groups = {(tuple(sorted(row['case_ids'])),tuple(sorted(row['control_ids'])))
              for row in primary['rows'] if row['status']=='RESULT'}
    fallback = unresolved_design(packet['study_id'],'No unique source-supported independent contrast for deterministic repair')
    if len(groups)!=1:return fallback
    case_ids,control_ids = (list(ids) for ids in next(iter(groups)))
    verified = verify_source_design(packet,case_ids,control_ids)
    if verified['status']!='PASS':return fallback
    return {**fallback,'status':'DESIGN_PROPOSED','case_ids':case_ids,'control_ids':control_ids,
            'excluded_ids':sorted({r['sample_id'] for r in packet['samples']}-set(case_ids+control_ids)),
            'biological_unit':'Explicit source-declared unit IDs','pairing':'independent','unresolved':[],
            'proposed_action':'recompute_welch'}


def compare_development_case(store, case, client_factory=None, *, preserved_reviews=None):
    """Run all four methods on identical preserved inputs; no reference labels."""
    from .agent_workflow_2026_09_10 import write_new, finding_alarm, finding_records, parse_design
    from collections import Counter
    from src.pipeline_2026_09_10 import assess_design, numerical_audit
    from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
    if case.get('preparation_status') == 'UNAVAILABLE_SHARED_PREREQUISITE':
        return compare_unavailable_case(store, case)
    study,packet,primary,clocks=(case[key] for key in ('study_id','packet','primary','clocks'))
    outputs=[]
    for method in COMPARISON_METHODS:
        prefix=study+'/'+method
        (store.root/prefix).mkdir(parents=True,exist_ok=True)
        started=time.monotonic()
        review=None
        if method in {'independent_audit','self_review'}:
            def review_stage(directory):
                try:
                    if preserved_reviews is None:
                        design=assess_design(packet,client_factory(study,method),directory,mode=method,
                                             primary=primary if method=='self_review' else None)
                    else:
                        saved=preserved_reviews[method]
                        if saved['source_packet'] != packet:
                            raise PipelineIntegrityError('Preserved review and comparison source packets differ')
                        if saved['result']['execution_status'] != 'PASS':
                            raise ValueError('Preserved model review failed: '+saved['result'].get('error','unspecified failure'))
                        design=parse_design(json.dumps(saved['model_design']),packet,
                                            primary=primary if method=='self_review' else None)
                        write_new(directory/'preserved-review.json',saved)
                    value={'execution_status':'PASS','design':design}
                except PipelineIntegrityError:
                    raise
                except Exception as error:
                    value={'execution_status':'CANNOT VERIFY','error':f'{type(error).__name__}: {error}',
                           'design':unresolved_design(study,'Model review did not yield a verified source design')}
                value['review_execution_scope']=('CURRENT_MODEL_CALL' if preserved_reviews is None
                                                  else 'REPARSE_PRESERVED_MODEL_RESPONSE')
                write_new(directory/'review-result.json',value)
                return value
            review=store.run(prefix+'/review',review_stage)

        def decision_stage(directory):
            if method=='existing_pipeline':
                retained=[{**row,'accepted_result':row['status']=='RESULT'} for row in primary['rows']]
                count=sum(row['accepted_result'] for row in retained)
                value={'study_decision_status':'PASS' if count else 'CANNOT VERIFY',
                       'accepted_results':count,'final_outputs':retained,'checks':[],
                       'repair_proposed':False,'flagged_material_error':False,
                       'historical_release_state':case.get('historical_release_state','CANNOT VERIFY'),
                       'decision_scope':'Prespecified no-audit retention of producer RESULT rows; not a reconstructed historical release verdict'}
                write_new(directory/'original-analysis.json',primary)
            else:
                design=review['design'] if review else deterministic_design(packet,primary)
                try:
                    value=numerical_audit(design,primary,clocks,directory,packet=packet)
                    value['repair_proposed']=any(check['repair_kind']=='CORRECTION_PROPOSAL' for check in value['checks'])
                    value['new_analysis_proposed']=sum(check['repair_kind'] in {
                        'NEW_ANALYSIS_AFTER_ABSTENTION', 'NEW_ANALYSIS_FOR_MISSING_COMPARISON'} for check in value['checks'])
                    value['finding_categories']=dict(Counter(item['category'] for item in finding_records(design)))
                    # Compare allegations with a separate reference. Source
                    # warnings and candidate failures have different subjects.
                    value['raw_model_alarm']=finding_alarm(design)
                    value['flagged_material_error']=any(check['producer_error_allegation'] for check in value['checks'])
                    value['failed_proposal_results']=sum(bool(check['proposal_errors']) for check in value['checks'])
                except PipelineIntegrityError:
                    raise
                except Exception as error:
                    value={'study_decision_status':'CANNOT VERIFY','accepted_results':0,'final_outputs':[],
                           'checks':[],'repair_proposed':False,'flagged_material_error':False,
                           'error':f'{type(error).__name__}: {error}'}
            value.update(study_id=study,family_id=case['family_id'],method=method,
                         independent_reference_status='CANNOT VERIFY',accepted_result_error_rate=None)
            value.setdefault('raw_model_alarm', False)
            value.setdefault('failed_proposal_results', 0)
            value['error_signal_definition']=ERROR_SIGNAL_DEFINITION
            write_new(directory/'method-result.json',value)
            return value
        decision=store.run(prefix+'/decision',decision_stage)
        outputs.append({'study_id':study,'family_id':case['family_id'],'method':method,
                        'decision_status':decision['study_decision_status'],
                        'accepted':decision['accepted_results']>0,'retained_clock_results':decision['accepted_results'],
                        'flagged_material_error':decision['flagged_material_error'],
                        'raw_model_alarm':decision['raw_model_alarm'],
                        'failed_proposal_results':decision['failed_proposal_results'],
                        'error_signal_definition':decision['error_signal_definition'],
                        'repair_proposed':decision['repair_proposed'],
                        'review_execution_status':review['execution_status'] if review else 'NOT_APPLICABLE',
                        'execution_error':decision.get('error',review.get('error') if review else None),
                        'independent_reference_status':'CANNOT VERIFY','accepted_result_error_rate':None,
                        'current_invocation_seconds':time.monotonic()-started,
                        'method_result':decision})
    return outputs


def compare_unavailable_case(store, case):
    """Retain all method slots without implying that their analyses executed."""
    from .agent_workflow_2026_09_10 import write_new
    outputs = []
    for method in COMPARISON_METHODS:
        prefix = case['study_id'] + '/' + method
        (store.root / prefix).mkdir(parents=True, exist_ok=True)
        def decision(directory):
            value = {'study_id': case['study_id'], 'family_id': case['family_id'], 'method': method,
                'execution_status': 'NOT_RUN_SHARED_PREREQUISITE', 'preparation_attempt': case['attempt'],
                'study_decision_status': 'CANNOT VERIFY', 'accepted_results': 0, 'final_outputs': [],
                'candidates': [], 'checks': [], 'repair_proposed': False, 'flagged_material_error': False,
                'raw_model_alarm': False, 'failed_proposal_results': 0,
                'error_signal_definition': ERROR_SIGNAL_DEFINITION,
                'decision_scope': 'Shared prerequisite unavailable; no numerical result or scientific reference is inferred'}
            write_new(directory / 'method-result.json', value)
            return value
        result = store.run(prefix + '/decision', decision)
        outputs.append({'study_id': case['study_id'], 'family_id': case['family_id'], 'method': method,
            'decision_status': 'CANNOT VERIFY', 'accepted': False, 'retained_clock_results': 0,
            'flagged_material_error': False, 'raw_model_alarm': False, 'failed_proposal_results': 0,
            'error_signal_definition': ERROR_SIGNAL_DEFINITION, 'repair_proposed': False,
            'review_execution_status': 'NOT_RUN_SHARED_PREREQUISITE',
            'execution_error': 'Recorded shared prerequisite is unavailable',
            'independent_reference_status': 'CANNOT VERIFY', 'accepted_result_error_rate': None,
            'current_invocation_seconds': None, 'method_result': result})
    return outputs


def replay_development_comparison(archive, output):
    """Reparse frozen reviews and recheck decisions, with no model/clock run."""
    from .agent_workflow_2026_09_10 import sha, write_new
    from src.pipeline_2026_09_10 import read_primary_outputs
    from src.pipeline_runtime_2026_09_10 import StageStore
    archive, output = Path(archive).resolve(), Path(output).resolve()
    if output == archive or output.is_relative_to(archive) or archive.is_relative_to(output):
        raise ValueError('Replay output must be separate from the preserved comparison')
    archived_files = {str(path): sha(path) for path in archive.rglob('*') if path.is_file()}
    if not archived_files:
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {archive}/config.json')
    old_config = json.loads((archive/'config.json').read_text())
    if old_config.get('development_only') is not True:
        raise ValueError('Only exposed development comparisons can be replayed with changed methods')

    def stage(name, filename):
        directory = (archive/name).resolve()
        if not directory.is_relative_to(archive):
            raise ValueError('Archived stage escapes the comparison directory')
        if sha(directory/'complete.json') != (directory/'complete.sha256').read_text().strip():
            raise ValueError('Archived stage receipt changed: '+name)
        receipt = json.loads((directory/'complete.json').read_text())
        for relative, digest in receipt['outputs'].items():
            path = (directory/relative).resolve()
            if not path.is_relative_to(directory) or archived_files.get(str(path)) != digest:
                raise ValueError('Archived stage output is missing, changed or outside its stage: '+relative)
        matches = [relative for relative in receipt['outputs'] if Path(relative).name == filename]
        if len(matches) != 1:
            raise ValueError('Missing or ambiguous preserved comparison artifact: '+filename)
        return json.loads((directory/matches[0]).read_text()), receipt['result']

    cases, receipt = stage('00_inputs', 'prepared-cases.json')
    if cases != receipt['cases'] or not cases or len({case['study_id'] for case in cases}) != len(cases):
        raise ValueError('Prepared cases differ from their receipt or have duplicate identities')
    external_inputs, saved_reviews, old_decisions = {}, {}, {}
    for case in cases:
        external_inputs.update(case['source_sha256'])
        for path, digest in case['source_sha256'].items():
            if sha(path) != digest:
                raise ValueError('Preserved comparison source changed: '+path)
        primary = case['primary']
        for path in (primary['path'], case['clocks']['path'], primary.get('manifest_path')):
            if path is not None and str(Path(path).resolve()) not in case['source_sha256']:
                raise ValueError('Comparison file lacks a preserved source hash: '+path)
        checked = read_primary_outputs(primary['path'], case['clocks'], case['packet'], primary.get('manifest_path'))
        if checked['rows'] != primary['rows']:
            raise ValueError('Prepared producer rows differ from their source CSV')
        if primary.get('manifest_path') and any(checked[key] != primary[key]
                for key in ('schema', 'comparisons', 'unassigned_ids')):
            raise ValueError('Prepared producer comparisons differ from their manifest')
        for method in COMPARISON_METHODS:
            prefix = case['study_id']+'/'+method
            decision, declared = stage(prefix+'/decision', 'method-result.json')
            if decision != declared:
                raise ValueError('Preserved method result differs from its receipt')
            old_decisions[case['study_id'], method] = decision
            if method in {'independent_audit', 'self_review'}:
                result, declared = stage(prefix+'/review', 'review-result.json')
                if result != declared:
                    raise ValueError('Preserved review differs from its receipt')
                source, raw, location = case['packet'], None, 'PRESERVED_FAILED_REVIEW'
                if result['execution_status'] == 'PASS':
                    receipt = json.loads((archive/prefix/'review/complete.json').read_text())
                    names = {Path(relative).name for relative in receipt['outputs']}
                    if {'source-packet.json', 'model-design.json'} <= names:
                        source = stage(prefix+'/review', 'source-packet.json')[0]
                        raw = stage(prefix+'/review', 'model-design.json')[0]
                        location = 'REVIEW_STAGE_ARTIFACTS'
                    elif names == {'design.json', 'review-result.json'}:
                        # Earlier assess_design saved the normalized design only.
                        # Recover its source and raw response from the bounded API
                        # record, and bind both to this exact study/method/result.
                        slots = list((archive/'api').glob(case['study_id']+'-'+method+'-*'))
                        if len(slots) != 1:
                            raise ValueError('Legacy review needs one unambiguous archived API response')
                        slot = slots[0]
                        request = json.loads((slot/'request.json').read_text())
                        usage = json.loads((slot/'receipt.json').read_text())
                        if sha(slot/'response.json') != usage['response_sha256']:
                            raise ValueError('Legacy review response differs from its API receipt')
                        if [message['role'] for message in request['messages']] != ['system', 'user']:
                            raise ValueError('Legacy review has an unsupported message layout')
                        body = json.loads(request['messages'][1]['content'])
                        source = body if method == 'independent_audit' else body['source_packet']
                        if method == 'self_review' and body['preceding_analysis']['rows'] != primary['rows']:
                            raise ValueError('Legacy self-review received different producer rows')
                        response = json.loads((slot/'response.json').read_text())
                        raw = json.loads(response['choices'][0]['message']['content'])
                        location = str(slot)
                    else:
                        raise ValueError('Incomplete or unsupported preserved review layout')
                    from .agent_workflow_2026_09_10 import parse_design
                    if source != case['packet'] or parse_design(json.dumps(raw), source,
                            primary=primary if method == 'self_review' else None) != result['design']:
                        raise ValueError('Raw preserved review differs from its source packet or normalized result')
                saved_reviews.setdefault(case['study_id'], {})[method] = {
                    'result': result, 'source_packet': source, 'model_design': raw,
                    'raw_response_location': location,
                    'source_review_stage': str(archive/prefix/'review')}
    root = Path(__file__).resolve().parents[1]
    # The producer is not executed. Pin every module used by replay decisions.
    modules = [root/'audit_system'/name for name in ('evaluation_2026_09_10.py',
        'agent_workflow_2026_09_10.py', 'source_packet.py', 'adapters/measurement_pilot.py',
        'claim_gate.py', 'policy.py', 'models.py', 'ledger.py', 'hashing.py', 'dependency_graph.py')]
    modules += [root/'src'/name for name in ('pipeline_2026_09_10.py', 'pipeline_runtime_2026_09_10.py', 'pipeline_inputs_2026_09_10.py')]
    config = {'archive': str(archive), 'archive_sha256': archived_files, 'source_sha256': external_inputs,
              'implementation_sha256': {str(path): sha(path) for path in modules},
              'error_signal_definition': ERROR_SIGNAL_DEFINITION, 'development_only': True,
              'generation_calls': 0, 'clock_execution_performed': False}
    output.mkdir(parents=True, exist_ok=True)
    with (output/'.run.lock').open('a') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('Another replay process owns this run')
        store = StageStore(output, config)
        rows = []
        for case in cases:
            rows.extend(compare_development_case(store, case, preserved_reviews=saved_reviews[case['study_id']]))
        def summarize(directory):
            records = []
            for row in rows:
                previous = old_decisions[row['study_id'], row['method']]
                value = {key: item for key, item in row.items() if key not in {'method_result', 'current_invocation_seconds'}}
                value.update(previous_method_flag=previous['flagged_material_error'],
                             previous_retained_results=previous['accepted_results'],
                             finding_attribution_counts=dict(Counter(
                                 item['attribution'] for check in row['method_result']['checks']
                                 for item in check['finding_attribution'])))
                records.append(value)
            result = {'records': records, 'study_count': len(cases),
                      'observed_family_count': len({case['family_id'] for case in cases}),
                      'error_signal_definition': ERROR_SIGNAL_DEFINITION,
                      'scope': 'Changed-code replay of preserved development responses; no fresh model or clock execution',
                      'previous_flags_scope': 'Historical method definitions retained; changes are not false-alarm or accuracy estimates',
                      'generation_calls': 0, 'estimated_incremental_cost_usd': 0,
                      'independent_reference_available': False, 'held_out_studies_evaluated': 0,
                      'accepted_result_error_rate': None, 'clock_execution_performed': False}
            write_new(directory/'replay-summary.json', result)
            with (directory/'replay-records.csv').open('x', newline='') as handle:
                writer = csv.DictWriter(handle, fieldnames=list(records[0]), lineterminator='\n'); writer.writeheader()
                writer.writerows({key: json.dumps(item) if isinstance(item, dict) else item
                                  for key, item in value.items()} for value in records)
            return result
        result = store.run('99_summary', summarize)
        if archived_files != {str(path): sha(path) for path in archive.rglob('*') if path.is_file()}:
            raise ValueError('Preserved comparison changed during replay')
        for path, digest in (external_inputs | config['implementation_sha256']).items():
            if sha(path) != digest:
                raise ValueError('Replay source or implementation changed: '+path)
        return result


def run_comparison_plan(plan_path, output, *, prepare_only=False, client_factory=None):
    """Execute a pinned, bounded development comparison with resumable stages."""
    from .agent_workflow_2026_09_10 import sha, write_new, review_messages
    from .report_2026_09_10 import load_archived_comparison
    from src.pipeline_runtime_2026_09_10 import BoundedOpenAI, StageStore
    plan_path,output=Path(plan_path).resolve(),Path(output).resolve()
    plan=json.loads(plan_path.read_text())
    if plan.get('schema') not in {'clockbase-development-comparison-v1', 'clockbase-development-comparison-v2'} or plan.get('split')!='development':
        raise ValueError('This adapter requires an explicit development comparison plan')
    if plan.get('provider')!='openai' or plan.get('model')!='gpt-5.6-sol':
        raise ValueError('This comparator adapter has no provider/model substitution')
    if plan.get('methods')!=list(COMPARISON_METHODS) or plan.get('review_calls_per_model_arm')!=1:
        raise ValueError('Matched comparison requires the declared four arms and one call per model arm')
    cases=plan['cases']; ids=[case['study_id'] for case in cases]
    if not ids or len(ids)!=len(set(ids)) or any(not re.fullmatch(r'GSE\d+',study) for study in ids):
        raise ValueError('Comparison requires unique enrolled GEO studies')
    if plan['max_calls']!=2*len(cases):raise ValueError('API allocation must match both model arms for every case')
    if plan.get('baseline_rule')!='retain_producer_RESULT_rows_without_audit':
        raise ValueError('Baseline retention rule must be declared before execution')
    limits={key:plan[key] for key in ('max_calls','max_cost','output_tokens','prompt_bytes','input_per_million','output_per_million')}
    if any(type(value) not in (int,float) or not math.isfinite(value) or value<=0 for value in limits.values()):
        raise ValueError('Positive finite API limits and explicit tariff bounds are required')
    if any(type(limits[key]) is not int for key in ('max_calls','output_tokens','prompt_bytes')):
        raise ValueError('Call, token and prompt-byte limits must be integers')
    family_path=Path(plan['family_manifest']).resolve()
    opener=gzip.open if family_path.suffix=='.gz' else open
    with opener(family_path,'rt') as handle:families=json.load(handle)
    family_by_study={row['study_id']:row['family_id'] for row in families['studies']}
    enrollment=[{'study_id':study,'family_id':family_by_study.get(study)} for study in ids]
    verify_family_enrollment(enrollment,families,'development')
    archive_manifest=Path(plan['archive_identity_manifest']).resolve()
    expected=json.loads(archive_manifest.read_text())['source_sha256']
    if plan['schema'] == 'clockbase-development-comparison-v2':
        return run_attempted_comparison(plan_path, output, plan, limits, family_path, archive_manifest,
            expected, enrollment, prepare_only=prepare_only, client_factory=client_factory)
    prepared=[{**load_archived_comparison(case,expected),'family_id':family_by_study[case['study_id']]} for case in cases]
    reservations=[]
    for case in prepared:
        for method in ('independent_audit','self_review'):
            messages=review_messages(case['packet'],mode=method,primary=case['primary'] if method=='self_review' else None)
            wire=[{'role':{'system':'system','human':'user'}[message.type],'content':message.content} for message in messages]
            request={'model':plan['model'],'messages':wire,'max_completion_tokens':limits['output_tokens'],
                     'reasoning_effort':'low','store':False}
            size=len(json.dumps(request,ensure_ascii=False).encode())
            if size>limits['prompt_bytes']:raise ValueError('Comparison prompt exceeds its declared byte limit before any call')
            reserve=((size+100*len(wire))*limits['input_per_million']+limits['output_tokens']*limits['output_per_million'])/1e6
            reservations.append({'study_id':case['study_id'],'method':method,'request_bytes':size,'upper_estimated_cost_usd':reserve})
    if sum(row['upper_estimated_cost_usd'] for row in reservations)>limits['max_cost']:
        raise ValueError('Declared budget cannot reserve both model arms for every selected study')
    inputs={str(path):sha(path) for path in (plan_path,family_path,archive_manifest)}
    for case in prepared:inputs.update(case['source_sha256'])
    root=Path(__file__).resolve().parents[1]
    code=[root/'audit_system'/name for name in ('evaluation_2026_09_10.py','report_2026_09_10.py',
        'agent_workflow_2026_09_10.py','source_packet.py','adapters/measurement_pilot.py','claim_gate.py',
        'policy.py','models.py','ledger.py','hashing.py','dependency_graph.py')]
    code += [root/'src'/name for name in ('pipeline_2026_09_10.py','pipeline_runtime_2026_09_10.py','pipeline_inputs_2026_09_10.py')]
    config={'plan_sha256':sha(plan_path),'input_sha256':inputs,
            'implementation_sha256':{str(path):sha(path) for path in code},
            'resource_limits':limits,'provider':plan['provider'],'model':plan['model'],
            'concurrency':1,'automatic_generation_retries':0,'development_only':True,'clock_execution_performed':False}
    output.mkdir(parents=True,exist_ok=True)
    with (output/'.run.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise RuntimeError('Another comparison process owns this run')
        store=StageStore(output,config)
        def freeze_inputs(directory):
            write_new(directory/'comparison-plan.json',plan)
            write_new(directory/'prepared-cases.json',prepared)
            write_new(directory/'matched-request-reservations.json',reservations)
            return {'cases':prepared,'frozen_before_review':True,'reference_labels_supplied':False}
        frozen=store.run('00_inputs',freeze_inputs)
        if prepare_only:return {'prepared':True,'studies':ids,'max_calls':plan['max_calls'],'max_estimated_cost_usd':plan['max_cost']}
        if client_factory is None:
            def client_factory(study,method):
                return BoundedOpenAI(output/'api',study+'-'+method,model=plan['model'],**limits)
        rows=[]
        for case in frozen['cases']:rows.extend(compare_development_case(store,case,client_factory))
        if len(rows)!=len(cases)*len(COMPARISON_METHODS):raise ValueError('A selected method/case is missing')
        for path,digest in inputs.items():
            if sha(path)!=digest:raise ValueError('Comparison input changed during execution: '+path)
        for path,digest in config['implementation_sha256'].items():
            if sha(path)!=digest:raise ValueError('Comparison implementation changed during execution')

        def summarize(directory):
            summaries=[comparison_resource_record(output, row) for row in rows]
            write_new(directory/'comparison-records.json',summaries)
            with (directory/'comparison-records.csv').open('x',newline='') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(summaries[0]),lineterminator='\n');writer.writeheader()
                writer.writerows({key:json.dumps(value) if isinstance(value,list) else value for key,value in row.items()} for row in summaries)
            methods={method:{'attempted_studies':len(cases),
                'retained_studies':sum(row['accepted'] for row in summaries if row['method']==method),
                'method_flags':sum(row['flagged_material_error'] for row in summaries if row['method']==method),
                'raw_model_alarm_studies':sum(row['raw_model_alarm'] for row in summaries if row['method']==method),
                'failed_proposal_results':sum(row['failed_proposal_results'] for row in summaries if row['method']==method),
                'proposed_repair_studies':sum(row['repair_proposed'] for row in summaries if row['method']==method),
                'unresolved_studies':sum(row['decision_status']=='CANNOT VERIFY' for row in summaries if row['method']==method),
                'generation_attempts':sum(row['generation_attempts'] for row in summaries if row['method']==method),
                'estimated_known_cost_usd':sum(row['estimated_known_cost_usd'] for row in summaries if row['method']==method),
                'accepted_result_error_rate':None,'false_alarms':None,'missed_material_errors':None,'repair_correctness':None}
                for method in COMPARISON_METHODS}
            result={'methods':methods,'records':summaries,'held_out_studies_evaluated':0,'independent_reference_available':False,
                    'error_signal_definition':ERROR_SIGNAL_DEFINITION,
                    'estimated_known_cost_usd':sum(row['estimated_known_cost_usd'] for row in summaries),
                    'generation_attempts':sum(row['generation_attempts'] for row in summaries),
                    'clock_execution_performed':False,'historical_release_state':'CANNOT VERIFY',
                    'resource_match':'Same model, one call, identical per-call ceilings and verifier; realized tokens/time are reported separately',
                    'api_sha256':{str(path.relative_to(output)):sha(path) for path in (output/'api').rglob('*') if path.is_file()}}
            write_new(directory/'comparison-summary.json',result)
            return result
        result=store.run('99_summary',summarize)
        observed={str(path.relative_to(output)):sha(path) for path in (output/'api').rglob('*') if path.is_file()}
        if observed!=result['api_sha256']:raise ValueError('API evidence changed after comparison completion')
        return result


def comparison_resource_record(output, row):
    """Read the same measured/unknown review accounting for either plan version."""
    actor = row['study_id'] + '-' + row['method']
    slots = sorted((output / 'api').glob(actor + '-*'))
    receipts = [json.loads((slot / 'receipt.json').read_text()) for slot in slots if (slot / 'receipt.json').is_file()]
    stages = [json.loads(path.read_text()) for path in (output / row['study_id'] / row['method']).glob('*/complete.json')]
    return {key: value for key, value in row.items() if key not in {'current_invocation_seconds', 'method_result'}} | {
        'runtime_seconds': sum(stage['elapsed_seconds'] for stage in stages),
        'generation_attempts': sum((slot / 'reservation.json').is_file() for slot in slots),
        'responses_with_usage': len(receipts), 'unknown_billing_attempts': sum(not (slot / 'receipt.json').is_file() for slot in slots),
        'estimated_known_cost_usd': sum(receipt['estimated_cost_usd'] for receipt in receipts),
        'billed_cost_usd': None if slots else 0,
        'prompt_tokens': sum(receipt['usage']['prompt_tokens'] for receipt in receipts),
        'completion_tokens': sum(receipt['usage']['completion_tokens'] for receipt in receipts),
        'returned_models': sorted({receipt['returned_model'] for receipt in receipts}),
        'api_evidence_scope': 'Incremental review only; original producer execution is shared and not rerun'}


def run_attempted_comparison(plan_path, output, plan, limits, family_path, archive_manifest,
                             expected, enrollment, *, prepare_only=False, client_factory=None):
    """Versioned development adapter with complete terminal prerequisite attempts."""
    from .agent_workflow_2026_09_10 import sha, write_new, review_messages
    from .report_2026_09_10 import load_archived_comparison, load_unavailable_comparison
    from src.pipeline_runtime_2026_09_10 import StageStore, BoundedOpenAI
    if plan.get('shared_prerequisite_failure_rule') != 'retain_all_methods_without_execution':
        raise ValueError('Attempted comparison requires a declared shared-prerequisite rule')
    if not isinstance(expected, dict) or not expected:
        raise ValueError('Attempted comparison needs pinned source files, including failed-attempt receipts')
    inputs = {str(path.resolve()): sha(path) for path in (plan_path, family_path, archive_manifest)}
    for name, digest in expected.items():
        if not Path(name).is_absolute() or not isinstance(digest, str) or sha(name) != digest:
            raise ValueError('Attempted comparison source identity changed: ' + name)
        inputs[name] = digest
    root = Path(__file__).resolve().parents[1]
    code = [root / 'audit_system' / name for name in ('evaluation_2026_09_10.py', 'report_2026_09_10.py',
        'agent_workflow_2026_09_10.py', 'source_packet.py', 'adapters/measurement_pilot.py', 'claim_gate.py',
        'policy.py', 'models.py', 'ledger.py', 'hashing.py', 'dependency_graph.py')]
    code += [root / 'src' / name for name in ('pipeline_2026_09_10.py', 'pipeline_runtime_2026_09_10.py',
                                            'pipeline_inputs_2026_09_10.py')]
    config = {'plan_schema': plan['schema'], 'plan_sha256': sha(plan_path), 'input_sha256': inputs,
        'implementation_sha256': {str(path): sha(path) for path in code}, 'resource_limits': limits,
        'provider': plan['provider'], 'model': plan['model'], 'concurrency': 1,
        'automatic_generation_retries': 0, 'development_only': True, 'clock_execution_performed': False}
    output.mkdir(parents=True, exist_ok=True)
    with (output / '.run.lock').open('a') as lock:
        try: fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError: raise RuntimeError('Another comparison process owns this run')
        store = StageStore(output, config)
        def enroll(directory):
            value = {'schema': 'clockbase-comparison-enrollment-v1', 'split': 'development',
                     'methods': list(COMPARISON_METHODS), 'studies': enrollment}
            path = directory / 'enrollment.json'; write_new(path, value)
            write_new(directory / 'comparison-plan.json', plan)
            return {'enrollment': value, 'sha256': sha(path), 'frozen_before_case_preparation': True}
        enrolled = store.run('00_enrollment', enroll)
        def prepare(directory):
            prepared, reservations = [], []
            by_study = {row['study_id']: row['family_id'] for row in enrollment}
            for declared in plan['cases']:
                if 'failure_receipt' in declared:
                    case = load_unavailable_comparison({**declared, 'family_id': by_study[declared['study_id']]},
                                                       expected, enrolled['sha256'])
                else:
                    case = {**load_archived_comparison(declared, expected),
                            'family_id': by_study[declared['study_id']]}
                if any(expected.get(path) != digest for path, digest in case['source_sha256'].items()):
                    raise ValueError('Prepared case used an undeclared source identity')
                prepared.append(case)
                write_new(directory / (case['study_id'] + '-attempt.json'), {
                    'study_id': case['study_id'], 'family_id': case['family_id'],
                    'enrollment_sha256': enrolled['sha256'],
                    'status': case.get('preparation_status', 'PASS'),
                    'source_sha256': case['source_sha256'], 'unavailable_attempt': case.get('attempt')})
                if case.get('preparation_status') == 'UNAVAILABLE_SHARED_PREREQUISITE': continue
                for method in ('independent_audit', 'self_review'):
                    messages = review_messages(case['packet'], mode=method,
                                              primary=case['primary'] if method == 'self_review' else None)
                    wire = [{'role': {'system': 'system', 'human': 'user'}[m.type], 'content': m.content} for m in messages]
                    request = {'model': plan['model'], 'messages': wire, 'max_completion_tokens': limits['output_tokens'],
                               'reasoning_effort': 'low', 'store': False}
                    size = len(json.dumps(request, ensure_ascii=False).encode())
                    if size > limits['prompt_bytes']: raise ValueError('Comparison prompt exceeds its declared byte limit before any call')
                    reserve = ((size + 100 * len(wire)) * limits['input_per_million']
                               + limits['output_tokens'] * limits['output_per_million']) / 1e6
                    reservations.append({'study_id': case['study_id'], 'method': method,
                                         'request_bytes': size, 'upper_estimated_cost_usd': reserve})
            if sum(row['upper_estimated_cost_usd'] for row in reservations) > limits['max_cost']:
                raise ValueError('Declared budget cannot reserve all runnable method calls')
            if ('expected_runnable_model_calls' in plan
                    and plan['expected_runnable_model_calls'] != len(reservations)):
                raise ValueError('Runnable calls differ from the explicitly frozen execution scope')
            write_new(directory / 'prepared-cases.json', prepared)
            write_new(directory / 'matched-request-reservations.json', reservations)
            return {'cases': prepared, 'enrollment_sha256': enrolled['sha256'],
                    'reference_labels_supplied': False, 'maximum_runnable_calls': len(reservations)}
        frozen = store.run('00_inputs', prepare)
        if prepare_only:
            return {'prepared': True, 'studies': [row['study_id'] for row in enrollment],
                    'unavailable_studies': sum(case.get('preparation_status') == 'UNAVAILABLE_SHARED_PREREQUISITE'
                                              for case in frozen['cases']),
                    'maximum_runnable_calls': frozen['maximum_runnable_calls'],
                    'max_calls': plan['max_calls'], 'max_estimated_cost_usd': plan['max_cost']}
        if client_factory is None:
            def client_factory(study, method):
                return BoundedOpenAI(output / 'api', study + '-' + method, model=plan['model'], **limits)
        rows = [row for case in frozen['cases'] for row in compare_development_case(store, case, client_factory)]
        validate_attempted_universe(rows, enrolled['enrollment'])
        for path, digest in {**inputs, **config['implementation_sha256']}.items():
            if sha(path) != digest: raise ValueError('Attempted comparison source or implementation changed: ' + path)
        def summarize(directory):
            summaries = [comparison_resource_record(output, row) for row in rows]
            write_new(directory / 'comparison-records.json', summaries)
            api_identities = {str(path.relative_to(output)): sha(path) for path in (output / 'api').rglob('*') if path.is_file()}
            result = {'records': summaries, 'attempted_studies': len(enrollment), 'attempted_method_records': len(rows),
                'unavailable_prerequisite_studies': sum(case.get('preparation_status') == 'UNAVAILABLE_SHARED_PREREQUISITE'
                                                      for case in frozen['cases']),
                'generation_attempts': len(list((output / 'api').glob('*/reservation.json'))),
                'estimated_known_cost_usd': sum(row['estimated_known_cost_usd'] for row in summaries),
                'unknown_billing_attempts': sum(row['unknown_billing_attempts'] for row in summaries),
                'reference_labels_completed': False, 'independent_reference_available': False,
                'held_out_studies_evaluated': 0, 'clock_execution_performed': False,
                'api_sha256': api_identities, 'resource_scope': 'Incremental review only; upstream cost remains separate'}
            write_new(directory / 'comparison-summary.json', result)
            return result
        result = store.run('99_summary', summarize)
        if {str(path.relative_to(output)): sha(path) for path in (output / 'api').rglob('*') if path.is_file()} != result['api_sha256']:
            raise ValueError('Attempted comparison API evidence changed after completion')
        return result


def build_families(records, exposed, seed=20260910):
    """Connect source-declared identities without reading outcomes or verdicts.

    Partial sources may establish links but cannot establish their absence.
    Linked accessions without records stay in the graph as unresolved nodes.
    Source completeness is relative to this input universe, not all GEO.
    """
    parent, nodes, owners, edges = {}, {}, {}, {}
    exposed = set(exposed)

    def node(accession):
        if not isinstance(accession, str) or not re.fullmatch(r'GSE\d+', accession):
            raise ValueError('Invalid family study identity')
        if accession not in parent:
            parent[accession] = accession
            nodes[accession] = {'sources': set(), 'samples': set(), 'primary_sample_sets': [],
                                'unresolved': [], 'primary_series': False}
        return nodes[accession]

    def root(accession):
        while parent[accession] != accession:
            parent[accession] = parent[parent[accession]]
            accession = parent[accession]
        return accession

    def join(a, b, kind, value, sources):
        node(a); node(b)
        if a == b:
            return
        left, right = sorted((root(a), root(b)))
        parent[right] = left
        key = (min(a,b), max(a,b), kind)
        witness = {'study_a': key[0], 'study_b': key[1], 'kind': kind,
                   'identifier': value, 'source_ids': sorted(set(sources))}
        if key not in edges or json.dumps(witness, sort_keys=True) < json.dumps(edges[key], sort_keys=True):
            edges[key] = witness

    for record in records:
        accession = record['accession']; current = node(accession)
        source = record['source_id']
        if not isinstance(source, str) or not source:
            raise ValueError('Family links require a source identity')
        current['sources'].add(source)
        current['unresolved'].extend(record.get('unresolved', []))
        samples = record.get('sample_ids', [])
        if any(not re.fullmatch(r'GSM\d+', sample) for sample in samples):
            raise ValueError('Invalid family sample identity')
        current['samples'].update(samples)
        if record.get('source_kind') == 'geo_series_primary' and not record.get('unresolved'):
            current['primary_series'] = True
            current['primary_sample_sets'].append(set(samples))
        groups = {'GSM': samples, 'BioProject': record.get('bioproject_ids', []),
                  'PubMed': record.get('pubmed_ids', [])}
        for kind, identifiers in groups.items():
            for identifier in sorted(set(identifiers)):
                pattern = {'GSM': r'GSM\d+', 'BioProject': r'PRJ(?:NA|EB|DB)\d+', 'PubMed': r'\d+'}[kind]
                if not isinstance(identifier, str) or not re.fullmatch(pattern, identifier):
                    raise ValueError('Invalid family linkage identifier: ' + kind)
                key = (kind, identifier)
                if key in owners:
                    previous, previous_source = owners[key]
                    join(previous, accession, kind, identifier, [source, previous_source])
                else:
                    owners[key] = (accession, source)
        for linked in record.get('linked_accessions', []):
            join(accession, linked, 'GEO_relation', accession + '|' + linked, [source])
    for accession in exposed:
        node(accession)
    components = defaultdict(list)
    for accession in nodes:
        components[root(accession)].append(accession)
    families, studies = [], []
    for members in sorted((sorted(values) for values in components.values()), key=lambda values: values[0]):
        digest = hashlib.sha256((str(seed) + '|' + ';'.join(members)).encode()).hexdigest()
        exposed_members = sorted(set(members) & exposed)
        gaps = []
        for accession in members:
            info = nodes[accession]
            if not info['primary_series']:
                info['unresolved'].append('Primary GEO series links and complete membership are unavailable')
            elif any(samples != info['samples'] for samples in info['primary_sample_sets']):
                info['unresolved'].append('Primary series membership differs from other source memberships')
            gaps.extend({'study': accession, 'reason': reason} for reason in sorted(set(info['unresolved'])))
        proposed = 'held_out' if int(digest, 16) % 5 == 0 else 'development'
        split = 'development' if exposed_members else 'UNRESOLVED' if gaps else proposed
        family_id = 'family:' + digest
        families.append({'family_id': family_id, 'members': members, 'exposed_members': exposed_members,
                         'split': split, 'provisional_hash_assignment': proposed,
                         'source_linkage_status': 'CANNOT VERIFY' if gaps else 'PASS', 'unresolved': gaps})
        for accession in members:
            info = nodes[accession]
            studies.append({'study_id': accession, 'family_id': family_id, 'split': split,
                            'directly_exposed': accession in exposed, 'family_exposed': bool(exposed_members),
                            'source_ids': sorted(info['sources']), 'known_sample_count': len(info['samples']),
                            'source_linkage_status': 'CANNOT VERIFY' if gaps else 'PASS'})
    # One witness per observed study-pair/link type is enough to establish an
    # edge. This is not a count or exhaustive list of shared identifiers.
    return {'schema': 'source-study-families-v1', 'seed': seed, 'families': families,
            'studies': sorted(studies, key=lambda row: row['study_id']),
            'edges': [edges[key] for key in sorted(edges)],
            'edge_witness_scope': 'One canonical witness per connected study pair and link type; not all shared identifiers',
            'exposed_studies': sorted(exposed), 'evaluation_release': 'CLOSED',
            'scope': 'Connected families in supplied source universe; external/unrecorded links are not ruled out',
            'held_out_studies_evaluated': 0}


def verify_family_enrollment(records, manifest, split):
    if manifest.get('schema') != 'source-study-families-v1' or split not in {'development', 'held_out'}:
        raise ValueError('Unrecognized source family manifest or cohort split')
    if split == 'held_out' and manifest.get('evaluation_release') != 'OPEN':
        raise ValueError('Held-out evaluation release is closed')
    members = {}
    for family in manifest['families']:
        expected_exposure = sorted(set(family['members']) & set(manifest['exposed_studies']))
        if family['exposed_members'] != expected_exposure:
            raise ValueError('Family exposure differs from the development registry')
        if family['split'] == 'held_out' and family['exposed_members']:
            raise ValueError('Exposed family cannot be held out')
        expected = hashlib.sha256((str(manifest['seed']) + '|' + ';'.join(sorted(family['members']))).encode()).hexdigest()
        if family['family_id'] != 'family:' + expected:
            raise ValueError('Family identity differs from its complete member list')
        for study in family['members']:
            if study in members:
                raise ValueError('Study is assigned to more than one family')
            members[study] = family
    declared = {row['study_id']: row for row in manifest['studies']}
    if len(declared) != len(manifest['studies']) or set(declared) != set(members):
        raise ValueError('Family and study manifests disagree')
    for edge in manifest['edges']:
        if edge['study_a'] not in members or edge['study_b'] not in members or (
                members[edge['study_a']]['family_id'] != members[edge['study_b']]['family_id']):
            raise ValueError('A source family link crosses declared families')
    for study, family in members.items():
        if declared[study]['family_id'] != family['family_id'] or declared[study]['split'] != family['split']:
            raise ValueError('Family and study split declarations disagree')
    for row in records:
        family = members.get(row['study_id'])
        if family is None or row['family_id'] != family['family_id']:
            raise ValueError('Evaluation family differs from source-linked enrollment')
        if family['split'] != split or declared[row['study_id']]['family_id'] != family['family_id']:
            raise ValueError('Evaluation study is outside the declared cohort split')
        if split == 'held_out' and family['source_linkage_status'] != 'PASS':
            raise ValueError('Unresolved source linkage cannot establish a held-out family')


def fraction(n, d):
    return n/d if d else None


def metrics(rows):
    accepted = [r for r in rows if r['accepted']]
    known = [r for r in rows if r['reference_status'] != 'CANNOT VERIFY']
    accepted_known = [r for r in accepted if r['final_reference_status'] != 'CANNOT VERIFY']
    errors = [r for r in rows if r['reference_status'] == 'FAIL']
    correct = [r for r in rows if r['reference_status'] == 'PASS']
    bad_accepted = sum(r['final_reference_status']=='FAIL' for r in accepted)
    unknown_accepted = sum(r['final_reference_status']=='CANNOT VERIFY' for r in accepted)
    repair = [r for r in rows if r['repair_proposed']]
    verified_repair = [r for r in repair if r['repair_status']!='CANNOT VERIFY']
    costs = [r['cost_usd'] for r in rows if r['cost_usd'] is not None]
    runtimes = [r['runtime_seconds'] for r in rows if r['runtime_seconds'] is not None]
    estimates = [r['estimated_cost_usd'] for r in rows if r.get('estimated_cost_usd') is not None]
    result = {'attempted_studies': len(rows), 'accepted_studies': len(accepted),
            'reference_resolved_studies': len(known),
            'accepted_reference_resolved_studies': len(accepted_known),
            'accepted_material_errors': bad_accepted,
            'accepted_reference_unresolved_studies': unknown_accepted,
            'accepted_result_error_rate': fraction(bad_accepted,len(accepted_known)),
            'accepted_error_lower_bound': fraction(bad_accepted,len(accepted)),
            'accepted_error_upper_bound': fraction(bad_accepted+unknown_accepted,len(accepted)),
            'coverage': fraction(len(accepted),len(rows)),
            'missed_material_errors': sum(not r['flagged_material_error'] for r in errors),
            'material_error_denominator': len(errors),
            'false_alarms': sum(r['flagged_material_error'] for r in correct),
            'no_material_error_denominator': len(correct),
            'proposed_repairs': len(repair), 'verified_repairs': len(verified_repair),
            'verified_correct_repairs': sum(r['repair_status']=='PASS' for r in verified_repair),
            'repair_correctness': fraction(sum(r['repair_status']=='PASS' for r in verified_repair),len(verified_repair)),
            'rejected_studies': sum(r['decision_status']=='FAIL' for r in rows),
            'unresolved_studies': sum(r['decision_status']=='CANNOT VERIFY' for r in rows),
            'unresolved_repairs': sum(r['repair_status']=='CANNOT VERIFY' for r in repair),
            'introduced_accepted_errors': sum(r['reference_status']=='PASS' and r['final_reference_status']=='FAIL' for r in accepted),
            'total_runtime_seconds': sum(runtimes) if len(runtimes) == len(rows) else None,
            'known_runtime_seconds': sum(runtimes), 'runtime_known_studies': len(runtimes),
            'known_cost_usd': sum(costs), 'cost_known_studies': len(costs),
            'estimated_known_cost_usd': sum(estimates), 'estimated_cost_known_studies': len(estimates)}
    if any(row.get('reference_role') == 'reference_unavailable' for row in rows):
        result['reference_unavailable_studies'] = sum(row.get('reference_role') == 'reference_unavailable' for row in rows)
        result['shared_prerequisite_unavailable_studies'] = sum(
            row.get('reference_availability') == 'UNAVAILABLE_SHARED_PREREQUISITE' for row in rows)
    if any(row.get('reference_role') == 'reference_incomplete' for row in rows):
        result['reference_incomplete_studies'] = sum(row.get('reference_role') == 'reference_incomplete' for row in rows)
        result['prospective_preparation_unavailable_studies'] = sum(
            row.get('preparation_status') == 'UNAVAILABLE_PROSPECTIVE_PREPARATION' for row in rows)
    return result


def evaluate(records=None, resamples=2000, seed=20260910, *, reference_freezes=None, family_manifest=None, split='development'):
    """Evaluate only frozen, source-linked independent reference submissions."""
    from .agent_workflow_2026_09_10 import sha
    if not isinstance(reference_freezes, list) or not reference_freezes:
        raise ValueError('Evaluation requires pinned frozen reference artifacts; provenance strings alone are insufficient')
    if split not in {'development', 'held_out'}: raise ValueError('Unknown evaluation split')
    derived, seen, archives = [], set(), set()
    for item in reference_freezes:
        path = str(Path(item['path']).resolve())
        if path in seen or sha(path) != item.get('sha256'): raise ValueError('Duplicate or changed frozen reference artifact')
        seen.add(path); derived.extend(read_reference_records(path))
        frozen = json.loads(Path(path).read_text())
        if frozen.get('schema') == 'clockbase-reference-unavailable-v1':
            archives.add(str(Path(frozen['source_identities']['archive']).resolve()))
        elif frozen.get('schema') == 'clockbase-reference-incomplete-v1':
            archives.add(str(Path(frozen['state']['archive']).resolve()))
        else:
            numerical = Path(frozen['source_identities']['numerical_archive'])
            archives.add(str(Path(json.loads((numerical / 'config.json').read_text())['archive']).resolve()))
    if len(archives) != 1:
        raise ValueError('A verified comparison requires one common enrolled method archive; independently selected archives cannot define its denominator')
    archive = Path(next(iter(archives)))
    context = reference_archive_context(archive)
    if context['split'] != split:
        raise ValueError('Reference artifacts cannot supply a different development/held-out evaluation split')
    if context.get('plan') is not None:
        analysis = context['plan'].get('analysis', {})
        if (set(analysis) != {'primary_endpoint', 'study_weighting', 'uncertainty', 'resamples', 'seed'}
                or analysis.get('primary_endpoint') != 'accepted_result_error_rate'
                or analysis.get('study_weighting') != 'equal'
                or analysis.get('uncertainty') != 'paired_family_cluster_percentile_bootstrap'
                or type(analysis.get('resamples')) is not int or analysis['resamples'] < 1
                or type(analysis.get('seed')) is not int
                or type(resamples) is not int or type(seed) is not int
                or analysis['resamples'] != resamples or analysis['seed'] != seed):
            raise ValueError('Evaluation analysis or bootstrap settings differ from the prospective frozen plan')
    enrolled_cases = reference_cases(archive)
    enrollment = context.get('enrollment') or {'schema': 'clockbase-comparison-enrollment-v1', 'split': split,
                  'methods': list(COMPARISON_METHODS),
                  'studies': [{'study_id': case['study_id'], 'family_id': case['family_id']} for case in enrolled_cases]}
    if context.get('enrollment') is not None:
        if not isinstance(context.get('family_manifest'), dict):
            raise ValueError('Prospective evaluation requires its jointly frozen source-family manifest')
        if family_manifest is not None and family_manifest != context['family_manifest']:
            raise ValueError('Evaluation family manifest differs from the prospective joint freeze')
        family_manifest = context['family_manifest']
    validate_attempted_universe(derived, enrollment)
    if records is not None:
        ordered = lambda values: sorted(values, key=lambda row: (row['study_id'], row['method']))
        if ordered(records) != ordered(derived): raise ValueError('Supplied labels or decisions differ from frozen reference records')
    result = summarize_records(derived, resamples, seed, family_manifest=family_manifest, split=split,
                               enrollment=enrollment)
    result.update(reference_artifacts_verified=True, reference_freezes=reference_freezes,
                  enrollment_archive=str(archive), enrollment=enrollment,
                  enrollment_archive_config_sha256=sha(archive / 'config.json'),
                  evidence_scope='Frozen independent assessor submissions with checked source links; interpretation and assessor independence are declarations, not authenticated by software')
    unavailable = {row['study_id'] for row in derived if row.get('reference_role') == 'reference_unavailable'}
    if unavailable:
        result.update(reference_unavailable_studies=sorted(unavailable),
            evidence_scope='Frozen assessor submissions and explicitly unavailable shared-prerequisite references; unavailable records contain no scientific judgment or invented adjudicator')
    incomplete = {row['study_id'] for row in derived if row.get('reference_role') == 'reference_incomplete'}
    if incomplete:
        result.update(reference_incomplete_studies=sorted(incomplete),
            evidence_scope='Frozen independent judgments where available, plus pinned administrative missingness; method acceptance and all enrolled attempts are retained without inventing scientific labels')
    upstream = {}
    for row in derived:
        if 'shared_upstream_resources' not in row: continue
        if row['study_id'] in upstream and upstream[row['study_id']] != row['shared_upstream_resources']:
            raise ValueError('Shared producer resources differ between method slots')
        upstream[row['study_id']] = row['shared_upstream_resources']
    if upstream:
        result['shared_upstream_resources_by_study'] = upstream
        result['shared_upstream_resource_scope'] = 'Count each study once; excluded from all incremental method-cost totals'
    return result


def summarize_records(records, resamples=2000, seed=20260910, *, family_manifest=None, split='development', enrollment=None):
    """Descriptive label arithmetic; direct calls do not verify reference evidence."""
    if not records or resamples<1:
        raise ValueError('Nonempty adjudicated study records and positive bootstrap resamples required')
    if family_manifest is not None:
        verify_family_enrollment(records, family_manifest, split)
    elif split == 'held_out':
        raise ValueError('Held-out evaluation requires source-linked family enrollment')
    if enrollment is not None:
        if enrollment.get('split') != split: raise ValueError('Evaluation split differs from enrollment')
        enrolled_counts = validate_attempted_universe(records, enrollment)
    elif split == 'held_out':
        raise ValueError('Held-out evaluation requires the frozen attempted study-by-method roster')
    else:
        enrolled_counts = None
    methods = defaultdict(list)
    identities = {}
    seen = set()
    for row in records:
        for field in ('study_id','family_id','method','reference_source'):
            if not isinstance(row.get(field),str) or not row[field].strip():
                raise ValueError('Missing study/family/method or independent reference provenance')
        if row.get('reference_role') == 'reference_unavailable':
            if (row.get('reference_adjudicator') is not None
                    or row.get('reference_availability') != 'UNAVAILABLE_SHARED_PREREQUISITE'
                    or any(row.get(key) != 'CANNOT VERIFY' for key in
                           ('reference_status', 'final_reference_status', 'repair_status', 'decision_status'))
                    or any(row.get(key) is not False for key in ('accepted', 'flagged_material_error', 'repair_proposed'))):
                raise ValueError('Unavailable reference cannot invent an assessor or scientific judgment, or erase an executed decision')
        elif row.get('reference_role') == 'reference_incomplete':
            known = row.get('known_reference_sources')
            adjudicator = row.get('reference_adjudicator')
            if (row.get('reference_availability') != 'MISSING_OR_PARTIAL_REFERENCE'
                    or row.get('reference_missing_phase') not in {'SOURCE_DESIGN', 'NUMERICAL_REVIEW'}
                    or not isinstance(known, list) or any(not isinstance(item, str) or not item.strip() for item in known)
                    or bool(known) != (isinstance(adjudicator, str) and bool(adjudicator.strip()))
                    or not known and (adjudicator is not None or any(row.get(key) != 'CANNOT VERIFY'
                            for key in ('reference_status', 'final_reference_status', 'repair_status')))):
                raise ValueError('Incomplete reference needs pinned missingness and cannot invent an assessor or known label')
        elif row.get('reference_role')!='independent_source_adjudication':
            raise ValueError('Model consensus or prior pipeline labels cannot supply reference judgments')
        elif not isinstance(row.get('reference_adjudicator'), str) or not row['reference_adjudicator'].strip():
            raise ValueError('Missing independent reference adjudicator')
        for field in ('accepted','flagged_material_error','repair_proposed'):
            if type(row.get(field)) is not bool:
                raise ValueError('Decision fields must be explicit booleans')
        for field in ('reference_status','final_reference_status','repair_status','decision_status'):
            if row.get(field) not in ('PASS','FAIL','CANNOT VERIFY'):
                raise ValueError('Unknown reference or repair status')
        if row['accepted'] != (row['decision_status']=='PASS'):
            raise ValueError('Acceptance must agree with the method decision, not the reference')
        if not row['repair_proposed'] and row['repair_status']!='CANNOT VERIFY':
            raise ValueError('Unproposed repair cannot be counted as verified')
        if row['accepted'] and (not isinstance(row.get('final_reference_source'),str) or not row['final_reference_source'].strip()):
            raise ValueError('Accepted outputs need their own independent reference provenance')
        for field in ('runtime_seconds','cost_usd','estimated_cost_usd'):
            if field in ('runtime_seconds', 'cost_usd') and field not in row:
                raise ValueError('Time/cost must be explicitly measured or unknown')
            value=row.get(field)
            if value is None:
                continue
            if type(value) not in (int,float) or not math.isfinite(value) or value<0:
                raise ValueError('Invalid time/cost measurement')
        key=(row['method'],row['study_id'])
        if key in seen:
            raise ValueError('Duplicate method-study row; clocks cannot be counted as studies')
        seen.add(key)
        identity=(row['family_id'],row['reference_status'],row['reference_source'])
        if row['study_id'] in identities and identities[row['study_id']]!=identity:
            raise ValueError('Study reference or family differs between methods')
        identities[row['study_id']]=identity
        methods[row['method']].append(row)
    universe={r['study_id'] for r in records}
    if any({r['study_id'] for r in rows}!=universe for rows in methods.values()):
        raise ValueError('Methods must cover the same attempted study universe, including failures')
    families=sorted({r['family_id'] for r in records})
    scores={method:metrics(rows) for method,rows in methods.items()}
    generator=random.Random(seed)
    bootstrap={method:{key:[] for key in ('accepted_result_error_rate','coverage')} for method in methods}
    pairs=list(combinations(sorted(methods),2))
    differences={pair:{key:[] for key in ('accepted_result_error_rate','coverage')} for pair in pairs}
    # All methods receive the same family resample; sibling studies remain together.
    for _ in range(resamples if len(families)>1 else 0):
        selected=generator.choices(families,k=len(families))
        sampled={}
        for method,rows in methods.items():
            values=metrics([row for f in selected for row in rows if row['family_id']==f])
            sampled[method]=values
            for key in bootstrap[method]:
                if values[key] is not None:
                    bootstrap[method][key].append(values[key])
        for a,b in pairs:
            for key in differences[a,b]:
                if sampled[a][key] is not None and sampled[b][key] is not None:
                    differences[a,b][key].append(sampled[b][key]-sampled[a][key])
    def interval(values):
        ordered=sorted(values)
        # A boundary-only bootstrap cannot establish absence of population errors.
        if len(ordered)<.95*resamples or len(set(ordered))<2:
            return None
        return [ordered[int(.025*(len(ordered)-1))],ordered[int(.975*(len(ordered)-1))]]
    for method in methods:
        scores[method]['family_count']=len(families)
        for key,values in bootstrap[method].items():
            # Undefined no-acceptance resamples are counted, never replaced with zero.
            scores[method][key+'_bootstrap_defined']=len(values)
            scores[method][key+'_ci95']=interval(values)
    contrasts=[]
    for a,b in pairs:
        comparison={'method_a':a,'method_b':b,'direction':'method_b minus method_a'}
        for key,values in differences[a,b].items():
            comparison[key+'_difference']=(None if scores[a][key] is None or scores[b][key] is None else scores[b][key]-scores[a][key])
            comparison[key+'_difference_ci95']=interval(values)
            comparison[key+'_bootstrap_defined']=len(values)
        contrasts.append(comparison)
    return {'unit':'study', 'uncertainty_unit':'study family', 'seed':seed,
            'weighting': 'Each enrolled study has equal weight; related studies are resampled together',
            'attempted_universe_verified': enrollment is not None, 'enrollment_counts': enrolled_counts,
            'reference_artifacts_verified': False,
            'evidence_scope': 'Descriptive arithmetic on supplied labels; independent reference evidence has not been checked',
            'family_enrollment_verified': family_manifest is not None, 'split': split,
            'resamples':resamples, 'methods':scores,'paired_method_comparisons':contrasts,
            'ci_method':'Paired family-cluster percentile bootstrap; absent with insufficient defined draws or no variation',
            'interpretation':'Error is conditional on reference-resolved acceptance; report coverage and unknown-reference bounds together. A missing interval does not imply certainty.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--records',type=Path)
    mode.add_argument('--compare-plan',type=Path)
    mode.add_argument('--replay-comparison',type=Path,help='Recheck preserved development reviews without model or clock calls')
    mode.add_argument('--reference-archive',type=Path,help='Export source-only reference design packets from a development archive')
    mode.add_argument('--freeze-reference-design',type=Path,help='Validate and freeze a completed independent design judgment')
    mode.add_argument('--reference-outcomes',type=Path,help='Release anonymized numerical outputs after source-design adjudication')
    mode.add_argument('--freeze-reference-outcomes',type=Path,help='Validate and freeze a completed numerical reference judgment')
    mode.add_argument('--freeze-reference-unavailable',type=Path,help='Method archive with an enrolled shared-prerequisite failure; no scientific judgment is created')
    mode.add_argument('--freeze-reference-incomplete',type=Path,help='Retain completed method decisions and prior judgments when independent reference assessment is missing or partial')
    mode.add_argument('--evaluate-references',type=Path,help='JSON list of frozen reference paths and SHA256 identities')
    mode.add_argument('--prepare-source-readiness',metavar='BATCH_ID',help='Project one fixed source-only readiness batch after parent freeze')
    mode.add_argument('--freeze-source-readiness',type=Path,help='Preserve and mechanically validate a completed source-readiness submission')
    mode.add_argument('--report-source-readiness',type=Path,help='JSON list of all pinned readiness versions in supersession order; no accuracy evaluation')
    parser.add_argument('--readiness-release',type=Path,help='Parent-frozen source-readiness protocol/runtime/code binding')
    parser.add_argument('--readiness-release-sha256',help='Expected parent release SHA256')
    parser.add_argument('--readiness-batch',help='Fixed batch for --freeze-source-readiness')
    parser.add_argument('--readiness-packet-sha256',help='Expected --source-design packet SHA256')
    parser.add_argument('--readiness-version',help='New named submission version')
    parser.add_argument('--readiness-supersedes',type=Path,help='JSON path/SHA256 binding of the previous readiness version to retain')
    parser.add_argument('--source-design',type=Path,help='Source-design packet named by a submitted reference judgment')
    parser.add_argument('--study',help='Enrolled study for an unavailable or incomplete reference freeze')
    parser.add_argument('--missingness-receipt',type=Path,help='Pinned administrative receipt required for --freeze-reference-incomplete')
    parser.add_argument('--prior-frozen-outcomes',type=Path,help='Earlier numerical reference judgment to retain in an incomplete reference freeze')
    parser.add_argument('--frozen-design',type=Path,help='Frozen reference design required for numerical review')
    parser.add_argument('--numerical-archive',type=Path,help='Anonymous numerical packet archive named by a reference judgment')
    parser.add_argument('--reference-freezes',type=Path,help='Pinned reference artifacts required when checking supplied --records')
    parser.add_argument('--prepare-comparison',action='store_true')
    parser.add_argument('--family-manifest',type=Path)
    parser.add_argument('--split',choices=['development','held_out'],default='development')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    try:
        if args.prepare_source_readiness or args.freeze_source_readiness or args.report_source_readiness:
            if args.readiness_release is None or args.readiness_release_sha256 is None:
                raise ValueError('Source readiness requires a parent-frozen release and its expected SHA256')
            if args.prepare_source_readiness:
                result=prepare_source_readiness_batch(args.readiness_release,args.readiness_release_sha256,
                    args.prepare_source_readiness,args.output)
            elif args.freeze_source_readiness:
                if not all((args.source_design,args.readiness_packet_sha256,args.readiness_batch,args.study,args.readiness_version)):
                    raise ValueError('Source readiness freeze requires the exact batch, study, packet hash and named version')
                result=freeze_source_readiness_submission(args.readiness_release,args.readiness_release_sha256,
                    args.readiness_batch,args.study,{'path':str(args.source_design),'sha256':args.readiness_packet_sha256},
                    args.freeze_source_readiness,args.readiness_version,args.output,
                    supersedes=_readiness_json(args.readiness_supersedes) if args.readiness_supersedes else None)
            else:
                result=report_source_readiness(args.readiness_release,args.readiness_release_sha256,
                    _readiness_json(args.report_source_readiness),args.output)
            print(json.dumps(result,indent=2))
        elif args.reference_archive:
            print(json.dumps(export_reference_designs(args.reference_archive,args.output),indent=2))
        elif args.freeze_reference_design:
            if args.source_design is None:raise ValueError('Reference freeze requires --source-design')
            print(json.dumps(freeze_reference_design(args.source_design,args.freeze_reference_design,args.output),indent=2))
        elif args.reference_outcomes:
            if args.frozen_design is None:raise ValueError('Numerical reference release requires --frozen-design')
            print(json.dumps(export_reference_outcomes(args.reference_outcomes,args.frozen_design,args.output),indent=2))
        elif args.freeze_reference_outcomes:
            if args.numerical_archive is None:raise ValueError('Numerical reference freeze requires --numerical-archive')
            print(json.dumps(freeze_reference_outcomes(args.numerical_archive,args.freeze_reference_outcomes,args.output),indent=2))
        elif args.freeze_reference_unavailable:
            if not args.study: raise ValueError('Unavailable reference freeze requires an enrolled --study')
            print(json.dumps(freeze_reference_unavailable(args.freeze_reference_unavailable,args.study,args.output),indent=2))
        elif args.freeze_reference_incomplete:
            if not args.study or args.missingness_receipt is None:
                raise ValueError('Incomplete reference freeze requires --study and --missingness-receipt')
            print(json.dumps(freeze_reference_incomplete(args.freeze_reference_incomplete,args.study,
                args.missingness_receipt,args.output,frozen_design=args.frozen_design,
                frozen_outcomes=args.prior_frozen_outcomes),indent=2))
        elif args.replay_comparison:
            result=replay_development_comparison(args.replay_comparison,args.output)
            print(json.dumps(result,indent=2))
        elif args.compare_plan:
            result=run_comparison_plan(args.compare_plan,args.output,prepare_only=args.prepare_comparison)
            print(json.dumps({key:value for key,value in result.items() if key not in {'records','api_sha256'}}))
        else:
            if args.family_manifest is None:raise ValueError('Evaluation requires --family-manifest')
            references = args.evaluate_references or args.reference_freezes
            if references is None:raise ValueError('Evaluation requires frozen reference artifacts')
            opener=gzip.open if args.family_manifest.suffix=='.gz' else open
            with opener(args.family_manifest,'rt') as source:families=json.load(source)
            result=evaluate(json.loads(args.records.read_text()) if args.records else None,
                            reference_freezes=json.loads(references.read_text()),family_manifest=families,split=args.split)
            with args.output.open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False)
    except Exception as error:
        raise SystemExit(f'Evaluation stopped: {type(error).__name__}: {error}')
