"""Extract a complete historical trace without executing pickle objects.

This inspects preserved records. It does not execute archived code, rerun clocks,
or establish correctness of the recorded analysis.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import pickle
import re
import threading


class SourceOnlyUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        raise pickle.UnpicklingError('Executable pickle object refused: ' + module + '.' + name)

    def persistent_load(self, pid):
        raise pickle.UnpicklingError('Persistent pickle reference refused')


CANONICAL_STUDY = re.compile(r'GSE[0-9]+')
HEX_SHA256 = re.compile(r'[0-9a-f]{64}')
MAX_ARCHIVE_BYTES = 2 * 1024**3
_SECRET = re.compile(
    r'(?i)(?:Authorization\s*:\s*Bearer\s+[^\s,;]+|'
    r'\bsk-(?:proj-)?[A-Za-z0-9_-]{8,}|\bgh[pousr]_[A-Za-z0-9]{20,}|'
    r'\b(?:api[_ -]?key|access[_ -]?token|secret)\s*[:=]\s*[^\s,;]+)')
_PRIVATE_PATH = re.compile(r"(?<![A-Za-z0-9_])/(?:Entropy|Lasagna|home/kying|Users|tmp|var/tmp)/[^\s\"'<>]*")
_EMAIL = re.compile(r'(?i)\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b')
_CONTACT_LINE = re.compile(
    r'(?im)^.*(?:!Series_contact_|!Sample_contact_|'
    r'\bcontact[_ ](?:name|email|phone|address)\s*[:=\t]).*$')


def _file_sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _canonical_trace_value(value):
    """Convert inert pickle primitives to a deterministic JSON representation."""
    if isinstance(value, tuple):
        return [_canonical_trace_value(item) for item in value]
    if isinstance(value, list):
        return [_canonical_trace_value(item) for item in value]
    if isinstance(value, dict):
        if any(not isinstance(key, str) for key in value):
            raise ValueError('Trace dictionaries require string keys')
        return {key: _canonical_trace_value(item) for key, item in value.items()}
    if value is None or type(value) in (str, bool, int):
        return value
    if type(value) is float and math.isfinite(value):
        return value
    raise ValueError('Unsupported trace value type: ' + type(value).__name__)


def _sanitized_trace_text(text):
    """Redact credentials, private paths and contact fields with hash mapping."""
    if not isinstance(text, str):
        raise ValueError('Trace evidence text must be a string')
    spans = []
    for kind, pattern in (('credential', _SECRET), ('private_path', _PRIVATE_PATH),
                          ('email', _EMAIL), ('contact_field', _CONTACT_LINE)):
        spans.extend((match.start(), match.end(), kind) for match in pattern.finditer(text))
    merged = []
    for start, end, kind in sorted(spans):
        if merged and start < merged[-1][1]:
            old_start, old_end, kinds = merged[-1]
            merged[-1] = (old_start, max(old_end, end), kinds | {kind})
        else:
            merged.append((start, end, {kind}))
    output, mapping, cursor = [], [], 0
    for index, (start, end, kinds) in enumerate(merged):
        original = text[start:end]
        digest = hashlib.sha256(original.encode()).hexdigest()
        placeholder = '[REDACTED_' + '_'.join(sorted(kinds)).upper() + f'_{index:03d}_{digest[:12]}]'
        output.extend((text[cursor:start], placeholder))
        mapping.append({'placeholder': placeholder, 'kinds': sorted(kinds),
                        'original_sha256': digest, 'original_characters': len(original)})
        cursor = end
    output.append(text[cursor:])
    safe = ''.join(output)
    if _SECRET.search(safe) or _PRIVATE_PATH.search(safe) or _EMAIL.search(safe) or _CONTACT_LINE.search(safe):
        raise ValueError('Trace redaction left recognized private content')
    return safe, mapping


class TraceArchiveLoader:
    """Thread-safe, source-only cache for post-design historical trace evidence.

    The loader never returns an archive object. Callers receive only a bounded,
    sanitized evidence record after supplying the already-frozen source-design
    hash. A request mismatch withholds all message text.
    """
    def __init__(self, archives, *, max_archive_bytes=MAX_ARCHIVE_BYTES,
                 max_evidence_chars=8 * 1024**2):
        if (not isinstance(archives, dict) or not 1 <= len(archives) <= 3 or
                type(max_archive_bytes) is not int or not 0 < max_archive_bytes <= MAX_ARCHIVE_BYTES or
                type(max_evidence_chars) is not int or max_evidence_chars <= 0):
            raise ValueError('Require one to three pinned archives and positive bounded size limits')
        records, paths = {}, set()
        for modality, record in archives.items():
            if (modality not in {'mouse_rna', 'human_rna', 'human_dnam'} or
                    not isinstance(record, dict) or set(record) != {'path', 'sha256'} or
                    not isinstance(record['sha256'], str) or not HEX_SHA256.fullmatch(record['sha256'])):
                raise ValueError('Archive declarations require a known modality, path and lowercase SHA-256')
            path = Path(record['path']).resolve()
            if path in paths:
                raise ValueError('Each modality requires a distinct trace archive')
            paths.add(path)
            records[modality] = {'path': path, 'sha256': record['sha256']}
        self._records = records
        self._max_archive_bytes = max_archive_bytes
        self._max_evidence_chars = max_evidence_chars
        self._cache = {}
        self._failures = {}
        self._load_counts = {modality: 0 for modality in records}
        self._lock = threading.RLock()

    def _load(self, modality):
        if modality not in self._records:
            raise KeyError('No pinned trace archive for modality: ' + str(modality))
        with self._lock:
            if modality in self._failures:
                raise ValueError(self._failures[modality])
            if modality in self._cache:
                return self._cache[modality]
            record = self._records[modality]
            try:
                stat = record['path'].stat()
                if not 0 < stat.st_size <= self._max_archive_bytes:
                    raise ValueError('Trace archive is empty or exceeds the frozen byte limit')
                before = _file_sha256(record['path'])
                if before != record['sha256']:
                    raise ValueError('Trace archive SHA-256 differs from its pinned declaration')
                with record['path'].open('rb') as handle:
                    archive = SourceOnlyUnpickler(handle).load()
                    if handle.read(1):
                        raise ValueError('Unparsed trailing content in trace archive')
                if _file_sha256(record['path']) != before:
                    raise ValueError('Trace archive changed while being loaded')
                if (not isinstance(archive, dict) or
                        any(not isinstance(key, str) or not CANONICAL_STUDY.fullmatch(key)
                            for key in archive)):
                    raise ValueError('Trace archive requires canonical ASCII GSE keys')
                self._cache[modality] = (archive, stat.st_size)
                self._load_counts[modality] += 1
                return self._cache[modality]
            except (OSError, ValueError, pickle.UnpicklingError, EOFError) as error:
                message = type(error).__name__ + ': ' + str(error)
                self._failures[modality] = message
                raise ValueError(message) from error

    def load_counts(self):
        """Return diagnostic counts without exposing loaded archive objects."""
        with self._lock:
            return dict(self._load_counts)

    def verify_sources_unchanged(self):
        """Recheck every declared archive after a bounded bulk review."""
        ledger = []
        with self._lock:
            for modality, record in sorted(self._records.items()):
                stat = record['path'].stat()
                if not 0 < stat.st_size <= self._max_archive_bytes:
                    raise ValueError('Trace archive is empty or exceeds the frozen byte limit')
                digest = _file_sha256(record['path'])
                if digest != record['sha256']:
                    raise ValueError('Trace archive changed after its pinned declaration')
                ledger.append({'modality': modality, 'source_sha256': digest,
                               'source_bytes': stat.st_size,
                               'loaded_once': self._load_counts[modality] == 1})
        return ledger

    def study_evidence(self, modality, study_id, *, frozen_source_design_sha256,
                       expected_prompt=None):
        if (not isinstance(study_id, str) or not CANONICAL_STUDY.fullmatch(study_id) or
                not isinstance(frozen_source_design_sha256, str) or
                not HEX_SHA256.fullmatch(frozen_source_design_sha256)):
            raise ValueError('Canonical study ID and frozen source-design SHA-256 are required')
        archive, source_bytes = self._load(modality)
        if study_id not in archive:
            return {'schema': 'clockbase-sanitized-trace-evidence-v1', 'study_id': study_id,
                    'modality': modality, 'availability': 'NOT_FOUND',
                    'source_sha256': self._records[modality]['sha256'],
                    'source_bytes': source_bytes, 'source_design_sha256': frozen_source_design_sha256,
                    'eligible_for_post_design_context': False, 'evidence': [],
                    'scientific_status': 'CANNOT VERIFY'}
        entry = archive[study_id]
        if not isinstance(entry, dict) or entry.get('study_id') != study_id:
            raise ValueError('Trace archive key and declared study identity differ')
        result = entry.get('result')
        messages = result[0] if isinstance(result, tuple) and result else result
        binding = verify_trace_request(messages, study_id, expected_prompt)
        base = {'schema': 'clockbase-sanitized-trace-evidence-v1', 'study_id': study_id,
                'modality': modality, 'availability': 'FOUND',
                'source_sha256': self._records[modality]['sha256'],
                'source_bytes': source_bytes, 'source_entry_count': len(archive),
                'source_design_sha256': frozen_source_design_sha256,
                'request_binding': binding, 'archive_status': entry.get('status'),
                'entry_canonical_sha256': hashlib.sha256(json.dumps(
                    _canonical_trace_value(entry), sort_keys=True, separators=(',', ':'),
                    ensure_ascii=False, allow_nan=False).encode()).hexdigest(),
                'scientific_status': 'CANNOT VERIFY',
                'scope': 'Sanitized historical trace evidence after source-design freeze; request identity does not establish input-to-result lineage or scientific correctness'}
        if binding['status'] == 'FAIL' or entry.get('status') != 'success':
            return {**base, 'eligible_for_post_design_context': False, 'evidence': [],
                    'withheld_reason': ('Trace request identity is FAIL' if binding['status'] == 'FAIL'
                                        else 'Archived trace status is not success')}
        if not isinstance(messages, list) or any(not isinstance(message, str) for message in messages):
            raise ValueError('Unsupported historical message representation')
        terminal = result[1] if isinstance(result, tuple) and len(result) > 1 else None
        values = [('message', index, message) for index, message in enumerate(messages)]
        if isinstance(terminal, str):
            values.append(('terminal_reply', 0, terminal))
        if sum(len(value) for _, _, value in values) > self._max_evidence_chars:
            raise ValueError('Sanitized trace evidence exceeds the frozen character limit')
        evidence = []
        for ordinal, (kind, index, original) in enumerate(values):
            safe, redactions = _sanitized_trace_text(original)
            evidence.append({'evidence_id': f'T{ordinal:05d}', 'kind': kind, 'message_index': index,
                'safe_text': safe, 'safe_sha256': hashlib.sha256(safe.encode()).hexdigest(),
                'original_sha256': hashlib.sha256(original.encode()).hexdigest(),
                'original_characters': len(original), 'redactions': redactions,
                'mentioned_studies': sorted(set(re.findall(r'\bGSE[0-9]+\b', safe)))})
        identity_scope = ('REQUEST_BOUND_CONTEXT' if binding['status'] == 'PASS'
                          else 'UNBOUND_CONTEXT')
        return {**base, 'eligible_for_post_design_context': True,
                'identity_scope': identity_scope,
                'identity_limit': ('The explicit initial request matches this study.'
                                   if binding['status'] == 'PASS' else
                                   'No contradictory request was found, but the initial request identity is absent; this trace cannot alone support a resolved source, grouping or scientific claim.'),
                'evidence': evidence,
                'evidence_count': len(evidence), 'evidence_characters': sum(len(row['safe_text']) for row in evidence),
                'raw_archive_object_exposed': False, 'archived_code_executed': False}


def verify_trace_request(messages, study, expected_prompt=None):
    """Check explicit initial requests, not incidental mentions of related studies."""
    if not isinstance(study, str) or not re.fullmatch(r'GSE[0-9]+', study):
        return {'status': 'FAIL', 'expected_study': study,
                'issues': ['Expected study must be a canonical ASCII GEO Series accession']}
    if not isinstance(messages, list) or any(not isinstance(message, str) for message in messages):
        return {'status': 'CANNOT VERIFY', 'expected_study': study, 'issues': ['Unsupported trace message representation']}
    human = []
    for index, message in enumerate(messages):
        match = re.fullmatch(r'=+\s*Human Message\s*=+\r?\n\r?\n([\s\S]*)', message)
        if match: human.append((index, match.group(1)))
    issues = []
    if len(human) != 1 or human[0][0] != 0:
        issues.append('Trace requires one initial Human Message belonging to this study')
    requests = set()
    for _, body in human:
        requests.update(value.strip() for value in re.findall(r'^Study identifier:([^\r\n]+)', body, re.MULTILINE))
        requests.update(re.findall(r'^Analyze the aging effects of conditions from study (GSE[0-9]+)\.', body, re.MULTILINE))
    if requests and requests != {study}: issues.append('Initial requested study differs from the trace record')
    if expected_prompt is not None and (len(human) != 1 or human[0][1].strip() != expected_prompt.strip()):
        issues.append('Initial prompt differs from the source-bound pipeline request')
    return {'status': 'FAIL' if issues else 'PASS' if requests == {study} else 'CANNOT VERIFY',
            'expected_study': study, 'declared_request_studies': sorted(requests), 'human_messages': len(human),
            'complete_prompt_checked': expected_prompt is not None, 'issues': issues or ([] if requests else ['No explicit study identity in initial request']),
            'scope': 'Trace/request identity; does not establish historical input-to-result lineage'}


def extract_trace(source, study, output, max_bytes=2 * 1024**3):
    from audit_system.agent_workflow_2026_09_10 import sha, write_new
    from src.pipeline_runtime_2026_09_10 import StageStore
    source, output = Path(source).resolve(), Path(output).resolve()
    if source.is_relative_to(output): raise ValueError('Trace output must be separate from its source file')
    stat = source.stat()
    if not 0 < stat.st_size <= max_bytes: raise ValueError('Trace source exceeds declared size limit or is empty')
    source_hash = sha(source)
    config = {'source_path': str(source), 'source_sha256': source_hash, 'source_bytes': stat.st_size,
              'source_mtime_utc': datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
              'study_id': study, 'max_bytes': max_bytes, 'implementation_sha256': sha(__file__),
              'scope': 'Complete selected historical trace extraction; no code or clock execution'}
    store = StageStore(output, config)

    def extract(directory):
        with source.open('rb') as handle:
            traces = SourceOnlyUnpickler(handle).load()
            if handle.read(1): raise ValueError('Unparsed trailing content in trace source')
        if not isinstance(traces, dict) or study not in traces: raise ValueError('Requested study is absent from trace source')
        entry = traces[study]
        if not isinstance(entry, dict) or entry.get('study_id') != study:
            raise ValueError('Trace key and declared study identity differ')
        result = entry.get('result')
        messages = result[0] if isinstance(result, tuple) and result else result
        if not isinstance(messages, list) or any(not isinstance(message, str) for message in messages):
            raise ValueError('Unsupported historical message representation; no implicit string conversion')
        tuple_paths = []
        def normalize(value, pointer=''):
            if isinstance(value, tuple):
                tuple_paths.append(pointer)
                return [normalize(item, pointer + '/' + str(index)) for index, item in enumerate(value)]
            if isinstance(value, list):
                return [normalize(item, pointer + '/' + str(index)) for index, item in enumerate(value)]
            if isinstance(value, dict):
                if any(not isinstance(key, str) for key in value): raise ValueError('Trace JSON requires string dictionary keys')
                return {key: normalize(item, pointer + '/' + key.replace('~', '~0').replace('/', '~1')) for key, item in value.items()}
            if value is None or type(value) in (str, bool, int): return value
            if type(value) is float and math.isfinite(value): return value
            raise ValueError('Unsupported trace value type: ' + type(value).__name__)
        normalized = normalize(entry)
        serialized = json.dumps(normalized, ensure_ascii=False, allow_nan=False)
        # Abort before saving selected content. Never echo the detected string.
        if re.search(r'\bsk-(?:proj-)?[A-Za-z0-9_-]{16,}|\bgh[pousr]_[A-Za-z0-9]{20,}|\bBearer\s+[A-Za-z0-9._-]{20,}', serialized):
            raise ValueError('Credential-like content detected in selected trace; extraction not saved')
        if sha(source) != source_hash: raise ValueError('Trace source changed during extraction')
        selected = {'schema': 'clockbase-preserved-trace-v1', 'study_id': study, 'entry': normalized,
                    'serialization': {'format': 'JSON', 'tuple_paths': tuple_paths, 'truncated': False},
                    'source_sha256': source_hash}
        write_new(directory / 'selected-trace.json', selected)
        index = [{'message_index': i, 'characters': len(message),
                  'sha256_utf8': hashlib.sha256(message.encode()).hexdigest(),
                  'mentioned_studies': sorted(set(re.findall(r'\bGSE\d+\b', message)))}
                 for i, message in enumerate(messages)]
        write_new(directory / 'message-index.json', index)
        binding = verify_trace_request(messages, study)
        write_new(directory / 'trace-request-check.json', binding)
        report = {'status': 'PASS', 'scope': config['scope'], 'study_id': study, 'source_entry_count': len(traces),
                  'selected_entry_fields': sorted(entry), 'message_count': len(messages),
                  'mentioned_studies': sorted({sid for row in index for sid in row['mentioned_studies']}),
                  'selected_trace': str(directory / 'selected-trace.json'), 'message_index': str(directory / 'message-index.json'),
                  'selected_trace_sha256': sha(directory / 'selected-trace.json'),
                  'pickle_executable_objects_permitted': False, 'trace_code_executed': False,
                  'trace_request_status': binding['status'],
                  'clock_execution_performed': False, 'scientific_status': 'CANNOT VERIFY',
                  'secret_scan_scope': 'Known token prefixes and bearer credentials; not a guarantee that all sensitive content is recognized'}
        write_new(directory / 'inspection.json', report)
        return report
    return store.run('00_trace_extraction', extract)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--study', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-bytes', type=int, default=2 * 1024**3)
    args = parser.parse_args()
    try:
        print(json.dumps(extract_trace(args.source, args.study, args.output, args.max_bytes), indent=2))
    except Exception as error:
        raise SystemExit(f'Trace inspection stopped: {type(error).__name__}: {error}')
