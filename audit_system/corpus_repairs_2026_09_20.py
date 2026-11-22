"""Versioned corpus correction candidates and independently checked withholding.

This stage never replaces an original or assigns a scientific repair. Historical
sample lineage and orientation are not inferred from matching numbers or agent
agreement. Literal typed-field defects can withhold a record in a separate CSV;
missing/undefined values and scientific disagreements remain unresolved.
"""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys

from .legacy_analysis_review_2026_09_19 import SCIENTIFIC_CHECKS
from .record_corrections import CorrectionError, patch_csv
from .statistical_rerun_2026_09_20 import WELCH, STUDENT, PAIRED, OLS


SCHEMA = 'clockbase-corpus-repairs-v1'
STATUSES = {'PASS', 'FAIL', 'CANNOT VERIFY'}
PROBABILITIES = {'p_value', 'pvalue', 'fdr', 'adjusted_p_value', 'q_value', 'qvalue'}
MISSING = {'', 'na', 'nan', 'n/a', '<na>', 'none', 'null', 'missing', 'undefined', 'not available', 'cannot verify'}
# This is an explicit whole-record release policy after a verified material
# typed-field defect. Other cleared fields are withheld, not individually called
# erroneous. Rubric components and unrelated records remain unchanged.
RELEASE_FIELDS = PROBABILITIES | {
    'effect_size', 'mean_difference', 'standard_error', 'confidence_interval', 'ci_low', 'ci_high',
    'neg_log10_pvalue', 'neg_log10_fdr', 'is_significant', 'direction', 'trust_score', 'final_score', 'rank', 'ranking',
}


def _bytes(value):
    return (json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, indent=2) + '\n').encode()


def _hash(raw):
    return hashlib.sha256(raw).hexdigest()


def _sha(path):
    try:
        return _hash(Path(path).read_bytes())
    except OSError as error:
        raise ValueError(f'Cannot read repair evidence {path}: {error}') from error


def _write(path, value):
    path = Path(path)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as handle:
            handle.write(value if isinstance(value, bytes) else _bytes(value))
    except OSError as error:
        raise ValueError(f'Cannot preserve corpus repair output {path}: {error}') from error


def _json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot parse corpus repair evidence {path}: {error}') from error


def _parse(raw):
    """Independent reader with literal record spans; ragged records stay visible."""
    text = raw.decode('utf-8')
    lines = text.splitlines(keepends=True)
    reader = csv.reader(io.StringIO(text, newline=''), strict=True)
    header = next(reader, [])
    previous = reader.line_num
    records, spans = [], []
    for row in reader:
        records.append(row)
        spans.append(''.join(lines[previous:reader.line_num]).encode())
        previous = reader.line_num
    # Re-read only the header to retain its physical newline/quoting bytes.
    header_reader = csv.reader(io.StringIO(text, newline=''), strict=True)
    next(header_reader, None)
    header_span = ''.join(lines[:header_reader.line_num]).encode()
    return header, records, spans, header_span


def _contract(artifact, header):
    declared = artifact.get('numeric_fields', [])
    supported = PROBABILITIES | {'effect_size', 'mean_difference', 'standard_error', 'ci_low', 'ci_high'}
    if not isinstance(declared, list) or any(not isinstance(field, str) or field not in header or field not in supported for field in declared):
        raise ValueError('Explicit numeric-field declaration is malformed or names an absent column')
    return {field: 'probability' if field in PROBABILITIES else 'finite_number'
            for field in (set(header) & PROBABILITIES) | set(declared)}


def _detect(fields, contract):
    checks = []
    for name, kind in sorted(contract.items()):
        value = fields[name]
        if value.strip().lower() in MISSING:
            status, reason = 'CANNOT VERIFY', 'Missing or explicitly undefined value; no material defect is inferred'
        else:
            try:
                if kind == 'probability':
                    number = Decimal(value)
                    valid = number.is_finite() and 0 <= number <= 1
                else:
                    number = float(value)
                    valid = math.isfinite(number)
            except (TypeError, ValueError, OverflowError, InvalidOperation):
                valid = False
            status = 'PASS' if valid else 'FAIL'
            reason = 'Literal value satisfies its declared numeric domain' if valid else 'Literal value contradicts its declared numeric domain'
        checks.append({'field': name, 'contract': kind, 'status': status, 'original_value': value, 'reason': reason})
    return checks


def _independent_literal_defect(value, kind):
    """Decimal readback; separate from the producer's binary-float detection."""
    if not isinstance(value, str) or value.strip().lower() in MISSING:
        return False
    try:
        number = Decimal(value)
    except InvalidOperation:
        return True
    if not number.is_finite():
        return True
    if kind == 'probability':
        return number < 0 or number > 1
    if kind == 'finite_number':
        return abs(number) > Decimal.from_float(sys.float_info.max)
    raise ValueError('Unknown independent numeric contract')


def _numerical_binding(receipt, path, calculations):
    """Require the on-disk R receipt, its checksum, request and every output."""
    if receipt is None or path is None:
        return {'status': 'CANNOT VERIFY', 'reason': 'File-bound independent numerical receipt was not supplied'}, {}
    try:
        path = Path(path).resolve()
        stored = _json(path)
        if stored != receipt or _sha(path) != path.with_name('receipt.sha256').read_text().strip():
            raise ValueError('Independent numerical receipt differs from its checksum or supplied record')
        if stored.get('schema') != 'clockbase-independent-statistics-v1':
            raise ValueError('Unsupported independent numerical receipt schema')
        for name, expected in stored['outputs_sha256'].items():
            target = (path.parent / name).resolve()
            if not target.is_relative_to(path.parent) or _sha(target) != expected:
                raise ValueError('Independent numerical output changed')
        request_path = path.parent / 'request.json'
        request = _json(request_path)
        if _sha(request_path) != stored['request_sha256'] or request['calculations'] != calculations:
            raise ValueError('Independent numerical request does not bind these calculations')
        prediction = request['predictions']
        if (not prediction.get('sha256') or _sha(prediction['path']) != prediction['sha256']
                or stored['predictions_sha256_before'] != prediction['sha256']
                or stored['predictions_sha256_after'] != prediction['sha256']
                or stored['source_hash_status'] != 'PASS'):
            raise ValueError('Independent numerical source hash does not bind the preserved calculation')
        results = stored['results']
        if len(results) != len(calculations) or {r['input_index'] for r in results} != set(range(len(calculations))):
            raise ValueError('Independent numerical results omit or duplicate input indices')
        by_index = {r['input_index']: r for r in results}
        for index, calculation in enumerate(calculations):
            if any(by_index[index].get(key) != calculation.get(key) for key in ('comparison_id', 'clock')):
                raise ValueError('Independent numerical result identity differs from its request')
        return {'status': 'PASS', 'receipt_sha256': _sha(path), 'request_sha256': stored['request_sha256'],
                'predictions_sha256': prediction['sha256'], 'scope': 'Receipt and source binding; individual numerical verdicts remain separate'}, by_index
    except (ValueError, KeyError, TypeError, OSError) as error:
        return {'status': 'FAIL', 'reason': str(error)}, {}


def _candidate(row, decision, design, calculations, numeric, numeric_binding):
    proposal = {'row_id': row['row_id'], 'status': 'CANNOT VERIFY', 'calculation': None,
                'historical_sample_lineage': 'CANNOT VERIFY', 'original_orientation': 'CANNOT VERIFY',
                'clock_applicability': 'CANNOT VERIFY', 'scientific_status': 'CANNOT VERIFY',
                'accepted_scientific_repair': False, 'original_replaced': False}
    refs = decision.get('source_comparison_ids')
    if not isinstance(refs, list) or len(refs) != 1:
        return proposal | {'reason': 'Exactly one explicit source comparison binding is required'}
    comparisons = [c for c in design.get('comparisons', []) if c.get('id') == refs[0]]
    fields = row.get('fields')
    if len(comparisons) != 1 or not fields or not isinstance(fields.get('clock'), str):
        return proposal | {'reason': 'Original clock or source comparison is absent or ambiguous'}
    binding = row.get('identity_binding', {})
    bound_key = binding.get('bound_study_key')
    if (binding.get('identity_status') != 'PASS' or not isinstance(bound_key, str)
            or not re.fullmatch(r'(human_rna|mouse_rna|human_dnam):GSE[0-9]+', bound_key)
            or bound_key.split(':')[1] != design.get('study_id')
            or binding.get('candidate_study_keys') != [bound_key] or binding.get('identity_issues')):
        return proposal | {'reason': 'Original row lacks an unambiguous verified source-study identity binding'}
    if 'gse_id' in fields and fields['gse_id'] != design.get('study_id'):
        return proposal | {'reason': 'Literal source-study identity contradicts the verified row binding'}
    proposal['identity_binding'] = binding
    comparison = comparisons[0]
    hits = [(i, c) for i, c in enumerate(calculations) if c.get('clock') == fields['clock'] and c.get('comparison_id') == refs[0]]
    if len(hits) != 1:
        return proposal | {'reason': 'Exact source comparison and literal clock do not identify one calculation'}
    index, calculation = hits[0]
    proposal.update(source_comparison_id=refs[0], clock=fields['clock'], calculation_index=index)
    if design.get('frozen_before_saved_outcomes') is not True or comparison.get('execution_design_status') != 'PASS':
        return proposal | {'reason': 'Source design was not frozen and executable before saved outcomes'}
    expected_design = _hash(json.dumps(design, sort_keys=True, allow_nan=False).encode())
    if calculation.get('source_design_sha256') != expected_design or any(
            calculation.get(key) != comparison.get(key) for key in ('case_ids', 'control_ids')):
        return proposal | {'reason': 'Calculation design hash or exact source groups differs from the frozen design'}
    result = numeric.get(index, {})
    if (calculation.get('status') != 'PASS' or not isinstance(calculation.get('calculation'), dict)
            or numeric_binding['status'] != 'PASS' or result.get('numerical_status') != 'PASS'
            or result.get('status') != 'PASS' or any(check.get('status') != 'PASS' for check in result.get('checks', []))):
        return proposal | {'reason': 'Source calculation lacks a file-bound independent R PASS'}
    if not result.get('checks'):
        return proposal | {'reason': 'Independent numerical checks are absent'}
    normalized = re.sub(r'[^a-z0-9]+', ' ', comparison.get('method', '').lower()).strip()
    expected_mode = ('descriptive' if comparison.get('mode') == 'descriptive' else
        'welch' if normalized in WELCH else 'student' if normalized in STUDENT else
        'paired' if normalized in PAIRED else
        ('ols_hc3' if normalized == 'ols hc3 adjusted contrast' else 'ols') if normalized in OLS else None)
    if expected_mode is None or result.get('method') != expected_mode:
        return proposal | {'reason': 'Independent method does not match the frozen source method'}
    return proposal | {'status': 'PASS', 'reason': 'Conditional calculation verified under the explicit proposed source comparison; historical repair is not accepted',
        'scope': 'Source-comparison candidate only; matching arithmetic does not certify historical membership, orientation or aging interpretation',
        'case_ids': calculation['case_ids'], 'control_ids': calculation['control_ids'],
        'calculation': calculation['calculation'], 'independent_receipt_sha256': numeric_binding['receipt_sha256'],
        'predictions_sha256': numeric_binding['predictions_sha256']}


def verify_corpus_repairs(output_dir):
    """Read outputs independently: exact field changes, unchanged records, sources.

    This function does not call patch_csv or the defect-detection function. It
    checks release mechanics and literal domains; it never adjudicates science.
    """
    output = Path(output_dir).resolve()
    manifest = _json(output / 'withholding-manifest.json')
    checks = []
    def check(passed, name, **detail):
        checks.append({'name': name, 'status': 'PASS' if passed else 'FAIL', **detail})
    for file in manifest['files']:
        source = Path(file['source_path'])
        original = output / file['preserved_original']
        check(_sha(source) == file['before_sha256'], 'source hash remains unchanged', artifact_id=file['artifact_id'])
        raw = original.read_bytes()
        check(_hash(raw) == file['before_sha256'], 'preserved original bytes match', artifact_id=file['artifact_id'])
        if not file.get('withheld_copy'):
            continue
        changed = (output / file['withheld_copy']).read_bytes()
        check(_hash(changed) == file['after_sha256'], 'versioned copy hash matches', artifact_id=file['artifact_id'])
        header, before, old_spans, old_header = _parse(raw)
        after_header, after, new_spans, new_header = _parse(changed)
        check(header == after_header and old_header == new_header and len(before) == len(after), 'literal header and row count preserved')
        operations = {op['record_number']: op for op in file['operations']}
        for number, (old, new) in enumerate(zip(before, after), 2):
            op = operations.get(number)
            if op is None:
                check(old == new and old_spans[number-2] == new_spans[number-2], 'unmodified record bytes preserved', record_number=number)
                continue
            check(dict(zip(header, old)) == op['before'], 'original target identity matches', record_number=number)
            expected = [op['changes'].get(field, value) for field, value in zip(header, old)]
            check(new == expected and all(value == '' for value in op['changes'].values()), 'only declared withholding fields changed', record_number=number)
            defects = op['literal_defects']
            check(bool(defects) and all(_independent_literal_defect(op['before'][d['field']], d['contract']) for d in defects),
                  'material typed-field defect independently confirmed', record_number=number)
    return {'schema': 'clockbase-corpus-withholding-readback-v1',
            'status': 'PASS' if checks and all(c['status'] == 'PASS' for c in checks) else 'CANNOT VERIFY' if not checks else 'FAIL',
            'checks': checks, 'check_counts': dict(Counter(c['status'] for c in checks)),
            'scientific_status': 'CANNOT VERIFY', 'accepted_scientific_repairs': 0,
            'scope': 'Independent source/copy/field readback and literal-domain withholding; no scientific replacement acceptance'}


def build_repair_candidates(original_rows, frozen_design, calculations, numerical_receipt,
                            scoped_checks, row_decisions, output_dir, *, artifact_sources,
                            source_root, evidence_index=None, numerical_receipt_path=None):
    """Preserve originals and build versioned candidates/withholding copies.

    Artifact sources are the corpus artifact list (artifact_id, relative path,
    sha256). Optional artifact ``numeric_fields`` explicitly establishes scalar
    numeric column semantics such as effect_size; no effect contract is guessed.
    The complete 16-check records are retained for every supplied source row.
    """
    root, output = Path(source_root).resolve(), Path(output_dir).resolve()
    artifacts = {a['artifact_id']: a for a in artifact_sources}
    if len(artifacts) != len(artifact_sources):
        raise ValueError('Duplicate repair source artifact identities')
    rows = {r['row_id']: r for r in original_rows}
    decisions = {r['row_id']: r for r in row_decisions}
    if len(rows) != len(original_rows) or len(decisions) != len(row_decisions) or set(rows) != set(decisions):
        raise ValueError('Every original row requires exactly one preserved decision')
    coverage = []
    for scoped in scoped_checks:
        if set(scoped.get('checks', {})) != set(SCIENTIFIC_CHECKS):
            raise ValueError('Every review scope requires all 16 scientific checks')
        if any(c.get('status') not in STATUSES for c in scoped['checks'].values()):
            raise ValueError('Unknown scoped scientific status')
        coverage.extend(scoped['row_ids'])
    if Counter(coverage) != Counter(rows.keys()):
        raise ValueError('Scoped checks omit, duplicate or invent original rows')
    if any(d.get('status') not in STATUSES or d.get('action') not in {'retain', 'withhold', 'propose_rerun', 'unresolved'} for d in decisions.values()):
        raise ValueError('Unknown original row decision or action')
    source_data = {}
    for artifact_id, artifact in artifacts.items():
        relative = Path(artifact['path'])
        path = (root / relative).resolve()
        if relative.is_absolute() or not path.is_relative_to(root) or path.is_relative_to(output):
            raise ValueError('Artifact source must be inside its source root and outside the repair output')
        if not path.is_file():
            raise FileNotFoundError('Cannot proceed - required data file not found: ' + str(path))
        raw = path.read_bytes()
        if _hash(raw) != artifact['sha256'] or raw.startswith(b'version https://git-lfs.github.com/spec/v1'):
            raise ValueError('Original source changed or lacks its LFS payload')
        try:
            header, records, _, _ = _parse(raw)
            error = None
        except (UnicodeError, csv.Error) as caught:
            header, records, error = [], [], str(caught)
        supplied = [r for r in original_rows if r['artifact_id'] == artifact_id]
        if error is None and (len(supplied) != len(records) or {r['logical_record'] for r in supplied} != set(range(2, len(records)+2))):
            raise ValueError('Original row inventory does not cover every source record')
        for row in supplied:
            if row['artifact_sha256'] != artifact['sha256']:
                raise ValueError('Original row source hash mismatch')
            if error is None:
                values = records[row['logical_record']-2]
                fields = dict(zip(header, values)) if len(header) == len(set(header)) == len(values) else None
                if row['header'] != header or row['values'] != values or row['fields'] != fields:
                    raise ValueError('Original row content differs from its immutable CSV source')
        source_data[artifact_id] = (path, raw, header, records, error)
    if set(r['artifact_id'] for r in original_rows) - set(artifacts):
        raise ValueError('Original row references an absent source artifact')
    numerical_binding, independent_rows = _numerical_binding(numerical_receipt, numerical_receipt_path, calculations)
    inputs = {'schema': SCHEMA, 'original_rows': original_rows, 'frozen_source_design': frozen_design,
        'calculations': calculations, 'numerical_receipt': numerical_receipt, 'scoped_checks': scoped_checks,
        'row_decisions': row_decisions, 'artifact_sources': artifact_sources, 'source_root': str(root),
        'source_evidence_index': evidence_index, 'numerical_binding': numerical_binding,
        'implementation_sha256': _sha(__file__),
        'dependency_sha256': {name: _sha(Path(__file__).with_name(name)) for name in
            ('record_corrections.py', 'statistical_rerun_2026_09_20.py', 'legacy_analysis_review_2026_09_19.py')}}
    input_hash = _hash(_bytes(inputs))
    output.mkdir(parents=True, exist_ok=True)
    receipt_path = output / 'receipt.json'
    if receipt_path.is_file():
        prior = _json(receipt_path)
        if _sha(receipt_path) != (output / 'receipt.sha256').read_text().strip() or prior['inputs_sha256'] != input_hash:
            raise ValueError('Repair receipt or inputs changed; use a new output directory')
        for name, expected in prior['outputs_sha256'].items():
            target = (output / name).resolve()
            if not target.is_relative_to(output) or _sha(target) != expected:
                raise ValueError('Repair output changed on resume')
        return prior
    if any(output.iterdir()):
        raise FileExistsError('Incomplete corpus repair output is preserved; use a new directory')
    _write(output / 'inputs.json', inputs)
    manifest = {'schema': SCHEMA, 'files': [], 'policy': 'Withhold record release fields after independently checkable literal typed-field defects; no numeric replacement'}
    row_outputs, candidates = [], []
    for index, (artifact_id, artifact) in enumerate(artifacts.items()):
        path, raw, header, records, parse_error = source_data[artifact_id]
        name = f'{index:04d}-{artifact["sha256"][:12]}.csv'
        preserved = 'originals/' + name
        _write(output / preserved, raw)
        file = {'artifact_id': artifact_id, 'source_path': str(path), 'before_sha256': artifact['sha256'],
                'preserved_original': preserved, 'row_count': len(records), 'operations': [],
                'source_parse_status': 'CANNOT VERIFY' if parse_error else 'PASS'}
        contract = _contract(artifact, header) if not parse_error else {}
        file['numeric_field_contracts'] = contract
        parseable_for_patch = not parse_error and bool(header) and len(header) == len(set(header)) and all(len(r) == len(header) for r in records)
        for row in [r for r in original_rows if r['artifact_id'] == artifact_id]:
            decision = decisions[row['row_id']]
            literal = _detect(row['fields'], contract) if row['fields'] is not None and not parse_error else []
            defects = [c for c in literal if c['status'] == 'FAIL']
            candidate = _candidate(row, decision, frozen_design, calculations, independent_rows, numerical_binding)
            candidates.append(candidate)
            item = {'row_id': row['row_id'], 'artifact_id': artifact_id, 'logical_record': row['logical_record'],
                'original_artifact_sha256': artifact['sha256'], 'model_decision': decision,
                'literal_checks': literal, 'candidate_status': candidate['status'], 'scientific_status': 'CANNOT VERIFY',
                'effective_action': 'CANDIDATE_ONLY' if candidate['status'] == 'PASS' else 'UNRESOLVED',
                'original_overwritten': False, 'accepted_scientific_repair': False}
            if defects and parseable_for_patch:
                changes = {field: '' for field in header if field in RELEASE_FIELDS and row['fields'][field] != ''}
                if changes:
                    file['operations'].append({'record_number': row['logical_record'], 'before': row['fields'],
                        'changes': changes, 'row_id': row['row_id'], 'literal_defects': defects})
                    item.update(effective_action='WITHHOLD_IN_VERSIONED_COPY', withheld_fields=list(changes),
                        withholding_scope='Whole-record release withheld after literal typed-field defect; other cleared fields are not individually declared erroneous')
            elif decision['action'] == 'withhold':
                item.update(effective_action='WITHHOLDING_PROPOSED',
                    unresolved_reason='Source/agent judgment alone does not satisfy this literal correction gate')
            if not parseable_for_patch:
                item['withholding_limitation'] = 'Malformed or ambiguous source CSV is preserved; no automatic file rewrite'
            row_outputs.append(item)
        if file['operations']:
            try:
                changed = patch_csv(raw, file)
            except CorrectionError as error:
                raise ValueError('Versioned withholding failed its existing correction guard: ' + str(error)) from error
            file['withheld_copy'] = 'withheld/' + name
            file['after_sha256'] = _hash(changed)
            _write(output / file['withheld_copy'], changed)
        manifest['files'].append(file)
    _write(output / 'withholding-manifest.json', manifest)
    readback = verify_corpus_repairs(output)
    _write(output / 'independent-readback.json', readback)
    for item in row_outputs:
        if item['effective_action'] == 'WITHHOLD_IN_VERSIONED_COPY':
            item['withholding_verification'] = readback['status']
            if readback['status'] != 'PASS':
                item['effective_action'] = 'WITHHOLDING_UNVERIFIED'
    _write(output / 'row-decisions.json', row_outputs)
    _write(output / 'correction-candidates.json', candidates)
    receipt = {'schema': SCHEMA, 'recorded_at_utc': datetime.now(timezone.utc).isoformat(),
        'status': readback['status'], 'inputs_sha256': input_hash, 'original_rows': len(original_rows),
        'original_files': len(artifacts), 'originals_preserved': readback['status'] == 'PASS',
        'row_action_counts': dict(Counter(r['effective_action'] for r in row_outputs)),
        'attempted_withholding_rows': sum(len(f['operations']) for f in manifest['files']),
        'attempted_withholding_files': sum(bool(f['operations']) for f in manifest['files']),
        'withheld_rows': sum(len(f['operations']) for f in manifest['files']) if readback['status'] == 'PASS' else 0,
        'withheld_files': sum(bool(f['operations']) for f in manifest['files']) if readback['status'] == 'PASS' else 0,
        'conditional_numeric_candidates': sum(c['status'] == 'PASS' for c in candidates),
        'unresolved_candidates': sum(c['status'] != 'PASS' for c in candidates),
        'numerical_binding': numerical_binding, 'independent_readback_counts': readback['check_counts'],
        'scientific_status': 'CANNOT VERIFY', 'accepted_scientific_repairs': 0, 'numeric_replacements': 0,
        'new_rank_assignments': 0, 'new_adjusted_p_values': 0, 'clock_runs': 0, 'model_calls': 0, 'api_cost_usd': 0,
        'outputs_sha256': {str(p.relative_to(output)): _sha(p) for p in sorted(output.rglob('*')) if p.is_file()},
        'limitations': ['Original sample lineage, comparison orientation and clock applicability remain unresolved.',
            'All 16 scoped model checks are preserved, but agreement is not independent scientific truth.',
            'Only declared literal numeric defects trigger automatic withholding; scientific corrections require a separate source-supported acceptance gate.']}
    _write(receipt_path, receipt)
    _write(output / 'receipt.sha256', (_sha(receipt_path) + '\n').encode())
    return receipt
