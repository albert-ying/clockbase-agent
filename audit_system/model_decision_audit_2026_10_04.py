"""Re-derive saved model decisions from raw ledger responses without sending requests.

  python -m audit_system.model_decision_audit_2026_10_04 --verifier DIR --out DIR
"""
import argparse
import csv
import hashlib
import json
import time
from collections import defaultdict
from pathlib import Path

from audit_system import design_audit_2026_09_30 as design_audit
from audit_system import design_reconstruction_2026_09_28 as reconstruction
from audit_system import individuals_audit_2026_10_10 as individuals
from audit_system import intervention_classification_2026_09_29 as classification
from audit_system import label_adjudication_2026_09_27 as adjudication
from audit_system.adjudication_files import decision_files
from audit_system.paths import GEO, VERIFIER_ROOT

SCHEMA = 'clockbase-model-decision-audit-v1'
KEY = ('modality', 'gse_id', 'comparison_name', 'condition_name')
FIELDS = {
    'adjudication': ('decision', 'judgments', 'issues_agreed'),
    'reconstruction': ('status', 'model', 'case_gsms', 'control_gsms', 'rules', 'sizes'),
    'classification': ('status', 'answers', 'contrast_type', 'intervention_class', 'manipulation',
                       'canonical_target', 'relevance_score', 'relevance_score_abs_diff', 'relevance_codes_agreed'),
    'design_audit': ('status', 'reasons', 'orientation', 'specimen', 'tumor', 'answers'),
    'individuals_audit': ('status', 'findings', 'individuals', 'findings_by_reader'),
    'rows': ('status', 'source', 'kind', 'case_gsms', 'control_gsms', 'n_case', 'n_control', 'design_rules'),
}
PARSERS = {'adjudication': adjudication.parse_judgment, 'reconstruction': reconstruction.parse_rules,
           'classification': classification.parse, 'design_audit': design_audit.parse,
           'individuals_audit': individuals.parse}
NON_MODEL = {'METADATA_UNAVAILABLE', 'GEO_INDIVIDUAL_FIELD', 'NOT_READ_ARMS_ABOVE_LIMIT'}


def read_records(path):
    """Read JSONL with a source line number in any error."""
    try:
        with Path(path).open() as handle:
            for number, line in enumerate(handle, 1):
                if line.strip():
                    try:
                        record = json.loads(line)
                        if not isinstance(record, dict):
                            raise ValueError('record must be an object')
                    except (ValueError, TypeError) as error:
                        raise ValueError(f'{path}:{number}: {error}') from error
                    yield number, record
    except OSError as error:
        raise OSError(f'Cannot read {path}: {error}') from error


def source(path):
    path = Path(path)
    digest = hashlib.sha256()
    try:
        with path.open('rb') as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b''):
                digest.update(block)
        stat = path.stat()
    except OSError as error:
        raise OSError(f'Cannot record source {path}: {error}') from error
    return {'path': str(path.resolve()), 'sha256': digest.hexdigest(), 'size': stat.st_size,
            'mtime_ns': stat.st_mtime_ns}


def json_value(value):
    """Use the same JSON representation as saved records (including tuple sizes and integer ids)."""
    return json.loads(json.dumps(value))


def differences(stored, derived, kind):
    derived = json_value(derived)
    return {field: {'stored': stored.get(field), 'derived': derived.get(field)}
            for field in FIELDS[kind] if stored.get(field) != derived.get(field)}


def example(report, name, context, detail):
    examples = report.setdefault(name + '_examples', [])
    if len(examples) < 20:
        examples.append({**context, **detail})


class Responses:
    def __init__(self, path, parser):
        self.records, self.parsed, self.parser = defaultdict(list), {}, parser
        if path.exists():
            for number, record in read_records(path):
                self.records[(record.get('request_sha256'), record.get('model_requested'))].append((number, record))

    def parse(self, number, record):
        if number not in self.parsed:
            text = record.get('response_text')
            try:
                parsed = self.parser(text) if isinstance(text, str) and text.strip() else None
                error = None if parsed else 'Parser returned no answer'
            except Exception as failure:
                parsed, error = None, f'{type(failure).__name__}: {failure}'
            self.parsed[number] = parsed, error
        return self.parsed[number]

    def receipt(self, receipt, report, context, saved=None):
        records = self.records.get((receipt.get('request_sha256'), receipt.get('model')), [])
        with_response = [(n, r) for n, r in records if isinstance(r.get('response_text'), str)]
        if not with_response:
            report['missing_receipts'] += 1
            example(report, 'missing_receipts', context, {'receipt': receipt,
                    'reason': 'No matching model/request with a stored response', 'ledger_entries': len(records)})
            return None, False, {}
        report['receipts_found'] += 1
        matching = [(n, r) for n, r in with_response if r.get('status') == receipt.get('status')]
        candidates = matching or with_response
        variants = {r['response_text'] for _, r in candidates}
        if len(variants) > 1:
            report['ambiguous_receipts'] += 1
            example(report, 'ambiguous_receipts', context, {'receipt': receipt, 'ledger_lines': [n for n, _ in candidates]})
        number, record = candidates[-1]
        # Receipts have no ledger line or time. A saved judgment can identify one of several raw responses;
        # only freshly parsed raw text is then used to derive the decision.
        if saved is not None and len(variants) > 1:
            for n, r in candidates:
                parsed, _ = self.parse(n, r)
                if parsed == saved:
                    number, record = n, r
                    break
        parsed, error = self.parse(number, record)
        delta = {}
        if not matching:
            delta['receipt_status'] = {'stored': receipt.get('status'), 'derived': record.get('status')}
        parsed_status = 'OK' if parsed else 'UNPARSEABLE'
        if record.get('status') != parsed_status:
            delta['response_status'] = {'stored': record.get('status'), 'derived': parsed_status}
        if error:
            report['unparseable_responses'] += 1
            example(report, 'unparseable_responses', context, {'receipt': receipt, 'ledger_line': number, 'error': error})
        return parsed, True, delta


def comparison_key(record):
    return tuple(record.get(k) for k in KEY)


def comparison_id(record, kind, index, candidates):
    if kind == 'design_audit':
        return candidates.get(comparison_key(record))
    if kind == 'classification':
        ids = {a.get('id') for a in record.get('answers', []) if isinstance(a, dict) and 'id' in a}
        return int(next(iter(ids))) if len(ids) == 1 else candidates.get(comparison_key(record))
    return index


def classification_candidates(folder, report):
    """Recover missing ids only from an atlas matching the run's recorded input hash."""
    summary = folder / 'summary.json'
    if not summary.exists():
        return {}
    try:
        saved = json.loads(summary.read_text())
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot read classification input provenance {summary}: {error}') from error
    digest = saved.get('atlas_sha256')
    if not digest:
        return {}
    root = folder.parent
    paths = sorted(set(root.glob('atlas_v3_combined.csv')) | set(root.glob('*/atlas_v3_combined.csv'))
                   | set(root.glob('final-*/atlas-v3/atlas_v3_combined.csv')))
    for path in paths:
        provenance = source(path)
        if provenance['sha256'] != digest:
            continue
        report['sources'].extend([source(summary), provenance])
        candidates = {}
        with path.open(newline='') as handle:
            # run() enumerates usable rows after dropping duplicate comparison keys, in input order.
            for record in csv.DictReader(handle):
                if record.get('v3_final_class') in ('VERIFIED', 'PROVISIONAL'):
                    key = comparison_key(record)
                    if key not in candidates:
                        candidates[key] = len(candidates)
        return candidates
    return {}


def rederive(record, kind, responses, report, context, cid=None, meta=None):
    receipts = record.get('receipts') or []
    if not receipts:
        if record.get('decision', record.get('status')) in NON_MODEL:
            report['non_model_records'] += 1
            return None, {}, 'non_model'
        report['missing_receipts'] += 1
        example(report, 'missing_receipts', context, {'reason': 'Decision has no receipts'})
        return None, {}, 'missing'
    answers, judgments, delta, complete = {}, [], {}, True
    saved = iter(record.get('judgments') or [])
    for i, receipt in enumerate(receipts):
        hint = next(saved, None) if kind == 'adjudication' and receipt.get('status') == 'OK' else None
        parsed, found, mismatch = responses.receipt(receipt, report, context, hint)
        complete &= found
        delta.update({f'receipt_{i}_{k}': v for k, v in mismatch.items()})
        model = receipt.get('model')
        if kind == 'adjudication':
            if parsed:
                judgments.append(parsed)
        elif kind == 'reconstruction':
            if receipt.get('batch') == 'whole':
                answers[model] = reconstruction.normalize_ids(parsed) if parsed else None
            elif model not in answers or answers[model] is not None:
                answers.setdefault(model, {})
                if parsed:
                    answers[model].update(reconstruction.normalize_ids(parsed))
                else:
                    answers[model] = None
        else:
            answers[model] = parsed
    if not complete:
        return None, delta, 'missing'
    if len({r.get('model') for r in receipts}) != 2:
        report['unresolved_records'] += 1
        example(report, 'unresolved_records', context, {'reason': 'Expected two independent model receipts'})
        return None, delta, 'unresolved'
    if kind == 'adjudication':
        derived = {'decision': adjudication.decide(judgments) if len(judgments) == 2 else 'ADJUDICATION_INCOMPLETE',
                   'judgments': judgments,
                   'issues_agreed': sorted(set.intersection(*[set(j.get('issues') or []) for j in judgments]))
                   if len(judgments) == 2 else []}
    elif kind == 'individuals_audit':  # one request per comparison and reader; the record carries its arms
        derived = individuals.reconcile(record['case_gsms'], record['control_gsms'], [answers.get(r.get('model')) for r in receipts])
    elif not all(answers.values()):
        derived = {'status': 'MODEL_CALL_FAILED'}
    elif cid is None:
        report['unresolved_records'] += 1
        example(report, 'unresolved_records', context, {'reason': 'Comparison id unavailable; no positional guess used'})
        return None, delta, 'unresolved'
    elif kind == 'reconstruction':
        derived = reconstruction.reconcile(meta, answers, cid)
    elif not all(cid in a for a in answers.values()):
        derived = {'status': 'MODEL_ANSWER_MISSING' if kind == 'design_audit' else 'MODEL_ANSWER_TRUNCATED'}
    else:
        reconcile = design_audit.reconcile if kind == 'design_audit' else classification.reconcile
        derived = reconcile([a[cid] for a in answers.values()])
    delta.update(differences(record, derived, kind))
    return derived, delta, 'differing' if delta else 'identical'


def audit_run(folder, kind):
    report = {'run': folder.name, 'kind': kind, 'records': 0, 'receipts_found': 0, 'rederived_identically': 0,
              'differing': 0, 'missing_receipts': 0, 'records_missing_receipts': 0, 'unparseable_responses': 0,
              'ambiguous_receipts': 0, 'unresolved_records': 0, 'non_model_records': 0, 'sources': [], 'files': {}}
    files = decision_files(folder) if kind == 'adjudication' else [folder / 'comparisons.jsonl']
    files = [f for f in files if f.exists()]
    ledger = folder / 'ledger.jsonl'
    if ledger.exists():
        report['sources'].append(source(ledger))
    responses = Responses(ledger, PARSERS[kind])
    candidates = {}
    if kind == 'design_audit' and (folder / 'candidates.csv').exists():
        path = folder / 'candidates.csv'
        report['sources'].append(source(path))
        with path.open(newline='') as handle:
            # run() saves the same sorted candidate order that build_items() enumerates.
            candidates = {comparison_key(r): i for i, r in enumerate(csv.DictReader(handle))}
    elif kind == 'classification':
        candidates = classification_candidates(folder, report)
    indices, parents, meta_gse, meta = defaultdict(int), {}, None, None

    def record_result(path, number, record, outcome, delta):
        report['records'] += 1
        stats = report['files'].setdefault(path.name, {'records': 0, 'identical': 0, 'differing': 0, 'missing': 0,
                                                     'unresolved': 0, 'non_model': 0})
        stats['records'] += 1
        stats[outcome] += 1
        if outcome in ('identical', 'differing', 'missing'):
            field = {'identical': 'rederived_identically', 'differing': 'differing',
                     'missing': 'records_missing_receipts'}[outcome]
            report[field] += 1
        if delta:
            example(report, 'differing', {'file': path.name, 'line': number,
                    **{k: record[k] for k in ('row_ordinal', 'gse_id', 'comparison_name') if k in record}}, {'fields': delta})

    for path in files:
        report['sources'].append(source(path))
        for number, record in read_records(path):
            context = {'file': path.name, 'line': number, **{k: record[k] for k in ('gse_id', 'row_ordinal') if k in record}}
            cid = comparison_id(record, kind, indices[record.get('gse_id')], candidates)
            if kind == 'reconstruction':
                gse = record['gse_id']
                indices[gse] += 1
                if gse != meta_gse:
                    modality = folder.name.split('-')[1]
                    for sub in ('processed_metadata_preharmonize', 'processed_metadata'):
                        metadata_path = GEO / reconstruction.MODALITIES[modality]['folder'] / sub / f'{gse}.csv'
                        if metadata_path.exists():
                            report['sources'].append(source(metadata_path))
                            break
                    meta, meta_gse = reconstruction.load_meta(modality, gse), gse
                    if meta is None and record.get('status') != 'METADATA_UNAVAILABLE':
                        raise FileNotFoundError(f'Cannot proceed - required data file not found: {metadata_path}')
            try:
                derived, delta, outcome = rederive(record, kind, responses, report, context, cid, meta)
            except (ValueError, TypeError, KeyError, AttributeError, IndexError) as error:
                derived, delta, outcome = None, {}, 'unresolved'
                report['unresolved_records'] += 1
                example(report, 'unresolved_records', context, {'reason': f'{type(error).__name__}: {error}'})
            record_result(path, number, record, outcome, delta)
            if kind == 'reconstruction':
                for ordinal in record.get('row_ordinals', []):
                    parents[ordinal] = (record, derived, outcome)
    rows = folder / 'rows.jsonl'
    if kind == 'reconstruction' and rows.exists():
        report['sources'].append(source(rows))
        for number, record in read_records(rows):
            parent, derived, outcome = parents.get(record.get('row_ordinal'), (None, None, 'missing'))
            delta = {}
            if parent is None:
                report['missing_receipts'] += 1
                example(report, 'missing_receipts', {'file': rows.name, 'line': number},
                        {'reason': 'No comparison record supplies this exported row\'s receipts'})
            elif derived is not None:
                if derived['status'] == 'AGREED':
                    expected = {'status': 'UNIQUE', 'source': 'design_reconstruction', 'kind': 'two_sample',
                                'case_gsms': derived['case_gsms'], 'control_gsms': derived['control_gsms'],
                                'n_case': len(derived['case_gsms']), 'n_control': len(derived['control_gsms']),
                                'design_rules': derived['rules']}
                    delta = differences(record, expected, 'rows')
                else:
                    delta = {'status': {'stored': record.get('status'), 'derived': 'ROW_NOT_EXPORTABLE_' + derived['status']}}
                outcome = 'differing' if delta else 'identical'
            elif outcome in ('unresolved', 'non_model'):
                report['unresolved_records' if outcome == 'unresolved' else 'non_model_records'] += 1
            record_result(rows, number, record, outcome, delta)
    return report


def run(verifier, out):
    verifier, out = Path(verifier), Path(out)
    if out.exists():
        raise FileExistsError(f'Output directory must not exist: {out}')
    if not verifier.is_dir():
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {verifier}')
    # The audited runs are read-only: the output may sit in a final run directory, never inside a decision run.
    audited = ('adjudication-', 'design-', 'classification-', 'individuals-', 'membership-', 'replay-', 'link-')
    resolved = out.resolve()
    if resolved.is_relative_to(verifier.resolve()):
        relative = resolved.relative_to(verifier.resolve())
        if not relative.parts or relative.parts[0].startswith(audited):
            raise ValueError('Output must not be inside an audited run directory')
    reports = []
    for folder in sorted(verifier.iterdir()):
        if not folder.is_dir():
            continue
        kind = ('adjudication' if folder.name.startswith('adjudication-') else
                'design_audit' if folder.name.startswith('design-audit-') else
                'reconstruction' if folder.name.startswith(('design-dnam-', 'design-hs_rna-', 'design-mm_rna-')) else
                'classification' if folder.name.startswith('classification-') else
                'individuals_audit' if folder.name.startswith('individuals-audit-') else None)
        if kind and (list(folder.glob('decisions-*.jsonl')) if kind == 'adjudication' else (folder / 'comparisons.jsonl').exists()):
            print(f'Auditing {folder.name}', flush=True)
            reports.append(audit_run(folder, kind))
    report = {'schema': SCHEMA, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
              'verifier': str(verifier.resolve()), 'runs': reports,
              'module_sources': [source(Path(m.__file__)) for m in (adjudication, reconstruction, classification, design_audit,
                                                                    individuals)]}
    lines = ['# Saved model decision audit', '',
             'Raw response text was re-parsed with the source modules and reconciled with their decision functions. '
             'No requests were sent. The verifier directory was read only.', '',
             'Reconstruction comparison ids follow each study\'s saved comparison order; selected arms use source GEO metadata. '
             'Exported rows inherit their comparison\'s receipts. Classification ids come from saved answer ids '
             'or an atlas whose SHA256 matches the run summary, using the module\'s input order and eligibility; '
             'design-audit ids come from candidates.csv. Missing ids are unresolved, never guessed. '
             'Individuals-audit records carry their arms and are reconciled from their two answers.', '',
             'Receipts match both request hash and requested model. Status distinguishes retries. For multiple raw responses '
             'with the same status, adjudication judgments identify a matching freshly parsed response when possible; '
             'otherwise the last matching ledger entry is used. Ambiguous receipts are counted separately. '
             'Identical means the fields match a stored response, not that the receipt uniquely identifies a ledger line. '
             'A parse failure can still reproduce an incomplete decision. Missing raw responses prevent re-derivation. '
             'Metadata-only records are counted separately. Counts of receipts refer to comparison/decision receipt occurrences, '
             'not unique requests; exported rows do not recount inherited receipts.', '',
             '| Run | Records | Identical | Differing | Missing records | Unparseable receipts | Unresolved | Non-model |',
             '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in reports:
        lines.append(f'| {r["run"]} | {r["records"]} | {r["rederived_identically"]} | {r["differing"]} | '
                     f'{r["records_missing_receipts"]} | {r["unparseable_responses"]} | {r["unresolved_records"]} | {r["non_model_records"]} |')
    lines += ['', 'audit.json records source paths, timestamps and SHA256 hashes, per-file counts, and up to 20 examples '
              'per finding category per run. Reconciliation uses the current source modules; historical logic changes '
              'may therefore appear as differences.', '']
    try:
        out.mkdir(parents=True, exist_ok=False)
        (out / 'audit.json').write_text(json.dumps(report, indent=2) + '\n')
        (out / 'README.md').write_text('\n'.join(lines))
    except OSError as error:
        raise OSError(f'Cannot write audit to {out}: {error}') from error
    print('\n'.join(lines))
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--verifier', type=Path, default=VERIFIER_ROOT)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        run(args.verifier, args.out)
    except (OSError, ValueError) as error:
        parser.exit(1, f'{error}\n')


if __name__ == '__main__':
    main()
