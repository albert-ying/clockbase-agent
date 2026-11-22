"""Execute S4 through the verifier's audited client and ledger."""
import hashlib
import json
import logging
import os
import threading
from pathlib import Path

from audit_system import label_adjudication_2026_09_27 as adjudication
from audit_system import paths
from .answers import parse_answer
from .prompts import StudyDescriptions, build_prompt
from .scorer import export_predictions, read_jsonl


SCHEMA = 'clockbase-s4-no-data-v1'
DEFAULT_MODELS = 'gp-6-as,cl-op-5,gp-4o,deepseek-v4-pro,qwen3.8-max,ge-3.1-pr-p'
MAX_TOKENS = 400
SAMPLING = {'gateway': f'temperature 0, max_tokens {MAX_TOKENS}',
            'openai:': ('no temperature sent (provider default; reasoning models accept none), '
                        f'max_completion_tokens {max(4 * MAX_TOKENS, 16000)} including reasoning')}
OUTPUT_FILES = ('manifest.json', 'requests.jsonl', 'ledger.jsonl', 'predictions.jsonl',
                'cost_summary.json', 's4_predictions.jsonl')


def request_body(model, prompt):
    return {'model': model, 'messages': [{'role': 'user', 'content': prompt}],
            'max_tokens': MAX_TOKENS, 'temperature': 0}


def request_hash(model, prompt):
    """Match the audited client's canonical request hash exactly."""
    return hashlib.sha256(json.dumps(request_body(model, prompt), sort_keys=True).encode()).hexdigest()


class S4Ledger(adjudication.Ledger):
    def __init__(self, path):
        self.path, self.lock, self.cache = Path(path), threading.Lock(), {}
        # Returned malformed answers stay cached, so a restart cannot buy a different answer.
        if self.path.exists():
            for record in read_jsonl(self.path):
                if record.get('status') in {'OK', 'UNPARSEABLE'}:
                    self.cache[record['request_sha256']] = record

    def append(self, record):
        super().append(record)
        if record.get('status') == 'UNPARSEABLE':
            self.cache[record['request_sha256']] = record


def _append(path, row):
    try:
        with Path(path).open('a') as handle:
            handle.write(json.dumps(row, allow_nan=False) + '\n')
    except (OSError, ValueError) as error:
        raise RuntimeError(f'Cannot append S4 record to {path}: {error}') from error


def _write_json(path, value):
    try:
        Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    except (OSError, ValueError) as error:
        raise RuntimeError(f'Cannot write S4 summary {path}: {error}') from error


def usage_summary(ledger_path):
    """Count ledger usage once per response, never once per cache replay."""
    per_model = {}
    if Path(ledger_path).exists():
        for record in read_jsonl(ledger_path):
            model = per_model.setdefault(record['model_requested'], {
                'requests': 0, 'attempts': 0, 'failures': 0, 'prompt_tokens': 0,
                'completion_tokens': 0, 'total_tokens': 0, 'usage_missing': 0, 'models_returned': []})
            model['requests'] += 1
            model['attempts'] += record.get('attempts', 0)
            model['failures'] += record.get('status') != 'OK'
            usage = record.get('usage') or {}
            model['usage_missing'] += not bool(usage)
            prompt = usage.get('prompt_tokens', usage.get('input_tokens')) or 0
            completion = usage.get('completion_tokens', usage.get('output_tokens')) or 0
            model['prompt_tokens'] += prompt
            model['completion_tokens'] += completion
            model['total_tokens'] += usage.get('total_tokens') or (prompt + completion)
            returned = record.get('model_returned')
            if returned and returned not in model['models_returned']:
                model['models_returned'].append(returned)
    return {'schema': SCHEMA, 'usage_by_model': per_model,
            'cost_currency': None, 'cost': None,
            'cost_note': 'Token usage only; no provider price schedule was supplied.'}


def _prepare_output(out, ground_truth):
    try:
        digest = hashlib.sha256(ground_truth.read_bytes()).hexdigest()
        out.mkdir(parents=True, exist_ok=True)
        manifest = {'schema': SCHEMA, 'ground_truth_sha256': digest,
                    'ground_truth': paths.portable(ground_truth),
                    'ground_truth_mtime_ns': ground_truth.stat().st_mtime_ns,
                    'max_completion_tokens': MAX_TOKENS, 'sampling': SAMPLING}
        path = out / 'manifest.json'
        if path.exists():
            saved = json.loads(path.read_text())
            if saved != manifest:
                raise ValueError(f'S4 output identity mismatch: {path}; use another output directory')
        else:
            occupied = [name for name in OUTPUT_FILES if (out / name).exists()]
            if occupied:
                raise FileExistsError(f'Refusing to replace unowned S4 outputs in {out}: {occupied}')
            _write_json(path, manifest)
    except OSError as error:
        raise RuntimeError(f'Cannot prepare S4 output {out}: {error}') from error


def run(ground_truth, out, models=None, limit=None, max_requests=None, descriptions=None):
    """Run or replay record/model requests; persist explicit incomplete predictions."""
    ground_truth, out = Path(ground_truth), Path(out)
    records = list(read_jsonl(ground_truth))
    ids = [row['id'] for row in records]
    if len(ids) != len(set(ids)):
        raise ValueError('Ground-truth record ids must be unique')
    if limit is not None:
        if limit < 0:
            raise ValueError('S4 limit must be nonnegative')
        records = records[:limit]
    if models is None:
        models = os.environ.get('BENCH_S4_MODELS', DEFAULT_MODELS)
    if isinstance(models, str):
        models = models.split(',')
    models = tuple(dict.fromkeys(model.strip() for model in models if model.strip()))
    if not models:
        raise ValueError('At least one S4 model is required')
    # Gateway names are sent at temperature 0. The shared client's direct OpenAI branch (openai: names) sends no
    # temperature, because OpenAI reasoning models accept none; those answers are sampled once and replayed from the
    # ledger, never re-requested. The manifest records both settings.
    max_requests = int(os.environ.get('BENCH_S4_MAX_REQUESTS', '3000')) if max_requests is None else max_requests
    if max_requests < 0:
        raise ValueError('S4 request budget must be nonnegative')
    _prepare_output(out, ground_truth)
    ledger_path, predictions = out / 'ledger.jsonl', out / 'predictions.jsonl'
    try:
        ledger = S4Ledger(ledger_path)
    except (OSError, ValueError) as error:
        raise RuntimeError(f'Cannot load S4 ledger {ledger_path}: {error}') from error
    previous = list(read_jsonl(ledger_path)) if ledger_path.exists() else []
    budget = {'sent': len(previous), 'max': max_requests, 'lock': threading.Lock()}
    request_path = out / 'requests.jsonl'
    known_requests = {row['request_hash'] for row in read_jsonl(request_path)} if request_path.exists() else set()
    descriptions = descriptions if descriptions is not None else StudyDescriptions()
    statuses, cached_count = {}, 0
    for record in records:
        # Missing empirical GEO inputs stop the run; no synthetic substitute is constructed.
        study = descriptions.get(record)
        try:
            prepared = build_prompt(record, study)
            preparation_error = None
        except ValueError as error:
            prepared = {'flags': sorted(set(getattr(error, 'flags', []) + ['prompt_preparation_failed'])),
                        'sources': study.get('sources', [])}
            preparation_error = str(error)
        for model in models:
            digest = None
            cached = False
            if preparation_error:
                receipt = {'status': 'PROMPT_FAILED', 'error': preparation_error}
            else:
                prompt = prepared['prompt']
                digest = request_hash(model, prompt)
                if digest not in known_requests:
                    _append(request_path, {'schema': SCHEMA, 'request_hash': digest,
                                          'request': request_body(model, prompt),
                                          'prompt_flags': prepared['flags'],
                                          'blinding_issues': prepared['blinding_issues'],
                                          'sources': prepared['sources'],
                                          'estimated_prompt_tokens': prepared['estimated_prompt_tokens']})
                    known_requests.add(digest)
                try:
                    receipt, cached = adjudication.call(model, prompt, ledger, budget,
                                                       parse=parse_answer, max_tokens=MAX_TOKENS)
                except Exception as error:
                    # The shared client can fail before transport (missing credentials) or on payload structure.
                    receipt = {'status': 'CLIENT_ERROR', 'request_sha256': digest,
                               'model_requested': model, 'error': f'{type(error).__name__}: {error}'[:300]}
                    ledger.append(receipt)
                    logging.error('S4 client failed for %s / %s: %s', record['id'], model, receipt['error'])
            parsed = parse_answer(receipt.get('response_text')) if receipt.get('status') == 'OK' else None
            status = receipt['status']
            if status == 'OK' and parsed is None:
                status = 'UNPARSEABLE'
            result = {'record_id': record['id'], 'model': model, 'status': status,
                      'prediction': None, 'significance': None, 'confidence': None,
                      'recognizes_publication': None, 'publication': None,
                      'request_hash': digest, 'prompt_flags': prepared['flags'], 'cached': cached,
                      'raw_response_reference': f'ledger.jsonl#request_sha256={digest}'
                      if digest and status not in {'BUDGET_STOP', 'NOT_CACHED'} else None}
            if parsed:
                result.update(parsed)
            if receipt.get('error'):
                result['error'] = receipt['error']
            _append(predictions, result)
            statuses[status] = statuses.get(status, 0) + 1
            cached_count += cached
            logging.info('%s / %s: %s%s', record['id'], model, status, ' (cached)' if cached else '')
    if not predictions.exists():
        predictions.touch()
    scored = export_predictions(predictions, out / 's4_predictions.jsonl')
    summary = usage_summary(ledger_path)
    summary.update(models=list(models), records=len(records), statuses=statuses, cached_responses=cached_count,
                   budget_max_requests=max_requests, requests_preexisting=len(previous),
                   requests_sent_this_run=budget['sent'] - len(previous),
                   cache_only=os.environ.get('AUDIT_CACHE_ONLY') == '1', engine_rows=len(scored))
    _write_json(out / 'cost_summary.json', summary)
    return summary
