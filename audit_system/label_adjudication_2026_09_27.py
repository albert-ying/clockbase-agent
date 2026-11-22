"""Dual-model adjudication of case/control labels for recovered atlas comparisons.

For each row with a recovered membership, two models from different families (NeoLink gp-6-as =
GPT-6 Astra, cl-op-5 = Claude Opus 5) independently read the literal GEO sample fields of the case
and control arms and judge whether the case arm is the labelled condition, whether the control is
appropriate, and whether the stored effect direction matches the label. Only agreement decides:

  both CORRECT  -> LABEL_SUPPORTED          both INVERTED -> LABEL_INVERTED_CONFIRMED (sign repair)
  otherwise     -> LABEL_UNCONFIRMED

Every request/response, token usage and failure is kept in a ledger. Identical requests are served
from the ledger instead of being re-sent.

  python -m audit_system.label_adjudication_2026_09_27 --modality dnam --membership DIR --out DIR \
      [--select inverted|all|sample] [--limit N] [--workers 8] [--max-requests N]

Run with the provider API key in the environment: python -m ...
"""
import argparse
import hashlib
import json
import os
import re
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from audit_system.adjudication_files import decision_files
from audit_system.membership_recovery_2026_09_27 import GEO, MODALITIES, sha256
from audit_system.row_verdicts_2026_09_27 import label_orientation

SCHEMA = 'clockbase-label-adjudication-v1'
# Two independent models; AUDIT_MODELS overrides the NeoLink pair (e.g. 'openai:gpt-6-astra,openai:gpt-6-sol').
MODELS = tuple(os.environ.get('AUDIT_MODELS', 'gp-6-as,cl-op-5').split(','))
MAX_SAMPLES_SHOWN = 12
FIELDS_SKIPPED = {'GSM'}
ISSUE_CODES = ['CASE_NOT_CONDITION', 'CONTROL_INAPPROPRIATE', 'TISSUE_OR_CELLTYPE_CONFOUNDED', 'TIMEPOINT_CONFOUNDED',
               'TECHNICAL_REPLICATES', 'NON_STANDARD_ASSAY', 'CELL_LINE_ONLY', 'MIXED_PLATFORM_OR_BATCH', 'OTHER']

PROMPT = """You audit one comparison from an aging-clock atlas built by an automated agent.
The atlas row claims: comparison "{comparison}", condition (case) "{condition}", clock {clock}.
The stored effect is {direction}: the case arm has a {higher} mean clock value than the control arm.

Below are the literal GEO sample fields of the samples the agent actually used. Arm A is the arm the
stored effect treats as the case; arm B is the control.

ARM A ({n_a} samples):
{arm_a}

ARM B ({n_b} samples):
{arm_b}

Judge only from these fields. Return JSON only, with keys:
  "arm_a_is_condition": "YES" | "NO" | "UNCLEAR"   (arm A matches the labelled case condition)
  "arm_b_is_appropriate_control": "YES" | "NO" | "UNCLEAR"
  "orientation": "CORRECT" | "INVERTED" | "UNCLEAR"   (INVERTED if arm B, not arm A, is the labelled condition)
  "issues": list drawn from {codes}
  "evidence": up to 3 short verbatim field values supporting your judgment
  "rationale": one sentence"""


def sample_fields(modality, gse):
    folder = MODALITIES[modality]['folder']
    path = GEO / folder / 'processed_metadata_preharmonize' / f'{gse}.csv'
    if not path.exists():
        path = GEO / folder / 'processed_metadata' / f'{gse}.csv'
    if not path.exists():
        return None
    frame = pd.read_csv(path, dtype=str)
    key = 'GSM' if 'GSM' in frame.columns else frame.columns[0]
    return frame.set_index(key)


def describe_arm(meta, gsms):
    rows = meta.reindex(gsms)
    columns = [c for c in rows.columns if c not in FIELDS_SKIPPED and rows[c].notna().any()]
    lines = []
    for gsm, row in rows.head(MAX_SAMPLES_SHOWN).iterrows():
        values = '; '.join(f'{c}={str(row[c])[:80]}' for c in columns if pd.notna(row[c]))
        lines.append(f'- {gsm}: {values}')
    if len(rows) > MAX_SAMPLES_SHOWN:
        lines.append(f'- ... {len(rows) - MAX_SAMPLES_SHOWN} more samples')
    return '\n'.join(lines)


class Ledger:
    def __init__(self, path):
        self.path, self.lock, self.cache = Path(path), threading.Lock(), {}
        if self.path.exists():
            for line in open(self.path):
                record = json.loads(line)
                if record.get('status') == 'OK':
                    self.cache[record['request_sha256']] = record

    def append(self, record):
        with self.lock, open(self.path, 'a') as handle:
            handle.write(json.dumps(record, default=str) + '\n')
            if record.get('status') == 'OK':
                self.cache[record['request_sha256']] = record


def parse_judgment(text):
    """Strict JSON first; otherwise read the enumerated fields literally (verbatim evidence may break JSON)."""
    match = re.search(r'\{.*\}', text, re.S)
    if match:
        try:
            return json.loads(match.group(0))
        except json.JSONDecodeError:
            pass
    fields = {}
    for key, options in (('arm_a_is_condition', 'YES|NO|UNCLEAR'), ('arm_b_is_appropriate_control', 'YES|NO|UNCLEAR'),
                         ('orientation', 'CORRECT|INVERTED|UNCLEAR')):
        found = re.search(rf'"{key}"\s*:\s*"({options})"', text)
        if not found:
            return None
        fields[key] = found.group(1)
    fields['issues'] = [c for c in ISSUE_CODES if f'"{c}"' in text]
    fields['parse'] = 'field_extraction'
    return fields


def call(model, prompt, ledger, budget, parse=None, max_tokens=800, cache_only=False):
    parse = parse or parse_judgment
    body = {'model': model, 'messages': [{'role': 'user', 'content': prompt}], 'max_tokens': max_tokens, 'temperature': 0}
    digest = hashlib.sha256(json.dumps(body, sort_keys=True).encode()).hexdigest()
    if digest in ledger.cache:
        return ledger.cache[digest], True
    if cache_only or os.environ.get('AUDIT_CACHE_ONLY') == '1':  # replay from the ledger only; never send
        return {'status': 'NOT_CACHED', 'request_sha256': digest}, False
    with budget['lock']:
        if budget['sent'] >= budget['max']:
            return {'status': 'BUDGET_STOP'}, False
        budget['sent'] += 1
    record = {'request_sha256': digest, 'model_requested': model, 'time': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    if model.startswith('openai:'):
        # OpenAI API directly: reasoning models take max_completion_tokens (reasoning counts against it) and no temperature.
        sent = {'model': model.split(':', 1)[1], 'messages': body['messages'], 'max_completion_tokens': max(4 * max_tokens, 16000)}
        url, key = 'https://api.openai.com/v1/chat/completions', os.environ['OPENAI_API_KEY']
    else:
        sent, url, key = body, os.environ['TOKENSWITCH_BASE_URL'].rstrip('/') + '/v1/chat/completions', os.environ['TOKENSWITCH_API_KEY']
    request = urllib.request.Request(url, data=json.dumps(sent).encode(), method='POST',
                                     headers={'Authorization': 'Bearer ' + key, 'Content-Type': 'application/json'})
    for attempt in range(3):  # retries cover transport failures only; a returned answer is never re-requested
        try:
            with urllib.request.urlopen(request, timeout=600) as response:
                payload = json.loads(response.read())
        except Exception as error:
            record.update(status='TRANSPORT_ERROR', error=f'{type(error).__name__}: {error}'[:300], attempts=attempt + 1)
            time.sleep(5 * (attempt + 1))
            continue
        text = payload['choices'][0]['message']['content']
        parsed = parse(text)
        record.update(status='OK' if parsed else 'UNPARSEABLE', model_returned=payload.get('model'),
                      usage=payload.get('usage'), response_text=text, parsed=parsed, attempts=attempt + 1)
        record.pop('error', None)
        break
    ledger.append(record)
    return record, False


def decide(judgments):
    orientations = [j.get('orientation') for j in judgments]
    if all(o == 'CORRECT' for o in orientations):
        return 'LABEL_SUPPORTED'
    if all(o == 'INVERTED' for o in orientations):
        return 'LABEL_INVERTED_CONFIRMED'
    return 'LABEL_UNCONFIRMED'


LABEL_CODES = {'LABEL_UNCONFIRMED', 'LABEL_ORIENTATION_INVERTED'}
HARD_ELSEWHERE = {'ARM_TOO_SMALL', 'EFFECT_NOT_REPRODUCED', 'INPUT_EMPTY', 'NON_STANDARD_ASSAY', 'FDR_BELOW_RAW_P',
                  'RNA_INPUT_ARTIFACT', 'COVERAGE_UNAVAILABLE', 'CLOCK_COVERAGE_LOW', 'CLOCK_UNREGISTERED'}


def run(modality, membership_dir, out_dir, select='inverted', limit=None, workers=8, max_requests=5000, verdicts_dir=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = pd.read_json(Path(membership_dir) / 'rows.jsonl', lines=True)
    rows = rows[rows.status == 'UNIQUE'].copy()
    for column in ('case_level', 'control_level', 'stratum'):
        if column not in rows.columns:
            rows[column] = None
    rows['heuristic'] = [label_orientation(r.comparison_name, r.condition_name, r.case_level, r.control_level, r.stratum)
                         if r.case_level is not None else 'UNCONFIRMED' for r in rows.itertuples()]
    if select == 'unconfirmed':
        # Usable rows still carrying LABEL_UNCONFIRMED: adjudication can confirm them or reveal an inversion.
        repairs = pd.read_json(Path(verdicts_dir) / 'repairs.jsonl', lines=True)
        wanted = repairs[repairs.final_class.isin(['VERIFIED', 'PROVISIONAL'])
                         & repairs.soft_flags.apply(lambda f: isinstance(f, list) and 'LABEL_UNCONFIRMED' in f)]
        rows = rows[rows.row_ordinal.isin(set(wanted.row_ordinal))]
    elif select == 'pending':
        # Rows the repair stage left waiting for an orientation decision.
        repairs = pd.read_json(Path(verdicts_dir) / 'repairs.jsonl', lines=True)
        rows = rows[rows.row_ordinal.isin(set(repairs[repairs.final_class == 'PENDING_ADJUDICATION'].row_ordinal))]
    elif select == 'decisive':
        # Only rows whose verdict the label decides: a label code present and no other hard failure.
        verdicts = pd.read_json(Path(verdicts_dir) / 'row_verdicts.jsonl', lines=True)
        decisive = verdicts[verdicts.reasons.apply(lambda r: bool(LABEL_CODES & set(r)) and not (HARD_ELSEWHERE & set(r)))]
        rows = rows[rows.row_ordinal.isin(set(decisive.row_ordinal))]
    elif select == 'inverted':
        rows = rows[rows.heuristic == 'INVERTED']
    elif select == 'sample':
        rows = rows.sample(min(len(rows), limit or 200), random_state=11)
    elif select == 'supported_sample':
        # Validation sample: usable rows whose orientation rests on the label-token heuristic alone (heuristic
        # SUPPORTED, never adjudicated). Both models' agreement estimates the heuristic's error rate.
        repairs = pd.read_json(Path(verdicts_dir) / 'repairs.jsonl', lines=True)
        usable = set(repairs[repairs.final_class.isin(['VERIFIED', 'PROVISIONAL'])].row_ordinal)
        done = {json.loads(line)['row_ordinal'] for path in decision_files(out_dir) for line in open(path)}
        pool = rows[(rows.heuristic == 'SUPPORTED') & rows.row_ordinal.isin(usable) & ~rows.row_ordinal.isin(done)
                    & (rows.get('source', pd.Series('', index=rows.index)) != 'design_reconstruction')]
        rows = pool.sample(min(len(pool), limit or 300), random_state=20260930)
        limit = None
    if limit:
        rows = rows.head(limit)
    ledger = Ledger(out_dir / 'ledger.jsonl')
    budget = {'max': max_requests, 'sent': 0, 'lock': threading.Lock()}
    metadata_cache, lock = {}, threading.Lock()

    def adjudicate(row):
        with lock:
            if row.gse_id not in metadata_cache:
                metadata_cache[row.gse_id] = sample_fields(modality, row.gse_id)
            meta = metadata_cache[row.gse_id]
        if meta is None:
            return {'row_ordinal': int(row.row_ordinal), 'decision': 'METADATA_UNAVAILABLE'}
        prompt = PROMPT.format(comparison=row.comparison_name, condition=row.condition_name, clock=row.clock,
                               direction='positive' if row.stored_effect > 0 else 'negative',
                               higher='higher' if row.stored_effect > 0 else 'lower',
                               n_a=len(row.case_gsms), n_b=len(row.control_gsms),
                               arm_a=describe_arm(meta, row.case_gsms), arm_b=describe_arm(meta, row.control_gsms),
                               codes=ISSUE_CODES)
        judgments, receipts = [], []
        for model in MODELS:
            record, cached = call(model, prompt, ledger, budget)
            receipts.append({'model': model, 'status': record.get('status'), 'cached': cached,
                             'request_sha256': record.get('request_sha256')})
            if record.get('status') == 'OK':
                judgments.append(record['parsed'])
        decision = decide(judgments) if len(judgments) == len(MODELS) else 'ADJUDICATION_INCOMPLETE'
        return {'row_ordinal': int(row.row_ordinal), 'gse_id': row.gse_id, 'comparison_name': row.comparison_name,
                'condition_name': row.condition_name, 'clock': row.clock, 'heuristic': row.heuristic,
                'decision': decision, 'judgments': judgments, 'receipts': receipts,
                'issues_agreed': sorted(set.intersection(*[set(j.get('issues') or []) for j in judgments])) if len(judgments) == 2 else []}

    results = []
    with ThreadPoolExecutor(workers) as pool:
        for result in pool.map(adjudicate, list(rows.itertuples())):
            results.append(result)
    stamp = time.strftime('%Y%m%dT%H%M%S')
    # One file per run: a later run with the same selection must not erase an earlier run's decisions.
    with open(out_dir / f'decisions-{select}-{stamp}.jsonl', 'w') as sink:
        for result in results:
            sink.write(json.dumps(result, default=str) + '\n')
    usage = {'prompt_tokens': 0, 'completion_tokens': 0}
    per_model = {}
    for line in open(out_dir / 'ledger.jsonl'):
        record = json.loads(line)
        u = record.get('usage') or {}
        m = per_model.setdefault(record['model_requested'], {'requests': 0, 'prompt_tokens': 0, 'completion_tokens': 0, 'failures': 0})
        m['requests'] += 1
        m['failures'] += record.get('status') != 'OK'
        m['prompt_tokens'] += u.get('prompt_tokens') or 0
        m['completion_tokens'] += u.get('completion_tokens') or 0
    summary = {'schema': SCHEMA, 'modality': modality, 'select': select, 'rows': len(results),
               'decisions': pd.Series([r['decision'] for r in results]).value_counts().to_dict(),
               'ledger_usage_by_model': per_model, 'module_sha256': sha256(__file__),
               'membership_sha256': sha256(Path(membership_dir) / 'rows.jsonl'), 'models': MODELS,
               'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / f'summary-{select}-{stamp}.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--modality', choices=sorted(MODALITIES), required=True)
    parser.add_argument('--membership', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--select', choices=['inverted', 'all', 'sample', 'decisive', 'pending', 'unconfirmed',
                                             'supported_sample'], default='inverted')
    parser.add_argument('--verdicts', type=Path, help='row verdicts dir (required for --select decisive)')
    parser.add_argument('--limit', type=int)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--max-requests', type=int, default=5000)
    args = parser.parse_args(argv)
    run(args.modality, args.membership, args.out, args.select, args.limit, args.workers, args.max_requests, args.verdicts)


if __name__ == '__main__':
    main()
