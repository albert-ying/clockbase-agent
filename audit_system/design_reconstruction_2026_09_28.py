"""Reconstruct case/control groups for atlas comparisons whose historical membership is unrecoverable.

One request per study and model lists that study's pending comparisons and a summary of its
sample metadata (every field with its values and sample counts). Each model returns, for every
comparison, a selection rule per arm: {field: [allowed values]}. The rules are applied to the
metadata deterministically to obtain GSM sets. A comparison is reconstructed only when GPT-6 Astra
and Claude Opus 5 produce identical case and control GSM sets; otherwise the disagreement is kept.
Output rows use the merged-membership schema (source = design_reconstruction), so row repairs can
recompute them with the declared estimand.

  with-key python -m audit_system.design_reconstruction_2026_09_28 --modality dnam \
      --repairs DIR --out DIR [--limit-studies N] [--workers 8] [--max-requests N]
"""
import argparse
import hashlib
import json
import os
import re
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from audit_system.label_adjudication_2026_09_27 import Ledger, call as ledger_call
from audit_system.membership_recovery_2026_09_27 import GEO, MODALITIES, sha256

SCHEMA = 'clockbase-design-reconstruction-v1'
MODELS = tuple(os.environ.get('AUDIT_MODELS', 'gp-6-as,cl-op-5').split(','))  # see label_adjudication_2026_09_27
MAX_LEVELS_LISTED, MAX_EXAMPLES = 60, 25
BATCH = 10

PROMPT = """You reconstruct the sample groups of comparisons that an automated agent reported for GEO study
{gse} ({n} samples). Its sample metadata fields, with each value and its sample count, are:

{fields}

For each comparison below, define the CASE arm (the labelled condition) and the CONTROL arm as selection
rules over these fields. A rule is an object mapping field names to lists of allowed values, copied
exactly as shown; a sample belongs to an arm when it matches every field in the rule. Use the most
specific appropriate control (same tissue, cell type, time point and other conditions as the case where
the metadata allow). If the comparison cannot be formed from these fields, set "feasible": false.

Comparisons:
{comparisons}

Return JSON only: {{"comparisons": [{{"id": <id>, "feasible": true|false, "case": {{...}}, "control": {{...}},
"note": "<one sentence>"}}]}}"""


def load_meta(modality, gse):
    folder = MODALITIES[modality]['folder']
    for sub in ('processed_metadata_preharmonize', 'processed_metadata'):
        path = GEO / folder / sub / f'{gse}.csv'
        if path.exists():
            frame = pd.read_csv(path, dtype=str)
            key = 'GSM' if 'GSM' in frame.columns else frame.columns[0]
            return frame.drop_duplicates(key).set_index(key)
    return None


def field_summary(meta):
    lines = []
    for column in meta.columns:
        counts = meta[column].fillna('NA').value_counts()
        if len(counts) <= 1:
            continue
        if len(counts) <= MAX_LEVELS_LISTED:
            shown = '; '.join(f'"{v}" ({n})' for v, n in counts.items())
        else:
            shown = '; '.join(f'"{v}"' for v in counts.index[:MAX_EXAMPLES]) + f'; ... ({len(counts)} distinct values)'
        lines.append(f'- {column}: {shown}')
    return '\n'.join(lines)


def apply_rule(meta, rule):
    if not isinstance(rule, dict) or not rule:
        return None
    mask = pd.Series(True, index=meta.index)
    for field, allowed in rule.items():
        if field not in meta.columns:
            return None
        allowed = allowed if isinstance(allowed, list) else [allowed]
        mask &= meta[field].fillna('NA').isin([str(v) for v in allowed])
    return sorted(meta.index[mask])


def parse_rules(text):
    match = re.search(r'\{.*\}', text or '', re.S)
    if not match:
        return None
    try:
        return {c['id']: c for c in json.loads(match.group(0)).get('comparisons', []) if 'id' in c}
    except (json.JSONDecodeError, TypeError, AttributeError):
        return None


def normalize_ids(rules):
    """Comparison ids come back as ints when parsed live and as strings when read from the JSON ledger."""
    if not rules:
        return rules
    return {int(k) if str(k).lstrip('-').isdigit() else k: v for k, v in rules.items()}


def reconcile(meta, rules_by_model, comparison_id):
    sets = []
    for model, rules in rules_by_model.items():
        rule = normalize_ids(rules or {}).get(comparison_id)
        if not rule or not rule.get('feasible', True):
            return {'status': 'INFEASIBLE_OR_MISSING', 'model': model}
        case, control = apply_rule(meta, rule.get('case')), apply_rule(meta, rule.get('control'))
        if not case or not control:
            return {'status': 'RULE_SELECTS_NO_SAMPLES', 'model': model}
        if set(case) & set(control):
            return {'status': 'ARMS_OVERLAP', 'model': model}
        sets.append((case, control, rule))
    if len(sets) == 2 and sets[0][:2] == sets[1][:2]:
        return {'status': 'AGREED', 'case_gsms': sets[0][0], 'control_gsms': sets[0][1],
                'rules': {m: s[2] for m, s in zip(rules_by_model, sets)}}
    return {'status': 'MODELS_DISAGREE', 'rules': {m: s[2] for m, s in zip(rules_by_model, sets)},
            'sizes': [(len(s[0]), len(s[1])) for s in sets]}


def run(modality, repairs_dir, out_dir, limit_studies=None, workers=8, max_requests=20000):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    repairs = pd.read_json(Path(repairs_dir) / 'repairs.jsonl', lines=True)
    pending = repairs[repairs.final_class == 'PENDING_DESIGN']
    comparisons = (pending.groupby(['gse_id', 'comparison_name', 'condition_name'], dropna=False)
                   .row_ordinal.apply(list).reset_index())
    studies = list(comparisons.groupby('gse_id'))
    if limit_studies:
        studies = studies[:limit_studies]
    ledger = Ledger(out_dir / 'ledger.jsonl')
    budget = {'max': max_requests, 'sent': 0, 'lock': threading.Lock()}

    def study(item):
        gse, frame = item
        meta = load_meta(modality, gse)
        if meta is None or len(meta) < 2:
            return [{'gse_id': gse, 'comparison_name': r.comparison_name, 'condition_name': r.condition_name,
                     'row_ordinals': r.row_ordinal, 'status': 'METADATA_UNAVAILABLE'} for r in frame.itertuples()]
        items = [{'id': i, 'comparison': r.comparison_name, 'condition': r.condition_name}
                 for i, r in enumerate(frame.itertuples())]
        fields = field_summary(meta)

        def prompt_for(batch):
            return PROMPT.format(gse=gse, n=len(meta), fields=fields,
                                 comparisons='\n'.join(f'{c["id"]}. "{c["comparison"]}" (condition: "{c["condition"]}")' for c in batch))

        rules_by_model, receipts = {}, []
        # A completed whole-study answer is reused from the ledger; otherwise request in batches so the
        # answer fits the output limit (whole-study answers for many comparisons were truncated).
        whole = {m: ledger_call(m, prompt_for(items), ledger, budget, parse=parse_rules, max_tokens=4000, cache_only=True)
                 for m in MODELS}
        for model in MODELS:
            record, cached = whole[model]
            if record.get('status') == 'OK':
                receipts.append({'model': model, 'status': 'OK', 'request_sha256': record.get('request_sha256'),
                                 'cached': True, 'batch': 'whole'})
                rules_by_model[model] = record.get('parsed')
                continue
            merged, ok = {}, True
            for start in range(0, len(items), BATCH):
                batch = items[start:start + BATCH]
                record, cached = ledger_call(model, prompt_for(batch), ledger, budget, parse=parse_rules, max_tokens=8000)
                receipts.append({'model': model, 'status': record.get('status'), 'request_sha256': record.get('request_sha256'),
                                 'cached': cached, 'batch': [c['id'] for c in batch]})
                if record.get('status') == 'OK':
                    merged.update(record.get('parsed') or {})
                else:
                    ok = False
            rules_by_model[model] = merged if ok else None
        out = []
        for c, r in zip(items, frame.itertuples()):
            result = reconcile(meta, rules_by_model, c['id']) if all(rules_by_model.values()) else {'status': 'MODEL_CALL_FAILED'}
            out.append({'gse_id': gse, 'comparison_name': r.comparison_name, 'condition_name': r.condition_name,
                        'row_ordinals': r.row_ordinal, 'receipts': receipts, **result})
        return out

    results = []
    with ThreadPoolExecutor(workers) as pool:
        for chunk in pool.map(study, studies):
            results.extend(chunk)
    with open(out_dir / 'comparisons.jsonl', 'w') as sink:
        for r in results:
            sink.write(json.dumps(r, default=str) + '\n')
    with open(out_dir / 'rows.jsonl', 'w') as sink:
        for r in results:
            if r['status'] != 'AGREED':
                continue
            for ordinal in r['row_ordinals']:
                sink.write(json.dumps({'row_ordinal': ordinal, 'status': 'UNIQUE', 'source': 'design_reconstruction',
                                       'kind': 'two_sample', 'case_gsms': r['case_gsms'], 'control_gsms': r['control_gsms'],
                                       'n_case': len(r['case_gsms']), 'n_control': len(r['control_gsms']),
                                       'design_rules': r['rules']}) + '\n')
    summary = {'schema': SCHEMA, 'modality': modality, 'studies': len(studies), 'comparisons': len(results),
               'rows': int(sum(len(r['row_ordinals']) for r in results)),
               'status': pd.Series([r['status'] for r in results]).value_counts().to_dict(),
               'rows_reconstructed': int(sum(len(r['row_ordinals']) for r in results if r['status'] == 'AGREED')),
               'module_sha256': sha256(__file__), 'repairs_sha256': sha256(Path(repairs_dir) / 'repairs.jsonl'),
               'models': MODELS, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--modality', choices=sorted(MODALITIES), required=True)
    parser.add_argument('--repairs', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--limit-studies', type=int)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--max-requests', type=int, default=20000)
    args = parser.parse_args(argv)
    run(args.modality, args.repairs, args.out, args.limit_studies, args.workers, args.max_requests)


if __name__ == '__main__':
    main()
