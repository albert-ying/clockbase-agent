"""Two audited models propose literal rules; only identical sample selections are accepted."""
import json
import os
import threading

import pandas as pd

from audit_system.label_adjudication_2026_09_27 import Ledger, call


def parse_json(text):
    try:
        value = json.loads(text)
        return value if isinstance(value, dict) else None
    except (TypeError, ValueError):
        return None


class Adjudicator:
    def __init__(self, ledger_path, models=None, max_requests=0, cache_only=True, transport=None):
        self.models = tuple(models or os.environ.get('BENCH_GROUP_MODELS', 'gp-6-as,cl-op-5').split(','))
        if len(self.models) != 2 or len(set(self.models)) != 2 or not all(self.models):
            raise ValueError('BENCH_GROUP_MODELS must contain two distinct model identifiers')
        self.ledger = Ledger(ledger_path)
        self.budget = {'max': max_requests, 'sent': 0, 'lock': threading.Lock()}
        self.cache_only, self.transport = cache_only, transport or call

    def ask(self, prompt):
        responses = []
        for model in self.models:
            try:
                response, cached = self.transport(model, prompt, self.ledger, self.budget, parse=parse_json,
                                                  max_tokens=2000, cache_only=self.cache_only)
            except Exception as error:
                response, cached = {'status': 'CALL_ERROR', 'error': f'{type(error).__name__}: {error}'}, False
            if response.get('status') in ('NOT_CACHED', 'BUDGET_STOP', 'CALL_ERROR'):
                self.ledger.append({**response, 'model_requested': model})
            responses.append({**response, 'model_requested': model, 'cached': cached})
        return responses


def apply_rules(meta, rules):
    """AND fields, OR literal allowed values within a field; no expression evaluation."""
    if not isinstance(rules, dict) or not rules:
        raise ValueError('Each arm needs at least one literal field rule')
    selected = pd.Series(True, index=meta.index)
    for field, allowed in rules.items():
        if field not in meta:
            raise ValueError(f'Unknown metadata field: {field}')
        if not isinstance(allowed, list) or not allowed or not all(isinstance(v, str) for v in allowed):
            raise ValueError(f'Allowed values for {field} must be a nonempty list of strings')
        unknown = set(allowed) - set(meta[field].dropna().astype(str))
        if unknown:
            raise ValueError(f'Unknown literal values for {field}: {sorted(unknown)}')
        selected &= meta[field].isin(allowed)
    return sorted(map(str, meta.index[selected]))


def validate_groups(case, control, meta):
    if not isinstance(case, list) or not isinstance(control, list):
        raise ValueError('Manual arms must be lists of GSM identifiers')
    if not all(isinstance(g, str) for g in case + control):
        raise ValueError('Sample identifiers must be strings')
    if not case or not control:
        raise ValueError('An arm is empty')
    if len(set(case)) != len(case) or len(set(control)) != len(control):
        raise ValueError('An arm contains duplicated sample identifiers')
    if set(case) & set(control):
        raise ValueError('Case and control overlap')
    missing = (set(case) | set(control)) - set(meta.index)
    if missing:
        raise ValueError(f'GSMs absent from metadata: {sorted(missing)}')
    return {'case_gsms': sorted(case), 'control_gsms': sorted(control)}


def assign_groups(record, meta, adjudicator, manual=None):
    if manual is not None:
        try:
            selected = validate_groups(manual.get('case_gsms'), manual.get('control_gsms'), meta)
            return {'status': 'MANUAL', **selected, 'provenance': manual.get('provenance')}
        except (ValueError, AttributeError) as error:
            return {'status': 'UNRESOLVED', 'error': str(error), 'manual': manual}
    if meta.empty:
        return {'status': 'UNRESOLVED', 'error': 'Sample metadata unavailable'}
    fields = {c: sorted(meta[c].dropna().astype(str).unique().tolist()) for c in sorted(meta.columns)}
    prompt = ('Assign the published contrast using only the literal GEO metadata below. '
              'Return JSON with feasible (boolean), case and control (objects mapping field names to '
              'nonempty lists of literal allowed values), and note. Rules use AND across fields and OR '
              'within each list. No ranges, substrings, inferred ages or invented fields. '
              'Sample counts in the definitions describe the published analysis and can differ from the '
              'deposited series; do not reject rules because of them, and report any difference in note. '
              'Return feasible=false only if the metadata cannot identify which samples belong to each group.\n'
              + json.dumps(
                  {'gse': record['gse'], 'case_definition': record['case_definition'],
                   'control_definition': record['control_definition'], 'fields': fields}, sort_keys=True))
    answers, selections = adjudicator.ask(prompt), []
    for answer in answers:
        rules = answer.get('parsed') or {}
        try:
            if answer.get('status') != 'OK' or rules.get('feasible') is not True:
                raise ValueError('No feasible rules returned')
            selected = validate_groups(apply_rules(meta, rules.get('case')), apply_rules(meta, rules.get('control')), meta)
            selections.append({'status': 'OK', 'rules': rules, **selected})
        except ValueError as error:
            selections.append({'status': 'UNRESOLVED', 'rules': rules, 'error': str(error)})
    agreed = all(s['status'] == 'OK' for s in selections) and all(
        selections[0][arm] == selections[1][arm] for arm in ('case_gsms', 'control_gsms'))
    return {'status': 'AGREED' if agreed else 'UNRESOLVED', 'answers': answers, 'selections': selections,
            **({arm: selections[0][arm] for arm in ('case_gsms', 'control_gsms')} if agreed else {})}
