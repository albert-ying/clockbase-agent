"""Run the published-clock benchmark without modifying verifier or GEO inputs."""
import argparse
import logging
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from .groups import Adjudicator, assign_groups
from .io import read_json, read_jsonl, receipt
from .match import atlas_result, rank_rows
from .metrics import metrics
from .predictions import PredictionReader, atlas_rows, export_path
from .report import write_report
from .stats import evaluate, threshold

from audit_system import membership_recovery_2026_09_27, row_repairs_2026_09_28, label_adjudication_2026_09_27

LOGGER = logging.getLogger(__name__)


def load_s4(path):
    """Accept JSONL predictions keyed by ground_truth_id (or id) and model."""
    predictions = {}
    for record in read_jsonl(path):
        identity = record.get('ground_truth_id', record.get('id'))
        model = record.get('model', record.get('model_requested'))
        if not identity or not model:
            raise ValueError(f'S4 prediction needs ground_truth_id and model: {path}')
        key = (identity, model)
        if key in predictions:
            raise ValueError(f'Duplicate S4 prediction: {key}')
        value = record.get('prediction', record.get('parsed', record))
        if not isinstance(value, dict):
            raise ValueError(f'S4 prediction must be an object: {key}')
        valid = record.get('status', 'OK') == 'OK' and value.get('direction') in ('higher', 'lower', 'none')
        predictions[key] = {**value, 'status': record.get('status', 'OK'), 'covered': valid,
                            'significant': value.get('significant') if isinstance(value.get('significant'), bool) else None}
    return predictions


def validate_records(records):
    seen = set()
    for record in records:
        required = ('id', 'tier', 'gse', 'modality', 'clock_clockbase', 'measure', 'case_definition',
                    'control_definition', 'direction', 'source')
        missing = [k for k in required if k not in record]
        if missing:
            raise ValueError(f'Ground-truth record missing {missing}')
        if record['id'] in seen:
            raise ValueError(f"Duplicate ground-truth id: {record['id']}")
        seen.add(record['id'])
        if record['tier'] not in ('T1', 'T2') or record['direction'] not in ('higher', 'lower', 'none'):
            raise ValueError(f"Invalid tier or direction: {record['id']}")
        if not isinstance(record['source'], dict) or not record['source'].get('quote'):
            raise ValueError(f"Ground-truth record needs its source quote: {record['id']}")


def run(ground_truth, out, groups_manual=None, cache_only=True, s4_predictions=None,
        bootstrap=1000, seed=20261005, max_requests=0):
    ground_truth, out = Path(ground_truth), Path(out)
    records = read_jsonl(ground_truth)
    validate_records(records)
    manual = read_json(groups_manual) if groups_manual else {}
    if not isinstance(manual, dict):
        raise ValueError('Manual groups must be an object keyed by ground-truth id')
    unknown = set(manual) - {r['id'] for r in records}
    if unknown:
        raise ValueError(f'Manual groups have unknown record ids: {sorted(unknown)}')
    if bootstrap < 0 or max_requests < 0:
        raise ValueError('Bootstrap replicates and request budget must be nonnegative')
    if s4_predictions is None:
        candidate = ground_truth.parent / 's4_predictions.jsonl'
        s4_predictions = candidate if candidate.exists() else None
    s4 = load_s4(s4_predictions) if s4_predictions else {}
    cache_only = cache_only or os.environ.get('AUDIT_CACHE_ONLY') == '1'
    try:
        out.mkdir(parents=True, exist_ok=False)
    except OSError as error:
        raise OSError(f'Cannot create fresh output directory {out}: {error}') from error
    adjudicator = Adjudicator(out / 'ledger.jsonl', max_requests=max_requests, cache_only=cache_only)
    # An existing ledger may be supplied read-only for request-hash cache replay.
    cache_path = os.environ.get('BENCH_GROUP_LEDGER')
    if cache_path:
        from audit_system.label_adjudication_2026_09_27 import Ledger
        adjudicator.ledger.cache.update(Ledger(cache_path).cache)
    reader, results = PredictionReader(), []
    for record in records:
        LOGGER.info('Comparing %s', record['id'])
        result = dict(record)
        try:
            rows = atlas_rows(record['modality'], record['gse'])
            prediction = reader.load(record['modality'], record['gse'], record['clock_clockbase'], rows)
            groups = assign_groups(record, prediction['metadata'], adjudicator, manual.get(record['id']))
            match = rank_rows(record, rows, groups, adjudicator)
            selected = next((c for c in match['candidates'] if c['row_ordinal'] == match['selected_row_ordinal']), None)
            systems = {s: atlas_result(record, selected, s) for s in ('S1', 'S2')}
            systems['S3'] = evaluate(prediction['samples'], groups, record['measure'], threshold(record)) if prediction['status'] == 'OK' else {
                'status': prediction['status'], 'covered': False}
            result.update(groups=groups, match=match, systems=systems,
                          predictions={k: v for k, v in prediction.items() if k not in ('samples', 'metadata')})
            path = export_path(record['modality'], record['gse'])
            result['atlas_source'] = receipt(path) if path.exists() else None
        except (OSError, ValueError, KeyError, TypeError) as error:
            LOGGER.warning('%s: %s', record['id'], error)
            result.update(error=f'{type(error).__name__}: {error}', groups={'status': 'UNRESOLVED'},
                          match={'status': 'UNRESOLVED', 'candidates': []},
                          systems={s: {'status': 'RECORD_ERROR', 'covered': False} for s in ('S1', 'S2', 'S3')})
        for (identity, model), prediction in s4.items():
            if identity == record['id']:
                result['systems'][f'S4:{model}'] = prediction
        results.append(result)
    summary = metrics(results, bootstrap, seed)
    manifest = {'schema': 'clockbase-published-clock-engine-v1', 'created': datetime.now(timezone.utc).isoformat(),
                'ground_truth': receipt(ground_truth), 'manual_groups': receipt(groups_manual) if groups_manual else None,
                's4_predictions': receipt(s4_predictions) if s4_predictions else None,
                'cache_ledger': receipt(cache_path) if cache_path else None, 'cache_only': cache_only,
                'models': adjudicator.models, 'request_budget': max_requests, 'requests_sent': adjudicator.budget['sent'],
                'records': len(results), 'bootstrap': summary['bootstrap'],
                'software': {'python': sys.version, 'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
                'verifier_code': [receipt(module.__file__) for module in (
                    membership_recovery_2026_09_27, row_repairs_2026_09_28, label_adjudication_2026_09_27)],
                'code': [receipt(p) for p in sorted(Path(__file__).parent.glob('*.py'))]}
    write_report(out, results, summary, manifest)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    command = commands.add_parser('run')
    command.add_argument('--ground-truth', type=Path, required=True)
    command.add_argument('--out', type=Path, required=True)
    command.add_argument('--groups-manual', type=Path)
    command.add_argument('--cache-only', action='store_true')
    command.add_argument('--s4-predictions', type=Path)
    command.add_argument('--bootstrap', type=int, default=1000)
    command.add_argument('--seed', type=int, default=20261005)
    command.add_argument('--max-requests', type=int, default=0, help='zero prevents live requests by default')
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    run(args.ground_truth, args.out, args.groups_manual, args.cache_only, args.s4_predictions,
        args.bootstrap, args.seed, args.max_requests)


if __name__ == '__main__':
    main()
