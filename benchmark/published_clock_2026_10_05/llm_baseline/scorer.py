"""Export normalized S4 predictions for the benchmark engine."""
import json
from pathlib import Path


def read_jsonl(path):
    path = Path(path)
    if not path.exists():
        raise RuntimeError(f'Cannot proceed - required data file not found: {path}')
    try:
        with path.open() as handle:
            for ordinal, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                try:
                    value = json.loads(line)
                except ValueError as error:
                    raise ValueError(f'Invalid JSON at {path}:{ordinal}: {error}') from error
                if not isinstance(value, dict):
                    raise ValueError(f'Expected an object at {path}:{ordinal}')
                yield value
    except OSError as error:
        raise RuntimeError(f'Cannot read {path}: {error}') from error


def export_predictions(predictions, output):
    """Keep the latest attempt per record/model, including explicit failed rows."""
    latest = {}
    for row in read_jsonl(predictions):
        latest[row['record_id'], row['model']] = row
    rows = []
    for row in latest.values():
        ok = row.get('status') == 'OK'
        rows.append({'ground_truth_id': row['record_id'], 'model': row['model'], 'system': 'S4',
                     'status': row['status'], 'direction': row['prediction'] if ok else None,
                     'significant': row['significance'] if ok else None,
                     'confidence': row['confidence'] if ok else None,
                     'recognizes_publication': row['recognizes_publication'] if ok else None,
                     'publication': row.get('publication') if ok else None,
                     'request_hash': row.get('request_hash'),
                     'raw_response_reference': row.get('raw_response_reference'),
                     'prompt_flags': row.get('prompt_flags', [])})
    output = Path(output)
    try:
        with output.open('w') as handle:
            for row in rows:
                handle.write(json.dumps(row, allow_nan=False) + '\n')
    except OSError as error:
        raise RuntimeError(f'Cannot write S4 engine input {output}: {error}') from error
    return rows
