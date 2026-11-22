"""JSON inputs and portable, hashed source receipts."""
import json
import math
from datetime import datetime, timezone
from pathlib import Path

from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import sha256


def receipt(path):
    path = Path(path)
    try:
        stat = path.stat()
        return {'path': paths.portable(path), 'sha256': sha256(path), 'bytes': stat.st_size,
                'modified': datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat()}
    except OSError as error:
        raise OSError(f'Cannot read source {path}: {error}') from error


def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot read JSON {path}: {error}') from error


def read_jsonl(path):
    records = []
    try:
        with open(path) as handle:
            for number, line in enumerate(handle, 1):
                if line.strip():
                    try:
                        records.append(json.loads(line))
                    except ValueError as error:
                        raise ValueError(f'{path}:{number}: {error}') from error
    except OSError as error:
        raise OSError(f'Cannot proceed - required data file not found: {path}') from error
    return records


def clean(value):
    """Use JSON null for undefined statistics, including numpy scalar NaNs."""
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean(v) for v in value]
    if hasattr(value, 'item'):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_json(path, value):
    try:
        Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + '\n')
    except OSError as error:
        raise OSError(f'Cannot write result {path}: {error}') from error
