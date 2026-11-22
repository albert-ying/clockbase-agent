"""Atomic, logically append-only storage for audit records."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from .hashing import canonical_json_bytes
from .models import (
    AuditRecord,
    Record,
    RecordValidationError,
    record_from_dict,
    record_to_dict,
)


class LedgerError(RuntimeError):
    """Base class for ledger storage failures."""


class DuplicateRecordError(LedgerError):
    """Raised when an immutable identifier is reused with new content."""


class CorruptLedgerError(LedgerError):
    """Raised when stored ledger bytes do not parse as valid records."""


@dataclass(frozen=True)
class LedgerSnapshot:
    records: Mapping[str, Record]
    payloads: Mapping[str, dict[str, object]]
    record_ids: tuple[str, ...]
    sha256: str


class Ledger:
    """Store immutable JSON Lines records with atomic whole-file replacement."""

    def __init__(self, path: Path):
        self.path = Path(path)
        self.lock_path = self.path.with_name(f"{self.path.name}.lock")

    def snapshot(self) -> LedgerSnapshot:
        try:
            raw = self.path.read_bytes() if self.path.exists() else b""
        except OSError as exc:
            raise LedgerError(f"Cannot read ledger {self.path}: {exc}") from exc
        records: dict[str, Record] = {}
        payloads: dict[str, dict[str, object]] = {}
        for line_number, line in enumerate(raw.splitlines(), start=1):
            if not line.strip():
                raise CorruptLedgerError(
                    f"Ledger {self.path} has a blank line at line {line_number}"
                )
            try:
                payload = json.loads(line)
                record = record_from_dict(payload)
                canonical = record_to_dict(record)
            except (json.JSONDecodeError, RecordValidationError) as exc:
                raise CorruptLedgerError(
                    f"Ledger {self.path} has an invalid record at line "
                    f"{line_number}: {exc}"
                ) from exc
            prior = payloads.get(record.record_id)
            if prior is not None and prior != canonical:
                raise CorruptLedgerError(
                    f"Ledger {self.path} reuses {record.record_id} with "
                    "different content"
                )
            records[record.record_id] = record
            payloads[record.record_id] = canonical
        record_ids = tuple(records)
        return LedgerSnapshot(
            records=MappingProxyType(records),
            payloads=MappingProxyType(payloads),
            record_ids=record_ids,
            sha256=hashlib.sha256(raw).hexdigest(),
        )

    def append(self, record: AuditRecord) -> str:
        try:
            payload = record_to_dict(record)
            encoded = canonical_json_bytes(payload) + b"\n"
        except (RecordValidationError, TypeError, ValueError):
            raise

        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.lock_path.touch(exist_ok=True)
        temporary_path: Path | None = None
        try:
            with self.lock_path.open("r+b") as lock_handle:
                fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
                snapshot = self.snapshot()
                prior = snapshot.payloads.get(record.record_id)
                if prior is not None:
                    if prior == payload:
                        return snapshot.sha256
                    raise DuplicateRecordError(
                        f"Record identifier {record.record_id} already exists "
                        "with different content"
                    )
                existing = self.path.read_bytes() if self.path.exists() else b""
                with tempfile.NamedTemporaryFile(
                    mode="wb",
                    dir=self.path.parent,
                    prefix=f".{self.path.name}.",
                    delete=False,
                ) as temporary:
                    temporary_path = Path(temporary.name)
                    temporary.write(existing)
                    temporary.write(encoded)
                    temporary.flush()
                    os.fsync(temporary.fileno())
                os.replace(temporary_path, self.path)
                temporary_path = None
                return hashlib.sha256(existing + encoded).hexdigest()
        except OSError as exc:
            raise LedgerError(f"Cannot append to ledger {self.path}: {exc}") from exc
        finally:
            if temporary_path is not None:
                try:
                    temporary_path.unlink(missing_ok=True)
                except OSError:
                    pass
