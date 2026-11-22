"""Canonical content identity for audit records and files."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path


class ArtifactIdentityError(RuntimeError):
    """Raised when artifact bytes cannot be inspected."""


@dataclass(frozen=True)
class FileIdentity:
    path: Path
    sha256: str
    byte_count: int


def canonical_json_bytes(value: object) -> bytes:
    """Serialize a JSON-compatible value to one canonical byte sequence."""

    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def sha256_value(value: object) -> str:
    """Return the SHA-256 digest of a canonical JSON value."""

    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def inspect_file(path: Path) -> FileIdentity:
    """Reopen a file and return its resolved path, digest, and byte count."""

    source = Path(path)
    digest = hashlib.sha256()
    byte_count = 0
    try:
        with source.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
                byte_count += len(block)
    except OSError as exc:
        raise ArtifactIdentityError(
            f"Cannot inspect artifact identity for {source}: {exc}"
        ) from exc
    return FileIdentity(source.resolve(), digest.hexdigest(), byte_count)
