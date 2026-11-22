"""Build a local W0 evidence ledger from selected existing reports."""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path

from .adapters import adapt_measurement_pilot, adapt_protected_repair
from .hashing import canonical_json_bytes, inspect_file
from .ledger import Ledger


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a read-only ClockBase W0 evidence ledger.",
    )
    parser.add_argument("--measurement-pilot", type=Path, required=True)
    parser.add_argument("--protected-manifest", type=Path, required=True)
    parser.add_argument("--protected-summary", type=Path, required=True)
    parser.add_argument("--release-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def build(args: argparse.Namespace) -> dict[str, object]:
    output = args.output
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {output}")

    source_paths = (
        args.measurement_pilot,
        args.protected_manifest,
        args.protected_summary,
    )
    source_identities = [inspect_file(path) for path in source_paths]
    records = (
        *adapt_measurement_pilot(args.measurement_pilot, args.release_id),
        *adapt_protected_repair(
            args.protected_manifest,
            args.protected_summary,
            args.release_id,
        ),
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=output.parent,
        prefix=f".{output.name}.",
        delete=True,
    ) as handle:
        temporary_path = Path(handle.name)
    temporary_lock = temporary_path.with_name(f"{temporary_path.name}.lock")
    try:
        ledger = Ledger(temporary_path)
        for record in records:
            ledger.append(record)
        snapshot = ledger.snapshot()
        if len(snapshot.record_ids) != len(records):
            raise RuntimeError(
                "Ledger record count differs from the validated adapter records"
            )
        os.replace(temporary_path, output)
    finally:
        temporary_path.unlink(missing_ok=True)
        temporary_lock.unlink(missing_ok=True)

    return {
        "release_id": args.release_id,
        "ledger_sha256": snapshot.sha256,
        "record_count": len(snapshot.record_ids),
        "source_report_sha256": {
            str(identity.path): identity.sha256 for identity in source_identities
        },
    }


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        summary = build(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(canonical_json_bytes(summary).decode("utf-8"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
