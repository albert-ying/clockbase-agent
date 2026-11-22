"""Read-only status and immutable in-progress snapshots for a corpus audit run.

The running controller is never imported. A study counts as complete only when
its terminal JSON and a valid matching SHA256 sidecar both exist. A terminal
without its sidecar, an empty/partial sidecar, or a progress event without a
terminal remains pending. Frozen output hashes are retained for the public
exporter's later byte-level readback; this status command checks their paths
and existence without rereading every live stage payload.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import sys


SHA256 = re.compile(r"[0-9a-f]{64}")
STUDY_ID = re.compile(r"GSE[0-9]+", re.ASCII)
MODALITIES = {"human_rna", "mouse_rna", "human_dnam"}
TERMINAL_STATUSES = {"REVIEWED", "PARTIALLY_REVIEWED", "UNRESOLVED",
                     "INTEGRITY_FAILURE", "BUDGET_DEFERRED", "CANNOT VERIFY"}
BUDGET_STATUSES = {"USAGE_ESTIMATE", "RESERVED_OR_UNCERTAIN"}


class StatusIntegrityError(RuntimeError):
    """A stable frozen or completed run artifact failed its integrity check."""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_json(path: Path):
    try:
        return json.loads(path.read_text())
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise StatusIntegrityError(f"Cannot read JSON {path}: {error}") from error


def _nonnegative_integer(value, name: str) -> int:
    if type(value) is not int or value < 0:
        raise StatusIntegrityError(f"{name} must be a nonnegative integer")
    return value


def _inside(root: Path, path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise StatusIntegrityError(f"Path escapes the run directory: {path}")
    return resolved


def _read_queue(run: Path, frozen: dict) -> tuple[list[dict], dict]:
    receipt = _read_json(run / "queue-receipt.json")
    queue_path = run / "queue.jsonl.gz"
    if set(receipt) != {"studies", "sha256"} or not SHA256.fullmatch(str(receipt["sha256"])):
        raise StatusIntegrityError("Malformed queue receipt")
    if _sha(queue_path) != receipt["sha256"]:
        raise StatusIntegrityError("Frozen queue hash differs from its receipt")
    try:
        with gzip.open(queue_path, "rt") as handle:
            queue = [json.loads(line) for line in handle]
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise StatusIntegrityError(f"Cannot read frozen queue: {error}") from error
    if _nonnegative_integer(receipt["studies"], "queue studies") != len(queue):
        raise StatusIntegrityError("Frozen queue count differs from its receipt")
    keys = []
    for task in queue:
        try:
            key, study, modality = task["study_key"], task["study_id"], task["modality"]
            artifact_ids = task["artifact_ids"]
        except (KeyError, TypeError) as error:
            raise StatusIntegrityError("Queue task lacks a required identity") from error
        if (modality not in MODALITIES or not STUDY_ID.fullmatch(str(study)) or
                key != f"{modality}:{study}" or not isinstance(artifact_ids, list) or
                any(not isinstance(item, str) or not item for item in artifact_ids) or
                len(artifact_ids) != len(set(artifact_ids))):
            raise StatusIntegrityError(f"Unsupported queue identity: {key!r}")
        keys.append(key)
    if len(keys) != len(set(keys)):
        raise StatusIntegrityError("Frozen queue has duplicate study keys")
    expected = frozen.get("denominators", {}).get("canonical_candidate_study_by_biological_modality")
    if expected is not None and expected != len(queue):
        raise StatusIntegrityError("Queue count differs from the frozen study-modality denominator")
    return queue, receipt


def _read_gzip_json_lines(path: Path):
    try:
        with gzip.open(path, "rt") as handle:
            for number, line in enumerate(handle, 1):
                yield number, json.loads(line)
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise StatusIntegrityError(f"Cannot read frozen worklist {path}: {error}") from error


def _load_worklist(frozen: dict) -> tuple[dict[str, dict], dict[str, str]]:
    """Load source row identities without importing the running producer."""
    worklist = Path(frozen["worklist"])
    summary = _read_json(worklist / "summary.json")
    outputs = summary.get("outputs", {})
    required = ("artifact-worklist.jsonl.gz", "clock-row-worklist.jsonl.gz")
    for name in required:
        expected = outputs.get(name)
        path = worklist / name
        if not path.is_file() or not SHA256.fullmatch(str(expected)) or _sha(path) != expected:
            raise StatusIntegrityError(f"Frozen worklist member changed: {name}")
    artifacts = {}
    for _, artifact in _read_gzip_json_lines(worklist / required[0]):
        identifier = artifact.get("artifact_id")
        if (not isinstance(identifier, str) or not identifier or identifier in artifacts or
                not isinstance(artifact.get("path"), str) or
                not SHA256.fullmatch(str(artifact.get("sha256"))) or
                type(artifact.get("parsed_rows")) is not int or artifact["parsed_rows"] < 0):
            raise StatusIntegrityError("Malformed or duplicate frozen artifact identity")
        artifacts[identifier] = artifact
    if len(artifacts) != frozen["denominators"]["artifacts"]:
        raise StatusIntegrityError("Frozen artifact count differs from its denominator")
    row_identity = {}
    for _, row in _read_gzip_json_lines(worklist / required[1]):
        artifact = artifacts.get(row.get("artifact_id"))
        if artifact is None or row.get("artifact_sha256") != artifact["sha256"]:
            raise StatusIntegrityError("Frozen row has no matching artifact binding")
        ordinal = row.get("logical_csv_record")
        if type(ordinal) is not int or ordinal < 2:
            raise StatusIntegrityError("Frozen row has an invalid logical record number")
        coordinates = [artifact["path"], artifact["sha256"], ordinal]
        full = "row:" + hashlib.sha256(json.dumps(
            coordinates, ensure_ascii=False, sort_keys=True, separators=(",", ":")
        ).encode()).hexdigest()
        short = "R" + hashlib.sha256(json.dumps(
            coordinates, sort_keys=True, allow_nan=False
        ).encode()).hexdigest()[:24]
        if row.get("row_id") != full or short in row_identity:
            raise StatusIntegrityError("Frozen row identity is malformed or duplicated")
        identity = row.get("identity_status")
        if identity not in {"PASS", "UNRESOLVED"}:
            raise StatusIntegrityError("Frozen row has an unsupported identity status")
        row_identity[short] = identity
    if len(row_identity) != frozen["denominators"]["parsed_logical_rows"]:
        raise StatusIntegrityError("Frozen row count differs from its denominator")
    return artifacts, row_identity


def _validate_freeze(run: Path) -> tuple[dict, list[dict], dict]:
    config_path, sidecar = run / "config.json", run / "config.sha256"
    try:
        expected = sidecar.read_text().strip()
    except OSError as error:
        raise StatusIntegrityError(f"Cannot read frozen configuration sidecar: {error}") from error
    if not SHA256.fullmatch(expected) or _sha(config_path) != expected:
        raise StatusIntegrityError("Frozen configuration changed or lacks its valid digest")
    frozen = _read_json(config_path)
    if frozen.get("schema") != "clockbase-corpus-audit-v1":
        raise StatusIntegrityError("Unsupported corpus configuration schema")
    config = frozen.get("config")
    if not isinstance(config, dict):
        raise StatusIntegrityError("Frozen configuration payload is missing")
    for key in ("max_total_cost_usd", "max_total_calls"):
        if key not in config:
            raise StatusIntegrityError(f"Frozen configuration lacks {key}")
    implementation = frozen.get("implementation_sha256")
    if not isinstance(implementation, dict) or not implementation:
        raise StatusIntegrityError("Frozen implementation ledger is missing")
    for relative, digest in implementation.items():
        if not isinstance(relative, str) or not SHA256.fullmatch(str(digest)):
            raise StatusIntegrityError("Malformed frozen implementation ledger")
        target = _inside(run / "implementation", run / "implementation" / relative)
        if not target.is_file() or _sha(target) != digest:
            raise StatusIntegrityError(f"Frozen implementation snapshot changed: {relative}")
    worklist = Path(frozen.get("worklist", ""))
    summary = worklist / "summary.json"
    expected_summary = frozen.get("worklist_summary_sha256")
    if not summary.is_file() or not SHA256.fullmatch(str(expected_summary)) or _sha(summary) != expected_summary:
        raise StatusIntegrityError("Frozen worklist summary changed or is unavailable")
    source = _read_json(summary)
    if source.get("denominators") != frozen.get("denominators"):
        raise StatusIntegrityError("Frozen worklist and run denominators differ")
    queue, receipt = _read_queue(run, frozen)
    return frozen, queue, receipt


def _read_progress(run: Path, queued: set[str]) -> tuple[dict[str, list[dict]], dict]:
    path = run / "progress.jsonl"
    if not path.exists():
        return {}, {"complete_events": 0, "ignored_partial_tail": False}
    try:
        raw = path.read_bytes()
    except OSError as error:
        raise StatusIntegrityError(f"Cannot read progress ledger: {error}") from error
    lines = raw.splitlines()
    ignored_tail = bool(raw and not raw.endswith(b"\n"))
    if ignored_tail:
        lines = lines[:-1]
    events: dict[str, list[dict]] = defaultdict(list)
    for number, raw_line in enumerate(lines, 1):
        try:
            event = json.loads(raw_line)
        except (UnicodeError, json.JSONDecodeError) as error:
            raise StatusIntegrityError(f"Malformed complete progress record {number}: {error}") from error
        if not isinstance(event, dict) or event.get("study_key") not in queued:
            raise StatusIntegrityError(f"Progress record {number} has an unknown study identity")
        if event.get("status") == "CONTROLLER_FAILURE":
            if event.get("terminal_sha256") is not None:
                raise StatusIntegrityError("Controller failure unexpectedly binds a terminal")
        else:
            digest = event.get("terminal_sha256")
            if digest is not None and not SHA256.fullmatch(str(digest)):
                raise StatusIntegrityError("Progress record has a malformed terminal hash")
        _nonnegative_integer(event.get("reviewed_rows", 0), "progress reviewed_rows")
        events[event["study_key"]].append(event)
    return dict(events), {"complete_events": len(lines), "ignored_partial_tail": ignored_tail}


def _partial_or_digest(sidecar: Path) -> str | None:
    try:
        text = sidecar.read_text().strip()
    except OSError:
        return None
    # write_text can expose an empty or prefix-only sidecar to a concurrent
    # observer. Only a complete digest is treated as stable evidence.
    return text if SHA256.fullmatch(text) else None


def _terminal(run: Path, task: dict) -> tuple[dict | None, str | None, str | None]:
    study = run / "studies" / task["study_key"].replace(":", "-")
    terminal_path, sidecar = study / "terminal.json", study / "terminal.sha256"
    if not terminal_path.exists() or not sidecar.exists():
        return None, None, "INCOMPLETE_PAIR" if terminal_path.exists() or sidecar.exists() else None
    expected = _partial_or_digest(sidecar)
    if expected is None:
        return None, None, "INCOMPLETE_SIDECAR"
    observed = _sha(terminal_path)
    if observed != expected:
        raise StatusIntegrityError(f"Terminal hash mismatch: {task['study_key']}")
    terminal = _read_json(terminal_path)
    if (terminal.get("study_key") != task["study_key"] or terminal.get("study_id") != task["study_id"] or
            terminal.get("modality") != task["modality"] or terminal.get("artifact_ids") != task["artifact_ids"]):
        raise StatusIntegrityError(f"Terminal identity differs from its queue task: {task['study_key']}")
    if terminal.get("status") not in TERMINAL_STATUSES:
        raise StatusIntegrityError(f"Unsupported terminal status: {task['study_key']}")
    if terminal.get("scientific_acceptance") != "CANNOT VERIFY" or terminal.get("accepted_scientific_repairs") != 0:
        raise StatusIntegrityError(f"Terminal inflates scientific acceptance: {task['study_key']}")
    reviewed = _nonnegative_integer(terminal.get("reviewed_rows"), "terminal reviewed_rows")
    calculations = _nonnegative_integer(terminal.get("calculations"), "terminal calculations")
    unresolved = _nonnegative_integer(terminal.get("unresolved_calculations"), "terminal unresolved_calculations")
    # PASS calculations and unresolved targets are disjoint counts.
    row_ids = terminal.get("row_ids", [])
    if not isinstance(row_ids, list) or len(row_ids) != reviewed or len(row_ids) != len(set(row_ids)):
        raise StatusIntegrityError(f"Terminal row identities disagree with reviewed_rows: {task['study_key']}")
    if any(not isinstance(row, str) or not re.fullmatch(r"R[0-9a-f]{24}", row) for row in row_ids):
        raise StatusIntegrityError(f"Malformed terminal row identity: {task['study_key']}")
    identity_counts = terminal.get("identity_status_counts", {})
    if not isinstance(identity_counts, dict) or any(type(value) is not int or value < 0 for value in identity_counts.values()):
        raise StatusIntegrityError(f"Malformed identity counts: {task['study_key']}")
    # Producer identity counts describe all input rows, including rows not
    # reached by a partial review. Reviewed-row identity is derived below from
    # the frozen worklist rather than inferred from this study-level summary.
    if "input_rows" in terminal:
        input_rows = _nonnegative_integer(terminal["input_rows"], "terminal input_rows")
        if sum(identity_counts.values()) != input_rows:
            raise StatusIntegrityError(f"Identity counts disagree with input_rows: {task['study_key']}")
    outputs = terminal.get("output_sha256")
    if not isinstance(outputs, dict) or not outputs:
        raise StatusIntegrityError(f"Terminal output ledger is missing: {task['study_key']}")
    for relative, digest in outputs.items():
        if not isinstance(relative, str) or not SHA256.fullmatch(str(digest)):
            raise StatusIntegrityError(f"Malformed terminal output binding: {task['study_key']}")
        target = _inside(study, study / relative)
        if not target.is_file():
            raise StatusIntegrityError(f"Terminal-bound output is missing: {task['study_key']}:{relative}")
    return terminal, observed, None


def _read_budget(run: Path, config: dict) -> dict:
    path = run / "budget.sqlite3"
    if not path.is_file():
        return {"max_cost_usd": config["max_total_cost_usd"], "max_calls": config["max_total_calls"],
                "call_reservations": 0, "estimated_or_reserved_usd": 0.0, "by_status": [],
                "billed_cost_usd": None, "read_scope": "READ_ONLY_SQLITE_TRANSACTION"}
    uri = "file:" + str(path.resolve()) + "?mode=ro"
    try:
        db = sqlite3.connect(uri, uri=True, timeout=60)
        db.execute("PRAGMA query_only=ON")
        db.execute("BEGIN")
        budget_config = db.execute("SELECT cost,calls FROM config WHERE id=1").fetchone()
        rows = db.execute(
            "SELECT id,reserved,charged,status,updated,binding_sha256 FROM reservations ORDER BY id"
        ).fetchall()
        db.rollback()
    except sqlite3.Error as error:
        raise StatusIntegrityError(f"Cannot read aggregate budget ledger: {error}") from error
    finally:
        if "db" in locals():
            db.close()
    expected = (config["max_total_cost_usd"], config["max_total_calls"])
    if budget_config != expected:
        raise StatusIntegrityError("Aggregate budget configuration differs from the frozen run")
    by_status: dict[str, list[float]] = defaultdict(lambda: [0, 0.0])
    for identifier, reserved, charged, status, updated, binding in rows:
        if (not SHA256.fullmatch(str(identifier)) or status not in BUDGET_STATUSES or
                not SHA256.fullmatch(str(binding)) or not isinstance(updated, str) or not updated or
                type(reserved) not in (int, float) or type(charged) not in (int, float) or
                not math.isfinite(reserved) or not math.isfinite(charged) or reserved <= 0 or charged < 0 or
                charged > reserved + 1e-10):
            raise StatusIntegrityError("Malformed aggregate budget reservation")
        if status == "RESERVED_OR_UNCERTAIN" and abs(charged - reserved) > 1e-10:
            raise StatusIntegrityError("Uncertain reservation does not retain its conservative charge")
        by_status[status][0] += 1
        by_status[status][1] += charged
    total = sum(value[1] for value in by_status.values())
    if len(rows) > config["max_total_calls"] or total > config["max_total_cost_usd"] + 1e-10:
        raise StatusIntegrityError("Aggregate budget ledger exceeds the frozen cap")
    return {"max_cost_usd": config["max_total_cost_usd"], "max_calls": config["max_total_calls"],
            "call_reservations": len(rows), "estimated_or_reserved_usd": total,
            "by_status": [{"status": status, "calls": values[0], "usd": values[1]}
                          for status, values in sorted(by_status.items())],
            "billed_cost_usd": None, "read_scope": "READ_ONLY_SQLITE_TRANSACTION"}


def inspect_run(run: Path) -> dict:
    """Return an exporter-compatible in-progress report without changing *run*."""
    run = Path(run).resolve()
    frozen, queue, queue_receipt = _validate_freeze(run)
    queued = {task["study_key"]: task for task in queue}
    artifacts, row_identity = _load_worklist(frozen)
    if any(artifact not in artifacts for task in queue for artifact in task["artifact_ids"]):
        raise StatusIntegrityError("Queue references an unknown frozen artifact")
    progress, progress_meta = _read_progress(run, set(queued))
    terminals: dict[str, dict] = {}
    terminal_hashes: dict[str, str] = {}
    partial: dict[str, str] = {}
    for task in queue:
        terminal, digest, partial_reason = _terminal(run, task)
        if terminal is not None:
            terminals[task["study_key"]] = terminal
            terminal_hashes[task["study_key"]] = digest
        elif partial_reason:
            partial[task["study_key"]] = partial_reason

    failures = {}
    duplicate_events = 0
    recovered_failures = 0
    progress_without_terminal = 0
    for key, events in progress.items():
        duplicate_events += max(0, len(events) - 1)
        hashes = {event["terminal_sha256"] for event in events if event.get("terminal_sha256") is not None}
        if len(hashes) > 1:
            raise StatusIntegrityError(f"Progress binds conflicting terminals: {key}")
        failures_for_key = [event for event in events if event.get("status") == "CONTROLLER_FAILURE"]
        if key in terminals:
            if hashes and hashes != {terminal_hashes[key]}:
                raise StatusIntegrityError(f"Progress and terminal hashes differ: {key}")
            if failures_for_key:
                recovered_failures += 1
        elif failures_for_key:
            failures[key] = {"study_key": key, "status": "CONTROLLER_FAILURE",
                             "error": "Controller failure recorded in progress.jsonl; detail is unavailable in the progress ledger."}
        elif events:
            progress_without_terminal += 1

    results = [terminals[key] for key in sorted(terminals)] + [failures[key] for key in sorted(failures)]
    result_keys = {result["study_key"] for result in results}
    status_counts = dict(sorted(Counter(result["status"] for result in results).items()))
    reviewed_assignments = sum(result.get("reviewed_rows", 0) for result in results)
    reviewed_ids = {row for result in results for row in result.get("row_ids", [])}
    calculations = sum(result.get("calculations", 0) for result in results)
    unresolved = sum(result.get("unresolved_calculations", 0) for result in results)
    if len(reviewed_ids) > reviewed_assignments:
        raise StatusIntegrityError("Aggregate terminal counts are internally inconsistent")
    denominators = frozen["denominators"]
    if reviewed_ids - set(row_identity):
        raise StatusIntegrityError("Terminal references a row outside the frozen worklist")
    assigned_artifacts = {artifact for task in queue for artifact in task["artifact_ids"]}
    attempted_artifacts = {artifact for result in results for artifact in result.get("artifact_ids", [])}
    bound_reviewed = sum(row_identity[row] == "PASS" for row in reviewed_ids)
    unassigned = [artifact for identifier, artifact in artifacts.items() if identifier not in assigned_artifacts]
    pending = len(queue) - len(results)
    if pending < 0:
        raise StatusIntegrityError("Result membership exceeds the frozen queue")
    budget = _read_budget(run, frozen["config"])
    coverage = {
        "saved_artifacts": denominators["artifacts"],
        "parsed_logical_rows": denominators["parsed_logical_rows"],
        "zero_row_artifacts": denominators["zero_row_artifacts"],
        "artifacts_assigned_to_candidate_studies": len(assigned_artifacts),
        "artifacts_without_candidate_study": len(unassigned),
        "rows_without_candidate_study": sum(item["parsed_rows"] for item in unassigned),
        "unassigned_artifact_ids": sorted(item["artifact_id"] for item in unassigned),
        "attempted_artifacts": len(attempted_artifacts),
        "reviewed_unique_rows": len(reviewed_ids),
        "reviewed_rows_with_bound_study_identity": bound_reviewed,
        "reviewed_rows_with_unresolved_study_identity": len(reviewed_ids) - bound_reviewed,
        "not_reviewed_rows": denominators["parsed_logical_rows"] - len(reviewed_ids),
        "pending_study_tasks": pending,
        "scope": "Live terminal coverage only. Running or partially written studies remain pending; source-less row counts require the completed worklist coverage pass.",
    }
    progress_summary = {
        "queued_study_tasks": len(queue),
        "completed_terminal_studies": len(terminals),
        "controller_failure_studies": len(failures),
        "selected_completed_or_failed_studies": len(results),
        "pending_or_running_studies": pending,
        "partial_terminal_pairs": len(partial),
        "partial_terminal_reasons": dict(sorted(Counter(partial.values()).items())),
        "progress_events": progress_meta["complete_events"],
        "ignored_partial_progress_tail": progress_meta["ignored_partial_tail"],
        "duplicate_or_resume_progress_events": duplicate_events,
        "recovered_controller_failure_studies": recovered_failures,
        "progress_without_complete_terminal": progress_without_terminal,
        "result_study_keys_sha256": hashlib.sha256(
            json.dumps(sorted(result_keys), separators=(",", ":")).encode()
        ).hexdigest(),
    }
    return {
        "schema": "clockbase-corpus-status-snapshot-v1",
        "status": "IN_PROGRESS",
        "created_at_utc": _now(),
        "selected_studies": len(results),
        "whole_corpus_studies": len(queue),
        "status_counts": status_counts,
        "reviewed_row_assignments": reviewed_assignments,
        "unique_reviewed_rows": len(reviewed_ids),
        "calculations": calculations,
        "unresolved_calculations": unresolved,
        "accepted_scientific_repairs": 0,
        "scientific_error_rate": None,
        "scientific_acceptance": "CANNOT VERIFY",
        "inventory_denominators": denominators,
        "coverage": coverage,
        "progress": progress_summary,
        "budget": budget,
        "source_bindings": {
            "config_sha256": _sha(run / "config.json"),
            "queue_sha256": queue_receipt["sha256"],
            "worklist_summary_sha256": frozen["worklist_summary_sha256"],
            "verified_terminal_sidecars": len(terminals),
            "terminal_output_files_present": sum(len(item["output_sha256"]) for item in terminals.values()),
            "terminal_output_file_hashes_rechecked": 0,
            "output_hash_scope": "Terminal ledgers are sealed and paths exist; export_corpus rechecks every payload hash before release.",
        },
        "results": results,
        "scope": "Read-only operational progress snapshot. Completion and calculations are terminal-record counts, not scientific acceptance or an error-rate estimate.",
    }


def write_snapshot(run: Path, destination: Path, report: dict) -> dict:
    run, destination = Path(run).resolve(), Path(destination)
    destination = (run / destination).resolve() if not destination.is_absolute() else destination.resolve()
    if not destination.is_relative_to(run):
        raise ValueError("Snapshot path must be inside the run directory")
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = (json.dumps(report, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode()
    try:
        descriptor = os.open(destination, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except FileExistsError as error:
        raise FileExistsError(f"Snapshot already exists: {destination}") from error
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        # Preserve a failed/partial attempt for diagnosis; never replace it.
        raise
    return {"path": str(destination), "sha256": hashlib.sha256(payload).hexdigest(), "bytes": len(payload)}


def _public_summary(report: dict) -> dict:
    return {key: value for key, value in report.items() if key != "results"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--snapshot", type=Path)
    args = parser.parse_args(argv)
    report = inspect_run(args.run)
    result = _public_summary(report)
    if args.snapshot is not None:
        result["snapshot"] = write_snapshot(args.run, args.snapshot, report)
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Corpus status stopped: {type(error).__name__}: {error}", file=sys.stderr)
        raise SystemExit(1)
