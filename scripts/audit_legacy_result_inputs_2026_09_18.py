#!/usr/bin/env python3
"""Read-only census of the legacy result CSV boundary; never execute plotting."""
import argparse
import ast
from collections import Counter
import csv
import gzip
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import struct
import sys
import time

import pandas as pd


ROOTS = ("result_mm_rna", "result_hs_rna", "hs_dnam_sonnet37", "result_hs_dnam")
SAVED_ROOTS = ("result_mm_rna", "hs_rna_sonnet37", "result_hs_rna", "hs_dnam_sonnet37", "result_hs_dnam")
REQUIRED = ("gse_id", "comparison_name", "condition_name", "condition_category",
            "clock", "effect_size", "p_value", "trust_score", "direction")
LIMITS = {"files": 30000, "bytes": 2000000000, "wall_seconds": 600}
SAVED_CORPUS_LIMITS = {"files": 100000, "bytes": 10000000000, "wall_seconds": 1800}
SAVED_SCORE_TABLES = (
    "data/detailed_scores_result_mm_rna.csv",
    "data/scores_summary_hs_rna.csv",
    "data/scores_summary_hs_dnam.csv",
)
SAVED_TRACE_ARCHIVES = {
    "successful_traces_mouse.pkl": "5b1271b726458e1e4cbdef00af34d0abbf93c33f5932d9bf50a6456bd66ba24e",
    "successful_traces_human_rna.pkl": "b1680f27ebfe9d912a11b02ca65cd854a52b97229c87131db76be9cbece4da14",
    "successful_traces_human_dnam.pkl": "7207f7368782d23416f63a15c9937fdbfc1d824afd14c73b01e86097b8176548",
}
SAVED_MODALITY = {
    "result_mm_rna": "mouse_rna",
    "hs_rna_sonnet37": "human_rna_sonnet37",
    "result_hs_rna": "human_rna",
    "hs_dnam_sonnet37": "human_dnam_sonnet37",
    "result_hs_dnam": "human_dnam",
}
WORKLIST_MODALITY = {name: modality.removesuffix("_sonnet37") for name, modality in SAVED_MODALITY.items()}
WORKLIST_GEO_FOLDERS = {"mouse_rna": "GEO_mm_rna", "human_rna": "GEO_hs_rna", "human_dnam": "GEO_hs_dnam"}
WORKLIST_SCHEMA = "clockbase-comparison-audit-worklist-v1"
WORKLIST_STAGES = ("NOT_YET_AUDITED", "REVIEWED", "RERUN", "VERIFIED", "UNRESOLVED")
WORKLIST_DEPENDENCIES = (
    "SOURCE_DESIGN_REVIEW_REQUIRED", "SOURCE_SUPPORTED_SAMPLE_GROUPS_REQUIRED",
    "STATISTICAL_METHOD_PLAN_REQUIRED", "CLOCK_APPLICABILITY_REVIEW_REQUIRED",
    "SCIENTIFIC_CLAIM_REVIEW_REQUIRED", "HISTORICAL_PREDICTION_BINDING_UNVERIFIED",
    "INDEPENDENT_VERIFICATION_REQUIRED",
)


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''): digest.update(block)
    return digest.hexdigest()


def numeric_stream_digest(frame):
    """Identity plus binary64 values; encoding is shared, parsing is separate."""
    import math
    digest = hashlib.sha256()
    counts = Counter()
    for feature, values in zip(frame.index, frame.to_numpy()):
        digest.update((json.dumps(str(feature), ensure_ascii=False, separators=(',', ':')) + '\n').encode())
        encoded = bytearray()
        for value in values:
            value = float(value)
            if math.isnan(value): encoded.extend(b'M'); counts['missing'] += 1
            elif math.isinf(value): encoded.extend(b'P' if value > 0 else b'N'); counts['nonfinite'] += 1
            else:
                encoded.extend(b'F' + struct.pack('<d', value)); counts['finite'] += 1
                counts['outside_beta_range'] += int(not 0 <= value <= 1)
        digest.update(encoded)
    return digest.hexdigest(), dict(counts)


def strict_csv_records(path):
    with path.open(newline='') as handle:
        reader = csv.reader(handle, strict=True)
        header = next(reader)
        rows = list(reader)
    if any(len(row) != len(header) for row in rows):
        raise ValueError('Ragged source result CSV: ' + str(path))
    return header, rows


def _saved_issue(code, status, row=None, detail=""):
    record = {"code": code, "status": status, "detail": detail}
    if row is not None:
        record["row"] = row
    return record


def _is_lfs_pointer(path):
    with path.open("rb") as handle:
        prefix = handle.read(256)
    return prefix.startswith(b"version https://git-lfs.github.com/spec/v1")


def _finite_number(value):
    import math
    if value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return number if math.isfinite(number) else False


def _saved_csv_audit(path, repo, parent_root, limits):
    """Stream one saved result CSV and retain all structural findings."""
    started = time.monotonic()
    stat = path.stat()
    record = {
        "path": str(path.relative_to(repo)), "root": parent_root,
        "modality": SAVED_MODALITY[parent_root], "parent_study_id": path.parent.name,
        "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns, "sha256": sha256(path),
        "status": "RUNNING", "rows": 0, "header": [], "missing_required_columns": [],
        "literal_accessions": [], "issues": [], "full_row_duplicate_groups": [],
        "conflicting_literal_key_groups": [], "elapsed_seconds": None,
    }
    if _is_lfs_pointer(path):
        record.update(status="CANNOT VERIFY", reason="LFS_POINTER",
                      issues=[_saved_issue("LFS_POINTER", "CANNOT VERIFY", detail="CSV payload is a Git LFS pointer")],
                      elapsed_seconds=time.monotonic() - started)
        return record
    seen_rows = {}
    key_rows = {}
    try:
        with path.open(newline="") as handle:
            reader = csv.reader(handle, strict=True)
            try:
                header = next(reader)
            except StopIteration:
                header = []
            record["header"] = header
            if not header:
                record["issues"].append(_saved_issue("EMPTY_HEADER", "FAIL"))
            if len(header) != len(set(header)):
                record["issues"].append(_saved_issue("DUPLICATE_HEADER", "FAIL"))
            record["missing_required_columns"] = [column for column in REQUIRED if column not in header]
            for column in record["missing_required_columns"]:
                record["issues"].append(_saved_issue("MISSING_REQUIRED_COLUMN", "FAIL", detail=column))
            positions = {}
            for index, column in enumerate(header):
                positions.setdefault(column, index)
            for row_number, row in enumerate(reader, 2):
                record["rows"] += 1
                if time.monotonic() - started > limits["wall_seconds"]:
                    raise TimeoutError("Saved-corpus wall-time limit exceeded")
                if len(row) != len(header):
                    record["issues"].append(_saved_issue("RAGGED_ROW", "FAIL", row_number,
                                                         f"fields={len(row)} expected={len(header)}"))
                    continue
                row_tuple = tuple(row)
                row_hash = hashlib.sha256(json.dumps(row, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()
                seen_rows.setdefault(row_hash, []).append(row_number)
                key = tuple(row[positions[column]] for column in ("gse_id", "comparison_name", "condition_name", "clock")
                            if column in positions)
                if len(key) == 4:
                    key_rows.setdefault(key, []).append((row_number, row_hash))
                for identity_column in ("gse_id", "comparison_name", "condition_name", "clock"):
                    if identity_column in positions and row[positions[identity_column]] == "":
                        record["issues"].append(_saved_issue("EMPTY_IDENTITY", "FAIL", row_number, identity_column))
                if "gse_id" in positions and row[positions["gse_id"]] != path.parent.name:
                    record["issues"].append(_saved_issue("PARENT_ID_MISMATCH", "FAIL", row_number,
                                                         f"literal={row[positions['gse_id']]} parent={path.parent.name}"))
                for numeric_column in ("effect_size", "p_value", "trust_score", "direction"):
                    if numeric_column not in positions:
                        continue
                    value = row[positions[numeric_column]]
                    numeric = _finite_number(value)
                    if numeric is None:
                        record["issues"].append(_saved_issue("MISSING_STATISTIC", "CANNOT VERIFY", row_number, numeric_column))
                    elif numeric is False:
                        record["issues"].append(_saved_issue("NONFINITE_OR_NONNUMERIC", "FAIL", row_number, numeric_column))
                    elif numeric_column == "p_value" and not 0 <= numeric <= 1:
                        record["issues"].append(_saved_issue("P_VALUE_OUT_OF_BOUNDS", "FAIL", row_number, value))
    except csv.Error as exc:
        record["issues"].append(_saved_issue("CSV_PARSE_ERROR", "FAIL", detail=str(exc)))
    except UnicodeError as exc:
        record["issues"].append(_saved_issue("TEXT_DECODE_ERROR", "FAIL", detail=str(exc)))
    record["literal_accessions"] = sorted({key[0] for key in key_rows if key})
    for row_hash, rows in sorted(seen_rows.items()):
        if len(rows) > 1:
            record["full_row_duplicate_groups"].append({"row_numbers": rows, "row_sha256": row_hash})
    for key, entries in sorted(key_rows.items()):
        distinct = sorted({row_hash for _, row_hash in entries})
        if len(distinct) > 1:
            record["conflicting_literal_key_groups"].append({
                "key": list(key), "row_numbers": [row_number for row_number, _ in entries],
                "distinct_row_sha256": distinct,
            })
    statuses = {issue["status"] for issue in record["issues"]}
    record["status"] = "FAIL" if "FAIL" in statuses else ("CANNOT VERIFY" if "CANNOT VERIFY" in statuses else "PASS")
    record["elapsed_seconds"] = time.monotonic() - started
    return record


def _saved_register(path, repo, kind, expected_sha256=None):
    record = {"path": str(path.relative_to(repo)), "kind": kind, "status": "CANNOT VERIFY",
              "available": path.exists(), "sha256": None, "bytes": None, "mtime_ns": None,
              "expected_sha256": expected_sha256, "parsed": False}
    if not path.exists():
        record["reason"] = "MISSING_PAYLOAD"
        return record
    stat = path.stat(); record.update(bytes=stat.st_size, mtime_ns=stat.st_mtime_ns, sha256=sha256(path))
    if _is_lfs_pointer(path):
        record["reason"] = "LFS_POINTER"
    elif expected_sha256 and record["sha256"] != expected_sha256:
        record.update(status="FAIL", reason="LFS_OBJECT_HASH_MISMATCH")
    else:
        record.update(status="PASS", reason="REGISTERED_ONLY")
    return record


def execute_saved_corpus(repo, plan_path, output, limits=None, roots=None, score_tables=None, trace_archives=None):
    """Bounded saved-output census; score tables and traces are registered, never parsed."""
    limits = limits or SAVED_CORPUS_LIMITS
    roots = tuple(roots or SAVED_ROOTS)
    score_tables = tuple(score_tables or SAVED_SCORE_TABLES)
    trace_archives = trace_archives or SAVED_TRACE_ARCHIVES
    started = time.monotonic()
    output.mkdir(parents=True, exist_ok=False)
    source_records = []
    root_records = []
    for root_name in roots:
        root = repo / root_name
        if not root.is_dir():
            root_records.append({"path": root_name, "status": "CANNOT VERIFY", "files": 0, "reason": "MISSING_ROOT"})
            continue
        found = sorted(path for path in root.rglob("aging_analysis*.csv") if path.is_file())
        root_records.append({"path": root_name, "status": "PRESENT", "files": len(found)})
        source_records.extend((path, root_name) for path in found)
    total_bytes = sum(path.stat().st_size for path, _ in source_records)
    if len(source_records) > limits["files"] or total_bytes > limits["bytes"]:
        raise ValueError("Saved-corpus inventory exceeds frozen file/byte limits")
    registrations = []
    for relative in score_tables:
        registrations.append(_saved_register(repo / relative, repo, "detailed_score_table"))
    for name, expected in sorted(trace_archives.items()):
        registrations.append(_saved_register(repo / name, repo, "trace_archive", expected))
    try:
        commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "CANNOT VERIFY: repository is not a Git checkout"
    manifest = {
        "schema": "clockbase-saved-analysis-census-v1", "status": "RUNNING",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(), "repository": str(repo),
        "commit": commit,
        "plan_path": str(plan_path), "plan_sha256": sha256(plan_path),
        "implementation_sha256": sha256(Path(__file__)), "runner_command": sys.argv,
        "python": sys.version, "platform": platform.platform(), "limits": limits,
        "root_records": root_records, "inventory_files": len(source_records), "inventory_bytes": total_bytes,
        "registration_count": len(registrations), "scientific_status": "CANNOT VERIFY",
        "network_requests": 0, "paid_model_api_calls": 0, "clock_executions": 0,
        "score_tables_parsed": False, "trace_archives_parsed": False,
        "scope": "All recursively matched aging_analysis*.csv artifacts; score tables and trace archives registered only.",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    counts = Counter(); rows_total = 0; accessions = set(); accession_modalities = set(); records = []
    try:
        with (output / "artifact-ledger.jsonl").open("x") as ledger:
            for path, root_name in source_records:
                record = _saved_csv_audit(path, repo, root_name, limits)
                records.append(record); rows_total += record["rows"]
                accessions.update(record["literal_accessions"])
                accession_modalities.update((acc, record["modality"]) for acc in record["literal_accessions"])
                counts[(root_name, record["status"])] += 1
                ledger.write(json.dumps(record, ensure_ascii=False) + "\n")
        with (output / "registration-ledger.jsonl").open("x") as ledger:
            for record in registrations:
                ledger.write(json.dumps(record, ensure_ascii=False) + "\n")
        with (output / "validation-table.csv").open("x", newline="") as handle:
            writer = csv.writer(handle); writer.writerow(["root", "status", "artifacts"])
            writer.writerows([root_name, status, count] for (root_name, status), count in sorted(counts.items()))
        for record in records:
            current = repo / record["path"]
            if not current.exists() or sha256(current) != record["sha256"]:
                raise ValueError("Source changed during saved-corpus census: " + record["path"])
        current_paths = [(str(path.relative_to(repo)), root_name) for root_name in roots if (repo / root_name).is_dir()
                         for path in sorted((repo / root_name).rglob("aging_analysis*.csv")) if path.is_file()]
        if current_paths != [(record["path"], record["root"]) for record in records]:
            raise ValueError("Saved-corpus inventory changed during census")
        if sha256(Path(__file__)) != manifest["implementation_sha256"]:
            raise ValueError("Auditor implementation changed during census")
        manifest.update(status="PASS", elapsed_seconds=time.monotonic() - started,
                        input_hashes_rechecked=len(records), input_inventory_rechecked=True,
                        counts=[{"root": root, "status": status, "artifacts": count}
                                for (root, status), count in sorted(counts.items())],
                        denominators={"artifacts": len(records), "rows": rows_total,
                                      "literal_accessions": len(accessions),
                                      "accession_by_modality": len(accession_modalities)},
                        registrations=registrations,
                        outputs={name: sha256(output / name) for name in
                                 ("artifact-ledger.jsonl", "registration-ledger.jsonl", "validation-table.csv")})
    except Exception as exc:
        manifest.update(status="FAIL", error=f"{type(exc).__name__}: {exc}",
                        elapsed_seconds=time.monotonic() - started, completed_records=len(records))
        raise
    finally:
        (output / "summary.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": manifest["status"], "artifacts": len(records),
                      "rows": rows_total, "output": str(output)}, indent=2))


def _worklist_id(kind, values):
    encoded = json.dumps(values, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    return kind + ":" + hashlib.sha256(encoded).hexdigest()


def _worklist_stages():
    # This queue does not import prior study-level reviews as row-level verification.
    return {stage: "NOT_YET_AUDITED" for stage in
            ("design_review", "statistics_review", "scientific_review", "repair", "rerun", "independent_verification")}


def _write_gz_jsonl(path, records):
    with path.open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="\n") as handle:
                for record in records:
                    handle.write(json.dumps(record, ensure_ascii=False, separators=(",", ":")) + "\n")


def _read_gz_jsonl(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            yield json.loads(line)


def _worklist_file(path, role):
    record = {"path": str(path), "role": role, "available": path.is_file(),
              "sha256": None, "bytes": None, "mtime_ns": None,
              "payload_status": "CANNOT VERIFY", "scientific_status": "CANNOT VERIFY"}
    if not record["available"]:
        record["reason"] = "MISSING_FILE"
        return record
    before = path.stat()
    record.update(sha256=sha256(path), bytes=before.st_size, mtime_ns=before.st_mtime_ns)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("Source changed while hashing: " + str(path))
    if _is_lfs_pointer(path):
        record["reason"] = "LFS_POINTER"
    elif before.st_size == 0:
        record["reason"] = "EMPTY_FILE"
    else:
        record.update(payload_status="PASS", reason="LOCATED_AND_HASHED_ONLY")
    return record


def _recheck_worklist_file(record):
    path = Path(record["path"])
    if path.is_file() != record["available"]:
        raise ValueError("Input availability changed: " + str(path))
    if record["available"] and sha256(path) != record["sha256"]:
        raise ValueError("Input payload changed: " + str(path))


def _worklist_family_map(manifest_path):
    if manifest_path is None:
        return {}, {"status": "UNRESOLVED", "reason": "NO_FROZEN_FAMILY_MANIFEST",
                    "leakage_safe_split_established": False}, []
    manifest = json.loads(manifest_path.read_text())
    enrollment = manifest_path.parent / "study-enrollment.csv"
    expected = manifest.get("output_sha256", {}).get("study-enrollment.csv")
    if not expected or sha256(enrollment) != expected:
        raise ValueError("Family enrollment is not bound to the frozen family manifest")
    families = {}
    with enrollment.open(newline="") as handle:
        reader = csv.DictReader(handle, strict=True)
        if not {"study_id", "family_id", "source_linkage_status"} <= set(reader.fieldnames or []):
            raise ValueError("Frozen family enrollment lacks source identity fields")
        for row in reader:
            study = row["study_id"]
            if study in families:
                raise ValueError("Duplicate study in frozen family enrollment: " + study)
            # Deliberately do not read split, outcome, exposure or reference labels.
            families[study] = {"family_id": row["family_id"],
                               "source_linkage_status": row["source_linkage_status"]}
    binding = {"status": "REUSED_SOURCE_FAMILY_GROUPS", "manifest_path": str(manifest_path),
               "manifest_sha256": sha256(manifest_path), "enrollment_sha256": expected,
               "source_coverage_complete": manifest.get("source_coverage_complete", False),
               "scope": manifest.get("scope", "Supplied source universe only"),
               "leakage_safe_split_established": False,
               "split_and_outcome_labels_used": False}
    return families, binding, [_worklist_file(manifest_path, "family_manifest"),
                               _worklist_file(enrollment, "source_family_enrollment")]


def _worklist_resume(output, census_ledger, geo_root, family_manifest):
    manifest = json.loads((output / "summary.json").read_text())
    expected = {"census_ledger": str(census_ledger), "geo_root": str(geo_root),
                "family_manifest": str(family_manifest) if family_manifest else None}
    if manifest.get("schema") != WORKLIST_SCHEMA or manifest.get("status") != "PASS":
        raise ValueError("Only a completed worklist can resume; preserve failed output and use a new directory")
    if manifest.get("arguments") != expected or manifest.get("implementation_sha256") != sha256(Path(__file__)):
        raise ValueError("Worklist resume arguments or implementation changed")
    for name, expected_hash in manifest["outputs"].items():
        if sha256(output / name) != expected_hash:
            raise ValueError("Worklist output changed: " + name)
    for record in _read_gz_jsonl(output / "input-files.jsonl.gz"):
        _recheck_worklist_file(record)
    print(json.dumps({"status": "PASS", "resumed": True, "output": str(output),
                      "denominators": manifest["denominators"], "model_calls": 0,
                      "scientific_acceptance": "CANNOT VERIFY"}, indent=2), flush=True)
    return manifest


def execute_worklist(repo, census_ledger, geo_root, output, family_manifest=None, resume=False):
    """Expand an immutable census into an exhaustive queue, never scientific acceptance.

    A literal comparison groups saved labels within one artifact, not verified
    biological membership. Source families only share design work; no new split
    or temporal eligibility is inferred. Completed resumes rehash every input.
    """
    if output.exists():
        if resume:
            return _worklist_resume(output, census_ledger, geo_root, family_manifest)
        raise FileExistsError("Refusing to overwrite worklist: " + str(output))
    if resume:
        raise ValueError("Cannot resume an absent worklist")
    started = time.monotonic()
    census_summary = census_ledger.parent / "summary.json"
    census = json.loads(census_summary.read_text())
    if census.get("schema") != "clockbase-saved-analysis-census-v1" or census.get("status") != "PASS":
        raise ValueError("A completed saved-analysis census is required")
    if census.get("outputs", {}).get(census_ledger.name) != sha256(census_ledger):
        raise ValueError("Census ledger is not bound to its completed summary")
    with census_ledger.open() as handle:
        artifacts = [json.loads(line) for line in handle]
    if len(artifacts) != census["denominators"]["artifacts"]:
        raise ValueError("Census artifact denominator differs from the ledger")
    if len({r["path"] for r in artifacts}) != len(artifacts):
        raise ValueError("Duplicate artifact paths in census")
    families, family_binding, inputs = _worklist_family_map(family_manifest)
    inputs += [_worklist_file(census_ledger, "census_ledger"), _worklist_file(census_summary, "census_summary")]
    output.mkdir(parents=True, exist_ok=False)
    manifest = {"schema": WORKLIST_SCHEMA, "status": "RUNNING",
                "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
                "implementation_sha256": sha256(Path(__file__)), "python": sys.version,
                "arguments": {"census_ledger": str(census_ledger), "geo_root": str(geo_root),
                              "family_manifest": str(family_manifest) if family_manifest else None},
                "runner_command": sys.argv, "census_ledger_sha256": sha256(census_ledger),
                "family_binding": family_binding, "stage_vocabulary": list(WORKLIST_STAGES),
                "scientific_acceptance": "CANNOT VERIFY", "model_calls": 0, "clock_executions": 0,
                "statistical_reruns": 0, "source_repairs": 0,
                "limits": {"max_artifacts": 100000, "max_input_bytes": 20000000000, "wall_seconds": 1800},
                "scope": "All census artifacts and parsed logical rows; no scientific review imported or performed.",
                "row_identity": "SHA256 of artifact path, artifact SHA256 and logical CSV record number (header=1)",
                "comparison_identity": "Within artifact: literal gse_id, comparison_name, condition_name, condition_category and covariate; clocks, methods and duplicate rows retained",
                "grouping_limit": "Literal labels are not sample groups, biological units, independent studies or a multiplicity family.",
                "source_input_limit": "Current companion files are located and hashed, not established as historical analysis inputs or valid clock predictions.",
                "prior_review_limit": "Earlier study-level pilots are not automatically imported as exact-row review, rerun or verification events.",
                "resume_scope": "Completed output only, with all source and output hashes checked; interrupted outputs are preserved and require a new directory."}
    try:
        manifest["commit"] = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True,
                                                      stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        manifest["commit"] = None
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    artifact_rows, comparisons, studies, row_counts = [], {}, {}, Counter()
    parsed_rows = 0

    def canonical(value):
        return bool(isinstance(value, str) and re.fullmatch(r"GSE[0-9]+", value, flags=re.ASCII))

    def register_study(modality, study):
        key = modality + ":" + study
        if key not in studies:
            studies[key] = {"study_key": key, "study_id": study, "modality": modality,
                            "source_family": families.get(study), "artifact_ids": set(),
                            "row_candidates": 0, "identity_bound_rows": 0, "comparison_ids": set()}
        return key

    def rows():
        nonlocal parsed_rows
        for artifact in artifacts:
            path = repo / artifact["path"]
            if not path.resolve().is_relative_to(repo.resolve()):
                raise ValueError("Census artifact escapes repository: " + artifact["path"])
            if not path.is_file():
                raise FileNotFoundError("Cannot proceed - required data file not found: " + str(path))
            input_record = _worklist_file(path, "saved_analysis")
            if input_record["sha256"] != artifact["sha256"]:
                raise ValueError("Saved CSV changed since the census: " + artifact["path"])
            inputs.append(input_record)
            modality = WORKLIST_MODALITY[artifact["root"]]
            artifact_id = _worklist_id("artifact", [artifact["path"], artifact["sha256"]])
            parent = artifact["parent_study_id"]
            parent_key = register_study(modality, parent) if canonical(parent) else None
            if parent_key:
                studies[parent_key]["artifact_ids"].add(artifact_id)
            item = {"artifact_id": artifact_id, "path": artifact["path"], "sha256": artifact["sha256"],
                    "root": artifact["root"], "modality": modality, "parent_study_id": parent,
                    "census_structural_status": artifact["status"], "census_issues": artifact["issues"],
                    "census_rows": artifact["rows"], "parsed_rows": 0, "comparison_ids": set(),
                    "header": [], "parse_error": None, "stages": _worklist_stages(),
                    "scientific_acceptance": "CANNOT VERIFY"}
            artifact_rows.append(item)
            if _is_lfs_pointer(path):
                item["parse_error"] = "LFS_POINTER"
                continue
            try:
                with path.open(newline="", encoding="utf-8") as handle:
                    reader = csv.reader(handle, strict=True)
                    item["header"] = header = next(reader, [])
                    positions = {name: i for i, name in enumerate(header) if header.count(name) == 1}
                    for logical, row in enumerate(reader, 2):
                        parsed_rows += 1
                        item["parsed_rows"] += 1
                        def cell(name):
                            i = positions.get(name)
                            return row[i] if i is not None and i < len(row) else None
                        row_id = _worklist_id("row", [artifact["path"], artifact["sha256"], logical])
                        literal = {name: cell(name) for name in
                                   ("gse_id", "comparison_name", "condition_name", "condition_category", "covariate")}
                        comparison_id = _worklist_id("comparison", [artifact_id, literal])
                        study = literal["gse_id"]
                        study_key = register_study(modality, study) if canonical(study) else None
                        candidates = sorted({key for key in (parent_key, study_key) if key})
                        identity_bound = bool(study_key and parent_key == study_key and len(row) == len(header))
                        identity_issues = []
                        if not canonical(parent): identity_issues.append("INVALID_PARENT_STUDY_ID")
                        if not canonical(study): identity_issues.append("INVALID_OR_AMBIGUOUS_ROW_STUDY_ID")
                        if study != parent: identity_issues.append("PARENT_ROW_STUDY_ID_CONTRADICTION")
                        if len(row) != len(header): identity_issues.append("RAGGED_ROW")
                        if len(header) != len(set(header)): identity_issues.append("DUPLICATE_HEADER")
                        if "DUPLICATE_HEADER" in identity_issues: identity_bound = False
                        for key in candidates:
                            studies[key]["artifact_ids"].add(artifact_id)
                            studies[key]["row_candidates"] += 1
                            studies[key]["comparison_ids"].add(comparison_id)
                        if identity_bound:
                            studies[study_key]["identity_bound_rows"] += 1
                        record = {"row_id": row_id, "artifact_id": artifact_id,
                                  "artifact_sha256": artifact["sha256"], "logical_csv_record": logical,
                                  "physical_line_end": reader.line_num,
                                  "row_sha256": hashlib.sha256(json.dumps(row, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest(),
                                  "comparison_id": comparison_id, "literal_comparison": literal,
                                  "clock": cell("clock"), "statistical_test": cell("statistical_test"),
                                  "model_version": cell("model_version"), "source_root": artifact["root"],
                                  "candidate_study_keys": candidates,
                                  "bound_study_key": study_key if identity_bound else None,
                                  "identity_status": "PASS" if identity_bound else "UNRESOLVED",
                                  "identity_issues": identity_issues, "stages": _worklist_stages(),
                                  "scientific_acceptance": "CANNOT VERIFY"}
                        comp = comparisons.setdefault(comparison_id, {
                            "comparison_id": comparison_id, "artifact_id": artifact_id,
                            "literal_comparison": literal, "modality": modality,
                            "source_root": artifact["root"], "row_ids": [], "clocks": set(),
                            "statistical_tests": set(), "model_versions": set(), "identity_issues": set(),
                            "candidate_study_keys": set(), "bound_study_keys": set(), "stages": _worklist_stages(),
                            "scientific_acceptance": "CANNOT VERIFY"})
                        comp["row_ids"].append(row_id)
                        for target, source in (("clocks", "clock"), ("statistical_tests", "statistical_test"), ("model_versions", "model_version")):
                            if record[source] is not None: comp[target].add(record[source])
                        comp["identity_issues"].update(identity_issues)
                        comp["candidate_study_keys"].update(candidates)
                        if identity_bound: comp["bound_study_keys"].add(study_key)
                        item["comparison_ids"].add(comparison_id)
                        row_counts["identity_bound" if identity_bound else "identity_unresolved"] += 1
                        row_counts["clock_label_present" if cell("clock") else "clock_label_missing"] += 1
                        yield record
            except (csv.Error, UnicodeError) as exc:
                item["parse_error"] = type(exc).__name__ + ": " + str(exc)
            if item["parsed_rows"] != artifact["rows"]:
                raise ValueError("Parsed logical row count differs from census: " + artifact["path"])
            if time.monotonic() - started > manifest["limits"]["wall_seconds"]:
                raise TimeoutError("Worklist wall-time limit exceeded")

    try:
        if len(artifacts) > manifest["limits"]["max_artifacts"]:
            raise ValueError("Worklist artifact cap exceeded")
        _write_gz_jsonl(output / "clock-row-worklist.jsonl.gz", rows())
        if parsed_rows != census["denominators"]["rows"]:
            raise ValueError("Worklist logical row denominator differs from census")
        source_counts = Counter()
        for study in studies.values():
            folder = geo_root / WORKLIST_GEO_FOLDERS[study["modality"]]
            source_files = {}
            for role, directory in (("source_metadata", "metadata"), ("processed_metadata", "processed_metadata"),
                                    ("companion_predictions", "clock_result")):
                record = _worklist_file(folder / directory / (study["study_id"] + ".csv"), role)
                inputs.append(record)
                source_files[role] = {key: record[key] for key in ("path", "available", "sha256", "bytes", "payload_status", "reason")}
                source_counts[(study["modality"], role, record["reason"])] += 1
            study.update(source_files=source_files, stages=_worklist_stages(),
                         ready_for_source_design_preparation=source_files["source_metadata"]["payload_status"] == "PASS",
                         ready_for_statistical_rerun=False, scientific_acceptance="CANNOT VERIFY")
            study["dependencies"] = list(WORKLIST_DEPENDENCIES)
            for role in ("source_metadata", "companion_predictions"):
                if source_files[role]["payload_status"] != "PASS":
                    study["dependencies"].append(role.upper() + "_" + source_files[role]["reason"])
            for key in ("artifact_ids", "comparison_ids"):
                study[key] = sorted(study[key])
            if time.monotonic() - started > manifest["limits"]["wall_seconds"]:
                raise TimeoutError("Worklist source inventory wall-time limit exceeded")
        if sum(record["bytes"] or 0 for record in inputs) > manifest["limits"]["max_input_bytes"]:
            raise ValueError("Worklist source input byte cap exceeded")
        for comparison in comparisons.values():
            for key in ("clocks", "statistical_tests", "model_versions", "identity_issues", "candidate_study_keys", "bound_study_keys"):
                comparison[key] = sorted(comparison[key])
            bound = comparison["bound_study_keys"]
            unambiguous = len(bound) == 1 and not comparison["identity_issues"]
            comparison["source_design_group_key"] = bound[0] if unambiguous else None
            comparison["source_family"] = studies[bound[0]]["source_family"] if unambiguous else None
            comparison["ready_for_source_design_preparation"] = bool(unambiguous and studies[bound[0]]["ready_for_source_design_preparation"])
            comparison["ready_for_statistical_rerun"] = False
            comparison["execution_disposition"] = "UNRESOLVED"
            comparison["dependencies"] = list(studies[bound[0]]["dependencies"]) if unambiguous else list(WORKLIST_DEPENDENCIES) + ["SOURCE_STUDY_BINDING_UNRESOLVED"]
        for item in artifact_rows:
            item["comparison_ids"] = sorted(item["comparison_ids"])
            item["execution_disposition"] = "UNRESOLVED"
            item["dependencies"] = (["NO_PARSED_RESULT_ROWS"] if item["parsed_rows"] == 0 else []) + list(WORKLIST_DEPENDENCIES)
        _write_gz_jsonl(output / "artifact-worklist.jsonl.gz", artifact_rows)
        _write_gz_jsonl(output / "comparison-worklist.jsonl.gz", (comparisons[key] for key in sorted(comparisons)))
        _write_gz_jsonl(output / "study-inputs.jsonl.gz", (studies[key] for key in sorted(studies)))
        for record in inputs:
            _recheck_worklist_file(record)
        if sha256(Path(__file__)) != manifest["implementation_sha256"]:
            raise ValueError("Worklist implementation changed during execution")
        _write_gz_jsonl(output / "input-files.jsonl.gz", inputs)
        denominator = {"artifacts": len(artifact_rows), "zero_row_artifacts": sum(r["parsed_rows"] == 0 for r in artifact_rows),
                       "parsed_logical_rows": parsed_rows, "clock_labelled_rows": row_counts["clock_label_present"],
                       "clock_label_missing_rows": row_counts["clock_label_missing"],
                       "literal_comparisons_within_artifacts": len(comparisons),
                       "canonical_candidate_study_accessions": len({r["study_id"] for r in studies.values()}),
                       "canonical_candidate_study_by_biological_modality": len(studies),
                       "identity_bound_rows": row_counts["identity_bound"],
                       "identity_unresolved_rows": row_counts["identity_unresolved"],
                       "source_family_groups_located": len({r["source_family"]["family_id"] for r in studies.values() if r["source_family"]}),
                       "study_modality_without_source_family": sum(r["source_family"] is None for r in studies.values()),
                       "source_design_preparation_ready_study_modalities": sum(r["ready_for_source_design_preparation"] for r in studies.values()),
                       "source_design_preparation_ready_comparisons": sum(r["ready_for_source_design_preparation"] for r in comparisons.values()),
                       "statistical_rerun_ready_comparisons": 0,
                       "reviewed_comparisons_in_this_worklist": 0, "rerun_comparisons_in_this_worklist": 0,
                       "independently_verified_comparisons_in_this_worklist": 0}
        names = ("clock-row-worklist.jsonl.gz", "artifact-worklist.jsonl.gz", "comparison-worklist.jsonl.gz", "study-inputs.jsonl.gz", "input-files.jsonl.gz")
        manifest.update(status="PASS", denominators=denominator,
                        input_files=len(inputs), input_bytes=sum(r["bytes"] or 0 for r in inputs),
                        input_hashes_rechecked=sum(r["available"] for r in inputs),
                        missing_input_paths_rechecked=sum(not r["available"] for r in inputs),
                        source_counts=[{"modality": k[0], "role": k[1], "status": k[2], "study_modalities": v} for k, v in sorted(source_counts.items())],
                        outputs={name: sha256(output / name) for name in names},
                        elapsed_seconds=time.monotonic() - started)
    except Exception as exc:
        manifest.update(status="FAIL", error=type(exc).__name__ + ": " + str(exc),
                        completed_artifacts=len(artifact_rows), completed_rows=parsed_rows,
                        elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        (output / "summary.json").write_text(json.dumps(manifest, indent=2) + "\n")
        (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": manifest["status"], "output": str(output), "denominators": denominator}, indent=2), flush=True)
    return manifest


def execute_reference_cache(repo, plan_path, cache, output):
    """Run the corrected reference reader and compare, without running clocks."""
    import numpy as np
    sys.path.insert(0, str(repo))
    from src.pipeline_inputs_2026_09_10 import reference_means, MATRIX_READER_CONTRACT
    plan = json.loads(plan_path.read_text())
    if plan.get('schema') != 'clockbase-reference-cache-readback-v1':
        raise ValueError('Expected a pinned reference-cache readback plan')
    for item in plan['inputs']:
        if sha256(Path(item['path'])) != item['sha256']:
            raise ValueError('Reference readback input changed: ' + item['path'])
    sources = {item['role']: item for item in plan['inputs']}
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    result = {'schema': plan['schema'], 'status': 'RUNNING', 'plan_sha256': sha256(plan_path),
              'reader_contract': MATRIX_READER_CONTRACT, 'recorded_at_utc': datetime.now(timezone.utc).isoformat(),
              'implementation_sha256': sha256(repo / 'src/pipeline_inputs_2026_09_10.py'),
              'runner_sha256': sha256(Path(__file__)), 'empirical_clock_executions': 0, 'paid_model_api_calls': 0,
              'scientific_status': 'CANNOT VERIFY'}
    try:
        values, record = reference_means(Path(sources['reference_source']['path']), cache)
        prior = pd.read_csv(sources['previous_means']['path'], index_col=0, converters={0: str},
                            float_precision='round_trip')['mean']
        same_ids = list(values.index) == list(prior.index)
        if not same_ids: raise ValueError('Reference mean feature membership/order changed')
        a, b = values.to_numpy(dtype=float), prior.to_numpy(dtype=float)
        difference = np.abs(a - b)
        parity = bool(np.allclose(a, b, rtol=plan['rtol'], atol=plan['atol'], equal_nan=False))
        target = cache / ('GSE40279-means-' + record['source_sha256'][:16] + '-' + MATRIX_READER_CONTRACT + '.csv')
        resumed, resumed_record = reference_means(Path(sources['reference_source']['path']), cache)
        if not resumed.equals(values) or record != resumed_record:
            raise ValueError('Versioned reference-cache resume differs from computation')
        for item in plan['inputs']:
            if sha256(Path(item['path'])) != item['sha256']:
                raise ValueError('Reference input changed during readback: ' + item['path'])
        result.update(status='PASS' if parity else 'FAIL', source_record=record,
                      prior_feature_order_equal=same_ids, prior_numeric_parity=parity,
                      reference_features=len(values), exactly_equal_values=int(np.sum(a == b)),
                      max_abs_difference=float(difference.max()), rtol=plan['rtol'], atol=plan['atol'],
                      versioned_cache_path=str(target), versioned_cache_sha256=sha256(target),
                      old_source_and_cache_unchanged=True, completed_resume_equal=True,
                      elapsed_seconds=time.monotonic() - started)
    except Exception as exc:
        result.update(status='FAIL', error=type(exc).__name__ + ': ' + str(exc), elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        with (output / 'receipt.json').open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps(result, indent=2), flush=True)


def verify_geo_reader(plan_path, root_summary_path, output):
    """Exercise the separate GEO producer parser without importing its jobs."""
    import numpy as np
    plan = json.loads(plan_path.read_text()); expected = json.loads(root_summary_path.read_text())
    if expected['plan_sha256'] != sha256(plan_path) or expected['status'] != 'PASS':
        raise ValueError('Root source audit is absent or not bound to this plan')
    source = Path(plan['geo_root']) / 'src/4_compute_age.py'
    before = sha256(source)
    tree = ast.parse(source.read_text())
    names = {'require_observed_beta_matrix', '_sha256_file', '_open_geo_matrix', '_scan_geo_matrix', 'read_gse_matrix'}
    constants = {'GEO_TABLE_BEGIN', 'GEO_TABLE_END', 'GEO_SAMPLE_ACCESSION'}
    nodes = [node for node in tree.body if
             (isinstance(node, ast.FunctionDef) and node.name in names)
             or (isinstance(node, ast.Assign) and len(node.targets) == 1
                 and isinstance(node.targets[0], ast.Name) and node.targets[0].id in constants)]
    if {node.name for node in nodes if isinstance(node, ast.FunctionDef)} != names:
        raise ValueError('GEO source parser functions changed; inspect before loading')
    namespace = {'np': np, 'pd': pd, 'gzip': gzip, 'csv': csv, 're': re, 'hashlib': hashlib}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), 'exec'), namespace)
    output.mkdir(parents=True, exist_ok=False)
    result = {'schema': 'clockbase-geo-reader-parity-v1', 'status': 'RUNNING',
              'plan_sha256': sha256(plan_path), 'root_summary_sha256': sha256(root_summary_path),
              'producer_source_sha256': before, 'runner_sha256': sha256(Path(__file__)),
              'recorded_at_utc': datetime.now(timezone.utc).isoformat(), 'checks': [],
              'empirical_clock_executions': 0, 'model_calls': 0, 'network_calls': 0,
              'scientific_status': 'CANNOT VERIFY'}
    inputs = {(r['study_id'], r['role']): r for r in plan['inputs']}
    started = time.monotonic()
    try:
        for target in expected['dnam']:
            study = target['study_id']; source_record = inputs[study, 'series_matrix']
            path = Path(source_record['path'])
            if sha256(path) != source_record['sha256']: raise ValueError('Frozen matrix changed')
            frame = namespace['read_gse_matrix'](str(path))
            digest, counts = numeric_stream_digest(frame)
            try:
                namespace['require_observed_beta_matrix'](frame); beta_status = 'PASS'
            except ValueError:
                beta_status = 'FAIL'
            checks = {'source_columns': list(frame.columns) == target['source_sample_ids'],
                      'source_features': len(frame) == target['source_features'],
                      'complete_numeric_and_identity_digest': digest == target['numeric_stream_sha256'],
                      'numeric_counts': counts == target['numeric_counts'],
                      'beta_input_contract': beta_status == target['beta_input_contract']['status'],
                      'source_unchanged': sha256(path) == source_record['sha256']}
            result['checks'].append({'study_id': study, 'checks': checks,
                                     'status': 'PASS' if all(checks.values()) else 'FAIL',
                                     'numeric_stream_sha256': digest, 'beta_input_status': beta_status})
            print('Verified GEO reader ' + study, flush=True)
        if sha256(source) != before: raise ValueError('GEO source code changed during verification')
        result.update(status='PASS' if all(r['status'] == 'PASS' for r in result['checks']) else 'FAIL',
                      source_code_unchanged=True, elapsed_seconds=time.monotonic() - started)
    except Exception as exc:
        result.update(status='FAIL', error=type(exc).__name__ + ': ' + str(exc), elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        with (output / 'receipt.json').open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    if result['status'] != 'PASS': raise ValueError('GEO reader and root source values differ')


def execute_identity_origin(repo, plan_path, output):
    """Post-hoc source diagnosis; preserve predictions and never infer IDs."""
    sys.path.insert(0, str(repo))
    from src.pipeline_inputs_2026_09_10 import read_matrix, _inspect_matrix_stream, MATRIX_READER_CONTRACT
    plan = json.loads(plan_path.read_text())
    if plan.get('schema') != 'clockbase-sample-identity-origin-v1' or plan.get('scope') != 'post_hoc_source_audit':
        raise ValueError('Expected pinned post-hoc identity-origin plan')
    started = time.monotonic()
    inputs = {(r['study_id'], r['role']): r for r in plan['inputs']}
    if len(inputs) != len(plan['inputs']): raise ValueError('Duplicate source roles')
    if sum(r.get('size_bytes', 0) for r in plan['inputs']) > plan['resources']['max_source_bytes']:
        raise ValueError('Source byte limit exceeded')
    for record in plan['inputs']:
        path = Path(record['path'])
        if record['available']:
            if sha256(path) != record['sha256']: raise ValueError('Source changed: ' + str(path))
        elif path.exists(): raise ValueError('Previously missing input now exists; revise the plan explicitly')
    output.mkdir(parents=True, exist_ok=False)
    result = {'schema': 'clockbase-sample-identity-origin-result-v1', 'status': 'RUNNING',
              'recorded_at_utc': datetime.now(timezone.utc).isoformat(),
              'scope': plan['scope'], 'plan_sha256': sha256(plan_path),
              'commit': subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip(),
              'reader_contract': MATRIX_READER_CONTRACT, 'python': sys.version, 'pandas': pd.__version__,
              'implementation_sha256': sha256(repo / 'src/pipeline_inputs_2026_09_10.py'),
              'runner_sha256': sha256(Path(__file__)), 'dnam': [], 'rna': [],
              'source_repairs': 0, 'clock_executions': 0, 'paid_model_api_calls': 0,
              'scientific_status': 'CANNOT VERIFY',
              'limits': ['Post-hoc selected source faults, not held-out accuracy or a scientific error rate.',
                         'Reproducing an error mechanism does not prove the historical producing runtime.',
                         'Renaming samples cannot restore dropped probes or truncated measurements.',
                         'No saved prediction, source measurement, figure coordinate or frozen benchmark is changed.']}
    try:
        for study in plan['dnam_studies']:
            print('Auditing DNAm source ' + study, flush=True)
            matrix = Path(inputs[study, 'series_matrix']['path'])
            clock = Path(inputs[study, 'predictions']['path'])
            header, rows = strict_csv_records(clock)
            saved_ids = [row[0] for row in rows]
            with gzip.open(matrix, 'rt', newline='') as handle:
                for physical, line in enumerate(handle, 1):
                    if line.strip() == '!series_matrix_table_begin': break
                else: raise ValueError('No physical begin marker in ' + str(matrix))
            with gzip.open(matrix, 'rt', newline='') as handle:
                source = _inspect_matrix_stream(handle, geo=True, delimiter='\t')
                handle.seek(source['data_start'])
                raw = pd.read_csv(handle, sep='\t', names=source['header'], header=None, index_col=0,
                                  converters={0: str}, nrows=len(source['features']), float_precision='round_trip')
            if list(raw.index) != source['features'] or list(raw.columns) != source['header'][1:]:
                raise ValueError('Corrected parsing differs from literal source identity')
            old_header = list(pd.read_csv(matrix, sep='\t', skiprows=physical, nrows=0).columns)
            with gzip.open(matrix, 'rt', newline='') as handle:
                reader = csv.reader(handle, delimiter='\t', strict=True)
                for logical, row in enumerate(reader, 1):
                    if row == ['!series_matrix_table_begin']: break
            lost = physical - logical
            if lost < 0: raise ValueError('Unexpected physical/logical source relation')
            digest, counts = numeric_stream_digest(raw)
            metadata_header, metadata = strict_csv_records(Path(inputs[study, 'metadata']['path']))
            column = metadata_header.index('GSM'); meta_ids = [row[column] for row in metadata]
            try:
                checked = read_matrix(matrix, beta=True)
                measurement = {'status': 'PASS', 'reason': 'Existing beta input contract passed'}
                if not raw.equals(checked): raise AssertionError('Active reader differs from inspected stream')
            except ValueError as exc:
                measurement = {'status': 'FAIL', 'reason': str(exc)}
            record = {'study_id': study, 'source_sha256': inputs[study, 'series_matrix']['sha256'],
                      'saved_prediction_sha256': inputs[study, 'predictions']['sha256'],
                      'source_samples': len(source['header']) - 1, 'source_features': len(source['features']),
                      'source_sample_ids': source['header'][1:], 'saved_sample_labels': saved_ids,
                      'source_header_sha256': hashlib.sha256((json.dumps(source['header'], ensure_ascii=False, separators=(',', ':')) + '\n').encode()).hexdigest(),
                      'source_feature_order_sha256': hashlib.sha256(''.join(json.dumps(v, ensure_ascii=False, separators=(',', ':')) + '\n' for v in source['features']).encode()).hexdigest(),
                      'numeric_stream_sha256': digest, 'numeric_counts': counts,
                      'preamble_physical_records': physical, 'preamble_logical_records': logical,
                      'old_reader_header': old_header, 'old_header_equals_saved_labels': old_header[1:] == saved_ids,
                      'old_reader_lost_feature_count': lost, 'old_reader_lost_feature_ids': source['features'][:lost],
                      'canonical_metadata_order_matches_header': meta_ids == source['header'][1:],
                      'canonical_metadata_membership_matches_header': set(meta_ids) == set(source['header'][1:]),
                      'source_parser_status': 'PASS', 'beta_input_contract': measurement,
                      'saved_prediction_identity_status': 'FAIL', 'prediction_repair_status': 'CANNOT VERIFY'}
            result['dnam'].append(record)
            (output / (study + '-source.json')).write_text(json.dumps(record, indent=2) + '\n')
            del raw
            if time.monotonic() - started > plan['resources']['max_wall_seconds_per_run']:
                raise TimeoutError('Pinned source-audit wall limit exceeded')
        for study in plan['rna_studies']:
            print('Auditing RNA source ' + study, flush=True)
            tables = []
            for role in ('scaled_data', 'scaled_diff_data', 'yugene_data', 'yugene_diff_data'):
                path = Path(inputs[study, role]['path']); ids = []; ragged = []
                with path.open(newline='') as handle:
                    reader = csv.reader(handle, strict=True); header = next(reader)
                    for n, row in enumerate(reader, 1):
                        ids.append(row[0] if row else '')
                        if len(row) != len(header):
                            ragged.append({'csv_data_record': n, 'sample_id': row[0] if row else '',
                                           'fields': len(row), 'expected_fields': len(header)})
                with path.open('rb') as handle:
                    handle.seek(-1, 2); final_newline = handle.read(1) in (b'\n', b'\r')
                tables.append({'role': role, 'source_sha256': inputs[study, role]['sha256'],
                               'rows': len(ids), 'fields': len(header), 'sample_ids': ids,
                               'invalid_ids': [v for v in ids if not re.fullmatch(r'GSM[0-9]+', v)],
                               'duplicate_id_excess': len(ids) - len(set(ids)), 'ragged_records': ragged,
                               'ends_in_newline': final_newline, 'csv_status': 'FAIL' if ragged else 'PASS'})
            header, rows = strict_csv_records(Path(inputs[study, 'predictions']['path']))
            saved_ids = [row[0] for row in rows]
            blank_records = [n for n, row in enumerate(rows, 1) if not row[0]]
            patterns = Counter(tuple(header[i] for i, value in enumerate(row[1:], 1) if value == '') for row in rows)
            compatible = []
            for table in tables:
                # Demonstrate the old concat's identity projection only; do not assign IDs to outputs.
                projected = table['sample_ids'] + [''] * (max(t['rows'] for t in tables) - table['rows'])
                if projected == saved_ids: compatible.append(table['role'])
            record = {'study_id': study, 'preprocessing': tables, 'prediction_rows': len(rows),
                      'saved_prediction_sha256': inputs[study, 'predictions']['sha256'],
                      'blank_id_records': blank_records,
                      'missing_clock_patterns': [{'records': count, 'missing_columns': list(names)} for names, count in patterns.items()],
                      'legacy_position_concat_compatible_first_inputs': compatible,
                      'metadata_available': inputs[study, 'metadata']['available'],
                      'input_integrity_status': 'FAIL', 'saved_prediction_identity_status': 'FAIL',
                      'prediction_repair_status': 'CANNOT VERIFY',
                      'historical_execution_binding': 'CANNOT VERIFY'}
            result['rna'].append(record)
            (output / (study + '-source.json')).write_text(json.dumps(record, indent=2) + '\n')
        for record in plan['inputs']:
            if record['available'] and sha256(Path(record['path'])) != record['sha256']:
                raise ValueError('Source changed during audit: ' + record['path'])
        result.update(status='PASS', source_hashes_rechecked=sum(r['available'] for r in plan['inputs']),
                      elapsed_seconds=time.monotonic() - started)
        with (output / 'validation-table.csv').open('x', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['study_id', 'modality', 'prediction_rows', 'invalid_identity_rows', 'observed_fault',
                             'parser_or_csv_status', 'prediction_repair_status'])
            for row in result['dnam']:
                writer.writerow([row['study_id'], 'human_dnam', len(row['saved_sample_labels']), len(row['saved_sample_labels']),
                                 'multiline_metadata_header_shift', row['beta_input_contract']['status'], 'CANNOT VERIFY'])
            for row in result['rna']:
                writer.writerow([row['study_id'], 'mouse_rna', row['prediction_rows'], len(row['blank_id_records']),
                                 'truncated_preprocessing_and_position_concat', 'FAIL', 'CANNOT VERIFY'])
    except Exception as exc:
        result.update(status='FAIL', error=type(exc).__name__ + ': ' + str(exc), elapsed_seconds=time.monotonic() - started)
        raise
    finally:
        with (output / 'summary.json').open('x') as handle: json.dump(result, handle, indent=2); handle.write('\n')
    print(json.dumps({'status': result['status'], 'studies': len(result['dnam']) + len(result['rna']),
                      'elapsed_seconds': result['elapsed_seconds'], 'output': str(output)}), flush=True)


def reader_from_source(path):
    """Only the source's boundary function executes, with explicit dependencies."""
    tree = ast.parse(path.read_text())
    names = {"read_result_csv", "read_cached_result_csv"}
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                 and node.name in names]
    if {node.name for node in functions} != names or len(functions) != len(names):
        raise ValueError("Expected both raw and cached result boundary definitions")
    namespace = {"csv": csv, "io": io, "Path": Path, "re": re, "pd": pd}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace["read_result_csv"], namespace["read_cached_result_csv"]


def source_paths(root):
    def raise_error(error):
        raise error
    for directory, subdirectories, files in os.walk(root, onerror=raise_error):
        subdirectories.sort()
        if "aging_analysis.csv" in files:
            yield Path(directory) / "aging_analysis.csv"


def execute(repo, output):
    started = time.monotonic()
    implementation = repo / "s01_visualize_results.py"
    reader, cache_reader = reader_from_source(implementation)
    roots = []
    paths = []
    for name in ROOTS:
        path = repo / name
        if not path.is_dir():
            roots.append({"path": name, "status": "CANNOT VERIFY", "files": 0,
                          "message": f"Cannot proceed - required data file not found: {path}"})
            continue
        entries = list(source_paths(path))
        paths.extend((entry, name, "study_result") for entry in entries)
        roots.append({"path": name, "status": "PRESENT", "files": len(entries)})
    # These are existing saved tables, not source-authenticated caches.
    caches = sorted((repo / "result").glob("combined_df_*.csv"))
    paths.extend((path, "result", "combined_cache") for path in caches)
    total_bytes = sum(path.stat().st_size for path, _, _ in paths)
    if len(paths) > LIMITS["files"] or total_bytes > LIMITS["bytes"]:
        raise ValueError("Input inventory exceeds the frozen file/byte limits")
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema": "clockbase-legacy-result-census-v1", "status": "RUNNING",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "repository": str(repo),
        "commit": subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        "implementation_sha256": sha256(implementation), "runner_sha256": sha256(Path(__file__)),
        "python": sys.version, "pandas": pd.__version__, "platform": platform.platform(),
        "limits": LIMITS, "roots": roots, "inventory_files": len(paths), "inventory_bytes": total_bytes,
        "scientific_status": "CANNOT VERIFY",
        "scope": "Artifact-level CSV structure, required columns and exact accession binding. First failure only; no numerical filtering, source authentication, biological approval or error-rate estimate.",
        "network_requests": 0, "paid_model_api_calls": 0, "clock_executions": 0,
        "top_level_plotting_executed": False, "source_repairs": 0,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    counts = Counter()
    records = []
    try:
        with (output / "source-ledger.jsonl").open("x") as handle:
            for path, root, kind in paths:
                if time.monotonic() - started > LIMITS["wall_seconds"]:
                    raise TimeoutError("Frozen census wall-time limit exceeded")
                stat = path.stat()
                record = {"path": str(path.relative_to(repo)), "root": root, "kind": kind,
                          "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns, "sha256": sha256(path),
                          "parent_study_id": path.parent.name if kind == "study_result" else None}
                kwargs = {"expected_study_id": path.parent.name, "required_columns": REQUIRED} if kind == "study_result" else {}
                try:
                    frame = reader(path, **kwargs) if kind == "study_result" else cache_reader(path)
                    record.update(status="PASS", reason="BOUNDARY_PASSED", message="",
                                  rows=len(frame), study_ids=sorted(frame["gse_id"].unique().tolist()))
                except ValueError as exc:
                    reason = str(exc).split(":", 1)[0]
                    record.update(status="CANNOT VERIFY" if reason in {"PARENT_ID_UNRESOLVED", "LFS_POINTER"} else "FAIL",
                                  reason=reason, message=str(exc), rows=None, study_ids=None)
                records.append(record)
                counts[(root, kind, record["status"], record["reason"])] += 1
                handle.write(json.dumps(record, ensure_ascii=False) + "\n")
        # Verify the complete inventory again; input mutation invalidates the run.
        for record in records:
            if sha256(repo / record["path"]) != record["sha256"]:
                raise ValueError(f"Source changed during census: {record['path']}")
        current_paths = [(str(path.relative_to(repo)), name, "study_result")
                         for name in ROOTS if (repo / name).is_dir() for path in source_paths(repo / name)]
        current_paths += [(str(path.relative_to(repo)), "result", "combined_cache")
                          for path in sorted((repo / "result").glob("combined_df_*.csv"))]
        if current_paths != [(r["path"], r["root"], r["kind"]) for r in records]:
            raise ValueError("Source inventory changed during census")
        if sha256(implementation) != manifest["implementation_sha256"]:
            raise ValueError("Reader implementation changed during census")
        with (output / "validation-table.csv").open("x", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["root", "artifact_kind", "boundary_status", "first_failure", "files"])
            writer.writerows([*key, count] for key, count in sorted(counts.items()))
        manifest.update(status="PASS", elapsed_seconds=time.monotonic() - started,
                        input_hashes_rechecked=len(records), input_inventory_rechecked=True,
                        counts=[dict(zip(("root", "kind", "status", "reason", "files"), (*key, count)))
                                for key, count in sorted(counts.items())],
                        outputs={name: sha256(output / name) for name in ("source-ledger.jsonl", "validation-table.csv")})
    except Exception as exc:
        manifest.update(status="FAIL", error=f"{type(exc).__name__}: {exc}",
                        elapsed_seconds=time.monotonic() - started, completed_records=len(records))
        raise
    finally:
        (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": manifest["status"], "files": len(records), "output": str(output),
                      "counts": manifest["counts"]}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--identity-origin-plan", type=Path,
                        help="Run the pinned post-hoc molecular/sample-identity source audit instead of the legacy result census")
    parser.add_argument('--reference-cache-plan', type=Path)
    parser.add_argument('--reference-cache', type=Path)
    parser.add_argument('--geo-reader-plan', type=Path)
    parser.add_argument('--root-source-summary', type=Path)
    parser.add_argument('--saved-corpus-plan', type=Path,
                        help='Run the bounded saved-output census against every aging_analysis*.csv under all frozen roots')
    parser.add_argument('--worklist-census', type=Path,
                        help='Expand a completed artifact-ledger.jsonl into an exhaustive scientific-audit queue')
    parser.add_argument('--worklist-geo-root', type=Path)
    parser.add_argument('--worklist-family-manifest', type=Path)
    parser.add_argument('--resume', action='store_true', help='Revalidate a completed worklist without overwriting it')
    args = parser.parse_args()
    try:
        if sum(bool(value) for value in (args.reference_cache_plan, args.geo_reader_plan,
                                         args.identity_origin_plan, args.saved_corpus_plan, args.worklist_census)) > 1:
            raise ValueError('Choose exactly one optional source-verification mode')
        if (args.resume or args.worklist_geo_root or args.worklist_family_manifest) and not args.worklist_census:
            raise ValueError('Worklist options require --worklist-census')
        if args.worklist_census:
            if args.worklist_geo_root is None:
                raise ValueError('Worklist requires --worklist-geo-root')
            execute_worklist(args.repo.resolve(), args.worklist_census.resolve(), args.worklist_geo_root.resolve(),
                             args.output.resolve(), args.worklist_family_manifest.resolve() if args.worklist_family_manifest else None,
                             args.resume)
        elif args.geo_reader_plan:
            if args.root_source_summary is None: raise ValueError('GEO reader check requires the root source summary')
            verify_geo_reader(args.geo_reader_plan.resolve(), args.root_source_summary.resolve(), args.output.resolve())
        elif args.reference_cache_plan:
            if args.reference_cache is None or args.identity_origin_plan:
                raise ValueError('Reference-cache mode requires its cache path and excludes the source-audit mode')
            execute_reference_cache(args.repo.resolve(), args.reference_cache_plan.resolve(),
                                    args.reference_cache.resolve(), args.output.resolve())
        elif args.identity_origin_plan:
            execute_identity_origin(args.repo.resolve(), args.identity_origin_plan.resolve(), args.output.resolve())
        elif args.saved_corpus_plan:
            execute_saved_corpus(args.repo.resolve(), args.saved_corpus_plan.resolve(), args.output.resolve())
        else:
            execute(args.repo.resolve(), args.output.resolve())
    except (OSError, ValueError, TimeoutError) as exc:
        parser.exit(1, f"Legacy input census failed: {exc}\n")


if __name__ == "__main__":
    main()
