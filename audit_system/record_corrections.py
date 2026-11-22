"""Apply source-adjudicated CSV corrections without changing unrelated bytes.

This module propagates explicit review decisions; it is not an automatic
scientific adjudicator. Dry-run is the default. Complete input/output hashes
make an interrupted multi-file application resumable and repeated runs idempotent.
"""
from __future__ import annotations

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
import logging
import os
from pathlib import Path
import tempfile

from .hashing import inspect_file

LOG = logging.getLogger(__name__)
IDENTITY_FIELDS = {"gse_id", "dataset", "case_id", "comparison_name", "clock"}
# Explicit tier values in questionaire.md; Q7 is context and is not scored.
SCORE_TIERS = {
    "q0_plausibility": {-50, 0, 30},
    "q1_experimental_model": {0, 15, 18, 20, 25},
    "q2_1_aging_pathways": {0, 5, 10, 15},
    "q2_2_age_diseases": {0, 5, 8, 10},
    "q3_1_experimental_rigor": {0, 10},
    "q4_1_translational_potential": {0, 2, 5},
    "q5_1_literature_saturation": {0, 3, 10, 15},
    "q5_2_mechanistic_novelty": {5, 10, 15, 20},
}


class CorrectionError(RuntimeError):
    """A correction cannot be applied to the supplied source version."""


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def read_bytes(path: Path) -> bytes:
    try:
        return path.read_bytes()
    except OSError as exc:
        raise CorrectionError(f"Cannot read correction input {path}: {exc}") from exc


def local_path(root: Path, value: str) -> Path:
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise CorrectionError(f"Expected a repository-relative path: {value}")
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise CorrectionError(f"Correction path escapes repository: {value}")
    return path


def csv_records(raw: bytes):
    """Return fields, rows and physical line spans, preserving quoted newlines."""
    text = raw.decode("utf-8")
    lines = text.splitlines(keepends=True)
    reader = csv.reader(io.StringIO(text, newline=""))
    try:
        fields = next(reader)
        if len(set(fields)) != len(fields):
            raise CorrectionError("Duplicate CSV column names")
        rows = []
        previous = reader.line_num
        for record_number, values in enumerate(reader, 2):
            if len(values) != len(fields):
                raise CorrectionError(f"Malformed CSV record {record_number}")
            rows.append((dict(zip(fields, values)), previous, reader.line_num))
            previous = reader.line_num
    except (StopIteration, csv.Error) as exc:
        raise CorrectionError(f"Cannot parse correction CSV: {exc}") from exc
    return fields, rows, lines


def validate_score_change(row: dict, changes: dict, number: int) -> None:
    """Check new assignments without endorsing unchanged historical scores.

    This checks tier membership and arithmetic only. Source support, novelty,
    applicability and biological efficacy still require separate adjudication.
    """
    touched = set(changes) & (set(SCORE_TIERS) | {"final_score", "trust_score"})
    if not touched:
        return

    def numeric(field, value):
        try:
            result = Decimal(value)
        except (InvalidOperation, ValueError) as exc:
            raise CorrectionError(f"Non-numeric score {field} in record {number}") from exc
        if not result.is_finite():
            raise CorrectionError(f"Non-finite score {field} in record {number}")
        return result

    for field in touched:
        value = changes[field]
        if value == "":  # Withholding a result is distinct from assigning zero.
            continue
        parsed = numeric(field, value)
        if field in SCORE_TIERS and parsed not in SCORE_TIERS[field]:
            raise CorrectionError(f"Score outside rubric tiers: {field}={value}, record {number}")
    after = row | changes
    if set(SCORE_TIERS) | {"final_score"} <= set(after):
        values = {field: after[field] for field in SCORE_TIERS}
        if after["final_score"] == "":
            return  # Preserve components as provenance if the total is withheld.
        if any(value == "" for value in values.values()):
            raise CorrectionError(f"Scored record has missing components: {number}")
        total = sum(numeric(field, value) for field, value in values.items())
        if numeric("final_score", after["final_score"]) != total:
            raise CorrectionError(f"Score total does not equal its components: {number}")
    if after.get("trust_score", "") and after.get("final_score", ""):
        if numeric("trust_score", after["trust_score"]) != numeric("final_score", after["final_score"]):
            raise CorrectionError(f"Linked score fields disagree: {number}")


def validate_probability_change(row: dict, changes: dict, number: int) -> None:
    """Reject contradictory new probability assignments, not historical rows.

    ClockBase's fdr field records BH-adjusted P values. A row check cannot
    establish complete-family multiplicity control or validate a raw P value.
    Withholding is an empty field, never an invented neutral measurement.
    """
    fields = {"pvalue", "p_value", "fdr"}
    if not fields.intersection(changes):
        return
    after = row | changes

    def probability(field, value):
        if value == "":
            return None
        try:
            parsed = Decimal(value)
        except (InvalidOperation, ValueError) as exc:
            raise CorrectionError(f"Non-numeric {field} in record {number}") from exc
        if not parsed.is_finite() or not 0 <= parsed <= 1:
            raise CorrectionError(f"Invalid probability {field} in record {number}")
        return parsed

    for field in fields.intersection(changes):
        probability(field, changes[field])
    dependencies = {"pvalue": ("neg_log10_pvalue",), "p_value": ("neg_log10_pvalue",),
                    "fdr": ("neg_log10_fdr", "is_significant")}
    for field in fields.intersection(changes):
        new_value = probability(field, after[field])
        try:
            old_value = probability(field, row[field])
        except CorrectionError:
            old_value = None
        if new_value == old_value and new_value is not None:
            continue
        for derived in dependencies[field]:
            if not after.get(derived, ""):
                continue
            if new_value is None or derived not in changes:
                raise CorrectionError(f"Stale derived probability field {derived}: record {number}")
            if derived == "is_significant":
                expected = "true" if new_value < Decimal("0.05") else "false"
                if after[derived].lower() != expected:
                    raise CorrectionError(f"Derived threshold flag disagrees with FDR: record {number}")
            else:
                try:
                    observed = Decimal(after[derived])
                except (InvalidOperation, ValueError) as exc:
                    raise CorrectionError(f"Invalid derived probability field {derived}: record {number}") from exc
                if new_value == 0 or not observed.is_finite():
                    raise CorrectionError(f"Undefined or nonfinite {derived}: record {number}")
                expected = -new_value.log10()
                allowed = max(Decimal("1e-12"), abs(expected) * Decimal("1e-10"))
                if abs(observed - expected) > allowed:
                    raise CorrectionError(f"Derived logarithm disagrees with probability: record {number}")
    if "fdr" not in row:
        return
    q = probability("fdr", after["fdr"])
    for field in ("pvalue", "p_value"):
        if field not in row:
            continue
        # Removing an unsupported q does not require adjudicating an unchanged
        # raw P or other unresolved historical fields.
        if q is None:
            continue
        p = probability(field, after[field])
        if p is None:
            raise CorrectionError(f"Raw P withheld but inherited FDR retained: record {number}")
        if field in changes and "fdr" not in changes:
            try:
                old_p = probability(field, row[field])
            except CorrectionError:
                old_p = None
            if old_p != p:
                raise CorrectionError(f"Raw P changed without explicit FDR disposition: record {number}")
        if q < p and p - q > p * Decimal("1e-10"):
            raise CorrectionError(f"BH-adjusted P is below raw P: record {number}")


def patch_csv(raw: bytes, specification: dict) -> bytes:
    """Check every expected row before constructing a changed byte sequence."""
    if digest(raw) != specification["before_sha256"]:
        raise CorrectionError("Input hash does not match the adjudicated version")
    fields, rows, lines = csv_records(raw)
    if len(rows) != specification["row_count"]:
        raise CorrectionError("CSV row count changed")
    operations = specification["operations"]
    seen = set()
    replacements = []
    for operation in operations:
        number = operation["record_number"]
        if number in seen or number < 2 or number > len(rows) + 1:
            raise CorrectionError(f"Duplicate or invalid target record: {number}")
        seen.add(number)
        row, start, end = rows[number - 2]
        if row != operation["before"]:
            raise CorrectionError(f"Record identity or expected fields changed: {number}")
        changes = operation["changes"]
        if not changes or not set(changes) <= set(fields):
            raise CorrectionError(f"Unknown or empty correction fields: {number}")
        if set(changes) & IDENTITY_FIELDS:
            raise CorrectionError(f"Identity changes require separate adjudication: {number}")
        if any(not isinstance(value, str) for value in changes.values()):
            raise CorrectionError("CSV corrections must preserve explicit string values")
        validate_probability_change(row, changes, number)
        validate_score_change(row, changes, number)
        after = row | changes
        tail = lines[end - 1]
        newline = "\r\n" if tail.endswith("\r\n") else "\n"
        buffer = io.StringIO(newline="")
        writer = csv.writer(buffer, lineterminator=newline)
        writer.writerow([after[field] for field in fields])
        replacement = buffer.getvalue()
        if not tail.endswith(("\n", "\r")):
            replacement = replacement.removesuffix(newline)
        replacements.append((start, end, replacement))
    for start, end, replacement in sorted(replacements, reverse=True):
        lines[start:end] = [replacement]
    return "".join(lines).encode("utf-8")


def prepare(root: Path, manifest: dict) -> list[dict]:
    """Validate the entire bundle before any file is written."""
    if manifest.get("schema_version") != "clockbase-record-corrections-v1":
        raise CorrectionError("Unsupported correction manifest")
    for evidence in manifest["evidence_files"]:
        path = local_path(root, evidence["path"])
        if inspect_file(path).sha256 != evidence["sha256"]:
            raise CorrectionError(f"Evidence hash changed: {path}")
    plans = []
    seen = set()
    for specification in manifest["files"]:
        output = local_path(root, specification["output"])
        if output in seen:
            raise CorrectionError(f"Duplicate output: {output}")
        seen.add(output)
        existing = read_bytes(output) if output.exists() else None
        if existing is not None and digest(existing) == specification["after_sha256"]:
            plans.append({"path": output, "bytes": existing, "existing": existing,
                          "state": "already_applied", "specification": specification})
            continue
        source = local_path(root, specification["source"])
        if output != source and existing is not None:
            raise CorrectionError(f"Refusing to replace a conflicting output: {output}")
        result = patch_csv(read_bytes(source), specification)
        if digest(result) != specification["after_sha256"]:
            raise CorrectionError(f"Planned output hash differs: {output}")
        plans.append({"path": output, "bytes": result, "existing": existing,
                      "state": "ready", "specification": specification})
    return plans


def atomic_write(path: Path, raw: bytes) -> None:
    temporary = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        mode = path.stat().st_mode & 0o777 if path.exists() else 0o644
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.",
                                         delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
        temporary.chmod(mode)
        os.replace(temporary, path)
        temporary = None
    except OSError as exc:
        raise CorrectionError(f"Cannot write corrected file {path}: {exc}") from exc
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def run(root: Path, manifest_path: Path, apply: bool = False) -> dict:
    manifest_raw = read_bytes(manifest_path)
    try:
        manifest = json.loads(manifest_raw)
    except ValueError as exc:
        raise CorrectionError(f"Invalid correction manifest: {exc}") from exc
    plans = prepare(root, manifest)
    report = {
        "schema_version": "clockbase-record-correction-application-v1",
        "manifest_sha256": digest(manifest_raw),
        "scope": manifest["scope"],
        "files": [{"path": str(plan["path"].relative_to(root.resolve())),
                   "state": plan["state"],
                   "rows": plan["specification"]["row_count"],
                   "changed_records": len(plan["specification"]["operations"]),
                   "before_sha256": plan["specification"]["before_sha256"],
                   "after_sha256": digest(plan["bytes"])}
                  for plan in plans],
        "scientific_acceptance": "CANNOT VERIFY: propagation does not adjudicate biology.",
    }
    receipt_path = local_path(root, manifest["application_receipt"])
    if receipt_path.exists():
        try:
            prior = json.loads(read_bytes(receipt_path))
        except ValueError as exc:
            raise CorrectionError(f"Invalid prior application receipt: {exc}") from exc
        if prior.get("manifest_sha256") != report["manifest_sha256"]:
            raise CorrectionError("Refusing to replace a receipt for another manifest")
    if apply:
        for plan in plans:
            current = read_bytes(plan["path"]) if plan["path"].exists() else None
            if current != plan["existing"]:
                raise CorrectionError(f"Output changed during preflight: {plan['path']}")
        for plan in plans:
            if plan["state"] == "ready":
                LOG.info("Applying %s", plan["path"])
                atomic_write(plan["path"], plan["bytes"])
            if inspect_file(plan["path"]).sha256 != digest(plan["bytes"]):
                raise CorrectionError(f"Post-write verification failed: {plan['path']}")
        if not receipt_path.exists():
            atomic_write(receipt_path, (json.dumps(report, indent=2) + "\n").encode())
    report["applied"] = apply
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    try:
        result = run(args.root.resolve(), args.manifest, args.apply)
    except (CorrectionError, OSError, KeyError, ValueError) as exc:
        LOG.error("Correction stopped: %s", exc)
        return 1
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
