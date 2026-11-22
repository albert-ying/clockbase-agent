"""Translate protected-repair manifests into immutable audit records."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from ..hashing import ArtifactIdentityError, inspect_file
from ..models import (
    ActorRecord,
    ActorRole,
    ActorType,
    ArtifactRecord,
    AuditRecord,
    CheckRecord,
    CheckStatus,
    EvidenceClassification,
    IdentityAssurance,
)


def _timestamp(path: Path) -> str:
    value = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc)
    return value.isoformat().replace("+00:00", "Z")


def _read_json(path: Path, label: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ArtifactIdentityError(f"Cannot parse {label} {path}: {exc}") from exc
    if not isinstance(value, dict) or not value.get("schema_version"):
        raise ArtifactIdentityError(f"{label} lacks schema_version: {path}")
    return value


def adapt_protected_repair(
    manifest_path: Path,
    summary_path: Path,
    release_id: str,
) -> tuple[AuditRecord, ...]:
    """Verify and translate a protected repair without approving its actors."""

    manifest_path = Path(manifest_path)
    summary_path = Path(summary_path)
    manifest_identity = inspect_file(manifest_path)
    summary_identity = inspect_file(summary_path)
    manifest = _read_json(manifest_path, "protected repair manifest")
    summary = _read_json(summary_path, "protected repair summary")
    created_at = _timestamp(manifest_path)
    suffix = manifest_identity.sha256[:16]
    adapter_id = f"actor:protected-repair-adapter:{suffix}"
    checker_id = f"actor:protected-repair-checker:{suffix}"
    manifest_id = f"artifact:protected-repair-manifest:{suffix}"
    summary_id = f"artifact:protected-repair-summary:{summary_identity.sha256[:16]}"
    records: list[AuditRecord] = [
        ActorRecord(
            record_id=adapter_id,
            release_id=release_id,
            created_at=created_at,
            actor_type=ActorType.PROGRAM,
            role=ActorRole.EVIDENCE_CUSTODIAN,
            organization="ClockBase independent audit",
            model_family=None,
            authority_scope=("artifact_registration",),
            identity_assurance=IdentityAssurance.VERIFIED,
        ),
        ActorRecord(
            record_id=checker_id,
            release_id=release_id,
            created_at=created_at,
            actor_type=ActorType.PROGRAM,
            role=ActorRole.DETERMINISTIC_CHECKER,
            organization="ClockBase protected repair verifier",
            model_family=None,
            authority_scope=("file_identity", "protected_repair_summary"),
            identity_assurance=IdentityAssurance.VERIFIED,
        ),
        ArtifactRecord(
            record_id=manifest_id,
            release_id=release_id,
            created_at=created_at,
            uri=str(manifest_identity.path),
            sha256=manifest_identity.sha256,
            byte_count=manifest_identity.byte_count,
            classification=EvidenceClassification.DOCUMENT,
            created_by_actor_id=adapter_id,
        ),
        ArtifactRecord(
            record_id=summary_id,
            release_id=release_id,
            created_at=created_at,
            uri=str(summary_identity.path),
            sha256=summary_identity.sha256,
            byte_count=summary_identity.byte_count,
            classification=EvidenceClassification.DOCUMENT,
            created_by_actor_id=adapter_id,
            parent_artifact_ids=(manifest_id,),
        ),
    ]
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, dict):
        raise ArtifactIdentityError("Protected repair manifest lacks artifacts")
    for name, declaration in sorted(artifacts.items()):
        if not isinstance(declaration, dict):
            raise ArtifactIdentityError(f"Invalid artifact declaration: {name}")
        artifact_path = manifest_path.parent / str(declaration.get("path", ""))
        actual = inspect_file(artifact_path)
        if actual.sha256 != declaration.get("sha256"):
            raise ArtifactIdentityError(
                f"Protected repair artifact digest mismatch: {artifact_path}"
            )
        if artifact_path.resolve() == summary_path.resolve() and (
            actual.sha256 != summary_identity.sha256
            or actual.byte_count != summary_identity.byte_count
        ):
            raise ArtifactIdentityError(
                f"Protected repair summary mismatch: {summary_path}"
            )
        records.append(
            ArtifactRecord(
                record_id=f"artifact:protected-repair:{name}:{actual.sha256[:16]}",
                release_id=release_id,
                created_at=created_at,
                uri=str(actual.path),
                sha256=actual.sha256,
                byte_count=actual.byte_count,
                classification=EvidenceClassification.DERIVED,
                created_by_actor_id=adapter_id,
                parent_artifact_ids=(manifest_id,),
            )
        )
    records.extend(
        (
            CheckRecord(
                record_id=f"check:protected-repair:summary:{suffix}",
                release_id=release_id,
                created_at=created_at,
                subject_id=summary_id,
                predicate="protected_repair_summary",
                status=CheckStatus.PASS,
                checker_actor_id=checker_id,
                evidence_artifact_ids=(manifest_id, summary_id),
                expected={"manifest_schema": manifest["schema_version"]},
                observed={"summary_schema": summary["schema_version"]},
            ),
            CheckRecord(
                record_id=f"check:protected-repair:repair-actor:{suffix}",
                release_id=release_id,
                created_at=created_at,
                subject_id=manifest_id,
                predicate="repair_actor_identity",
                status=CheckStatus.BLOCKED,
                checker_actor_id=checker_id,
                evidence_artifact_ids=(manifest_id,),
                expected={"identity_assurance": "VERIFIED"},
                observed={"identity_assurance": "UNVERIFIED"},
            ),
            CheckRecord(
                record_id=f"check:protected-repair:authority:{suffix}",
                release_id=release_id,
                created_at=created_at,
                subject_id=manifest_id,
                predicate="verification_authority",
                status=CheckStatus.BLOCKED,
                checker_actor_id=checker_id,
                evidence_artifact_ids=(manifest_id,),
                expected={"authority": "INDEPENDENT_VERIFIED"},
                observed={"authority": "UNRESOLVED"},
            ),
        )
    )
    return tuple(records)
