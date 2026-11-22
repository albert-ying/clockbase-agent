"""Immutable record types for the ClockBase audit ledger."""

from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
from pathlib import Path
from typing import Any, ClassVar, Mapping


class RecordValidationError(ValueError):
    """Raised when a record violates the version 1 exchange contract."""


class StringEnum(str, Enum):
    """An enum that serializes to its stable string value."""


class EvidenceClassification(StringEnum):
    EMPIRICAL = "EMPIRICAL"
    DERIVED = "DERIVED"
    SIMULATION = "SIMULATION"
    DOCUMENT = "DOCUMENT"
    SOFTWARE_FIXTURE = "SOFTWARE_FIXTURE"


class ActorType(StringEnum):
    HUMAN = "HUMAN"
    PROGRAM = "PROGRAM"
    MODEL = "MODEL"


class ActorRole(StringEnum):
    EVIDENCE_CUSTODIAN = "EVIDENCE_CUSTODIAN"
    PRIMARY_ANALYST = "PRIMARY_ANALYST"
    REPAIRER = "REPAIRER"
    AUDITOR = "AUDITOR"
    DETERMINISTIC_CHECKER = "DETERMINISTIC_CHECKER"
    HUMAN_ADJUDICATOR = "HUMAN_ADJUDICATOR"
    GATE = "GATE"


class IdentityAssurance(StringEnum):
    UNVERIFIED = "UNVERIFIED"
    VERIFIED = "VERIFIED"


class CheckStatus(StringEnum):
    PASS = "PASS"
    FAIL = "FAIL"
    BLOCKED = "BLOCKED"
    NEEDS_REVIEW = "NEEDS_REVIEW"
    NOT_APPLICABLE = "NOT_APPLICABLE"


class FindingStatus(StringEnum):
    PROPOSED = "PROPOSED"
    CONFIRMED = "CONFIRMED"
    REJECTED = "REJECTED"
    RESOLVED = "RESOLVED"


class Severity(StringEnum):
    S0 = "S0"
    S1 = "S1"
    S2 = "S2"
    S3 = "S3"


class DependencyRequirement(StringEnum):
    REQUIRED_EMPIRICAL = "REQUIRED_EMPIRICAL"
    REQUIRED_DERIVED = "REQUIRED_DERIVED"
    REQUIRED_CHECK = "REQUIRED_CHECK"
    SUPPORTING = "SUPPORTING"
    CONTRADICTING = "CONTRADICTING"


class ClaimDisposition(StringEnum):
    PROPOSED = "PROPOSED"
    WITHDRAWN = "WITHDRAWN"


class DecisionStatus(StringEnum):
    APPROVED = "APPROVED"
    REJECTED = "REJECTED"
    BLOCKED = "BLOCKED"
    WITHDRAWN = "WITHDRAWN"


def _require_text(name: str, value: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise RecordValidationError(f"{name} must be a nonempty string")


@dataclass(frozen=True)
class AuditRecord:
    record_id: str
    release_id: str
    created_at: str

    record_type: ClassVar[str] = "audit_record"

    def __post_init__(self) -> None:
        _require_text("record_id", self.record_id)
        _require_text("release_id", self.release_id)
        _require_text("created_at", self.created_at)


@dataclass(frozen=True)
class ArtifactRecord(AuditRecord):
    uri: str
    sha256: str
    byte_count: int
    classification: EvidenceClassification
    created_by_actor_id: str
    parent_artifact_ids: tuple[str, ...] = ()
    quarantined: bool = False

    record_type: ClassVar[str] = "artifact"

    def __post_init__(self) -> None:
        super().__post_init__()
        _require_text("uri", self.uri)
        _require_text("sha256", self.sha256)
        _require_text("created_by_actor_id", self.created_by_actor_id)
        if self.byte_count < 0:
            raise RecordValidationError("byte_count must be nonnegative")


@dataclass(frozen=True)
class ActorRecord(AuditRecord):
    actor_type: ActorType
    role: ActorRole
    organization: str
    model_family: str | None
    authority_scope: tuple[str, ...]
    identity_assurance: IdentityAssurance = IdentityAssurance.UNVERIFIED

    record_type: ClassVar[str] = "actor"

    def __post_init__(self) -> None:
        super().__post_init__()
        _require_text("organization", self.organization)
        if self.actor_type is ActorType.MODEL and not self.model_family:
            raise RecordValidationError("model_family is required for model actors")


@dataclass(frozen=True)
class CheckRecord(AuditRecord):
    subject_id: str
    predicate: str
    status: CheckStatus
    checker_actor_id: str
    evidence_artifact_ids: tuple[str, ...]
    expected: Mapping[str, Any]
    observed: Mapping[str, Any]

    record_type: ClassVar[str] = "check"

    def __post_init__(self) -> None:
        super().__post_init__()
        _require_text("subject_id", self.subject_id)
        _require_text("predicate", self.predicate)
        _require_text("checker_actor_id", self.checker_actor_id)


@dataclass(frozen=True)
class FindingRecord(AuditRecord):
    subject_id: str
    severity: Severity
    status: FindingStatus
    author_actor_id: str
    evidence_artifact_ids: tuple[str, ...]
    affected_subject_ids: tuple[str, ...]
    allegation: str
    requested_test: str

    record_type: ClassVar[str] = "finding"

    def __post_init__(self) -> None:
        super().__post_init__()
        _require_text("subject_id", self.subject_id)
        _require_text("author_actor_id", self.author_actor_id)
        _require_text("allegation", self.allegation)
        _require_text("requested_test", self.requested_test)


@dataclass(frozen=True)
class DependencyRecord(AuditRecord):
    parent_id: str
    child_id: str
    requirement: DependencyRequirement

    record_type: ClassVar[str] = "dependency"

    def __post_init__(self) -> None:
        super().__post_init__()
        _require_text("parent_id", self.parent_id)
        _require_text("child_id", self.child_id)


@dataclass(frozen=True)
class ClaimRecord(AuditRecord):
    exact_text: str
    claim_type: str
    scope: str
    required_dependency_ids: tuple[str, ...]
    falsifier: str
    disposition: ClaimDisposition = ClaimDisposition.PROPOSED

    record_type: ClassVar[str] = "claim"

    def __post_init__(self) -> None:
        super().__post_init__()
        _require_text("exact_text", self.exact_text)
        _require_text("claim_type", self.claim_type)
        _require_text("scope", self.scope)
        _require_text("falsifier", self.falsifier)


@dataclass(frozen=True)
class DecisionRecord(AuditRecord):
    claim_id: str
    status: DecisionStatus
    reason_codes: tuple[str, ...]
    evaluated_record_ids: tuple[str, ...]
    gate_version: str
    ledger_sha256: str

    record_type: ClassVar[str] = "decision"

    def __post_init__(self) -> None:
        super().__post_init__()
        _require_text("claim_id", self.claim_id)
        _require_text("gate_version", self.gate_version)
        _require_text("ledger_sha256", self.ledger_sha256)


Record = (
    ArtifactRecord
    | ActorRecord
    | CheckRecord
    | FindingRecord
    | DependencyRecord
    | ClaimRecord
    | DecisionRecord
)


def _json_value(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return {field.name: _json_value(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def record_to_dict(record: AuditRecord) -> dict[str, Any]:
    """Return the stable JSON representation of one record."""

    if not isinstance(record, AuditRecord):
        raise RecordValidationError("record must inherit AuditRecord")
    payload = _json_value(record)
    payload["record_type"] = record.record_type
    return payload


def record_from_dict(payload: Mapping[str, Any]) -> Record:
    """Parse one JSON record and enforce its typed constructor contract."""

    if not isinstance(payload, Mapping):
        raise RecordValidationError("record payload must be an object")
    values = dict(payload)
    record_type = values.pop("record_type", None)
    try:
        if record_type == "artifact":
            values["classification"] = EvidenceClassification(
                values["classification"]
            )
            values["parent_artifact_ids"] = tuple(
                values.get("parent_artifact_ids", ())
            )
            return ArtifactRecord(**values)
        if record_type == "actor":
            values["actor_type"] = ActorType(values["actor_type"])
            values["role"] = ActorRole(values["role"])
            values["authority_scope"] = tuple(values.get("authority_scope", ()))
            values["identity_assurance"] = IdentityAssurance(
                values.get("identity_assurance", IdentityAssurance.UNVERIFIED)
            )
            return ActorRecord(**values)
        if record_type == "check":
            values["status"] = CheckStatus(values["status"])
            values["evidence_artifact_ids"] = tuple(
                values.get("evidence_artifact_ids", ())
            )
            return CheckRecord(**values)
        if record_type == "finding":
            values["severity"] = Severity(values["severity"])
            values["status"] = FindingStatus(values["status"])
            values["evidence_artifact_ids"] = tuple(
                values.get("evidence_artifact_ids", ())
            )
            values["affected_subject_ids"] = tuple(
                values.get("affected_subject_ids", ())
            )
            return FindingRecord(**values)
        if record_type == "dependency":
            values["requirement"] = DependencyRequirement(values["requirement"])
            return DependencyRecord(**values)
        if record_type == "claim":
            values["required_dependency_ids"] = tuple(
                values.get("required_dependency_ids", ())
            )
            values["disposition"] = ClaimDisposition(
                values.get("disposition", ClaimDisposition.PROPOSED)
            )
            return ClaimRecord(**values)
        if record_type == "decision":
            values["status"] = DecisionStatus(values["status"])
            values["reason_codes"] = tuple(values.get("reason_codes", ()))
            values["evaluated_record_ids"] = tuple(
                values.get("evaluated_record_ids", ())
            )
            return DecisionRecord(**values)
    except (KeyError, TypeError, ValueError) as exc:
        raise RecordValidationError(
            f"Invalid {record_type or 'unknown'} record: {exc}"
        ) from exc
    raise RecordValidationError(f"Unsupported record_type: {record_type!r}")
