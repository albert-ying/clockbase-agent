"""Fail-closed evaluation of evidence dependencies for one claim."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote, urlparse

from .dependency_graph import DependencyGraph
from .hashing import ArtifactIdentityError, inspect_file, sha256_value
from .ledger import LedgerSnapshot
from .models import (
    ActorRecord,
    ActorRole,
    ArtifactRecord,
    CheckRecord,
    CheckStatus,
    ClaimDisposition,
    ClaimRecord,
    DecisionRecord,
    DecisionStatus,
    DependencyRecord,
    DependencyRequirement,
    EvidenceClassification,
    FindingRecord,
    FindingStatus,
    Severity,
)
from .policy import AuthorityPolicy


GATE_VERSION = "clockbase_claim_gate_v2"


class ClaimGateError(RuntimeError):
    """Raised when the requested claim cannot be evaluated as a record."""


def _artifact_path(uri: str) -> Path | None:
    parsed = urlparse(uri)
    if parsed.scheme == "file":
        return Path(unquote(parsed.path))
    if parsed.scheme:
        return None
    return Path(uri)


def evaluate_claim(
    snapshot: LedgerSnapshot,
    claim_id: str,
    policy: AuthorityPolicy,
) -> DecisionRecord:
    """Evaluate all declared requirements and retain every decision reason."""

    claim = snapshot.records.get(claim_id)
    if not isinstance(claim, ClaimRecord):
        raise ClaimGateError(f"Claim record not found: {claim_id}")

    rejection: set[str] = set()
    blocking: set[str] = set()
    evaluated: set[str] = {claim_id}

    if claim.disposition is ClaimDisposition.WITHDRAWN:
        return _decision(
            claim,
            snapshot,
            DecisionStatus.WITHDRAWN,
            (),
            evaluated,
        )

    for dependency_id in claim.required_dependency_ids:
        record = snapshot.records.get(dependency_id)
        if not isinstance(record, DependencyRecord):
            blocking.add("MISSING_REQUIRED_DEPENDENCY")
        else:
            evaluated.add(dependency_id)

    graph = DependencyGraph.from_snapshot(snapshot, claim.release_id)
    dependencies = graph.required_ancestors(claim_id)
    if not dependencies:
        blocking.add("NO_REQUIRED_EVIDENCE")
    connected = {record.record_id for record in dependencies}
    for dependency_id in claim.required_dependency_ids:
        if dependency_id in snapshot.records and dependency_id not in connected:
            blocking.add("DISCONNECTED_REQUIRED_DEPENDENCY")
    relevant_subjects = {claim_id}
    for dependency in dependencies:
        evaluated.add(dependency.record_id)
        evaluated.add(dependency.parent_id)
        relevant_subjects.add(dependency.parent_id)
        parent = snapshot.records[dependency.parent_id]
        if isinstance(parent, ArtifactRecord):
            if dependency.requirement is DependencyRequirement.REQUIRED_CHECK:
                blocking.add("DEPENDENCY_TYPE_MISMATCH")
                continue
            _evaluate_artifact(
                parent,
                dependency.requirement,
                rejection,
                blocking,
            )
        elif isinstance(parent, CheckRecord):
            if dependency.requirement is not DependencyRequirement.REQUIRED_CHECK:
                blocking.add("DEPENDENCY_TYPE_MISMATCH")
                continue
            _evaluate_check(
                parent,
                snapshot,
                policy,
                rejection,
                blocking,
                evaluated,
            )
        else:
            blocking.add("DEPENDENCY_TYPE_MISMATCH")

    for record in snapshot.records.values():
        if record.release_id != claim.release_id:
            continue
        if not isinstance(record, FindingRecord):
            continue
        affected = {record.subject_id, *record.affected_subject_ids}
        if not (affected & relevant_subjects):
            continue
        evaluated.add(record.record_id)
        if record.status is FindingStatus.CONFIRMED and record.severity in {
            Severity.S0,
            Severity.S1,
        }:
            author = snapshot.records.get(record.author_actor_id)
            finding_blocking: set[str] = set()
            if (not isinstance(author, ActorRecord) or author.release_id != record.release_id or
                    not policy.check_can_pass(author, "finding_confirmation").allowed):
                finding_blocking.add("UNAUTHORIZED_FINDING_CONFIRMATION")
            else:
                evaluated.add(author.record_id)
            if not record.evidence_artifact_ids:
                finding_blocking.add("MISSING_FINDING_EVIDENCE")
            for evidence_id in record.evidence_artifact_ids:
                evidence = snapshot.records.get(evidence_id)
                if not isinstance(evidence, ArtifactRecord) or evidence.release_id != record.release_id:
                    finding_blocking.add("MISSING_FINDING_EVIDENCE")
                    continue
                evaluated.add(evidence_id)
                _evaluate_artifact(evidence, DependencyRequirement.SUPPORTING, rejection, finding_blocking)
            blocking.update(finding_blocking)
            if not finding_blocking:
                rejection.add("CONFIRMED_MATERIAL_FINDING")
        elif record.status is FindingStatus.PROPOSED and record.severity in {
            Severity.S0,
            Severity.S1,
            Severity.S2,
        }:
            blocking.add("UNRESOLVED_MATERIAL_FINDING")

    if rejection:
        status = DecisionStatus.REJECTED
        reasons = rejection | blocking
    elif blocking:
        status = DecisionStatus.BLOCKED
        reasons = blocking
    else:
        status = DecisionStatus.APPROVED
        reasons = set()
    return _decision(claim, snapshot, status, reasons, evaluated)


def _evaluate_artifact(
    artifact: ArtifactRecord,
    requirement: DependencyRequirement,
    rejection: set[str],
    blocking: set[str],
) -> None:
    if artifact.quarantined:
        rejection.add("QUARANTINED_ANCESTOR")
    if (
        requirement is DependencyRequirement.REQUIRED_EMPIRICAL
        and artifact.classification is EvidenceClassification.SIMULATION
    ):
        rejection.add("SIMULATION_AS_EMPIRICAL")
    elif (
        requirement is DependencyRequirement.REQUIRED_EMPIRICAL
        and artifact.classification is not EvidenceClassification.EMPIRICAL
    ):
        blocking.add("EMPIRICAL_EVIDENCE_REQUIRED")
    elif (
        requirement is DependencyRequirement.REQUIRED_DERIVED
        and artifact.classification is not EvidenceClassification.DERIVED
    ):
        blocking.add("DERIVED_EVIDENCE_REQUIRED")

    path = _artifact_path(artifact.uri)
    if path is None:
        if artifact.classification is not EvidenceClassification.SOFTWARE_FIXTURE:
            blocking.add("UNVERIFIABLE_ARTIFACT_URI")
        return
    try:
        identity = inspect_file(path)
    except ArtifactIdentityError:
        blocking.add("ARTIFACT_UNAVAILABLE")
        return
    if identity.sha256 != artifact.sha256 or identity.byte_count != artifact.byte_count:
        blocking.add("SOURCE_IDENTITY_DRIFT")


def _evaluate_check(
    check: CheckRecord,
    snapshot: LedgerSnapshot,
    policy: AuthorityPolicy,
    rejection: set[str],
    blocking: set[str],
    evaluated: set[str],
) -> None:
    if check.status in {
        CheckStatus.BLOCKED,
        CheckStatus.NEEDS_REVIEW,
        CheckStatus.NOT_APPLICABLE,
    }:
        blocking.add(f"REQUIRED_CHECK_{check.status.value}")
        return
    check_blocking: set[str] = set()
    checker = snapshot.records.get(check.checker_actor_id)
    if not isinstance(checker, ActorRecord) or checker.release_id != check.release_id:
        check_blocking.add("MISSING_CHECKER_ACTOR")
    else:
        evaluated.add(checker.record_id)
        authority = policy.check_can_pass(checker, check.predicate)
        if not authority.allowed:
            check_blocking.add(f"UNAUTHORIZED_{check.status.value}")

    subject = snapshot.records.get(check.subject_id)
    if subject is None or subject.release_id != check.release_id:
        check_blocking.add("MISSING_CHECK_SUBJECT")
    else:
        evaluated.add(subject.record_id)
    if not check.evidence_artifact_ids:
        check_blocking.add("MISSING_CHECK_EVIDENCE")
    # A cited source does not need a second dependency edge to retain its
    # integrity requirements. Check every cited artifact before trusting a verdict.
    for evidence_id in check.evidence_artifact_ids:
        evidence = snapshot.records.get(evidence_id)
        if not isinstance(evidence, ArtifactRecord) or evidence.release_id != check.release_id:
            check_blocking.add("MISSING_CHECK_EVIDENCE")
            continue
        evaluated.add(evidence_id)
        _evaluate_artifact(evidence, DependencyRequirement.SUPPORTING, rejection, check_blocking)
    if isinstance(subject, ArtifactRecord):
        _evaluate_artifact(subject, DependencyRequirement.SUPPORTING, rejection, check_blocking)
        creator = snapshot.records.get(subject.created_by_actor_id)
        if isinstance(creator, ActorRecord) and isinstance(checker, ActorRecord):
            evaluated.add(creator.record_id)
            if (
                creator.role is ActorRole.REPAIRER
                and creator.record_id == checker.record_id
            ):
                check_blocking.add("SELF_VERIFICATION")
    blocking.update(check_blocking)
    if check.status is CheckStatus.FAIL and not check_blocking:
        rejection.add("REQUIRED_CHECK_FAILED")


def _decision(
    claim: ClaimRecord,
    snapshot: LedgerSnapshot,
    status: DecisionStatus,
    reasons: set[str] | tuple[str, ...],
    evaluated: set[str],
) -> DecisionRecord:
    basis = {
        "claim_id": claim.record_id,
        "ledger_sha256": snapshot.sha256,
        "gate_version": GATE_VERSION,
    }
    timestamp = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    return DecisionRecord(
        record_id=f"decision:{sha256_value(basis)[:24]}",
        release_id=claim.release_id,
        created_at=timestamp,
        claim_id=claim.record_id,
        status=status,
        reason_codes=tuple(sorted(reasons)),
        evaluated_record_ids=tuple(sorted(evaluated)),
        gate_version=GATE_VERSION,
        ledger_sha256=snapshot.sha256,
    )
