"""Authority and independence rules for audit decisions."""

from __future__ import annotations

from dataclasses import dataclass

from .models import ActorRecord, ActorRole, ActorType, IdentityAssurance


_HARD_CHECKS = {
    "effect_recalculation",
    "file_identity",
    "signature_verification",
    "source_integrity",
    "sample_identity",
    "source_design",
    "measurement_applicability",
    "finding_confirmation",
}


@dataclass(frozen=True)
class PolicyResult:
    allowed: bool
    reason_code: str


class AuthorityPolicy:
    """Evaluate externally enforceable role and model separation rules."""

    def check_can_pass(
        self,
        actor: ActorRecord,
        predicate: str,
    ) -> PolicyResult:
        if predicate not in actor.authority_scope:
            return PolicyResult(False, "PREDICATE_OUTSIDE_AUTHORITY")
        if actor.actor_type is ActorType.MODEL and predicate in _HARD_CHECKS:
            return PolicyResult(False, "MODEL_CANNOT_AUTHORIZE_HARD_CHECK")
        if actor.role not in {
            ActorRole.DETERMINISTIC_CHECKER,
            ActorRole.HUMAN_ADJUDICATOR,
            ActorRole.EVIDENCE_CUSTODIAN,
        }:
            return PolicyResult(False, "ROLE_CANNOT_AUTHORIZE_CHECK")
        return PolicyResult(True, "AUTHORIZED")

    def verify_separation(
        self,
        producer: ActorRecord,
        repairer: ActorRecord | None,
        verifier: ActorRecord,
    ) -> PolicyResult:
        if repairer is not None and repairer.record_id == verifier.record_id:
            return PolicyResult(False, "SELF_VERIFICATION")
        model_actors = [
            actor
            for actor in (producer, verifier)
            if actor.actor_type is ActorType.MODEL
        ]
        if any(
            actor.identity_assurance is not IdentityAssurance.VERIFIED
            for actor in model_actors
        ):
            return PolicyResult(False, "UNVERIFIED_MODEL_IDENTITY")
        if (
            producer.actor_type is ActorType.MODEL
            and verifier.actor_type is ActorType.MODEL
            and producer.model_family == verifier.model_family
        ):
            return PolicyResult(False, "SAME_MODEL_FAMILY")
        return PolicyResult(True, "INDEPENDENT")
