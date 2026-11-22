import tempfile
import unittest
from pathlib import Path

from audit_system.claim_gate import evaluate_claim
from audit_system.hashing import inspect_file
from audit_system.ledger import Ledger
from audit_system.models import (
    ActorRecord,
    ActorRole,
    ActorType,
    ArtifactRecord,
    CheckRecord,
    CheckStatus,
    ClaimDisposition,
    ClaimRecord,
    DecisionStatus,
    DependencyRecord,
    DependencyRequirement,
    EvidenceClassification,
    FindingRecord,
    FindingStatus,
    IdentityAssurance,
    Severity,
)
from audit_system.policy import AuthorityPolicy


class ClaimGateTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)
        self.ledger = Ledger(self.directory / "records.jsonl")
        self.release_id = "release:test"
        self.created_at = "2026-09-07T00:00:00Z"
        self.add_actor(
            "actor:custodian",
            ActorType.PROGRAM,
            ActorRole.EVIDENCE_CUSTODIAN,
            ("artifact_registration",),
        )
        self.add_actor(
            "actor:checker",
            ActorType.PROGRAM,
            ActorRole.DETERMINISTIC_CHECKER,
            ("file_identity", "effect_recalculation"),
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def add_actor(self, record_id, actor_type, role, scope, family=None):
        self.ledger.append(
            ActorRecord(
                record_id=record_id,
                release_id=self.release_id,
                created_at=self.created_at,
                actor_type=actor_type,
                role=role,
                organization="ClockBase audit",
                model_family=family,
                authority_scope=scope,
                identity_assurance=IdentityAssurance.VERIFIED,
            )
        )

    def add_artifact(
        self,
        record_id="artifact:source",
        classification=EvidenceClassification.EMPIRICAL,
        quarantined=False,
        path=None,
        creator="actor:custodian",
    ):
        if path is None:
            uri = f"fixture://{record_id}"
            digest = "0" * 64
            byte_count = 1
        else:
            identity = inspect_file(path)
            uri = str(identity.path)
            digest = identity.sha256
            byte_count = identity.byte_count
        self.ledger.append(
            ArtifactRecord(
                record_id=record_id,
                release_id=self.release_id,
                created_at=self.created_at,
                uri=uri,
                sha256=digest,
                byte_count=byte_count,
                classification=classification,
                created_by_actor_id=creator,
                quarantined=quarantined,
            )
        )

    def add_check(self, status=CheckStatus.PASS, checker="actor:checker"):
        self.ledger.append(
            CheckRecord(
                record_id="check:effect",
                release_id=self.release_id,
                created_at=self.created_at,
                subject_id="artifact:source",
                predicate="effect_recalculation",
                status=status,
                checker_actor_id=checker,
                evidence_artifact_ids=("artifact:source",),
                expected={"status": "valid"},
                observed={"status": status.value},
            )
        )

    def add_claim(self, disposition=ClaimDisposition.PROPOSED, required_ids=()):
        self.ledger.append(
            ClaimRecord(
                record_id="claim:central",
                release_id=self.release_id,
                created_at=self.created_at,
                exact_text="The central ClockBase result is supported.",
                claim_type="scientific",
                scope="submitted manuscript",
                required_dependency_ids=required_ids,
                falsifier="The source-linked calculation fails.",
                disposition=disposition,
            )
        )

    def add_dependency(self, record_id, parent, requirement):
        self.ledger.append(
            DependencyRecord(
                record_id=record_id,
                release_id=self.release_id,
                created_at=self.created_at,
                parent_id=parent,
                child_id="claim:central",
                requirement=requirement,
            )
        )

    def evaluate(self):
        return evaluate_claim(
            self.ledger.snapshot(),
            "claim:central",
            AuthorityPolicy(),
        )

    def test_complete_authorized_evidence_can_be_approved(self):
        source = self.directory / "approved-source.csv"
        source.write_text("sample,value\nS1,1\n", encoding="utf-8")
        self.add_artifact(path=source)
        self.add_check()
        self.add_claim()
        self.add_dependency(
            "dependency:source",
            "artifact:source",
            DependencyRequirement.REQUIRED_EMPIRICAL,
        )
        self.add_dependency(
            "dependency:check",
            "check:effect",
            DependencyRequirement.REQUIRED_CHECK,
        )

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.APPROVED)
        self.assertEqual(decision.reason_codes, ())

    def test_claim_without_required_evidence_cannot_be_approved(self):
        self.add_claim()
        decision = self.evaluate()
        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("NO_REQUIRED_EVIDENCE", decision.reason_codes)

    def test_declared_requirement_must_connect_to_the_claim(self):
        source = self.directory / "source.csv"
        source.write_text("software fixture\n", encoding="utf-8")
        self.add_artifact(path=source)
        self.add_check()
        self.add_claim(required_ids=("dependency:disconnected",))
        self.ledger.append(DependencyRecord(
            record_id="dependency:disconnected", release_id=self.release_id,
            created_at=self.created_at, parent_id="artifact:source",
            child_id="check:effect", requirement=DependencyRequirement.REQUIRED_EMPIRICAL))
        decision = self.evaluate()
        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("DISCONNECTED_REQUIRED_DEPENDENCY", decision.reason_codes)

    def test_pass_check_cannot_hide_changed_cited_evidence(self):
        source = self.directory / "source.csv"
        source.write_text("software fixture\n", encoding="utf-8")
        self.add_artifact(path=source)
        self.add_check()
        self.add_claim()
        self.add_dependency("dependency:check", "check:effect", DependencyRequirement.REQUIRED_CHECK)
        source.write_text("changed software fixture\n", encoding="utf-8")
        decision = self.evaluate()
        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("SOURCE_IDENTITY_DRIFT", decision.reason_codes)

    def test_pass_check_requires_available_subject_and_evidence(self):
        self.add_check()
        self.add_claim()
        self.add_dependency("dependency:check", "check:effect", DependencyRequirement.REQUIRED_CHECK)
        decision = self.evaluate()
        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("MISSING_CHECK_SUBJECT", decision.reason_codes)
        self.assertIn("MISSING_CHECK_EVIDENCE", decision.reason_codes)

    def test_artifact_cannot_satisfy_a_required_check(self):
        source = self.directory / "source.csv"
        source.write_text("software fixture\n", encoding="utf-8")
        self.add_artifact(path=source)
        self.add_claim()
        self.add_dependency("dependency:wrong-type", "artifact:source", DependencyRequirement.REQUIRED_CHECK)
        decision = self.evaluate()
        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("DEPENDENCY_TYPE_MISMATCH", decision.reason_codes)

    def test_missing_required_dependency_blocks_claim(self):
        self.add_claim(required_ids=("dependency:missing",))

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("MISSING_REQUIRED_DEPENDENCY", decision.reason_codes)

    def test_failed_required_check_rejects_claim(self):
        source = self.directory / "failed-check-source.csv"
        source.write_text("software fixture\n", encoding="utf-8")
        self.add_artifact(path=source)
        self.add_check(CheckStatus.FAIL)
        self.add_claim()
        self.add_dependency(
            "dependency:check",
            "check:effect",
            DependencyRequirement.REQUIRED_CHECK,
        )

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.REJECTED)
        self.assertIn("REQUIRED_CHECK_FAILED", decision.reason_codes)

    def test_model_failure_does_not_become_confirmed_rejection(self):
        source = self.directory / "model-failure-source.csv"
        source.write_text("software fixture\n", encoding="utf-8")
        self.add_actor("actor:model", ActorType.MODEL, ActorRole.AUDITOR,
                       ("effect_recalculation",), family="gpt")
        self.add_artifact(path=source)
        self.add_check(CheckStatus.FAIL, checker="actor:model")
        self.add_claim()
        self.add_dependency("dependency:check", "check:effect", DependencyRequirement.REQUIRED_CHECK)
        decision = self.evaluate()
        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("UNAUTHORIZED_FAIL", decision.reason_codes)

    def test_simulation_cannot_satisfy_empirical_requirement(self):
        self.add_artifact(classification=EvidenceClassification.SIMULATION)
        self.add_claim()
        self.add_dependency(
            "dependency:source",
            "artifact:source",
            DependencyRequirement.REQUIRED_EMPIRICAL,
        )

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.REJECTED)
        self.assertIn("SIMULATION_AS_EMPIRICAL", decision.reason_codes)

    def test_quarantined_required_artifact_rejects_claim(self):
        self.add_artifact(quarantined=True)
        self.add_claim()
        self.add_dependency(
            "dependency:source",
            "artifact:source",
            DependencyRequirement.REQUIRED_EMPIRICAL,
        )

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.REJECTED)
        self.assertIn("QUARANTINED_ANCESTOR", decision.reason_codes)

    def test_proposed_material_finding_blocks_claim(self):
        self.add_claim()
        self.ledger.append(
            FindingRecord(
                record_id="finding:one",
                release_id=self.release_id,
                created_at=self.created_at,
                subject_id="claim:central",
                severity=Severity.S1,
                status=FindingStatus.PROPOSED,
                author_actor_id="actor:checker",
                evidence_artifact_ids=(),
                affected_subject_ids=("claim:central",),
                allegation="The source assignment may be incorrect.",
                requested_test="Reconstruct sample assignments.",
            )
        )

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("UNRESOLVED_MATERIAL_FINDING", decision.reason_codes)

    def test_model_cannot_promote_its_finding_to_confirmed(self):
        self.add_actor("actor:model", ActorType.MODEL, ActorRole.AUDITOR,
                       ("finding_confirmation",), family="gpt")
        self.add_claim()
        self.ledger.append(FindingRecord(
            record_id="finding:unsupported-confirmation", release_id=self.release_id,
            created_at=self.created_at, subject_id="claim:central", severity=Severity.S1,
            status=FindingStatus.CONFIRMED, author_actor_id="actor:model",
            evidence_artifact_ids=(), affected_subject_ids=(),
            allegation="Model-declared error", requested_test="Verify from source evidence"))
        decision=self.evaluate()
        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("UNAUTHORIZED_FINDING_CONFIRMATION", decision.reason_codes)
        self.assertIn("MISSING_FINDING_EVIDENCE", decision.reason_codes)

    def test_authorized_source_supported_finding_can_reject(self):
        source=self.directory / "finding-source.txt"
        source.write_text("software fixture\n", encoding="utf-8")
        self.add_artifact(path=source)
        self.add_actor("actor:adjudicator", ActorType.HUMAN, ActorRole.HUMAN_ADJUDICATOR,
                       ("finding_confirmation",))
        self.add_claim()
        self.add_dependency("dependency:source", "artifact:source", DependencyRequirement.REQUIRED_EMPIRICAL)
        self.ledger.append(FindingRecord(
            record_id="finding:source-supported", release_id=self.release_id,
            created_at=self.created_at, subject_id="claim:central", severity=Severity.S1,
            status=FindingStatus.CONFIRMED, author_actor_id="actor:adjudicator",
            evidence_artifact_ids=("artifact:source",), affected_subject_ids=(),
            allegation="Adjudicated fixture error", requested_test="Source review completed in fixture"))
        decision=self.evaluate()
        self.assertEqual(decision.status, DecisionStatus.REJECTED)
        self.assertIn("CONFIRMED_MATERIAL_FINDING", decision.reason_codes)

    def test_changed_source_bytes_block_prior_identity(self):
        path = self.directory / "source.csv"
        path.write_text("sample,value\nS1,1\n", encoding="utf-8")
        self.add_artifact(path=path)
        self.add_claim()
        self.add_dependency(
            "dependency:source",
            "artifact:source",
            DependencyRequirement.REQUIRED_EMPIRICAL,
        )
        path.write_text("sample,value\nS1,2\n", encoding="utf-8")

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("SOURCE_IDENTITY_DRIFT", decision.reason_codes)

    def test_model_pass_is_advisory_and_blocks_claim(self):
        self.add_actor(
            "actor:model",
            ActorType.MODEL,
            ActorRole.AUDITOR,
            ("effect_recalculation",),
            family="gpt",
        )
        self.add_artifact()
        self.add_check(checker="actor:model")
        self.add_claim()
        self.add_dependency(
            "dependency:check",
            "check:effect",
            DependencyRequirement.REQUIRED_CHECK,
        )

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("UNAUTHORIZED_PASS", decision.reason_codes)

    def test_repairer_checking_own_artifact_blocks_claim(self):
        self.add_actor(
            "actor:repairer",
            ActorType.HUMAN,
            ActorRole.REPAIRER,
            ("effect_recalculation",),
        )
        self.add_artifact(creator="actor:repairer")
        self.add_check(checker="actor:repairer")
        self.add_claim()
        self.add_dependency(
            "dependency:check",
            "check:effect",
            DependencyRequirement.REQUIRED_CHECK,
        )

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.BLOCKED)
        self.assertIn("SELF_VERIFICATION", decision.reason_codes)

    def test_withdrawn_claim_does_not_delete_evidence(self):
        self.add_artifact()
        self.add_claim(disposition=ClaimDisposition.WITHDRAWN)

        decision = self.evaluate()

        self.assertEqual(decision.status, DecisionStatus.WITHDRAWN)
        self.assertIn("artifact:source", self.ledger.snapshot().records)


if __name__ == "__main__":
    unittest.main()
