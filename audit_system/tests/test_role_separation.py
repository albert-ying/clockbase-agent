import unittest

from audit_system.models import (
    ActorRecord,
    ActorRole,
    ActorType,
    IdentityAssurance,
)
from audit_system.policy import AuthorityPolicy


class RoleSeparationTests(unittest.TestCase):
    def setUp(self):
        self.policy = AuthorityPolicy()
        self.base = {
            "release_id": "release:test",
            "created_at": "2026-09-07T00:00:00Z",
            "organization": "ClockBase audit",
            "identity_assurance": IdentityAssurance.VERIFIED,
        }

    def actor(self, record_id, actor_type, role, family=None, scope=()):
        return ActorRecord(
            record_id=record_id,
            actor_type=actor_type,
            role=role,
            model_family=family,
            authority_scope=scope,
            **self.base,
        )

    def test_repairer_cannot_verify_own_repair(self):
        primary = self.actor(
            "actor:primary",
            ActorType.MODEL,
            ActorRole.PRIMARY_ANALYST,
            "deepseek",
        )
        repairer = self.actor(
            "actor:repairer",
            ActorType.HUMAN,
            ActorRole.REPAIRER,
            scope=("effect_recalculation",),
        )

        result = self.policy.verify_separation(
            producer=primary,
            repairer=repairer,
            verifier=repairer,
        )

        self.assertFalse(result.allowed)
        self.assertEqual(result.reason_code, "SELF_VERIFICATION")

    def test_same_model_family_is_not_independent(self):
        producer = self.actor(
            "actor:primary",
            ActorType.MODEL,
            ActorRole.PRIMARY_ANALYST,
            "gpt",
        )
        verifier = self.actor(
            "actor:auditor",
            ActorType.MODEL,
            ActorRole.AUDITOR,
            "gpt",
        )

        result = self.policy.verify_separation(
            producer=producer,
            repairer=None,
            verifier=verifier,
        )

        self.assertFalse(result.allowed)
        self.assertEqual(result.reason_code, "SAME_MODEL_FAMILY")

    def test_unverified_model_identity_is_not_independent(self):
        producer = self.actor(
            "actor:primary",
            ActorType.MODEL,
            ActorRole.PRIMARY_ANALYST,
            "deepseek",
        )
        verifier = ActorRecord(
            record_id="actor:auditor",
            release_id="release:test",
            created_at="2026-09-07T00:00:00Z",
            actor_type=ActorType.MODEL,
            role=ActorRole.AUDITOR,
            organization="External audit",
            model_family="gpt",
            authority_scope=(),
        )

        result = self.policy.verify_separation(
            producer=producer,
            repairer=None,
            verifier=verifier,
        )

        self.assertFalse(result.allowed)
        self.assertEqual(result.reason_code, "UNVERIFIED_MODEL_IDENTITY")

    def test_model_cannot_authoritatively_pass_numerical_check(self):
        auditor = self.actor(
            "actor:auditor",
            ActorType.MODEL,
            ActorRole.AUDITOR,
            "gpt",
            scope=("effect_recalculation",),
        )

        result = self.policy.check_can_pass(
            auditor,
            "effect_recalculation",
        )

        self.assertFalse(result.allowed)
        self.assertEqual(result.reason_code, "MODEL_CANNOT_AUTHORIZE_HARD_CHECK")

    def test_deterministic_checker_can_pass_declared_predicate(self):
        checker = self.actor(
            "actor:checker",
            ActorType.PROGRAM,
            ActorRole.DETERMINISTIC_CHECKER,
            scope=("file_identity",),
        )

        result = self.policy.check_can_pass(checker, "file_identity")

        self.assertTrue(result.allowed)
        self.assertEqual(result.reason_code, "AUTHORIZED")


if __name__ == "__main__":
    unittest.main()
