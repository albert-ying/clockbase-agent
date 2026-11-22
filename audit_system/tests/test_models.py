import unittest
from dataclasses import FrozenInstanceError

from audit_system.models import (
    ActorRecord,
    ActorRole,
    ActorType,
    CheckRecord,
    CheckStatus,
    RecordValidationError,
    record_to_dict,
)


class ModelTests(unittest.TestCase):
    def test_records_are_immutable(self):
        actor = ActorRecord(
            record_id="actor:checker",
            release_id="release:test",
            created_at="2026-09-07T00:00:00Z",
            actor_type=ActorType.PROGRAM,
            role=ActorRole.DETERMINISTIC_CHECKER,
            organization="ClockBase audit",
            model_family=None,
            authority_scope=("file_identity",),
        )

        with self.assertRaises(FrozenInstanceError):
            actor.role = ActorRole.REPAIRER

    def test_record_serialization_uses_enum_values_and_json_arrays(self):
        check = CheckRecord(
            record_id="check:file-identity",
            release_id="release:test",
            created_at="2026-09-07T00:00:00Z",
            subject_id="artifact:source",
            predicate="file_identity",
            status=CheckStatus.PASS,
            checker_actor_id="actor:checker",
            evidence_artifact_ids=("artifact:source",),
            expected={"sha256": "abc"},
            observed={"sha256": "abc"},
        )

        payload = record_to_dict(check)

        self.assertEqual(payload["record_type"], "check")
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(payload["evidence_artifact_ids"], ["artifact:source"])

    def test_empty_record_identifier_is_rejected(self):
        with self.assertRaisesRegex(RecordValidationError, "record_id"):
            ActorRecord(
                record_id="",
                release_id="release:test",
                created_at="2026-09-07T00:00:00Z",
                actor_type=ActorType.PROGRAM,
                role=ActorRole.DETERMINISTIC_CHECKER,
                organization="ClockBase audit",
                model_family=None,
                authority_scope=("file_identity",),
            )

    def test_model_actor_requires_model_family(self):
        with self.assertRaisesRegex(RecordValidationError, "model_family"):
            ActorRecord(
                record_id="actor:model",
                release_id="release:test",
                created_at="2026-09-07T00:00:00Z",
                actor_type=ActorType.MODEL,
                role=ActorRole.AUDITOR,
                organization="Independent audit",
                model_family=None,
                authority_scope=(),
            )

    def test_unknown_enum_value_is_rejected(self):
        with self.assertRaises(ValueError):
            CheckStatus("CONFIDENT")


if __name__ == "__main__":
    unittest.main()
