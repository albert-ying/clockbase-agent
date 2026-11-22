import tempfile
import unittest
from pathlib import Path

from audit_system.ledger import CorruptLedgerError, DuplicateRecordError, Ledger
from audit_system.models import (
    ActorRecord,
    ActorRole,
    ActorType,
    RecordValidationError,
)


class LedgerTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)
        self.ledger = Ledger(self.directory / "records.jsonl")

    def tearDown(self):
        self.temporary_directory.cleanup()

    @staticmethod
    def actor(record_id: str, organization: str = "ClockBase audit"):
        return ActorRecord(
            record_id=record_id,
            release_id="release:test",
            created_at="2026-09-07T00:00:00Z",
            actor_type=ActorType.PROGRAM,
            role=ActorRole.DETERMINISTIC_CHECKER,
            organization=organization,
            model_family=None,
            authority_scope=("file_identity",),
        )

    def test_append_round_trips_a_validated_record(self):
        digest = self.ledger.append(self.actor("actor:one"))
        snapshot = self.ledger.snapshot()

        self.assertEqual(snapshot.sha256, digest)
        self.assertEqual(snapshot.record_ids, ("actor:one",))
        self.assertEqual(
            snapshot.records["actor:one"].organization,
            "ClockBase audit",
        )
        self.assertTrue(self.ledger.path.read_bytes().endswith(b"\n"))

    def test_byte_identical_duplicate_is_idempotent(self):
        record = self.actor("actor:one")
        first = self.ledger.append(record)
        bytes_after_first = self.ledger.path.read_bytes()

        second = self.ledger.append(record)

        self.assertEqual(first, second)
        self.assertEqual(bytes_after_first, self.ledger.path.read_bytes())
        self.assertEqual(self.ledger.snapshot().record_ids, ("actor:one",))

    def test_conflicting_record_identifier_is_rejected(self):
        self.ledger.append(self.actor("actor:one", organization="A"))
        before = self.ledger.path.read_bytes()

        with self.assertRaisesRegex(DuplicateRecordError, "actor:one"):
            self.ledger.append(self.actor("actor:one", organization="B"))

        self.assertEqual(before, self.ledger.path.read_bytes())

    def test_invalid_record_leaves_ledger_unchanged(self):
        self.ledger.append(self.actor("actor:one"))
        before = self.ledger.path.read_bytes()

        with self.assertRaises(RecordValidationError):
            self.ledger.append(object())

        self.assertEqual(before, self.ledger.path.read_bytes())

    def test_malformed_existing_line_blocks_append_without_repair(self):
        self.ledger.path.write_bytes(b'{"record_id":"broken"\n')
        before = self.ledger.path.read_bytes()

        with self.assertRaisesRegex(CorruptLedgerError, "line 1"):
            self.ledger.append(self.actor("actor:one"))

        self.assertEqual(before, self.ledger.path.read_bytes())

    def test_snapshot_rejects_two_payloads_with_the_same_identifier(self):
        first = self.actor("actor:one", organization="A")
        second = self.actor("actor:one", organization="B")
        from audit_system.hashing import canonical_json_bytes
        from audit_system.models import record_to_dict

        self.ledger.path.write_bytes(
            canonical_json_bytes(record_to_dict(first))
            + b"\n"
            + canonical_json_bytes(record_to_dict(second))
            + b"\n"
        )

        with self.assertRaisesRegex(CorruptLedgerError, "actor:one"):
            self.ledger.snapshot()


if __name__ == "__main__":
    unittest.main()
