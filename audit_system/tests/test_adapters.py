import json
import tempfile
import unittest
from pathlib import Path

from audit_system.adapters import adapt_measurement_pilot, adapt_protected_repair
from audit_system.hashing import ArtifactIdentityError, inspect_file
from audit_system.models import (
    ArtifactRecord,
    CheckRecord,
    CheckStatus,
    record_to_dict,
)


class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)
        self.source = self.directory / "archived-export.csv"
        self.source.write_text("sample,value\nS1,1\n", encoding="utf-8")
        self.module = self.directory / "measurement.py"
        self.module.write_text("VERSION = 'test'\n", encoding="utf-8")
        source_identity = inspect_file(self.source)
        module_identity = inspect_file(self.module)
        self.pilot_path = self.directory / "pilot.json"
        self.pilot_path.write_text(
            json.dumps(
                {
                    "scientific_status": "NOT_EVALUATED",
                    "release_ready": False,
                    "input_kind": "real",
                    "input_source": {
                        "path": str(self.source),
                        "sha256": source_identity.sha256,
                        "size_bytes": source_identity.byte_count,
                    },
                    "measurement_module": {
                        "path": str(self.module),
                        "sha256": module_identity.sha256,
                        "size_bytes": module_identity.byte_count,
                        "expected_sha256": module_identity.sha256,
                    },
                    "rows_scanned": 3,
                    "rows_selected": 2,
                    "study": "GSE_TEST",
                    "audit": {
                        "version": "measurement-preflight-test",
                        "scope": "single_study_single_clock_preflight",
                        "value_column": "clock_a",
                        "row_count": 2,
                        "preflight_status": "FAIL",
                        "checks": [
                            {
                                "name": "finite_predictions",
                                "status": "PASS",
                                "message": "Values are finite.",
                                "details": {"invalid_count": 0},
                            },
                            {
                                "name": "reference_parity",
                                "status": "BLOCKED",
                                "message": "Reference is absent.",
                                "details": {},
                            },
                        ],
                        "release_ready": False,
                    },
                },
                indent=2,
            ),
            encoding="utf-8",
        )

    def tearDown(self):
        self.temporary_directory.cleanup()

    def test_measurement_adapter_preserves_distinct_statuses_and_bytes(self):
        before = self.pilot_path.read_bytes()

        records = adapt_measurement_pilot(
            self.pilot_path,
            "release:pilot",
        )

        self.assertEqual(before, self.pilot_path.read_bytes())
        checks = {
            record.predicate: record
            for record in records
            if isinstance(record, CheckRecord)
        }
        self.assertEqual(checks["preflight_status"].status, CheckStatus.FAIL)
        self.assertEqual(
            checks["scientific_status"].observed["raw_status"],
            "NOT_EVALUATED",
        )
        self.assertEqual(checks["release_ready"].status, CheckStatus.BLOCKED)
        self.assertEqual(
            checks["reference_parity"].status,
            CheckStatus.BLOCKED,
        )
        payload = json.dumps([record_to_dict(record) for record in records])
        self.assertIn(inspect_file(self.pilot_path).sha256, payload)

    def test_adapter_rejects_tampered_report_against_expected_identity(self):
        expected = inspect_file(self.pilot_path)
        self.pilot_path.write_text("{}", encoding="utf-8")

        with self.assertRaises(ArtifactIdentityError):
            adapt_measurement_pilot(
                self.pilot_path,
                "release:pilot",
                expected_identity=expected,
            )

    def test_protected_repair_adapter_registers_artifacts_and_blocks_authority(self):
        protected = self.directory / "protected"
        protected.mkdir()
        original = protected / "original.csv"
        corrected = protected / "corrected.csv"
        original.write_text("id,value\n1,1\n", encoding="utf-8")
        corrected.write_text("id,value\n1,2\n", encoding="utf-8")
        summary = protected / "summary.json"
        summary.write_text(
            json.dumps({"schema_version": "repair-summary-v1"}),
            encoding="utf-8",
        )
        manifest = protected / "manifest.json"
        manifest.write_text(
            json.dumps(
                {
                    "schema_version": "repair-manifest-v1",
                    "artifacts": {
                        "original_spine": {
                            "path": original.name,
                            "sha256": inspect_file(original).sha256,
                            "row_count": 1,
                        },
                        "corrected_spine": {
                            "path": corrected.name,
                            "sha256": inspect_file(corrected).sha256,
                            "row_count": 1,
                        },
                        "validation_summary": {
                            "path": summary.name,
                            "sha256": inspect_file(summary).sha256,
                            "row_count": 1,
                        },
                    },
                }
            ),
            encoding="utf-8",
        )
        before = {path: path.read_bytes() for path in protected.iterdir()}

        records = adapt_protected_repair(
            manifest,
            summary,
            "release:repair",
        )

        self.assertEqual(
            before,
            {path: path.read_bytes() for path in protected.iterdir()},
        )
        artifacts = {
            record.record_id: record
            for record in records
            if isinstance(record, ArtifactRecord)
        }
        self.assertTrue(
            any("original_spine" in record_id for record_id in artifacts)
        )
        self.assertTrue(
            any("corrected_spine" in record_id for record_id in artifacts)
        )
        checks = {
            record.predicate: record
            for record in records
            if isinstance(record, CheckRecord)
        }
        self.assertEqual(
            checks["repair_actor_identity"].status,
            CheckStatus.BLOCKED,
        )
        self.assertEqual(
            checks["verification_authority"].status,
            CheckStatus.BLOCKED,
        )


if __name__ == "__main__":
    unittest.main()
