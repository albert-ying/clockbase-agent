import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from audit_system.hashing import inspect_file
from audit_system.ledger import Ledger


class BuildW0LedgerTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary_directory.name)
        self.pilot, self.manifest, self.summary = self.write_reports()

    def tearDown(self):
        self.temporary_directory.cleanup()

    def write_reports(self):
        source = self.directory / "source.csv"
        source.write_text("sample,value\nS1,1\n", encoding="utf-8")
        module = self.directory / "measurement.py"
        module.write_text("VERSION = 'test'\n", encoding="utf-8")
        pilot = self.directory / "pilot.json"
        pilot.write_text(
            json.dumps(
                {
                    "scientific_status": "NOT_EVALUATED",
                    "release_ready": False,
                    "input_kind": "real",
                    "input_source": {
                        "path": str(source),
                        "sha256": inspect_file(source).sha256,
                        "size_bytes": inspect_file(source).byte_count,
                    },
                    "measurement_module": {
                        "path": str(module),
                        "sha256": inspect_file(module).sha256,
                        "size_bytes": inspect_file(module).byte_count,
                    },
                    "study": "GSE_TEST",
                    "audit": {
                        "version": "measurement-preflight-test",
                        "preflight_status": "FAIL",
                        "checks": [],
                    },
                }
            ),
            encoding="utf-8",
        )
        original = self.directory / "original.csv"
        corrected = self.directory / "corrected.csv"
        original.write_text("id,value\n1,1\n", encoding="utf-8")
        corrected.write_text("id,value\n1,2\n", encoding="utf-8")
        summary = self.directory / "summary.json"
        summary.write_text(
            json.dumps({"schema_version": "summary-v1"}),
            encoding="utf-8",
        )
        manifest = self.directory / "manifest.json"
        manifest.write_text(
            json.dumps(
                {
                    "schema_version": "manifest-v1",
                    "artifacts": {
                        "original_spine": {
                            "path": original.name,
                            "sha256": inspect_file(original).sha256,
                        },
                        "corrected_spine": {
                            "path": corrected.name,
                            "sha256": inspect_file(corrected).sha256,
                        },
                        "validation_summary": {
                            "path": summary.name,
                            "sha256": inspect_file(summary).sha256,
                        },
                    },
                }
            ),
            encoding="utf-8",
        )
        return pilot, manifest, summary

    def run_build(self, output, pilot=None):
        environment = os.environ.copy()
        project_root = Path(__file__).resolve().parents[2]
        environment["PYTHONPATH"] = str(project_root)
        return subprocess.run(
            [
                sys.executable,
                "-m",
                "audit_system.build_w0_ledger",
                "--measurement-pilot",
                str(pilot or self.pilot),
                "--protected-manifest",
                str(self.manifest),
                "--protected-summary",
                str(self.summary),
                "--release-id",
                "release:test-build",
                "--output",
                str(output),
            ],
            cwd=project_root,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )

    def test_build_creates_valid_ledger_and_machine_readable_summary(self):
        output = self.directory / "ledger.jsonl"
        inputs = {
            path: path.read_bytes()
            for path in (self.pilot, self.manifest, self.summary)
        }

        completed = self.run_build(output)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        summary = json.loads(completed.stdout)
        snapshot = Ledger(output).snapshot()
        self.assertEqual(summary["release_id"], "release:test-build")
        self.assertEqual(summary["ledger_sha256"], snapshot.sha256)
        self.assertEqual(summary["record_count"], len(snapshot.record_ids))
        self.assertGreater(summary["record_count"], 10)
        self.assertEqual(
            inputs,
            {path: path.read_bytes() for path in inputs},
        )

    def test_existing_output_is_not_overwritten(self):
        output = self.directory / "ledger.jsonl"
        output.write_text("protected\n", encoding="utf-8")

        completed = self.run_build(output)

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("refusing to overwrite", completed.stderr.lower())
        self.assertEqual(output.read_text(encoding="utf-8"), "protected\n")

    def test_invalid_input_leaves_no_partial_output(self):
        output = self.directory / "ledger.jsonl"
        invalid = self.directory / "invalid-pilot.json"
        invalid.write_text("{}", encoding="utf-8")

        completed = self.run_build(output, pilot=invalid)

        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("lacks an audit version", completed.stderr.lower())
        self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
