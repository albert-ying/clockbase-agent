"""Software-fixture tests for Luna batch resume state and integrity gates."""

from __future__ import annotations

import json
import importlib.util
import io
from contextlib import redirect_stdout
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

from audit_system.tests import require_development_artifact

_ROOT = Path(__file__).resolve().parents[2]
_RUNNER = _ROOT / "paper/revision_2026-09-07/evidence/saved-analysis-corpus-audit-2026-09-19/run_legacy_review_batch_2026_09_19.py"


class BatchResumeFixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        global batch
        require_development_artifact(_RUNNER)
        spec = importlib.util.spec_from_file_location("legacy_batch_runner_fixture", _RUNNER)
        batch = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = batch
        assert spec.loader is not None
        spec.loader.exec_module(batch)

    def test_scientific_batch_completes_three_case_preflight_without_transport(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            plans = []
            for case_id in ("human_rna-GSE185657-001", "human_rna-GSE120159-001", "human_dnam-GSE71837-001"):
                path = root / "cases" / case_id / "plan.json"
                (path.parent / "prepared").mkdir(parents=True)
                path.write_text(json.dumps({"model": "gpt-5.6-luna", "review_profile": "scientific_checks_v2", "limits": {"max_calls": 1, "max_cost": 0.2}}))
                (path.parent / "prepared" / "prepared-request.json").write_text("software fixture request")
                plans.append({"case_id": case_id, "plan": {"sha256": batch.sha256(path)}})
            (root / "preparation-manifest.json").write_text(json.dumps({"review_profile": "scientific_checks_v2", "aggregate_budget": {"model": "gpt-5.6-luna"}, "plans": plans}))
            output = root / "receipt.json"
            validated = {"study_id": "software-fixture", "source_design_status": "COMPLETED", "trace_request_status": "CANNOT VERIFY", "analysis_preflight": {"status": "PASS"}, "prediction_preflight": {"status": "CANNOT VERIFY"}}
            with mock.patch.object(sys, "argv", ["runner", "--root", str(root), "--output", str(output), "--check"]), mock.patch.object(batch, "validate_plan", return_value={"validated": validated}), mock.patch.object(batch, "run") as transport, redirect_stdout(io.StringIO()):
                self.assertEqual(batch.main(), 0)
            transport.assert_not_called()
            receipt = json.loads(output.read_text())
            self.assertEqual(receipt["status"], "PASS")
            self.assertEqual(receipt["plan_count"], 3)
            self.assertEqual(receipt["counts"], {"pass": 3, "blocked": 0, "fail": 0, "pending": 0})

    def test_pending_dominates_intermediate_failure(self):
        status, execution, counts = batch._batch_state(
            [{"status": "FAIL"}, {"status": "PASS"}], 3, False)
        self.assertEqual(status, "INCOMPLETE")
        self.assertEqual(execution, "RUN_IN_PROGRESS")
        self.assertEqual(counts, {"pass": 1, "blocked": 0, "fail": 1, "pending": 1})

    def _prepared(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        plan = root / "plan.json"
        plan.write_text("{}")
        prepared = root / "prepared"
        prepared.mkdir()
        request = prepared / "prepared-request.json"
        result = prepared / "prepare-result.json"
        request.write_text("request")
        result.write_text("result")
        import hashlib
        digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
        manifest = prepared / "manifest.json"
        manifest.write_text(json.dumps({"plan_sha256": digest(plan), "prepared_request_sha256": digest(request), "prepare_result_sha256": digest(result)}))
        return root, plan, request, digest

    def test_modified_input_rejects_recorded_pass(self):
        root, plan, request, digest = self._prepared()
        row = {"plan_sha256": digest(plan), "prepared_request_sha256": digest(request)}
        plan.write_text('{"changed":true}')
        with mock.patch.object(batch, "validate_plan"):
            with self.assertRaises(ValueError):
                batch._validate_recorded_pass(plan, row, run_mode=False)

    def test_run_resume_requires_completed_output_and_validates_cached_output(self):
        root, plan, request, digest = self._prepared()
        row = {"plan_sha256": digest(plan), "prepared_request_sha256": digest(request)}
        with mock.patch.object(batch, "validate_plan"):
            with self.assertRaises(ValueError):
                batch._validate_recorded_pass(plan, row, run_mode=True)
        (plan.parent / "prepared" / "review-result.json").write_text("cached")
        with mock.patch.object(batch, "validate_plan"), mock.patch.object(batch, "run", return_value={"status": "PASS"}) as cached:
            batch._validate_recorded_pass(plan, row, run_mode=True)
        cached.assert_called_once()


if __name__ == "__main__":
    unittest.main()
