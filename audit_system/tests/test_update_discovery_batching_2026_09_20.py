"""Offline regression fixtures for bounded GEO ESummary discovery and tariffs."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from src import update_pipeline_2026_09_20 as update


class DiscoveryBatchingTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)

    def fixture(self, count=205, malformed_batch=None, max_requests=10):
        seen = self.root / "seen.jsonl"
        seen.write_text(json.dumps({"study_id": "GSE1"}) + "\n")
        ids = [str(value) for value in range(1, count + 1)]
        plan = {
            "schema": "clockbase-update-discovery-plan-v1", "query": "senescence",
            "modality": "human_rna", "date_from": "2025-01-01", "date_to": "2025-12-31",
            "seen_worklist": {"path": str(seen), "sha256": update.sha(seen)},
            "limits": {"max_requests": max_requests, "max_response_bytes": 1000000,
                       "max_total_response_bytes": 10000000, "max_seconds": 30,
                       "max_query_ids": count, "max_candidates": count,
                       "minimum_interval_seconds": 0.35}}
        path = self.root / "plan.json"
        path.write_text(json.dumps(plan))

        def transport(url, cap):
            if "esearch" in url:
                return json.dumps({"esearchresult": {"idlist": ids, "count": str(count)}}).encode()
            requested = update.urllib.parse.parse_qs(update.urllib.parse.urlparse(url).query)["id"][0].split(",")
            returned = list(requested)
            if malformed_batch is not None and requested[0] == str(malformed_batch):
                returned = requested[1:] + ["999999"]
            result = {uid: {"accession": "GSE1", "pdat": "2025/01/01"} for uid in returned}
            result["uids"] = returned
            return json.dumps({"result": result}).encode()
        return path, transport

    def test_more_than_one_hundred_uids_are_batched_and_completed_resume_is_offline(self):
        path, transport = self.fixture()
        output = self.root / "output"
        with patch.object(update.time, "sleep"):
            result = update.discover(path, output, transport=transport)
        batches = json.loads((output / "esummary-batches.json").read_text())
        self.assertEqual([len(row["requested_uids"]) for row in batches["batches"]], [100, 100, 5])
        self.assertEqual([row["status"] for row in batches["batches"]], ["PASS", "PASS", "PASS"])
        self.assertEqual(result["requests_used"], 4)
        self.assertEqual(result["candidate_status_counts"], {"SEEN_ACCESSION": 205})
        def no_network(*unused): raise AssertionError("Completed discovery resume made a request")
        self.assertEqual(result, update.discover(path, output, transport=no_network))

    def test_request_cap_marks_unrequested_batch_uids_unresolved(self):
        path, transport = self.fixture(count=150, max_requests=2)
        output = self.root / "output"
        with patch.object(update.time, "sleep"):
            result = update.discover(path, output, transport=transport)
        batches = json.loads((output / "esummary-batches.json").read_text())
        self.assertEqual([row["status"] for row in batches["batches"]], ["PASS", "FAIL"])
        self.assertIsNone(batches["batches"][1]["receipt_path"])
        self.assertEqual(result["candidate_status_counts"], {"SEEN_ACCESSION": 100, "UNRESOLVED": 50})
        self.assertEqual(sum(result["candidate_status_counts"].values()), 150)

    def test_missing_and_unexpected_summary_uids_fail_only_the_affected_batch(self):
        path, transport = self.fixture(count=150, malformed_batch=101)
        output = self.root / "output"
        with patch.object(update.time, "sleep"):
            result = update.discover(path, output, transport=transport)
        batches = json.loads((output / "esummary-batches.json").read_text())
        self.assertEqual([row["status"] for row in batches["batches"]], ["PASS", "FAIL"])
        self.assertIn("UID mismatch", batches["batches"][1]["error"])
        self.assertEqual(result["candidate_status_counts"], {"SEEN_ACCESSION": 100, "UNRESOLVED": 50})

    def test_batch_manifest_tampering_stops_completed_resume(self):
        path, transport = self.fixture(count=101)
        output = self.root / "output"
        with patch.object(update.time, "sleep"):
            update.discover(path, output, transport=transport)
        (output / "esummary-batches.json").write_text("{}\n")
        with self.assertRaisesRegex(update.UpdateIntegrityError, "Completed child output changed"):
            update.discover(path, output, transport=transport)


class UpdateTariffContractTests(unittest.TestCase):
    def plan(self, pipeline=None):
        return {"schema": update.SCHEMA, "studies": [], "pipeline": pipeline or {},
                "limits": {"max_studies": 1, "max_input_bytes": 1, "max_samples_per_study": 1,
                           "max_child_seconds": 1, "concurrency": 1, "allow_network": False,
                           "max_total_model_calls": 1, "max_total_cost_usd": 1.0},
                "comprehensive_audit": {"model": "gpt-5.6-luna", "max_calls": 1,
                                        "max_cost_usd": 1.0, "input_usd_per_million": 0.5,
                                        "output_usd_per_million": 1.8}}

    def test_default_luna_tariffs_are_current_conservative_bounds(self):
        _, options, _ = update.validate_plan(self.plan())
        self.assertEqual((options["input_usd_per_million"], options["output_usd_per_million"]), (0.5, 1.8))
        self.assertEqual((options["audit_input_usd_per_million"], options["audit_output_usd_per_million"]), (0.5, 1.8))

    def test_nondefault_producer_model_requires_both_explicit_tariffs(self):
        with self.assertRaisesRegex(ValueError, "non-default model"):
            update.validate_plan(self.plan({"model": "gpt-5.6-sol"}))
        _, options, _ = update.validate_plan(self.plan({"model": "gpt-5.6-sol",
            "input_usd_per_million": 2.0, "output_usd_per_million": 8.0}))
        self.assertEqual(options["model"], "gpt-5.6-sol")

    def test_nondefault_audit_model_requires_both_explicit_tariffs(self):
        with self.assertRaisesRegex(ValueError, "non-default audit_model"):
            update.validate_plan(self.plan({"audit_model": "gpt-5.6-sol"}))

    def test_partial_role_tariff_pairs_are_rejected(self):
        for pipeline in ({"input_usd_per_million": 1.0}, {"audit_output_usd_per_million": 2.0},
                         {"model": "gpt-5.6-sol", "input_usd_per_million": 2.0},
                         {"audit_model": "gpt-5.6-sol", "audit_output_usd_per_million": 8.0}):
            with self.subTest(pipeline=pipeline), self.assertRaisesRegex(ValueError, "tariff"):
                update.validate_plan(self.plan(pipeline))


if __name__ == "__main__":
    unittest.main()
