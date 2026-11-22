"""Focused software-fixture tests for the preserved-analysis review adapter.

These fixtures exercise contracts only; they are not empirical study results.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from langchain_core.messages import AIMessage

from audit_system import legacy_analysis_review_2026_09_19 as review


class LegacyReviewFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.files = {}
        self.study = "GSE123"
        self._write("metadata.csv", 'GSM,title\nGSM1,"case\nwith context"\nGSM2,control\n')
        self._write("metadata-prompt.csv", 'GSM,title\nGSM1,"case\nwith context"\nGSM2,control\n')
        packet = {"study_id": self.study, "source_sha256": "a" * 64,
                  "samples": [{"sample_id": "GSM1", "title": "case GSM1"},
                              {"sample_id": "GSM2", "title": "control GSM2"}]}
        self._write_json("source-packet.json", packet)
        self._write("predictions.csv", "sample_id,ClockA\nGSM1,1.0\nGSM2,2.0\n")
        self._write("analysis.csv", "gse_id,clock,case_ids,control_ids,comparison\nGSE123,ClockA,GSM1,GSM2,legacy comparison\n")
        message = "================================ Human Message =================================\n\nStudy identifier:GSE123\nReview preserved output."
        trace = {"schema": "clockbase-preserved-trace-v1", "study_id": self.study,
                 "entry": {"study_id": self.study, "status": "success",
                            "result": [[message], "<solution>done</solution>"]}}
        self._write_json("selected-trace.json", trace)
        self._write_json("trace-check.json", {"trace_request_status": "PASS"})
        self._write_json("source-design.json", {"schema": "source-design", "study_id": self.study,
                                                "scientific_status": "CANNOT VERIFY"})
        self._write_json("source-design-config.json", {"schema": "source-design-config", "plan": {
            "study_id": self.study, "source_sha256": {"metadata.csv": "PLACEHOLDER"}},
            "source_packet_sha256": "PLACEHOLDER"})
        self._write_json("source-preflight-packet.json", packet)
        self._write("reviewed_code.py", "def preserved_fixture():\n    return 'software fixture'\n")
        self.plan_path = self.root / "plan.json"
        self._write_json("plan.json", self._plan())

    def tearDown(self):
        self.tmp.cleanup()

    def _write(self, name, text):
        path = self.root / name
        if isinstance(text, bytes):
            path.write_bytes(text)
        else:
            path.write_text(text, encoding="utf-8", newline="")
        self.files[name] = path
        return path

    def _write_json(self, name, value):
        return self._write(name, json.dumps(value, ensure_ascii=False, indent=2) + "\n")

    def _record(self, name):
        path = self.files[name]
        return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}

    def _plan(self):
        packet_record = self._record("source-packet.json")
        metadata_record = self._record("metadata.csv")
        metadata_prompt_record = self._record("metadata-prompt.csv")
        config_record = self._record("source-design-config.json")
        packet_canonical = review.canonical_json_sha256(json.loads(self.files["source-packet.json"].read_text()))
        config_value = json.loads(self.files["source-design-config.json"].read_text())
        config_value["plan"]["source_sha256"]["metadata.csv"] = metadata_record["sha256"]
        config_value["source_packet_sha256"] = packet_canonical
        self._write_json("source-design-config.json", config_value)
        config_record = self._record("source-design-config.json")
        source_design_record = self._record("source-design.json")
        source_design_record.update({"status": "COMPLETED", "metadata_sha256": metadata_record["sha256"],
                                     "source_packet_sha256": packet_canonical,
                                     "source_packet_canonical_sha256": packet_canonical,
                                     "completion_sha256": source_design_record["sha256"],
                                     "run_config_sha256": config_record["sha256"],
                                     "preflight_source_packet": self._record("source-preflight-packet.json"),
                                     "completion_output_hashes": {"source-design.json": source_design_record["sha256"]}})
        return {
            "schema": review.SCHEMA, "split": "development", "provider": "openai",
            "model": "gpt-5.6-luna", "study_id": self.study, "modality": "RNA",
            "context_scope": {
                "saved_analysis": "software fixture CSV",
                "predictions": "software fixture companion CSV",
                "historical_lineage": "not established in fixture",
                "source_design_agent": "software fixture",
                "scientific_acceptance": "CANNOT VERIFY",
            },
            "limits": {"max_calls": 1, "max_cost": 1.0, "output_tokens": 500,
                       "prompt_bytes": 200000, "input_per_million": 1.0,
                       "output_per_million": 1.0},
            "artifacts": {
                "metadata": metadata_record,
                "metadata_prompt": {**metadata_prompt_record, "original_sha256": metadata_record["sha256"],
                                     "redaction_profile": "software-fixture allowlist"},
                "source_packet": packet_record,
                "analysis_csv": self._record("analysis.csv"),
                "predictions": self._record("predictions.csv"),
                "selected_trace": self._record("selected-trace.json"),
                "trace_request": self._record("trace-check.json"),
                "source_design": source_design_record,
                "source_design_config": config_record,
                "code": [self._record("reviewed_code.py")],
            },
        }

    def _response(self, **overrides):
        domains = {}
        for domain in review.DOMAINS:
            domains[domain] = {"status": "CANNOT VERIFY", "summary": "fixture review",
                               "evidence": [{"artifact_id": "analysis_csv", "quote": "GSM1"}],
                               "unresolved": ["software fixture has no scientific acceptance"],
                               "proposed_repair": None, "verification_needed": ["independent verification"]}
        value = {"schema": review.REVIEW_SCHEMA, "study_id": self.study,
                 "scientific_status": "CANNOT VERIFY", "domains": domains,
                 "original_membership": [{"target": "fixture", "status": "CANNOT VERIFY",
                    "case_ids": [], "control_ids": [], "excluded_ids": [],
                    "evidence": [{"artifact_id": "analysis_csv", "quote": "GSM1"},
                                 {"artifact_id": "analysis_csv", "quote": "GSM2"}]}]}
        value.update(overrides)
        return json.dumps(value, ensure_ascii=False)

    def _prepare(self):
        out = self.root / "review"
        review.prepare(self.plan_path, out)
        return out

    def _scientific_response(self):
        plan = self._plan()
        plan["review_profile"] = review.SCIENTIFIC_PROFILE
        self._write_json("plan.json", plan)
        value = json.loads(self._response())
        value["scientific_checks"] = {
            key: {"status": "CANNOT VERIFY", "scope": "Software fixture comparison", "summary": "Software fixture lacks design evidence",
                  "evidence": [{"artifact_id": "analysis_csv", "quote": "GSM1"}],
                  "unresolved": ["Independent biological unit not supplied"],
                  "proposed_repair": None, "verification_needed": ["Source verification required"]}
            for key in review.SCIENTIFIC_CHECKS
        }
        return value

    def test_scientific_profile_requires_explicit_questions_and_remains_advisory(self):
        response = self._scientific_response()
        out = self._prepare()
        request = json.loads((out / "prepared-request.json").read_text())
        prompt = json.loads(request["messages"][1]["content"])
        self.assertEqual(set(prompt["scientific_check_questions"]), set(review.SCIENTIFIC_CHECKS))
        self.assertEqual(len(review.SCIENTIFIC_CHECKS), 16)
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=json.dumps(response))):
            result = review.run(self.plan_path, out)
        saved = json.loads((out / "review.json").read_text())
        self.assertEqual(saved["scientific_checks_summary"]["counts"]["CANNOT VERIFY"], 16)
        self.assertEqual(result["scientific_status"], "CANNOT VERIFY")
        self.assertEqual(result["accepted_repairs"], [])
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke") as invoke:
            self.assertEqual(review.run(self.plan_path, out), result)
            invoke.assert_not_called()

    def test_scientific_profile_rejects_missing_method_grouping_or_claim_checks(self):
        for key in ("group_membership", "pairing_repeated_measures", "estimand_and_test_choice", "scientific_claim_scope"):
            with self.subTest(key=key):
                response = self._scientific_response()
                del response["scientific_checks"][key]
                out = self.root / key
                review.prepare(self.plan_path, out)
                with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke", return_value=AIMessage(content=json.dumps(response))):
                    with self.assertRaisesRegex(ValueError, "sixteen scientific checks"):
                        review.run(self.plan_path, out)

    def test_scientific_pass_rejects_unresolved_or_uncited_assertions(self):
        for violation in ("unresolved", "uncited", "bad_quote", "no_reason", "no_scope", "uncited_unknown"):
            with self.subTest(violation=violation):
                response = self._scientific_response()
                check = response["scientific_checks"]["distribution_and_assumptions"]
                if violation == "unresolved":
                    check["status"] = "PASS"
                elif violation == "uncited":
                    check.update(status="PASS", evidence=[], unresolved=[])
                elif violation == "bad_quote":
                    check["evidence"][0]["quote"] = "The assumptions all passed"
                elif violation == "no_reason":
                    check["unresolved"] = []
                elif violation == "no_scope":
                    check["scope"] = ""
                else:
                    check["evidence"] = []
                out = self.root / violation
                review.prepare(self.plan_path, out)
                with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke", return_value=AIMessage(content=json.dumps(response))):
                    with self.assertRaises(ValueError):
                        review.run(self.plan_path, out)

    def test_scientific_check_can_fail_without_accepting_or_repairing_result(self):
        response = self._scientific_response()
        response["scientific_checks"]["scientific_claim_scope"].update(status="FAIL", unresolved=[])
        out = self._prepare()
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke", return_value=AIMessage(content=json.dumps(response))):
            result = review.run(self.plan_path, out)
        saved = json.loads((out / "review.json").read_text())
        self.assertEqual(saved["scientific_checks"]["scientific_claim_scope"]["status"], "FAIL")
        self.assertEqual(result["scientific_status"], "CANNOT VERIFY")
        self.assertEqual(result["accepted_repairs"], [])

    def test_old_profile_has_no_new_scientific_check_requirement(self):
        out = self._prepare()
        request = json.loads((out / "prepared-request.json").read_text())
        payload = json.loads(request["messages"][1]["content"])
        self.assertNotIn("scientific_check_questions", payload)
        self.assertNotIn("review_profile", payload)

    def test_unknown_scientific_profile_is_rejected_before_transport(self):
        plan = self._plan()
        plan["review_profile"] = "unknown"
        self._write_json("plan.json", plan)
        with self.assertRaisesRegex(ValueError, "review_profile"):
            self._prepare()

    def test_prepare_preserves_multiline_csv_and_does_not_call_transport(self):
        out = self._prepare()
        request = json.loads((out / "prepared-request.json").read_text())
        self.assertIn("case\\nwith context", request["messages"][1]["content"])
        self.assertTrue((out / "manifest.json").is_file())

    def test_run_contract_and_resume_integrity(self):
        out = self._prepare()
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=self._response())) as invoke:
            result = review.run(self.plan_path, out)
        self.assertEqual(result["scientific_status"], "CANNOT VERIFY")
        self.assertEqual(result["response_contract_status"], "PASS")
        self.assertEqual(invoke.call_count, 1)
        resumed = review.run(self.plan_path, out)
        self.assertEqual(resumed["review_response_sha256"], result["review_response_sha256"])
        (out / "review.json").write_text("{}")
        with self.assertRaises(ValueError):
            review.run(self.plan_path, out)

    def test_prepared_role_order_and_input_mutations_are_rejected(self):
        out = self._prepare()
        request_path = out / "prepared-request.json"
        request = json.loads(request_path.read_text())
        request["messages"].reverse()
        request_path.write_text(json.dumps(request))
        with self.assertRaises(ValueError):
            review.run(self.plan_path, out)

    def test_wrong_gsm_and_invented_quote_are_rejected(self):
        out = self._prepare()
        wrong = json.loads(self._response())
        wrong["original_membership"][0]["case_ids"] = ["GSM999"]
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=json.dumps(wrong))):
            with self.assertRaises(ValueError):
                review.run(self.plan_path, out)
        out2 = self.root / "review2"
        review.prepare(self.plan_path, out2)
        invented = json.loads(self._response())
        invented["domains"][review.DOMAINS[0]]["evidence"][0]["quote"] = "not in source"
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=json.dumps(invented))):
                with self.assertRaises(ValueError):
                    review.run(self.plan_path, out2)

    def test_original_membership_pass_is_rejected_without_source_roster_artifact(self):
        out = self._prepare()
        response = json.loads(self._response())
        response["original_membership"][0].update({
            "status": "PASS", "case_ids": ["GSM1"], "control_ids": ["GSM2"]})
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=json.dumps(response))):
            with self.assertRaises(ValueError):
                review.run(self.plan_path, out)

    def test_unbound_trace_cannot_resolve_domain_without_non_trace_evidence(self):
        plan = self._plan()
        trace = json.loads(self.files["selected-trace.json"].read_text())
        trace["entry"]["result"][0][0] = trace["entry"]["result"][0][0].replace("GSE123", "GSE999")
        self._write_json("selected-trace.json", trace)
        self._write_json("trace-check.json", {"trace_request_status": "FAIL"})
        self._write_json("plan.json", self._plan())
        out = self._prepare()
        response = json.loads(self._response())
        response["domains"][review.DOMAINS[0]].update({
            "status": "FAIL", "evidence": [{"artifact_id": "selected_trace", "quote": "GSE999"}],
            "unresolved": []})
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=json.dumps(response))):
            with self.assertRaises(ValueError):
                review.run(self.plan_path, out)

    def test_absent_analysis_identity_is_cannot_verify(self):
        self._write("analysis.csv", "clock,case_ids,control_ids\nClockA,GSM1,GSM2\n")
        self._write_json("plan.json", self._plan())
        bundle = review.validate_plan(self.plan_path)
        self.assertEqual(bundle["validated"]["analysis_identity"]["status"], "CANNOT VERIFY")

    def test_unknown_declared_clock_is_rejected(self):
        out = self._prepare()
        response = json.loads(self._response())
        response["original_membership"][0]["clock"] = "Clock-not-in-predictions"
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=json.dumps(response))):
            with self.assertRaises(ValueError):
                review.run(self.plan_path, out)

    def test_source_mutation_after_prepare_is_rejected(self):
        out = self._prepare()
        self.files["metadata.csv"].write_text(self.files["metadata.csv"].read_text() + "\n")
        with self.assertRaises(ValueError):
            review.run(self.plan_path, out)

    def test_uncertain_billing_and_prompt_cap_stop_without_transport(self):
        out = self._prepare()
        slot = out / "api" / "legacy_review-000"
        slot.mkdir(parents=True)
        (slot / "reservation.json").write_text(json.dumps({"reserved_usd": 0.01}))
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke") as invoke:
            with self.assertRaises(RuntimeError):
                review.run(self.plan_path, out)
            invoke.assert_not_called()

        plan = self._plan()
        plan["limits"]["prompt_bytes"] = 10
        self._write_json("plan-small-cap.json", plan)
        with self.assertRaises(ValueError):
            review.prepare(self.root / "plan-small-cap.json", self.root / "small-cap")

    def test_trace_identity_failure_is_retained_as_failure(self):
        trace = json.loads(self.files["selected-trace.json"].read_text())
        trace["entry"]["result"][0][0] = trace["entry"]["result"][0][0].replace("GSE123", "GSE999")
        self._write_json("selected-trace.json", trace)
        self._write_json("trace-check.json", {"trace_request_status": "FAIL"})
        self._write_json("plan.json", self._plan())
        out = self._prepare()
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=self._response())):
            result = review.run(self.plan_path, out)
        self.assertEqual(result["trace_identity_status"], "FAIL")
        self.assertEqual(result["scientific_status"], "CANNOT VERIFY")

    def test_missing_trace_allows_source_review_without_synthetic_transcript(self):
        plan = self._plan()
        plan["artifacts"]["selected_trace"] = {"status": "MISSING", "reason": "fixture trace was not extracted"}
        self._write_json("trace-check.json", {"trace_request_status": "CANNOT VERIFY"})
        plan["artifacts"]["trace_request"] = self._record("trace-check.json")
        self._write_json("plan.json", plan)
        out = self._prepare()
        with mock.patch("src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke",
                        return_value=AIMessage(content=self._response())):
            result = review.run(self.plan_path, out)
        self.assertEqual(result["trace_identity_status"], "CANNOT VERIFY")

    def test_private_trace_requires_explicit_prompt_derivative(self):
        trace = json.loads(self.files["selected-trace.json"].read_text())
        trace["entry"]["result"][0][0] += "\nPath: /Entropy/private/input.csv"
        self._write_json("selected-trace.json", trace)
        self._write_json("plan.json", self._plan())
        with self.assertRaises(ValueError):
            self._prepare()

    def test_unsafe_pickle_and_unknown_clock_are_rejected(self):
        self._write("selected-trace.json", b"not-json")
        self._write_json("plan.json", self._plan())
        with self.assertRaises(ValueError):
            review.validate_plan(self.plan_path)

    def test_pending_source_design_and_oversized_request_are_rejected(self):
        plan = self._plan()
        plan["artifacts"]["source_design"]["status"] = "PENDING"
        self._write_json("plan.json", plan)
        with self.assertRaises(ValueError):
            review.validate_plan(self.plan_path)

    def test_analysis_identity_contradiction_is_retained_for_review(self):
        self._write("analysis.csv", "gse_id,clock,case_ids,control_ids,comparison\nGSE999,ClockA,GSM1,GSM2,contradictory fixture\n")
        self._write_json("plan.json", self._plan())
        bundle = review.validate_plan(self.plan_path)
        self.assertEqual(bundle["validated"]["analysis_identity"]["status"], "FAIL")

    def test_empty_analysis_artifact_is_preflight_blocked_not_silently_repaired(self):
        self._write("analysis.csv", "gse_id,clock,case_ids,control_ids,comparison\n")
        self._write_json("plan.json", self._plan())
        bundle = review.validate_plan(self.plan_path)
        self.assertEqual(bundle["validated"]["analysis_preflight"]["status"], "CANNOT VERIFY")


if __name__ == "__main__":
    unittest.main()
