"""Software fixtures for update orchestration; not empirical scientific results."""
import copy
import json
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from src import update_pipeline_2026_09_20 as update
from audit_system.tests import require_development_artifact


class UpdatePipelineTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.calls = []

    def declaration(self, path):
        return {"path": str(path), "sha256": update.sha(path)}

    def study(self, study="GSE1", sample="GSM1"):
        source = self.root / study; source.mkdir(exist_ok=True)
        series = source / "series.soft"
        series.write_text(f"^SERIES = {study}\n!Series_status = Public on Jan 03 2025\n!Series_sample_id = {sample}\n")
        samples = source / "samples.soft"
        samples.write_text(f"^SAMPLE = {sample}\n!Sample_title = Software fixture\n!Sample_organism_ch1 = Homo sapiens\n")
        matrix = source / "counts.csv"
        matrix.write_text(f"gene,{sample}\nfixture,1\n")
        return {"study_id": study, "modality": "human_rna", "sources": {
            "series": self.declaration(series), "samples": self.declaration(samples), "matrix": self.declaration(matrix)}}

    def plan(self, studies):
        return {"schema": update.SCHEMA, "studies": studies,
                "limits": {"max_studies": 5, "max_input_bytes": 1000000, "max_samples_per_study": 20,
                           "max_child_seconds": 10, "concurrency": 1, "allow_network": False,
                           "max_total_model_calls": 20, "max_total_cost_usd": 3.0},
                "pipeline": {"max_calls": 8, "max_cost_usd": 1.0},
                "comprehensive_audit": {"model": "gpt-5.6-luna", "max_calls": 2, "max_cost_usd": 0.5,
                                        "input_usd_per_million": 0.25, "output_usd_per_million": 1.2}}

    def execute(self, plan, until="input", runner=None):
        name = str(len(list(self.root.glob("plan-*.json"))))
        path = self.root / ("plan-" + name + ".json")
        path.write_text(json.dumps(plan))
        return update.run(path, self.root / "registry", until, executor=runner or self.fake_executor)

    def fake_executor(self, request_path, log, timeout):
        request = json.loads(request_path.read_text())
        self.calls.append(request)
        log.write_text("SOFTWARE FIXTURE: no pipeline, network, clock or model execution\n")
        if request["phase"] == "pipeline":
            pipeline = Path(request["pipeline_root"]); pipeline.mkdir(exist_ok=True)
            (pipeline / (request["until"] + ".txt")).write_text("SOFTWARE FIXTURE\n")
        update.write_new(Path(request["result_path"]), {"status": "SOFTWARE_FIXTURE", "execution_status": "COMPLETED"})
        return 0

    def test_duplicate_resume_skips_execution_and_preserves_quarantine(self):
        plan = self.plan([self.study()])
        first = self.execute(plan); second = self.execute(plan)
        self.assertEqual(len(self.calls), 1)
        self.assertEqual(first["studies"][0]["version_id"], second["studies"][0]["version_id"])
        self.assertEqual(second["studies"][0]["version_event"], "EXISTING_VERSION")
        self.assertEqual(second["reserved_model_calls"], 0)
        self.assertEqual(second["studies"][0]["first_public_date"], "2025-01-03")
        self.assertEqual(second["original_paper_2025_eligibility"], "QUARANTINED")
        self.assertEqual(second["scientific_status"], "CANNOT VERIFY")

    def test_changed_input_requires_explicit_version_and_preserves_original(self):
        entry = self.study(); plan = self.plan([entry])
        first = self.execute(plan); previous = first["studies"][0]["version_id"]
        changed = self.root / "GSE1/changed-counts.csv"; changed.write_text("gene,GSM1\nfixture,2\n")
        entry["sources"]["matrix"] = self.declaration(changed)
        blocked = self.execute(plan)
        self.assertEqual(blocked["studies"][0]["status"], "UNRESOLVED")
        self.assertIn("supersedes_version", blocked["studies"][0]["error"])
        self.assertEqual(len(self.calls), 1)
        entry["supersedes_version"] = previous
        resumed = self.execute(plan)
        self.assertEqual(resumed["studies"][0]["version_event"], "NEW_VERSION")
        self.assertNotEqual(resumed["studies"][0]["version_id"], previous)
        self.assertEqual(len(self.calls), 2)
        self.assertTrue((self.root / "registry/studies/human_rna/GSE1" / previous / "version.json").exists())

    def test_source_snapshot_is_preserved_and_tampering_stops_resume(self):
        plan = self.plan([self.study()])
        first = self.execute(plan)
        version_root = Path(first["studies"][0]["version_root"])
        version = json.loads((version_root / "version.json").read_text())
        for name, digest in version["implementation"].items():
            self.assertEqual(update.sha(version_root / "implementation" / name), digest)
        changed = version_root / "implementation/src/update_pipeline_2026_09_20.py"
        changed.write_text("Software fixture: snapshot changed\n")
        resumed = self.execute(plan)
        self.assertEqual(resumed["studies"][0]["status"], "UNRESOLVED")
        self.assertIn("Frozen implementation snapshot changed", resumed["studies"][0]["error"])
        self.assertEqual(len(self.calls), 1)

    def test_failure_retained_and_other_study_continues_without_retry(self):
        plan = self.plan([self.study(), self.study("GSE2", "GSM2")])
        def runner(request_path, log, timeout):
            if "GSE1" in str(request_path):
                self.calls.append("failed"); log.write_text("Software fixture failure\n"); return 1
            return self.fake_executor(request_path, log, timeout)
        first = self.execute(plan, runner=runner)
        self.assertEqual([r["status"] for r in first["studies"]], ["FAILED", "PIPELINE_COMPLETED_THROUGH_INPUT"])
        self.execute(plan, runner=runner)
        self.assertEqual(len(self.calls), 2)
        progressed = self.execute(self.plan([plan["studies"][0]]), until="report", runner=runner)
        self.assertIn("Earlier child failed", progressed["studies"][0]["error"])
        self.assertEqual(len(self.calls), 2)

    def test_input_can_progress_to_report_without_duplicate_budget_reservation(self):
        plan = self.plan([self.study()])
        self.execute(plan)
        result = self.execute(plan, until="report")
        self.assertEqual(len(self.calls), 3)
        self.assertEqual(result["reserved_model_calls"], 10)
        self.assertEqual(result["reserved_cost_usd"], 1.5)
        self.execute(plan, until="report")
        self.assertEqual(len(self.calls), 3)

    def test_unsupported_comprehensive_model_stops_before_producer_or_reservation(self):
        plan=self.plan([self.study()])
        plan['comprehensive_audit']['model']='unsupported-model-software-fixture'
        with patch.object(update,'reserve',side_effect=AssertionError('No paid reservation allowed')):
            with self.assertRaisesRegex(ValueError,'unsupported model rejected before execution'):
                self.execute(plan,until='report')
        self.assertEqual(self.calls,[])
        self.assertFalse((self.root/'registry').exists())

    def test_explicit_gpt6_switch_preserves_versions_and_lifetime_reservations(self):
        plan=self.plan([self.study()])
        first=self.execute(plan,until="report")
        prior=Path(first["studies"][0]["version_root"])
        prior_bytes=(prior/"version.json").read_bytes()
        old_id=first["studies"][0]["version_id"]
        plan["pipeline"].update(model="gpt-6-luna",audit_model="gpt-6-luna",
            input_usd_per_million=.20,output_usd_per_million=.75,
            audit_input_usd_per_million=.20,audit_output_usd_per_million=.75)
        plan["comprehensive_audit"].update(model="gpt-6-luna",
            input_usd_per_million=.20,output_usd_per_million=.75)
        rejected=self.execute(plan,until="report")
        self.assertIn("supersedes_version",rejected["studies"][0]["error"])
        self.assertEqual(len(self.calls),2)
        self.assertEqual(rejected["reserved_model_calls"],10)
        plan["studies"][0]["supersedes_version"]=old_id
        revised=self.execute(plan,until="report")
        new_root=Path(revised["studies"][0]["version_root"])
        version=update.read_sealed(new_root/"version.json")
        self.assertNotEqual(version["version_id"],old_id)
        self.assertEqual(version["supersedes_version"],old_id)
        self.assertEqual(version["pipeline_options"]["model"],"gpt-6-luna")
        self.assertEqual(version["pipeline_options"]["audit_model"],"gpt-6-luna")
        self.assertEqual(version["comprehensive_audit"],plan["comprehensive_audit"])
        self.assertEqual(self.calls[-1]["audit_request"]["budget"],plan["comprehensive_audit"])
        self.assertEqual(revised["reserved_model_calls"],20)
        self.assertEqual(revised["reserved_cost_usd"],3.0)
        self.execute(plan,until="report")
        self.assertEqual(len(self.calls),4)
        self.assertEqual((prior/"version.json").read_bytes(),prior_bytes)
        self.assertEqual(update.DEFAULT_OPTIONS["model"],"gpt-5.6-luna")
        self.assertEqual(update.DEFAULT_OPTIONS["audit_model"],"gpt-5.6-luna")

    def test_gpt6_requires_role_tariffs_before_any_execution(self):
        for role in ("model","audit_model"):
            with self.subTest(role=role):
                plan=self.plan([self.study()])
                plan["pipeline"][role]="gpt-6-luna"
                with self.assertRaisesRegex(ValueError,"both explicit role-specific tariffs"):
                    self.execute(plan,until="report")
                self.assertEqual(self.calls,[])
                self.assertFalse((self.root/"registry").exists())

    def test_optional_service_tier_changes_version_and_forwards_both_children(self):
        plan=self.plan([self.study()])
        original=self.execute(plan,until='report')
        old_root=Path(original['studies'][0]['version_root'])
        old_version=(old_root/'version.json').read_bytes()
        prior=update.read_sealed(old_root/'version.json')
        self.assertNotIn('service_tier',update.DEFAULT_OPTIONS)
        self.assertNotIn('service_tier',prior['pipeline_options'])
        self.assertNotIn('service_tier',prior['comprehensive_audit'])
        plan['pipeline']['service_tier']='default'
        plan['comprehensive_audit']['service_tier']='default'
        blocked=self.execute(plan,until='report')
        self.assertIn('supersedes_version',blocked['studies'][0]['error'])
        self.assertEqual(len(self.calls),2)
        plan['studies'][0]['supersedes_version']=prior['version_id']
        current=self.execute(plan,until='report')
        version=update.read_sealed(Path(current['studies'][0]['version_root'])/'version.json')
        self.assertNotEqual(version['version_id'],prior['version_id'])
        for record in (version,version['identity']):
            self.assertEqual(record['pipeline_options']['service_tier'],'default')
            self.assertEqual(record['comprehensive_audit']['service_tier'],'default')
        self.assertEqual(self.calls[-1]['audit_request']['budget']['service_tier'],'default')
        request=copy.deepcopy(self.calls[-2]);request['result_path']=str(self.root/'tier-child-fixture.json')
        with patch('src.pipeline_2026_09_10.main',return_value={'status':'SOFTWARE_FIXTURE'}) as child:
            update.child(request)
        argv=child.call_args.args[0]
        self.assertEqual(argv[argv.index('--service-tier')+1],'default')
        self.assertEqual(current['reserved_model_calls'],20)
        self.assertEqual(current['reserved_cost_usd'],3)
        self.assertEqual((old_root/'version.json').read_bytes(),old_version)

    def test_invalid_service_tier_stops_before_controller_or_budget(self):
        for role in ('pipeline','comprehensive_audit'):
            for value in (None,False,1,[],{},'', 'auto','priority'):
                plan=self.plan([self.study()]);plan[role]['service_tier']=value
                with self.subTest(role=role,value=value), \
                        patch.object(update,'reserve',side_effect=AssertionError('No reservation')):
                    with self.assertRaisesRegex(ValueError,'service_tier'):
                        self.execute(plan,until='report')
                self.assertEqual(self.calls,[])
                self.assertFalse((self.root/'registry').exists())

    def test_comprehensive_pair_policy_is_separate_from_budget_and_bound_to_version(self):
        plan=self.plan([self.study()])
        self.assertNotIn('comprehensive_audit_policy',plan)
        self.assertEqual(len(update.validate_plan(plan)[2]),5)
        plan['comprehensive_audit_policy']='source_label_pairs_v1'
        first=self.execute(plan,until='report')
        version_root=Path(first['studies'][0]['version_root'])
        version=update.read_sealed(version_root/'version.json')
        self.assertEqual(version['comprehensive_audit_policy'],'source_label_pairs_v1')
        self.assertEqual(version['identity']['comprehensive_audit_policy'],'source_label_pairs_v1')
        audit=next(call for call in self.calls if call['phase']=='comprehensive_audit')
        request=audit['audit_request']
        self.assertEqual(request['comprehensive_audit_policy'],'source_label_pairs_v1')
        self.assertEqual(request['budget'],plan['comprehensive_audit'])
        self.assertEqual(len(request['budget']),5)
        self.execute(plan,until='report')
        self.assertEqual(len(self.calls),2)
        plan['comprehensive_audit_policy']='model_selected_policy'
        with self.assertRaisesRegex(ValueError,'controller-selected comprehensive'):
            update.validate_plan(plan)
        del plan['comprehensive_audit_policy']
        plan['pipeline']['source_design_policy']='source_label_pairs_v1'
        with self.assertRaisesRegex(ValueError,'Unsupported child pipeline option'):
            update.validate_plan(plan)

    def test_comprehensive_pair_policy_change_requires_new_version(self):
        plan=self.plan([self.study()])
        first=self.execute(plan)
        old_id=first['studies'][0]['version_id']
        plan['comprehensive_audit_policy']='source_label_pairs_v1'
        stopped=self.execute(plan)
        self.assertEqual(stopped['studies'][0]['status'],'UNRESOLVED')
        self.assertIn('supersedes_version',stopped['studies'][0]['error'])
        self.assertEqual(len(self.calls),1)
        plan['studies'][0]['supersedes_version']=old_id
        revised=self.execute(plan)
        self.assertNotEqual(revised['studies'][0]['version_id'],old_id)
        self.assertEqual(revised['studies'][0]['version_event'],'NEW_VERSION')
        self.assertEqual(len(self.calls),2)

    def test_aggregate_pipeline_and_post_audit_budget_is_enforced(self):
        plan = self.plan([self.study(), self.study("GSE2", "GSM2")])
        plan["limits"].update(max_total_model_calls=10, max_total_cost_usd=1.5)
        result = self.execute(plan, until="report")
        self.assertEqual(len(self.calls), 2)
        self.assertEqual(result["reserved_model_calls"], 10)
        self.assertEqual(result["reserved_cost_usd"], 1.5)
        self.assertEqual(result["studies"][1]["status"], "UNRESOLVED")
        self.assertIn("Aggregate model-call", result["studies"][1]["error"])

    def test_unresolved_comprehensive_audit_stays_unresolved_and_is_not_retried(self):
        plan = self.plan([self.study()])
        def runner(request_path, log, timeout):
            request = json.loads(request_path.read_text())
            if request["phase"] == "pipeline": return self.fake_executor(request_path, log, timeout)
            self.calls.append(request)
            log.write_text("SOFTWARE FIXTURE: unresolved audit\n")
            update.write_new(Path(request["result_path"]), {"status": "UNRESOLVED", "execution_status": "UNRESOLVED"})
            return 0
        first = self.execute(plan, until="report", runner=runner)
        self.assertEqual(first["studies"][0]["status"], "PENDING_COMPREHENSIVE_AUDIT")
        self.assertEqual(first["studies"][0]["comprehensive_audit"]["status"], "UNRESOLVED")
        self.assertEqual(first["status"], "COMPLETED_WITH_UNRESOLVED_CASES")
        self.execute(plan, until="report", runner=runner)
        self.assertEqual(len(self.calls), 2)

    def test_comprehensive_dependencies_are_bound_and_cannot_escape_repository(self):
        hook = self.root / "audit_system/corpus_audit_2026_09_20.py"
        hook.parent.mkdir()
        hook.write_text("CODE_FILES = ('verifier.py',)\nraise RuntimeError('Do not import this fixture')\n")
        dependency = self.root / "verifier.py"; dependency.write_text("fixture = 1\n")
        with patch.object(update, "ROOT", self.root), patch.object(update, "IMPLEMENTATIONS", ()):
            before = update.implementation_hashes()
            self.assertEqual(before["verifier.py"], update.sha(dependency))
            dependency.write_text("fixture = 2\n")
            with self.assertRaisesRegex(update.UpdateIntegrityError, "implementation changed"):
                update.verify_implementation(before)
            self.assertNotEqual(before, update.implementation_hashes())
            hook.write_text("CODE_FILES = ('../escape.py',)\n")
            with self.assertRaisesRegex(update.UpdateIntegrityError, "escapes repository"):
                update.implementation_hashes()

    def test_source_and_output_tampering_are_not_resumed(self):
        plan = self.plan([self.study()])
        result = self.execute(plan)
        pipeline = Path(result["studies"][0]["version_root"]) / "pipeline"
        (pipeline / "input.txt").write_text("changed")
        resumed = self.execute(plan)
        self.assertIn("Completed child output changed", resumed["studies"][0]["error"])
        self.assertEqual(len(self.calls), 1)
        Path(plan["studies"][0]["sources"]["matrix"]["path"]).write_text("changed source")
        changed = self.execute(plan)
        self.assertIn("Source SHA256 differs", changed["studies"][0]["error"])
        self.assertEqual(len(self.calls), 1)

    def test_discovery_snapshot_is_bound_and_duplicate_studies_rejected(self):
        study = self.study()
        discovery = self.root / "discovered.json"
        discovery.write_text(json.dumps({"schema": "clockbase-update-candidates-v1", "studies": [study]}))
        plan = self.plan([]); plan["discovery"] = self.declaration(discovery)
        result = self.execute(plan)
        self.assertEqual(len(result["studies"]), 1)
        plan["studies"] = [study]
        with self.assertRaisesRegex(ValueError, "Duplicate study/modality"):
            self.execute(plan)

    def test_series_dates_and_membership_require_source_support(self):
        study = self.study()
        series = Path(study["sources"]["series"]["path"])
        series.write_text("^SERIES = GSE1\n!Series_last_update_date = Jan 03 2025\n!Series_sample_id = GSM1\n")
        study["sources"]["series"] = self.declaration(series)
        result = self.execute(self.plan([study]))
        self.assertEqual(result["studies"][0]["status"], "UNRESOLVED")
        self.assertEqual(len(self.calls), 0)

    def test_missing_sample_metadata_and_wrong_organism_stop_before_execution(self):
        study = self.study()
        sample_path = Path(study["sources"]["samples"]["path"])
        sample_path.write_text("^SAMPLE = GSM2\n!Sample_organism_ch1 = Homo sapiens\n")
        study["sources"]["samples"] = self.declaration(sample_path)
        result = self.execute(self.plan([study]))
        self.assertIn("metadata disagree", result["studies"][0]["error"])
        sample_path.write_text("^SAMPLE = GSM1\n!Sample_organism_ch1 = Mus musculus\n")
        study["sources"]["samples"] = self.declaration(sample_path)
        result = self.execute(self.plan([study]))
        self.assertIn("organism differs", result["studies"][0]["error"])
        self.assertEqual(len(self.calls), 0)

    def test_byte_limit_network_off_and_lfs_pointer(self):
        entry = self.study(); plan = self.plan([entry]); plan["limits"]["max_input_bytes"] = 1
        result = self.execute(plan)
        self.assertIn("byte cap", result["studies"][0]["error"])
        plan["limits"]["max_input_bytes"] = 1000000
        entry["sources"]["matrix"] = {"url": "https://ftp.ncbi.nlm.nih.gov/test.csv", "max_bytes": 100}
        result = self.execute(plan)
        self.assertIn("Network acquisition disabled", result["studies"][0]["error"])
        path = self.root / "pointer.csv"
        path.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 123\n")
        entry["sources"]["matrix"] = self.declaration(path)
        result = self.execute(plan)
        self.assertIn("unmaterialized LFS", result["studies"][0]["error"])
        self.assertEqual(len(self.calls), 0)

    def test_frozen_metadata_uses_source_parser_without_network(self):
        entry = self.study(); limits = self.plan([entry])["limits"]
        sources = update.collect_sources(entry, self.root / "acq", limits, {})
        output = self.root / "metadata"; output.mkdir()
        result = update.frozen_metadata("GSE1", output, sources)
        self.assertEqual(result["sample_ids"], ["GSM1"])
        self.assertEqual(result["series_first_public_date"], "2025-01-03")
        context = json.loads(Path(result["study_context"]).read_text())
        self.assertEqual(context["acquisition"], [{"step": "publication-methods", "status": "NOT_REQUESTED"}])

    def discovery_fixture(self):
        seen = self.root / "seen.jsonl"
        seen.write_text(json.dumps({"study_id": "GSE1"}) + "\n")
        sample_sources = self.root / "seen-samples.json"
        sample_sources.write_text(json.dumps([{"accession": "GSE1", "sample_ids": ["GSM1"]}]))
        plan = {"schema": "clockbase-update-discovery-plan-v1", "query": "senescence", "modality": "human_rna",
                "date_from": "2025-01-01", "date_to": "2025-12-31",
                "seen_worklist": self.declaration(seen), "seen_sample_sources": self.declaration(sample_sources),
                "limits": {"max_requests": 10, "max_response_bytes": 100000, "max_total_response_bytes": 1000000,
                           "max_seconds": 30, "max_query_ids": 3, "max_candidates": 3, "minimum_interval_seconds": 0.35}}
        path = self.root / "discovery-plan.json"; path.write_text(json.dumps(plan))
        def transport(url, cap):
            if "esearch" in url: return json.dumps({"esearchresult": {"idlist": ["1", "2", "3"], "count": "4"}}).encode()
            if "esummary" in url:
                requested = update.urllib.parse.parse_qs(update.urllib.parse.urlparse(url).query)["id"][0].split(",")
                result = {uid: {"accession": "GSE" + str(int(uid) + 10), "pdat": "2025/01/01"} for uid in requested}
                result["uids"] = requested
                return json.dumps({"result": result}).encode()
            if "GSE11" in url: study, sample, date = "GSE11", "GSM1", "Jan 03 2025"
            elif "GSE12" in url: study, sample, date = "GSE12", "GSM2", "Jan 03 2026"
            else: study, sample, date = "GSE13", "GSM3", "Jan 03 2025"
            return f"^SERIES = {study}\n!Series_status = Public on {date}\n!Series_sample_id = {sample}\n".encode()
        return path, transport

    def test_discovery_retains_every_uid_and_distinguishes_source_dates_overlap_and_missing_assay(self):
        plan, transport = self.discovery_fixture()
        output = self.root / "discovery"
        with patch.object(update.time, "sleep"):
            result = update.discover(plan, output, transport=transport)
        self.assertEqual(result["processed_uids"], 3)
        self.assertEqual(result["unretrieved_uid_count"], 1)
        self.assertFalse(result["query_fully_enumerated"])
        self.assertEqual(result["candidate_status_counts"], {"SHARED_SEEN_SAMPLES": 1, "OUTSIDE_FIRST_PUBLIC_WINDOW": 1, "UNRESOLVED_NO_ASSAY_INPUT": 1})
        self.assertEqual(result["paid_model_calls"], 0)
        self.assertEqual(result["molecular_downloads"], 0)
        def no_network(*args): raise AssertionError("Discovery resume made a network request")
        resumed = update.discover(plan, output, transport=no_network)
        self.assertEqual(result, resumed)

    def test_discovery_request_budget_leaves_unresolved_candidates_visible(self):
        path, transport = self.discovery_fixture()
        plan = json.loads(path.read_text()); plan["limits"]["max_requests"] = 3; path.write_text(json.dumps(plan))
        with patch.object(update.time, "sleep"):
            result = update.discover(path, self.root / "discovery", transport=transport)
        self.assertEqual(result["requests_used"], 3)
        self.assertEqual(sum(result["candidate_status_counts"].values()), 3)
        self.assertEqual(result["candidate_status_counts"]["UNRESOLVED"], 2)

    def test_discovery_resume_rejects_changed_archived_response(self):
        plan, transport = self.discovery_fixture(); output = self.root / "discovery"
        with patch.object(update.time, "sleep"):
            update.discover(plan, output, transport=transport)
        (output / "requests/esearch/response.bin").write_text("changed")
        with self.assertRaisesRegex(update.UpdateIntegrityError, "Completed child output changed"):
            update.discover(plan, output, transport=transport)

    def test_discovery_leaves_unprocessed_enumerated_candidates(self):
        path, transport = self.discovery_fixture()
        plan = json.loads(path.read_text()); plan["limits"]["max_candidates"] = 1; path.write_text(json.dumps(plan))
        with patch.object(update.time, "sleep"):
            result = update.discover(path, self.root / "discovery", transport=transport)
        self.assertEqual(result["unprocessed_enumerated_uids"], 2)
        self.assertEqual(result["candidate_status_counts"]["NOT_YET_PROCESSED"], 2)

    def test_discovery_continues_frozen_uid_snapshot_without_reissuing_search(self):
        path, transport = self.discovery_fixture()
        plan = json.loads(path.read_text()); plan["limits"]["max_candidates"] = 1; path.write_text(json.dumps(plan))
        first = self.root / "first-discovery"
        with patch.object(update.time, "sleep"):
            update.discover(path, first, transport=transport)
        plan["query_snapshot"] = self.declaration(first / "summary.json")
        plan["candidate_offset"] = 1
        path.write_text(json.dumps(plan))
        def no_search(url, cap):
            if "esearch" in url: raise AssertionError("Repeated query may shift the candidate order")
            return transport(url, cap)
        with patch.object(update.time, "sleep"):
            result = update.discover(path, self.root / "next-discovery", transport=no_search)
        self.assertEqual(result["enumerated_uids"], 3)
        self.assertEqual(result["requests_used"], 2)
        self.assertEqual(result["candidate_status_counts"]["OUTSIDE_FIRST_PUBLIC_WINDOW"], 1)
        self.assertEqual(result["query_snapshot_source"]["summary_sha256"], update.sha(first / "summary.json"))

    def test_dnam_reference_must_be_declared_before_paid_stage_or_preprocessing(self):
        entry = self.study(); entry["modality"] = "human_dnam"
        result = self.execute(self.plan([entry]), until="report")
        self.assertIn("pinned DNAm reference", result["studies"][0]["error"])
        self.assertEqual(len(self.calls), 0)

    def test_coding_format_retry_is_opt_in_bounded_and_changes_version_identity(self):
        default_plan = self.plan([self.study()])
        _, defaults, _ = update.validate_plan(default_plan)
        self.assertEqual(defaults["coding_validation_retries"], 0)
        retry_plan = copy.deepcopy(default_plan)
        retry_plan["pipeline"]["coding_validation_retries"] = 1
        _, retry, _ = update.validate_plan(retry_plan)
        self.assertEqual(retry["coding_validation_retries"], 1)
        self.assertNotEqual(update.identity(defaults), update.identity(retry))
        retry_plan["pipeline"]["coding_validation_retries"] = 2
        with self.assertRaisesRegex(ValueError, "coding_validation_retries must be 0 or 1"):
            update.validate_plan(retry_plan)

    def test_coding_completion_policy_preserves_default_options_and_is_explicit(self):
        plan = self.plan([self.study()])
        _, defaults, _ = update.validate_plan(plan)
        self.assertNotIn("coding_completion_policy", defaults)
        plan["pipeline"]["coding_completion_policy"] = "legacy"
        _, legacy, _ = update.validate_plan(plan)
        self.assertEqual(defaults, legacy)
        plan["pipeline"]["coding_completion_policy"] = "outputs_v1"
        _, enabled, _ = update.validate_plan(plan)
        self.assertEqual(enabled["coding_completion_policy"], "outputs_v1")
        self.assertEqual(enabled["coding_validation_retries"], 0)
        self.assertNotEqual(update.identity(defaults), update.identity(enabled))
        for invalid in (None, True, 1, [], "automatically_retry"):
            with self.subTest(invalid=invalid):
                plan["pipeline"]["coding_completion_policy"] = invalid
                with self.assertRaisesRegex(ValueError, "Unknown coding completion policy"):
                    update.validate_plan(plan)

    def test_coding_completion_policy_reaches_the_existing_child_cli(self):
        plan = self.plan([self.study()])
        plan["pipeline"].update(coding_completion_policy="outputs_v1", coding_validation_retries=1)
        self.execute(plan)
        request = copy.deepcopy(self.calls[-1])
        request["result_path"] = str(self.root / "child-cli-result.json")
        with patch("src.pipeline_2026_09_10.main", return_value={"status": "SOFTWARE_FIXTURE"}) as main:
            update.child(request)
        argv = main.call_args.args[0]
        self.assertEqual(argv[argv.index("--coding-completion-policy") + 1], "outputs_v1")
        self.assertEqual(argv[argv.index("--coding-validation-retries") + 1], "1")


class ConservativeSourceContractTests(unittest.TestCase):
    """Archived-response regression and mutated software fixtures; no API calls."""
    def setUp(self):
        from src import pipeline_2026_09_10 as pipeline
        self.pipeline = pipeline
        base = update.ROOT / 'paper/revision_2026-09-07/evidence/continuous-audit-pipeline-2026-09-20/update-controller'
        prior = base / 'registry-003/studies/human_rna/GSE127814/81c7b90e1f432ae1b2887bd8b69ee1296383ad8a95502c41bfa26bd63092212f/pipeline'
        require_development_artifact(prior / '03_blinded_design/attempt-000/source-packet.json')
        self.packet = json.loads((prior / '03_blinded_design/attempt-000/source-packet.json').read_text())
        self.responses = []
        for number in (0, 1):
            directory = prior / f'api/design-{number:03d}'
            receipt = json.loads((directory / 'receipt.json').read_text())
            self.assertEqual(update.sha(directory / 'response.json'), receipt['response_sha256'])
            self.responses.append(json.loads((directory / 'response.json').read_text())['choices'][0]['message']['content'])

    def parse_conservative(self, raw):
        candidate, changes = self.pipeline.conservative_source_status(raw)
        checked = self.pipeline.parse_design(json.dumps(candidate), self.packet, review_profile='source_domains_v4')
        return candidate, changes, checked

    def test_strict_default_preserves_prior_rejection(self):
        args = self.pipeline.parser().parse_args(['--study', 'GSE127814', '--modality', 'human_rna', '--output', 'unused'])
        self.assertEqual(args.source_contract_policy, 'strict')
        self.assertEqual(update.DEFAULT_OPTIONS['source_contract_policy'], 'strict')
        with self.assertRaisesRegex(ValueError, 'Domain PASS conflicts'):
            self.pipeline.parse_design(self.responses[1], self.packet, review_profile='source_domains_v4')

    def test_json_response_mode_is_explicit_and_only_applies_to_source_design_actor(self):
        args = self.pipeline.parser().parse_args(['--study', 'GSE127814', '--modality', 'human_rna', '--output', 'unused',
            '--input-usd-per-million', '0.5', '--output-usd-per-million', '1.8'])
        with TemporaryDirectory() as folder:
            root = Path(folder)
            self.assertFalse(self.pipeline.build_client(args, root, 'design').json_response)
            args.source_contract_policy = 'conservative_status_v1'
            self.assertTrue(self.pipeline.build_client(args, root, 'design').json_response)
            self.assertFalse(self.pipeline.build_client(args, root, 'primary').json_response)

    def test_unsupported_comparison_status_is_rejected_with_precise_feedback(self):
        from langchain_core.messages import AIMessage
        from types import SimpleNamespace
        base = update.ROOT / 'paper/revision_2026-09-07/evidence/continuous-audit-pipeline-2026-09-20/update-controller'
        path = base / 'registry-005/studies/human_rna/GSE127814/b2030a53e4342b273303a6bdb97885f25a61dd19a34450d4ad11f7b40b221bce/pipeline/api/design-001/response.json'
        raw = json.loads(path.read_text())['choices'][0]['message']['content']
        candidate, changes = self.pipeline.conservative_source_status(raw)
        self.assertEqual(candidate, json.loads(raw))
        self.assertEqual(changes, [])
        with self.assertRaisesRegex(ValueError, 'Model cannot approve'):
            self.parse_conservative(raw)
        captured = []
        def invoke(messages):
            captured.append(copy.deepcopy(messages))
            return AIMessage(content=raw)
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            with self.assertRaisesRegex(ValueError, 'Model cannot approve'):
                self.pipeline.assess_design(self.packet, SimpleNamespace(invoke=invoke), directory,
                    review_profile='source_domains_v4', validation_retries=1,
                    source_binding_profile='literal_source_v1', source_contract_policy='conservative_status_v1')
            self.assertEqual(len(captured), 2)
            self.assertIn(self.pipeline.SOURCE_STATUS_CONTRACT, captured[0][0].content)
            feedback = json.loads((directory / 'response-validation/000-validation.json').read_text())
            self.assertEqual(len(feedback['status_issues']), 4)
            self.assertEqual(feedback['status_issues'][0], {'path': '$.comparisons[0].status',
                'received': 'SOURCE_DEFINED', 'allowed': ['DESIGN_PROPOSED', 'CANNOT VERIFY']})
            self.assertIn('Unsupported source-design status enum', feedback['error'])
            self.assertEqual((directory / 'contract-adjustments/001-original-response.txt').read_text(), raw)

    def test_only_status_changes_and_every_original_limitation_group_and_citation_remains(self):
        original = json.loads(self.responses[1])
        candidate, changes, checked = self.parse_conservative(self.responses[1])
        self.assertEqual(len(changes), 1)
        self.assertEqual(changes[0]['path'], '$.audit_domains.grouping.status')
        self.assertEqual(candidate['audit_domains']['grouping']['status'], 'CANNOT VERIFY')
        restored = copy.deepcopy(candidate)
        restored['audit_domains']['grouping']['status'] = 'PASS'
        self.assertEqual(restored, original)
        self.assertEqual(checked['scientific_status'], 'CANNOT VERIFY')
        self.assertTrue(all(len(c['case_ids']) == len(c['control_ids']) == 1 for c in checked['comparisons']))
        self.assertTrue(all(c['proposed_action'] == 'withhold_inference' and c['unresolved'] for c in checked['comparisons']))
        unchanged, next_changes, _ = self.parse_conservative(json.dumps(candidate))
        self.assertEqual(unchanged, candidate)
        self.assertEqual(next_changes, [])

    def test_invalid_or_missing_quotations_and_duplicate_keys_still_fail(self):
        with self.assertRaisesRegex(ValueError, 'Unsupported source quotation'):
            self.parse_conservative(self.responses[0])
        missing = json.loads(self.responses[1]); missing['audit_domains']['grouping']['evidence'] = []
        with self.assertRaisesRegex(ValueError, 'lacks source evidence'):
            self.parse_conservative(json.dumps(missing))
        duplicate = self.responses[1].replace('"schema":', '"schema":"invalid","schema":', 1)
        with self.assertRaisesRegex(ValueError, 'Duplicate model response key'):
            self.parse_conservative(duplicate)

    def test_invented_or_overlapping_ids_and_required_fields_still_fail(self):
        original = json.loads(self.responses[1])
        invalid = copy.deepcopy(original); invalid['comparisons'][0]['case_ids'] = ['GSM999999999']
        overlap = copy.deepcopy(original); overlap['comparisons'][0]['case_ids'] = overlap['comparisons'][0]['control_ids']
        missing = copy.deepcopy(original); del missing['audit_domains']['grouping']['verification_needed']
        for candidate in (invalid, overlap, missing):
            with self.subTest(candidate=candidate), self.assertRaises(ValueError):
                self.parse_conservative(json.dumps(candidate))

    def test_assessment_archives_original_and_partial_admissibility_without_scientific_acceptance(self):
        from langchain_core.messages import AIMessage
        from types import SimpleNamespace
        client = SimpleNamespace(invoke=lambda unused: AIMessage(content=self.responses[1]))
        with TemporaryDirectory() as folder:
            directory = Path(folder)
            result = self.pipeline.assess_design(self.packet, client, directory, review_profile='source_domains_v4',
                source_binding_profile='literal_source_v1', source_contract_policy='conservative_status_v1')
            self.assertEqual((directory / 'contract-adjustments/000-original-response.txt').read_text(), self.responses[1])
            self.assertEqual(json.loads((directory / 'model-design.json').read_text()), json.loads(self.responses[1]))
            receipt = json.loads((directory / 'source-contract-admissibility.json').read_text())
            self.assertEqual(receipt['contract_status'], 'PARTIALLY_ADMISSIBLE')
            self.assertFalse(receipt['accepted_result'])
            self.assertEqual(receipt['scientific_status'], 'CANNOT VERIFY')
            self.assertEqual(result['audit_domains']['grouping']['status'], 'CANNOT VERIFY')
            self.assertEqual(receipt['original_response_sha256'], update.sha(directory / 'contract-adjustments/000-original-response.txt'))


class SourceDesignReplayTests(unittest.TestCase):
    """Archived-source regressions and labeled mutation fixtures; no model calls."""
    def setUp(self):
        from src import pipeline_2026_09_10 as pipeline
        self.pipeline = pipeline
        self.temporary = TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.base = update.ROOT / 'paper/revision_2026-09-07/evidence/continuous-audit-pipeline-2026-09-20/update-controller'
        origin = self.base / 'registry-007/studies/human_dnam/GSE247655/2fee535cd9f9932d4905b84c068a7af673fbf22a3a89c6800db47de15f96de68'
        self.origin = self.copy_source_archive(origin)
        self.manifest = self.root / 'replay.json'
        self.freeze()
        self.packet = json.loads((self.origin / 'pipeline/03_blinded_design/attempt-000/source-packet.json').read_text())
        self.parameters = dict(review_profile='source_domains_v4', source_binding_profile='literal_source_v1',
                               source_contract_policy='conservative_status_v1', source_model='gpt-5.6-luna')

    def copy_source_archive(self, origin):
        """Build a path-relocated SOFTWARE FIXTURE from sealed source-only bytes.

        Git preserves the historical absolute locators. They cannot identify
        files in another checkout, so do not ask the production replay binder
        to accept the unmodified relocated archive. Verify its local members,
        then rebase only the private fixture's explicit metadata roles and
        corresponding source packet/transport/request locators. Model replies,
        dates, source bytes, groups and limitations remain unchanged.
        """
        origin = Path(origin).resolve()
        require_development_artifact(origin / 'version.json')
        version = update.read_sealed(origin / 'version.json')
        self.assertEqual(version['version_id'], update.identity(version['identity']))
        names = {'version.json', 'version.json.sha256', 'pipeline/config.json'}
        stages = {}
        for stage_name in ('01_metadata', '03_blinded_design'):
            stage = origin / 'pipeline' / stage_name
            complete = stage / 'complete.json'
            self.assertEqual(update.sha(complete), (stage / 'complete.sha256').read_text().strip())
            record = json.loads(complete.read_text()); stages[stage_name] = record
            self.assertIsInstance(record['outputs'], dict)
            self.assertTrue(record['outputs'])
            names.update(str(path.relative_to(origin)) for path in (complete, stage / 'complete.sha256'))
            for relative, digest in record['outputs'].items():
                self.assertFalse(Path(relative).is_absolute())
                self.assertNotIn('..', Path(relative).parts)
                member = (stage / relative).resolve()
                self.assertTrue(member.is_relative_to(stage))
                self.assertEqual(update.sha(member), digest)
                names.add(str(member.relative_to(origin)))
        slots = sorted((origin / 'pipeline/api').glob('design-*'))
        self.assertIn(len(slots), (1, 2))
        self.assertEqual([slot.name for slot in slots], [f'design-{i:03d}' for i in range(len(slots))])
        for slot in slots:
            names.update(str((slot / name).relative_to(origin)) for name in
                         ('request.json', 'response.json', 'reservation.json', 'receipt.json'))
        target = self.root / ('source-copy-' + str(getattr(self, 'source_copy_count', 0)))
        self.source_copy_count = getattr(self, 'source_copy_count', 0) + 1
        original_hashes = {name: update.sha(origin / name) for name in names}
        for name in sorted(names):
            destination = target / name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(origin / name, destination)
            self.assertEqual(update.sha(destination), original_hashes[name])
        receipt = target / 'pipeline/01_metadata/complete.json'
        metadata = stages['01_metadata']
        roles = {'metadata': 'metadata.csv', 'study_context': 'study-context.json',
                 'source_sample_records': 'samples.soft', 'series_source': 'series.soft'}
        historical_roots, relocations = set(), {}
        for role, basename in roles.items():
            members = [name for name in metadata['outputs'] if Path(name).name == basename]
            self.assertEqual(len(members), 1)
            relative = Path('pipeline/01_metadata') / members[0]
            previous = Path(metadata['result'][role])
            self.assertTrue(previous.is_absolute())
            self.assertEqual(previous.parts[-len(relative.parts):], relative.parts)
            historical_roots.add(str(Path(*previous.parts[:-len(relative.parts)])))
            relocated = target / relative
            relocations[role] = {'original': str(previous), 'fixture': str(relocated),
                                 'sha256': metadata['outputs'][members[0]]}
            metadata['result'][role] = str(relocated)
        self.assertEqual(len(historical_roots), 1)
        receipt.write_text(json.dumps(metadata))
        (receipt.parent / 'complete.sha256').write_text(update.sha(receipt) + '\n')
        stage = target / 'pipeline/03_blinded_design'
        completed = stages['03_blinded_design']
        packets = [stage / name for name in completed['outputs'] if Path(name).name == 'source-packet.json']
        self.assertEqual(len(packets), 1)
        packet_path = packets[0]
        original_packet = json.loads(packet_path.read_text())
        packet = copy.deepcopy(original_packet)
        for field, role in (('metadata_source', 'metadata'), ('sample_source', 'source_sample_records')):
            self.assertEqual(packet[field]['path'], relocations[role]['original'])
            self.assertEqual(packet[field]['sha256'], relocations[role]['sha256'])
            packet[field]['path'] = relocations[role]['fixture']
        config = json.loads((target / 'pipeline/config.json').read_text())
        def source_wire(value):
            messages = self.pipeline.review_messages(value, mode='independent_audit', primary=None,
                review_profile=config['audit_profile'], source_binding_profile=config['source_binding_profile'])
            wire = [{'role': 'system', 'content': messages[0].content}, {'role': 'user', 'content': messages[1].content}]
            if config.get('source_contract_policy') == 'conservative_status_v1':
                wire[0]['content'] += self.pipeline.SOURCE_STATUS_CONTRACT
            return wire
        old_wire, new_wire = source_wire(original_packet), source_wire(packet)
        transport = packet_path.parent / 'source-transport.json'
        self.assertEqual(json.loads(transport.read_text()), json.loads(old_wire[1]['content']))
        packet_path.write_text(json.dumps(packet))
        transport.write_text(json.dumps(json.loads(new_wire[1]['content'])))
        for slot in slots:
            request_path = target / 'pipeline/api' / slot.name / 'request.json'
            request = json.loads(request_path.read_text())
            self.assertEqual(request['messages'][:2], old_wire)
            request['messages'][:2] = new_wire
            request_path.write_text(json.dumps(request))
        complete = stage / 'complete.json'
        completed['outputs'] = {name: update.sha(stage / name) for name in completed['outputs']}
        complete.write_text(json.dumps(completed))
        (stage / 'complete.sha256').write_text(update.sha(complete) + '\n')
        # This marker is fixture documentation, outside the replay input set.
        (target / 'SOFTWARE_FIXTURE.json').write_text(json.dumps({
            'scope': 'Path-relocated software regression fixture, not a new empirical assessment',
            'archive_path': str(origin), 'archive_file_sha256': original_hashes, 'locator_relocations': relocations,
            'new_model_calls': 0, 'response_bytes_changed': False}))
        self.assertEqual(original_hashes, {name: update.sha(origin / name) for name in names})
        update.source_design_replay_binding(target)
        return target

    def freeze(self):
        self.manifest.write_text(json.dumps(update.source_design_replay_binding(self.origin)))
        self.digest = update.sha(self.manifest)

    def rebind_source_fixture(self):
        """Rebind intentionally changed software fixtures, never original evidence."""
        stage = self.origin / 'pipeline/03_blinded_design'
        receipt = stage / 'complete.json'
        record = json.loads(receipt.read_text())
        record['outputs'] = {name: update.sha(stage / name) for name in record['outputs']}
        receipt.write_text(json.dumps(record))
        (stage / 'complete.sha256').write_text(update.sha(receipt) + '\n')
        for slot in (self.origin / 'pipeline/api').glob('design-*'):
            path = slot / 'receipt.json'; value = json.loads(path.read_text())
            value['response_sha256'] = update.sha(slot / 'response.json')
            path.write_text(json.dumps(value))
        self.freeze()

    def replay(self, **overrides):
        destination = self.root / ('replay-output-' + str(len(list(self.root.glob('replay-output-*')))))
        destination.mkdir()
        with patch.object(self.pipeline.BoundedOpenAI, 'invoke', side_effect=AssertionError('Replay invoked a model')):
            result = self.pipeline.replay_design(self.packet, self.manifest, self.digest, destination,
                                                 **{**self.parameters, **overrides})
        return result, destination

    def relocated_metadata(self):
        """Copy completed acquisition outputs; never alter the source archive."""
        self.freeze()
        original = self.origin / 'pipeline/01_metadata'
        current = self.root / ('new-pipeline-' + str(getattr(self, 'metadata_copy_count', 0))) / '01_metadata'
        self.metadata_copy_count = getattr(self, 'metadata_copy_count', 0) + 1
        shutil.copytree(original, current)
        receipt = current / 'complete.json'
        record = json.loads(receipt.read_text())
        for key, value in record['result'].items():
            if isinstance(value, str) and value.startswith(str(original)):
                record['result'][key] = str(current) + value[len(str(original)):]
        receipt.write_text(json.dumps(record))
        (current / 'complete.sha256').write_text(update.sha(receipt) + '\n')
        metadata = record['result']
        self.packet = self.pipeline.source_packet(Path(metadata['metadata']), self.packet['study_id'],
            self.packet['target'], context=metadata['study_context'],
            comparison_scope='source_defined_interventions', sample_source=metadata['source_sample_records'])
        self.parameters['current_metadata_stage'] = current
        return current

    def relocate_context_time(self, acquired=1893456000.25):
        """Regenerate only the acquisition timestamp in a sealed copied stage."""
        current = self.relocated_metadata()
        receipt = current / 'complete.json'
        record = json.loads(receipt.read_text())
        context = Path(record['result']['study_context'])
        value = json.loads(context.read_text())
        value['source_acquired_unix'] = acquired
        context.write_text(json.dumps(value))
        members = [name for name in record['outputs'] if Path(name).name == 'study-context.json']
        self.assertEqual(len(members), 1)
        record['outputs'][members[0]] = update.sha(context)
        receipt.write_text(json.dumps(record))
        (current / 'complete.sha256').write_text(update.sha(receipt) + '\n')
        metadata = record['result']
        self.packet = self.pipeline.source_packet(Path(metadata['metadata']), self.packet['study_id'],
            self.packet['target'], context=metadata['study_context'],
            comparison_scope='source_defined_interventions', sample_source=metadata['source_sample_records'])
        return current

    def test_fixture_replays_after_archive_copy_with_original_location_unavailable(self):
        relocated_archive = self.root / 'portable-checkout/copied-source'
        shutil.copytree(self.origin, relocated_archive)
        original_location = self.origin
        original_location.rename(self.root / 'unavailable-original-location')
        self.assertFalse(original_location.exists())
        before = {str(path.relative_to(relocated_archive)): update.sha(path)
                  for path in relocated_archive.rglob('*') if path.is_file()}
        self.origin = self.copy_source_archive(relocated_archive)
        self.packet = json.loads((self.origin / 'pipeline/03_blinded_design/attempt-000/source-packet.json').read_text())
        self.freeze()
        result, output = self.replay()
        self.assertFalse(original_location.exists())
        self.assertEqual(before, {str(path.relative_to(relocated_archive)): update.sha(path)
                                for path in relocated_archive.rglob('*') if path.is_file()})
        self.assertEqual(len(result['comparisons']), 4)
        self.assertEqual(json.loads((output / 'source-design-replay.json').read_text())['new_model_calls'], 0)
        self.assertTrue((self.origin / 'SOFTWARE_FIXTURE.json').is_file())

    def test_fixture_import_rejects_changed_source_bytes_before_rebasing(self):
        changed_archive = self.root / 'copied-invalid-source'
        shutil.copytree(self.origin, changed_archive)
        source = changed_archive / 'pipeline/01_metadata/attempt-000/metadata.csv'
        source.write_text(source.read_text() + '\n')
        with self.assertRaises(AssertionError): self.copy_source_archive(changed_archive)

    def test_archived_dnam_design_replays_without_outcomes_calls_or_changed_date(self):
        self.assertFalse((self.origin / 'pipeline/04_preprocessing').exists())
        result, output = self.replay()
        previous = json.loads((self.origin / 'pipeline/03_blinded_design/complete.json').read_text())['result']
        self.assertEqual(result, previous)
        self.assertEqual(len(result['comparisons']), 4)
        receipt = json.loads((output / 'source-design-replay.json').read_text())
        self.assertEqual(receipt['new_model_calls'], 0)
        self.assertEqual(receipt['new_estimated_cost_usd'], 0)
        self.assertFalse(receipt['new_independent_assessment'])
        original = self.origin / 'pipeline/api/design-000/response.json'
        self.assertEqual(receipt['original_source_calls'][0]['response_created_unix'], json.loads(original.read_text())['created'])
        archived = receipt['archived_files']['pipeline/api/design-000/response.json']
        self.assertEqual(update.sha(output / archived['path']), update.sha(original))
        self.assertEqual(len(list(output.rglob('source-packet.json'))), 1)
        self.assertEqual(len(list(output.rglob('design.json'))), 1)

    def default_tier_fixture(self):
        """Opt-in SOFTWARE MUTATION of a private copy; no new model evidence."""
        config_path=self.origin/'pipeline/config.json'
        config=json.loads(config_path.read_text());config['service_tier']='default'
        config_path.write_text(json.dumps(config))
        version_path=self.origin/'version.json';version=json.loads(version_path.read_text())
        for record in (version,version['identity']):
            record['pipeline_options']['service_tier']='default'
        version['version_id']=update.identity(version['identity'])
        version_path.write_text(json.dumps(version))
        (self.origin/'version.json.sha256').write_text(update.sha(version_path)+'\n')
        for slot in (self.origin/'pipeline/api').glob('design-*'):
            for filename in ('request.json','response.json'):
                path=slot/filename;value=json.loads(path.read_text());value['service_tier']='default'
                path.write_text(json.dumps(value))
        marker=self.origin/'SOFTWARE_FIXTURE.json';value=json.loads(marker.read_text())
        value.update(response_bytes_changed=True,mutation_scope='Explicit service-tier software fixture only')
        marker.write_text(json.dumps(value))
        self.rebind_source_fixture()

    def test_explicit_tier_cannot_adopt_unbound_legacy_source_assessment(self):
        from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
        with self.assertRaisesRegex(PipelineIntegrityError,'Archived source service tier'):
            self.replay(service_tier='default')
        version=update.read_sealed(self.origin/'version.json')
        sources={**version['sources'],'source_design_replay':{'path':str(self.manifest),'sha256':self.digest}}
        options={**version['pipeline_options'],'service_tier':'default'}
        entry={key:version[key] for key in ('study_id','modality','target')}
        with self.assertRaisesRegex(update.UpdateIntegrityError,'source-review profile differs'):
            update.validate_replay_sources(entry,sources,options)

    def test_source_replay_binds_explicit_tier_across_config_version_request_and_response(self):
        from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
        self.default_tier_fixture()
        result,output=self.replay(service_tier='default')
        self.assertEqual(result,json.loads((self.origin/'pipeline/03_blinded_design/complete.json').read_text())['result'])
        self.assertEqual(json.loads((output/'source-design-replay.json').read_text())['new_model_calls'],0)
        with self.assertRaisesRegex(PipelineIntegrityError,'Archived source service tier'):
            self.replay()
        config=self.origin/'pipeline/config.json';version=self.origin/'version.json'
        request=self.origin/'pipeline/api/design-000/request.json'
        response=self.origin/'pipeline/api/design-000/response.json'
        originals={path:path.read_bytes() for path in (config,version,request,response)}
        mutations=((config,(),None),(version,('pipeline_options',),None),
                   (version,('identity','pipeline_options'),None),(request,(),None),
                   (response,(),None),(response,(),'priority'))
        for path,parents,tier in mutations:
            for original,data in originals.items():original.write_bytes(data)
            value=json.loads(path.read_text());target=value
            for key in parents:target=target[key]
            if tier is None:del target['service_tier']
            else:target['service_tier']=tier
            path.write_text(json.dumps(value))
            version_value=json.loads(version.read_text())
            version_value['version_id']=update.identity(version_value['identity'])
            version.write_text(json.dumps(version_value))
            (version.parent/'version.json.sha256').write_text(update.sha(version)+'\n')
            self.rebind_source_fixture()
            with self.subTest(path=path.name,parents=parents,tier=tier), self.assertRaises(PipelineIntegrityError):
                self.replay(service_tier='default')

    def test_controller_checks_tier_in_both_origin_version_profiles(self):
        self.default_tier_fixture()
        version_path=self.origin/'version.json';original=version_path.read_bytes()
        version=json.loads(original)
        entry={key:version[key] for key in ('study_id','modality','target')}
        options=copy.deepcopy(version['pipeline_options'])
        def validate():
            sources={**version['sources'],'source_design_replay':{'path':str(self.manifest),'sha256':self.digest}}
            update.validate_replay_sources(entry,sources,options)
        validate()
        for nested in (False,True):
            version=json.loads(original);profile=version['identity'] if nested else version
            del profile['pipeline_options']['service_tier']
            version['version_id']=update.identity(version['identity'])
            version_path.write_text(json.dumps(version))
            (version_path.parent/'version.json.sha256').write_text(update.sha(version_path)+'\n')
            self.freeze()
            with self.subTest(identity=nested), self.assertRaisesRegex(update.UpdateIntegrityError,'profile differs'):
                validate()

    def test_replay_preserves_a_rejected_source_reply_and_its_supported_retry(self):
        origin = self.base / 'registry-006/studies/human_rna/GSE127814/0624025b3deb6ac935f86a064728c463ad1bd828cd05f44c34a7df4a3cbca079'
        self.origin = self.copy_source_archive(origin)
        self.packet = json.loads((self.origin / 'pipeline/03_blinded_design/attempt-000/source-packet.json').read_text())
        self.freeze()
        result, output = self.replay()
        receipt = json.loads((output / 'source-design-replay.json').read_text())
        self.assertEqual(len(receipt['original_source_calls']), 2)
        self.assertEqual([item['slot'] for item in receipt['original_source_calls']], ['design-000', 'design-001'])
        self.assertTrue(all(c['proposed_action'] == 'withhold_inference' for c in result['comparisons']))

    def test_replay_checks_v2_feedback_and_rejects_forged_or_unknown_feedback(self):
        """Versioned software fixture; original responses and evidence remain unchanged."""
        from audit_system.agent_workflow_2026_09_10 import design_validation_feedback
        origin = self.base / 'registry-006/studies/human_rna/GSE127814/0624025b3deb6ac935f86a064728c463ad1bd828cd05f44c34a7df4a3cbca079'
        self.origin = self.copy_source_archive(origin)
        stage = self.origin / 'pipeline/03_blinded_design/attempt-000'
        self.packet = json.loads((stage/'source-packet.json').read_text())
        path = stage/'response-validation/000-validation.json'
        legacy = json.loads(path.read_text())
        raw = json.loads((self.origin/'pipeline/api/design-000/response.json').read_text())['choices'][0]['message']['content']
        feedback = design_validation_feedback(raw, self.packet, ValueError(legacy['error']))
        for key in ('status_issues','status_contract'):
            if key in legacy: feedback[key] = legacy[key]
        next_request = self.origin/'pipeline/api/design-001/request.json'
        request = json.loads(next_request.read_text())
        for mutation in ('valid_v2','invented_coverage','unknown_schema'):
            with self.subTest(mutation=mutation):
                candidate = copy.deepcopy(feedback)
                if mutation == 'invented_coverage': candidate['coverage_issue_count'] += 1
                if mutation == 'unknown_schema': candidate['schema'] = 'unavailable-feedback-v99'
                path.write_text(json.dumps(candidate))
                request['messages'][-1]['content'] = json.dumps(candidate, ensure_ascii=False)
                next_request.write_text(json.dumps(request))
                self.rebind_source_fixture()
                if mutation == 'valid_v2':
                    result, output = self.replay()
                    self.assertTrue(all(c['proposed_action'] == 'withhold_inference' for c in result['comparisons']))
                    self.assertEqual(json.loads((output/'source-design-replay.json').read_text())['new_model_calls'], 0)
                else:
                    with self.assertRaisesRegex(ValueError, 'feedback'):
                        self.replay()

    def test_relocated_metadata_requires_equal_completed_source_bytes_and_preserves_request(self):
        current = self.relocated_metadata()
        original_request = self.origin / 'pipeline/api/design-000/request.json'
        request_hash = update.sha(original_request)
        result, output = self.replay()
        self.assertEqual(result, json.loads((self.origin / 'pipeline/03_blinded_design/complete.json').read_text())['result'])
        receipt = json.loads((output / 'source-design-replay.json').read_text())
        self.assertEqual(len(receipt['source_locator_relocations']), 2)
        self.assertEqual({r['role'] for r in receipt['source_locator_relocations']}, {'metadata', 'source_sample_records'})
        for relocation in receipt['source_locator_relocations']:
            self.assertEqual(relocation['original']['sha256'], relocation['current']['sha256'])
            self.assertEqual(relocation['current']['metadata_stage_sha256'], update.sha(current / 'complete.json'))
        self.assertEqual(json.loads((output / 'source-packet.json').read_text()), self.packet)
        self.assertEqual(update.sha(original_request), request_hash)
        archived_request = output / receipt['archived_files']['pipeline/api/design-000/request.json']['path']
        self.assertEqual(update.sha(archived_request), request_hash)
        self.assertEqual(receipt['new_model_calls'], 0)

    def test_timestamp_only_context_relocation_is_stage_bound_and_preserves_request(self):
        current = self.relocate_context_time()
        original_request = self.origin / 'pipeline/api/design-000/request.json'
        request_hash = update.sha(original_request)
        result, output = self.replay()
        previous = json.loads((self.origin / 'pipeline/03_blinded_design/complete.json').read_text())['result']
        self.assertEqual(result, previous)
        receipt = json.loads((output / 'source-design-replay.json').read_text())
        relocations = {row['role']: row for row in receipt['source_locator_relocations']}
        self.assertEqual(set(relocations), {'metadata', 'source_sample_records', 'study_context'})
        context = relocations['study_context']
        self.assertEqual(context['field'], '$.context_sha256')
        self.assertNotEqual(context['original']['sha256'], context['current']['sha256'])
        self.assertNotEqual(context['original']['source_acquired_unix'], context['current']['source_acquired_unix'])
        self.assertEqual(context['current']['metadata_stage_sha256'], update.sha(current / 'complete.json'))
        self.assertEqual(update.sha(original_request), request_hash)
        self.assertEqual(receipt['new_model_calls'], 0)

    def test_timestamp_context_relocation_rejects_other_content_or_document_change(self):
        current = self.relocate_context_time()
        receipt = current / 'complete.json'
        original_receipt = receipt.read_text()
        record = json.loads(original_receipt)
        path = Path(record['result']['study_context'])
        original_context = path.read_text()
        original_packet = copy.deepcopy(self.packet)
        for mutation in ('interpretation', 'document', 'source_date'):
            with self.subTest(mutation=mutation):
                value = json.loads(original_context)
                if mutation == 'interpretation':
                    value['interpretation'] += ' changed software fixture'
                elif mutation == 'document':
                    value['documents'][0]['fields']['title'] += ' changed software fixture'
                else:
                    value['documents'][0]['source_date'] = '2099-01-01'
                path.write_text(json.dumps(value))
                record = json.loads(original_receipt)
                member = next(name for name in record['outputs'] if Path(name).name == 'study-context.json')
                record['outputs'][member] = update.sha(path); receipt.write_text(json.dumps(record))
                (current / 'complete.sha256').write_text(update.sha(receipt) + '\n')
                self.packet['context_sha256'] = update.sha(path)
                with self.assertRaisesRegex(ValueError, 'changed beyond its acquisition time'):
                    self.replay()
                path.write_text(original_context); receipt.write_text(original_receipt)
                (current / 'complete.sha256').write_text(update.sha(receipt) + '\n')
                self.packet = copy.deepcopy(original_packet)

    def test_timestamp_context_relocation_requires_finite_positive_nonboolean_times(self):
        current = self.relocate_context_time()
        receipt = current / 'complete.json'; original_receipt = receipt.read_text()
        record = json.loads(original_receipt); path = Path(record['result']['study_context'])
        original_context = path.read_text(); original_packet = copy.deepcopy(self.packet)
        for invalid in ('MISSING', None, True, 0, -1, float('inf')):
            with self.subTest(invalid=invalid):
                value = json.loads(original_context)
                if invalid == 'MISSING':
                    value.pop('source_acquired_unix')
                else:
                    value['source_acquired_unix'] = invalid
                path.write_text(json.dumps(value))
                record = json.loads(original_receipt)
                member = next(name for name in record['outputs'] if Path(name).name == 'study-context.json')
                record['outputs'][member] = update.sha(path); receipt.write_text(json.dumps(record))
                (current / 'complete.sha256').write_text(update.sha(receipt) + '\n')
                self.packet['context_sha256'] = update.sha(path)
                with self.assertRaisesRegex(ValueError, 'finite and positive'):
                    self.replay()
                path.write_text(original_context); receipt.write_text(original_receipt)
                (current / 'complete.sha256').write_text(update.sha(receipt) + '\n')
                self.packet = copy.deepcopy(original_packet)

    def test_context_hash_and_stage_role_must_bind_exact_payload(self):
        current = self.relocate_context_time()
        digest = self.packet['context_sha256']
        self.packet['context_sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'not bound to its completed metadata role'):
            self.replay()
        self.packet['context_sha256'] = digest
        receipt = current / 'complete.json'; record = json.loads(receipt.read_text())
        unbound = current / 'attempt-000/unbound-context.json'
        source = Path(record['result']['study_context']); shutil.copyfile(source, unbound)
        record['result']['study_context'] = str(unbound)
        receipt.write_text(json.dumps(record)); (current / 'complete.sha256').write_text(update.sha(receipt) + '\n')
        with self.assertRaisesRegex(ValueError, 'not bound to its completed metadata role'):
            self.replay()

    def test_empirical_rna_archives_replay_with_timestamp_only_context_change_and_no_calls(self):
        archives = {
            'GSE203635': '94fe164f634bc133a1391ff7b2a41b9cebec008268bab7fa829160963c3c27b4',
            'GSE227798': 'b04349fb2721b56cdb7e005613c824762f59cac0cabae2f2ec5a99a6a0abfb0c',
            'GSE240194': '30ef24bc614946dc0555ab34ecd0c971d74887f1f7d5aebb9a51a7d4e8200950',
        }
        for index, (study, version) in enumerate(archives.items()):
            with self.subTest(study=study):
                origin = self.base / 'registry-012/studies/human_rna' / study / version
                self.origin = self.copy_source_archive(origin)
                self.packet = json.loads((self.origin / 'pipeline/03_blinded_design/attempt-000/source-packet.json').read_text())
                self.freeze()
                current = self.relocate_context_time(1893456000.25 + index)
                result, output = self.replay()
                previous = json.loads((self.origin / 'pipeline/03_blinded_design/complete.json').read_text())['result']
                self.assertEqual(result, previous)
                receipt = json.loads((output / 'source-design-replay.json').read_text())
                self.assertEqual(receipt['new_model_calls'], 0)
                self.assertEqual(receipt['new_estimated_cost_usd'], 0)
                self.assertEqual({row['role'] for row in receipt['source_locator_relocations']},
                                 {'metadata', 'source_sample_records', 'study_context'})
                self.assertEqual(receipt['source_locator_relocations'][-1]['current']['metadata_stage_sha256'],
                                 update.sha(current / 'complete.json'))

    def test_locator_relocation_without_bound_current_stage_is_rejected(self):
        self.relocated_metadata()
        with self.assertRaisesRegex(ValueError, 'current completed metadata stage'):
            self.replay(current_metadata_stage=None)

    def test_unbound_or_outside_stage_equal_byte_locator_is_rejected(self):
        current = self.relocated_metadata()
        for name in ('unbound', 'outside'):
            with self.subTest(name=name):
                copied = (current / 'attempt-000/unbound.csv' if name == 'unbound' else self.root / 'outside.csv')
                shutil.copyfile(current / 'attempt-000/metadata.csv', copied)
                self.packet['metadata_source']['path'] = str(copied)
                with self.assertRaisesRegex(ValueError, 'not bound to its completed metadata role'):
                    self.replay()

    def test_relocation_cannot_switch_roles_even_when_file_is_a_stage_member(self):
        current = self.relocated_metadata()
        self.packet['metadata_source']['path'] = str(current / 'attempt-000/samples.soft')
        with self.assertRaisesRegex(ValueError, 'not bound to its completed metadata role'):
            self.replay()

    def test_one_relocated_path_still_requires_both_roles_in_current_stage(self):
        self.relocated_metadata()
        original = json.loads((self.origin / 'pipeline/03_blinded_design/attempt-000/source-packet.json').read_text())
        self.packet['sample_source']['path'] = original['sample_source']['path']
        with self.assertRaisesRegex(ValueError, 'completed metadata role: sample_source'):
            self.replay()

    def test_relocated_source_hash_or_content_change_is_rejected(self):
        current = self.relocated_metadata()
        digest = self.packet['metadata_source']['sha256']
        self.packet['metadata_source']['sha256'] = '0' * 64
        with self.assertRaisesRegex(ValueError, 'identical blinded source packet'):
            self.replay()
        self.packet['metadata_source']['sha256'] = digest
        path = current / 'attempt-000/metadata.csv'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'inventory changed'):
            self.replay()

    def test_relocation_preserves_every_nonpath_packet_field(self):
        self.relocated_metadata()
        self.packet['samples'][0]['title'] += ' changed software fixture'
        with self.assertRaisesRegex(ValueError, 'identical blinded source packet'):
            self.replay()

    def test_relocation_requires_valid_receipt_seal_and_no_stage_escape(self):
        current = self.relocated_metadata()
        receipt = current / 'complete.json'
        record = json.loads(receipt.read_text())
        receipt.write_text(json.dumps(record) + '\n')
        with self.assertRaisesRegex(ValueError, 'sealed completed metadata stage'):
            self.replay()
        outside = current.parent / 'outside.csv'
        shutil.copyfile(current / 'attempt-000/metadata.csv', outside)
        record['outputs']['../outside.csv'] = update.sha(outside)
        receipt.write_text(json.dumps(record))
        (current / 'complete.sha256').write_text(update.sha(receipt) + '\n')
        with self.assertRaisesRegex(ValueError, 'escapes its stage'):
            self.replay()

    def test_changed_packet_profile_or_source_model_rejected(self):
        for field, value in [('source_model', 'gpt-5.6-sol'), ('source_contract_policy', 'strict'),
                             ('source_binding_profile', 'legacy')]:
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.replay(**{field: value})
        self.packet['samples'][0]['sample_id'] = 'GSM999999999'
        with self.assertRaisesRegex(ValueError, 'identical blinded source packet'):
            self.replay()

    def test_nested_tampering_and_omitted_manifest_member_rejected(self):
        value = json.loads(self.manifest.read_text())
        value['files'].pop('pipeline/api/design-000/request.json')
        self.manifest.write_text(json.dumps(value)); self.digest = update.sha(self.manifest)
        with self.assertRaisesRegex(ValueError, 'inventory changed'):
            self.replay()
        self.freeze()
        path = self.origin / 'pipeline/api/design-000/request.json'
        path.write_text(path.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'inventory changed'):
            self.replay()

    def test_outcome_contaminated_source_request_rejected_after_fixture_rebinding(self):
        path = self.origin / 'pipeline/api/design-000/request.json'
        request = json.loads(path.read_text())
        request['messages'].append({'role': 'user', 'content': 'SOFTWARE FIXTURE: preceding outcome exposed'})
        path.write_text(json.dumps(request)); self.freeze()
        with self.assertRaisesRegex(ValueError, 'Original source request'):
            self.replay()

    def test_current_validator_cannot_change_archived_groups_or_limitations(self):
        original = self.pipeline.parse_design
        def changed(*args, **kwargs):
            value = original(*args, **kwargs)
            value['comparisons'][0]['unresolved'] = []
            return value
        with patch.object(self.pipeline, 'parse_design', side_effect=changed):
            with self.assertRaisesRegex(ValueError, 'would change the archived design'):
                self.replay()

    def test_unsupported_quotation_still_fails_current_source_parser(self):
        api = self.origin / 'pipeline/api/design-000/response.json'
        response = json.loads(api.read_text())
        value = json.loads(response['choices'][0]['message']['content'])
        value['evidence'][0]['quote'] = 'UNSUPPORTED SOFTWARE FIXTURE QUOTATION'
        raw = json.dumps(value)
        response['choices'][0]['message']['content'] = raw
        api.write_text(json.dumps(response))
        stage = self.origin / 'pipeline/03_blinded_design/attempt-000'
        (stage / 'contract-adjustments/000-original-response.txt').write_text(raw)
        (stage / 'contract-adjustments/000-candidate.json').write_text(raw)
        self.rebind_source_fixture()
        with self.assertRaisesRegex(ValueError, 'Unsupported source quotation'):
            self.replay()

    def test_nested_sources_enter_controller_hashes_and_resume_checks(self):
        helper = UpdatePipelineTests(); helper.setUp(); self.addCleanup(helper.doCleanups)
        entry = helper.study()
        entry['source_design_replay'] = {'path': str(self.manifest), 'sha256': self.digest}
        sources = update.collect_sources(entry, self.root / 'acquisition', helper.plan([entry])['limits'], {})
        roles = [key for key in sources if key.startswith('source_design_replay:')]
        self.assertEqual(len(roles), len(json.loads(self.manifest.read_text())['files']))
        identity = update.identity({key: value['sha256'] for key, value in sources.items()})
        self.assertNotEqual(identity, update.identity({key: value['sha256'] for key, value in sources.items() if key not in roles}))
        nested = Path(sources[roles[-1]]['path']); nested.write_text(nested.read_text() + '\n')
        with self.assertRaisesRegex(ValueError, 'Frozen source changed'):
            update.verify_sources(sources)

    def test_stage_escape_is_rejected_without_reading_outcomes(self):
        path = self.origin / 'pipeline/03_blinded_design/complete.json'
        record = json.loads(path.read_text())
        record['outputs']['../05_clock_predictions/forbidden.csv'] = 'a' * 64
        path.write_text(json.dumps(record)); (path.parent / 'complete.sha256').write_text(update.sha(path) + '\n')
        with self.assertRaisesRegex(ValueError, 'escapes its stage'):
            self.freeze()

    def test_default_has_no_replay_and_unbound_cli_replay_fails(self):
        argv = ['--study', 'GSE1', '--modality', 'human_rna', '--output', str(self.root / 'unused')]
        self.assertIsNone(self.pipeline.parser().parse_args(argv).source_design_replay)
        with self.assertRaisesRegex(ValueError, 'requires both'):
            self.pipeline.main(argv + ['--source-design-replay', str(self.manifest)])


class SourceDiagnosticVersionTests(unittest.TestCase):
    """Resealed software copies test dispatch; archived assessments stay untouched."""

    setUp = SourceDesignReplayTests.setUp
    copy_source_archive = SourceDesignReplayTests.copy_source_archive
    freeze = SourceDesignReplayTests.freeze
    rebind_source_fixture = SourceDesignReplayTests.rebind_source_fixture
    replay = SourceDesignReplayTests.replay

    def current_fixture(self):
        from audit_system import source_packet as source
        from audit_system import agent_workflow_2026_09_10 as workflow
        self.attempt = self.origin / 'pipeline/03_blinded_design/attempt-000'
        self.original_diagnostic_implementation = {
            name: json.loads((self.origin / 'pipeline/config.json').read_text())['implementation'][name]
            for name in ('audit_system/source_packet.py', 'audit_system/agent_workflow_2026_09_10.py')}
        design = json.loads((self.attempt / 'design.json').read_text())
        self.original_diagnostics = {}
        def replace(name, value):
            path = self.attempt / name
            self.original_diagnostics[name] = path.read_text()
            path.write_text(json.dumps(value))
        replace('domain-verification.json', workflow.verify_domain_review(design, self.packet))
        for comparison in design['comparisons']:
            packet = {**self.packet, 'target': comparison['target']}
            for name, reader in (('comparison-binding.json', source.verify_requested_comparison),
                                 ('comparison-proposals.json', source.propose_source_comparisons),
                                 ('grouping-diagnostics.json', source.source_grouping_diagnostics)):
                replace(comparison['comparison_id'] + '/' + name,
                        reader(packet, comparison['case_ids'], comparison['control_ids']))
        self.set_implementation({
            'audit_system/source_packet.py': update.sha(Path(source.__file__)),
            'audit_system/agent_workflow_2026_09_10.py': update.sha(Path(workflow.__file__)),
        })
        self.rebind_source_fixture()
        return design

    def set_implementation(self, declarations, *, config=True, version=True, identity=True):
        for relative, active in (('pipeline/config.json', config), ('version.json', version)):
            if not active:
                continue
            path = self.origin / relative
            value = json.loads(path.read_text())
            value['implementation'].update(declarations)
            if relative == 'version.json' and identity:
                value['identity']['implementation'].update(declarations)
                value['version_id'] = update.identity(value['identity'])
            path.write_text(json.dumps(value))
            if relative == 'version.json':
                (self.origin / 'version.json.sha256').write_text(update.sha(path) + '\n')

    def test_legacy_replay_records_recognized_profile_and_preserves_result(self):
        from audit_system.source_packet import SOURCE_DIAGNOSTIC_LEGACY
        result, output = self.replay()
        receipt = json.loads((output / 'source-design-replay.json').read_text())
        self.assertEqual(receipt['source_diagnostic_profile'], SOURCE_DIAGNOSTIC_LEGACY)
        self.assertEqual([item['slot'] for item in receipt['original_source_calls']], ['design-000'])
        self.assertEqual(result, json.loads((self.origin / 'pipeline/03_blinded_design/complete.json').read_text())['result'])
        self.assertEqual(receipt['new_model_calls'], 0)

    def test_current_marker_and_implementation_replay_exactly(self):
        from audit_system.source_packet import SOURCE_DIAGNOSTIC_CURRENT
        expected = self.current_fixture()
        result, output = self.replay()
        self.assertEqual(result, expected)
        receipt = json.loads((output / 'source-design-replay.json').read_text())
        self.assertEqual(receipt['source_diagnostic_profile'], SOURCE_DIAGNOSTIC_CURRENT)
        self.assertEqual([item['slot'] for item in receipt['original_source_calls']], ['design-000'])
        self.assertEqual(receipt['new_model_calls'], 0)

    def test_current_marker_removal_or_relabeling_cannot_select_legacy(self):
        self.current_fixture()
        path = self.attempt / 'domain-verification.json'
        original = path.read_text()
        for marker in (None, 'source-consistency-v1', 'source-consistency-v99', True):
            with self.subTest(marker=marker):
                value = json.loads(original)
                if marker is None:
                    value.pop('source_diagnostic_profile')
                else:
                    value['source_diagnostic_profile'] = marker
                path.write_text(json.dumps(value)); self.rebind_source_fixture()
                with self.assertRaisesRegex(ValueError, 'diagnostic'):
                    self.replay()

    def test_coordinated_legacy_diagnostics_cannot_downgrade_current_implementation(self):
        self.current_fixture()
        for name, original in self.original_diagnostics.items():
            (self.attempt / name).write_text(original)
        self.rebind_source_fixture()
        with self.assertRaisesRegex(ValueError, 'recognized legacy implementation'):
            self.replay()

    def test_unknown_or_conflicting_implementation_is_rejected(self):
        self.current_fixture()
        self.set_implementation({'audit_system/source_packet.py': 'f' * 64}, version=False)
        self.rebind_source_fixture()
        with self.assertRaisesRegex(ValueError, 'implementation differs'):
            self.replay()
        self.set_implementation({'audit_system/source_packet.py': 'f' * 64})
        self.rebind_source_fixture()
        with self.assertRaisesRegex(ValueError, 'profile/implementation'):
            self.replay()

    def test_downgraded_top_level_implementation_cannot_override_bound_current_identity(self):
        self.current_fixture()
        for name, original in self.original_diagnostics.items():
            (self.attempt / name).write_text(original)
        self.set_implementation(self.original_diagnostic_implementation, identity=False)
        self.rebind_source_fixture()
        with self.assertRaisesRegex(ValueError, 'implementation differs between config, version and identity'):
            self.replay()

    def test_mixed_nested_profiles_and_changed_domain_result_reject(self):
        self.current_fixture()
        path = self.attempt / 'domain-verification.json'; original = path.read_text()
        for mutation in ('target', 'design', 'domain'):
            with self.subTest(mutation=mutation):
                value = json.loads(original)
                if mutation == 'target':
                    value['comparisons'][0]['target_binding']['profile'] = 'exact-source-comparison-v1'
                elif mutation == 'design':
                    value['comparisons'][0]['source_design_profile']['profile'] = 'unavailable-profile-v99'
                else:
                    value['domains']['grouping']['status'] = 'PASS'
                path.write_text(json.dumps(value)); self.rebind_source_fixture()
                with self.assertRaisesRegex(ValueError, 'profiles|archived design'):
                    self.replay()

    def test_resealed_changed_comparison_descriptors_reject(self):
        self.current_fixture()
        comparison = json.loads((self.attempt / 'design.json').read_text())['comparisons'][0]
        base = self.attempt / comparison['comparison_id']
        for filename in ('comparison-binding.json', 'comparison-proposals.json', 'grouping-diagnostics.json'):
            with self.subTest(filename=filename):
                path = base / filename; original = path.read_text(); value = json.loads(original)
                if isinstance(value, list):
                    value.append({'software_fixture_unbound_proposal': True})
                else:
                    value['software_fixture_unbound_claim'] = True
                path.write_text(json.dumps(value)); self.rebind_source_fixture()
                with self.assertRaisesRegex(ValueError, 'comparison diagnostics differ'):
                    self.replay()
                path.write_text(original); self.rebind_source_fixture()

    def test_current_default_retains_unknown_and_contradiction_checks(self):
        from audit_system import source_packet as source
        from audit_system import agent_workflow_2026_09_10 as workflow
        packet = {'study_id': 'SOFTWARE_FIXTURE', 'source_sha256': 'f' * 64,
                  'target': 'treatment: drug versus treatment: vehicle',
                  'samples': [{'sample_id': sid, 'characteristics_ch1': 'treatment: ' + arm}
                              for sid, arm in [('A', 'unknown'), ('B', 'drug'), ('C', 'vehicle'), ('D', 'vehicle')]]}
        self.assertEqual(source.verify_requested_comparison(packet, ['A', 'B'], ['C', 'D'])['status'], 'CANNOT VERIFY')
        self.assertEqual(source.verify_requested_comparison(packet, ['A', 'B'], ['C', 'D'],
            diagnostic_profile=source.SOURCE_DIAGNOSTIC_LEGACY)['status'], 'FAIL')
        packet['samples'][0]['characteristics_ch1'] = 'treatment: vehicle'
        self.assertEqual(source.verify_requested_comparison(packet, ['A', 'B'], ['C', 'D'])['status'], 'FAIL')
        for reader in (source.verify_requested_comparison, source.verify_source_design,
                       source.propose_source_comparisons, source.source_grouping_diagnostics):
            with self.subTest(reader=reader.__name__), self.assertRaisesRegex(ValueError, 'Unknown source diagnostic profile'):
                reader(packet, ['A', 'B'], ['C', 'D'], diagnostic_profile='unknown-v99')


class DnamIdatInterfaceTests(unittest.TestCase):
    """Receipt/interface software fixtures; no IDAT decoding or empirical clocks."""

    def setUp(self):
        self.helper = UpdatePipelineTests()
        self.helper.setUp()
        self.addCleanup(self.helper.doCleanups)
        self.root = self.helper.root
        self.entry = self.helper.study()
        self.entry['modality'] = 'human_dnam'
        self.wrapper = self.root / 'idat-provenance.json'
        self.wrapper.write_text('{"software_fixture": true}\n')
        raw = self.root / 'channel.idat'
        raw.write_bytes(b'SOFTWARE FIXTURE: no IDAT data')
        self.files = {role: dict(value) for role, value in self.entry['sources'].items()}
        self.files['raw:GSM1:Grn'] = self.helper.declaration(raw)
        self.verified = {'study_id': 'GSE1', 'sample_ids': ['GSM1'], 'files': self.files,
                         'measurement': {'profile': 'SOFTWARE_FIXTURE', 'scientific_status': 'CANNOT VERIFY'},
                         'scientific_status': 'CANNOT VERIFY'}
        self.entry['sources']['dnam_idat_provenance'] = self.helper.declaration(self.wrapper)
        verifier = patch('audit_system.report_2026_09_10.verify_dnam_idat_provenance',
                         return_value=self.verified, create=True)
        self.verifier = verifier.start()
        self.addCleanup(verifier.stop)

    def collect(self):
        return update.collect_sources(self.entry, self.root / 'acquisition',
                                      self.helper.plan([self.entry])['limits'], {})

    def test_transitive_idat_files_enter_version_and_changed_raw_stops_resume(self):
        result = self.helper.execute(self.helper.plan([self.entry]))
        version = json.loads((Path(result['studies'][0]['version_root']) / 'version.json').read_text())
        for role, value in self.files.items():
            key = 'dnam_idat_provenance:' + role
            self.assertEqual(version['identity']['source_sha256'][key], value['sha256'])
            self.assertEqual(version['sources'][key]['path'], value['path'])
        Path(self.files['raw:GSM1:Grn']['path']).write_bytes(b'changed software fixture')
        resumed = self.helper.execute(self.helper.plan([self.entry]))
        self.assertEqual(resumed['studies'][0]['status'], 'UNRESOLVED')
        self.assertIn('nested file differs', resumed['studies'][0]['error'])
        self.assertEqual(len(self.helper.calls), 1)

    def test_idat_requires_distinct_dnam_source_role_and_primary_bindings(self):
        for change in ('rna', 'supplement', 'foreign_series'):
            with self.subTest(change=change):
                saved = copy.deepcopy(self.entry)
                if change == 'rna':
                    self.entry['modality'] = 'human_rna'
                elif change == 'supplement':
                    self.entry['sources']['dnam_supplement'] = self.entry['sources']['matrix']
                else:
                    original = Path(self.entry['sources']['series']['path'])
                    duplicate = self.root / 'same-bytes-different-series.soft'
                    duplicate.write_bytes(original.read_bytes())
                    self.entry['sources']['series'] = self.helper.declaration(duplicate)
                with self.assertRaisesRegex(ValueError, 'DNAm IDAT provenance'):
                    self.collect()
                self.entry = saved

    def test_cached_refresh_rechecks_external_idat_graph_without_rerunning_update(self):
        result = self.helper.execute(self.helper.plan([self.entry]))
        external = update._refresh_update_external(result)
        raw = Path(self.files['raw:GSM1:Grn']['path'])
        self.assertEqual(external[str(raw)], self.files['raw:GSM1:Grn']['sha256'])
        refresh = self.root / 'refresh'
        update._refresh_phase(refresh, 'update-input', {}, lambda unused: ({'software_fixture': True}, external))
        raw.write_bytes(b'changed after refresh software fixture')
        with self.assertRaisesRegex(ValueError, 'Refresh external output changed'):
            update._refresh_phase(refresh, 'update-input', {}, lambda unused: self.fail('cached update was rerun'))
        self.assertEqual(len(self.helper.calls), 1)

    def test_idat_graph_rejects_wrong_scope_and_ambiguous_samples(self):
        for field, value in [('study_id', 'GSE2'), ('scientific_status', 'PASS'),
                             ('sample_ids', ['GSM1', 'GSM1']), ('files', {'matrix': self.files['matrix']})]:
            with self.subTest(field=field):
                changed = {**self.verified, field: value}
                self.verifier.return_value = changed
                with self.assertRaisesRegex(ValueError, 'DNAm IDAT provenance'):
                    self.collect()

    def test_controller_forwards_idat_role_to_child_cli(self):
        from src import pipeline_2026_09_10 as pipeline
        sources = self.collect()
        version_path = self.root / 'child-version.json'
        version = {'study_id': 'GSE1', 'modality': 'human_dnam', 'target': 'Software fixture',
                   'sources': sources, 'implementation': {}, 'pipeline_options': {}}
        update.seal(version_path, version)
        request = {'version_path': str(version_path), 'phase': 'pipeline', 'until': 'input',
                   'pipeline_root': str(self.root / 'pipeline'), 'cache': str(self.root / 'cache'),
                   'result_path': str(self.root / 'child-result.json')}
        with patch.object(pipeline, 'main', return_value={'status': 'SOFTWARE_FIXTURE'}) as called:
            update.child(request)
        argv = called.call_args.args[0]
        index = argv.index('--dnam-idat-provenance')
        self.assertEqual(argv[index + 1], str(self.wrapper))

    def test_direct_pipeline_rechecks_transitive_sources_on_every_exit(self):
        from src import pipeline_2026_09_10 as pipeline
        args = pipeline.parser().parse_args(['--study', 'GSE1', '--modality', 'human_dnam',
                    '--matrix', self.entry['sources']['matrix']['path'], '--output', str(self.root / 'pipeline'),
                    '--dnam-idat-provenance', str(self.wrapper)])
        def changed(*unused):
            Path(self.files['raw:GSM1:Grn']['path']).write_bytes(b'changed during pipeline software fixture')
            return {'status': 'SOFTWARE_FIXTURE'}
        with patch.object(pipeline, '_run', side_effect=changed):
            with self.assertRaisesRegex(ValueError, 'source changed during pipeline execution'):
                pipeline.run(args)

    def test_direct_pipeline_records_and_passes_upstream_lineage(self):
        from src import pipeline_2026_09_10 as pipeline
        output = self.root / 'pipeline'
        args = pipeline.parser().parse_args(['--study', 'GSE1', '--modality', 'human_dnam',
                    '--matrix', self.entry['sources']['matrix']['path'], '--output', str(output),
                    '--dnam-idat-provenance', str(self.wrapper), '--until', 'preprocess'])
        declared = self.helper.declaration(self.wrapper)
        metadata = {'sample_ids': ['GSM1'], 'metadata': 'software-fixture.csv',
                    'series_source': self.files['series']['path'], 'study_context': {},
                    'source_sample_records': self.files['samples']['path']}
        raw = {'path': self.files['matrix']['path'], 'idat_provenance': declared}
        configs = []
        class Store:
            def __init__(self, folder, config):
                self.folder = folder
                configs.append(config)
            def run(self, stage, function):
                folder = self.folder / stage
                folder.mkdir(parents=True, exist_ok=True)
                return function(folder)
        with patch.object(pipeline, 'StageStore', Store), patch.object(pipeline, 'checked_assets'), \
                patch.object(pipeline, 'acquire_metadata', return_value=metadata), \
                patch.object(pipeline, 'source_packet', return_value={}), \
                patch.object(pipeline, 'build_client', return_value=object()), \
                patch.object(pipeline, 'assess_design', return_value={}), \
                patch.object(pipeline, 'acquire_dnam', return_value=raw) as acquire, \
                patch.object(pipeline, 'preprocess_dnam', return_value={'status': 'SOFTWARE_FIXTURE'}) as preprocess, \
                patch.object(pipeline.subprocess, 'check_output', return_value='R SOFTWARE FIXTURE'), \
                patch('importlib.metadata.version', return_value='SOFTWARE_FIXTURE'):
            self.assertEqual(pipeline.run(args)['status'], 'SOFTWARE_FIXTURE')
        self.assertEqual(acquire.call_args.kwargs['idat_provenance'], self.wrapper)
        self.assertEqual(preprocess.call_args.kwargs['idat_provenance'], declared)
        self.assertEqual(configs[0]['dnam_idat_provenance_sources']['raw:GSM1:Grn'], self.files['raw:GSM1:Grn'])


if __name__ == "__main__":
    unittest.main()
