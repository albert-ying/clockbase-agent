"""Software fixtures for read-only corpus status snapshots."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import tempfile
import unittest

from audit_system.corpus_status_2026_09_20 import (
    StatusIntegrityError, inspect_run, write_snapshot,
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class CorpusStatusTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name); self.run = self.root / "run"; self.run.mkdir()
        self.worklist = self.root / "worklist"; self.worklist.mkdir()
        self.tasks = [self.task("GSE1"), self.task("GSE2")]
        artifact_rows = [
            {"artifact_id": "artifact:GSE1", "path": "one.csv", "sha256": "1" * 64, "parsed_rows": 3},
            {"artifact_id": "artifact:GSE2", "path": "two.csv", "sha256": "2" * 64, "parsed_rows": 1},
            {"artifact_id": "artifact:orphan", "path": "orphan.csv", "sha256": "3" * 64, "parsed_rows": 1},
        ]
        rows = []
        for artifact in artifact_rows:
            for ordinal in range(2, artifact["parsed_rows"] + 2):
                coordinates = [artifact["path"], artifact["sha256"], ordinal]
                rows.append({"artifact_id": artifact["artifact_id"], "artifact_sha256": artifact["sha256"],
                    "logical_csv_record": ordinal, "row_id": "row:" + hashlib.sha256(json.dumps(
                        coordinates, ensure_ascii=False, sort_keys=True, separators=(",", ":")
                    ).encode()).hexdigest(), "identity_status": "PASS" if artifact["artifact_id"] != "artifact:orphan" else "UNRESOLVED"})
        outputs = {}
        for name, values in (("artifact-worklist.jsonl.gz", artifact_rows), ("clock-row-worklist.jsonl.gz", rows)):
            with gzip.open(self.worklist / name, "wt") as handle:
                for value in values:
                    handle.write(json.dumps(value) + "\n")
            outputs[name] = sha(self.worklist / name)
        summary = {"denominators": self.denominators(), "outputs": outputs}
        self.write(self.worklist / "summary.json", summary)
        implementation = self.run / "implementation" / "fixture.py"
        implementation.parent.mkdir(); implementation.write_text("frozen\n")
        config = {
            "schema": "clockbase-corpus-audit-v1",
            "config": {"max_total_cost_usd": 4.0, "max_total_calls": 20},
            "implementation_sha256": {"fixture.py": sha(implementation)},
            "worklist": str(self.worklist),
            "worklist_summary_sha256": sha(self.worklist / "summary.json"),
            "denominators": self.denominators(),
        }
        self.write(self.run / "config.json", config)
        (self.run / "config.sha256").write_text(sha(self.run / "config.json") + "\n")
        with gzip.open(self.run / "queue.jsonl.gz", "wt") as handle:
            for task in self.tasks:
                handle.write(json.dumps(task) + "\n")
        self.write(self.run / "queue-receipt.json", {"studies": 2, "sha256": sha(self.run / "queue.jsonl.gz")})
        self.database()

    @staticmethod
    def denominators():
        return {"artifacts": 3, "zero_row_artifacts": 1, "parsed_logical_rows": 5,
                "canonical_candidate_study_by_biological_modality": 2,
                "identity_bound_rows": 4, "identity_unresolved_rows": 1}

    @staticmethod
    def task(study):
        return {"study_key": f"human_rna:{study}", "study_id": study, "modality": "human_rna",
                "artifact_ids": [f"artifact:{study}"], "source_files": {}}

    @staticmethod
    def write(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value, sort_keys=True))

    def database(self):
        db = sqlite3.connect(self.run / "budget.sqlite3")
        db.execute("CREATE TABLE config (id INTEGER PRIMARY KEY, cost REAL, calls INTEGER)")
        db.execute("INSERT INTO config VALUES (1,4.0,20)")
        db.execute("CREATE TABLE reservations (id TEXT PRIMARY KEY, reserved REAL, charged REAL, status TEXT, updated TEXT, binding_sha256 TEXT)")
        db.commit(); db.close()

    def add_budget(self):
        db = sqlite3.connect(self.run / "budget.sqlite3")
        db.execute("INSERT INTO reservations VALUES (?,?,?,?,?,?)",
                   ("a" * 64, .25, .10, "USAGE_ESTIMATE", "2026-09-20T00:00:00Z", "b" * 64))
        db.execute("INSERT INTO reservations VALUES (?,?,?,?,?,?)",
                   ("c" * 64, .30, .30, "RESERVED_OR_UNCERTAIN", "2026-09-20T00:00:01Z", "d" * 64))
        db.commit(); db.close()

    def terminal(self, task, *, reviewed=1, calculations=2, unresolved=1):
        study = self.run / "studies" / task["study_key"].replace(":", "-")
        study.mkdir(parents=True, exist_ok=True)
        output = study / "stage.json"; output.write_text("{}\n")
        artifact = {"path": "one.csv" if task["study_id"] == "GSE1" else "two.csv",
                    "sha256": "1" * 64 if task["study_id"] == "GSE1" else "2" * 64}
        rows = ["R" + hashlib.sha256(json.dumps(
                    [artifact["path"], artifact["sha256"], i + 2], sort_keys=True, allow_nan=False
                ).encode()).hexdigest()[:24] for i in range(reviewed)]
        value = {"study_key": task["study_key"], "study_id": task["study_id"], "modality": task["modality"],
                 "artifact_ids": task["artifact_ids"], "status": "REVIEWED",
                 "scientific_acceptance": "CANNOT VERIFY", "reviewed_rows": reviewed,
                 "calculations": calculations, "unresolved_calculations": unresolved,
                 "accepted_scientific_repairs": 0, "row_ids": rows,
                 "identity_status_counts": {"PASS": reviewed}, "output_sha256": {"stage.json": sha(output)}}
        self.write(study / "terminal.json", value)
        (study / "terminal.sha256").write_text(sha(study / "terminal.json") + "\n")
        return value

    def progress(self, *events, partial_tail=None):
        data = b"".join((json.dumps(event) + "\n").encode() for event in events)
        if partial_tail is not None:
            data += partial_tail
        (self.run / "progress.jsonl").write_bytes(data)

    def event(self, task, status="REVIEWED", terminal=None):
        return {"at_utc": "2026-09-20T00:00:00Z", "study_key": task["study_key"], "status": status,
                "reviewed_rows": terminal.get("reviewed_rows", 0) if terminal else 0,
                "terminal_sha256": sha(self.run / "studies" / task["study_key"].replace(":", "-") / "terminal.json")
                if terminal else None}

    def test_empty_pending_queue_and_immutable_snapshot(self):
        report = inspect_run(self.run)
        self.assertEqual(report["progress"]["completed_terminal_studies"], 0)
        self.assertEqual(report["coverage"]["pending_study_tasks"], 2)
        self.assertEqual(report["selected_studies"], 0)
        path = self.run / "reports" / "status-0001.json"
        binding = write_snapshot(self.run, path, report)
        self.assertEqual(binding["sha256"], sha(path))
        self.assertEqual(json.loads(path.read_text())["status"], "IN_PROGRESS")
        with self.assertRaises(FileExistsError):
            write_snapshot(self.run, path, report)
        relative = write_snapshot(self.run, Path("reports/status-0002.json"), report)
        self.assertTrue(Path(relative["path"]).is_relative_to(self.run))
        with self.assertRaises(ValueError):
            write_snapshot(self.run, self.root / "outside.json", report)

    def test_consistent_subset_counts_without_inflation(self):
        terminal = self.terminal(self.tasks[0], reviewed=2, calculations=3, unresolved=1)
        self.progress(self.event(self.tasks[0], terminal=terminal), partial_tail=b'{"study_key"')
        report = inspect_run(self.run)
        self.assertEqual((report["selected_studies"], report["unique_reviewed_rows"]), (1, 2))
        self.assertEqual((report["calculations"], report["unresolved_calculations"]), (3, 1))
        self.assertEqual(report["coverage"]["pending_study_tasks"], 1)
        self.assertEqual(report["coverage"]["rows_without_candidate_study"], 1)
        self.assertTrue(report["progress"]["ignored_partial_progress_tail"])

    def test_partial_terminal_pair_remains_pending(self):
        study = self.run / "studies" / "human_rna-GSE1"; study.mkdir(parents=True)
        (study / "terminal.json").write_text("{partially written")
        report = inspect_run(self.run)
        self.assertEqual(report["selected_studies"], 0)
        self.assertEqual(report["progress"]["partial_terminal_pairs"], 1)
        (study / "terminal.sha256").write_text("abc")
        report = inspect_run(self.run)
        self.assertEqual(report["selected_studies"], 0)
        self.assertEqual(report["progress"]["partial_terminal_reasons"], {"INCOMPLETE_SIDECAR": 1})

    def test_terminal_or_queue_tamper_fails_closed(self):
        self.terminal(self.tasks[0])
        path = self.run / "studies" / "human_rna-GSE1" / "terminal.json"
        path.write_text(path.read_text() + " ")
        with self.assertRaisesRegex(StatusIntegrityError, "Terminal hash mismatch"):
            inspect_run(self.run)
        path = self.run / "queue.jsonl.gz"
        path.write_bytes(path.read_bytes() + b"x")
        with self.assertRaisesRegex(StatusIntegrityError, "queue hash"):
            inspect_run(self.run)

    def test_duplicate_resume_and_controller_failure_do_not_inflate(self):
        terminal = self.terminal(self.tasks[0])
        failure = self.event(self.tasks[0], status="CONTROLLER_FAILURE")
        success = self.event(self.tasks[0], terminal=terminal)
        unresolved_failure = self.event(self.tasks[1], status="CONTROLLER_FAILURE")
        self.progress(failure, success, success, unresolved_failure)
        report = inspect_run(self.run)
        self.assertEqual(report["selected_studies"], 2)
        self.assertEqual(report["progress"]["completed_terminal_studies"], 1)
        self.assertEqual(report["progress"]["controller_failure_studies"], 1)
        self.assertEqual(report["progress"]["recovered_controller_failure_studies"], 1)
        self.assertEqual(report["progress"]["duplicate_or_resume_progress_events"], 2)
        self.assertEqual(report["unique_reviewed_rows"], 1)
        self.assertEqual(report["status_counts"], {"CONTROLLER_FAILURE": 1, "REVIEWED": 1})

    def test_conflicting_resume_hashes_fail_closed(self):
        terminal = self.terminal(self.tasks[0])
        event = self.event(self.tasks[0], terminal=terminal)
        other = dict(event, terminal_sha256="f" * 64)
        self.progress(event, other)
        with self.assertRaisesRegex(StatusIntegrityError, "conflicting terminals"):
            inspect_run(self.run)

    def test_budget_is_read_only_and_conservative(self):
        self.add_budget()
        before = (sha(self.run / "budget.sqlite3"), (self.run / "budget.sqlite3").stat().st_mtime_ns)
        report = inspect_run(self.run)
        after = (sha(self.run / "budget.sqlite3"), (self.run / "budget.sqlite3").stat().st_mtime_ns)
        self.assertEqual(before, after)
        self.assertEqual(report["budget"]["call_reservations"], 2)
        self.assertAlmostEqual(report["budget"]["estimated_or_reserved_usd"], .40)
        self.assertIsNone(report["budget"]["billed_cost_usd"])

    def test_unresolved_targets_are_disjoint_from_pass_calculations(self):
        terminal = self.terminal(self.tasks[0], reviewed=1, calculations=0, unresolved=2)
        self.progress(self.event(self.tasks[0], terminal=terminal))
        report = inspect_run(self.run)
        self.assertEqual((report["calculations"], report["unresolved_calculations"]), (0, 2))

    def test_early_unresolved_terminal_can_have_no_row_identity_summary(self):
        value = self.terminal(self.tasks[0], reviewed=0, calculations=0, unresolved=0)
        study = self.run / "studies" / "human_rna-GSE1"
        value["status"] = "UNRESOLVED"; value.pop("row_ids"); value.pop("identity_status_counts")
        self.write(study / "terminal.json", value)
        (study / "terminal.sha256").write_text(sha(study / "terminal.json"))
        report = inspect_run(self.run)
        self.assertEqual(report["status_counts"], {"UNRESOLVED": 1})


if __name__ == "__main__":
    unittest.main()
