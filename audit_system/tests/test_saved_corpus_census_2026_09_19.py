import csv
import gzip
import importlib.util
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "audit_legacy_result_inputs_2026_09_18.py"
SPEC = importlib.util.spec_from_file_location("saved_census", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)

HEADER = [
    "gse_id", "comparison_name", "condition_name", "condition_category", "clock",
    "effect_size", "p_value", "trust_score", "direction",
]


def write_csv(path, rows, header=HEADER):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def plan_file(tmp_path):
    path = tmp_path / "PLAN.md"
    path.write_text("frozen saved-corpus fixture plan\n")
    return path


def check_saved_csv_diagnostics_keep_all_findings(tmp_path):
    path = tmp_path / "result_mm_rna" / "GSE1" / "aging_analysis_sonnet37.csv"
    rows = [
        ["GSE1", "cmp", "condition", "cat", "clock", "1", "0.2", "0.5", "1"],
        ["GSE1", "cmp", "condition", "cat", "clock", "1", "0.2", "0.5", "1"],
        ["GSE1", "cmp", "condition", "cat", "clock", "2", "", "0.5", "-1"],
        ["GSE2", "other", "condition", "cat", "clock2", "nan", "1.2", "inf", "0"],
    ]
    write_csv(path, rows)
    record = AUDIT._saved_csv_audit(path, tmp_path, "result_mm_rna", {"wall_seconds": 60})
    codes = {issue["code"] for issue in record["issues"]}
    assert record["status"] == "FAIL"
    assert {"PARENT_ID_MISMATCH", "NONFINITE_OR_NONNUMERIC", "P_VALUE_OUT_OF_BOUNDS", "MISSING_STATISTIC"} <= codes
    assert len(record["full_row_duplicate_groups"]) == 1
    assert len(record["conflicting_literal_key_groups"]) == 1
    assert record["literal_accessions"] == ["GSE1", "GSE2"]


def check_saved_corpus_census_includes_all_roots_and_both_filename_forms(tmp_path):
    roots = tuple(AUDIT.SAVED_ROOTS)
    for index, root in enumerate(roots):
        for suffix in ("", "_sonnet37"):
            write_csv(tmp_path / root / f"GSE{index + 1}" / f"aging_analysis{suffix}.csv", [
                [f"GSE{index + 1}", "cmp", "condition", "cat", "clock", "1", "0.2", "0.5", "1"],
            ])
    output = tmp_path / "out"
    AUDIT.execute_saved_corpus(tmp_path, plan_file(tmp_path), output, limits={"files": 100, "bytes": 10_000, "wall_seconds": 60})
    summary = json.loads((output / "summary.json").read_text())
    assert summary["status"] == "PASS"
    assert summary["denominators"] == {"artifacts": 10, "rows": 10, "literal_accessions": 5, "accession_by_modality": 5}
    ledger = [json.loads(line) for line in (output / "artifact-ledger.jsonl").read_text().splitlines()]
    assert {Path(item["path"]).name for item in ledger} == {"aging_analysis.csv", "aging_analysis_sonnet37.csv"}
    assert all(item["status"] == "PASS" for item in ledger)


def check_lfs_pointer_is_registered_as_unavailable_without_parsing(tmp_path):
    path = tmp_path / "result_hs_rna" / "GSE9" / "aging_analysis.csv"
    path.parent.mkdir(parents=True)
    path.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 10\n")
    record = AUDIT._saved_csv_audit(path, tmp_path, "result_hs_rna", {"wall_seconds": 60})
    assert record["status"] == "CANNOT VERIFY"
    assert record["reason"] == "LFS_POINTER"
    assert record["rows"] == 0


def check_census_refuses_overwrite_and_detects_input_mutation(tmp_path):
    path = tmp_path / "result_mm_rna" / "GSE1" / "aging_analysis.csv"
    write_csv(path, [["GSE1", "cmp", "condition", "cat", "clock", "1", "0.2", "0.5", "1"]])
    plan = plan_file(tmp_path)
    output = tmp_path / "out"
    AUDIT.execute_saved_corpus(tmp_path, plan, output, limits={"files": 10, "bytes": 10_000, "wall_seconds": 60})
    with unittest.TestCase().assertRaises(FileExistsError):
        AUDIT.execute_saved_corpus(tmp_path, plan, output, limits={"files": 10, "bytes": 10_000, "wall_seconds": 60})
    original = AUDIT._saved_csv_audit
    mutated = {"done": False}

    def mutate_after_audit(*args, **kwargs):
        record = original(*args, **kwargs)
        if not mutated["done"]:
            with path.open("a") as handle:
                handle.write("GSE1,cmp2,condition,cat,clock,2,0.3,0.6,-1\n")
            mutated["done"] = True
        return record

    with patch.object(AUDIT, "_saved_csv_audit", mutate_after_audit):
        with unittest.TestCase().assertRaisesRegex(ValueError, "Source changed during saved-corpus census"):
            AUDIT.execute_saved_corpus(tmp_path, plan, tmp_path / "mutated-out",
                                       limits={"files": 10, "bytes": 10_000, "wall_seconds": 60})


class SavedCorpusTests(unittest.TestCase):
    """Run the census fixtures through the repository's standard test runner."""

    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.tmp_path = Path(self.directory.name)

    def test_diagnostics(self):
        check_saved_csv_diagnostics_keep_all_findings(self.tmp_path)

    def test_roots_and_filenames(self):
        check_saved_corpus_census_includes_all_roots_and_both_filename_forms(self.tmp_path)

    def test_lfs_pointer(self):
        check_lfs_pointer_is_registered_as_unavailable_without_parsing(self.tmp_path)

    def test_overwrite_and_mutation(self):
        check_census_refuses_overwrite_and_detects_input_mutation(self.tmp_path)

    def worklist_fixture(self):
        root = self.tmp_path
        row = ["GSE1", "cmp", "condition", "cat", "clock", "1", "0.2", "0.5", "1"]
        write_csv(root / "result_mm_rna/GSE1/aging_analysis.csv", [row, row])
        write_csv(root / "result_mm_rna/GSE1/aging_analysis_sonnet37.csv", [row])
        write_csv(root / "result_mm_rna/GSE2/aging_analysis.csv", [row])
        empty = root / "result_mm_rna/GSE3/aging_analysis.csv"
        empty.parent.mkdir(parents=True)
        empty.write_text("")
        ragged = root / "result_mm_rna/GSE4/aging_analysis.csv"
        write_csv(ragged, [["GSE4", "ragged"], ["GSE4", "cmp", "condition", "cat", "clock", "2", "0.3", "1", "1"]])
        malformed = root / "result_mm_rna/GSE5/aging_analysis.csv"
        write_csv(malformed, [["GSE5", "cmp", "condition", "cat", "clock", "2", "0.3", "1", "1"]])
        with malformed.open("a") as handle:
            handle.write('"unterminated\n')
        pointer = root / "result_mm_rna/GSE6/aging_analysis.csv"
        pointer.parent.mkdir(parents=True)
        pointer.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 12\n")
        geo = root / "geo"
        for role in ("metadata", "clock_result"):
            write_csv(geo / "GEO_mm_rna" / role / "GSE1.csv", [["GSM1", "control"]], ["GSM", "group"])
        census = root / "census"
        AUDIT.execute_saved_corpus(root, plan_file(root), census, roots=("result_mm_rna",),
                                   limits={"files": 100, "bytes": 100000, "wall_seconds": 60})
        return census / "artifact-ledger.jsonl", geo

    def test_worklist_preserves_all_variants_rows_empty_and_malformed_artifacts(self):
        ledger, geo = self.worklist_fixture()
        result = AUDIT.execute_worklist(self.tmp_path, ledger, geo, self.tmp_path / "worklist")
        out = self.tmp_path / "worklist"
        artifacts = list(AUDIT._read_gz_jsonl(out / "artifact-worklist.jsonl.gz"))
        rows = list(AUDIT._read_gz_jsonl(out / "clock-row-worklist.jsonl.gz"))
        comparisons = list(AUDIT._read_gz_jsonl(out / "comparison-worklist.jsonl.gz"))
        self.assertEqual(result["denominators"]["artifacts"], 7)
        self.assertEqual(result["denominators"]["zero_row_artifacts"], 2)
        self.assertEqual(result["denominators"]["parsed_logical_rows"], 7)
        self.assertEqual(result["denominators"]["literal_comparisons_within_artifacts"], 6)
        self.assertEqual(len({r["row_id"] for r in rows}), 7)
        # Identical rows in one file retain different IDs and one shared row digest.
        duplicate = [r for r in rows if r["artifact_id"] == artifacts[0]["artifact_id"]]
        self.assertEqual(len({r["row_sha256"] for r in duplicate}), 1)
        self.assertEqual(len({r["row_id"] for r in duplicate}), 2)
        self.assertTrue(any(r["parse_error"] and "CSV" not in r["parse_error"] for r in artifacts))
        self.assertEqual(sum(r["parse_error"] == "LFS_POINTER" for r in artifacts), 1)
        self.assertTrue(any("PARENT_ROW_STUDY_ID_CONTRADICTION" in r["identity_issues"] for r in rows))
        self.assertTrue(any("RAGGED_ROW" in r["identity_issues"] for r in rows))
        self.assertTrue(all(r["scientific_acceptance"] == "CANNOT VERIFY" for r in rows + comparisons + artifacts))
        self.assertTrue(all(set(r["stages"].values()) == {"NOT_YET_AUDITED"} for r in rows + comparisons))
        self.assertEqual(result["denominators"]["statistical_rerun_ready_comparisons"], 0)
        self.assertEqual(result["denominators"]["source_design_preparation_ready_comparisons"], 2)
        self.assertTrue(all("p_value" not in row and "effect_size" not in row for row in rows))

    def test_worklist_refuses_changed_census_csv_and_duplicate_census_paths(self):
        ledger, geo = self.worklist_fixture()
        source = self.tmp_path / "result_mm_rna/GSE1/aging_analysis.csv"
        original = source.read_bytes()
        source.write_bytes(original + b"\n")
        with self.assertRaisesRegex(ValueError, "changed since the census"):
            AUDIT.execute_worklist(self.tmp_path, ledger, geo, self.tmp_path / "changed")
        source.write_bytes(original)
        lines = ledger.read_text().splitlines()
        ledger.write_text("\n".join(lines + [lines[0]]) + "\n")
        summary = json.loads((ledger.parent / "summary.json").read_text())
        summary["outputs"][ledger.name] = AUDIT.sha256(ledger)
        summary["denominators"]["artifacts"] += 1
        (ledger.parent / "summary.json").write_text(json.dumps(summary))
        with self.assertRaisesRegex(ValueError, "Duplicate artifact paths"):
            AUDIT.execute_worklist(self.tmp_path, ledger, geo, self.tmp_path / "duplicates")

    def test_worklist_resume_verifies_inputs_outputs_and_missing_path_state(self):
        ledger, geo = self.worklist_fixture()
        output = self.tmp_path / "worklist"
        AUDIT.execute_worklist(self.tmp_path, ledger, geo, output)
        before = {p.name: p.read_bytes() for p in output.iterdir()}
        AUDIT.execute_worklist(self.tmp_path, ledger, geo, output, resume=True)
        self.assertEqual(before, {p.name: p.read_bytes() for p in output.iterdir()})
        with self.assertRaises(FileExistsError):
            AUDIT.execute_worklist(self.tmp_path, ledger, geo, output)
        missing = geo / "GEO_mm_rna/metadata/GSE2.csv"
        write_csv(missing, [["GSM2", "control"]], ["GSM", "group"])
        with self.assertRaisesRegex(ValueError, "availability changed"):
            AUDIT.execute_worklist(self.tmp_path, ledger, geo, output, resume=True)
        missing.unlink()
        metadata = geo / "GEO_mm_rna/metadata/GSE1.csv"
        original = metadata.read_bytes()
        metadata.write_bytes(original + b"\n")
        with self.assertRaisesRegex(ValueError, "Input payload changed"):
            AUDIT.execute_worklist(self.tmp_path, ledger, geo, output, resume=True)
        metadata.write_bytes(original)
        (output / "comparison-worklist.jsonl.gz").write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError, "Worklist output changed"):
            AUDIT.execute_worklist(self.tmp_path, ledger, geo, output, resume=True)

    def test_worklist_reuses_only_bound_family_source_identity(self):
        ledger, geo = self.worklist_fixture()
        family = self.tmp_path / "families"
        enrollment = family / "study-enrollment.csv"
        write_csv(enrollment, [["GSE1", "family:test", "CANNOT VERIFY", "DO_NOT_IMPORT"]],
                  ["study_id", "family_id", "source_linkage_status", "reference_label"])
        manifest_path = family / "manifest.json"
        manifest_path.write_text(json.dumps({"source_coverage_complete": False,
            "scope": "Source fixture universe", "output_sha256": {enrollment.name: AUDIT.sha256(enrollment)}}))
        output = self.tmp_path / "worklist"
        result = AUDIT.execute_worklist(self.tmp_path, ledger, geo, output, family_manifest=manifest_path)
        studies = list(AUDIT._read_gz_jsonl(output / "study-inputs.jsonl.gz"))
        self.assertEqual(studies[0]["source_family"]["family_id"], "family:test")
        self.assertFalse(result["family_binding"]["leakage_safe_split_established"])
        self.assertNotIn("DO_NOT_IMPORT", str(studies))
        enrollment.write_text(enrollment.read_text() + "\n")
        with self.assertRaisesRegex(ValueError, "not bound"):
            AUDIT.execute_worklist(self.tmp_path, ledger, geo, self.tmp_path / "bad-family", family_manifest=manifest_path)

    def test_worklist_rechecks_every_source_after_preparation(self):
        ledger, geo = self.worklist_fixture()
        original = AUDIT._write_gz_jsonl
        def mutate_after_rows(path, rows):
            original(path, rows)
            if path.name == "study-inputs.jsonl.gz":
                with (geo / "GEO_mm_rna/metadata/GSE1.csv").open("a") as handle:
                    handle.write("GSM2,treatment\n")
        with patch.object(AUDIT, "_write_gz_jsonl", mutate_after_rows):
            with self.assertRaisesRegex(ValueError, "Input payload changed"):
                AUDIT.execute_worklist(self.tmp_path, ledger, geo, self.tmp_path / "mutated")
        self.assertEqual(json.loads((self.tmp_path / "mutated/summary.json").read_text())["status"], "FAIL")


if __name__ == "__main__":
    unittest.main()
