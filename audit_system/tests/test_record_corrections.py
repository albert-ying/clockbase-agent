"""Source-linked regression checks and explicitly constructed boundary fixtures."""
import copy
import json
from pathlib import Path
import subprocess
import unittest

from audit_system.record_corrections import (
    CorrectionError, csv_records, digest, local_path, patch_csv, prepare,
    validate_probability_change,
)
from audit_system.tests import require_development_artifact

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "paper/revision_2026-09-07/source/corrections/2026-09-10-PLX5622-record-corrections.json"
NOVELTY_MANIFEST = ROOT / "paper/revision_2026-09-07/source/corrections/2026-09-15-PLX5622-novelty-priority.json"
CURRENT_MANIFEST = ROOT / "paper/revision_2026-09-07/source/corrections/2026-09-16-PLX5622-component-review.json"


class ProbabilityCorrectionTests(unittest.TestCase):
    """Adversarial software fixtures; these are not empirical observations."""

    def test_p_change_cannot_silently_inherit_adjusted_p(self):
        row = {"pvalue": "0.01", "fdr": "0.5"}
        with self.assertRaisesRegex(CorrectionError, "explicit FDR disposition"):
            validate_probability_change(row, {"pvalue": "0.1"}, 2)
        # Lexical-only change preserves the numeric input identity.
        validate_probability_change(row, {"pvalue": "0.010"}, 2)

    def test_impossible_q_and_withheld_p_are_rejected(self):
        for field in ("pvalue", "p_value"):
            row = {field: "0.1", "fdr": "0.2"}
            with self.subTest(field=field):
                with self.assertRaisesRegex(CorrectionError, "below raw P"):
                    validate_probability_change(row, {"fdr": "0.01"}, 2)
                with self.assertRaisesRegex(CorrectionError, "Raw P withheld"):
                    validate_probability_change(row, {field: ""}, 2)

    def test_withholding_keeps_unknown_source_fields(self):
        row = {"pvalue": "unresolved original field", "fdr": "0.01"}
        validate_probability_change(row, {"fdr": ""}, 2)
        validate_probability_change({"pvalue": "0.7", "fdr": "0.00001"},
                                    {"final_score": "78"}, 2)

    def test_finite_probability_limits_and_exact_zeros(self):
        for field in ("pvalue", "p_value", "fdr"):
            for value in ("NaN", "Infinity", "-0.1", "1.1", "invalid"):
                with self.subTest(field=field, value=value), self.assertRaises(CorrectionError):
                    validate_probability_change({field: "0.1"}, {field: value}, 2)
        validate_probability_change({"pvalue": "0.1", "fdr": "0.2"},
                                    {"pvalue": "0", "fdr": "0"}, 2)

    def test_empirical_inherited_contradiction_can_only_be_withheld_or_explicitly_replaced(self):
        source = ROOT / "result/main_results_result_mm_rna_outlier_marked.csv"
        require_development_artifact(source)
        rows = csv_records(source.read_bytes())[1]
        row = next(row for row, _, _ in rows if row["gse_id"] == "GSE65770"
                   and row["comparison_name"] == "Q140 vs WT at 6 month")
        with self.assertRaisesRegex(CorrectionError, "below raw P"):
            validate_probability_change(row, {"fdr": row["fdr"]}, 2)
        validate_probability_change(row, {"fdr": ""}, 2)

    def test_withholding_cannot_leave_logarithms_or_selection_flags(self):
        row = {"pvalue": "0.01", "fdr": "0.04", "neg_log10_pvalue": "2",
               "neg_log10_fdr": "1.3979400086720376", "is_significant": "True"}
        for changes in ({"fdr": ""}, {"fdr": "", "neg_log10_fdr": ""},
                        {"pvalue": "", "fdr": "", "neg_log10_fdr": "", "is_significant": ""}):
            with self.subTest(changes=changes), self.assertRaisesRegex(CorrectionError, "Stale derived"):
                validate_probability_change(row, changes, 2)
        validate_probability_change(row, {key: "" for key in row}, 2)

    def test_explicit_derived_fields_must_match_the_changed_probability(self):
        row = {"pvalue": "0.01", "fdr": "0.04", "neg_log10_pvalue": "2",
               "neg_log10_fdr": "1.3979400086720376", "is_significant": "True"}
        valid = {"pvalue": "0.1", "fdr": "0.2", "neg_log10_pvalue": "1",
                 "neg_log10_fdr": "0.6989700043360189", "is_significant": "False"}
        validate_probability_change(row, valid, 2)
        for changed in ({"is_significant": "True"}, {"neg_log10_pvalue": "2"},
                        {"neg_log10_fdr": "NaN"}, {"neg_log10_fdr": "-1"},
                        {"pvalue": "0", "neg_log10_pvalue": "300"}):
            with self.subTest(changed=changed), self.assertRaises(CorrectionError):
                validate_probability_change(row, valid | changed, 2)


class RecordCorrectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        for path in (MANIFEST, NOVELTY_MANIFEST, CURRENT_MANIFEST):
            require_development_artifact(path)
        cls.manifest = json.loads(MANIFEST.read_text())
        cls.inputs = []
        cls.outputs = []
        cls.rows = []
        for declaration in cls.manifest["files"]:
            source = ROOT / declaration["source"]
            require_development_artifact(source)
            raw = source.read_bytes()
            if digest(raw) != declaration["before_sha256"]:
                raw = subprocess.check_output(
                    ["git", "show", cls.manifest["source_commit"] + ":" + declaration["source"]],
                    cwd=ROOT,
                )
            cls.inputs.append(raw)
            after = patch_csv(raw, declaration)
            cls.outputs.append(after)
            cls.rows.append((csv_records(raw)[1], csv_records(after)[1]))

    def test_declared_source_and_output_hashes_match(self):
        for declaration, raw, after in zip(self.manifest["files"], self.inputs, self.outputs):
            self.assertEqual(digest(raw), declaration["before_sha256"])
            self.assertEqual(digest(after), declaration["after_sha256"])

    def test_all_unlisted_rows_and_fields_are_preserved(self):
        for declaration, (before, after) in zip(self.manifest["files"], self.rows):
            operations = {o["record_number"]: o for o in declaration["operations"]}
            self.assertEqual(len(before), len(after))
            for number, (old, new) in enumerate(zip(before, after), 2):
                expected = old[0] | operations.get(number, {}).get("changes", {})
                self.assertEqual(expected, new[0])

    def test_ad_score_is_neutralized_without_manufacturing_top_rank(self):
        for (_, after) in self.rows:
            field = "case_id" if "case_id" in after[0][0] else "dataset"
            targets = [r for r, _, _ in after if r[field] == "GSE134151-10407"]
            self.assertEqual(len(targets), 1)
            self.assertEqual(float(targets[0]["final_score"]), 61.0)
            self.assertLess(float(targets[0]["final_score"]), 100.0)
            if "q0_plausibility" in targets[0]:
                self.assertEqual(targets[0]["q0_plausibility"], "0")

    def test_both_age_confounded_comparisons_withheld_in_both_tables(self):
        for _, after in self.rows[:2]:
            targets = [r for r, _, _ in after if r["gse_id"] == "GSE166236"
                       and r["comparison_name"].startswith("PLX5622 effect on")]
            self.assertEqual(len(targets), 4)
            for row in targets:
                for field in ["pvalue", "effect_size", "direction", "log2FoldChange", "fdr", "final_score"]:
                    self.assertEqual(row[field], "")
                self.assertEqual(row["audit_decision"], "FAIL_DESIGN_CONFOUNDING")

    def test_duplicates_are_not_reactivated(self):
        for before, after in self.rows[:2]:
            for (old, _, _), (new, _, _) in zip(before, after):
                if old["outlier_status"].startswith("Excluded"):
                    self.assertEqual(new["outlier_status"], old["outlier_status"])
                    for field in ["pvalue", "effect_size", "direction", "final_score"]:
                        self.assertEqual(new[field], old[field])

    def test_positive_and_unsigned_effect_encodings_are_preserved(self):
        for before, after in self.rows[:2]:
            for (old, _, _), (new, _, _) in zip(before, after):
                if old["gse_id"] in {"GSE171321", "GSE132516"}:
                    for field in ["effect_size", "direction", "log2FoldChange", "pvalue"]:
                        self.assertEqual(new[field], old[field])
                    if old["gse_id"] == "GSE132516":
                        self.assertGreater(float(new["effect_size"]), 0)
                        self.assertLess(float(new["log2FoldChange"]), 0)

    def test_modified_input_and_wrong_record_identity_are_rejected(self):
        declaration = self.manifest["files"][0]
        with self.assertRaises(CorrectionError):
            patch_csv(self.inputs[0] + b"\n", declaration)
        changed = copy.deepcopy(declaration)
        changed["operations"][0]["before"]["gse_id"] = "not-the-source-identifier"
        with self.assertRaises(CorrectionError):
            patch_csv(self.inputs[0], changed)

    def test_duplicate_targets_and_identity_changes_are_rejected(self):
        declaration = copy.deepcopy(self.manifest["files"][0])
        declaration["operations"].append(copy.deepcopy(declaration["operations"][0]))
        with self.assertRaises(CorrectionError):
            patch_csv(self.inputs[0], declaration)
        declaration = copy.deepcopy(self.manifest["files"][0])
        declaration["operations"][0]["changes"]["gse_id"] = "different-identifier"
        with self.assertRaises(CorrectionError):
            patch_csv(self.inputs[0], declaration)

    def test_applied_bundle_is_idempotent_and_paths_stay_in_repository(self):
        current = json.loads(CURRENT_MANIFEST.read_text())
        plans = prepare(ROOT, current)
        self.assertTrue(all(p["state"] == "already_applied" for p in plans))
        with self.assertRaises(CorrectionError):
            local_path(ROOT, "../outside.csv")
        with self.assertRaises(CorrectionError):
            local_path(ROOT, "/absolute.csv")

    def test_priority_reassessment_preserves_statistics_and_records_uncertainty(self):
        current = json.loads(NOVELTY_MANIFEST.read_text())
        self.assertEqual(current["linked_record_count"], 3)
        for declaration in current["files"]:
            original = subprocess.check_output(
                ["git", "show", current["source_commit"] + ":" + declaration["source"]],
                cwd=ROOT,
            )
            after = patch_csv(original, declaration)
            self.assertEqual(digest(after), declaration["after_sha256"])
            checkpoint = json.loads(CURRENT_MANIFEST.read_text())["source_commit"]
            self.assertEqual(after, subprocess.check_output(
                ["git", "show", checkpoint + ":" + declaration["output"]], cwd=ROOT))
            fields, old_rows, _ = csv_records(original)
            _, new_rows, _ = csv_records(after)
            identity = "dataset" if "dataset" in fields else "case_id"
            self.assertEqual(len(old_rows), len(new_rows))
            for (old, _, _), (new, _, _) in zip(old_rows, new_rows):
                if old[identity] != "GSE134151-10407":
                    self.assertEqual(old, new)
                    continue
                allowed = {"final_score", "trust_score", "q5_2_mechanistic_novelty",
                           "audit_decision", "audit_rationale"}
                self.assertEqual({k: v for k, v in old.items() if k not in allowed},
                                 {k: v for k, v in new.items() if k not in allowed})
                self.assertEqual(old["final_score"], "61.0")
                self.assertEqual(new["final_score"], "71.0")
                self.assertEqual(new["audit_decision"], "AUTHOR_PRIORITY_NOVELTY_REASSESSMENT")
                self.assertIn("Independent novelty and biological aging CANNOT VERIFY",
                              new["audit_rationale"])
                if "q5_2_mechanistic_novelty" in fields:
                    self.assertEqual(new["q5_2_mechanistic_novelty"], "15")
                    self.assertEqual(new["q5_1_literature_saturation"], "3")
                    components = ["q0_plausibility", "q1_experimental_model", "q2_1_aging_pathways",
                                  "q2_2_age_diseases", "q3_1_experimental_rigor",
                                  "q4_1_translational_potential", "q5_1_literature_saturation",
                                  "q5_2_mechanistic_novelty"]
                    self.assertEqual(sum(float(new[k]) for k in components), 71.0)

    def test_superseded_application_cannot_overwrite_priority_reassessment(self):
        for obsolete in (self.manifest, json.loads(NOVELTY_MANIFEST.read_text())):
            with self.assertRaises(CorrectionError):
                prepare(ROOT, obsolete)

    def test_component_review_changes_only_adjudicated_fields(self):
        current = json.loads(CURRENT_MANIFEST.read_text())
        for declaration in current["files"]:
            raw = subprocess.check_output(
                ["git", "show", current["source_commit"] + ":" + declaration["source"]], cwd=ROOT)
            after = patch_csv(raw, declaration)
            output = ROOT / declaration["output"]
            require_development_artifact(output)
            self.assertEqual(after, output.read_bytes())
            before_rows = csv_records(raw)[1]
            after_rows = csv_records(after)[1]
            self.assertEqual(len(before_rows), len(after_rows))
            changed = []
            for i, ((old, _, _), (new, _, _)) in enumerate(zip(before_rows, after_rows), 2):
                if old == new:
                    continue
                changed.append(i)
                allowed = {"final_score", "trust_score", "q2_1_aging_pathways",
                           "q3_1_experimental_rigor", "audit_decision", "audit_rationale"}
                self.assertEqual({k: v for k, v in old.items() if k not in allowed},
                                 {k: v for k, v in new.items() if k not in allowed})
                self.assertEqual((old["final_score"], new["final_score"]), ("71.0", "78.0"))
                self.assertIn("Independent novelty and biological aging CANNOT VERIFY",
                              new["audit_rationale"])
                if "q3_1_experimental_rigor" in new:
                    self.assertEqual(new["q3_1_experimental_rigor"], "10")
                    self.assertEqual(new["q2_1_aging_pathways"], "15")
                    self.assertEqual(new["q0_plausibility"], "0")
                    self.assertEqual(new["q5_2_mechanistic_novelty"], "15")
            self.assertEqual(changed, [op["record_number"] for op in declaration["operations"]])

    def test_new_score_assignments_reject_invalid_tiers_and_inconsistent_totals(self):
        current = json.loads(CURRENT_MANIFEST.read_text())
        declaration = current["files"][-1]
        raw = subprocess.check_output(
            ["git", "show", current["source_commit"] + ":" + declaration["source"]], cwd=ROOT)
        for changes in ({"q3_1_experimental_rigor": "8"},
                        {"q4_1_translational_potential": "3"},
                        {"final_score": "108.0"}, {"final_score": "NaN"},
                        {"q3_1_experimental_rigor": "Infinity"}):
            invalid = copy.deepcopy(declaration)
            invalid["operations"][0]["changes"].update(changes)
            with self.subTest(changes=changes), self.assertRaises(CorrectionError):
                patch_csv(raw, invalid)

    def test_partial_component_withholding_cannot_keep_numeric_total(self):
        current = json.loads(CURRENT_MANIFEST.read_text())
        declaration = current["files"][-1]
        raw = subprocess.check_output(
            ["git", "show", current["source_commit"] + ":" + declaration["source"]], cwd=ROOT)
        declaration["operations"][0]["changes"]["q3_1_experimental_rigor"] = ""
        with self.assertRaises(CorrectionError):
            patch_csv(raw, declaration)

    def test_linked_score_fields_cannot_disagree(self):
        current = json.loads(CURRENT_MANIFEST.read_text())
        declaration = current["files"][0]
        raw = subprocess.check_output(
            ["git", "show", current["source_commit"] + ":" + declaration["source"]], cwd=ROOT)
        declaration["operations"][0]["changes"]["trust_score"] = "108.0"
        with self.assertRaises(CorrectionError):
            patch_csv(raw, declaration)


if __name__ == "__main__":
    unittest.main()
