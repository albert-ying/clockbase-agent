"""Empirical regression checks for comparison-specific score reconciliation."""
import importlib.util
import json
from pathlib import Path
import unittest
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("merge", ROOT / "s03_merge_scores_human.py")
merge = importlib.util.module_from_spec(spec)
spec.loader.exec_module(merge)


class IdentityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.master = pd.read_csv(ROOT / "result/main_results_hs_dnam_outlier_marked.csv")
        cls.sample = json.loads((ROOT / "minipaper_fork/output_hs_dnam/GSE140344-40/sample.json").read_text())

    def resolve(self, master, sample=None):
        self.assertTrue(hasattr(merge, "resolve_identity"), "Exact comparison resolver is missing")
        return merge.resolve_identity(master, sample or self.sample, "GrimAgeV2")

    def test_saved_iqplus_identity_selects_grim_not_lowest_p_or_all_t21(self):
        result = self.resolve(self.master)
        self.assertEqual(result["status"], "unique_identity_match")
        row = self.master.loc[result["master_index"]]
        self.assertEqual(row.comparison_name, "Trisomy21_IQplus_vs_Control")
        self.assertEqual(row.clock, "GrimAgeV2")
        self.assertAlmostEqual(row.pvalue, 0.003339363409602, places=14)

    def test_wrong_comparison_does_not_fall_back_to_same_study(self):
        master = self.master[self.master.comparison_name != self.sample["comparison_name"]]
        self.assertEqual(self.resolve(master)["status"], "no_identity_match")

    def test_duplicate_identical_rows_are_rejected(self):
        duplicated = pd.concat([self.master, self.master], ignore_index=True)
        self.assertEqual(self.resolve(duplicated)["status"], "ambiguous_identity")

    def test_nonfinite_statistic_is_rejected(self):
        sample = dict(self.sample, pvalue=float("nan"))
        self.assertEqual(self.resolve(self.master, sample)["status"], "invalid_scoring_identity")


if __name__ == "__main__":
    unittest.main()
