"""Integration regression checks using the archived empirical export, not fixtures."""
import json
from pathlib import Path
import runpy
import unittest

ROOT = Path(__file__).resolve().parents[1]


class ReconstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = runpy.run_path(str(ROOT / "scripts/reconstruct_down_syndrome_models_2026-09-07.py"))["run"]()

    def test_primary_results_remain_unchanged(self):
        previous = json.loads((ROOT / "docs/2026-09-07-down-syndrome-model-reconstruction.json").read_text())
        for old, new in zip(previous["GSE52588_results"], self.result["GSE52588_results"]):
            self.assertEqual(old["primary"], new["primary"])

    def test_diagnostics_preserve_family_identity_and_robust_inference(self):
        for row in self.result["GSE52588_results"]:
            self.assertIn("diagnostics", row)
            diagnostic = row["diagnostics"]
            self.assertEqual(len(diagnostic["leave_one_family_out"]), 29)
            self.assertEqual({x["family"] for x in diagnostic["leave_one_family_out"]},
                             {x["family"] for x in self.result["family_pairs"]})
            self.assertAlmostEqual(diagnostic["hc3_secondary"]["coefficient"], row["primary"]["coefficient"])
            if row["clock"] == "GrimAgeV2":
                self.assertGreater(diagnostic["hc3_secondary"]["pvalue"], 0.05)

    def test_replay_comparisons_are_persisted(self):
        self.assertIn("source_pvalue_comparisons", self.result)
        self.assertEqual(len(self.result["source_pvalue_comparisons"]), 10)
        self.assertTrue(all(r["matches"] for r in self.result["source_pvalue_comparisons"]))


if __name__ == "__main__":
    unittest.main()
