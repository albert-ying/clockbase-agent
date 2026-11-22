"""Adversarial CSV/identity fixtures for the supplementary roster boundary."""
import csv
from pathlib import Path
import tempfile
import unittest

from scripts.withdraw_supplementary_table_s8 import source_roster


class SupplementaryRosterTests(unittest.TestCase):
    def write_fixture(self, path, header, rows):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    def test_unquoted_comma_row_width_is_rejected_before_identity_use(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.csv"
            path.write_text("gse_id,comparison_name,covariate\nGSE168846,muDX400 versus control,cell line, response\n")
            with self.assertRaisesRegex(ValueError, "4 fields; expected 3"):
                source_roster(path, "Mouse RNA-seq")

    def test_blank_joint_and_whitespace_identifiers_remain_unassigned(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.csv"
            self.write_fixture(path, ["gse_id", "comparison_name"], [
                ["GSE1", "A"], ["GSE1/GSE2", "B"], ["", "C"], [" GSE3", "D"]])
            groups, unresolved, count = source_roster(path, "Mouse RNA-seq")
            self.assertEqual(groups, {("GSE1", "Mouse RNA-seq"): [2]})
            self.assertEqual([r["source_identifier"] for r in unresolved], ["GSE1/GSE2", "", " GSE3"])
            self.assertEqual(count, 4)
            self.assertTrue(all(r["status"] == "CANNOT VERIFY" for r in unresolved))

    def test_duplicate_headers_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.csv"
            self.write_fixture(path, ["gse_id", "gse_id"], [["GSE1", "GSE2"]])
            with self.assertRaisesRegex(ValueError, "Unique source headers"):
                source_roster(path, "Mouse RNA-seq")

    def test_repeated_source_rows_keep_every_record_witness(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fixture.csv"
            self.write_fixture(path, ["gse_id", "comparison_name"], [["GSE1", "A"], ["GSE1", "B"]])
            groups, unresolved, count = source_roster(path, "Human RNA-seq")
            self.assertEqual(groups, {("GSE1", "Human RNA-seq"): [2, 3]})
            self.assertEqual(unresolved, [])
            self.assertEqual(count, 2)


if __name__ == "__main__":
    unittest.main()
