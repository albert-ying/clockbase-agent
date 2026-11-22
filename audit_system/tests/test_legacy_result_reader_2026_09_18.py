"""CSV boundary tests; fabricated inputs below are software fixtures only."""
import ast
import contextlib
import csv
import io
import os
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from audit_system.tests import require_development_artifact


REPO = Path(__file__).resolve().parents[2]
IMPLEMENTATION = REPO / "s01_visualize_results.py"


def load_reader_functions():
    """Load only pure boundary functions; never execute legacy plot generation."""
    tree = ast.parse(IMPLEMENTATION.read_text())
    names = {"read_result_csv", "read_cached_result_csv", "load_data", "process_directory"}
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)
                 and node.name in names]
    if {node.name for node in functions} != names:
        raise AssertionError("Legacy boundary functions are missing")
    namespace = {"csv": csv, "io": io, "os": os, "Path": Path,
                 "re": re, "np": np, "pd": pd}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(IMPLEMENTATION), "exec"), namespace)
    return namespace


class LegacyResultReaderTests(unittest.TestCase):
    def setUp(self):
        self.functions = load_reader_functions()
        self.read = self.functions["read_result_csv"]
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def write_source(self, text, study="GSE1"):
        path = self.root / study / "aging_analysis.csv"
        path.parent.mkdir(exist_ok=True)
        path.write_bytes(text if isinstance(text, bytes) else text.encode())
        return path

    def test_preserved_gse168846_is_rejected_before_column_shift(self):
        source = REPO / "result_mm_rna/GSE168846/aging_analysis.csv"
        require_development_artifact(source)
        before = source.read_bytes()
        with self.assertRaisesRegex(ValueError, "CSV_RAGGED.*12 fields; expected 11"):
            self.read(source, expected_study_id="GSE168846")
        self.assertEqual(source.read_bytes(), before)

    def test_preserved_valid_source_keeps_pandas_values_and_types(self):
        source = REPO / "result_mm_rna/GSE121759/aging_analysis.csv"
        require_development_artifact(source)
        pd.testing.assert_frame_equal(
            self.read(source, expected_study_id="GSE121759"), pd.read_csv(source),
            check_exact=True,
        )

    def test_quoted_comma_and_bom_are_valid_csv(self):
        path = self.write_source('\ufeffgse_id,covariate,value\nGSE1,"age, sex",0.25\n')
        result = self.read(path, expected_study_id="GSE1")
        self.assertEqual(result.loc[0, "covariate"], "age, sex")
        pd.testing.assert_frame_equal(result, pd.read_csv(path), check_exact=True)

    def test_ambiguous_csv_and_missing_payloads_are_rejected(self):
        fixtures = [
            ("gse_id,value\n", "EMPTY_SOURCE"),
            ("gse_id,value\nGSE1,1,extra\n", "CSV_RAGGED"),
            ('gse_id,value\nGSE1,"unterminated\n', "CSV_READ_ERROR"),
            ("gse_id,value,value\nGSE1,1,2\n", "AMBIGUOUS_HEADER"),
            ("gse_id,\nGSE1,1\n", "AMBIGUOUS_HEADER"),
            ("other,value\nGSE1,1\n", "MISSING_COLUMNS"),
            ("version https://git-lfs.github.com/spec/v1\noid sha256:fixture\n", "LFS_POINTER"),
            (b"gse_id,value\nGSE1,\xff\n", "CSV_READ_ERROR"),
            (b"gse_id,value\nGSE1,1\x002\n", "CSV_READ_ERROR"),
        ]
        for content, reason in fixtures:
            with self.subTest(reason=reason, content=content):
                with self.assertRaisesRegex(ValueError, reason):
                    self.read(self.write_source(content), expected_study_id="GSE1")

    def test_identifiers_are_not_trimmed_split_or_coerced(self):
        for identity in ("", " GSE1", "GSE1 ", "GSE1/GSE2", "GSE１", "comparison label", "NA"):
            with self.subTest(identity=identity):
                path = self.write_source(f"gse_id,value\n{identity},1\n")
                with self.assertRaisesRegex(ValueError, "INVALID_GSE_ID"):
                    self.read(path)

    def test_parent_mismatch_and_unknown_parent_are_distinct(self):
        path = self.write_source("gse_id,value\nGSE2,1\n")
        with self.assertRaisesRegex(ValueError, "PARENT_ID_MISMATCH"):
            self.read(path, expected_study_id="GSE1")
        with self.assertRaisesRegex(ValueError, "PARENT_ID_UNRESOLVED"):
            self.read(path, expected_study_id="unresolved folder")

    def test_cache_accepts_multiple_ids_but_not_invalid_identity(self):
        path = self.write_source("gse_id,pvalue\nGSE1,0.2\nGSE2,0.3\n")
        pd.testing.assert_frame_equal(self.read(path), pd.read_csv(path), check_exact=True)
        path.write_text("gse_id,pvalue\nGSE1,0.2\ncomparison label,0.3\n")
        with self.assertRaisesRegex(ValueError, "INVALID_GSE_ID"):
            self.read(path)

    def test_cache_requires_columns_used_by_the_plotting_branch(self):
        # Small boundary fixture derived from a preserved cache; no scientific result.
        source = REPO / "result/combined_df_hs_rna.csv"
        require_development_artifact(source)
        frame = pd.read_csv(source, nrows=2)
        path = self.root / "cache.csv"
        frame.to_csv(path, index=False)
        reader = self.functions["read_cached_result_csv"]
        pd.testing.assert_frame_equal(reader(path), pd.read_csv(path), check_exact=True)
        for missing in ("pvalue", "log2FoldChange", "comparison_name", "trust_score"):
            with self.subTest(missing=missing):
                frame.drop(columns=[missing]).to_csv(path, index=False)
                with self.assertRaisesRegex(ValueError, f"MISSING_COLUMNS.*{missing}"):
                    reader(path)

    def test_directory_failure_prevents_partial_aggregation(self):
        # One preserved valid source plus one malformed software fixture.
        source = REPO / "result_mm_rna/GSE121759/aging_analysis.csv"
        require_development_artifact(source)
        self.write_source(source.read_bytes(), "GSE121759")
        bad = self.write_source("gse_id,value\nGSE2,1,extra\n", "GSE2")
        stream = io.StringIO()
        with contextlib.redirect_stdout(stream), self.assertRaisesRegex(ValueError, "1 rejected"):
            self.functions["process_directory"](self.root)
        self.assertIn("CSV_RAGGED", stream.getvalue())
        self.assertIn(str(bad), stream.getvalue())

    def test_directory_errors_and_rejection_order_are_explicit(self):
        self.write_source("gse_id,value\nGSE9,1,extra\n", "GSE9")
        self.write_source("gse_id,value\nGSE2,1,extra\n", "GSE2")
        with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(ValueError) as captured:
            self.functions["process_directory"](self.root)
        self.assertLess(str(captured.exception).index("GSE2"), str(captured.exception).index("GSE9"))
        with self.assertRaises(FileNotFoundError):
            self.functions["process_directory"](self.root / "absent")
        (self.root / "empty").mkdir()
        with self.assertRaisesRegex(ValueError, "No accepted"):
            self.functions["process_directory"](self.root / "empty")
        def unreadable_walk(base, *, onerror):
            onerror(PermissionError("fixture directory unreadable"))
            return []
        with patch.object(os, "walk", unreadable_walk), self.assertRaises(PermissionError):
            self.functions["process_directory"](self.root)


if __name__ == "__main__":
    unittest.main()
