import math
import tempfile
import unittest
from pathlib import Path

from audit_system.hashing import (
    ArtifactIdentityError,
    canonical_json_bytes,
    inspect_file,
    sha256_value,
)


class HashingTests(unittest.TestCase):
    def test_canonical_hash_ignores_dictionary_insertion_order(self):
        left = {"b": 2, "a": 1}
        right = {"a": 1, "b": 2}

        self.assertEqual(canonical_json_bytes(left), b'{"a":1,"b":2}')
        self.assertEqual(
            sha256_value(left),
            "43258cff783fe7036d8a43033f830adfc60ec037382473548ac742b888292777",
        )
        self.assertEqual(sha256_value(left), sha256_value(right))

    def test_canonical_json_rejects_nonfinite_numbers(self):
        for value in (math.nan, math.inf, -math.inf):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    canonical_json_bytes({"value": value})

    def test_file_identity_records_exact_bytes_and_changes_with_content(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "source.csv"
            path.write_bytes(b"sample,value\nS1,1\n")

            before = inspect_file(path)

            self.assertEqual(before.path, path.resolve())
            self.assertEqual(before.byte_count, 18)
            self.assertEqual(
                before.sha256,
                "f1c7a285edd9a29c9e4d865bf19ac778c617e11fc57c54b5f6d7c85b8254830b",
            )

            path.write_bytes(b"sample,value\nS1,2\n")
            after = inspect_file(path)

            self.assertNotEqual(before.sha256, after.sha256)

    def test_missing_file_reports_identity_failure(self):
        missing = Path("/definitely-missing-clockbase-audit-source")

        with self.assertRaisesRegex(
            ArtifactIdentityError,
            "definitely-missing-clockbase-audit-source",
        ):
            inspect_file(missing)


if __name__ == "__main__":
    unittest.main()
