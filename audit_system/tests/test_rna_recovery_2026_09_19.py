"""Small CSV software fixtures; these are not empirical measurements."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from audit_system.rna_recovery_2026_09_19 import (
    compare_saved_prefix, file_record, r_median_differences, read_saved_table,
    validate_plan_scope,
)


class RnaRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.source = Path(self.temp.name) / 'software-fixture.csv'

    def write(self, text):
        self.source.write_text(text)
        return self.source

    def test_r_median_propagates_a_single_missing_cell_and_uses_even_midpoint(self):
        values = np.array([[.1, .2, np.nan], [.9, .4, .7]])
        differences, medians = r_median_differences(values)
        np.testing.assert_allclose(medians, [.5, .3, np.nan], rtol=0, atol=1e-15)
        np.testing.assert_allclose(differences, [[-.4, -.1, np.nan], [.4, .1, np.nan]], rtol=0, atol=1e-15)

    def test_odd_median_and_all_na_column(self):
        values = np.array([[.1, np.nan], [.8, np.nan], [.3, np.nan]])
        differences, medians = r_median_differences(values)
        np.testing.assert_allclose(medians, [.3, np.nan], rtol=0, atol=0)
        self.assertEqual(differences[2, 0], 0)

    def test_r_header_and_identity_are_preserved(self):
        result = read_saved_table(self.write('"","100","101"\n"GSM1",NA,0.2\n"GSM2",NA,0.4\n'), yugene=True)
        self.assertEqual(result['samples'], ['GSM1', 'GSM2'])
        self.assertEqual(result['genes'], ['100', '101'])
        np.testing.assert_allclose(result['values'], [[np.nan, .2], [np.nan, .4]], rtol=0, atol=0)

    def test_rejects_ragged_unterminated_duplicate_and_invalid_identities(self):
        fixtures = [
            'sample_id,g,g\nGSM1,0,0\n', 'sample_id,g,\nGSM1,0,0\n',
            'sample_id,g\nGSM1,0\nGSM1,0\n', 'sample_id,g\nGSM1,0',
            'sample_id,g\nGSM1,0,1\n', 'sample_id,g\nGSM１,0\n',
            'sample_id,g\n GSM1,0\n', 'sample_id,g\n',
        ]
        for fixture in fixtures:
            with self.subTest(fixture=fixture), self.assertRaises(ValueError):
                read_saved_table(self.write(fixture))

    def test_rejects_undocumented_missingness_nonfinite_and_domain_values(self):
        for token in ['', 'NaN', 'Na', 'null', 'Inf', '1e999', '-0.1', '1.1', ' 0.1', '0\x00']:
            with self.subTest(token=token), self.assertRaises(ValueError):
                read_saved_table(self.write('sample_id,g\nGSM1,' + token + '\n'), yugene=True)

    def test_expected_roster_and_feature_order_must_match(self):
        path = self.write('sample_id,g,h\nGSM1,0,1\nGSM2,0,1\n')
        with self.assertRaisesRegex(ValueError, 'Sample identity/order mismatch'):
            read_saved_table(path, expected_samples=['GSM2', 'GSM1'])
        with self.assertRaisesRegex(ValueError, 'Feature identity/order mismatch'):
            read_saved_table(path, expected_genes=['h', 'g'])

    def test_detects_source_mutation(self):
        path = self.write('sample_id,g\nGSM1,0\n')
        before = file_record(path)
        after = dict(before, sha256='0' * 64)
        with patch('audit_system.rna_recovery_2026_09_19.file_record', side_effect=[before, after]):
            with self.assertRaisesRegex(ValueError, 'Source changed'):
                read_saved_table(path)

    def test_all_terminated_partial_cells_checked_but_last_parseable_token_excluded(self):
        path = self.write('sample_id,g,h,i\nGSM1,NA,0.1,0.2\nGSM2,NA,0.3,987')
        matrix = np.array([[np.nan, .1, .2], [np.nan, .3, .4]])
        report = compare_saved_prefix(path, ['GSM1', 'GSM2'], ['g', 'h', 'i'], matrix, atol=1e-12)
        self.assertEqual(report['checked_cells'], 5)
        self.assertEqual(report['checked_na_cells'], 2)
        self.assertEqual(report['checked_numeric_cells'], 3)
        self.assertEqual(report['partial_terminated_cells'], 2)
        self.assertEqual(report['unverified_final_token'], '987')
        self.assertEqual(report['unverified_final_token_feature'], 'i')

    def test_prefix_fails_on_na_mask_or_intact_value_disagreement(self):
        matrix = np.array([[np.nan, .1], [np.nan, .3]])
        for fixture in ['sample_id,g,h\nGSM1,0,0.1\nGSM2,NA,0',
                        'sample_id,g,h\nGSM1,NA,0.2\nGSM2,NA,0',
                        'sample_id,g,h\nGSM2,NA,0.1\nGSM1,NA,0']:
            with self.subTest(fixture=fixture), self.assertRaises(ValueError):
                compare_saved_prefix(self.write(fixture), ['GSM1', 'GSM2'], ['g', 'h'], matrix, atol=1e-12)

    def test_prefix_tolerance_has_no_relative_component(self):
        matrix = np.array([[.1], [.2]])
        path = self.write('sample_id,g\nGSM1,0.1000000000005\nGSM2,0')
        self.assertEqual(compare_saved_prefix(path, ['GSM1', 'GSM2'], ['g'], matrix, atol=1e-12)['status'], 'PASS')
        path = self.write('sample_id,g\nGSM1,0.100000000002\nGSM2,0')
        with self.assertRaisesRegex(ValueError, 'differ beyond'):
            compare_saved_prefix(path, ['GSM1', 'GSM2'], ['g'], matrix, atol=1e-12)

    def test_arithmetic_mode_retains_negative_tokens_without_domain_approval(self):
        path = self.write('sample_id,g\nGSM1,-2.22044604925031e-16\nGSM2,0\n')
        with self.assertRaisesRegex(ValueError, 'outside'):
            read_saved_table(path, yugene=True)
        table = read_saved_table(path, yugene=False)
        self.assertEqual(table['values'][0, 0], -2.22044604925031e-16)
        differences, median = r_median_differences(table['values'])
        self.assertLess(median[0], 0)
        self.assertLess(differences[0, 0], 0)

    def test_arithmetic_mode_requires_preserved_plan_and_prohibits_repair_certification(self):
        original = {'schema': 'clockbase-rna-truncation-recovery-plan-v1', 'inputs': [],
                    'definition_sources': [], 'studies': ['GSE1'],
                    'reconstruction': {'comparison_atol': 1e-12, 'comparison_rtol': 0,
                                       'missingness': 'NA', 'definition': 'median', 'comparison_scope': 'prefix'}}
        self.source.write_text(json.dumps(original))
        plan = dict(original, schema='clockbase-rna-truncation-arithmetic-plan-v2',
                    arithmetic_only=True, source_value_policy='preserve',
                    repair_must_remain_unverified=True, v1_plan=file_record(self.source))
        self.assertFalse(validate_plan_scope(original))
        self.assertTrue(validate_plan_scope(plan))
        for key, value in [('arithmetic_only', False), ('source_value_policy', 'clamp'),
                           ('repair_must_remain_unverified', False), ('studies', ['GSE2'])]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                validate_plan_scope(dict(plan, **{key: value}))
        with self.assertRaisesRegex(ValueError, 'comparison_atol'):
            validate_plan_scope(dict(plan, reconstruction=dict(plan['reconstruction'], comparison_atol=1e-3)))


if __name__ == '__main__':
    unittest.main()
