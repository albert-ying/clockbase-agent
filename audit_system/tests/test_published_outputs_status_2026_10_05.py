"""Tests of published row matching. All fixtures are synthetic."""
import contextlib
import csv
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from audit_system import published_outputs_status_2026_09_29 as ps


class PublishedStatusTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR'))
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        self.repairs = self.root / 'repairs'
        self.repairs.mkdir()
        self.published = self.root / 'bars.csv'
        self.out = self.root / 'status'
        self.atlas_row = {'gse_id': 'GSE1', 'comparison_name': 'Ouabain vs vehicle', 'clock': 'RNA',
                          'effect_size': 5.737, 'pvalue': 8.59e-6, 'v3_final_class': 'VERIFIED',
                          'v3_actions': 'SWAP_ARMS_CONFIRMED_INVERSION;RECOMPUTE_DECLARED_ESTIMAND',
                          'v3_cohen_d': -5.737, 'v3_p': 8.59e-6, 'v3_fdr': .001,
                          'v3_n_case': 3, 'v3_n_control': 3}
        self.item = {'modality': 'human RNA-seq', 'gse_id': 'GSE1', 'condition_name': 'ouabain',
                     'effect': -5.737, 'p': 8.59e-6, 'direction': 'anti-aging'}
        published = {'Figure3_bars': (self.published, {'effect': 'effect', 'p': 'p', 'claim': 'direction'})}
        patch = mock.patch.dict(ps.PUBLISHED, published, clear=True)
        patch.start()
        self.addCleanup(patch.stop)

    @staticmethod
    def write_csv(path, rows):
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def run_status(self, rows=None, items=None):
        self.write_csv(self.repairs / 'atlas_repaired.csv', rows or [self.atlas_row])
        self.write_csv(self.published, items or [self.item])
        with contextlib.redirect_stdout(io.StringIO()), mock.patch.object(ps.paths, 'portable', return_value='inputs:fixture.csv'):
            ps.run({'hs_rna': self.repairs}, self.out)
        with (self.out / 'Figure3_bars_v3_status.csv').open() as handle:
            return list(csv.DictReader(handle))

    def test_sign_reversed_bar_keeps_row_identity_and_confirmed_swap_fate(self):
        row = self.run_status()[0]
        self.assertEqual(row['match'], 'SIGN_REVERSED')
        self.assertEqual(row['atlas_modality'], 'hs_rna')
        self.assertEqual(row['atlas_row_ordinal'], '0')
        self.assertEqual(row['comparison_name'], self.atlas_row['comparison_name'])
        self.assertEqual(row['v3_final_class'], 'VERIFIED')
        self.assertEqual(row['outcome'], 'CONFIRMED_CLAIMED_DIRECTION')
        self.assertEqual(row['direction_flipped_by_repair'], 'True')
        self.assertEqual(float(row['atlas_historical_effect']), 5.737)
        summary = json.loads((self.out / 'summary.json').read_text(), parse_constant=lambda value: self.fail(value))
        self.assertEqual(summary['Figure3_bars']['match'], {'SIGN_REVERSED': 1})
        self.assertEqual(summary['Figure3_bars']['direction_flipped_by_repair'], 1)
        self.assertEqual(summary['Figure3_bars']['outcomes'], {'CONFIRMED_CLAIMED_DIRECTION': 1})
        self.assertNotIn('/Entropy', (self.out / 'summary.json').read_text())

    def test_exact_match_has_priority_over_sign_reversed(self):
        exact = {**self.atlas_row, 'effect_size': -5.737, 'v3_cohen_d': 5.737, 'v3_actions': ''}
        row = self.run_status([self.atlas_row, exact])[0]
        self.assertEqual(row['match'], 'EXACT')
        self.assertEqual(row['atlas_row_ordinal'], '1')
        self.assertEqual(row['outcome'], 'CONTRADICTS_CLAIMED_DIRECTION')
        self.assertEqual(row['direction_flipped_by_repair'], 'False')

    def test_sign_reversed_requires_matching_p_and_study(self):
        items = [{**self.item, 'p': .02}, {**self.item, 'gse_id': 'GSE2'}]
        rows = self.run_status(items=items)
        self.assertEqual([r['match'] for r in rows], ['NO_EXACT_ROW', 'STUDY_NOT_IN_ATLAS'])

    def test_sign_reversed_candidates_with_different_fates_are_ambiguous(self):
        other = {**self.atlas_row, 'v3_cohen_d': 5.737}
        row = self.run_status([self.atlas_row, other])[0]
        self.assertEqual(row['match'], 'AMBIGUOUS')
        self.assertEqual(row['candidates'], '2')

    def test_duplicate_match_uses_representative_fate(self):
        duplicate = {**self.atlas_row, 'v3_final_class': 'DUPLICATE', 'v3_representative_class': 'PROVISIONAL'}
        row = self.run_status([duplicate])[0]
        self.assertEqual(row['match'], 'SIGN_REVERSED')
        self.assertEqual(row['v3_final_class'], 'DUPLICATE')
        self.assertEqual(row['outcome'], 'CONFIRMED_CLAIMED_DIRECTION')


if __name__ == '__main__':
    unittest.main()
