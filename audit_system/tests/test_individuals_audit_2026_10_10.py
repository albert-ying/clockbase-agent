import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import pandas as pd

from audit_system import individuals_audit_2026_10_10 as ia

CASE, CONTROL = ['G0', 'G1', 'G2', 'G3'], ['G4', 'G5', 'G6', 'G7']


def answer(individuals, shared='YES', material='PRIMARY_OR_IN_VIVO'):
    return {'material': material, 'shared_individuals': shared, 'individual_evidence': 'titles', 'individuals': individuals}


REPLICATES = {'G0': 'm1', 'G1': 'm1', 'G2': 'm2', 'G3': 'm2', 'G4': 'm3', 'G5': 'm3', 'G6': 'm4', 'G7': 'm4'}
PAIRS = {'G0': 'a', 'G1': 'b', 'G2': 'c', 'G3': 'd', 'G4': 'a', 'G5': 'b', 'G6': 'c', 'G7': 'd'}
NONE = {g: None for g in CASE + CONTROL}


class Reconcile(unittest.TestCase):
    def test_parse_accepts_an_object_with_surrounding_text(self):
        self.assertEqual(ia.parse('{"a": 1}'), {'a': 1})
        self.assertEqual(ia.parse('Here it is:\n{"a": 1}\nDone.'), {'a': 1})
        self.assertIsNone(ia.parse('[1, 2]'))
        self.assertIsNone(ia.parse('no json'))

    def test_agreed_technical_replicates_save_one_partition_with_reader_independent_names(self):
        other = {g: v.replace('m', 'mouse_') for g, v in REPLICATES.items()}
        out = ia.reconcile(CASE, CONTROL, [answer(REPLICATES), answer(other)])
        self.assertEqual((out['status'], out['findings']), ('SHARED_AGREED', ['REPEATED_INDIVIDUAL']))
        self.assertEqual(out['individuals'], {'G0': 'I001', 'G1': 'I001', 'G2': 'I002', 'G3': 'I002',
                                              'G4': 'I003', 'G5': 'I003', 'G6': 'I004', 'G7': 'I004'})

    def test_agreed_pairs(self):
        out = ia.reconcile(CASE, CONTROL, [answer(PAIRS), answer(PAIRS)])
        self.assertEqual((out['status'], out['findings']), ('SHARED_AGREED', ['PAIRED_ARMS']))

    def test_no_shared_individuals_and_cell_lines_are_independent(self):
        distinct = {g: g for g in CASE + CONTROL}
        self.assertEqual(ia.reconcile(CASE, CONTROL, [answer(NONE, 'NO'), answer(distinct, 'NO')])['status'], 'INDEPENDENT')
        self.assertEqual(ia.reconcile(CASE, CONTROL, [answer(REPLICATES, material='CELL_LINE'), answer(NONE, 'NO')])['status'],
                         'INDEPENDENT')

    def test_different_partitions_change_nothing(self):
        self.assertEqual(ia.reconcile(CASE, CONTROL, [answer(REPLICATES), answer(PAIRS)])['status'], 'SHARED_DISAGREE')

    def test_a_reader_who_cannot_tell_leaves_the_comparison_undetermined(self):
        self.assertEqual(ia.reconcile(CASE, CONTROL, [answer(REPLICATES), answer(NONE, 'CANNOT_TELL')])['status'], 'UNDETERMINED')
        partial = {'G0': 'm1', 'G1': 'm1'}
        self.assertEqual(ia.reconcile(CASE, CONTROL, [answer(REPLICATES), answer(partial)])['status'], 'UNDETERMINED')

    def test_a_missing_answer_fails(self):
        self.assertEqual(ia.reconcile(CASE, CONTROL, [answer(REPLICATES), None])['status'], 'MODEL_CALL_FAILED')
        self.assertEqual(ia.reconcile(CASE, CONTROL, [answer(REPLICATES), {'material': 'UNCLEAR'}])['status'], 'MODEL_CALL_FAILED')

    def test_individuals_outside_the_arms_are_ignored(self):
        extra = {**REPLICATES, 'G99': 'm1'}
        self.assertEqual(ia.reconcile(CASE, CONTROL, [answer(extra), answer(REPLICATES)])['status'], 'SHARED_AGREED')


class Items(unittest.TestCase):
    def build(self, raw, case=CASE, control=CONTROL):
        frame = pd.DataFrame([{'modality': 'mouse RNA-seq', 'gse_id': 'GSE1', 'comparison_name': 'KO vs WT',
                               'condition_name': 'KO', 'condition_category': 'Genetic', 'atlas_row_ordinal': 7}])
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp) / 'GEO_mm_rna' / 'processed_metadata_preharmonize'
            folder.mkdir(parents=True)
            if raw is not None:
                raw.rename_axis('GSM').to_csv(folder / 'GSE1.csv')
            with mock.patch.object(ia, 'GEO', Path(tmp)), \
                    mock.patch.object(ia.design_audit, 'arms_from_repairs', return_value={('mm_rna', 7): (case, control)}):
                return ia.build_items(frame, {})[0]

    def test_prompt_shows_geo_fields_of_both_arms_and_no_result(self):
        raw = pd.DataFrame({'title': [f'KO_rep{i}' for i in range(4)] + [f'WT_rep{i}' for i in range(4)], 'tissue': ['liver'] * 8},
                           index=CASE + CONTROL)
        item = self.build(raw)
        self.assertIn('G0: title=KO_rep0', item['prompt'])
        self.assertIn('tissue=liver', item['prompt'])
        self.assertEqual(len(item['prompt_sha256']), 64)
        for word in ('P value =', 'cohen', 'mean_diff'):
            self.assertNotIn(word, item['prompt'])

    def test_series_with_a_geo_individual_field_missing_table_or_large_arms_are_not_read(self):
        raw = pd.DataFrame({'title': list('abcdefgh'), 'mouse id': ['1', '1', '2', '2', '3', '3', '4', '4']}, index=CASE + CONTROL)
        self.assertEqual(self.build(raw)['skip'], 'GEO_INDIVIDUAL_FIELD')
        self.assertEqual(self.build(None)['skip'], 'METADATA_UNAVAILABLE')
        many = [f'G{i}' for i in range(130)]
        big = pd.DataFrame({'title': [f't{i}' for i in range(130)]}, index=many)
        self.assertEqual(self.build(big, many[:65], many[65:])['skip'], 'NOT_READ_ARMS_ABOVE_LIMIT')


if __name__ == '__main__':
    unittest.main()
