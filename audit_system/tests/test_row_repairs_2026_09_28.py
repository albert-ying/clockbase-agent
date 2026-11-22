"""Unit tests for row repairs. All fixtures are synthetic."""
import unittest

import numpy as np
import pandas as pd

from audit_system import row_repairs_2026_09_28 as rr

VALUES = pd.Series([1.0, 1.2, 0.9, 1.1, 3.0, 3.2, 2.9, 3.1], index=[f'GSM{i}' for i in range(8)])
ROW = {'gse_id': 'GSE1', 'comparison_name': 'Drug vs Vehicle', 'condition_name': 'Drug', 'clock': 'c',
       'statistical_test': 't-test', 'pvalue': 0.01, 'effect_size': -2.0, 'direction': -1, 'fdr': 0.02}
MEMBER = {'status': 'UNIQUE', 'source': 'replay', 'case_gsms': ['GSM0', 'GSM1', 'GSM2', 'GSM3'],
          'control_gsms': ['GSM4', 'GSM5', 'GSM6', 'GSM7']}


def verdict(*reasons, verdict='WITHHELD'):
    return {'row_ordinal': 0, 'verdict': verdict, 'reasons': list(reasons)}


class RowRepairTests(unittest.TestCase):
    def test_single_unit_arms_use_minimum_reason_before_constancy(self):
        # Synthetic software fixture; one sample per arm.
        member = {**MEMBER, **{key: [MEMBER[key][0]]
                              for key in ('case_gsms', 'control_gsms')}}
        record = rr.repair_row(ROW, member, verdict(), None, VALUES,
                               {'path': 'synthetic'}, {})
        self.assertEqual(record['final_class'], 'DESCRIPTIVE')
        self.assertEqual(record['recomputed']['n_case'], len(member['case_gsms']))
        self.assertEqual(record['recomputed']['n_control'], len(member['control_gsms']))
        self.assertNotIn('undefined_reason', record['recomputed'])
        self.assertNotIn(rr.design_estimand.CONSTANT_OUTPUT_CODE, record['reasons'])
        self.assertIn('DESCRIPTIVE_ONLY_ARM_BELOW_3', [a['action'] for a in record['actions']])
        self.assertEqual(rr.assign_fdr([record]), [])

    def test_constant_outputs_are_descriptive_and_leave_bh(self):
        # Synthetic fixtures cover equal and unequal arm constants and a
        # near-constant output at the declared round-off threshold.
        for a, b in (([0.0510765934346411] * 4, [0.0510765934346411] * 4),
                     ([2.] * 4, [5.] * 4),
                     ([1., 1. + 1e-13, 1., 1.], [2.] * 4)):
            with self.subTest(a=a, b=b):
                values = pd.Series(a + b, index=VALUES.index)
                record = rr.repair_row(ROW, MEMBER, verdict(), None, values, {'path': 'synthetic'}, {},
                                       series={'ages': {g: 20 + i for i, g in enumerate(VALUES.index)}})
                self.assertEqual(record['final_class'], 'DESCRIPTIVE')
                self.assertEqual(record['pending_reason'], ['outputs constant within both arms'])
                self.assertTrue(np.isfinite(record['recomputed']['mean_diff']))
                for key in ('cohen_d', 't', 'p', 'fdr_family'):
                    self.assertNotIn(key, record['recomputed'])
                self.assertEqual(rr.assign_fdr([record]), [])
                self.assertNotIn('design', record)

    def test_variation_in_one_arm_retains_defined_test(self):
        values = VALUES.copy()
        values.iloc[:4] = 2.
        record = rr.repair_row(ROW, MEMBER, verdict(), None, values, {'path': 'synthetic'}, {})
        self.assertEqual(record['final_class'], 'VERIFIED')
        self.assertIsNotNone(record['recomputed']['p'])

    def repair(self, reasons, adjudicated=None, member=MEMBER, text=None):
        return rr.repair_row(ROW, member, verdict(*reasons), adjudicated, VALUES, {'path': 'x'}, text or {})

    def test_confirmed_inversion_swaps_arms_and_recomputes(self):
        record = self.repair(['MEMBERSHIP_RECOVERED', 'LABEL_ORIENTATION_INVERTED'], 'LABEL_INVERTED_CONFIRMED')
        self.assertEqual(record['final_class'], 'VERIFIED')
        self.assertGreater(record['recomputed']['mean_diff'], 0)
        self.assertIn('SWAP_ARMS_CONFIRMED_INVERSION', [a['action'] for a in record['actions']])

    def test_unadjudicated_inversion_waits(self):
        record = self.repair(['LABEL_ORIENTATION_INVERTED'])
        self.assertEqual(record['final_class'], 'PENDING_ADJUDICATION')

    def test_assay_samples_removed_then_small_arm_becomes_descriptive(self):
        text = {'GSM0': 'MeDIP DNA', 'GSM1': 'MeDIP DNA'}
        record = self.repair(['NON_STANDARD_ASSAY'], text=text)
        self.assertEqual(record['actions'][0]['removed_gsms'], ['GSM0', 'GSM1'])
        self.assertEqual(record['final_class'], 'DESCRIPTIVE')
        self.assertIsNone(record['recomputed'].get('p'))
        self.assertEqual(record['recomputed']['n_case'], 2)

    def test_soft_flags_keep_row_as_provisional(self):
        record = self.repair(['MEMBERSHIP_RECOVERED', 'TISSUE_OUT_OF_DOMAIN'])
        self.assertEqual((record['final_class'], record['soft_flags']), ('PROVISIONAL', ['TISSUE_OUT_OF_DOMAIN']))
        supported = self.repair(['LABEL_UNCONFIRMED'], 'LABEL_SUPPORTED')
        self.assertEqual(supported['final_class'], 'VERIFIED')

    def test_unconfirmed_adjudication_of_a_label_match_leaves_a_soft_flag(self):
        record = self.repair(['MEMBERSHIP_RECOVERED'], 'LABEL_UNCONFIRMED')
        self.assertEqual((record['final_class'], record['soft_flags']), ('PROVISIONAL', ['LABEL_UNCONFIRMED']))

    def test_unrecovered_and_input_rows_are_routed_not_dropped(self):
        self.assertEqual(self.repair(['MEMBERSHIP_UNRECOVERED'], member={'status': 'NONE'})['final_class'], 'PENDING_DESIGN')
        self.assertEqual(self.repair(['INPUT_EMPTY'])['final_class'], 'PENDING_INPUT')

    def test_bh_matches_reference(self):
        p = [0.01, 0.04, 0.03, 0.2]
        expected = [0.04, 0.0533333, 0.0533333, 0.2]
        np.testing.assert_allclose(rr.bh(p), expected, rtol=1e-5)


    def test_adjudicated_inversion_swaps_even_without_heuristic_flag(self):
        record = self.repair(['MEMBERSHIP_RECOVERED', 'LABEL_UNCONFIRMED'], 'LABEL_INVERTED_CONFIRMED')
        self.assertEqual(record['final_class'], 'VERIFIED')
        self.assertIn('SWAP_ARMS_CONFIRMED_INVERSION', [a['action'] for a in record['actions']])

    def test_adjudicated_design_issues_become_soft_flags(self):
        record = rr.repair_row(ROW, MEMBER, verdict('LABEL_UNCONFIRMED'), 'LABEL_SUPPORTED', VALUES, {'path': 'x'}, {},
                               ['TECHNICAL_REPLICATES', 'CASE_NOT_CONDITION'])
        self.assertEqual((record['final_class'], record['soft_flags']), ('PROVISIONAL', ['ADJUDICATED_TECHNICAL_REPLICATES']))

    def test_detached_samples_dropped_or_routed_and_used_arms_recorded(self):
        evidence = {'detached_gsms': ['GSM0', 'GSM4']}
        record = rr.repair_row(ROW, MEMBER, {**verdict('INPUT_PROFILE_DETACHED'), 'evidence': evidence}, None, VALUES,
                               {'path': 'x'}, {})
        self.assertEqual(record['actions'][0], {'action': 'DROP_DETACHED_INPUT_SAMPLES', 'removed_gsms': ['GSM0', 'GSM4']})
        self.assertEqual((record['final_class'], record['recomputed']['n_case']), ('VERIFIED', 3))
        self.assertEqual(record['inputs']['case_gsms'], ['GSM1', 'GSM2', 'GSM3'])
        self.assertEqual(record['inputs']['stored_case_gsms'], MEMBER['case_gsms'])
        emptied = rr.repair_row(ROW, MEMBER, {**verdict('INPUT_PROFILE_DETACHED'),
                                              'evidence': {'detached_gsms': MEMBER['control_gsms']}}, None, VALUES, {'path': 'x'}, {})
        self.assertEqual((emptied['final_class'], emptied['pending_reason']), ('PENDING_INPUT', ['INPUT_PROFILE_DETACHED']))
        repaired = rr.repair_row(ROW, MEMBER, {**verdict('INPUT_PROFILE_DETACHED'),
                                               'evidence': {'detached_gsms': MEMBER['control_gsms']}}, None, VALUES,
                                 {'path': 'x'}, {}, (), (VALUES, {'path': 'repaired'}, {'verdict': 'VERIFIED', 'reasons': []}))
        self.assertEqual(repaired['final_class'], 'VERIFIED')
        self.assertEqual(repaired['actions'][0]['action'], 'USE_REPAIRED_INPUT')

    def test_swapped_rows_record_swapped_arms(self):
        record = self.repair(['LABEL_ORIENTATION_INVERTED'], 'LABEL_INVERTED_CONFIRMED')
        self.assertEqual(record['inputs']['case_gsms'], MEMBER['control_gsms'])
        self.assertEqual(record['inputs']['stored_case_gsms'], MEMBER['case_gsms'])

    def test_repeated_tests_point_to_one_representative(self):
        def rec(ordinal, cls, case, control, clock='c'):
            return {'row_ordinal': ordinal, 'gse_id': 'GSE1', 'clock': clock, 'final_class': cls, 'actions': [],
                    'inputs': {'case_gsms': case, 'control_gsms': control}, 'recomputed': {'p': 0.01}}
        records = [rec(0, 'PROVISIONAL', ['a', 'b'], ['c', 'd']), rec(1, 'VERIFIED', ['b', 'a'], ['d', 'c']),
                   rec(2, 'VERIFIED', ['c', 'd'], ['a', 'b']), rec(3, 'VERIFIED', ['a', 'b'], ['c', 'd'], clock='other'),
                   rec(4, 'PENDING_INPUT', ['a', 'b'], ['c', 'd'])]
        self.assertEqual(rr.mark_duplicates(records), 1)
        self.assertEqual([r['final_class'] for r in records], ['DUPLICATE', 'VERIFIED', 'VERIFIED', 'VERIFIED', 'PENDING_INPUT'])
        self.assertEqual((records[0]['duplicate_of'], records[0]['representative_class']), (1, 'VERIFIED'))

    def test_inappropriate_control_restates_a_confirmed_inversion(self):
        swapped = rr.repair_row(ROW, MEMBER, verdict('LABEL_ORIENTATION_INVERTED'), 'LABEL_INVERTED_CONFIRMED', VALUES,
                                {'path': 'x'}, {}, ['CASE_NOT_CONDITION', 'CONTROL_INAPPROPRIATE'])
        self.assertEqual((swapped['final_class'], swapped['soft_flags']), ('VERIFIED', []))
        kept = rr.repair_row(ROW, MEMBER, verdict('LABEL_UNCONFIRMED'), 'LABEL_SUPPORTED', VALUES, {'path': 'x'}, {},
                             ['CONTROL_INAPPROPRIATE'])
        self.assertEqual(kept['soft_flags'], ['ADJUDICATED_CONTROL_INAPPROPRIATE'])


    def test_design_rows_are_label_supported_and_repaired_inputs_used(self):
        design = {**MEMBER, 'source': 'design_reconstruction', 'design_rules': {'m': {}}}
        record = rr.repair_row(ROW, design, verdict('MEMBERSHIP_UNRECOVERED', 'INPUT_EMPTY'), None, None, None, {}, (),
                               (VALUES, {'path': 'repaired'}, {'verdict': 'PROVISIONAL', 'reasons': []}))
        actions = [a['action'] for a in record['actions']]
        self.assertEqual(actions[:2], ['DESIGN_RECONSTRUCTED_DUAL_MODEL', 'USE_REPAIRED_INPUT'])
        self.assertEqual((record['final_class'], record['soft_flags']), ('PROVISIONAL', ['REPAIRED_CLOCK_PROVISIONAL']))
        withheld = rr.repair_row(ROW, design, verdict('INPUT_EMPTY'), None, None, None, {}, (),
                                 (VALUES, {'path': 'repaired'}, {'verdict': 'WITHHELD', 'reasons': ['MISSING_REQUIRED_SITES']}))
        self.assertEqual(withheld['final_class'], 'PENDING_INPUT')


    def test_repaired_input_replaces_historical_coverage_reasons(self):
        design = {**MEMBER, 'source': 'design_reconstruction'}
        record = rr.repair_row(ROW, design, verdict('INPUT_EMPTY', 'CLOCK_COVERAGE_LOW'), None, None, None, {}, (),
                               (VALUES, {'path': 'repaired'}, {'verdict': 'VERIFIED', 'reasons': []}))
        self.assertEqual(record['final_class'], 'VERIFIED')


def series(subjects=None, ages=None):
    return {'subject_field': 'individual' if subjects else None, 'subjects': subjects or {},
            'subject_source': 'GEO sample field' if subjects else None, 'age_field': 'harmonized.age' if ages else None,
            'ages': ages or {}, 'files': []}


class DesignEstimandTests(unittest.TestCase):
    def repair(self, values=VALUES, member=MEMBER, row=ROW, modality='mm_rna', contrast_types=(), **design):
        return rr.repair_row(row, member, verdict('MEMBERSHIP_RECOVERED'), None, values, {'path': 'x'}, {},
                             series=series(**design), modality=modality, contrast_types=contrast_types)

    def test_independent_samples_keep_the_declared_estimand(self):
        record = self.repair()
        plain = rr.repair_row(ROW, MEMBER, verdict('MEMBERSHIP_RECOVERED'), None, VALUES, {'path': 'x'}, {})
        self.assertEqual(record['recomputed'], plain['recomputed'])
        self.assertEqual(record['estimand'], rr.ESTIMAND)
        self.assertNotIn('design', record)
        self.assertEqual([a['action'] for a in record['actions']], ['RECOMPUTE_DECLARED_ESTIMAND'])

    def test_paired_arms_replace_the_welch_result_and_keep_it_beside(self):
        subjects = {f'GSM{i}': f's{i % 4}' for i in range(8)}
        values = VALUES + pd.Series([0.0, 0.05, -0.02, 0.01, 0.3, -0.1, 0.2, 0.0], index=VALUES.index)
        record = self.repair(values=values, subjects=subjects)
        self.assertEqual(record['design']['method'], 'PAIRED')
        self.assertEqual([a['action'] for a in record['actions']], ['RECOMPUTE_DECLARED_ESTIMAND', 'RECOMPUTE_DESIGN_ESTIMAND'])
        self.assertEqual((record['recomputed']['n_case'], record['recomputed_declared']['n_case']), (4, 4))
        self.assertNotEqual(record['recomputed']['p'], record['recomputed_declared']['p'])
        self.assertIn('Paired t-test', record['estimand'])
        self.assertEqual(record['final_class'], 'VERIFIED')

    def test_technical_replicates_of_too_few_individuals_become_descriptive(self):
        subjects = {f'GSM{i}': f's{i // 2}' for i in range(8)}  # two individuals per arm, two samples each
        record = self.repair(subjects=subjects)
        self.assertEqual(record['design']['method'], 'UNIT_MEANS')
        self.assertEqual((record['recomputed']['n_case'], record['recomputed']['n_control']), (2, 2))
        self.assertEqual(record['final_class'], 'DESCRIPTIVE')
        self.assertEqual(record['actions'][-1]['action'], 'DESCRIPTIVE_ONLY_UNITS_BELOW_3')
        self.assertIsNotNone(record['recomputed_declared']['p'])

    def test_recorded_ages_enter_for_human_age_clocks_only(self):
        ages = {f'GSM{i}': 40.0 + 3 * i for i in range(4)} | {f'GSM{i}': 42.0 + 3 * (i - 4) for i in range(4, 8)}
        row = {**ROW, 'clock': 'Horvathv1'}
        record = self.repair(row=row, modality='dnam', ages=ages)
        self.assertEqual(record['design']['method'], 'AGE_ADJUSTED')
        self.assertIn('Least squares', record['estimand'])
        self.assertNotIn('design', self.repair(row=row, modality='mm_rna', ages=ages))
        self.assertNotIn('design', self.repair(row={**ROW, 'clock': 'DunedinPACE'}, modality='dnam', ages=ages))
        self.assertNotIn('design', self.repair(row=row, modality='dnam', ages=ages, contrast_types=['AGE']))

    def test_age_separating_the_arms_is_a_soft_flag(self):
        ages = {f'GSM{i}': 20.0 + i for i in range(4)} | {f'GSM{i}': 60.0 + i for i in range(4, 8)}
        record = self.repair(row={**ROW, 'clock': 'Horvathv1'}, modality='dnam', ages=ages)
        self.assertEqual((record['final_class'], record['soft_flags']), ('PROVISIONAL', ['AGE_SEPARATES_ARMS']))
        self.assertNotIn('recomputed_declared', record)
        self.assertEqual(record['estimand'], rr.ESTIMAND)

    def test_overflow_clock_values_are_dropped_from_the_arms(self):
        values = VALUES.copy()
        values['GSM0'] = 1.05e39
        record = self.repair(values=values)
        self.assertEqual(record['actions'][0], {'action': 'DROP_OVERFLOW_CLOCK_VALUES', 'removed_gsms': ['GSM0']})
        self.assertEqual(record['inputs']['case_gsms'], ['GSM1', 'GSM2', 'GSM3'])
        self.assertLess(abs(record['recomputed']['mean_case']), 10)
        kept = rr.repair_row(ROW, MEMBER, verdict('MEMBERSHIP_RECOVERED'), None, values, {'path': 'x'}, {})
        self.assertGreater(kept['recomputed']['mean_case'], 1e30)  # without the design step nothing is dropped

    def test_an_arm_of_overflow_values_waits_for_input(self):
        values = VALUES.copy()
        values[['GSM0', 'GSM1', 'GSM2', 'GSM3']] = 1.0e39
        record = self.repair(values=values)
        self.assertEqual((record['final_class'], record['pending_reason']), ('PENDING_INPUT', ['CLOCK_VALUES_OVERFLOW']))


if __name__ == '__main__':
    unittest.main()
