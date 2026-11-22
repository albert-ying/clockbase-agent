import unittest

import pandas as pd

from benchmark.published_clock_2026_10_05.groupbench import design_checks as dc


def frame(**columns):
    n = len(next(iter(columns.values())))
    return pd.DataFrame(columns, index=[f'G{i}' for i in range(n)])


class SubjectTable(unittest.TestCase):
    def test_named_geo_field_with_repeats_is_used(self):
        raw = frame(title=list('abcdef'), individual=['A02', 'A05', 'A12', 'A02', 'A05', 'A12'])
        table = dc.subject_table(raw, None)
        self.assertEqual((table['field'], table['source']), ('individual', 'geo'))
        self.assertEqual(table['subjects']['G3'], 'A02')

    def test_descriptive_field_is_not_a_subject(self):
        raw = frame(donor=['healthy male', 'MDD diagnosed male', 'mix', 'healthy male', 'mix', 'mix'])
        self.assertIsNone(dc.subject_table(raw, None))

    def test_status_field_name_is_not_a_subject(self):
        raw = frame(**{'subject status': ['1', '2', '1', '2']})
        self.assertIsNone(dc.subject_table(raw, None))

    def test_unique_ids_do_not_qualify(self):
        raw = frame(**{'patient id': ['P1', 'P2', 'P3', 'P4']})
        self.assertIsNone(dc.subject_table(raw, None))

    def test_harmonized_field_is_used_only_in_variant_v2(self):
        raw = frame(title=list('abcd'))
        processed = frame(**{'harmonized.subject': ['1', '2', '1', '2']})
        self.assertIsNone(dc.subject_table(raw, processed))
        self.assertEqual(dc.subject_table(raw, processed, harmonized=True)['source'], 'harmonized')

    def test_finest_named_field_wins(self):
        raw = frame(donor=['1', '1', '1', '2', '2', '2'], **{'patient id': ['a1', 'a1', 'a2', 'b1', 'b1', 'b2']})
        self.assertEqual(dc.subject_table(raw, None)['field'], 'patient id')


class SubjectFindings(unittest.TestCase):
    def test_repeated_subject_within_an_arm(self):
        subjects = {'G0': 's1', 'G1': 's1', 'G2': 's2', 'G3': 's3', 'G4': 's4'}
        out = dc.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4'], subjects)
        self.assertEqual(out['findings'], ['REPEATED_SUBJECT'])

    def test_paired_arms(self):
        subjects = {'G0': 's1', 'G1': 's2', 'G2': 's3', 'G3': 's1', 'G4': 's2', 'G5': 's3'}
        out = dc.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], subjects)
        self.assertEqual(out['findings'], ['PAIRED_AS_INDEPENDENT'])

    def test_one_subject_is_not_flagged(self):
        subjects = {f'G{i}': 's1' for i in range(6)}
        self.assertEqual(dc.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], subjects)['findings'], [])

    def test_independent_arms(self):
        subjects = {f'G{i}': f's{i}' for i in range(6)}
        self.assertEqual(dc.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], subjects)['findings'], [])

    def test_unknown_subjects_are_not_assessable(self):
        self.assertFalse(dc.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], {'G0': 's1'})['assessable'])


class Age(unittest.TestCase):
    def test_human_median_gap(self):
        ages = {'G0': 30, 'G1': 32, 'G2': 34, 'G3': 60, 'G4': 62, 'G5': 64}
        self.assertEqual(dc.age_finding(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], ages, True)['criterion'], 'median_gap')

    def test_balanced_arms(self):
        ages = {'G0': 30, 'G1': 40, 'G2': 50, 'G3': 31, 'G4': 40, 'G5': 49}
        self.assertIsNone(dc.age_finding(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], ages, True)['criterion'])

    def test_mouse_relative_gap(self):
        ages = {'G0': 8, 'G1': 8, 'G2': 12, 'G3': 12}
        self.assertEqual(dc.age_finding(['G0', 'G1'], ['G2', 'G3'], ages, False)['criterion'], 'median_gap')
        self.assertIsNone(dc.age_finding(['G0', 'G2'], ['G1', 'G3'], ages, False)['criterion'])

    def test_missing_age_is_not_assessed(self):
        self.assertIsNone(dc.age_finding(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], {'G0': 30, 'G3': 40}, True))

    def test_age_declared_by_name_or_classification(self):
        self.assertEqual(dc.age_declared({'comparison_name': 'Old vs young liver'}), 'name')
        self.assertIsNone(dc.age_declared({'comparison_name': 'Fast vs Late rebound post-romidepsin'}))
        self.assertIsNone(dc.age_declared({'comparison_name': 'Disease vs age-matched controls'}))
        self.assertIsNone(dc.age_declared({'comparison_name': 'p53 knockout vs wild type'}))
        row = {'comparison_name': 'Week 6 vs baseline', 'ranking_evidence': {'classification': {'answers': [
            {'contrast_type': 'INTERVENTION'}, {'contrast_type': 'DEVELOPMENT_OR_TIME'}]}}}
        self.assertEqual(dc.age_declared(row), 'classification')

    def test_harmonized_age_comes_first(self):
        processed = frame(**{'harmonized.age': ['30', '40', None]})
        raw = frame(age=['31', '41', '51'])
        self.assertEqual(dc.age_table(raw, processed, True), ({'G0': 30.0, 'G1': 40.0}, 'harmonized.age'))

    def test_geo_age_at_sampling_is_the_human_fallback(self):
        raw = frame(**{'age_sampling': ['35', '37 years', '48', '31'], 'age at diagnosis': ['20', '21', '22', '23']})
        processed = frame(**{'harmonized.age': [None, None, None, None]})
        ages, field = dc.age_table(raw, processed, True)
        self.assertEqual((field, ages['G1']), ('age_sampling', 37.0))
        self.assertEqual(dc.age_table(raw, processed, False), ({}, None))

    def test_other_ages_are_not_used(self):
        raw = frame(**{'age at diagnosis': ['20', '21', '22', '23'], 'gestational age': ['38', '39', '40', '41'],
                       'age': ['young adult', 'elderly adult', 'young adult', 'elderly adult']})
        self.assertEqual(dc.age_table(raw, None, True), ({}, None))


class Formula(unittest.TestCase):
    def test_age_term_on_the_right_hand_side(self):
        self.assertTrue(dc.formula_has_age('YingDamAge ~ condition + age + gender'))
        self.assertTrue(dc.formula_has_age('Hazard2_All_All_All_scaleddiff_EN_tAge ~ C(condition) + age_numeric'))
        self.assertTrue(dc.formula_has_age('clock ~ group + patient_age_numeric'))

    def test_clock_names_and_look_alike_terms_are_not_age(self):
        self.assertFalse(dc.formula_has_age('PhenoAge ~ condition + gender'))
        self.assertFalse(dc.formula_has_age('Hazard_Mouse_All_All_scaleddiff_EN_tAge ~ C(genotype) + C(agent)'))
        self.assertFalse(dc.formula_has_age('GrimAgeV2 ~ group + passage + stage'))
        self.assertFalse(dc.formula_has_age('clock ~ group + gestational_age + age_at_diagnosis'))
        self.assertFalse(dc.formula_has_age(None))


class Recompute(unittest.TestCase):
    def test_paired_design_uses_subject_differences(self):
        subjects = {f'A{i}': f's{i}' for i in range(5)} | {f'B{i}': f's{i}' for i in range(5)}
        base = [10, 30, 50, 70, 90]
        values = {f'A{i}': base[i] + 2 + 0.1 * i for i in range(5)} | {f'B{i}': base[i] for i in range(5)}
        case, control = [f'A{i}' for i in range(5)], [f'B{i}' for i in range(5)]
        unpaired = dc.welch([values[g] for g in case], [values[g] for g in control])
        paired = dc.recompute(case, control, values, subjects, {}, ['PAIRED_AS_INDEPENDENT'])
        self.assertGreater(unpaired['p'], 0.5)
        self.assertLess(paired['p'], 0.001)
        self.assertEqual(dc.change(dc.call(unpaired['estimate'], unpaired['p']), dc.call(paired['estimate'], paired['p'])),
                         'SIGNIFICANCE_GAINED')

    def test_repeated_subjects_collapse_to_means(self):
        subjects = {'A0': 's1', 'A1': 's1', 'A2': 's1', 'A3': 's2', 'B0': 's3', 'B1': 's4'}
        values = {'A0': 1.0, 'A1': 2.0, 'A2': 3.0, 'A3': 6.0, 'B0': 0.0, 'B1': 1.0}
        out = dc.recompute(['A0', 'A1', 'A2', 'A3'], ['B0', 'B1'], values, subjects, {}, ['REPEATED_SUBJECT'])
        self.assertEqual(out['n'], [2, 2])
        self.assertAlmostEqual(out['estimate'], 4.0 - 0.5)

    def test_one_subject_arm_is_not_testable(self):
        subjects = {'A0': 's1', 'A1': 's1', 'B0': 's2', 'B1': 's3'}
        values = {'A0': 1.0, 'A1': 2.0, 'B0': 0.0, 'B1': 1.0}
        out = dc.recompute(['A0', 'A1'], ['B0', 'B1'], values, subjects, {}, ['REPEATED_SUBJECT'])
        self.assertEqual(out['status'], 'ONE_SUBJECT_ARM')

    def test_age_adjustment_removes_an_age_driven_difference(self):
        case = [f'A{i}' for i in range(6)]
        control = [f'B{i}' for i in range(6)]
        noise = [0.3, -0.2, 0.1, -0.3, 0.2, -0.1]
        ages = {g: 50 + 4 * i for i, g in enumerate(case)} | {g: 40 + 4 * i for i, g in enumerate(control)}
        values = {g: ages[g] + noise[i % 6] for i, g in enumerate(case + control)}
        unadjusted = dc.welch([values[g] for g in case], [values[g] for g in control])
        out = dc.recompute(case, control, values, {}, ages, ['AGE_IMBALANCE'])
        self.assertAlmostEqual(unadjusted['estimate'], 10.0, places=6)
        self.assertEqual(out['status'], 'OK')
        self.assertLess(abs(out['estimate']), 0.5)
        self.assertAlmostEqual(out['age_slope'], 1.0, places=1)

    def test_age_separating_the_arms_is_not_estimable(self):
        ages = {'A0': 20, 'A1': 22, 'A2': 24, 'B0': 60, 'B1': 62, 'B2': 64}
        values = {g: a + 1.0 for g, a in ages.items()}
        out = dc.recompute(['A0', 'A1', 'A2'], ['B0', 'B1', 'B2'], values, {}, ages, ['AGE_IMBALANCE'])
        self.assertEqual(out['status'], 'AGE_SEPARATES_ARMS')


class Row(unittest.TestCase):
    def test_unflagged_row_does_not_load_predictions(self):
        row = {'gse': 'GSE1', 'row_ordinal': 1, 'clock': 'c', 'final_class': 'VERIFIED', 'case_gsms': ['G0', 'G1'],
               'control_gsms': ['G2', 'G3'], 'recomputed': {'mean_diff': 1.0, 'p': 0.01, 'fdr_family': 0.02}}
        record = dc.check_row(row, 'dnam', None, {}, lambda clock: self.fail('predictions loaded'))
        self.assertEqual(record['findings'], [])
        self.assertEqual(record['released']['call'], 'higher')

    def test_named_age_contrast_is_exempt(self):
        ages = {'G0': 20, 'G1': 22, 'G2': 60, 'G3': 62}
        row = {'gse': 'GSE1', 'row_ordinal': 1, 'clock': 'c', 'final_class': 'VERIFIED', 'case_gsms': ['G2', 'G3'],
               'control_gsms': ['G0', 'G1'], 'comparison_name': 'Old vs young', 'recomputed': {'mean_diff': 1.0, 'p': 0.01}}
        record = dc.check_row(row, 'dnam', None, ages, lambda clock: {})
        self.assertEqual(record['findings'], [])
        self.assertEqual(record['age_declared'], 'name')

    def test_pace_clock_is_exempt(self):
        ages = {'G0': 20, 'G1': 22, 'G2': 60, 'G3': 62}
        row = {'gse': 'GSE1', 'row_ordinal': 1, 'clock': 'DunedinPACE', 'final_class': 'VERIFIED',
               'case_gsms': ['G2', 'G3'], 'control_gsms': ['G0', 'G1'], 'comparison_name': 'Disease vs control',
               'recomputed': {'mean_diff': 1.0, 'p': 0.01}, 'evidence': {'clock_target': 'pace of aging (years per year)'}}
        self.assertEqual(dc.check_row(row, 'dnam', None, ages, lambda clock: {})['findings'], [])
        row['evidence'] = {'clock_target': 'chronological age (years)'}
        self.assertEqual(dc.check_row(row, 'dnam', None, ages, lambda clock: {})['findings'], ['AGE_IMBALANCE'])

    def test_mouse_age_is_checked_only_in_variant_v2(self):
        ages = {'G0': 8, 'G1': 8, 'G2': 20, 'G3': 20}
        row = {'gse': 'GSE1', 'row_ordinal': 1, 'clock': 'c', 'final_class': 'VERIFIED', 'case_gsms': ['G2', 'G3'],
               'control_gsms': ['G0', 'G1'], 'comparison_name': 'Knockout vs wild type', 'recomputed': {'mean_diff': 1.0, 'p': 0.01}}
        self.assertEqual(dc.check_row(row, 'mm_rna', None, ages, lambda clock: {})['findings'], [])
        self.assertEqual(dc.check_row(row, 'mm_rna', None, ages, lambda clock: {}, mouse_age=True)['findings'], ['AGE_IMBALANCE'])

    def test_candidate_keeps_the_released_test_without_design_information(self):
        row = {'gse': 'GSE1', 'row_ordinal': 1, 'clock': 'c', 'final_class': 'VERIFIED', 'case_gsms': ['G0', 'G1'],
               'control_gsms': ['G2', 'G3'], 'recomputed': {'mean_diff': 1.0, 'p': 0.01}}
        record = dc.check_row(row, 'mm_rna', None, {}, lambda clock: self.fail('predictions loaded'))
        self.assertEqual(record['candidate']['method'], 'Welch t-test (released)')
        self.assertEqual(record['candidate']['call'], 'higher')

    def test_candidate_adjusts_for_age_when_ages_are_recorded(self):
        case, control = [f'A{i}' for i in range(6)], [f'B{i}' for i in range(6)]
        noise = [0.3, -0.2, 0.1, -0.3, 0.2, -0.1]
        ages = {g: 50 + 4 * i for i, g in enumerate(case)} | {g: 40 + 4 * i for i, g in enumerate(control)}
        values = {g: ages[g] + noise[i % 6] for i, g in enumerate(case + control)}
        row = {'gse': 'GSE1', 'row_ordinal': 1, 'clock': 'Horvathv1', 'final_class': 'VERIFIED', 'case_gsms': case,
               'control_gsms': control, 'comparison_name': 'Disease vs control', 'recomputed': {'mean_diff': 10.0, 'p': 0.01}}
        record = dc.check_row(row, 'dnam', None, ages, lambda clock: values)
        self.assertIn('least squares', record['candidate']['method'])
        self.assertEqual(record['candidate']['call'], 'ns')
        self.assertEqual(dc.check_row(row, 'mm_rna', None, ages, lambda clock: values)['candidate']['method'],
                         'Welch t-test (released)')

    def test_candidate_falls_back_when_all_samples_have_one_age(self):
        case, control = ['A0', 'A1', 'A2'], ['B0', 'B1', 'B2']
        ages = {g: 30.0 for g in case + control}
        values = {'A0': 5.0, 'A1': 6.0, 'A2': 7.0, 'B0': 1.0, 'B1': 2.0, 'B2': 3.0}
        row = {'gse': 'GSE1', 'row_ordinal': 1, 'clock': 'Horvathv1', 'final_class': 'VERIFIED', 'case_gsms': case,
               'control_gsms': control, 'comparison_name': 'Drug vs vehicle', 'recomputed': {'mean_diff': 4.0, 'p': 0.008}}
        record = dc.check_row(row, 'dnam', None, ages, lambda clock: values)
        self.assertEqual(record['candidate']['method'], 'Welch t-test (released)')
        self.assertEqual(record['candidate']['call'], 'higher')

    def test_candidate_pairs_subjects(self):
        subjects = {f'A{i}': f's{i}' for i in range(5)} | {f'B{i}': f's{i}' for i in range(5)}
        base = [10, 30, 50, 70, 90]
        values = {f'A{i}': base[i] + 2 + 0.1 * i for i in range(5)} | {f'B{i}': base[i] for i in range(5)}
        row = {'gse': 'GSE1', 'row_ordinal': 1, 'clock': 'c', 'final_class': 'VERIFIED',
               'case_gsms': [f'A{i}' for i in range(5)], 'control_gsms': [f'B{i}' for i in range(5)],
               'recomputed': {'mean_diff': 2.2, 'p': 0.9}}
        table = {'field': 'individual', 'source': 'geo', 'subjects': subjects}
        record = dc.check_row(row, 'hs_rna', table, {}, lambda clock: values)
        self.assertIn('paired', record['candidate']['method'])
        self.assertEqual(record['candidate']['call'], 'higher')

    def test_unit_key_ignores_arm_order(self):
        self.assertEqual(dc.unit_key('GSE1', ['a', 'b'], ['c']), dc.unit_key('GSE1', ['c'], ['b', 'a']))


if __name__ == '__main__':
    unittest.main()
