import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system import design_estimand_2026_10_10 as de


def frame(**columns):
    n = len(next(iter(columns.values())))
    return pd.DataFrame(columns, index=[f'G{i}' for i in range(n)])


def series(subjects=None, ages=None):
    return {'subject_field': 'individual' if subjects else None, 'subjects': subjects or {},
            'subject_source': 'GEO sample field' if subjects else None, 'age_field': 'harmonized.age' if ages else None,
            'ages': ages or {}, 'files': []}


class Tables(unittest.TestCase):
    def test_small_constant_individual_and_paired_arms_keep_minimum_reason(self):
        # Synthetic software fixtures cover each arm below the minimum.
        for n in range(1, de.MIN_UNITS):
            for left, right in ((n, n), (n, de.MIN_UNITS), (de.MIN_UNITS, n)):
                a, b = [float(left)] * left, [float(right)] * right
                units_a = {str(i): (v, None, []) for i, v in enumerate(a)}
                units_b = {str(i): (v, None, []) for i, v in enumerate(b)}
                for result in (de.two_sample(a, b), de.paired(units_a, units_b)[0]):
                    self.assertNotIn('undefined_reason', result)
                    self.assertFalse(set(result) & {'p', 't', 'cohen_d'})

    def test_constant_guard_covers_each_estimand(self):
        # All inputs in this test are synthetic.
        units_a = {str(i): (2., 20. + i, []) for i in range(4)}
        units_b = {str(i): (5., 20. + i, []) for i in range(4)}
        results = [de.two_sample([2.] * 4, [5.] * 4),
                   de.paired(units_a, units_b)[0],
                   de.age_adjusted(units_a, units_b)[1]]
        for result in results:
            self.assertEqual(result['undefined_reason'], 'outputs constant within both arms')
            self.assertEqual(result['mean_diff'], -3.)
            self.assertFalse(set(result) & {'cohen_d', 't', 'p'})
        self.assertFalse(de.constant_output([], [1., 1.]))
        self.assertTrue(de.constant_output([1.], [2.]))
        self.assertTrue(de.constant_output([1., 1. + 1e-13], [1., 1.]))
        self.assertFalse(de.constant_output([1., 1. + 1e-10], [1., 1.]))

    def test_near_constant_age_fit_is_invariant_to_response_baseline(self):
        # Exactly representable synthetic variation isolates intercept cancellation.
        a = {str(i): (float(2 ** -30 * (i * i + 2)), 20. + i, []) for i in range(5)}
        b = {str(i): (float(2 ** -30 * (i * i)), 20. + i, []) for i in range(5)}
        shifted = lambda units: {k: (v - 1., age, gsms) for k, (v, age, gsms) in units.items()}
        first = de.age_adjusted(a, b)[1]
        second = de.age_adjusted(shifted(a), shifted(b))[1]
        for field in ('mean_diff', 'cohen_d', 'p', 't'):
            self.assertEqual(first[field], second[field])

    def test_named_geo_field_with_repeats_is_the_individual(self):
        raw = frame(title=list('abcdef'), individual=['A02', 'A05', 'A12', 'A02', 'A05', 'A12'])
        field, subjects = de.subject_table(raw)
        self.assertEqual((field, subjects['G3']), ('individual', 'A02'))

    def test_descriptive_status_and_unique_fields_do_not_qualify(self):
        self.assertEqual(de.subject_table(frame(donor=['healthy male', 'MDD diagnosed male', 'mix', 'healthy male', 'mix', 'mix'])),
                         (None, {}))
        self.assertEqual(de.subject_table(frame(**{'subject status': ['1', '2', '1', '2']})), (None, {}))
        self.assertEqual(de.subject_table(frame(**{'patient id': ['P1', 'P2', 'P3', 'P4']})), (None, {}))
        self.assertEqual(de.subject_table(None), (None, {}))

    def test_finest_named_field_wins(self):
        raw = frame(donor=['1', '1', '1', '2', '2', '2'], **{'patient id': ['a1', 'a1', 'a2', 'b1', 'b1', 'b2']})
        self.assertEqual(de.subject_table(raw)[0], 'patient id')

    def test_harmonized_age_first_then_geo_age_at_sampling_for_humans(self):
        processed = frame(**{'harmonized.age': ['30', '40', None]})
        self.assertEqual(de.age_table(frame(age=['31', '41', '51']), processed, True),
                         ('harmonized.age', {'G0': 30.0, 'G1': 40.0}))
        raw = frame(**{'age_sampling': ['35', '37 years', '48', '31'], 'age at diagnosis': ['20', '21', '22', '23']})
        empty = frame(**{'harmonized.age': [None, None, None, None]})
        field, ages = de.age_table(raw, empty, True)
        self.assertEqual((field, ages['G1']), ('age_sampling', 37.0))
        self.assertEqual(de.age_table(raw, empty, False), (None, {}))
        other = frame(**{'age at diagnosis': ['20', '21', '22', '23'], 'gestational age': ['38', '39', '40', '41'],
                         'age': ['young adult', 'elderly adult', 'young adult', 'elderly adult']})
        self.assertEqual(de.age_table(other, None, True), (None, {}))

    def test_load_series_reads_both_sample_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            geo = Path(tmp)
            for sub, table in (('processed_metadata', frame(**{'harmonized.age': ['30', '31', '60', '61']})),
                               ('processed_metadata_preharmonize', frame(individual=['a', 'b', 'a', 'b']))):
                (geo / 'F' / sub).mkdir(parents=True)
                table.rename_axis('GSM').to_csv(geo / 'F' / sub / 'GSE1.csv')
            loaded = de.load_series(geo, 'F', 'GSE1', True)
            self.assertEqual((loaded['subject_field'], loaded['age_field'], len(loaded['files'])), ('individual', 'harmonized.age', 2))
            self.assertEqual(de.load_series(geo, 'F', 'GSE2', True)['subjects'], {})


class Findings(unittest.TestCase):
    def test_repeated_paired_and_independent(self):
        self.assertEqual(de.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4'],
                                             {'G0': 's1', 'G1': 's1', 'G2': 's2', 'G3': 's3', 'G4': 's4'}), ['REPEATED_INDIVIDUAL'])
        self.assertEqual(de.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'],
                                             {'G0': 's1', 'G1': 's2', 'G2': 's3', 'G3': 's1', 'G4': 's2', 'G5': 's3'}), ['PAIRED_ARMS'])
        self.assertEqual(de.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], {f'G{i}': f's{i}' for i in range(6)}), [])

    def test_one_individual_or_unknown_individuals_give_no_finding(self):
        self.assertEqual(de.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], {f'G{i}': 's1' for i in range(6)}), [])
        self.assertEqual(de.subject_findings(['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], {'G0': 's1'}), [])

    def test_age_declared_by_names_or_classification(self):
        self.assertTrue(de.age_declared(['Old vs young liver']))
        self.assertFalse(de.age_declared(['Fast vs Late rebound post-romidepsin', 'Disease vs age-matched controls']))
        self.assertFalse(de.age_declared(['p53 knockout vs wild type', None, float('nan')]))
        self.assertTrue(de.age_declared(['Week 6 vs baseline'], ['INTERVENTION', 'DEVELOPMENT_OR_TIME']))

    def test_overflow_values(self):
        values = pd.Series({'G0': 34.0, 'G1': 1.05e39, 'G2': -2000.0, 'G3': float('nan')})
        self.assertEqual(de.overflow(values, ['G0', 'G1', 'G2', 'G3', 'G9']), ['G1', 'G2'])


class Estimate(unittest.TestCase):
    def test_independent_samples_without_ages_keep_the_declared_estimand(self):
        values = pd.Series({f'G{i}': float(i) for i in range(6)})
        self.assertIsNone(de.estimate(values, ['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], series(), 'mm_rna', 'c'))
        ages = {f'G{i}': 8.0 + i for i in range(6)}
        self.assertIsNone(de.estimate(values, ['G0', 'G1', 'G2'], ['G3', 'G4', 'G5'], series(ages=ages), 'mm_rna', 'c'))

    def test_paired_arms_use_within_individual_differences(self):
        subjects = {f'A{i}': f's{i}' for i in range(5)} | {f'B{i}': f's{i}' for i in range(5)}
        base = [10, 30, 50, 70, 90]
        values = pd.Series({f'A{i}': base[i] + 2 + 0.1 * i for i in range(5)} | {f'B{i}': float(base[i]) for i in range(5)})
        case, control = [f'A{i}' for i in range(5)], [f'B{i}' for i in range(5)]
        out = de.estimate(values, case, control, series(subjects), 'hs_rna', 'c')
        expected = stats.ttest_rel([values[g] for g in case], [values[g] for g in control])
        self.assertEqual(out['design']['method'], 'PAIRED')
        self.assertEqual((out['recomputed']['n_case'], out['recomputed']['n_control']), (5, 5))
        self.assertAlmostEqual(out['recomputed']['p'], expected.pvalue)
        self.assertAlmostEqual(out['recomputed']['mean_diff'], 2.2)
        self.assertGreater(stats.ttest_ind([values[g] for g in case], [values[g] for g in control], equal_var=False).pvalue, 0.5)

    def test_fewer_than_three_pairs_give_no_p(self):
        subjects = {'A0': 's0', 'A1': 's1', 'B0': 's0', 'B1': 's1'}
        values = pd.Series({'A0': 3.0, 'A1': 5.0, 'B0': 1.0, 'B1': 2.0})
        out = de.estimate(values, ['A0', 'A1'], ['B0', 'B1'], series(subjects), 'dnam', 'Horvathv1')
        self.assertEqual(out['design']['pairs'], 2)
        self.assertNotIn('p', out['recomputed'])

    def test_repeated_individuals_collapse_to_means(self):
        subjects = {'A0': 's1', 'A1': 's1', 'A2': 's1', 'A3': 's2', 'A4': 's5', 'B0': 's3', 'B1': 's4', 'B2': 's6'}
        values = pd.Series({'A0': 1.0, 'A1': 2.0, 'A2': 3.0, 'A3': 6.0, 'A4': 7.0, 'B0': 0.0, 'B1': 1.0, 'B2': 2.0})
        out = de.estimate(values, ['A0', 'A1', 'A2', 'A3', 'A4'], ['B0', 'B1', 'B2'], series(subjects), 'mm_rna', 'c')
        self.assertEqual(out['design']['method'], 'UNIT_MEANS')
        self.assertEqual((out['recomputed']['n_case'], out['recomputed']['n_control']), (3, 3))
        self.assertAlmostEqual(out['recomputed']['mean_diff'], 5.0 - 1.0)
        self.assertAlmostEqual(out['recomputed']['p'], stats.ttest_ind([2.0, 6.0, 7.0], [0.0, 1.0, 2.0], equal_var=False).pvalue)

    def test_one_individual_per_arm_gives_no_p(self):
        subjects = {'A0': 's1', 'A1': 's1', 'A2': 's1', 'B0': 's2', 'B1': 's2', 'B2': 's2'}
        values = pd.Series({'A0': 1.0, 'A1': 2.0, 'A2': 3.0, 'B0': 0.0, 'B1': 1.0, 'B2': 2.0})
        out = de.estimate(values, ['A0', 'A1', 'A2'], ['B0', 'B1', 'B2'], series(subjects), 'mm_rna', 'c')
        self.assertEqual((out['recomputed']['n_case'], out['recomputed']['n_control']), (1, 1))
        self.assertNotIn('p', out['recomputed'])

    def test_age_enters_for_human_age_clocks(self):
        case, control = [f'A{i}' for i in range(6)], [f'B{i}' for i in range(6)]
        noise = [0.3, -0.2, 0.1, -0.3, 0.2, -0.1]
        ages = {g: 50.0 + 4 * i for i, g in enumerate(case)} | {g: 40.0 + 4 * i for i, g in enumerate(control)}
        values = pd.Series({g: ages[g] + noise[i % 6] for i, g in enumerate(case + control)})
        out = de.estimate(values, case, control, series(ages=ages), 'dnam', 'Horvathv1', ['Disease vs control'])
        self.assertEqual(out['design']['method'], 'AGE_ADJUSTED')
        self.assertLess(abs(out['recomputed']['mean_diff']), 0.5)
        self.assertAlmostEqual(out['design']['age_slope'], 1.0, places=1)
        self.assertGreater(out['recomputed']['p'], 0.05)
        for exempt in ({'modality': 'mm_rna'}, {'clock': 'DunedinPACE'}, {'names': ['Old vs young']},
                       {'contrast_types': ['AGE', 'OTHER']}):
            args = {'modality': 'dnam', 'clock': 'Horvathv1', 'names': ['Disease vs control'], 'contrast_types': [], **exempt}
            self.assertIsNone(de.estimate(values, case, control, series(ages=ages), **args))

    def test_age_adjustment_matches_ordinary_least_squares(self):
        rng = np.random.default_rng(1)
        case, control = [f'A{i}' for i in range(8)], [f'B{i}' for i in range(9)]
        ages = {g: float(a) for g, a in zip(case + control, rng.uniform(30, 70, 17))}
        values = pd.Series({g: 0.9 * ages[g] + (3.0 if g in case else 0.0) + rng.normal(0, 2) for g in case + control})
        out = de.estimate(values, case, control, series(ages=ages), 'hs_rna', 'score')['recomputed']
        x = np.array([[1.0, float(g in case), ages[g]] for g in case + control])
        y = values[case + control].to_numpy()
        beta, residuals, *_ = np.linalg.lstsq(x, y, rcond=None)
        sigma2 = residuals[0] / (len(y) - 3)
        t = beta[1] / np.sqrt(sigma2 * np.linalg.inv(x.T @ x)[1, 1])
        self.assertAlmostEqual(out['mean_diff'], beta[1])
        self.assertAlmostEqual(out['t'], t)
        self.assertAlmostEqual(out['p'], 2 * stats.t.sf(abs(t), len(y) - 3))
        self.assertAlmostEqual(out['cohen_d'], beta[1] / np.sqrt(sigma2))

    def test_age_separating_the_arms_keeps_the_statistic_with_a_flag(self):
        ages = {'A0': 20.0, 'A1': 22.0, 'A2': 24.0, 'B0': 60.0, 'B1': 62.0, 'B2': 64.0}
        values = pd.Series({g: a + 1.0 for g, a in ages.items()})
        out = de.estimate(values, ['A0', 'A1', 'A2'], ['B0', 'B1', 'B2'], series(ages=ages), 'dnam', 'Horvathv1')
        self.assertIsNone(out['recomputed'])
        self.assertEqual(out['soft'], ['AGE_SEPARATES_ARMS'])

    def test_one_age_for_all_samples_keeps_the_declared_estimand(self):
        ages = {g: 30.0 for g in ('A0', 'A1', 'A2', 'B0', 'B1', 'B2')}
        values = pd.Series({'A0': 5.0, 'A1': 6.0, 'A2': 7.0, 'B0': 1.0, 'B1': 2.0, 'B2': 3.0})
        self.assertIsNone(de.estimate(values, ['A0', 'A1', 'A2'], ['B0', 'B1', 'B2'], series(ages=ages), 'dnam', 'Horvathv1'))

    def test_repeated_individuals_with_ages_fit_age_on_individual_means(self):
        subjects = {f'A{i}': f'a{i // 2}' for i in range(8)} | {f'B{i}': f'b{i // 2}' for i in range(8)}
        ages = {f'A{i}': 40.0 + 5 * (i // 2) for i in range(8)} | {f'B{i}': 45.0 + 5 * (i // 2) for i in range(8)}
        values = pd.Series({g: a + (2.0 if g.startswith('A') else 0.0) + 0.1 * int(g[1]) for g, a in ages.items()})
        out = de.estimate(values, [f'A{i}' for i in range(8)], [f'B{i}' for i in range(8)], series(subjects, ages), 'dnam', 'PhenoAge')
        self.assertEqual(out['design']['method'], 'AGE_ADJUSTED_UNIT_MEANS')
        self.assertEqual((out['recomputed']['n_case'], out['recomputed']['n_control']), (4, 4))

    def test_samples_without_a_value_are_ignored(self):
        subjects = {'A0': 's0', 'A1': 's1', 'A2': 's2', 'B0': 's0', 'B1': 's1', 'B2': 's2', 'B3': 's3'}
        values = pd.Series({'A0': 3.0, 'A1': 5.0, 'A2': 9.0, 'B0': 1.0, 'B1': 2.0, 'B2': 4.0, 'B3': float('nan')})
        out = de.estimate(values, ['A0', 'A1', 'A2', 'A9'], ['B0', 'B1', 'B2', 'B3'], series(subjects), 'mm_rna', 'c')
        self.assertEqual((out['design']['samples_case'], out['design']['samples_control'], out['design']['pairs']), (3, 3, 3))


if __name__ == '__main__':
    unittest.main()
