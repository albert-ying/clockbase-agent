"""Unit tests for numerical membership recovery. All fixtures are synthetic."""
import unittest
import unittest.mock
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system import membership_recovery_2026_09_27 as mr


def fixture():
    ids = [f'GSM{i}' for i in range(12)]
    meta = pd.DataFrame({'treatment': ['drug'] * 4 + ['vehicle'] * 4 + ['naive'] * 4,
                         'tissue': ['liver', 'liver', 'lung', 'lung'] * 3}, index=ids)
    rng = np.random.default_rng(3)
    values = pd.Series(rng.normal(50, 5, 12), index=ids)
    return meta, values


class MembershipRecoveryTests(unittest.TestCase):
    def test_saved_predictions_recover_without_copying_reference_membership(self):
        meta, values = fixture()
        clock = 'fixture_clock'
        case, control = values.iloc[:4], values.iloc[4:8]
        p = stats.ttest_ind(case, control, equal_var=False).pvalue
        effect = mr.effect_sizes(case, control)['cohen_d_pooled']
        row = {'row_ordinal': 0, 'comparison_name': 'drug vs vehicle', 'condition_name': 'drug',
               'clock': clock, 'statistical_test': 'Welch t-test', 'pvalue': p, 'effect_size': effect}
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / 'replay/inputs').mkdir(parents=True)
            source = root / 'replay/inputs/GSE1.csv'
            pd.DataFrame({'gsm_id': values.index, clock + '_mean': values.values}).to_csv(source, index=False)
            state = {'modality': 'hs_rna', 'config': mr.MODALITIES['hs_rna'],
                     'prediction_replay': root / 'replay', 'prediction_inputs': {}}
            with unittest.mock.patch.object(mr, 'GEO', root), unittest.mock.patch.dict(mr._STATE, state, clear=True), \
                    unittest.mock.patch.object(mr, '_metadata', return_value=meta):
                _, records, _ = mr._study_rows('GSE1', [row, {**row, 'pvalue': p * 1.37}])
                self.assertEqual(records[0]['status'], 'UNIQUE')
                self.assertEqual(records[0]['prediction_source'], 'replay_inputs')
                self.assertEqual(records[0]['case_gsms'], sorted(case.index))
                self.assertEqual(records[0]['prediction_input']['sha256'], mr.sha256(source))
                self.assertEqual(records[1]['status'], 'NONE')
                # An available canonical source still takes precedence.
                canonical = root / 'GEO_hs_rna/clock_result/GSE1.csv'
                canonical.parent.mkdir(parents=True)
                pd.DataFrame({'sample_id': values.index, clock: values.values + 2}).to_csv(canonical, index=False)
                sources = list(mr._prediction_sources('GSE1', clock))
                self.assertEqual([name for name, _ in sources], ['clock_result'])
                np.testing.assert_allclose(sources[0][1], values + 2)

    def test_replay_predictions_require_explicit_source(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            state = {'modality': 'hs_rna', 'config': mr.MODALITIES['hs_rna'], 'prediction_replay': None}
            with unittest.mock.patch.object(mr, 'GEO', root), unittest.mock.patch.dict(mr._STATE, state, clear=True):
                self.assertEqual(list(mr._prediction_sources('GSE1', 'clock')), [])
                mr._STATE['prediction_replay'] = root / 'replay'
                with self.assertRaisesRegex(FileNotFoundError, 'required data file not found'):
                    list(mr._prediction_sources('GSE1', 'clock'))

    def test_unique_recovery_with_orientation_from_effect_sign(self):
        meta, values = fixture()
        case, control = values.iloc[:4], values.iloc[4:8]
        p = stats.ttest_ind(case, control, equal_var=False).pvalue
        d = mr.effect_sizes(case, control)['cohen_d_pooled']
        result = mr.recover_row(values, meta, p, d)
        self.assertEqual(result['status'], 'UNIQUE')
        self.assertEqual(result['case_gsms'], sorted(case.index))
        self.assertEqual(result['test'], 'welch')
        self.assertIn('cohen_d_pooled', result['effect_definitions_matched'])
        flipped = mr.recover_row(values, meta, p, -d)
        self.assertEqual(flipped['case_gsms'], sorted(control.index))

    def test_stratified_split_and_student_test(self):
        meta, values = fixture()
        case = values[(meta.treatment == 'drug') & (meta.tissue == 'liver')]
        control = values[(meta.treatment == 'naive') & (meta.tissue == 'liver')]
        p = stats.ttest_ind(case, control, equal_var=True).pvalue
        result = mr.recover_row(values, meta, p, mr.effect_sizes(case, control)['cohen_d_pooled'])
        self.assertIn(result['status'], {'UNIQUE', 'AMBIGUOUS'})
        self.assertEqual(result['stratum'], ['tissue', 'liver'])
        self.assertEqual(sorted(result['case_gsms'] + result['control_gsms']), sorted(case.index.union(control.index)))

    def test_no_candidate_reproduces_arbitrary_p(self):
        meta, values = fixture()
        self.assertEqual(mr.recover_row(values, meta, 0.123456789, 0.5)['status'], 'NONE')

    def test_distinct_memberships_with_equal_p_are_ambiguous(self):
        ids = [f'GSM{i}' for i in range(8)]
        meta = pd.DataFrame({'a': ['x'] * 2 + ['y'] * 2 + ['z'] * 4, 'b': ['u', 'v'] * 4}, index=ids)
        values = pd.Series([1., 2., 1., 2., 5., 6., 5., 6.], index=ids)  # GSM0/1 and GSM2/3 are exchangeable
        p = stats.ttest_ind(values.iloc[:2], values.iloc[4:], equal_var=False).pvalue
        result = mr.recover_row(values, meta, p, -1.0)
        self.assertEqual(result['status'], 'AMBIGUOUS')
        self.assertGreaterEqual(result['matching_memberships'], 2)

    def test_two_sample_t_detection(self):
        for text, expected in [('Welch t-test', True), ("Welch's t-test", True), ('Two-sample t-test', True),
                               ('paired t-test', False), ('One-sample t-test', False), ('linear model', False),
                               ('Mann-Whitney U', False)]:
            self.assertEqual(mr.is_two_sample_t(text), expected, text)


    def test_duplicated_sample_ids_are_dropped_not_fatal(self):
        meta, values = fixture()
        case, control = values.iloc[:4], values.iloc[4:8]
        p = stats.ttest_ind(case, control, equal_var=False).pvalue
        doubled = pd.concat([values, pd.Series([99.0], index=['GSM11'])])  # GSM11 appears twice
        result = mr.recover_row(doubled, meta, p, mr.effect_sizes(case, control)['cohen_d_pooled'])
        self.assertEqual(result['status'], 'UNIQUE')
        self.assertEqual(result['duplicated_sample_ids'], ['GSM11'])

    def test_study_errors_are_recorded_per_row(self):
        rows = [{'row_ordinal': 0, 'comparison_name': 'a vs b', 'condition_name': 'a', 'clock': 'c',
                 'statistical_test': 't-test', 'pvalue': 0.01, 'effect_size': 1.0}]
        with unittest.mock.patch.object(mr, '_study_rows', side_effect=RuntimeError('boom')):
            gse, records, _ = mr._study(('GSE1', rows))
        self.assertEqual(records[0]['status'], 'STUDY_ERROR')
        self.assertIn('boom', records[0]['error'])


if __name__ == '__main__':
    unittest.main()
