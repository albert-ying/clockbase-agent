"""Unit tests for linking atlas rows to replayed tests. All fixtures are synthetic."""
import unittest

from audit_system import replay_link_2026_09_27 as link

INDEX = {'Horvathv1': {1.0: {'GSM1'}, 2.0: {'GSM2'}, 3.0: {'GSM3'}, 7.0: {'GSM4'}, 8.0: {'GSM5'}, 9.0: {'GSM6'}},
         'PhenoAge': {11.0: {'GSM1'}, 12.0: {'GSM2'}}}
STUDY = {'GSM1', 'GSM2', 'GSM3', 'GSM4', 'GSM5', 'GSM6'}
ROW = {'pvalue': 0.0123, 'effect_size': 2.5, 'clock': 'Horvathv1'}


def two_sample(first, second, p=0.0123, clock='Horvathv1'):
    return {'fn': 'ttest_ind', 'p': p, 'groups': [{'clock': clock, 'gsms': first}, {'clock': clock, 'gsms': second}]}


class ReplayLinkTests(unittest.TestCase):
    def test_two_sample_orientation_follows_stored_effect(self):
        record = link.explain(ROW, two_sample(['GSM1', 'GSM2', 'GSM3'], ['GSM4', 'GSM5', 'GSM6']), INDEX, STUDY)
        self.assertEqual(record['case_gsms'], ['GSM4', 'GSM5', 'GSM6'])  # higher mean carries the positive effect
        self.assertFalse(record['first_argument_is_case'])

    def test_rejects_wrong_p_clock_or_foreign_samples(self):
        self.assertIsNone(link.explain(ROW, two_sample(['GSM1'], ['GSM4'], p=0.5), INDEX, STUDY))
        self.assertIsNone(link.explain(ROW, two_sample(['GSM1'], ['GSM2'], clock='PhenoAge'), INDEX, STUDY))
        self.assertIsNone(link.explain(ROW, two_sample(['GSM1'], ['GSM99']), INDEX, STUDY))

    def test_model_fit_uses_matched_term_design_column(self):
        test = {'fn': 'sm.OLS', 'pvalues': {'Intercept': 1e-5, 'group[T.treated]': 0.0123, 'age': 0.4},
                'endog': [1.0, 2.0, 3.0, 7.0, 8.0, 9.0], 'exog_names': ['Intercept', 'group[T.treated]', 'age'],
                'exog': [[1, 0, 30], [1, 0, 31], [1, 0, 32], [1, 1, 33], [1, 1, 34], [1, 1, 35]], 'formula': 'y ~ group + age'}
        record = link.explain(ROW, test, INDEX, STUDY)
        self.assertEqual(record['kind'], 'model')
        self.assertEqual(record['case_gsms'], ['GSM4', 'GSM5', 'GSM6'])
        self.assertEqual(record['covariates'], ['age'])

    def test_bucket_candidates_are_found_near_boundaries(self):
        by_bucket = {link.bucket(0.0123): [{'id': 1}]}
        self.assertEqual(link.candidates(by_bucket, 0.0123 * (1 + 1e-10)), [{'id': 1}])


    def test_model_fit_uses_analyzed_members_when_endog_is_absent(self):
        test = {'fn': 'sm.OLS', 'pvalues': {'Intercept': 1e-5, 'group': 0.0123}, 'formula': 'Horvathv1 ~ group',
                'members': ['GSM1', 'GSM2', 'GSM4'], 'clock': None, 'exog_names': ['Intercept', 'group'],
                'exog': [[1, 0], [1, 0], [1, 1]]}
        record = link.explain(ROW, test, INDEX, STUDY)
        self.assertEqual(record['case_gsms'], ['GSM4'])
        self.assertIsNone(link.explain({**ROW, 'clock': 'PhenoAge'}, test, INDEX, STUDY))


if __name__ == '__main__':
    unittest.main()
