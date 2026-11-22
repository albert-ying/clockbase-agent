import unittest

import pandas as pd

from audit_system.followup_judgment_2026_10_18 import admitted, judged_tables, parse, reconcile


def answer(age='OLD', background='HEALTHY', control='UNTREATED'):
    return {'age_at_start': age, 'background': background, 'control': control}


KEY = {'modality': 'mouse RNA-seq', 'gse_id': 'GSE1', 'comparison_name': 'A vs B', 'condition_name': 'A'}


class ReconcileTest(unittest.TestCase):
    def test_agreement_keeps_the_answer_and_disagreement_drops_it(self):
        result = reconcile([answer(), answer(age='ADULT')])
        self.assertIsNone(result['age_at_start'])
        self.assertEqual((result['background'], result['control']), ('HEALTHY', 'UNTREATED'))

    def test_unknown_values_are_not_stated(self):
        self.assertEqual(reconcile([answer(age='AGED'), answer(age='AGED')])['age_at_start'], 'NOT_STATED')


class AdmissionTest(unittest.TestCase):
    def test_lowering_needs_old_healthy_untreated(self):
        self.assertTrue(admitted('age-lowering', answer()))
        self.assertFalse(admitted('age-lowering', answer(age='ADULT')))
        self.assertFalse(admitted('age-lowering', answer(background='DISEASE_MODEL')))
        self.assertFalse(admitted('age-lowering', answer(control='ACTIVE_COMPARATOR')))
        self.assertFalse(admitted('age-lowering', {**answer(), 'age_at_start': None}))

    def test_raising_does_not_ask_for_age(self):
        self.assertTrue(admitted('age-raising', answer(age='YOUNG')))
        self.assertFalse(admitted('age-raising', answer(background='NOT_STATED')))

    def test_judgement_version_keeps_statistical_order(self):
        ranked = pd.DataFrame([{**KEY, 'comparison_name': name, 'direction_v3': 'age-lowering', 'statistical_rank': rank,
                                'condition_category': 'Drug'} for rank, name in enumerate(['x', 'y', 'z'], 1)])
        results = [{**KEY, 'comparison_name': 'x', 'status': 'JUDGED', **answer(age='ADULT')},
                   {**KEY, 'comparison_name': 'y', 'status': 'JUDGED', **answer()},
                   {**KEY, 'comparison_name': 'z', 'status': 'MODEL_CALL_FAILED'}]
        _, versions = judged_tables(ranked, results)
        kept = versions['age-lowering']
        self.assertEqual(kept.comparison_name.tolist(), ['y'])
        self.assertEqual((kept.statistical_rank.iloc[0], kept.judgement_rank.iloc[0]), (2, 1))


class ParseTest(unittest.TestCase):
    def test_salvages_fields_when_evidence_breaks_json(self):
        text = ('{"comparisons": [{"id": 3, "age_at_start": "OLD", "background": "HEALTHY", "control": "UNTREATED", '
                '"evidence": ["age: 22 months]}, {"id": 4, "age_at_start": "YOUNG", "background": "DISEASE_MODEL"')
        parsed = parse(text)
        self.assertEqual(list(parsed), [3])
        self.assertEqual(parsed[3]['age_at_start'], 'OLD')


if __name__ == '__main__':
    unittest.main()
