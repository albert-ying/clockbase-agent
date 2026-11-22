"""Unit tests for intervention classification reconciliation. All fixtures are synthetic."""
import unittest

import pandas as pd

from audit_system import intervention_classification_2026_09_29 as ic


def answer(**overrides):
    base = {'contrast_type': 'INTERVENTION', 'intervention_class': 'GENE_LOSS', 'canonical_target': 'SIRT6',
            'manipulation': 'DECREASE', 'relevance_score': 7, 'relevance_codes': ['IN_VIVO', 'ESTABLISHED_AGING_PATHWAY']}
    base.update(overrides)
    return base


class ClassificationTests(unittest.TestCase):
    def test_target_normalization(self):
        self.assertEqual(ic.normalize_target('Sirt6'), ic.normalize_target('SIRT6'))
        self.assertEqual(ic.normalize_target('mouse Sirt6'), 'sirt6')
        self.assertEqual(ic.normalize_target('caloric-restriction'), 'caloric restriction')
        self.assertEqual(ic.normalize_target(None), 'NA')

    def test_agreement_and_relevance(self):
        out = ic.reconcile([answer(), answer(canonical_target='Sirt6', relevance_score=5, relevance_codes=['IN_VIVO'])])
        self.assertEqual(out['status'], 'AGREED')
        self.assertEqual((out['relevance_score'], out['relevance_score_abs_diff']), (6.0, 2.0))
        self.assertEqual(out['relevance_codes_agreed'], ['IN_VIVO'])
        self.assertEqual(ic.reconcile([answer(), answer(manipulation='INCREASE')])['status'], 'MODELS_DISAGREE')
        self.assertEqual(ic.reconcile([answer(), None])['status'], 'MODEL_ANSWER_MISSING')

    def test_arm_contrast_lists_differing_fields_only(self):
        meta = pd.DataFrame({'genotype': ['WT', 'WT', 'KO', 'KO'], 'tissue': ['liver'] * 4},
                            index=['GSM1', 'GSM2', 'GSM3', 'GSM4'])
        text = ic.arm_contrast(meta, ['GSM3', 'GSM4'], ['GSM1', 'GSM2'])
        self.assertIn("genotype: case=['KO'] control=['WT']", text)
        self.assertNotIn('tissue', text)

    def test_parse(self):
        self.assertEqual(ic.parse('{"comparisons": [{"id": "3", "contrast_type": "AGE"}]}'),
                         {3: {'id': '3', 'contrast_type': 'AGE'}})


    def test_parse_salvages_complete_objects_from_truncated_answer(self):
        text = ('{"comparisons": [\n  {"id": 1, "contrast_type": "AGE", "relevance_codes": ["IN_VIVO"]},\n'
                '  {"id": 2, "contrast_type": "INTERVENTION", "relevance_codes": ["CELL_LINE", "TUM')
        self.assertEqual(list(ic.parse(text)), [1])


    def test_grouping_key_tolerates_granularity_not_substance(self):
        drug_a = answer(intervention_class='DRUG', canonical_target='Decitabine', manipulation='DECREASE')
        drug_b = answer(intervention_class='DRUG', canonical_target='decitabine', manipulation='EXPOSURE')
        self.assertEqual(ic.reconcile([drug_a, drug_b])['status'], 'AGREED')
        combo = answer(canonical_target='DNMT1;DNMT3B')
        self.assertEqual(ic.reconcile([combo, answer(canonical_target='DNMT3B + DNMT1')])['status'], 'AGREED')
        self.assertEqual(ic.reconcile([combo, answer(canonical_target='DNMT1')])['status'], 'MODELS_DISAGREE')
        gain = answer(intervention_class='GENE_MODULATION', manipulation='INCREASE')
        self.assertEqual(ic.reconcile([answer(), gain])['status'], 'MODELS_DISAGREE')


if __name__ == '__main__':
    unittest.main()
