import unittest

import pandas as pd

from audit_system.design_audit_2026_09_30 import describe_arms, pack, parse, reconcile


def answer(specimen='ORGANISM_DIRECT', tumor='NO', confounds=(), orientation='CASE_IS_CONDITION'):
    return {'orientation': orientation, 'specimen': specimen, 'tumor': tumor, 'confounds': list(confounds)}


class ReconcileTest(unittest.TestCase):
    def test_rankable_only_when_both_admit(self):
        self.assertEqual(reconcile([answer(), answer()])['status'], 'RANKABLE')

    def test_culture_disagreement_and_tumor_block(self):
        self.assertEqual(reconcile([answer(), answer('EX_VIVO_CULTURE')])['reasons'], ['SPECIMEN_DISAGREE'])
        self.assertEqual(reconcile([answer('EX_VIVO_CULTURE'), answer('EX_VIVO_CULTURE')])['reasons'],
                         ['SPECIMEN_EX_VIVO_CULTURE'])
        self.assertEqual(reconcile([answer(tumor='YES'), answer(tumor='YES')])['reasons'], ['TUMOR_YES'])

    def test_orientation_must_be_confirmed_by_both(self):
        self.assertEqual(reconcile([answer(orientation='REVERSED'), answer(orientation='REVERSED')])['reasons'],
                         ['ORIENTATION_REVERSED'])
        self.assertEqual(reconcile([answer(), answer(orientation='UNCLEAR')])['reasons'], ['ORIENTATION_DISAGREE'])

    def test_a_single_model_confound_blocks(self):
        result = reconcile([answer(confounds=['TIME_POINT', 'SEX']), answer(confounds=['SEX'])])
        self.assertEqual(result['reasons'], ['CONFOUND_AGREED_SEX', 'CONFOUND_ONE_MODEL_TIME_POINT'])

    def test_unknown_values_are_unclear(self):
        result = reconcile([answer('IN_VIVO'), answer('IN_VIVO')])
        self.assertEqual(result['reasons'], ['SPECIMEN_UNCLEAR'])
        self.assertIn('CONFOUNDS_MISSING', reconcile([{'orientation': 'CASE_IS_CONDITION', 'specimen': 'ORGANISM_DIRECT',
                                                        'tumor': 'NO'}, answer()])['reasons'])


class ParseTest(unittest.TestCase):
    def test_salvages_fields_when_evidence_breaks_json(self):
        text = ('{"comparisons": [{"id": 3, "orientation": "CASE_IS_CONDITION", "specimen": "EX_VIVO_CULTURE", "tumor": "NO", '
                '"confounds": [], "evidence": ["cultured in "GM-CSF" media"]}, {"id": 4, "orientation": "REVERSED", '
                '"specimen": "ORGANISM_DIRECT", "tumor": "NO", "confounds": ["TIME_POINT"], "evidence": ["4h"]}')
        parsed = parse(text)
        self.assertEqual(parsed[3]['specimen'], 'EX_VIVO_CULTURE')
        self.assertEqual((parsed[4]['confounds'], parsed[4]['orientation']), (['TIME_POINT'], 'REVERSED'))


class FormatTest(unittest.TestCase):
    def test_shared_protocol_shown_once_and_varying_fields_per_sample(self):
        records = pd.DataFrame({'title': ['WT 4h', 'KO 4h', 'WT 0h'],
                                'growth_protocol_ch1': ['RPMI + GM-CSF'] * 3,
                                'characteristics_ch1': ['genotype: WT', 'genotype: KO', 'genotype: WT']},
                               index=['GSM1', 'GSM2', 'GSM3'])
        text = describe_arms(records, ['GSM2'], ['GSM1'])
        self.assertIn('growth_protocol_ch1: RPMI + GM-CSF', text.split('Case arm')[0])
        self.assertIn('GSM2: title=KO 4h; characteristics_ch1=genotype: KO', text)

    def test_pack_limits(self):
        items = [{'text': 'x' * 10000} for _ in range(5)]
        self.assertEqual([len(b) for b in pack(items)], [2, 2, 1])


if __name__ == '__main__':
    unittest.main()
