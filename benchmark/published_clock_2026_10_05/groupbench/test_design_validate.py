import unittest

import pandas as pd

from benchmark.published_clock_2026_10_05.groupbench import design_validate as dv

ITEM = {'modality': 'dnam', 'case': ['G0', 'G1', 'G2'], 'control': ['G3', 'G4', 'G5'], 'nonpace': True}


def answer(**fields):
    base = {'material': 'PRIMARY_OR_IN_VIVO', 'shared_individuals': 'NO', 'individuals': {}, 'age_unit': None, 'ages': {},
            'age_or_time_is_the_contrast': False}
    return {**base, **fields}


class ReaderFindings(unittest.TestCase):
    def test_paired_individuals(self):
        individuals = {'G0': 'a', 'G1': 'b', 'G2': 'c', 'G3': 'a', 'G4': 'b', 'G5': 'c'}
        out = dv.reader_findings(ITEM, answer(shared_individuals='YES', individuals=individuals))
        self.assertTrue(out['PAIRED_AS_INDEPENDENT'])
        self.assertFalse(out['REPEATED_SUBJECT'])
        self.assertTrue(out['subject_determined'])

    def test_no_shared_individuals_is_determined_without_identifiers(self):
        out = dv.reader_findings(ITEM, answer(individuals={g: None for g in ITEM['case'] + ITEM['control']}))
        self.assertTrue(out['subject_determined'])
        self.assertFalse(out['REPEATED_SUBJECT'] or out['PAIRED_AS_INDEPENDENT'])

    def test_cannot_tell_is_not_determined(self):
        out = dv.reader_findings(ITEM, answer(shared_individuals='CANNOT_TELL', individuals={'G0': None}))
        self.assertFalse(out['subject_determined'])

    def test_cell_lines_are_not_individuals(self):
        individuals = {'G0': 'HeLa', 'G1': 'HeLa', 'G2': 'HeLa', 'G3': 'HeLa', 'G4': 'HeLa', 'G5': 'MCF7'}
        out = dv.reader_findings(ITEM, answer(material='CELL_LINE', shared_individuals='YES', individuals=individuals))
        self.assertTrue(out['subject_determined'])
        self.assertFalse(out['REPEATED_SUBJECT'])

    def test_ages_are_converted_to_years_for_human_series(self):
        ages = {'G0': 240, 'G1': 252, 'G2': 264, 'G3': 600, 'G4': 612, 'G5': 624}
        out = dv.reader_findings(ITEM, answer(age_unit='months', ages=ages))
        self.assertAlmostEqual(out['age_gap'], -30.0)
        self.assertTrue(out['AGE_IMBALANCE'])

    def test_age_contrast_and_pace_units_are_exempt(self):
        ages = {'G0': 20, 'G1': 21, 'G2': 22, 'G3': 60, 'G4': 61, 'G5': 62}
        self.assertFalse(dv.reader_findings(ITEM, answer(age_unit='years', ages=ages, age_or_time_is_the_contrast=True))['AGE_IMBALANCE'])
        self.assertFalse(dv.reader_findings({**ITEM, 'nonpace': False}, answer(age_unit='years', ages=ages))['AGE_IMBALANCE'])

    def test_unusable_answer(self):
        self.assertIsNone(dv.reader_findings(ITEM, {'material': 'UNCLEAR'}))
        self.assertIsNone(dv.reader_findings(ITEM, None))


class Reference(unittest.TestCase):
    def reader(self, **fields):
        return {'subject_determined': True, 'REPEATED_SUBJECT': False, 'PAIRED_AS_INDEPENDENT': False,
                'AGE_IMBALANCE': False, **fields}

    def test_agreement_gives_the_reference(self):
        readers = [self.reader(PAIRED_AS_INDEPENDENT=True), self.reader(PAIRED_AS_INDEPENDENT=True)]
        self.assertTrue(dv.reference(readers, 'PAIRED_AS_INDEPENDENT'))
        self.assertFalse(dv.reference(readers, 'REPEATED_SUBJECT'))

    def test_disagreement_or_missing_reader_is_unresolved(self):
        self.assertIsNone(dv.reference([self.reader(AGE_IMBALANCE=True), self.reader()], 'AGE_IMBALANCE'))
        self.assertIsNone(dv.reference([self.reader(), None], 'REPEATED_SUBJECT'))

    def test_undetermined_subjects_leave_subject_findings_unresolved(self):
        readers = [self.reader(subject_determined=False), self.reader()]
        self.assertIsNone(dv.reference(readers, 'REPEATED_SUBJECT'))
        self.assertFalse(dv.reference(readers, 'AGE_IMBALANCE'))


class Packages(unittest.TestCase):
    def test_strata(self):
        self.assertEqual(dv.stratum({'PAIRED_AS_INDEPENDENT', 'AGE_IMBALANCE'}), 'SUBJECT')
        self.assertEqual(dv.stratum({'AGE_IMBALANCE'}), 'AGE')
        self.assertEqual(dv.stratum(set()), 'NONE')

    def test_package_shows_geo_fields_and_no_check_output(self):
        raw = pd.DataFrame({'title': ['P1_w0', 'P2_w0', 'P1_w6', 'P2_w6'], 'tissue': ['blood'] * 4,
                            'individual': ['P1', 'P2', 'P1', 'P2']}, index=['G0', 'G1', 'G2', 'G3'])
        unit = {'gse': 'GSE1', 'comparison_name': 'Week 6 vs week 0', 'condition_name': 'treatment',
                'case': ['G2', 'G3'], 'control': ['G0', 'G1']}
        prompt = dv.package(unit, raw, {'title': 'T', 'summary': 'S', 'design': 'D'})
        self.assertIn('tissue=blood', prompt.split('Group A')[0])
        self.assertIn('G2: title=P1_w6; individual=P1', prompt)
        for word in ('REPEATED_SUBJECT', 'PAIRED_AS_INDEPENDENT', 'AGE_IMBALANCE', 'harmonized'):
            self.assertNotIn(word, prompt)


if __name__ == '__main__':
    unittest.main()
