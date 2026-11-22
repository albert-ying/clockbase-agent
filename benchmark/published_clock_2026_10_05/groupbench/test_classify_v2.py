import unittest

from .classify_v2 import classify_unit

FACTORS = ['compound', 'tissue', 'sex', 'background', 'age']


def mouse(gsm, compound='control', sex='Male', tissue='Liver', age=600, subject=None, category='labelled'):
    expert = {'Intervention': 'Control' if compound == 'control' else f'{compound}; 14ppm', 'Sex': sex}
    return gsm, {'expert': expert, 'geo': {'treatment': compound}, 'category': category, 'subject': subject,
                 'coded': {'compound': compound, 'tissue': tissue, 'sex': sex, 'background': 'B6', 'age': age}}


def declared(case='rapamycin', reference='control'):
    return {'declared_factor': 'intervention_compound', 'status': 'agreed',
            'case_rule': [[{'source': 'geo', 'field': 'treatment', 'values': [case]}]],
            'reference_rule': [[{'source': 'geo', 'field': 'treatment', 'values': [reference]}]]}


def unit(arm_a, arm_b, case_is_arm_a=True, direction=1.0, effect=1.0, p=0.01, test='Welch t-test', values=None, **kw):
    clock = {'clock': 'tAge', 'historical': {'pvalue': p, 'direction': direction, 'effect_size': effect,
                                             'statistical_test': test}, 'values': values or {}}
    variant = {'index': 0, 'declared': kw.get('decl', declared()), 'case_is_arm_a': case_is_arm_a, 'clocks': [clock],
               'adjusted': kw.get('adjusted', False), 'series_case_level': kw.get('series_case', []),
               'series_reference_level': kw.get('series_reference', [])}
    return {'unit_id': 'U', 'arm_a': arm_a, 'arm_b': arm_b, 'variants': [variant]}


def factor_of(declared_contrast):
    return 'compound'


class ClassifyV2Tests(unittest.TestCase):
    def setUp(self):
        self.samples = dict([mouse(f'R{i}', 'rapamycin') for i in range(4)] + [mouse(f'C{i}') for i in range(4)] +
                            [mouse(f'F{i}', 'rapamycin', sex='Female') for i in range(4)] +
                            [mouse(f'M{i}', 'metformin') for i in range(4)] +
                            [mouse(f'K{i}', 'rapamycin', tissue='Kidney') for i in range(4)])
        self.values = {**{f'R{i}': 1.0 + i / 10 for i in range(4)}, **{f'C{i}': 0.0 + i / 10 for i in range(4)},
                       **{f'F{i}': 1.0 for i in range(4)}}

    def run_unit(self, u, set_name='A1', factors=FACTORS):
        return classify_unit(u, self.samples, set_name, factors, factor_of)

    def ids(self, prefix):
        return [f'{prefix}{i}' for i in range(4)]

    def test_correct_and_reported_inversion(self):
        good = self.run_unit(unit(self.ids('R'), self.ids('C'), values=self.values))
        self.assertEqual((good['class'], good['error_types']), ('CORRECT', []))
        inverted = self.run_unit(unit(self.ids('R'), self.ids('C'), direction=-1.0, effect=-1.0, values=self.values))
        self.assertEqual(inverted['error_types'], ['REPORTED_RESULT_INVERSION'])

    def test_arm_order_with_right_direction_is_a_record_error_only(self):
        result = self.run_unit(unit(self.ids('C'), self.ids('R'), case_is_arm_a=True, direction=1.0, effect=-1.0,
                                    values=self.values))
        self.assertEqual(result['class'], 'CORRECT')
        self.assertEqual(result['variants'][0]['record_errors'][0]['type'], 'EFFECT_SIGN_INCONSISTENCY')

    def test_impure_mismatch_and_wrong_factor(self):
        impure = self.run_unit(unit(['R0', 'R1', 'M0', 'M1'], self.ids('C'), values=self.values))
        self.assertEqual(impure['error_types'], ['IMPURE_ARM'])
        self.samples.update(dict([mouse(f'G{i}', 'glucose') for i in range(4)]))
        mismatch = self.run_unit(unit(self.ids('M'), self.ids('G'), values=self.values, decl=declared()))
        self.assertEqual(mismatch['error_types'], ['LABEL_MISMATCH'])
        wrong = self.run_unit(unit(self.ids('R'), self.ids('F'), values=self.values))
        self.assertEqual(wrong['error_types'], ['WRONG_FACTOR'])

    def test_confounding_and_study_design_exemption(self):
        confounded = self.run_unit(unit(self.ids('F'), self.ids('C'), values=self.values))
        self.assertEqual(confounded['error_types'], ['CONFOUNDED'])
        design = self.run_unit(unit(self.ids('F'), self.ids('C'), values=self.values,
                                    series_case=self.ids('F'), series_reference=self.ids('C')))
        self.assertEqual(design['class'], 'CORRECT')
        self.assertIn('study-design imbalance in sex', design['variants'][0]['flags'])
        partial = self.run_unit(unit(['R0', 'R1'], self.ids('C'), values=self.values, p=None,
                                     series_case=self.ids('R'), series_reference=self.ids('C')))
        self.assertIn('incomplete case arm', partial['variants'][0]['flags'])
        tissue = self.run_unit(unit(['R0', 'R1', 'K0', 'K1'], ['C0', 'C1', 'C2', 'C3'], values=self.values))
        self.assertEqual(tissue['error_types'], ['CONFOUNDED'])
        self.samples.update(dict([mouse('X0', 'rapamycin', sex='Female')]))
        chance = self.run_unit(unit(['R0', 'R1', 'X0'], ['C0', 'C1', 'C2'], values=self.values))
        self.assertEqual(chance['class'], 'CORRECT')

    def test_subjects_small_arms_and_degenerate(self):
        self.samples.update(dict([mouse('S1', 'rapamycin', subject='a'), mouse('S2', 'rapamycin', subject='a')]))
        repeated = self.run_unit(unit(['S1', 'S2', 'R0'], self.ids('C'), values=self.values))
        self.assertIn('REPEATED_SUBJECT', repeated['error_types'])
        self.samples.update(dict([mouse(f'P{i}', 'rapamycin', subject=f's{i}') for i in range(3)] +
                                 [mouse(f'Q{i}', subject=f's{i}') for i in range(3)]))
        paired = self.run_unit(unit(['P0', 'P1', 'P2'], ['Q0', 'Q1', 'Q2'], values=self.values))
        self.assertIn('PAIRED_AS_INDEPENDENT', paired['error_types'])
        small = self.run_unit(unit(['R0', 'R1'], ['C0', 'C1'], values=self.values))
        self.assertIn('SMALL_ARM_INFERENCE', small['error_types'])
        degenerate = self.run_unit(unit(self.ids('R'), self.ids('C'), p=1, effect=0, values=self.values))
        self.assertIn('DEGENERATE_TEST', degenerate['error_types'])

    def test_unencoded_ambiguous_and_unassignable(self):
        unencoded = self.run_unit(unit(self.ids('R'), self.ids('C'), values=self.values, decl={
            'declared_factor': 'geo_field', 'status': 'agreed',
            'case_rule': [[{'source': 'geo', 'field': 'batch', 'values': ['1']}]],
            'reference_rule': [[{'source': 'geo', 'field': 'batch', 'values': ['2']}]]}))
        self.assertEqual(unencoded['class'], 'NOT_ASSESSABLE_UNENCODED')
        ambiguous = self.run_unit(unit(self.ids('R'), self.ids('C'), values=self.values,
                                       decl={'declared_factor': 'unclear', 'status': 'ambiguous'}))
        self.assertEqual(ambiguous['class'], 'AMBIGUOUS_LABEL')
        self.samples.update(dict([mouse(f'U{i}', 'rapamycin', category='unassignable') for i in range(2)]))
        lost = self.run_unit(unit(['U0', 'U1', 'R0'], self.ids('C'), values=self.values))
        self.assertEqual(lost['class'], 'NOT_ASSESSABLE')

    def test_a2_age_confounding_spares_rate_clocks(self):
        human = {}
        for i in range(6):
            human[f'H{i}'] = {'expert': {'condition': 'HIV'}, 'geo': {'disease': 'HIV'}, 'category': 'labelled',
                              'subject': None, 'coded': {'condition': 'HIV', 'tissue': 'Blood', 'sex': 'M', 'age': 60 + i}}
            human[f'N{i}'] = {'expert': {'condition': 'HC'}, 'geo': {'disease': 'none'}, 'category': 'labelled',
                              'subject': None, 'coded': {'condition': 'HC', 'tissue': 'Blood', 'sex': 'M', 'age': 30 + i}}
        decl = {'declared_factor': 'condition', 'status': 'agreed',
                'case_rule': [[{'source': 'expert', 'field': 'condition', 'values': ['HIV']}]],
                'reference_rule': [[{'source': 'expert', 'field': 'condition', 'values': ['HC']}]]}
        u = unit([f'H{i}' for i in range(6)], [f'N{i}' for i in range(6)], decl=decl)
        result = classify_unit(u, human, 'A2', ['condition', 'tissue', 'sex', 'age'], lambda d: 'condition')
        self.assertEqual(result['error_types'], ['CONFOUNDED_AGE'])
        u['variants'][0]['clocks'][0]['clock'] = 'DunedinPACE'
        result = classify_unit(u, human, 'A2', ['condition', 'tissue', 'sex', 'age'], lambda d: 'condition')
        self.assertEqual(result['class'], 'CORRECT')


if __name__ == '__main__':
    unittest.main()
