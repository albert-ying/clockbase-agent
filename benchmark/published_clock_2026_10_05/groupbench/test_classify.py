import unittest

from .classify import classify, completeness, total_variation

FACTORS = ['intervention', 'tissue', 'sex', 'strain', 'age']


def sample(intervention='Control', tissue='Liver', sex='Male', strain='B6', age=600, kind=None):
    return {'intervention': intervention, 'tissue': tissue, 'sex': sex, 'strain': strain, 'age': age,
            'type': kind or ('Control' if intervention == 'Control' else 'Compound')}


def run(case, control, labels):
    return classify(case, control, labels, FACTORS, 'type', {'Control'}, age_factor='age')


class ClassifyTests(unittest.TestCase):
    def setUp(self):
        self.labels = {}
        for i in range(4):
            self.labels[f'C{i}'] = sample()
            self.labels[f'R{i}'] = sample('Rapamycin; 14ppm')
            self.labels[f'R2{i}'] = sample('Rapamycin; 42ppm')
            self.labels[f'A{i}'] = sample('Acarbose; 1000ppm')
            self.labels[f'Y{i}'] = sample(age=90)
            self.labels[f'K{i}'] = sample(tissue='Kidney')
            self.labels[f'F{i}'] = sample('Rapamycin; 14ppm', sex='Female')

    def ids(self, prefix, n=4):
        return [f'{prefix}{i}' for i in range(n)]

    def test_correct_intervention_and_age_contrasts(self):
        self.assertEqual(run(self.ids('R'), self.ids('C'), self.labels)['class'], 'CORRECT')
        result = run(self.ids('C'), self.ids('Y'), self.labels)
        self.assertEqual((result['class'], result['contrast_factor']), ('CORRECT', 'age'))

    def test_inversions(self):
        self.assertEqual(run(self.ids('C'), self.ids('R'), self.labels)['class'], 'POSSIBLE_INVERSION')
        self.assertEqual(run(self.ids('Y'), self.ids('C'), self.labels)['class'], 'POSSIBLE_INVERSION')

    def test_confounded_mixed_and_no_contrast(self):
        self.assertEqual(run(self.ids('F'), self.ids('C'), self.labels)['class'], 'CONFOUNDED')
        self.assertEqual(run(self.ids('R', 2) + self.ids('A', 2), self.ids('C'), self.labels)['class'], 'MIXED_ARM')
        self.assertEqual(run(self.ids('C', 2), ['C2', 'C3'], self.labels)['class'], 'NO_CONTRAST')

    def test_pooled_doses_and_tissue_balance(self):
        result = run(self.ids('R', 2) + self.ids('R2', 2), self.ids('C'), self.labels)
        self.assertEqual(result['class'], 'CORRECT')
        self.assertIn('pooled_doses_case', result['flags'])
        result = run(self.ids('R', 2) + ['K0', 'K1'], self.ids('C', 2) + ['K2', 'K3'], self.labels)
        self.assertEqual(result['class'], 'MIXED_ARM')

    def test_coverage_and_completeness(self):
        result = run(self.ids('R') + ['X1', 'X2'], self.ids('C'), self.labels)
        self.assertEqual(result['class'], 'NOT_ASSESSABLE')
        self.assertEqual(completeness(self.ids('R', 2), self.labels, self.labels, FACTORS), 0.5)
        self.assertAlmostEqual(total_variation(['a', 'a', 'b', 'b'], ['a', 'a', 'a', 'a']), 0.5)


if __name__ == '__main__':
    unittest.main()
