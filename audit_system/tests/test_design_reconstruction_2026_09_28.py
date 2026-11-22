"""Unit tests for design reconstruction rule handling. All fixtures are synthetic."""
import unittest

import pandas as pd

from audit_system import design_reconstruction_2026_09_28 as dr

META = pd.DataFrame({'genotype': ['WT', 'WT', 'WT', 'KO', 'KO', 'KO'], 'tissue': ['liver'] * 3 + ['liver', 'lung', 'liver']},
                    index=[f'GSM{i}' for i in range(6)])


class DesignReconstructionTests(unittest.TestCase):
    def test_rules_apply_literally(self):
        self.assertEqual(dr.apply_rule(META, {'genotype': ['KO'], 'tissue': ['liver']}), ['GSM3', 'GSM5'])
        self.assertIsNone(dr.apply_rule(META, {'unknown': ['x']}))

    def test_agreement_requires_identical_sets(self):
        a = {0: {'id': 0, 'case': {'genotype': ['KO']}, 'control': {'genotype': ['WT']}}}
        b = {0: {'id': 0, 'case': {'genotype': 'KO'}, 'control': {'genotype': ['WT'], 'tissue': ['liver']}}}
        self.assertEqual(dr.reconcile(META, {'m1': a, 'm2': b}, 0)['status'], 'AGREED')
        c = {0: {'id': 0, 'case': {'genotype': ['KO'], 'tissue': ['liver']}, 'control': {'genotype': ['WT']}}}
        self.assertEqual(dr.reconcile(META, {'m1': a, 'm2': c}, 0)['status'], 'MODELS_DISAGREE')

    def test_infeasible_overlap_and_parse(self):
        self.assertEqual(dr.reconcile(META, {'m1': {0: {'id': 0, 'feasible': False}}, 'm2': {}}, 0)['status'],
                         'INFEASIBLE_OR_MISSING')
        overlap = {0: {'id': 0, 'case': {'tissue': ['liver']}, 'control': {'genotype': ['WT']}}}
        self.assertEqual(dr.reconcile(META, {'m1': overlap, 'm2': overlap}, 0)['status'], 'ARMS_OVERLAP')
        self.assertEqual(dr.parse_rules('x {"comparisons": [{"id": 1, "feasible": false}]} y'), {1: {'id': 1, 'feasible': False}})


    def test_ids_from_json_ledger_are_normalized(self):
        a = {'0': {'id': 0, 'case': {'genotype': ['KO']}, 'control': {'genotype': ['WT']}}}
        b = {0: {'id': 0, 'case': {'genotype': ['KO']}, 'control': {'genotype': ['WT']}}}
        self.assertEqual(dr.reconcile(META, {'m1': a, 'm2': b}, 0)['status'], 'AGREED')


if __name__ == '__main__':
    unittest.main()
