"""Unit tests for the error census. All fixtures are synthetic."""
import unittest

from audit_system import error_census_2026_09_29 as ec

ROW = {'pvalue': 0.01, 'effect_size': 1.2, 'direction': 1}


def repair(reasons, actions=()):
    return {'reasons': list(reasons), 'actions': [{'action': a} for a in actions]}


class CensusTests(unittest.TestCase):
    def test_reproduced_row_with_detached_inputs_is_material(self):
        defects, reproduced = ec.classify(repair(['MEMBERSHIP_RECOVERED', 'INPUT_PROFILE_DETACHED'],
                                                 ['DROP_DETACHED_INPUT_SAMPLES']), ROW)
        self.assertTrue(reproduced)
        self.assertEqual(defects, {'DETACHED_INPUT_PROFILES'})
        self.assertIn('DETACHED_INPUT_PROFILES', ec.MATERIAL)

    def test_inversion_and_direction_field(self):
        defects, _ = ec.classify(repair(['MEMBERSHIP_RECOVERED'], ['SWAP_ARMS_CONFIRMED_INVERSION']),
                                 {**ROW, 'direction': -1})
        self.assertEqual(defects, {'ARMS_INVERTED_CONFIRMED', 'DIRECTION_FIELD_CONTRADICTS_EFFECT'})

    def test_unreproduced_rows_are_counted_as_such(self):
        defects, reproduced = ec.classify(repair(['MEMBERSHIP_UNRECOVERED']), ROW)
        self.assertFalse(reproduced)
        self.assertEqual(defects, {'HISTORICAL_ANALYSIS_NOT_REPRODUCIBLE'})


if __name__ == '__main__':
    unittest.main()
