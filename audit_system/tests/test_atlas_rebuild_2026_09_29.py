"""Unit tests for the v3 atlas ranking. All fixtures are synthetic."""
import unittest

import numpy as np
import pandas as pd

from audit_system import atlas_rebuild_2026_09_29 as ar


class AtlasRebuildTests(unittest.TestCase):
    def test_hedges_olkin_se(self):
        self.assertAlmostEqual(ar.hedges_olkin_se(0.0, 10, 10), np.sqrt(0.2))
        self.assertAlmostEqual(ar.hedges_olkin_se(1.0, 10, 10), np.sqrt(0.2 + 1 / 40))

    def test_robust_flag_and_rank_bound(self):
        frame = pd.DataFrame({'v3_cohen_d': [-2.0, 2.0, 0.3, -2.0], 'v3_n_case': [10, 10, 10, 10],
                              'v3_n_control': [10, 10, 10, 10], 'v3_fdr': [0.01, 0.01, 0.01, 0.2],
                              'v3_final_class': ['VERIFIED', 'PROVISIONAL', 'VERIFIED', 'VERIFIED'],
                              'intervention': ['a', 'b', 'c', 'd'], 'condition_name': ['a', 'b', 'c', 'd']})
        out = ar.annotate(frame)
        self.assertEqual(out.v3_robust.tolist(), [True, True, False, False])
        self.assertAlmostEqual(out.v3_rank_bound[0], out.v3_rank_bound[1])  # symmetric magnitude
        self.assertGreater(out.v3_rank_bound[0], 0)
        self.assertLess(out.v3_d_ci_high[0], 0)


if __name__ == '__main__':
    unittest.main()
