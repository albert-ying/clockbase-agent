"""Unit tests for the blood clock-structure analysis. All fixtures are synthetic."""
import unittest

import numpy as np
import pandas as pd

from audit_system import blood_clock_structure_2026_09_30 as bs


def synthetic(seed=0, studies=4, n=60):
    # Two latent factors beyond age: A drives Horvathv1 and Hannum, B drives the GrimAge versions. Studies
    # differ by large offsets, which dominate pooled raw correlations but not within-study residuals.
    rng = np.random.default_rng(seed)
    parts = []
    for s in range(studies):
        age = rng.uniform(20, 80, n)
        a, b = rng.normal(size=n), rng.normal(size=n)
        values = {c: age + rng.normal(0, 3, n) + 30 * s for c in bs.CLOCKS}
        for c in ('Horvathv1', 'Hannum'):
            values[c] = values[c] + 5 * a
        for c in ('GrimAgeV1', 'GrimAgeV2'):
            values[c] = values[c] + 5 * b
        parts.append(pd.DataFrame({'gse_id': f'GSE{s}', 'age': age, **values}))
    return pd.concat(parts, ignore_index=True)


class BloodStructureTests(unittest.TestCase):
    def test_constant_clock_in_a_study_contributes_no_correlation(self):
        frame = synthetic()
        frame.loc[frame.gse_id == 'GSE0', 'Horvathv1'] = 61.702  # imputed constant, as from a failed re-derivation
        residuals = bs.age_residuals(frame)
        self.assertTrue(residuals.loc[frame.gse_id == 'GSE0', 'Horvathv1'].isna().all())
        pooled = bs.pooled_within_study(residuals, frame.gse_id, 'spearman')
        self.assertTrue(np.isfinite(pooled.loc['Horvathv1', 'Hannum']))  # the other studies still inform it

    def test_residuals_remove_age_within_each_study(self):
        frame = synthetic()
        residuals = bs.age_residuals(frame)
        for _, rows in frame.groupby('gse_id'):
            for clock in ('Horvathv1', 'GrimAgeV1'):
                self.assertAlmostEqual(np.corrcoef(residuals.loc[rows.index, clock], rows.age)[0, 1], 0, places=8)

    def test_within_study_pooling_recovers_latent_groups(self):
        frame = synthetic()
        pooled = bs.pooled_within_study(bs.age_residuals(frame), frame.gse_id, 'pearson')
        info = bs.structure(pooled)
        self.assertEqual(info['nearest_neighbour']['Horvathv1'], 'Hannum')
        self.assertEqual(info['nearest_neighbour']['GrimAgeV1'], 'GrimAgeV2')
        self.assertNotEqual(info['clusters_k3']['Horvathv1'], info['clusters_k3']['GrimAgeV1'])
        self.assertTrue(np.allclose(np.diag(pooled), 1))


if __name__ == '__main__':
    unittest.main()
