import unittest

import numpy as np
import pandas as pd

from audit_system.input_integrity_2026_09_30 import Correlation, screen_dnam


class CorrelationTest(unittest.TestCase):
    def test_chunked_pairwise_complete_matches_numpy(self):
        rng = np.random.default_rng(0)
        x = rng.random((500, 3))
        y = x[:, 0] * 0.5 + rng.random(500) * 0.5
        x[rng.random((500, 3)) < 0.1] = np.nan
        acc = Correlation(3)
        for start in range(0, 500, 70):
            acc.add(x[start:start + 70], y[start:start + 70])
        for j in range(3):
            keep = np.isfinite(x[:, j])
            self.assertAlmostEqual(acc.r()[j], np.corrcoef(x[keep, j], y[keep])[0, 1], places=10)


class ScreenTest(unittest.TestCase):
    def test_probe_detached_sample_loses_reference_correlation(self):
        # Synthetic fixture: a bimodal reference methylome, two samples near it and one sample whose
        # values were shuffled across probes (same distribution, detached from probe identities).
        rng = np.random.default_rng(1)
        probes = [f'cg{i:08d}' for i in range(4000)]
        reference = pd.Series(np.where(rng.random(4000) < 0.5, 0.1, 0.85) + rng.normal(0, 0.03, 4000), index=probes)
        good = np.clip(reference.to_numpy()[:, None] + rng.normal(0, 0.05, (4000, 2)), 0, 1)
        detached = rng.permutation(good[:, 0])
        frame = pd.DataFrame(np.column_stack([good, detached]), index=pd.Index(probes, name='probe_id'),
                             columns=['GSM1', 'GSM2', 'GSM3'])
        import tempfile, os
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as folder:
            path = f'{folder}/GSE0.csv.gz'
            frame.to_csv(path)
            result = screen_dnam(path, reference).set_index('gsm')
        self.assertGreater(result.loc['GSM1', 'r_reference'], 0.9)
        self.assertGreater(result.loc['GSM2', 'r_reference'], 0.9)
        self.assertLess(abs(result.loc['GSM3', 'r_reference']), 0.1)
        self.assertEqual(result.loc['GSM3', 'features'], 4000)
        self.assertEqual(result.loc['GSM3', 'share_outside_unit'], 0)
        self.assertFalse(result.loc['GSM3', 'input_truncated'])

    def test_truncated_stream_keeps_chunks_read(self):
        import gzip, os, tempfile
        from audit_system import input_integrity_2026_09_30 as ii
        n = 200000
        probes = [f'cg{i:08d}' for i in range(n)]
        frame = pd.DataFrame({'GSM1': np.linspace(0, 1, n)}, index=pd.Index(probes, name='probe_id'))
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as folder:
            whole = gzip.compress(frame.to_csv().encode())
            path = f'{folder}/GSE0.csv.gz'
            with open(path, 'wb') as handle:
                handle.write(whole[:len(whole) * 2 // 3])
            original = ii.CELLS_PER_CHUNK
            ii.CELLS_PER_CHUNK = 1000
            try:
                result = screen_dnam(path, pd.Series(np.linspace(0, 1, n), index=probes))
            finally:
                ii.CELLS_PER_CHUNK = original
        self.assertTrue(result.input_truncated[0])
        self.assertGreater(result.features[0], 0)
        self.assertLess(result.features[0], n)
        self.assertAlmostEqual(result.r_reference[0], 1.0)


if __name__ == '__main__':
    unittest.main()
