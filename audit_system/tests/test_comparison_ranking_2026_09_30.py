"""Unit tests for comparison ranking. All fixtures are synthetic."""
import json
import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from audit_system import comparison_ranking_2026_09_30 as cr


class ComparisonRankingTests(unittest.TestCase):
    def test_exact_ci_contains_estimate_and_is_wider_for_small_n(self):
        lo, hi = cr.exact_d_ci(1.0, 20, 20)
        self.assertLess(lo, 1.0)
        self.assertGreater(hi, 1.0)
        small = cr.exact_d_ci(10.0, 3, 3)
        self.assertGreater(small[1] - small[0], hi - lo)

    def test_exact_ci_coverage_matches_t_test_at_zero_boundary(self):
        # The lower bound is 0 exactly when the one-sided p of the pooled t-test is 0.025.
        n1 = n2 = 10
        t = stats.t.ppf(0.975, n1 + n2 - 2)
        d = t * np.sqrt(1 / n1 + 1 / n2)
        lo, _ = cr.exact_d_ci(d, n1, n2)
        self.assertAlmostEqual(lo, 0.0, places=6)
        adjusted, _ = cr.exact_d_ci(d, n1, n2, 1 - 0.05 / 5)
        self.assertLess(adjusted, 0.0)


def row(gse, comparison, clock, d, category='Drug', flags='[]', context='in_vivo', fdr=0.001, n=8):
    return {'modality': 'human DNAm', 'gse_id': gse, 'comparison_name': comparison, 'condition_name': comparison,
            'condition_category': category, 'clock': clock, 'intervention': comparison, 'v3_final_class': 'VERIFIED',
            'v3_soft_flags': flags, 'v3_n_case': n, 'v3_n_control': n, 'v3_cohen_d': d, 'v3_mean_diff': d,
            'v3_fdr': fdr, 'v3_p': fdr / 10, 'v3_context': context}


class RankingRunTests(unittest.TestCase):
    def run_ranking(self, rows, classes, audit=None):
        folder = tempfile.mkdtemp(dir=os.environ.get('TMPDIR'))
        pd.DataFrame(rows).to_csv(f'{folder}/atlas.csv', index=False)
        Path(f'{folder}/cls').mkdir()
        with open(f'{folder}/cls/comparisons.jsonl', 'w') as sink:
            for c in classes:
                sink.write(json.dumps(c) + '\n')
        if audit is not None:
            Path(f'{folder}/audit').mkdir()
            with open(f'{folder}/audit/comparisons.jsonl', 'w') as sink:
                for a in audit:
                    sink.write(json.dumps(a) + '\n')
        cr.run(f'{folder}/atlas.csv', f'{folder}/cls', f'{folder}/out', f'{folder}/audit' if audit is not None else None)
        return pd.read_csv(f'{folder}/out/ranked_comparisons_all.csv'), json.loads(Path(f'{folder}/out/summary.json').read_text())

    @staticmethod
    def agreed(gse, comparison, contrast='INTERVENTION', relevance=6):
        return {'modality': 'human DNAm', 'gse_id': gse, 'comparison_name': comparison, 'condition_name': comparison,
                'status': 'AGREED', 'contrast_type': contrast, 'canonical_target': 'x', 'relevance_score': relevance,
                'relevance_score_abs_diff': 0}

    def test_eligibility_selection_adjustment_and_audit(self):
        rows = [row('GSE1', 'A', 'Horvathv1', -1.5), row('GSE1', 'A', 'PhenoAge', -0.2),  # best of k = 2 clocks
                row('GSE2', 'B', 'Horvathv1', 2.0, flags="['ADJUDICATED_TIMEPOINT_CONFOUNDED']"),
                row('GSE3', 'C', 'Horvathv1', 2.0, category='Disease'),
                row('GSE4', 'D', 'Horvathv1', 2.0, context='in_vitro')]
        classes = [self.agreed('GSE1', 'A'), self.agreed('GSE2', 'B'), self.agreed('GSE3', 'C'), self.agreed('GSE4', 'D')]
        ranked, summary = self.run_ranking(rows, classes)
        self.assertEqual(ranked.gse_id.tolist(), ['GSE1'])  # B: design flag; C: contrast mismatch; D: context
        self.assertEqual((ranked.clock[0], ranked.selection_k[0]), ('Horvathv1', 2))
        self.assertAlmostEqual(ranked.ci_level[0], 0.975)
        self.assertIn('GPT-4', summary['context_source'])
        audit = [{'modality': 'human DNAm', 'gse_id': 'GSE1', 'comparison_name': 'A', 'condition_name': 'A',
                  'status': 'NOT_RANKABLE', 'reasons': ['SPECIMEN_EX_VIVO_CULTURE']},
                 {'modality': 'human DNAm', 'gse_id': 'GSE4', 'comparison_name': 'D', 'condition_name': 'D',
                  'status': 'RANKABLE', 'reasons': []}]
        ranked, summary = self.run_ranking(rows, classes, audit)
        self.assertEqual(ranked.gse_id.tolist(), ['GSE4'])  # the audit replaces the GPT-4 label
        self.assertEqual(summary['excluded_by_design_audit'], {'SPECIMEN_EX_VIVO_CULTURE': 1})

    def test_contrast_type_agreement_and_unranked_clock(self):
        rows = [row('GSE5', 'E', 'Horvathv1', -1.5), row('GSE6', 'F', 'Horvathv1', -1.5),
                row('GSE7', 'G', 'GrimAgeV2', -3.0), row('GSE7', 'G', 'Horvathv1', -1.5)]
        # E: models agree on the contrast type but not the target -> eligible; F: contrast types differ -> not.
        disagree = {**self.agreed('GSE5', 'E'), 'status': 'MODELS_DISAGREE', 'canonical_target': None}
        differ = {**self.agreed('GSE6', 'F'), 'status': 'MODELS_DISAGREE', 'contrast_type': None}
        ranked, summary = self.run_ranking(rows, [disagree, differ, self.agreed('GSE7', 'G')])
        self.assertEqual(sorted(ranked.gse_id), ['GSE5', 'GSE7'])
        g = ranked[ranked.gse_id == 'GSE7'].iloc[0]
        self.assertEqual((g.clock, g.selection_k), ('Horvathv1', 1))  # GrimAgeV2 neither ranked nor counted in k
        self.assertEqual(summary['unranked_clock_rows_verified_or_provisional'], 1)


if __name__ == '__main__':
    unittest.main()
