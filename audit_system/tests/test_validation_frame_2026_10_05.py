"""Unit tests for validation-frame recovery. All fixtures are synthetic."""
import inspect
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from audit_system import validation_frame_2026_10_05 as vf


def write_jsonl(path, records):
    path.write_text(''.join(json.dumps(r) + '\n' for r in records))


def member(ordinal, **overrides):
    record = {'row_ordinal': ordinal, 'gse_id': f'GSE{ordinal}', 'comparison_name': 'Resveratrol vs DMSO',
              'condition_name': 'Resveratrol', 'clock': 'c', 'status': 'UNIQUE', 'source': 'replay',
              'case_level': 'resveratrol', 'control_level': 'DMSO', 'stratum': None}
    record.update(overrides)
    return record


class ValidationFrameTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='w18-validation-tests-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.membership, self.repairs, self.decisions = [self.root / n for n in ('membership', 'repairs', 'decisions')]
        for p in (self.membership, self.repairs, self.decisions):
            p.mkdir()

    def test_original_selection_excludes_unusable_adjudicated_and_designed(self):
        records = [member(i) for i in range(12)]
        records[2]['source'] = 'design_reconstruction'
        records[3]['status'] = 'AMBIGUOUS'
        records[4]['case_level'], records[4]['control_level'] = 'DMSO', 'resveratrol'
        records[5]['case_level'] = None
        write_jsonl(self.membership / 'rows.jsonl', records)
        write_jsonl(self.repairs / 'repairs.jsonl', [{'row_ordinal': i,
                    'final_class': 'DUPLICATE' if i == 6 else 'PROVISIONAL' if i % 2 else 'VERIFIED'} for i in range(12)])
        write_jsonl(self.decisions / 'decisions-pending-20260929T231958.jsonl',
                    [{'row_ordinal': 1, 'decision': 'ADJUDICATION_INCOMPLETE'}])
        frame, sample, digest = vf.redraw(self.membership, self.repairs, self.decisions, limit=3)
        self.assertEqual(list(frame.row_ordinal), [0, 7, 8, 9, 10, 11])
        self.assertEqual(list(sample.row_ordinal), list(frame.sample(3, random_state=20260930).row_ordinal))
        self.assertEqual(digest, vf.SELECTION_SHA256)
        self.assertEqual(sorted(p.name for p in self.decisions.iterdir()), ['decisions-pending-20260929T231958.jsonl'])

    def test_small_frame_and_missing_source_column(self):
        records = [member(i) for i in range(2)]
        for record in records:
            del record['source']
        write_jsonl(self.membership / 'rows.jsonl', records)
        write_jsonl(self.repairs / 'repairs.jsonl', [{'row_ordinal': i, 'final_class': 'VERIFIED'} for i in range(2)])
        frame, sample, _ = vf.redraw(self.membership, self.repairs, self.decisions)
        self.assertEqual(len(frame), 2)
        self.assertEqual(len(sample), 2)

    def test_cutoff_excludes_validation_and_later_runs_regardless_of_mtime(self):
        names = ['decisions-pending-20260929T231958.jsonl', 'decisions-supported_sample-20261001T162223.jsonl',
                 'decisions-unconfirmed-20261003T010000.jsonl']
        for name in names:
            write_jsonl(self.decisions / name, [])
        self.assertEqual([p.name for p in vf.before_sample(self.decisions)], names[:1])
        write_jsonl(self.decisions / 'decisions-undated.jsonl', [])
        with self.assertRaisesRegex(ValueError, 'Cannot date'):
            vf.before_sample(self.decisions)

    def test_sample_comparison_checks_identity_multiplicity_and_order(self):
        records = [{k: member(i)[k] for k in vf.IDENTITY} for i in range(3)]
        sample = pd.DataFrame(records)
        self.assertTrue(vf.compare_sample(sample, records)['ordered_equal'])
        result = vf.compare_sample(sample, records[::-1])
        self.assertTrue(result['rows_equal'])
        self.assertFalse(result['ordered_equal'])
        records[1]['clock'] = 'other_clock'
        result = vf.compare_sample(sample, records)
        self.assertFalse(result['rows_equal'])
        self.assertEqual(len(result['recorded_only']), 1)
        self.assertEqual(len(result['redrawn_only']), 1)
        self.assertFalse(vf.compare_sample(sample, [records[0]] * 3)['rows_equal'])

    def test_modified_selection_fails_closed(self):
        source = inspect.getsource(vf.labels).replace('random_state=20260930', 'random_state=1')
        with self.assertRaisesRegex(ValueError, 'differs from the recorded'):
            vf.selection_function(source)

    def test_weighted_rate_uses_frame_weights_and_all_decisions(self):
        results = {'a': {'frame_size': 100, 'sample_size': 10, 'decisions': {'LABEL_INVERTED_CONFIRMED': 1,
                    'LABEL_SUPPORTED': 7, 'LABEL_UNCONFIRMED': 2}},
                   'b': {'frame_size': 300, 'sample_size': 10, 'decisions': {'LABEL_INVERTED_CONFIRMED': 2,
                    'LABEL_SUPPORTED': 5, 'LABEL_UNCONFIRMED': 3}}}
        result = vf.weighted_rate(results)
        self.assertAlmostEqual(result['rate'], 0.175)
        self.assertEqual(result['frame_total'], 400)
        self.assertAlmostEqual(result['expected_inversions'], 70)
        self.assertAlmostEqual(result['variance'], (0.25 ** 2 * .1 * .9 + .75 ** 2 * .2 * .8) / 9)

    def test_candidate_selection_ignores_post_sample_chain_and_rejects_overwrites(self):
        for name, created in [('before', '2026-09-30T17:05:21-0400'), ('after', '2026-10-01T23:23:36-0400')]:
            folder = self.root / name
            for sub in ('merged-dnam', 'repairs-dnam', 'verdicts-dnam'):
                (folder / sub).mkdir(parents=True)
            for sub, file in [('merged-dnam', 'rows.jsonl'), ('repairs-dnam', 'repairs.jsonl'),
                              ('verdicts-dnam', 'row_verdicts.jsonl')]:
                write_jsonl(folder / sub / file, [{'synthetic': True}])
            manifest = {'created': created, 'membership_sha256': vf.sha256(folder / 'merged-dnam/rows.jsonl'),
                        'repairs_sha256': vf.sha256(folder / 'repairs-dnam/repairs.jsonl'),
                        'verdicts_sha256': vf.sha256(folder / 'verdicts-dnam/row_verdicts.jsonl')}
            (folder / 'repairs-dnam/manifest.json').write_text(json.dumps(manifest))
        candidates, rejected = vf.saved_chains(self.root, 'dnam', manifest['membership_sha256'])
        self.assertEqual([r[1].name for r in candidates], ['before'])
        self.assertEqual(rejected, [])
        (self.root / 'before/repairs-dnam/repairs.jsonl').write_text('overwritten\n')
        candidates, rejected = vf.saved_chains(self.root, 'dnam', manifest['membership_sha256'])
        self.assertEqual(candidates, [])
        self.assertEqual(len(rejected), 1)


if __name__ == '__main__':
    unittest.main()
