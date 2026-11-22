"""Disclosed constructed software fixtures, not empirical scientific evidence."""
import copy
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from audit_system.corpus_repairs_2026_09_20 import build_repair_candidates, verify_corpus_repairs
from audit_system.independent_statistics_2026_09_20 import verify_calculations
from audit_system.legacy_analysis_review_2026_09_19 import SCIENTIFIC_CHECKS
from audit_system.statistical_rerun_2026_09_20 import recompute_contrast, record_hash


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class CorpusRepairTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='.corpus-repairs-test-')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root / 'original.csv'
        self.header = ['gse_id', 'clock', 'statistical_test', 'p_value', 'effect_size', 'fdr',
                       'direction', 'trust_score', 'neg_log10_pvalue', 'is_significant', 'notes']
        self.values = [
            ['GSE1', 'clock', 'Paired t-test', '1.2', '2', '0.03', 'lower', '61', '2', 'True', 'original, quoted\nmultiline'],
            ['GSE1', 'clock', 'Paired t-test', 'NaN', 'undefined', '', '', '', '', '', 'missing values are not fabricated numbers'],
            ['GSE1', 'clock', 'Paired t-test', '0.04', '<function result>', '0.08', 'lower', '61', '', '', 'unchanged literal source'],
        ]
        self.write_source()
        self.case = ['GSM1', 'GSM2', 'GSM3']; self.control = ['GSM4', 'GSM5', 'GSM6']
        self.design = {'study_id': 'GSE1', 'frozen_before_saved_outcomes': True,
            'comparisons': [{'id': 'paired', 'mode': 'inferential', 'method': 'Paired t-test',
                'execution_design_status': 'PASS', 'case_ids': self.case, 'control_ids': self.control,
                'independence': {'status': 'PASS'}, 'limitations': ['Constructed software fixture']}]}
        frame = pd.DataFrame({'clock': [4., 8., 7., 1., 4., 5.]}, index=self.case+self.control)
        self.predictions = self.root / 'predictions.csv'; frame.to_csv(self.predictions, index_label='sample_id')
        pairs = {'pairs': [{'unit_id': str(i), 'case_id': a, 'control_id': b}
                           for i, (a, b) in enumerate(zip(self.case, self.control))]}
        candidate = recompute_contrast(frame, self.case, self.control, 'clock', 'Paired t-test', design=pairs)
        self.calculations = [{'comparison_id': 'paired', 'clock': 'clock', 'case_ids': self.case,
            'control_ids': self.control, 'source_design_sha256': record_hash(self.design),
            'status': 'PASS', 'calculation': candidate}]
        self.numeric = None
        self.numeric_path = None

    def write_source(self):
        buffer = io.StringIO(newline='')
        writer = csv.writer(buffer, lineterminator='\r\n')
        writer.writerow(self.header); writer.writerows(self.values)
        self.source.write_bytes(buffer.getvalue().encode())
        self.artifacts = [{'artifact_id': 'A1', 'path': 'original.csv', 'sha256': digest(self.source)}]
        self.rows = [{'row_id': f'R{i}', 'artifact_id': 'A1', 'artifact_sha256': digest(self.source),
            'logical_record': i+1, 'header': self.header[:], 'values': values[:],
            'identity_binding': {'identity_status': 'PASS', 'bound_study_key': 'human_rna:GSE1',
                'candidate_study_keys': ['human_rna:GSE1'], 'identity_issues': []},
            'fields': dict(zip(self.header, values)) if len(values) == len(self.header) and len(set(self.header)) == len(self.header) else None}
            for i, values in enumerate(self.values, 1)]
        self.decisions = [{'row_id': row['row_id'], 'status': 'PASS', 'action': 'retain',
            'source_comparison_ids': ['paired'], 'reason': 'Software fixture model decision',
            'evidence_ids': ['S0'], 'unresolved': [], 'proposed_repair': None} for row in self.rows]
        self.checks = [{'review_stage': 'fixture', 'row_ids': [r['row_id'] for r in self.rows],
            'checks': {name: {'status': 'CANNOT VERIFY', 'reason': 'Software fixture only',
                'evidence_ids': ['S0'], 'unresolved': ['No scientific reference'], 'scope': 'fixture'} for name in SCIENTIFIC_CHECKS}}]

    def run_numeric(self):
        rscript = shutil.which('Rscript')
        if not rscript:
            self.skipTest('Independent base R execution is unavailable')
        self.numeric = verify_calculations(self.calculations, self.predictions, self.root / 'r-verification', rscript=rscript)
        self.numeric_path = self.root / 'r-verification/receipt.json'
        self.assertEqual(self.numeric['counts']['PASS'], 1)

    def repair(self, name='repairs', **kwargs):
        return build_repair_candidates(self.rows, self.design, self.calculations, self.numeric,
            self.checks, self.decisions, self.root / name, artifact_sources=self.artifacts,
            source_root=self.root, numerical_receipt_path=self.numeric_path, **kwargs)

    def output(self, name):
        return json.loads((self.root / 'repairs' / name).read_text())

    def test_literal_defect_withholds_versioned_record_and_preserves_other_bytes(self):
        original = self.source.read_bytes()
        receipt = self.repair()
        self.assertEqual(receipt['status'], 'PASS')
        self.assertEqual(receipt['withheld_rows'], 1)
        self.assertEqual(receipt['numeric_replacements'], 0)
        self.assertEqual(receipt['accepted_scientific_repairs'], 0)
        self.assertEqual(self.source.read_bytes(), original)
        manifest = self.output('withholding-manifest.json')
        file = manifest['files'][0]
        self.assertEqual((self.root/'repairs'/file['preserved_original']).read_bytes(), original)
        with (self.root/'repairs'/file['withheld_copy']).open(newline='') as handle:
            reader = csv.DictReader(handle); after = list(reader)
        self.assertEqual(after[0]['p_value'], '')
        self.assertEqual(after[0]['fdr'], '')
        self.assertEqual(after[0]['neg_log10_pvalue'], '')
        self.assertEqual(after[0]['is_significant'], '')
        self.assertEqual(after[0]['notes'], self.rows[0]['fields']['notes'])
        self.assertEqual(after[1:], [r['fields'] for r in self.rows[1:]])
        self.assertEqual(verify_corpus_repairs(self.root/'repairs')['status'], 'PASS')

    def test_missing_nan_undefined_markers_do_not_become_material_defects(self):
        self.values = [['GSE1','clock','Paired t-test',mark,mark,'','','','','','missing']
                       for mark in ('', 'NA', 'NaN', 'undefined', 'None')]
        self.write_source(); self.artifacts[0]['numeric_fields'] = ['effect_size']
        receipt = self.repair()
        self.assertEqual(receipt['withheld_rows'], 0)
        self.assertEqual(receipt['withheld_files'], 0)
        for row in self.output('row-decisions.json'):
            self.assertTrue(all(c['status'] == 'CANNOT VERIFY' for c in row['literal_checks']))

    def test_effect_field_requires_explicit_numeric_contract(self):
        self.assertEqual(self.repair('untyped')['withheld_rows'], 1)
        self.artifacts[0]['numeric_fields'] = ['effect_size']
        receipt = self.repair()
        self.assertEqual(receipt['withheld_rows'], 2)
        checks = self.output('row-decisions.json')[2]['literal_checks']
        self.assertEqual(next(c['status'] for c in checks if c['field'] == 'effect_size'), 'FAIL')

    def test_exact_r_verified_candidate_retains_historical_and_scientific_uncertainty(self):
        self.run_numeric()
        receipt = self.repair()
        self.assertEqual(receipt['conditional_numeric_candidates'], 3)
        for row in self.output('correction-candidates.json'):
            self.assertEqual(row['status'], 'PASS')
            self.assertEqual(row['historical_sample_lineage'], 'CANNOT VERIFY')
            self.assertEqual(row['original_orientation'], 'CANNOT VERIFY')
            self.assertFalse(row['accepted_scientific_repair'])
            self.assertEqual(row['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(receipt['accepted_scientific_repairs'], 0)

    def test_clock_alias_and_multiple_proposed_comparisons_are_not_silently_resolved(self):
        self.run_numeric()
        self.values[1][1] = 'clock alias'; self.write_source()
        self.decisions[2]['source_comparison_ids'] = ['paired', 'other']
        receipt = self.repair()
        self.assertEqual(receipt['conditional_numeric_candidates'], 1)
        candidates = self.output('correction-candidates.json')
        self.assertIn('literal clock', candidates[1]['reason'])
        self.assertIn('Exactly one', candidates[2]['reason'])

    def test_other_study_in_mixed_source_file_cannot_receive_this_study_candidate(self):
        self.run_numeric()
        self.values[1][0] = 'GSE2'; self.write_source()
        receipt = self.repair()
        self.assertEqual(receipt['conditional_numeric_candidates'], 2)
        self.assertIn('source-study identity', self.output('correction-candidates.json')[1]['reason'])

    def test_bound_future_producer_without_gse_column_is_supported_but_unbound_row_stops(self):
        self.run_numeric()
        self.header.pop(0)
        self.values = [values[1:] for values in self.values]
        self.write_source()
        for row in self.rows:
            row['identity_binding']['producer_stage_sha256'] = 'a'*64
        self.rows[1]['identity_binding']['identity_status'] = 'CANNOT VERIFY'
        self.rows[2].pop('identity_binding')
        receipt = self.repair()
        self.assertEqual(receipt['conditional_numeric_candidates'], 1)
        self.assertEqual(self.output('correction-candidates.json')[0]['status'], 'PASS')

    def test_groups_design_freeze_and_independent_request_must_match(self):
        self.run_numeric()
        self.calculations[0]['case_ids'] = self.case[::-1]
        receipt = self.repair()
        self.assertEqual(receipt['conditional_numeric_candidates'], 0)
        self.assertEqual(receipt['numerical_binding']['status'], 'FAIL')
        self.assertIn('exact source groups', self.output('correction-candidates.json')[0]['reason'])

    def test_in_memory_pass_without_file_bound_r_receipt_is_not_accepted(self):
        self.numeric = {'schema': 'clockbase-independent-statistics-v1', 'status': 'PASS',
                        'results': [{'input_index': 0, 'status': 'PASS', 'numerical_status': 'PASS'}]}
        receipt = self.repair()
        self.assertEqual(receipt['conditional_numeric_candidates'], 0)
        self.assertEqual(receipt['numerical_binding']['status'], 'CANNOT VERIFY')

    def test_tampered_r_output_blocks_candidate_but_not_independent_literal_withholding(self):
        self.run_numeric()
        (self.root/'r-verification/comparison-000000-reference.csv').write_text('tampered')
        receipt = self.repair()
        self.assertEqual(receipt['conditional_numeric_candidates'], 0)
        self.assertEqual(receipt['numerical_binding']['status'], 'FAIL')
        self.assertEqual(receipt['withheld_rows'], 1)

    def test_all_agent_checks_pass_still_cannot_accept_a_scientific_repair(self):
        self.run_numeric()
        for check in self.checks[0]['checks'].values():
            check.update(status='PASS', unresolved=[])
        receipt = self.repair()
        self.assertEqual(receipt['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(receipt['accepted_scientific_repairs'], 0)
        self.assertEqual(receipt['withheld_rows'], 1)  # Mechanical defect overrides model retain.

    def test_scientific_only_withholding_is_retained_as_a_proposal(self):
        self.decisions[1].update(status='FAIL', action='withhold', reason='Model alleges an unsupported biological claim')
        receipt = self.repair()
        row = self.output('row-decisions.json')[1]
        self.assertEqual(row['effective_action'], 'WITHHOLDING_PROPOSED')
        self.assertEqual(receipt['withheld_rows'], 1)

    def test_all_sixteen_checks_and_exact_decision_coverage_are_required(self):
        original = copy.deepcopy(self.checks)
        self.checks[0]['checks'].pop(next(iter(SCIENTIFIC_CHECKS)))
        with self.assertRaisesRegex(ValueError, '16 scientific checks'):
            self.repair()
        self.checks = original
        self.checks[0]['row_ids'].pop()
        with self.assertRaisesRegex(ValueError, 'omit, duplicate or invent'):
            self.repair()
        self.checks = original; self.decisions.pop()
        with self.assertRaisesRegex(ValueError, 'exactly one'):
            self.repair()

    def test_missing_row_changed_source_and_row_content_tampering_stop(self):
        original = self.source.read_bytes()
        self.source.write_bytes(original + b'changed')
        with self.assertRaisesRegex(ValueError, 'Original source changed'):
            self.repair()
        self.source.write_bytes(original)
        self.rows[0]['fields']['p_value'] = '0.1'
        with self.assertRaisesRegex(ValueError, 'row content differs'):
            self.repair()

    def test_ragged_source_is_preserved_without_automatic_rewrite(self):
        self.values[-1] = self.values[-1][:-2]; self.write_source()
        receipt = self.repair()
        self.assertEqual(receipt['withheld_files'], 0)
        self.assertEqual(receipt['withheld_rows'], 0)
        self.assertIn('Malformed or ambiguous', self.output('row-decisions.json')[0]['withholding_limitation'])
        self.assertEqual(self.output('row-decisions.json')[0]['literal_checks'][1]['status'], 'FAIL')

    def test_probability_domain_checks_do_not_round_tiny_negative_or_above_one_into_validity(self):
        for index, value in enumerate(('-1e-400', '1.00000000000000000000000001', 'Infinity')):
            self.values[0][3] = value; self.write_source()
            with self.subTest(value=value):
                receipt = self.repair(f'out-{index}')
                self.assertEqual(receipt['withheld_rows'], 1)
                self.assertEqual(receipt['status'], 'PASS')

    def test_resume_checks_all_output_hashes_and_independent_readback_detects_new_changes(self):
        receipt = self.repair()
        self.assertEqual(self.repair(), receipt)
        manifest = self.output('withholding-manifest.json')
        changed = self.root/'repairs'/manifest['files'][0]['withheld_copy']
        changed.write_bytes(changed.read_bytes().replace(b'unchanged literal source', b'changed literal source'))
        self.assertEqual(verify_corpus_repairs(self.root/'repairs')['status'], 'FAIL')
        with self.assertRaisesRegex(ValueError, 'output changed on resume'):
            self.repair()

    def test_partial_output_and_source_hash_changes_are_never_overwritten(self):
        output = self.root/'repairs'; output.mkdir()
        (output/'inputs.json').write_text('preserved partial artifact')
        with self.assertRaisesRegex(FileExistsError, 'preserved'):
            self.repair()
        self.assertEqual((output/'inputs.json').read_text(), 'preserved partial artifact')

    def test_independent_readback_rejects_a_faulty_writer_without_claiming_withholding(self):
        # Deliberate software fault injection: the writer fails to apply fields.
        # Independent source/copy readback must catch it and preserve the attempt.
        with patch('audit_system.corpus_repairs_2026_09_20.patch_csv', side_effect=lambda raw, spec: raw):
            receipt = self.repair()
        self.assertEqual(receipt['status'], 'FAIL')
        self.assertEqual(receipt['attempted_withholding_rows'], 1)
        self.assertEqual(receipt['withheld_rows'], 0)
        self.assertEqual(receipt['accepted_scientific_repairs'], 0)
        self.assertEqual(self.output('row-decisions.json')[0]['effective_action'], 'WITHHOLDING_UNVERIFIED')
        self.assertTrue((self.root/'repairs/independent-readback.json').is_file())


if __name__ == '__main__':
    unittest.main()
