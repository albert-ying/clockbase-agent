"""Constructed software fixtures; these are not empirical audit-accuracy evidence.

The producer is imported here solely as the system under test. The independent
verifier imports no producer or numerical Python library and runs base R against
the prediction CSV. Fixture values below are disclosed, deterministic inputs.
"""
import csv
import copy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from audit_system.independent_statistics_2026_09_20 import _csv, _prepared, model_verification_view, verify_calculations
from audit_system.statistical_rerun_2026_09_20 import recompute_contrast


class IndependentStatisticsTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix='independent-statistics-test-')
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)
        self.source = self.root / 'predictions.csv'
        self.ids = [f'GSM{i}' for i in range(1, 7)]
        self.frame = pd.DataFrame({'clock': [4., 8., 7., 1., 4., 5.]}, index=self.ids)
        self.case, self.control = self.ids[:3], self.ids[3:]
        # Leading zero unit labels must survive the Python-to-R boundary.
        self.pairs = {'pairs': [dict(case_id=a, control_id=b, unit_id=unit)
            for unit, a, b in zip(['01', '1', '02'], self.case, self.control)]}
        self.write_source()

    def write_source(self, header='sample_id'):
        self.frame.to_csv(self.source, index_label=header)

    def record(self, method='Welch t-test', design=None, case=None):
        selected = self.case if case is None else case
        candidate = recompute_contrast(self.frame, selected, self.control, 'clock', method, design=design)
        return {'comparison_id': 'constructed-software-fixture', 'clock': 'clock',
                'case_ids': selected, 'control_ids': self.control, 'status': 'PASS',
                'scientific_status': 'CANNOT VERIFY', 'calculation': candidate}

    def verify(self, records, directory='out', **kwargs):
        if 'rscript' not in kwargs:
            executable = shutil.which('Rscript')
            if not executable:
                self.skipTest('Rscript unavailable; independent numerical execution cannot be tested')
            kwargs['rscript'] = executable
        return verify_calculations(records, self.source, self.root / directory, **kwargs)

    def assert_pass(self, receipt, count=1):
        self.assertEqual(receipt['status'], 'PASS', receipt['results'])
        self.assertEqual(receipt['counts'], {'PASS': count, 'FAIL': 0, 'CANNOT VERIFY': 0})
        self.assertEqual(receipt['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(receipt['accepted_scientific_repairs'], 0)
        self.assertEqual(receipt['model_calls'], 0)
        self.assertEqual(receipt['clock_runs'], 0)
        self.assertTrue(receipt['r_process']['executed'])

    def test_paired_explicit_map_and_literal_unit_ids_survive_arm_reordering(self):
        receipt = self.verify([self.record('Paired t-test', self.pairs, self.case[::-1])])
        self.assert_pass(receipt)
        reference = receipt['results'][0]['reference_values']
        self.assertEqual(dict(reference['within_unit_difference']), {'01': 3., '1': 4., '02': 2.})
        self.assertEqual(reference['n_biological_units'][0][1], 3)
        self.assertEqual(reference['effect_size'][0][1], 3.)

    def test_welch_and_student_match_independent_r_with_unequal_group_sizes(self):
        self.frame = pd.DataFrame({'clock': [4., 8., 20., 1., 4., 5., 2.]},
                                  index=[f'GSM{i}' for i in range(1, 8)])
        self.control = list(self.frame.index[3:])
        self.write_source()
        receipt = self.verify([self.record('Welch t-test'), self.record('Student t-test')])
        self.assert_pass(receipt, 2)
        a, b = receipt['results']
        self.assertNotEqual(a['reference_values']['p_value'], b['reference_values']['p_value'])
        self.assertNotEqual(a['reference_values']['degrees_of_freedom'], b['reference_values']['degrees_of_freedom'])

    def test_zero_mean_difference_has_finite_independently_computed_standard_error(self):
        self.frame['clock'] = [0., 1., 2., -1., 1., 3.]
        self.write_source()
        receipt = self.verify([self.record()])
        self.assert_pass(receipt)
        reference = receipt['results'][0]['reference_values']
        self.assertEqual(reference['p_value'][0][1], 1.)
        self.assertGreater(reference['standard_error'][0][1], 0.)

    def descriptive(self):
        return {'comparison_id': 'constructed-assay-fixture', 'clock': 'clock',
                'case_ids': self.case, 'control_ids': self.control, 'status': 'PASS',
                'calculation': {'case_mean': 19/3, 'control_mean': 10/3,
                    'mean_difference': 3., 'n_case_assays': 3, 'n_control_assays': 3,
                    'p_value': None, 'effect_size': None,
                    'scope': 'Constructed technical assays; no population inference'}}

    def test_descriptive_assays_verify_means_but_have_no_inferential_result(self):
        receipt = self.verify([self.descriptive()])
        self.assert_pass(receipt)
        self.assertNotIn('p_value', receipt['results'][0]['reference_values'])
        self.assertNotIn('ci_low', receipt['results'][0]['reference_values'])

    def test_descriptive_output_cannot_smuggle_in_an_inferential_p_value(self):
        record = self.descriptive()
        record['calculation']['p_value'] = .001
        receipt = self.verify([record])
        self.assertEqual(receipt['counts']['FAIL'], 1)
        self.assertFalse(receipt['r_process']['executed'])
        self.assertIn('inferential', receipt['results'][0]['reason'])

    def ols_design(self):
        self.ids = [f'GSM{i}' for i in range(1, 9)]
        self.case, self.control = self.ids[:4], self.ids[4:]
        # y=10+2*case+3*baseline+orthogonal disclosed residuals.
        self.frame = pd.DataFrame({'clock': [13., 14., 17., 22., 8., 15., 18., 17.]}, index=self.ids)
        self.write_source()
        return {'covariates': [{'name': 'baseline', 'values': dict(zip(self.ids, [0., 1., 2., 3., 0., 1., 2., 3.]))}]}

    def test_classical_and_hc3_ols_verify_coefficients_covariance_ci_and_residual_df(self):
        design = self.ols_design()
        receipt = self.verify([self.record('OLS adjusted contrast', design),
                               self.record('OLS HC3 adjusted contrast', design)])
        self.assert_pass(receipt, 2)
        a, b = [row['reference_values'] for row in receipt['results']]
        self.assertAlmostEqual(a['adjusted_mean_difference'][0][1], 2.)
        self.assertAlmostEqual(a['coefficients'][2][1], 3.)
        self.assertEqual(a['degrees_of_freedom'][0][1], 5)
        self.assertNotEqual(a['standard_error'], b['standard_error'])

    def test_all_t_interval_methods_match_r_and_reject_wrong_intervals(self):
        covariates = self.ols_design()
        pairs = {'pairs': [dict(unit_id=str(i), case_id=a, control_id=b)
                          for i, (a, b) in enumerate(zip(self.case, self.control))]}
        records = []
        for method, design in (('Paired t-test', pairs), ('Welch t-test', None),
                ('Student t-test', None), ('OLS adjusted contrast', covariates),
                ('OLS HC3 adjusted contrast', covariates)):
            for confidence in (.8, .95, .99):
                record = self.record(method, design)
                record['comparison_id'] = method + ':' + str(confidence)
                record['calculation'] = recompute_contrast(self.frame, self.case,
                    self.control, 'clock', method, design=design, confidence_level=confidence)
                records.append(record)
        self.assert_pass(self.verify(records, 'all-methods'), 15)
        mutations = copy.deepcopy(records)
        for record in mutations:
            record['calculation']['confidence_interval'][0] += .001
            record['calculation']['confidence_interval'][1] += .001
        receipt = self.verify(mutations, 'wrong-intervals')
        self.assertEqual(receipt['counts'], {'PASS': 0, 'FAIL': 15, 'CANNOT VERIFY': 0})
        for result in receipt['results']:
            self.assertEqual([check['name'] for check in result['checks'] if check['status'] == 'FAIL'],
                             ['confidence_interval'])

    def test_p_ci_effect_definition_and_internal_source_binding_mutations_fail(self):
        changes = [lambda c: c.update(p_value=.123456789),
                   lambda c: c.update(confidence_interval=[-100, 100]),
                   lambda c: c.update(effect_size_definition='different denominator'),
                   lambda c: c.update(clock='other clock'),
                   lambda c: c.update(case_ids=self.case[::-1])]
        for index, change in enumerate(changes):
            record = self.record(); change(record['calculation'])
            with self.subTest(index=index):
                receipt = self.verify([record], f'out-{index}')
                self.assertEqual(receipt['counts']['FAIL'], 1)

    def test_tiny_p_values_use_relative_comparison_without_absolute_tolerance_loophole(self):
        # This disclosed software fixture produces a small computed P value;
        # the mutation is 100-fold and remains below the ordinary absolute tol.
        self.ids = [f'GSM{i}' for i in range(1, 41)]
        self.case, self.control = self.ids[:20], self.ids[20:]
        self.frame = pd.DataFrame({'clock': [100.+i/10 for i in range(20)] + [i/10 for i in range(20)]}, index=self.ids)
        self.write_source()
        record = self.record()
        self.assertLess(record['calculation']['p_value'] * 100, 1e-10)
        self.assert_pass(self.verify([record], 'good'))
        record['calculation']['p_value'] *= 100
        receipt = self.verify([record], 'mutated')
        self.assertEqual(receipt['counts']['FAIL'], 1)
        self.assertIn('p_value', [c['name'] for c in receipt['results'][0]['checks'] if c['status'] == 'FAIL'])

    def test_pair_map_and_covariate_matrix_mutations_are_detected(self):
        paired = self.record('Paired t-test', self.pairs)
        pairs = paired['calculation']['pairs']
        pairs[0]['control_id'], pairs[1]['control_id'] = pairs[1]['control_id'], pairs[0]['control_id']
        self.assertEqual(self.verify([paired], 'pairs')['counts']['FAIL'], 1)
        design = self.ols_design()
        adjusted = self.record('OLS adjusted contrast', design)
        adjusted['calculation']['design_matrix'][0][2] += 1
        receipt = self.verify([adjusted], 'covariates')
        self.assertEqual(receipt['counts']['FAIL'], 1)
        self.assertIn('design_matrix binds declared sample covariates',
            [c['name'] for c in receipt['results'][0]['checks'] if c['status'] == 'FAIL'])

    def test_invalid_pair_shape_and_incomplete_covariates_fail_contract_without_crash(self):
        paired = self.record('Paired t-test', self.pairs)
        paired['calculation']['pairs'][0] = 'not a pairing object'
        self.assertEqual(self.verify([paired], 'pairs')['counts']['FAIL'], 1)
        design = self.ols_design()
        adjusted = self.record('OLS adjusted contrast', design)
        adjusted['calculation']['covariates'][0]['values'].pop(self.ids[0])
        self.assertEqual(self.verify([adjusted], 'covariates')['counts']['FAIL'], 1)

    def test_ambiguous_clock_does_not_get_an_automatic_alias(self):
        record = self.record()
        self.source.write_text('sample_id,clock,clock\n' + ''.join(f'{sid},{i},{i+1}\n' for i, sid in enumerate(self.ids)))
        record.update(status='CANNOT VERIFY'); record.pop('calculation')
        receipt = self.verify([record])
        self.assertEqual(receipt['counts']['CANNOT VERIFY'], 1)
        self.assertEqual(receipt['verified_input_stops'], 1)
        self.assertIn('ambiguous', receipt['results'][0]['reason'])

    def test_duplicate_and_conflicting_source_identities_are_unverifiable(self):
        record = self.record()
        sources = [self.source.read_text() + 'GSM1,10\n',
                   'sample_id,gsm_id,clock\nGSM1,GSM2,1\n',
                   'sample_id,sample_id,clock\nGSM1,GSM2,1\n']
        for index, source in enumerate(sources):
            self.source.write_text(source)
            with self.subTest(index=index):
                receipt = self.verify([record], f'out-{index}')
                self.assertEqual(receipt['counts']['CANNOT VERIFY'], 1)
                self.assertFalse(receipt['r_process']['executed'])

    def test_repeated_identity_schema_finding_is_retained_with_numerical_pass(self):
        record = self.record()
        rows = [[sid, str(value), sid, sid] for sid, value in zip(self.ids, self.frame['clock'])]
        with self.source.open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['sample_id', 'clock', 'sample_id', 'gsm_id'])
            writer.writerows(rows)
        original_bytes = self.source.read_bytes()
        receipt = self.verify([record])
        self.assert_pass(receipt)
        identity = json.loads((self.root / 'out/source-identity.json').read_text())
        self.assertEqual(receipt['source_identity'], identity)
        self.assertEqual(identity['columns'], [0, 2, 3])
        finding = identity['schema_findings'][0]
        self.assertEqual(finding['header'], 'sample_id')
        self.assertEqual(finding['positions'], [0, 2])
        self.assertEqual(finding['provenance_status'], 'UNRESOLVED')
        self.assertIn('every source row', finding['scope'])
        self.assertEqual(self.source.read_bytes(), original_bytes)
        expected_hash = hashlib.sha256(original_bytes).hexdigest()
        self.assertEqual(receipt['predictions_sha256_before'], expected_hash)
        self.assertEqual(receipt['predictions_sha256_after'], expected_hash)
        self.assertEqual(receipt, self.verify([record]))

    def test_blank_index_header_requires_full_column_canonical_gsm_binding(self):
        record = self.record()
        self.write_source(header='')
        receipt = self.verify([record], 'canonical')
        self.assert_pass(receipt)
        identity = json.loads((self.root / 'canonical/source-identity.json').read_text())
        self.assertEqual(identity['blank_header_binding']['format'], 'GSM accession')
        self.source.write_text(self.source.read_text().replace('GSM1,', '0,'))
        invalid = self.verify([record], 'row-numbers')
        self.assertEqual(invalid['counts']['CANNOT VERIFY'], 1)
        self.assertIn('canonical', invalid['results'][0]['reason'])

    def test_source_missing_nonfinite_ragged_empty_and_lfs_cases_retain_stops(self):
        record = self.record()
        sources = [self.source.read_text().replace('GSM6,5.0\n', ''),
                   self.source.read_text().replace('GSM1,4.0', 'GSM1,nan'),
                   'sample_id,clock\nGSM1,1,extra\n', '\n', '',
                   'version https://git-lfs.github.com/spec/v1\noid sha256:abc\nsize 123\n']
        for index, source in enumerate(sources):
            self.source.write_text(source)
            with self.subTest(index=index):
                receipt = self.verify([record], f'out-{index}')
                self.assertEqual(receipt['counts']['CANNOT VERIFY'], 1)
                self.assertFalse(receipt['r_process']['executed'])
        self.source.unlink()
        missing = self.verify([record], 'missing-file')
        self.assertIn('Cannot proceed - required data file not found:', missing['results'][0]['reason'])

    def test_unresolved_design_is_retained_and_not_treated_as_independently_adjudicated(self):
        record = self.record()
        record.pop('calculation'); record.update(status='CANNOT VERIFY', reason='Unknown biological units')
        receipt = self.verify([record])
        self.assertEqual(receipt['counts']['CANNOT VERIFY'], 1)
        self.assertEqual(receipt['verified_input_stops'], 0)
        self.assertFalse(receipt['r_process']['executed'])
        self.assertIn('not adjudicated', receipt['results'][0]['reason'])

    def test_missing_r_and_timeout_never_fabricate_a_numerical_pass(self):
        record = self.record()
        missing = self.verify([record], 'missing-r', rscript=str(self.root / 'no-Rscript'))
        self.assertEqual(missing['counts']['CANNOT VERIFY'], 1)
        self.assertIn('unavailable', missing['r_process']['error'])
        with patch('audit_system.independent_statistics_2026_09_20.subprocess.run',
                   side_effect=subprocess.TimeoutExpired('disclosed process-control fixture', .1, output=b'partial stdout')):
            timed = self.verify([record], 'timeout', timeout_seconds=.1)
        self.assertEqual(timed['counts']['CANNOT VERIFY'], 1)
        self.assertTrue(timed['r_process']['timed_out'])
        self.assertEqual(timed['r_process']['stdout'], 'partial stdout')
        self.assertTrue((self.root / 'timeout/request.json').is_file())

    def test_degenerate_or_confounded_reference_design_is_explicitly_unverifiable(self):
        record = self.record()
        self.frame['clock'] = [1., 1., 1., 2., 2., 2.]
        self.write_source()
        receipt = self.verify([record], 'degenerate')
        self.assertEqual(receipt['counts']['CANNOT VERIFY'], 1)
        self.assertIn('Independent R stopped', receipt['results'][0]['reason'])
        design = self.ols_design()
        record = self.record('OLS adjusted contrast', design)
        record['calculation']['covariates'][0]['values'] = dict(zip(self.ids, [1, 1, 1, 1, 0, 0, 0, 0]))
        receipt = self.verify([record], 'confounded')
        self.assertEqual(receipt['counts']['CANNOT VERIFY'], 1)
        self.assertIn('rank-deficient', receipt['results'][0]['reason'])

    def test_resume_binds_source_and_every_generated_output_hash(self):
        record = self.record()
        receipt = self.verify([record])
        self.assert_pass(receipt)
        self.assertEqual(receipt, self.verify([record]))
        for relative, expected in receipt['outputs_sha256'].items():
            self.assertEqual(hashlib.sha256((self.root / 'out' / relative).read_bytes()).hexdigest(), expected)
        (self.root / 'out/comparison-000000-reference.csv').write_text('tampered')
        with self.assertRaisesRegex(ValueError, 'output changed on resume'):
            self.verify([record])
        self.assertTrue((self.root / 'out/receipt.json').is_file())

    def test_changed_source_and_partial_previous_output_require_new_directory(self):
        record = self.record()
        self.verify([record])
        self.source.write_text(self.source.read_text().replace('GSM1,4.0', 'GSM1,4.1'))
        with self.assertRaisesRegex(ValueError, 'inputs or implementation changed'):
            self.verify([record])
        partial = self.root / 'partial'; partial.mkdir()
        (partial / 'request.json').write_text('preserved partial output')
        with self.assertRaisesRegex(FileExistsError, 'preserved'):
            self.verify([record], 'partial')
        self.assertEqual((partial / 'request.json').read_text(), 'preserved partial output')

    def test_receipt_mutation_is_detected_and_empty_input_cannot_pass(self):
        record = self.record(); self.verify([record])
        (self.root / 'out/receipt.json').write_text('{}')
        with self.assertRaisesRegex(ValueError, 'receipt changed'):
            self.verify([record])
        empty = self.verify([], 'empty')
        self.assertEqual(empty['status'], 'CANNOT VERIFY')
        self.assertFalse(empty['r_process']['executed'])

    def test_model_view_retains_safe_numeric_evidence_and_excludes_local_paths_and_errors(self):
        receipt = self.verify([self.record(), self.descriptive()])
        receipt['results'][0].update(reason='/Entropy/private/path: raw process error',
                                    scope='arbitrary local diagnostic /home/private')
        receipt['r_process']['command'] = ['/private/command']
        receipt['results'][0]['checks'].append({'name': '/private/untrusted-check', 'status': 'FAIL',
                                              'observed': '/private/candidate', 'expected': '/private/reference'})
        view = model_verification_view(receipt)
        serialized = json.dumps(view)
        for fragment in ('/Entropy', '/home', '/private', 'raw process error', 'command', 'outputs_sha256'):
            self.assertNotIn(fragment, serialized)
        self.assertEqual(view['counts']['PASS'], 2)
        self.assertEqual(view['results'][0]['check_counts']['FAIL'], 1)
        self.assertEqual(view['results'][0]['unlisted_check_count'], 1)
        self.assertEqual(view['results'][0]['clock'], 'clock')
        self.assertEqual(view['results'][0]['reference_metrics']['p_value'],
                         receipt['results'][0]['reference_values']['p_value'][0][1])
        self.assertEqual(view['results'][1]['reference_metrics']['n_case_assays'], 3.)
        self.assertNotIn('n_case', view['results'][1]['reference_metrics'])
        self.assertEqual(view['scientific_status'], 'CANNOT VERIFY')

    def test_model_view_hashes_unsafe_labels_and_omits_nonnumeric_metric_values(self):
        receipt = self.verify([self.record()])
        row = receipt['results'][0]
        row['clock'] = '/home/owner/unpublished.csv'
        row['comparison_id'] = 'arbitrary\nfree text /Entropy/private'
        row['method'] = 'unexpected /private/model'
        row['reference_values']['p_value'] = [['index /private', '/private/value']]
        view = model_verification_view(receipt)
        result = view['results'][0]
        self.assertIsNone(result['clock']); self.assertIsNone(result['comparison_id'])
        self.assertEqual(len(result['clock_sha256']), 64)
        self.assertEqual(len(result['comparison_id_sha256']), 64)
        self.assertIsNone(result['method'])
        self.assertNotIn('p_value', result['reference_metrics'])
        self.assertNotIn('/private', json.dumps(view))
        self.assertNotIn('/home', json.dumps(view))
        self.assertEqual(view['scientific_status'], 'CANNOT VERIFY')


class IndependentIdentityReaderTests(unittest.TestCase):
    """Deterministic CSV contract fixtures; no R or scientific inference."""

    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory(prefix='independent-identity-test-')
        self.addCleanup(self.scratch.cleanup)
        self.source = Path(self.scratch.name) / 'predictions.csv'

    def write(self, header, rows):
        with self.source.open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(header)
            writer.writerows(rows)

    def test_all_allowlisted_repeated_labels_preserve_every_header_and_cell(self):
        for label in ('sample_id', 'gsm_id', 'GSM', 'geo_accession'):
            with self.subTest(label=label):
                header = [label, 'clock', 'annotation', label, label]
                rows = [['GSM1', '01.2300', 'literal, quoted cell', 'GSM1', 'GSM1'],
                        ['GSM2', '-0.0000', 'unchanged text', 'GSM2', 'GSM2']]
                self.write(header, rows)
                before = self.source.read_bytes()
                parsed_header, records, binding = _csv(self.source)
                self.assertEqual(parsed_header, header)
                self.assertEqual(list(records), ['GSM1', 'GSM2'])
                self.assertEqual(list(records.values()), rows)
                self.assertEqual(binding['columns'], [0, 3, 4])
                self.assertEqual(binding['rows'], 2)
                self.assertEqual(binding['schema_findings'][0]['header'], label)
                self.assertEqual(binding['schema_findings'][0]['positions'], [0, 3, 4])
                self.assertEqual(binding['schema_findings'][0]['provenance_status'], 'UNRESOLVED')
                self.assertEqual(self.source.read_bytes(), before)

    def test_mixed_aliases_and_blank_binding_check_conflicts_outside_selected_arms(self):
        header = ['', 'sample_id', 'clock', 'gsm_id', 'sample_id', 'GSM', 'geo_accession']
        rows = [[f'GSM{i}', f'GSM{i}', str(i), f'GSM{i}', f'GSM{i}', f'GSM{i}', f'GSM{i}']
                for i in range(1, 8)]
        self.write(header, rows)
        parsed_header, records, binding = _csv(self.source)
        selection = {'case_ids': ['GSM1', 'GSM2', 'GSM3'],
                     'control_ids': ['GSM4', 'GSM5', 'GSM6'], 'clock': 'clock'}
        self.assertEqual(len(_prepared(selection, parsed_header, records)[2]), 6)
        self.assertEqual(binding['blank_header_binding']['index'], 0)
        self.assertEqual(binding['schema_findings'][0]['positions'], [1, 4])
        # The seventh source row is not used by the contrast. Every identity
        # position, including the blank binding and each alias, must agree.
        for position in binding['columns']:
            changed = [list(row) for row in rows]
            changed[-1][position] = 'GSM999'
            self.write(header, changed)
            with self.subTest(position=position), self.assertRaisesRegex(ValueError, 'conflict'):
                _csv(self.source)

    def test_repeated_labels_do_not_relax_malformed_duplicate_or_ragged_stops(self):
        header = ['sample_id', 'clock', 'sample_id']
        cases = [([['GSM1', '1', 'GSM1'], ['GSM1', '2', 'GSM1']], 'duplicated'),
                 ([['GSM1', '1', 'GSM1'], ['GSM２', '2', 'GSM２']], 'canonical'),
                 ([['GSM1', '1', 'GSM1'], [' GSM2', '2', ' GSM2']], 'canonical'),
                 ([['GSM1', '1', 'GSM1'], ['', '2', '']], 'canonical'),
                 ([['GSM1', '1', 'GSM1'], ['GSM2', '2']], 'ragged')]
        for rows, reason in cases:
            self.write(header, rows)
            with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                _csv(self.source)
        self.write(['sampleId', 'clock', 'sampleId'], [['GSM1', '1', 'GSM1']])
        with self.assertRaisesRegex(ValueError, 'absent'):
            _csv(self.source)

    def test_repeated_identity_does_not_relax_selected_clock_or_nonfinite_stops(self):
        selection = {'case_ids': ['GSM1'], 'control_ids': ['GSM2'], 'clock': 'clock'}
        cases = [(['sample_id', 'clock', 'sample_id', 'clock'],
                  [['GSM1', '1', 'GSM1', '1'], ['GSM2', '2', 'GSM2', '2']], 'ambiguous'),
                 (['sample_id', 'clock', 'sample_id'],
                  [['GSM1', 'nan', 'GSM1'], ['GSM2', '2', 'GSM2']], 'nonfinite'),
                 (['sample_id', 'clock', 'sample_id'],
                  [['GSM1', 'inf', 'GSM1'], ['GSM2', '2', 'GSM2']], 'nonfinite')]
        for header, rows, reason in cases:
            self.write(header, rows)
            parsed_header, records, _ = _csv(self.source)
            with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                _prepared(selection, parsed_header, records)


if __name__ == '__main__':
    unittest.main()
