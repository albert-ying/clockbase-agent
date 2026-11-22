"""Constructed software fixtures, not empirical audit-accuracy evidence."""
import copy
import json
import math
import sys
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import stdtr

from audit_system.statistical_rerun_2026_09_20 import (
    _critical_t, recompute_contrast, bind_design, select_groups, clock_frame, run, sha, SCHEMA, POOLED, record_hash)
from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError


class StatisticalRerunTests(unittest.TestCase):
    def setUp(self):
        self.frame = pd.DataFrame({'clock': [4., 8., 7., 1., 4., 5.]}, index=list('ABCDEF'))
        self.case, self.control = list('ABC'), list('DEF')
        self.pairs = {'pairs': [dict(case_id=a, control_id=b, unit_id=str(i))
                               for i, (a, b) in enumerate(zip(self.case, self.control))]}

    def test_critical_t_matches_analytic_df_one_and_two(self):
        for confidence in (.5, .8, .9, .95, .99, .999, 1-1e-6, 1-1e-12):
            probability = (1+confidence)/2
            tail = 1-probability
            expected = {1: 1/math.tan(math.pi*tail),
                        2: (2*probability-1)/math.sqrt(2*probability*tail)}
            for df, reference in expected.items():
                with self.subTest(df=df, confidence=confidence):
                    self.assertTrue(math.isclose(_critical_t(df, confidence), reference,
                                                 rel_tol=1e-12, abs_tol=5e-14))

    def test_critical_t_invalid_and_unrepresentable_inputs_stop(self):
        for df in (True, np.bool_(True), 0, -1, math.inf, math.nan, '5', None):
            with self.subTest(df=df), self.assertRaises(ValueError):
                _critical_t(df, .95)
        for confidence in (True, 0, -1, 1, 2, math.inf, math.nan, '.95', None,
                           math.nextafter(0., 1.), math.nextafter(1., 0.)):
            with self.subTest(confidence=confidence), self.assertRaises(ValueError):
                _critical_t(5, confidence)

    def test_critical_t_cdf_and_bracket_failures_stop(self):
        module = 'audit_system.statistical_rerun_2026_09_20'
        for value in (math.nan, math.inf, -.1, .6):
            with self.subTest(cdf=value), patch(module+'.stdtr', return_value=value), \
                    self.assertRaisesRegex(ValueError, 'CDF'):
                _critical_t(5, .95)
        with patch(module+'.stdtr', return_value=.5) as cdf, \
                self.assertRaisesRegex(ValueError, 'bracket'):
            _critical_t(5, .95)
        self.assertLessEqual(cdf.call_count, 1024)

    def test_critical_t_solver_failure_or_bad_root_is_not_accepted(self):
        target = 'audit_system.statistical_rerun_2026_09_20.optimize.brentq'
        for failure in (RuntimeError('software nonconvergence'), ValueError('software bad bracket')):
            with patch(target, side_effect=failure), self.assertRaisesRegex(ValueError, 'inversion failed'):
                _critical_t(5, .95)
        for root, converged in ((1., True), (math.nan, True), (math.inf, True),
                                (0., True), (2.57, False)):
            with self.subTest(root=root, converged=converged), \
                    patch(target, return_value=(root, SimpleNamespace(converged=converged))), \
                    self.assertRaisesRegex(ValueError, 'accuracy'):
                _critical_t(5, .95)

    def test_critical_t_adjacent_float_refinement_keeps_tail_accuracy_requirement(self):
        confidence = 1-1e-12
        tail = 1-(1+confidence)/2
        critical = _critical_t(100, confidence)
        self.assertLessEqual(abs(float(stdtr(100,-critical))-tail),
                             64*sys.float_info.epsilon*tail)

    def test_all_t_interval_methods_preserve_arm_reversal_symmetry(self):
        ids = [f'GSM{i}' for i in range(1, 9)]
        frame = pd.DataFrame({'clock': [13., 14., 17., 22., 8., 15., 18., 17.]}, index=ids)
        covariates = {'covariates': [{'name': 'baseline',
            'values': dict(zip(ids, [0., 1., 2., 3., 0., 1., 2., 3.]))}]}
        pairs = {'pairs': [dict(unit_id=str(i), case_id=a, control_id=b)
                          for i, (a, b) in enumerate(zip(ids[:4], ids[4:]))]}
        reverse_pairs = {'pairs': [dict(unit_id=p['unit_id'], case_id=p['control_id'],
                                       control_id=p['case_id']) for p in pairs['pairs']]}
        for method, design, reverse in (
                ('Paired t-test', pairs, reverse_pairs), ('Welch t-test', None, None),
                ('Student t-test', None, None), ('OLS adjusted contrast', covariates, covariates),
                ('OLS HC3 adjusted contrast', covariates, covariates)):
            for confidence in (.8, .95, .99):
                with self.subTest(method=method, confidence=confidence):
                    a = recompute_contrast(frame, ids[:4], ids[4:], 'clock', method,
                                           design=design, confidence_level=confidence)
                    b = recompute_contrast(frame, ids[4:], ids[:4], 'clock', method,
                                           design=reverse, confidence_level=confidence)
                    for observed, expected in zip(a['confidence_interval'],
                                                  [-v for v in b['confidence_interval'][::-1]]):
                        self.assertTrue(math.isclose(observed, expected, rel_tol=1e-12, abs_tol=5e-14))
                    self.assertAlmostEqual(a.get('adjusted_mean_difference', a['mean_difference']),
                        -b.get('adjusted_mean_difference', b['mean_difference']), places=12)
                    self.assertEqual(a['scientific_status'], 'CANNOT VERIFY')
                    with patch('audit_system.statistical_rerun_2026_09_20._critical_t',
                               side_effect=ValueError('software inversion failure')), \
                            self.assertRaisesRegex(ValueError, 'inversion failure'):
                        recompute_contrast(frame, ids[:4], ids[4:], 'clock', method,
                                           design=design, confidence_level=confidence)

    def test_paired_map_not_input_order_controls_results_and_ci(self):
        a = recompute_contrast(self.frame, self.case, self.control, 'clock', 'paired t-test', design=self.pairs)
        b = recompute_contrast(self.frame, self.case[::-1], self.control, 'clock', 'paired t-test', design=self.pairs)
        expected = stats.ttest_1samp([3., 4., 2.], 0)
        self.assertAlmostEqual(a['p_value'], expected.pvalue, places=14)
        self.assertEqual(a['confidence_interval'], b['confidence_interval'])
        self.assertEqual(a['n_biological_units'], 3)
        self.assertAlmostEqual(a['effect_size'], 3.)
        self.assertNotEqual(a['effect_size'], a['pooled_standardized_difference'])
        self.assertLess(a['confidence_interval'][0], 3.)
        self.assertGreater(a['confidence_interval'][1], 3.)
        self.assertEqual(a['scientific_status'], 'CANNOT VERIFY')

    def test_legacy_callable_does_not_backfill_original_estimator_declarations(self):
        from src.pipeline_2026_09_10 import recompute_declared_contrast
        legacy = recompute_declared_contrast(self.frame, self.case, self.control, 'clock', 'Welch')
        explicit = recompute_declared_contrast(self.frame, self.case, self.control, 'clock', 'Welch', design={})
        self.assertNotIn('effect_size_definition', legacy)
        self.assertNotIn('scientific_status', legacy)
        self.assertEqual(explicit['effect_size_definition'], POOLED)
        self.assertEqual(explicit['p_value'], legacy['p_value'])

    def test_paired_missing_duplicated_or_invented_units_stop(self):
        for mutate in (
            lambda d: d['pairs'].pop(),
            lambda d: d['pairs'][1].update(unit_id='0'),
            lambda d: d['pairs'][1].update(control_id='D'),
            lambda d: d['pairs'][1].update(case_id='unavailable'),
            lambda d: d['pairs'][1].update(unit_id=''),
        ):
            design = copy.deepcopy(self.pairs); mutate(design)
            with self.subTest(design=design), self.assertRaises(ValueError):
                recompute_contrast(self.frame, self.case, self.control, 'clock', 'paired t-test', design=design)
        with self.assertRaises(NotImplementedError):
            recompute_contrast(self.frame, self.case, self.control, 'clock', 'paired t-test')

    def test_ols_recovers_prespecified_coefficient_and_hc3_is_explicit(self):
        ids = [f's{i}' for i in range(8)]
        g = np.array([1., 1., 1., 1., 0., 0., 0., 0.])
        x = np.array([0., 1., 2., 3., 0., 1., 2., 3.])
        # Residual fixture is orthogonal to intercept, group and covariate.
        e = np.array([1., -1., -1., 1., -2., 2., 2., -2.])
        frame = pd.DataFrame({'clock': 10+2*g+3*x+e}, index=ids)
        design = {'covariates': [{'name': 'baseline', 'values': dict(zip(ids, x))}]}
        result = recompute_contrast(frame, ids[:4], ids[4:], 'clock', 'OLS adjusted contrast', design=design)
        robust = recompute_contrast(frame, ids[:4], ids[4:], 'clock', 'OLS HC3 adjusted contrast', design=design)
        self.assertAlmostEqual(result['adjusted_mean_difference'], 2.)
        self.assertAlmostEqual(result['coefficients'][2], 3.)
        self.assertAlmostEqual(result['standard_error'], 2**.5)
        self.assertAlmostEqual(result['p_value'], 2*stats.t.sf(2/(2**.5), 5))
        self.assertNotEqual(result['standard_error'], robust['standard_error'])
        self.assertIsNone(result['effect_size'])
        self.assertEqual(robust['variance_estimator'], 'HC3 with residual-df t reference')

    def test_covariate_confounding_missingness_and_duplicate_names_stop(self):
        for covariates in (
            [{'name': 'confounded', 'values': dict(zip(list('ABCDEF'), [1, 1, 1, 0, 0, 0]))}],
            [{'name': 'missing', 'values': {'A': 2}}],
            [{'name': 'case', 'values': dict.fromkeys(list('ABCDEF'), 2)}],
            [{'name': 'nonfinite', 'values': dict.fromkeys(list('ABCDEF'), float('nan'))}],
        ):
            with self.subTest(covariates=covariates), self.assertRaises(ValueError):
                recompute_contrast(self.frame, self.case, self.control, 'clock', 'OLS adjusted contrast',
                                   design={'covariates': covariates})

    def test_exact_source_filters_pairs_and_covariates_preserve_membership(self):
        metadata = [{'GSM': s, 'arm': 'case' if s in self.case else 'control',
                     'donor': str(i % 3), 'sex': 'F' if i % 2 else 'M'} for i, s in enumerate('ABCDEF')]
        proposal = {'case_filter': {'arm': 'case'}, 'control_filter': {'arm': 'control'},
                    'pair_field': 'donor', 'covariates': [{'field': 'sex', 'type': 'categorical',
                    'levels': ['M', 'F'], 'role': 'baseline_confounder', 'justification': 'software fixture'}]}
        case, control = select_groups(metadata, proposal)
        design = bind_design(metadata, proposal, case, control)
        self.assertEqual(design['pairs'], self.pairs['pairs'])
        self.assertEqual(len(design['covariates']), 1)
        metadata[0]['donor'] = 'unmatched'
        with self.assertRaisesRegex(ValueError, 'Unmatched'):
            bind_design(metadata, proposal, case, control)

    def test_duplicate_clock_and_sample_identity_stop_without_alias_guessing(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'predictions.csv'
            path.write_text('sample_id,clock,clock,other\nA,1,2,3\nB,2,4,5\n')
            with self.assertRaisesRegex(ValueError, 'ambiguous'):
                clock_frame(path, 'clock')
            self.assertEqual(clock_frame(path, 'other').loc['A', 'other'], '3')
            path.write_text('sample_id,gsm_id,clock\nA,B,3\n')
            with self.assertRaisesRegex(ValueError, 'conflicting'):
                clock_frame(path, 'clock')

    def label_fixture(self):
        """Constructed label fixture; source-backed study checks are separate."""
        metadata = [{'GSM': sid, 'title': f'1HR_{arm}_ANIMAL_{unit}', 'source': 'fixture'}
                    for arm, ids in [('LTP', self.case), ('CTRL', self.control)]
                    for unit, sid in enumerate(ids, 1)]
        metadata += [{'GSM': sid+'later', 'title': row['title'].replace('1HR', '3HR'), 'source': 'fixture'}
                     for sid, row in zip(self.case+self.control, list(metadata))]
        rule = {'schema': 'clockbase-delimited-pair-v1', 'delimiter': '_', 'token_count': 4,
                'fixed_tokens': {'2': 'ANIMAL'}, 'stratum_tokens': {'0': '1HR'},
                'unit_tokens': [2, 3], 'arm_token': 1, 'case_token': 'LTP', 'control_token': 'CTRL',
                'scope_filters': {'source': 'fixture'}}
        return metadata, {'pair_field': 'title', 'pair_identity': rule}

    def test_derived_pairs_preserve_raw_labels_complete_membership_and_time(self):
        metadata, proposal = self.label_fixture()
        before = copy.deepcopy((metadata, proposal))
        result = bind_design(metadata, proposal, self.case, self.control)
        shuffled = bind_design(metadata[::-1], proposal, self.case[::-1], self.control[::-1])
        self.assertEqual(result, shuffled)
        self.assertEqual((metadata, proposal), before)
        self.assertEqual([(p['case_id'], p['control_id']) for p in result['pairs']], list(zip(self.case, self.control)))
        self.assertEqual(result['pair_identity']['counts'], {'source_records': 12, 'scoped_records': 12,
                         'paired_label_units': 3, 'included_samples': 6, 'excluded_samples': 6})
        self.assertEqual({e['reason'] for e in result['pair_identity']['excluded_samples']}, {'OTHER_DECLARED_STRATUM'})
        self.assertEqual(result['pair_identity']['scientific_status'], 'CANNOT VERIFY')
        for witness in result['pair_identity']['raw_field_witnesses']:
            self.assertEqual(witness['value'], next(r['title'] for r in metadata if r['GSM'] == witness['sample_id']))
            self.assertEqual(json.loads(witness['derived_unit_id'])['stratum'], ['1HR'])

    def test_derived_pairs_never_pool_times_or_drop_unselected_matching_animals(self):
        metadata, proposal = self.label_fixture()
        with self.assertRaisesRegex(ValueError, 'complete source groups'):
            bind_design(metadata, proposal, self.case[:2], self.control[:2])
        with self.assertRaisesRegex(ValueError, 'complete source groups'):
            bind_design(metadata, proposal, self.case+[s+'later' for s in self.case],
                        self.control+[s+'later' for s in self.control])
        later = copy.deepcopy(proposal); later['pair_identity']['stratum_tokens']['0'] = '3HR'
        a = bind_design(metadata, proposal, self.case, self.control)
        b = bind_design(metadata, later, [s+'later' for s in self.case], [s+'later' for s in self.control])
        self.assertFalse({p['unit_id'] for p in a['pairs']} & {p['unit_id'] for p in b['pairs']})

    def test_derived_pair_rule_rejects_ambiguous_roles_or_undeclared_tokens(self):
        metadata, proposal = self.label_fixture()
        for mutation in (
            {'schema': 'unknown'}, {'delimiter': ''}, {'token_count': True}, {'token_count': 33},
            {'fixed_tokens': {'02': 'ANIMAL'}}, {'fixed_tokens': {'-1': 'ANIMAL'}},
            {'stratum_tokens': {}}, {'unit_tokens': [True, 3]}, {'unit_tokens': [3, 3]},
            {'unit_tokens': [0, 3]}, {'unit_tokens': [1, 3]}, {'unit_tokens': []},
            {'fixed_tokens': {}, 'unit_tokens': [3]}, {'arm_token': True},
            {'control_token': 'LTP'}, {'scope_filters': {'source': None}}, {'regex': '.*'},
        ):
            changed = copy.deepcopy(proposal); changed['pair_identity'].update(mutation)
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                bind_design(metadata, changed, self.case, self.control)

    def test_derived_pair_metadata_mismatch_never_aggregates_or_deletes(self):
        metadata, proposal = self.label_fixture()
        for mutate in (
            lambda rows: rows[0].update(title=''),
            lambda rows: rows[0].update(title='1HR_LTP_ANIMAL_1_extra'),
            lambda rows: rows[0].update(title='1HR_LTP__1'),
            lambda rows: rows[0].update(title='1HR_LTP_REPLICATE_1'),
            lambda rows: rows[0].update(title='1HR_CTRL_ANIMAL_1'),
            lambda rows: rows[0].update(title='1HR_OTHER_ANIMAL_1'),
            lambda rows: rows[1].update(title='1HR_LTP_ANIMAL_1'),
            lambda rows: rows[3].update(title='1HR_CTRL_ANIMAL_4'),
            lambda rows: rows[1].update(GSM='A'),
            lambda rows: rows.pop(3),
        ):
            changed = copy.deepcopy(metadata); mutate(changed)
            with self.subTest(metadata=changed), self.assertRaises(ValueError):
                bind_design(changed, proposal, self.case, self.control)

    def test_derived_pair_scope_exclusions_are_preserved_and_cannot_hide_selected_rows(self):
        metadata, proposal = self.label_fixture()
        metadata.append({'GSM': 'outside', 'source': 'another assay', 'title': 'different label format'})
        result = bind_design(metadata, proposal, self.case, self.control)
        excluded = {row['sample_id']: row for row in result['pair_identity']['excluded_samples']}
        self.assertEqual(excluded['outside']['reason'], 'OUTSIDE_LITERAL_SCOPE')
        self.assertEqual(excluded['outside']['scope_field_values'], {'source': 'another assay'})
        metadata[0]['source'] = 'another assay'
        with self.assertRaisesRegex(ValueError, 'complete source groups'):
            bind_design(metadata, proposal, self.case, self.control)

    def test_derived_pairs_require_literal_source_field_and_disjoint_ids(self):
        metadata, proposal = self.label_fixture()
        for change in ({'pair_field': None}, {'pair_field': 'not_a_field'}):
            with self.subTest(change=change), self.assertRaises(ValueError):
                bind_design(metadata, proposal | change, self.case, self.control)
        for case, control in ((self.case, self.case), (self.case+[self.case[0]], self.control), ([], self.control)):
            with self.subTest(case=case, control=control), self.assertRaises(ValueError):
                bind_design(metadata, proposal, case, control)

    def test_derived_pairs_use_existing_conditional_arithmetic_without_acceptance(self):
        metadata, proposal = self.label_fixture()
        derived = bind_design(metadata, proposal, self.case, self.control)
        result = recompute_contrast(self.frame, self.case, self.control, 'clock', 'paired t-test', design=derived)
        explicit = recompute_contrast(self.frame, self.case, self.control, 'clock', 'paired t-test', design=self.pairs)
        for key in ('p_value', 'mean_difference', 'effect_size', 'confidence_interval'):
            self.assertEqual(result[key], explicit[key])
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')

    def test_unnamed_gsm_index_requires_exact_declared_source_column(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'predictions.csv'
            path.write_text(',clock\nGSM2,2\nGSM1,4\n')
            with self.assertRaisesRegex(ValueError, 'explicit sample identity'):
                clock_frame(path, 'clock')
            binding = {'index': 0, 'header': '', 'format': 'GSM accession'}
            self.assertEqual(clock_frame(path, 'clock', sample_id_column=binding).loc['GSM1', 'clock'], '4')
            for mutation in ({'index': 1}, {'header': 'sample_id'}, {'format': 'row number'}):
                with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                    clock_frame(path, 'clock', sample_id_column=binding|mutation)
            path.write_text(',clock\n0,2\n1,4\n')
            with self.assertRaisesRegex(ValueError, 'literal GSM'):
                clock_frame(path, 'clock', sample_id_column=binding)

    def make_plan(self, root):
        for name, text in {
            'metadata.csv': 'GSM,arm,donor\nA,case,1\nB,case,2\nC,case,3\nD,control,1\nE,control,2\nF,control,3\n',
            'predictions.csv': 'sample_id,clock\nA,4\nB,8\nC,7\nD,1\nE,4\nF,5\n',
            'original.csv': 'gse_id,clock,p_value,effect_size\nGSE1,clock,0.5,invalid_function\n',
        }.items():
            (root/name).write_text(text)
        (root/'reference.json').write_text(json.dumps({'target': {'case_ids': self.case, 'control_ids': self.control,
            'pair_map': [{'donor': str(i+1), 'case_gsm': a, 'control_gsm': b}
                         for i, (a, b) in enumerate(zip(self.case, self.control))]}, 'method': 'paired_t_by_donor'}))
        identity = lambda name: {'path': str(root/name), 'sha256': sha(root/name)}
        plan = {'schema': SCHEMA, 'evaluation_role': 'exposed_development',
                'multiplicity': {'action': 'withhold_adjusted_inference', 'reason': 'fixture family unresolved'},
                'confidence_level': .95, 'relative_tolerance': 1e-8, 'absolute_tolerance': 1e-12,
                'independent_design_reference': identity('reference.json'),
                'sources': {'GSE1': {k: identity(k+'.csv') for k in ('metadata', 'predictions', 'original')}},
                'proposals': [{'proposal_id': 'paired-fixture', 'study_id': 'GSE1', 'reference_target': ['target'],
                'reference_method': ['method'],
                'clock': 'clock', 'original_row_ordinal': 1, 'case_filter': {'arm': 'case'},
                'original_record_sha256': record_hash({'gse_id': 'GSE1', 'clock': 'clock', 'p_value': '0.5', 'effect_size': 'invalid_function'}),
                'control_filter': {'arm': 'control'}, 'mode': 'inferential_conditional', 'pair_field': 'donor',
                'statistical_test': 'paired t-test', 'scientific_limitations': ['Software fixture; no scientific evidence'],
                'saved_effect_definition': POOLED}]}
        path = root/'plan.json'; path.write_text(json.dumps(plan)); return path, plan

    def test_runner_preserves_invalid_original_and_verifies_resume(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path, plan = self.make_plan(root)
            summary = run(path, root/'out')
            self.assertEqual(summary['rerun_targets'], 1)
            self.assertEqual(summary['accepted_scientific_repairs'], 0)
            self.assertEqual(run(path, root/'out'), summary)
            result_path = root/'out/paired-fixture/attempt-000/result.json'
            result = json.loads(result_path.read_text())
            self.assertEqual(result['original']['effect_size'], 'invalid_function')
            self.assertEqual(result['detections'][0]['field'], 'effect_size')
            self.assertEqual(result['release_status'], 'WITHHELD')
            result_path.write_text('{}')
            with self.assertRaisesRegex(PipelineIntegrityError, 'output changed'):
                run(path, root/'out')

    def test_runner_missing_sample_stops_without_imputation(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path, plan = self.make_plan(root)
            prediction = root/'predictions.csv'; prediction.write_text(prediction.read_text().replace('F,5\n', ''))
            reference = root/'reference.json'; ref = json.loads(reference.read_text())
            ref['target']['missing_prediction_ids'] = ['F']; reference.write_text(json.dumps(ref))
            plan['independent_design_reference']['sha256'] = sha(reference)
            plan['sources']['GSE1']['predictions']['sha256'] = sha(prediction); path.write_text(json.dumps(plan))
            summary = run(path, root/'out')
            result = json.loads((root/'out/paired-fixture/attempt-000/result.json').read_text())
            self.assertEqual(summary['unresolved_targets'], 1)
            self.assertEqual(result['missing_prediction_ids'], ['F'])
            self.assertIsNone(result['candidate'])

    def test_runner_refuses_changed_inputs_and_reference_groups(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path, plan = self.make_plan(root)
            (root/'metadata.csv').write_text('altered')
            with self.assertRaisesRegex(PipelineIntegrityError, 'Frozen source changed'):
                run(path, root/'out')
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path, plan = self.make_plan(root)
            reference = root/'reference.json'; reference.write_text(json.dumps({'target': {'case_ids': ['A'], 'control_ids': ['D']}}))
            plan['independent_design_reference']['sha256'] = sha(reference); path.write_text(json.dumps(plan))
            with self.assertRaisesRegex(PipelineIntegrityError, 'reference'):
                run(path, root/'out')

    def test_runner_cannot_change_reviewed_method_to_make_result_pass(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path, plan = self.make_plan(root)
            plan['proposals'][0]['statistical_test'] = 'Welch'
            path.write_text(json.dumps(plan))
            with self.assertRaisesRegex(PipelineIntegrityError, 'mode/test'):
                run(path, root/'out')

    def test_runner_checks_exact_reference_pairing_not_only_arm_sets(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); path, plan = self.make_plan(root)
            reference = root/'reference.json'; ref = json.loads(reference.read_text())
            ref['target']['pair_map'][0]['control_gsm'] = 'E'
            ref['target']['pair_map'][1]['control_gsm'] = 'D'
            reference.write_text(json.dumps(ref)); plan['independent_design_reference']['sha256'] = sha(reference)
            path.write_text(json.dumps(plan))
            with self.assertRaisesRegex(PipelineIntegrityError, 'pair map'):
                run(path, root/'out')

    def test_runner_checks_reference_missing_set_and_original_row_identity(self):
        for change in ('missing', 'original'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); path, plan = self.make_plan(root)
                if change == 'missing':
                    reference = root/'reference.json'; ref = json.loads(reference.read_text())
                    ref['target']['missing_prediction_ids'] = ['A']; reference.write_text(json.dumps(ref))
                    plan['independent_design_reference']['sha256'] = sha(reference)
                else:
                    plan['proposals'][0]['original_record_sha256'] = '0'*64
                path.write_text(json.dumps(plan))
                with self.assertRaises(PipelineIntegrityError):
                    run(path, root/'out')


if __name__ == '__main__':
    unittest.main()
