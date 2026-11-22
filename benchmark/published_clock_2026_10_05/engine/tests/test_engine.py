"""Synthetic fixtures for source selection, grouping, matching, statistics and reports.

All numeric arrays in this module are synthetic test inputs, not empirical results.
Fixtures are explicit constants; the only random process is the seeded paper bootstrap.
"""
import hashlib
import json
import os
import tempfile
import unittest
import warnings
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

from audit_system import paths
from audit_system import membership_recovery_2026_09_27 as membership
from benchmark.published_clock_2026_10_05.engine import cli, groups, match, metrics, predictions, stats
from benchmark.published_clock_2026_10_05.engine.io import clean, receipt, write_json


def record(identity='synthetic-1', **kwargs):
    return {'id': identity, 'tier': 'T1', 'gse': 'GSE1', 'modality': 'human DNAm',
            'source': {'doi': 'synthetic-paper-1', 'quote': 'Synthetic fixture; no published result.',
                       'citation': 'synthetic collection'}, 'clock_reported': 'Horvath', 'clock_clockbase': 'Horvathv1',
            'measure': 'age', 'case_label': 'case', 'control_label': 'control',
            'case_definition': 'group=case', 'control_definition': 'group=control', 'direction': 'higher',
            'significant': True, 'threshold': 'P < 0.05', 'p_value': 0.01, 'effect': None, **kwargs}


def atlas_row(ordinal=1, **kwargs):
    return {'row_ordinal': ordinal, 'gse_id': 'GSE1', 'clock': 'Horvathv1', 'comparison_name': 'case vs control',
            'condition_name': 'case', 'case_gsms': ['GSM1', 'GSM2', 'GSM3'],
            'control_gsms': ['GSM4', 'GSM5', 'GSM6'], 'final_class': 'VERIFIED',
            'historical': {'pvalue': 0.01, 'effect_size': 2.0, 'direction': 'higher'},
            'recomputed': {'p': 0.02, 'cohen_d': 1.0, 'mean_diff': 3.0, 'n_case': 3, 'n_control': 3}, **kwargs}


def arms(**kwargs):
    return {'status': 'MANUAL', 'case_gsms': ['GSM1', 'GSM2', 'GSM3'],
            'control_gsms': ['GSM4', 'GSM5', 'GSM6'], **kwargs}


class FixtureCase(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='synthetic-fixture-', dir=Path(__file__).parent)
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        for name in ('GEO', 'VERIFIER_ROOT', 'STREAM_EMPTY', 'STREAM_DETACHED'):
            context = patch.object(paths, name, self.root / name)
            context.start()
            self.addCleanup(context.stop)
        context = patch.dict(os.environ, {'AUDIT_CACHE_ONLY': '1', 'BENCH_GROUP_MODELS': 'synthetic-a,synthetic-b'})
        context.start()
        self.addCleanup(context.stop)

    def csv(self, path, frame, index=False):
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=index)
        return path

    def json(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        write_json(path, value)
        return path

    def metadata(self):
        meta = pd.DataFrame({'GSM': [f'GSM{i}' for i in range(1, 7)], 'group': ['case'] * 3 + ['control'] * 3,
                             'harmonized.age': ['10', '20', '30', '10', '20', '30'],
                             'harmonized.sex': ['male', 'female', 'female', 'male', 'female', 'female']})
        self.csv(paths.GEO / 'GEO_hs_dnam/processed_metadata/GSE1.csv', meta)
        return meta.set_index('GSM')

    def adjudicator(self, answers):
        iterator = iter(answers)
        def transport(model, prompt, ledger, budget, **kwargs):
            return {'status': 'OK', 'parsed': next(iterator), 'response_text': 'synthetic'}, False
        return groups.Adjudicator(self.root / 'ledger.jsonl', transport=transport)


class PredictionTests(FixtureCase):
    def test_historical_qced_and_raw_and_grimage(self):
        self.metadata()
        raw = self.csv(paths.GEO / 'result/hs_dnam_age.csv.gz', pd.DataFrame(
            {'GSE_ID': ['GSE1'], 'GSM_ID': ['GSM1'], 'Horvathv1': [11], 'DNAmGrimAge': [12], 'DNAmGrimAge.1': [13]}))
        qced = self.csv(paths.GEO / 'result/qced_human_dnam_age.csv', pd.DataFrame(
            {'gse_id': ['GSE1'], 'gsm_id': ['GSM1'], 'Horvathv1': [21]}))
        reader = predictions.PredictionReader()
        previous = membership._STATE.copy()
        value = reader.load('human DNAm', 'GSE1', 'Horvathv1')
        self.assertEqual(value['samples'].loc['GSM1', 'clock'], 21)
        self.assertEqual(value['samples'].loc['GSM1', 'age'], 10)
        self.assertEqual(value['samples'].loc['GSM1', 'sex'], 'male')
        self.assertEqual(value['source']['sha256'], hashlib.sha256(qced.read_bytes()).hexdigest())
        self.assertEqual(membership._STATE, previous)
        value = reader.load('dnam', 'GSE1', 'GrimAgeV2')
        self.assertEqual(value['samples'].loc['GSM1', 'clock'], 13)
        self.assertEqual(value['source']['sha256'], receipt(raw)['sha256'])
        preferred = reader.load('dnam', 'GSE1', 'Horvathv1', [atlas_row(evidence={'prediction_source': 'raw'})])
        self.assertEqual(preferred['samples'].loc['GSM1', 'clock'], 11)

    def test_rna_mean_and_duplicate_exclusion(self):
        self.csv(paths.GEO / 'GEO_hs_rna/clock_result/GSE2.csv', pd.DataFrame(
            {'sample_id': ['GSM1', 'GSM2', 'GSM2'], 'tAge_mean': [10, 20, 30]}))
        value = predictions.PredictionReader().load('human RNA-seq', 'GSE2', 'tAge')
        self.assertEqual(value['source_kind'], 'clock_result')
        self.assertEqual(list(value['samples'].index), ['GSM1'])
        self.assertEqual(value['duplicated_sample_ids'], ['GSM2'])
        self.assertTrue(np.isnan(value['samples'].loc['GSM1', 'age']))

    def test_repaired_retained_only_and_root_precedence(self):
        for root, output in [(paths.STREAM_EMPTY, 30), (paths.STREAM_DETACHED, 40)]:
            study = root / 'studies/GSE1'
            self.json(study / 'status.json', {'state': 'REPAIRED', 'clock_verdicts': [{'clock': 'Horvathv1', 'verdict': 'VERIFIED'}]})
            self.csv(study / 'clock_predictions.csv', pd.DataFrame({'Horvathv1': [output, 99]}, index=['GSM1', 'GSM2']), True)
            self.csv(study / 'sample_status.csv', pd.DataFrame({'status': ['RETAINED', 'EXCLUDED']}, index=['GSM1', 'GSM2']), True)
        value = predictions.PredictionReader().load('dnam', 'GSE1', 'Horvathv1')
        self.assertEqual(value['samples'].loc['GSM1', 'clock'], 30)
        self.assertNotIn('GSM2', value['samples'].index)
        self.assertIn('sample_status.csv', value['source'])

    def test_recorded_replay_overrides_historical(self):
        self.metadata()
        path = self.csv(paths.VERIFIER_ROOT / 'replay/inputs/GSE1.csv', pd.DataFrame(
            {'gsm_id': ['GSM1', 'GSM2'], 'Horvathv1': [50, 60]}))
        rows = [atlas_row(predictions=receipt(path))]
        value = predictions.PredictionReader().load('dnam', 'GSE1', 'Horvathv1', rows)
        self.assertEqual(value['samples'].loc['GSM1', 'clock'], 50)
        self.assertEqual(value['source_kind'], 'verifier_replay')

    def test_missing_and_conflicting_and_mismatched_sources(self):
        reader = predictions.PredictionReader()
        value = reader.load('dnam', 'GSE1', 'Horvathv1', [atlas_row(predictions={'path': str(self.root / 'missing.csv')})])
        self.assertIn('Cannot proceed - required data file not found:', value['error'])
        rows = [atlas_row(predictions={'path': 'a'}), atlas_row(2, predictions={'path': 'b'})]
        self.assertEqual(reader.load('dnam', 'GSE1', 'Horvathv1', rows)['status'], 'AMBIGUOUS_PREDICTION_SOURCE')
        path = self.csv(paths.VERIFIER_ROOT / 'replay/inputs/GSE1.csv', pd.DataFrame({'gsm_id': ['GSM1'], 'Horvathv1': [1]}))
        rows = [atlas_row(predictions={'path': str(path), 'sha256': 'bad'})]
        self.assertEqual(reader.load('dnam', 'GSE1', 'Horvathv1', rows)['status'], 'SOURCE_HASH_MISMATCH')
        self.assertEqual(reader.load('dnam', 'GSE1', None)['status'], 'NO_CLOCK_MAPPING')

    def test_withheld_repaired_clock(self):
        study = paths.STREAM_DETACHED / 'studies/GSE1'
        self.json(study / 'status.json', {'state': 'REPAIRED', 'clock_verdicts': [{'clock': 'Horvathv1', 'verdict': 'WITHHELD'}]})
        self.csv(study / 'clock_predictions.csv', pd.DataFrame({'Horvathv1': [1]}, index=['GSM1']), True)
        self.csv(study / 'sample_status.csv', pd.DataFrame({'status': ['RETAINED']}, index=['GSM1']), True)
        row = atlas_row(predictions=receipt(study / 'clock_predictions.csv'))
        value = predictions.PredictionReader().load('dnam', 'GSE1', 'Horvathv1', [row])
        self.assertEqual(value['status'], 'CLOCK_WITHHELD')

    def test_conflicting_hashes_at_same_path(self):
        reader = predictions.PredictionReader()
        rows = [atlas_row(predictions={'path': 'same.csv', 'sha256': 'one'}),
                atlas_row(2, predictions={'path': 'same.csv', 'sha256': 'two'})]
        value = reader.load('dnam', 'GSE1', 'Horvathv1', rows)
        self.assertEqual(value['status'], 'AMBIGUOUS_PREDICTION_SOURCE')
        self.assertEqual(len(value['sources']), 2)


class GroupTests(FixtureCase):
    def test_different_rules_same_selections_agree(self):
        meta = self.metadata()
        meta['diagnosis'] = meta['group']
        a = {'feasible': True, 'case': {'group': ['case']}, 'control': {'group': ['control']}}
        b = {'feasible': True, 'case': {'diagnosis': ['case']}, 'control': {'diagnosis': ['control']}}
        result = groups.assign_groups(record(), meta, self.adjudicator([a, b]))
        self.assertEqual(result['status'], 'AGREED')
        self.assertEqual(result['case_gsms'], arms()['case_gsms'])
        self.assertEqual(len(result['answers']), 2)

    def test_disagreement_keeps_both_selections(self):
        a = {'feasible': True, 'case': {'group': ['case']}, 'control': {'group': ['control']}}
        b = {'feasible': True, 'case': {'group': ['control']}, 'control': {'group': ['case']}}
        result = groups.assign_groups(record(), self.metadata(), self.adjudicator([a, b]))
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertNotEqual(result['selections'][0]['case_gsms'], result['selections'][1]['case_gsms'])

    def test_bad_rules_overlap_and_missing_ids(self):
        meta = self.metadata()
        for rules in ({}, {'missing': ['case']}, {'group': ['invented']}, {'group': 'case'}, {'group': []}):
            with self.subTest(rules=rules), self.assertRaises(ValueError):
                groups.apply_rules(meta, rules)
        for manual in ({'case_gsms': ['GSM1'], 'control_gsms': ['GSM1']},
                       {'case_gsms': ['GSM404'], 'control_gsms': ['GSM4']},
                       {'case_gsms': ['GSM1', 'GSM1'], 'control_gsms': ['GSM4']}):
            self.assertEqual(groups.assign_groups(record(), meta, None, manual)['status'], 'UNRESOLVED')

    def test_manual_and_rules_intersection(self):
        meta = self.metadata()
        manual = {**arms(), 'provenance': 'synthetic fixture'}
        self.assertEqual(groups.assign_groups(record(), meta, None, manual)['status'], 'MANUAL')
        self.assertEqual(groups.apply_rules(meta, {'group': ['case'], 'harmonized.sex': ['female']}), ['GSM2', 'GSM3'])

    def test_cache_only_has_no_transport(self):
        adjudicator = groups.Adjudicator(self.root / 'ledger.jsonl')
        with patch('urllib.request.urlopen', side_effect=AssertionError('network forbidden')):
            result = groups.assign_groups(record(), self.metadata(), adjudicator)
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertEqual([a['status'] for a in result['answers']], ['NOT_CACHED', 'NOT_CACHED'])
        self.assertEqual(adjudicator.budget['sent'], 0)
        self.assertEqual(len((self.root / 'ledger.jsonl').read_text().splitlines()), 2)

    def test_fake_http_transport_uses_verifier_ledger_and_cache(self):
        rules = {'feasible': True, 'case': {'group': ['case']}, 'control': {'group': ['control']}}
        payload = {'model': 'synthetic-returned', 'usage': {'total_tokens': 5},
                   'choices': [{'message': {'content': json.dumps(rules)}}]}
        class Response:
            def __enter__(self):
                return self
            def __exit__(self, *args):
                return False
            def read(self):
                return json.dumps(payload).encode()
        with patch.dict(os.environ, {'AUDIT_CACHE_ONLY': '0', 'TOKENSWITCH_BASE_URL': 'https://synthetic.invalid',
                                    'TOKENSWITCH_API_KEY': 'synthetic-test-token'}), patch('urllib.request.urlopen', return_value=Response()) as http:
            adjudicator = groups.Adjudicator(self.root / 'ledger.jsonl', max_requests=2, cache_only=False)
            result = groups.assign_groups(record(), self.metadata(), adjudicator)
            self.assertEqual(result['status'], 'AGREED')
            adjudicator.cache_only = True
            replay = groups.assign_groups(record(), self.metadata(), adjudicator)
            self.assertEqual(http.call_count, 2)
            self.assertTrue(all(a['cached'] for a in replay['answers']))
            self.assertEqual(adjudicator.budget['sent'], 2)
            self.assertTrue(all(a['request_sha256'] and a['response_text'] for a in result['answers']))

    def test_budget_and_model_pair(self):
        with self.assertRaises(ValueError):
            groups.Adjudicator(self.root / 'ledger.jsonl', models=['same', 'same'])
        with patch.dict(os.environ, {'AUDIT_CACHE_ONLY': '0'}), patch('urllib.request.urlopen', side_effect=AssertionError()):
            adjudicator = groups.Adjudicator(self.root / 'ledger.jsonl', cache_only=False, max_requests=0)
            self.assertEqual(adjudicator.ask('synthetic')[0]['status'], 'BUDGET_STOP')


class StatisticsTests(unittest.TestCase):
    def setUp(self):
        self.samples = pd.DataFrame({'clock': [14, 26, 35, 10, 20, 31, 45], 'age': [10, 20, 30, 10, 20, 30, 40]},
                                    index=[f'GSM{i}' for i in range(1, 8)])

    def test_welch_pooled_d_and_mean_difference(self):
        result = stats.evaluate(self.samples, arms(), 'age', 0.05)
        a, b = self.samples.clock.iloc[:3], self.samples.clock.iloc[3:6]
        self.assertAlmostEqual(result['p'], ttest_ind(a, b, equal_var=False).pvalue)
        self.assertAlmostEqual(result['mean_diff'], a.mean() - b.mean())
        self.assertAlmostEqual(result['cohen_d'], (a.mean() - b.mean()) / np.sqrt((a.var() + b.var()) / 2))
        self.assertEqual(result['n_case'], 3)
        self.assertEqual(result['direction'], 'higher')

    def test_residual_fit_includes_ungrouped_sample(self):
        result = stats.evaluate(self.samples, arms(), 'acceleration_residual', 0.05)
        design = np.column_stack([np.ones(7), self.samples.age])
        beta = np.linalg.lstsq(design, self.samples.clock, rcond=None)[0]
        residual = self.samples.clock - design @ beta
        self.assertAlmostEqual(result['mean_diff'], residual.iloc[:3].mean() - residual.iloc[3:6].mean())
        self.assertEqual(result['fit']['n'], 7)
        self.assertIn('GSM7', result['fit']['gsms'])
        self.assertAlmostEqual(result['fit']['slope'], beta[1])

    def test_difference_and_rate(self):
        result = stats.evaluate(self.samples, arms(), 'acceleration_difference')
        values = self.samples.clock - self.samples.age
        self.assertAlmostEqual(result['mean_diff'], values.iloc[:3].mean() - values.iloc[3:6].mean())
        self.assertEqual(stats.evaluate(self.samples, arms(), 'rate')['p'], stats.evaluate(self.samples, arms(), 'age')['p'])
        self.assertIsNone(result['significant'])

    def test_missing_ages_and_descriptive_arm(self):
        self.samples.loc['GSM1', 'age'] = np.nan
        result = stats.evaluate(self.samples, arms(), 'acceleration_residual', 0.05)
        self.assertEqual(result['status'], 'DESCRIPTIVE')
        self.assertEqual(result['excluded_case_gsms'], ['GSM1'])
        self.assertIsNone(result.get('p'))
        self.assertIsNone(result['significant'])
        self.assertEqual(result['n_case'], 2)

    def test_unusable_fit_and_unsupported_measure(self):
        self.samples['age'] = 5
        self.assertEqual(stats.evaluate(self.samples, arms(), 'acceleration_residual')['status'], 'INSUFFICIENT_RECORDED_AGE')
        self.assertEqual(stats.evaluate(self.samples, arms(), 'other')['status'], 'UNSUPPORTED_MEASURE')
        self.assertFalse(stats.evaluate(self.samples, arms(status='UNRESOLVED'), 'age')['covered'])

    def test_zero_variance_and_empty_arms(self):
        self.samples['clock'] = 5
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)
            result = stats.evaluate(self.samples, arms(), 'age')
        self.assertIsNone(result['cohen_d'])
        self.assertIsNone(result['p'])
        result = stats.evaluate(self.samples, arms(case_gsms=['GSM404']), 'age')
        self.assertFalse(result['covered'])
        self.assertEqual(result['status'], 'EMPTY_ARM')

    def test_threshold(self):
        for text in ('P < 0.05', 'p ≤ .05', 'P<5e-2', 0.05):
            self.assertEqual(stats.threshold({'threshold': text}), 0.05)
        for text in (None, 'FDR < 0.05', 'P < 0.05 adjusted', 'not stated', True):
            self.assertIsNone(stats.threshold({'threshold': text}))


class MatchTests(FixtureCase):
    def test_swapped_orientation_and_all_candidates(self):
        row = atlas_row(case_gsms=arms()['control_gsms'], control_gsms=arms()['case_gsms'])
        result = match.rank_rows(record(), [row, atlas_row(2, case_gsms=['GSM404'])], arms())
        self.assertEqual(result['selected_row_ordinal'], 1)
        self.assertEqual(len(result['candidates']), 2)
        self.assertTrue(result['candidates'][0]['swapped'])
        self.assertEqual(result['candidates'][0]['jaccard_case'], 1)
        self.assertEqual(match.atlas_result(record(), result['candidates'][0], 'S2')['cohen_d'], -1)

    def test_boundary_and_both_arms_required(self):
        result = match.rank_rows(record(), [atlas_row(control_gsms=['GSM4', 'GSM5', 'GSM6', 'GSM7', 'GSM8', 'GSM9'])], arms())
        self.assertEqual(result['status'], 'MATCHED')
        self.assertEqual(result['candidates'][0]['jaccard_control'], 0.5)
        result = match.rank_rows(record(), [atlas_row(control_gsms=['GSM404'])], arms())
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertEqual(match.rank_rows(record(), [atlas_row(clock='other')], arms())['candidates'], [])

    def test_rows_without_recovered_samples_do_not_match(self):
        result = match.rank_rows(record(), [atlas_row(case_gsms=None, control_gsms=None)], arms())
        self.assertEqual(result['status'], 'UNRESOLVED')
        self.assertEqual(result['candidates'][0]['jaccard_case'], 0)

    def test_text_model_agreement_and_disagreement(self):
        answer = {'candidates': [{'row_ordinal': 1, 'orientation': 'SWAPPED'}]}
        result = match.rank_rows(record(), [atlas_row()], {'status': 'UNRESOLVED'}, self.adjudicator([answer, answer]))
        self.assertEqual(result['status'], 'MATCHED')
        self.assertTrue(result['candidates'][0]['swapped'])
        other = {'candidates': [{'row_ordinal': 1, 'orientation': 'SAME'}]}
        result = match.rank_rows(record(), [atlas_row()], {'status': 'UNRESOLVED'}, self.adjudicator([answer, other]))
        self.assertEqual(result['status'], 'UNRESOLVED')

    def test_pending_descriptive_and_historical_swap(self):
        candidate = match.rank_rows(record(), [atlas_row(final_class='PENDING_INPUT')], arms())['candidates'][0]
        self.assertTrue(match.atlas_result(record(), candidate, 'S1')['covered'])
        self.assertFalse(match.atlas_result(record(), candidate, 'S2')['covered'])
        candidate['row']['final_class'] = 'DESCRIPTIVE'
        self.assertIsNone(match.atlas_result(record(), candidate, 'S2')['p'])
        candidate['row']['actions'] = [{'action': 'SWAP_ARMS_CONFIRMED_INVERSION'}]
        self.assertEqual(match.atlas_result(record(), candidate, 'S1')['direction'], 'lower')
        self.assertEqual(match.atlas_result(record(), candidate, 'S2')['direction'], 'higher')


class MetricTests(unittest.TestCase):
    def items(self):
        items = []
        for i in range(6):
            item = record(str(i), source={'doi': f'synthetic-paper-{i // 2}', 'quote': 'synthetic'},
                          significant=i < 4, direction='higher' if i < 4 else 'none')
            item['systems'] = {'S3': {'covered': i != 0, 'direction': 'higher', 'significant': i < 3,
                                     'cohen_d': float(i), 'measure': 'age'}}
            items.append(item)
        return items

    def test_denominators_and_missing_coverage(self):
        point = metrics.point_metrics(self.items(), 'S3')
        self.assertEqual(point['coverage'], 5 / 6)
        self.assertEqual(point['direction_denominator'], 3)
        self.assertEqual(point['direction_agreement'], 1)
        self.assertEqual(point['significance_denominator'], 5)
        self.assertEqual(point['significance_agreement'], 4 / 5)
        self.assertEqual(point['nonsignificant_agreement'], 1)

    def test_paper_bootstrap_seed_and_single_paper(self):
        items = self.items()
        a = metrics.summarize(items, 'S3', 50, 123)
        self.assertEqual(a, metrics.summarize(items, 'S3', 50, 123))
        self.assertEqual(a['papers'], 3)
        self.assertEqual(a['intervals']['coverage']['valid_resamples'], 50)
        low = a['intervals']['coverage']['low']
        self.assertLessEqual(low, a['coverage'])
        single = metrics.summarize(items[:2], 'S3', 50)
        self.assertIsNone(single['intervals']['coverage']['low'])

    def test_unit_compatible_effect_correlation(self):
        items = self.items()
        for i, item in enumerate(items):
            item['effect'] = {'value': float(i), 'unit': 'SD', 'type': "Cohen's d"}
        point = metrics.point_metrics(items, 'S3')
        self.assertAlmostEqual(point['effect_spearman'], 1)
        self.assertEqual(point['effect_denominator'], 5)
        items[2]['systems']['S3']['measure'] = 'rate'
        items[3]['effect']['type'] = 'beta'
        self.assertEqual(metrics.point_metrics(items, 'S3')['effect_denominator'], 3)
        self.assertIsNone(metrics.effect_pair(items[1], items[1]['systems']['S3'], 'S1'))

    def test_effect_unit_groups_keep_separate_bootstrap_intervals(self):
        items = self.items()
        for i, item in enumerate(items):
            item['systems']['S3']['covered'] = True
            if i < 3:
                item['effect'] = {'value': i, 'unit': 'SD', 'type': "Cohen's d"}
            else:
                item['effect'] = {'value': i, 'unit': 'years', 'type': 'mean difference'}
                item['systems']['S3']['mean_diff'] = i
        point = metrics.summarize(items, 'S3', 20, 123)
        self.assertIsNone(point['effect_spearman'])
        self.assertEqual(len(point['effect_groups']), 2)
        self.assertTrue(all('interval' in group for group in point['effect_groups']))
        self.assertTrue(any(group['interval']['valid_resamples'] for group in point['effect_groups']))

    def test_breakdowns_tiers_and_calibration(self):
        items = self.items()
        items[5]['tier'] = 'T2'
        items[5]['significant'] = None
        items[1]['post_cutoff'] = True
        items[1]['systems']['S4:synthetic'] = {'covered': True, 'direction': 'higher', 'confidence': 0.8}
        output = metrics.metrics(items, 10)
        self.assertEqual(set(output['tiers']), {'T1', 'T2'})
        self.assertIn('horvath', output['tiers']['T1']['clock_family'])
        self.assertIn('post_cutoff', output['tiers']['T1'])
        point = output['tiers']['T1']['overall']['all']['S4:synthetic']
        self.assertAlmostEqual(point['confidence_brier'], 0.04)
        self.assertEqual(point['confidence_denominator'], 1)
        self.assertEqual(output['tiers']['T2']['overall']['all']['S3']['direction_denominator'], 0)


class CliTests(FixtureCase):
    def inputs(self):
        meta = self.metadata()
        path = self.csv(paths.VERIFIER_ROOT / 'replay/inputs/GSE1.csv', pd.DataFrame(
            {'gsm_id': list(meta.index), 'Horvathv1': [30, 32, 33, 20, 21, 23]}))
        export = predictions.export_path('dnam', 'GSE1')
        self.json(export, {'rows': [atlas_row(predictions=receipt(path))]})
        truth = self.root / 'ground_truth.jsonl'
        truth.write_text(json.dumps(record()) + '\n')
        manual = self.json(self.root / 'manual.json', {'synthetic-1': arms()})
        return truth, manual

    def test_full_run_manual_no_network_and_report_files(self):
        truth, manual = self.inputs()
        s4 = self.root / 's4_predictions.jsonl'
        s4.write_text(json.dumps({'ground_truth_id': 'synthetic-1', 'model': 'synthetic',
                                 'direction': 'higher', 'significant': True, 'confidence': 0.9}) + '\n')
        out = self.root / 'result'
        with patch('urllib.request.urlopen', side_effect=AssertionError('network forbidden')):
            cli.main(['run', '--ground-truth', str(truth), '--out', str(out), '--groups-manual', str(manual),
                      '--cache-only', '--bootstrap', '10'])
        self.assertEqual(set(p.name for p in out.iterdir()), {'summary.json', 'manifest.json', 'records.jsonl', 'per-record.csv', 'report.md'})
        manifest = json.loads((out / 'manifest.json').read_text())
        self.assertEqual(manifest['requests_sent'], 0)
        result = json.loads((out / 'records.jsonl').read_text())
        self.assertEqual(result['systems']['S3']['status'], 'INFERENTIAL')
        self.assertIn('S4:synthetic', result['systems'])
        self.assertIn('published_direction', (out / 'per-record.csv').read_text())
        self.assertIn('## Limits', (out / 'report.md').read_text())
        self.assertNotIn('NaN', (out / 'records.jsonl').read_text())

    def test_cache_only_unresolved_is_retained(self):
        truth, _ = self.inputs()
        out = self.root / 'cache'
        with patch('urllib.request.urlopen', side_effect=AssertionError('network forbidden')):
            cli.run(truth, out, cache_only=True, bootstrap=0)
        result = json.loads((out / 'records.jsonl').read_text())
        self.assertEqual(result['groups']['status'], 'UNRESOLVED')
        self.assertEqual(result['systems']['S3']['status'], 'UNRESOLVED_GROUPS')
        self.assertEqual(len((out / 'ledger.jsonl').read_text().splitlines()), 4)

    def test_bad_inputs_and_no_overwrite(self):
        truth, manual = self.inputs()
        with self.assertRaises(OSError):
            cli.run(truth, self.root, manual)
        truth.write_text(json.dumps(record()) + '\n' + json.dumps(record()) + '\n')
        with self.assertRaises(ValueError):
            cli.run(truth, self.root / 'bad')
        self.assertEqual(clean({'nan': np.nan, 'inf': np.inf}), {'nan': None, 'inf': None})

    def test_s4_status_duplicate_and_alias(self):
        path = self.root / 's4.jsonl'
        value = {'id': 'synthetic-1', 'model_requested': 'synthetic', 'status': 'NOT_CACHED',
                 'parsed': {'direction': 'higher'}}
        path.write_text(json.dumps(value) + '\n')
        self.assertFalse(cli.load_s4(path)[('synthetic-1', 'synthetic')]['covered'])
        path.write_text((json.dumps(value) + '\n') * 2)
        with self.assertRaises(ValueError):
            cli.load_s4(path)


if __name__ == '__main__':
    unittest.main()
