"""Focused direct-Series clock-plan integrity tests; no clocks or network calls."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

from scripts import run_recovered_dnam_clocks_2026_09_16 as runner


class DirectSeriesClockPlanTests(unittest.TestCase):
    def fixture(self, root):
        ids = ['GSM1001', 'GSM1002']
        series = root/'series.txt.gz'
        with gzip.open(series, 'wt', encoding='utf-8', newline='') as handle:
            handle.write('!Series_title\t"line one\nline two"\n')
            handle.write('!series_matrix_table_begin\n')
            handle.write('ID_REF\t' + '\t'.join(ids) + '\n')
            handle.write('cg0001\t0.2\t0.3\n')
            handle.write('!series_matrix_table_end\n')
        metadata = root/'metadata.csv'
        with metadata.open('w', newline='') as handle:
            writer = csv.writer(handle); writer.writerow(['GSM', 'title'])
            writer.writerows([[ids[0], 'a'], [ids[1], 'b']])
        origin_plan = root/'origin-plan.json'; origin_plan.write_text('{}\n')
        origin_summary = root/'origin-summary.json'; origin_summary.write_text('{}\n')
        origin_table = root/'origin-table.csv'; origin_table.write_text('study,status\n')
        reference = root/'reference.txt.gz'; reference.write_bytes(b'reference-source')
        cache = root/'GSE40279-means-fixture-geo_stream_round_trip_v2.csv'
        cache.write_text('id,mean\ncg0001,0.25\n')
        reference_receipt = root/'reference-receipt.json'
        reference_receipt.write_text(json.dumps({
            'status': 'PASS', 'reader_contract': 'geo_stream_round_trip_v2',
            'source_record': {'source_sha256': runner.sha(reference),
                              'means_sha256': runner.sha(cache)},
            'versioned_cache_path': str(cache),
            'versioned_cache_sha256': runner.sha(cache)}) + '\n')
        def record(study, status):
            readback = root/(study + '-source.json')
            beta = {'status': 'PASS', 'reason': 'Existing beta input contract passed'}
            if status == 'STOPPED':
                beta = {'status': 'FAIL', 'reason': 'AUDIT_FAIL: measurement range conflicts with input contract'}
            payload = {'study_id': study, 'source_sha256': runner.sha(series),
                'source_samples': 2, 'source_features': 1, 'source_sample_ids': ids,
                'source_header_sha256': 'header', 'source_feature_order_sha256': 'features',
                'numeric_stream_sha256': 'numeric', 'source_parser_status': 'PASS',
                'beta_input_contract': beta}
            readback.write_text(json.dumps(payload) + '\n')
            return {'study_id': study, 'execution_status': status,
                'stop_reason': None if status == 'ELIGIBLE' else beta['reason'],
                'series_matrix': {'path': str(series), 'sha256': runner.sha(series),
                                  'size_bytes': series.stat().st_size},
                'metadata': {'path': str(metadata), 'sha256': runner.sha(metadata),
                             'size_bytes': metadata.stat().st_size},
                'source_readback': {'path': str(readback), 'sha256': runner.sha(readback)},
                'source_samples': 2, 'source_features': 1,
                'ordered_source_sample_ids': ids,
                'source_header_sha256': 'header',
                'source_feature_order_sha256': 'features',
                'numeric_stream_sha256': 'numeric',
                'numeric_counts': {'finite': 2, 'outside_beta_range': 0},
                'source_parser_status': 'PASS', 'beta_input_contract': beta,
                'saved_prediction_identity_status': 'FAIL',
                'prior_prediction_use': 'PROHIBITED'}

        studies = [record('GSE900001', 'ELIGIBLE'), record('GSE900002', 'STOPPED')]
        independent_source = root/'independent-source.json'
        independent_source.write_text(json.dumps({
            'source_only': True, 'clock_executions': 0, 'model_api_calls': 0,
            'records': [{
                'role': 'series_matrix', 'study_id': entry['study_id'],
                'sha256_before': entry['series_matrix']['sha256'],
                'sha256_after': entry['series_matrix']['sha256'],
                'source_mutated_during_read': False,
                'sample_ids': entry['ordered_source_sample_ids'],
                'sample_id_count': entry['source_samples'],
                'feature_count': entry['source_features'],
                'header_identity_sha256': entry['source_header_sha256'],
                'feature_order_sha256': entry['source_feature_order_sha256'],
                'value_stream_sha256': entry['numeric_stream_sha256'],
                'finite_value_count': entry['numeric_counts']['finite'],
                'missing_value_count': 0,
                'out_of_beta_count': 0,
                'source_status': entry['beta_input_contract']['status'],
                'arithmetic_status': 'PASS', 'scientific_status': 'CANNOT VERIFY'}
                for entry in studies]}) + '\n')
        required_checks = {key: True for key in (
            'source_sha256', 'sample_count', 'feature_count', 'sample_ids',
            'full_header_sha256', 'feature_order_sha256', 'numeric_stream_sha256',
            'finite_count', 'missing_count', 'out_of_beta_count',
            'source_contract_status')}
        source_partition = {'total': 2, 'PASS': 1, 'FAIL': 1, 'CANNOT VERIFY': 0}
        independent = root/'independent.json'
        independent.write_text(json.dumps({
            'schema': 'fixture-source-cross-readback-v2', 'status': 'PASS',
            'status_scope': 'saved-source-cross-readback-only', 'source_only': True,
            'clock_executions': 0, 'model_calls': 0, 'model_api_calls': 0,
            'scientific_status': 'CANNOT VERIFY', 'failures': [],
            'source_status_partition': source_partition,
            'plan_path': str(origin_plan), 'plan_sha256': runner.sha(origin_plan),
            'root_summary_path': str(origin_summary),
            'root_summary_sha256': runner.sha(origin_summary),
            'run_readback_path': str(independent_source),
            'run_readback_sha256': runner.sha(independent_source),
            'dnam_study_count': len(studies),
            'dnam_comparisons': [{
                'study_id': entry['study_id'], 'checks': required_checks,
                'all_checks_pass': True,
                'ours_source_status': entry['beta_input_contract']['status'],
                'root_contract_status': entry['beta_input_contract']['status']}
                for entry in studies]}) + '\n')
        plan = {'schema': runner.DIRECT_PLAN_SCHEMA, 'status': 'READY_FOR_EXECUTION',
            'source_origin_plan': {'path': str(origin_plan), 'sha256': runner.sha(origin_plan)},
            'source_origin_summary': {'path': str(origin_summary), 'sha256': runner.sha(origin_summary)},
            'source_origin_validation_table': {'path': str(origin_table), 'sha256': runner.sha(origin_table)},
            'independent_source_readback': {
                'status': 'PASS', 'schema': 'fixture-source-cross-readback-v2',
                'status_scope': 'saved-source-cross-readback-only',
                'path': str(independent), 'sha256': runner.sha(independent),
                'source_run': {'path': str(independent_source),
                               'sha256': runner.sha(independent_source)},
                'source_origin_plan_sha256': runner.sha(origin_plan),
                'root_summary_sha256': runner.sha(origin_summary),
                'dnam_study_count': len(studies),
                'required_checks': list(required_checks),
                'source_status_partition_retained': source_partition,
                'scientific_status': 'CANNOT VERIFY'},
            'reference': {'source': str(reference), 'source_sha256': runner.sha(reference),
                          'reader_contract': 'geo_stream_round_trip_v2',
                          'cache_path': str(cache), 'cache_sha256': runner.sha(cache),
                          'receipt_path': str(reference_receipt),
                          'receipt_sha256': runner.sha(reference_receipt),
                          'receipt_status': 'PASS'},
            'clock_execution': {'clocks': ['Horvathv1', 'Hannum', 'PhenoAge'],
                                'coefficient_arithmetic_tolerance': {'rtol': 1e-10, 'atol': 1e-10},
                                'historical_prediction_join': False,
                                'historical_prediction_comparison': False},
            'resource_limits': {'study_concurrency': 1, 'threads_per_process': 1,
                                'automatic_retries': 0, 'network_calls': 0,
                                'paid_model_api_calls': 0},
            'denominator': {'selected_studies': 2, 'eligible_studies': 1,
                            'stopped_studies': 1, 'eligible_source_samples': 2,
                            'stopped_source_samples': 2},
            'studies': studies}
        plan_path = root/'plan.json'; plan_path.write_text(json.dumps(plan, indent=2) + '\n')
        return plan_path, plan, series, metadata, reference, cache

    def args(self, plan, reference, cache, output, study='GSE900001', prepare=True,
             historical=None):
        return argparse.Namespace(recovery=None, direct_series_plan=plan, study=study,
            reference=reference, reference_cache=cache.parent, historical=historical,
            output=output, prepare=prepare)

    def test_prepare_validates_bindings_without_output_or_clocks(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); plan, _, _, _, reference, cache=self.fixture(root)
            output=root/'output'
            result=runner.run_direct_series(self.args(plan,reference,cache,output))
            self.assertEqual(result['execution_status'],'ELIGIBLE')
            self.assertEqual(result['source_samples'],2)
            self.assertFalse(output.exists())

    def test_stopped_study_is_terminal_nonexecution_record(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); plan, _, _, _, reference, cache=self.fixture(root)
            output=root/'output'
            result=runner.run_direct_series(
                self.args(plan,reference,cache,output,study='GSE900002',prepare=False)
            )
            self.assertEqual(result['status'],'STOPPED')
            self.assertEqual(result['candidate_status'],'NOT_EXECUTED')
            self.assertFalse(result['execution_started'])
            self.assertFalse(output.exists())

    def test_hash_metadata_order_and_historical_use_are_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); plan_path, plan, _, metadata, reference, cache=self.fixture(root)
            output=root/'output'; metadata.write_text('GSM,title\nGSM1001,a\nGSM1002,changed\n')
            with self.assertRaisesRegex(ValueError,'bound file changed'):
                runner.run_direct_series(self.args(plan_path,reference,cache,output))
            self.assertFalse(output.exists())

            metadata.write_text('GSM,title\nGSM1002,b\nGSM1001,a\n')
            for entry in plan['studies']:
                entry['metadata']['sha256']=runner.sha(metadata)
                entry['metadata']['size_bytes']=metadata.stat().st_size
            plan_path.write_text(json.dumps(plan,indent=2)+'\n')
            binding=runner.load_direct_series_plan(plan_path,'GSE900001')
            self.assertTrue(binding['metadata_reindexed_to_primary_header'])
            reordered=root/'reordered-snapshot';reordered.mkdir()
            manifest=runner.snapshot_direct_sources(binding,reordered)
            self.assertTrue(manifest['metadata_reindexed_to_primary_header'])
            with Path(manifest['metadata_primary_header_order_view']['path']).open(newline='') as handle:
                self.assertEqual([row['GSM'] for row in csv.DictReader(handle)],
                                 ['GSM1001','GSM1002'])
            self.assertFalse(output.exists())

            metadata.write_text('GSM,title\nGSM1002,b\nGSM1003,c\n')
            for entry in plan['studies']:
                entry['metadata']['sha256']=runner.sha(metadata)
                entry['metadata']['size_bytes']=metadata.stat().st_size
            plan_path.write_text(json.dumps(plan,indent=2)+'\n')
            with self.assertRaisesRegex(ValueError,'metadata membership differs'):
                runner.run_direct_series(self.args(plan_path,reference,cache,output))
            self.assertFalse(output.exists())

            with self.assertRaisesRegex(ValueError,'historical is prohibited'):
                runner.run_direct_series(
                    self.args(plan_path,reference,cache,output,historical=root/'old.csv')
                )
            self.assertFalse(output.exists())

    def test_snapshot_is_byte_identical_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); plan, _, series, metadata, _, _=self.fixture(root)
            binding=runner.load_direct_series_plan(plan,'GSE900001',for_execution=True)
            output=root/'snapshot';output.mkdir()
            result=runner.snapshot_direct_sources(binding,output)
            self.assertEqual(Path(result['files']['series_matrix']['path']).read_bytes(),
                             series.read_bytes())
            self.assertEqual(Path(result['files']['metadata']['path']).read_bytes(),
                             metadata.read_bytes())
            with Path(result['metadata_primary_header_order_view']['path']).open(newline='') as handle:
                self.assertEqual([row['GSM'] for row in csv.DictReader(handle)],
                                 ['GSM1001','GSM1002'])
            with self.assertRaises(FileExistsError):
                runner.snapshot_direct_sources(binding,output)

    def test_scoped_pass_rejects_stale_root_and_mismatched_source_bindings(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); plan_path, plan, _, _, _, _=self.fixture(root)
            receipt_path=Path(plan['independent_source_readback']['path'])
            receipt=json.loads(receipt_path.read_text())
            receipt['plan_sha256']='0'*64
            receipt_path.write_text(json.dumps(receipt)+'\n')
            plan['independent_source_readback']['sha256']=runner.sha(receipt_path)
            plan_path.write_text(json.dumps(plan,indent=2)+'\n')
            with self.assertRaisesRegex(ValueError,'stale root bindings'):
                runner.load_direct_series_plan(plan_path,'GSE900001')

        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); plan_path, plan, _, _, _, _=self.fixture(root)
            source_path=Path(plan['independent_source_readback']['source_run']['path'])
            source=json.loads(source_path.read_text())
            source['records'][0]['value_stream_sha256']='0'*64
            source_path.write_text(json.dumps(source)+'\n')
            source_sha=runner.sha(source_path)
            plan['independent_source_readback']['source_run']['sha256']=source_sha
            receipt_path=Path(plan['independent_source_readback']['path'])
            receipt=json.loads(receipt_path.read_text())
            receipt['run_readback_sha256']=source_sha
            receipt_path.write_text(json.dumps(receipt)+'\n')
            plan['independent_source_readback']['sha256']=runner.sha(receipt_path)
            plan_path.write_text(json.dumps(plan,indent=2)+'\n')
            with self.assertRaisesRegex(ValueError,'values differ for GSE900001'):
                runner.load_direct_series_plan(plan_path,'GSE900001')

    def test_pending_independent_readback_allows_prepare_but_blocks_execution(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); plan_path, plan, _, _, reference, cache=self.fixture(root)
            plan['status']='HOLD_PENDING_INDEPENDENT_SOURCE_READBACK'
            plan['independent_source_readback']={'status':'PENDING'}
            plan_path.write_text(json.dumps(plan,indent=2)+'\n')
            prepared=runner.run_direct_series(self.args(plan_path,reference,cache,root/'out'))
            self.assertTrue(prepared['prepared'])
            with self.assertRaisesRegex(ValueError,'pending independent'):
                runner.run_direct_series(
                    self.args(plan_path,reference,cache,root/'out',prepare=False)
                )
            self.assertFalse((root/'out').exists())

    def test_recovery_cli_shape_remains_available(self):
        args=runner.parse_args(['--recovery','recovery','--reference','reference',
                                '--output','output','--prepare'])
        self.assertEqual(args.recovery,Path('recovery'))
        self.assertIsNone(args.direct_series_plan)
        self.assertIsNone(args.study)

    def test_direct_summary_is_labeled_as_numerical_candidate(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder); predictions=root/'predictions.csv'
            predictions.write_text('sample_id,Hannum\nGSM1001,10\nGSM1002,20\n')
            receipt={'path':str(predictions),'clocks':['Hannum'],
                     'checks':[{'clock':'Hannum','required_features':71,
                                'observed_features_min':60,'numerical_status':'PASS',
                                'max_abs_independent_difference':0.0}]}
            recovery={'source_samples':2,'study_id':'GSE900001'}
            result=runner.summarize(receipt,recovery,None,root,source_mode='direct_series')
            self.assertEqual(result['candidate_status'],'NUMERICAL_CANDIDATE')
            self.assertEqual(result['measurement_unit_applicability'],'CANNOT VERIFY')
            self.assertEqual(result['biological_sample_unit_status'],'CANNOT VERIFY')
            self.assertNotIn('max_absolute_change',result['rows'][0])


if __name__ == '__main__':
    unittest.main()
