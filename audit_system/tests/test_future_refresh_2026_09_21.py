"""Software fixtures for one-config refresh; no empirical or model performance claims."""
import copy
import gzip
import io
import json
import shutil
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import h5py

from src import update_pipeline_2026_09_20 as update
from audit_system import archs4_update_2026_09_10 as intake
from audit_system.tests import test_archs4_batch_extraction_2026_09_21 as batch_tests
from audit_system.tests import require_development_artifact


class FutureRefreshTests(unittest.TestCase):
    def setUp(self):
        self.fixture = batch_tests.Archs4BatchExtractionTest()
        self.fixture.setUp(); self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root
        self.output = self.root / 'refresh-output'
        self.calls = []
        self.child_calls = []
        self.ids = ['1', '2', '3']
        self.query_total = 3
        self.seen = self.root / 'saved-analyses.jsonl'; self.seen.write_text('')
        self.family = self.root / 'sample-sources.json'
        # A future metadata receipt must not count as prior analysis membership.
        self.family.write_text(json.dumps([{'accession': 'GSE1', 'sample_ids': ['GSM1', 'GSM2'], 'source_kind': 'geo_series_receipt'}]))
        self.plan = {'schema': update.REFRESH_SCHEMA, 'revision': 'software-fixture-v1',
            'modality': 'human_rna', 'query': 'SOFTWARE_FIXTURE', 'date_from': '2025-01-01', 'date_to': '2025-12-31',
            'archs4': {'path': str(self.fixture.release), 'sha256': update.sha(self.fixture.release)}, 'max_release_bytes': 1000000,
            'seen_worklist': {'path': str(self.seen), 'sha256': update.sha(self.seen)},
            'seen_sample_sources': {'path': str(self.family), 'sha256': update.sha(self.family)},
            'discovery_limits': {'max_requests': 20, 'max_response_bytes': 100000, 'max_total_response_bytes': 1000000,
                'max_seconds': 30, 'max_query_ids': 10, 'max_candidates': 3, 'minimum_interval_seconds': .35},
            'network_budget': {'max_requests': 60, 'max_bytes': 3000000},
            'extraction_limits': {'max_studies': 3, 'max_samples_per_study': 4, 'max_extracted_bytes_per_study': 100000,
                'max_total_extracted_bytes': 300000, 'concurrency': 1},
            'update_limits': {'max_studies': 3, 'max_input_bytes': 20000000, 'max_samples_per_study': 4,
                'max_child_seconds': 60, 'max_total_model_calls': 30, 'max_total_cost_usd': 10., 'concurrency': 1, 'allow_network': False},
            'pipeline': {'max_calls': 2, 'max_cost_usd': 1.},
            'comprehensive_audit': {'model': 'gpt-5.6-luna', 'max_calls': 2, 'max_cost_usd': 1.,
                'input_usd_per_million': .5, 'output_usd_per_million': 1.8}}
        self.path = self.root / 'refresh.json'
        self.write_plan()
        self.addCleanup(patch.stopall)
        patch('urllib.request.urlopen', side_effect=AssertionError('Network forbidden in software fixture')).start()
        patch.object(update.time, 'sleep').start()

    def write_plan(self, path=None):
        (path or self.path).write_text(json.dumps(self.plan, sort_keys=True))

    def transport(self, url, cap):
        self.calls.append(url)
        if 'esearch.fcgi' in url:
            return json.dumps({'esearchresult': {'idlist': self.ids, 'count': str(self.query_total)}}).encode()
        if 'esummary.fcgi' in url:
            ids = update.urllib.parse.parse_qs(update.urllib.parse.urlparse(url).query)['id'][0].split(',')
            return json.dumps({'result': {'uids': ids, **{uid: {'accession': 'GSE' + uid, 'pdat': '2025/01/01'} for uid in ids}}}).encode()
        query = update.urllib.parse.parse_qs(update.urllib.parse.urlparse(url).query)
        role = 'series' if query['targ'] == ['self'] else 'samples'
        return Path(self.fixture.sources[query['acc'][0]][role]['path']).read_bytes()

    def software_child(self, request_path, log, timeout):
        """Consumer boundary fixture; it does not validate clocks or science."""
        request = json.loads(request_path.read_text()); self.child_calls.append(request)
        version = update.read_sealed(Path(request['version_path']))
        update.verify_sources(version['sources'])
        if request['phase'] == 'pipeline':
            directory = Path(request['pipeline_root']); directory.mkdir(parents=True, exist_ok=True)
            (directory / ('software-' + request['until'] + '.json')).write_text(json.dumps({'software_fixture': True, 'study': version['study_id']}))
            result = {'software_fixture': True, 'scientific_status': 'CANNOT VERIFY'}
        else:
            directory = Path(request['audit_request']['output']); directory.mkdir(parents=True)
            (directory / 'software-review.json').write_text('{"software_fixture":true,"scientific_status":"CANNOT VERIFY"}')
            result = {'execution_status': 'COMPLETED', 'scientific_status': 'CANNOT VERIFY', 'software_fixture': True}
        Path(request['result_path']).write_text(json.dumps(result))
        return 0

    def execute(self, until='report', executor=None):
        self.assertTrue((update.ROOT / update.REFRESH_VERIFIER).is_file())
        return update.refresh(self.path, self.output, until, transport=self.transport, executor=executor or self.software_child)

    def run_root(self):
        return self.output / 'refreshes' / update.identity(self.plan)

    def test_one_config_source_graph_and_pipeline_boundary_retains_all_studies(self):
        result = self.execute()
        self.assertEqual(result['enumerated_uids'], 3)
        self.assertEqual(result['dispositions'], {'AUDIT_EXECUTION_COMPLETED': 2, 'SOURCE_PROFILE_UNRESOLVED': 1})
        self.assertEqual(len(self.calls), 8)
        self.assertEqual(len(self.child_calls), 4)
        self.assertEqual(result['reserved_cost_usd'], 4.)
        self.assertEqual(result['reserved_model_calls'], 8)
        self.assertTrue(all(r['scientific_status'] == 'CANNOT VERIFY' for r in result['studies']))
        handoff = update.read_sealed(self.run_root() / 'phases/handoff/handoff/handoff.json')
        self.assertEqual(handoff['counts']['candidate_cells_compared'], 6)
        self.assertEqual(handoff['counts']['enrolled_studies'], 3)
        self.assertEqual(len(list(self.run_root().glob('phases/*/plan.json'))), 4)
        seen = update.read_sealed(self.run_root() / 'phases/seen-membership/result.json')
        self.assertEqual(seen['membership_categories'], {})
        for request in self.child_calls:
            version = update.read_sealed(Path(request['version_path']))
            self.assertIn('matrix_provenance', version['sources'])
            self.assertNotIn('archs4', version['sources'])

    def test_completed_resume_has_zero_network_or_child_calls_and_rejects_tamper(self):
        original = self.execute()
        with patch.object(update, 'run_child', side_effect=AssertionError('No child on resume')):
            result = update.refresh(self.path, self.output, transport=lambda *x: self.fail('Network on resume'),
                                    executor=lambda *x: self.fail('Child on resume'))
        self.assertEqual(result, original)
        matrix = next((self.run_root() / 'phases/count-extraction').rglob('counts.csv'))
        matrix.write_text(matrix.read_text() + 'SOFTWARE_CORRUPTION\n')
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'phase output changed'):
            self.execute()

    def test_unsupported_audit_model_stops_before_discovery_or_registry(self):
        self.plan['comprehensive_audit']['model']='unsupported-model-software-fixture'
        self.write_plan()
        with self.assertRaisesRegex(ValueError,'unsupported model rejected before execution'):
            update.refresh(self.path,self.output,transport=self.transport,executor=self.software_child)
        self.assertEqual(self.calls,[])
        self.assertEqual(self.child_calls,[])
        self.assertFalse(self.output.exists())

    def test_input_then_report_reuses_sources_and_reserves_models_only_when_needed(self):
        first = self.execute('input')
        self.assertEqual(first['reserved_model_calls'], 0)
        network_count = len(self.calls)
        second = self.execute('report')
        self.assertEqual(len(self.calls), network_count)
        self.assertEqual(second['reserved_model_calls'], 8)
        self.assertEqual({r['version_id'] for r in (v['pipeline'] for v in first['studies'] if v['pipeline'])},
                         {r['version_id'] for r in (v['pipeline'] for v in second['studies'] if v['pipeline'])})
        self.assertEqual(self.execute('input'), first)

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_actual_input_stage_consumes_verified_software_counts_without_network_or_models(self):
        def input_child(request_path, log, timeout):
            request = json.loads(request_path.read_text())
            self.assertEqual(request['until'], 'input')
            self.assertEqual(request['phase'], 'pipeline')
            update.child(request)
            return 0
        result = self.execute('input', input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 2, 'SOURCE_PROFILE_UNRESOLVED': 1})
        self.assertEqual(result['reserved_model_calls'], 0)

    def test_new_revision_keeps_lifetime_reservations_and_attempted_membership(self):
        first = self.execute()
        first_path = self.root / 'first.json'; first_path.write_bytes(self.path.read_bytes())
        self.plan['revision'] = 'software-fixture-v2'; self.write_plan()
        second = self.execute()
        self.assertEqual(second['reserved_network_requests'], 40)
        self.assertEqual(second['reserved_cost_usd'], 4.)
        self.assertEqual(second['dispositions']['SEEN_ACCESSION'], 2)
        self.assertEqual(len(self.child_calls), 4)
        memberships = [json.loads(line) for line in (self.run_root() / 'phases/seen-membership/seen/studies.jsonl').read_text().splitlines()]
        self.assertEqual({row['membership'] for row in memberships}, {'ANALYSIS_EXECUTED'})
        self.plan['revision'] = 'software-fixture-v3'; self.write_plan(); self.execute()
        self.plan['revision'] = 'software-fixture-v4'; self.write_plan()
        call_count = len(self.calls)
        with self.assertRaisesRegex(RuntimeError, 'Lifetime discovery reservation'):
            self.execute()
        self.assertEqual(len(self.calls), call_count)
        self.path.write_bytes(first_path.read_bytes()); self.plan = json.loads(self.path.read_text())
        self.assertEqual(self.execute(), first)

    def test_seen_execution_evidence_remains_bound_when_later_refresh_resumes(self):
        self.execute()
        self.plan['revision'] = 'software-fixture-v2'; self.write_plan()
        self.execute()
        receipt = update.read_sealed(self.run_root() / 'phases/seen-membership/phase.json')
        outputs = [Path(path) for path in receipt['external_sha256'] if Path(path).name == 'software-report.json']
        self.assertEqual(len(outputs), 2)
        outputs[0].write_text('{"software_fixture":"changed"}')
        request_count, child_count = len(self.calls), len(self.child_calls)
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'Refresh external output changed'):
            self.execute()
        self.assertEqual((len(self.calls), len(self.child_calls)), (request_count, child_count))

    def test_budget_limits_cannot_be_reset_in_existing_output(self):
        self.execute('input')
        self.plan['revision'] = 'changed-budget'; self.plan['network_budget']['max_requests'] += 1; self.write_plan()
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'lifetime limits changed'):
            self.execute()

    def test_catalogue_only_seen_declaration_is_not_prior_analysis(self):
        self.seen.write_text('{"study_id":"GSE1","source_kind":"geo_series_receipt"}\n')
        self.plan['seen_worklist']['sha256'] = update.sha(self.seen); self.write_plan()
        result = self.execute()
        self.assertIn('seen-membership', result['failures'])
        self.assertEqual(self.calls, [])
        self.assertIsNone(result['query_total_uids'])

    def test_no_candidates_and_unenumerated_query_have_distinct_denominators(self):
        self.ids = []; self.query_total = 0
        empty = self.execute()
        self.assertEqual(empty['status'], 'COMPLETED')
        self.assertEqual(empty['query_total_uids'], 0)
        self.assertFalse(self.child_calls)
        self.ids = ['1']; self.query_total = 9
        self.plan['revision'] = 'incomplete-query'; self.write_plan()
        limited = self.execute('input')
        self.assertEqual(limited['unretrieved_uid_count'], 8)
        self.assertEqual(limited['status'], 'COMPLETED_WITH_UNRESOLVED_CASES')

    def test_dnam_without_reference_retains_source_ledger_without_assay_or_models(self):
        self.plan['modality'] = 'human_dnam'; del self.plan['archs4']; del self.plan['max_release_bytes']; self.write_plan()
        result = self.execute()
        self.assertEqual(result['dispositions'], {'DNAM_REFERENCE_UNRESOLVED': 3})
        self.assertEqual(len(self.calls), 5)
        self.assertEqual(self.child_calls, [])
        self.assertEqual(result['reserved_cost_usd'], 0)
        self.assertNotIn('source-profile', result['phase_receipts'])

    def test_incomplete_discovery_does_not_reset_request_reservation_on_resume(self):
        self.plan['discovery_limits']['max_requests'] = 3; self.write_plan()
        result = self.execute()
        self.assertEqual(result['enumerated_uids'], 3)
        self.assertEqual(len(self.calls), 3)
        self.assertEqual(result['reserved_network_requests'], 3)
        self.assertEqual(self.execute(), result)
        self.assertEqual(len(self.calls), 3)
        self.assertEqual(self.child_calls, [])

    def test_changed_origin_and_changed_phase_inventory_are_rejected(self):
        self.execute('input')
        extra = self.run_root() / 'phases/source-profile/unbound.json'; extra.write_text('{}')
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'phase inventory changed'):
            self.execute('input')
        extra.unlink()
        with self.fixture.release.open('ab') as handle: handle.write(b'CHANGED SOFTWARE FIXTURE')
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'Source SHA256 differs'):
            self.execute('input')

    def test_absent_release_study_is_retained_while_other_inputs_continue(self):
        with h5py.File(self.fixture.release, 'a') as handle:
            handle['meta/samples/series_id'][2] = b'GSE9'
        self.plan['archs4']['sha256'] = update.sha(self.fixture.release); self.write_plan()
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1, 'ASSAY_INPUT_UNRESOLVED': 1, 'SOURCE_PROFILE_UNRESOLVED': 1})
        full_profile = update.read_sealed(self.run_root() / 'phases/source-profile/source-profiles.json')
        selected = update.read_sealed(self.run_root() / 'phases/extraction-subset/source-profiles.json')
        self.assertEqual(len(full_profile['studies']), 3)
        self.assertEqual([r['study_id'] for r in selected['studies']], ['GSE1', 'GSE3'])
        self.assertIn('GSE2', selected['excluded_from_extraction'])
        self.assertEqual(selected['verification_script_sha256'], update.sha(update.__file__))
        self.assertEqual(self.execute('input'), result)

    def test_oversized_source_study_does_not_change_its_group_or_block_others(self):
        self.plan['extraction_limits']['max_samples_per_study'] = 1; self.write_plan()
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'SOURCE_SAMPLE_LIMIT_UNRESOLVED': 1, 'PIPELINE_COMPLETED_THROUGH_INPUT': 1, 'SOURCE_PROFILE_UNRESOLVED': 1})
        original = update.read_sealed(self.run_root() / 'phases/source-profile/source-profiles.json')
        self.assertEqual(original['studies'][0]['sample_ids'], ['GSM1', 'GSM2'])
        self.assertEqual(len(self.child_calls), 1)

    def test_source_readback_failure_blocks_downstream_and_retains_all_cases(self):
        original = intake.run_archs4_batch_extraction
        def corrupted_software_batch(path):
            result = original(path)
            matrix = next((self.run_root() / 'phases/count-extraction').rglob('counts.csv'))
            matrix.write_text(matrix.read_text() + 'SOFTWARE_CORRUPTION\n')
            return result
        with patch.object(intake, 'run_archs4_batch_extraction', side_effect=corrupted_software_batch):
            result = self.execute()
        self.assertIn('independent-source-readback', result['failures'])
        self.assertEqual(result['dispositions'], {'INDEPENDENT_SOURCE_VERIFICATION_UNRESOLVED': 2, 'SOURCE_PROFILE_UNRESOLVED': 1})
        self.assertEqual(self.child_calls, [])
        self.assertEqual(result['reserved_cost_usd'], 0)
        self.assertEqual(self.execute(), result)

    def test_failed_producer_does_not_hide_case_or_retry_and_others_continue(self):
        def first_fails(path, log, timeout):
            request = json.loads(path.read_text())
            version = update.read_sealed(Path(request['version_path']))
            if version['study_id'] == 'GSE1':
                self.child_calls.append(request)
                return 1
            return self.software_child(path, log, timeout)
        result = self.execute(executor=first_fails)
        self.assertEqual(result['dispositions'], {'FAILED': 1, 'AUDIT_EXECUTION_COMPLETED': 1, 'SOURCE_PROFILE_UNRESOLVED': 1})
        self.assertEqual(result['reserved_cost_usd'], 3.)
        count = len(self.child_calls)
        self.assertEqual(self.execute(executor=first_fails), result)
        self.assertEqual(len(self.child_calls), count)

    def test_interrupted_phase_is_preserved_without_an_implicit_retry(self):
        directory = self.root / 'interrupted/phases/discovery'; directory.mkdir(parents=True)
        inputs = {'software_fixture': '0' * 64}
        update.seal(directory / 'started.json', {'phase': 'discovery', 'input_sha256': inputs})
        callback = lambda *args: self.fail('Interrupted phase repeated')
        first = update._refresh_phase(self.root / 'interrupted', 'discovery', inputs, callback)
        self.assertEqual(first['status'], 'UNRESOLVED')
        self.assertIn('Interrupted phase retained', first['error'])
        self.assertEqual(update._refresh_phase(self.root / 'interrupted', 'discovery', inputs, callback), first)


class MouseFutureRefreshTests(unittest.TestCase):
    """The same bounded source graph must retain mouse identity end to end."""

    setUp = FutureRefreshTests.setUp
    write_plan = FutureRefreshTests.write_plan
    transport = FutureRefreshTests.transport
    software_child = FutureRefreshTests.software_child
    execute = FutureRefreshTests.execute
    run_root = FutureRefreshTests.run_root

    def mouse_sources(self):
        self.fixture.mouse_fixture()
        self.plan['modality'] = 'mouse_rna'
        self.plan['archs4'] = {'path': str(self.fixture.release),
                               'sha256': update.sha(self.fixture.release)}
        self.write_plan()

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_mouse_input_binds_species_and_preserves_profile_stop(self):
        self.mouse_sources()

        def input_child(request_path, log, timeout):
            request = json.loads(request_path.read_text())
            self.assertEqual(request['until'], 'input')
            update.child(request)
            return 0

        result = self.execute('input', input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 2,
                                                  'SOURCE_PROFILE_UNRESOLVED': 1})
        self.assertEqual(result['reserved_model_calls'], 0)
        base = self.run_root() / 'phases'
        profile = update.read_sealed(base / 'source-profile/source-profiles.json')
        batch = update.read_sealed(base / 'count-extraction/batch/batch-summary.json')
        checked = update.read_sealed(base / 'independent-source-readback/result.json')
        handoff = update.read_sealed(base / 'handoff/handoff/handoff.json')
        self.assertEqual([document['modality'] for document in (profile, batch, checked, handoff)],
                         ['mouse_rna'] * 4)
        self.assertEqual([document['species'] for document in (checked, handoff)], ['mouse'] * 2)
        self.assertEqual(batch['source_hdf5']['species'], 'mouse')
        self.assertEqual(handoff['counts']['candidate_cells_compared'], 6)
        self.assertEqual(handoff['counts']['candidate_studies'], 2)
        self.assertEqual(self.execute('input', lambda *args: self.fail('child on resume')), result)

    def test_human_release_cannot_enter_mouse_refresh_with_valid_hash(self):
        self.plan['modality'] = 'mouse_rna'
        self.write_plan()
        result = self.execute('input')
        self.assertEqual(result['reserved_model_calls'], 0)
        self.assertEqual(result['dispositions'], {'ASSAY_INPUT_UNRESOLVED': 3})
        self.assertFalse((self.run_root() / 'phases/count-extraction/batch/batch-summary.json').exists())


class ArchivedHumanMatrixReplayTests(unittest.TestCase):
    """The pre-species human receipt remains valid; a foreign species does not."""

    def test_archived_human_wrapper_and_cross_species_rejection(self):
        wrapper_path = (Path(__file__).resolve().parents[2] /
            'paper/revision_2026-09-07/evidence/continuous-audit-pipeline-2026-09-20/'
            'update-controller/archs4-batch-handoff-001/wrappers/human_rna-GSE127814.json')
        require_development_artifact(wrapper_path)
        wrapper = json.loads(wrapper_path.read_text())
        self.assertNotIn('species', wrapper['source_hdf5'])
        sources = {role: wrapper[role] for role in ('matrix', 'series', 'samples')}
        # The receipt names the matrices it was built from by their paths on the machine that ran the batch; they are
        # development artefacts and are not part of the code release or the rebuild bundle.
        for role, source in sources.items():
            location = source.get('path') if isinstance(source, dict) else source
            if location and not Path(str(location)).exists():
                self.skipTest(f'archived {role} source of the receipt is not available here: {location}')
        result = update.read_matrix_provenance(wrapper_path, update.sha(wrapper_path),
                                               sources, 'GSE127814', 'human_rna')
        self.assertEqual(result['measurement']['cells_compared'], 403116)
        self.assertFalse(result['measurement']['hdf5_rehashed_for_this_child'])
        with tempfile.TemporaryDirectory() as name:
            changed = Path(name) / 'foreign-species.json'
            wrapper['source_hdf5']['species'] = 'mouse'
            changed.write_text(json.dumps(wrapper))
            with self.assertRaisesRegex(update.UpdateIntegrityError, 'HDF5 origin scope'):
                update.read_matrix_provenance(changed, update.sha(changed), sources,
                                              'GSE127814', 'human_rna')


class DnamFutureRefreshTests(unittest.TestCase):
    """Source-bound software fixtures; they do not validate clock applicability."""

    def setUp(self):
        FutureRefreshTests.setUp(self)
        self.ids = ['286412']; self.query_total = 1
        self.study = 'GSE286412'; self.sample_ids = ['GSM111', 'GSM112']
        self.reference = self.root / 'reference.tsv'
        self.reference.write_text('SOFTWARE FIXTURE REFERENCE; INPUT STAGE DOES NOT CONSUME VALUES\n')
        del self.plan['archs4']; del self.plan['max_release_bytes']
        self.plan.update(modality='human_dnam', dnam_reference={
            'path': str(self.reference), 'sha256': update.sha(self.reference)})
        self.plan['discovery_limits'].update(max_candidates=1, max_requests=12)
        self.plan['extraction_limits']['max_studies'] = 1
        self.plan['update_limits']['max_studies'] = 1
        self.plan['network_budget'] = {'max_requests': 24, 'max_bytes': 3000000}
        self.series_beta_text = 'Values are reported as methylation beta values.'
        self.supplement_name = self.study + '_Betas_software.tsv.gz'
        self.matrix_name = self.study + '_series_matrix.txt.gz'
        self.empty_matrix = True
        self.matrix_value = .1
        self.extra_matrix_name = None
        self.extra_supplement_name = None
        self.sample_molecule = 'genomic DNA'
        self.sample_platforms = ['GPL13534', 'GPL13534']
        self.mismatched_membership = False
        self.matrix_study = self.study
        self.matrix_identity_row = None
        self.supplement_aliases = self.sample_ids.copy()
        self.write_plan()

    write_plan = FutureRefreshTests.write_plan
    execute = FutureRefreshTests.execute
    run_root = FutureRefreshTests.run_root
    software_child = FutureRefreshTests.software_child

    @staticmethod
    def gz(payload):
        target = io.BytesIO()
        with gzip.GzipFile(fileobj=target, mode='wb', filename='', mtime=0) as handle:
            handle.write(payload.encode())
        return target.getvalue()

    def series_payload(self):
        prefix = self.study[:-3] + 'nnn'
        supplement = f'ftp://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{self.study}/suppl/{self.supplement_name}'
        lines = [f'^SERIES = {self.study}', f'!Series_geo_accession = {self.study}',
            '!Series_status = Public on Jan 01 2025', '!Series_type = Methylation profiling by array',
            '!Series_platform_id = GPL13534', f'!Series_data_processing = {self.series_beta_text}',
            *(f'!Series_sample_id = {sample}' for sample in self.sample_ids),
            f'!Series_supplementary_file = {supplement}']
        if self.extra_supplement_name:
            lines.append(f'!Series_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{self.study}/suppl/{self.extra_supplement_name}')
        return ('\n'.join(lines) + '\n').encode()

    def samples_payload(self):
        rows = []
        ids = self.sample_ids if not self.mismatched_membership else [self.sample_ids[0], 'GSM99']
        for number, (sample, platform) in enumerate(zip(ids, self.sample_platforms), 1):
            rows.extend([f'^SAMPLE = {sample}', f'!Sample_title = Software sample [{sample}]',
                '!Sample_source_name_ch1 = software source', '!Sample_organism_ch1 = Homo sapiens',
                f'!Sample_molecule_ch1 = {self.sample_molecule}', f'!Sample_platform_id = {platform}',
                f'!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/{sample}/suppl/{sample}_R0{number}C0{number}_Grn.idat.gz',
                f'!Sample_supplementary_file = ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSMnnn/{sample}/suppl/{sample}_R0{number}C0{number}_Red.idat.gz'])
        return ('\n'.join(rows) + '\n').encode()

    def matrix_payload(self):
        identity_row = self.matrix_identity_row or f'!Series_geo_accession = {self.matrix_study}'
        lines = [f'^SERIES = {self.matrix_study}', identity_row,
            '!Sample_geo_accession\t"GSM111"\t"GSM112"', '!series_matrix_table_begin',
            '"ID_REF"\t"GSM111"\t"GSM112"']
        if not self.empty_matrix:
            lines.append(f'cg00000001\t{self.matrix_value}\t0.2')
        lines.append('!series_matrix_table_end')
        return self.gz('\n'.join(lines) + '\n')

    def supplement_payload(self):
        return self.gz(' '.join('"' + value + '"' for value in self.supplement_aliases)
                       + '\ncg00000001 0.1 0.2\nch.1.2R 0.3 NA\n')

    def transport(self, url, cap):
        self.calls.append(url)
        if 'esearch.fcgi' in url:
            payload = json.dumps({'esearchresult': {'idlist': self.ids, 'count': str(self.query_total)}}).encode()
        elif 'esummary.fcgi' in url:
            payload = json.dumps({'result': {'uids': self.ids, self.ids[0]: {'accession': self.study, 'pdat': '2025/01/01'}}}).encode()
        elif 'targ=self' in url:
            payload = self.series_payload()
        elif 'targ=gsm' in url:
            payload = self.samples_payload()
        elif url.endswith('/matrix/'):
            names = [self.matrix_name] + ([self.extra_matrix_name] if self.extra_matrix_name else [])
            payload = ''.join(f'<a href="{name}">{name}</a>\n' for name in names).encode()
        elif url.endswith('/' + self.matrix_name):
            payload = self.matrix_payload()
        elif url.endswith('/' + self.supplement_name):
            payload = self.supplement_payload()
        else:
            raise AssertionError('Unexpected software-fixture URL: ' + url)
        if len(payload) > cap:
            return payload
        return payload

    def actual_input_child(self, request_path, log, timeout):
        request = json.loads(request_path.read_text())
        self.child_calls.append(request)
        self.assertEqual((request['phase'], request['until']), ('pipeline', 'input'))
        update.child(request)
        return 0

    def resolution(self):
        return json.loads((self.run_root() / 'phases/dnam-input-resolution/source/resolution.json').read_text())

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_empty_matrix_unique_beta_supplement_runs_existing_recovery_and_input_stage(self):
        result = self.execute('input', self.actual_input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1})
        self.assertEqual((result['reserved_model_calls'], result['reserved_cost_usd']), (0, 0))
        resolution = self.resolution(); row = resolution['studies'][0]
        self.assertEqual((resolution['ready_studies'], resolution['requests_used']), (1, 4))
        self.assertEqual(row['sample_ids'], self.sample_ids)
        self.assertEqual(row['numeric_validation']['independent_verification']['status'], 'PASS')
        self.assertEqual(row['numeric_validation']['reconstruction']['feature_count'], 2)
        self.assertEqual(set(row['sources']), {'series', 'samples', 'matrix', 'dnam_supplement', 'dnam_reference'})
        version = update.read_sealed(Path(self.child_calls[0]['version_path']))
        self.assertEqual(set(version['sources']), {'series', 'samples', 'matrix', 'dnam_supplement', 'dnam_reference'})
        molecular = json.loads((Path(self.child_calls[0]['pipeline_root']) / '02_molecular_input/complete.json').read_text())['result']
        self.assertEqual((molecular['feature_count'], molecular['sample_count'], molecular['repair_status']), (2, 2, 'PASS'))

    def xlsx_supplement_payload(self, header=None):
        """Source transport fixture; no empirical measurement or model claims."""
        from openpyxl import Workbook
        workbook = Workbook(); sheet = workbook.active
        sheet.append(header or ['ID_REF', *self.sample_ids])
        sheet.append(['cg00000001', .1, .2]); sheet.append(['ch.1.2R', .3, None])
        output = io.BytesIO(); workbook.save(output); workbook.close()
        return output.getvalue()

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_xlsx_source_supplement_completes_input_with_independent_cell_verification(self):
        self.supplement_name = self.study + '_Betas_software.xlsx'
        with patch.object(self, 'supplement_payload', side_effect=self.xlsx_supplement_payload):
            result = self.execute('input', self.actual_input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1})
        verification = self.resolution()['studies'][0]['numeric_validation']['independent_verification']
        self.assertEqual((verification['values_compared'], verification['missing_values']), (4, 1))
        self.assertIn('openpyxl', verification['parser'])
        self.assertEqual((result['reserved_model_calls'], result['reserved_cost_usd']), (0, 0))

    def test_xlsx_decompressed_limit_stops_before_reconstruction_or_child(self):
        self.supplement_name = self.study + '_Betas_software.xlsx'
        self.plan['extraction_limits']['max_extracted_bytes_per_study'] = 1000
        self.write_plan()
        with patch.object(self, 'supplement_payload', side_effect=self.xlsx_supplement_payload):
            result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('workbook exceeds decompressed input byte limit', result['studies'][0]['reason'])
        self.assertFalse(self.child_calls)
        self.assertEqual(self.resolution()['requests_used'], 4)

    def test_xlsx_aggregate_columns_stop_before_model_or_child(self):
        self.supplement_name = self.study + '_Betas_software.xlsx'
        with patch.object(self, 'supplement_payload', return_value=self.xlsx_supplement_payload(['ID_REF', 'Control', 'Treated'])):
            result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('no unique source binding', result['studies'][0]['reason'])
        self.assertFalse(self.child_calls)
        self.assertEqual(self.resolution()['requests_used'], 4)

    def test_mixed_assay_reason_reports_declared_types_without_claiming_methylation_absent(self):
        payload = self.series_payload() + b'!Series_type = Expression profiling by high throughput sequencing\n'
        with patch.object(self, 'series_payload', return_value=payload):
            result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('single declared Series assay type', result['studies'][0]['reason'])
        self.assertIn('Expression profiling by high throughput sequencing', result['studies'][0]['reason'])
        self.assertFalse(self.child_calls)

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_archived_geo_tab_quoted_study_identity_reaches_existing_recovery(self):
        # Literal identity records from GSE286412 and GSE282439 archived
        # Series Matrix files; this is a software fixture, not an assay read.
        for study in ('GSE286412', 'GSE282439'):
            with self.subTest(study=study):
                self.study = study
                self.ids = [study.removeprefix('GSE')]
                self.matrix_study = study
                self.matrix_name = study + '_series_matrix.txt.gz'
                self.supplement_name = study + '_Betas_software.tsv.gz'
                self.matrix_identity_row = f'!Series_geo_accession\t"{study}"'
                self.plan['revision'] = study + '-quoted-source-fixture'
                self.output = self.root / (study + '-quoted-output')
                self.write_plan()
                result = self.execute('input', self.actual_input_child)
                self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1})
                self.assertEqual((self.resolution()['requests_used'],
                                  self.resolution()['studies'][0]['numeric_validation']['independent_verification']['status']),
                                 (4, 'PASS'))

    def test_ambiguous_tab_study_identity_stops_with_source_receipts(self):
        self.matrix_identity_row = '!Series_geo_accession\t"GSE286412"\t"GSE286412"'
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('ambiguous GEO record', result['studies'][0]['reason'])
        self.assertFalse(self.child_calls)
        self.assertEqual(self.resolution()['requests_used'], 3)

    def test_nonempty_matrix_needs_source_beta_declaration_and_range_guard(self):
        self.empty_matrix = False; self.series_beta_text = 'Processed array signal values.'
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('lacks a literal source beta-unit declaration', result['studies'][0]['reason'])
        self.assertFalse(self.child_calls)

        self.output = self.root / 'range-output'; self.plan['revision'] = 'range-v2'; self.series_beta_text = 'Beta values.'
        self.matrix_value = 2.; self.write_plan()
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('measurement range conflicts', result['studies'][0]['reason'])
        self.assertFalse(self.child_calls)

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_source_declared_nonempty_beta_matrix_runs_input_without_supplement(self):
        self.empty_matrix = False
        result = self.execute('input', self.actual_input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1})
        self.assertEqual(self.resolution()['requests_used'], 3)
        self.assertEqual(set(self.resolution()['studies'][0]['sources']),
                         {'series', 'samples', 'matrix', 'dnam_reference'})
        self.assertEqual((result['reserved_model_calls'], result['reserved_cost_usd']), (0, 0))
        molecular = json.loads((Path(self.child_calls[0]['pipeline_root']) / '02_molecular_input/complete.json').read_text())['result']
        self.assertEqual((molecular['feature_count'], molecular['sample_count']), (1, 2))

    def test_source_profile_membership_and_unique_source_routes_fail_closed(self):
        self.mismatched_membership = True
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('membership differs', result['studies'][0]['reason'])
        self.assertEqual(len(self.calls), 4)  # Discovery plus one Sample request; no assay download.

        self.output = self.root / 'ambiguous-output'; self.plan['revision'] = 'ambiguous-v2'; self.mismatched_membership = False
        self.extra_matrix_name = self.study + '_other_series_matrix.txt.gz'; self.write_plan()
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('exactly one unambiguous', result['studies'][0]['reason'])

    def test_supplement_ambiguity_and_alias_without_source_owner_stop_before_child(self):
        self.extra_supplement_name = self.study + '_Betas_second.tsv.gz'
        first = self.execute('input')
        self.assertIn('exactly one source-declared', first['studies'][0]['reason'])
        self.assertFalse(self.child_calls)
        self.output = self.root / 'alias-output'; self.plan['revision'] = 'alias-v2'
        self.extra_supplement_name = None; self.supplement_aliases[1] = 'UNDECLARED_R02C02'; self.write_plan()
        second = self.execute('input')
        self.assertIn('no unique source binding', second['studies'][0]['reason'])
        self.assertFalse(self.child_calls)

    def test_matrix_identity_and_compressed_expansion_limits_stop_before_child(self):
        self.matrix_study = 'GSE999999'
        first = self.execute('input')
        self.assertIn('matrix study identity differs', first['studies'][0]['reason'])
        self.assertFalse(self.child_calls)
        self.output = self.root / 'expansion-output'; self.plan['revision'] = 'expansion-v2'
        self.matrix_study = self.study
        self.plan['extraction_limits']['max_extracted_bytes_per_study'] = 10; self.write_plan()
        second = self.execute('input')
        self.assertIn('decompressed input byte limit', second['studies'][0]['reason'])
        self.assertFalse(self.child_calls)

    def test_request_cap_exhaustion_retains_case_without_unrecorded_network_or_child(self):
        self.plan['discovery_limits']['max_requests'] = 5; self.write_plan()
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('Insufficient request remainder', result['studies'][0]['reason'])
        self.assertEqual(len(self.calls), 3)  # ESearch, ESummary and Series discovery only.
        self.assertFalse(self.child_calls)
        self.assertEqual(result['reserved_network_requests'], 5)

    def test_three_remaining_requests_reject_fourth_before_network(self):
        self.plan['discovery_limits']['max_requests'] = 6; self.write_plan()
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('Insufficient request remainder for source-declared DNAm beta supplement',
                      result['studies'][0]['reason'])
        self.assertEqual(len(self.calls), 6)  # Three discovery plus three DNAm requests.
        self.assertEqual(self.resolution()['requests_used'], 3)
        self.assertFalse((self.run_root() / ('phases/dnam-input-resolution/source/requests/'
                                            + self.study + '-beta-supplement')).exists())
        self.assertFalse(self.child_calls)
        self.assertEqual(result['reserved_network_requests'], 6)

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_three_remaining_requests_complete_nonempty_declared_beta(self):
        self.empty_matrix = False
        self.plan['discovery_limits']['max_requests'] = 6; self.write_plan()
        result = self.execute('input', self.actual_input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1})
        self.assertEqual(len(self.calls), 6)
        self.assertEqual(self.resolution()['requests_used'], 3)
        self.assertEqual(result['reserved_network_requests'], 6)

    def test_oversized_download_is_sealed_as_failed_request_without_retry(self):
        original_transport = self.transport
        def over_limit(url, cap):
            if url.endswith('/' + self.matrix_name):
                self.calls.append(url)
                return b'X' * (cap + 1)
            return original_transport(url, cap)
        self.transport = over_limit
        result = self.execute('input')
        self.assertEqual(result['dispositions'], {'DNAM_INPUT_UNRESOLVED': 1})
        self.assertIn('oversized DNAm source response', result['studies'][0]['reason'])
        receipt = update.read_sealed(self.run_root() / (
            'phases/dnam-input-resolution/source/requests/' + self.study + '-matrix/receipt.json'))
        self.assertEqual((receipt['status'], receipt['accounted_bytes']), ('FAIL', receipt['byte_limit'] + 1))
        self.assertFalse(self.child_calls)
        request_count = len(self.calls)
        self.assertEqual(self.execute('input'), result)
        self.assertEqual(len(self.calls), request_count)

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_completed_resume_and_later_revision_preserve_terminal_without_duplicate_dnam_requests(self):
        first = self.execute('input', self.actual_input_child)
        call_count, child_count = len(self.calls), len(self.child_calls)
        self.assertEqual(update.refresh(self.path, self.output, 'input',
                         transport=lambda *args: self.fail('network on completed resume'),
                         executor=lambda *args: self.fail('child on completed resume')), first)
        self.plan['revision'] = 'software-fixture-v2'; self.write_plan()
        second = self.execute('input', self.actual_input_child)
        self.assertEqual(second['dispositions'], {'SEEN_ACCESSION': 1})
        self.assertEqual(len(self.child_calls), child_count)
        # The new revision rechecks query/Series metadata, but performs no Sample or matrix request.
        self.assertEqual(len(self.calls), call_count + 2)
        self.assertNotIn('update-input', second['phase_receipts'])

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_tampered_prior_dnam_phase_cannot_be_used_as_duplicate_suppression(self):
        self.execute('input', self.actual_input_child)
        response = self.run_root() / ('phases/dnam-input-resolution/source/requests/' + self.study + '-matrix/response.bin')
        with response.open('ab') as handle: handle.write(b'changed')
        self.plan['revision'] = 'tamper-v2'; self.write_plan()
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'Prior DNAm phase output changed'):
            self.execute('input')


class DnamIdatFutureRefreshTests(unittest.TestCase):
    """IDAT transport/producer boundary fixtures, never empirical normalization."""
    def setUp(self):
        DnamFutureRefreshTests.setUp(self)
        from src import pipeline_inputs_2026_09_10 as adapter
        self.adapter = adapter
        self.plan['dnam_input_route'] = update.DNAM_IDAT_ROUTE
        self.plan['extraction_limits'].update(max_extracted_bytes_per_study=500000, max_total_extracted_bytes=1000000)
        self.sample_platforms = ['GPL23976', 'GPL23976']
        self.missing_channel = self.foreign_owner = self.processing_failure = self.verification_failure = False
        self.external_evidence = self.root / 'software-external-idat-evidence.txt'
        self.external_evidence.write_text('SOFTWARE FIXTURE; not numerical validation\n')
        self.process_calls = []
        self.fallback_fetch = patch.object(adapter, 'fetch', side_effect=AssertionError('Hidden IDAT network is forbidden')).start()
        self.processor = patch.object(adapter, 'process_dnam_idat', side_effect=self.software_processor).start()
        patch('audit_system.report_2026_09_10.verify_dnam_idat_provenance', side_effect=self.software_verifier).start()
        self.write_plan()

    write_plan = DnamFutureRefreshTests.write_plan
    execute = DnamFutureRefreshTests.execute
    run_root = DnamFutureRefreshTests.run_root
    resolution = DnamFutureRefreshTests.resolution
    actual_input_child = DnamFutureRefreshTests.actual_input_child
    software_child = DnamFutureRefreshTests.software_child
    matrix_payload = DnamFutureRefreshTests.matrix_payload
    supplement_payload = DnamFutureRefreshTests.supplement_payload
    gz = staticmethod(DnamFutureRefreshTests.gz)

    def series_payload(self):
        return DnamFutureRefreshTests.series_payload(self).replace(b'GPL13534', b'GPL23976')

    def samples_payload(self):
        value = DnamFutureRefreshTests.samples_payload(self).decode()
        if self.missing_channel:
            value = '\n'.join(line for line in value.splitlines() if 'GSM112_R02C02_Red.idat.gz' not in line) + '\n'
        if self.foreign_owner: value = value.replace('/GSM111/suppl/', '/GSM999/suppl/')
        return value.encode()

    def transport(self, url, cap):
        if url.endswith('.idat.gz'):
            self.calls.append(url)
            # Magic-only transport fixture; never passed to an IDAT decoder.
            return gzip.compress(('IDAT software fixture: '+url).encode(), mtime=0)
        return DnamFutureRefreshTests.transport(self, url, cap)

    def software_processor(self, plan_path, acquisition_path, output, *, max_output_bytes):
        self.process_calls.append((plan_path, acquisition_path, output, max_output_bytes))
        output.mkdir()
        (output / 'software-processing-attempt.txt').write_text('SOFTWARE FIXTURE; no R or normalization\n')
        if self.processing_failure: raise ValueError('software fixture processing stop')
        plan = json.loads(plan_path.read_text()); acquisition = json.loads(acquisition_path.read_text())
        self.assertEqual(plan['profile'], self.adapter.DNAM_IDAT_SUPPORTED_PROFILE)
        self.assertEqual(set(plan['limits']), set(self.adapter.DNAM_IDAT_LIMITS))
        self.assertLessEqual(plan['limits']['processing_seconds'], self.plan['update_limits']['max_child_seconds'])
        matrix = output / 'software-beta.csv'; matrix.write_text('probe_id,GSM111,GSM112\ncg00000001,0.1,0.2\n')
        def declared(path): return {'path': str(path.resolve()), 'sha256': update.sha(path)}
        files = {name: {'path': item['path'], 'sha256': item['sha256']}
                 for name, item in {**plan['source_files'], **acquisition['files']}.items()}
        files.update(matrix=declared(matrix), protocol=declared(plan_path), acquisition=declared(acquisition_path),
                     software_external_evidence=declared(self.external_evidence))
        provenance = output / 'software-provenance.json'
        provenance.write_text(json.dumps({'software_fixture': True, 'study_id': self.study,
            'sample_ids': self.sample_ids, 'files': files, 'scientific_status': 'CANNOT VERIFY'}))
        return {'matrix': declared(matrix), 'dnam_idat_provenance': declared(provenance)}

    def software_verifier(self, path, matrix_path=None):
        if self.verification_failure: raise ValueError('software fixture independent verification stop')
        value = json.loads(Path(path).read_text()); self.assertTrue(value['software_fixture'])
        for item in value['files'].values():
            if update.sha(item['path']) != item['sha256']: raise ValueError('Software fixture transitive bytes changed')
        if matrix_path is not None:
            self.assertEqual(Path(matrix_path).resolve(), Path(value['files']['matrix']['path']))
        return {**value, 'measurement': {'software_fixture': True, 'scientific_status': 'CANNOT VERIFY'}}

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_fresh_route_binds_all_candidates_and_requests_without_fallback_fetch(self):
        result = self.execute('input', self.actual_input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1})
        row = self.resolution()['studies'][0]
        self.assertEqual(set(row['sources']), {'series', 'samples', 'matrix', 'dnam_reference', 'dnam_idat_provenance'})
        self.assertEqual(row['input_route'], update.DNAM_IDAT_ROUTE); self.assertFalse(row['author_method_reproduction'])
        self.assertEqual(self.resolution()['requests_used'], 5)
        self.assertEqual(len([url for url in self.calls if url.endswith('.idat.gz')]), 4)
        self.assertFalse(any('/matrix/' in url for url in self.calls)); self.fallback_fetch.assert_not_called()
        version = update.read_sealed(Path(self.child_calls[0]['version_path']))
        self.assertIn('dnam_idat_provenance:raw:GSM111:Grn', version['sources'])
        self.assertIn('dnam_idat_provenance:software_external_evidence', version['sources'])
        frozen = update.read_sealed(self.run_root() / 'version.json')
        self.assertIn('src/dnam_idat_preprocess_2026_09_24.R', frozen['implementation'])
        phase = update.read_sealed(self.run_root() / 'phases/dnam-input-resolution/phase.json')
        self.assertEqual(phase['external_sha256'][str(self.external_evidence)], update.sha(self.external_evidence))
        self.assertEqual(result['reserved_model_calls'], 0)

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_default_processed_beta_route_does_not_call_idat_adapters(self):
        del self.plan['dnam_input_route']; self.empty_matrix = False; self.write_plan()
        result = self.execute('input', self.actual_input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1})
        self.processor.assert_not_called(); self.fallback_fetch.assert_not_called()
        self.assertNotIn('input_route', self.resolution())
        self.assertEqual(set(self.resolution()['studies'][0]['sources']), {'series', 'samples', 'matrix', 'dnam_reference'})

    def test_route_policy_rejects_other_modalities_and_unknown_values_before_network(self):
        for modality, route in [('human_rna', update.DNAM_IDAT_ROUTE), ('human_dnam', 'invented'), ('human_dnam', None)]:
            with self.subTest(modality=modality, route=route):
                self.plan.update(modality=modality, dnam_input_route=route); self.write_plan()
                with self.assertRaisesRegex(ValueError, 'route'): self.execute('input')
                self.assertEqual(self.calls, [])

    def test_missing_channel_and_foreign_owner_stop_before_raw_download(self):
        for label in ('missing_channel', 'foreign_owner'):
            with self.subTest(label=label):
                self.output = self.root / label; self.plan['revision'] = label; self.write_plan()
                self.missing_channel = label == 'missing_channel'; self.foreign_owner = label == 'foreign_owner'
                self.calls.clear(); self.execute('input')
                self.assertEqual(self.resolution()['ready_studies'], 0)
                self.assertFalse(any(url.endswith('.idat.gz') for url in self.calls))
                self.processor.assert_not_called(); self.fallback_fetch.assert_not_called()

    def test_source_species_molecule_platform_and_sample_identity_gate_raw_downloads(self):
        original = self.samples_payload()
        changes = [(b'Homo sapiens', b'Mus musculus'), (b'genomic DNA', b'total RNA'),
                   (b'GPL23976', b'GPL13534'), (b'^SAMPLE = GSM112', b'^SAMPLE = GSM999')]
        for index, (before, after) in enumerate(changes):
            with self.subTest(source_change=before):
                self.output = self.root / ('source-gate-'+str(index)); self.calls.clear()
                with patch.object(self, 'samples_payload', return_value=original.replace(before, after)):
                    self.execute('input')
                self.assertEqual(self.resolution()['ready_studies'], 0)
                self.assertFalse(any(url.endswith('.idat.gz') for url in self.calls))
                self.processor.assert_not_called(); self.fallback_fetch.assert_not_called()

    def test_returned_matrix_must_match_verified_provenance_role(self):
        def conflicting_result(*args, **kwargs):
            value = self.software_processor(*args, **kwargs)
            different = Path(args[2])/'different-software-matrix.csv'
            different.write_text('SOFTWARE FIXTURE; wrong matrix binding\n')
            value['matrix'] = {'path':str(different), 'sha256':update.sha(different)}
            return value
        self.processor.side_effect = conflicting_result
        self.execute('input')
        self.assertEqual(self.resolution()['ready_studies'], 0)
        self.assertEqual(self.child_calls, [])

    def test_contradictory_repeated_source_fields_stop_before_raw_downloads(self):
        original = self.samples_payload()
        for field, value, other in [('organism_ch1', 'Homo sapiens', 'Mus musculus'),
                                    ('molecule_ch1', 'genomic DNA', 'total RNA'),
                                    ('platform_id', 'GPL23976', 'GPL13534')]:
            with self.subTest(field=field):
                self.output = self.root / ('contradictory-' + field); self.calls.clear()
                line = ('!Sample_' + field + ' = ' + value).encode()
                repeated = line + b'\n' + ('!Sample_' + field + ' = ' + other).encode()
                # The parser stores the contradictory physical field as .1.
                with patch.object(self, 'samples_payload', return_value=original.replace(line, repeated)):
                    self.execute('input')
                self.assertEqual(self.resolution()['ready_studies'], 0)
                self.assertIn('contradictory ' + field, self.resolution()['studies'][0]['reason'])
                self.assertFalse(any(url.endswith('.idat.gz') for url in self.calls))
                self.processor.assert_not_called(); self.fallback_fetch.assert_not_called()
                self.assertEqual(self.child_calls, [])

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_equivalent_repeated_organism_and_molecule_values_remain_supported(self):
        repeated = self.samples_payload().replace(b'!Sample_organism_ch1 = Homo sapiens',
            b'!Sample_organism_ch1 = Homo sapiens\n!Sample_organism_ch1 =  HOMO  sapiens ')
        repeated = repeated.replace(b'!Sample_molecule_ch1 = genomic DNA',
            b'!Sample_molecule_ch1 = genomic DNA\n!Sample_molecule_ch1 =  Genomic  dna ')
        with patch.object(self, 'samples_payload', return_value=repeated):
            result = self.execute('input', self.actual_input_child)
        self.assertEqual(result['dispositions'], {'PIPELINE_COMPLETED_THROUGH_INPUT': 1})
        self.assertEqual(self.resolution()['ready_studies'], 1)
        self.assertEqual(len([url for url in self.calls if url.endswith('.idat.gz')]), 4)
        self.fallback_fetch.assert_not_called()

    def test_raw_route_network_count_cap_stops_before_partial_pair_download(self):
        self.plan['discovery_limits']['max_requests'] = 7; self.write_plan(); self.execute('input')
        self.assertIn('every declared IDAT channel', self.resolution()['studies'][0]['reason'])
        self.assertFalse(any(url.endswith('.idat.gz') for url in self.calls)); self.processor.assert_not_called()

    def test_raw_transport_byte_failure_is_accounted_and_never_retried(self):
        original = self.transport
        def oversized(url, cap):
            if url.endswith('.idat.gz'):
                self.calls.append(url); return b'x' * (cap + 1)
            return original(url, cap)
        self.transport = oversized; self.execute('input')
        self.assertEqual(self.resolution()['ready_studies'], 0)
        failures = [update.read_sealed(path) for path in (self.run_root() / 'phases/dnam-input-resolution/source/requests').glob('*idat*/receipt.json')]
        self.assertEqual(len(failures), 1); self.assertEqual(failures[0]['status'], 'FAIL')
        self.assertEqual(failures[0]['accounted_bytes'], failures[0]['byte_limit'] + 1)
        before = len(self.calls); self.execute('input'); self.assertEqual(len(self.calls), before)
        self.processor.assert_not_called(); self.fallback_fetch.assert_not_called()

    def test_raw_transport_timeout_preserves_failure_accounting_and_stops(self):
        original = self.transport
        def timed_out(url, cap):
            if url.endswith('.idat.gz'):
                self.calls.append(url); raise TimeoutError('software fixture network time cap')
            return original(url, cap)
        self.transport = timed_out; self.execute('input')
        self.assertEqual(self.resolution()['ready_studies'], 0)
        failures = [update.read_sealed(path) for path in (self.run_root() / 'phases/dnam-input-resolution/source/requests').glob('*idat*/receipt.json')]
        self.assertEqual(len(failures), 1); self.assertEqual(failures[0]['status'], 'FAIL')
        self.assertEqual(failures[0]['accounted_bytes'], failures[0]['byte_limit'] + 1)
        self.assertIn('time cap', failures[0]['error']); self.processor.assert_not_called()
        before = len(self.calls); self.execute('input'); self.assertEqual(len(self.calls), before)

    def test_controller_rejects_oversized_retained_processing_outputs(self):
        def oversized_result(*args, **kwargs):
            value = self.software_processor(*args, **kwargs)
            (Path(args[2])/'oversized-software-only.bin').write_bytes(b'x' * 500001)
            return value
        self.processor.side_effect = oversized_result; self.execute('input')
        resolution = self.resolution()
        self.assertEqual(resolution['ready_studies'], 0)
        self.assertIn('retained outputs exceed', resolution['studies'][0]['reason'])
        self.assertGreater(resolution['reconstructed_bytes'], 500000)
        self.assertEqual(self.child_calls, [])
        before = len(self.process_calls); self.execute('input'); self.assertEqual(len(self.process_calls), before)

    def test_processing_and_verification_failures_remain_terminal(self):
        for label in ('processing_failure', 'verification_failure'):
            with self.subTest(label=label):
                self.output = self.root / label; self.plan['revision'] = label; self.write_plan()
                self.processing_failure = label == 'processing_failure'; self.verification_failure = label == 'verification_failure'
                self.execute('input'); self.assertEqual(self.resolution()['ready_studies'], 0)
                self.assertIn('software fixture', self.resolution()['studies'][0]['reason'])
                attempt = self.run_root() / ('phases/dnam-input-resolution/source/'+self.study+'/processing/software-processing-attempt.txt')
                self.assertTrue(attempt.is_file())
                before_calls, before_process = len(self.calls), len(self.process_calls); self.execute('input')
                self.assertEqual((len(self.calls), len(self.process_calls)), (before_calls, before_process))
                self.assertEqual(self.child_calls, [])

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_completed_resume_uses_no_network_acquisition_or_processing(self):
        self.execute('input', self.actual_input_child)
        with patch.object(self.adapter, 'acquire_dnam_idat', side_effect=AssertionError('Acquisition on resume')):
            with patch.object(self.adapter, 'process_dnam_idat', side_effect=AssertionError('Processing on resume')):
                update.refresh(self.path, self.output, 'input', transport=lambda *x: self.fail('Network on resume'),
                               executor=lambda *x: self.fail('Child on resume'))

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_changed_external_transitive_bytes_reject_completed_resume(self):
        self.execute('input', self.actual_input_child); self.external_evidence.write_text('Changed software source\n')
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'external'): self.execute('input')

    @unittest.skipUnless(shutil.which('Rscript'), 'Input pipeline runtime metadata requires Rscript on PATH')
    def test_changed_raw_bytes_reject_completed_resume(self):
        self.execute('input', self.actual_input_child)
        raw = next((self.run_root() / 'phases/dnam-input-resolution/source').glob('*/acquisition/uncompressed/*_Grn.idat'))
        raw.write_bytes(raw.read_bytes() + b'changed software fixture')
        with self.assertRaisesRegex(update.UpdateIntegrityError, 'phase output changed'): self.execute('input')

    def test_new_revision_preserves_prior_idat_terminal_stop(self):
        self.missing_channel = True; self.execute('input'); self.missing_channel = False
        self.plan['revision'] = 'explicit-next-revision'; self.write_plan(); self.calls.clear(); self.execute('input')
        self.assertEqual(self.resolution()['studies'][0]['status'], 'PRIOR_DNAM_REFRESH_TERMINAL')
        self.assertFalse(any(url.endswith('.idat.gz') or 'targ=gsm' in url for url in self.calls)); self.processor.assert_not_called()


if __name__ == '__main__':
    unittest.main()
