"""Software fixtures for release/order/accounting; no empirical performance data."""
import copy
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from langchain_core.messages import HumanMessage
from audit_system import heldout_2026_09_17 as controller
from audit_system.agent_workflow_2026_09_10 import sha, write_new
from audit_system.evaluation_2026_09_10 import COMPARISON_METHODS, build_families, reference_stage
from src.pipeline_inputs_2026_09_10 import acquire_study_context
from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError


class ProspectiveTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='clockbase-prospective-fixtures-',
                                                dir=os.environ.get('TMPDIR'))
        self.root = Path(self.temp.name)
        template = self.root / 'git-template'; template.mkdir()
        self.git_environment = patch.dict(os.environ, {
            'GIT_CONFIG_NOSYSTEM': '1', 'GIT_CONFIG_GLOBAL': os.devnull,
            'GIT_TEMPLATE_DIR': str(template)})
        self.git_environment.start(); self.addCleanup(self.git_environment.stop)
        self.patch_root = patch.object(controller, 'ROOT', self.root)
        self.patch_root.start()
        self.addCleanup(self.patch_root.stop); self.addCleanup(self.temp.cleanup)

    def git(self, *arguments):
        result = subprocess.run(['git', *arguments], cwd=self.root, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(f'Git fixture command {arguments!r} failed: {result.stderr.strip()}')
        return result.stdout.strip()

    def commit(self):
        self.git('add', '.')
        self.git('-c', 'user.name=ClockBase software fixture', '-c', 'user.email=fixture@example.test',
                 'commit', '-qm', 'Software fixture freeze')
        return self.git('rev-parse', 'HEAD')

    def fixture(self):
        """A tiny committed repository tests the gate without real GEO access."""
        self.git('init', '-q')
        for relative in controller.IMPLEMENTATIONS:
            path = self.root / relative; path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('# Software fixture identity; never executed.\n')
        assets = self.root / 'pipeline_assets/rna_4_6_2026_09_10/manifest.json'
        assets.parent.mkdir(parents=True); write_new(assets, {'files': {}})
        primary = self.root / 'source'; primary.mkdir()
        series = primary / 'series.soft'
        series.write_text('^SERIES = GSE1\n!Series_status = Public on Jan 01 2020\n'
            '!Series_title = Software fixture\n!Series_summary = Software fixture groups\n'
            '!Series_overall_design = A treated group and its control.\n'
            '!Series_sample_id = GSM1\n!Series_sample_id = GSM2\n')
        samples = primary / 'samples.soft'
        samples.write_text(''.join('^SAMPLE = '+sid+'\n!Sample_title = Software fixture '+sid+
            '\n!Sample_organism_ch1 = Mus musculus\n!Sample_molecule_ch1 = total RNA\n'
            for sid in ('GSM1', 'GSM2')))
        metadata = primary / 'metadata.csv'
        with metadata.open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=['GSM', 'gse_id', 'title', 'organism_ch1', 'molecule_ch1'])
            writer.writeheader(); writer.writerows({'GSM': sid, 'gse_id': 'GSE1',
                'title': 'Software fixture '+sid, 'organism_ch1': 'Mus musculus',
                'molecule_ch1': 'total RNA'} for sid in ('GSM1', 'GSM2'))
        context = acquire_study_context('GSE1', series, primary, publication_methods=False)
        record = {'accession': 'GSE1', 'sample_ids': ['GSM1','GSM2'], 'pubmed_ids': [], 'relations': [],
            'metadata': str(metadata), 'source_sample_records': str(samples),
            'series_source': str(series), 'study_context': str(context)}
        family = build_families([{'accession': 'GSE1', 'sample_ids': ['GSM1','GSM2'],
            'source_kind': 'geo_series_primary', 'source_id': 'software fixture', 'unresolved': []}], {'GSE1'})
        raw = self.root / 'software-fixture-counts.csv'; raw.write_text('gene,GSM1,GSM2\nE1,10,12\n')
        case = {'study_id': 'GSE1', 'family_id': family['families'][0]['family_id'],
            'modality': 'mouse_rna', 'metadata_source': {'kind': 'pinned', 'record': record,
                'source_sha256': {str(p): sha(p) for p in (series, samples, metadata, context)}},
            'molecular_source': {'kind': 'raw_counts', 'path': str(raw), 'sha256': sha(raw)},
            'expected_linkage': {k: record[k] for k in ('sample_ids','pubmed_ids','relations')}}
        paths = {'family_manifest': family, 'exposure_manifest': {'exposed_studies': ['GSE1']},
                 'cohort_manifest': {'schema': 'clockbase-prospective-cohort-v1',
                                     'split': 'development', 'cases': [case]}}
        plan = {'schema': controller.PLAN_SCHEMA, 'split': 'development', 'methods': list(COMPARISON_METHODS),
            'baseline_rule': 'retain_producer_RESULT_rows_without_audit', 'preprocessing': 'clockbase_legacy',
            'producer_rule': 'one_fresh_original_CodingAgent_per_ready_study',
            'comparison_scope': 'source_defined_interventions', 'target': 'Source-defined interventions',
            'provider': 'openai', 'model': 'gpt-5.6-sol', 'concurrency': 1, 'automatic_retries': 0,
            'max_samples_per_study': 2000,
            'implementation_sha256': {p: sha(self.root / p) for p in controller.IMPLEMENTATIONS},
            'reference_procedure': {'kind': 'separate_source_first_agent'},
            'analysis': copy.deepcopy(controller.ANALYSIS),
            'primary_limits': {'max_calls': 6, 'max_cost': 6, 'output_tokens': 4000,
                'prompt_bytes': 48000, 'input_per_million': 5, 'output_per_million': 20},
            'review_limits': {'max_calls': 1, 'max_cost': 1, 'output_tokens': 6000,
                'prompt_bytes': 48000, 'input_per_million': 5, 'output_per_million': 20},
            'budget': {'initial_authorized_cap_usd': 5000, 'previously_allocated_usd': 2,
                       'this_allocation_usd': 8},
            'asset_manifest': {'path': str(assets), 'sha256': sha(assets)}}
        for key, value in paths.items():
            path = self.root / (key + '.json'); write_new(path, value)
            plan[key] = {'path': str(path), 'sha256': sha(path)}
        protocol = self.root / 'protocol.md'; protocol.write_text('Software fixture protocol, no empirical claims.\n')
        plan['protocol'] = {'path': str(protocol), 'sha256': sha(protocol)}
        self.plan_path = self.root / 'plan.json'; write_new(self.plan_path, plan)
        commit = self.commit(); self.release_path = self.root / 'release.json'
        write_new(self.release_path, {'schema': controller.RELEASE_SCHEMA, 'state': 'OPEN',
            'freeze_commit': commit, 'plan': {'path': str(self.plan_path), 'sha256': sha(self.plan_path)}})
        self.commit()
        return plan, case

    def refreeze(self, plan):
        self.plan_path.write_text(json.dumps(plan))
        commit = self.commit(); release = controller.read_json(self.release_path)
        release.update(freeze_commit=commit, plan={'path': str(self.plan_path), 'sha256': sha(self.plan_path)})
        self.release_path.write_text(json.dumps(release)); self.commit()

    def missing_state(self, archive, case):
        evidence = self.root / 'software-reference-timeout.txt'
        evidence.write_text('Software fixture: adjudicator did not return a response.\n')
        missing = self.root / 'source-reference-missingness.json'
        write_new(missing, {'schema': 'clockbase-reference-missingness-v1', 'study_id': case['study_id'],
            'family_id': case['family_id'], 'archive_config_sha256': sha(archive / 'config.json'),
            'phase': 'SOURCE_DESIGN', 'reason': 'Software fixture: no independent response',
            'evidence': [{'path': str(evidence), 'sha256': sha(evidence)}]})
        states = self.root / 'source-reference-states.json'
        write_new(states, {'schema': 'clockbase-source-reference-states-v1',
            'archive_config_sha256': sha(archive / 'config.json'), 'studies': [
                {'study_id': case['study_id'], 'kind': 'MISSING',
                 'artifact': {'path': str(missing), 'sha256': sha(missing)}}]})
        return states

    def test_joint_gate_blocks_uncommitted_closed_and_changed_release_before_sources(self):
        self.fixture(); release = controller.read_json(self.release_path)
        for change in ({**release, 'state': 'CLOSED'}, {**release, 'extra': 'uncommitted'}):
            self.release_path.write_text(json.dumps(change))
            with patch('src.pipeline_inputs_2026_09_10.acquire_metadata') as acquire:
                with self.assertRaises(PipelineIntegrityError):
                    controller.run(self.release_path, self.root / 'run')
                acquire.assert_not_called(); self.assertFalse((self.root / 'run').exists())
        self.release_path.write_bytes(controller.git_payload(self.release_path, 'HEAD'))
        path = self.root / controller.IMPLEMENTATIONS[0]; path.write_text('changed implementation')
        with self.assertRaisesRegex(PipelineIntegrityError, 'changed'):
            controller.validate_release(self.release_path)

    def test_committed_plan_cannot_replace_legacy_or_exceed_allocation(self):
        plan, _ = self.fixture()
        for key, value in [('preprocessing','published_mouse_rle'), ('budget', {
                'initial_authorized_cap_usd':5000, 'previously_allocated_usd':4999, 'this_allocation_usd':8})]:
            changed = copy.deepcopy(plan); changed[key] = value; self.refreeze(changed)
            with self.assertRaises(PipelineIntegrityError): controller.validate_release(self.release_path)

    def test_release_rejects_missing_and_conflicting_input_descriptors_before_sources(self):
        plan,case=self.fixture();path=Path(plan['cohort_manifest']['path'])
        original=controller.read_json(path)
        variants=[]
        missing=copy.deepcopy(original);missing['cases'][0]['molecular_source']={'kind':'raw_counts'};variants.append(missing)
        changed=copy.deepcopy(original);other=copy.deepcopy(case);other['study_id']='GSE2'
        other['molecular_source']['sha256']='f'*64;changed['cases'].append(other);variants.append(changed)
        for cohort in variants:
            path.write_text(json.dumps(cohort));plan['cohort_manifest']['sha256']=sha(path);self.refreeze(plan)
            with patch('src.pipeline_inputs_2026_09_10.acquire_metadata') as acquire:
                with self.assertRaises(PipelineIntegrityError): controller.run(self.release_path,self.root/'run')
                acquire.assert_not_called()
        self.assertFalse((self.root/'run').exists())

    def test_nonfinite_budget_is_rejected(self):
        plan,_=self.fixture();plan['budget']['previously_allocated_usd']=float('nan');self.refreeze(plan)
        with self.assertRaisesRegex(PipelineIntegrityError,'Nonfinite'):
            controller.validate_release(self.release_path)

    def test_resumed_producer_and_review_verify_external_api_ledgers(self):
        for name,slot in [('GSE1/coding_agent','baseline-api/GSE1/primary-000'),
                          ('GSE1/independent_audit/review','api/GSE1-independent_audit-000'),
                          ('GSE1/self_review/review','api/GSE1-self_review-000')]:
            with self.subTest(stage=name):
                archive=self.root/name.replace('/','-');store=controller.ProspectiveStore(archive,{'software_fixture':True})
                (archive/name).parent.mkdir(parents=True,exist_ok=True)
                api=archive/slot;api.mkdir(parents=True)
                reservation=api/'reservation.json';write_new(reservation,{'reserved_usd':.5})
                def action(directory):write_new(directory/'software-result.json',{'ok':True});return {'ok':True}
                store.run(name,action)
                self.assertEqual(store.run(name,action),{'ok':True})
                reservation.write_text('{"reserved_usd":0}')
                with self.assertRaisesRegex(PipelineIntegrityError,'API request'):
                    store.run(name,action)

    def test_source_packets_exist_before_any_numerical_read_and_resume_unchanged(self):
        self.fixture(); archive = self.root / 'run'
        with patch.object(controller, '_input_path', side_effect=AssertionError('Numerical read before reference')):
            result = controller.run(self.release_path, archive)
            self.assertEqual(result['phase'], 'SOURCE_ONLY')
            packet = controller.read_json(result['sources'][0]['source_packet'])
            controller.validate_reference_packet_context(packet)
            self.assertNotIn('primary', json.dumps(packet)); self.assertNotIn('clock_predictions', json.dumps(packet))
            before = {str(p.relative_to(archive)):sha(p) for p in archive.rglob('*') if p.is_file()}
            controller.run(self.release_path, archive)
            self.assertEqual(before, {str(p.relative_to(archive)):sha(p) for p in archive.rglob('*') if p.is_file()})
            with self.assertRaisesRegex(PipelineIntegrityError, 'reference state'):
                controller.run(self.release_path, archive, until='producer')

    def test_reference_packet_mutation_cannot_enter_a_source_judgment(self):
        self.fixture(); result = controller.run(self.release_path, self.root / 'run')
        packet = controller.read_json(result['sources'][0]['source_packet'])
        packet['source_packet']['target'] = 'Changed after source export'
        with self.assertRaisesRegex(PipelineIntegrityError, 'differs'):
            controller.validate_reference_packet_context(packet)

    def test_terminal_failure_is_retained_without_retry_and_integrity_is_not_missingness(self):
        store = controller.ProspectiveStore(self.root / 'run', {'software_fixture': True}); calls = []
        def missing(_directory): calls.append(1); raise TimeoutError('Software fixture timeout')
        result = controller.terminal_stage(store, 'source', missing)
        self.assertEqual(result['status'], 'CANNOT VERIFY')
        self.assertEqual(controller.terminal_stage(store, 'source', missing), result)
        self.assertEqual(len(calls), 1)
        def corrupted(_directory): raise PipelineIntegrityError('Software fixture mutation')
        with self.assertRaises(PipelineIntegrityError): controller.terminal_stage(store, 'changed', corrupted)
        self.assertFalse((store.root / 'changed/complete.json').exists())
        with self.assertRaisesRegex(PipelineIntegrityError, 'Interrupted'):
            controller.terminal_stage(store, 'changed', missing)

    def test_actor_ceiling_is_enforced_in_code_and_unknown_spending_stays_reserved(self):
        limits = {'max_calls':6, 'max_cost':6, 'output_tokens':4000,'prompt_bytes':48000,
                  'input_per_million':5,'output_per_million':20}
        with patch('src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke', return_value='software fixture') as transport:
            client = controller.ActorClient(self.root/'api', 'primary', limits=limits,
                                            total_calls=100, total_cost=100, model='gpt-5.6-sol')
            for _ in range(6): client.invoke([HumanMessage(content='Software fixture')])
            with self.assertRaisesRegex(RuntimeError, 'turn limit'):
                client.invoke([HumanMessage(content='Software fixture')])
            self.assertEqual(transport.call_count,6)
            slot=self.root/'api/primary-000';slot.mkdir();write_new(slot/'reservation.json',{'reserved_usd':5.99})
            another=controller.ActorClient(self.root/'api','primary',limits=limits,total_calls=100,total_cost=100,model='gpt-5.6-sol')
            another.calls=1
            with self.assertRaisesRegex(RuntimeError,'allocation exhausted'):
                another.invoke([HumanMessage(content='Software fixture')])
            self.assertEqual(transport.call_count,6)

    def test_missing_source_reference_still_permits_legacy_producer_and_all_methods(self):
        _, case = self.fixture(); archive=self.root/'run';controller.run(self.release_path,archive)
        states=self.missing_state(archive,case); seen=[]
        def preprocess(raw,species,directory,**kwargs):
            self.assertEqual(kwargs['normalization'],'clockbase_legacy')
            self.assertTrue((archive/'GSE1/reference_state/complete.json').is_file())
            self.assertNotIn('references',kwargs['source']);seen.append('preprocess')
            path=directory/'software-processed.csv';path.write_text('software fixture only\n');return {'path':str(path)}
        def predict(processed,species,directory):
            seen.append('predict');path=directory/'software-predictions.csv'
            path.write_text('sample_id,fixture_clock\nGSM1,0\nGSM2,1\n')
            return {'path':str(path),'clocks':['fixture_clock']}
        def producer(packet,clocks,client,directory):
            seen.append('producer');self.assertNotIn('references',packet)
            value={'rows':[{'clock':'fixture_clock','case_ids':[],'control_ids':[],
                'statistical_test':'','p_value':'','mean_difference':'','effect_size':'',
                'status':'CANNOT VERIFY','reason':'Software fixture abstention'}]}
            write_new(directory/'software-primary.json',value);return value
        with patch('src.pipeline_inputs_2026_09_10.checked_assets'), \
                patch('src.pipeline_inputs_2026_09_10.preprocess_rna',side_effect=preprocess), \
                patch('src.pipeline_inputs_2026_09_10.predict_rna',side_effect=predict), \
                patch('src.pipeline_2026_09_10.primary_analysis',side_effect=producer), \
                patch('src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke',side_effect=TimeoutError('Software fixture no response')):
            result=controller.run(self.release_path,archive,until='methods',reference_states=states)
            self.assertEqual(seen,['preprocess','predict','producer'])
            self.assertEqual({r['method'] for r in result['records']},set(COMPARISON_METHODS))
            self.assertEqual(result['attempted_studies'],1)
            before={str(p.relative_to(archive)):sha(p) for p in archive.rglob('*') if p.is_file()}
            controller.run(self.release_path,archive,until='methods',reference_states=states)
            self.assertEqual(seen,['preprocess','predict','producer'])
            self.assertEqual(before,{str(p.relative_to(archive)):sha(p) for p in archive.rglob('*') if p.is_file()})

    def test_preparation_failure_retains_four_unexecuted_method_slots(self):
        _,case=self.fixture();archive=self.root/'run';controller.run(self.release_path,archive)
        states=self.missing_state(archive,case)
        with patch('src.pipeline_inputs_2026_09_10.checked_assets'), \
                patch('src.pipeline_inputs_2026_09_10.preprocess_rna',side_effect=ValueError('Software fixture unavailable input')), \
                patch('src.pipeline_2026_09_10.primary_analysis') as producer:
            result=controller.run(self.release_path,archive,until='methods',reference_states=states)
            producer.assert_not_called();self.assertEqual(len(result['records']),4)
            self.assertTrue(all(not r['accepted'] and r['review_execution_status']=='NOT_RUN_SHARED_PREREQUISITE'
                                for r in result['records']))
            cases,_=reference_stage(archive,'00_inputs','prepared-cases.json')
            self.assertEqual(cases[0]['preparation_status'],controller.UNAVAILABLE)


if __name__ == '__main__':
    unittest.main()
