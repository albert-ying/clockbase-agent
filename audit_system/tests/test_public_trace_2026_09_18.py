"""Software fixtures for public trace integrity; these are not study results."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy import stats

from audit_system.agent_workflow_2026_09_10 import sha
from audit_system.public_trace_2026_09_18 import (
    EFFECT, STAGES, _failed_retry_request, _ids, _matrix_qc, _safe, _snapshot,
    export, export_failed, numerical_record, prepare, prepare_failed, verify_numbers,
    grouping_readback,
)
from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError


class SourceReferencePublicTraceTests(unittest.TestCase):
    """Declared software fixtures only; no empirical labels or outcomes."""

    def setUp(self):
        from audit_system.tests.test_source_readiness_2026_09_24 import SourceReadinessTests
        from audit_system import public_trace_2026_09_18 as trace
        self.trace = trace
        self.fixture = SourceReadinessTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.root = self.fixture.root

    def plan(self, versions, name='export-plan.json'):
        registry = self.root/(name+'.registry.json')
        registry.write_text(json.dumps([self.fixture.frozen_binding(v) for v in versions])+'\n')
        plan = self.root/name
        result = self.trace.prepare_source_reference(
            self.fixture.release['path'], self.fixture.release['sha256'], registry, sha(registry), plan)
        self.assertEqual(result['new_model_calls'], 0)
        return plan, registry

    def accepted(self, name='review-001', supersedes=None, mutate=None):
        value = self.fixture.submission(comparison=True)
        if mutate is not None:
            mutate(value)
        result = self.fixture.freeze(value, name, supersedes=supersedes)
        self.assertEqual(result['validation_status'], 'ACCEPTED_SOURCE_ATTRIBUTION')
        return result

    def test_source_reference_preserves_domains_denominators_and_unknown_model_cost(self):
        version = self.accepted()
        plan, _ = self.plan([version])
        output = self.root/'public'
        result = self.trace.export_source_reference(plan, output)
        trace = json.loads((output/self.fixture.sid/'audit_trace.json').read_text())
        cohort = json.loads((output/'readiness.json').read_text())
        history = json.loads((output/self.fixture.sid/'downloads/source-reference-history.json').read_text())
        self.assertEqual((result['traces'], result['enrolled_studies']), (1, 80))
        self.assertEqual((cohort['allowed_studies'], cohort['administrative_holds']), (78, 2))
        self.assertEqual(len(cohort['studies']), 80)
        held = [r for r in cohort['studies'] if not r['source_access_allowed']]
        self.assertEqual(len(held), 2)
        self.assertTrue(all(set(r['domains'].values()) == {'NOT_ASSESSED'} for r in held))
        self.assertEqual(len(history['versions'][0]['validation']['domains']), 9)
        self.assertEqual(len(trace['domains']), 7)
        self.assertEqual(trace['models'], [])
        self.assertTrue(all(v is None for v in trace['resource_usage'].values()))
        self.assertEqual(trace['overall']['status'], 'CANNOT VERIFY')
        self.assertTrue(all(c['status'] == 'CANNOT VERIFY' for c in trace['comparisons']))
        self.assertIsNone(cohort['scientific_error_rate'])
        self.assertEqual(cohort['accepted_scientific_repairs'], 0)
        with (output/'readiness.csv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 80)
        with (output/'domain-counts.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 9)
        self.assertTrue(all(sum(int(r[k]) for k in ('PASS','FAIL','CANNOT VERIFY','NOT_ASSESSED')) == 78 for r in rows))
        for artifact in trace['artifacts']:
            self.assertEqual(sha(output/self.fixture.sid/artifact['download_path']), artifact['sha256'])
        by_id = {a['id']: a for a in trace['artifacts']}
        self.assertEqual(len(by_id), len(trace['artifacts']))
        self.assertEqual(by_id['readiness-json']['download_path'], 'downloads/readiness.json')
        self.assertEqual(by_id['readiness-csv']['download_path'], 'downloads/readiness.csv')
        self.assertNotEqual(by_id['readiness-json']['sha256'], by_id['readiness-csv']['sha256'])

    def test_long_inline_prose_retains_full_source_fields_without_shortening_quotes(self):
        rationale = 'SOFTWARE FIXTURE long source rationale ' * 500
        version = self.accepted(mutate=lambda v: v['domains']['sample_identity'].update(rationale=rationale))
        plan, _ = self.plan([version])
        output = self.root/'long-prose-public'
        self.trace.export_source_reference(plan, output)
        trace = json.loads((output/self.fixture.sid/'audit_trace.json').read_text())
        summary = next(d for d in trace['domains'] if d['id'] == 'sample_identity')['model_assessment']['summary']
        self.assertLessEqual(len(summary), 12000)
        self.assertIn('Inline text omitted', summary)
        history = json.loads((output/self.fixture.sid/'downloads/source-reference-history.json').read_text())
        self.assertEqual(history['versions'][0]['judgment']['domains']['sample_identity']['rationale'], rationale)
        evidence = next(d for d in trace['domains'] if d['id'] == 'sample_identity')['model_assessment']['evidence']
        original_quotes = {v['quote'] for v in history['versions'][0]['judgment']['evidence'].values()}
        self.assertTrue(all(e['quote'] in original_quotes for e in evidence))
        self.assertIn('Inline text omitted', self.trace._source_reference_inline('   '))

    def test_rejected_current_version_does_not_fall_back_to_earlier_acceptance(self):
        first = self.accepted()
        invalid = self.fixture.submission(comparison=True)
        invalid['comparisons'][0]['control_ids'] = invalid['comparisons'][0]['case_ids'][:]
        second = self.fixture.freeze(invalid, 'review-002', supersedes=self.fixture.frozen_binding(first))
        self.assertEqual(second['validation_status'], 'REJECTED')
        plan, _ = self.plan([first, second])
        output = self.root/'rejected-public'
        result = self.trace.export_source_reference(plan, output)
        trace = json.loads((output/self.fixture.sid/'audit_trace.json').read_text())
        cohort = json.loads((output/'readiness.json').read_text())
        row = next(r for r in cohort['studies'] if r['study_id'] == self.fixture.sid)
        self.assertEqual(result['traces'], 1)
        self.assertEqual(result['source_validated_studies'], 0)
        self.assertEqual(set(row['domains'].values()), {'NOT_ASSESSED'})
        self.assertEqual([e['status'] for e in trace['events']], ['PASS', 'FAIL'])
        self.assertTrue(all(d['model_assessment']['status'] == 'CANNOT VERIFY'
                            and d['verification']['status'] == 'FAIL' for d in trace['domains']))
        self.assertEqual(trace['comparisons'], [])
        self.assertTrue(all('NOT_ASSESSED' in d['model_assessment']['summary'] for d in trace['domains']))

    def test_omitting_a_superseded_version_is_rejected(self):
        first = self.accepted()
        second = self.accepted('review-002', self.fixture.frozen_binding(first))
        with self.assertRaisesRegex(ValueError, 'complete ordered supersession chain'):
            self.plan([second])

    def test_changed_raw_source_stops_before_output_creation(self):
        version = self.accepted()
        plan, _ = self.plan([version])
        source = Path(self.fixture.allowed['studies'][0]['sources']['samples']['path'])
        source.write_bytes(source.read_bytes()+b'\nSOFTWARE FIXTURE CHANGED SOURCE\n')
        with self.assertRaises(PipelineIntegrityError):
            self.trace.export_source_reference(plan, self.root/'changed-raw-public')
        self.assertFalse((self.root/'changed-raw-public').exists())

    def test_changed_frozen_submission_stops_before_output_creation(self):
        version = self.accepted()
        plan, _ = self.plan([version])
        submission = Path(version['frozen_submission']).with_name('submission.json')
        submission.write_bytes(submission.read_bytes()+b' ')
        with self.assertRaisesRegex(ValueError, 'output changed'):
            self.trace.export_source_reference(plan, self.root/'changed-submission-public')
        self.assertFalse((self.root/'changed-submission-public').exists())

    def test_changed_registry_or_export_identity_cannot_silently_replan(self):
        version = self.accepted()
        plan, registry = self.plan([version])
        registry.write_text('[]\n')
        with self.assertRaises(PipelineIntegrityError):
            self.trace.export_source_reference(plan, self.root/'changed-registry-public')
        data = json.loads(plan.read_text())
        data['export_identity']['implementation_sha256']['audit_system/public_trace_2026_09_18.py'] = '0'*64
        plan.write_text(json.dumps(data))
        with self.assertRaisesRegex(ValueError, 'implementation changed'):
            self.trace.export_source_reference(plan, self.root/'changed-code-public')

    def test_private_prose_is_explicitly_withheld_and_original_is_unchanged(self):
        secret = 'SOFTWARE FIXTURE private source /home/fixture/private/source.tsv'
        version = self.accepted(mutate=lambda v: v['domains']['sample_identity'].update(rationale=secret))
        original = Path(version['frozen_submission']).with_name('submission.json')
        before = sha(original)
        plan, _ = self.plan([version])
        output = self.root/'withheld-public'
        self.trace.export_source_reference(plan, output)
        for path in output.rglob('*.json'):
            text = path.read_text()
            self.assertNotIn(secret, text)
            self.assertNotIn('/home/fixture/', text)
        trace = (output/self.fixture.sid/'audit_trace.json').read_text()
        self.assertIn('Text withheld from public export; SHA256', trace)
        self.assertEqual(sha(original), before)
        self.assertIn(secret, original.read_text())

    def test_unparseable_rejected_submission_is_preserved_by_hash(self):
        from audit_system import evaluation_2026_09_10 as readiness
        packet, _ = self.fixture.project()
        source = self.root/'malformed-software-fixture.json'
        source.write_text('{"declared_software_fixture": ')
        version = readiness.freeze_source_readiness_submission(
            self.fixture.release['path'], self.fixture.release['sha256'], self.fixture.batch,
            self.fixture.sid, packet, source, 'malformed-001', self.root/'malformed-store')
        self.assertEqual(version['validation_status'], 'REJECTED')
        plan, _ = self.plan([version])
        output = self.root/'malformed-public'
        self.trace.export_source_reference(plan, output)
        history = json.loads((output/self.fixture.sid/'downloads/source-reference-history.json').read_text())
        self.assertTrue(history['versions'][0]['unreadable_submission_preserved_by_hash'])
        self.assertEqual(history['versions'][0]['submission_sha256'], sha(source))
        self.assertIsNone(history['versions'][0]['judgment'])

    def test_inline_evidence_cap_preserves_quotes_and_discloses_omissions(self):
        evidence = {}
        for index in range(8):
            evidence[str(index)] = {'entity': 'GSM1', 'field': 'SOFTWARE_FIXTURE_FIELD',
                'line_start': index+1, 'line_end': index+1,
                'quote': 'SOFTWARE FIXTURE quote '+str(index)}
        evidence['private'] = {**evidence['0'], 'quote': '/home/fixture/private/text'}
        shown, note = self.trace._source_reference_evidence({'evidence': evidence}, list(evidence))
        self.assertEqual(len(shown), 5)
        self.assertEqual([r['quote'] for r in shown], [evidence[str(i)]['quote'] for i in range(5)])
        self.assertIn('5 of 9', note)
        self.assertIn('1 omitted', note)

    def test_export_refuses_to_overwrite_an_existing_bundle(self):
        version = self.accepted()
        plan, _ = self.plan([version])
        output = self.root/'no-overwrite'
        self.trace.export_source_reference(plan, output)
        before = {str(p): sha(p) for p in output.rglob('*') if p.is_file()}
        with self.assertRaisesRegex(ValueError, 'new directory'):
            self.trace.export_source_reference(plan, output)
        self.assertEqual(before, {str(p): sha(p) for p in output.rglob('*') if p.is_file()})


class PublicTraceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.frame = pd.DataFrame({'clock': [1., 2., 4., 8.]},
                                  index=['GSM1','GSM2','GSM3','GSM4'])
        self.row = {'clock':'clock','status':'RESULT','case_ids':['GSM1','GSM2'],
                    'control_ids':['GSM3','GSM4'],'statistical_test':'Welch unpaired two-sided t-test',
                    'mean_difference':-4.5,
                    'p_value':float(stats.ttest_ind([1,2],[4,8],equal_var=False).pvalue),
                    'effect_size':-4.5/np.sqrt(4.25),'reason':'Software fixture'}

    def test_conditional_arithmetic_does_not_invent_effect_definition(self):
        checks=verify_numbers(self.row,self.frame)
        self.assertEqual([c['status'] for c in checks],['PASS','PASS','CANNOT VERIFY'])
        self.assertIsNone(numerical_record(self.row)['effect_size_definition'])
        self.assertIn('not supplied',checks[-1]['reason'])

    def test_public_grouping_readback_retains_prose_and_reversed_target_failures(self):
        packet = {'study_id':'GSE1', 'source_sha256':'a'*64, 'samples':[
            {'sample_id':sid, 'characteristics_ch1':'treatment: '+('drug' if sid in ['GSM1','GSM2'] else 'none')}
            for sid in self.frame.index]}
        comparison = {'comparison_id':'fixture','comparison_name':'Software fixture',
                      'target':'Drug-treated samples versus untreated samples',
                      'case_ids':['GSM1','GSM2'], 'control_ids':['GSM3','GSM4']}
        record, checks = grouping_readback(packet,[comparison])
        self.assertEqual([c['status'] for c in checks], ['CANNOT VERIFY','PASS'])
        self.assertEqual(record['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(record['repairs_applied'], [])
        witness = record['comparisons'][0]['source_membership_witnesses'][0]
        self.assertFalse(witness['applied'])
        self.assertEqual(witness['case_ids'], comparison['case_ids'])
        comparison['target'] = 'none versus drug'
        record, checks = grouping_readback(packet,[comparison])
        self.assertEqual(checks[0]['status'], 'FAIL')
        self.assertEqual(record['scientific_status'], 'CANNOT VERIFY')

    def test_explicit_effect_and_changed_estimates(self):
        row={**self.row,'effect_size_definition':EFFECT}
        self.assertEqual([c['status'] for c in verify_numbers(row,self.frame)],['PASS']*3)
        row['effect_size'] *= -1
        self.assertEqual(verify_numbers(row,self.frame)[-1]['status'],'FAIL')
        row['p_value']=0.8
        self.assertEqual(verify_numbers(row,self.frame)[1]['status'],'FAIL')

    def test_unperformed_and_zero_estimates_are_distinct(self):
        row={**self.row,'status':'CANNOT VERIFY','p_value':'','mean_difference':'','effect_size':''}
        result=numerical_record(row)
        self.assertEqual(result['state'],'ABSTAINED')
        self.assertIsNone(result['p_value'])
        self.assertEqual(len(verify_numbers(row,self.frame)),1)
        row.update(status='RESULT',p_value=0,mean_difference=0,effect_size=0)
        self.assertEqual(numerical_record(row)['state'],'COMPUTED')
        self.assertEqual(numerical_record(row)['p_value'],0.)
        row['status']='CANNOT VERIFY'
        self.assertEqual(numerical_record(row)['state'],'INCONSISTENT')

    def test_paired_one_sided_and_unknown_methods_not_substituted(self):
        for method in ('paired t-test','one-sided Welch t-test','ANCOVA','invented method'):
            with self.subTest(method=method):
                checks=verify_numbers({**self.row,'statistical_test':method},self.frame)
                self.assertEqual(next(c for c in checks if c['name']=='P value')['status'],'CANNOT VERIFY')

    def test_missing_clock_unknown_or_overlapping_membership(self):
        for row in ({**self.row,'clock':'absent'}, {**self.row,'case_ids':['GSM1','GSM3']},
                    {**self.row,'case_ids':['GSM1','GSM8']}):
            self.assertEqual(verify_numbers(row,self.frame)[0]['status'],'FAIL')
        for ids in (['GSM1','GSM1'],['GSM8'],[{}]):
            with self.assertRaises(ValueError): _ids(ids,set(self.frame.index))

    def test_invalid_numeric_and_private_text_refused(self):
        for value in (True,float('nan'),float('inf'),-0.1,1.1):
            with self.assertRaises(ValueError): numerical_record({**self.row,'p_value':value})
        for text in ('/Entropy/private/source.csv','Bearer sensitive-token','/home/user/key',
                     'from /mnt/private/file.csv', 'Contact example@example.org'):
            with self.assertRaises(ValueError): _safe({'quote':text})

    def test_missing_beta_is_preserved_and_out_of_range_is_failure(self):
        path=self.root/'beta.csv'
        path.write_text('probe,GSM1,GSM2\ncg1,0.5,NA\ncg2,0,1\n')
        features,samples,check=_matrix_qc(path,{'GSM1','GSM2'},'human_dnam')
        self.assertEqual((features,samples),(2,2));self.assertEqual(check['status'],'CANNOT VERIFY')
        self.assertIn('1 missing',check['reason'])
        path.write_text('probe,GSM1,GSM2\ncg1,0.5,1.1\n')
        self.assertEqual(_matrix_qc(path,{'GSM1','GSM2'},'human_dnam')[2]['status'],'FAIL')

    def test_rna_domain_and_identity(self):
        path=self.root/'counts.csv'
        path.write_text('gene,GSM1,GSM2\ngene1,5,0\n')
        self.assertEqual(_matrix_qc(path,{'GSM1','GSM2'},'mouse_rna')[2]['status'],'PASS')
        path.write_text('gene,GSM1,GSM2\ngene1,-1,0.5\n')
        self.assertEqual(_matrix_qc(path,{'GSM1','GSM2'},'mouse_rna')[2]['status'],'FAIL')
        with self.assertRaises(ValueError): _matrix_qc(path,{'GSM1','GSM3'},'mouse_rna')
        path.write_text('')
        with self.assertRaisesRegex(ValueError,'empty'): _matrix_qc(path,{'GSM1','GSM2'},'mouse_rna')

    def _stage_fixture(self):
        run=self.root/'run';run.mkdir()
        (run/'config.json').write_text(json.dumps({'study':'GSE1'}))
        for stage in STAGES:
            directory=run/stage;directory.mkdir()
            (directory/'evidence.txt').write_text('software integrity fixture')
            receipt=directory/'complete.json'
            receipt.write_text(json.dumps({'outputs':{'evidence.txt':sha(directory/'evidence.txt')}}))
            receipt.with_suffix('.sha256').write_text(sha(receipt))
        return run

    def test_missing_stage_payload_and_receipt_changes_rejected(self):
        run=self._stage_fixture()
        self.assertEqual(len(_snapshot(run)),25)
        path=run/STAGES[0]/'evidence.txt';path.write_text('changed')
        with self.assertRaises(PipelineIntegrityError): _snapshot(run)
        path.write_text('software integrity fixture')
        (run/STAGES[0]/'complete.json').write_text('{}')
        with self.assertRaises(PipelineIntegrityError): _snapshot(run)
        (run/STAGES[0]/'complete.json').unlink()
        with self.assertRaises(ValueError): _snapshot(run)

    def test_recomputed_archive_receipt_cannot_override_frozen_plan(self):
        run=self._stage_fixture();plan=self.root/'plan.json'
        prepare(run,plan)
        directory=run/STAGES[0]
        (directory/'evidence.txt').write_text('changed with a new receipt')
        receipt=directory/'complete.json'
        receipt.write_text(json.dumps({'outputs':{'evidence.txt':sha(directory/'evidence.txt')}}))
        receipt.with_suffix('.sha256').write_text(sha(receipt))
        with self.assertRaises(PipelineIntegrityError): export(plan,self.root/'public')

    def test_stage_output_path_cannot_escape_stage(self):
        run=self._stage_fixture();directory=run/STAGES[0]
        receipt=directory/'complete.json'
        receipt.write_text(json.dumps({'outputs':{'../config.json':sha(run/'config.json')}}))
        receipt.with_suffix('.sha256').write_text(sha(receipt))
        with self.assertRaises(PipelineIntegrityError): _snapshot(run)

    def _failed_fixture(self, *, responses=1):
        run = self.root / 'failed-run'; run.mkdir()
        (run / 'config.json').write_text(json.dumps({'study': 'GSE999001', 'max_prompt_bytes': 48000}))
        packet = {'study_id': 'GSE999001', 'target': 'Software fixture target',
                  'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE999001',
                  'source_sha256': 'a' * 64,
                  'samples': [{'sample_id': 'GSM999001', 'treatment': 'control'}],
                  'study_documents': [{'source_id': 'GEO:GSE999001',
                      'source_url': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE999001',
                      'source_sha256': 'b' * 64, 'fields': {'title': 'Software fixture'}}]}
        for index, name in enumerate(STAGES[:2]):
            stage = run / name; (stage / 'attempt-000').mkdir(parents=True)
            payload = stage / 'attempt-000' / ('metadata.csv' if index == 0 else 'input.csv')
            payload.write_text('sample_id\nGSM999001\n')
            receipt = stage / 'complete.json'
            receipt.write_text(json.dumps({'result': {'sample_ids': ['GSM999001'],
                'metadata': str(payload)}, 'outputs': {
                str(payload.relative_to(stage)): sha(payload)}}))
            (stage / 'complete.sha256').write_text(sha(receipt) + '\n')
        attempt = run / '03_blinded_design' / 'attempt-000'; (attempt / 'response-validation').mkdir(parents=True)
        packet['source_sha256'] = sha(run / '01_metadata/attempt-000/metadata.csv')
        (attempt / 'source-packet.json').write_text(json.dumps(packet))
        (attempt / 'source-transport.json').write_text(json.dumps(packet))
        (attempt / 'failure.json').write_text(json.dumps({'status': 'CANNOT VERIFY',
            'error': 'ValueError: software fixture validation failure', 'elapsed_seconds': 1.5}))
        api = run / 'api'; api.mkdir()
        request = {'model': 'gpt-5.6-luna', 'messages': [
            {'role': 'system', 'content': 'Software fixture'}, {'role': 'user', 'content': 'Source'}],
            'max_completion_tokens': 100, 'reasoning_effort': 'low', 'store': False}
        feedback = {'schema': 'clockbase-design-validation-feedback-v1', 'validation_status': 'FAIL',
                    'error': 'Unsupported source quotation', 'citation_issue_count': 1,
                    'citation_issues': [{'path': '$.evidence[0]', 'error': 'Unsupported source quotation'}],
                    'citation_issues_omitted': 0, 'instruction': 'Software feedback'}
        previous = None
        for index in range(responses):
            slot = api / f'design-{index:03d}'; slot.mkdir()
            if index == 0:
                current = request
            else:
                current = _failed_retry_request(previous[0], previous[1], previous[2])
            response = {'model': 'gpt-5.6-luna', 'choices': [{'message': {
                'content': json.dumps({'invalid': True, 'index': index})}}]}
            reservation = {'model': 'gpt-5.6-luna', 'reserved_usd': .2, 'started_unix': 1, 'automatic_retries': 0, 'concurrency': 1}
            receipt = {'response_sha256': None,
                       'requested_model': 'gpt-5.6-luna', 'returned_model': 'gpt-5.6-luna',
                       'estimated_cost_usd': .1, 'billed_cost_usd': None, 'elapsed_seconds': 2.}
            (slot / 'request.json').write_text(json.dumps(current, ensure_ascii=False))
            (slot / 'reservation.json').write_text(json.dumps(reservation))
            (slot / 'response.json').write_text(json.dumps(response))
            receipt['response_sha256'] = sha(slot / 'response.json')
            (slot / 'receipt.json').write_text(json.dumps(receipt))
            (attempt / 'response-validation' / f'{index:03d}-response.txt').write_text(response['choices'][0]['message']['content'])
            (attempt / 'response-validation' / f'{index:03d}-validation.json').write_text(json.dumps(feedback))
            previous = (current, response, feedback)
        return run

    def test_failed_trace_exports_without_downstream_claims_and_reconstructs_retry(self):
        run = self._failed_fixture(responses=1); plan = self.root / 'failed-plan.json'; output = self.root / 'failed-public'
        self.assertEqual(prepare_failed(run, plan)['api_responses'], 1)
        result = export_failed(plan, output)
        self.assertEqual(result['next_retry_transmitted'], False)
        trace = json.loads((output / 'audit_trace.json').read_text())
        self.assertEqual(trace['schema'], 'clockbase-audit-trace-v1')
        self.assertEqual({row['model_assessment']['status'] for row in trace['domains']}, {'CANNOT VERIFY'})
        self.assertNotIn('execution', trace); self.assertNotIn('data', trace); self.assertNotIn('analysis', trace)
        self.assertNotIn('DESIGN_PROPOSED', json.dumps(trace))
        history = json.loads((output / 'public/failure-history.json').read_text())
        self.assertEqual(history['attempts'][0]['request_bytes'], 202)
        self.assertEqual(history['next_retry']['prompt_byte_limit'], 48000)
        self.assertFalse('/' + str(self.root).lstrip('/') in json.dumps(trace))

    def test_failed_trace_supports_sent_retry_and_resume_integrity(self):
        run = self._failed_fixture(responses=2); plan = self.root / 'retry-plan.json'; output = self.root / 'retry-public'
        prepare_failed(run, plan); first = export_failed(plan, output)
        self.assertTrue(first['next_retry_transmitted'])
        history = json.loads((output / 'public/failure-history.json').read_text())
        self.assertEqual(len(history['attempts']), 2); self.assertIsNone(history['next_retry'])
        resumed = export_failed(plan, output); self.assertTrue(resumed['resume'])
        (output / 'public/failure-history.json').write_text('changed')
        with self.assertRaises(PipelineIntegrityError): export_failed(plan, output)

    @staticmethod
    def _rewrite_failed_feedback(run, feedback):
        """Update only labeled software fixture transport; no empirical response repair."""
        previous = None
        for index, slot in enumerate(sorted((run/'api').iterdir())):
            request_path = slot/'request.json'
            request = json.loads(request_path.read_text())
            response = json.loads((slot/'response.json').read_text())
            if previous is not None:
                request = _failed_retry_request(previous[0], previous[1], feedback)
                request_path.write_text(json.dumps(request, ensure_ascii=False))
            path = run/f'03_blinded_design/attempt-000/response-validation/{index:03d}-validation.json'
            path.write_text(json.dumps(feedback))
            previous = request, response

    @staticmethod
    def _v3_coverage_feedback():
        return {'schema':'clockbase-design-validation-feedback-v3', 'validation_status':'FAIL',
            'error':'Assigned samples lack source quotations', 'citation_issue_count':0,
            'citation_issues':[], 'citation_issues_omitted':0, 'instruction':'Software fixture feedback',
            'coverage_issue_count':1, 'coverage_issues_omitted':0, 'coverage_issues':[
                {'path':'$.comparisons[0].case_ids', 'kind':'ASSIGNED_SAMPLE_QUOTATION_MISSING',
                 'missing_sample_ids':['GSM999001'], 'missing_count':1, 'omitted_count':0}]}

    def test_v3_failed_history_projects_bounded_coverage_and_preserves_retry_binding(self):
        run = self._failed_fixture(responses=2)
        feedback = self._v3_coverage_feedback()
        self._rewrite_failed_feedback(run, feedback)
        plan, output = self.root/'v3-plan.json', self.root/'v3-public'
        prepare_failed(run, plan); result = export_failed(plan, output)
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        history = json.loads((output/'public/failure-history.json').read_text())
        for index, attempt in enumerate(history['attempts']):
            self.assertEqual(attempt['validation_feedback_schema'], feedback['schema'])
            for key in ('coverage_issue_count', 'coverage_issues', 'coverage_issues_omitted'):
                self.assertEqual(attempt[key], feedback[key])
            self.assertNotIn('instruction', attempt)
            self.assertNotIn('response', attempt)
            self.assertNotIn('invalid', json.dumps(attempt))
            source_path = f'03_blinded_design/attempt-000/response-validation/{index:03d}-validation.json'
            self.assertEqual(history['source_sha256'][source_path], sha(run/source_path))
        self.assertTrue(history['attempts'][0]['retry_transmitted_after_this_response'])
        self.assertFalse(history['attempts'][1]['retry_transmitted_after_this_response'])
        self.assertIsNone(history['next_retry'])
        self.assertTrue(export_failed(plan, output)['resume'])
        feedback['coverage_issues'][0]['path'] = '$.comparisons[0].control_ids'
        self._rewrite_failed_feedback(run, feedback)
        with self.assertRaises(PipelineIntegrityError):
            export_failed(plan, self.root/'changed-public')

    def test_v3_failed_history_rejects_unbounded_or_malformed_coverage_before_and_after_plan(self):
        mutations = ('count_bool', 'count_mismatch', 'omitted_negative', 'too_many_issues',
                     'missing_bool', 'missing_mismatch', 'missing_omitted', 'duplicate_ids',
                     'too_many_ids', 'unknown_kind', 'nonobject', 'quote_payload', 'unowned_sample',
                     'unowned_field', 'sample_id_field', 'nonstring_field', 'arbitrary_path', 'wrong_kind_path')
        for phase in ('before', 'after'):
            for mutation in mutations:
                with self.subTest(phase=phase, mutation=mutation):
                    old = self.root; self.root = old/(phase+'-'+mutation); self.root.mkdir()
                    try:
                        run = self._failed_fixture(); plan = self.root/'plan.json'
                        feedback = self._v3_coverage_feedback(); self._rewrite_failed_feedback(run, feedback)
                        if phase == 'after': prepare_failed(run, plan)
                        item = feedback['coverage_issues'][0]
                        if mutation == 'count_bool': feedback['coverage_issue_count'] = True
                        elif mutation == 'count_mismatch': feedback['coverage_issue_count'] = 2
                        elif mutation == 'omitted_negative': feedback['coverage_issues_omitted'] = -1
                        elif mutation == 'too_many_issues':
                            feedback['coverage_issues'] *= 21; feedback['coverage_issue_count'] = 21
                        elif mutation == 'missing_bool': item['missing_count'] = True
                        elif mutation == 'missing_mismatch': item['missing_count'] = 2
                        elif mutation == 'missing_omitted': item['omitted_count'] = 1
                        elif mutation == 'duplicate_ids':
                            item['missing_sample_ids'] *= 2; item['missing_count'] = 2
                        elif mutation == 'too_many_ids':
                            item['missing_sample_ids'] = [f'GSM{i:06d}' for i in range(21)]; item['missing_count'] = 21
                        elif mutation == 'unknown_kind': item['kind'] = 'SCIENTIFIC_ACCEPTANCE'
                        elif mutation == 'nonobject': feedback['coverage_issues'][0] = None
                        elif mutation == 'quote_payload': item['quote'] = 'Raw content must not be projected'
                        elif mutation == 'unowned_sample': item['missing_sample_ids'] = ['SOFTWARE_RESPONSE_TEXT_SENTINEL']
                        elif mutation in ('unowned_field', 'sample_id_field', 'nonstring_field'):
                            item.update(kind='COVARIATE_FIELD_QUOTATION_MISSING', path='$.covariate_plan[0].source_fields')
                            item.pop('missing_sample_ids')
                            item['missing_source_fields'] = [{'unowned_field':'SOFTWARE_RESPONSE_TEXT_SENTINEL',
                                'sample_id_field':'sample_id', 'nonstring_field':'nonstring'}[mutation]]
                            if mutation == 'nonstring_field':
                                packet_path = run/'03_blinded_design/attempt-000/source-packet.json'
                                packet = json.loads(packet_path.read_text()); packet['samples'][0]['nonstring'] = {'nested':'value'}
                                packet_path.write_text(json.dumps(packet))
                        elif mutation == 'arbitrary_path': item['path'] = '$.SOFTWARE_RESPONSE_TEXT_SENTINEL'
                        elif mutation == 'wrong_kind_path': item['path'] = '$.covariate_plan[0].source_fields'
                        self._rewrite_failed_feedback(run, feedback)
                        with self.assertRaisesRegex(ValueError, 'coverage issue'):
                            if phase == 'before': prepare_failed(run, plan)
                            else: export_failed(plan, self.root/'public')
                    finally:
                        self.root = old

    def test_v3_failed_history_retains_source_owned_condition_and_field_diagnostics(self):
        run = self._failed_fixture()
        feedback = self._v3_coverage_feedback()
        feedback['coverage_issues'] += [
            {'path':'$.comparisons[0].findings[0].condition', 'kind':'CONDITION_SAMPLE_QUOTATION_MISSING',
             'missing_sample_ids':['GSM999001'], 'missing_count':1, 'omitted_count':0},
            {'path':'$.covariate_plan[0].source_fields', 'kind':'COVARIATE_FIELD_QUOTATION_MISSING',
             'missing_source_fields':['treatment'], 'missing_count':1, 'omitted_count':0}]
        feedback['coverage_issue_count'] = 3
        self._rewrite_failed_feedback(run, feedback)
        plan, output = self.root/'plan.json', self.root/'public'
        prepare_failed(run, plan); export_failed(plan, output)
        history = json.loads((output/'public/failure-history.json').read_text())
        self.assertEqual(history['attempts'][0]['coverage_issues'], feedback['coverage_issues'])

    def test_v1_v2_failed_history_keep_legacy_projection_shape(self):
        for version in ('v1', 'v2'):
            with self.subTest(version=version):
                old = self.root; self.root = old/version; self.root.mkdir()
                try:
                    run = self._failed_fixture(responses=2)
                    feedback = self._v3_coverage_feedback()
                    feedback['schema'] = 'clockbase-design-validation-feedback-' + version
                    self._rewrite_failed_feedback(run, feedback)
                    plan, output = self.root/'plan.json', self.root/'public'
                    prepare_failed(run, plan); export_failed(plan, output)
                    history = json.loads((output/'public/failure-history.json').read_text())
                    for attempt in history['attempts']:
                        self.assertEqual(set(attempt), {'index', 'api_slot', 'request_bytes', 'prompt_byte_limit',
                            'prompt_limit_check', 'response_sha256', 'requested_model', 'returned_model',
                            'estimated_cost_usd', 'billed_cost_usd', 'elapsed_seconds', 'validation_status',
                            'validation_error', 'citation_issue_count', 'citation_issues',
                            'citation_issues_omitted', 'retry_transmitted_after_this_response'})
                finally:
                    self.root = old

    def test_failed_trace_rejects_changed_source_api_receipt_private_text_and_downstream_stage(self):
        for mutation in ('source', 'api', 'receipt', 'downstream'):
            with self.subTest(mutation=mutation):
                fixture_root = self.root / mutation; fixture_root.mkdir()
                original_root = self.root
                self.root = fixture_root
                run = self._failed_fixture(); plan = self.root / (mutation + '-plan.json')
                self.root = original_root
                prepare_failed(run, plan)
                if mutation == 'source':
                    p = run / '03_blinded_design/attempt-000/source-packet.json'
                    packet = json.loads(p.read_text()); packet['samples'][0]['private'] = '/Entropy/private'
                    p.write_text(json.dumps(packet));
                    with self.assertRaises(PipelineIntegrityError): export_failed(plan, self.root / (mutation + '-out'))
                elif mutation == 'api':
                    p = run / 'api/design-000/response.json'; p.write_text(p.read_text() + '\n')
                    with self.assertRaises(PipelineIntegrityError): export_failed(plan, self.root / (mutation + '-out'))
                elif mutation == 'receipt':
                    p = run / 'api/design-000/receipt.json'; receipt = json.loads(p.read_text()); receipt['estimated_cost_usd'] = .3
                    p.write_text(json.dumps(receipt))
                    with self.assertRaises(PipelineIntegrityError): export_failed(plan, self.root / (mutation + '-out'))
                else:
                    stage = run / '04_preprocessing'; stage.mkdir(); (stage / 'complete.json').write_text('{}')
                    with self.assertRaises(ValueError): prepare_failed(run, plan)

    def test_failed_trace_hardens_cap_identity_private_source_roster_and_retry(self):
        # These mutations happen before plan preparation so they exercise the
        # source/contract guards rather than only the frozen hash comparison.
        fixture_root = self.root / 'hardening'; fixture_root.mkdir()
        original_root = self.root; self.root = fixture_root
        run = self._failed_fixture(responses=1)
        self.root = original_root
        packet_path = run / '03_blinded_design/attempt-000/source-packet.json'
        packet = json.loads(packet_path.read_text())
        packet['samples'][0]['private_note'] = '/Entropy/private-before-plan'
        packet_path.write_text(json.dumps(packet))
        plan = fixture_root / 'private-plan.json'
        prepare_failed(run, plan)
        with self.assertRaises(ValueError): export_failed(plan, fixture_root / 'private-public')

        fixture_root = self.root / 'identity'; fixture_root.mkdir()
        original_root = self.root; self.root = fixture_root
        run = self._failed_fixture(responses=1)
        self.root = original_root
        packet_path = run / '03_blinded_design/attempt-000/source-packet.json'
        packet = json.loads(packet_path.read_text()); packet['source_sha256'] = 'c' * 64
        (run / '03_blinded_design/attempt-000/source-transport.json').write_text(json.dumps(packet))
        packet_path.write_text(json.dumps(packet))
        with self.assertRaises(PipelineIntegrityError): prepare_failed(run, fixture_root / 'source-mismatch-plan.json')

        fixture_root = self.root / 'roster'; fixture_root.mkdir()
        original_root = self.root; self.root = fixture_root
        run = self._failed_fixture(responses=1)
        self.root = original_root
        packet_path = run / '03_blinded_design/attempt-000/source-packet.json'
        packet = json.loads(packet_path.read_text()); packet['samples'].append({'sample_id': 'GSM999002'})
        packet_path.write_text(json.dumps(packet))
        with self.assertRaises(ValueError): prepare_failed(run, fixture_root / 'roster-plan.json')

        fixture_root = self.root / 'cap'; fixture_root.mkdir()
        original_root = self.root; self.root = fixture_root
        run = self._failed_fixture(responses=1)
        self.root = original_root
        with self.assertRaises(ValueError): prepare_failed(run, fixture_root / 'cap-plan.json', prompt_bytes=96000)
        plan = fixture_root / 'normal-plan.json'; prepare_failed(run, plan)
        saved = json.loads(plan.read_text()); saved['prompt_byte_limit'] = 96000
        plan.write_text(json.dumps(saved))
        with self.assertRaises(PipelineIntegrityError): export_failed(plan, fixture_root / 'cap-public')

        fixture_root = self.root / 'retry-tamper'; fixture_root.mkdir()
        original_root = self.root; self.root = fixture_root
        run = self._failed_fixture(responses=2)
        self.root = original_root
        retry_path = run / 'api/design-001/request.json'
        retry = json.loads(retry_path.read_text()); retry['messages'][-1]['content'] += ' tampered'
        retry_path.write_text(json.dumps(retry, ensure_ascii=False))
        with self.assertRaises(PipelineIntegrityError): prepare_failed(run, fixture_root / 'retry-plan.json')

        fixture_root = self.root / 'over-cap'; fixture_root.mkdir()
        original_root = self.root; self.root = fixture_root
        run = self._failed_fixture(responses=1)
        self.root = original_root
        config_path = run / 'config.json'; config = json.loads(config_path.read_text()); config['max_prompt_bytes'] = 100
        config_path.write_text(json.dumps(config))
        with self.assertRaises(ValueError): prepare_failed(run, fixture_root / 'over-cap-plan.json')


class LiteralFieldFailureTests(unittest.TestCase):
    """Constructed source/API files exercise integrity, not model or study accuracy."""
    _failed_fixture = PublicTraceTests._failed_fixture

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    @staticmethod
    def _write(path, value):
        path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

    @staticmethod
    def _design(kind='parse'):
        from audit_system.agent_workflow_2026_09_10 import AUDIT_DOMAINS, DOMAIN_DESIGN_SCHEMA
        value = {'schema': DOMAIN_DESIGN_SCHEMA, 'study_id': 'GSE999001',
            'status': 'CANNOT VERIFY', 'comparisons': [], 'unassigned_ids': ['GSM999001'],
            'findings': [], 'unresolved': ['UNPUBLISHED_SOFTWARE_RESPONSE_SENTINEL'],
            'evidence': [{'sample_id': 'GSM999001', 'field': 'title', 'quote_from_field': True}],
            'audit_domains': {name: {'status': 'CANNOT VERIFY', 'summary': 'Software fixture only',
                'evidence': [], 'unresolved': ['No empirical evidence'], 'proposed_repair': None,
                'verification_needed': ['Independent evidence']} for name in AUDIT_DOMAINS},
            'covariate_plan': [{'comparison_name': '*', 'name': 'Software metadata',
                'source_fields': ['title'] if kind != 'valid' else [], 'decision': 'unknown',
                'rationale': 'Software evidence-coverage fixture', 'evidence': [],
                'missingness_action': 'Remain unresolved'}]}
        if kind == 'materialization':
            # The first field is valid; failure of the second must leave no saved
            # candidate and no partially returned change list.
            value['evidence'].append({'sample_id': 'GSM999001', 'field': 'title',
                                      'quote_from_field': True, 'quote': 'ambiguous'})
        return value

    def _literal_fixture(self, folder, kind='parse'):
        import shutil
        from types import SimpleNamespace
        from langchain_core.messages import AIMessage
        from audit_system.agent_workflow_2026_09_10 import source_packet
        from audit_system.public_trace_2026_09_18 import _export_identity
        from src.pipeline_2026_09_10 import assess_design
        folder.mkdir()
        old = self.root
        try:
            self.root = folder
            run = self._failed_fixture()
        finally:
            self.root = old
        metadata = run / '01_metadata/attempt-000/metadata.csv'
        metadata.write_text('GSM,title,source_name_ch1\nGSM999001,Software title,Software source\n')
        complete = run / '01_metadata/complete.json'
        record = json.loads(complete.read_text())
        record['outputs']['attempt-000/metadata.csv'] = sha(metadata)
        self._write(complete, record)
        (complete.parent / 'complete.sha256').write_text(sha(complete) + '\n')
        config = {'study':'GSE999001', 'target':'Software fixture target', 'max_prompt_bytes':48000,
            'model':'gpt-6-luna', 'audit_model':'gpt-6-luna', 'max_output_tokens':1000,
            'audit_profile':'source_domains_v4', 'source_binding_profile':'legacy',
            'source_contract_policy':'literal_field_refs_v1', 'audit_validation_retries':0,
            'audit_input_usd_per_million':.2, 'audit_output_usd_per_million':.75,
            'implementation':_export_identity()['implementation_sha256']}
        self._write(run/'config.json', config)
        packet = source_packet(metadata, config['study'], config['target'])
        attempt = run/'03_blinded_design/attempt-000'
        shutil.rmtree(attempt); attempt.mkdir()
        shutil.rmtree(run/'api'); (run/'api').mkdir()
        raw = json.dumps(self._design(kind), ensure_ascii=False)
        calls = []
        def invoke(messages):
            calls.append(messages)
            slot = run/'api/design-000'; slot.mkdir()
            wire = [{'role':'system','content':messages[0].content}, {'role':'user','content':messages[1].content}]
            request = {'model':config['audit_model'],'messages':wire,'max_completion_tokens':1000,
                       'reasoning_effort':'low','store':False,'response_format':{'type':'json_object'}}
            usage = {'prompt_tokens':10,'completion_tokens':5,'total_tokens':15}
            response = {'model':config['audit_model'],'choices':[{'message':{'role':'assistant','content':raw},
                                                                  'finish_reason':'stop'}],'usage':usage}
            self._write(slot/'request.json', request); self._write(slot/'response.json', response)
            reserved = ((len(json.dumps(request, ensure_ascii=False).encode())+200)*.2+1000*.75)/1e6
            self._write(slot/'reservation.json', {'reserved_usd':reserved,'started_unix':1.,
                'model':config['audit_model'],'automatic_retries':0,'concurrency':1})
            self._write(slot/'receipt.json', {'usage':usage,'estimated_cost_usd':(10*.2+5*.75)/1e6,
                'billed_cost_usd':None,'response_sha256':sha(slot/'response.json'),
                'requested_model':config['audit_model'],'returned_model':config['audit_model'],'elapsed_seconds':.5})
            return AIMessage(content=raw)
        try:
            assess_design(packet, SimpleNamespace(invoke=invoke), attempt,
                review_profile='source_domains_v4', source_binding_profile='legacy',
                source_contract_policy='literal_field_refs_v1', validation_retries=0)
        except (ValueError, TypeError, KeyError) as error:
            self._write(attempt/'failure.json', {'status':'CANNOT VERIFY',
                'error':f'{type(error).__name__}: {error}','elapsed_seconds':1.})
        else:
            self.fail('Rejection fixture unexpectedly passed its source contract')
        self.assertEqual(len(calls), 1)
        return run

    def _reseal_adjustment(self, run, *, raw=None):
        """Deliberately adversarial software reseal, not an empirical repair."""
        from src.pipeline_2026_09_10 import source_contract_candidate, source_contract_visible_packet
        attempt = run/'03_blinded_design/attempt-000'; slot = run/'api/design-000'
        adjust = attempt/'contract-adjustments'; response = json.loads((slot/'response.json').read_text())
        if raw is not None:
            response['choices'][0]['message']['content'] = raw
            self._write(slot/'response.json', response)
            (adjust/'000-original-response.txt').write_text(raw)
        receipt = json.loads((slot/'receipt.json').read_text()); receipt['response_sha256']=sha(slot/'response.json')
        self._write(slot/'receipt.json', receipt)
        packet, transport = [json.loads((attempt/name).read_text()) for name in ['source-packet.json','source-transport.json']]
        value, changes = source_contract_candidate(response['choices'][0]['message']['content'],
            source_contract_visible_packet(packet, transport), 'literal_field_refs_v1')
        self._write(adjust/'000-candidate.json', value)
        record=json.loads((adjust/'000-receipt.json').read_text())
        record.update(changes=changes,candidate_sha256=sha(adjust/'000-candidate.json'),
            original_response_sha256=sha(adjust/'000-original-response.txt'),
            source_packet_sha256=sha(attempt/'source-packet.json'),source_transport_sha256=sha(attempt/'source-transport.json'))
        self._write(adjust/'000-receipt.json', record)
        request=json.loads((slot/'request.json').read_text())
        reservation=json.loads((slot/'reservation.json').read_text())
        reservation['reserved_usd']=((len(json.dumps(request,ensure_ascii=False).encode())+200)*.2+1000*.75)/1e6
        self._write(slot/'reservation.json', reservation)

    def _mutate(self, run, mutation):
        from audit_system.agent_workflow_2026_09_10 import review_messages
        from src.pipeline_2026_09_10 import source_contract_instructions
        attempt=run/'03_blinded_design/attempt-000'; adjust=attempt/'contract-adjustments'; slot=run/'api/design-000'
        def change(path, update):
            value=json.loads(path.read_text()); update(value); self._write(path,value)
        receipt=adjust/'000-receipt.json'
        if mutation.startswith('missing_'):
            name={'missing_original':'000-original-response.txt','missing_candidate':'000-candidate.json',
                  'missing_adjustment':'000-receipt.json'}[mutation]
            (adjust/name).unlink()
        elif mutation == 'candidate_resealed':
            change(adjust/'000-candidate.json', lambda v:v['unresolved'].append('Adversarial changed candidate'))
            change(receipt, lambda v:v.update(candidate_sha256=sha(adjust/'000-candidate.json')))
        elif mutation == 'changes_omitted':
            change(receipt,lambda v:v.update(changes=[]))
        elif mutation in {'receipt_extra','schema','policy','receipt_error','accepted_numeric','attempt_boolean','normalization_limits','implementation_hash'}:
            key,value={'receipt_extra':('extra',True),'schema':('schema','unknown'),'policy':('policy','strict'),
                'receipt_error':('error','ValueError: fabricated failure'),'accepted_numeric':('accepted_result',0),
                'attempt_boolean':('attempt',False),'normalization_limits':('normalization_limits',{}),
                'implementation_hash':('implementation_sha256','0'*64)}[mutation]
            change(receipt,lambda v:v.update({key:value}))
        elif mutation == 'failure_error':
            change(attempt/'failure.json',lambda v:v.update(error='ValueError: fabricated failure'))
        elif mutation == 'failure_status':
            change(attempt/'failure.json',lambda v:v.update(status='PASS'))
        elif mutation == 'config_code':
            change(run/'config.json',lambda v:v['implementation'].update({'src/pipeline_2026_09_10.py':'0'*64}))
        elif mutation == 'config_retry':
            change(run/'config.json',lambda v:v.update(audit_validation_retries=1))
        elif mutation == 'config_policy':
            change(run/'config.json',lambda v:v.update(source_contract_policy='unknown'))
        elif mutation == 'request_role':
            change(slot/'request.json',lambda v:v['messages'][1].update(role='assistant'))
        elif mutation == 'request_extra':
            change(slot/'request.json',lambda v:v.update(tools=[]))
        elif mutation == 'request_option':
            change(slot/'request.json',lambda v:v.update(store=True))
        elif mutation == 'request_model':
            change(run/'config.json',lambda v:v.update(audit_model='gpt-5.6-luna'))
        elif mutation == 'request_payload':
            change(slot/'request.json',lambda v:v['messages'][1].update(content='{}'))
        elif mutation == 'response_finish':
            change(slot/'response.json',lambda v:v['choices'][0].update(finish_reason='length'))
            change(slot/'receipt.json',lambda v:v.update(response_sha256=sha(slot/'response.json')))
        elif mutation == 'usage':
            change(slot/'receipt.json',lambda v:v['usage'].update(prompt_tokens=11))
        elif mutation == 'cost':
            change(slot/'receipt.json',lambda v:v.update(estimated_cost_usd=.1))
        elif mutation == 'original_resealed':
            (adjust/'000-original-response.txt').write_text('changed')
            change(receipt,lambda v:v.update(original_response_sha256=sha(adjust/'000-original-response.txt')))
        elif mutation == 'mixed_validation':
            (attempt/'response-validation').mkdir()
        elif mutation == 'success_result':
            (attempt/'design.json').write_text('{}')
        elif mutation == 'stage_success':
            (attempt.parent/'complete.json').write_text('{}')
        elif mutation == 'downstream':
            (run/'04_preprocessing').mkdir()
        elif mutation == 'extra_upstream':
            (run/'01_metadata/extra.txt').write_text('Unlisted software fixture')
        elif mutation == 'extra_api_directory':
            (slot/'extra').mkdir()
        elif mutation == 'symlink':
            (adjust/'extra-link').symlink_to(receipt)
        elif mutation == 'typed_duplicate_key':
            text=receipt.read_text();receipt.write_text(text.replace('"attempt": 0','"attempt": 0, "attempt": 0'))
        elif mutation == 'source_substitution':
            change(attempt/'source-packet.json',lambda v:v['samples'][0].update(title='Substituted software title'))
            packet=json.loads((attempt/'source-packet.json').read_text())
            messages=review_messages(packet,review_profile='source_domains_v4',source_binding_profile='legacy')
            self._write(attempt/'source-transport.json',json.loads(messages[1].content))
            change(slot/'request.json',lambda v:v.update(messages=[
                {'role':'system','content':messages[0].content+source_contract_instructions('literal_field_refs_v1')},
                {'role':'user','content':messages[1].content}]))
            self._reseal_adjustment(run)
        elif mutation == 'valid_relabel':
            self._reseal_adjustment(run,raw=json.dumps(self._design('valid'),ensure_ascii=False))
        else:
            raise AssertionError('Unknown software mutation: '+mutation)

    def test_zero_retry_parser_rejection_has_separate_diagnostic_and_no_invented_retry(self):
        run=self._literal_fixture(self.root/'parse');plan=self.root/'plan.json';output=self.root/'public'
        prepared=prepare_failed(run,plan);self.assertEqual(prepared['validation_records'],0)
        result=export_failed(plan,output);self.assertFalse(result['next_retry_transmitted'])
        history=json.loads((output/'public/failure-history.json').read_text());attempt=history['attempts'][0]
        self.assertEqual(attempt['recorded_rejection']['contract_status'],'REJECTED')
        self.assertEqual(attempt['export_time_diagnostic']['rejection_step'],'strict_parse')
        self.assertEqual(attempt['recorded_rejection']['recorded_change_count'],1)
        self.assertIsNone(history['next_retry']);self.assertEqual(history['response_validation_records'],0)
        self.assertNotIn('validation_status',attempt);self.assertNotIn('citation_issue_count',attempt)
        for path in output.rglob('*.json'):
            self.assertNotIn('UNPUBLISHED_SOFTWARE_RESPONSE_SENTINEL',path.read_text())
        trace=json.loads((output/'audit_trace.json').read_text())
        self.assertEqual(trace['overall']['status'],'CANNOT VERIFY');self.assertEqual(trace['comparisons'],[])
        self.assertEqual({x['model_assessment']['status'] for x in trace['domains']},{'CANNOT VERIFY'})
        self.assertTrue(export_failed(plan,output)['resume'])

    def test_zero_retry_materialization_rejection_has_no_candidate_or_partial_changes(self):
        run=self._literal_fixture(self.root/'materialization','materialization')
        plan=self.root/'plan.json';output=self.root/'public';prepare_failed(run,plan);export_failed(plan,output)
        attempt=json.loads((output/'public/failure-history.json').read_text())['attempts'][0]
        self.assertEqual(attempt['export_time_diagnostic']['rejection_step'],'materialization')
        self.assertIsNone(attempt['recorded_rejection']['candidate_sha256'])
        self.assertEqual(attempt['recorded_rejection']['recorded_change_count'],0)
        self.assertFalse((run/'03_blinded_design/attempt-000/contract-adjustments/000-candidate.json').exists())

    def test_rejection_mutations_fail_before_and_after_preparation(self):
        mutations=['missing_original','missing_candidate','missing_adjustment','candidate_resealed','changes_omitted',
            'receipt_extra','schema','policy','receipt_error','failure_error','failure_status','accepted_numeric',
            'attempt_boolean','normalization_limits','implementation_hash','config_code','config_retry','config_policy',
            'request_role','request_extra','request_option','request_model','request_payload','response_finish','usage',
            'cost','original_resealed','mixed_validation','success_result','stage_success','downstream',
            'extra_upstream','extra_api_directory','symlink','typed_duplicate_key','source_substitution','valid_relabel']
        for phase in ['before','after']:
            for mutation in mutations:
                with self.subTest(phase=phase,mutation=mutation):
                    root=self.root/(phase+'-'+mutation);run=self._literal_fixture(root)
                    plan=root/'plan.json';output=root/'public'
                    if phase=='after':prepare_failed(run,plan)
                    self._mutate(run,mutation)
                    with self.assertRaises((ValueError,PipelineIntegrityError)):
                        if phase=='before':prepare_failed(run,plan)
                        else:export_failed(plan,output)
                    self.assertFalse(output.exists())

    def test_materialization_failure_refuses_spurious_candidate_or_changes(self):
        for name in ['candidate','changes']:
            with self.subTest(mutation=name):
                root=self.root/name;run=self._literal_fixture(root,'materialization')
                adjust=run/'03_blinded_design/attempt-000/contract-adjustments'
                if name=='candidate':self._write(adjust/'000-candidate.json',{})
                else:
                    value=json.loads((adjust/'000-receipt.json').read_text());value['changes']=[{}]
                    self._write(adjust/'000-receipt.json',value)
                with self.assertRaises((ValueError,PipelineIntegrityError)):
                    prepare_failed(run,root/'plan.json')

    def test_plan_binds_current_normalizer_and_final_output_bytes(self):
        run=self._literal_fixture(self.root/'binding');plan=self.root/'plan.json';output=self.root/'public'
        prepare_failed(run,plan);value=json.loads(plan.read_text())
        self.assertIn('src/pipeline_2026_09_10.py',value['export_identity']['implementation_sha256'])
        export_failed(plan,output)
        (output/'public/failure-history.json').write_text('{}')
        with self.assertRaises(PipelineIntegrityError):export_failed(plan,output)


class CorpusProjectionTests(unittest.TestCase):
    """Constructed receipt/CSV fixtures exercise export, not empirical science."""
    def setUp(self):
        import gzip, hashlib
        from audit_system.legacy_analysis_review_2026_09_19 import SCIENTIFIC_CHECKS
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name); self.run = self.root/'run'; self.run.mkdir()
        self.study = self.run/'studies/human_rna-GSE1'; self.study.mkdir(parents=True)
        self.csv = self.root/'source.csv'; self.csv.write_text('gse_id,clock,p_value\nGSE1,clock,0.5\n')
        self.row_id = 'R'+hashlib.sha256(json.dumps(['source.csv',sha(self.csv),2],sort_keys=True,allow_nan=False).encode()).hexdigest()[:24]
        self.task = {'study_key':'human_rna:GSE1','study_id':'GSE1','modality':'human_rna','source_files':{}}
        queue=[self.task,{'study_key':'human_rna:GSE2','study_id':'GSE2','modality':'human_rna','source_files':{}}]
        with gzip.open(self.run/'queue.jsonl.gz','wt') as h:
            h.write(''.join(json.dumps(r)+'\n' for r in queue))
        self.write(self.run/'queue-receipt.json',{'sha256':sha(self.run/'queue.jsonl.gz'),'studies':2})
        self.write(self.run/'config.json',{'denominators':{'parsed_logical_rows':4},'implementation_sha256':{'fixture.py':'a'*64},
            'config':{'provider':'software_fixture','api':{'model':'none'},'concurrency':1},'source_scope':'Disclosed software fixture'})
        (self.run/'config.sha256').write_text(sha(self.run/'config.json'))
        self.write(self.study/'config.json',{'task':self.task,'artifacts':[{'path':'source.csv','sha256':sha(self.csv),'artifact_id':'fixture'}]})
        self.stage('01_source_index',{'evidence':{'S00000':{'artifact':'GEO','field':'treatment','value':'fixture','sample_ids':['GSM1']}}})
        self.stage('02_blinded_design',{'study_id':'GSE1','evidence_ids':['S00000'],'comparisons':[],
            'limitations':['Private source at /Entropy/secret/input.csv'],'scientific_acceptance':'CANNOT VERIFY'})
        self.stage('03_conditional_reruns',[{'comparison_id':'fixture','clock':'clock','case_ids':['GSM1'],'control_ids':['GSM2'],
            'status':'CANNOT VERIFY','reason':'Software fixture has no measurements'}])
        checks={k:{'status':'CANNOT VERIFY','reason':'Fixture lacks evidence','evidence_ids':[],'unresolved':['Fixture'],'scope':'fixture'} for k in SCIENTIFIC_CHECKS}
        self.aggregate={'rows':[{'row_id':self.row_id,'status':'CANNOT VERIFY','action':'unresolved','evidence_ids':['N00000']}],
            'scoped_scientific_checks':[{'review_stage':'04_review_0000','row_ids':[self.row_id],'checks':checks}], 'design_disagreements':[]}
        self.stage('05_aggregation',self.aggregate)
        self.terminal={'study_key':self.task['study_key'],'status':'REVIEWED','reviewed_rows':1,'output_sha256':{}}
        self.seal()
        self.report={'selected_studies':1,'whole_corpus_studies':2,'results':[self.terminal],
            'accepted_scientific_repairs':0,'scientific_error_rate':None,'unique_reviewed_rows':1,'calculations':0,
            'created_at_utc':'2026-09-20T00:00:00Z','scope':'Disclosed software fixture','budget':{}}
        self.report_path=self.run/'reports/report-0000.json'; self.write(self.report_path,self.report)
    def write(self,path,value):
        path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps(value,sort_keys=True))
    def stage(self,name,value):
        d=self.study/name; self.write(d/'result.json',value)
        self.write(d/'complete.json',{'result':value,'outputs':{'result.json':sha(d/'result.json')}})
        (d/'complete.sha256').write_text(sha(d/'complete.json'))
    def seal(self):
        self.terminal['output_sha256']={str(p.relative_to(self.study)):sha(p) for p in self.study.rglob('*') if p.is_file() and p.name not in ('terminal.json','terminal.sha256')}
        self.write(self.study/'terminal.json',self.terminal); (self.study/'terminal.sha256').write_text(sha(self.study/'terminal.json'))
    def export(self):
        from audit_system.public_trace_2026_09_18 import export_corpus
        return export_corpus(self.run,self.root/'public',source_root=self.root)
    def test_exact_denominators_scopes_identity_and_sanitization(self):
        self.export(); a=json.loads((self.root/'public/assessment.json').read_text())
        self.assertEqual(a['summary']['attempted_study_tasks'],1); self.assertEqual(a['summary']['pending_study_tasks'],1)
        self.assertEqual(a['summary']['rows_without_completed_review'],3)
        p=json.loads((self.root/'public/studies/human_rna-GSE1.json').read_text())
        self.assertEqual(len(p['review']['scoped_scientific_checks'][0]['checks']),16)
        self.assertEqual(p['original_rows'][0]['row_id'],self.row_id)
        self.assertEqual(p['review_evidence_bindings'][0]['conditional_evidence']['N00000']['calculation_index'],0)
        self.assertIn('Text withheld',p['frozen_source_design']['limitations'][0])
        self.assertNotIn('/Entropy/',json.dumps(p)); self.assertEqual(p['accepted_scientific_repairs'],0)
        manifest=json.loads((self.root/'public/artifact-manifest.json').read_text())
        for relative,digest in manifest['files'].items(): self.assertEqual(sha(self.root/'public'/relative),digest)
        with self.assertRaises(ValueError): self.export()

    def interrupt_fixture(self):
        """A software lifecycle fixture; it supplies no empirical measurements."""
        import shutil
        (self.study/'terminal.json').unlink()
        (self.study/'terminal.sha256').unlink()
        shutil.rmtree(self.study/'05_aggregation')
        failure = {'study_key': self.task['study_key'], 'status': 'CONTROLLER_FAILURE',
                   'error': 'Disclosed software transport-interruption fixture'}
        self.write(self.study/'04_review_0000/attempt-000/failure.json',
                   {'status': 'FAIL', 'error': failure['error'], 'elapsed_seconds': 1})
        self.report.update(results=[failure], unique_reviewed_rows=0, calculations=0,
                           unresolved_calculations=0)
        self.write(self.report_path, self.report)
        return failure

    def test_interrupted_stops_remain_outside_terminal_totals(self):
        self.interrupt_fixture()
        self.export()
        a=json.loads((self.root/'public/assessment.json').read_text())
        p=json.loads((self.root/'public/studies/human_rna-GSE1.json').read_text())
        self.assertEqual((a['summary']['conditional_calculations'], a['summary']['unresolved_calculations']), (0, 0))
        self.assertEqual((p['calculations'], p['independent_numerics'], p['review']['rows']), ([], None, []))
        evidence=p['intermediate_evidence']
        self.assertEqual(evidence['source_report_sha256'], sha(self.report_path))
        self.assertEqual(evidence['calculations'][0]['reason'], 'Software fixture has no measurements')
        self.assertEqual(a['intermediate_evidence_summary']['unresolved_targets'], 1)
        self.assertEqual(a['intermediate_evidence_summary']['failed_attempts'], 1)
        self.assertEqual(evidence['completed_stages']['03_conditional_reruns']['receipt_sha256'],
                         sha(self.study/'03_conditional_reruns/complete.json'))

    def test_interrupted_task_without_calculation_stage_is_explicit(self):
        import shutil
        self.interrupt_fixture()
        shutil.rmtree(self.study/'03_conditional_reruns')
        self.export()
        a=json.loads((self.root/'public/assessment.json').read_text())
        p=json.loads((self.root/'public/studies/human_rna-GSE1.json').read_text())
        self.assertEqual(p['intermediate_evidence']['calculations'], [])
        self.assertEqual(a['intermediate_evidence_summary']['study_tasks'], 1)
        self.assertEqual(a['intermediate_evidence_summary']['unresolved_targets'], 0)

    def test_interrupted_completed_empty_stage_is_preserved(self):
        self.interrupt_fixture()
        self.stage('03_conditional_reruns', [])
        self.export()
        p=json.loads((self.root/'public/studies/human_rna-GSE1.json').read_text())
        self.assertEqual(p['intermediate_evidence']['calculations'], [])
        self.assertIn('03_conditional_reruns', p['intermediate_evidence']['completed_stages'])

    def test_mixed_completed_and_interrupted_tasks_keep_terminal_coverage(self):
        import shutil
        completed=self.study
        interrupted=self.run/'studies/human_rna-GSE2'
        shutil.copytree(completed, interrupted)
        original_task=self.task
        self.study=interrupted
        self.task={**original_task, 'study_key': 'human_rna:GSE2', 'study_id': 'GSE2'}
        config=json.loads((interrupted/'config.json').read_text());config['task']=self.task
        self.write(interrupted/'config.json', config)
        # PASS is only a lifecycle/counting flag in this software fixture.
        self.stage('03_conditional_reruns', [{'comparison_id':'fixture', 'clock':'clock',
            'case_ids':['GSM1'], 'control_ids':['GSM2'], 'status':'PASS',
            'scientific_status':'CANNOT VERIFY', 'meaning':'Software counting fixture; no measurements'}])
        failure=self.interrupt_fixture()
        self.report.update(selected_studies=2, results=[self.terminal, failure], unique_reviewed_rows=1,
                           calculations=0, unresolved_calculations=1)
        self.write(self.report_path,self.report)
        self.export()
        a=json.loads((self.root/'public/assessment.json').read_text())
        self.assertEqual(a['summary']['reviewed_unique_rows'], 1)
        self.assertEqual((a['summary']['conditional_calculations'], a['summary']['unresolved_calculations']), (0, 1))
        self.assertEqual(a['intermediate_evidence_summary']['calculated_targets'], 1)
        self.assertEqual(a['summary']['status_counts'], {'REVIEWED':1, 'CONTROLLER_FAILURE':1})

    def test_unmapped_interrupted_review_is_refused_without_discarding_it(self):
        self.interrupt_fixture()
        self.stage('04_review_0000', {'rows': [], 'scope': 'Disclosed unsupported software fixture'})
        before=sha(self.study/'04_review_0000/complete.json')
        with self.assertRaisesRegex(ValueError, 'Unsupported interrupted corpus stage'):
            self.export()
        self.assertEqual(sha(self.study/'04_review_0000/complete.json'), before)
    def test_independent_scientific_hold_remains_distinct_from_integrity_pass(self):
        from audit_system.public_trace_2026_09_18 import export_corpus
        review = {'run': {'report_sha256': sha(self.report_path), 'config_sha256': sha(self.run/'config.json'),
            'queue_sha256': sha(self.run/'queue.jsonl.gz')}, 'status': 'PASS', 'decision':'HOLD',
            'scientific_findings':['Disclosed software fixture: unsupported model FAIL remains unresolved.'],
            'source_ledger':{'path':'/Entropy/private'}}
        path=self.root/'independent.json';self.write(path,review)
        export_corpus(self.run,self.root/'public',source_root=self.root,independent_review_path=path)
        p=json.loads((self.root/'public/independent-review-public.json').read_text())
        self.assertEqual((p['status'],p['decision']),('PASS','HOLD')); self.assertNotIn('source_ledger',p)
        review['run']['report_sha256']='a'*64;self.write(path,review)
        with self.assertRaises(ValueError):
            export_corpus(self.run,self.root/'another',source_root=self.root,independent_review_path=path)

    def test_malformed_tail_retains_parsed_rows_and_explicit_issue(self):
        import hashlib
        self.csv.write_text(self.csv.read_text()+'"unterminated')
        config=json.loads((self.study/'config.json').read_text()); config['artifacts'][0]['sha256']=sha(self.csv)
        self.write(self.study/'config.json',config)
        rid='R'+hashlib.sha256(json.dumps(['source.csv',sha(self.csv),2],sort_keys=True,allow_nan=False).encode()).hexdigest()[:24]
        self.aggregate['rows'][0]['row_id']=rid
        self.aggregate['scoped_scientific_checks'][0]['row_ids']=[rid]
        self.stage('05_aggregation',self.aggregate);self.seal();self.write(self.report_path,self.report)
        self.export();p=json.loads((self.root/'public/studies/human_rna-GSE1.json').read_text())
        self.assertEqual(len(p['original_rows']),1)
        self.assertEqual(p['original_parse_issues'][0]['status'],'CANNOT VERIFY')
        self.assertEqual(p['original_rows'][0]['row_id'],rid)
        self.assertEqual(p['accepted_scientific_repairs'],0)

    def test_changed_original_source_is_refused(self):
        self.csv.write_text(self.csv.read_text()+'GSE1,clock,0.1\n')
        with self.assertRaises(PipelineIntegrityError): self.export()
    def test_changed_stage_or_terminal_is_refused(self):
        (self.study/'03_conditional_reruns/result.json').write_text('[]')
        with self.assertRaises(PipelineIntegrityError): self.export()
    def test_summary_and_scientific_inflation_are_refused(self):
        for key,value in [('unique_reviewed_rows',2),('calculations',1),('accepted_scientific_repairs',1)]:
            changed={**self.report,key:value}; self.write(self.report_path,changed)
            with self.assertRaises(ValueError): self.export()
    def test_missing_check_rejected_even_with_valid_hashes(self):
        self.aggregate['scoped_scientific_checks'][0]['checks'].pop('group_membership')
        self.stage('05_aggregation',self.aggregate); self.seal(); self.write(self.report_path,self.report)
        with self.assertRaises(ValueError): self.export()
    def test_partial_review_retains_unexecuted_repair_status(self):
        self.terminal['status']='PARTIALLY_REVIEWED'
        self.terminal['repair_summary']={'status':'CANNOT VERIFY','reason':'Incomplete model review'}
        self.seal(); self.write(self.report_path,self.report); self.export()
        p=json.loads((self.root/'public/studies/human_rna-GSE1.json').read_text())
        self.assertIsNone(p['repairs']); self.assertEqual(p['repair_stage_limit']['status'],'CANNOT VERIFY')

    def stopped_numeric_fixture(self, *, clock='N/A', defect=None):
        """Receipt-contract fixture for the observed stopped N/A export case."""
        calculation={'comparison_id':'source_target','clock':clock,'case_ids':['GSM1'],
            'control_ids':['GSM2'],'status':'CANNOT VERIFY','reason':'Clock unavailable; no calculation'}
        original={'comparison_id':'source_target','clock':clock,'input_index':0,
            'status':'CANNOT VERIFY','numerical_status':'CANNOT VERIFY',
            'scientific_status':'CANNOT VERIFY','checks':[],'reference_values':{}}
        if defect is not None:
            field,value=defect;original[field]=value
        self.stage('03_conditional_reruns',[calculation])
        self.stage('03b_independent_r',{'schema':'clockbase-independent-statistics-v1',
            'status':'CANNOT VERIFY','scientific_status':'CANNOT VERIFY','results':[original]})
        self.seal();self.write(self.report_path,self.report)
        return calculation,original

    def test_source_bound_stopped_na_clock_retains_exact_public_identity(self):
        self.stopped_numeric_fixture();self.export()
        p=json.loads((self.root/'public/studies/human_rna-GSE1.json').read_text())
        calculation=p['calculations'][0];numeric=p['independent_numerics']['results'][0]
        self.assertEqual((numeric['input_index'],numeric['comparison_id'],numeric['clock']),
                         (0,calculation['comparison_id'],calculation['clock']))
        self.assertEqual(numeric['clock'],'N/A')
        self.assertEqual(numeric['clock_label_status'],'SOURCE_BOUND_LITERAL_STOP_LABEL')
        self.assertEqual(numeric['numerical_status'],'CANNOT VERIFY')
        self.assertIsNone(calculation['calculation'])
        self.assertEqual(p['scientific_status'],'CANNOT VERIFY')

    def test_resealed_stopped_target_identity_mismatch_is_rejected(self):
        for defect in [('input_index',1),('input_index',False),('comparison_id','other_target'),('clock','other_clock')]:
            with self.subTest(defect=defect):
                self.stopped_numeric_fixture(defect=defect)
                with self.assertRaisesRegex(ValueError,'identity differs'):
                    self.export()

    def test_na_clock_cannot_be_promoted_to_numerical_pass(self):
        self.stopped_numeric_fixture(defect=('numerical_status','PASS'))
        with self.assertRaisesRegex(ValueError,'Unresolved clock label'):
            self.export()

    def test_other_unsupported_clock_labels_stay_withheld(self):
        self.stopped_numeric_fixture(clock='/Entropy/private/clock');self.export()
        p=json.loads((self.root/'public/studies/human_rna-GSE1.json').read_text())
        numeric=p['independent_numerics']['results'][0]
        self.assertIsNone(numeric['clock'])
        self.assertEqual(numeric['clock_label_status'],'WITHHELD_UNSAFE_OR_UNSUPPORTED_LABEL')
        self.assertNotIn('/Entropy/private/clock',json.dumps(p))


class FutureProjectionBindingTests(unittest.TestCase):
    """Constructed binding fixtures; no empirical clock or scientific result."""
    def setUp(self):
        from audit_system.public_trace_2026_09_18 import _future_terminal_result, _future_public_date, _future_producer_binding
        self.terminal_check = _future_terminal_result
        self.date_check = _future_public_date
        self.producer_check = _future_producer_binding
        self.terminal = {'scientific_acceptance':'CANNOT VERIFY','accepted_scientific_repairs':0,'input_rows':1,'scientific_error_rate':None}
        self.envelope = {**self.terminal,'execution_status':'COMPLETED','budget':{},'pipeline_root_sha256':'a'*64,'original_paper_2025_eligibility':'QUARANTINED'}
        self.packet = {'study_id':'GSE1','study_key':'human_rna:GSE1','independent_numerics':{'predictions_sha256':'b'*64},
            'original_rows':[{'identity_binding':{'identity_status':'PASS','bound_study_key':'human_rna:GSE1','producer_stage_sha256':'c'*64}}]}
        self.summary = {'reviewed_rows':1,'original_rows':1}

    def test_terminal_accepts_only_accounting_envelope(self):
        self.terminal_check(self.terminal,self.envelope)
        for field in ('input_rows','scientific_acceptance','accepted_scientific_repairs','original_paper_2025_eligibility','execution_status'):
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.terminal_check(self.terminal,{**self.envelope,field:'changed'})

    def test_missing_null_field_and_unknown_envelope_rejected(self):
        for remove, extra in [('scientific_error_rate',{}),(None,{'unbound_acceptance':True})]:
            broken={**self.envelope,**extra}
            if remove: broken.pop(remove)
            with self.assertRaises(ValueError): self.terminal_check(self.terminal,broken)

    def test_first_public_date_uses_source_quote_not_modification_time(self):
        self.date_check('2025-02-20','Public on Feb 20 2025','Status: Public on Feb 20 2025','2025-02-20')
        for args in [
            ('2025-02-21','Public on Feb 20 2025','Public on Feb 20 2025',None),
            ('2025-02-20','Public on Feb 20 2025','Modified on Feb 20 2025',None),
            ('2025-02-20','Public on Feb 20 2025','Public on Feb 20 2025','2024-02-20'),
            ('2025-02-20','Modified on Feb 20 2025','Modified on Feb 20 2025',None)]:
            with self.subTest(args=args), self.assertRaises(ValueError): self.date_check(*args)

    def test_exact_producer_identity_and_r_prediction_binding(self):
        self.producer_check(self.packet,self.summary,self.terminal,'GSE1','b'*64,'c'*64)
        for study,prediction,producer in [('GSE2','b'*64,'c'*64),('GSE1','d'*64,'c'*64),('GSE1','b'*64,'d'*64)]:
            with self.subTest(study=study,prediction=prediction,producer=producer),self.assertRaises(ValueError):
                self.producer_check(self.packet,self.summary,self.terminal,study,prediction,producer)

    def test_missing_or_ambiguous_identity_is_not_inferred(self):
        for binding in ({},{'identity_status':'CANNOT VERIFY','bound_study_key':'human_rna:GSE1','producer_stage_sha256':'c'*64},
                        {'identity_status':'PASS','bound_study_key':'human_rna:GSE2','producer_stage_sha256':'c'*64}):
            packet={**self.packet,'original_rows':[{'identity_binding':binding}]}
            with self.subTest(binding=binding),self.assertRaises(ValueError):
                self.producer_check(packet,self.summary,self.terminal,'GSE1','b'*64,'c'*64)

    def test_complete_review_does_not_mask_unreviewed_producer_rows(self):
        for summary in ({'reviewed_rows':0,'original_rows':1},{'reviewed_rows':1,'original_rows':2}):
            with self.subTest(summary=summary),self.assertRaises(ValueError):
                self.producer_check(self.packet,summary,self.terminal,'GSE1','b'*64,'c'*64)


class FuturePredictionMembershipTests(unittest.TestCase):
    """Disclosed file-contract fixtures, not empirical clock predictions."""
    def setUp(self):
        from audit_system.public_trace_2026_09_18 import _future_prediction_source
        import tempfile
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.pipeline = Path(self.temporary.name)/'pipeline'
        self.path = self.pipeline/'05_clock_predictions/attempt-000/predictions.csv'
        self.path.parent.mkdir(parents=True)
        self.path.write_text('sample_id,fixture_clock_a,fixture_clock_b,fixture_clock_c\nGSM1,1,2,3\n')
        self.stage = {'result':{'path':str(self.path),'clocks':['fixture_clock_a','fixture_clock_b','fixture_clock_c']},
                      'outputs':{'attempt-000/predictions.csv':sha(self.path)}}
        self.bind = _future_prediction_source
    def checked(self,path,expected):
        if sha(path) != expected: raise PipelineIntegrityError('Software fixture changed')
        return path
    def test_dnam_multiple_clock_file_without_rna_result_hash_is_bound(self):
        self.assertEqual(self.bind(self.pipeline,self.stage,self.checked),self.path)
        self.assertEqual(len(self.stage['result']['clocks']),3)
    def test_optional_result_hash_must_agree_with_member(self):
        self.stage['result']['prediction_sha256']=sha(self.path)
        self.bind(self.pipeline,self.stage,self.checked)
        self.stage['result']['prediction_sha256']='0'*64
        with self.assertRaises(ValueError):self.bind(self.pipeline,self.stage,self.checked)
    def test_missing_member_or_changed_bytes_are_rejected(self):
        original=self.stage['outputs'];self.stage['outputs']={}
        with self.assertRaises(ValueError):self.bind(self.pipeline,self.stage,self.checked)
        self.stage['outputs']=original;self.path.write_text('modified software fixture')
        with self.assertRaises(PipelineIntegrityError):self.bind(self.pipeline,self.stage,self.checked)
    def test_external_path_and_duplicate_alias_are_rejected(self):
        self.stage['result']['path']=str(self.pipeline/'elsewhere.csv')
        with self.assertRaises((ValueError,PipelineIntegrityError)):self.bind(self.pipeline,self.stage,self.checked)
        self.stage['result']['path']=str(self.path)
        self.stage['outputs']['attempt-000/../attempt-000/predictions.csv']=sha(self.path)
        with self.assertRaises(ValueError):self.bind(self.pipeline,self.stage,self.checked)


class MeasurementProjectionTests(unittest.TestCase):
    """Disclosed source-contract fixtures, not molecular or scientific results."""
    def setUp(self):
        self.fixture=CorpusProjectionTests('test_exact_denominators_scopes_identity_and_sanitization')
        self.fixture.setUp();self.addCleanup(self.fixture.doCleanups)
        f=self.fixture;self.root=f.root;self.study=f.study;self.pipeline=self.root/'pipeline'
        from audit_system.public_trace_2026_09_18 import _corpus_record_hash
        self.hash=_corpus_record_hash
        self.index={'study_id':'GSE1','sample_ids':['GSM1'],'evidence':{'S00000':{'artifact':'GEO','value':'Software fixture'},
            'L00000':{'artifact':'availability','field':'evidence_scope','value':'Original source-only software fixture','sample_ids':[]}}}
        f.stage('01_source_index',self.index)
        self.design=json.loads((self.study/'02_blinded_design/complete.json').read_text())['result']
        receipts={}
        for name in ('02_molecular_input','04_preprocessing','05_clock_predictions'):
            directory=self.pipeline/name;f.write(directory/'fixture.json',{'disclosed_software_fixture':True})
            f.write(directory/'complete.json',{'result':{},'outputs':{'fixture.json':sha(directory/'fixture.json')}})
            (directory/'complete.sha256').write_text(sha(directory/'complete.json'))
            receipts[name]={'path':str(directory/'complete.json'),'sha256':sha(directory/'complete.json')}
        self.measurement={'schema':'clockbase-pipeline-measurement-evidence-v1','study_id':'GSE1','modality':'human_rna',
            'evidence':{'M00000':{'kind':'fixture','value':{'source_matrix_present':True},'sample_ids':['GSM1'],
                'source':{'role':'software_fixture','sha256':sha(f.csv)}}},
            'source_files':{'software_fixture':{'path':str(f.csv),'sha256':sha(f.csv)}},'stage_receipts':receipts,
            'scientific_status':'CANNOT VERIFY','limitations':['Software fixture only']}
        f.task['pipeline_measurement_context']={'root':str(self.pipeline)}
        cfg=json.loads((self.study/'config.json').read_text());cfg['task']=f.task;f.write(self.study/'config.json',cfg)
        f.stage('02c_measurement_evidence',self.measurement)
        f.terminal['measurement_evidence_sha256']=self.hash(self.measurement);f.terminal['measurement_evidence_records']=1
        f.aggregate['scoped_scientific_checks'][0]['checks']['measurement_applicability']['evidence_ids']=['M00000','L00001']
        f.stage('05_aggregation',f.aggregate)
        self.scope=f.aggregate['scoped_scientific_checks'][0]
        f.stage('04_review_0000',{'checks':self.scope['checks'],'rows':f.aggregate['rows']})
        contract_path=self.study/'04_review_0000/contract-00.json'
        f.write(contract_path,{'status':'PASS','attempt':0})
        complete=json.loads((self.study/'04_review_0000/complete.json').read_text());complete['outputs']['contract-00.json']=sha(contract_path)
        f.write(self.study/'04_review_0000/complete.json',complete);(self.study/'04_review_0000/complete.sha256').write_text(sha(self.study/'04_review_0000/complete.json'))
        self.later_scope={'artifact':'availability','field':'post_design_measurement_scope','value':'Exact fixture scope, not reconstructed by the exporter.','sample_ids':[]}
        self.payload={'study_id':'GSE1','source':{**self.index,'evidence':{**self.index['evidence'],**self.measurement['evidence'],
            'L00001':self.later_scope,f.row_id:{'artifact':'saved_analysis'}}},'frozen_source_design':self.design,
            'measurement_evidence_sha256':self.hash(self.measurement)}
        self.bind_request()
    def bind_request(self):
        f=self.fixture;self.request_path=self.study/'api/review_0000_00-000/request.json'
        request={'messages':[{'role':'system','content':'Disclosed software fixture'},{'role':'user','content':json.dumps(self.payload)}]}
        f.write(self.request_path,request)
        binding={'schema':'clockbase-aggregate-call-binding-v1','call_slot':'review_0000_00-000','request_file_sha256':sha(self.request_path),
            'aggregate_key':self.hash({'study':str(self.study.resolve()),'actor':'review_0000_00','system':'Disclosed software fixture','payload':self.payload})}
        f.write(self.study/'api/review_0000_00-000.aggregate-binding.json',binding)
        f.seal();f.write(f.report_path,f.report)
    def project(self):
        from audit_system.public_trace_2026_09_18 import _corpus_projection_helpers,_project_corpus_study
        return _project_corpus_study(self.study,self.fixture.terminal,self.fixture.task,self.root,_corpus_projection_helpers({}))[0]
    def test_bound_measurement_records_and_exact_later_scope_are_public(self):
        p=self.project()
        self.assertEqual(p['source_evidence']['L00001'],self.later_scope)
        self.assertEqual(p['source_evidence']['L00000'],self.index['evidence']['L00000'])
        self.assertEqual(p['measurement_evidence']['evidence'],self.measurement['evidence'])
        self.assertEqual(p['measurement_evidence']['review_request_bindings'][0]['request_sha256'],sha(self.request_path))
        self.assertEqual(p['measurement_evidence']['registry_sha256'],self.hash(self.measurement))
        self.assertNotIn(str(self.root),json.dumps(p));self.assertNotIn('source_files',p['measurement_evidence'])
        self.assertEqual(p['scientific_status'],'CANNOT VERIFY')
    def test_changed_measurement_source_fails_even_when_public_values_look_plausible(self):
        self.fixture.csv.write_text('gse_id,clock,p_value\nGSE1,clock,0.9\n')
        with self.assertRaises(PipelineIntegrityError):self.project()
    def test_changed_registry_hash_is_not_accepted(self):
        self.fixture.terminal['measurement_evidence_sha256']='a'*64;self.bind_request()
        with self.assertRaises(ValueError):self.project()
    def test_request_with_different_m_values_is_rejected_after_resealing(self):
        self.payload['source']['evidence']['M00000']={'kind':'fixture','value':{'source_matrix_present':False}}
        self.bind_request()
        with self.assertRaises(ValueError):self.project()
    def test_unbound_api_request_is_rejected(self):
        self.request_path.write_text(self.request_path.read_text()+' ')
        self.fixture.seal()
        with self.assertRaises(ValueError):self.project()
    def test_missing_exact_post_design_scope_is_rejected(self):
        self.payload['source']['evidence'].pop('L00001');self.bind_request()
        with self.assertRaises(ValueError):self.project()
    def test_changed_frozen_design_in_review_payload_is_rejected(self):
        self.payload['frozen_source_design']={**self.design,'comparisons':[{'id':'invented'}]};self.bind_request()
        with self.assertRaises(ValueError):self.project()
    def test_absent_measurement_stage_cannot_satisfy_claimed_registry(self):
        import shutil
        shutil.rmtree(self.study/'02c_measurement_evidence');self.fixture.seal()
        with self.assertRaises(ValueError):self.project()

    def test_recovered_matrix_requires_bound_metadata_in_addition_to_numeric_stages(self):
        self.measurement['evidence']['M00000']['kind']='source_linked_dnam_input_recovery'
        def refresh_binding():
            self.fixture.stage('02c_measurement_evidence',self.measurement)
            self.fixture.terminal['measurement_evidence_sha256']=self.hash(self.measurement)
            self.payload['measurement_evidence_sha256']=self.hash(self.measurement)
            self.bind_request()
        refresh_binding()
        with self.assertRaisesRegex(ValueError,'Measurement registry'):self.project()
        directory=self.pipeline/'01_metadata';self.fixture.write(directory/'complete.json',{'result':{},'outputs':{}})
        (directory/'complete.sha256').write_text(sha(directory/'complete.json'))
        self.measurement['stage_receipts']['01_metadata']={'path':str(directory/'complete.json'),'sha256':sha(directory/'complete.json')}
        refresh_binding()
        self.assertIn('01_metadata',self.project()['measurement_evidence']['stage_sha256_by_name'])
        (directory/'complete.json').write_text('{}')
        with self.assertRaises(PipelineIntegrityError):self.project()


class GeoSourceBodyTests(unittest.TestCase):
    def test_download_format_keeps_gzip_tsv_and_plain_csv_distinct(self):
        import gzip
        from audit_system.public_trace_2026_09_18 import _tabular_format
        with tempfile.TemporaryDirectory() as folder:
            for name,payload,expected in (
                ('matrix.tsv.gz',gzip.compress(b'probe\tGSM1\n'),(True,'.tsv')),
                ('predictions.csv',b'sample_id,clock\n',(False,'.csv')),
                ('compressed.csv',gzip.compress(b'sample_id,clock\n'),(True,'.csv'))):
                path=Path(folder)/name;path.write_bytes(payload)
                self.assertEqual(_tabular_format(path),expected)

    def test_recovered_gzip_tsv_preserves_sample_and_value_checks(self):
        import gzip
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'recovered-beta.tsv.gz'
            path.write_bytes(gzip.compress(b'probe\tGSM1\tGSM2\ncg1\t0.1\t0.5\ncg2\tNA\t0.2\n'))
            features,samples,check=_matrix_qc(path,{'GSM1','GSM2'},'human_dnam')
            self.assertEqual((features,samples,check['status']),(2,2,'CANNOT VERIFY'))
            self.assertIn('4 cells inspected; 1 missing',check['reason'])
            with self.assertRaises(ValueError):_matrix_qc(path,{'GSM1','GSM3'},'human_dnam')

    def test_first_public_quote_uses_gzip_magic_without_filename_guessing(self):
        import gzip
        from audit_system.public_trace_2026_09_18 import _geo_source_text,_future_public_date
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'response.bin';text='!Series_status = Public on Jun 03 2025\n'
            for payload in (text.encode(),gzip.compress(text.encode())):
                path.write_bytes(payload)
                self.assertEqual(_geo_source_text(path),text)
                _future_public_date('2025-06-03','Public on Jun 03 2025',_geo_source_text(path))
                with self.assertRaises(ValueError):
                    _future_public_date('2025-06-04','Public on Jun 03 2025',_geo_source_text(path))


class AuditOnlyFollowupExportTests(unittest.TestCase):
    """Constructed filesystem contracts; projected science is a disclosed stub.

    MeasurementProjectionTests separately checks the unmocked shared projector.
    These cases isolate replay, immutable lineage, accounting and output safety.
    No model or molecular/statistical computation is performed by this fixture.
    """
    def setUp(self):
        import hashlib
        from audit_system.public_trace_2026_09_18 import _digest,_corpus_record_hash
        self.hash=_corpus_record_hash;self.digest=_digest
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        identity={'disclosed_software_fixture':True};vid=_digest(identity)
        self.version=self.root/vid;self.pipeline=self.version/'pipeline';self.prior=self.version/'comprehensive-audit/study'
        self.follow=self.root/'followup';self.study=self.follow/'audit/study';self.parent=self.root/'prior-public/GSE1'
        self.output=self.root/'next-public/GSE1'
        self.version.mkdir();self.follow.mkdir();self.parent.mkdir(parents=True)
        self.source=self.root/'fixture-source.csv';self.source.write_text('label\nDisclosed software fixture\n')
        sources={'matrix':{'path':str(self.source),'sha256':sha(self.source)}}
        version={'version_id':vid,'identity':identity,'study_id':'GSE1','modality':'human_rna','sources':sources,'implementation':{}}
        self.sealed(self.version/'version.json',version)
        self.write(self.pipeline/'config.json',{'study':'GSE1','software_fixture':True})
        prediction=self.pipeline/'05_clock_predictions/software-fixture-predictions.csv'
        prediction.parent.mkdir(parents=True);prediction.write_text('id,fixture\nGSM1,0\n')
        self.stage(self.pipeline,'05_clock_predictions',{'path':str(prediction),'prediction_sha256':sha(prediction)})
        self.stage(self.pipeline,'06_coding_agent',{'software_fixture':True})
        self.index={'study_id':'GSE1','sample_ids':['GSM1'],'evidence':{}}
        self.design={'source_index_sha256':self.hash(self.index),'frozen_before_saved_outcomes':True,'comparisons':[]}
        for directory in (self.prior,self.study):
            self.stage(directory,'01_source_index',self.index);self.stage(directory,'02_blinded_design',self.design)
        self.write(self.prior/'terminal.json',{'status':'REVIEWED','output_sha256':{}})
        (self.prior/'terminal.sha256').write_text(sha(self.prior/'terminal.json'))
        replay={'root':str(self.prior),'files':{str(p):{'path':str(p),'sha256':sha(p)} for p in self.prior.rglob('*') if p.is_file()}}
        allowance={'model':'gpt-5.6-luna','max_cost_usd':1.,'max_calls':2,'input_usd_per_million':0.5,'output_usd_per_million':1.8}
        request={'study_id':'GSE1','modality':'human_rna','version_id':vid,'source_bindings':sources,'pipeline_root':str(self.pipeline),
                 'output':str(self.follow/'audit'),'budget':allowance}
        self.frozen={'schema':'clockbase-pipeline-audit-followup-v1','revision':'software-fixture','request':request,
            'source_design_replay':replay,'input_files':{str(self.prior/'terminal.json'):{'path':str(self.prior/'terminal.json'),'sha256':sha(self.prior/'terminal.json')}},
            'implementation_sha256':{},'new_producer_or_clock_execution':False,'created_at_utc':'2026-09-21T00:00:00Z'}
        self.sealed(self.follow/'followup.json',self.frozen)
        task={'study_id':'GSE1','study_key':'human_rna:GSE1','modality':'human_rna','producer_version_id':vid,
            'source_design_replay':replay,'audit_followup':{'revision':'software-fixture','manifest':{'path':str(self.follow/'followup.json'),'sha256':sha(self.follow/'followup.json')}}}
        self.write(self.study/'config.json',{'task':task})
        self.stage(self.study,'02_blinded_design',self.design,{'replay-provenance.json':{'binding':replay,'new_model_calls':0,'new_independent_design_assessment':False}})
        self.write(self.study/'api/review_0000_00-000/request.json',{'software_fixture':True})
        self.write(self.study/'api/review_0000_00-000/receipt.json',{'requested_model':'gpt-5.6-luna','returned_model':'gpt-5.6-luna',
            'usage':{'prompt_tokens':10,'completion_tokens':5},'estimated_cost_usd':0.000014,'billed_cost_usd':None,'elapsed_seconds':0.0})
        self.terminal={'study_id':'GSE1','study_key':'human_rna:GSE1','modality':'human_rna','status':'REVIEWED','input_rows':1,
            'scientific_acceptance':'CANNOT VERIFY','accepted_scientific_repairs':0,'finished_at_utc':'2026-09-21T00:01:00Z','output_sha256':{}}
        self.result={**self.terminal,'execution_status':'COMPLETED','budget':{'call_reservations':1,'estimated_or_reserved_usd':0.000014,'billed_cost_usd':None},
            'pipeline_root_sha256':sha(self.pipeline/'config.json'),'original_paper_2025_eligibility':'QUARANTINED'}
        self.reseal()
        rows=[{'row_id':'R'+'a'*24,'fields':{'fixture':'software only'},'identity_binding':{'identity_status':'PASS','bound_study_key':'human_rna:GSE1',
            'producer_stage_sha256':sha(self.pipeline/'06_coding_agent/complete.json')}}]
        self.packet={'study_id':'GSE1','study_key':'human_rna:GSE1','original_rows':rows,'frozen_source_design':self.design,
            'independent_numerics':{'predictions_sha256':sha(prediction)},'measurement_evidence':{'registry_sha256':'a'*64,'evidence':{'M00000':{'fixture':True}},'review_request_bindings':[{'fixture':True}]}}
        self.summary={'study_id':'GSE1','study_key':'human_rna:GSE1','original_rows':1,'reviewed_rows':1,'artifact_id':'human_rna-GSE1','terminal_sha256':sha(self.study/'terminal.json')}
        self.write(self.parent/'public/comprehensive-audit.json',{'original_rows':rows})
        self.trace={'schema':'software-fixture','study_id':'GSE1','future_ingestion':{'version_id':vid,'source_version_sha256':sha(self.version/'version.json')},
            'comprehensive_audit':{'terminal_sha256':sha(self.prior/'terminal.json')},'future_independent_review':{'prior_fixture':True},
            'artifacts':[{'id':'comprehensive-audit','label':'Software fixture','sha256':sha(self.parent/'public/comprehensive-audit.json'),'download_path':'public/comprehensive-audit.json'}]}
        self.parent_receipt()
    def write(self,path,value):
        path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,sort_keys=True))
    def sealed(self,path,value):
        self.write(path,value);path.with_suffix(path.suffix+'.sha256').write_text(sha(path))
    def stage(self,root,name,result,extras=None):
        directory=root/name;self.write(directory/'result.json',result)
        outputs={'result.json':sha(directory/'result.json')}
        if name=='05_clock_predictions' and 'path' in result:
            prediction=Path(result['path']);outputs[str(prediction.relative_to(directory))]=sha(prediction)
        for rel,value in (extras or {}).items():self.write(directory/rel,value);outputs[rel]=sha(directory/rel)
        self.write(directory/'complete.json',{'result':result,'outputs':outputs});(directory/'complete.sha256').write_text(sha(directory/'complete.json'))
    def reseal(self):
        self.terminal['output_sha256']={str(p.relative_to(self.study)):sha(p) for p in self.study.rglob('*') if p.is_file() and p.name not in ('terminal.json','terminal.sha256')}
        self.write(self.study/'terminal.json',self.terminal);(self.study/'terminal.sha256').write_text(sha(self.study/'terminal.json'))
        self.result.update(self.terminal);self.sealed(self.follow/'result.json',self.result)
    def parent_receipt(self):
        self.write(self.parent/'audit_trace.json',self.trace)
        self.write(self.parent/'export-receipt.json',{'output_sha256':{str(p.relative_to(self.parent)):sha(p) for p in self.parent.rglob('*') if p.is_file() and p.name!='export-receipt.json'}})
    def export(self,**kwargs):
        from unittest.mock import patch
        from audit_system.public_trace_2026_09_18 import export_future_followup
        with patch('audit_system.public_trace_2026_09_18._project_corpus_study',return_value=(self.packet,self.summary)):
            return export_future_followup(self.follow,self.parent,self.output,**kwargs)
    def test_followup_keeps_original_files_replay_and_cost_scopes_separate(self):
        before={str(p.relative_to(self.parent)):p.read_bytes() for p in self.parent.rglob('*') if p.is_file()}
        result=self.export();self.assertEqual(result['status'],'PASS');self.assertEqual(result['new_model_calls'],0)
        trace=json.loads((self.output/'audit_trace.json').read_text());follow=trace['audit_followup']
        self.assertTrue(follow['grouping_design_reused']);self.assertFalse(follow['clocks_rerun']);self.assertFalse(follow['producer_rerun'])
        self.assertEqual(follow['additional_audit_budget']['call_reservations'],1);self.assertIsNone(follow['independent_review'])
        self.assertNotIn('future_independent_review',trace)
        self.assertEqual(json.loads((self.output/'public/superseded-audit-trace.json').read_text()),self.trace)
        self.assertEqual(before,{str(p.relative_to(self.parent)):p.read_bytes() for p in self.parent.rglob('*') if p.is_file()})
        receipt=json.loads((self.output/'export-receipt.json').read_text())
        for relative,digest in receipt['output_sha256'].items():self.assertEqual(sha(self.output/relative),digest)
        with self.assertRaises(ValueError):self.export()
    def test_changed_parent_bundle_is_refused(self):
        (self.parent/'public/comprehensive-audit.json').write_text('{}')
        with self.assertRaises(PipelineIntegrityError):self.export()
    def test_followup_time_is_canonical_utc_without_changing_the_instant(self):
        from datetime import datetime
        self.terminal['finished_at_utc']='2026-09-21T00:01:00.538731+00:00';self.reseal()
        self.export();recorded=json.loads((self.output/'audit_trace.json').read_text())['recorded_at_utc']
        self.assertEqual(recorded,'2026-09-21T00:01:00.538731Z')
        self.assertEqual(datetime.fromisoformat(recorded),datetime.fromisoformat(self.result['finished_at_utc']))
    def test_followup_time_without_explicit_utc_is_refused(self):
        self.terminal['finished_at_utc']='2026-09-21T00:01:00';self.reseal()
        with self.assertRaisesRegex(ValueError,'identify UTC explicitly'):self.export()
        self.assertFalse(self.output.exists())
    def test_changed_design_is_refused_even_after_new_stage_sealing(self):
        self.stage(self.study,'02_blinded_design',{**self.design,'comparisons':[{'id':'changed'}]})
        with self.assertRaises(ValueError):self.export()
    def test_changed_original_rows_cannot_be_promoted(self):
        self.packet['original_rows'][0]['fields']={'fixture':'changed'}
        with self.assertRaises(ValueError):self.export()
    def test_inconsistent_additional_cost_is_refused(self):
        self.result['budget']['estimated_or_reserved_usd']=0.5;self.reseal()
        with self.assertRaises(ValueError):self.export()
    def test_independent_readback_requires_all_three_execution_hashes(self):
        review={'status':'PASS','scientific_status':'CANNOT VERIFY','scope':'Disclosed software fixture',
            'followup_manifest_sha256':sha(self.follow/'followup.json'),'result_sha256':sha(self.follow/'result.json'),
            'terminal_sha256':'f'*64}
        path=self.root/'independent.json';self.write(path,review)
        with self.assertRaises(ValueError):self.export(independent_review_path=path)
        review['terminal_sha256']=sha(self.study/'terminal.json');self.write(path,review)
        self.export(independent_review_path=path)
        trace=json.loads((self.output/'audit_trace.json').read_text());self.assertEqual(trace['audit_followup']['independent_review']['source_receipt_sha256'],sha(path))
    def test_missing_request_receipt_cannot_appear_as_zero_additional_cost(self):
        (self.study/'api/review_0000_00-000/receipt.json').unlink();self.result['budget']={'call_reservations':0,'estimated_or_reserved_usd':0.0};self.reseal()
        with self.assertRaises(ValueError):self.export()
    def test_usage_tokens_must_reproduce_frozen_tariff_cost(self):
        path=self.study/'api/review_0000_00-000/receipt.json';receipt=json.loads(path.read_text())
        receipt['usage']['prompt_tokens']=100;self.write(path,receipt);self.reseal()
        with self.assertRaises(ValueError):self.export()


class FreshSourceReassessmentExportTests(unittest.TestCase):
    """Software-only lineage fixtures; science and producer export are stubbed."""
    write = AuditOnlyFollowupExportTests.write
    sealed = AuditOnlyFollowupExportTests.sealed
    stage = AuditOnlyFollowupExportTests.stage
    reseal = AuditOnlyFollowupExportTests.reseal
    parent_receipt = AuditOnlyFollowupExportTests.parent_receipt

    def setUp(self):
        import shutil
        AuditOnlyFollowupExportTests.setUp(self)
        shutil.rmtree(self.prior/'02_blinded_design')
        prior = {'status':'UNRESOLVED','study_id':'GSE1','output_sha256':{}}
        self.write(self.prior/'terminal.json',prior)
        (self.prior/'terminal.sha256').write_text(sha(self.prior/'terminal.json'))
        self.sealed(self.version/'comprehensive-audit-result.json',
            {**prior,'execution_status':'UNRESOLVED','budget':{'call_reservations':2}})
        series = self.root/'series.txt';series.write_text('!Series_status = Public on Jun 03 2025\n')
        version = json.loads((self.version/'version.json').read_text())
        version['sources']['series'] = {'path':str(series),'sha256':sha(series)}
        version['source_design'] = {'series_first_public_date':'2025-06-03','source_status_quote':'Public on Jun 03 2025'}
        self.sealed(self.version/'version.json',version)
        self.frozen.update(schema='clockbase-pipeline-audit-followup-v2',mode='fresh_source_v1')
        self.frozen.pop('source_design_replay')
        self.frozen['request']['source_bindings'] = version['sources']
        self.frozen['input_files'] = {str(p):{'path':str(p),'sha256':sha(p)} for p in
            (self.prior/'terminal.json',self.version/'comprehensive-audit-result.json')}
        self.source_request = {'model':'gpt-5.6-luna','messages':[
            {'role':'system','content':'Disclosed software-only source fixture'},
            {'role':'user','content':json.dumps({'source_index':self.index})}]}
        self.frozen['fresh_source'] = {'prior_audit':str(self.prior.parent),
            'source_index_sha256':self.hash(self.index),
            'initial_request_sha256':self.hash({'system':self.source_request['messages'][0]['content'],
                                              'payload':{'source_index':self.index}})}
        self.sealed(self.follow/'followup.json',self.frozen)
        task = json.loads((self.study/'config.json').read_text())['task']
        task.pop('source_design_replay')
        task['audit_followup'] = {'revision':self.frozen['revision'],'mode':'fresh_source_v1',
            'manifest':{'path':str(self.follow/'followup.json'),'sha256':sha(self.follow/'followup.json')}}
        self.write(self.study/'config.json',{'task':task})
        shutil.rmtree(self.study/'02_blinded_design')
        self.stage(self.study,'02_blinded_design',self.design)
        self.write(self.study/'api/source_design_00-000/request.json',self.source_request)
        self.write(self.study/'api/source_design_00-000/receipt.json',
            json.loads((self.study/'api/review_0000_00-000/receipt.json').read_text()))
        self.result['budget'].update(call_reservations=2,estimated_or_reserved_usd=0.000028,max_calls=2,max_cost_usd=1.)
        self.reseal()
        self.review_path=self.root/'independent-readback.json'
        self.review={'status':'PASS','scientific_status':'CANNOT VERIFY','scope':'Disclosed software-only readback fixture',
            'followup_manifest_sha256':sha(self.follow/'followup.json'),'result_sha256':sha(self.follow/'result.json'),
            'terminal_sha256':sha(self.study/'terminal.json'),
            'source_claim_findings':[{'claim':'Example label','assessment':'Unsupported by this software fixture',
                'source_evidence':[],'consequence':'Retain unresolved; do not accept a repair.'}]}
        self.write(self.review_path,self.review)

    def export(self, review=True):
        from unittest.mock import patch
        from audit_system.public_trace_2026_09_18 import export_audit_reassessment
        plan=self.root/'plan.json'
        def prepare_stub(run,path):
            self.write(path,{'software_fixture':True})
        def export_stub(path,output):
            output.mkdir(parents=True);(output/'public').mkdir()
            self.write(output/'audit_trace.json',{'study_id':'GSE1','models':[],'artifacts':[]})
            self.write(output/'export-receipt.json',{'software_fixture':True})
        with patch('audit_system.public_trace_2026_09_18._project_corpus_study',return_value=(self.packet,self.summary)), \
                patch('audit_system.public_trace_2026_09_18.prepare',side_effect=prepare_stub), \
                patch('audit_system.public_trace_2026_09_18.export',side_effect=export_stub):
            return export_audit_reassessment(self.follow,self.output,plan_path=plan,
                independent_review_path=self.review_path if review else None)

    def test_fresh_source_trace_retains_independent_findings_and_prior_failure(self):
        protected={str(p):p.read_bytes() for d in (self.follow,self.version) for p in d.rglob('*') if p.is_file()}
        result=self.export();self.assertEqual(result['status'],'PASS');self.assertEqual(result['new_model_calls'],0)
        trace=json.loads((self.output/'audit_trace.json').read_text());audit=trace['audit_reassessment']
        self.assertEqual(audit['schema'],'clockbase-public-audit-reassessment-v2')
        self.assertEqual(audit['prior_audit']['status'],'UNRESOLVED')
        self.assertEqual(audit['independent_review']['source_claim_findings'],self.review['source_claim_findings'])
        self.assertEqual(audit['independent_readback_sha256'],sha(self.review_path))
        self.assertNotIn('execution_receipt_sha256',audit)
        self.assertTrue(audit['new_independent_design_assessment']);self.assertFalse(audit['grouping_design_reused'])
        self.assertEqual(audit['additional_audit_budget']['call_reservations'],2)
        self.assertIn('independent-source-readback',{a['id'] for a in trace['artifacts']})
        self.assertEqual(protected,{p:Path(p).read_bytes() for p in protected})

    def test_fresh_source_export_requires_bound_independent_readback(self):
        with self.assertRaisesRegex(ValueError,'independent readback'):self.export(review=False)
        self.review['terminal_sha256']='0'*64;self.write(self.review_path,self.review)
        with self.assertRaisesRegex(ValueError,'does not bind'):self.export()
        self.assertFalse(self.output.exists())

    def test_first_source_request_cannot_gain_prior_outcomes(self):
        self.source_request['messages'][1]['content']=json.dumps({'prior_verdict':'PASS'})
        self.write(self.study/'api/source_design_00-000/request.json',self.source_request);self.reseal()
        with self.assertRaisesRegex(ValueError,'prelaunch blinded packet'):self.export()

    def test_replayed_design_cannot_be_described_as_fresh(self):
        path=self.study/'config.json';cfg=json.loads(path.read_text());cfg['task']['source_design_replay']={}
        self.write(path,cfg);self.reseal()
        with self.assertRaisesRegex(ValueError,'unchanged producer'):self.export()

    def test_changed_source_index_is_refused(self):
        index={**self.index,'sample_ids':['GSM2']}
        self.stage(self.study,'01_source_index',index)
        self.stage(self.study,'02_blinded_design',{**self.design,'source_index_sha256':self.hash(index)})
        self.reseal()
        with self.assertRaisesRegex(ValueError,'frozen source index'):self.export()

    def test_changed_original_failure_is_refused(self):
        self.write(self.prior/'terminal.json',{'status':'REVIEWED'})
        with self.assertRaises(PipelineIntegrityError):self.export()




class PartialFutureBindingTests(unittest.TestCase):
    """Constructed lifecycle fixtures; no empirical measurements or API calls."""
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.version = self.root/('a'*64)
        self.pipeline = self.version/'pipeline'
        self.pipeline.mkdir(parents=True)
        from audit_system.public_trace_2026_09_18 import _partial_future_snapshot
        self.validate = _partial_future_snapshot
        self.write(self.pipeline/'config.json', {'study': 'GSE999001'})
        self.write(self.version/'version.json', {'version_id': self.version.name, 'study_id': 'GSE999001',
            'sources': {}, 'source_sha256': {}, 'implementation': {}})
        self.seal(self.version/'version.json', '.json.sha256')
        for name in STAGES[:5]: self.stage(self.pipeline/name)
        self.write(self.version/'pipeline-report.json', {'status': 'FAILED', 'output_sha256': {}})
        self.seal(self.version/'pipeline-report.json', '.json.sha256')
        self.failure = self.pipeline/'06_coding_agent/attempt-000/failure.json'
        self.write(self.failure, {'status': 'CANNOT VERIFY', 'error': 'Disclosed software contract fixture', 'elapsed_seconds': 1.0})
        self.api(self.pipeline/'api/primary-000')

    def write(self, path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))

    def seal(self, path, suffix):
        path.with_suffix(suffix).write_text(sha(path))

    def stage(self, directory):
        self.write(directory/'attempt-000/fixture.json', {'disclosed_software_fixture': True})
        self.write(directory/'complete.json', {'result': {}, 'outputs': {'attempt-000/fixture.json': sha(directory/'attempt-000/fixture.json')}})
        self.seal(directory/'complete.json', '.sha256')

    def api(self, directory):
        usage = {'prompt_tokens': 1, 'completion_tokens': 1}
        self.write(directory/'request.json', {'model': 'software-fixture'})
        self.write(directory/'reservation.json', {'model': 'software-fixture'})
        self.write(directory/'response.json', {'model': 'software-fixture', 'usage': usage})
        self.write(directory/'receipt.json', {'requested_model': 'software-fixture', 'returned_model': 'software-fixture',
            'usage': usage, 'response_sha256': sha(directory/'response.json'), 'estimated_cost_usd': 0.0})

    def binding(self, complete=False):
        return {'schema': 'clockbase-partial-update-binding-v1', 'integrity_status': 'PASS',
            'study_id': 'GSE999001', 'version_id': self.version.name, 'version_sha256': sha(self.version/'version.json'),
            'artifact_sha256': {str(p.relative_to(self.version)): sha(p) for p in self.version.rglob('*') if p.is_file()},
            'completed_stage_count': 8 if complete else 5, 'pipeline_status': 'COMPLETED' if complete else 'FAILED',
            'comprehensive_audit_status': 'UNRESOLVED' if complete else 'MISSING',
            'api_call_count': 3 if complete else 1, 'estimated_cost_usd': 0.0,
            'scientific_status': 'CANNOT VERIFY', 'original_paper_2025_eligibility': 'QUARANTINED'}

    def unresolved(self):
        self.failure.unlink()
        for name in STAGES[5:]: self.stage(self.pipeline/name)
        self.write(self.version/'pipeline-report.json', {'status': 'COMPLETED', 'output_sha256': {}})
        self.seal(self.version/'pipeline-report.json', '.json.sha256')
        self.write(self.version/'pipeline-report-result.json', {'execution_status': 'PASS', 'accepted_results': 0, 'scientific_status': 'CANNOT VERIFY'})
        audit = self.version/'comprehensive-audit/study'
        self.write(audit/'config.json', {'task': {'producer_version_id': self.version.name}})
        self.stage(audit/'01_source_index')
        self.failure = audit/'02_blinded_design/attempt-000/failure.json'
        self.write(self.failure, {'status': 'CANNOT VERIFY', 'error': 'Disclosed software source contract fixture', 'elapsed_seconds': 1.0})
        for index in range(2):
            self.write(self.failure.parent/f'contract-{index:02d}.json', {'status': 'REJECTED_RESPONSE_CONTRACT',
                'attempt': index, 'error': 'Disclosed software fixture', 'scientific_status': 'CANNOT VERIFY'})
            self.api(audit/f'api/source_design_{index:02d}-000')
        terminal = {'study_id': 'GSE999001', 'status': 'UNRESOLVED', 'scientific_acceptance': 'CANNOT VERIFY',
            'accepted_scientific_repairs': 0, 'reviewed_rows': 0, 'calculations': 0, 'output_sha256': {}}
        result = {**terminal, 'execution_status': 'UNRESOLVED', 'budget': {},
            'pipeline_root_sha256': sha(self.pipeline/'config.json'), 'original_paper_2025_eligibility': 'QUARANTINED'}
        self.write(audit/'terminal.json', terminal); self.seal(audit/'terminal.json', '.sha256')
        self.write(self.version/'comprehensive-audit-result.json', result)
        self.write(self.version/'comprehensive-audit.json', {'status': 'UNRESOLVED', 'result': result, 'output_sha256': {}})
        self.seal(self.version/'comprehensive-audit.json', '.json.sha256')
        return audit

    def test_coding_stop_keeps_five_stages_and_full_exporter_rejects_it(self):
        value = self.validate(self.version, self.binding())
        self.assertEqual(list(value[3]), list(STAGES[:5]))
        self.assertEqual(value[4]['outcome'], 'coding_stop')
        self.assertEqual(value[4]['failure']['stage'], '06_coding_agent')
        self.assertEqual(len(value[4]['api_calls']), 1)
        with self.assertRaises(ValueError): _snapshot(self.pipeline)

    def test_changed_bytes_or_added_files_cannot_use_original_binding(self):
        bound = self.binding(); self.write(self.failure, {'status': 'PASS'})
        with self.assertRaises(PipelineIntegrityError): self.validate(self.version, bound)

    def test_coding_stop_rejects_resealed_later_outputs_or_missing_failure(self):
        for name in ('pipeline/07_audit_and_repair', 'comprehensive-audit'):
            path = self.version/name; path.mkdir()
            with self.assertRaises(ValueError): self.validate(self.version, self.binding())
            path.rmdir()
        self.failure.unlink()
        with self.assertRaises(ValueError): self.validate(self.version, self.binding())

    def test_partial_identity_cost_and_stage_count_are_not_inferred(self):
        for key, value in [('study_id', 'GSE2'), ('version_sha256', 'b'*64), ('completed_stage_count', 8),
                           ('api_call_count', 2), ('estimated_cost_usd', 1.0), ('scientific_status', 'PASS')]:
            bound = self.binding(); bound[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError): self.validate(self.version, bound)

    def test_resealed_stage_receipt_does_not_hide_a_changed_output(self):
        self.write(self.pipeline/STAGES[0]/'attempt-000/fixture.json', {'changed_fixture': True})
        with self.assertRaises(PipelineIntegrityError): self.validate(self.version, self.binding())

    def test_unresolved_comprehensive_keeps_two_rejected_contracts(self):
        self.unresolved(); value = self.validate(self.version, self.binding(True))
        self.assertEqual(len(value[3]), 8)
        self.assertEqual(value[4]['comprehensive_status'], 'UNRESOLVED')
        self.assertEqual(len(value[4]['contract_attempts']), 2)

    def test_unresolved_rejects_resealed_terminal_binding_and_later_stage(self):
        audit = self.unresolved()
        path = audit/'config.json'; original = json.loads(path.read_text())
        self.write(path, {'task': {'producer_version_id': 'b'*64}})
        with self.assertRaises(ValueError): self.validate(self.version, self.binding(True))
        self.write(path, original); self.stage(audit/'03_conditional_reruns')
        with self.assertRaises(ValueError): self.validate(self.version, self.binding(True))

    def test_rejected_contract_cannot_be_accepted_or_omitted(self):
        self.unresolved(); path = self.failure.parent/'contract-01.json'
        self.write(path, {'status': 'PASS', 'attempt': 1})
        with self.assertRaises(ValueError): self.validate(self.version, self.binding(True))
        path.unlink()
        with self.assertRaises(ValueError): self.validate(self.version, self.binding(True))



class PlannedClockCheckpointTests(unittest.TestCase):
    """Constructed software lifecycle only; no downloaded data, models or clocks."""
    def setUp(self):
        from audit_system.public_trace_2026_09_18 import _planned_clock_snapshot, export_planned_clock_checkpoint
        from audit_system.agent_workflow_2026_09_10 import source_packet, review_messages
        from src.pipeline_inputs_2026_09_10 import acquire_study_context, parse_samples
        from decimal import Decimal
        import hashlib
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.validate=_planned_clock_snapshot;self.export=export_planned_clock_checkpoint
        self.study='GSE999001';self.ids=['GSM999001','GSM999002']
        source=self.root/'source';source.mkdir()
        series=source/'series.soft';samples=source/'samples.soft'
        series.write_text('^SERIES = GSE999001\n!Series_status = Public on Jan 01 2026\n'
            '!Series_title = Disclosed software lifecycle fixture\n!Series_overall_design = Constructed contract test only\n'
            '!Series_sample_id = GSM999001\n!Series_sample_id = GSM999002\n')
        samples.write_text(''.join('^SAMPLE = '+sid+'\n!Sample_title = Disclosed software fixture '+str(i)+'\n'
            '!Sample_source_name_ch1 = Software fixture\n!Sample_organism_ch1 = Homo sapiens\n'
            '!Sample_molecule_ch1 = genomic DNA\n!Sample_platform_id = GPL13534\n' for i,sid in enumerate(self.ids)))
        self.options={'model':'gpt-6-luna','audit_model':'gpt-6-luna','max_calls':2,'max_cost_usd':.4,
            'max_output_tokens':100,'max_prompt_bytes':50000,'audit_input_usd_per_million':.25,
            'audit_output_usd_per_million':.75,'audit_validation_retries':getattr(self,'fixture_retries',1),'source_contract_policy':'literal_field_refs_v1',
            'audit_profile':'source_domains_v4','source_binding_profile':'literal_source_v1',
            'comparison_scope':'source_defined_interventions'}
        from src.pipeline_2026_09_10 import assess_design
        import inspect
        implementation=Path(inspect.getsourcefile(assess_design)).read_bytes()
        implementation_name='src/pipeline_2026_09_10.py'
        sources={role:{'path':str(p),'sha256':sha(p)} for role,p in [('series',series),('samples',samples)]}
        identity={'study_id':self.study,'modality':'human_dnam','target':'Disclosed software fixture',
            'sources':sources,'source_sha256':{k:v['sha256'] for k,v in sources.items()},
            'pipeline_options':self.options,'comprehensive_audit':{'model':'gpt-6-luna'},
            'implementation':{implementation_name:hashlib.sha256(implementation).hexdigest()},
            'source_design':{'sample_ids':self.ids,'series_first_public_date':'2026-01-01','source_status_quote':'Public on Jan 01 2026'},
            'explicit_analysis_revision':'software_checkpoint_fixture','supersedes_version':'a'*64}
        from src.update_pipeline_2026_09_20 import identity as producer_version_identity
        version_id=producer_version_identity(identity)
        registry=self.root/'registry';self.version=registry/'studies/human_dnam'/self.study/version_id
        self.pipeline=self.version/'pipeline';self.pipeline.mkdir(parents=True)
        self.write(self.version/'version.json',{'identity':identity,'version_id':version_id,**identity});self.seal(self.version/'version.json')
        implementation_path=self.version/'implementation'/implementation_name
        implementation_path.parent.mkdir(parents=True);implementation_path.write_bytes(implementation)
        self.write(self.pipeline/'config.json',{'study':self.study,'modality':'human_dnam','target':identity['target'],
            'implementation':identity['implementation'],**self.options})
        self.plan=self.root/'update-plan.json';self.scope=self.root/'planned-scope.json';self.freeze=self.root/'protocol-freeze.json'
        self.write(self.plan,{'schema':'clockbase-continuous-update-plan-v1','pipeline':self.options,
            'comprehensive_audit':identity['comprehensive_audit'],'limits':{'allow_network':False,'concurrency':1},
            'studies':[{'study_id':self.study,'modality':'human_dnam','target':identity['target'],'sources':sources,
                'supersedes_version':identity['supersedes_version'],'analysis_revision':identity['explicit_analysis_revision']}]})
        self.write(self.scope,{'allowed_until_sequence':['input','design','clocks'],'prohibited_until':['agent','audit','report'],
            'input_plan':{'sha256':sha(self.plan)},'registry':str(registry),'study_id':self.study,
            'analysis_revision':identity['explicit_analysis_revision'],'supersedes_version':identity['supersedes_version']})
        meta=self.pipeline/STAGES[0]/'attempt-000';meta.mkdir(parents=True)
        (meta/'samples.soft').write_bytes(samples.read_bytes());(meta/'series.soft').write_bytes(series.read_bytes())
        rows=parse_samples(samples.read_text(),self.study)
        with (meta/'metadata.csv').open('w',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=sorted(rows[0]));writer.writeheader();writer.writerows(rows)
        context=acquire_study_context(self.study,meta/'series.soft',meta,publication_methods=False)
        metadata={'metadata':str(meta/'metadata.csv'),'sample_ids':self.ids,'study_context':str(context),
            'source_sample_records':str(meta/'samples.soft'),'series_source':str(meta/'series.soft')}
        self.stage(STAGES[0],metadata)
        packet=source_packet(meta/'metadata.csv',self.study,identity['target'],context=context,
            sample_source=meta/'samples.soft',comparison_scope=self.options['comparison_scope'])
        messages=review_messages(packet,review_profile='source_domains_v4',source_binding_profile='literal_source_v1')
        transport=json.loads(messages[1].content)
        self.attempt=self.pipeline/STAGES[2]/'attempt-000'
        from types import SimpleNamespace
        from langchain_core.messages import AIMessage
        import copy
        self.attempt.mkdir(parents=True)
        final=LiteralFieldFailureTests._design('valid');final['unassigned_ids']=self.ids
        first=copy.deepcopy(final);first['covariate_plan'][0]['source_fields']=['title']
        responses=[first,final] if self.options['audit_validation_retries'] else [final];self.cost=0.;calls=[]
        self.first=self.root/'first-source-request.json'
        def invoke(messages):
            index=len(calls);calls.append(index);slot=self.pipeline/'api'/('design-%03d'%index)
            request={'model':'gpt-6-luna','max_completion_tokens':100,'reasoning_effort':'low','store':False,
                'response_format':{'type':'json_object'},'messages':[{'role':{'human':'user','ai':'assistant'}.get(m.type,m.type),
                    'content':m.content} for m in messages]}
            if index==0:self.write(self.first,request)
            content=json.dumps(responses[index],ensure_ascii=False)
            usage={'prompt_tokens':100+index,'completion_tokens':10}
            response={'model':'gpt-6-luna','usage':usage,'choices':[{'finish_reason':'stop','message':{'content':content}}]}
            self.write(slot/'request.json',request);self.write(slot/'response.json',response)
            raw=json.dumps(request,ensure_ascii=False).encode()
            reserved=float((Decimal(len(raw)+100*len(request['messages']))*Decimal('.25')+Decimal(100)*Decimal('.75'))/Decimal(1000000))
            cost=float((Decimal(usage['prompt_tokens'])*Decimal('.25')+Decimal(10)*Decimal('.75'))/Decimal(1000000))
            self.cost+=cost
            self.write(slot/'reservation.json',{'model':'gpt-6-luna','automatic_retries':0,'concurrency':1,'reserved_usd':reserved,'started_unix':1767312000.+index})
            self.write(slot/'receipt.json',{'requested_model':'gpt-6-luna','returned_model':'gpt-6-luna',
                'usage':usage,'response_sha256':sha(slot/'response.json'),'estimated_cost_usd':cost,'billed_cost_usd':None,'elapsed_seconds':1.})
            return AIMessage(content=content)
        self.design=assess_design(packet,SimpleNamespace(invoke=invoke),self.attempt,
            review_profile='source_domains_v4',source_binding_profile='literal_source_v1',
            source_contract_policy='literal_field_refs_v1',validation_retries=self.options['audit_validation_retries'])
        self.call_count=len(calls)
        self.stage(STAGES[2],self.design)
        for name in (STAGES[1],STAGES[3]):
            self.write(self.pipeline/name/'attempt-000/fixture.json',{'disclosed_software_fixture':True});self.stage(name,{})
        self.predictions=self.pipeline/STAGES[4]/'attempt-000/clock_predictions.csv';self.predictions.parent.mkdir(parents=True)
        self.predictions.write_text('sample_id,SoftwareFixtureClock\nGSM999001,1\nGSM999002,2\n')
        self.clocks={'path':str(self.predictions),'clocks':['SoftwareFixtureClock'],'sample_count':2,'clock_predictions_recomputed':True}
        self.stage(STAGES[4],self.clocks)
        self.write(self.version/'pipeline-clocks-request.json',{'phase':'pipeline','version_path':str(self.version/'version.json'),
            'pipeline_root':str(self.pipeline),'cache':str(registry/'cache'),'until':'clocks','result_path':str(self.version/'pipeline-clocks-result.json')})
        self.write(self.version/'pipeline-clocks-result.json',self.clocks);self.wrapper()
        self.write(self.freeze,{'schema':'clockbase-planned-clock-protocol-freeze-v1','status':'FROZEN','study_id':self.study,
            'version_id':version_id,'utc':'2026-01-01T00:00:00Z','plan_sha256':sha(self.plan),'first_source_request_sha256':sha(self.first),
            'allowed_until_sequence':['input','design','clocks'],'prohibited_until':['agent','audit','report'],
            'provider':'openai','model':'gpt-6-luna','serial_concurrency':1,'max_new_reserved_calls':2,'max_new_reserved_usd':'.4',
            'scientific_status':'CANNOT VERIFY','accepted_scientific_repairs':0,
            'model_calls_at_freeze':0,'new_reservations_at_freeze':0,'clock_runs_at_freeze':0})

    def write(self,path,value):
        path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

    def seal(self,path):
        path.with_name(path.name+'.sha256').write_text(sha(path)+'\n')

    def stage(self,name,result=None):
        folder=self.pipeline/name
        if result is None:result=json.loads((folder/'complete.json').read_text())['result']
        outputs={str(p.relative_to(folder)):sha(p) for p in folder.rglob('*') if p.is_file() and p.name not in {'complete.json','complete.sha256'}}
        self.write(folder/'complete.json',{'result':result,'outputs':outputs})
        (folder/'complete.sha256').write_text(sha(folder/'complete.json')+'\n')

    def wrapper(self):
        outputs={str(p):sha(p) for p in self.pipeline.rglob('*') if p.is_file() and p.name!='.run.lock'}
        outputs[str(self.version/'pipeline-clocks-result.json')]=sha(self.version/'pipeline-clocks-result.json')
        self.write(self.version/'pipeline-clocks.json',{'status':'COMPLETED','until':'clocks','output_sha256':outputs})
        self.seal(self.version/'pipeline-clocks.json')

    def binding(self):
        return {'schema':'clockbase-planned-clock-checkpoint-binding-v1','integrity_status':'PASS','scientific_status':'CANNOT VERIFY',
            'accepted_scientific_repairs':0,'original_paper_2025_eligibility':'QUARANTINED','completed_stage_count':5,
            'pipeline_status':'COMPLETED','comprehensive_audit_status':'NOT_RUN','study_id':self.study,
            'version_id':self.version.name,'version_sha256':sha(self.version/'version.json'),'api_call_count':self.call_count,'estimated_cost_usd':self.cost,
            **{name:{'path':str(path),'sha256':sha(path)} for name,path in [('update_plan',self.plan),('planned_scope',self.scope),
                ('first_source_request',self.first),('protocol_freeze',self.freeze)]},
            'artifact_sha256':{str(p.relative_to(self.version)):sha(p) for p in sorted(self.version.rglob('*')) if p.is_file()}}

    def test_planned_export_keeps_contract_history_without_analysis_failure(self):
        value=self.validate(self.version,self.binding())[4]
        self.assertEqual(value['scientific_status'],'CANNOT VERIFY');self.assertEqual(value['accepted_scientific_repairs'],0)
        self.assertEqual([x['recorded_validation_status'] for x in value['source_response_history']],['FAIL','PASS'])
        self.assertEqual([x['retry_transmitted'] for x in value['source_response_history']],[True,False])
        binding=self.root/'binding.json';self.write(binding,self.binding())
        self.export(self.version,self.root/'export',readback_path=binding,plan_path=self.root/'export-plan.json')
        trace=json.loads((self.root/'export/audit_trace.json').read_text())
        self.assertNotIn('partial_update',trace);self.assertEqual(trace['comparisons'],[])
        self.assertEqual(trace['resource_usage']['responses_received'],2)
        self.assertEqual((self.root/'export/public/sample-predictions.csv').read_bytes(),self.predictions.read_bytes())

    def test_self_consistent_spaced_identity_is_not_a_producer_version(self):
        import hashlib
        from unittest.mock import patch
        fixture = PlannedClockCheckpointTests()
        try:
            # Construct the old, internally consistent fixture without changing
            # the exporter. A producer-derived fixture is exercised above.
            with patch('src.update_pipeline_2026_09_20.identity', side_effect=lambda value:
                       hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()):
                fixture.setUp()
            with self.assertRaisesRegex(ValueError, '^Planned checkpoint version/configuration identity differs$'):
                fixture.validate(fixture.version, fixture.binding())
        finally:
            fixture.doCleanups()

    def test_scope_or_freeze_cannot_be_retroactive_or_allow_later_stages(self):
        for path,key,value in [(self.scope,'prohibited_until',[]),(self.scope,'supersedes_version','b'*64),
            (self.freeze,'utc','2027-01-01T00:00:00Z'),(self.freeze,'first_source_request_sha256','b'*64),
            (self.freeze,'model_calls_at_freeze',1),(self.freeze,'max_new_reserved_calls',3)]:
            original=path.read_bytes();payload=json.loads(original);payload[key]=value;self.write(path,payload)
            with self.subTest(key=key),self.assertRaises(ValueError):self.validate(self.version,self.binding())
            path.write_bytes(original)

    def test_failure_or_later_output_cannot_become_planned_stop(self):
        for relative in ['pipeline/06_coding_agent/evidence.json','pipeline/03_blinded_design/attempt-000/failure.json',
                         'pipeline-report.json','comprehensive-audit/terminal.json','pipeline/api/primary-000/request.json']:
            path=self.version/relative;self.write(path,{'disclosed_software_fixture':True});self.wrapper()
            with self.subTest(path=relative),self.assertRaises(ValueError):self.validate(self.version,self.binding())
            path.unlink()
            while path.parent not in (self.pipeline,self.version,self.attempt) and not list(path.parent.iterdir()):
                path=path.parent;path.rmdir()
            self.wrapper()

    def test_resealed_source_cost_membership_and_history_changes_rejected(self):
        cases=[('pipeline/api/design-001/receipt.json','estimated_cost_usd',1.),
            ('pipeline/api/design-001/receipt.json','returned_model','different-model'),
            ('pipeline/03_blinded_design/attempt-000/response-validation/000-validation.json','validation_status','PASS'),
            ('pipeline/03_blinded_design/attempt-000/response-validation/001-validation.json','validation_status','FAIL'),
            ('pipeline/03_blinded_design/attempt-000/source-packet.json','target','Changed source packet')]
        for relative,key,value in cases:
            path=self.version/relative;original=path.read_bytes();payload=json.loads(original);payload[key]=value;self.write(path,payload)
            self.stage(STAGES[2]);self.wrapper()
            with self.subTest(path=relative,key=key),self.assertRaises(ValueError):self.validate(self.version,self.binding())
            path.write_bytes(original);self.stage(STAGES[2]);self.wrapper()
        original=self.predictions.read_bytes()
        for data in ['sample_id,SoftwareFixtureClock\nGSM999002,1\nGSM999001,2\n',
                     'sample_id,SoftwareFixtureClock\nGSM999001,NaN\nGSM999002,2\n']:
            self.predictions.write_text(data);self.stage(STAGES[4]);self.wrapper()
            with self.assertRaises(ValueError):self.validate(self.version,self.binding())
        self.predictions.write_bytes(original)

    def test_snapshot_change_and_resealed_change_after_export_plan_are_rejected(self):
        binding=self.root/'binding.json';self.write(binding,self.binding())
        self.export(self.version,self.root/'export',readback_path=binding,plan_path=self.root/'export-plan.json')
        self.predictions.write_text('sample_id,SoftwareFixtureClock\nGSM999001,3\nGSM999002,2\n')
        with self.assertRaises(PipelineIntegrityError):self.validate(self.version,json.loads(binding.read_text()))
        self.stage(STAGES[4]);self.wrapper();self.write(binding,self.binding())
        with self.assertRaises(PipelineIntegrityError):
            self.export(self.version,self.root/'export-2',readback_path=binding,plan_path=self.root/'export-plan.json')

    def test_old_partial_profile_cannot_accept_planned_completion(self):
        from audit_system.public_trace_2026_09_18 import _partial_future_snapshot
        binding=self.binding();binding['schema']='clockbase-partial-update-binding-v1'
        with self.assertRaises((ValueError,FileNotFoundError)):_partial_future_snapshot(self.version,binding)

    def test_materialized_candidate_adjustment_and_saved_design_cannot_diverge(self):
        for basename,key,value in [('contract-adjustments/001-candidate.json','status','PASS'),
            ('contract-adjustments/001-receipt.json','accepted_result',True),
            ('contract-adjustments/000-receipt.json','contract_status','UNCHANGED_VALID'),
            ('source-contract-admissibility.json','candidate_sha256','a'*64),
            ('design.json','unresolved',['Changed advisory scientific claim'])]:
            path=self.attempt/basename;original=path.read_bytes();record=json.loads(original);record[key]=value;self.write(path,record)
            self.stage(STAGES[2],record if basename=='design.json' else self.design);self.wrapper()
            with self.subTest(file=basename),self.assertRaises(ValueError):self.validate(self.version,self.binding())
            path.write_bytes(original);self.stage(STAGES[2],self.design);self.wrapper()

    def test_resealed_source_instruction_substitution_is_not_source_only(self):
        path=self.pipeline/'api/design-000/request.json';record=json.loads(path.read_text())
        record['messages'][0]['content']='Changed instructions with the same public metadata'
        self.write(path,record);self.write(self.first,record)
        freeze=json.loads(self.freeze.read_text());freeze['first_source_request_sha256']=sha(self.first);self.write(self.freeze,freeze)
        self.wrapper()
        with self.assertRaises(ValueError):self.validate(self.version,self.binding())

    def test_independent_readbacks_are_byte_bound_downloads_without_acceptance(self):
        binding=self.binding();binding['independent_readbacks']={}
        for name,key,digest in [('independent-source-readback','source_design_sha256',sha(self.attempt/'design.json')),
                               ('independent-clock-readback','source_predictions_sha256',sha(self.predictions))]:
            original=self.root/(name+'-original.json');self.write(original,{'disclosed_software_fixture':True})
            public=self.root/(name+'.json');record={'schema':'disclosed-software-readback-v1','study_id':self.study,
                'version_id':self.version.name,'source_version_sha256':sha(self.version/'version.json'),
                'source_receipt_sha256':sha(original),key:digest,'scientific_status':'CANNOT VERIFY','accepted_scientific_repairs':0}
            self.write(public,record);binding['independent_readbacks'][name]={'path':str(public),'sha256':sha(public),
                'source_receipt':{'path':str(original),'sha256':sha(original)}}
        binding_path=self.root/'readback-binding.json';self.write(binding_path,binding)
        self.export(self.version,self.root/'with-readbacks',readback_path=binding_path,plan_path=self.root/'readbacks-plan.json')
        for name,record in binding['independent_readbacks'].items():
            self.assertEqual((self.root/'with-readbacks/public'/(name+'.json')).read_bytes(),Path(record['path']).read_bytes())
        name='independent-clock-readback';record=binding['independent_readbacks'][name];path=Path(record['path'])
        for key,value in [('source_predictions_sha256','a'*64),('accepted_scientific_repairs',1),
                          ('scientific_status','PASS'),('scope','Private file /Entropy/secret/fixture.csv')]:
            original=path.read_bytes();payload=json.loads(original);payload[key]=value;self.write(path,payload);record['sha256']=sha(path)
            with self.subTest(key=key),self.assertRaises(ValueError):self.validate(self.version,binding)
            path.write_bytes(original)

    def test_zero_retry_success_does_not_invent_historical_validation(self):
        fixture=PlannedClockCheckpointTests();fixture.fixture_retries=0
        try:
            fixture.setUp();value=fixture.validate(fixture.version,fixture.binding())[4]
            self.assertEqual(len(value['api_calls']),1)
            self.assertEqual(value['source_response_history'][0]['recorded_validation_status'],'NOT_RECORDED')
            self.assertIsNone(value['source_response_history'][0]['validation_sha256'])
            self.assertFalse(value['source_response_history'][0]['retry_transmitted'])
        finally:fixture.doCleanups()


if __name__=='__main__': unittest.main()
