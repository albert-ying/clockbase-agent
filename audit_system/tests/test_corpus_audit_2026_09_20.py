"""Software fixtures for corpus orchestration; not empirical audit evidence."""
import copy
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from audit_system.corpus_audit_2026_09_20 import (
    Budget, evidence_index, selected, validate_check, validate_design,
    conditional_calculations, review_request, validate_review, historical_rows,
    DOMAINS, SCIENTIFIC_CHECKS, declaration, StudyClient,
    admissible_check, completed_pipeline_stage, validate_configuration,
    coverage_report, design_request, replay_source_design, verify_pipeline_followup,
    source_selector_feedback, checked_response, audit_pipeline_result, audit_study,
    review_row_scope_feedback, REVIEW_ROW_SCOPE_RULE, REVIEW_ROW_ID_INSTRUCTION,
    prepare_pipeline_followup, run_pipeline_followup, record_hash, TARGET, DESIGN_SYSTEM,
    prepare_selected_stage_recovery, run_selected_stage_recovery, REVIEW_SYSTEM,
)
from audit_system.agent_workflow_2026_09_10 import sha, source_packet
from src.pipeline_runtime_2026_09_10 import StageStore, PipelineIntegrityError


class CorpusTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.packet={'study_id':'GSE1','source_url':'https://example.invalid/software-fixture',
          'source_sha256':'a'*64,'samples':[
            {'sample_id':f'GSM{i+1}','condition':'treated' if i<3 else 'control',
             'donor':str(i%3), 'assay':'RNA','protocol':'same\nquoted "protocol"'} for i in range(6)]}
        self.index=evidence_index(self.packet)

    def tearDown(self):
        self.tmp.cleanup()

    def ref(self,field,value):
        return next(k for k,v in self.index['evidence'].items() if v['field']==field and v['value']==value)

    def check(self,status='CANNOT VERIFY'):
        return {'status':status,'reason':'Software fixture only','evidence_ids':[self.ref('assay','RNA')],
                'unresolved':[] if status=='PASS' else ['Unknown in this fixture']}

    def design(self):
        return {'study_id':'GSE1','checks':{k:self.check() for k in DOMAINS},
            'unresolved_targets':[], 'comparisons':[{'id':'paired','target':'fixture paired contrast',
             'case_all':[self.ref('condition','treated')], 'control_all':[self.ref('condition','control')],
             'mode':'inferential','method':'Paired t-test','unit_field':'donor','pair_field':'donor',
             'independence':self.check('PASS'),'covariates':[],
             'evidence_ids':[self.ref('assay','RNA')],'limitations':['Software fixture']} ]}

    def test_lossless_evidence_index_and_exact_membership(self):
        restored={sid:{} for sid in self.index['sample_ids']}
        for key,record in self.index['evidence'].items():
            if key.startswith('S'):
                for sid in record['sample_ids']:
                    restored[sid][record['field']]=record['value']
        self.assertEqual(restored,{s['sample_id']:{k:v for k,v in s.items() if k!='sample_id'} for s in self.packet['samples']})
        self.assertEqual(selected(self.index,[self.ref('condition','treated'),self.ref('donor','1')]),['GSM2'])
        with self.assertRaises(ValueError): selected(self.index,['L00000'])

    def test_source_identity_and_unknown_references_rejected(self):
        value=self.design();value['study_id']='GSE2'
        with self.assertRaises(ValueError):validate_design(value,self.index)
        value=self.design();value['comparisons'][0]['case_all']=['S99999']
        with self.assertRaises(ValueError):validate_design(value,self.index)
        value=self.design();value['comparisons'][0]['control_all']=value['comparisons'][0]['case_all']
        with self.assertRaises(ValueError):validate_design(value,self.index)

    def test_selector_feedback_reports_disjoint_constraints_without_repair(self):
        value=self.design()
        value['comparisons'][0]['case_all'] += [self.ref('donor','0'),self.ref('donor','1')]
        before=copy.deepcopy((value,self.index))
        feedback=source_selector_feedback(value,self.index)
        self.assertEqual(feedback['issue_count'],1)
        issue=feedback['issues'][0]
        self.assertEqual(issue['empty_arms'],['case'])
        self.assertEqual(issue['arms']['case']['selected_members']['count'],0)
        self.assertEqual(issue['arms']['control']['selected_members']['values'],['GSM4','GSM5','GSM6'])
        cells=issue['arms']['case']['all']['selectors']
        self.assertEqual(cells[-2]['source_members']['values'],['GSM1','GSM4'])
        self.assertEqual(cells[-1]['source_members']['values'],['GSM2','GSM5'])
        self.assertEqual((value,self.index),before)
        with self.assertRaisesRegex(ValueError,'Source groups are empty or overlap'):
            validate_design(value,self.index)

    def test_selector_feedback_accepts_valid_union_and_reports_full_overlap(self):
        value=self.design();comparison=value['comparisons'][0]
        comparison['case_any']=[self.ref('donor','0'),self.ref('donor','1')]
        self.assertEqual(selected(self.index,comparison['case_all'],comparison['case_any']),['GSM1','GSM2'])
        self.assertEqual(source_selector_feedback(value,self.index)['issue_count'],0)
        repeated=self.design();repeated['comparisons'][0]['case_any']=repeated['comparisons'][0]['case_all'][:]
        self.assertEqual(source_selector_feedback(repeated,self.index)['issue_count'],0)
        comparison['control_all']=[self.ref('donor','0')]
        issue=source_selector_feedback(value,self.index)['issues'][0]
        self.assertEqual(issue['overlap_members']['values'],['GSM1'])
        self.assertEqual(issue['arms']['control']['selected_members']['values'],['GSM1','GSM4'])

    def test_selector_feedback_handles_malformed_duplicate_and_unknown_refs(self):
        for refs in (None,'S00000',[{}],['S99999'],['L00000'],[self.ref('donor','0')]*2):
            with self.subTest(refs=refs):
                value=self.design();value['comparisons'][0]['case_all']=refs
                feedback=source_selector_feedback(value,self.index)
                self.assertEqual(feedback['issue_count'],1)
                self.assertIsNotNone(feedback['issues'][0]['arms']['case']['selection_error'])
        for value in (None,[],{}, {'comparisons':'bad'}, {'comparisons':[None]}):
            with self.subTest(value=value):
                feedback=source_selector_feedback(value,self.index)
                self.assertEqual(feedback['issue_count'],1)
                self.assertIn('structure_error',feedback['issues'][0])

    def test_selector_feedback_bounds_details_after_full_membership_selection(self):
        packet=copy.deepcopy(self.packet)
        packet['samples']=[{'sample_id':'GSM'+str(i+1),'condition':'same','long':'x'*2000}
                           for i in range(50)]
        index=evidence_index(packet)
        ref=next(k for k,v in index['evidence'].items() if v.get('field')=='condition')
        value={'comparisons':[{'id':str(i),'case_all':[ref],'control_all':[ref]} for i in range(25)]}
        feedback=source_selector_feedback(value,index)
        self.assertEqual(feedback['issue_count'],25)
        self.assertEqual(len(feedback['issues'])+feedback['issues_omitted'],25)
        self.assertLessEqual(len(feedback['issues']),20)
        self.assertLessEqual(len(json.dumps(feedback,ensure_ascii=False).encode()),feedback['display_byte_limit'])
        overlap=feedback['issues'][0]['overlap_members']
        self.assertEqual(overlap['count'],50)
        self.assertEqual(overlap['omitted_count'],38)
        self.assertEqual(len(overlap['values']),12)

    def test_selector_retry_retains_exact_feedback_and_stops_after_one_retry(self):
        bad=self.design();bad['comparisons'][0]['case_all'] += [self.ref('donor','0'),self.ref('donor','1')]
        wire=json.dumps(bad);calls=[]
        class Client:
            def ask(self,actor,system,payload):
                calls.append(copy.deepcopy(payload));return wire
        directory=self.root/'selector-retry';directory.mkdir()
        payload=design_request(self.index);before=copy.deepcopy(payload)
        with self.assertRaisesRegex(ValueError,'Source groups are empty or overlap'):
            checked_response(Client(),'source_design','unchanged initial system',payload,
                lambda value:validate_design(value,self.index),directory,1,
                feedback_builder=lambda value:source_selector_feedback(value,self.index))
        self.assertEqual(len(calls),2)
        self.assertEqual(calls[0],before)
        self.assertEqual(payload,before)
        first=json.loads((directory/'contract-00.json').read_text())
        self.assertEqual(calls[1]['contract_validation_feedback'],first)
        self.assertEqual(first['source_selector_diagnostics']['issue_count'],1)
        self.assertEqual(calls[1]['previous_response'],wire)
        self.assertTrue((directory/'contract-01.json').is_file())

    def test_selector_feedback_does_not_retry_transport_or_parse_failures_with_builder(self):
        from unittest.mock import Mock
        builder=Mock()
        client=Mock();client.ask.side_effect=RuntimeError('transport unavailable')
        directory=self.root/'selector-transport';directory.mkdir()
        with self.assertRaisesRegex(RuntimeError,'transport unavailable'):
            checked_response(client,'source_design','system',{},lambda value:value,directory,1,
                             feedback_builder=builder)
        self.assertEqual(client.ask.call_count,1);builder.assert_not_called()
        self.assertEqual(list(directory.iterdir()),[])
        client=Mock();client.ask.return_value='not json'
        with self.assertRaises(ValueError):
            checked_response(client,'source_design','system',{},lambda value:value,directory,1,
                             feedback_builder=builder)
        self.assertEqual(client.ask.call_count,2);builder.assert_not_called()

    def test_selector_diagnostic_failure_preserves_original_contract_error(self):
        from unittest.mock import Mock
        client=Mock();client.ask.return_value='{}'
        directory=self.root/'selector-diagnostic-failure';directory.mkdir()
        def reject(value):raise ValueError('Original source contract failure')
        def diagnostic(value):raise TypeError('Malformed diagnostic source')
        with self.assertRaisesRegex(ValueError,'Original source contract failure'):
            checked_response(client,'source_design','system',{},reject,directory,0,feedback_builder=diagnostic)
        record=json.loads((directory/'contract-00.json').read_text())
        self.assertEqual(record['error'],'ValueError: Original source contract failure')
        self.assertEqual(record['source_selector_diagnostics']['status'],'UNAVAILABLE')
        self.assertEqual(client.ask.call_count,1)

    def test_default_contract_feedback_keeps_existing_non_source_shape(self):
        from unittest.mock import Mock
        client=Mock();client.ask.return_value='{}'
        directory=self.root/'default-feedback';directory.mkdir()
        def reject(value):raise ValueError('Original review contract failure')
        with self.assertRaisesRegex(ValueError,'Original review contract failure'):
            checked_response(client,'review','system',{},reject,directory,0)
        self.assertEqual(json.loads((directory/'contract-00.json').read_text()),{
            'status':'REJECTED_RESPONSE_CONTRACT','attempt':0,
            'error':'ValueError: Original review contract failure','scientific_status':'CANNOT VERIFY'})

    def test_limitations_cannot_be_promoted_to_pass(self):
        value=self.check('PASS');value['evidence_ids']=['L00000']
        with self.assertRaises(ValueError):validate_check(value,self.index['evidence'],'test')
        value=self.check('PASS');value['unresolved']=['unknown units']
        with self.assertRaises(ValueError):validate_check(value,self.index['evidence'],'test')

    def test_aggregate_budget_atomic_and_uncertain_cost_retained(self):
        budget=Budget(self.root/'budget.sqlite3',1.0,100)
        def reserve(i):
            try: budget.reserve(str(i),.3);return True
            except RuntimeError:return False
        with ThreadPoolExecutor(max_workers=8) as pool:
            self.assertEqual(sum(pool.map(reserve,range(20))),3)
        self.assertAlmostEqual(budget.report()['estimated_or_reserved_usd'],.9)
        with budget.connect() as db:
            key=db.execute('SELECT id FROM reservations LIMIT 1').fetchone()[0]
        budget.settle(key,.1)
        self.assertAlmostEqual(budget.report()['estimated_or_reserved_usd'],.7)
        budget.reserve(key,.3)
        self.assertEqual(budget.report()['call_reservations'],3)
        with self.assertRaises(ValueError):Budget(self.root/'budget.sqlite3',2.0,100)

    def predictions(self,duplicate=False):
        path=self.root/'predictions.csv'
        with path.open('w',newline='') as handle:
            writer=csv.writer(handle);writer.writerow(['','clock','clock'] if duplicate else ['','clock'])
            for i,v in enumerate([5,7,11,2,4,6]):
                writer.writerow([f'GSM{i+1}',v,v] if duplicate else [f'GSM{i+1}',v])
        return path

    def test_paired_rerun_uses_units_not_row_order(self):
        design=validate_design(self.design(),self.index)
        result=conditional_calculations(design,self.packet,self.predictions(),['clock'])[0]
        self.assertEqual(result['status'],'PASS')
        self.assertEqual(result['calculation']['n_biological_units'],3)
        self.assertEqual(result['calculation']['within_unit_differences'],[3,3,5])
        self.assertEqual(result['scientific_status'],'CANNOT VERIFY')

    def test_repeated_units_and_ambiguous_clock_stop(self):
        value=self.design();value['comparisons'][0]['method']='Welch t-test'
        value['comparisons'][0]['pair_field']=None
        design=validate_design(value,self.index)
        result=conditional_calculations(design,self.packet,self.predictions(),['clock'])[0]
        self.assertEqual(result['status'],'CANNOT VERIFY')
        self.assertIn('Repeated biological units',result['reason'])
        result=conditional_calculations(validate_design(self.design(),self.index),self.packet,self.predictions(True),['clock'])[0]
        self.assertEqual(result['status'],'CANNOT VERIFY')
        self.assertIn('ambiguous',result['reason'])

    def test_review_must_account_for_every_row(self):
        design=validate_design(self.design(),self.index)
        rows=[{'row_id':'R1','fields':{'p_value':'.01'},'artifact_sha256':'b'*64,'logical_record':2,
               'identity_binding':{'identity_status':'PASS','bound_study_key':'human_rna:GSE1'}}]
        payload=review_request(self.index,design,rows,[],[])
        value={'study_id':'GSE1','checks':{k:{**self.check(),'scope':'fixture'} for k in SCIENTIFIC_CHECKS},
               'rows':[], 'design_disagreements':[]}
        with self.assertRaises(ValueError):validate_review(value,payload,rows)
        value['rows']=[{**self.check(),'row_id':'R1','action':'unresolved','source_comparison_ids':[],'proposed_repair':None}]
        self.assertEqual(validate_review(value,payload,rows)['contract_status'],'PASS')
        value['rows'][0]['source_comparison_ids']=['invented']
        with self.assertRaises(ValueError):validate_review(value,payload,rows)

    def review_scope_fixture(self, count=2):
        rows=[{'row_id':f'R{i:05d}','fields':{'p_value':'.01','clock':'fixture'},
               'artifact_sha256':'b'*64,'logical_record':i+2,
               'identity_binding':{'identity_status':'PASS','bound_study_key':'human_rna:GSE1'}}
              for i in range(count)]
        inventory=[{'row_id':f'R{i:05d}','artifact_id':'fixture','target':{'clock':'fixture'}}
                   for i in range(count+11)]
        payload=review_request(self.index,validate_design(self.design(),self.index),rows,[],[],
                               whole_study_inventory=inventory)
        value={'study_id':'GSE1','checks':{k:{**self.check(),'scope':'fixture'} for k in SCIENTIFIC_CHECKS},
               'rows':[{**self.check(),'row_id':r['row_id'],'action':'unresolved',
                        'source_comparison_ids':[],'proposed_repair':None} for r in rows],
               'design_disagreements':[]}
        return rows,inventory,payload,value

    def test_review_chunk_scope_keeps_whole_inventory_and_exact_row_bindings(self):
        rows,inventory,payload,value=self.review_scope_fixture(24)
        before=copy.deepcopy((self.index,rows,inventory,payload,value))
        self.assertEqual(payload['review_row_ids'],[r['row_id'] for r in rows])
        self.assertEqual(payload['review_row_scope'],REVIEW_ROW_SCOPE_RULE)
        self.assertEqual(payload['required_output']['rows'][0]['row_id'],REVIEW_ROW_ID_INSTRUCTION)
        self.assertEqual(payload['whole_study_saved_output_inventory'],inventory)
        self.assertEqual(len(inventory),35)
        evidence=payload['source']['evidence']
        self.assertEqual({k for k,v in evidence.items() if v['artifact']=='saved_analysis'},set(payload['review_row_ids']))
        for row in rows:
            self.assertEqual(evidence[row['row_id']],{'artifact':'saved_analysis','value':row['fields'],
                'artifact_sha256':row['artifact_sha256'],'logical_record':row['logical_record'],
                'identity_binding':row['identity_binding'],'sample_ids':[]})
        result=validate_review(value,payload,rows)
        self.assertEqual(result['scientific_acceptance'],'CANNOT VERIFY')
        self.assertTrue(all(r['status']=='CANNOT VERIFY' and r['action']=='unresolved' for r in result['rows']))
        self.assertEqual((self.index,rows,inventory,payload,value),before)

    def test_review_empty_scope_requires_empty_rows_and_all_scientific_checks(self):
        rows,inventory,payload,value=self.review_scope_fixture(0)
        self.assertEqual(payload['review_row_ids'],[])
        self.assertEqual(validate_review(value,payload,rows)['rows'],[])
        invented=copy.deepcopy(value);invented['rows']=[{'row_id':inventory[0]['row_id']}]
        with self.assertRaisesRegex(ValueError,'omitted, duplicated or invented'):
            validate_review(invented,payload,rows)
        del value['checks'][next(iter(SCIENTIFIC_CHECKS))]
        with self.assertRaisesRegex(ValueError,'every scientific check'):validate_review(value,payload,rows)

    def test_review_row_coverage_rejects_missing_extra_duplicate_and_malformed_without_filtering(self):
        rows,inventory,payload,value=self.review_scope_fixture()
        mutations=[value['rows'][:1],value['rows']+[{'row_id':inventory[-1]['row_id']}],
                   [value['rows'][0],value['rows'][0]],[value['rows'][0],{'row_id':'unknown'}],
                   [value['rows'][0],None],[value['rows'][0],['R00001']],
                   [value['rows'][0],{'row_id':['R00001']}],[value['rows'][0],{}],
                   [value['rows'][0],{'row_id':None}],[value['rows'][0],{'row_id':''}],None,{}]
        for entries in mutations:
            with self.subTest(entries=entries):
                bad={**value,'rows':copy.deepcopy(entries)};before=copy.deepcopy(bad)
                with self.assertRaisesRegex(ValueError,'omitted, duplicated or invented'):
                    validate_review(bad,payload,rows)
                self.assertEqual(bad,before)

    def test_review_controller_scope_cannot_rebind_original_rows(self):
        rows,inventory,payload,value=self.review_scope_fixture()
        for scope in (list(reversed(payload['review_row_ids'])),[],[inventory[-1]['row_id']]):
            with self.subTest(scope=scope),self.assertRaises(PipelineIntegrityError):
                validate_review(value,{**payload,'review_row_ids':scope},rows)
        with self.assertRaises(PipelineIntegrityError):validate_review(value,payload,rows+rows[:1])
        historical=copy.deepcopy(payload);del historical['review_row_ids'];del historical['review_row_scope']
        self.assertEqual(validate_review(value,historical,rows),validate_review(value,payload,rows))

    def test_review_row_scope_diagnostics_separate_malformed_and_duplicate_ids(self):
        value={'rows':[{'row_id':'R1'},{'row_id':'R1'},{'row_id':'R3'},None,{},
                       {'row_id':['R2']},{'row_id':''}]}
        before=copy.deepcopy(value)
        result=review_row_scope_feedback(value,['R1','R2'])
        self.assertEqual(result['missing_row_ids']['values'],['R2'])
        self.assertEqual(result['extra_row_ids']['values'],['R3'])
        self.assertEqual(result['duplicate_row_ids']['values'],['R1'])
        self.assertEqual(result['duplicate_entry_count'],1)
        self.assertEqual(result['malformed_entry_count'],4)
        self.assertEqual([v['index'] for v in result['malformed_entries']],[3,4,5,6])
        self.assertFalse(result['exact_row_coverage']);self.assertEqual(value,before)
        self.assertEqual(result['scientific_status'],'CANNOT VERIFY')
        self.assertTrue(review_row_scope_feedback({'rows':[]},[])['exact_row_coverage'])
        for malformed in (None,[],{}, {'rows':None},{'rows':{}}):
            with self.subTest(malformed=malformed):
                self.assertIsNotNone(review_row_scope_feedback(malformed,[])['rows_structure_error'])

    def test_review_row_scope_diagnostics_count_full_input_before_bounding_display(self):
        value={'rows':[{'row_id':'x'*2000+str(i)} for i in range(100)]+[None]*100}
        result=review_row_scope_feedback(value,[f'R{i}' for i in range(100)])
        self.assertEqual(result['missing_row_ids']['count'],100)
        self.assertEqual(result['extra_row_ids']['count'],100)
        self.assertEqual(result['extra_row_ids']['omitted_count'],76)
        self.assertEqual(result['malformed_entry_count'],100)
        self.assertEqual(result['malformed_entries_omitted'],76)
        self.assertLess(len(json.dumps(result).encode()),16000)
        for bad in (None,'R1',['R1','R1'],[None],[[]],['']):
            with self.subTest(bad=bad),self.assertRaises(ValueError):review_row_scope_feedback({},bad)

    def test_review_retry_reports_scope_only_and_retains_rejected_text(self):
        rows,inventory,payload,value=self.review_scope_fixture()
        bad=copy.deepcopy(value);bad['rows'].append({**bad['rows'][0],'row_id':inventory[-1]['row_id']})
        wire=json.dumps(bad);calls=[]
        class Client:
            def ask(self,actor,system,request):calls.append(copy.deepcopy(request));return wire
        directory=self.root/'review-retry';directory.mkdir()
        with self.assertRaisesRegex(ValueError,'omitted, duplicated or invented'):
            checked_response(Client(),'review','unchanged system',payload,
                lambda v:validate_review(v,payload,rows),directory,1,
                feedback_builder=lambda v:review_row_scope_feedback(v,payload['review_row_ids']),
                feedback_key='review_row_scope_diagnostics')
        self.assertEqual(len(calls),2);self.assertEqual(calls[0],payload)
        record=json.loads((directory/'contract-00.json').read_text())
        self.assertNotIn('source_selector_diagnostics',record)
        self.assertEqual(record['review_row_scope_diagnostics']['extra_row_ids']['values'],[inventory[-1]['row_id']])
        self.assertEqual(calls[1]['contract_validation_feedback'],record)
        self.assertEqual(calls[1]['previous_response'],wire)
        self.assertEqual(calls[1]['source'],payload['source'])
        self.assertEqual(calls[1]['frozen_source_design'],payload['frozen_source_design'])
        self.assertTrue((directory/'contract-01.json').is_file())

    def test_review_diagnostic_failure_and_transport_keep_original_errors(self):
        from unittest.mock import Mock
        client=Mock();client.ask.return_value='{}'
        directory=self.root/'review-diagnostic-failure';directory.mkdir()
        def reject(value):raise ValueError('Original review contract failure')
        def diagnostic(value):raise TypeError('Malformed diagnostic controller scope')
        with self.assertRaisesRegex(ValueError,'Original review contract failure'):
            checked_response(client,'review','system',{},reject,directory,0,feedback_builder=diagnostic,
                             feedback_key='review_row_scope_diagnostics')
        record=json.loads((directory/'contract-00.json').read_text())
        self.assertEqual(record['error'],'ValueError: Original review contract failure')
        self.assertEqual(record['review_row_scope_diagnostics']['status'],'UNAVAILABLE')
        directory=self.root/'review-transport';directory.mkdir()
        client=Mock();client.ask.side_effect=RuntimeError('TRANSPORT_TIMEOUT')
        builder=Mock()
        with self.assertRaisesRegex(RuntimeError,'TRANSPORT_TIMEOUT'):
            checked_response(client,'review','system',{},reject,directory,1,feedback_builder=builder,
                             feedback_key='review_row_scope_diagnostics')
        self.assertEqual(client.ask.call_count,1);builder.assert_not_called()
        self.assertEqual(list(directory.iterdir()),[])

    def test_review_controller_integrity_error_never_enters_contract_retry(self):
        from unittest.mock import Mock
        rows,inventory,payload,value=self.review_scope_fixture()
        payload['review_row_ids']=[]
        client=Mock();client.ask.return_value=json.dumps(value)
        builder=Mock();directory=self.root/'review-integrity';directory.mkdir()
        with self.assertRaisesRegex(PipelineIntegrityError,'original row bindings'):
            checked_response(client,'review','system',payload,
                lambda v:validate_review(v,payload,rows),directory,1,
                feedback_builder=builder,feedback_key='review_row_scope_diagnostics')
        self.assertEqual(client.ask.call_count,1);builder.assert_not_called()
        self.assertEqual(list(directory.iterdir()),[])

    def test_ragged_saved_row_is_retained(self):
        path=self.root/'result.csv';path.write_text('gse_id,clock,p_value\nGSE1,clock\nGSE1,clock,0.2\n')
        artifact={**declaration(path),'path':'result.csv','artifact_id':'artifact1'}
        rows,issues=historical_rows([artifact],self.root)
        self.assertEqual(len(rows),2);self.assertIsNone(rows[0]['fields'])
        self.assertEqual(rows[0]['logical_record'],2)
        path.write_text('changed')
        with self.assertRaises(ValueError):historical_rows([artifact],self.root)

    def test_prompt_reservation_includes_nested_json_escaping(self):
        config={'api':{'model':'gpt-5.6-luna','input_per_million':.2,'output_per_million':1.2,
          'output_tokens':10,'prompt_bytes':10000,'max_calls_per_study':2,'max_cost_per_study':1.0}}
        budget=Budget(self.root/'budget.sqlite3',1.,2)
        output=self.root/'study';output.mkdir()
        body={'model':'gpt-5.6-luna','usage':{'prompt_tokens':10,'completion_tokens':3},
              'choices':[{'finish_reason':'stop','message':{'content':'{}','refusal':None}}]}
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
             patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                   return_value=SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')):
            StudyClient(output,config,budget).ask('design','fixture',{'x':'"\n'*100})
        raw=json.dumps({'model':'gpt-5.6-luna','messages':[{'role':'system','content':'fixture'},
            {'role':'user','content':json.dumps({'x':'"\n'*100},ensure_ascii=False,separators=(',',':'))}],
            'max_completion_tokens':10,'reasoning_effort':'low','store':False},ensure_ascii=False).encode()
        expected=((len(raw)+200)*.2+12)/1e6
        reservation=json.loads((output/'api/design-000/reservation.json').read_text())
        aggregate=json.loads((output/'api/design-000.aggregate-binding.json').read_text())
        self.assertAlmostEqual(reservation['reserved_usd'],expected)
        self.assertAlmostEqual(aggregate['reserved_usd'],expected)
        self.assertAlmostEqual(budget.report()['estimated_or_reserved_usd'],(10*.2+3*1.2)/1e6)

    def test_gpt6_configuration_keeps_provider_and_budget_validation(self):
        config={'provider':'openai','concurrency':1,'review_rows_per_call':1,
            'max_total_calls':2,'max_total_cost_usd':.10,'validation_retries':0,
            'api':{'model':'gpt-6-luna','input_per_million':.20,'output_per_million':.75,
                   'output_tokens':100,'prompt_bytes':10000,'max_calls_per_study':2,
                   'max_cost_per_study':.10,'json_response':True}}
        validate_configuration(config)
        for provider,model in [('openrouter','gpt-6-luna'),('openai','gpt-6-astra'),
                               ('openai','gpt-6-luna-alias')]:
            with self.subTest(provider=provider,model=model):
                changed=copy.deepcopy(config);changed['provider']=provider;changed['api']['model']=model
                with self.assertRaisesRegex(ValueError,'no provider or model fallback'):
                    validate_configuration(changed)
        for field in ('input_per_million','output_per_million','max_cost_per_study'):
            changed=copy.deepcopy(config);changed['api'][field]=0
            with self.assertRaisesRegex(ValueError,'Positive finite API budget and tariff'):
                validate_configuration(changed)

    def test_gpt6_transport_records_requested_returned_model_and_rejects_mismatch(self):
        for returned in ('gpt-6-luna','gpt-5.6-luna'):
            with self.subTest(returned=returned):
                output=self.root/returned;output.mkdir()
                config={'api':{'model':'gpt-6-luna','input_per_million':.20,
                    'output_per_million':.75,'output_tokens':100,'prompt_bytes':10000,
                    'max_calls_per_study':1,'max_cost_per_study':.10,'json_response':True}}
                budget=Budget(output/'budget.sqlite3',.10,1)
                body={'model':returned,'usage':{'prompt_tokens':10,'completion_tokens':3},
                    'choices':[{'finish_reason':'stop','message':{'content':'{}'}}]}
                client=StudyClient(output,config,budget)
                with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                     patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                         return_value=SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')) as transport:
                    for _ in range(2):
                        if returned=='gpt-6-luna':
                            self.assertEqual(client.ask('design','Software fixture JSON',{}),'{}')
                        else:
                            with self.assertRaisesRegex(PipelineIntegrityError,'Returned model differs'):
                                client.ask('design','Software fixture JSON',{})
                    self.assertEqual(transport.call_count,1)
                slot=output/'api/design-000'
                request=json.loads((slot/'request.json').read_text())
                receipt=json.loads((slot/'receipt.json').read_text())
                self.assertEqual(request['model'],'gpt-6-luna')
                self.assertEqual(request['reasoning_effort'],'low')
                self.assertEqual(request['response_format'],{'type':'json_object'})
                self.assertNotIn('tools',request)
                self.assertEqual(receipt['requested_model'],'gpt-6-luna')
                self.assertEqual(receipt['returned_model'],returned)
                self.assertAlmostEqual(receipt['estimated_cost_usd'],(10*.20+3*.75)/1e6)
                self.assertIsNone(receipt['billed_cost_usd'])
                self.assertEqual((output/'design-response.txt').exists(),returned=='gpt-6-luna')

    def test_partial_rejection_retains_other_judgments_and_raw_response(self):
        value=self.design();value['checks']['data_provenance']=self.check('PASS')
        value['checks']['data_provenance']['unresolved']=['Missing acquisition history']
        original=copy.deepcopy(value)
        result=validate_design(value,self.index)
        self.assertEqual(value,original)
        self.assertEqual(result['contract_status'],'PARTIAL_WITH_CONSERVATIVE_REJECTIONS')
        self.assertEqual(result['checks']['data_provenance']['status'],'CANNOT VERIFY')
        self.assertEqual(result['comparisons'][0]['execution_design_status'],'PASS')
        row={**self.check('PASS'),'action':'retain','source_comparison_ids':['paired'],
             'proposed_repair':'inadmissible proposal','unresolved':['Unknown lineage']}
        result=admissible_check(row,self.index['evidence'],'row')
        self.assertEqual(result['model_status'],'PASS')
        self.assertEqual(result['model_action'],'retain')
        self.assertEqual(result['action'],'unresolved')
        self.assertEqual(result['source_comparison_ids'],[])
        self.assertIsNone(result['proposed_repair'])
        self.assertEqual(row['action'],'retain')

    def test_unknown_citation_rejects_only_affected_check(self):
        value=self.design();value['checks']['data_provenance']['evidence_ids']=['S99999']
        result=validate_design(value,self.index)
        check=result['checks']['data_provenance']
        self.assertEqual(check['admissibility'],'REJECTED_EVIDENCE')
        self.assertEqual(check['evidence_ids'],['L00000'])
        self.assertEqual(result['comparisons'][0]['execution_design_status'],'PASS')

    def test_unbound_trace_never_resolves_scientific_claim(self):
        trace={'identity_scope':'UNBOUND_CONTEXT','evidence':[{'evidence_id':'T00000','safe_text':'Software fixture'}]}
        payload=review_request(self.index,validate_design(self.design(),self.index),[],[],[],trace_evidence=trace)
        check={**self.check('PASS'),'evidence_ids':['T00000'],'scope':'causal interpretation'}
        value={'study_id':'GSE1','checks':{k:{**self.check(),'scope':'fixture'} for k in SCIENTIFIC_CHECKS},
               'rows':[],'design_disagreements':[]}
        value['checks']['scientific_claim_scope']=check
        self.assertIn('scientific_claim_scope',SCIENTIFIC_CHECKS)
        result=validate_review(value,payload,[])
        self.assertEqual(result['checks']['scientific_claim_scope']['status'],'CANNOT VERIFY')
        self.assertIn('trace',result['checks']['scientific_claim_scope']['contract_findings'][0])
        # Combining an unbound trace with a valid source citation cannot hide the dependency.
        check['evidence_ids'].append(self.ref('assay','RNA'))
        self.assertEqual(admissible_check(check,payload['source']['evidence'],'test')['status'],'CANNOT VERIFY')

    def test_independent_verification_has_citable_source_bindings(self):
        view={'results':[{'comparison_id':'paired','clock':'clock','status':'PASS'}],
              'request_sha256':'b'*64,'predictions_sha256':'c'*64,'source_hash_status':'PASS'}
        payload=review_request(self.index,validate_design(self.design(),self.index),[],[],[],independent_numerics=view)
        self.assertEqual(payload['source']['evidence']['V00000']['predictions_sha256'],'c'*64)
        self.assertEqual(payload['source']['evidence']['V00000']['scientific_status'],'CANNOT VERIFY')

    def measurement(self):
        return {'study_id':'GSE1','modality':'human_rna','evidence':{
            'M00000':{'artifact':'molecular_input','sample_ids':self.index['sample_ids'],
                      'value':{'feature_count':3,'source_sha256':'f'*64}},
            'M00001':{'artifact':'clock_execution','sample_ids':self.index['sample_ids'],
                      'value':{'clock':'SOFTWARE_FIXTURE_CLOCK','units':'CANNOT VERIFY'}}}}

    def test_measurement_review_preserves_the_blinded_index(self):
        before=copy.deepcopy(self.index)
        design_before=design_request(self.index)
        payload=review_request(self.index,validate_design(self.design(),self.index),[],[],[],
                               measurement_evidence=self.measurement())
        self.assertEqual(self.index,before)
        self.assertEqual(design_request(self.index),design_before)
        self.assertEqual(payload['source']['evidence']['L00000'],before['evidence']['L00000'])
        self.assertIn('L00001',payload['source']['evidence'])
        self.assertEqual(len(payload['measurement_evidence_sha256']),64)
        self.assertNotIn('SOFTWARE_FIXTURE_CLOCK',json.dumps(design_before))
        self.assertFalse(any(key.startswith('M') for key in design_before['source']['evidence']))

    def test_measurement_registry_rejects_identity_and_identifier_collisions(self):
        malformed=self.measurement();malformed['evidence']=list(malformed['evidence'].values())
        with self.assertRaisesRegex(ValueError,'identifier-to-record mapping'):
            review_request(self.index,{},[],[],[],measurement_evidence=malformed)
        value=self.measurement();value['study_id']='GSE2'
        with self.assertRaises(PipelineIntegrityError):
            review_request(self.index,{},[],[],[],measurement_evidence=value)
        value=self.measurement();value['evidence']['S00000']={}
        with self.assertRaises(PipelineIntegrityError):
            review_request(self.index,{},[],[],[],measurement_evidence=value)

    def test_future_measurement_judgments_require_method_and_context(self):
        payload=review_request(self.index,validate_design(self.design(),self.index),[],[],[],
                               measurement_evidence=self.measurement())
        value={'study_id':'GSE1','checks':{k:{**self.check(),'scope':'fixture'} for k in SCIENTIFIC_CHECKS},
               'rows':[],'design_disagreements':[]}
        check={**self.check('PASS'),'scope':'fixture measurement only'}
        value['checks']['measurement_applicability']=check
        with self.assertRaisesRegex(ValueError,'lacks supplied method evidence'):
            validate_review(value,payload,[])
        check['evidence_ids']=['M00000']
        with self.assertRaisesRegex(ValueError,'lacks source evidence'):
            validate_review(value,payload,[])
        check['evidence_ids'].append(self.ref('assay','RNA'))
        self.assertEqual(validate_review(value,payload,[])['contract_status'],'PASS')
        value['checks']['clock_target_and_scale']={**self.check('PASS'),'scope':'fixture clock'}
        with self.assertRaisesRegex(ValueError,'lacks supplied method evidence'):
            validate_review(value,payload,[])

    def test_measurement_presence_never_establishes_biological_units(self):
        payload=review_request(self.index,validate_design(self.design(),self.index),[],[],[],
                               measurement_evidence=self.measurement())
        value={'study_id':'GSE1','checks':{k:{**self.check(),'scope':'fixture'} for k in SCIENTIFIC_CHECKS},
               'rows':[],'design_disagreements':[]}
        value['checks']['biological_replication']={**self.check('PASS'),'scope':'fixture',
                                                  'evidence_ids':['M00000']}
        with self.assertRaisesRegex(ValueError,'lacks source evidence'):
            validate_review(value,payload,[])

    def replay_binding(self):
        root=self.root/'prior';store=StageStore(root,{'software_fixture':True})
        store.run('01_source_index',lambda _:self.index)
        design=validate_design(self.design(),self.index)
        store.run('02_blinded_design',lambda _:design)
        return {'root':str(root),'files':{str(p):declaration(p) for p in root.rglob('*') if p.is_file()}},design

    def test_identical_source_design_replay_keeps_original_timestamp_and_groups(self):
        binding,design=self.replay_binding()
        with patch('audit_system.corpus_audit_2026_09_20.StudyClient.ask',side_effect=AssertionError('No new model assessment')):
            self.assertEqual(replay_source_design(binding,self.index),design)
        changed=copy.deepcopy(self.index);changed['source_sha256']='b'*64
        with self.assertRaisesRegex(PipelineIntegrityError,'identical blinded source index'):
            replay_source_design(binding,changed)

    def test_source_label_pairs_require_controller_policy_and_source_owned_context(self):
        packet={'study_id':'GSE1','source_url':'https://example.invalid/fixture',
                'source_sha256':'a'*64,'samples':[]}
        for time in ('1HR','3HR'):
            for arm in ('CTRL','LTP'):
                for animal in (1,2):
                    packet['samples'].append({'sample_id':f'GSM{len(packet["samples"])+1}',
                        'title':f'{time}_{arm}_ANIMAL_{animal}', 'source_name_ch1':'hippocampus',
                        'data_processing':'animal and treatment modeled within each time point'})
        index=evidence_index(packet)
        def ref(field,value):
            return next(key for key,cell in index['evidence'].items()
                        if cell['field']==field and cell['value']==value)
        context=ref('data_processing','animal and treatment modeled within each time point')
        value={'study_id':'GSE1','checks':{name:{'status':'CANNOT VERIFY','reason':'fixture',
                    'evidence_ids':[context],'unresolved':['Scientific adequacy unresolved']}
                    for name in DOMAINS},'unresolved_targets':[],
               'comparisons':[{'id':'one_hour','target':'1HR LTP minus CTRL',
                    'case_all':[ref('source_name_ch1','hippocampus')],
                    'case_any':[ref('title',f'1HR_LTP_ANIMAL_{i}') for i in (1,2)],
                    'control_all':[ref('source_name_ch1','hippocampus')],
                    'control_any':[ref('title',f'1HR_CTRL_ANIMAL_{i}') for i in (1,2)],
                    'mode':'inferential','method':'Paired t-test','unit_field':'title','pair_field':'title',
                    'independence':{'status':'PASS','reason':'software fixture only',
                                    'evidence_ids':[context],'unresolved':[]},
                    'covariates':[],'evidence_ids':[context],
                    'limitations':['Source label alone cannot establish physical pairing'],
                    'pair_identity':{'schema':'clockbase-delimited-pair-v1','delimiter':'_',
                        'token_count':4,'fixed_tokens':{'2':'ANIMAL'},'stratum_tokens':{'0':'1HR'},
                        'unit_tokens':[2,3],'arm_token':1,'case_token':'LTP','control_token':'CTRL',
                        'scope_filters':{'source_name_ch1':'hippocampus'}},
                    'pair_identity_evidence_ids':[context]}]}
        before=copy.deepcopy(value)
        with self.assertRaisesRegex(ValueError,'not enabled by the controller'):
            validate_design(value,index)
        self.assertEqual(value,before)
        self.assertNotIn('source_design_policy',design_request(index))
        self.assertIn('pair_identity',design_request(index,source_design_policy='source_label_pairs_v1')
                     ['required_output']['comparisons'][0])
        result=validate_design(value,index,source_design_policy='source_label_pairs_v1')
        self.assertEqual(result['source_design_policy'],'source_label_pairs_v1')
        comparison=result['comparisons'][0]
        self.assertEqual(comparison['execution_design_status'],'PASS')
        self.assertEqual(len(comparison['numeric_design']['pairs']),2)
        self.assertEqual(len(comparison['numeric_design']['pair_identity']['excluded_samples']),4)
        self.assertEqual(comparison['numeric_design']['pair_identity']['scientific_status'],'CANNOT VERIFY')
        response_policy=copy.deepcopy(value);response_policy['source_design_policy']='source_label_pairs_v1'
        with self.assertRaisesRegex(ValueError,'Model response cannot declare'):
            validate_design(response_policy,index,source_design_policy='source_label_pairs_v1')
        document_index=copy.deepcopy(index)
        document_index['evidence']['D00000']={
            'artifact':'source_document', 'sample_ids':index['sample_ids'],
            'fields':{'overall_design':'Software fixture for animal-matched slices'}}
        document_context=copy.deepcopy(value)
        document_context['comparisons'][0]['pair_identity_evidence_ids']=[context,'D00000']
        with self.assertRaisesRegex(ValueError,'source-owned context evidence'):
            validate_design(document_context,document_index,source_design_policy='source_label_pairs_v1')
        wildcard_scope=copy.deepcopy(value)
        wildcard_scope['comparisons'][0]['pair_identity']['scope_filters']={'title':'1HR_*_ANIMAL_*'}
        original_wildcard=copy.deepcopy(wildcard_scope)
        stopped=validate_design(wildcard_scope,index,source_design_policy='source_label_pairs_v1')
        self.assertEqual(stopped['comparisons'][0]['execution_design_status'],'CANNOT VERIFY')
        self.assertIn('Derived complete source groups differ',stopped['comparisons'][0]['execution_design_reason'])
        self.assertEqual(wildcard_scope,original_wildcard)
        no_extra_scope=copy.deepcopy(value)
        no_extra_scope['comparisons'][0]['pair_identity']['scope_filters']={}
        unfiltered=validate_design(no_extra_scope,index,source_design_policy='source_label_pairs_v1')
        self.assertEqual(unfiltered['comparisons'][0]['numeric_design']['pairs'],
                         comparison['numeric_design']['pairs'])
        wrong_context=copy.deepcopy(value)
        wrong_context['comparisons'][0]['pair_identity_evidence_ids']=[ref('title','1HR_LTP_ANIMAL_1')]
        with self.assertRaisesRegex(ValueError,'does not cover every selected sample'):
            validate_design(wrong_context,index,source_design_policy='source_label_pairs_v1')
        titles_only=copy.deepcopy(value)
        titles_only['comparisons'][0]['pair_identity_evidence_ids']=[
            ref('title',f'1HR_{arm}_ANIMAL_{i}') for arm in ('CTRL','LTP') for i in (1,2)]
        with self.assertRaisesRegex(ValueError,'contextual fields beyond'):
            validate_design(titles_only,index,source_design_policy='source_label_pairs_v1')
        wrong_method=copy.deepcopy(value)
        wrong_method['comparisons'][0]['method']='Welch t-test'
        with self.assertRaisesRegex(ValueError,'declared paired test'):
            validate_design(wrong_method,index,source_design_policy='source_label_pairs_v1')
        omitted=copy.deepcopy(value)
        omitted['comparisons'][0]['case_any']=omitted['comparisons'][0]['case_any'][:1]
        incomplete=validate_design(omitted,index,source_design_policy='source_label_pairs_v1')
        self.assertEqual(incomplete['comparisons'][0]['execution_design_status'],'CANNOT VERIFY')
        self.assertIn('Derived complete source groups differ',incomplete['comparisons'][0]['execution_design_reason'])

    def test_empty_pair_evidence_does_not_claim_a_derived_pair_rule(self):
        # Regression for archived outputs with pair_identity:null and [].
        # These synthetic source fixtures do not establish a biological design.
        for policy in ('literal_fields_v1', 'source_label_pairs_v1'):
            for mode, method in (('inferential', 'Paired t-test'),
                                 ('descriptive', 'descriptive'), ('unresolved', 'unsupported')):
                with self.subTest(policy=policy, mode=mode):
                    value=self.design();item=value['comparisons'][0]
                    item.update(mode=mode, method=method, pair_identity=None,
                                pair_identity_evidence_ids=[])
                    if mode!='inferential':
                        item.update(pair_field=None, unit_field=None, independence=self.check())
                    before=copy.deepcopy(value)
                    result=validate_design(value,self.index,source_design_policy=policy)
                    comparison=result['comparisons'][0]
                    self.assertEqual(comparison['case_ids'],['GSM1','GSM2','GSM3'])
                    self.assertEqual(comparison['control_ids'],['GSM4','GSM5','GSM6'])
                    self.assertEqual(value,before)
                    self.assertEqual(result['checks']['biological_units']['status'],'CANNOT VERIFY')
                    if mode=='inferential':
                        self.assertEqual(len(comparison['numeric_design']['pairs']),3)
                    else:
                        self.assertNotIn('numeric_design',comparison)
                    if mode=='unresolved':
                        self.assertEqual(comparison['execution_design_status'],'CANNOT VERIFY')

    def test_nonempty_or_malformed_pair_evidence_still_requires_a_rule(self):
        for evidence in ([self.ref('assay','RNA')], '', {}, False, 0, 'S00000'):
            with self.subTest(evidence=evidence):
                value=self.design()
                value['comparisons'][0].update(pair_identity=None, pair_identity_evidence_ids=evidence)
                with self.assertRaisesRegex(ValueError,'without a pairing rule'):
                    validate_design(value,self.index,source_design_policy='source_label_pairs_v1')
        for evidence in ([], None):
            with self.subTest(rule_evidence=evidence):
                value=self.design()
                value['comparisons'][0].update(pair_identity={'schema':'clockbase-delimited-pair-v1'},
                                              pair_identity_evidence_ids=evidence)
                with self.assertRaisesRegex(ValueError,'source-owned context evidence'):
                    validate_design(value,self.index,source_design_policy='source_label_pairs_v1')

    def test_source_label_replay_cannot_change_controller_policy(self):
        root=self.root/'policy-replay'
        store=StageStore(root,{'config':{'source_design_policy':'source_label_pairs_v1'}})
        store.run('01_source_index',lambda _:self.index)
        design=validate_design(self.design(),self.index,source_design_policy='source_label_pairs_v1')
        store.run('02_blinded_design',lambda _:design)
        binding={'root':str(root),'files':{str(path):declaration(path)
                 for path in root.rglob('*') if path.is_file()}}
        self.assertEqual(replay_source_design(binding,self.index,
            source_design_policy='source_label_pairs_v1'),design)
        with self.assertRaisesRegex(PipelineIntegrityError,'controller policy differs'):
            replay_source_design(binding,self.index)
        unbound=copy.deepcopy(binding)
        del unbound['files'][str(root/'config.json')]
        with self.assertRaisesRegex(PipelineIntegrityError,'lacks its archived controller policy binding'):
            replay_source_design(unbound,self.index,source_design_policy='source_label_pairs_v1')
        redirected=copy.deepcopy(binding)
        redirected['files'][str(root/'config.json')]=declaration(root/'01_source_index/complete.json')
        with self.assertRaisesRegex(PipelineIntegrityError,'exact archived controller policy binding'):
            replay_source_design(redirected,self.index,source_design_policy='source_label_pairs_v1')

    def test_future_audit_rejects_request_policy_different_from_sealed_version(self):
        from src.update_pipeline_2026_09_20 import seal
        version_root=self.root/'future-version'
        version_root.mkdir()
        seal(version_root/'version.json',{'comprehensive_audit_policy':'source_label_pairs_v1'})
        with self.assertRaisesRegex(PipelineIntegrityError,'policy differs from the sealed update version'):
            audit_pipeline_result({'study_id':'GSE1','modality':'human_rna',
                'pipeline_root':str(version_root/'pipeline'),
                'comprehensive_audit_policy':'literal_fields_v1'})

    def test_changed_design_receipt_cannot_be_reused(self):
        binding,_=self.replay_binding()
        (Path(binding['root'])/'02_blinded_design/complete.json').write_text('{}')
        with self.assertRaises(PipelineIntegrityError):replay_source_design(binding,self.index)

    def test_followup_freeze_rejects_implementation_or_input_drift(self):
        from src.update_pipeline_2026_09_20 import seal
        code=self.root/'code';code.mkdir();(code/'fixture.py').write_text('# software fixture')
        output=self.root/'followup';(output/'implementation').mkdir(parents=True)
        (output/'implementation/fixture.py').write_bytes((code/'fixture.py').read_bytes())
        source=self.root/'source.csv';source.write_text('software fixture')
        frozen={'schema':'clockbase-pipeline-audit-followup-v1',
                'request':{'output':str(output/'audit')},
                'implementation_sha256':{'fixture.py':declaration(code/'fixture.py')['sha256']},
                'input_files':{'source':declaration(source)}}
        seal(output/'followup.json',frozen)
        with patch('audit_system.corpus_audit_2026_09_20.ROOT',code), patch('audit_system.corpus_audit_2026_09_20.CODE_FILES',('fixture.py',)):
            self.assertEqual(verify_pipeline_followup(output/'followup.json'),frozen)
            source.write_text('changed fixture')
            with self.assertRaises(PipelineIntegrityError):verify_pipeline_followup(output/'followup.json')
            source.write_text('software fixture');(code/'fixture.py').write_text('# changed')
            with self.assertRaises(PipelineIntegrityError):verify_pipeline_followup(output/'followup.json')

    def test_pipeline_consumed_file_must_be_in_verified_stage(self):
        root=self.root/'pipeline';store=StageStore(root,{'software_fixture':True})
        external=self.root/'elsewhere.csv';external.write_text('not a stage output')
        def wrong(directory):
            (directory/'metadata.csv').write_text('fixture')
            return {'metadata':str(external)}
        store.run('01_metadata',wrong)
        with self.assertRaisesRegex(PipelineIntegrityError,'not a verified stage output'):
            completed_pipeline_stage(root,'01_metadata',('metadata',))
        def valid(directory):
            path=directory/'clock.csv';path.write_text('fixture')
            return {'path':str(path)}
        store.run('05_clock_predictions',valid)
        self.assertTrue(Path(completed_pipeline_stage(root,'05_clock_predictions',('path',))['path']).is_file())

    def test_source_document_evidence_excludes_local_paths(self):
        packet=copy.deepcopy(self.packet)
        packet['study_documents']=[{'source_id':'fixture','source_kind':'geo_series',
            'source_url':'https://example.invalid/fixture','source_sha256':'e'*64,
            'fields':{'overall_design':'matched donors'},'sample_ids':['GSM1'],
            'local_path':'/private/software_fixture'}]
        index=evidence_index(packet)
        self.assertEqual(index['evidence']['D00000']['fields']['overall_design'],'matched donors')
        self.assertNotIn('local_path',index['evidence']['D00000'])

    def test_row_with_unresolved_study_identity_cannot_be_retained(self):
        design=validate_design(self.design(),self.index)
        rows=[{'row_id':'R1','fields':{'p_value':'.01'},'artifact_sha256':'b'*64,'logical_record':2,
               'identity_binding':{'identity_status':'CANNOT VERIFY','bound_study_key':None}}]
        payload=review_request(self.index,design,rows,[],[])
        value={'study_id':'GSE1','checks':{k:{**self.check(),'scope':'fixture'} for k in SCIENTIFIC_CHECKS},
               'rows':[{**self.check('PASS'),'row_id':'R1','action':'retain','source_comparison_ids':['paired'],
                        'proposed_repair':None}],'design_disagreements':[]}
        result=validate_review(value,payload,rows)
        self.assertEqual(result['rows'][0]['status'],'CANNOT VERIFY')
        self.assertEqual(result['rows'][0]['action'],'unresolved')
        self.assertEqual(result['rows'][0]['model_status'],'PASS')
        rows[0].pop('identity_binding')
        payload=review_request(self.index,design,rows,[],[])
        result=validate_review(value,payload,rows)
        self.assertEqual(result['rows'][0]['status'],'CANNOT VERIFY')
        self.assertEqual(result['rows'][0]['action'],'unresolved')

    def test_coverage_keeps_unassigned_artifacts_and_unreviewed_rows(self):
        artifacts={'a':{'artifact_id':'a','parsed_rows':2},'b':{'artifact_id':'b','parsed_rows':0},
                   'c':{'artifact_id':'c','parsed_rows':1}}
        queue=[{'study_key':'human_rna:GSE1','artifact_ids':['a']}]
        results=[{'study_key':'human_rna:GSE1','artifact_ids':['a'],'row_ids':['R1','R2']}]
        rows={'R1':{'identity_status':'PASS'},'R2':{'identity_status':'CANNOT VERIFY'},'R3':{'identity_status':'CANNOT VERIFY'}}
        receipt=coverage_report(artifacts,queue,results,rows)
        self.assertEqual(receipt['saved_artifacts'],3)
        self.assertEqual(receipt['parsed_logical_rows'],3)
        self.assertEqual(receipt['artifacts_without_candidate_study'],2)
        self.assertEqual(receipt['rows_without_candidate_study'],1)
        self.assertEqual(receipt['not_reviewed_rows'],1)
        self.assertEqual(receipt['reviewed_rows_with_bound_study_identity'],1)

    def test_omitted_reconstruction_is_not_a_material_error(self):
        design=validate_design(self.design(),self.index)
        rows=[{'row_id':'R1','fields':{'comparison_name':'original scientific target'},
               'artifact_sha256':'b'*64,'logical_record':2,
               'identity_binding':{'identity_status':'PASS','bound_study_key':'human_rna:GSE1'}}]
        payload=review_request(self.index,design,rows,[],[])
        value={'study_id':'GSE1','checks':{k:{**self.check(),'scope':'fixture'} for k in SCIENTIFIC_CHECKS},
               'rows':[{**self.check('FAIL'),'row_id':'R1','action':'withhold',
                        'failure_basis':'design_reconstruction_disagreement',
                        'source_comparison_ids':['paired'],'proposed_repair':'Analyze a different contrast'}],
               'design_disagreements':['The frozen source reconstruction omitted an original contrast']}
        result=validate_review(value,payload,rows)
        self.assertEqual(result['rows'][0]['status'],'CANNOT VERIFY')
        self.assertEqual(result['rows'][0]['model_status'],'FAIL')
        self.assertIsNone(result['rows'][0]['proposed_repair'])
        self.assertEqual(value['rows'][0]['status'],'FAIL')
        value['rows'][0]['failure_basis']='missing_evidence'
        self.assertEqual(validate_review(value,payload,rows)['rows'][0]['status'],'CANNOT VERIFY')

    def test_different_target_cannot_be_a_repair(self):
        design=validate_design(self.design(),self.index)
        rows=[{'row_id':'R1','fields':{'p_value':'1.1'},'artifact_sha256':'b'*64,'logical_record':2,
               'identity_binding':{'identity_status':'PASS','bound_study_key':'human_rna:GSE1'}}]
        payload=review_request(self.index,design,rows,[],[])
        value={'study_id':'GSE1','checks':{k:{**self.check(),'scope':'fixture'} for k in SCIENTIFIC_CHECKS},
               'rows':[{**self.check('FAIL'),'row_id':'R1','action':'propose_rerun',
                        'failure_basis':'invalid_numeric_or_record','source_comparison_ids':['paired'],
                        'repair_target_scope':'different_target','proposed_repair':'Separate target proposal'}],
               'design_disagreements':[]}
        result=validate_review(value,payload,rows)['rows'][0]
        self.assertEqual(result['status'],'FAIL')
        self.assertEqual(result['action'],'unresolved')
        self.assertEqual(result['source_comparison_ids'],[])
        self.assertEqual(result['separate_analysis_proposal']['text'],'Separate target proposal')
        self.assertIsNone(result['proposed_repair'])

    def timeout_task(self, name):
        """Local software rows and absent outcomes isolate transport behavior."""
        root = self.root / name
        root.mkdir()
        metadata = root / 'software-metadata.csv'
        metadata.write_text('GSM,title\nGSM1,Software fixture\n')
        saved = root / 'software-saved-rows.csv'
        saved.write_text('clock,comparison_name\nsoftware_clock,first\nsoftware_clock,second\n')
        artifact = {'artifact_id': 'software-artifact', 'path': str(saved.relative_to(root)),
                    'sha256': sha(saved)}
        task = {'study_key': 'human_rna:GSE1', 'study_id': 'GSE1', 'modality': 'human_rna',
                'source_files': {'source_metadata': {**declaration(metadata), 'available': True},
                    'companion_predictions': {'path': str(root / 'absent-software-outcomes.csv'),
                        'sha256': '0' * 64, 'available': False}}}
        config = {'validation_retries': 1, 'review_rows_per_call': 1,
            'api': {'model': 'gpt-5.6-luna', 'input_per_million': 0.2,
                'output_per_million': 1.2, 'output_tokens': 1000, 'prompt_bytes': 100000,
                'max_calls_per_study': 8, 'max_cost_per_study': 1.0, 'json_response': True}}
        budget = Budget(root / 'budget.sqlite', 1.0, 8)
        return root, task, [artifact], config, budget

    def fixture_transport(self, timeout_actor):
        """Mock curl replies; no empirical clock values or network are used."""
        def transport(command, **kwargs):
            request_path = Path(command[command.index('--data-binary') + 1][1:])
            actor = request_path.parent.name
            if actor == timeout_actor:
                raise subprocess.TimeoutExpired(['software-only-private-command'], 190,
                    output=b'software-only-private-output', stderr=b'software-only-private-stderr')
            request = json.loads(request_path.read_text())
            payload = json.loads(request['messages'][1]['content'])
            if actor.startswith('source_design'):
                value = self.design()
            else:
                value = {'study_id': 'GSE1',
                    'checks': {k: {**self.check(), 'scope': 'Software transport fixture'} for k in SCIENTIFIC_CHECKS},
                    'rows': [{**self.check(), 'row_id': row_id, 'action': 'unresolved',
                        'source_comparison_ids': [], 'proposed_repair': None}
                        for row_id, record in payload['source']['evidence'].items()
                        if record['artifact'] == 'saved_analysis'],
                    'design_disagreements': []}
            body = {'model': 'gpt-5.6-luna', 'usage': {'prompt_tokens': 10, 'completion_tokens': 3},
                'choices': [{'finish_reason': 'stop', 'message': {'content': json.dumps(value)}}]}
            return SimpleNamespace(returncode=0, stdout=json.dumps(body).encode() + b'\n200')
        return transport

    def assert_uncertain_terminal_resume(self, root, task, artifacts, config, budget, result, transport):
        output = root / 'study'
        calls = transport.call_count
        pending = budget.report()
        uncertain = [r for r in pending['by_status'] if r['status'] == 'RESERVED_OR_UNCERTAIN']
        self.assertEqual(len(uncertain), 1)
        self.assertEqual(uncertain[0]['calls'], 1)
        self.assertGreater(uncertain[0]['usd'], 0)
        before = {str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()}
        self.assertEqual((output / 'terminal.sha256').read_text().strip(), sha(output / 'terminal.json'))
        reconciled = budget.reconcile(root)
        self.assertEqual(reconciled['uncertain_reservations'], 1)
        self.assertAlmostEqual(budget.report()['estimated_or_reserved_usd'], pending['estimated_or_reserved_usd'])
        self.assertEqual(audit_study(task, artifacts, root, output, config, budget), result)
        self.assertEqual(transport.call_count, calls)
        self.assertEqual(before, {str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()})
        self.assertEqual(result['accepted_scientific_repairs'], 0)
        self.assertEqual(result['scientific_acceptance'], 'CANNOT VERIFY')
        diagnostic = next(output.rglob('transport-error.json'))
        self.assertIn(str(diagnostic.relative_to(output)), result['output_sha256'])
        self.assertFalse((diagnostic.parent / 'receipt.json').exists())
        self.assertFalse((diagnostic.parent / 'response.json').exists())
        for prohibited in ('software-only-private-command', 'software-only-private-output', 'software-only-private-stderr'):
            self.assertNotIn(prohibited, (output / 'terminal.json').read_text())
            self.assertNotIn(prohibited, diagnostic.read_text())

    def test_source_design_timeout_seals_unresolved_and_preserves_uncertain_resume(self):
        root, task, artifacts, config, budget = self.timeout_task('source-timeout')
        with patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}), \
                patch('audit_system.corpus_audit_2026_09_20.source_packet', return_value=self.packet), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                      side_effect=self.fixture_transport('source_design_00-000')) as transport:
            result = audit_study(task, artifacts, root, root / 'study', config, budget)
            self.assertEqual(result['status'], 'UNRESOLVED')
            self.assertIn('OPENAI_TRANSPORT_TIMEOUT', result['model_error'])
            self.assertEqual(result['reviewed_rows'], 0)
            self.assertEqual(transport.call_count, 1)
            self.assertFalse((root / 'study/api/source_design_01-000').exists())
            self.assert_uncertain_terminal_resume(root, task, artifacts, config, budget, result, transport)

    def test_review_timeout_retains_partial_judgments_and_unreviewed_rows(self):
        for chunk in (0, 1):
            with self.subTest(timeout_chunk=chunk):
                root, task, artifacts, config, budget = self.timeout_task('review-timeout-' + str(chunk))
                actor = f'review_{chunk:04d}_00-000'
                with patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}), \
                        patch('audit_system.corpus_audit_2026_09_20.source_packet', return_value=self.packet), \
                        patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                              side_effect=self.fixture_transport(actor)) as transport:
                    result = audit_study(task, artifacts, root, root / 'study', config, budget)
                    self.assertEqual(result['status'], 'PARTIALLY_REVIEWED')
                    self.assertEqual(result['reviewed_rows'], 1)
                    self.assertEqual(result['unreviewed_rows'], 1)
                    self.assertEqual(transport.call_count, 3)
                    failure = result['review_failures']
                    self.assertEqual(len(failure), 1)
                    self.assertEqual(failure[0]['stage'], f'04_review_{chunk:04d}')
                    self.assertIn('OPENAI_TRANSPORT_TIMEOUT', failure[0]['error'])
                    aggregation = completed_pipeline_stage(root / 'study', '05_aggregation')
                    self.assertEqual(aggregation['unreviewed_row_ids'], failure[0]['row_ids'])
                    self.assertEqual(len(aggregation['rows']), 1)
                    self.assertTrue(set(aggregation['unreviewed_row_ids']).isdisjoint(result['row_ids']))
                    self.assertFalse((root / f'study/api/review_{chunk:04d}_01-000').exists())
                    self.assert_uncertain_terminal_resume(root, task, artifacts, config, budget, result, transport)

    def test_transport_integrity_error_still_seals_integrity_failure(self):
        root, task, artifacts, config, budget = self.timeout_task('integrity')
        with patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}), \
                patch('audit_system.corpus_audit_2026_09_20.source_packet', return_value=self.packet), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                      side_effect=PipelineIntegrityError('software integrity failure')) as transport:
            result = audit_study(task, artifacts, root, root / 'study', config, budget)
            self.assertEqual(result['status'], 'INTEGRITY_FAILURE')
            self.assertIn('PipelineIntegrityError: software integrity failure', result['model_error'])
            self.assertEqual(transport.call_count, 1)
            self.assertEqual(list((root / 'study').rglob('transport-error.json')), [])


class FreshSourceFollowupTests(unittest.TestCase):
    """Offline software fixtures exercise lifecycle, never scientific accuracy."""
    module='audit_system.corpus_audit_2026_09_20'

    def setUp(self):
        from src.update_pipeline_2026_09_20 import seal,identity
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.version=self.root/'original';self.version.mkdir()
        self.pipeline=self.version/'pipeline';self.prior=self.version/'comprehensive-audit'
        self.allowance={'model':'gpt-5.6-luna','max_cost_usd':0.60,'max_calls':4,
                        'input_usd_per_million':0.5,'output_usd_per_million':1.8}
        self.allowance_path=self.root/'new-allowance.json';self.write(self.allowance_path,self.allowance)
        raw=self.root/'source.txt';raw.write_text('Software-only source declaration')
        old={**self.allowance,'max_calls':24,'max_cost_usd':5.0}
        value={'identity':{'software_fixture':'source-stop','comprehensive_audit':old},
               'study_id':'GSE1','modality':'human_rna',
               'sources':{'raw':declaration(raw)},'comprehensive_audit':old}
        value['version_id']=identity(value['identity']);seal(self.version/'version.json',value)
        store=StageStore(self.pipeline,{'software_fixture':True})
        def output(directory,name,text,field='path'):
            path=directory/name;path.write_text(text);return {field:str(path)}
        self.metadata=store.run('01_metadata',lambda d:output(d,'metadata.csv',
            'GSM,title,characteristics_ch1,p_value,prior_verdict\n'
            'GSM1,Software case,arm: case,PRIVATE_P_VALUE,PRIVATE_VERDICT\n'
            'GSM2,Software control,arm: control,PRIVATE_P_VALUE,PRIVATE_VERDICT\n','metadata'))['metadata']
        store.run('02_molecular_input',lambda d:output(d,'matrix.csv','gene,GSM1,GSM2\nfixture,1,2\n'))
        store.run('04_preprocessing',lambda d:output(d,'processed.csv','gene,GSM1,GSM2\nfixture,1,2\n'))
        store.run('05_clock_predictions',lambda d:output(d,'predictions.csv','GSM,software_clock\nGSM1,1\nGSM2,2\n'))
        store.run('06_coding_agent',lambda d:output(d,'saved.csv',
            'gse_id,clock,comparison_name,mean_difference,effect_size\n'
            'GSE1,software_clock,PRIVATE_PREVIOUS_TARGET,,\n'))
        store.run('08_report',lambda d:output(d,'report.html','Software fixture report','report'))
        self.request={'study_id':'GSE1','modality':'human_rna','version_id':value['version_id'],
            'pipeline_root':str(self.pipeline),'source_bindings':value['sources'],'output':str(self.prior),
            'budget':old,'scientific_status':'CANNOT VERIFY','original_paper_2025_eligibility':'QUARANTINED'}
        self.write(self.version/'comprehensive-audit-request.json',{'audit_request':self.request})
        self.index=evidence_index(source_packet(self.metadata,'GSE1',TARGET,comparison_scope='source_defined_interventions'))
        self.seen=[]
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=self.transport('reject')) as intercepted:
            result=audit_pipeline_result(self.request)
        self.assertEqual(intercepted.call_count,2);self.assertEqual(result['status'],'UNRESOLVED')
        seal(self.version/'pipeline-report.json',{'status':'COMPLETED','until':'report',
             'output_sha256':self.files(self.pipeline),'result':{'software_fixture':True}})
        seal(self.version/'comprehensive-audit.json',{'status':'UNRESOLVED','result':result,
             'output_sha256':self.files(self.prior)})
        self.output=self.root/'followup';self.seen=[]

    def tearDown(self):self.tmp.cleanup()

    @staticmethod
    def write(path,value):
        Path(path).parent.mkdir(parents=True,exist_ok=True)
        Path(path).write_text(json.dumps(value,indent=2)+'\n')

    @staticmethod
    def files(root):return {str(p.resolve()):sha(p) for p in Path(root).rglob('*') if p.is_file()}

    def reseal(self,path,value):
        self.write(path,value);path.with_suffix(path.suffix+'.sha256').write_text(sha(path)+'\n')

    def reseal_prior(self):
        """Counterfactual software mutations are resealed to test semantic guards."""
        study=self.prior/'study';terminal=study/'terminal.json';value=json.loads(terminal.read_text())
        value['output_sha256']={str(p.relative_to(study)):sha(p) for p in study.rglob('*')
            if p.is_file() and p not in (terminal,study/'terminal.sha256')}
        self.write(terminal,value);(study/'terminal.sha256').write_text(sha(terminal)+'\n')
        wrapper=self.version/'comprehensive-audit.json';record=json.loads(wrapper.read_text())
        record['result'].update(value);record['output_sha256']=self.files(self.prior)
        self.reseal(wrapper,record)

    def design(self,reject=False):
        eid=next(k for k in self.index['evidence'] if k.startswith('S'))
        check={'status':'CANNOT VERIFY','reason':'Software fixture: no scientific reference',
               'evidence_ids':[eid],'unresolved':['No biological-unit claim in software fixture']}
        result={'study_id':'GSE1','checks':{k:copy.deepcopy(check) for k in DOMAINS},
                'comparisons':[],'unresolved_targets':['Software fixture lacks a supported contrast']}
        if reject:
            result['comparisons']=[{'id':'software','target':'software overlap',
                'case_all':[eid],'control_all':[eid],'mode':'unresolved','method':'unsupported',
                'unit_field':None,'pair_field':None,'independence':copy.deepcopy(check),'covariates':[],
                'evidence_ids':[eid],'limitations':['Software fixture']}]
        return result

    def transport(self,mode,*,after_first=None):
        def intercept(command,**kwargs):
            self.assertIn('--data-binary',command,'Only intercepted model requests are expected')
            path=Path(command[command.index('--data-binary')+1][1:]);request=json.loads(path.read_text())
            actor=path.parent.name;payload=json.loads(request['messages'][1]['content'])
            self.seen.append((actor,request))
            if mode=='timeout':raise subprocess.TimeoutExpired(['offline-fixture'],190)
            if mode=='integrity':raise PipelineIntegrityError('Software integrity failure')
            if actor.startswith('source_design'):
                value=self.design(reject=mode=='reject' or (mode=='retry' and actor=='source_design_00-000'))
            else:
                check={'status':'CANNOT VERIFY','reason':'Software fixture only','evidence_ids':[],
                       'unresolved':['No scientific reference in software fixture']}
                value={'study_id':'GSE1','checks':{k:{**check,'scope':'Software fixture'} for k in SCIENTIFIC_CHECKS},
                       'rows':[{**check,'row_id':rid,'action':'unresolved','source_comparison_ids':[],
                                'proposed_repair':None} for rid in payload['review_row_ids']],
                       'design_disagreements':[]}
            body={'model':'gpt-5.6-luna','usage':{'prompt_tokens':10,'completion_tokens':3},
                  'choices':[{'finish_reason':'content_filter' if mode=='policy' else 'stop',
                              'message':{'content':json.dumps(value)}}]}
            if after_first is not None and len(self.seen)==1:after_first()
            return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')
        return intercept

    def prepare(self):
        return prepare_pipeline_followup(self.version,self.prior,self.allowance_path,self.output,
                                         'software-fresh',fresh_source_design=True)

    def run_offline(self,mode='valid',**kwargs):
        measurement={'study_id':'GSE1','evidence':{'M00000':{'artifact':'software_fixture',
            'value':'No empirical measurement applicability is asserted','sample_ids':[]}},
            'source_files':{},'stage_receipts':{}}
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=self.transport(mode,**kwargs)) as transport, \
                patch('audit_system.pipeline_measurement_evidence_2026_09_20.build_measurement_evidence',return_value=measurement):
            result=run_pipeline_followup(self.output)
        return result,transport.call_count

    def test_normal_adapter_gpt6_budget_is_bound_before_model_transport(self):
        from src.update_pipeline_2026_09_20 import identity
        request=copy.deepcopy(self.request)
        request['budget']={**self.allowance,'model':'gpt-6-luna',
            'input_usd_per_million':.20,'output_usd_per_million':.75}
        request['output']=str(self.root/'gpt6-normal-audit')
        with patch(self.module+'.StudyClient.ask',side_effect=AssertionError('Unbound request must stop')):
            with self.assertRaisesRegex(PipelineIntegrityError,'sealed update version'):
                audit_pipeline_result(request)
        # Counterfactual software-only version binds the new budget explicitly.
        version_path=self.version/'version.json';version=json.loads(version_path.read_text())
        version['comprehensive_audit']=request['budget']
        self.reseal(version_path,version)
        with patch(self.module+'.StudyClient.ask',side_effect=AssertionError('Identity mismatch must stop')):
            with self.assertRaisesRegex(PipelineIntegrityError,'sealed update version'):
                audit_pipeline_result(request)
        version['identity']['comprehensive_audit']=request['budget']
        version['version_id']=identity(version['identity'])
        request['version_id']=version['version_id'];self.reseal(version_path,version)
        seen=[]
        def inspect(task,artifacts,root,output,config,budget):
            seen.append(copy.deepcopy(config))
            self.assertEqual(config['api']['model'],'gpt-6-luna')
            self.assertEqual(config['api']['input_per_million'],.20)
            self.assertEqual(config['api']['output_per_million'],.75)
            self.assertEqual(config['max_total_calls'],request['budget']['max_calls'])
            self.assertEqual(config['max_total_cost_usd'],request['budget']['max_cost_usd'])
            return {'status':'UNRESOLVED','scientific_acceptance':'CANNOT VERIFY',
                    'accepted_scientific_repairs':0}
        with patch(self.module+'.audit_study',side_effect=inspect):
            result=audit_pipeline_result(request)
        self.assertEqual(len(seen),1)
        self.assertEqual(result['scientific_acceptance'],'CANNOT VERIFY')
        self.assertEqual(result['accepted_scientific_repairs'],0)
        for changed_field,value in [('model','gpt-5.6-luna'),('max_calls',100),('input_usd_per_million',.01)]:
            changed=copy.deepcopy(request);changed['budget'][changed_field]=value
            with patch(self.module+'.audit_study',side_effect=AssertionError('Budget tampering must stop')):
                with self.assertRaisesRegex(PipelineIntegrityError,'sealed update version'):
                    audit_pipeline_result(changed)

    def test_normal_adapter_binds_and_forwards_explicit_service_tier(self):
        from src.update_pipeline_2026_09_20 import identity
        request=copy.deepcopy(self.request)
        request['budget']={**self.allowance,'model':'gpt-6-luna','service_tier':'default'}
        request['output']=str(self.root/'tier-normal-audit')
        version_path=self.version/'version.json';version=json.loads(version_path.read_text())
        for record in (version,version['identity']):record['comprehensive_audit']=request['budget']
        version['version_id']=identity(version['identity']);request['version_id']=version['version_id']
        self.reseal(version_path,version)
        def inspect(task,artifacts,root,output,config,budget):
            self.assertEqual(config['api']['service_tier'],'default')
            self.assertEqual(config['api']['model'],'gpt-6-luna')
            self.assertNotIn('service_tier',config)
            return {'status':'UNRESOLVED','scientific_acceptance':'CANNOT VERIFY'}
        with patch(self.module+'.audit_study',side_effect=inspect) as study:
            audit_pipeline_result(request)
            self.assertEqual(study.call_count,1)
        for value in (None,False,[],{},'auto','priority'):
            request['budget']['service_tier']=value
            for record in (version,version['identity']):record['comprehensive_audit']=request['budget']
            version['version_id']=identity(version['identity']);request['version_id']=version['version_id']
            self.reseal(version_path,version)
            with self.subTest(value=value), patch(self.module+'.Budget',side_effect=AssertionError('No budget initialization')):
                with self.assertRaisesRegex(ValueError,'service_tier'):
                    audit_pipeline_result(request)

    def test_fresh_followup_optional_tier_preserves_old_requests_and_usage(self):
        before=self.files(self.version)
        self.allowance.update(model='gpt-6-luna',service_tier='default')
        self.write(self.allowance_path,self.allowance)
        with patch(self.module+'.StudyClient.ask',side_effect=AssertionError('No model call in preparation')):
            self.prepare()
        frozen=verify_pipeline_followup(self.output/'followup.json')
        self.assertEqual(frozen['request']['budget']['service_tier'],'default')
        self.assertEqual(frozen['request']['budget']['model'],'gpt-6-luna')
        self.assertEqual(before,self.files(self.version))
        self.assertEqual(frozen['fresh_source']['prior_budget']['call_reservations'],2)
        self.assertTrue(all('service_tier' not in json.loads(path.read_text())
                            for path in (self.prior/'study/api').glob('*/request.json')))

    def test_opt_in_source_packet_and_complete_prior_inputs_are_frozen_without_calls(self):
        before=self.files(self.version)
        with patch(self.module+'.StudyClient.ask',side_effect=AssertionError('No calls during preparation')):
            result=self.prepare();frozen=verify_pipeline_followup(self.output/'followup.json')
        self.assertEqual(result['model_calls'],0);self.assertEqual(before,self.files(self.version))
        self.assertNotIn('source_design_replay',frozen)
        self.assertTrue(set(before)<=set(frozen['input_files']))
        self.assertIn(str(self.prior/'budget.sqlite3'),frozen['input_files'])
        self.assertEqual(frozen['request'],{**self.request,'output':str(self.output/'audit'),'budget':self.allowance})
        payload=design_request(self.index)
        self.assertEqual(frozen['fresh_source']['initial_request_sha256'],record_hash({'system':DESIGN_SYSTEM,'payload':payload}))
        for withheld in ('PRIVATE_P_VALUE','PRIVATE_VERDICT','PRIVATE_PREVIOUS_TARGET','previous_response'):
            self.assertNotIn(withheld,json.dumps(payload))

    def reviewed_prior_for_allowance_tests(self):
        """Create a separate reviewed software fixture; no empirical output is used."""
        self.prior=self.root/'reviewed-software-audit'
        self.request={**self.request,'output':str(self.prior)}
        self.write(self.version/'comprehensive-audit-request.json',{'audit_request':self.request})
        measurement={'study_id':'GSE1','evidence':{'M00000':{'artifact':'software_fixture',
            'value':'No empirical applicability is asserted','sample_ids':[]}},
            'source_files':{},'stage_receipts':{}}
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=self.transport('valid')), \
                patch('audit_system.pipeline_measurement_evidence_2026_09_20.build_measurement_evidence',return_value=measurement):
            result=audit_pipeline_result(self.request)
        self.assertEqual(result['status'],'REVIEWED')
        self.reseal(self.version/'comprehensive-audit.json',{'status':'COMPLETED','result':result,
                    'output_sha256':self.files(self.prior)})

    def assert_allowance_prepares_without_execution(self,fresh):
        before=self.files(self.version);prior_before=self.files(self.prior)
        with patch(self.module+'.StudyClient.ask',side_effect=AssertionError('No model in prepare')), \
                patch(self.module+'.Budget',side_effect=AssertionError('No new budget in prepare')), \
                patch(self.module+'.audit_pipeline_result',side_effect=AssertionError('No audit in prepare')), \
                patch('src.update_pipeline_2026_09_20.run_child',side_effect=AssertionError('No producer/clock child')):
            result=prepare_pipeline_followup(self.version,self.prior,self.allowance_path,self.output,
                                            'software-model-allowance',fresh_source_design=fresh)
            frozen=verify_pipeline_followup(self.output/'followup.json')
        self.assertEqual(result['model_calls'],0)
        self.assertEqual(self.files(self.version),before)
        self.assertEqual(self.files(self.prior),prior_before)
        self.assertFalse((self.output/'audit').exists())
        self.assertIs(frozen['new_producer_or_clock_execution'],False)
        self.assertEqual(frozen['request'],{**self.request,'output':str(self.output/'audit'),
                                          'budget':json.loads(self.allowance_path.read_text())})
        self.assertEqual(frozen['implementation_sha256']['audit_system/corpus_audit_2026_09_20.py'],
                         sha(Path(__import__(self.module,fromlist=['__file__']).__file__)))
        self.assertEqual(frozen['schema'],'clockbase-pipeline-audit-followup-v2' if fresh
                         else 'clockbase-pipeline-audit-followup-v1')
        self.assertEqual('source_design_replay' in frozen,not fresh)
        if fresh:
            payload=design_request(self.index)
            self.assertEqual(frozen['fresh_source']['initial_request_sha256'],
                             record_hash({'system':DESIGN_SYSTEM,'payload':payload}))
            for hidden in ('PRIVATE_P_VALUE','PRIVATE_VERDICT','PRIVATE_PREVIOUS_TARGET','previous_response'):
                self.assertNotIn(hidden,json.dumps(payload))
        else:
            self.assertEqual(frozen['source_design_replay']['root'],str(self.prior/'study'))

    def test_explicit_luna_models_prepare_fresh_without_calls_or_prior_changes(self):
        for model,input_rate,output_rate in [('gpt-5.6-luna',.5,1.8),('gpt-6-luna',.25,.75)]:
            with self.subTest(model=model):
                self.output=self.root/('fresh-'+model)
                self.write(self.allowance_path,{**self.allowance,'model':model,
                    'input_usd_per_million':input_rate,'output_usd_per_million':output_rate})
                self.assert_allowance_prepares_without_execution(True)

    def test_explicit_luna_models_prepare_replay_without_calls_or_prior_changes(self):
        self.reviewed_prior_for_allowance_tests()
        for model,input_rate,output_rate in [('gpt-5.6-luna',.5,1.8),('gpt-6-luna',.25,.75)]:
            with self.subTest(model=model):
                self.output=self.root/('replay-'+model)
                self.write(self.allowance_path,{**self.allowance,'model':model,
                    'input_usd_per_million':input_rate,'output_usd_per_million':output_rate})
                self.assert_allowance_prepares_without_execution(False)

    def invalid_followup_allowances(self):
        good={**self.allowance,'model':'gpt-6-luna','input_usd_per_million':.25,'output_usd_per_million':.75}
        cases=[]
        for key in good:
            missing=copy.deepcopy(good);missing.pop(key);cases.append(('missing-'+key,missing))
        for value in ('gpt-6-astra','gpt-6-luna-alias','',None,True,['gpt-6-luna']):
            cases.append(('model-'+repr(value),{**good,'model':value}))
        for key in ('max_cost_usd','input_usd_per_million','output_usd_per_million'):
            for value in (0,-1,True,float('nan'),float('inf')):
                cases.append((key+'-'+repr(value),{**good,key:value}))
        for value in (0,-1,True,4.0):
            cases.append(('calls-'+repr(value),{**good,'max_calls':value}))
        cases.append(('extra-field',{**good,'model_fallback':'gpt-5.6-luna'}))
        return cases

    def assert_invalid_allowances_stop_before_preparation(self,fresh):
        before=self.files(self.version);prior_before=self.files(self.prior)
        for i,(name,value) in enumerate(self.invalid_followup_allowances()):
            with self.subTest(fresh=fresh,invalid=name):
                self.output=self.root/('invalid-'+str(fresh)+'-'+str(i))
                self.write(self.allowance_path,value)
                with patch(self.module+'.StudyClient.ask',side_effect=AssertionError('No model in prepare')), \
                        patch(self.module+'.Budget',side_effect=AssertionError('No new budget in prepare')):
                    with self.assertRaises((ValueError,PipelineIntegrityError)):
                        prepare_pipeline_followup(self.version,self.prior,self.allowance_path,self.output,
                                                  'software-invalid',fresh_source_design=fresh)
                self.assertFalse(self.output.exists())
        self.assertEqual(self.files(self.version),before)
        self.assertEqual(self.files(self.prior),prior_before)

    def test_fresh_allowance_rejects_missing_unknown_models_and_invalid_tariffs_caps(self):
        self.assert_invalid_allowances_stop_before_preparation(True)

    def test_replay_allowance_rejects_missing_unknown_models_and_invalid_tariffs_caps(self):
        self.reviewed_prior_for_allowance_tests()
        self.assert_invalid_allowances_stop_before_preparation(False)

    def test_gpt6_allowance_does_not_relax_original_gpt56_response_binding(self):
        self.write(self.allowance_path,{**self.allowance,'model':'gpt-6-luna',
                   'input_usd_per_million':.25,'output_usd_per_million':.75})
        slot=self.prior/'study/api/source_design_00-000'
        response_path=slot/'response.json';response=json.loads(response_path.read_text())
        response['model']='gpt-6-luna';self.write(response_path,response)
        receipt_path=slot/'receipt.json';receipt=json.loads(receipt_path.read_text())
        receipt['returned_model']='gpt-6-luna';receipt['response_sha256']=sha(response_path)
        self.write(receipt_path,receipt);self.reseal_prior()
        before=self.files(self.version)
        with patch(self.module+'.StudyClient.ask',side_effect=AssertionError('No model in prepare')):
            with self.assertRaisesRegex(PipelineIntegrityError,'completed matching model response'):
                self.prepare()
        self.assertFalse(self.output.exists());self.assertEqual(self.files(self.version),before)

    def test_legacy_default_still_requires_completed_prior_review(self):
        with self.assertRaisesRegex(ValueError,'requires a completed prior review'):
            prepare_pipeline_followup(self.version,self.prior,self.allowance_path,self.output,'legacy')
        self.assertFalse(self.output.exists())

    def test_resealed_mode_request_index_and_inventory_mutations_fail_before_budget(self):
        self.prepare();path=self.output/'followup.json';original=json.loads(path.read_text())
        mutations=[lambda x:x.update(source_design_replay={'root':str(self.prior/'study')}),
                   lambda x:x.update(schema='clockbase-pipeline-audit-followup-v1'),
                   lambda x:x.update(mode='replay'),
                   lambda x:x['request'].update(scientific_status='PASS'),
                   lambda x:x['fresh_source'].update(source_index_sha256='0'*64),
                   lambda x:x['input_files'].pop(str(self.prior/'budget.sqlite3'))]
        for mutate in mutations:
            with self.subTest(mutation=mutate):
                value=copy.deepcopy(original);mutate(value);self.reseal(path,value)
                with patch(self.module+'.Budget',side_effect=AssertionError('Budget must not be opened')):
                    with self.assertRaises(PipelineIntegrityError):run_pipeline_followup(self.output)

    def test_changed_prior_stage_and_missing_report_are_ineligible(self):
        path=self.pipeline/'08_report/complete.sha256';path.unlink()
        with self.assertRaises(PipelineIntegrityError):self.prepare()
        self.assertFalse(self.output.exists())

    def test_resealed_contract_followed_by_transport_failure_is_ineligible(self):
        failure=self.prior/'study/02_blinded_design/attempt-000/failure.json'
        value=json.loads(failure.read_text());value['error']='RuntimeError: OPENAI_TRANSPORT_TIMEOUT';self.write(failure,value)
        terminal=self.prior/'study/terminal.json';value=json.loads(terminal.read_text())
        value['model_error']='RuntimeError: OPENAI_TRANSPORT_TIMEOUT';self.write(terminal,value);self.reseal_prior()
        with self.assertRaisesRegex(PipelineIntegrityError,'final source failure'):self.prepare()

    def test_resealed_request_only_successor_is_ineligible(self):
        self.write(self.prior/'study/api/source_design_02-000/request.json',{'software_fixture':True})
        self.reseal_prior()
        with self.assertRaisesRegex(PipelineIntegrityError,'uncompleted or unrelated API'):self.prepare()

    def test_resealed_reservation_only_successor_is_ineligible(self):
        with sqlite3.connect(self.prior/'budget.sqlite3') as db:
            db.execute('INSERT INTO reservations VALUES (?,?,?,?,?,?)',
                       ('extra',0.1,0.1,'RESERVED_OR_UNCERTAIN','software-fixture','0'*64))
        self.reseal_prior()
        with self.assertRaisesRegex(PipelineIntegrityError,'unmatched calls'):self.prepare()

    def test_runtime_rechecks_prior_eligibility_after_recomputed_outer_seals(self):
        self.prepare();path=self.output/'followup.json';frozen=json.loads(path.read_text())
        terminal=self.prior/'study/terminal.json';value=json.loads(terminal.read_text());value['status']='REVIEWED'
        self.write(terminal,value);self.reseal_prior()
        for name in frozen['input_files']:frozen['input_files'][name]=declaration(name)
        self.reseal(path,frozen)
        with patch(self.module+'.Budget',side_effect=AssertionError('No budget after ineligible terminal')):
            with self.assertRaisesRegex(PipelineIntegrityError,'unresolved wrapper'):run_pipeline_followup(self.output)

    def test_successful_fresh_run_preserves_sources_and_completed_resume_has_zero_calls(self):
        self.prepare();before=self.files(self.version);result,calls=self.run_offline()
        self.assertEqual((result['status'],calls),('REVIEWED',2));self.assertEqual(result['accepted_scientific_repairs'],0)
        self.assertEqual(result['scientific_acceptance'],'CANNOT VERIFY')
        task=json.loads((self.output/'audit/study/config.json').read_text())['task'];self.assertNotIn('source_design_replay',task)
        initial=json.loads(self.seen[0][1]['messages'][1]['content']);self.assertEqual(initial,design_request(self.index))
        self.assertNotIn('previous_response',initial);saved=self.files(self.output/'audit/study')
        resumed,calls=self.run_offline('timeout');self.assertEqual(calls,0);self.assertEqual(resumed,result)
        self.assertEqual(saved,self.files(self.output/'audit/study'));self.assertEqual(before,self.files(self.version))

    def test_completed_contract_retry_and_unresolved_resume_preserve_each_response(self):
        self.prepare();result,calls=self.run_offline('reject')
        self.assertEqual((result['status'],calls),('UNRESOLVED',2));self.assertEqual(result['budget']['call_reservations'],2)
        second=json.loads(self.seen[1][1]['messages'][1]['content'])
        self.assertIn('source_selector_diagnostics',second['contract_validation_feedback'])
        self.assertEqual(second['source'],self.index)
        saved=self.files(self.output/'audit/study');resumed,calls=self.run_offline('valid')
        self.assertEqual(calls,0);self.assertEqual(resumed,result);self.assertEqual(saved,self.files(self.output/'audit/study'))

    def test_one_source_contract_retry_can_continue_to_the_same_reviewer(self):
        self.prepare();result,calls=self.run_offline('retry')
        self.assertEqual((result['status'],calls),('REVIEWED',3));self.assertEqual(result['budget']['max_calls'],4)
        self.assertEqual(result['budget']['max_cost_usd'],0.6)

    def test_transport_policy_integrity_and_budget_failures_are_not_contract_retries(self):
        for mode in ('timeout','policy','integrity','budget'):
            with self.subTest(mode=mode):
                self.output=self.root/('followup-'+mode)
                if mode=='budget':
                    self.write(self.allowance_path,{**self.allowance,'max_cost_usd':0.000001})
                self.prepare();self.seen=[];result,calls=self.run_offline(mode)
                self.assertEqual(calls,0 if mode=='budget' else 1)
                self.assertIn(result['status'],('UNRESOLVED','INTEGRITY_FAILURE','BUDGET_DEFERRED'))
                self.assertFalse((self.output/'audit/study/api/source_design_01-000').exists())

    def test_source_drift_after_first_response_stops_before_retry_reservation(self):
        self.prepare();self.seen=[]
        with self.assertRaises(PipelineIntegrityError):
            self.run_offline('retry',after_first=lambda:Path(self.metadata).write_text('changed software source'))
        self.assertEqual(len(self.seen),1)
        self.assertFalse((self.output/'audit/study/api/source_design_01-000').exists())
        terminal=json.loads((self.output/'audit/study/terminal.json').read_text())
        self.assertEqual(terminal['status'],'INTEGRITY_FAILURE')


class SelectedStageRecoveryTests(unittest.TestCase):
    """Fabricated software packets and intercepted transport, never empirical results."""
    module='audit_system.corpus_audit_2026_09_20'

    def write(self,path,value):
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(value,indent=2)+'\n')

    def gz(self,path,values):
        path.parent.mkdir(parents=True,exist_ok=True)
        with gzip.open(path,'wt') as handle:
            for value in values:handle.write(json.dumps(value)+'\n')

    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.old=self.root/'historical';self.repo=self.root/'repo';self.output=self.root/'recovery'
        self.feas=self.root/'feasibility';self.diagnosis=self.root/'diagnosis';self.seen=[]
        self.config={'provider':'openai','concurrency':1,'ordering_seed':'software-fixture',
            'review_rows_per_call':24,'validation_retries':1,'max_total_calls':100,'max_total_cost_usd':10,
            'api':{'model':'gpt-5.6-luna','input_per_million':.5,'output_per_million':1.8,
                   'output_tokens':16000,'prompt_bytes':900000,'max_calls_per_study':64,
                   'max_cost_per_study':10,'json_response':True},'freeze_inputs':[]}
        self.write(self.old/'config.json',{'config':self.config})
        self.stage_list=[];self.requests=[];self.tasks=[];self.queue=[];self.artifacts=[];self.policy=[];self.targets={}
        for number,stage in getattr(self,'selected_cases',[(1,'02_blinded_design'),(2,'04_review_0000')]):
            study=f'GSE{number}';key='human_rna:'+study;rid=f'R{number:05d}';event='event:software'+str(number)
            metadata=self.repo/(study+'-metadata.json');predictions=self.repo/(study+'-predictions.csv')
            self.write(metadata,{'software_fixture':True,'study_id':study});predictions.write_text('software_fixture\n1\n')
            artifact=self.repo/(study+'-saved.csv');artifact.write_text('clock,p_value\nsoftware,NA\n')
            artifact_record={'artifact_id':'artifact'+str(number),'path':artifact.name,'sha256':sha(artifact),'parsed_rows':1}
            self.artifacts.append(artifact_record)
            source={'study_id':study,'source_url':'https://example.invalid/software-fixture','source_sha256':sha(metadata),
                    'samples':[{'sample_id':'GSM'+str(i),'condition':'software','unit':'unknown'} for i in range(1,4)]}
            index=evidence_index(source);ref=next(k for k in index['evidence'] if k.startswith('S'))
            check={'status':'CANNOT VERIFY','reason':'Fabricated software test packet','evidence_ids':[ref],
                   'unresolved':['No empirical scientific claim in software fixture']}
            design={'study_id':study,'checks':{k:copy.deepcopy(check) for k in DOMAINS},
                    'comparisons':[],'unresolved_targets':['Software fixture']}
            row={'row_id':rid,'fields':{'clock':'software','p_value':'NA'},'artifact_sha256':sha(artifact),
                 'logical_record':2,'identity_binding':{'identity_status':'PASS','bound_study_key':key}}
            numerical={'results':[{'clock':'software','status':'FAIL','reason':'Preserved software numerical failure'}],
                       'request_sha256':'b'*64,'predictions_sha256':sha(predictions),'source_hash_status':'PASS'}
            calculation={'comparison_id':'software','clock':'software','case_ids':['GSM1'],'control_ids':['GSM2'],
                         'status':'CANNOT VERIFY','reason':'Preserved software numerical stop'}
            payload=design_request(index) if stage=='02_blinded_design' else review_request(index,design,[row],[calculation],[],
                whole_study_inventory=[{'row_id':rid,'artifact_id':artifact_record['artifact_id'],'target':{'clock':'software'}}],
                independent_numerics=numerical)
            system=DESIGN_SYSTEM if stage=='02_blinded_design' else REVIEW_SYSTEM
            request={'model':'gpt-5.6-luna','messages':[{'role':'system','content':system},
                     {'role':'user','content':json.dumps(payload,ensure_ascii=False,separators=(',',':'))}],
                     'max_completion_tokens':16000,'reasoning_effort':'low','store':False,
                     'response_format':{'type':'json_object'}}
            task={'study_key':key,'study_id':study,'modality':'human_rna','artifact_ids':[artifact_record['artifact_id']],
                  'source_files':{'source_metadata':{**declaration(metadata),'available':True},
                                  'companion_predictions':{**declaration(predictions),'available':True}}}
            self.queue.append(task);child=self.old/'studies'/key.replace(':','-')
            actor='source_design' if stage=='02_blinded_design' else 'review_0000';request_path=child/'api'/(actor+'_00-000')/'request.json'
            self.write(request_path,request)
            self.write(child/'01_source_index/attempt-000/source-index.json',index)
            if stage!='02_blinded_design':
                self.write(child/'02_blinded_design/attempt-000/source-design.json',design)
                self.write(child/'03_conditional_reruns/attempt-000/calculations.json',[calculation])
                self.write(child/'03b_independent_r/attempt-000/numerical.json',numerical)
            error='ValueError: Source groups are empty or overlap' if stage=='02_blinded_design' else 'ValueError: Review omitted, duplicated or invented a saved row'
            self.write(child/stage/'attempt-000/failure.json',{'error':error})
            terminal={'study_key':key,'status':'UNRESOLVED' if stage=='02_blinded_design' else 'PARTIALLY_REVIEWED',
                      'model_error':error if stage=='02_blinded_design' else None,'review_failures':[] if stage=='02_blinded_design' else
                      [{'stage':stage,'row_ids':[rid],'error':error}],
                      'output_sha256':{str(p.relative_to(child)):sha(p) for p in child.rglob('*') if p.is_file()}}
            self.write(child/'terminal.json',terminal)
            self.stage_list.append({'study_key':key,'stage':stage,'event_id':event})
            self.requests.append({'study_key':key,'stage':stage,'event_id':event,'initial_request_planning_eligible':True,
                'planning_hold_reasons':[],'original_block_row_ids':[rid],'original_request_path':str(request_path),
                'original_request_sha256':sha(request_path),
                'current_accounting_request_sha256':hashlib.sha256(json.dumps(request,ensure_ascii=False).encode()).hexdigest()})
            self.tasks.append({'study_key':key,'original_queue_record':task,'original_terminal_sha256':sha(child/'terminal.json'),
                               'planning_hold_reasons':[]})
            self.policy.append({'short_row_id':rid,'attached_task_keys':[key],
                'next_action':'response_contract_development_candidate','provider_policy_hold':False,
                'completed_recovery_overlays':[],'failure_contexts':[{'study_key':key,'stage':stage,'event_id':event}]})
            self.targets[event]=[rid]
        self.gz(self.old/'queue.jsonl.gz',self.queue);self.gz(self.root/'artifact-worklist.jsonl.gz',self.artifacts)
        self.write(self.root/'summary.json',{'outputs':{'artifact-worklist.jsonl.gz':sha(self.root/'artifact-worklist.jsonl.gz')}})
        self.write(self.old/'config.json',{'config':self.config,'worklist':str(self.root),
                                        'worklist_summary_sha256':sha(self.root/'summary.json')})
        self.write(self.root/'old-cost.json',{'conservative_usd':100.0})
        self.pool={'schema':'clockbase-selected-stage-cost-pool-v1','allocation_status':'FROZEN','pool_cap_usd':4000,
                   'max_new_cost_usd':.6,'max_new_calls':4,'prior_conservative_usd':100.0,
                   'charges':[{'run_id':'software-prior','receipt':declaration(self.root/'old-cost.json'),
                      'amount_pointer':'/conservative_usd','conservative_usd':100.0,'evidence_kind':'USAGE_PLUS_UNCERTAIN'}]}
        self.plan={'schema':'clockbase-selected-stage-recovery-plan-v1','phase':'selected_failed_stages_only',
                   'ordering_seed':'software-deterministic-order','expected_stages':2,'expected_tasks':2,
                   'historical_root':str(self.old),'historical_config':declaration(self.old/'config.json'),
                   'historical_queue_sha256':sha(self.old/'queue.jsonl.gz'),'repo':str(self.repo),
                   'artifact_worklist':declaration(self.root/'artifact-worklist.jsonl.gz')}
        self.sync_inputs()

    def tearDown(self):self.tmp.cleanup()

    def sync_inputs(self):
        self.write(self.feas/'planning-stage-list.json',self.stage_list)
        self.gz(self.feas/'request-size-roster.jsonl.gz',self.requests)
        self.gz(self.feas/'task-source-manifest.jsonl.gz',self.tasks)
        self.write(self.feas/'receipt.json',{'status':'PASS','outputs_sha256':{
            name:sha(self.feas/name) for name in ('planning-stage-list.json','request-size-roster.jsonl.gz','task-source-manifest.jsonl.gz')}})
        self.gz(self.diagnosis/'unreviewed-row-roster.jsonl.gz',self.policy)
        self.write(self.diagnosis/'receipt.json',{'status':'PASS','outputs_sha256':{
            'unreviewed-row-roster.jsonl.gz':sha(self.diagnosis/'unreviewed-row-roster.jsonl.gz')}})
        self.write(self.root/'row-targets.json',{'status':'PASS','per_stage_target_short_row_ids':self.targets})
        self.write(self.root/'pool.json',self.pool)
        self.plan.update(feasibility_receipt=declaration(self.feas/'receipt.json'),
                         diagnosis_receipt=declaration(self.diagnosis/'receipt.json'),
                         row_scope=declaration(self.root/'row-targets.json'),cost_pool=declaration(self.root/'pool.json'))
        self.write(self.root/'plan.json',self.plan)

    def prepare(self):return prepare_selected_stage_recovery(self.root/'plan.json',self.output)

    def transport(self,mode='valid',after_first=None):
        def call(command,**kwargs):
            self.assertIn('--data-binary',command,'No numerical/producer subprocess is allowed')
            path=Path(command[command.index('--data-binary')+1][1:]);request=json.loads(path.read_text())
            payload=json.loads(request['messages'][1]['content']);self.seen.append(request)
            if mode=='interrupt':raise KeyboardInterrupt('Software hard interruption')
            if mode=='timeout':raise subprocess.TimeoutExpired(['software-transport'],190)
            if mode=='policy':
                body={'error':{'code':'bio_policy','type':'policy','message':'Software fixture refusal'}}
                return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n400')
            ref=next(k for k in payload['source']['evidence'] if k.startswith('S'))
            check={'status':'CANNOT VERIFY','reason':'Fabricated software test response','evidence_ids':[ref],
                   'unresolved':['No biological reference in software fixture']}
            if 'review_row_ids' in payload:
                value={'study_id':payload['study_id'],'checks':{k:{**check,'scope':'Software fixture'} for k in SCIENTIFIC_CHECKS},
                       'rows':[{**check,'row_id':rid,'action':'unresolved','source_comparison_ids':[],
                                'proposed_repair':None} for rid in payload['review_row_ids']],'design_disagreements':[]}
            else:value={'study_id':payload['source']['study_id'],'checks':{k:copy.deepcopy(check) for k in DOMAINS},
                        'comparisons':[],'unresolved_targets':['Software fixture']}
            if mode=='reject' or (mode in ('retry','retry_incomplete') and path.parent.name.endswith('_00-000')):value={}
            if after_first is not None and len(self.seen)==1:after_first()
            body={'model':'gpt-5.6-luna','usage':{'prompt_tokens':10,'completion_tokens':3},
                  'choices':[{'finish_reason':'stop','message':{'content':json.dumps(value)}}]}
            if mode=='incomplete' or (mode=='retry_incomplete' and path.parent.name.endswith('_01-000')):
                body['choices'][0]['finish_reason']='length'
            if mode=='mismatch':body['model']='unexpected-model'
            return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')
        return call

    def run_offline(self,mode='valid',after_first=None):
        self.seen=[]
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
             patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=self.transport(mode,after_first)) as transport, \
             patch(self.module+'.conditional_calculations',side_effect=AssertionError('No numerical rerun')), \
             patch(self.module+'.audit_study',side_effect=AssertionError('No full study rerun')):
            result=run_selected_stage_recovery(self.output)
        return result,transport.call_count

    def reseal(self,path,value):
        self.write(path,value);Path(str(path)+'.sha256').write_text(sha(path)+'\n')

    def test_preparation_is_offline_and_freezes_exact_stages(self):
        with patch(self.module+'.Budget',side_effect=AssertionError('No Budget during preparation')), \
             patch(self.module+'.StudyClient.ask',side_effect=AssertionError('No calls during preparation')):
            result=self.prepare()
        self.assertEqual((result['stages'],result['tasks'],result['new_model_calls']),(2,2,0))
        self.assertFalse((self.output/'budget.sqlite3').exists())

    def test_complete_resume_zero_calls_and_old_numerical_verdicts_unchanged(self):
        self.prepare();before={str(p):sha(p) for p in self.old.rglob('*') if p.is_file()}
        result,calls=self.run_offline();self.assertEqual(calls,2)
        self.assertEqual(result['status_counts'],{'SOURCE_CONTRACT_COMPLETE':1,'ROW_REVIEW_COMPLETE':1})
        self.assertEqual(result['unique_primary_missing_rows_reviewed'],1)
        review=next(json.loads(r['messages'][1]['content']) for r in self.seen if 'review_row_ids' in json.loads(r['messages'][1]['content']))
        self.assertEqual(review['independent_numerical_verification']['results'][0]['status'],'FAIL')
        self.assertEqual(review['source']['evidence']['N00000']['value']['status'],'CANNOT VERIFY')
        resumed,calls=self.run_offline('timeout');self.assertEqual(calls,0)
        self.assertEqual(resumed['status_counts'],result['status_counts'])
        self.assertEqual(before,{str(p):sha(p) for p in self.old.rglob('*') if p.is_file()})

    def test_single_contract_correction_and_terminal_contract_stop(self):
        self.prepare();result,calls=self.run_offline('reject')
        self.assertEqual((calls,result['status_counts']),(4,{'CONTRACT_STOP':2}))
        self.assertTrue(any('contract_validation_feedback' in json.loads(r['messages'][1]['content']) for r in self.seen))
        resumed,calls=self.run_offline();self.assertEqual(calls,0);self.assertEqual(resumed['status_counts'],result['status_counts'])

    def test_transport_uncertainty_and_policy_stop_do_not_retry(self):
        for mode in ('timeout','policy'):
            with self.subTest(mode=mode):
                self.output=self.root/('recovery-'+mode);self.prepare();result,calls=self.run_offline(mode)
                self.assertEqual((calls,result['unfinished_stages']),(1,1));self.assertEqual(result['status_counts'],{'STOPPED':1})
                self.assertGreater(result['budget']['estimated_or_reserved_usd'],0)
                resumed,calls=self.run_offline();self.assertEqual(calls,0);self.assertEqual(resumed['unfinished_stages'],1)

    def test_policy_overlap_and_completed_overlay_are_rejected_before_output(self):
        for kind in ('policy','overlay'):
            with self.subTest(kind=kind):
                before=copy.deepcopy(self.policy)
                self.policy[0]['next_action']='provider_policy_hold_no_retry' if kind=='policy' else 'read_existing_completed_recovery_before_new_work'
                self.sync_inputs()
                with self.assertRaisesRegex(PipelineIntegrityError,'overlaps'):self.prepare()
                self.assertFalse(self.output.exists());self.policy=before

    def test_source_drift_and_symlink_fail_before_budget(self):
        self.prepare();metadata=Path(self.queue[0]['source_files']['source_metadata']['path'])
        metadata.write_text('changed software source')
        with patch(self.module+'.Budget',side_effect=AssertionError('No Budget before source checks')):
            with self.assertRaises(PipelineIntegrityError):self.run_offline()
        self.assertFalse((self.output/'budget.sqlite3').exists())

    def test_source_drift_after_response_cannot_trigger_correction(self):
        self.prepare();metadata=Path(self.queue[0]['source_files']['source_metadata']['path'])
        result,calls=self.run_offline('reject',lambda:metadata.write_text('changed software source'))
        self.assertEqual((calls,result['status_counts']),(1,{'STOPPED':1}))

    def test_resealed_spending_increase_fails_before_budget(self):
        self.prepare();path=self.output/'selected-stage-recovery.json';record=json.loads(path.read_text())
        record['config']['max_total_cost_usd']=100;self.reseal(path,record)
        with patch(self.module+'.Budget',side_effect=AssertionError('No Budget after changed allocation')):
            with self.assertRaises(PipelineIntegrityError):self.run_offline()

    def test_task_row_and_stage_substitution_fail_in_preparation(self):
        for kind in ('task','row','stage'):
            with self.subTest(kind=kind):
                before=copy.deepcopy(self.stage_list);target_before=copy.deepcopy(self.targets)
                if kind=='task':self.stage_list[0]['study_key']='human_rna:GSE2'
                elif kind=='stage':self.stage_list[0]['stage']='04_review_9999'
                else:self.targets[self.stage_list[0]['event_id']]=['R99999']
                self.sync_inputs()
                with self.assertRaises(PipelineIntegrityError):self.prepare()
                self.assertFalse(self.output.exists());self.stage_list=before;self.targets=target_before


    def test_completed_contract_correction_preserves_both_attempts(self):
        self.prepare();result,calls=self.run_offline('retry')
        self.assertEqual(calls,4)
        self.assertEqual(result['status_counts'],{'SOURCE_CONTRACT_COMPLETE':1,'ROW_REVIEW_COMPLETE':1})
        self.assertEqual(len(list((self.output/'studies').rglob('contract-00.json'))),2)
        self.assertEqual(self.run_offline()[1],0)

    def test_incomplete_and_model_mismatch_stop_without_contract_retry(self):
        for mode,count in [('incomplete',1),('mismatch',1),('retry_incomplete',2)]:
            with self.subTest(mode=mode):
                self.output=self.root/('recovery-'+mode);self.prepare();result,calls=self.run_offline(mode)
                self.assertEqual(calls,count);self.assertEqual(result['status_counts'],{'STOPPED':1})
                self.assertEqual(self.run_offline()[1],0)

    def test_hard_interruption_resumes_without_uncertain_repeat(self):
        self.prepare()
        with self.assertRaises(KeyboardInterrupt):self.run_offline('interrupt')
        self.assertEqual(len(self.seen),1)
        result,calls=self.run_offline()
        self.assertEqual(calls,0);self.assertEqual(result['status_counts'],{'STOPPED':1})
        self.assertEqual(result['budget']['by_status'][0]['status'],'RESERVED_OR_UNCERTAIN')
        self.assertEqual(self.run_offline()[1],0)

    def test_budget_stop_retains_unfinished_and_sends_zero_calls(self):
        self.pool['max_new_cost_usd']=.000001;self.sync_inputs();self.prepare()
        result,calls=self.run_offline();self.assertEqual(calls,0)
        self.assertEqual(result['unfinished_stages'],1);self.assertEqual(result['status_counts'],{'STOPPED':1})
        self.assertEqual(result['budget']['call_reservations'],0)

    def test_symlinks_live_dangling_and_output_target_unchanged(self):
        self.prepare();source=Path(self.queue[0]['source_files']['source_metadata']['path'])
        original=source.read_bytes();target=self.root/'outside-marker';target.write_bytes(original)
        source.unlink();source.symlink_to(target)
        with self.assertRaises(PipelineIntegrityError):self.run_offline()
        self.assertEqual(target.read_bytes(),original);source.unlink();source.symlink_to(self.root/'missing')
        with self.assertRaises(PipelineIntegrityError):self.run_offline()
        source.unlink();source.write_bytes(original)
        (self.output/'budget.sqlite3').symlink_to(target)
        with self.assertRaises(PipelineIntegrityError):self.run_offline()
        self.assertEqual(target.read_bytes(),original)

    def test_resealed_packet_evidence_and_identity_substitution_fails_before_budget(self):
        self.prepare();manifest_path=self.output/'selected-stage-recovery.json'
        original=json.loads(manifest_path.read_text());entry=original['stages'][0]
        packet_path=Path(entry['packet_binding']['path']);packet=json.loads(packet_path.read_text())
        for field in ('payload','dependencies','target_row_ids','stage'):
            with self.subTest(field=field):
                altered=copy.deepcopy(packet);manifest=copy.deepcopy(original)
                if field=='payload':altered[field]['source']['source_sha256']='f'*64
                elif field=='dependencies':altered[field]=[]
                elif field=='target_row_ids':altered[field]=['R99999']
                else:altered[field]='04_review_9999'
                self.write(packet_path,altered)
                manifest['stages'][0]={**{k:v for k,v in altered.items() if k not in {'system','payload','rows','dependencies'}},
                                      'packet_binding':declaration(packet_path)}
                self.reseal(manifest_path,manifest)
                with patch(self.module+'.Budget',side_effect=AssertionError('No Budget for resealed packet')):
                    with self.assertRaises(PipelineIntegrityError):self.run_offline()
        self.write(packet_path,packet);self.reseal(manifest_path,original)

    def test_completed_stage_uncertain_original_and_size_stop_reject_preparation(self):
        task=self.tasks[0];child=self.old/'studies'/task['study_key'].replace(':','-')
        terminal_path=child/'terminal.json';original=json.loads(terminal_path.read_text())
        for kind in ('complete','uncertain'):
            with self.subTest(kind=kind):
                terminal=copy.deepcopy(original)
                rel='02_blinded_design/complete.json' if kind=='complete' else 'api/source_design_00-000/reservation.json'
                self.write(child/rel,{'software_fixture':True});terminal['output_sha256'][rel]=sha(child/rel)
                self.write(terminal_path,terminal);task['original_terminal_sha256']=sha(terminal_path);self.sync_inputs()
                with self.assertRaises(PipelineIntegrityError):self.prepare()
                self.assertFalse(self.output.exists());(child/rel).unlink()
        self.write(terminal_path,original);task['original_terminal_sha256']=sha(terminal_path)
        config=json.loads((self.old/'config.json').read_text());config['config']['api']['prompt_bytes']=100
        self.write(self.old/'config.json',config);self.plan['historical_config']=declaration(self.old/'config.json');self.sync_inputs()
        with self.assertRaisesRegex(PipelineIntegrityError,'envelope'):self.prepare()
        self.assertFalse(self.output.exists())

    def test_deterministic_enrollment_and_implementation_drift(self):
        self.prepare();manifest=json.loads((self.output/'selected-stage-recovery.json').read_text())
        self.stage_list.reverse();self.sync_inputs();self.output=self.root/'second';self.prepare()
        second=json.loads((self.output/'selected-stage-recovery.json').read_text())
        self.assertEqual([r['stage_id'] for r in manifest['stages']],[r['stage_id'] for r in second['stages']])
        code=self.output/'implementation'/'audit_system/corpus_audit_2026_09_20.py';code.write_text('changed fixture code')
        with patch(self.module+'.Budget',side_effect=AssertionError('No Budget for source drift')):
            with self.assertRaises(PipelineIntegrityError):self.run_offline()

    def test_resealed_resume_coverage_cannot_invent_reviewed_rows(self):
        self.prepare();self.run_offline()
        path=next(self.output.glob('results/*.json'));result=json.loads(path.read_text())
        result['primary_missing_rows_reviewed']=['R99999'];self.reseal(path,result)
        with patch(self.module+'.Budget',side_effect=AssertionError('No Budget for changed resume result')):
            with self.assertRaises(PipelineIntegrityError):self.run_offline()

    def test_source_only_completion_has_no_saved_row_review(self):
        self.stage_list=self.stage_list[:1];self.targets={self.stage_list[0]['event_id']:['R00001']}
        self.plan.update(expected_stages=1,expected_tasks=1);self.sync_inputs();self.prepare()
        result,calls=self.run_offline();self.assertEqual(calls,1)
        self.assertEqual(result['unique_primary_missing_rows_reviewed'],0)
        self.assertEqual(result['reviewed_row_assignments'],0)
        self.assertEqual(result['status_counts'],{'SOURCE_CONTRACT_COMPLETE':1})
        self.assertEqual(self.run_offline()[1],0)

    def test_original_r_csv_and_rejected_text_are_bound(self):
        child=self.old/'studies/human_rna-GSE2';term=child/'terminal.json'
        value=json.loads(term.read_text())
        for rel in ('03b_independent_r/attempt-000/base-r-results.csv','review_0000_00-response.txt'):
            p=child/rel;p.write_text('Original fabricated software witness')
            value['output_sha256'][rel]=sha(p)
        self.write(term,value);self.tasks[1]['original_terminal_sha256']=sha(term);self.sync_inputs();self.prepare()
        for rel in ('03b_independent_r/attempt-000/base-r-results.csv','review_0000_00-response.txt'):
            p=child/rel;original=p.read_bytes();p.write_text('Changed fabricated software witness')
            with self.assertRaises(PipelineIntegrityError):self.run_offline()
            p.write_bytes(original)

    def test_missing_terminal_completed_stage_still_binds_provider_response(self):
        self.prepare();self.run_offline()
        manifest=json.loads((self.output/'selected-stage-recovery.json').read_text())
        packet=next(r for r in manifest['stages'] if r['stage']=='04_review_0000')
        result=self.output/'results'/(packet['stage_id']+'.json')
        result.unlink();Path(str(result)+'.sha256').unlink()
        complete=self.output/'studies/human_rna-GSE2/04_review_0000/complete.json'
        record=json.loads(complete.read_text());record['result']['rows'][0]['row_id']='R99999'
        self.write(complete,record);complete.with_suffix('.sha256').write_text(sha(complete)+'\n')
        with patch(self.module+'.Budget',side_effect=AssertionError('No Budget for changed cached stage')):
            with self.assertRaises(PipelineIntegrityError):self.run_offline()

    def test_timeout_terminal_cannot_be_downgraded_to_contract_stop(self):
        self.prepare();self.run_offline('timeout')
        path=next(self.output.glob('results/*.json'));value=json.loads(path.read_text())
        value.update(status='CONTRACT_STOP',stop_scope='stage');self.reseal(path,value)
        with patch(self.module+'.Budget',side_effect=AssertionError('No Budget for altered stop kind')):
            with self.assertRaises(PipelineIntegrityError):self.run_offline()

    def test_request_argument_mutation_rejected_before_reservation(self):
        import audit_system.corpus_audit_2026_09_20 as module
        for kind in ('actor','system','payload'):
            with self.subTest(kind=kind):
                self.output=self.root/('request-'+kind);self.prepare()
                original=module.StudyClient.ask
                def altered(client,actor,system,payload):
                    if kind=='actor':actor='different_actor_00'
                    elif kind=='system':system+=' changed'
                    else:payload={**payload,'study_id':'GSE999999'}
                    return original(client,actor,system,payload)
                with patch.object(module.StudyClient,'ask',altered):result,calls=self.run_offline()
                self.assertEqual(calls,0);self.assertEqual(result['budget']['call_reservations'],0)
                self.assertEqual(result['status_counts'],{'STOPPED':1})

    def test_sorted_membership_census_preserves_multirow_request_order(self):
        child=self.old/'studies/human_rna-GSE2';request_path=child/'api/review_0000_00-000/request.json'
        request=json.loads(request_path.read_text());payload=json.loads(request['messages'][1]['content'])
        payload['source']['evidence']['R00000']=copy.deepcopy(payload['source']['evidence']['R00002'])
        payload['review_row_ids']=['R00002','R00000']
        request['messages'][1]['content']=json.dumps(payload,ensure_ascii=False,separators=(',',':'))
        self.write(request_path,request)
        terminal_path=child/'terminal.json';terminal=json.loads(terminal_path.read_text())
        terminal['review_failures'][0]['row_ids']=['R00002','R00000']
        terminal['output_sha256']['api/review_0000_00-000/request.json']=sha(request_path)
        self.write(terminal_path,terminal);self.tasks[1]['original_terminal_sha256']=sha(terminal_path)
        self.requests[1].update(original_block_row_ids=['R00000','R00002'],original_request_sha256=sha(request_path),
            current_accounting_request_sha256=hashlib.sha256(json.dumps(request,ensure_ascii=False).encode()).hexdigest())
        self.sync_inputs();self.prepare();result,calls=self.run_offline()
        self.assertEqual(result['reviewed_row_assignments'],2);self.assertEqual(result['unique_primary_missing_rows_reviewed'],1)
        sent=next(json.loads(r['messages'][1]['content']) for r in self.seen if 'review_row_ids' in json.loads(r['messages'][1]['content']))
        self.assertEqual(sent['review_row_ids'],['R00002','R00000'])
        # Historical packets had no explicit row-list field. Add only the new
        # measured sorted field while retaining the unsorted old evidence.
        expected_payload=copy.deepcopy(payload);expected_payload['review_row_ids']=['R00000','R00002']
        old_payload=copy.deepcopy(payload);old_payload.pop('review_row_ids')
        request['messages'][1]['content']=json.dumps(old_payload,ensure_ascii=False,separators=(',',':'))
        self.write(request_path,request);terminal['output_sha256']['api/review_0000_00-000/request.json']=sha(request_path)
        self.write(terminal_path,terminal);self.tasks[1]['original_terminal_sha256']=sha(terminal_path)
        expected_request=copy.deepcopy(request)
        expected_request['messages'][1]['content']=json.dumps(expected_payload,ensure_ascii=False,separators=(',',':'))
        self.requests[1].update(original_request_sha256=sha(request_path),
            current_accounting_request_sha256=hashlib.sha256(json.dumps(expected_request,ensure_ascii=False).encode()).hexdigest())
        self.sync_inputs();self.output=self.root/'legacy-no-explicit-list';self.prepare();result,calls=self.run_offline()
        sent=next(json.loads(r['messages'][1]['content']) for r in self.seen if 'review_row_ids' in json.loads(r['messages'][1]['content']))
        self.assertEqual(sent['review_row_ids'],['R00000','R00002'])
        self.assertEqual(sent['source'],old_payload['source'])
        manifest=json.loads((self.output/'selected-stage-recovery.json').read_text())
        descriptor=next(r for r in manifest['stages'] if r['stage']=='04_review_0000')
        self.assertEqual(descriptor['original_row_ids'],['R00002','R00000'])
        for rows in (['R00002'],['R00000','R00002','R00002']):
            self.requests[1]['original_block_row_ids']=rows;self.sync_inputs()
            self.output=self.root/('negative-'+str(len(rows)))
            with self.assertRaises(PipelineIntegrityError):self.prepare()
            self.assertFalse(self.output.exists())


class SelectedStageContinuationTests(unittest.TestCase):
    """Fabricated transport fixtures; no empirical scientific result or API call."""

    def setUp(self):
        from audit_system import corpus_audit_2026_09_20 as module
        self.m=module
        self.f=SelectedStageRecoveryTests('runTest');self.f.setUp()
        self.addCleanup(self.f.tearDown)
        self.parent=self.f.output
        self.f.prepare()
        original=self.f.transport('reject')
        def oversized_software_reply(command,**kwargs):
            reply=original(command,**kwargs)
            body=json.loads(reply.stdout.rsplit(b'\n',1)[0])
            body['choices'][0]['message']['content']=json.dumps({'fabricated_software_output':'x'*910000})
            return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')
        # Reproduce the predecessor's old cohort-wide classification explicitly.
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=oversized_software_reply) as call, \
                patch(self.m.__name__+'._selected_context_stop',return_value=False):
            stopped=self.m.run_selected_stage_recovery(self.parent)
        self.assertEqual((stopped['terminal_stages'],stopped['unfinished_stages'],call.call_count),(1,1,1))
        self.assertEqual(stopped['status_counts'],{'STOPPED':1})
        self.output=self.f.root/'continuation'
        self.f.write(self.f.root/'terminal-readback.json',{'status':'PASS','scope':'Fabricated software fixture',
            'observation_kind':'TERMINAL_SELECTED_STAGE_READBACK',
            'manifest_sha256':sha(self.parent/'selected-stage-recovery.json'),
            'report_sha256':sha(self.parent/'reports/report-0000.json'),
            'budget_database_sha256':sha(self.parent/'budget.sqlite3'),
            'stages_attempted':1,'stages_remaining':1,'budget':{'call_reservations':1,'uncertain_reservations':0}})
        self.f.write(self.f.root/'protocol.json',{'scope':'Fabricated continuation protocol; not empirical evidence'})
        self.plan={'schema':self.m.CONTINUATION_SCHEMA,'phase':'selected_failed_stages_only',
            'parent_output':str(self.parent),'parent_manifest':declaration(self.parent/'selected-stage-recovery.json'),
            'parent_report':declaration(self.parent/'reports/report-0000.json'),
            'parent_budget':declaration(self.parent/'budget.sqlite3'),
            'terminal_readback':declaration(self.f.root/'terminal-readback.json'),
            'protocol':declaration(self.f.root/'protocol.json'),'output':str(self.output),
            'model':'gpt-6-luna','input_per_million':.25,'output_per_million':.75,
            'expected_stages':1,'expected_tasks':1,'expected_prior_calls':1,
            'expected_prior_cost_usd':stopped['budget']['estimated_or_reserved_usd'],'excluded_study_ids':[]}
        self.plan_path=self.f.root/'continuation-plan.json'
        self.before={str(p.relative_to(self.parent)):sha(p) for p in self.parent.rglob('*') if p.is_file()}

    def prepare(self):
        self.f.write(self.plan_path,self.plan)
        return self.m.prepare_selected_stage_recovery(self.plan_path,self.output)

    def run_offline(self):
        self.f.seen=[];original=self.f.transport()
        def model_reply(command,**kwargs):
            reply=original(command,**kwargs)
            request=json.loads(Path(command[command.index('--data-binary')+1][1:]).read_text())
            body=json.loads(reply.stdout.rsplit(b'\n',1)[0]);body['model']=request['model']
            return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=model_reply) as call:
            result=self.m.run_selected_stage_recovery(self.output)
        return result,call.call_count

    def test_only_untouched_stages_run_with_prior_costs_and_original_bytes_preserved(self):
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('No ledger during preparation')):
            prepared=self.prepare()
        self.assertEqual(prepared['stages'],1)
        self.assertFalse((self.output/'budget.sqlite3').exists())
        self.assertFalse(self.parent.with_name(self.parent.name+'.continuation-owner.json').exists())
        result,calls=self.run_offline()
        self.assertEqual(calls,1)
        self.assertEqual(result['budget']['call_reservations'],2)
        self.assertEqual(result['budget']['max_calls'],4)
        self.assertEqual(result['budget']['max_cost_usd'],.6)
        self.assertEqual(result['allocation_carry_forward']['new_call_reservations'],1)
        self.assertEqual(result['terminal_stages'],1)
        self.assertEqual(result['accepted_scientific_repairs'],0)
        self.assertTrue(all(r['model']=='gpt-6-luna' for r in self.f.seen))
        self.assertEqual(self.before,{str(p.relative_to(self.parent)):sha(p) for p in self.parent.rglob('*') if p.is_file()})
        resumed,calls=self.run_offline()
        self.assertEqual(calls,0)
        self.assertEqual(resumed['budget']['call_reservations'],2)

    def test_sibling_cannot_claim_the_same_allocation(self):
        self.prepare();self.run_offline()
        self.output=self.f.root/'sibling';self.plan['output']=str(self.output)
        self.plan_path=self.f.root/'sibling-plan.json';self.prepare()
        with self.assertRaisesRegex(PipelineIntegrityError,'another continuation'):
            self.run_offline()
        self.assertFalse((self.output/'budget.sqlite3').exists())

    def test_predecessor_lock_prevents_budget_and_transport(self):
        self.prepare()
        with self.m.locked(self.parent/'.selected-stage.lock'):
            with self.assertRaisesRegex(RuntimeError,'already owned'):
                self.run_offline()
        self.assertFalse((self.output/'budget.sqlite3').exists())

    def test_budget_drift_cannot_be_repaired_by_updating_the_report_alone(self):
        self.prepare()
        with sqlite3.connect(self.parent/'budget.sqlite3') as db:
            db.execute('UPDATE reservations SET charged=0')
        with self.assertRaisesRegex(PipelineIntegrityError,'changed'):
            self.run_offline()
        self.assertFalse((self.output/'budget.sqlite3').exists())

    def test_excluded_study_or_changed_counts_refuse_preparation(self):
        manifest=json.loads((self.parent/'selected-stage-recovery.json').read_text())
        terminal=json.loads((self.parent/'reports/report-0000.json').read_text())['results'][0]['stage_id']
        remaining=next(p for p in manifest['stages'] if p['stage_id']!=terminal)
        self.plan['excluded_study_ids']=[remaining['study_key'].split(':')[1]]
        with self.assertRaisesRegex(PipelineIntegrityError,'excluded study'):self.prepare()
        self.assertFalse(self.output.exists())
        self.plan['excluded_study_ids']=[];self.plan['expected_stages']=2
        with self.assertRaisesRegex(PipelineIntegrityError,'count changed'):self.prepare()
        self.assertFalse(self.output.exists())

    def test_a_partial_unsealed_stage_is_never_treated_as_untouched(self):
        manifest=json.loads((self.parent/'selected-stage-recovery.json').read_text())
        terminal=json.loads((self.parent/'reports/report-0000.json').read_text())['results'][0]['stage_id']
        remaining=next(p for p in manifest['stages'] if p['stage_id']!=terminal)
        child=self.parent/'studies'/remaining['study_key'].replace(':','-')
        self.f.write(child/remaining['stage']/'attempt-000/failure.json',{'scope':'Fabricated partial fixture'})
        with self.assertRaisesRegex(PipelineIntegrityError,'partially attempted'):self.prepare()
        self.assertFalse(self.output.exists())

    def test_model_substitution_and_budget_enlargement_fail_before_budget(self):
        self.plan['model']='gpt-5.6-luna'
        with self.assertRaisesRegex(PipelineIntegrityError,'Explicit bounded'):self.prepare()
        self.plan['model']='gpt-6-luna';self.plan['output_per_million']=100000
        with self.assertRaisesRegex(PipelineIntegrityError,'retained allocation'):self.prepare()
        self.assertFalse(self.output.exists())

    def test_sqlite_sidecar_and_extra_old_call_binding_refuse_preparation(self):
        sidecar=Path(str(self.parent/'budget.sqlite3')+'-wal');sidecar.write_bytes(b'fabricated software WAL')
        with self.assertRaisesRegex(PipelineIntegrityError,'sidecar'):self.prepare()
        sidecar.unlink()
        binding=next((self.parent/'studies').rglob('api/*.aggregate-binding.json'))
        value=json.loads(binding.read_text());value['aggregate_key']='f'*64
        self.f.write(binding.with_name('orphan_00-000.aggregate-binding.json'),value)
        with self.assertRaisesRegex(PipelineIntegrityError,'call bindings'):self.prepare()
        self.assertFalse(self.output.exists())

    def test_oversized_correction_is_unresolved_while_the_next_independent_stage_runs(self):
        self.f.output=self.f.root/'prospective';self.f.prepare();self.f.seen=[]
        original=self.f.transport()
        def first_reply_too_large(command,**kwargs):
            reply=original(command,**kwargs)
            if len(self.f.seen)>1:return reply
            body=json.loads(reply.stdout.rsplit(b'\n',1)[0])
            body['choices'][0]['message']['content']=json.dumps({'fabricated_software_output':'x'*910000})
            return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=first_reply_too_large) as call:
            result=self.m.run_selected_stage_recovery(self.f.output)
            resumed=self.m.run_selected_stage_recovery(self.f.output)
        self.assertEqual(call.call_count,2)
        self.assertEqual((result['terminal_stages'],result['unfinished_stages']),(2,0))
        self.assertEqual(result['status_counts']['CONTEXT_STOP'],1)
        self.assertEqual(result['budget']['call_reservations'],2)
        self.assertEqual(resumed['budget']['call_reservations'],2)
        self.assertEqual(result['accepted_scientific_repairs'],0)
        report=json.loads((self.f.output/'reports/report-0000.json').read_text())
        first,second=report['results']
        self.assertEqual(first['status'],'CONTEXT_STOP');self.assertEqual(first['stop_scope'],'stage')
        self.assertEqual(first['primary_missing_rows_reviewed'],[])
        self.assertNotEqual(first['study_key'],second['study_key'])
        self.assertIn(second['status'],{'SOURCE_CONTRACT_COMPLETE','ROW_REVIEW_COMPLETE'})
        child=self.f.output/'studies'/first['study_key'].replace(':','-')
        self.assertFalse(any((child/'api').glob('*_01*')))

    def test_context_stop_requires_a_completed_rejected_reply_and_no_correction_hold(self):
        manifest=json.loads((self.parent/'selected-stage-recovery.json').read_text())
        descriptor=manifest['stages'][0];packet=json.loads(Path(descriptor['packet_binding']['path']).read_text())
        child=self.parent/'studies'/packet['study_key'].replace(':','-')
        self.assertTrue(self.m._selected_context_stop(child,packet,manifest['config'],self.m.SELECTED_CONTEXT_ERROR))
        marker=child/'api'/(packet['actor']+'_01-000.aggregate-binding.json')
        self.f.write(marker,{'scope':'Fabricated uncertain correction marker'})
        self.assertFalse(self.m._selected_context_stop(child,packet,manifest['config'],self.m.SELECTED_CONTEXT_ERROR))
        marker.unlink()
        self.assertFalse(self.m._selected_context_stop(child,packet,manifest['config'],'RuntimeError: provider failure'))
        slot=child/'api'/(packet['actor']+'_00-000')
        response=json.loads((slot/'response.json').read_text())
        response['choices'][0]['message']['refusal']='Fabricated software provider refusal'
        self.f.write(slot/'response.json',response)
        receipt=json.loads((slot/'receipt.json').read_text());receipt['response_sha256']=sha(slot/'response.json')
        self.f.write(slot/'receipt.json',receipt)
        self.assertFalse(self.m._selected_context_stop(child,packet,manifest['config'],self.m.SELECTED_CONTEXT_ERROR))

    def test_explicit_provider_refusal_keeps_cohort_stopped_without_a_correction(self):
        self.f.output=self.f.root/'prospective-refusal';self.f.prepare();self.f.seen=[]
        original=self.f.transport()
        def policy_reply(command,**kwargs):
            reply=original(command,**kwargs)
            body=json.loads(reply.stdout.rsplit(b'\n',1)[0])
            body['choices'][0]['message'].update(
                content=json.dumps({'fabricated_software_output':'x'*910000}),
                refusal='Fabricated software provider refusal')
            return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=policy_reply) as call:
            result=self.m.run_selected_stage_recovery(self.f.output)
            resumed=self.m.run_selected_stage_recovery(self.f.output)
        self.assertEqual(call.call_count,1)
        self.assertEqual((result['terminal_stages'],result['unfinished_stages']),(1,1))
        self.assertEqual(result['status_counts'],{'STOPPED':1})
        self.assertEqual(resumed['status_counts'],{'STOPPED':1})
        self.assertEqual(result['budget']['call_reservations'],1)
        report=json.loads((self.f.output/'reports/report-0000.json').read_text())
        stop=report['results'][0]
        self.assertEqual(stop['stop_scope'],'cohort')
        self.assertIn('OPENAI_PROVIDER_REFUSAL',stop['error'])
        child=self.f.output/'studies'/stop['study_key'].replace(':','-')
        self.assertFalse(any((child/'api').glob('*_01*')))
        self.assertFalse(any((child/stop['stage']).glob('attempt-*/contract-*.json')))
        self.assertTrue(next((child/'api').glob('*/receipt.json')).is_file())

    def test_unbound_reference_receipt_cannot_authorize_another_terminal_state(self):
        receipt=json.loads((self.f.root/'terminal-readback.json').read_text())
        receipt['manifest_sha256']='a'*64;self.f.write(self.f.root/'terminal-readback.json',receipt)
        self.plan['terminal_readback']=declaration(self.f.root/'terminal-readback.json')
        with self.assertRaisesRegex(PipelineIntegrityError,'Independent terminal readback'):self.prepare()
        self.assertFalse(self.output.exists())


class ServiceTierCorpusTests(unittest.TestCase):
    """Offline request/accounting fixtures; dollar values are test guard inputs."""
    module='audit_system.corpus_audit_2026_09_20'

    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.config={'provider':'openai','concurrency':1,'review_rows_per_call':1,
            'validation_retries':1,'max_total_calls':4,'max_total_cost_usd':.1,
            'api':{'model':'gpt-6-luna','service_tier':'default','input_per_million':.2,
                   'output_per_million':.75,'output_tokens':20,'prompt_bytes':10000,
                   'max_calls_per_study':4,'max_cost_per_study':.1,'json_response':True}}

    @staticmethod
    def response(**options):
        body={'model':'gpt-6-luna','usage':{'prompt_tokens':10,'completion_tokens':3},
              'choices':[{'finish_reason':'stop','message':{'content':'{}'}}],**options}
        return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')

    def test_request_builder_matches_wire_binding_and_tiered_reconciliation(self):
        from src.pipeline_runtime_2026_09_10 import openai_request
        budget=Budget(self.root/'budget.sqlite3',.1,4)
        client=StudyClient(self.root,self.config,budget)
        payload={'software_fixture':'unicode μ and "quoted" content'}
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'},clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                      return_value=self.response(service_tier='default')) as transport:
            self.assertEqual(client.ask('source_design_00','Software fixture',payload),'{}')
            self.assertEqual(client.ask('source_design_00','Software fixture',payload),'{}')
            self.assertEqual(transport.call_count,1)
        slot=self.root/'api/source_design_00-000'
        request=json.loads((slot/'request.json').read_text())
        wire=[{'role':'system','content':'Software fixture'},
              {'role':'user','content':json.dumps(payload,ensure_ascii=False,separators=(',',':'))}]
        expected=openai_request('gpt-6-luna',wire,20,json_response=True,service_tier='default')
        self.assertEqual((slot/'request.json').read_bytes(),(json.dumps(expected,indent=2,allow_nan=False)+'\n').encode())
        binding=json.loads((self.root/'api/source_design_00-000.aggregate-binding.json').read_text())
        self.assertEqual(binding['request_file_sha256'],sha(slot/'request.json'))
        reserved=((len(json.dumps(request,ensure_ascii=False).encode())+200)*.2+20*.75)/1e6
        self.assertEqual(binding['reserved_usd'],reserved)
        with budget.connect() as db:
            self.assertEqual(db.execute('SELECT binding_sha256 FROM reservations').fetchone()[0],record_hash(binding))
        result=budget.reconcile(self.root)
        self.assertEqual(result['usage_estimates'],1)
        self.assertEqual(result['uncertain_reservations'],0)
        self.assertAlmostEqual(result['estimated_or_reserved_usd'],(10*.2+3*.75)/1e6)

    def test_legacy_top_level_marker_does_not_change_request_or_aggregate_binding(self):
        from audit_system.corpus_audit_2026_09_20 import _selected_reenvelope
        config=copy.deepcopy(self.config);del config['api']['service_tier']
        config['service_tier']='default'
        validate_configuration(config)
        budget=Budget(self.root/'budget.sqlite3',.1,4)
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'},clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=self.response()):
            StudyClient(self.root,config,budget).ask('review_00','Software fixture',{})
        slot=self.root/'api/review_00-000'
        expected={'model':'gpt-6-luna','messages':[{'role':'system','content':'Software fixture'},
                  {'role':'user','content':'{}'}],'max_completion_tokens':20,
                  'reasoning_effort':'low','store':False,'response_format':{'type':'json_object'}}
        wire_bytes=(json.dumps(expected,indent=2,allow_nan=False)+'\n').encode()
        self.assertEqual((slot/'request.json').read_bytes(),wire_bytes)
        binding=json.loads((self.root/'api/review_00-000.aggregate-binding.json').read_text())
        self.assertEqual(binding['request_file_sha256'],hashlib.sha256(wire_bytes).hexdigest())
        packet={'system':'Software fixture','payload':{}}
        before=_selected_reenvelope([packet],config)
        del config['service_tier']
        self.assertEqual(_selected_reenvelope([packet],config),before)
        self.assertEqual(before[0]['request_sha256'],hashlib.sha256(json.dumps(expected,ensure_ascii=False).encode()).hexdigest())
        self.assertEqual(budget.reconcile(self.root)['usage_estimates'],1)

    def test_invalid_explicit_tiers_stop_before_reservation_and_transport(self):
        for value in (None,False,1,[],{},'','auto','priority'):
            config=copy.deepcopy(self.config);config['api']['service_tier']=value
            with self.subTest(value=value), self.assertRaisesRegex(ValueError,'service_tier'):
                validate_configuration(config)
            with patch(self.module+'.Budget.reserve') as reserve, \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport:
                client=StudyClient(self.root,config,SimpleNamespace(reserve=reserve))
                with self.assertRaisesRegex(ValueError,'service_tier'):
                    client.ask('review_00','Software fixture',{})
                reserve.assert_not_called();transport.assert_not_called()
            self.assertFalse((self.root/'api').exists())

    def test_tier_failure_keeps_full_aggregate_hold_and_never_contract_retries(self):
        for i,options in enumerate(({}, {'service_tier':None}, {'service_tier':'priority'})):
            root=self.root/str(i);budget=Budget(root/'budget.sqlite3',.1,4)
            client=StudyClient(root,self.config,budget)
            with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'},clear=True), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=self.response(**options)) as transport:
                with self.assertRaisesRegex(PipelineIntegrityError,'Returned service tier'):
                    checked_response(client,'review','Software fixture',{},
                        lambda _:self.fail('Rejected tier reached contract validation'),root/'contracts',1)
                with self.assertRaisesRegex(PipelineIntegrityError,'Prior service-tier rejection'):
                    client.ask('review_00','Software fixture',{})
                self.assertEqual(transport.call_count,1)
            slot=root/'api/review_00-000'
            reserve=json.loads((slot/'reservation.json').read_text())['reserved_usd']
            self.assertFalse((slot/'receipt.json').exists())
            self.assertFalse((root/'review_00-response.txt').exists())
            self.assertEqual(budget.report()['by_status'],[{'status':'RESERVED_OR_UNCERTAIN','calls':1,'usd':reserve}])
            reconciled=budget.reconcile(root)
            self.assertEqual(reconciled['uncertain_reservations'],1)
            self.assertEqual(reconciled['usage_estimates'],0)
            self.assertEqual(reconciled['estimated_or_reserved_usd'],reserve)
            self.assertEqual(list(root.glob('contracts/contract-*.json')),[])

    def test_cached_rejected_tier_cannot_settle_aggregate_or_reconcile_as_success(self):
        for i,tier in enumerate((None,'priority')):
            root=self.root/str(i);budget=Budget(root/'initial.sqlite3',.1,4)
            with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'},clear=True), \
                    patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=self.response(service_tier='default')):
                StudyClient(root,self.config,budget).ask('review_00','Software fixture',{})
            slot=root/'api/review_00-000';response=json.loads((slot/'response.json').read_text())
            if tier is None:del response['service_tier']
            else:response['service_tier']=tier
            (slot/'response.json').write_text(json.dumps(response))
            receipt=json.loads((slot/'receipt.json').read_text());receipt['response_sha256']=sha(slot/'response.json')
            (slot/'receipt.json').write_text(json.dumps(receipt))
            binding=json.loads((root/'api/review_00-000.aggregate-binding.json').read_text())
            # Software interruption fixture: a cache exists before aggregate settlement.
            resumed=Budget(root/'resumed.sqlite3',.1,4)
            resumed.reserve(binding['aggregate_key'],binding['reserved_usd'],record_hash(binding))
            with patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport, \
                    patch.object(resumed,'settle',side_effect=AssertionError('No rejected-tier settlement')):
                with self.assertRaisesRegex(PipelineIntegrityError,'Returned service tier'):
                    StudyClient(root,self.config,resumed).ask('review_00','Software fixture',{})
                transport.assert_not_called()
            self.assertEqual(resumed.report()['by_status'][0]['status'],'RESERVED_OR_UNCERTAIN')
            self.assertEqual(resumed.report()['estimated_or_reserved_usd'],binding['reserved_usd'])
            with self.assertRaisesRegex(PipelineIntegrityError,'Rejected service tier cannot have settled usage'):
                resumed.reconcile(root)
            self.assertEqual(resumed.report()['estimated_or_reserved_usd'],binding['reserved_usd'])

    def test_reconciliation_rejects_malformed_tier_even_with_rehashed_binding(self):
        budget=Budget(self.root/'budget.sqlite3',.1,4)
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'},clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=self.response(service_tier='default')):
            StudyClient(self.root,self.config,budget).ask('review_00','Software fixture',{})
        slot=self.root/'api/review_00-000';request=json.loads((slot/'request.json').read_text())
        binding_path=self.root/'api/review_00-000.aggregate-binding.json'
        binding=json.loads(binding_path.read_text())
        for value in (None,False,[],{},'priority'):
            request['service_tier']=value;(slot/'request.json').write_text(json.dumps(request))
            binding['request_file_sha256']=sha(slot/'request.json');binding_path.write_text(json.dumps(binding))
            with self.subTest(value=value), self.assertRaisesRegex(PipelineIntegrityError,'invalid service tier'):
                budget.reconcile(self.root)

    def test_corpus_source_replay_requires_exact_tier_config_and_response_bindings(self):
        fixture=CorpusTests(methodName='runTest');fixture.setUp();self.addCleanup(fixture.tearDown)
        binding,design=fixture.replay_binding();root=Path(binding['root'])
        self.assertEqual(replay_source_design(binding,fixture.index),design)
        with self.assertRaisesRegex(PipelineIntegrityError,'service tier differs'):
            replay_source_design(binding,fixture.index,service_tier='default')
        config_path=root/'config.json';config=json.loads(config_path.read_text())
        config['config']={'api':{'service_tier':'default'}};config_path.write_text(json.dumps(config))
        budget=Budget(fixture.root/'budget.sqlite3',.1,4)
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'},clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=self.response(service_tier='default')):
            StudyClient(root,self.config,budget).ask('source_design_00','Software fixture',{})
        binding['files']={str(p):declaration(p) for p in root.rglob('*') if p.is_file()}
        self.assertEqual(replay_source_design(binding,fixture.index,service_tier='default'),design)
        for filename in ('config.json','api/source_design_00-000/request.json','api/source_design_00-000/response.json'):
            changed=copy.deepcopy(binding);del changed['files'][str(root/filename)]
            with self.subTest(file=filename), self.assertRaisesRegex(PipelineIntegrityError,'binding'):
                replay_source_design(changed,fixture.index,service_tier='default')
        slot=root/'api/source_design_00-000';response=json.loads((slot/'response.json').read_text())
        response['service_tier']='priority';(slot/'response.json').write_text(json.dumps(response))
        receipt=json.loads((slot/'receipt.json').read_text());receipt['response_sha256']=sha(slot/'response.json')
        (slot/'receipt.json').write_text(json.dumps(receipt))
        binding['files']={str(p):declaration(p) for p in root.rglob('*') if p.is_file()}
        with self.assertRaisesRegex(PipelineIntegrityError,'Returned service tier'):
            replay_source_design(binding,fixture.index,service_tier='default')


if __name__=='__main__':unittest.main()


class SelectedStageContinuationV2Tests(unittest.TestCase):
    """Fabricated software fixtures and mocked HTTP; never empirical audit data."""

    def setUp(self):
        from audit_system import corpus_audit_2026_09_20 as module
        self.m=module
        self.f=SelectedStageRecoveryTests('runTest')
        self.f.selected_cases=getattr(self,'fixture_cases',[(1,'02_blinded_design')]+[(n,'04_review_0000') for n in range(2,5)])
        self.f.setUp();self.addCleanup(self.f.tearDown)
        self.f.pool.update(max_new_calls=32,max_new_cost_usd=8)
        self.f.plan.update(expected_stages=4,expected_tasks=4);self.f.sync_inputs()
        self.root=self.f.output;self.parent=self.f.root/'continuation-v1';self.output=self.f.root/'continuation-v2'
        self.f.prepare();self.f.seen=[];original=self.f.transport('reject')
        def oversized_reply(command,**kwargs):
            reply=original(command,**kwargs);body=json.loads(reply.stdout.rsplit(b'\n',1)[0])
            body['choices'][0]['message']['content']=json.dumps({'software_fixture':'x'*910000})
            return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=oversized_reply) as call, \
                patch(self.m.__name__+'._selected_context_stop',return_value=False):
            root_report=self.m.run_selected_stage_recovery(self.root)
        self.assertEqual((call.call_count,root_report['unfinished_stages']),(1,3))
        self.root_readback=self.f.root/'root-readback.json';self.make_readback(self.root,self.root_readback)
        self.protocol=self.f.root/'software-protocol.json'
        self.f.write(self.protocol,{'scope':'Fabricated software-only continuation protocol; no empirical data'})
        self.parent_plan=self.f.root/'continuation-v1-plan.json'
        self.f.write(self.parent_plan,{'schema':self.m.CONTINUATION_SCHEMA,'phase':'selected_failed_stages_only',
            'parent_output':str(self.root),'parent_manifest':declaration(self.root/'selected-stage-recovery.json'),
            'parent_report':declaration(self.root/'reports/report-0000.json'),'parent_budget':declaration(self.root/'budget.sqlite3'),
            'terminal_readback':declaration(self.root_readback),'protocol':declaration(self.protocol),'output':str(self.parent),
            'model':'gpt-6-luna','input_per_million':.25,'output_per_million':.75,
            'expected_stages':3,'expected_tasks':3,'expected_prior_calls':1,
            'expected_prior_cost_usd':root_report['budget']['estimated_or_reserved_usd'],'excluded_study_ids':[]})
        self.m.prepare_selected_stage_recovery(self.parent_plan,self.parent)
        prior,calls=self.run_at(self.parent,mode='length')
        self.assertEqual((calls,prior['status_counts'],prior['unfinished_stages']),(1,{'STOPPED':1},2))
        self.parent_readback=self.f.root/'closed-census.json';self.make_readback(self.parent,self.parent_readback,closed=True)
        self.plan_path=self.f.root/'continuation-v2-plan.json'
        self.plan={'schema':self.m.CONTINUATION_V2_SCHEMA,'phase':'selected_failed_stages_only',
            'parent_output':str(self.parent),'parent_manifest':declaration(self.parent/'selected-stage-recovery.json'),
            'parent_report':declaration(self.parent/'reports/report-0000.json'),'parent_budget':declaration(self.parent/'budget.sqlite3'),
            'terminal_readback':declaration(self.parent_readback),'protocol':declaration(self.protocol),'output':str(self.output),
            'model':'gpt-6-luna','input_per_million':.25,'output_per_million':.75,
            'expected_stages':2,'expected_tasks':2,'expected_prior_calls':2,
            'expected_prior_cost_usd':prior['budget']['estimated_or_reserved_usd'],
            'excluded_study_ids':['GSE298220','GSE295890'],
            'allocation_root_manifest':declaration(self.root/'selected-stage-recovery.json'),
            'ancestors':[self.ancestor(self.root,self.f.root/'plan.json',self.root_readback,
                                     declaration(self.root.with_name(self.root.name+'.continuation-owner.json'))),
                         self.ancestor(self.parent,self.parent_plan,self.parent_readback,None)],
            'output_limit_policy':self.m.OUTPUT_LIMIT_POLICY,'service_tier':'default'}
        self.ancestor_bytes={str(p):sha(p) for directory in (self.root,self.parent) for p in directory.rglob('*') if p.is_file()}
        self.root_owner=self.root.with_name(self.root.name+'.continuation-owner.json')
        self.root_owner_bytes=self.root_owner.read_bytes()

    def ancestor(self,directory,plan,readback,owner):
        return {'output':str(directory),'manifest':declaration(directory/'selected-stage-recovery.json'),
            'plan':declaration(plan),'report':declaration(directory/'reports/report-0000.json'),
            'budget':declaration(directory/'budget.sqlite3'),'terminal_readback':declaration(readback),'successor_owner':owner}

    def make_readback(self,directory,path,closed=False):
        report=json.loads((directory/'reports/report-0000.json').read_text())
        if closed:
            value={'schema':'clockbase-historical-gpt6-stopped-readback-v1','status':'PASS','failed_checks':[],'error':None,
                'scope':'Fabricated software census; not empirical evidence',
                'observation_kind':'CLOSED_RESULT_CENSUS_WITH_UNVERIFIED_PROCESS_EXIT',
                'process_returncode':None,'exact_process_exit_time':None,
                'manifest_sha256':sha(directory/'selected-stage-recovery.json'),
                'report_sha256':sha(directory/'reports/report-0000.json'),'budget_sha256':sha(directory/'budget.sqlite3'),
                'summaries':{'process':{'process_returncode':None,'exact_process_exit_time':None,'exit_cause':'UNVERIFIED',
                                      'terminal_process_receipt_present':False},
                    'stages':{'unsealed_stage_activity':[],'enrolled_stages':report['enrolled_stages'],
                        'status_counts':report['status_counts'],'sealed_stages':report['terminal_stages'],
                        'unfinished_stages':report['unfinished_stages']},
                    'budget':{'combined_calls':report['budget']['call_reservations'],'unsettled_reservations':0,
                        'combined_configured_usage_estimate_usd':str(report['budget']['estimated_or_reserved_usd']),
                        'ledger_statuses':{'USAGE_ESTIMATE':report['budget']['call_reservations']}}}}
        else:
            value={'status':'PASS','scope':'Fabricated software readback; not empirical evidence',
                'observation_kind':'TERMINAL_SELECTED_STAGE_READBACK',
                'manifest_sha256':sha(directory/'selected-stage-recovery.json'),
                'report_sha256':sha(directory/'reports/report-0000.json'),'budget_database_sha256':sha(directory/'budget.sqlite3'),
                'stages_attempted':report['terminal_stages'],'stages_remaining':report['unfinished_stages'],
                'budget':{'call_reservations':report['budget']['call_reservations'],'uncertain_reservations':0}}
        self.f.write(path,value)

    def prepare(self):
        self.f.write(self.plan_path,self.plan)
        return self.m.prepare_selected_stage_recovery(self.plan_path,self.output)

    def run_at(self,directory,mode='valid',mutate=None,after=None):
        self.f.seen=[];original=self.f.transport();first_actor=[None]
        def reply(command,**kwargs):
            if mode=='timeout':raise subprocess.TimeoutExpired(['software-fixture-transport'],190)
            response=original(command,**kwargs)
            path=Path(command[command.index('--data-binary')+1][1:]);request=json.loads(path.read_text())
            if first_actor[0] is None:first_actor[0]=path.parent.parent.parent
            first=path.parent.parent.parent==first_actor[0]
            body=json.loads(response.stdout.rsplit(b'\n',1)[0]);body.update(model=request['model'],service_tier='default')
            body['choices'][0]['message']['refusal']=None
            if first and mode=='correction_length' and path.parent.name.endswith('_00-000'):
                body['choices'][0]['message']['content']='{}'
            if first and (mode=='length' or mode=='correction_length' and path.parent.name.endswith('_01-000')):
                body['choices'][0].update(finish_reason='length',message={'content':'','refusal':None})
                body['usage']={'prompt_tokens':10,'completion_tokens':16000,'total_tokens':16010,
                               'completion_tokens_details':{'reasoning_tokens':16000}}
                if mutate is not None:mutate(body,request,path)
            if after is not None:after(body,request,path)
            return SimpleNamespace(returncode=0,stdout=json.dumps(body).encode()+b'\n200')
        with patch.dict(os.environ,{'OPENAI_API_KEY':'software-fixture-key'}), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=reply) as call, \
                patch(self.m.__name__+'.conditional_calculations',side_effect=AssertionError('No numerical execution in software fixture')), \
                patch(self.m.__name__+'.audit_study',side_effect=AssertionError('No full-study execution in software fixture')):
            result=self.m.run_selected_stage_recovery(directory)
        return result,call.call_count

    def run_offline(self,mode='valid',mutate=None,after=None):return self.run_at(self.output,mode,mutate,after)

    def first_packet(self,directory=None):
        directory=directory or self.output
        manifest=json.loads((directory/'selected-stage-recovery.json').read_text())
        packet=json.loads(Path(manifest['stages'][0]['packet_binding']['path']).read_text())
        return manifest,packet,directory/'studies'/packet['study_key'].replace(':','-')

    def assert_ancestors_unchanged(self):
        self.assertEqual(self.root_owner.read_bytes(),self.root_owner_bytes)
        self.assertEqual(self.ancestor_bytes,{str(p):sha(p) for directory in (self.root,self.parent)
                                             for p in directory.rglob('*') if p.is_file()})

    def assert_prepare_refused(self,pattern=None):
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('No Budget during invalid preparation')), \
                patch(self.m.__name__+'.StudyClient',side_effect=AssertionError('No client during invalid preparation')):
            with self.assertRaisesRegex(PipelineIntegrityError,pattern or '.'):
                self.prepare()
        self.assertFalse(self.output.exists())
        self.assertFalse(self.parent.with_name(self.parent.name+'.continuation-owner.json').exists())

    def test_exact_untouched_union_offline_prepare_and_unchanged_ancestors(self):
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('No Budget during preparation')), \
                patch(self.m.__name__+'.StudyClient',side_effect=AssertionError('No client during preparation')):
            prepared=self.prepare()
        self.assertEqual((prepared['stages'],prepared['tasks'],prepared['new_model_calls']),(2,2,0))
        self.assertFalse((self.output/'budget.sqlite3').exists())
        self.assertFalse(self.parent.with_name(self.parent.name+'.continuation-owner.json').exists())
        root_manifest=json.loads((self.root/'selected-stage-recovery.json').read_text())
        root_results=json.loads((self.root/'reports/report-0000.json').read_text())['results']
        parent_results=json.loads((self.parent/'reports/report-0000.json').read_text())['results']
        attempted={r['stage_id'] for r in root_results+parent_results}
        manifest=json.loads((self.output/'selected-stage-recovery.json').read_text())
        self.assertEqual([p['stage_id'] for p in manifest['stages']],
                         [p['stage_id'] for p in root_manifest['stages'] if p['stage_id'] not in attempted])
        for descriptor in manifest['stages']:
            inherited=next(p for p in root_manifest['stages'] if p['stage_id']==descriptor['stage_id'])
            for field in ('study_key','stage','original_row_ids','target_row_ids'):
                self.assertEqual(descriptor[field],inherited[field])
        result,calls=self.run_offline()
        self.assertEqual(calls,2)
        self.assertEqual(result['budget']['call_reservations'],4) # Two physical ancestors; not 1+2 cumulative snapshots.
        self.assertEqual((result['budget']['max_calls'],result['budget']['max_cost_usd']),(32,8))
        self.assertEqual(result['allocation_carry_forward']['prior_calls'],2)
        self.assertEqual(result['allocation_carry_forward']['new_call_reservations'],2)
        self.assertEqual(result['allocation_carry_forward']['ancestors'],self.plan['ancestors'])
        self.assertEqual(self.run_offline('timeout')[1],0)
        self.assert_ancestors_unchanged()

    def test_chain_metadata_shape_and_policy_reject_before_budget(self):
        original=copy.deepcopy(self.plan)
        changes={
            'missing_root':lambda p:p['ancestors'].pop(0),
            'duplicate':lambda p:p['ancestors'].insert(0,copy.deepcopy(p['ancestors'][0])),
            'reordered':lambda p:p['ancestors'].reverse(),
            'cycle':lambda p:p['ancestors'][-1].update(output=p['output']),
            'missing_owner':lambda p:p['ancestors'][0].update(successor_owner=None),
            'last_owner':lambda p:p['ancestors'][-1].update(successor_owner=p['ancestors'][0]['successor_owner']),
            'wrong_parent':lambda p:p.update(parent_output=p['ancestors'][0]['output']),
            'missing_policy':lambda p:p.pop('output_limit_policy'),
            'unknown_policy':lambda p:p.update(output_limit_policy='arbitrary_incomplete_v1'),
            'tier':lambda p:p.update(service_tier='priority'),
            'held_exclusion_missing':lambda p:p.update(excluded_study_ids=[]),
            'stage_override':lambda p:p.update(stage_ids=[]),
            'enlarged_price':lambda p:p.update(output_per_million=100),
            'wrong_model':lambda p:p.update(model='gpt-5.6-luna'),
            'wrong_count':lambda p:p.update(expected_stages=3),
            'double_count':lambda p:p.update(expected_prior_calls=3)}
        for label,change in changes.items():
            with self.subTest(label=label):
                self.plan=copy.deepcopy(original);change(self.plan);self.assert_prepare_refused()
        self.plan=original

    def test_symlink_alias_and_overlap_are_refused(self):
        original=copy.deepcopy(self.plan);alias=self.f.root/'root-alias';alias.symlink_to(self.root,target_is_directory=True)
        for path in (alias,self.root/'nested',self.root.parent):
            with self.subTest(path=str(path)):
                self.plan=copy.deepcopy(original);self.plan['ancestors'][0]['output']=str(path)
                self.assert_prepare_refused()
        self.plan=original

    def test_all_ancestor_locks_are_required_and_released(self):
        self.prepare()
        for directory in (self.root,self.parent):
            with self.subTest(ancestor=directory.name),self.m.locked(directory/'.selected-stage.lock'), \
                    patch(self.m.__name__+'.Budget',side_effect=AssertionError('No Budget behind another owner')):
                with self.assertRaisesRegex(RuntimeError,'already owned'):self.run_offline()
            self.assertFalse((self.output/'budget.sqlite3').exists())
            self.assertFalse(self.parent.with_name(self.parent.name+'.continuation-owner.json').exists())
        self.assertEqual(self.run_offline()[1],2)

    def test_existing_edges_and_sibling_claims_are_immutable(self):
        original=self.root_owner.read_bytes();self.root_owner.write_text('{}\n')
        self.assert_prepare_refused();self.root_owner.write_bytes(original)
        self.prepare();self.run_offline()
        self.assertEqual(self.root_owner.read_bytes(),original)
        self.output=self.f.root/'sibling';self.plan['output']=str(self.output);self.plan_path=self.f.root/'sibling-plan.json'
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('Sibling must not reserve')):
            with self.assertRaisesRegex(PipelineIntegrityError,'another continuation'):self.prepare()
        self.assertFalse(self.output.exists())

    def test_claim_crash_allows_only_exact_restart(self):
        self.prepare()
        with patch(self.m.__name__+'.Budget',side_effect=RuntimeError('Software interruption after claim')):
            with self.assertRaisesRegex(RuntimeError,'Software interruption'):self.run_offline()
        owner=self.parent.with_name(self.parent.name+'.continuation-owner.json');before=owner.read_bytes()
        self.assertFalse((self.output/'budget.sqlite3').exists())
        manifest_path=self.output/'selected-stage-recovery.json';saved=manifest_path.read_bytes()
        manifest=json.loads(saved);manifest['software_tamper']=True;self.f.reseal(manifest_path,manifest)
        with self.assertRaisesRegex(PipelineIntegrityError,'another continuation'):self.run_offline()
        manifest_path.write_bytes(saved);Path(str(manifest_path)+'.sha256').write_text(sha(manifest_path)+'\n')
        self.assertEqual(self.run_offline()[1],2);self.assertEqual(owner.read_bytes(),before)

    def test_ancestor_sidecars_and_partial_activity_are_never_untouched(self):
        for directory in (self.root,self.parent):
            for suffix in ('-wal','-shm','-journal'):
                with self.subTest(ancestor=directory.name,sidecar=suffix):
                    sidecar=Path(str(directory/'budget.sqlite3')+suffix);sidecar.write_bytes(b'software fixture')
                    self.assert_prepare_refused('sidecar');sidecar.unlink()
        manifest=json.loads((self.parent/'selected-stage-recovery.json').read_text())
        packet=json.loads(Path(manifest['stages'][1]['packet_binding']['path']).read_text())
        child=self.parent/'studies'/packet['study_key'].replace(':','-')
        marker=child/packet['stage']/'attempt-000/failure.json';self.f.write(marker,{'scope':'Unsealed software fixture'})
        self.assert_prepare_refused('partially attempted')

    def test_held_metadata_denies_packet_reads(self):
        original=(self.root/'selected-stage-recovery.json').read_bytes();plan=copy.deepcopy(self.plan)
        reader=self.m._selected_json
        def metadata_only(binding):
            if '/packets/' in binding['path']:raise AssertionError('Held packet body was opened')
            return reader(binding)
        for study in ('GSE298220','GSE295890'):
            with self.subTest(study=study):
                manifest=json.loads(original);manifest['stages'][0]['study_key']='human_rna:'+study
                self.f.reseal(self.root/'selected-stage-recovery.json',manifest)
                self.plan=copy.deepcopy(plan)
                binding=declaration(self.root/'selected-stage-recovery.json')
                self.plan['ancestors'][0]['manifest']=binding;self.plan['allocation_root_manifest']=binding
                with patch(self.m.__name__+'._selected_json',side_effect=metadata_only):
                    self.assert_prepare_refused('held or excluded')

    def test_closed_census_retains_null_exit_and_rejects_pid_only_evidence(self):
        before=self.parent_readback.read_bytes();self.prepare();self.run_offline()
        self.assertEqual(self.parent_readback.read_bytes(),before)
        census=json.loads(before)
        self.assertIsNone(census['process_returncode']);self.assertIsNone(census['exact_process_exit_time'])
        ancestor=self.plan['ancestors'][-1];report=json.loads(Path(ancestor['report']['path']).read_text())
        rows=self.m._selected_prior_budget(ancestor['budget'],report['budget'])
        changes=[lambda r:r.update(process_returncode=0),
                 lambda r:r['summaries']['process'].update(process_returncode=0),
                 lambda r:r['summaries']['process'].update(exit_cause='PID_ABSENT'),
                 lambda r:r['summaries']['stages'].update(unsealed_stage_activity=['unverified']),
                 lambda r:r['summaries']['stages'].update(sealed_stages=0),
                 lambda r:r.update(observation_kind='PID_ABSENT')]
        for index,change in enumerate(changes):
            with self.subTest(index=index):
                value=json.loads(before);change(value);path=self.f.root/f'invalid-readback-{index}.json';self.f.write(path,value)
                with self.assertRaises(PipelineIntegrityError):
                    self.m._selected_readback_adapter(declaration(path),ancestor,report,1,2,rows)

    def test_initial_output_limit_is_zero_coverage_and_next_stage_runs(self):
        self.prepare();result,calls=self.run_offline('length')
        self.assertEqual(calls,2);self.assertEqual(result['output_limit_stages'],1)
        self.assertEqual((result['terminal_stages'],result['unfinished_stages']),(2,0))
        self.assertEqual(result['budget']['call_reservations'],4)
        report=json.loads((self.output/'reports/report-0000.json').read_text());first,second=report['results']
        self.assertEqual((first['status'],first['stop_scope']),('OUTPUT_LIMIT_STOP','stage'))
        self.assertEqual(first['reviewed_row_ids'],[]);self.assertEqual(first['primary_missing_rows_reviewed'],[])
        self.assertIn(second['status'],{'ROW_REVIEW_COMPLETE','SOURCE_CONTRACT_COMPLETE'})
        self.assertEqual(result['reviewed_row_assignments'],len(second['reviewed_row_ids']))
        self.assertEqual((result['accepted_scientific_repairs'],result['scientific_error_rate'],result['scientific_acceptance']),
                         (0,None,'CANNOT VERIFY'))
        manifest,packet,child=self.first_packet()
        self.assertFalse(any((child/'api').glob(packet['actor']+'_01*')))
        self.assertFalse(any(child.glob(packet['actor']+'_*response.txt')))
        self.assertFalse(any((child/packet['stage']).glob('attempt-*/result.json')))
        self.assertFalse((child/packet['stage']/'complete.json').exists())
        saved={str(p):sha(p) for p in (self.output/'results').rglob('*') if p.is_file()}
        resumed,calls=self.run_offline('timeout')
        self.assertEqual(calls,0);self.assertEqual(resumed['budget'],result['budget'])
        self.assertEqual(saved,{str(p):sha(p) for p in (self.output/'results').rglob('*') if p.is_file()})
        self.assert_ancestors_unchanged()

    def test_correction_output_limit_retains_both_calls_and_exact_feedback(self):
        self.prepare();result,calls=self.run_offline('correction_length')
        self.assertEqual((calls,result['output_limit_stages'],result['unfinished_stages']),(3,1,0))
        self.assertEqual(result['budget']['call_reservations'],5)
        manifest,packet,child=self.first_packet()
        self.assertFalse(any((child/'api').glob(packet['actor']+'_02*')))
        self.assertEqual(len(list((child/packet['stage']).glob('attempt-*/contract-00.json'))),1)
        expected=self.m._selected_correction_payload(child,packet,manifest['config'])
        request=json.loads((child/'api'/(packet['actor']+'_01-000')/'request.json').read_text())
        self.assertEqual(json.loads(request['messages'][1]['content']),expected)
        self.assertEqual(self.run_offline()[1],0)
        feedback=next((child/packet['stage']).glob('attempt-*/contract-00.json'))
        saved=json.loads(feedback.read_text());saved['error']='ValueError: Fabricated different rejection';self.f.write(feedback,saved)
        self.assertFalse(self.m._selected_output_limit_stop(child,packet,manifest['config'],
            ValueError(self.m.SELECTED_OUTPUT_LIMIT_MESSAGE),self.output/'budget.sqlite3'))

    def test_refusal_length_keeps_cohort_stopped(self):
        self.prepare()
        result,calls=self.run_offline('length',lambda body,*_:body['choices'][0]['message'].update(refusal='software refusal'))
        self.assertEqual((calls,result['status_counts'],result['unfinished_stages']),(1,{'STOPPED':1},1))
        self.assertEqual(self.run_offline()[1],0)

    def test_missing_tier_length_keeps_cohort_stopped(self):
        self.prepare();result,calls=self.run_offline('length',lambda body,*_:body.pop('service_tier'))
        self.assertEqual((calls,result['status_counts']),(1,{'STOPPED':1}))

    def test_nondefault_tier_length_keeps_cohort_stopped(self):
        self.prepare();result,calls=self.run_offline('length',lambda body,*_:body.update(service_tier='priority'))
        self.assertEqual((calls,result['status_counts']),(1,{'STOPPED':1}))

    def test_model_substitution_length_keeps_cohort_stopped(self):
        self.prepare();result,calls=self.run_offline('length',lambda body,*_:body.update(model='gpt-5.6-luna'))
        self.assertEqual((calls,result['status_counts']),(1,{'STOPPED':1}))

    def test_other_finish_reason_keeps_cohort_stopped(self):
        self.prepare();result,calls=self.run_offline('length',lambda body,*_:body['choices'][0].update(finish_reason='content_filter'))
        self.assertEqual((calls,result['status_counts']),(1,{'STOPPED':1}))

    def test_timeout_never_retries_or_advances_the_cohort(self):
        self.prepare();result,calls=self.run_offline('timeout')
        self.assertEqual((calls,result['status_counts'],result['unfinished_stages']),(1,{'STOPPED':1},1))
        self.assertEqual(self.run_offline()[1],0)
        self.assertEqual(next(r['calls'] for r in result['budget']['by_status'] if r['status']=='RESERVED_OR_UNCERTAIN'),1)

    def test_narrow_witness_denies_ambiguous_raw_metadata(self):
        self.prepare();self.run_offline('length');manifest,packet,child=self.first_packet()
        slot=child/'api'/(packet['actor']+'_00-000');response_path=slot/'response.json';receipt_path=slot/'receipt.json'
        raw=response_path.read_bytes();saved_receipt=receipt_path.read_bytes()
        def witness():return self.m._selected_output_limit_stop(child,packet,manifest['config'],
            ValueError(self.m.SELECTED_OUTPUT_LIMIT_MESSAGE),self.output/'budget.sqlite3')
        self.assertTrue(witness())
        changes={
            'no_choices':lambda r:r.update(choices=[]),
            'extra_choice':lambda r:r['choices'].append(copy.deepcopy(r['choices'][0])),
            'missing_refusal':lambda r:r['choices'][0]['message'].pop('refusal'),
            'refusal':lambda r:r['choices'][0]['message'].update(refusal='software refusal'),
            'empty_refusal_string':lambda r:r['choices'][0]['message'].update(refusal=''),
            'nonempty_content':lambda r:r['choices'][0]['message'].update(content='{}'),
            'null_content':lambda r:r['choices'][0]['message'].update(content=None),
            'tool_call':lambda r:r['choices'][0]['message'].update(tool_calls=[{'software':True}]),
            'function_call':lambda r:r['choices'][0]['message'].update(function_call={'name':'software'}),
            'below_cap':lambda r:r['usage'].update(completion_tokens=15999),
            'missing_reasoning':lambda r:r['usage'].pop('completion_tokens_details'),
            'reasoning_below_cap':lambda r:r['usage']['completion_tokens_details'].update(reasoning_tokens=15999),
            'float_reasoning':lambda r:r['usage']['completion_tokens_details'].update(reasoning_tokens=16000.0),
            'bool_reasoning':lambda r:r['usage']['completion_tokens_details'].update(reasoning_tokens=True),
            'negative_usage':lambda r:r['usage'].update(prompt_tokens=-1),
            'extra_negative_usage':lambda r:r['usage'].update(prompt_tokens_details={'cached_tokens':-1}),
            'bad_total':lambda r:r['usage'].update(total_tokens=0),
            'missing_tier':lambda r:r.pop('service_tier'),
            'wrong_tier':lambda r:r.update(service_tier='flex'),
            'null_finish':lambda r:r['choices'][0].update(finish_reason=None),
            'tool_finish':lambda r:r['choices'][0].update(finish_reason='tool_calls'),
            'returned_model':lambda r:r.update(model='other-model')}
        for name,change in changes.items():
            with self.subTest(name=name):
                value=json.loads(raw);change(value);self.f.write(response_path,value)
                receipt=json.loads(saved_receipt);receipt.update(response_sha256=sha(response_path),usage=value['usage'])
                self.f.write(receipt_path,receipt);self.assertFalse(witness())
        response_path.write_bytes(raw);receipt_path.write_bytes(saved_receipt);self.assertTrue(witness())

    def test_exact_exception_and_failure_tag_take_precedence(self):
        self.prepare();self.run_offline('length');manifest,packet,child=self.first_packet()
        message=self.m.SELECTED_OUTPUT_LIMIT_MESSAGE
        class DerivedValueError(ValueError):pass
        for error in (RuntimeError(message),PipelineIntegrityError(message),DerivedValueError(message),
                      ValueError('prefix '+message),ValueError(message+' suffix')):
            with self.subTest(type=type(error).__name__,text=str(error)):
                self.assertFalse(self.m._selected_output_limit_stop(child,packet,manifest['config'],error,self.output/'budget.sqlite3'))
        for tag in (message,'RuntimeError: '+message,self.m.SELECTED_OUTPUT_LIMIT_ERROR+' suffix'):
            self.assertFalse(self.m._selected_output_limit_stop(child,packet,manifest['config'],tag,
                self.output/'budget.sqlite3',cached=True))
        failure=next((child/packet['stage']).glob('attempt-*/failure.json'));saved=json.loads(failure.read_text())
        saved['error']='RuntimeError: '+message;self.f.write(failure,saved)
        self.assertFalse(self.m._selected_output_limit_stop(child,packet,manifest['config'],ValueError(message),self.output/'budget.sqlite3'))

    def test_saved_transport_integrity_and_completion_artifacts_deny_witness(self):
        self.prepare();self.run_offline('length');manifest,packet,child=self.first_packet()
        slot=child/'api'/(packet['actor']+'_00-000')
        def witness():return self.m._selected_output_limit_stop(child,packet,manifest['config'],
            ValueError(self.m.SELECTED_OUTPUT_LIMIT_MESSAGE),self.output/'budget.sqlite3')
        for path in (slot/'provider-error.json',slot/'transport-error.json',
                     child/(packet['actor']+'_00-response.txt'),child/packet['stage']/'complete.json',
                     child/packet['stage']/'attempt-000/result.json'):
            with self.subTest(path=path.name):
                path.write_text('{}\n');self.assertFalse(witness());path.unlink()
        receipt_path=slot/'receipt.json';raw=receipt_path.read_bytes()
        for field,value in (('response_sha256','0'*64),('estimated_cost_usd',0),('billed_cost_usd',.01),
                            ('requested_model','other-model'),('returned_model','other-model')):
            with self.subTest(field=field):
                receipt=json.loads(raw);receipt[field]=value;self.f.write(receipt_path,receipt);self.assertFalse(witness())
        receipt_path.write_bytes(raw)
        request=slot/'request.json';saved=request.read_bytes();request.write_text('{}\n');self.assertFalse(witness());request.write_bytes(saved)
        self.assertTrue(witness())

    def test_partial_new_stage_is_never_automatically_repeated(self):
        self.prepare();manifest,packet,child=self.first_packet()
        marker=child/packet['stage']/'attempt-000/failure.json';self.f.write(marker,{'scope':'Software interrupted stage'})
        result,calls=self.run_offline()
        self.assertEqual((calls,result['status_counts'],result['unfinished_stages']),(0,{'STOPPED':1},1))
        self.assertEqual(self.run_offline()[1],0)

    def test_initial_witness_recovery_after_result_seal_interruption_has_no_repeat(self):
        self.prepare();original=self.m._selected_seal;interrupted=[False]
        def seal(path,value):
            if path.parent.name=='results' and not interrupted[0]:
                interrupted[0]=True;raise KeyboardInterrupt('Software interruption before result seal')
            return original(path,value)
        with patch(self.m.__name__+'._selected_seal',side_effect=seal):
            with self.assertRaises(KeyboardInterrupt):self.run_offline('length')
        self.assertEqual(len(self.f.seen),1)
        result,calls=self.run_offline()
        self.assertEqual((calls,result['output_limit_stages'],result['unfinished_stages']),(1,1,0))
        self.assertEqual(self.run_offline()[1],0)

    def test_ancestor_drift_after_prepare_and_during_call_cannot_advance(self):
        self.prepare();budget=self.parent/'budget.sqlite3';saved=budget.read_bytes()
        with sqlite3.connect(budget) as db:db.execute('UPDATE reservations SET charged=0')
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('No Budget after ancestor drift')):
            with self.assertRaises(PipelineIntegrityError):self.run_offline()
        budget.write_bytes(saved)
        def drift_after_first(*_):
            if len(self.f.seen)==1:
                with sqlite3.connect(budget) as db:db.execute('UPDATE reservations SET charged=0')
        result,calls=self.run_offline(after=drift_after_first)
        self.assertEqual((calls,result['status_counts']),(1,{'STOPPED':1}))
        self.assertEqual(result['reviewed_row_assignments'],0)

    def test_joint_study_and_aggregate_caps_are_not_reset(self):
        config={'api':{'prompt_bytes':900000,'input_per_million':.25,'output_tokens':16000,'output_per_million':.75,
                       'max_calls_per_study':64,'max_cost_per_study':10},'max_total_calls':7326,'max_total_cost_usd':1500,
                'validation_retries':1}
        # Explicit arithmetic fixture of the declared envelope; not measured calls or charges.
        packets=[{'study_key':'human_rna:GSE'+str(i//2+1)} for i in range(1721)]
        rows=[(str(i),0,91.11490345/2719,'USAGE_ESTIMATE','fixture') for i in range(2719)]
        self.m._selected_chain_envelope(config,packets,rows,{})
        total=91.11490345+1721*2*((900000+200)*.25+16000*.75)/1e6
        self.assertAlmostEqual(total,907.04100345,places=8)
        for changed in ({'max_total_calls':6160},{'max_total_cost_usd':907}):
            with self.subTest(changed=changed):
                with self.assertRaises(PipelineIntegrityError):self.m._selected_chain_envelope({**config,**changed},packets,rows,{})
        for prior in ((63,0),(0,9.9)):
            with self.subTest(prior=prior):
                with self.assertRaisesRegex(PipelineIntegrityError,'per-study'):
                    self.m._selected_chain_envelope(config,[{'study_key':'human_rna:GSE1'}],[],{'human_rna:GSE1':prior})

    def test_resume_gap_before_sealed_cohort_stop_refuses_before_budget(self):
        self.prepare()
        def second_stop(body,*_):
            if len(self.f.seen)==2:body['choices'][0]['finish_reason']='content_filter'
        result,calls=self.run_offline(after=second_stop)
        self.assertEqual((calls,result['status_counts'].get('STOPPED')),(2,1))
        manifest,packet,child=self.first_packet();saved=self.f.root/'preserved-gap-fixture';saved.mkdir()
        child.rename(saved/'first-stage-files')
        for suffix in ('.json','.json.sha256'):
            path=self.output/'results'/(packet['stage_id']+suffix);path.rename(saved/path.name)
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('Inconsistent resume reached Budget')):
            with self.assertRaisesRegex(PipelineIntegrityError,'prefix'):self.run_offline()

    def test_resume_activity_beyond_next_unsealed_stage_refuses_before_budget(self):
        self.prepare();manifest=json.loads((self.output/'selected-stage-recovery.json').read_text())
        descriptor=manifest['stages'][1]
        path=self.output/'studies'/descriptor['study_key'].replace(':','-')/descriptor['stage']/'attempt-000/failure.json'
        self.f.write(path,{'scope':'Fabricated later partial stage'})
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('Non-prefix activity reached Budget')):
            with self.assertRaisesRegex(PipelineIntegrityError,'prefix'):self.run_offline()

    def refresh_mutated_parent_budget_fixture(self):
        """Reseal a deliberately corrupted software fixture, never empirical evidence."""
        report_path=self.parent/'reports/report-0000.json';report=json.loads(report_path.read_text())
        with sqlite3.connect(self.parent/'budget.sqlite3') as db:
            entries=db.execute('SELECT status,charged FROM reservations').fetchall()
        statuses=sorted({status for status,_ in entries})
        report['budget'].update(call_reservations=len(entries),estimated_or_reserved_usd=sum(cost for _,cost in entries),
            by_status=[{'status':status,'calls':sum(s==status for s,_ in entries),
                        'usd':sum(cost for s,cost in entries if s==status)} for status in statuses])
        self.f.reseal(report_path,report);self.make_readback(self.parent,self.parent_readback,closed=True)
        for key,path in [('budget',self.parent/'budget.sqlite3'),('report',report_path),('terminal_readback',self.parent_readback)]:
            self.plan['ancestors'][-1][key]=declaration(path)
            self.plan[{'budget':'parent_budget','report':'parent_report','terminal_readback':'terminal_readback'}[key]]=declaration(path)
        self.plan['expected_prior_cost_usd']=report['budget']['estimated_or_reserved_usd']

    def test_uncertain_ancestor_cannot_be_resealed_into_eligibility(self):
        with sqlite3.connect(self.root/'budget.sqlite3') as db:root_key=db.execute('SELECT id FROM reservations').fetchone()[0]
        with sqlite3.connect(self.parent/'budget.sqlite3') as db:
            db.execute("UPDATE reservations SET status='RESERVED_OR_UNCERTAIN',charged=reserved WHERE id!=?",(root_key,))
        self.refresh_mutated_parent_budget_fixture()
        self.assert_prepare_refused('uncertain inherited call')

    def test_changed_inherited_charge_cannot_be_resealed_into_eligibility(self):
        with sqlite3.connect(self.root/'budget.sqlite3') as db:root_key=db.execute('SELECT id FROM reservations').fetchone()[0]
        with sqlite3.connect(self.parent/'budget.sqlite3') as db:
            db.execute('UPDATE reservations SET charged=0 WHERE id=?',(root_key,))
        self.refresh_mutated_parent_budget_fixture()
        self.assert_prepare_refused('lost, changed')

    def test_missing_or_orphaned_physical_calls_do_not_reset_the_ledger(self):
        before={directory:sha(directory/'budget.sqlite3') for directory in (self.root,self.parent)}
        root_call=next((self.root/'studies').rglob('api/*.aggregate-binding.json'));saved=root_call.read_bytes()
        root_call.unlink();self.assert_prepare_refused();root_call.write_bytes(saved)
        parent_call=next((self.parent/'studies').rglob('api/*.aggregate-binding.json'))
        orphan=parent_call.with_name('orphan_00-000.aggregate-binding.json');orphan.write_bytes(parent_call.read_bytes())
        self.assert_prepare_refused('outside the attempted stage census');orphan.unlink()
        receipt=next((self.parent/'studies').rglob('api/*/receipt.json'));saved=receipt.read_bytes()
        receipt.unlink();self.assert_prepare_refused();receipt.write_bytes(saved)
        self.assertEqual(before,{directory:sha(directory/'budget.sqlite3') for directory in (self.root,self.parent)})

    def test_later_v2_owner_promotion_passes_prefix_but_unsupported_stop_stays_closed(self):
        self.prepare();result,calls=self.run_offline('length',lambda body,*_:body['choices'][0].update(finish_reason='content_filter'))
        self.assertEqual((calls,result['unfinished_stages']),(1,1))
        current=self.output;current_plan=self.plan_path;readback=self.f.root/'later-v2-readback.json'
        self.make_readback(current,readback,closed=True)
        prior_plan=copy.deepcopy(self.plan)
        ancestors=copy.deepcopy(prior_plan['ancestors'])
        self.assertIsNone(ancestors[-1]['successor_owner'])
        ancestors[-1]['successor_owner']=declaration(self.parent.with_name(self.parent.name+'.continuation-owner.json'))
        ancestors.append(self.ancestor(current,current_plan,readback,None))
        self.output=self.f.root/'later-v2';self.plan_path=self.f.root/'later-v2-plan.json'
        self.plan.update(output=str(self.output),parent_output=str(current),
            parent_manifest=declaration(current/'selected-stage-recovery.json'),
            parent_report=declaration(current/'reports/report-0000.json'),parent_budget=declaration(current/'budget.sqlite3'),
            terminal_readback=declaration(readback),ancestors=ancestors,
            expected_stages=1,expected_tasks=1,expected_prior_calls=3,
            expected_prior_cost_usd=result['budget']['estimated_or_reserved_usd'])
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('Unsupported terminal state reached Budget')):
            with self.assertRaisesRegex(PipelineIntegrityError,'lacks separately admissible context or output-limit evidence'):
                self.prepare()
        self.assertFalse(self.output.exists())
        self.assertFalse(current.with_name(current.name+'.continuation-owner.json').exists())
        # Reaching this terminal-evidence rejection proves that the preceding null-to-bound
        # promotion and the complete immutable prefix passed their checks. No new run is admitted.
        self.assertEqual(json.loads(current_plan.read_text()),prior_plan)

    def test_settled_usage_over_reservation_keeps_cohort_stopped(self):
        self.prepare()
        def excessive_usage(body,*_):body['usage'].update(prompt_tokens=100000000,total_tokens=100016000)
        result,calls=self.run_offline('length',excessive_usage)
        self.assertEqual((calls,result['status_counts'],result['unfinished_stages']),(1,{'STOPPED':1},1))
        report=json.loads((self.output/'reports/report-0000.json').read_text())
        self.assertIn('PipelineIntegrityError',report['results'][0]['error'])
        self.assertEqual(result['reviewed_row_assignments'],0)

    def test_malformed_length_choices_do_not_trigger_contract_correction(self):
        self.prepare();result,calls=self.run_offline('length',lambda body,*_:body.update(choices=[]))
        self.assertEqual((calls,result['status_counts'],result['unfinished_stages']),(1,{'STOPPED':1},1))
        _,packet,child=self.first_packet()
        self.assertFalse(any((child/'api').glob(packet['actor']+'_01*')))

    def test_nonfinite_closed_census_cost_is_not_independent_evidence(self):
        value=json.loads(self.parent_readback.read_text())
        value['summaries']['budget']['combined_configured_usage_estimate_usd']='NaN'
        self.f.write(self.parent_readback,value)
        binding=declaration(self.parent_readback)
        self.plan['terminal_readback']=binding;self.plan['ancestors'][-1]['terminal_readback']=binding
        self.assert_prepare_refused('Closed census')

    def test_changed_strict_completion_envelope_is_rejected_before_budget(self):
        self.prepare();path=self.output/'selected-stage-recovery.json';saved=path.read_bytes()
        for key,value in (('output_tokens',16001),('max_calls_per_study',65),('max_cost_per_study',11)):
            with self.subTest(key=key):
                manifest=json.loads(saved);manifest['config']['api'][key]=value;self.f.reseal(path,manifest)
                with patch(self.m.__name__+'.Budget',side_effect=AssertionError('Changed envelope reached Budget')):
                    with self.assertRaisesRegex(ValueError,'output-limit policy envelope'):self.run_offline()
        path.write_bytes(saved);Path(str(path)+'.sha256').write_text(sha(path)+'\n')

    def test_attempted_packet_cannot_be_reintroduced_by_resealing_new_manifest(self):
        self.prepare();manifest_path=self.output/'selected-stage-recovery.json';manifest=json.loads(manifest_path.read_text())
        prior=json.loads((self.parent/'selected-stage-recovery.json').read_text())
        manifest['stages'][0]=prior['stages'][0];self.f.reseal(manifest_path,manifest)
        with patch(self.m.__name__+'.Budget',side_effect=AssertionError('An attempted stage reached Budget')):
            with self.assertRaisesRegex(PipelineIntegrityError,'planning records'):self.run_offline()

    def test_source_output_limits_preserve_zero_row_semantics(self):
        for mode in ('length','correction_length'):
            with self.subTest(mode=mode):
                fixture=SelectedStageContinuationV2Tests('runTest')
                fixture.fixture_cases=[(number,'02_blinded_design') for number in range(1,5)]
                fixture.setUp();self.addCleanup(fixture.doCleanups)
                fixture.prepare();report,calls=fixture.run_offline(mode)
                self.assertEqual(calls,2 if mode=='length' else 3)
                self.assertEqual(report['status_counts'],{'OUTPUT_LIMIT_STOP':1,'SOURCE_CONTRACT_COMPLETE':1})
                self.assertEqual((report['reviewed_row_assignments'],report['unique_primary_missing_rows_reviewed']),(0,0))
                self.assertEqual((report['accepted_scientific_repairs'],report['scientific_error_rate']),(0,None))
                self.assertEqual(fixture.run_offline()[1],0)

    def test_correction_requests_cannot_gain_third_call_or_changed_payload(self):
        self.prepare();self.run_offline('correction_length');manifest,packet,child=self.first_packet()
        slot=child/'api'/(packet['actor']+'_01-000');request=slot/'request.json';raw=request.read_bytes()
        def witness():return self.m._selected_output_limit_stop(child,packet,manifest['config'],
            ValueError(self.m.SELECTED_OUTPUT_LIMIT_MESSAGE),self.output/'budget.sqlite3')
        for field in ('previous_response','contract_validation_feedback'):
            with self.subTest(field=field):
                value=json.loads(raw);payload=json.loads(value['messages'][1]['content'])
                payload[field]='software-mutated';value['messages'][1]['content']=json.dumps(payload)
                self.f.write(request,value);self.assertFalse(witness())
        request.write_bytes(raw)
        binding=slot.with_name(slot.name+'.aggregate-binding.json')
        third=binding.with_name(packet['actor']+'_02-000.aggregate-binding.json');third.write_bytes(binding.read_bytes())
        self.assertFalse(witness());third.unlink();self.assertTrue(witness())

    def test_ancestor_drift_during_length_reply_is_a_cohort_integrity_stop(self):
        self.prepare();budget=self.parent/'budget.sqlite3';original_budget=budget.read_bytes()
        def drift_after_first(*_):
            if len(self.f.seen)==1:
                with sqlite3.connect(budget) as db:db.execute('UPDATE reservations SET charged=0')
        result,calls=self.run_offline('length',after=drift_after_first)
        self.assertEqual((calls,result['status_counts'],result['unfinished_stages']),(1,{'STOPPED':1},1))
        report=json.loads((self.output/'reports/report-0000.json').read_text())
        self.assertEqual(report['results'][0]['stop_scope'],'cohort')
        self.assertIn('PipelineIntegrityError',report['results'][0]['error'])
        self.assertEqual(result['reviewed_row_assignments'],0)
        manifest,packet,child=self.first_packet()
        failure=next((child/packet['stage']).glob('attempt-*/failure.json'))
        self.assertEqual(json.loads(failure.read_text())['error'],self.m.SELECTED_OUTPUT_LIMIT_ERROR)
        raw_response=child/'api'/(packet['actor']+'_00-000')/'response.json'
        self.assertEqual(json.loads(raw_response.read_text())['choices'][0]['finish_reason'],'length')
        retained={str(path):sha(path) for path in child.rglob('*') if path.is_file()}
        budget.write_bytes(original_budget)
        resumed,calls=self.run_offline()
        self.assertEqual((calls,resumed['status_counts'],resumed['unfinished_stages']),(0,{'STOPPED':1},1))
        self.assertEqual(retained,{str(path):sha(path) for path in child.rglob('*') if path.is_file()})
