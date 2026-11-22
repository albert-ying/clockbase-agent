"""Software fixtures for source-readback failure patterns; no empirical results.

Named regressions reproduce group/abstention structure only. Model responses and
statistics below are test inputs, not new model calls or biological findings.
"""
import copy
import json
import unittest

from audit_system.corpus_audit_2026_09_20 import (
    CODE_FILES, SCIENTIFIC_CHECKS, evidence_index, review_request, validate_review,
)
from audit_system.corpus_repairs_2026_09_20 import _candidate


class AttributionGuardTests(unittest.TestCase):
    def fixture(self, *, abstain=False, cases=None, controls=None, proposed_cases=None,
                proposed_controls=None, basis='observed_method_error', status='FAIL',
                action='propose_rerun', proposal='Software-only proposed rerun'):
        cases=['GSM1','GSM2','GSM3'] if cases is None else cases
        controls=['GSM4','GSM5','GSM6'] if controls is None else controls
        packet={'study_id':'GSE1','source_url':'https://example.invalid/software-fixture',
            'source_sha256':'a'*64,'samples':[
                {'sample_id':sid,'condition':'treated' if sid in cases else 'control'}
                for sid in cases+controls]}
        index=evidence_index(packet)
        source_ref=next(k for k in index['evidence'] if k.startswith('S'))
        design={'study_id':'GSE1','comparisons':[{'id':'declared',
            'case_ids':cases if proposed_cases is None else proposed_cases,
            'control_ids':controls if proposed_controls is None else proposed_controls}]}
        fields={'status':'CANNOT VERIFY' if abstain else 'RESULT',
            'case_ids':json.dumps(cases),'control_ids':json.dumps(controls),
            'statistical_test':'' if abstain else 'Welch two-sample t-test',
            'p_value':'' if abstain else '0.4','effect_size':'' if abstain else '0.1',
            'mean_difference':'' if abstain else '0.2','clock':'software-fixture-clock'}
        rows=[{'row_id':'R1','fields':fields,'artifact_sha256':'b'*64,'logical_record':2,
            'identity_binding':{'identity_status':'PASS','bound_study_key':'human_rna:GSE1'}}]
        payload=review_request(index,design,rows,[],[])
        check={'status':'CANNOT VERIFY','reason':'Software fixture scope',
            'evidence_ids':[source_ref],'unresolved':['No scientific validation in this fixture']}
        model={'study_id':'GSE1','checks':{k:{**copy.deepcopy(check),'scope':'Software fixture'} for k in SCIENTIFIC_CHECKS},
            'rows':[{'row_id':'R1','status':status,'reason':'Opaque allegation without diagnostic keywords',
                'evidence_ids':['R1',source_ref],'unresolved':[] if status=='PASS' else ['Software fixture uncertainty'],
                'action':action,'failure_basis':basis,'repair_target_scope':'same_original_target',
                'source_comparison_ids':['declared'],'proposed_repair':proposal}],
            'design_disagreements':[]}
        return model,payload,rows

    def validate(self, model, payload, rows):
        before=copy.deepcopy((model,payload,rows))
        result=validate_review(model,payload,rows)
        self.assertEqual((model,payload,rows),before,'Neither source rows nor raw model data may change')
        self.assertEqual(result['scientific_acceptance'],'CANNOT VERIFY')
        return result['rows'][0]

    def test_gse253002_abstention_cannot_be_an_observed_method_error(self):
        model,payload,rows=self.fixture(abstain=True,
            cases=['GSM8012239','GSM8012240','GSM8012241'],
            controls=['GSM8012236','GSM8012237','GSM8012238'],
            proposed_cases=['GSM8012239'],proposed_controls=['GSM8012236'])
        result=self.validate(model,payload,rows)
        self.assertEqual(result['status'],'CANNOT VERIFY')
        self.assertEqual(result['model_status'],'FAIL')
        self.assertEqual(result['model_action'],'propose_rerun')
        self.assertEqual(result['model_row'],model['rows'][0])
        self.assertEqual(result['action'],'unresolved')
        self.assertIsNone(result['proposed_repair'])
        self.assertEqual(result['source_comparison_ids'],[])
        self.assertEqual(result['original_group_binding']['status'],'FAIL')
        self.assertEqual(result['separate_analysis_proposal']['comparison_ids'],['declared'])
        candidate=_candidate(rows[0],result,payload['frozen_source_design'],[],{}, {'status':'PASS'})
        self.assertEqual(candidate['status'],'CANNOT VERIFY')
        self.assertIsNone(candidate['calculation'])
        self.assertFalse(candidate['accepted_scientific_repair'])

    def test_abstention_blocks_method_error_even_when_groups_match(self):
        args=self.fixture(abstain=True,action='withhold')
        result=self.validate(*args)
        self.assertEqual(result['original_group_binding']['status'],'PASS')
        self.assertEqual(result['status'],'CANNOT VERIFY')
        self.assertEqual(result['source_comparison_ids'],[])
        self.assertEqual(result['repair_target_scope'],'no_repair_proposed')
        self.assertEqual(result['model_action'],'withhold')

    def test_abstention_allows_explicitly_missing_or_blank_result_fields(self):
        for replacement in ('absent',None,'  '):
            with self.subTest(replacement=replacement):
                model,payload,rows=self.fixture(abstain=True)
                for key in ('statistical_test','p_value','effect_size','mean_difference'):
                    if replacement=='absent':rows[0]['fields'].pop(key)
                    else:rows[0]['fields'][key]=replacement
                self.assertEqual(self.validate(model,payload,rows)['status'],'CANNOT VERIFY')

    def test_supplied_statistic_or_method_is_not_silently_blank(self):
        for key,value in (('p_value','0'),('effect_size','0'),('effect',0),('mean_difference','0'),
                          ('statistical_test','Welch t-test'),('pvalue','0.2'),('p_value','NA')):
            with self.subTest(key=key,value=value):
                model,payload,rows=self.fixture(abstain=True)
                rows[0]['fields'][key]=value
                result=self.validate(model,payload,rows)
                self.assertEqual(result['status'],'FAIL')
                self.assertEqual(result['action'],'propose_rerun')

    def test_non_abstention_observed_error_is_not_dismissed(self):
        for status in ('RESULT','FAIL','unknown',None):
            with self.subTest(status=status):
                model,payload,rows=self.fixture(abstain=True)
                rows[0]['fields']['status']=status
                self.assertEqual(self.validate(model,payload,rows)['status'],'FAIL')

    def test_other_data_record_errors_remain_separately_possible(self):
        for basis in ('source_contradiction','invalid_numeric_or_record','claim_exceeds_evidence'):
            with self.subTest(basis=basis):
                result=self.validate(*self.fixture(abstain=True,basis=basis))
                self.assertEqual(result['status'],'FAIL')
                self.assertEqual(result['failure_basis'],basis)
        model,payload,rows=self.fixture(abstain=True,basis='invalid_numeric_or_record')
        rows[0]['fields']['p_value']='1.1'
        self.assertEqual(self.validate(model,payload,rows)['status'],'FAIL')

    def test_guard_does_not_read_rationale_keywords(self):
        for reason in ('No test was run', 'Paired t-test was definitely wrong', 'arbitrary text'):
            with self.subTest(reason=reason):
                model,payload,rows=self.fixture(abstain=True)
                model['rows'][0]['reason']=reason
                self.assertEqual(self.validate(model,payload,rows)['status'],'CANNOT VERIFY')
                rows[0]['fields']['status']='RESULT'
                self.assertEqual(self.validate(model,payload,rows)['status'],'FAIL')

    def test_gse293225_subset_is_not_original_verification_or_repair(self):
        model,payload,rows=self.fixture(cases=['GSM8878902','GSM8878903','GSM8878904'],
            controls=['GSM8878899','GSM8878900','GSM8878901'],proposed_cases=['GSM8878903'],
            basis='source_contradiction')
        result=self.validate(model,payload,rows)
        self.assertEqual(result['status'],'FAIL','A group mismatch alone cannot dismiss another allegation')
        self.assertEqual(result['model_row'],model['rows'][0])
        binding=result['original_group_binding']
        self.assertEqual(binding['status'],'FAIL')
        self.assertEqual(binding['comparisons'][0]['membership_difference']['case_ids'],
                         {'original_only':['GSM8878902','GSM8878904'],'proposed_only':[]})
        self.assertEqual(result['source_comparison_ids'],[])
        proposal=result['separate_analysis_proposal']
        self.assertEqual(proposal['scope'],'different_target')
        self.assertFalse(proposal['original_verification'])
        self.assertFalse(proposal['original_result_repair'])

    def test_unmatched_group_cannot_support_pass_retention(self):
        args=self.fixture(proposed_cases=['GSM1'],status='PASS',action='retain',proposal=None,basis=None)
        result=self.validate(*args)
        self.assertEqual(result['status'],'CANNOT VERIFY')
        self.assertEqual(result['model_status'],'PASS')
        self.assertEqual(result['model_action'],'retain')
        self.assertEqual(result['action'],'unresolved')
        self.assertTrue(result['unresolved'])

    def test_exact_membership_allows_reordered_groups_without_scientific_approval(self):
        args=self.fixture(proposed_cases=['GSM3','GSM1','GSM2'],proposed_controls=['GSM6','GSM4','GSM5'])
        result=self.validate(*args)
        self.assertEqual(result['original_group_binding']['status'],'PASS')
        self.assertEqual(result['status'],'FAIL')
        self.assertEqual(result['action'],'propose_rerun')
        self.assertEqual(result['source_comparison_ids'],['declared'])
        self.assertNotIn('separate_analysis_proposal',result)
        self.assertEqual(result['repair_acceptance'],'CANNOT VERIFY')

    def test_matched_pass_row_keeps_its_narrow_review_without_scientific_approval(self):
        args=self.fixture(status='PASS',action='retain',proposal=None,basis=None)
        result=self.validate(*args)
        self.assertEqual(result['original_group_binding']['status'],'PASS')
        self.assertEqual(result['status'],'PASS')
        self.assertEqual(result['action'],'retain')
        self.assertEqual(result['repair_acceptance'],'CANNOT VERIFY')

    def test_invalid_record_allegation_is_not_auto_cleared_for_abstaining_subset(self):
        # The empirical GSE253002 reviewer used this broader basis. This narrow
        # guard blocks subset attribution, but cannot adjudicate the allegation.
        args=self.fixture(abstain=True,basis='invalid_numeric_or_record',
                          proposed_cases=['GSM1'],proposed_controls=['GSM4'])
        result=self.validate(*args)
        self.assertEqual(result['status'],'FAIL')
        self.assertEqual(result['original_group_binding']['status'],'FAIL')
        self.assertEqual(result['source_comparison_ids'],[])
        self.assertEqual(result['action'],'unresolved')
        self.assertIsNone(result['proposed_repair'])

    def test_case_control_orientation_must_match(self):
        result=self.validate(*self.fixture(proposed_cases=['GSM4','GSM5','GSM6'],
                                          proposed_controls=['GSM1','GSM2','GSM3']))
        self.assertEqual(result['original_group_binding']['status'],'FAIL')
        self.assertEqual(result['action'],'unresolved')

    def test_missing_or_malformed_original_groups_do_not_establish_matching(self):
        for value in (None,'GSM1,GSM2,GSM3','["GSM1",','[]','["GSM1","GSM1"]',
                      '[1,2,3]','[" GSM1","GSM2","GSM3"]','{"case":"GSM1"}'):
            with self.subTest(value=value):
                model,payload,rows=self.fixture()
                rows[0]['fields']['case_ids']=value
                result=self.validate(model,payload,rows)
                self.assertEqual(result['original_group_binding']['status'],'CANNOT VERIFY')
                self.assertEqual(result['separate_analysis_proposal']['scope'],'cannot_bind')
                self.assertEqual(result['source_comparison_ids'],[])
        model,payload,rows=self.fixture();rows[0]['fields'].pop('control_ids')
        self.assertEqual(self.validate(model,payload,rows)['original_group_binding']['status'],'CANNOT VERIFY')
        model,payload,rows=self.fixture();rows[0]['fields']=None
        self.assertEqual(self.validate(model,payload,rows)['original_group_binding']['status'],'CANNOT VERIFY')

    def test_overlap_and_malformed_frozen_comparison_do_not_bind(self):
        model,payload,rows=self.fixture();rows[0]['fields']['control_ids']='["GSM1","GSM4"]'
        self.assertEqual(self.validate(model,payload,rows)['original_group_binding']['status'],'CANNOT VERIFY')
        model,payload,rows=self.fixture();payload['frozen_source_design']['comparisons'][0]['case_ids']=['GSM1','GSM1']
        self.assertEqual(self.validate(model,payload,rows)['original_group_binding']['status'],'CANNOT VERIFY')

    def test_every_referenced_comparison_must_have_exact_original_groups(self):
        model,payload,rows=self.fixture()
        payload['frozen_source_design']['comparisons'].append({'id':'subset','case_ids':['GSM1'],
            'control_ids':['GSM4','GSM5','GSM6']})
        model['rows'][0]['source_comparison_ids']=['declared','subset']
        result=self.validate(model,payload,rows)
        self.assertEqual([c['status'] for c in result['original_group_binding']['comparisons']],['PASS','FAIL'])
        self.assertEqual(result['source_comparison_ids'],[])

    def test_absent_or_duplicate_comparison_binding_is_not_an_original_repair(self):
        for refs in ([],['declared','declared']):
            with self.subTest(refs=refs):
                model,payload,rows=self.fixture();model['rows'][0]['source_comparison_ids']=refs
                result=self.validate(model,payload,rows)
                self.assertEqual(result['original_group_binding']['status'],'CANNOT VERIFY')
                self.assertEqual(result['source_comparison_ids'],[])

    def test_no_requested_comparison_does_not_add_group_gate(self):
        model,payload,rows=self.fixture(status='CANNOT VERIFY',action='unresolved',proposal=None,basis=None)
        model['rows'][0]['source_comparison_ids']=[];rows[0]['fields'].pop('case_ids')
        self.assertNotIn('original_group_binding',self.validate(model,payload,rows))

    def test_original_model_action_survives_a_later_identity_restriction(self):
        model,payload,rows=self.fixture(abstain=True)
        payload['source']['evidence']['R1']['identity_binding']=None
        result=self.validate(model,payload,rows)
        self.assertEqual(result['model_action'],'propose_rerun')
        self.assertEqual(result['model_row'],model['rows'][0])

    def test_recovered_dnam_verifier_is_in_new_runtime_inventory(self):
        self.assertIn('audit_system/report_2026_09_10.py',CODE_FILES)


if __name__=='__main__':unittest.main()
