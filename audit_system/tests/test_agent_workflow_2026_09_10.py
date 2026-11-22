"""Software fixtures exercise audit control flow; no paid calls or accuracy claims."""
import copy
import csv
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest
from unittest.mock import patch
from audit_system.tests import require_development_artifact

from audit_system import agent_workflow_2026_09_10 as workflow
from audit_system.source_packet import (replay_saved_outputs, article_methods,
    context_documents, CONTEXT_VERSION, read_study_context, compact_source_packet, expand_source_packet,
    source_metadata_summary, verify_requested_comparison, propose_source_comparisons,
    source_grouping_diagnostics)


class WorkflowTests(unittest.TestCase):
    def test_empty_molecular_input_prevents_model_call_and_records_withholding(self):
        module_path=Path(workflow.__file__).resolve().parents[1]/'reliability_audit_2026-06-16/measurement_audit_v0.py'
        require_development_artifact(module_path)
        matrix=self.root/'empty.csv'
        matrix.write_text(',A,B,C,D\n')
        with patch.object(workflow,'request_once') as call:
            result=workflow.run(self.path,'FIXTURE','drug versus none','tokenswitch',self.root/'run',matrix=matrix)
        call.assert_not_called()
        self.assertEqual(result['status'],'FAIL')
        self.assertEqual(result['generation_attempts'],0)
        self.assertEqual(result['execution']['action'],'withhold_inference')

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.path = self.root/'metadata.csv'
        self.path.write_text('GSM,title,characteristics_ch1,pvalue,harmonized.treatment,contact_email\n'
                             'A,case1,treatment: drug,.01,old verdict,private@example.test\n'
                             'B,case2,treatment: drug,.02,old verdict,private@example.test\n'
                             'C,control1,treatment: none,.3,old verdict,private@example.test\n'
                             'D,control2,treatment: none,.4,old verdict,private@example.test\n')
        self.packet = workflow.source_packet(self.path, 'FIXTURE', 'drug versus none')
        self.design = {'study_id': 'FIXTURE', 'status': 'DESIGN_PROPOSED',
                       'case_ids': ['A','B'], 'control_ids': ['C','D'], 'excluded_ids': [],
                       'biological_unit': 'software fixture', 'pairing': 'independent',
                       'findings': [], 'unresolved': [], 'proposed_action': 'recompute_welch',
                       'evidence': [{'sample_id': s['sample_id'], 'field': 'title', 'quote': s['title']}
                                    for s in self.packet['samples']]}

    def test_raw_only_packet_has_no_previous_annotations_outcomes_or_contacts(self):
        payload = json.dumps(self.packet)
        for excluded in ['.01','private@example','old verdict','harmonized.treatment']:
            self.assertNotIn(excluded, payload)
        self.assertEqual(len(self.packet['samples']), 4)

    def test_source_packet_preserves_all_numbered_protocols_without_outcome_columns(self):
        row = {'GSM': 'A', 'title': 'software fixture', 'extract_protocol_ch1': 'First extraction step',
               'extract_protocol_ch1.1': 'Second extraction step', 'treatment_protocol_ch1.1': 'Second treatment step',
               'growth_protocol_ch1.1': 'Second culture step', 'pvalue': '.123456',
               'data_processing': 'First processing step\nExact second line', 'data_processing.1': 'Second processing step',
               'label_ch1': 'Cy3', 'label_protocol_ch1.1': 'Second label step',
               'hyb_protocol': 'Hybridization method', 'scan_protocol.1': 'Second scan step',
               'data_processing.pvalue': 'PRIOR_OUTCOME_SENTINEL',
               'extract_protocol_ch1.pvalue': 'PRIOR_OUTCOME_SENTINEL', 'contact_email': 'private@example.test'}
        path = self.root / 'protocols.csv'
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row)); writer.writeheader(); writer.writerow(row)
        packet = workflow.source_packet(path, 'FIXTURE', 'Reconstruct source design')
        observed = packet['samples'][0]
        for key in ('extract_protocol_ch1', 'extract_protocol_ch1.1', 'treatment_protocol_ch1.1', 'growth_protocol_ch1.1',
                    'data_processing', 'data_processing.1', 'label_ch1', 'label_protocol_ch1.1',
                    'hyb_protocol', 'scan_protocol.1'):
            self.assertEqual(observed[key], row[key])
        for withheld in ('pvalue', 'PRIOR_OUTCOME_SENTINEL', 'private@example.test'):
            self.assertNotIn(withheld, json.dumps(packet))

    def test_source_descriptions_retain_replicate_class_without_derived_outcome_fields(self):
        row = {'GSM': 'GSM1', 'title': 'software fixture',
               'description': 'Technical_replicate', 'description.1': 'donor_1',
               'description.2': 'paired aliquot', 'description.pvalue': 'HIDDEN_P_VALUE',
               'description.verdict': 'HIDDEN_VERDICT', 'description.１': 'NON_ASCII_SUFFIX',
               'pvalue': 'HIDDEN_P_VALUE', 'contact_email': 'private@example.test'}
        path = self.root/'replicate-description.csv'
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row)); writer.writeheader(); writer.writerow(row)
        packet = workflow.source_packet(path, 'GSE1', 'source design')
        for field in ('description', 'description.1', 'description.2'):
            self.assertEqual(packet['samples'][0][field], row[field])
        serialized = json.dumps(packet)
        for hidden in ('HIDDEN_P_VALUE', 'HIDDEN_VERDICT', 'private@example.test', 'NON_ASCII_SUFFIX'):
            self.assertNotIn(hidden, serialized)

    def design_set_fixture(self):
        packet = copy.deepcopy(self.packet)
        packet['samples'].extend([{'sample_id': sid, 'title': 'second ' + sid,
                                   'characteristics_ch1': 'treatment: second'} for sid in ('E', 'F')])
        children = []
        for name, case in [('drug', ['A', 'B']), ('second', ['E', 'F'])]:
            child = {key: copy.deepcopy(value) for key, value in self.design.items() if key != 'study_id'}
            child.update(comparison_name=name, target=name + ' versus none', case_ids=case,
                         excluded_ids=sorted(set('ABCDEF') - set(case + ['C', 'D'])),
                         evidence=[{'sample_id': row['sample_id'], 'field': 'title', 'quote': row['title']}
                                   for row in packet['samples'] if row['sample_id'] in case + ['C', 'D']])
            children.append(child)
        return packet, {'schema': workflow.DESIGN_SET_SCHEMA, 'study_id': 'FIXTURE',
                        'status': 'DESIGN_PROPOSED', 'comparisons': children,
                        'unassigned_ids': {'remaining': True}, 'findings': [], 'unresolved': [], 'evidence': []}

    def test_multiple_comparisons_preserve_shared_controls_and_stable_identity(self):
        packet, raw = self.design_set_fixture()
        result = workflow.parse_design(json.dumps(raw), packet)
        self.assertEqual(len(result['comparisons']), 2)
        self.assertEqual(result['unassigned_ids'], [])
        self.assertEqual([child['control_ids'] for child in result['comparisons']], [['C', 'D']] * 2)
        self.assertEqual(workflow.parse_design(json.dumps(result), packet), result)
        reordered = copy.deepcopy(raw)
        reordered['comparisons'][0]['comparison_name'] = 'renamed'
        reordered['comparisons'][0]['case_ids'].reverse()
        changed = workflow.parse_design(json.dumps(reordered), packet)
        self.assertEqual(changed['comparisons'][0]['comparison_id'], result['comparisons'][0]['comparison_id'])
        self.assertNotEqual(workflow.comparison_identity('FIXTURE', ['A', 'B'], ['C', 'D']),
                            workflow.comparison_identity('FIXTURE', ['C', 'D'], ['A', 'B']))

    def test_design_set_rejects_duplicate_membership_identity_and_lost_samples(self):
        packet, raw = self.design_set_fixture()
        for mutation in ('duplicate_name', 'duplicate_groups', 'lost_sample', 'overlap', 'identity', 'unassigned'):
            with self.subTest(mutation=mutation):
                item = copy.deepcopy(raw)
                if mutation == 'duplicate_name': item['comparisons'][1]['comparison_name'] = ' DRUG '
                if mutation == 'duplicate_groups':
                    item['comparisons'][1] = {**item['comparisons'][0], 'comparison_name': 'other name'}
                if mutation == 'lost_sample': item['comparisons'][0]['excluded_ids'] = ['E']
                if mutation == 'overlap': item['comparisons'][0]['case_ids'].append('C')
                if mutation == 'identity': item['comparisons'][0]['comparison_id'] = 'forged'
                if mutation == 'unassigned': item['unassigned_ids'] = ['A']
                with self.assertRaises(ValueError): workflow.parse_design(json.dumps(item), packet)

    def test_typed_numerical_allegation_requires_visible_computed_result(self):
        packet, raw = self.design_set_fixture()
        primary = {'rows': [{'clock': 'clock', 'case_ids': ['A', 'B'], 'control_ids': ['C', 'D'],
                            'status': 'RESULT', 'statistical_test': 'Welch', 'p_value': .1,
                            'mean_difference': 1., 'effect_size': .5}]}
        result_id = workflow.producer_states(primary, packet)[0]['result_id']
        raw['findings'] = [{'category': 'NUMERICAL_RESULT_ERROR', 'statement': 'Software-test allegation',
                            'evidence': [], 'affected_result': result_id}]
        with self.assertRaisesRegex(ValueError, 'unavailable'): workflow.parse_design(json.dumps(raw), packet)
        checked = workflow.parse_design(json.dumps(raw), packet, primary=primary)
        self.assertTrue(workflow.finding_alarm(checked))
        primary['rows'][0].update(status='CANNOT VERIFY', p_value='', mean_difference='', effect_size='')
        with self.assertRaisesRegex(ValueError, 'COMPUTED'): workflow.parse_design(json.dumps(raw), packet, primary=primary)

    def test_observations_and_missing_evidence_are_not_material_error_flags(self):
        packet, raw = self.design_set_fixture()
        raw['findings'] = [{'category': 'OBSERVATION', 'statement': 'Controls are shared',
                            'evidence': [], 'affected_result': None}]
        self.assertFalse(workflow.finding_alarm(workflow.parse_design(json.dumps(raw), packet)))
        raw['findings'][0].update(category='SOURCE_AMBIGUITY', statement='Unresolved source interpretation')
        with self.assertRaisesRegex(ValueError, 'quotation'): workflow.parse_design(json.dumps(raw), packet)
        raw['findings'][0]['evidence'] = [raw['comparisons'][0]['evidence'][0]]
        checked = workflow.parse_design(json.dumps(raw), packet)
        self.assertFalse(workflow.finding_alarm(checked))
        self.assertEqual(len(workflow.blocking_findings(checked)), 1)
        raw['findings'][0]['category'] = 'SOURCE_DESIGN_CONCERN'
        self.assertTrue(workflow.finding_alarm(workflow.parse_design(json.dumps(raw), packet)))

    def test_comparison_binding_requires_source_labels_and_complete_directional_membership(self):
        result=verify_requested_comparison(self.packet,['A','B'],['C','D'])
        self.assertEqual(result['status'],'PASS')
        self.assertEqual(result['bindings'][0]['source_label'],'treatment')
        for case,control in [(['C','D'],['A','B']),(['A'],['C','D']),(['A','B','C'],['D'])]:
            with self.subTest(case=case,control=control):
                self.assertEqual(verify_requested_comparison(self.packet,case,control)['status'],'FAIL')

    def test_grouping_witness_does_not_relabel_prose_or_supply_reference_truth(self):
        packet = {**self.packet, 'target': 'Drug-treated specimens versus untreated specimens'}
        original = copy.deepcopy(packet)
        value = source_grouping_diagnostics(packet, ['A','B'], ['C','D'])
        self.assertEqual(value['target_binding']['status'], 'CANNOT VERIFY')
        self.assertEqual(value['scientific_status'], 'CANNOT VERIFY')
        self.assertFalse(value['reference_judgment'])
        self.assertEqual(len(value['source_membership_witnesses']), 1)
        witness = value['source_membership_witnesses'][0]
        self.assertEqual(witness['proposed_target'], 'treatment: drug versus treatment: none')
        self.assertEqual(witness['source_binding']['status'], 'PASS')
        self.assertEqual(witness['case_ids'], ['A','B'])
        self.assertEqual(witness['source_metadata_sha256'], packet['source_sha256'])
        self.assertFalse(witness['applied'])
        self.assertFalse(witness['reference_judgment'])
        self.assertEqual(witness['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(packet, original)

    def test_conditional_witness_cannot_hide_reversed_target_direction(self):
        value = source_grouping_diagnostics(self.packet, ['C','D'], ['A','B'])
        self.assertEqual(value['target_binding']['status'], 'FAIL')
        self.assertEqual(value['source_membership_witnesses'][0]['proposed_target'],
                         'treatment: none versus treatment: drug')
        self.assertFalse(value['source_membership_witnesses'][0]['applied'])

    def test_incomplete_or_invalid_membership_has_no_source_witness(self):
        for case, control in [(['A'],['C','D']), (['A','A'],['C','D']),
                              (['A','UNKNOWN'],['C','D']), ([],['C','D'])]:
            with self.subTest(case=case, control=control):
                value = source_grouping_diagnostics(self.packet, case, control)
                self.assertEqual(value['source_membership_witnesses'], [])
                self.assertEqual(value['scientific_status'], 'CANNOT VERIFY')

    def test_conflicting_source_declarations_have_no_membership_witness(self):
        packet = copy.deepcopy(self.packet)
        packet['samples'][0]['characteristics_ch1.2'] = 'treatment: other'
        self.assertEqual(source_grouping_diagnostics(packet, ['A','B'], ['C','D'])['source_membership_witnesses'], [])

    def test_all_source_witnesses_are_preserved_without_choosing_a_causal_label(self):
        packet = copy.deepcopy(self.packet)
        for row in packet['samples']:
            row['characteristics_ch1.2'] = 'batch: ' + ('one' if row['sample_id'] in ['A','B'] else 'two')
        value = source_grouping_diagnostics(packet, ['A','B'], ['C','D'])
        self.assertEqual({w['source_binding']['bindings'][0]['source_label']
                          for w in value['source_membership_witnesses']}, {'treatment','batch'})
        self.assertTrue(all(not w['applied'] for w in value['source_membership_witnesses']))

    def test_literal_source_prompt_is_opt_in_and_does_not_change_contract_or_source(self):
        packet, raw = self.domain_design_fixture()
        packet['comparison_scope'] = 'source_defined_interventions'
        old = workflow.review_messages(packet, review_profile='source_domains_v4')
        default = workflow.review_messages(packet, review_profile='source_domains_v4', source_binding_profile='legacy')
        new = workflow.review_messages(packet, review_profile='source_domains_v4', source_binding_profile='literal_source_v1')
        self.assertEqual([m.content for m in old], [m.content for m in default])
        self.assertEqual(old[1].content, new[1].content)
        self.assertIn('at the same age/context', old[0].content)
        self.assertNotIn('at the same age/context', new[0].content)
        self.assertIn(workflow.SOURCE_BINDING_INSTRUCTIONS, new[0].content)
        self.assertIn('evidence, audit_domains, covariate_plan.', new[0].content)
        with self.assertRaisesRegex(ValueError, 'Unknown source binding'):
            workflow.review_messages(packet, source_binding_profile='infer_controls')
        packet['comparison_scope'] = 'exact_target'
        with self.assertRaisesRegex(ValueError, 'requires source_defined'):
            workflow.review_messages(packet, source_binding_profile='literal_source_v1')

    def test_literal_source_assessment_preserves_raw_prose_and_emits_unapplied_witnesses(self):
        from src.pipeline_2026_09_10 import assess_design
        from langchain_core.messages import AIMessage
        packet, raw = self.domain_design_fixture()
        packet['comparison_scope'] = 'source_defined_interventions'
        for child in raw['comparisons']:
            child['target'] = 'Unbound software-fixture prose'
        client = unittest.mock.Mock()
        client.invoke.return_value = AIMessage(content=json.dumps(raw))
        path = self.root/'literal-assessment'; path.mkdir()
        value = assess_design(packet, client, path, review_profile='source_domains_v4', source_binding_profile='literal_source_v1')
        self.assertEqual(json.loads((path/'model-design.json').read_text()), raw)
        for child in value['comparisons']:
            witness = json.loads((path/child['comparison_id']/'grouping-diagnostics.json').read_text())
            self.assertEqual(witness['target_binding']['status'], 'CANNOT VERIFY')
            self.assertEqual(witness['scientific_status'], 'CANNOT VERIFY')
            self.assertTrue(witness['source_membership_witnesses'])

    def test_literal_source_binding_rejects_older_review_profiles_before_acquisition(self):
        from src.pipeline_2026_09_10 import run
        packet = {**self.packet, 'comparison_scope':'source_defined_interventions'}
        for profile in ('legacy','source_fields_v2','source_domains_v3'):
            with self.subTest(profile=profile):
                with self.assertRaisesRegex(ValueError, 'requires source_domains_v4'):
                    workflow.review_messages(packet, review_profile=profile, source_binding_profile='literal_source_v1')
                # Incomplete arguments deliberately prove this stops before data access.
                with self.assertRaisesRegex(ValueError, 'requires source_domains_v4'):
                    run(SimpleNamespace(comparison_scope='source_defined_interventions',
                        source_binding_profile='literal_source_v1', audit_profile=profile))

    def test_literal_source_producer_keeps_effect_schema_and_default_prompt(self):
        from src.pipeline_2026_09_10 import primary_analysis, EFFECT_DEFINITION
        packet = {**self.packet, 'comparison_scope': 'source_defined_interventions'}
        clocks = self.root/'clocks.csv'; clocks.write_text('sample_id,fixture\nA,1\n')
        prompts = []
        for index, kwargs in enumerate([{}, {'source_binding_profile':'legacy'},
                                       {'source_binding_profile':'literal_source_v1'}]):
            path = self.root/f'producer-{index}'; path.mkdir()
            with patch('src.pipeline_2026_09_10.CodingAgent') as agent:
                agent.return_value.go.side_effect = RuntimeError('Capture prompt before generation')
                with self.assertRaisesRegex(RuntimeError, 'Capture prompt'):
                    primary_analysis(packet, {'path':str(clocks),'clocks':['fixture']}, None, path,
                                     analysis_contract='explicit_effect_v2', **kwargs)
                prompts.append(agent.return_value.go.call_args.args[0])
        self.assertEqual(prompts[0], prompts[1])
        self.assertIn('and matched control at the same age.', prompts[0])
        self.assertNotIn('and matched control at the same age.', prompts[2])
        self.assertIn(workflow.SOURCE_BINDING_INSTRUCTIONS, prompts[2])
        self.assertIn('effect_size,effect_size_definition,status,reason', prompts[2])
        self.assertIn(EFFECT_DEFINITION, prompts[2])

    def test_source_binding_plan_identity_changes_cannot_resume(self):
        from audit_system.report_2026_09_10 import run_source_design_plan
        from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
        plan = {**self.source_design_plan_v2(), 'schema':'clockbase-source-design-run-v4',
                'review_profile':'source_domains_v4', 'comparison_scope':'source_defined_interventions',
                'source_binding_profile':'literal_source_v1'}
        path = self.root/'binding-plan.json'; path.write_text(json.dumps(plan))
        out = self.root/'binding-plan-run'
        with patch('src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke') as call:
            result = run_source_design_plan(path, out, prepare_only=True)
            self.assertEqual(result['generation_calls'], 0)
            request = json.loads((out/'00_source_preflight/attempt-000/prepared-request.json').read_text())
            self.assertIn(workflow.SOURCE_BINDING_INSTRUCTIONS, request['messages'][0]['content'])
            plan['source_binding_profile'] = 'legacy'; path.write_text(json.dumps(plan))
            with self.assertRaises(PipelineIntegrityError):
                run_source_design_plan(path, out, prepare_only=True)
            call.assert_not_called()

    def test_source_condition_matches_observed_groups_but_does_not_supply_reference_truth(self):
        packet, raw = self.design_set_fixture()
        finding = {'category': 'SOURCE_DESIGN_CONCERN', 'statement': 'Pooling distinct fixture treatments needs justification',
                   'affected_result': None, 'condition': {'operation': 'pools_groups', 'arm': 'case',
                   'groups': [['A', 'B'], ['E', 'F']]},
                   'evidence': [entry for child in raw['comparisons'] for entry in child['evidence']]}
        raw['findings'] = [finding]
        design = workflow.parse_design(json.dumps(raw), packet)
        original = {'clock': 'fixture-clock', 'case_ids': ['A', 'B', 'E', 'F'], 'control_ids': ['C', 'D'],
                    'status': 'RESULT', 'statistical_test': 'Welch', 'p_value': .1, 'mean_difference': 1., 'effect_size': .5}
        assessed = workflow.attribute_findings(design, original, packet)[0]
        self.assertTrue(assessed['producer_error_allegation'])
        self.assertEqual(assessed['attribution'], 'OBSERVED_ACTION_ALLEGATION')
        self.assertEqual(assessed['condition_verification']['witness']['observed_ids_by_source_group'], [['A', 'B'], ['E', 'F']])
        self.assertEqual(assessed['material_error_status'], 'CANNOT VERIFY')
        original['case_ids'] = ['A', 'B']
        separate = workflow.attribute_findings(design, original, packet)[0]
        self.assertFalse(separate['producer_error_allegation'])
        self.assertEqual(separate['attribution'], 'ACTION_NOT_PRESENT')
        self.assertEqual(workflow.finding_records(design)[0], finding)

    def test_condition_match_on_abstention_or_missing_row_is_not_a_computed_error(self):
        packet, raw = self.design_set_fixture()
        raw['findings'] = [{'category': 'SOURCE_DESIGN_CONCERN', 'statement': 'Fixture conditional concern',
                            'affected_result': None, 'condition': {'operation': 'uses_any', 'arm': 'case', 'sample_ids': ['A']},
                            'evidence': [raw['comparisons'][0]['evidence'][0]]}]
        design = workflow.parse_design(json.dumps(raw), packet)
        for status in ('CANNOT VERIFY', 'NOT_REPORTED'):
            original = {'clock': 'fixture-clock', 'case_ids': ['A', 'B'], 'control_ids': ['C', 'D'],
                        'status': status, 'p_value': '', 'mean_difference': '', 'effect_size': '', 'absence_record': True}
            item = workflow.attribute_findings(design, original, packet)[0]
            self.assertTrue(item['condition_verification']['matches'])
            self.assertFalse(item['producer_error_allegation'])
            self.assertEqual(item['attribution'], 'NO_COMPUTED_RESULT')
        self.assertTrue(workflow.finding_alarm(design))

    def test_condition_requires_exact_supported_schema_and_source_witnesses(self):
        packet, raw = self.design_set_fixture()
        good = {'category': 'SOURCE_DESIGN_CONCERN', 'statement': 'Fixture concern', 'affected_result': None,
                'condition': {'operation': 'pools_groups', 'arm': 'case', 'groups': [['A'], ['E']]},
                'evidence': [entry for child in raw['comparisons'] for entry in child['evidence']]}
        for change in ('unknown_operation', 'either_pool', 'overlap', 'empty', 'invented', 'unquoted', 'non_concern', 'extra'):
            with self.subTest(change=change):
                item = copy.deepcopy(good)
                if change == 'unknown_operation': item['condition']['operation'] = 'execute_python'
                if change == 'either_pool': item['condition']['arm'] = 'either'
                if change == 'overlap': item['condition']['groups'] = [['A', 'B'], ['B', 'E']]
                if change == 'empty': item['condition']['groups'] = [[], ['E']]
                if change == 'invented': item['condition']['groups'] = [['A'], ['INVENTED']]
                if change == 'unquoted': item['evidence'] = [raw['comparisons'][0]['evidence'][0]]
                if change == 'non_concern': item['category'] = 'OBSERVATION'
                if change == 'extra': item['condition']['code'] = 'untrusted'
                raw['findings'] = [item]
                with self.assertRaises(ValueError): workflow.parse_design(json.dumps(raw), packet)

    def test_only_a_condition_absent_from_proposed_groups_can_stop_blocking_recomputation(self):
        packet, raw = self.design_set_fixture()
        child = raw['comparisons'][0]
        child['findings'] = [{'category': 'SOURCE_DESIGN_CONCERN', 'statement': 'Fixture concern about other treatment',
                             'affected_result': None, 'condition': {'operation': 'uses_any', 'arm': 'case', 'sample_ids': ['E']},
                             'evidence': [entry for entry in raw['comparisons'][1]['evidence'] if entry['sample_id'] == 'E']}]
        result = workflow.parse_design(json.dumps(raw), packet)
        self.assertEqual(result['comparisons'][0]['proposed_action'], 'recompute_welch')
        child['findings'][0]['condition'] = None
        with self.assertRaisesRegex(ValueError, 'Welch proposal'): workflow.parse_design(json.dumps(raw), packet)
        child['findings'][0]['condition'] = {'operation': 'uses_any', 'arm': 'control', 'sample_ids': ['C']}
        child['findings'][0]['evidence'] = [entry for entry in child['evidence'] if entry['sample_id'] == 'C']
        with self.assertRaisesRegex(ValueError, 'Welch proposal'): workflow.parse_design(json.dumps(raw), packet)

    def test_treatment_substrings_unrelated_phenotypes_and_unrequested_strata_cannot_bind(self):
        for target in ['Does drug treat aging?', 'radiation resistant versus parental',
                       'drug versus none in a different tissue', 'drug-like versus none',
                       'Reconstruct study design', 'none versus none']:
            with self.subTest(target=target):
                packet={**self.packet,'target':target}
                self.assertEqual(verify_requested_comparison(packet,['A','B'],['C','D'])['status'],'CANNOT VERIFY')
        packet={**self.packet,'target':'Compare treatment: drug versus treatment: none'}
        self.assertEqual(verify_requested_comparison(packet,['A','B'],['C','D'])['status'],'PASS')

    def test_ambiguous_source_dimensions_and_conflicting_declarations_remain_unresolved(self):
        packet=copy.deepcopy(self.packet)
        for i,sample in enumerate(packet['samples']):
            sample['characteristics_ch1.1']='condition: '+('drug' if i in [0,2] else 'none')
        self.assertEqual(verify_requested_comparison(packet,['A','B'],['C','D'])['status'],'CANNOT VERIFY')
        packet=copy.deepcopy(self.packet)
        packet['samples'][0]['characteristics_ch1.1']='treatment: none'
        self.assertEqual(verify_requested_comparison(packet,['A','B'],['C','D'])['status'],'CANNOT VERIFY')

    def test_metadata_summaries_preserve_units_missingness_and_conflicts(self):
        packet=copy.deepcopy(self.packet)
        for sample,value in zip(packet['samples'],['50 years','50 years','10 months','']):
            sample['characteristics_ch1.1']='age: '+value
        before=copy.deepcopy(packet)
        result=source_metadata_summary(packet,detailed=True)['fields']['age']
        self.assertEqual(result['observed'],3);self.assertEqual(result['declared_missing'],1)
        self.assertEqual(result['numeric_by_literal_unit'],[
            {'unit':'months','n':1,'minimum':10.,'median':10.,'maximum':10.},
            {'unit':'years','n':2,'minimum':50.,'median':50.,'maximum':50.}])
        self.assertEqual(packet,before)
        packet['samples'][0]['characteristics_ch1.2']='age: 60 years'
        result=source_metadata_summary(packet)['fields']['age']
        self.assertEqual(result['conflicting'],1);self.assertEqual(result['observed'],2)
        with self.assertRaises(ValueError):source_metadata_summary(packet,['A','A'])
        with self.assertRaises(ValueError):source_metadata_summary(packet,['INVENTED'])

    def test_stratified_source_proposal_preserves_original_target_and_complete_membership(self):
        packet=copy.deepcopy(self.packet)
        packet['target']='An unsupported phenotype label'
        for sample in packet['samples']:
            sample['characteristics_ch1.1']='cell line: first line'
        packet['samples'] += [{**packet['samples'][0],'sample_id':'E','characteristics_ch1.1':'cell line: second line'},
                              {**packet['samples'][2],'sample_id':'F','characteristics_ch1.1':'cell line: second line'}]
        before=copy.deepcopy(packet)
        proposals=propose_source_comparisons(packet,['A','B'],['C','D'])
        self.assertEqual(len(proposals),1)
        self.assertEqual(proposals[0]['proposed_target'],'treatment: drug versus treatment: none within cell line: first line')
        self.assertEqual(proposals[0]['source_binding']['status'],'PASS')
        self.assertFalse(proposals[0]['applied'])
        self.assertEqual(proposals[0]['scientific_status'],'CANNOT VERIFY')
        self.assertEqual(packet,before)
        self.assertEqual(verify_requested_comparison(packet,['A','B'],['C','D'])['status'],'CANNOT VERIFY')

    def test_missing_stratum_information_cannot_silently_exclude_an_eligible_sample(self):
        packet=copy.deepcopy(self.packet)
        packet['target']='drug versus none within cell line: first line'
        for sample in packet['samples']:
            sample['characteristics_ch1.1']='cell line: first line'
        packet['samples'].append({**packet['samples'][0],'sample_id':'E'})
        del packet['samples'][-1]['characteristics_ch1.1']
        self.assertEqual(verify_requested_comparison(packet,['A','B'],['C','D'])['status'],'CANNOT VERIFY')
        self.assertEqual(propose_source_comparisons(packet,['A','B'],['C','D']),[])

    def test_proposals_cannot_redefine_groups_from_numeric_sample_indices(self):
        packet=copy.deepcopy(self.packet)
        packet['samples'][0]['characteristics_ch1']='treatment: unrelated'
        self.assertEqual(propose_source_comparisons(packet,['A','B'],['C','D']),[])
        self.assertEqual(propose_source_comparisons(packet,['A','A'],['C','D']),[])
        self.assertEqual(propose_source_comparisons(packet,[],[]),[])

    def test_metadata_summary_denominators_follow_source_selection_without_primary_results(self):
        packet={**self.packet,'primary':{'p_value':'PRIMARY_NUMERIC_SENTINEL'},'prior_verdict':'VERDICT_SENTINEL'}
        result=source_metadata_summary(packet,['A','B'])
        self.assertEqual(result['assay_samples'],2)
        self.assertEqual(result['fields']['treatment']['counts'],{'drug':2})
        self.assertNotIn('SENTINEL',json.dumps(result))
        self.assertEqual(source_metadata_summary(packet,[])['assay_samples'],0)

    def test_producer_states_separate_abstention_estimates_and_unknown_release(self):
        row={'clock':'fixture','status':'CANNOT VERIFY','case_ids':['A'],'control_ids':['B'],
             'p_value':'','mean_difference':'','effect_size':''}
        state=workflow.primary_row_state(row)
        self.assertEqual(state['state'],'ABSTAINED')
        self.assertEqual(state['historical_release_status'],'CANNOT VERIFY')
        self.assertEqual(workflow.primary_row_state({**row,'p_value':.1})['state'],'INCONSISTENT')
        computed={**row,'status':'RESULT','p_value':.1,'mean_difference':1.,'effect_size':.2,'statistical_test':'Welch'}
        self.assertEqual(workflow.primary_row_state(computed)['state'],'COMPUTED')
        for bad in ({'p_value':float('nan')},{'p_value':2},{'effect_size':''},{'statistical_test':''}):
            self.assertEqual(workflow.primary_row_state({**computed,**bad})['state'],'INCONSISTENT')

    def test_source_summaries_match_review_arms_but_producer_states_are_self_review_only(self):
        primary={'rows':[{'clock':'fixture','status':'CANNOT VERIFY','p_value':'','effect_size':'','mean_difference':''}]}
        independent=workflow.review_messages(self.packet)
        control=workflow.review_messages(self.packet,mode='self_review',primary=primary)
        own=json.loads(independent[1].content);other=json.loads(control[1].content)
        self.assertEqual(own['metadata_summary'],other['source_packet']['metadata_summary'])
        self.assertNotIn('row_states',own)
        self.assertEqual(other['preceding_analysis']['row_states'][0]['state'],'ABSTAINED')

    def test_self_review_access_does_not_leak_into_independent_design(self):
        primary={'rows':[{'reason':'PRECEDING_ANALYSIS_SENTINEL','p_value':.123}]}
        independent=workflow.review_messages(self.packet)
        control=workflow.review_messages(self.packet,mode='self_review',primary=primary)
        self.assertNotIn('PRECEDING_ANALYSIS_SENTINEL',independent[1].content)
        self.assertIn('PRECEDING_ANALYSIS_SENTINEL',control[1].content)
        self.assertEqual(json.loads(control[1].content)['source_packet'],json.loads(independent[1].content))
        self.assertNotIn('judgments are withheld',control[0].content)
        with self.assertRaises(ValueError):workflow.review_messages(self.packet,primary=primary)
        with self.assertRaises(ValueError):workflow.review_messages(self.packet,mode='self_review')

    def test_source_compaction_preserves_sparse_fields_empty_strings_unicode_and_context(self):
        packet=copy.deepcopy(self.packet)
        packet['samples'][0]['title']='Prefix 001 λ '
        packet['samples'][1]['title']='Prefix 002 λ '
        packet['samples'][0]['optional']=''
        packet['samples'][1]['optional']='  kept exactly  '
        packet['study_documents']=[{'source_id':'SOFTWARE_FIXTURE','fields':{'text':'Source quotation remains unchanged'}}]
        encoded=compact_source_packet(packet)
        self.assertEqual(expand_source_packet(json.loads(json.dumps(encoded))),packet)
        self.assertNotIn('samples',encoded)
        self.assertNotIn('optional',expand_source_packet(encoded)['samples'][2])
        self.assertEqual(expand_source_packet(encoded)['samples'][0]['optional'],'')
        self.assertEqual(self.packet['samples'][0]['title'],'case1')

    def test_source_decoder_rejects_changed_values_ids_dimensions_and_unknown_schema(self):
        encoded=compact_source_packet(self.packet)
        for mutation in ('value','identity','drop_row','unknown','boolean_index','drop_field'):
            with self.subTest(mutation=mutation):
                bad=copy.deepcopy(encoded);table=bad['sample_table']
                if mutation=='value':table['columns']['title']={'values':['changed']*4}
                if mutation=='identity':table['sample_ids'][0]='IMAGINED_SAMPLE'
                if mutation=='drop_row':table['sample_ids'].pop()
                if mutation=='unknown':table['extra']='unexpected'
                if mutation=='boolean_index':table['columns']['title']={'dictionary':['x'],'indices':[False]*4}
                if mutation=='drop_field':table['columns'].pop('title')
                with self.assertRaises(ValueError):expand_source_packet(bad)

    def test_compaction_rejects_duplicate_ids_and_nonstring_source_fields(self):
        for mutation in ('duplicate','nonstring','encoded_collision'):
            with self.subTest(mutation=mutation):
                packet=copy.deepcopy(self.packet)
                if mutation=='duplicate':packet['samples'][1]['sample_id']='A'
                if mutation=='nonstring':packet['samples'][0]['title']=123
                if mutation=='encoded_collision':packet['sample_table']={}
                with self.assertRaises(ValueError):compact_source_packet(packet)

    def test_compact_review_keeps_same_source_for_both_information_conditions(self):
        primary={'rows':[{'reason':'PRECEDING_ANALYSIS_SENTINEL'}]}
        independent=workflow.review_messages(self.packet,representation='compact')
        control=workflow.review_messages(self.packet,mode='self_review',primary=primary,representation='compact')
        source=json.loads(independent[1].content)
        self.assertEqual(source,json.loads(control[1].content)['source_packet'])
        self.assertEqual(expand_source_packet(source),self.packet)
        self.assertNotIn('PRECEDING_ANALYSIS_SENTINEL',independent[1].content)
        self.assertIn('PRECEDING_ANALYSIS_SENTINEL',control[1].content)
        self.assertIn('null direct values or index',independent[0].content)
        self.assertEqual(workflow.review_messages(self.packet)[0].content,workflow.SYSTEM)
        self.assertEqual(workflow.review_messages(self.packet)[1].content,json.dumps(self.packet,ensure_ascii=False))

    def test_oversized_unique_source_is_still_rejected_before_model_transport(self):
        from src.pipeline_runtime_2026_09_10 import BoundedOpenAI
        packet=copy.deepcopy(self.packet)
        for i,row in enumerate(packet['samples']):row['title']=str(i)+('x'*50000)
        client=BoundedOpenAI(self.root/'api','software_fixture',prompt_bytes=48000,max_calls=1)
        with patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport:
            with self.assertRaisesRegex(ValueError,'nothing truncated or sent'):
                client.invoke(workflow.review_messages(packet))
        transport.assert_not_called()
        self.assertFalse(list((self.root/'api').glob('*/reservation.json')))

    def selector_design(self):
        design=copy.deepcopy(self.design)
        select=lambda value:{'where':[{'field':'characteristics_ch1','equals':value}]}
        design['case_ids']=select('treatment: drug')
        design['control_ids']=select('treatment: none')
        design['excluded_ids']={'remaining':True}
        design['evidence']=[{'sample_selector':select(value),'field':'characteristics_ch1','quote':value}
            for value in ['treatment: drug','treatment: none']]
        return design

    def test_source_selectors_expand_complete_membership_and_check_each_quote(self):
        raw=self.selector_design();result=workflow.parse_design(json.dumps(raw),self.packet)
        self.assertEqual(result['case_ids'],['A','B']);self.assertEqual(result['control_ids'],['C','D'])
        self.assertEqual(result['excluded_ids'],[])
        self.assertEqual([q['sample_id'] for q in result['evidence']],['A','B','C','D'])
        self.assertIn('where',raw['case_ids'])
        raw['evidence'][0].update(field='title',quote='case1')
        with self.assertRaisesRegex(ValueError,'Unsupported source quotation'):
            workflow.parse_design(json.dumps(raw),self.packet)

    def test_selectors_reject_unknown_fuzzy_overlapping_or_incomplete_assignments(self):
        for mutation in ('unknown_field','substring','numeric','unknown_operator','empty_conditions','overlap','lost_sample','numeric_remainder'):
            with self.subTest(mutation=mutation):
                raw=self.selector_design();condition=raw['case_ids']['where'][0]
                if mutation=='unknown_field':condition['field']='imagined'
                if mutation=='substring':condition['equals']='drug'
                if mutation=='numeric':condition['equals']=1
                if mutation=='unknown_operator':condition['contains']=condition.pop('equals')
                if mutation=='empty_conditions':raw['case_ids']['where']=[]
                if mutation=='overlap':raw['control_ids']=copy.deepcopy(raw['case_ids'])
                if mutation=='lost_sample':raw['case_ids']=['A'];raw['excluded_ids']=[]
                if mutation=='numeric_remainder':raw['excluded_ids']={'remaining':1}
                with self.assertRaises(ValueError):workflow.parse_design(json.dumps(raw),self.packet)

    def test_unresolved_selector_design_retains_unresolved_scientific_assumptions(self):
        raw=self.selector_design();raw.update(status='CANNOT VERIFY',pairing='unknown',
            unresolved=['Biological unit identity is missing'],proposed_action='obtain_source_evidence')
        result=workflow.parse_design(json.dumps(raw),self.packet)
        self.assertEqual(result['status'],'CANNOT VERIFY');self.assertEqual(result['case_ids'],['A','B'])
        raw['proposed_action']='recompute_welch'
        with self.assertRaisesRegex(ValueError,'Welch proposal conflicts'):
            workflow.parse_design(json.dumps(raw),self.packet)

    def source_design_plan(self):
        return {'schema':'clockbase-source-design-run-v1','split':'development','provider':'openai',
            'model':'gpt-5.6-sol','study_id':'FIXTURE','metadata':str(self.path),'target':'drug versus none',
            'source_sha256':{str(self.path):workflow.sha(self.path)},'max_calls':1,'max_cost':.5,
            'output_tokens':3500,'prompt_bytes':48000,'input_per_million':5.,'output_per_million':20.}

    def test_source_design_preflight_pins_sources_and_calls_no_model(self):
        from audit_system.report_2026_09_10 import run_source_design_plan
        plan=self.root/'source-review-plan.json';plan.write_text(json.dumps(self.source_design_plan()))
        with patch('src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke') as invoke:
            result=run_source_design_plan(plan,self.root/'review',prepare_only=True)
            self.assertEqual(result['generation_calls'],0);self.assertEqual(result['samples'],4)
            self.assertEqual(run_source_design_plan(plan,self.root/'review',prepare_only=True),result)
        invoke.assert_not_called()
        self.path.write_text(self.path.read_text().replace('case1','changed'))
        with self.assertRaisesRegex(ValueError,'Frozen design source changed'):
            run_source_design_plan(plan,self.root/'review',prepare_only=True)

    def test_source_design_rejects_invalid_budgets_before_transport_or_output(self):
        from audit_system.report_2026_09_10 import run_source_design_plan
        for changes in ({'max_cost':.00001},{'max_cost':float('nan')},{'max_calls':2},
                        {'prompt_bytes':1},{'output_tokens':True},{'provider':'unapproved'}):
            with self.subTest(changes=changes):
                plan=self.root/'bad-plan.json';plan.write_text(json.dumps({**self.source_design_plan(),**changes}))
                with patch('src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke') as invoke:
                    with self.assertRaises(ValueError):run_source_design_plan(plan,self.root/'not-created')
                invoke.assert_not_called();self.assertFalse((self.root/'not-created').exists())

    def source_design_plan_v2(self):
        return {**self.source_design_plan(), 'schema':'clockbase-source-design-run-v2',
                'model':'gpt-5.6-luna', 'review_profile':'source_fields_v2', 'validation_retries':1,
                'max_calls':2, 'max_cost':.1, 'input_per_million':.5, 'output_per_million':2.4}

    def fixture_http_response(self, design=None, *, model='gpt-5.6-luna', finish='stop'):
        """HTTP software fixture, never an empirical or paid model response."""
        body = {'model':model, 'usage':{'prompt_tokens':300, 'completion_tokens':100},
                'choices':[{'finish_reason':finish, 'message':{'content':json.dumps(design or self.design)}}]}
        return SimpleNamespace(returncode=0, stdout=(json.dumps(body)+'\n200').encode())

    def test_source_fields_view_separates_provenance_and_preserves_original_and_legacy(self):
        path, _ = self.make_context()
        packet = workflow.source_packet(self.path,'FIXTURE','drug versus none',context=path)
        for document in packet['study_documents']:
            document['cohort_correspondence'] = 'INTERNAL_VERDICT_SENTINEL'
        original = copy.deepcopy(packet)
        old = workflow.review_messages(packet, review_profile='legacy')
        self.assertEqual([m.content for m in old], [m.content for m in workflow.review_messages(packet)])
        messages = workflow.review_messages(packet, review_profile='source_fields_v2')
        self.assertNotIn('Numerical outcomes and historical model\njudgments are withheld.', messages[0].content)
        self.assertIn('Published prose may contain findings.', messages[0].content)
        view = json.loads(messages[1].content)
        self.assertNotIn('INTERNAL_VERDICT_SENTINEL', messages[1].content)
        self.assertIn('INTERNAL_VERDICT_SENTINEL', old[1].content)
        self.assertEqual(packet, original)
        for before, after, provenance in zip(original['study_documents'], view['study_documents'],
                                             view['source_provenance']['documents']):
            self.assertEqual(set(after), {'source_id','fields'})
            self.assertEqual(after['fields'], before['fields'])
            self.assertEqual(after['source_id'], provenance['source_id'])
            self.assertEqual(provenance['source_sha256'], before['source_sha256'])
            self.assertNotIn('fields', provenance)
        self.assertEqual(view['source_provenance']['acquisition'], original.get('source_acquisition',[]))
        with self.assertRaises(ValueError):workflow.review_messages(packet, review_profile='unknown')

    def test_validation_feedback_reports_all_bad_citations_without_rewriting(self):
        raw = copy.deepcopy(self.design)
        raw['evidence'][0]['quote'] = 'case 1'
        raw['evidence'][1] = {'source_id':'B','field':'title','quote':'case2'}
        original = copy.deepcopy(raw)
        feedback = workflow.design_validation_feedback(json.dumps(raw), self.packet, ValueError('Invalid citations'))
        self.assertEqual(feedback['citation_issue_count'], 2)
        self.assertEqual(raw, original)
        self.assertEqual({item['path'] for item in feedback['citation_issues']}, {'$.evidence[0]','$.evidence[1]'})
        with self.assertRaises(ValueError):workflow.parse_design(json.dumps(raw),self.packet)

    def test_feedback_reports_covariate_coverage_hidden_by_an_earlier_bad_quote(self):
        packet, raw = self.domain_design_fixture()
        raw['audit_domains']['sample_identity']['evidence'] = [
            {'sample_id':'A', 'field':'sample_id', 'quote':'A'}]
        raw['covariate_plan'][0].update(source_fields=['title', 'characteristics_ch1'],
            evidence=[{'sample_id':'A', 'field':'characteristics_ch1', 'quote':'treatment: drug'}])
        original, source = copy.deepcopy(raw), copy.deepcopy(packet)
        feedback = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Unsupported source quotation'))
        self.assertEqual(feedback['citation_issue_count'], 1)
        self.assertEqual(feedback['coverage_issues'], [{
            'path':'$.covariate_plan[0].source_fields', 'kind':'COVARIATE_FIELD_QUOTATION_MISSING',
            'missing_source_fields':['title'], 'missing_count':1, 'omitted_count':0}])
        self.assertEqual((raw, packet), (original, source))
        with self.assertRaisesRegex(ValueError, 'Unsupported source quotation'):
            workflow.parse_design(json.dumps(raw), packet, review_profile='source_domains_v4')
        raw['audit_domains']['sample_identity']['evidence'] = [{'sample_id':'A','field':'title','quote':'case1'}]
        with self.assertRaisesRegex(ValueError, 'Covariate source fields'):
            workflow.parse_design(json.dumps(raw), packet, review_profile='source_domains_v4')

    def test_covariate_coverage_requires_its_own_valid_sample_quote(self):
        packet, raw = self.domain_design_fixture()
        packet['study_documents'] = [{'source_id':'DOC', 'fields':{'title':'Source document fixture'}}]
        raw['covariate_plan'][0]['source_fields'] = ['title']
        for evidence in ([], [{'source_id':'DOC','field':'title','quote':'Source document fixture'}],
                         [{'sample_id':'A','field':'title','quote':'invented quote'}]):
            with self.subTest(evidence=evidence):
                raw['covariate_plan'][0]['evidence'] = evidence
                feedback = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Coverage'))
                self.assertEqual(feedback['coverage_issues'][0]['missing_source_fields'], ['title'])
        raw['covariate_plan'][0].update(source_fields=[],
            evidence=[{'source_id':'DOC','field':'title','quote':'Source document fixture'}])
        self.assertEqual(workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Other'))['coverage_issue_count'], 0)

    def test_condition_coverage_checks_both_operations_without_borrowing_comparison_quotes(self):
        packet, raw = self.domain_design_fixture()
        finding = {'category':'SOURCE_DESIGN_CONCERN', 'statement':'Software fixture conditional concern',
            'affected_result':None, 'evidence':[{'sample_id':'A','field':'title','quote':'case1'}]}
        for condition in ({'operation':'uses_any','arm':'either','sample_ids':['A','B']},
                          {'operation':'pools_groups','arm':'case','groups':[['A'],['B']]}):
            with self.subTest(condition=condition):
                raw['findings'] = [{**finding, 'condition':condition}]
                feedback = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Other quote first'))
                self.assertEqual(feedback['coverage_issues'][0]['missing_sample_ids'], ['B'])
                with self.assertRaisesRegex(ValueError, 'condition sample lacks'):
                    workflow.parse_design(json.dumps(raw), packet, review_profile='source_domains_v4')

    def test_coverage_expands_selectors_and_never_uses_an_invalid_shared_quote(self):
        packet, raw = self.domain_design_fixture()
        evidence = [{'sample_selector':{'where':[{'field':'characteristics_ch1','equals':'treatment: drug'}]},
                     'field':'characteristics_ch1', 'quote':'treatment: drug'}]
        raw['covariate_plan'][0].update(source_fields=['characteristics_ch1'], evidence=copy.deepcopy(evidence))
        raw['findings'] = [{'category':'SOURCE_DESIGN_CONCERN', 'statement':'Software fixture',
            'affected_result':None, 'condition':{'operation':'uses_any','arm':'either','sample_ids':['A','B']},
            'evidence':copy.deepcopy(evidence)}]
        feedback = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Other'))
        self.assertEqual((feedback['citation_issue_count'], feedback['coverage_issue_count']), (0,0))
        raw['findings'][0]['evidence'][0].update(field='title', quote='case1')
        feedback = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Quotation'))
        self.assertEqual(feedback['citation_issue_count'], 1)
        self.assertEqual(feedback['coverage_issues'][0]['missing_sample_ids'], ['A','B'])

    def test_feedback_bounds_coverage_details_and_handles_malformed_responses(self):
        packet, raw = self.domain_design_fixture()
        item = raw['covariate_plan'][0]
        raw['covariate_plan'] = [{**copy.deepcopy(item), 'name':str(i),
                                 'source_fields':['field'+str(n) for n in range(25)]} for i in range(25)]
        feedback = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Unknown source field'))
        self.assertEqual(feedback['coverage_issue_count'], 25)
        self.assertEqual(len(feedback['coverage_issues']), 20)
        self.assertEqual(feedback['coverage_issues_omitted'], 5)
        self.assertEqual(feedback['coverage_issues'][0]['missing_count'], 25)
        self.assertEqual(feedback['coverage_issues'][0]['omitted_count'], 5)
        self.assertEqual(len(feedback['coverage_issues'][0]['missing_source_fields']), 20)
        for text in ('not json', 'null', '[]', '{"covariate_plan":[null,{"source_fields":[{}]}]}',
                     '{"findings":[null,{"condition":{"operation":"uses_any","sample_ids":3}}]}'):
            with self.subTest(text=text):
                self.assertEqual(workflow.design_validation_feedback(text, packet, ValueError('Invalid'))['coverage_issue_count'], 0)

    def test_feedback_schema_selection_preserves_legacy_and_rejects_unknown_versions(self):
        packet, raw = self.domain_design_fixture()
        raw['covariate_plan'][0]['source_fields'] = ['title']
        legacy = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Coverage'),
            schema='clockbase-design-validation-feedback-v1')
        self.assertEqual(set(legacy), {'schema','validation_status','error','citation_issue_count',
                                     'citation_issues','citation_issues_omitted','instruction'})
        current = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Coverage'))
        self.assertEqual(current['coverage_issues'][0]['missing_source_fields'], ['title'])
        self.assertEqual(current['citation_issues'], legacy['citation_issues'])
        for schema in (None, 'clockbase-design-validation-feedback-v0', 'anything', 2):
            with self.subTest(schema=schema), self.assertRaisesRegex(ValueError, 'feedback schema'):
                workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Coverage'), schema=schema)

    def test_first_retry_receives_all_coverage_issues_and_a_second_failure_still_stops(self):
        from langchain_core.messages import AIMessage
        from src.pipeline_2026_09_10 import assess_design
        packet, raw = self.domain_design_fixture()
        raw['audit_domains']['sample_identity']['evidence'] = [{'sample_id':'A','field':'sample_id','quote':'A'}]
        raw['covariate_plan'][0]['source_fields'] = ['title']
        response = json.dumps(raw)
        calls = []
        def invoke(messages):
            calls.append(copy.deepcopy(messages))
            return AIMessage(content=response)
        directory = self.root/'bounded-feedback'; directory.mkdir()
        with self.assertRaisesRegex(ValueError, 'Unsupported source quotation'):
            assess_design(packet, SimpleNamespace(invoke=invoke), directory,
                          review_profile='source_domains_v4', validation_retries=1)
        self.assertEqual(len(calls), 2)
        self.assertNotIn(workflow.SOURCE_COVERAGE_INSTRUCTIONS, calls[0][0].content)
        feedback = json.loads(calls[1][-1].content)
        self.assertEqual((feedback['citation_issue_count'], feedback['coverage_issue_count']), (1,1))
        self.assertIn(workflow.SOURCE_COVERAGE_INSTRUCTIONS, feedback['instruction'])
        self.assertEqual((directory/'response-validation/000-response.txt').read_text(), response)
        self.assertEqual((directory/'response-validation/001-response.txt').read_text(), response)
        self.assertFalse((directory/'design.json').exists())

    def test_v3_feedback_enumerates_missing_samples_in_both_comparison_arms(self):
        packet, raw = self.design_set_fixture()
        raw['comparisons'][0]['evidence'] = [entry for entry in raw['comparisons'][0]['evidence']
                                            if entry['sample_id'] in ('A', 'C')]
        original, source = copy.deepcopy(raw), copy.deepcopy(packet)
        feedback = workflow.design_validation_feedback(json.dumps(raw), packet,
            ValueError('Assigned samples lack source quotations'), schema='clockbase-design-validation-feedback-v3')
        self.assertEqual(feedback['coverage_issues'], [
            {'path':'$.comparisons[0].case_ids', 'kind':'ASSIGNED_SAMPLE_QUOTATION_MISSING',
             'missing_sample_ids':['B'], 'missing_count':1, 'omitted_count':0},
            {'path':'$.comparisons[0].control_ids', 'kind':'ASSIGNED_SAMPLE_QUOTATION_MISSING',
             'missing_sample_ids':['D'], 'missing_count':1, 'omitted_count':0}])
        self.assertEqual(feedback['citation_issue_count'], 0)
        self.assertEqual((raw, packet), (original, source))
        self.assertIn("comparison's own evidence list", feedback['instruction'])
        self.assertNotIn('case2', json.dumps(feedback))
        with self.assertRaisesRegex(ValueError, 'Assigned samples lack source quotations'):
            workflow.parse_design(json.dumps(raw), packet)

    def test_v3_feedback_valid_coverage_and_single_arm_omissions(self):
        for missing in (set(), {'A', 'B'}, {'C', 'D'}):
            with self.subTest(missing=missing):
                raw = copy.deepcopy(self.design)
                raw['evidence'] = [entry for entry in raw['evidence'] if entry['sample_id'] not in missing]
                feedback = workflow.design_validation_feedback(json.dumps(raw), self.packet, ValueError('Fixture'),
                    schema='clockbase-design-validation-feedback-v3')
                observed = {sid for item in feedback['coverage_issues'] for sid in item['missing_sample_ids']}
                self.assertEqual(observed, missing)
                if missing:
                    with self.assertRaisesRegex(ValueError, 'Assigned samples lack source quotations'):
                        workflow.parse_design(json.dumps(raw), self.packet)
                else:
                    self.assertEqual(workflow.parse_design(json.dumps(raw), self.packet)['case_ids'], ['A', 'B'])

    def test_v3_assigned_coverage_requires_local_valid_sample_quotes(self):
        packet, raw = self.design_set_fixture()
        packet['study_documents'] = [{'source_id':'DOC', 'fields':{'title':'Software source document'}}]
        comparison = raw['comparisons'][0]
        complete = copy.deepcopy(comparison['evidence'])
        raw['evidence'] = copy.deepcopy(complete)
        raw['findings'] = [{'category':'OBSERVATION', 'statement':'Software fixture',
                            'evidence':copy.deepcopy(complete), 'affected_result':None}]
        for evidence in ([], [complete[0], complete[0]],
                         [{'sample_id':'E', 'field':'title', 'quote':'second E'}],
                         [{'source_id':'DOC', 'field':'title', 'quote':'Software source document'}],
                         [{'sample_id':'A', 'field':'title', 'quote':'invented quotation'}]):
            with self.subTest(evidence=evidence):
                comparison['evidence'] = copy.deepcopy(evidence)
                feedback = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Fixture'),
                    schema='clockbase-design-validation-feedback-v3')
                gaps = {item['path']: item['missing_sample_ids'] for item in feedback['coverage_issues']}
                self.assertEqual(gaps['$.comparisons[0].control_ids'], ['C', 'D'])
                self.assertEqual(gaps['$.comparisons[0].case_ids'], ['B'] if evidence == [complete[0]] * 2 else ['A', 'B'])
                self.assertTrue(all(path.startswith('$.comparisons[0].') for path in gaps))
                with self.assertRaises(ValueError):
                    workflow.parse_design(json.dumps(raw), packet)

    def test_v3_assigned_coverage_resolves_exact_selectors_without_mutation(self):
        raw = copy.deepcopy(self.design)
        raw['case_ids'] = {'where':[{'field':'characteristics_ch1', 'equals':'treatment: drug'}]}
        raw['control_ids'] = {'where':[{'field':'characteristics_ch1', 'equals':'treatment: none'}]}
        raw['evidence'] = [{'sample_selector':copy.deepcopy(raw[key]), 'field':'characteristics_ch1', 'quote':value}
            for key, value in (('case_ids','treatment: drug'), ('control_ids','treatment: none'))]
        original = copy.deepcopy(raw)
        feedback = workflow.design_validation_feedback(json.dumps(raw), self.packet, ValueError('Fixture'),
            schema='clockbase-design-validation-feedback-v3')
        self.assertEqual(feedback['coverage_issue_count'], 0)
        self.assertEqual(workflow.parse_design(json.dumps(raw), self.packet)['control_ids'], ['C', 'D'])
        raw['evidence'][0].update(field='title', quote='case1')
        feedback = workflow.design_validation_feedback(json.dumps(raw), self.packet, ValueError('Fixture'),
            schema='clockbase-design-validation-feedback-v3')
        self.assertEqual(feedback['citation_issue_count'], 1)
        self.assertEqual(feedback['coverage_issues'][0]['missing_sample_ids'], ['A', 'B'])
        self.assertEqual(raw['case_ids'], original['case_ids'])

    def test_v3_assigned_coverage_is_bounded_and_handles_malformed_arms(self):
        packet = {'samples':[{'sample_id':f'S{i:02d}', 'title':'Software fixture'} for i in range(25)]}
        raw = {'comparisons':[{'case_ids':[row['sample_id'] for row in reversed(packet['samples'])],
                              'control_ids':[], 'evidence':[]} for _ in range(25)]}
        feedback = workflow.design_validation_feedback(json.dumps(raw), packet, ValueError('Fixture'),
            schema='clockbase-design-validation-feedback-v3')
        self.assertEqual((feedback['coverage_issue_count'], len(feedback['coverage_issues']),
                          feedback['coverage_issues_omitted']), (25, 20, 5))
        self.assertEqual(feedback['coverage_issues'][0]['missing_sample_ids'], [f'S{i:02d}' for i in range(20)])
        self.assertEqual((feedback['coverage_issues'][0]['missing_count'],
                          feedback['coverage_issues'][0]['omitted_count']), (25, 5))
        for arm in (None, True, 3, 'A', [None], [{}], [''], {'where':[]}, {'where':[{'field':'absent','equals':'x'}]}):
            with self.subTest(arm=arm):
                malformed = {'comparisons':[None, {'case_ids':arm, 'control_ids':[], 'evidence':[]}]}
                feedback = workflow.design_validation_feedback(json.dumps(malformed), self.packet, ValueError('Original error'),
                    schema='clockbase-design-validation-feedback-v3')
                self.assertEqual(feedback['coverage_issue_count'], 0)
                self.assertEqual(feedback['error'], 'Original error')

    def test_v3_feedback_does_not_weaken_duplicate_or_overlapping_membership_checks(self):
        for case in (['A', 'A', 'B'], ['A', 'B', 'C']):
            with self.subTest(case=case):
                raw = copy.deepcopy(self.design); raw['case_ids'] = case
                feedback = workflow.design_validation_feedback(json.dumps(raw), self.packet, ValueError('Membership'),
                    schema='clockbase-design-validation-feedback-v3')
                self.assertEqual(feedback['coverage_issue_count'], 0)
                with self.assertRaises(ValueError):
                    workflow.parse_design(json.dumps(raw), self.packet)

    def test_v3_assigned_feedback_never_echoes_unowned_member_text(self):
        raw = copy.deepcopy(self.design)
        sentinel = 'SOFTWARE_UNOWNED_RESPONSE_TEXT_SENTINEL' + 'x' * 1000
        raw['case_ids'].append(sentinel)
        raw['evidence'] = [entry for entry in raw['evidence'] if entry['sample_id'] != 'B']
        original = copy.deepcopy(raw)
        feedback = workflow.design_validation_feedback(json.dumps(raw), self.packet, ValueError('Invalid membership'),
            schema='clockbase-design-validation-feedback-v3')
        self.assertEqual(feedback['coverage_issues'][0]['missing_sample_ids'], ['B'])
        self.assertNotIn(sentinel, json.dumps(feedback))
        self.assertEqual(raw, original)
        with self.assertRaises(ValueError):
            workflow.parse_design(json.dumps(raw), self.packet)

    def test_v3_source_feedback_materializes_only_explicit_references_and_keeps_legacy_schemas(self):
        from src.pipeline_2026_09_10 import source_contract_candidate, source_contract_feedback, SOURCE_FIELD_REFERENCE_POLICY
        raw = copy.deepcopy(self.design)
        raw['evidence'] = [{'sample_id':'A', 'field':'title', 'quote_from_field':True},
                           {'sample_id':'C', 'field':'title', 'quote_from_field':True}]
        original = copy.deepcopy(raw)
        text = json.dumps(raw)
        candidate, changes = source_contract_candidate(text, self.packet, SOURCE_FIELD_REFERENCE_POLICY)
        feedback = source_contract_feedback(text, self.packet, ValueError('Assigned samples lack source quotations'),
                                            SOURCE_FIELD_REFERENCE_POLICY)
        self.assertEqual(feedback['schema'], 'clockbase-design-validation-feedback-v3')
        self.assertEqual(feedback['citation_issue_count'], 0)
        self.assertEqual([item['missing_sample_ids'] for item in feedback['coverage_issues']], [['B'], ['D']])
        self.assertEqual([entry['sample_id'] for entry in candidate['evidence']], ['A', 'C'])
        self.assertEqual(len(changes), 2)
        self.assertEqual(raw, original)
        with self.assertRaisesRegex(ValueError, 'Assigned samples lack source quotations'):
            workflow.parse_design(json.dumps(candidate), self.packet)
        for version in ('v1', 'v2'):
            legacy = source_contract_feedback(text, self.packet, ValueError('Assigned samples lack source quotations'),
                SOURCE_FIELD_REFERENCE_POLICY, schema='clockbase-design-validation-feedback-' + version)
            self.assertNotIn('ASSIGNED_SAMPLE_QUOTATION_MISSING', json.dumps(legacy))
            self.assertNotIn("comparison's own evidence list", legacy['instruction'])
            if version == 'v1': self.assertNotIn('coverage_issues', legacy)
            else: self.assertEqual(legacy['coverage_issue_count'], 0)
        raw['evidence'][0]['quote_from_field'] = False
        malformed = source_contract_feedback(json.dumps(raw), self.packet, ValueError('Original materializer error'),
                                              SOURCE_FIELD_REFERENCE_POLICY)
        self.assertEqual(malformed['error'], 'Original materializer error')
        self.assertEqual(malformed['citation_issue_count'], 2)
        self.assertEqual([item['missing_sample_ids'] for item in malformed['coverage_issues']], [['A', 'B'], ['C', 'D']])
        with self.assertRaises(ValueError):
            source_contract_candidate(json.dumps(raw), self.packet, SOURCE_FIELD_REFERENCE_POLICY)

    def test_v3_retry_receives_assigned_sample_gaps_and_second_omission_still_stops(self):
        from langchain_core.messages import AIMessage
        from src.pipeline_2026_09_10 import assess_design, SOURCE_FIELD_REFERENCE_POLICY
        packet, raw = self.domain_design_fixture()
        raw['comparisons'][0]['evidence'] = [
            {'sample_id':sid, 'field':'title', 'quote_from_field':True} for sid in ('A', 'C')]
        response, original = json.dumps(raw), copy.deepcopy(raw)
        calls = []
        def invoke(messages):
            calls.append(copy.deepcopy(messages)); return AIMessage(content=response)
        directory = self.root/'assigned-feedback'; directory.mkdir()
        with self.assertRaisesRegex(ValueError, 'Assigned samples lack source quotations'):
            assess_design(packet, SimpleNamespace(invoke=invoke), directory,
                          validation_retries=1, review_profile='source_domains_v4',
                          source_contract_policy=SOURCE_FIELD_REFERENCE_POLICY)
        self.assertEqual(len(calls), 2)
        feedback = json.loads(calls[1][-1].content)
        self.assertEqual(feedback['schema'], 'clockbase-design-validation-feedback-v3')
        self.assertEqual([item['missing_sample_ids'] for item in feedback['coverage_issues']], [['B'], ['D']])
        for i in (0, 1):
            self.assertEqual((directory/f'response-validation/{i:03d}-response.txt').read_text(), response)
            receipt = json.loads((directory/f'contract-adjustments/{i:03d}-receipt.json').read_text())
            self.assertEqual(receipt['contract_status'], 'REJECTED')
            self.assertFalse(receipt['accepted_result'])
        self.assertFalse((directory/'design.json').exists())
        self.assertEqual(raw, original)

    def domain_design_fixture(self):
        packet, raw = self.design_set_fixture()
        raw['schema'] = workflow.DOMAIN_DESIGN_SCHEMA
        raw['audit_domains'] = {name: {'status':'CANNOT VERIFY', 'summary':'Software fixture: source scope only',
            'evidence':[{'sample_id':'A','field':'characteristics_ch1','quote':'treatment: drug'}],
            'unresolved':['Independent evidence unavailable in software fixture'],
            'proposed_repair':None, 'verification_needed':['Obtain source evidence']}
            for name in workflow.AUDIT_DOMAINS}
        raw['covariate_plan'] = [{'comparison_name':'*', 'name':'age', 'source_fields':[],
            'decision':'unknown', 'rationale':'Not supplied', 'evidence':[],
            'missingness_action':'Retain unavailable; do not impute'}]
        return packet, raw

    def test_domain_review_preserves_groups_and_never_promotes_model_to_reference(self):
        packet, raw = self.domain_design_fixture()
        value = workflow.parse_design(json.dumps(raw), packet, review_profile='source_domains_v3')
        self.assertEqual(value['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(value['comparisons'][0]['case_ids'], ['A','B'])
        check = workflow.verify_domain_review(value, packet)
        self.assertFalse(check['independent_reference'])
        self.assertEqual(check['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(check['repairs_applied'], [])
        self.assertEqual(check['unique_study_sample_ids'], 6)
        self.assertNotIn('PASS', [row['status'] for row in check['domains'].values()])
        messages = workflow.review_messages(packet, review_profile='source_domains_v3')
        self.assertIn(workflow.DOMAIN_DESIGN_SCHEMA, messages[0].content)
        with self.assertRaises(ValueError): workflow.parse_design(json.dumps(raw), packet)

    def test_domain_review_rejects_missing_domains_and_claimed_unseen_measurement_approval(self):
        packet, raw = self.domain_design_fixture()
        for name in workflow.AUDIT_DOMAINS:
            bad = copy.deepcopy(raw); del bad['audit_domains'][name]
            with self.subTest(missing=name), self.assertRaises(ValueError):
                workflow.parse_design(json.dumps(bad), packet, review_profile='source_domains_v3')

    def test_version4_has_one_root_schema_and_preserves_frozen_version3_messages(self):
        packet, raw = self.domain_design_fixture()
        before = workflow.review_messages(packet, review_profile='source_domains_v3')
        messages = workflow.review_messages(packet, review_profile='source_domains_v4')
        self.assertIn('evidence, audit_domains, covariate_plan.', messages[0].content)
        self.assertNotIn('schema must be clockbase-source-design-set-v1', messages[0].content)
        self.assertNotIn('Keep every earlier design-set key', messages[0].content)
        self.assertEqual(before[1].content, messages[1].content)
        self.assertEqual([m.content for m in before], [m.content for m in workflow.review_messages(packet, review_profile='source_domains_v3')])
        self.assertEqual(workflow.parse_design(json.dumps(raw),packet,review_profile='source_domains_v3'),
                         workflow.parse_design(json.dumps(raw),packet,review_profile='source_domains_v4'))
        for name in ('data_provenance','measurement_applicability'):
            bad = copy.deepcopy(raw); bad['audit_domains'][name].update(status='PASS', unresolved=[])
            with self.subTest(unseen=name), self.assertRaisesRegex(ValueError, 'unseen'):
                workflow.parse_design(json.dumps(bad), packet, review_profile='source_domains_v3')

    def test_domain_evidence_and_covariates_cannot_invent_fields_or_supported_actions(self):
        packet, raw = self.domain_design_fixture()
        cases = []
        bad = copy.deepcopy(raw); bad['audit_domains']['grouping']['evidence'][0]['quote'] = 'not a source quotation'; cases.append(bad)
        bad = copy.deepcopy(raw); bad['covariate_plan'][0]['source_fields'] = ['age_invented']; cases.append(bad)
        bad = copy.deepcopy(raw); bad['covariate_plan'][0]['decision'] = 'include'; cases.append(bad)
        bad = copy.deepcopy(raw); bad['covariate_plan'][0]['comparison_name'] = 'absent'; cases.append(bad)
        bad = copy.deepcopy(raw); bad['audit_domains']['grouping']['status'] = 'PASS'; cases.append(bad)
        bad = copy.deepcopy(raw); bad['covariate_plan'][0]['source_fields'] = ['characteristics_ch1']; cases.append(bad)
        for index, bad in enumerate(cases):
            with self.subTest(case=index), self.assertRaises(ValueError):
                workflow.parse_design(json.dumps(bad), packet, review_profile='source_domains_v3')

    def test_version3_plan_prepares_without_model_and_retains_contract_failures(self):
        from audit_system.report_2026_09_10 import run_source_design_plan
        plan = {**self.source_design_plan_v2(), 'schema':'clockbase-source-design-run-v3',
                'review_profile':'source_domains_v3'}
        path = self.root/'domain-plan.json'; path.write_text(json.dumps(plan))
        with patch('src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke') as call:
            result = run_source_design_plan(path, self.root/'domain-run', prepare_only=True)
            call.assert_not_called()
        self.assertEqual(result['generation_calls'], 0)
        self.assertEqual(result['status'], 'PASS')

    def test_completed_response_feedback_is_bounded_and_preserves_every_attempt(self):
        from src.pipeline_2026_09_10 import assess_design
        from langchain_core.messages import AIMessage
        raw = copy.deepcopy(self.design); raw['evidence'][0]['quote'] = 'case 1'
        for sequence, expected_calls, succeeds in [([self.design],1,True),
                ([raw,self.design],2,True), ([raw,raw],2,False)]:
            with self.subTest(sequence=len(sequence), succeeds=succeeds):
                directory = self.root / ('response-'+str(expected_calls)+'-'+str(succeeds)); directory.mkdir()
                client = unittest.mock.Mock()
                client.invoke.side_effect = [AIMessage(content=json.dumps(item)) for item in sequence]
                if succeeds:
                    result = assess_design(self.packet,client,directory,review_profile='source_fields_v2',validation_retries=1)
                    self.assertEqual(result['case_ids'], ['A','B'])
                else:
                    with self.assertRaisesRegex(ValueError,'Unsupported source quotation'):
                        assess_design(self.packet,client,directory,review_profile='source_fields_v2',validation_retries=1)
                self.assertEqual(client.invoke.call_count,expected_calls)
                self.assertEqual(len(list((directory/'response-validation').glob('*-response.txt'))),expected_calls)
                self.assertEqual(json.loads((directory/'response-validation/000-response.txt').read_text()),sequence[0])
                if expected_calls == 2:
                    sent = client.invoke.call_args_list[1].args[0]
                    self.assertEqual([m.type for m in sent],['system','human','ai','human'])
                    self.assertEqual(json.loads(sent[2].content),raw)
                    self.assertEqual(json.loads(sent[3].content)['validation_status'],'FAIL')
                self.assertEqual(json.loads((directory/'source-packet.json').read_text()),self.packet)

    def test_invoke_failures_do_not_trigger_validation_retry(self):
        from src.pipeline_2026_09_10 import assess_design
        from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
        for index, error in enumerate([TimeoutError('uncertain transport'), ValueError('Incomplete response'),
                RuntimeError('Run API call/cost limit reached'), PipelineIntegrityError('Model mismatch')]):
            with self.subTest(error=error):
                directory=self.root/str(index); directory.mkdir()
                client=unittest.mock.Mock();client.invoke.side_effect=error
                with self.assertRaises(type(error)):
                    assess_design(self.packet,client,directory,review_profile='source_fields_v2',validation_retries=1)
                self.assertEqual(client.invoke.call_count,1)
                self.assertFalse((directory/'response-validation').exists())

    def test_versioned_plan_rejects_implicit_retries_and_invalid_options(self):
        from audit_system.report_2026_09_10 import run_source_design_plan
        changes = [{'schema':'clockbase-source-design-run-v1'}, {'review_profile':'legacy'},
                   {'validation_retries':True}, {'validation_retries':2}, {'max_calls':3}, {'max_calls':True}]
        for change in changes:
            with self.subTest(change=change):
                path=self.root/'invalid.json';path.write_text(json.dumps({**self.source_design_plan_v2(),**change}))
                with patch('src.pipeline_runtime_2026_09_10.subprocess.run') as transport:
                    with self.assertRaises(ValueError):run_source_design_plan(path,self.root/'not-created')
                transport.assert_not_called();self.assertFalse((self.root/'not-created').exists())

    def test_version2_transport_enforces_model_identity_and_completion_without_retry(self):
        from audit_system.report_2026_09_10 import run_source_design_plan
        path=self.root/'plan.json';path.write_text(json.dumps(self.source_design_plan_v2()))
        for index, response in enumerate([self.fixture_http_response(model='different-model'),
                                         self.fixture_http_response(finish='length')]):
            with self.subTest(index=index):
                output=self.root/('transport-'+str(index))
                with patch.dict('os.environ',{'OPENAI_API_KEY':'software-fixture-key'}), \
                     patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=response) as transport:
                    with self.assertRaises(ValueError):run_source_design_plan(path,output)
                self.assertEqual(transport.call_count,1)
                self.assertTrue((output/'api/source_design-000/receipt.json').exists())
                self.assertFalse((output/'api/source_design-001').exists())

    def test_completed_version2_rejects_changed_api_artifacts_without_transport(self):
        from audit_system.report_2026_09_10 import run_source_design_plan
        path=self.root/'plan.json';path.write_text(json.dumps(self.source_design_plan_v2()))
        for name in ['request.json','response.json','receipt.json','reservation.json']:
            with self.subTest(name=name):
                output=self.root/name.replace('.json','')
                with patch.dict('os.environ',{'OPENAI_API_KEY':'software-fixture-key'}), \
                     patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=self.fixture_http_response()) as transport:
                    first=run_source_design_plan(path,output)
                    self.assertEqual(run_source_design_plan(path,output),first)
                    self.assertEqual(transport.call_count,1)
                    artifact=output/'api/source_design-000'/name
                    artifact.write_text(artifact.read_text()+'\n')
                    with self.assertRaisesRegex(ValueError,'API artifacts changed'):
                        run_source_design_plan(path,output)
                    self.assertEqual(transport.call_count,1)

    def test_feedback_uses_same_persistent_budget_and_unknown_billing_is_not_repeated(self):
        from audit_system.report_2026_09_10 import run_source_design_plan
        path=self.root/'plan.json';plan=self.source_design_plan_v2();path.write_text(json.dumps(plan))
        prepared=run_source_design_plan(path,self.root/'preflight',prepare_only=True)
        plan['max_cost']=prepared['reserved_cost_upper_estimate_usd'];path.write_text(json.dumps(plan))
        bad=copy.deepcopy(self.design);bad['evidence'][0]['quote']='case 1'
        with patch.dict('os.environ',{'OPENAI_API_KEY':'software-fixture-key'}), \
             patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=self.fixture_http_response(bad)) as transport:
            with self.assertRaisesRegex(RuntimeError,'cost limit'):
                run_source_design_plan(path,self.root/'limited')
            self.assertEqual(transport.call_count,1)
            self.assertFalse((self.root/'limited/api/source_design-001').exists())
        with patch.dict('os.environ',{'OPENAI_API_KEY':'software-fixture-key'}), \
             patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=SimpleNamespace(returncode=1)) as transport:
            with self.assertRaisesRegex(RuntimeError,'transport failed'):
                run_source_design_plan(path,self.root/'uncertain')
            with self.assertRaisesRegex(RuntimeError,'uncertain billing'):
                run_source_design_plan(path,self.root/'uncertain')
            self.assertEqual(transport.call_count,1)
            self.assertFalse((self.root/'uncertain/api/source_design-001').exists())

    def source_design_batch_fixture(self, *, studies=3):
        plans=[]
        for index in range(studies):
            child={**self.source_design_plan_v2(),'study_id':'FIXTURE'+str(index)}
            path=self.root/f'child-{index}.json';path.write_text(json.dumps(child))
            plans.append({'path':str(path),'sha256':workflow.sha(path)})
        batch={'schema':'clockbase-source-design-batch-v1','split':'development', 'provider':'openai',
               'model':'gpt-5.6-luna','concurrency':1,'max_total_cost_usd':.3,'max_total_calls':6,'plans':plans}
        path=self.root/'batch.json';path.write_text(json.dumps(batch))
        return path,batch

    def test_batch_aggregate_budget_and_all_child_sources_checked_before_transport(self):
        from audit_system.report_2026_09_10 import run_source_design_batch
        for mutation in ['cost','calls','last_source','provider','child_hash']:
            with self.subTest(mutation=mutation):
                path,batch=self.source_design_batch_fixture()
                if mutation=='cost':batch['max_total_cost_usd']=.29
                if mutation=='calls':batch['max_total_calls']=5
                if mutation in ('last_source','provider'):
                    child_path=Path(batch['plans'][-1]['path']);child=json.loads(child_path.read_text())
                    if mutation=='last_source':child['source_sha256'][str(self.path)]='0'*64
                    else:child['provider']='different-provider'
                    child_path.write_text(json.dumps(child));batch['plans'][-1]['sha256']=workflow.sha(child_path)
                if mutation=='child_hash':batch['plans'][-1]['sha256']='0'*64
                path.write_text(json.dumps(batch))
                with patch('src.pipeline_runtime_2026_09_10.BoundedOpenAI.invoke') as transport:
                    with self.assertRaises(ValueError):run_source_design_batch(path,self.root/mutation)
                transport.assert_not_called()

    def test_public_trace_exports_attested_source_evidence_and_failed_runs_without_private_fields(self):
        """Software HTTP/source fixtures test packaging, not scientific performance."""
        # Import numerical libraries before replacing the shared subprocess
        # function with an HTTP fixture (NumPy probes CPU features at import).
        from src.pipeline_2026_09_10 import assess_design
        from audit_system.report_2026_09_10 import run_source_design_batch, export_source_audit_trace
        packet, raw = self.domain_design_fixture()
        ids = {old:'GSM'+str(i+1) for i, old in enumerate('ABCDEF')}
        def replace(value):
            if isinstance(value, dict): return {key:replace(item) for key,item in value.items()}
            if isinstance(value, list): return [replace(item) for item in value]
            return ids.get(value,value) if isinstance(value,str) else value
        packet, raw = replace(packet), replace(raw)
        for index, failed in enumerate((False, True)):
            study='GSE'+str(900001+index)
            metadata=self.root/(study+'.csv')
            with metadata.open('w',newline='') as handle:
                rows=[{'GSM':row['sample_id'],**{key:value for key,value in row.items() if key!='sample_id'}} for row in packet['samples']]
                writer=csv.DictWriter(handle,fieldnames=sorted(set().union(*(row.keys() for row in rows))))
                writer.writeheader();writer.writerows(rows)
            plan={**self.source_design_plan_v2(), 'schema':'clockbase-source-design-run-v3',
                  'review_profile':'source_domains_v3','study_id':study,'metadata':str(metadata),
                  'source_sha256':{str(metadata):workflow.sha(metadata)}}
            child=self.root/(study+'-plan.json');child.write_text(json.dumps(plan))
            batch={'schema':'clockbase-source-design-batch-v1','split':'development','provider':'openai',
                   'model':'gpt-5.6-luna','concurrency':1,'max_total_cost_usd':.1,'max_total_calls':2,
                   'plans':[{'path':str(child),'sha256':workflow.sha(child)}]}
            path=self.root/(study+'-batch.json');path.write_text(json.dumps(batch))
            response=copy.deepcopy(raw);response['study_id']=study
            if failed: response['audit_domains']['grouping']['evidence'][0]['quote']='not source evidence'
            run=self.root/(study+'-run')
            with patch.dict('os.environ',{'OPENAI_API_KEY':'software-fixture-key'}), \
                 patch('src.pipeline_runtime_2026_09_10.subprocess.run',return_value=self.fixture_http_response(response)):
                result=run_source_design_batch(path,run)
            self.assertEqual(result['rows'][0]['response_contract_status'],'FAIL' if failed else 'PASS')
            output=self.root/(study+'-public')
            exported=export_source_audit_trace(run/'studies'/study,output)
            self.assertEqual(export_source_audit_trace(run/'studies'/study,output),exported)
            trace=json.loads((output/'audit_trace.json').read_text())
            self.assertEqual(trace['overall']['status'],'CANNOT VERIFY')
            self.assertEqual(len(trace['domains']),7)
            self.assertNotIn(str(self.root),json.dumps(trace))
            self.assertNotIn('software-fixture-key',json.dumps(trace))
            self.assertEqual(trace['events'][1]['status'],'FAIL' if failed else 'PASS')
            self.assertEqual(len(trace['comparisons']),0 if failed else 2)
            (output/'audit_trace.json').write_text('{}')
            with self.assertRaisesRegex(ValueError,'export changed'):
                export_source_audit_trace(run/'studies'/study,output)

    def test_batch_resume_preserves_failed_response_and_enforces_child_artifact_integrity(self):
        from audit_system.report_2026_09_10 import run_source_design_batch
        path,_=self.source_design_batch_fixture(studies=2)
        good=copy.deepcopy(self.design);good['study_id']='FIXTURE0'
        bad=copy.deepcopy(self.design);bad['study_id']='FIXTURE1';bad['evidence'][0]['quote']='case 1'
        replies=[self.fixture_http_response(item) for item in [good,bad,bad]]
        output=self.root/'batch-run'
        with patch.dict('os.environ',{'OPENAI_API_KEY':'software-fixture-key'}), \
             patch('src.pipeline_runtime_2026_09_10.subprocess.run',side_effect=replies) as transport:
            preflight=run_source_design_batch(path,output,prepare_only=True)
            transport.assert_not_called();self.assertEqual(preflight['allocated_upper_limit_usd'],.2)
            first=run_source_design_batch(path,output)
            self.assertEqual(first['generation_attempts'],3)
            self.assertEqual([r['response_contract_status'] for r in first['rows']],['PASS','FAIL'])
            self.assertEqual(first['accepted_scientific_repairs'],0)
            snapshot={str(p):workflow.sha(p) for p in output.rglob('*') if p.is_file()}
            self.assertEqual(run_source_design_batch(path,output),first)
            self.assertEqual(snapshot,{str(p):workflow.sha(p) for p in output.rglob('*') if p.is_file()})
            self.assertEqual(transport.call_count,3)
            artifact=output/'studies/FIXTURE1/api/source_design-001/response.json'
            artifact.write_text(artifact.read_text()+'\n')
            with self.assertRaisesRegex(ValueError,'child artifacts changed'):
                run_source_design_batch(path,output)
            self.assertEqual(transport.call_count,3)

    def make_context(self):
        """Document fixtures test filtering/identity; they are not study evidence."""
        series = self.root/'series.soft'
        series.write_text('^SERIES = FIXTURE\n!Series_title = Design fixture\n'
            '!Series_overall_design = Separate experimental units; assigned drug or none.\n'
            '!Series_summary = RESULT_MUST_BE_WITHHELD\n!Series_pubmed_id = 1\n'
            + ''.join('!Series_sample_id = '+s+'\n' for s in ['A','B','C','D']))
        article = self.root/'PMC1.xml'
        article.write_text('<article><front><article-meta><article-id pub-id-type="pmcid">PMC1</article-id>'
            '<article-id pub-id-type="pmid">1</article-id><abstract><p>ABSTRACT_WITHHELD</p></abstract>'
            '</article-meta></front><body><sec><title>Results</title><p>RESULTS_WITHHELD</p>'
            '<sec sec-type="methods"><p>HIDDEN_METHODS_IN_RESULTS</p></sec></sec>'
            '<sec><title>Materials and Methods</title><sec><title>Study design</title>'
            '<p>Animals were assigned to separate intervention groups.</p></sec>'
            '<sec><title>Statistical analysis</title><p>PRIOR_STATISTICAL_RESULT</p></sec>'
            '<sec><title>Results</title><p>NESTED_RESULT_WITHHELD</p></sec></sec></body></article>')
        sources = [{'kind': kind, 'id': id_, 'path': path.name, 'sha256': workflow.sha(path),
                    'url': 'https://source.example.test/'+path.name} for kind,id_,path in [
                        ('geo_series','FIXTURE',series),('publication_jats','PMC1',article)]]
        context = {'version':CONTEXT_VERSION,'study_id':'FIXTURE','sources':sources,
                   'documents':context_documents('FIXTURE',sources,self.root),'acquisition':[]}
        path = self.root/'study-context.json'
        workflow.write_new(path,context)
        return path, context

    def test_study_context_excludes_abstract_results_and_mixed_sections(self):
        path,_ = self.make_context()
        packet = workflow.source_packet(self.path,'FIXTURE','drug versus none',context=path)
        serialized = json.dumps(packet)
        for withheld in ['RESULT_MUST_BE_WITHHELD','ABSTRACT_WITHHELD','RESULTS_WITHHELD',
                         'NESTED_RESULT_WITHHELD','HIDDEN_METHODS_IN_RESULTS','PRIOR_STATISTICAL_RESULT']:
            self.assertNotIn(withheld,serialized)
        self.assertIn('Separate experimental units',serialized)
        self.assertIn('Animals were assigned',serialized)

    def test_methods_quotes_checked_and_cannot_replace_sample_assignment_quotes(self):
        path,_ = self.make_context()
        packet = workflow.source_packet(self.path,'FIXTURE','drug versus none',context=path)
        quote = {'source_id':'PMC1:methods:1','field':'text','quote':'Animals were assigned'}
        design = {**self.design,'evidence':self.design['evidence']+[quote]}
        self.assertEqual(workflow.parse_design(json.dumps(design),packet),design)
        for evidence in [[quote], self.design['evidence']+[{**quote,'quote':'Independent animals confirmed'}],
                         self.design['evidence']+[{**quote,'source_id':'PMC999:methods:1'}]]:
            with self.subTest(evidence=evidence), self.assertRaises(ValueError):
                workflow.parse_design(json.dumps({**design,'evidence':evidence}),packet)

    def test_statistical_section_descendants_are_withheld_without_changing_source_ids(self):
        path,context=self.make_context()
        article=self.root/'PMC1.xml'
        article.write_text(article.read_text().replace(
            '<p>PRIOR_STATISTICAL_RESULT</p></sec>',
            '<sec><title>Estimation</title><p>NESTED_PRIOR_RESULT</p></sec></sec>'
            '<sec><title>Animal housing</title><p>Housing source fixture.</p></sec>'))
        context['sources'][1]['sha256']=workflow.sha(article)
        context['documents']=context_documents('FIXTURE',context['sources'],self.root)
        path.write_text(json.dumps(context))
        preserved=path.read_bytes()
        packet=workflow.source_packet(self.path,'FIXTURE','drug versus none',context=path)
        self.assertNotIn('NESTED_PRIOR_RESULT',json.dumps(packet))
        self.assertIn('NESTED_PRIOR_RESULT',json.dumps(read_study_context(path,'FIXTURE')))
        visible={doc['source_id']:doc for doc in packet['study_documents']}
        self.assertEqual(visible['PMC1:methods:3']['fields']['title'],'Animal housing')
        self.assertNotIn('PMC1:methods:2',visible)
        self.assertEqual(path.read_bytes(),preserved)

    def test_context_rejects_modified_raw_sources_and_modified_extractions(self):
        path,context = self.make_context()
        context['documents'][0]['fields']['overall_design']='Fabricated design'
        path.write_text(json.dumps(context))
        with self.assertRaisesRegex(ValueError,'differ from source'):
            read_study_context(path,'FIXTURE')
        (self.root/'series.soft').write_text('changed')
        with self.assertRaisesRegex(ValueError,'hash changed'):
            read_study_context(path,'FIXTURE')

    def test_context_rejects_cross_study_and_sample_mismatch(self):
        path,_ = self.make_context()
        with self.assertRaisesRegex(ValueError,'identity mismatch'):
            read_study_context(path,'ANOTHER')
        self.path.write_text(self.path.read_text().replace('D,control2','E,control2'))
        with self.assertRaisesRegex(ValueError,'sample membership differ'):
            workflow.source_packet(self.path,'FIXTURE','drug versus none',context=path)

    def test_publication_requires_matching_identity_and_body_or_geo_link(self):
        self.make_context()
        payload = (self.root/'PMC1.xml').read_bytes()
        with self.assertRaisesRegex(ValueError,'identity differs'):
            article_methods(payload,'GSE1','PMC2', ['1'])
        with self.assertRaisesRegex(ValueError,'no verified'):
            article_methods(payload,'GSE1','PMC1', [])
        with self.assertRaisesRegex(ValueError,'no verified'):
            article_methods(payload.replace(b'</article>',b'<back><ref-list><p>GSE1</p></ref-list></back></article>'),
                            'GSE1','PMC1', [])
        linked=payload.replace(b'Animals were assigned',b'GSE1: Animals were assigned')
        self.assertEqual(article_methods(linked,'GSE1','PMC1',[])[0]['link_basis'],
                         'EXACT_ACCESSION_IN_ARTICLE_BODY')
        with self.assertRaisesRegex(ValueError,'secondary reuse'):
            article_methods(linked,'GSE1','PMC1',['2'])
        with self.assertRaisesRegex(ValueError,'no verified'):
            article_methods(linked.replace(b'GSE1:',b'GSE10:'),'GSE1','PMC1',[])

    def test_context_cannot_read_outside_its_directory(self):
        path,context = self.make_context()
        context['sources'][0]['path']=str((self.root/'series.soft').resolve())
        path.write_text(json.dumps(context))
        with self.assertRaisesRegex(ValueError,'inside the context directory'):
            read_study_context(path,'FIXTURE')

    def test_source_quotations_and_complete_sample_partition_are_required(self):
        changes = [dict(case_ids=['A']), dict(control_ids=['A','C']),
                   dict(case_ids=['A','INVENTED']), dict(evidence=[]),
                   dict(evidence=[{'sample_id':'A','field':'title','quote':'invented quote'}])]
        for change in changes:
            with self.subTest(change=change), self.assertRaises(ValueError):
                workflow.parse_design(json.dumps({**self.design, **change}), self.packet)

    def test_unsupported_or_uncertain_design_cannot_execute_welch(self):
        for change in [dict(pairing='paired'), dict(pairing='unknown'),
                       dict(unresolved=['biological units unknown']), dict(findings=['confounded']),
                       dict(case_ids=['A'],control_ids=['B','C'],excluded_ids=['D']),
                       dict(status='PASS')]:
            with self.subTest(change=change), self.assertRaises(ValueError):
                workflow.parse_design(json.dumps({**self.design, **change}), self.packet)

    def test_duplicate_response_fields_rejected(self):
        raw=json.dumps(self.design).replace('"status":', '"status":"CANNOT VERIFY","status":')
        with self.assertRaises(ValueError):
            workflow.parse_design(raw,self.packet)

    def test_transport_failure_is_retained_and_resume_makes_no_second_call(self):
        with patch.object(workflow,'request_once',side_effect=RuntimeError('transport failure')) as call:
            a=workflow.run(self.path,'FIXTURE','drug versus none','tokenswitch',self.root/'run')
            b=workflow.run(self.path,'FIXTURE','drug versus none','tokenswitch',self.root/'run',resume=True)
        self.assertEqual(a,b)
        self.assertEqual(a['scientific_status'],'CANNOT VERIFY')
        self.assertEqual(call.call_count,1)
        self.assertIsNone(a['execution'])
        (self.root/'run/result.json').write_text('{"status":"PASS"}')
        with self.assertRaises(ValueError):
            workflow.run(self.path,'FIXTURE','drug versus none','tokenswitch',self.root/'run',resume=True)

    def test_prompt_limit_prevents_call_and_truncation(self):
        with patch.object(workflow,'request_once') as call, self.assertRaises(ValueError):
            workflow.run(self.path,'FIXTURE','drug versus none','tokenswitch',self.root/'run',max_prompt_bytes=10)
        call.assert_not_called()

    def test_complete_design_without_reference_is_not_scientific_acceptance(self):
        response={'model':'fixture-model','choices':[{'finish_reason':'stop','message':{'content':json.dumps(self.design)}}]}
        with patch.object(workflow,'request_once',return_value=response):
            result=workflow.run(self.path,'FIXTURE','drug versus none','openai',self.root/'run')
        self.assertEqual(result['design_status'],'DESIGN_PROPOSED')
        self.assertEqual(result['status'],'CANNOT VERIFY')
        self.assertIsNone(result['execution'])

    def test_second_implementation_checks_saved_fixture_and_detects_changed_result(self):
        root=Path(__file__).resolve().parents[2]/'paper/revision_2026-09-07/source/replay_fixture'
        samples=root/'GSE149866_samples.csv'; contract=root/'GSE149866_expected.json'
        require_development_artifact(samples)
        require_development_artifact(contract)
        replay=replay_saved_outputs(samples,contract)
        self.assertEqual(workflow.verify_replay(samples,contract,replay)['status'],'PASS')
        broken=copy.deepcopy(replay)
        broken['results'][0]['p_value']=.5
        self.assertEqual(workflow.verify_replay(samples,contract,broken)['status'],'FAIL')

    def test_cascade_escalates_complex_design_once_without_approving_consensus(self):
        calls=[]
        def child(provider,output,inputs,study,target,resume):
            calls.append((provider,inputs,target))
            output.mkdir()
            result={'status':'CANNOT VERIFY','design':{**self.design,'pairing':'paired'},
                    'generation_attempts':1,'billed_cost_usd':None}
            workflow.write_new(output/'result.json',result)
            return result
        result=workflow.run_cascade(self.path,'FIXTURE','drug versus none',self.root/'cascade',child=child)
        self.assertEqual([c[0] for c in calls],['tokenswitch','openai'])
        self.assertEqual(calls[0][1:],calls[1][1:])
        self.assertTrue(result['design_agreement'])
        self.assertEqual(result['generation_attempts'],2)
        self.assertFalse(result['accepted_result'])
        resumed=workflow.run_cascade(self.path,'FIXTURE','drug versus none',self.root/'cascade',resume=True,child=child)
        self.assertEqual(resumed,result)
        self.assertEqual(len(calls),2)

    def test_cascade_does_not_spend_second_call_on_missing_source_evidence_alone(self):
        def child(provider,output,*args):
            self.assertEqual(provider,'tokenswitch')
            output.mkdir()
            result={'status':'CANNOT VERIFY','generation_attempts':1,'billed_cost_usd':None,
                    'design':{**self.design,'pairing':'unknown','unresolved':['No donor identity']}}
            workflow.write_new(output/'result.json',result)
            return result
        result=workflow.run_cascade(self.path,'FIXTURE','drug versus none',self.root/'cascade',child=child)
        self.assertFalse(result['second_model_executed'])
        self.assertEqual(result['generation_attempts'],1)
        self.assertTrue(result['source_evidence_missing'])


if __name__=='__main__':
    unittest.main()
