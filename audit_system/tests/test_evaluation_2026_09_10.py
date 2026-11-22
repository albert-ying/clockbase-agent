"""Software fixtures only; these tests are not empirical audit error estimates."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from audit_system.evaluation_2026_09_10 import (summarize_records as evaluate, evaluate as verified_evaluate,
    build_families, verify_family_enrollment, export_reference_designs, freeze_reference_design,
    export_reference_outcomes, freeze_reference_outcomes, read_reference_records, COMPARISON_METHODS)


def fixtures():
    return [dict(study_id='FIXTURE'+str(i),family_id='FAMILY'+str(i//2),method=m,
                 reference_source='software fixture',reference_adjudicator='fixture authority',
                 reference_role='independent_source_adjudication',reference_status=s,
                 final_reference_status=s,final_reference_source='software fixture final calculation',
                 accepted=True,decision_status='PASS',flagged_material_error=False,repair_proposed=False,
                 repair_status='CANNOT VERIFY',runtime_seconds=1,cost_usd=None)
            for m in ['pipeline','audit'] for i,s in enumerate(['FAIL','PASS','CANNOT VERIFY','PASS'])]


class EvaluationTests(unittest.TestCase):
    def test_unavailable_case_remains_in_mixed_verified_denominator(self):
        from unittest.mock import patch
        from audit_system.tests.test_unified_pipeline_2026_09_10 import PipelineTests
        from audit_system.evaluation_2026_09_10 import (run_comparison_plan, reference_cases,
                                                       freeze_reference_unavailable)
        from audit_system.agent_workflow_2026_09_10 import sha
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            plan, prepared, _, _ = PipelineTests().attempted_plan_fixture(root, ready=True)
            class FailedClient:
                def invoke(self, messages): raise TimeoutError('Software fixture: no model response')
            archive = root / 'comparison'
            with patch('audit_system.report_2026_09_10.load_archived_comparison', return_value=prepared):
                run_comparison_plan(plan, archive, client_factory=lambda *args: FailedClient())
            case = next(row for row in reference_cases(archive) if row['study_id'] == 'GSE1')
            numerical, path, _, _ = ReferenceWorkflowTests().completed_numerical_fixture(root, archive, case)
            completed = freeze_reference_outcomes(numerical, path, root / 'completed-reference')
            missing = freeze_reference_unavailable(archive, 'GSE2', root / 'missing-reference')
            references = [{'path': item['frozen_outcomes'], 'sha256': item['sha256']} for item in (completed, missing)]
            result = verified_evaluate(reference_freezes=references, resamples=20)
            baseline = result['methods']['existing_pipeline']
            self.assertEqual(result['enrollment_counts']['expected_records'], 8)
            self.assertEqual(result['reference_unavailable_studies'], ['GSE2'])
            self.assertEqual(baseline['attempted_studies'], 2)
            self.assertEqual(baseline['accepted_studies'], 1)
            self.assertEqual(baseline['coverage'], .5)
            self.assertEqual(baseline['accepted_reference_unresolved_studies'], 1)
            self.assertIsNone(baseline['accepted_result_error_rate'])
            with self.assertRaisesRegex(ValueError, 'roster'):
                verified_evaluate(reference_freezes=references[:1], resamples=20)
            rows = read_reference_records(missing['frozen_outcomes'])
            altered = copy.deepcopy(rows); altered[0]['reference_status'] = 'PASS'
            with self.assertRaisesRegex(ValueError, 'Unavailable reference'):
                evaluate(altered, resamples=20)
            self.assertEqual(sha(missing['frozen_outcomes']), missing['sha256'])

    def test_frozen_roster_prevents_whole_study_or_method_disappearance(self):
        rows = fixtures()
        enrollment = {'schema': 'clockbase-comparison-enrollment-v1', 'split': 'development',
                      'methods': ['pipeline', 'audit'], 'studies': [
                          {'study_id': 'FIXTURE' + str(i), 'family_id': 'FAMILY' + str(i // 2)} for i in range(4)]}
        good = evaluate(rows, resamples=10, enrollment=enrollment)
        self.assertTrue(good['attempted_universe_verified'])
        self.assertEqual(good['enrollment_counts']['expected_records'], 8)
        for reduced in ([r for r in rows if r['study_id'] != 'FIXTURE0'],
                        [r for r in rows if r['method'] != 'audit'], rows[:-1]):
            with self.assertRaisesRegex(ValueError, 'roster'):
                evaluate(reduced, resamples=10, enrollment=enrollment)
        changed = copy.deepcopy(rows); changed[0]['family_id'] = 'another family'
        with self.assertRaisesRegex(ValueError, 'enrolled study family'):
            evaluate(changed, resamples=10, enrollment=enrollment)

    def source_fixture(self, accession, samples, **kwargs):
        return {'accession':accession,'sample_ids':samples,'source_kind':'geo_series_primary',
                'source_id':'software-fixture:'+accession,'unresolved':[],**kwargs}

    def test_transitive_family_exposure_spans_each_source_link_type(self):
        records=[self.source_fixture('GSE1',['GSM1']),
                 self.source_fixture('GSE2',['GSM1','GSM2'],bioproject_ids=['PRJNA10']),
                 self.source_fixture('GSE3',['GSM3'],bioproject_ids=['PRJNA10'],pubmed_ids=['12345']),
                 self.source_fixture('GSE4',['GSM4'],pubmed_ids=['12345'],linked_accessions=['GSE5']),
                 self.source_fixture('GSE5',['GSM5']),self.source_fixture('GSE6',['GSM6'])]
        result=build_families(records,{'GSE1'})
        family=next(f for f in result['families'] if 'GSE1' in f['members'])
        self.assertEqual(family['members'],['GSE1','GSE2','GSE3','GSE4','GSE5'])
        self.assertEqual(family['split'],'development')
        self.assertEqual({e['kind'] for e in result['edges']},{'GSM','BioProject','PubMed','GEO_relation'})
        reverse=build_families(list(reversed(records)),{'GSE1'})
        self.assertEqual(result['families'],reverse['families'])
        self.assertEqual(result['evaluation_release'],'CLOSED')

    def test_missing_linked_study_or_partial_membership_cannot_supply_held_out_family(self):
        records=[self.source_fixture('GSE1',['GSM1'],linked_accessions=['GSE2']),
                 self.source_fixture('GSE3',['GSM3'],source_kind='sample_metadata')]
        result=build_families(records,set())
        self.assertEqual({s['study_id'] for s in result['studies']},{'GSE1','GSE2','GSE3'})
        self.assertTrue(all(f['split']=='UNRESOLVED' for f in result['families']))
        self.assertTrue(all(f['source_linkage_status']=='CANNOT VERIFY' for f in result['families']))

    def test_source_membership_conflicts_and_namespace_collisions_are_not_erased(self):
        records=[self.source_fixture('GSE1',['GSM1']),
                 self.source_fixture('GSE1',['GSM1','GSM2'],source_kind='sample_metadata',source_id='fixture:second'),
                 self.source_fixture('GSE2',['GSM3'],pubmed_ids=['1'],bioproject_ids=['PRJNA1'])]
        result=build_families(records,set())
        self.assertEqual(len(result['families']),2)
        family=next(f for f in result['families'] if 'GSE1' in f['members'])
        self.assertEqual(family['split'],'UNRESOLVED')
        self.assertTrue(any('membership differs' in row['reason'] for row in family['unresolved']))

    def test_source_enrollment_rejects_relabeling_and_closed_held_out_release(self):
        result=build_families([self.source_fixture('GSE1',['GSM1']),self.source_fixture('GSE2',['GSM1'])],{'GSE1'})
        row={'study_id':'GSE2','family_id':result['families'][0]['family_id']}
        verify_family_enrollment([row],result,'development')
        for modified in ({**row,'family_id':'invented'}, {**row,'study_id':'GSE99'}):
            with self.assertRaises(ValueError):verify_family_enrollment([modified],result,'development')
        with self.assertRaisesRegex(ValueError,'release is closed'):
            verify_family_enrollment([row],result,'held_out')
        changed=copy.deepcopy(result); changed['families'][0]['exposed_members']=[]
        with self.assertRaisesRegex(ValueError,'exposure differs'):
            verify_family_enrollment([row],changed,'development')
        with self.assertRaisesRegex(ValueError,'source-linked family'):
            evaluate(fixtures(),resamples=10,split='held_out')

    def test_verified_repair_changes_final_error_without_erasing_original_detection_target(self):
        rows=fixtures()
        for row in rows:
            if row['method']=='audit' and row['reference_status']=='FAIL':
                row.update(flagged_material_error=True,repair_proposed=True,repair_status='PASS',final_reference_status='PASS')
        report=evaluate(rows,resamples=100)['methods']['audit']
        self.assertEqual(report['accepted_material_errors'],0)
        self.assertEqual(report['material_error_denominator'],1)
        self.assertEqual(report['repair_correctness'],1)
    def test_unknown_reference_is_not_correct_and_cost_unknown_is_retained(self):
        report=evaluate(fixtures(),resamples=100)['methods']['audit']
        self.assertEqual(report['accepted_result_error_rate'],1/3)
        self.assertEqual(report['accepted_error_lower_bound'],.25)
        self.assertEqual(report['accepted_error_upper_bound'],.5)
        self.assertEqual(report['coverage'],1)
        self.assertEqual(report['cost_known_studies'],0)
        self.assertEqual(report['family_count'],2)

    def test_zero_acceptance_has_no_zero_error_claim(self):
        rows=fixtures()
        for r in rows:
            if r['method']=='audit':r.update(accepted=False,decision_status='CANNOT VERIFY')
        report=evaluate(rows,resamples=100)['methods']['audit']
        self.assertIsNone(report['accepted_result_error_rate'])
        self.assertIsNone(report['accepted_result_error_rate_ci95'])
        self.assertEqual(report['coverage'],0)
        self.assertEqual(report['unresolved_studies'],4)
        self.assertEqual(report['rejected_studies'],0)

    def test_rejection_is_distinct_from_missing_evidence_and_comparison_is_paired(self):
        rows=fixtures()
        rows[4].update(accepted=False,decision_status='FAIL')
        report=evaluate(rows,resamples=100)
        self.assertEqual(report['methods']['audit']['rejected_studies'],1)
        self.assertEqual(report['methods']['audit']['unresolved_studies'],0)
        comparison=report['paired_method_comparisons'][0]
        self.assertEqual(comparison['method_b'],'pipeline')
        self.assertEqual(comparison['coverage_difference'],.25)
        self.assertIsNotNone(comparison['coverage_difference_ci95'])
        # Identical methods do not justify an interval with zero uncertainty.
        identical=evaluate(fixtures(),resamples=100)['paired_method_comparisons'][0]
        self.assertIsNone(identical['accepted_result_error_rate_difference_ci95'])

    def test_duplicate_studies_missing_attempts_or_agent_reference_rejected(self):
        a=fixtures();a.append(copy.deepcopy(a[0]))
        b=fixtures()[:-1]
        c=fixtures();c[0]['reference_role']='agent_agreement'
        d=fixtures();d[-1]['reference_status']='FAIL'
        for rows in [a,b,c,d]:
            with self.assertRaises(ValueError):evaluate(rows,resamples=100)


class ReferenceWorkflowTests(unittest.TestCase):
    """Software packets and assessor declarations, not empirical reference labels."""

    def archive_fixture(self, root, change_method=None, change_primary=None):
        from audit_system.tests.test_unified_pipeline_2026_09_10 import PipelineTests
        from audit_system.agent_workflow_2026_09_10 import sha, write_new
        from src.pipeline_runtime_2026_09_10 import StageStore
        packet, clocks, _, primary, _ = PipelineTests().measured_fixture(root)
        primary['rows'][0]['reason'] = 'NUMERICAL_OUTCOME_SENTINEL'
        if change_primary: change_primary(primary)
        source = root / 'software-fixture-primary.json'; write_new(source, primary)
        case = {'study_id': packet['study_id'], 'family_id': 'software-fixture-family',
                'packet': packet, 'primary': primary, 'clocks': clocks,
                'source_sha256': {str(source): sha(source), clocks['path']: sha(clocks['path'])}}
        archive = root / 'archive'; store = StageStore(archive, {'development_only': True})
        def inputs(directory):
            write_new(directory / 'prepared-cases.json', [case]); return {'cases': [case]}
        store.run('00_inputs', inputs)
        for method in COMPARISON_METHODS:
            (archive / case['study_id'] / method).mkdir(parents=True)
            def decision(directory):
                result = {'final_outputs': [{**primary['rows'][0], 'scientific_status': 'FAIL',
                    'accepted_result': False, 'selection_reason': 'METHOD_DECISION_SENTINEL'}], 'candidates': [],
                    'method': method, 'accepted_results': 0, 'decision_scope': 'Software fixture',
                    'study_decision_status': 'CANNOT VERIFY', 'checks': [], 'repair_proposed': False,
                    'flagged_material_error': False}
                if change_method: change_method(method, result, primary)
                write_new(directory / 'method-result.json', result); return result
            store.run(case['study_id'] + '/' + method + '/decision', decision)
        return archive, case

    def submitted_fixture(self, root, archive):
        prepared = export_reference_designs(archive, root / 'source-packets')
        packet = Path(prepared['packets'][0]['source_packet'])
        form = json.loads(Path(prepared['packets'][0]['form']).read_text())
        form['completed'] = True
        form['reviewer'].update(name='Independent software fixture assessor', independent_from_methods=True)
        for domain in form['domains'].values():
            domain['rationale'] = 'Software fixture: evidence remains unresolved'
        path = root / 'submitted-software-fixture.json'; path.write_text(json.dumps(form))
        return packet, path, form

    def test_source_reference_export_has_no_producer_outcomes_or_decisions_and_resumes(self):
        from audit_system.agent_workflow_2026_09_10 import sha
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            output = root / 'source-packets'; result = export_reference_designs(archive, output)
            packet = json.loads(Path(result['packets'][0]['source_packet']).read_text())
            self.assertEqual(packet['source_packet']['samples'], case['packet']['samples'])
            for path in output.rglob('*'):
                if path.suffix in {'.html', '.json'}:
                    self.assertNotIn('NUMERICAL_OUTCOME_SENTINEL', path.read_text())
                    self.assertNotIn('METHOD_DECISION_SENTINEL', path.read_text())
            serialized = json.dumps(packet)
            for field in ('"primary"', '"clocks"', '"p_value"', '"effect_size"', '"accepted_result"'):
                self.assertNotIn(field, serialized)
            self.assertFalse(result['reference_labels_completed'])
            before = {str(p): sha(p) for p in output.rglob('*') if p.is_file()}
            self.assertEqual(export_reference_designs(archive, output), result)
            self.assertEqual(before, {str(p): sha(p) for p in output.rglob('*') if p.is_file()})

    def test_uncompleted_form_cannot_release_numerical_reference_material(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, _ = self.archive_fixture(root)
            prepared = export_reference_designs(archive, root / 'source-packets')['packets'][0]
            with self.assertRaisesRegex(ValueError, 'uncompleted'):
                freeze_reference_design(prepared['source_packet'], prepared['form'], root / 'frozen')
            self.assertFalse((root / 'frozen').exists())

    def test_reference_freeze_requires_identity_rationale_and_supported_quotations(self):
        for change in ('name', 'role', 'independence', 'rationale', 'missing_quote', 'invented_quote', 'packet'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); archive, _ = self.archive_fixture(root)
                packet, judgment, form = self.submitted_fixture(root, archive)
                if change == 'name': form['reviewer']['name'] = ''
                if change == 'role': form['reviewer']['role'] = 'agent_agreement'
                if change == 'independence': form['reviewer']['independent_from_methods'] = False
                if change == 'rationale': form['domains']['sample_identity']['rationale'] = ''
                if change == 'missing_quote': form['domains']['sample_identity']['status'] = 'PASS'
                if change == 'invented_quote': form['domains']['sample_identity']['evidence'] = [
                    {'sample_id': 'GSM1', 'field': 'title', 'quote': 'Invented software source'}]
                if change == 'packet': form['packet_sha256'] = '0' * 64
                judgment.write_text(json.dumps(form))
                with self.assertRaises(ValueError): freeze_reference_design(packet, judgment, root / 'frozen')
                self.assertFalse((root / 'frozen').exists())

    def test_resolved_comparisons_require_exact_nonoverlapping_source_groups(self):
        for change in ('empty', 'overlap', 'invented', 'no_quote', 'duplicate'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); archive, case = self.archive_fixture(root)
                packet, judgment, form = self.submitted_fixture(root, archive)
                evidence = [{'sample_id': sample['sample_id'], 'field': 'characteristics_ch1.1',
                             'quote': sample['characteristics_ch1.1']} for sample in case['packet']['samples']]
                form['domains']['comparison_definition'].update(status='PASS', evidence=evidence)
                item = {'target': 'drugX versus vehicle', 'case_ids': ['GSM1', 'GSM2'],
                        'control_ids': ['GSM3', 'GSM4'], 'evidence': evidence}
                if change != 'empty': form['comparisons'] = [item]
                if change == 'overlap': item['control_ids'][0] = 'GSM1'
                if change == 'invented': item['control_ids'][0] = 'GSM999'
                if change == 'no_quote': item['evidence'] = []
                if change == 'duplicate': form['comparisons'].append(copy.deepcopy(item))
                judgment.write_text(json.dumps(form))
                with self.assertRaises(ValueError): freeze_reference_design(packet, judgment, root / 'frozen')

    def test_frozen_design_releases_anonymized_outputs_without_fabricating_labels(self):
        from audit_system.agent_workflow_2026_09_10 import sha
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, _ = self.archive_fixture(root)
            packet, judgment, _ = self.submitted_fixture(root, archive)
            frozen = freeze_reference_design(packet, judgment, root / 'frozen')
            self.assertFalse(frozen['scientific_reference_complete'])
            self.assertEqual(freeze_reference_design(packet, judgment, root / 'frozen'), frozen)
            output = root / 'numerical-reference'
            result = export_reference_outcomes(archive, frozen['frozen_design'], output)
            self.assertEqual(result['unique_outputs'], 1)
            self.assertEqual(result['mapped_output_occurrences'], 5)
            self.assertFalse(result['reference_labels_completed'])
            text = Path(result['reviewer_packet']).read_text()
            for phrase in ('NUMERICAL_OUTCOME_SENTINEL', 'METHOD_DECISION_SENTINEL', 'independent_audit',
                           'self_review', 'scientific_status', 'accepted_result', 'selection_reason'):
                self.assertNotIn(phrase, text)
            self.assertIn('p_value', text)
            mapping = json.loads(next(output.rglob('controller-output-mapping.json')).read_text())
            self.assertEqual({row['method'] for row in mapping}, set(COMPARISON_METHODS) | {'producer'})
            before = {str(p): sha(p) for p in output.rglob('*') if p.is_file()}
            self.assertEqual(export_reference_outcomes(archive, frozen['frozen_design'], output), result)
            self.assertEqual(before, {str(p): sha(p) for p in output.rglob('*') if p.is_file()})

    def test_numerical_export_preserves_effect_definitions_without_filling_legacy_rows(self):
        """The reviewer must receive the estimator declaration, including disagreement."""
        for definition in (None, 'cohen_d_pooled_sample_sd_case_minus_control'):
            with self.subTest(definition=definition), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                def change_primary(primary):
                    if definition is not None:
                        primary['rows'][0]['effect_size_definition'] = definition
                def change_method(method, result, primary):
                    if definition is not None and method == 'self_review':
                        # This conflicting estimator is a software fixture, not an approved repair.
                        result['final_outputs'][0]['effect_size_definition'] = 'glass_delta_control_sd'
                archive, _ = self.archive_fixture(root, change_method, change_primary)
                packet, judgment, _ = self.submitted_fixture(root, archive)
                frozen = freeze_reference_design(packet, judgment, root / 'frozen')
                released = export_reference_outcomes(archive, frozen['frozen_design'], root / 'numerical')
                content = json.loads(Path(released['reviewer_packet']).read_text())
                outputs = content['outcomes']
                self.assertEqual(released['mapped_output_occurrences'], 5)
                if definition is None:
                    self.assertEqual(len(outputs), 1)
                    self.assertNotIn('effect_size_definition', outputs[0])
                else:
                    self.assertEqual(len(outputs), 2)
                    self.assertEqual({row['effect_size_definition'] for row in outputs},
                                     {definition, 'glass_delta_control_sd'})
                    self.assertEqual(len({row['output_id'] for row in outputs}), 2)
                self.assertFalse(content['reference_labels_completed'])
                self.assertTrue(all('accepted_result' not in row and 'scientific_status' not in row
                                    for row in outputs))

    def test_changed_reference_input_or_freeze_cannot_open_numerical_review(self):
        from audit_system.agent_workflow_2026_09_10 import sha
        for change in ('judgment', 'packet', 'frozen', 'source'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); archive, case = self.archive_fixture(root)
                packet, judgment, form = self.submitted_fixture(root, archive)
                frozen = freeze_reference_design(packet, judgment, root / 'frozen')
                if change == 'judgment':
                    form['reviewer']['name'] = 'Changed assessor'; judgment.write_text(json.dumps(form))
                if change == 'packet': packet.write_text(packet.read_text() + '\n')
                if change == 'frozen':
                    path = Path(frozen['frozen_design']); data = json.loads(path.read_text())
                    data['judgment']['reviewer']['name'] = 'Forged local freeze'
                    path.write_text(json.dumps(data)); path.with_suffix('.sha256').write_text(sha(path))
                if change == 'source': Path(case['clocks']['path']).write_text('changed empirical source fixture\n')
                with self.assertRaises(ValueError):
                    export_reference_outcomes(archive, frozen['frozen_design'], root / 'numerical-reference')
                self.assertFalse((root / 'numerical-reference').exists())

    def test_source_export_rejects_unpinned_changes_and_held_out_archive(self):
        for change in ('metadata', 'source', 'held_out'):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); archive, case = self.archive_fixture(root)
                if change == 'metadata':
                    path = next(archive.glob('00_inputs/attempt-*/prepared-cases.json'))
                    path.write_text(path.read_text() + '\n')
                if change == 'source': Path(case['clocks']['path']).write_text('changed source fixture\n')
                if change == 'held_out':
                    path = archive / 'config.json'; data = json.loads(path.read_text())
                    data['development_only'] = False; path.write_text(json.dumps(data))
                with self.assertRaises(ValueError): export_reference_designs(archive, root / 'source-packets')
                self.assertFalse((root / 'source-packets').exists())

    def test_replay_reference_packets_follow_pinned_ancestry_and_detect_parent_changes(self):
        from audit_system.agent_workflow_2026_09_10 import sha, write_new
        from src.pipeline_runtime_2026_09_10 import StageStore
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            replay = root / 'replay'
            StageStore(replay, {'development_only': True, 'archive': str(archive),
                'archive_sha256': {str(path): sha(path) for path in archive.rglob('*') if path.is_file()}})
            result = export_reference_designs(replay, root / 'reference-packets')
            packet = json.loads(Path(result['packets'][0]['source_packet']).read_text())
            self.assertEqual(packet['source_packet']['samples'], case['packet']['samples'])
            write_new(archive / 'changed-parent.json', {'scope': 'Software source tamper fixture'})
            with self.assertRaisesRegex(ValueError, 'source archive changed'):
                export_reference_designs(replay, root / 'new-reference-packets')

    def completed_numerical_fixture(self, root, archive, case, resolved=False, frozen_target=None,
                                    design_change=None):
        from audit_system.agent_workflow_2026_09_10 import sha
        packet, submission, form = self.submitted_fixture(root, archive)
        evidence = [{'sample_id': row['sample_id'], 'field': 'characteristics_ch1.1',
                     'quote': row['characteristics_ch1.1']} for row in case['packet']['samples']]
        if resolved:
            for domain in form['domains'].values():
                domain.update(status='PASS', evidence=evidence,
                              rationale='Software fixture source-design judgment; not empirical evidence')
            row = case['primary']['rows'][0]
            form['comparisons'] = [{'target': frozen_target or case['packet']['target'], 'case_ids': row['case_ids'],
                                   'control_ids': row['control_ids'], 'evidence': evidence}]
        if design_change: design_change(form)
        submission.write_text(json.dumps(form))
        frozen = freeze_reference_design(packet, submission, root / 'design-freeze')
        numerical = root / 'numerical-reference'
        released = export_reference_outcomes(archive, frozen['frozen_design'], numerical)
        packet_path = Path(released['reviewer_packet']); packet = json.loads(packet_path.read_text())
        form = json.loads((packet_path.parent / 'outcome-form.json').read_text()); form['completed'] = True
        for item in form['assessments']:
            for domain in item['domains'].values(): domain['rationale'] = 'Unresolved software fixture judgment'
        if resolved:
            # These declarations are software fixtures, not empirical adjudications.
            proof = root / 'separate-software-fixture-calculation.json'
            proof.write_text(json.dumps({'scope': 'SOFTWARE FIXTURE', 'checked': True}))
            output_ids = [row['output_id'] for row in packet['outcomes']]
            inputs = {case['clocks']['path']: sha(case['clocks']['path'])}
            for source in case['clocks']['checks'][0]['measurement_sources'].values():
                inputs[source['path']] = source['sha256']
            for role in ('independent_calculation', 'measurement_evidence'):
                form['supporting_artifacts'].append({'artifact_id': role, 'path': str(proof), 'sha256': sha(proof),
                    'role': role, 'inputs': inputs, 'independent_from_methods': True, 'output_ids': output_ids})
            for item in form['assessments']:
                for name, domain in item['domains'].items():
                    role = {'numerical_accuracy': 'independent_calculation', 'measurement_applicability': 'measurement_evidence'}.get(name)
                    witness = ([{'artifact_id': role, 'pointer': '/checked', 'value': True}] if role else evidence)
                    domain.update(status='PASS', rationale='Software fixture resolved judgment', evidence=witness)
        path = root / 'submitted-numerical-software-fixture.json'; path.write_text(json.dumps(form))
        return numerical, path, form, packet

    def test_uncompleted_numerical_forms_and_bare_labels_cannot_supply_evaluation(self):
        with self.assertRaisesRegex(ValueError, 'provenance strings'):
            verified_evaluate(fixtures(), resamples=10)
        self.assertFalse(evaluate(fixtures(), resamples=10)['reference_artifacts_verified'])
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            numerical, path, form, _ = self.completed_numerical_fixture(root, archive, case)
            form['completed'] = False; path.write_text(json.dumps(form))
            with self.assertRaisesRegex(ValueError, 'uncompleted'):
                freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')

    def test_numerical_pass_cannot_erase_unresolved_or_failed_source_design(self):
        for source_status, material, expected in [('CANNOT VERIFY', None, 'CANNOT VERIFY'),
                                                  ('FAIL', None, 'CANNOT VERIFY'),
                                                  ('FAIL', True, 'FAIL'), ('FAIL', False, 'PASS')]:
            with self.subTest(status=source_status, material=material), tempfile.TemporaryDirectory() as folder:
                root = Path(folder); archive, case = self.archive_fixture(root)
                def change(form):
                    form['domains']['biological_units'].update(status=source_status, material=material)
                numerical, path, _, _ = self.completed_numerical_fixture(
                    root, archive, case, resolved=True, design_change=change)
                frozen = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
                value = json.loads(Path(frozen['frozen_outcomes']).read_text())
                self.assertEqual(set(value['labels'].values()), {expected})
                self.assertEqual({row['status'] for row in value['source_design_constraints'].values()}, {expected})

    def test_comparison_design_is_frozen_before_numerical_release_and_has_its_own_scope(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            def change(form):
                form['comparisons'][0]['domains'] = copy.deepcopy(form['domains'])
                form['domains']['pairing'].update(status='CANNOT VERIFY', rationale='Other comparisons remain unresolved')
            numerical, path, _, _ = self.completed_numerical_fixture(
                root, archive, case, resolved=True, design_change=change)
            frozen = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
            value = json.loads(Path(frozen['frozen_outcomes']).read_text())
            self.assertEqual(set(value['labels'].values()), {'PASS'})
            self.assertEqual({row['scope'] for row in value['source_design_constraints'].values()},
                             {'Comparison-specific source domains'})

    def test_numerical_freeze_preserves_unknown_labels_and_resumes(self):
        from audit_system.agent_workflow_2026_09_10 import sha
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            numerical, path, _, _ = self.completed_numerical_fixture(root, archive, case)
            destination = root / 'outcome-freeze'; result = freeze_reference_outcomes(numerical, path, destination)
            before = {str(p): sha(p) for p in destination.rglob('*') if p.is_file()}
            self.assertEqual(result['reference_resolved_outputs'], 0)
            records = read_reference_records(result['frozen_outcomes'])
            self.assertEqual(len(records), 4)
            report = verified_evaluate(reference_freezes=[{'path': result['frozen_outcomes'], 'sha256': result['sha256']}], resamples=10)
            self.assertTrue(report['reference_artifacts_verified'])
            for method in report['methods'].values():
                self.assertEqual(method['attempted_studies'], 1)
                self.assertEqual(method['reference_resolved_studies'], 0)
                self.assertIsNone(method['accepted_result_error_rate'])
            self.assertEqual(freeze_reference_outcomes(numerical, path, destination), result)
            self.assertEqual(before, {str(p): sha(p) for p in destination.rglob('*') if p.is_file()})

    def test_verified_evaluation_rejects_missing_whole_study_reference(self):
        from audit_system.agent_workflow_2026_09_10 import write_new
        from src.pipeline_runtime_2026_09_10 import StageStore
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); source_root = root / 'source'; source_root.mkdir()
            archive, case = self.archive_fixture(source_root)
            second = copy.deepcopy(case); second.update(study_id='SECOND_SOFTWARE_FIXTURE', family_id='second fixture family')
            second['packet']['study_id'] = second['study_id']
            combined = root / 'combined'; store = StageStore(combined, {'development_only': True})
            def inputs(directory):
                write_new(directory / 'prepared-cases.json', [case, second]); return {'cases': [case, second]}
            store.run('00_inputs', inputs)
            for study in (case['study_id'], second['study_id']):
                for method in COMPARISON_METHODS:
                    original = json.loads((archive / case['study_id'] / method / 'decision/complete.json').read_text())['result']
                    (combined / study / method).mkdir(parents=True)
                    def decision(directory, result=original):
                        write_new(directory / 'method-result.json', result); return result
                    store.run(study + '/' + method + '/decision', decision)
            numerical, path, _, _ = self.completed_numerical_fixture(root, combined, case)
            frozen = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
            with self.assertRaisesRegex(ValueError, 'roster'):
                verified_evaluate(reference_freezes=[{'path': frozen['frozen_outcomes'], 'sha256': frozen['sha256']}], resamples=10)

    def test_verified_evaluation_rejects_mixed_independently_selected_archives(self):
        with tempfile.TemporaryDirectory() as folder:
            refs = []
            for index in range(2):
                root = Path(folder) / str(index); root.mkdir()
                archive, case = self.archive_fixture(root)
                numerical, path, _, _ = self.completed_numerical_fixture(root, archive, case)
                frozen = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
                refs.append({'path': frozen['frozen_outcomes'], 'sha256': frozen['sha256']})
            with self.assertRaisesRegex(ValueError, 'common enrolled method archive'):
                verified_evaluate(reference_freezes=refs, resamples=10)

    def test_numerical_judgment_requires_full_coverage_source_witnesses_and_materiality(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            numerical, path, original, _ = self.completed_numerical_fixture(root, archive, case, resolved=True)
            for change in ('missing', 'duplicate', 'reviewer', 'rationale', 'quote', 'pointer', 'materiality', 'calculation', 'measurement', 'wrong_output'):
                with self.subTest(change=change):
                    form = copy.deepcopy(original); domain = form['assessments'][0]['domains']['numerical_accuracy']
                    if change == 'missing': form['assessments'] = []
                    if change == 'duplicate': form['assessments'].append(copy.deepcopy(form['assessments'][0]))
                    if change == 'reviewer': form['reviewer']['name'] = 'Another fixture name'
                    if change == 'rationale': domain['rationale'] = ''
                    if change == 'quote': domain['evidence'] = [{'sample_id': 'GSM1', 'field': 'characteristics_ch1.1', 'quote': 'invented'}]
                    if change == 'pointer': domain['evidence'][0]['value'] = False
                    if change == 'materiality': domain.update(status='FAIL', material=None)
                    if change == 'calculation': domain['evidence'] = form['assessments'][0]['domains']['statistical_validity']['evidence']
                    if change == 'measurement': form['assessments'][0]['domains']['measurement_applicability']['evidence'] = domain['evidence']
                    if change == 'wrong_output': form['supporting_artifacts'][0]['output_ids'] = ['missing-output']
                    path.write_text(json.dumps(form))
                    with self.assertRaises((ValueError, KeyError)):
                        freeze_reference_outcomes(numerical, path, root / ('bad-' + change))
                    self.assertFalse((root / ('bad-' + change)).exists())

    def test_original_error_accepted_result_and_repair_are_separately_referenced(self):
        from audit_system.agent_workflow_2026_09_10 import sha
        corrected = {}
        def original(primary):
            corrected.update(copy.deepcopy(primary['rows'][0])); primary['rows'][0]['p_value'] = 1e-40
        def methods(method, result, primary):
            result['final_outputs'] = [{**primary['rows'][0], 'accepted_result': True}]
            result.update(accepted_results=1, study_decision_status='PASS')
            if method == 'independent_audit':
                result['candidates'] = [copy.deepcopy(corrected)]
                result['final_outputs'] = [{**corrected, 'accepted_result': True}]
                result['checks'] = [{'clock': corrected['clock'], 'repair_kind': 'CORRECTION_PROPOSAL'}]
                result.update(repair_proposed=True, flagged_material_error=True)
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root, methods, original)
            numerical, path, form, packet = self.completed_numerical_fixture(root, archive, case, resolved=True)
            bad = {row['output_id'] for row in packet['outcomes'] if row['p_value'] == 1e-40}
            for item in form['assessments']:
                if item['output_id'] in bad:
                    item['domains']['numerical_accuracy'].update(status='FAIL', material=True,
                        rationale='Software fixture: original P value was deliberately corrupted')
            path.write_text(json.dumps(form)); frozen = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
            refs = [{'path': frozen['frozen_outcomes'], 'sha256': frozen['sha256']}]
            report = verified_evaluate(reference_freezes=refs, resamples=10)['methods']
            self.assertEqual(report['existing_pipeline']['accepted_material_errors'], 1)
            self.assertEqual(report['independent_audit']['accepted_material_errors'], 0)
            self.assertEqual(report['independent_audit']['verified_correct_repairs'], 1)
            self.assertEqual(report['independent_audit']['material_error_denominator'], 1)
            self.assertEqual(report['existing_pipeline']['missed_material_errors'], 1)
            rows = read_reference_records(frozen['frozen_outcomes']); rows[0]['reference_status'] = 'PASS'
            with self.assertRaisesRegex(ValueError, 'differ from frozen'):
                verified_evaluate(rows, reference_freezes=refs, resamples=10)
            with self.assertRaises(ValueError): verified_evaluate(reference_freezes=refs * 2, resamples=10)
            with self.assertRaisesRegex(ValueError, 'held-out'):
                verified_evaluate(reference_freezes=refs, split='held_out', resamples=10)
            before = {str(p): sha(p) for p in numerical.rglob('*') if p.is_file()}
            verified_evaluate(reference_freezes=refs, resamples=10)
            self.assertEqual(before, {str(p): sha(p) for p in numerical.rglob('*') if p.is_file()})

    def test_changed_reference_proof_and_prediction_invalidate_frozen_metrics(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            numerical, path, form, _ = self.completed_numerical_fixture(root, archive, case, resolved=True)
            result = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
            for source in (path, Path(form['supporting_artifacts'][0]['path']), Path(case['clocks']['path']),
                           Path(result['frozen_outcomes'])):
                with self.subTest(source=source.name):
                    before = source.read_bytes(); source.write_bytes(before + b'\n')
                    try:
                        with self.assertRaises(ValueError): read_reference_records(result['frozen_outcomes'])
                    finally: source.write_bytes(before)

    def test_unnecessary_and_new_analysis_candidates_are_not_verified_repairs(self):
        def methods(method, result, primary):
            result['candidates'] = [copy.deepcopy(primary['rows'][0])]
            result['final_outputs'] = [{**primary['rows'][0], 'accepted_result': True}]
            kind = 'CORRECTION_PROPOSAL' if method == 'independent_audit' else 'NEW_ANALYSIS_AFTER_ABSTENTION'
            result['checks'] = [{'clock': primary['rows'][0]['clock'], 'repair_kind': kind}]
            result.update(accepted_results=1, study_decision_status='PASS', repair_proposed=kind == 'CORRECTION_PROPOSAL')
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root, methods)
            numerical, path, _, _ = self.completed_numerical_fixture(root, archive, case, resolved=True)
            frozen = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
            rows = read_reference_records(frozen['frozen_outcomes'])
            self.assertTrue(all(row['repair_status'] == 'CANNOT VERIFY' for row in rows))
            self.assertEqual(sum(row['repair_proposed'] for row in rows), 1)

    def test_failed_candidate_does_not_contaminate_retained_original_reference(self):
        def methods(method, result, primary):
            result['final_outputs'] = [{**primary['rows'][0], 'accepted_result': True}]
            result.update(accepted_results=1, study_decision_status='PASS')
            if method == 'independent_audit':
                result['candidates'] = [{**primary['rows'][0], 'p_value': 1e-40}]
                result['checks'] = [{'clock': primary['rows'][0]['clock'], 'repair_kind': 'CORRECTION_PROPOSAL'}]
                result.update(repair_proposed=True, flagged_material_error=True)
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root, methods)
            numerical, path, form, packet = self.completed_numerical_fixture(root, archive, case, resolved=True)
            wrong = {row['output_id'] for row in packet['outcomes'] if row['p_value'] == 1e-40}
            for item in form['assessments']:
                if item['output_id'] in wrong:
                    item['domains']['numerical_accuracy'].update(status='FAIL', material=True,
                        rationale='Software fixture: proposed calculation is incorrect')
            path.write_text(json.dumps(form)); frozen = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
            report = verified_evaluate(reference_freezes=[{'path': frozen['frozen_outcomes'], 'sha256': frozen['sha256']}], resamples=10)
            audit = report['methods']['independent_audit']
            self.assertEqual(audit['accepted_material_errors'], 0)
            self.assertEqual(audit['verified_correct_repairs'], 0)
            self.assertEqual(audit['verified_repairs'], 1)
            self.assertEqual(audit['false_alarms'], 1)

    def test_nonmaterial_failure_does_not_become_material_error(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            numerical, path, form, _ = self.completed_numerical_fixture(root, archive, case, resolved=True)
            form['assessments'][0]['domains']['numerical_accuracy'].update(status='FAIL', material=False,
                rationale='Software fixture: disclosed rounding does not affect interpretation')
            path.write_text(json.dumps(form)); frozen = freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')
            self.assertTrue(all(row['reference_status'] == 'PASS' for row in read_reference_records(frozen['frozen_outcomes'])))

    def test_same_sample_groups_do_not_certify_another_declared_target(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); archive, case = self.archive_fixture(root)
            numerical, path, _, _ = self.completed_numerical_fixture(root, archive, case, resolved=True,
                frozen_target='Different software fixture hypothesis')
            with self.assertRaisesRegex(ValueError, 'comparison and target'):
                freeze_reference_outcomes(numerical, path, root / 'outcome-freeze')

    def test_unknown_runtime_is_not_reported_as_zero_total(self):
        rows = fixtures(); rows[0]['runtime_seconds'] = None
        report = evaluate(rows, resamples=10)['methods']['pipeline']
        self.assertIsNone(report['total_runtime_seconds'])
        self.assertEqual(report['runtime_known_studies'], 3)


if __name__=='__main__':unittest.main()
