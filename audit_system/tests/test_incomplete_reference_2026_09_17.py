"""Software fixtures only; no empirical labels, held-out sources or model calls."""
import contextlib
import copy
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from audit_system import evaluation_2026_09_10 as evaluation
from audit_system.agent_workflow_2026_09_10 import sha, write_new
from audit_system.tests import test_evaluation_2026_09_10 as legacy
from src.pipeline_runtime_2026_09_10 import StageStore


class IncompleteReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        # Stage progress is retained by production runs; fixture output is noise.
        self.stdout = contextlib.redirect_stdout(io.StringIO())
        self.stdout.__enter__()
        self.addCleanup(self.stdout.__exit__, None, None, None)

    @staticmethod
    def accepted(method, result, primary):
        if method == 'existing_pipeline':
            result['final_outputs'][0].update(accepted_result=True, scientific_status='PASS')
            result.update(accepted_results=1, study_decision_status='PASS')

    def archive(self, change_method=None):
        return legacy.ReferenceWorkflowTests().archive_fixture(self.root, change_method or self.accepted)

    def missing(self, archive, case, phase='NUMERICAL_REVIEW', name='missing'):
        evidence = self.root / (name + '-attempt.json')
        write_new(evidence, {'scope': 'SOFTWARE FIXTURE', 'status': 'UNAVAILABLE', 'phase': phase})
        path = self.root / (name + '.json')
        write_new(path, {'schema': 'clockbase-reference-missingness-v1',
            'study_id': case['study_id'], 'family_id': case['family_id'], 'phase': phase,
            'archive_config_sha256': sha(archive / 'config.json'),
            'reason': 'Software fixture: bounded reference attempt returned no usable judgment',
            'evidence': [{'path': str(evidence), 'sha256': sha(evidence)}]})
        return path

    def freeze(self, archive, case, missing=None, name='incomplete', **kwargs):
        missing = missing or self.missing(archive, case)
        result = evaluation.freeze_reference_incomplete(archive, case['study_id'], missing,
                                                        self.root / name, **kwargs)
        return result, evaluation.read_reference_records(result['frozen_outcomes'])

    @staticmethod
    def reference(result):
        return {'path': result['frozen_outcomes'], 'sha256': result['sha256']}

    @staticmethod
    def stage(store, name, filename, value):
        (store.root / name).parent.mkdir(parents=True, exist_ok=True)
        def save(directory):
            write_new(directory / filename, value)
            return value
        return store.run(name, save)

    def source_design(self, archive, case, status='PASS', material=None):
        packet, submitted, form = legacy.ReferenceWorkflowTests().submitted_fixture(self.root, archive)
        evidence = [{'sample_id': row['sample_id'], 'field': 'characteristics_ch1.1',
                     'quote': row['characteristics_ch1.1']} for row in case['packet']['samples']]
        for domain in form['domains'].values(): domain.update(status='PASS', evidence=evidence)
        form['domains']['biological_units'].update(status=status, material=material)
        row = case['primary']['rows'][0]
        form['comparisons'] = [{'target': case['packet']['target'], 'case_ids': row['case_ids'],
                               'control_ids': row['control_ids'], 'evidence': evidence}]
        submitted.write_text(json.dumps(form))
        return evaluation.freeze_reference_design(packet, submitted, self.root / 'design-freeze')['frozen_design']

    def test_executed_acceptance_is_not_reset_by_missing_reference(self):
        archive, case = self.archive()
        missing = self.missing(archive, case)
        result, rows = self.freeze(archive, case, missing)
        baseline = next(row for row in rows if row['method'] == 'existing_pipeline')
        self.assertTrue(baseline['accepted'])
        self.assertEqual(baseline['decision_status'], 'PASS')
        self.assertIsNone(baseline['reference_adjudicator'])
        self.assertEqual(set(baseline['reference_output_ids']), {'original', 'accepted', 'final', 'candidates'})
        self.assertEqual(baseline['accepted_output_count'], 1)
        self.assertEqual(baseline['reference_status'], 'CANNOT VERIFY')
        report = evaluation.evaluate(reference_freezes=[self.reference(result)], resamples=20)
        stats = report['methods']['existing_pipeline']
        self.assertEqual(report['enrollment_counts']['expected_records'], 4)
        self.assertEqual(stats['coverage'], 1)
        self.assertEqual(stats['accepted_reference_unresolved_studies'], 1)
        self.assertIsNone(stats['accepted_result_error_rate'])
        self.assertEqual((stats['accepted_error_lower_bound'], stats['accepted_error_upper_bound']), (0, 1))
        before = {str(path): sha(path) for path in (self.root / 'incomplete').rglob('*') if path.is_file()}
        self.assertEqual(result, evaluation.freeze_reference_incomplete(archive, case['study_id'], missing,
                                                                       self.root / 'incomplete'))
        self.assertEqual(before, {str(path): sha(path) for path in (self.root / 'incomplete').rglob('*') if path.is_file()})
        changed = copy.deepcopy(rows); changed[0]['accepted'] = False
        with self.assertRaisesRegex(ValueError, 'differ from frozen'):
            evaluation.evaluate(changed, reference_freezes=[self.reference(result)], resamples=10)
        with self.assertRaisesRegex(ValueError, 'shared-prerequisite'):
            evaluation.freeze_reference_unavailable(archive, case['study_id'], self.root / 'wrong-adapter')

    def test_known_usage_and_unknown_billing_survive_missingness(self):
        archive, case = self.archive()
        for number in (0, 1):
            slot = archive / 'api' / (case['study_id'] + '-independent_audit-%03d' % number)
            slot.mkdir(parents=True)
            write_new(slot / 'reservation.json', {'reserved_usd': 2})
            if number == 0:
                write_new(slot / 'receipt.json', {'estimated_cost_usd': .03,
                    'usage': {'prompt_tokens': 100, 'completion_tokens': 20}, 'returned_model': 'software-fixture'})
        _, rows = self.freeze(archive, case)
        audit = next(row for row in rows if row['method'] == 'independent_audit')
        self.assertEqual((audit['generation_attempts'], audit['responses_with_usage'], audit['unknown_billing_attempts']), (2, 1, 1))
        self.assertEqual(audit['estimated_known_cost_usd'], .03)
        self.assertEqual((audit['prompt_tokens'], audit['completion_tokens']), (100, 20))
        self.assertIsNone(audit['estimated_cost_usd'])
        self.assertIsNone(audit['billed_cost_usd'])
        self.assertIsNone(audit['cost_usd'])

    def test_source_pass_does_not_invent_numerical_pass_and_material_failure_survives(self):
        for status, material, expected in [('PASS', None, 'CANNOT VERIFY'), ('FAIL', True, 'FAIL'),
                                            ('FAIL', None, 'CANNOT VERIFY')]:
            with self.subTest(status=status, material=material):
                previous = self.root; self.root = previous / (status + str(material)); self.root.mkdir()
                try:
                    archive, case = self.archive()
                    design = self.source_design(archive, case, status, material)
                    result, rows = self.freeze(archive, case, frozen_design=design)
                    self.assertEqual({row['reference_status'] for row in rows}, {expected})
                    baseline = next(row for row in rows if row['accepted'])
                    self.assertEqual(baseline['final_reference_status'], expected)
                    self.assertEqual(baseline['known_reference_sources'], [design])
                    self.assertEqual(evaluation.evaluate(reference_freezes=[self.reference(result)], resamples=10)
                        ['methods']['existing_pipeline']['accepted_material_errors'], int(expected == 'FAIL'))
                finally: self.root = previous

    def test_partial_numerical_labels_and_verified_repair_are_not_discarded(self):
        def repaired(method, result, primary):
            self.accepted(method, result, primary)
            if method == 'independent_audit':
                candidate = {**primary['rows'][0], 'p_value': .25, 'mean_difference': .5,
                             'accepted_result': True, 'scientific_status': 'PASS'}
                result.update(final_outputs=[candidate], candidates=[candidate], accepted_results=1,
                    study_decision_status='PASS', repair_proposed=True, flagged_material_error=True,
                    checks=[{'repair_kind': 'CORRECTION_PROPOSAL', 'clock': candidate['clock'],
                             'primary_analysis_state': {}}])
            if method == 'self_review':
                result['candidates'] = [{**primary['rows'][0], 'p_value': .314159, 'mean_difference': .75,
                                         'accepted_result': False, 'scientific_status': 'CANNOT VERIFY'}]
        archive, case = self.archive(repaired)
        numerical, judgment, form, packet = legacy.ReferenceWorkflowTests().completed_numerical_fixture(
            self.root, archive, case, resolved=True)
        original = case['primary']['rows'][0]
        original_id = next(row['output_id'] for row in packet['outcomes'] if row['p_value'] == original['p_value'])
        unresolved_id = next(row['output_id'] for row in packet['outcomes'] if row['p_value'] == .314159)
        for item in form['assessments']:
            if item['output_id'] == original_id:
                item['domains']['numerical_accuracy'].update(status='FAIL', material=True)
            if item['output_id'] == unresolved_id:
                for domain in item['domains'].values(): domain.update(status='CANNOT VERIFY', material=None)
        judgment.write_text(json.dumps(form))
        complete = evaluation.freeze_reference_outcomes(numerical, judgment, self.root / 'complete')
        legacy_rows = evaluation.read_reference_records(complete['frozen_outcomes'])
        result, rows = self.freeze(archive, case, frozen_outcomes=complete['frozen_outcomes'])
        retained = json.loads(Path(result['frozen_outcomes']).read_text())['state']['labels']
        self.assertEqual(set(retained.values()), {'PASS', 'FAIL', 'CANNOT VERIFY'})
        self.assertEqual(retained, json.loads(Path(complete['frozen_outcomes']).read_text())['labels'])
        for old, new in zip(legacy_rows, rows):
            for key in ('accepted', 'decision_status', 'reference_status', 'final_reference_status',
                        'repair_proposed', 'repair_status', 'flagged_material_error'):
                self.assertEqual(old[key], new[key])
        audit = next(row for row in rows if row['method'] == 'independent_audit')
        self.assertEqual(audit['repair_status'], 'PASS')
        self.assertEqual(len(audit['reference_output_ids']['candidates']), 1)
        self.assertEqual(evaluation.evaluate(reference_freezes=[self.reference(result)], resamples=10)
                         ['methods']['independent_audit']['verified_correct_repairs'], 1)
        self.assertNotIn('reference_incomplete_studies', evaluation.evaluate(
            reference_freezes=[self.reference(complete)], resamples=10)['methods']['independent_audit'])
        self.assertEqual(set(legacy_rows[0]['reference_output_ids']), {'original', 'accepted'})

    def test_sealed_prior_design_is_loaded_without_optional_argument(self):
        archive, case = self.archive()
        design = self.source_design(archive, case, 'FAIL', True)
        store = StageStore(archive, json.loads((archive / 'config.json').read_text()))
        state = {'study_id': case['study_id'], 'references': [
            {'phase': 'SOURCE_DESIGN', 'path': design, 'sha256': sha(design)}], 'missingness': []}
        self.stage(store, case['study_id'] + '/reference_state', 'reference-state.json', state)
        missing = self.missing(archive, case)
        _, rows = self.freeze(archive, case, missing)
        self.assertEqual({row['reference_status'] for row in rows}, {'FAIL'})
        with self.assertRaisesRegex(ValueError, 'replace a sealed reference'):
            self.freeze(archive, case, missing, name='replacement', frozen_design=self.root / 'invented')
        with self.assertRaisesRegex(ValueError, 'sealed pre-numerical'):
            self.freeze(archive, case, self.missing(archive, case, 'SOURCE_DESIGN', 'late-source'), name='late')

    def test_evidence_method_and_freeze_tampering_are_rejected(self):
        for target in ('evidence', 'method', 'freeze', 'source'):
            with self.subTest(target=target):
                previous = self.root; self.root = previous / target; self.root.mkdir()
                try:
                    archive, case = self.archive()
                    result, _ = self.freeze(archive, case)
                    path = {'evidence': self.root / 'missing-attempt.json',
                            'method': next(archive.rglob('method-result.json')),
                            'freeze': Path(result['frozen_outcomes']),
                            'source': Path(case['clocks']['path'])}[target]
                    path.write_text(path.read_text() + '\n')
                    with self.assertRaises(ValueError): evaluation.read_reference_records(result['frozen_outcomes'])
                finally: self.root = previous

    def test_missing_method_and_mixed_archives_cannot_define_denominator(self):
        archive, case = self.archive()
        missing_stage = archive / case['study_id'] / 'self_review' / 'decision' / 'complete.json'
        saved = missing_stage.read_bytes(); missing_stage.unlink()
        missing = self.missing(archive, case)
        with self.assertRaises(FileNotFoundError): self.freeze(archive, case, missing)
        missing_stage.write_bytes(saved)
        first, _ = self.freeze(archive, case, missing)
        previous = self.root; self.root = previous / 'other'; self.root.mkdir()
        try:
            other, other_case = self.archive()
            second, _ = self.freeze(other, other_case)
        finally: self.root = previous
        with self.assertRaisesRegex(ValueError, 'common enrolled method archive'):
            evaluation.evaluate(reference_freezes=[self.reference(first), self.reference(second)], resamples=10)

    def prospective(self, *, unavailable=False, split='development', extra_study=False):
        """Contract fixture; the controller joint-release validator is mocked separately."""
        original, case = self.archive()
        original_study, study = case['study_id'], 'GSE1'
        case = copy.deepcopy(case); case['study_id'] = study; case['packet']['study_id'] = study
        studies = ['GSE1', 'GSE2'] if extra_study else ['GSE1']
        manifest = evaluation.build_families([{'accession': item, 'sample_ids': ['GSM1'],
            'source_kind': 'geo_series_primary', 'source_id': 'SOFTWARE FIXTURE ' + item,
            'unresolved': []} for item in studies], set(studies))
        # This is explicitly a software split fixture, never a scientific enrollment.
        if split == 'held_out':
            manifest['evaluation_release'] = 'OPEN'; manifest['exposed_studies'] = []
            manifest['families'][0].update(split=split, exposed_members=[])
            for item in manifest['studies']: item.update(split=split, family_exposed=False, directly_exposed=False)
        case = copy.deepcopy(case); case['family_id'] = manifest['families'][0]['family_id']
        enrollment = {'schema': 'clockbase-comparison-enrollment-v1', 'split': split,
            'methods': list(evaluation.COMPARISON_METHODS),
            'studies': [{'study_id': study, 'family_id': case['family_id']}]}
        if extra_study:
            enrollment['studies'].append({'study_id': 'GSE2', 'family_id': case['family_id']})
        archive = self.root / 'prospective'
        store = StageStore(archive, {'schema': 'clockbase-prospective-comparison-v1',
                                    'split': split, 'development_only': split == 'development'})
        def enroll(directory):
            path = directory / 'enrollment.json'; write_new(path, enrollment)
            return {'enrollment': enrollment, 'sha256': sha(path)}
        enrolled = store.run('00_enrollment', enroll)
        source_case = copy.deepcopy(case)
        if unavailable:
            name = study + '/coding_agent'; (archive / name).parent.mkdir(parents=True)
            def failed(directory):
                value = {'phase': name, 'status': 'CANNOT VERIFY', 'automatic_retry': False,
                         'error': 'SOFTWARE FIXTURE: producer attempt failed'}
                path = directory / 'terminal-attempt.json'; write_new(path, value)
                return {**value, 'terminal_receipt': str(path), 'terminal_sha256': sha(path)}
            terminal = store.run(name, failed); receipt = terminal['terminal_receipt']
            attempt = {'schema': 'clockbase-prospective-unavailable-attempt-v1',
                'enrollment_sha256': enrolled['sha256'], 'phase': name, 'failure_receipt': receipt,
                'failure_sha256': sha(receipt), 'recorded_failure': json.loads(Path(receipt).read_text()),
                'upstream_config': str(archive / 'config.json'), 'upstream_config_sha256': sha(archive / 'config.json')}
            case = {'study_id': study, 'family_id': case['family_id'],
                'preparation_status': 'UNAVAILABLE_PROSPECTIVE_PREPARATION', 'attempt': attempt,
                'source_sha256': {receipt: sha(receipt), str(archive / 'config.json'): sha(archive / 'config.json')}}
        cases = [case]
        if extra_study:
            second = copy.deepcopy(case); second['study_id'] = 'GSE2'
            second['packet']['study_id'] = second['study_id']; cases.append(second)
        def inputs(directory):
            write_new(directory / 'prepared-cases.json', cases); return {'cases': cases}
        store.run('00_inputs', inputs)
        for prepared in cases:
            if unavailable:
                evaluation.compare_unavailable_case(store, prepared)
            else:
                for method in evaluation.COMPARISON_METHODS:
                    result, _ = evaluation.reference_stage(original, original_study + '/' + method + '/decision', 'method-result.json')
                    result.update(study_id=prepared['study_id'], family_id=prepared['family_id'])
                    self.stage(store, prepared['study_id'] + '/' + method + '/decision', 'method-result.json', result)
        context = {'split': split, 'enrollment': enrollment, 'family_manifest': manifest}
        self.enterContext(patch('audit_system.heldout_2026_09_17.validate_archive_context', return_value=context))
        self.enterContext(patch('audit_system.heldout_2026_09_17.validate_reference_packet_context', return_value=context))
        return archive, case, source_case, store, context

    def test_prospective_preparation_failure_keeps_producer_cost_without_error_labels(self):
        archive, case, _, store, _ = self.prospective(unavailable=True)
        missing = self.missing(archive, case, 'SOURCE_DESIGN')
        self.stage(store, case['study_id'] + '/reference_state', 'reference-state.json',
            {'study_id': case['study_id'], 'references': [], 'missingness': [
                {'phase': 'SOURCE_DESIGN', 'path': str(missing), 'sha256': sha(missing)}]})
        for number in (0, 1):
            slot = archive / 'baseline-api' / case['study_id'] / ('primary-%03d' % number)
            slot.mkdir(parents=True); write_new(slot / 'reservation.json', {'reserved_usd': 1})
            if number == 0:
                write_new(slot / 'receipt.json', {'estimated_cost_usd': .04,
                    'usage': {'prompt_tokens': 100, 'completion_tokens': 50}, 'returned_model': 'software-fixture'})
        result, rows = self.freeze(archive, case, missing)
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(not row['accepted'] and row['reference_status'] == 'CANNOT VERIFY' for row in rows))
        self.assertEqual(rows[0]['reference_output_ids'], dict.fromkeys(('original', 'accepted', 'final', 'candidates'), []))
        shared = rows[0]['shared_upstream_resources']
        self.assertEqual((shared['generation_attempts'], shared['unknown_billing_attempts']), (2, 1))
        self.assertEqual(shared['estimated_known_cost_usd'], .04)
        self.assertIsNone(shared['estimated_cost_usd'])
        self.assertEqual(rows[0]['generation_attempts'], 0)
        report = evaluation.evaluate(reference_freezes=[self.reference(result)], resamples=10)
        self.assertEqual(len(report['shared_upstream_resources_by_study']), 1)
        self.assertEqual(report['methods']['existing_pipeline']['material_error_denominator'], 0)
        self.assertIsNone(report['methods']['existing_pipeline']['accepted_result_error_rate'])
        with self.assertRaisesRegex(ValueError, 'Legacy unavailable'):
            evaluation.freeze_reference_unavailable(archive, case['study_id'], self.root / 'legacy')

    def test_prospective_held_out_missing_reference_requires_joint_context_and_complete_roster(self):
        archive, case, _, store, context = self.prospective(split='held_out', extra_study=True)
        missing = self.missing(archive, case, 'SOURCE_DESIGN')
        self.stage(store, case['study_id'] + '/reference_state', 'reference-state.json',
            {'study_id': case['study_id'], 'references': [], 'missingness': [
                {'phase': 'SOURCE_DESIGN', 'path': str(missing), 'sha256': sha(missing)}]})
        second_case = next(row for row in evaluation.reference_cases(archive) if row['study_id'] == 'GSE2')
        second_missing = self.missing(archive, second_case, 'SOURCE_DESIGN', 'second-missing')
        self.stage(store, second_case['study_id'] + '/reference_state', 'reference-state.json',
            {'study_id': second_case['study_id'], 'references': [], 'missingness': [
                {'phase': 'SOURCE_DESIGN', 'path': str(second_missing), 'sha256': sha(second_missing)}]})
        first, _ = self.freeze(archive, case, missing)
        with self.assertRaisesRegex(ValueError, 'roster'):
            evaluation.evaluate(reference_freezes=[self.reference(first)], split='held_out', resamples=10)
        with self.assertRaisesRegex(ValueError, 'controller source-first'):
            evaluation.export_reference_designs(archive, self.root / 'post-hoc-packets')
        with self.assertRaisesRegex(ValueError, 'different development/held-out'):
            evaluation.evaluate(reference_freezes=[self.reference(first)], split='development', resamples=10)
        with self.assertRaisesRegex(ValueError, 'joint freeze'):
            evaluation.evaluate(reference_freezes=[self.reference(first)], split='held_out', resamples=10,
                                family_manifest={**context['family_manifest'], 'evaluation_release': 'CLOSED'})
        second, _ = self.freeze(archive, second_case, second_missing, name='second-incomplete')
        report = evaluation.evaluate(reference_freezes=[self.reference(first), self.reference(second)],
                                     split='held_out', resamples=10)
        self.assertEqual(report['enrollment_counts']['expected_records'], 8)
        self.assertEqual(report['methods']['existing_pipeline']['accepted_studies'], 2)
        self.assertEqual(report['methods']['existing_pipeline']['coverage'], 1)
        self.assertIsNone(report['methods']['existing_pipeline']['accepted_result_error_rate'])
        self.assertEqual(report['uncertainty_unit'], 'study family')

    def prospective_design(self, archive, case, store, *, failure=False):
        packet = {'schema': 'clockbase-reference-design-packet-v2', 'phase': 'SOURCE_ONLY',
            'study_id': case['study_id'], 'family_id': case['family_id'],
            'split': json.loads((archive / 'config.json').read_text())['split'],
            'evaluation_context': {'archive': str(archive), 'archive_config_sha256': sha(archive / 'config.json'),
                                   'source_stage': case['study_id'] + '/source_packet'},
            'source_packet': evaluation.reference_source_view(case['packet'])}
        source_directory = self.root / 'prospective-source'; source_directory.mkdir()
        source_path = source_directory / 'source-design.json'; write_new(source_path, packet)
        evidence = [{'sample_id': row['sample_id'], 'field': 'characteristics_ch1.1',
                     'quote': row['characteristics_ch1.1']} for row in packet['source_packet']['samples']]
        form = {'schema': 'clockbase-reference-design-judgment-v1', 'completed': True,
            'packet_sha256': sha(source_path), 'study_id': case['study_id'], 'family_id': case['family_id'],
            'reviewer': {'name': 'Software fixture assessor', 'role': 'independent_source_adjudication',
                         'independent_from_methods': True},
            'domains': {name: {'status': 'PASS', 'material': None, 'rationale': 'SOFTWARE FIXTURE',
                                'evidence': evidence} for name in evaluation.REFERENCE_DESIGN_DOMAINS},
            'comparisons': [{'target': case['packet']['target'], 'case_ids': case['primary']['rows'][0]['case_ids'],
                'control_ids': case['primary']['rows'][0]['control_ids'], 'evidence': evidence}],
            'additional_sources_needed': []}
        if failure: form['domains']['biological_units'].update(status='FAIL', material=True)
        submission = self.root / 'prospective-source-submission.json'; write_new(submission, form)
        frozen = evaluation.freeze_reference_design(source_path, submission, self.root / 'prospective-design')
        state = {'study_id': case['study_id'], 'references': [{'phase': 'SOURCE_DESIGN',
            'path': frozen['frozen_design'], 'sha256': frozen['sha256']}], 'missingness': []}
        self.stage(store, case['study_id'] + '/reference_state', 'reference-state.json', state)
        return frozen

    def test_prospective_complete_reference_uses_registered_source_design_and_new_versions(self):
        archive, case, _, store, _ = self.prospective(split='held_out')
        design = self.prospective_design(archive, case, store)
        self.assertEqual(evaluation.read_reference_design(design['frozen_design'])['schema'],
                         'clockbase-frozen-reference-design-v2')
        numerical = self.root / 'prospective-numerical'
        released = evaluation.export_reference_outcomes(archive, design['frozen_design'], numerical)
        packet = json.loads(Path(released['reviewer_packet']).read_text())
        self.assertEqual((packet['schema'], packet['split']), ('clockbase-reference-outcome-packet-v3', 'held_out'))
        form = json.loads((Path(released['reviewer_packet']).parent / 'outcome-form.json').read_text())
        form['completed'] = True
        for item in form['assessments']:
            for domain in item['domains'].values(): domain['rationale'] = 'SOFTWARE FIXTURE: unresolved reference domain'
        judgment = self.root / 'prospective-outcome-submission.json'; write_new(judgment, form)
        complete = evaluation.freeze_reference_outcomes(numerical, judgment, self.root / 'prospective-complete')
        self.assertEqual(json.loads(Path(complete['frozen_outcomes']).read_text())['schema'],
                         'clockbase-frozen-reference-outcomes-v2')
        report = evaluation.evaluate(reference_freezes=[self.reference(complete)], split='held_out', resamples=10)
        self.assertEqual(report['methods']['existing_pipeline']['coverage'], 1)
        self.assertEqual(report['enrollment_counts']['expected_records'], 4)
        wrapped, rows = self.freeze(archive, case, frozen_outcomes=complete['frozen_outcomes'])
        self.assertTrue(next(row['accepted'] for row in rows if row['method'] == 'existing_pipeline'))
        self.assertEqual(len(json.loads(Path(wrapped['frozen_outcomes']).read_text())['state']['known_reference_sources']), 2)

    def test_failed_preparation_preserves_known_source_failure_without_inventing_output_errors(self):
        archive, case, source_case, store, _ = self.prospective(unavailable=True)
        design = self.prospective_design(archive, source_case, store, failure=True)
        result, rows = self.freeze(archive, case)
        state = json.loads(Path(result['frozen_outcomes']).read_text())['state']
        self.assertEqual(state['known_source_design']['domains']['biological_units']['status'], 'FAIL')
        self.assertEqual(state['labels'], {})
        self.assertEqual({row['reference_status'] for row in rows}, {'CANNOT VERIFY'})
        self.assertEqual(state['known_reference_sources'], [design['frozen_design']])

    def test_integrity_failures_propagate_and_are_not_terminal_method_abstentions(self):
        from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError
        archive, case = self.archive()
        for kind in ('review', 'decision'):
            with self.subTest(kind=kind):
                output = self.root / ('integrity-' + kind); store = StageStore(output, {'development_only': True})
                target = 'src.pipeline_2026_09_10.' + ('assess_design' if kind == 'review' else 'numerical_audit')
                with patch(target, side_effect=PipelineIntegrityError('Software fixture source hash changed')):
                    with self.assertRaises(PipelineIntegrityError):
                        evaluation.compare_development_case(store, case, lambda *args: object())
                method = 'independent_audit' if kind == 'review' else 'deterministic'
                self.assertFalse((output / case['study_id'] / method / kind / 'complete.json').exists())
        output = self.root / 'ordinary-review-failure'; store = StageStore(output, {'development_only': True})
        with patch('src.pipeline_2026_09_10.assess_design', side_effect=ValueError('Software fixture unparseable response')):
            evaluation.compare_development_case(store, case, lambda *args: object())
        result, _ = evaluation.reference_stage(output, case['study_id'] + '/independent_audit/review', 'review-result.json')
        self.assertEqual(result['execution_status'], 'CANNOT VERIFY')

    def test_prospective_missing_source_state_cannot_be_repaired_after_execution(self):
        archive, case, source_case, store, _ = self.prospective()
        missing = self.missing(archive, case, 'SOURCE_DESIGN')
        with self.assertRaisesRegex(ValueError, 'sealed prior reference-state'):
            self.freeze(archive, case, missing)
        self.stage(store, case['study_id'] + '/reference_state', 'reference-state.json',
            {'study_id': case['study_id'], 'references': [], 'missingness': [
                {'phase': 'SOURCE_DESIGN', 'path': str(missing), 'sha256': sha(missing)}]})
        with self.assertRaisesRegex(ValueError, 'late source-design'):
            self.freeze(archive, case, missing, frozen_design=self.root / 'late-design.json')

    def test_missingness_receipt_cannot_supply_labels_or_wrong_context(self):
        archive, case = self.archive(); missing = self.missing(archive, case)
        original = json.loads(missing.read_text())
        for update in ({'labels': {'output': 'PASS'}}, {'study_id': 'other'}, {'phase': 'afterthought'},
                       {'evidence': []}, {'archive_config_sha256': '0' * 64}):
            with self.subTest(update=update):
                missing.write_text(json.dumps({**original, **update}))
                with self.assertRaisesRegex(ValueError, 'administrative'):
                    self.freeze(archive, case, missing)

    def test_prospective_analysis_settings_cannot_override_the_frozen_plan(self):
        archive, case, _, store, context = self.prospective()
        context['plan'] = {'analysis': {'primary_endpoint': 'accepted_result_error_rate',
            'study_weighting': 'equal', 'uncertainty': 'paired_family_cluster_percentile_bootstrap',
            'resamples': 10, 'seed': 20260910}}
        missing = self.missing(archive, case, 'SOURCE_DESIGN')
        self.stage(store, case['study_id'] + '/reference_state', 'reference-state.json',
            {'study_id': case['study_id'], 'references': [], 'missingness': [
                {'phase': 'SOURCE_DESIGN', 'path': str(missing), 'sha256': sha(missing)}]})
        result, _ = self.freeze(archive, case, missing)
        references = [self.reference(result)]
        self.assertEqual(evaluation.evaluate(reference_freezes=references, resamples=10)['resamples'], 10)
        for settings in ({'resamples': 9}, {'resamples': 10, 'seed': 20260911}):
            with self.subTest(settings=settings), self.assertRaisesRegex(ValueError, 'frozen plan'):
                evaluation.evaluate(reference_freezes=references, **settings)
        for field, value in [('primary_endpoint', 'coverage'), ('study_weighting', 'equal_family'),
                             ('uncertainty', 'unpaired_bootstrap')]:
            previous = context['plan']['analysis'][field]; context['plan']['analysis'][field] = value
            try:
                with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'frozen plan'):
                    evaluation.evaluate(reference_freezes=references, resamples=10)
            finally: context['plan']['analysis'][field] = previous

    def test_incomplete_reference_cli_dispatch_preserves_acceptance(self):
        archive, case = self.archive(); missing = self.missing(archive, case)
        command = [sys.executable, '-B', '-W', 'ignore', '-m', 'audit_system.evaluation_2026_09_10',
            '--freeze-reference-incomplete', str(archive), '--study', case['study_id'],
            '--missingness-receipt', str(missing), '--output', str(self.root / 'cli-freeze')]
        process = subprocess.run(command, cwd=Path(__file__).resolve().parents[2],
                                 capture_output=True, text=True, timeout=30)
        self.assertEqual(process.returncode, 0, process.stderr)
        paths = list((self.root / 'cli-freeze').rglob('frozen-outcomes.json'))
        self.assertEqual(len(paths), 1)
        rows = evaluation.read_reference_records(paths[0])
        self.assertTrue(next(row['accepted'] for row in rows if row['method'] == 'existing_pipeline'))


if __name__ == '__main__':
    unittest.main()
