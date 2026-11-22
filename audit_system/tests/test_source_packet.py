"""Software fixtures only; no scientific performance claims."""
import importlib
import importlib.util
import unittest
import json
import csv
from pathlib import Path
import tempfile
from audit_system.tests import require_development_artifact


class LocalEntryPointTests(unittest.TestCase):
    """Metadata examples are software fixtures; replay uses the preserved source fixture."""

    def setUp(self):
        self.m = importlib.import_module('audit_system.source_packet')
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.metadata = self.root/'metadata.csv'
        self.metadata.write_text('GSM,title,batch,p_value,Hazard_clock,gpt_verdict\n'
                                 'GSM1,case,run1,.001,999,PASS\n'
                                 'GSM2,control,run2,.9,-999,FAIL\n')
        self.fixture = Path(__file__).resolve().parents[2]/'paper/revision_2026-09-07/source/replay_fixture'

    def test_adapter_retains_batch_but_withholds_numerical_and_prior_verdict_fields(self):
        packet = self.m.metadata_packet(self.metadata, 'GSE1')
        payload = json.dumps(packet)
        for excluded in ['p_value', 'Hazard_clock', 'gpt_verdict', '999']:
            self.assertNotIn(excluded, payload)
        self.assertIn('run1', payload)
        self.assertEqual(packet['comparisons'], [])
        self.assertEqual(len(packet['samples']), 2)

    def test_method_alias_conflict_fails_without_choosing_a_source(self):
        self.metadata.write_text('GSM,title,data_processing,raw.data_processing\n'
                                 'GSM1,software fixture,beta,fluorescence\n')
        with self.assertRaisesRegex(ValueError, 'Conflicting native and raw method'):
            self.m.metadata_packet(self.metadata, 'GSE1')

    def test_method_narratives_preserve_multiline_values_and_exclude_suffix_lookalikes(self):
        row = {'GSM': 'GSM1', 'title': 'software fixture', 'data_processing': 'First line\nSecond line',
               'data_processing.1': 'Other step', 'raw.scan_protocol.1': 'Exact raw method',
               'data_processing.p_value': 'HIDDEN_RESULT', 'raw.scan_protocol.p_value': 'HIDDEN_RESULT',
               'data_processing_extra': 'HIDDEN_RESULT', 'contact_email': 'HIDDEN_CONTACT'}
        with self.metadata.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row)); writer.writeheader(); writer.writerow(row)
        packet = self.m.metadata_packet(self.metadata, 'GSE1')
        self.assertEqual(packet['samples'][0]['raw.data_processing'], row['data_processing'])
        self.assertEqual(packet['samples'][0]['raw.data_processing.1'], row['data_processing.1'])
        self.assertEqual(packet['samples'][0]['raw.scan_protocol.1'], row['raw.scan_protocol.1'])
        self.assertNotIn('HIDDEN_', json.dumps(packet))

    def test_duplicate_clock_names_fail_even_when_values_equal(self):
        self.metadata.write_text('GSM,clock,clock\nGSM1,1,1\n')
        with self.assertRaisesRegex(ValueError, 'Ambiguous duplicate'):
            self.m.read_source_csv(self.metadata)

    def test_repeated_sample_columns_require_rowwise_equality(self):
        self.metadata.write_text('sample_id,sample_id,title\nGSM1,GSM1,A\n')
        packet = self.m.metadata_packet(self.metadata, 'GSE1')
        self.assertEqual(packet['adapter']['collapsed_identical_id_columns'], ['sample_id'])
        self.metadata.write_text('sample_id,sample_id,title\nGSM1,GSM2,A\n')
        with self.assertRaisesRegex(ValueError, 'Conflicting duplicate'):
            self.m.metadata_packet(self.metadata, 'GSE1')

    def test_empty_ragged_and_pointer_inputs_fail(self):
        for content in ['GSM,title\n', 'GSM,title\nGSM1,A,extra\n',
                        'version https://git-lfs.github.com/spec/v1\noid sha256:abc\n']:
            self.metadata.write_text(content)
            with self.assertRaises(ValueError):
                self.m.metadata_packet(self.metadata, 'GSE1')

    def test_resume_preserves_output_and_rejects_changed_input(self):
        output = self.root/'run'
        inputs = {'metadata.csv': self.metadata}
        self.m.run_local('prepare', inputs, output, study_id='GSE1')
        before = (output/'manifest.json').read_bytes()
        self.m.run_local('prepare', inputs, output, study_id='GSE1', resume=True)
        self.assertEqual(before, (output/'manifest.json').read_bytes())
        with self.assertRaisesRegex(ValueError, 'Output exists'):
            self.m.run_local('prepare', inputs, output, study_id='GSE1')
        self.metadata.write_text(self.metadata.read_text().replace('run1', 'run3'))
        with self.assertRaisesRegex(ValueError, 'Resume rejected'):
            self.m.run_local('prepare', inputs, output, study_id='GSE1', resume=True)
        self.assertEqual(before, (output/'manifest.json').read_bytes())

    def test_invalid_input_does_not_publish_completion(self):
        output = self.root/'run'
        self.metadata.write_text('GSM,title\nGSM1,A\nGSM1,B\n')
        with self.assertRaises(ValueError):
            self.m.run_local('prepare', {'metadata.csv': self.metadata}, output, study_id='GSE1')
        self.assertFalse(output.exists())

    def test_empirical_saved_fixture_is_not_scientific_approval(self):
        samples = self.fixture/'GSE149866_samples.csv'
        contract = self.fixture/'GSE149866_expected.json'
        require_development_artifact(samples)
        require_development_artifact(contract)
        result = self.m.replay_saved_outputs(samples, contract)
        self.assertEqual(result['status'], 'PASS')
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(len(result['results']), 2)

    def test_partition_loss_overlap_singletons_and_unsupported_test_rejected(self):
        contract_path = self.fixture/'GSE149866_expected.json'
        require_development_artifact(contract_path)
        original = json.loads(contract_path.read_text())
        for modification in ['loss', 'overlap', 'singleton', 'paired']:
            contract = json.loads(json.dumps(original))
            if modification == 'loss':
                contract['background_gsm_ids_not_in_drug_test'].pop()
            elif modification == 'overlap':
                contract['control_gsm_ids'][0] = contract['case_gsm_ids'][0]
            elif modification == 'singleton':
                contract['background_gsm_ids_not_in_drug_test'] += contract['case_gsm_ids'][1:]
                contract['case_gsm_ids'] = contract['case_gsm_ids'][:1]
            else:
                contract['test'] = 'paired t-test'
            path = self.root/'contract.json'
            path.write_text(json.dumps(contract))
            with self.subTest(modification=modification), self.assertRaises(ValueError):
                samples = self.fixture/'GSE149866_samples.csv'
                require_development_artifact(samples)
                self.m.replay_saved_outputs(samples, path)

    def test_saved_numerical_mismatch_is_fail_not_repaired(self):
        contract_path = self.fixture/'GSE149866_expected.json'
        samples = self.fixture/'GSE149866_samples.csv'
        require_development_artifact(contract_path)
        require_development_artifact(samples)
        contract = json.loads(contract_path.read_text())
        first = next(iter(contract['expected_from_saved_analysis']))
        contract['expected_from_saved_analysis'][first]['p_value'] = 0.9
        path = self.root/'contract.json'
        path.write_text(json.dumps(contract))
        result = self.m.replay_saved_outputs(samples, path)
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(json.loads(path.read_text()), contract)


class SourcePacketTests(unittest.TestCase):
    def module(self):
        self.assertIsNotNone(importlib.util.find_spec('audit_system.source_packet'))
        return importlib.import_module('audit_system.source_packet')

    def test_preserves_every_sample_without_outcome_leakage(self):
        m = self.module()
        rows = [{'gse_id': 'GSE1', 'gsm_id': 'GSM1', 'title': 'A', 'raw_characteristics': 'time: 1', 'Hazard_clock': '9', 'gpt_case_control': 'case'},
                {'gse_id': 'GSE1', 'gsm_id': 'GSM2', 'title': 'B', 'raw_characteristics': 'time: 2', 'Hazard_clock': '-9'}]
        result = m.build_packet('GSE1', rows, [{'comparison_name': 'B vs A', 'p_value': '.001', 'effect_size': '9'}])
        self.assertEqual([r['gsm_id'] for r in result['samples']], ['GSM1', 'GSM2'])
        self.assertEqual(result['samples'][1]['raw_characteristics'], 'time: 2')
        self.assertNotIn('Hazard_clock', result['samples'][0])
        self.assertNotIn('gpt_case_control', result['samples'][0])
        self.assertNotIn('p_value', result['comparisons'][0])
        self.assertNotIn('effect_size', result['comparisons'][0])

    def test_duplicate_sample_keys_rejected(self):
        m = self.module()
        with self.assertRaises(ValueError):
            m.build_packet('GSE1', [{'gse_id': 'GSE1', 'gsm_id': 'GSM1'}] * 2, [])

    def test_mixed_study_samples_rejected(self):
        m = self.module()
        with self.assertRaises(ValueError):
            m.build_packet('GSE1', [{'gse_id': 'GSE2', 'gsm_id': 'GSM1'}], [])

    def test_incomplete_or_unknown_reply_rejected(self):
        m = self.module()
        for reply in ['{"status":"PASS"', '{}', '{"status":"PASS"}', '{"status":"PASS","findings":[],"unresolved":[]}']:
            with self.assertRaises(ValueError):
                m.parse_review(reply)
        review = m.parse_review('{"status":"NEEDS_REVIEW","findings":[],"unresolved":["replicate identity"]}')
        self.assertEqual(review['status'], 'NEEDS_REVIEW')

    def test_design_matrix_reports_rank_and_estimable_contrast_with_nuisance_alias(self):
        m = self.module()
        # The third column duplicates the intercept, so the nuisance block is
        # rank-deficient while the treatment contrast remains in the row space.
        result = m.design_matrix_diagnostic(
            [[1, 0, 1], [1, 1, 1], [1, 0, 1], [1, 1, 1]],
            ['intercept', 'treatment', 'batch_alias'],
            ['S1', 'S2', 'S3', 'S4'],
            {'treatment': [0, 1, 0], 'zero': [0, 0, 0]},
        )
        self.assertEqual((result['row_count'], result['column_count'], result['rank']), (4, 3, 2))
        self.assertEqual(result['residual_df'], 2)
        self.assertEqual(result['contrasts']['treatment']['status'], 'ESTIMABLE')
        self.assertEqual(result['contrasts']['zero']['status'], 'ESTIMABLE')
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(result['biological_independence_status'], 'CANNOT VERIFY')
        self.assertEqual(result['measurement_applicability_status'], 'CANNOT VERIFY')

    def test_fully_confounded_contrast_is_not_estimable(self):
        m = self.module()
        result = m.design_matrix_diagnostic(
            [[1, 0, 0], [1, 0, 0], [1, 1, 1], [1, 1, 1]],
            ['intercept', 'treatment', 'batch'],
            ['S1', 'S2', 'S3', 'S4'],
            {'treatment': [0, 1, 0]},
        )
        self.assertEqual(result['rank'], 2)
        self.assertEqual(result['contrasts']['treatment']['status'], 'NOT_ESTIMABLE')
        self.assertGreater(result['contrasts']['treatment']['row_space_residual_norm'], 1e-10)
        self.assertEqual(result['scientific_status'], 'CANNOT VERIFY')

    def test_design_matrix_rejects_identity_shape_finiteness_and_tolerance_errors(self):
        m = self.module()
        valid = ([[1, 0], [1, 1]], ['intercept', 'treatment'], ['S1', 'S2'], {'treatment': [0, 1]})
        cases = [
            (valid[0], ['intercept', 'intercept'], valid[2], valid[3]),
            (valid[0], ['', 'treatment'], valid[2], valid[3]),
            (valid[0], valid[1], ['S1', 'S1'], valid[3]),
            (valid[0], valid[1], ['', 'S2'], valid[3]),
            ([[1, 0]], valid[1], valid[2], valid[3]),
            ([[1, float('nan')], [1, 1]], valid[1], valid[2], valid[3]),
            (valid[0], valid[1], valid[2], {'treatment': [0]}),
        ]
        for design, names, samples, contrasts in cases:
            with self.subTest(design=design, names=names, samples=samples, contrasts=contrasts), self.assertRaises(ValueError):
                m.design_matrix_diagnostic(design, names, samples, contrasts)
        with self.assertRaises(ValueError):
            m.design_matrix_diagnostic(*valid, rank_rtol=0)
        with self.assertRaises(ValueError):
            m.design_matrix_diagnostic(*valid, estimability_tolerance=-1)


class SourceConsistencyTests(unittest.TestCase):
    """Constructed software fixtures; no empirical performance estimates."""

    def setUp(self):
        self.m = importlib.import_module('audit_system.source_packet')
        self.packet = {'study_id': 'FIXTURE', 'source_sha256': 'fixture',
            'target': 'treatment: exposed versus treatment: vehicle within time: late',
            'samples': [{'sample_id': sid, 'characteristics_ch1': 'treatment: ' + value,
                         'characteristics_ch1.1': 'time: late', 'platform_id': platform}
                        for sid, value, platform in [('A', 'exposed', 'P1'), ('B', 'exposed', 'P1'),
                                                     ('C', 'vehicle', 'P2'), ('D', 'vehicle', 'P2')]]}

    def check(self):
        return self.m.verify_requested_comparison(self.packet, ['A', 'B'], ['C', 'D'])

    def test_absent_category_does_not_hide_explicit_record_contradiction(self):
        for row in self.packet['samples']:
            row['characteristics_ch1'] = 'treatment: vehicle'
        result = self.check()
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual({r['sample_id'] for r in result['explicit_source_contradictions']}, {'A', 'B'})
        self.assertEqual(result['explicit_source_contradictions'][0]['evidence'][0]['quote'], 'treatment: vehicle')
        self.assertFalse(result['reference_judgment'])

    def test_wrong_requested_stratum_is_explicit_contradiction(self):
        self.packet['samples'][0]['characteristics_ch1.1'] = 'time: early'
        result = self.check()
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual(result['explicit_source_contradictions'][0]['source_label'], 'time')

    def test_missing_unknown_or_conflicting_declarations_are_not_contradictions(self):
        for value in (None, 'treatment: unknown', 'treatment: NA', 'treatment: missing'):
            row = self.packet['samples'][0]
            row.pop('characteristics_ch1', None)
            if value is not None:
                row['characteristics_ch1'] = value
            with self.subTest(value=value):
                self.assertEqual(self.check()['status'], 'CANNOT VERIFY')
        self.packet['samples'][0]['characteristics_ch1'] = 'treatment: other'
        self.packet['samples'][0]['characteristics_ch1.2'] = 'treatment: another value'
        self.assertEqual(self.check()['status'], 'CANNOT VERIFY')
        self.packet['samples'][0]['characteristics_ch1'] = 'treatment: exposed'
        self.packet['samples'][0]['characteristics_ch1.2'] = 'treatment: other'
        self.assertEqual(self.check()['status'], 'CANNOT VERIFY')

    def test_unknown_requested_arm_or_stratum_cannot_pass_as_a_category(self):
        for token in ('unknown', 'NA', 'not reported'):
            packet = json.loads(json.dumps(self.packet))
            for row in packet['samples'][:2]:
                row['characteristics_ch1'] = 'treatment: ' + token
            packet['target'] = 'treatment: ' + token + ' versus treatment: vehicle within time: late'
            with self.subTest(arm=token):
                self.assertEqual(self.m.verify_requested_comparison(packet, ['A', 'B'], ['C', 'D'])['status'], 'CANNOT VERIFY')
            packet = json.loads(json.dumps(self.packet))
            for row in packet['samples']:
                row['characteristics_ch1.1'] = 'time: ' + token
            packet['target'] = 'treatment: exposed versus treatment: vehicle within time: ' + token
            with self.subTest(stratum=token):
                self.assertEqual(self.m.verify_requested_comparison(packet, ['A', 'B'], ['C', 'D'])['status'], 'CANNOT VERIFY')

    def test_unknown_assigned_stratum_remains_unresolved(self):
        self.packet['samples'][0]['characteristics_ch1.1'] = 'time: unknown'
        self.assertEqual(self.check()['status'], 'CANNOT VERIFY')

    def test_omitted_possible_member_prevents_certifying_complete_membership(self):
        extra = {'sample_id': 'E', 'characteristics_ch1.1': 'time: late'}
        self.packet['samples'].append(extra)
        self.assertEqual(self.check()['status'], 'CANNOT VERIFY')
        extra['characteristics_ch1'] = 'treatment: exposed'
        extra['characteristics_ch1.1'] = 'time: unknown'
        self.assertEqual(self.check()['status'], 'CANNOT VERIFY')

    def test_proved_ineligible_omitted_record_does_not_create_uncertainty(self):
        extra = {'sample_id': 'E', 'characteristics_ch1.1': 'time: early'}
        self.packet['samples'].append(extra)
        self.assertEqual(self.check()['status'], 'PASS')
        extra['characteristics_ch1.1'] = 'time: unknown'
        extra['characteristics_ch1'] = 'treatment: other'
        self.assertEqual(self.check()['status'], 'PASS')

    def test_known_contradiction_is_retained_alongside_unresolved_declarations(self):
        self.packet['samples'][0]['characteristics_ch1'] = 'treatment: vehicle'
        self.packet['samples'][1]['characteristics_ch1'] = 'treatment: unknown'
        self.assertEqual(self.check()['status'], 'FAIL')
        self.packet['samples'][0]['characteristics_ch1'] = 'treatment: unknown'
        self.packet['samples'][0]['characteristics_ch1.1'] = 'time: early'
        self.assertEqual(self.check()['status'], 'FAIL')

    def test_known_eligible_omission_remains_a_failure_with_other_unknown_records(self):
        self.packet['samples'].append({'sample_id': 'E', 'characteristics_ch1.1': 'time: late'})
        result = self.m.verify_requested_comparison(self.packet, ['A'], ['C', 'D'])
        self.assertEqual(result['status'], 'FAIL')
        self.assertEqual([r['sample_id'] for r in result['known_membership_disagreements']], ['B'])

    def test_unknown_platform_observation_propagates_to_final_design_status(self):
        for row in self.packet['samples']:
            row.update(platform_id='P1', organism_ch1='Homo sapiens')
            for index, (label, value) in enumerate([('donor id', row['sample_id']), ('age', '40'),
                    ('sex', 'female'), ('genotype', 'wild type'), ('tissue', 'blood'), ('batch', 'one')], start=2):
                row['characteristics_ch1.' + str(index)] = label + ': ' + value
        self.assertEqual(self.m.verify_source_design(self.packet, ['A', 'B'], ['C', 'D'])['status'], 'PASS')
        for token in (None, '', 'unknown', 'NA', 'not reported'):
            for row in self.packet['samples']:
                row['platform_id'] = token
            with self.subTest(platform=token):
                result = self.m.verify_source_design(self.packet, ['A', 'B'], ['C', 'D'])
                self.assertEqual(result['platform_by_arm']['status'], 'CANNOT VERIFY')
                self.assertEqual(result['status'], 'CANNOT VERIFY')

    def test_unlabelled_prose_and_different_keys_do_not_create_semantic_errors(self):
        for target in ('Exposed specimens versus controls', 'treatment: exposed versus condition: vehicle'):
            self.packet['target'] = target
            with self.subTest(target=target):
                self.assertEqual(self.check()['status'], 'CANNOT VERIFY')
                self.assertNotIn('explicit_source_contradictions', self.check())

    def test_valid_membership_pass_does_not_certify_biological_independence(self):
        self.assertEqual(self.check()['status'], 'PASS')
        observation = self.m.source_field_by_arm(self.packet, ['A', 'B'], ['C', 'D'])
        self.assertEqual(observation['scientific_status'], 'CANNOT VERIFY')
        self.assertEqual(observation['biological_independence_status'], 'CANNOT VERIFY')

    def test_separation_has_exact_evidence_and_assay_record_counts(self):
        result = self.m.source_field_by_arm(self.packet, ['A', 'B'], ['C', 'D'])
        self.assertEqual((result['status'], result['pattern']), ('PASS', 'PERFECT_ARM_SEPARATION'))
        self.assertEqual(result['counts'], {'case': {'P1': 2}, 'control': {'P2': 2}})
        self.assertEqual(len(result['evidence']), 4)
        self.assertEqual(self.m.source_grouping_diagnostics(self.packet, ['A', 'B'], ['C', 'D'])['platform_by_arm'], result)

    def test_mixed_platforms_are_not_automatically_perfect_separation(self):
        self.packet['samples'][1]['platform_id'] = 'P2'
        result = self.m.source_field_by_arm(self.packet, ['A', 'B'], ['C', 'D'])
        self.assertEqual(result['pattern'], 'VARIES_WITHIN_ARM')

    def test_same_platform_is_a_constant_not_an_independence_certificate(self):
        for row in self.packet['samples']:
            row['platform_id'] = 'P1'
        result = self.m.source_field_by_arm(self.packet, ['A', 'B'], ['C', 'D'])
        self.assertEqual(result['pattern'], 'SAME_CONSTANT')
        self.assertEqual(result['biological_independence_status'], 'CANNOT VERIFY')

    def test_missing_platform_prevents_complete_separation_claim(self):
        for value in (None, '', 'unknown', 'NA'):
            self.packet['samples'][0]['platform_id'] = value
            with self.subTest(value=value):
                result = self.m.source_field_by_arm(self.packet, ['A', 'B'], ['C', 'D'])
                self.assertEqual((result['status'], result['pattern']), ('CANNOT VERIFY', 'UNRESOLVED'))
                self.assertEqual(result['unresolved_ids'], ['A'])

    def test_invalid_identity_and_empty_arms_cannot_yield_a_separation_claim(self):
        for case, control in ((['A', 'A'], ['C', 'D']), (['MISSING'], ['C', 'D']), (['A'], ['A'])):
            with self.subTest(case=case, control=control):
                self.assertEqual(self.m.source_field_by_arm(self.packet, case, control)['status'], 'FAIL')
        self.assertEqual(self.m.source_field_by_arm(self.packet, [], ['C', 'D'])['status'], 'CANNOT VERIFY')

    def test_repeated_source_units_do_not_suppress_a_record_level_observation(self):
        before = json.dumps(self.packet, sort_keys=True)
        result = self.m.source_field_by_arm(self.packet, ['A', 'B'], ['C', 'D'])
        self.assertEqual(before, json.dumps(self.packet, sort_keys=True))
        for row in self.packet['samples']:
            row['characteristics_ch1.2'] = 'donor id: shared donor'
        repeated = self.m.source_field_by_arm(self.packet, ['A', 'B'], ['C', 'D'])
        self.assertEqual(repeated['pattern'], result['pattern'])
        self.assertEqual(repeated['biological_independence_status'], 'CANNOT VERIFY')


if __name__ == '__main__':
    unittest.main()
