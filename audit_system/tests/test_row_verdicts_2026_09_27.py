"""Unit tests for verifier-v3 row verdicts. All fixtures are synthetic."""
import unittest
from unittest import mock

from audit_system import row_verdicts_2026_09_27 as rv


def recovered(**overrides):
    record = {'status': 'UNIQUE', 'case_gsms': ['GSM1', 'GSM2', 'GSM3'], 'control_gsms': ['GSM4', 'GSM5', 'GSM6'],
              'n_case': 3, 'n_control': 3, 'test': 'welch', 'field': 'treatment',
              'case_level': 'resveratrol', 'control_level': 'DMSO', 'stratum': None,
              'effect_definitions_matched': ['cohen_d_pooled'], 'prediction_source': 'clock_result'}
    record.update(overrides)
    return record


ROW = {'gse_id': 'GSE1', 'comparison_name': 'Resveratrol vs DMSO', 'condition_name': 'Resveratrol',
       'clock': 'Hazard_Mouse_All_All_scaleddiff_EN_tAge', 'fdr_below_raw_p': False}
RNA = {'modality': 'mm_rna', 'input_artifact': {}, 'sample_text': lambda gse: {}}


class RowVerdictTests(unittest.TestCase):
    def test_label_orientation(self):
        self.assertEqual(rv.label_orientation('Resveratrol vs DMSO', 'Resveratrol', 'resveratrol 10uM', 'DMSO', None),
                         'SUPPORTED')
        self.assertEqual(rv.label_orientation('Resveratrol vs DMSO', 'Resveratrol', 'DMSO', 'resveratrol', None),
                         'INVERTED')
        self.assertEqual(rv.label_orientation('Group A vs Group B', 'A', 'x1', 'x2', None), 'UNCONFIRMED')

    def test_clean_recovered_row_is_verified(self):
        verdict, reasons, evidence = rv.verdict_for(ROW, recovered(), RNA)
        self.assertEqual(verdict, 'VERIFIED')
        self.assertEqual(reasons, ['MEMBERSHIP_RECOVERED'])
        self.assertEqual(evidence['label_orientation'], 'SUPPORTED')

    def test_hard_failures_withhold(self):
        cases = [(recovered(n_case=2), 'ARM_TOO_SMALL'),
                 (recovered(effect_definitions_matched=[]), 'EFFECT_NOT_REPRODUCED'),
                 (recovered(case_level='DMSO', control_level='resveratrol'), 'LABEL_ORIENTATION_INVERTED'),
                 ({'status': 'NONE'}, 'MEMBERSHIP_UNRECOVERED'),
                 ({'status': 'NOT_TWO_SAMPLE_T'}, 'TEST_NOT_RECONSTRUCTIBLE'),
                 (None, 'MEMBERSHIP_UNRECOVERED')]
        for record, code in cases:
            verdict, reasons, _ = rv.verdict_for(ROW, record, RNA)
            self.assertEqual(verdict, 'WITHHELD', code)
            self.assertIn(code, reasons)
        verdict, reasons, _ = rv.verdict_for({**ROW, 'fdr_below_raw_p': True}, recovered(), RNA)
        self.assertEqual((verdict, 'FDR_BELOW_RAW_P' in reasons), ('WITHHELD', True))
        verdict, reasons, _ = rv.verdict_for(ROW, recovered(), {**RNA, 'input_artifact': {'GSE1': 'RNA_INPUT_ARTIFACT'}})
        self.assertIn('RNA_INPUT_ARTIFACT', reasons)

    def test_dnam_assay_coverage_and_tissue(self):
        row = {**ROW, 'clock': 'PhenoAge'}
        texts = {'GSM1': 'MeDIP DNA from liver', 'GSM4': 'liver'}
        context = {'modality': 'dnam', 'input_artifact': {}, 'sample_text': lambda gse: texts}
        with mock.patch.object(rv, 'dnam_coverage', return_value=0.95), \
                mock.patch.object(rv, 'sample_tissue', return_value='liver'):
            verdict, reasons, _ = rv.verdict_for(row, recovered(), context)
        self.assertEqual(verdict, 'WITHHELD')
        self.assertIn('NON_STANDARD_ASSAY', reasons)
        self.assertIn('TISSUE_OUT_OF_DOMAIN', reasons)
        context['sample_text'] = lambda gse: {}
        with mock.patch.object(rv, 'dnam_coverage', return_value=0.7), \
                mock.patch.object(rv, 'sample_tissue', return_value='whole blood'):
            verdict, reasons, _ = rv.verdict_for(row, recovered(), context)
        self.assertEqual((verdict, 'CLOCK_COVERAGE_MID' in reasons), ('PROVISIONAL', True))
        with mock.patch.object(rv, 'dnam_coverage', return_value=None), \
                mock.patch.object(rv, 'sample_tissue', return_value='whole blood'):
            verdict, reasons, _ = rv.verdict_for(row, recovered(), context)
        self.assertEqual((verdict, 'COVERAGE_UNAVAILABLE' in reasons), ('WITHHELD', True))


    def test_grimage_rows_carry_the_age_input_flag(self):
        context = {'modality': 'dnam', 'input_artifact': {}, 'sample_text': lambda gse: {}}
        with mock.patch.object(rv, 'dnam_coverage', return_value=0.95), \
                mock.patch.object(rv, 'sample_tissue', return_value='whole blood'):
            verdict, reasons, _ = rv.verdict_for({**ROW, 'clock': 'GrimAgeV2'}, recovered(), context)
            self.assertEqual((verdict, 'CLOCK_AGE_INPUT_FROM_HORVATH' in reasons), ('PROVISIONAL', True))
            verdict, reasons, _ = rv.verdict_for({**ROW, 'clock': 'PhenoAge'}, recovered(), context)
            self.assertEqual(verdict, 'VERIFIED')

    def test_detached_profiles_withhold(self):
        row = {**ROW, 'clock': 'Horvathv1'}
        context = {'modality': 'dnam', 'input_artifact': {}, 'sample_text': lambda gse: {},
                   'integrity': {'GSM1': 0.9, 'GSM2': 0.35, 'GSM3': 0.93, 'GSM4': -0.002, 'GSM5': 0.9}}
        with mock.patch.object(rv, 'dnam_coverage', return_value=0.95), \
                mock.patch.object(rv, 'sample_tissue', return_value='liver'):
            verdict, reasons, evidence = rv.verdict_for(row, recovered(), context)
        self.assertEqual(verdict, 'WITHHELD')
        self.assertIn('INPUT_PROFILE_DETACHED', reasons)
        self.assertIn('PROFILE_UNSCREENED', reasons)  # GSM6 absent from the screen
        self.assertIn('PROFILE_ATYPICAL', reasons)  # GSM2 at 0.35
        self.assertEqual(evidence['atypical_profile_gsms'], ['GSM2'])
        self.assertEqual((evidence['detached_gsms'], evidence['profile_r_reference_min']), (['GSM4'], -0.002))

    def test_assay_regex_treats_underscore_as_separator(self):
        for text, expected in [('tissue_1_oxBS', True), ('DMSO_7d_oxBS_Rep1', True), ('H3K4me3_ChIP_input', True),
                               ('5hmCIP_450K', True), ('tissue_1_BS', False), ('ChIP: 6264509100', False),
                               ('noIP_450K', False), ('whole blood', False)]:
            self.assertEqual(bool(rv.NON_STANDARD_ASSAY.search(text)), expected, text)


    def test_design_reconstructed_rows_get_member_checks_but_no_effect_check(self):
        design = recovered(source='design_reconstruction', effect_definitions_matched=[], case_level=None,
                           control_level=None)
        verdict, reasons, _ = rv.verdict_for(ROW, design, RNA)
        self.assertEqual(verdict, 'VERIFIED')
        self.assertIn('DESIGN_RECONSTRUCTED', reasons)
        verdict, reasons, _ = rv.verdict_for(ROW, {**design, 'n_case': 2}, RNA)
        self.assertEqual((verdict, 'ARM_TOO_SMALL' in reasons), ('WITHHELD', True))


if __name__ == '__main__':
    unittest.main()
