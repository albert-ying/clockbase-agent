"""Unit tests for the public release export. All fixtures are synthetic."""
import contextlib
import csv
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from audit_system import v3_public_export_2026_09_30 as pe
from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import sha256


def reject_constant(value):
    raise ValueError(f'Non-finite JSON constant: {value}')


def strict_json(path):
    return json.loads(path.read_text(), parse_constant=reject_constant)


class PublicExportTests(unittest.TestCase):
    def test_constant_output_reason_and_action_are_exportable(self):
        record = {'final_class': 'DESCRIPTIVE',
                  'reasons': ['OUTPUTS_CONSTANT_WITHIN_BOTH_ARMS'],
                  'pending_reason': ['outputs constant within both arms'],
                  'actions': [{'action': 'DESCRIPTIVE_ONLY_CONSTANT_OUTPUTS'}]}
        self.assertEqual(pe.undefined_codes(record), set())

    def setUp(self):
        self.folder = tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR'))
        self.addCleanup(self.folder.cleanup)
        self.root = Path(self.folder.name)
        self.final = self.root / 'verifier' / 'final'
        self.verifier = self.root / 'verifier'
        self.out = self.root / 'public-export'
        roots = {label: self.root / label for label in pe.ROOTS}
        for path in roots.values():
            path.mkdir()
        self.root_patch = mock.patch.dict(pe.ROOTS, roots, clear=True)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)
        self.predictions = roots['geo'] / 'GSE1' / 'predictions.csv'
        self.key = {'modality': 'human DNAm', 'gse_id': 'GSE1', 'comparison_name': 'A', 'condition_name': 'Drug A'}
        self.answers = [{'contrast_type': 'INTERVENTION', 'relevance_score': 7},
                        {'contrast_type': 'INTERVENTION', 'relevance_score': 9}]
        self.classification = {**self.key, 'status': 'AGREED', 'contrast_type': 'INTERVENTION',
                               'intervention_class': 'DRUG', 'manipulation': 'TREATMENT', 'canonical_target': 'drug a',
                               'relevance_score': 8.0, 'relevance_score_abs_diff': 2.0, 'answers': self.answers,
                               'receipts': [{'model': 'gp-6-as'}, {'model': 'cl-op-5'}]}
        self.audit_answers = [{'orientation': 'CASE_IS_CONDITION', 'specimen': 'ORGANISM_DIRECT', 'tumor': 'NO'}] * 2
        self.audit = {**self.key, 'status': 'RANKABLE', 'reasons': [], 'orientation': 'CASE_IS_CONDITION',
                      'specimen': 'ORGANISM_DIRECT', 'tumor': 'NO', 'answers': self.audit_answers}
        self.write_jsonl(self.verifier / 'classification-001' / 'comparisons.jsonl', [self.classification])
        self.write_jsonl(self.verifier / 'design-audit-001' / 'comparisons.jsonl', [self.audit])
        self.write_jsonl(self.verifier / 'design-audit-002' / 'comparisons.jsonl',
                         [{**self.audit, 'status': 'NOT_RANKABLE', 'reasons': ['TUMOR_YES']},
                          {**self.audit, 'comparison_name': 'B', 'status': 'NOT_RANKABLE', 'tumor': 'YES',
                           'reasons': ['TUMOR_YES']}])
        repairs = [self.repair(0), self.repair(1, clock='PhenoAge'), self.repair(2, comparison='B'),
                   self.repair(3, comparison='C')]
        repairs.extend([{**self.repair(4), 'final_class': 'DUPLICATE', 'duplicate_of': 0},
                        {**self.repair(5), 'recomputed': {'n_case': 3, 'n_control': 3, 'mean_diff': 2.0, 'p': 0.01}},
                        {**self.repair(6), 'final_class': 'PENDING_DESIGN'}, self.repair(7)])
        self.write_jsonl(self.final / 'repairs-dnam' / 'repairs.jsonl', repairs)
        self.write_jsonl(self.final / 'repairs-mm_rna' / 'repairs.jsonl',
                         [self.repair(i, gse=gse) for i, gse in enumerate([float('nan'), None, '', 'nan'])])
        self.write_jsonl(self.final / 'merged-dnam' / 'rows.jsonl',
                         [{'row_ordinal': 0, 'source': 'design_reconstruction',
                           'design_rules': {'path': str(roots['inputs'] / 'rules.json')}, 'historical_status': 'UNRECOVERED'}])
        self.write_jsonl(self.final / 'verdicts-dnam' / 'row_verdicts.jsonl',
                         [{'row_ordinal': 0, 'evidence': {'prediction_source': str(self.predictions),
                                                        'score': float('-inf')}}])
        self.write_jsonl(self.verifier / 'link-dnam-001' / 'row_links.jsonl',
                         [{'row_ordinal': 0, 'link_status': 'LINKED',
                           'link_traces': [str(roots['repo'] / 'trace.json')],
                           'link': {'trace': str(roots['verifier'] / 'trace.json'), 'test_line': 12}}])
        ranking = self.final / 'comparison-ranking'
        ranking.mkdir(parents=True)
        ranking_rows = [{**self.key, 'clock': 'Horvathv1', 'hedges_g': -2.0, 'g_ci_low': -3.0,
                         'g_ci_high': -1.0, 'ci_level': 0.975, 'selection_k': 2.0, 'direction_v3': 'age-lowering',
                         'v3_final_class': 'VERIFIED', 'v3_n_case': 3.0, 'v3_n_control': 3.0,
                         'v3_mean_diff': -2.000000000000001, 'v3_p': 0.01},
                        {**self.key, 'comparison_name': 'D', 'hedges_g': 1.5, 'direction_v3': 'age-raising'},
                        {**self.key, 'comparison_name': 'E', 'hedges_g': -1.2, 'direction_v3': 'age-lowering'}]
        with (ranking / 'ranked_comparisons_all.csv').open('w', newline='') as sink:
            writer = csv.DictWriter(sink, fieldnames=list(ranking_rows[0]))
            writer.writeheader()
            writer.writerows(ranking_rows)
        self.write_json(ranking / 'summary.json', {'score': float('nan')})
        self.write_json(ranking / 'nested' / 'detail.json', {'score': float('inf')})
        (ranking / 'table1_pro_aging.csv').write_text('comparison_name\nD\n')
        self.write_json(self.final / 'error-census' / 'census.json', {'missing': float('inf')})
        self.write_json(self.final / 'published-status' / 'summary.json', {'missing': float('-inf')})
        (self.final / 'published-status' / 'example_v3_status.csv').write_text('status\nVERIFIED\n')
        self.write_jsonl(self.verifier / 'design-audit-003' / 'comparisons.jsonl',
                         [{**self.key, 'comparison_name': 'Z', 'status': 'INCOMPLETE', 'score': float('nan')}])
        self.write_jsonl(self.verifier / 'classification-001' / 'comparisons.jsonl',
                         [self.classification, {**self.key, 'comparison_name': 'Z', 'relevance_score': float('inf')}])

    @staticmethod
    def write_json(path, value):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(value))  # source fixtures deliberately include legacy NaN/Infinity

    @staticmethod
    def write_jsonl(path, rows):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(''.join(json.dumps(row) + '\n' for row in rows))

    def repair(self, ordinal, gse='GSE1', comparison='A', clock='Horvathv1'):
        return {'row_ordinal': ordinal, 'gse_id': gse, 'comparison_name': comparison, 'condition_name': 'Drug A',
                'clock': clock, 'original': {'pvalue': float('nan'), 'effect_size': float('inf')},
                'original_verdict': 'PASS', 'reasons': [],
                'actions': [{'action': 'RECOMPUTE_DECLARED_ESTIMAND', 'source': str(self.predictions)}], 'final_class': 'VERIFIED',
                'recomputed': {'n_case': 3, 'n_control': 3, 'mean_diff': -2.0, 'p': 0.01, 'score': float('-inf')},
                'inputs': {'predictions': {'path': str(self.predictions)},
                           'case_gsms': ['GSM1', 'GSM2', 'GSM3'], 'control_gsms': ['GSM4', 'GSM5', 'GSM6']}}

    def run_export(self):
        with contextlib.redirect_stdout(io.StringIO()):
            pe.run(self.final, self.verifier, self.out)

    def test_every_json_and_jsonl_is_strict(self):
        self.run_export()
        for path in self.out.rglob('*'):
            if path.suffix == '.json':
                strict_json(path)
            elif path.suffix == '.jsonl':
                for line in path.read_text().splitlines():
                    json.loads(line, parse_constant=reject_constant)
        row = strict_json(self.out / 'studies/dnam/GSE1.json')['rows'][0]
        self.assertIsNone(row['historical']['pvalue'])
        self.assertIsNone(row['historical']['effect_size'])
        self.assertIsNone(row['recomputed']['score'])
        self.assertIsNone(row['evidence']['score'])
        self.assertEqual(strict_json(self.out / 'ranking/nested/detail.json'), {'score': None})

    def test_missing_study_id_is_grouped_under_no_study_id(self):
        self.run_export()
        study = strict_json(self.out / 'studies/mm_rna/NO_STUDY_ID.json')
        self.assertIsNone(study['gse_id'])
        self.assertEqual(len(study['rows']), 4)
        self.assertFalse(list(self.out.rglob('nan.json')))

    def test_paths_are_portable_through_nested_evidence(self):
        self.run_export()
        row = strict_json(self.out / 'studies/dnam/GSE1.json')['rows'][0]
        self.assertEqual(row['predictions']['path'], 'geo:GSE1/predictions.csv')
        self.assertEqual(row['actions'][0]['source'], 'geo:GSE1/predictions.csv')
        self.assertEqual(row['evidence']['prediction_source'], 'geo:GSE1/predictions.csv')
        self.assertEqual(row['replay']['traces'], ['repo:trace.json'])
        self.assertEqual(row['replay']['test']['trace'], 'verifier:trace.json')
        self.assertEqual(row['design_rules']['path'], 'inputs:rules.json')
        self.assertNotIn(str(self.root), (self.out / 'studies/dnam/GSE1.json').read_text())
        for label, root in pe.ROOTS.items():
            self.assertEqual(pe.release_value(root / 'example.csv'), f'{label}:example.csv')
        alias = self.root / 'geo-alias'
        alias.symlink_to(pe.ROOTS['geo'], target_is_directory=True)
        self.assertEqual(pe.release_value(str(alias / 'example.csv')), 'geo:example.csv')

    def test_ranking_evidence_and_first_design_audit_wins(self):
        self.run_export()
        rows = strict_json(self.out / 'studies/dnam/GSE1.json')['rows']
        evidence = rows[0]['ranking_evidence']
        self.assertEqual(evidence['classification'], {'status': 'AGREED', 'contrast_type': 'INTERVENTION',
                                                      'intervention_class': 'DRUG', 'manipulation': 'TREATMENT',
                                                      'canonical_target': 'drug a', 'relevance_score': 8.0,
                                                      'relevance_score_abs_diff': 2.0, 'answers': self.answers,
                                                      'receipts': [{'model': 'gp-6-as'}, {'model': 'cl-op-5'}],
                                                      'relevance_codes_agreed': None})
        self.assertEqual(evidence['design_audit'], {'run': 'design-audit-001',
                                                  **{k: self.audit.get(k) for k in pe.DESIGN_AUDIT_FIELDS}})
        self.assertEqual(evidence['ranking'], {'ranked': True, 'rank': 1, 'rank_within_direction': 1,
                                               'hedges_g': -2.0, 'g_ci_low': -3.0, 'g_ci_high': -1.0,
                                               'ci_level': 0.975, 'selection_k': 2.0,
                                               'direction_v3': 'age-lowering', 'clock': 'Horvathv1'})
        unranked = {'ranked': False, 'rank': None, 'rank_within_direction': None, 'hedges_g': None,
                    'g_ci_low': None, 'g_ci_high': None, 'ci_level': None, 'selection_k': None,
                    'direction_v3': None, 'clock': None}
        self.assertEqual(rows[1]['ranking_evidence']['ranking'], unranked)  # unselected clock
        self.assertEqual(rows[2]['ranking_evidence']['ranking'], unranked)
        self.assertIsNone(rows[2]['ranking_evidence']['classification'])
        self.assertEqual(rows[2]['ranking_evidence']['design_audit']['reasons'], ['TUMOR_YES'])
        self.assertEqual(rows[2]['ranking_evidence']['design_audit']['run'], 'design-audit-002')
        self.assertIsNone(rows[3]['ranking_evidence']['design_audit'])
        for row in rows[4:]:
            self.assertEqual(row['ranking_evidence']['ranking'], unranked)
        self.assertEqual(sum(row['ranking_evidence']['ranking']['ranked'] for row in rows), 1)
        ranked = pe.ranking_evidence(self.final, self.verifier)[2]
        third = ranked[(pe.comparison_key({**self.key, 'comparison_name': 'E'}), None)]
        self.assertEqual((third['rank'], third['rank_within_direction']), (3, 2))

    def test_missing_ranking_table_is_null(self):
        (self.final / 'comparison-ranking/ranked_comparisons_all.csv').unlink()
        self.run_export()
        row = strict_json(self.out / 'studies/dnam/GSE1.json')['rows'][0]
        self.assertIsNone(row['ranking_evidence']['ranking'])

    def test_manifest_covers_every_export_and_codes_document_roots(self):
        self.run_export()
        manifest = strict_json(self.out / 'manifest.json')
        files = {p.relative_to(self.out).as_posix() for p in self.out.rglob('*')
                 if p.is_file() and p != self.out / 'manifest.json'}
        self.assertEqual(set(manifest['files']), files)
        for name, digest in manifest['files'].items():
            self.assertEqual(digest, sha256(self.out / name))
        self.assertTrue({'README.md', 'ranking/ranked_comparisons_all.csv', 'ranking/table1_pro_aging.csv',
                         'classification/comparisons.jsonl', 'design-audit/design-audit-001.jsonl',
                         'design-audit/design-audit-002.jsonl', 'design-audit/design-audit-003.jsonl'} <= files)
        codes = strict_json(self.out / 'codes.json')
        self.assertEqual(codes['model_aliases'], pe.MODEL_ALIASES)
        self.assertEqual(codes['portable_paths']['roots'], list(pe.ROOTS))
        self.assertIn('root:relative/path', codes['portable_paths']['format'])
        self.assertIn('rank_within_direction', (self.out / 'README.md').read_text())

    def test_records_use_defined_codes_and_explain_every_check(self):
        self.run_export()
        codes = strict_json(self.out / 'codes.json')
        definitions = {name: codes[name] for name in pe.code_definitions()}
        for path in self.out.glob('studies/*/*.json'):
            study = strict_json(path)
            self.assertFalse(pe.undefined_codes(study, definitions), str(path))
            for row in study['rows']:
                self.assertEqual(set(row['checks']), set(pe.CHECKS))
                for name, check in row['checks'].items():
                    self.assertIn(check['outcome'], pe.CHECK_OUTCOMES)
                    self.assertIn('evidence', check)
                    self.assertIn('reason_codes', check)
                    self.assertIn(name, codes['checks'])
        self.assertEqual(pe.undefined_codes({'soft_flags': ['UNDEFINED_SOFT_CODE'],
                                            'actions': [{'action': 'UNDEFINED_ACTION'}],
                                            'adjudication': {'issues_agreed': ['UNDEFINED_ISSUE']}}),
                         {'UNDEFINED_SOFT_CODE', 'UNDEFINED_ACTION', 'UNDEFINED_ISSUE'})

    def test_export_rejects_an_undefined_record_code(self):
        path = self.final / 'repairs-dnam/repairs.jsonl'
        rows = pe.load_jsonl(path)
        rows[0]['pending_reason'] = ['UNDOCUMENTED_PENDING_REASON']
        self.write_jsonl(path, rows)
        with self.assertRaisesRegex(ValueError, 'UNDOCUMENTED_PENDING_REASON'):
            self.run_export()

    def check_row(self, evidence=None, reasons=(), modality='dnam', member=None, **updates):
        row = {**self.repair(0), 'historical': {'pvalue': .01, 'fdr': .02},
               'evidence': evidence, 'reasons': list(reasons), 'ranking_evidence': {'design_audit': None}}
        row.update(updates)
        return pe.row_checks(row, member or {}, modality)

    def test_successes_and_stopped_checks_are_distinguished(self):
        evidence = {'membership_status': 'UNIQUE', 'membership_source': 'recovery', 'test': 'Welch',
                    'n_case': 3, 'n_control': 3, 'effect_definitions_matched': ['cohen_d'],
                    'label_orientation': 'SUPPORTED', 'profile_r_reference_min': .8,
                    'clock_coverage_min': .9, 'clock_target': 'chronological age (years)'}
        checks = self.check_row(evidence, ('MEMBERSHIP_RECOVERED',), member={'status': 'UNIQUE', 'source': 'recovery'})
        for name in ('membership_recovery', 'test_type', 'stored_p_value', 'effect_reproduction',
                     'arm_orientation', 'arm_size', 'assay', 'dnam_input_integrity', 'dnam_clock_coverage',
                     'clock_target_age_input', 'fdr_vs_p', 'duplicate'):
            self.assertEqual(checks[name]['outcome'], 'PASS', name)
        self.assertEqual(checks['tissue_domain']['outcome'], 'NOT_APPLICABLE')
        self.assertEqual(checks['rna_input']['outcome'], 'NOT_APPLICABLE')
        self.assertEqual(checks['design_concerns']['outcome'], 'NOT_EVALUATED')
        stopped = self.check_row({'membership_status': 'NONE'}, ('MEMBERSHIP_UNRECOVERED',),
                                 final_class='PENDING_DESIGN', recomputed=None, actions=[])
        self.assertEqual(stopped['membership_recovery']['outcome'], 'FAIL')
        for name in ('effect_reproduction', 'arm_orientation', 'arm_size', 'assay', 'dnam_input_integrity',
                     'dnam_clock_coverage', 'clock_target_age_input', 'duplicate'):
            self.assertEqual(stopped[name]['outcome'], 'NOT_EVALUATED', name)

    def test_failure_soft_flag_and_recorded_repair_remain_visible(self):
        actions = [{'action': 'SWAP_ARMS_CONFIRMED_INVERSION'},
                   {'action': 'USE_REPAIRED_INPUT', 'clock_verdict': {'verdict': 'PROVISIONAL',
                                                                    'reasons': ['COVERAGE_0.60_TO_0.80']}},
                   {'action': 'RECOMPUTE_DECLARED_ESTIMAND'}]
        checks = self.check_row({'membership_status': 'UNIQUE', 'label_orientation': 'INVERTED',
                                 'clock_coverage_min': .3, 'clock_target': 'phenotypic age (years)'},
                                ('LABEL_ORIENTATION_INVERTED', 'INPUT_EMPTY', 'CLOCK_COVERAGE_LOW',
                                 'TISSUE_OUT_OF_DOMAIN', 'EFFECT_NOT_REPRODUCED'),
                                actions=actions, clock='PhenoAge', soft_flags=['ADJUDICATED_TECHNICAL_REPLICATES'])
        for name in ('arm_orientation', 'dnam_input_integrity', 'dnam_clock_coverage', 'effect_reproduction'):
            self.assertEqual(checks[name]['outcome'], 'FAIL', name)
        self.assertEqual(checks['arm_orientation']['after_repair']['outcome'], 'PASS')
        self.assertEqual(checks['dnam_clock_coverage']['after_repair']['outcome'], 'FLAG')
        self.assertEqual(checks['tissue_domain']['outcome'], 'FLAG')
        self.assertEqual(checks['design_concerns']['outcome'], 'FLAG')
        self.assertNotIn('after_repair', checks['effect_reproduction'])
        self.assertIn('replacement', checks['effect_reproduction'])

    def test_design_reconstruction_does_not_validate_historical_statistics(self):
        member = {'status': 'UNIQUE', 'source': 'design_reconstruction', 'historical_status': 'STORED_P_UNUSABLE'}
        checks = self.check_row({'membership_status': 'UNIQUE'}, ('DESIGN_RECONSTRUCTED',), member=member)
        self.assertEqual(checks['membership_recovery']['outcome'], 'PASS')
        self.assertEqual(checks['stored_p_value']['outcome'], 'FAIL')
        self.assertEqual(checks['effect_reproduction']['outcome'], 'NOT_APPLICABLE')
        member['historical_status'] = 'NOT_TWO_SAMPLE_T'
        checks = self.check_row({'membership_status': 'UNIQUE'}, ('DESIGN_RECONSTRUCTED',), member=member)
        self.assertEqual(checks['test_type']['outcome'], 'FAIL')
        self.assertEqual(checks['stored_p_value']['outcome'], 'NOT_EVALUATED')

    def test_profile_unscreened_and_rna_applicability(self):
        checks = self.check_row({'membership_status': 'UNIQUE', 'profile_r_reference_min': .8}, ('PROFILE_UNSCREENED',))
        self.assertEqual(checks['dnam_input_integrity']['outcome'], 'NOT_EVALUATED')
        checks = self.check_row({'membership_status': 'NONE'}, ('RNA_INPUT_ARTIFACT',), modality='hs_rna')
        self.assertEqual(checks['rna_input']['outcome'], 'FAIL')
        for name in ('assay', 'dnam_input_integrity', 'dnam_clock_coverage', 'tissue_domain'):
            self.assertEqual(checks[name]['outcome'], 'NOT_APPLICABLE')

    def test_atlas_prefixes_historical_columns_and_documents_all_columns(self):
        path = self.final / 'atlas-v3' / 'atlas_v3_combined.csv'
        path.parent.mkdir()
        source = 'gse_id,condition_name,clock,atlas_row_ordinal,pvalue,fdr,effect_size,direction,log2FoldChange,v3_final_class,v3_p\nGSE1,A,Horvathv1,0,.01,.02,-2,-1,-2,DESCRIPTIVE,\n'
        path.write_text(source)
        self.run_export()
        self.assertEqual(path.read_text(), source)
        with (self.out / 'atlas_v3_combined.csv').open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]['historical_pvalue'], '.01')
        self.assertEqual(rows[0]['v3_p'], '')
        self.assertEqual(rows[0]['v3_final_class'], 'DESCRIPTIVE')
        self.assertEqual(rows[0]['gse_id'], 'GSE1')
        for column in rows[0]:
            self.assertIn(f'| {column} |', (self.out / 'README.md').read_text())
        self.assertTrue({'historical_pvalue', 'historical_fdr', 'historical_effect_size',
                         'historical_direction', 'historical_log2FoldChange'} <= rows[0].keys())
        ranking = (self.out / 'ranking/ranked_comparisons_all.csv').read_bytes()
        self.assertEqual(ranking, (self.final / 'comparison-ranking/ranked_comparisons_all.csv').read_bytes())


class FinalRunCodeTests(unittest.TestCase):
    def test_final_run_records_have_definitions_for_every_emitted_code(self):
        final = paths.VERIFIER_ROOT / 'final-2026-10-04' / 'public-export'
        if not final.exists():
            self.skipTest('final-2026-10-04 public export is not installed')
        for path in final.glob('studies/*/*.json'):
            missing = pe.undefined_codes(strict_json(path))
            if missing:
                self.fail(f'{path.name}: undefined codes {sorted(missing)}')
        for path in final.rglob('*.jsonl'):
            for line in path.read_text().splitlines():
                missing = pe.undefined_codes(json.loads(line, parse_constant=reject_constant))
                if missing:
                    self.fail(f'{path.name}: undefined codes {sorted(missing)}')


if __name__ == '__main__':
    unittest.main()
