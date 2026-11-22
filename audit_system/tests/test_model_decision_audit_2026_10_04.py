"""Saved decision audits with small synthetic fixtures; no model requests or empirical measurements."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from audit_system import model_decision_audit_2026_10_04 as audit


class ModelDecisionAuditTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.network = patch('urllib.request.urlopen', side_effect=AssertionError('Audit must not send requests'))
        self.network.start()
        self.addCleanup(self.network.stop)

    def write(self, folder, name, records):
        folder.mkdir(parents=True, exist_ok=True)
        (folder / name).write_text(''.join(json.dumps(r) + '\n' for r in records))

    def receipts(self):
        return [{'model': m, 'status': 'OK', 'request_sha256': m + '-hash', 'cached': False} for m in ('m1', 'm2')]

    def ledger(self, answers, wrapped=False):
        return [{'model_requested': m, 'status': 'OK', 'request_sha256': m + '-hash',
                 'response_text': json.dumps({'comparisons': [a]} if wrapped else a),
                 'parsed': {'ignored': 'synthetic deliberately incorrect cached parse'}}
                for m, a in zip(('m1', 'm2'), answers)]

    def label_answer(self, orientation='CORRECT'):
        return {'arm_a_is_condition': 'YES', 'arm_b_is_appropriate_control': 'YES',
                'orientation': orientation, 'issues': ['OTHER']}

    def label_record(self, answers):
        return {'row_ordinal': 7, 'decision': audit.adjudication.decide(answers), 'judgments': answers,
                'issues_agreed': ['OTHER'], 'receipts': self.receipts()}

    def design_answer(self, cid=9, specimen='ORGANISM_DIRECT'):
        return {'id': cid, 'orientation': 'CASE_IS_CONDITION', 'specimen': specimen, 'tumor': 'NO',
                'confounds': [], 'evidence': ['synthetic specimen field']}

    def design_record(self, answers):
        return {'modality': 'human DNAm', 'gse_id': 'synthetic-study', 'comparison_name': 'synthetic case vs control',
                'condition_name': 'synthetic case', 'receipts': self.receipts(), **audit.design_audit.reconcile(answers)}

    def candidates(self, folder, record, cid=9):
        lines = [','.join(audit.KEY)]
        lines += [f'human DNAm,synthetic-other-{i},unused,unused' for i in range(cid)]
        lines += [','.join(record[k] for k in audit.KEY)]
        (folder / 'candidates.csv').write_text('\n'.join(lines) + '\n')

    def test_adjudication_reparses_raw_and_detects_tampered_decision(self):
        folder = self.root / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        record = self.label_record(answers)
        tampered = {**record, 'row_ordinal': 8, 'decision': 'LABEL_INVERTED_CONFIRMED', 'issues_agreed': []}
        self.write(folder, 'ledger.jsonl', self.ledger(answers))
        self.write(folder, 'decisions-all.jsonl', [record, tampered])
        report = audit.audit_run(folder, 'adjudication')
        self.assertEqual((report['records'], report['receipts_found'], report['rederived_identically'], report['differing']),
                         (2, 4, 1, 1))
        self.assertEqual(set(report['differing_examples'][0]['fields']), {'decision', 'issues_agreed'})

    def test_missing_and_wrong_model_receipts_are_not_identical(self):
        folder = self.root / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        ledger = self.ledger(answers)
        ledger[1]['model_requested'] = 'wrong-model'
        self.write(folder, 'ledger.jsonl', ledger)
        self.write(folder, 'decisions-all.jsonl', [self.label_record(answers)])
        report = audit.audit_run(folder, 'adjudication')
        self.assertEqual((report['receipts_found'], report['missing_receipts'], report['records_missing_receipts']), (1, 1, 1))
        self.assertEqual(report['rederived_identically'], 0)

    def test_unparseable_response_reproduces_incomplete_decision(self):
        folder = self.root / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        ledger = self.ledger(answers)
        ledger[1].update(status='UNPARSEABLE', response_text='synthetic broken response')
        record = self.label_record(answers)
        record.update(decision='ADJUDICATION_INCOMPLETE', judgments=answers[:1], issues_agreed=[])
        record['receipts'][1]['status'] = 'UNPARSEABLE'
        self.write(folder, 'ledger.jsonl', ledger)
        self.write(folder, 'decisions-all.jsonl', [record])
        report = audit.audit_run(folder, 'adjudication')
        self.assertEqual((report['unparseable_responses'], report['rederived_identically']), (1, 1))

    def test_retry_does_not_replace_historical_response(self):
        folder = self.root / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        ledger = self.ledger(answers)
        retry = {**ledger[0], 'response_text': json.dumps(self.label_answer('INVERTED'))}
        self.write(folder, 'ledger.jsonl', ledger + [retry])
        self.write(folder, 'decisions-all.jsonl', [self.label_record(answers)])
        report = audit.audit_run(folder, 'adjudication')
        self.assertEqual((report['ambiguous_receipts'], report['rederived_identically']), (1, 1))

    def test_missing_raw_response_is_not_replaced_by_cached_parse(self):
        folder = self.root / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        ledger = self.ledger(answers)
        ledger[1].pop('response_text')
        ledger[1]['parsed'] = answers[1]
        self.write(folder, 'ledger.jsonl', ledger)
        self.write(folder, 'decisions-all.jsonl', [self.label_record(answers)])
        report = audit.audit_run(folder, 'adjudication')
        self.assertEqual((report['missing_receipts'], report['records_missing_receipts']), (1, 1))
        self.assertEqual(report['rederived_identically'], 0)

    def test_stored_empty_response_is_found_but_unparseable(self):
        folder = self.root / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        ledger = self.ledger(answers)
        ledger[1].update(status='UNPARSEABLE', response_text='')
        record = self.label_record(answers)
        record.update(decision='ADJUDICATION_INCOMPLETE', judgments=answers[:1], issues_agreed=[])
        record['receipts'][1]['status'] = 'UNPARSEABLE'
        self.write(folder, 'ledger.jsonl', ledger)
        self.write(folder, 'decisions-all.jsonl', [record])
        report = audit.audit_run(folder, 'adjudication')
        self.assertEqual((report['receipts_found'], report['missing_receipts'], report['unparseable_responses'],
                          report['rederived_identically']), (2, 0, 1, 1))

    def test_receipt_status_tampering_is_reported(self):
        folder = self.root / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        record = self.label_record(answers)
        record['receipts'][0]['status'] = 'UNPARSEABLE'
        self.write(folder, 'ledger.jsonl', self.ledger(answers))
        self.write(folder, 'decisions-all.jsonl', [record])
        report = audit.audit_run(folder, 'adjudication')
        self.assertEqual(report['differing'], 1)
        self.assertIn('receipt_0_receipt_status', report['differing_examples'][0]['fields'])

    def test_design_audit_uses_candidate_ids_and_detects_tampering(self):
        folder = self.root / 'design-audit-001'
        answers = [self.design_answer(), self.design_answer()]
        record = self.design_record(answers)
        tampered = {**record, 'reasons': ['TUMOR_YES'], 'orientation': 'REVERSED'}
        self.write(folder, 'ledger.jsonl', self.ledger(answers, wrapped=True))
        self.write(folder, 'comparisons.jsonl', [record, tampered])
        self.candidates(folder, record)
        report = audit.audit_run(folder, 'design_audit')
        self.assertEqual((report['rederived_identically'], report['differing']), (1, 1))
        self.assertEqual(set(report['differing_examples'][0]['fields']), {'reasons', 'orientation'})

    def test_design_audit_salvages_raw_fields(self):
        folder = self.root / 'design-audit-001'
        text = ('{"comparisons": [{"id": 9, "orientation": "CASE_IS_CONDITION", "specimen": "EX_VIVO_CULTURE", '
                '"tumor": "NO", "confounds": [], "evidence": ["synthetic "broken" quote"]}')
        parsed = audit.design_audit.parse(text)[9]
        answers = [parsed, parsed]
        ledger = self.ledger(answers, wrapped=True)
        for r in ledger:
            r['response_text'] = text
        record = self.design_record(answers)
        self.write(folder, 'ledger.jsonl', ledger)
        self.write(folder, 'comparisons.jsonl', [record])
        self.candidates(folder, record)
        report = audit.audit_run(folder, 'design_audit')
        self.assertEqual(report['rederived_identically'], 1)
        self.assertEqual(report['unparseable_responses'], 0)

    def test_missing_candidates_are_unresolved(self):
        folder = self.root / 'design-audit-001'
        answers = [self.design_answer(), self.design_answer()]
        self.write(folder, 'ledger.jsonl', self.ledger(answers, wrapped=True))
        self.write(folder, 'comparisons.jsonl', [self.design_record(answers)])
        report = audit.audit_run(folder, 'design_audit')
        self.assertEqual((report['unresolved_records'], report['rederived_identically']), (1, 0))

    def test_design_audit_missing_answer_is_rederived_from_candidate_id(self):
        folder = self.root / 'design-audit-001'
        answers = [self.design_answer(), self.design_answer(cid=10)]
        record = {**self.design_record(answers), 'status': 'MODEL_ANSWER_MISSING'}
        for field in ('reasons', 'orientation', 'specimen', 'tumor', 'answers'):
            record.pop(field)
        self.write(folder, 'ledger.jsonl', self.ledger(answers, wrapped=True))
        self.write(folder, 'comparisons.jsonl', [record])
        self.candidates(folder, record)
        report = audit.audit_run(folder, 'design_audit')
        self.assertEqual(report['rederived_identically'], 1)

    def test_difference_examples_are_capped_at_twenty(self):
        folder = self.root / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        record = {**self.label_record(answers), 'decision': 'synthetic-tampered'}
        self.write(folder, 'ledger.jsonl', self.ledger(answers))
        self.write(folder, 'decisions-all.jsonl', [record] * 25)
        report = audit.audit_run(folder, 'adjudication')
        self.assertEqual(report['differing'], 25)
        self.assertEqual(len(report['differing_examples']), 20)

    def test_classification_reconciles_scores_and_targets_from_raw(self):
        folder = self.root / 'classification-001'
        a = {'id': 27, 'contrast_type': 'INTERVENTION', 'intervention_class': 'GENE_LOSS',
             'canonical_target': 'synthetic_gene', 'manipulation': 'DECREASE', 'relevance_score': 4,
             'relevance_codes': ['IN_VIVO']}
        b = {**a, 'relevance_score': 8}
        record = {**audit.classification.reconcile([a, b]), 'receipts': self.receipts()}
        self.write(folder, 'ledger.jsonl', self.ledger([a, b], wrapped=True))
        self.write(folder, 'comparisons.jsonl', [record, {**record, 'relevance_score': 9}])
        report = audit.audit_run(folder, 'classification')
        self.assertEqual((report['rederived_identically'], report['differing']), (1, 1))
        self.assertEqual(report['differing_examples'][0]['fields']['relevance_score']['derived'], 6)

    def test_classification_recovers_truncated_id_only_from_hash_matched_atlas(self):
        folder = self.root / 'classification-001'
        a = {'id': 1, 'contrast_type': 'OTHER', 'relevance_score': 1}
        record = {'modality': 'human DNAm', 'gse_id': 'synthetic-study', 'comparison_name': 'synthetic case',
                  'condition_name': 'synthetic condition', 'status': 'MODEL_ANSWER_TRUNCATED',
                  'receipts': self.receipts()}
        self.write(folder, 'ledger.jsonl', self.ledger([a, {**a, 'id': 2}], wrapped=True))
        self.write(folder, 'comparisons.jsonl', [record])
        atlas = self.root / 'atlas_v3_combined.csv'
        atlas.write_text(','.join((*audit.KEY, 'v3_final_class')) + '\n'
                         'human DNAm,synthetic-excluded,unused,unused,INVALID\n'
                         'human DNAm,synthetic-first,unused,unused,VERIFIED\n'
                         'human DNAm,synthetic-first,unused,unused,VERIFIED\n'
                         'human DNAm,synthetic-study,synthetic case,synthetic condition,PROVISIONAL\n')
        (folder / 'summary.json').write_text(json.dumps({'atlas_sha256': audit.source(atlas)['sha256']}))
        report = audit.audit_run(folder, 'classification')
        self.assertEqual((report['rederived_identically'], report['unresolved_records']), (1, 0))
        (folder / 'summary.json').write_text(json.dumps({'atlas_sha256': 'synthetic-wrong-hash'}))
        report = audit.audit_run(folder, 'classification')
        self.assertEqual((report['rederived_identically'], report['unresolved_records']), (0, 1))

    def test_reconstruction_merges_batches_and_checks_exported_arms(self):
        folder = self.root / 'design-dnam-001'
        rule = {'id': 0, 'case': {'group': ['case']}, 'control': {'group': ['control']}}
        other = {'id': 1, 'feasible': False}
        meta = pd.DataFrame({'group': ['case', 'control']}, index=['synthetic-GSM1', 'synthetic-GSM2'])
        receipts, ledger = [], []
        for model in ('m1', 'm2'):
            for cid, answer in enumerate((rule, other)):
                digest = f'{model}-{cid}'
                receipts.append({'model': model, 'status': 'OK', 'request_sha256': digest, 'batch': [cid]})
                ledger.append({'model_requested': model, 'status': 'OK', 'request_sha256': digest,
                               'response_text': json.dumps({'comparisons': [answer]})})
        result = audit.reconstruction.reconcile(meta, {'m1': {0: rule}, 'm2': {0: rule}}, 0)
        records = [{'gse_id': 'synthetic-study', 'row_ordinals': [1], 'receipts': receipts, **result},
                   {'gse_id': 'synthetic-study', 'row_ordinals': [2], 'receipts': receipts,
                    'status': 'INFEASIBLE_OR_MISSING', 'model': 'm1'}]
        exported = {'row_ordinal': 1, 'status': 'UNIQUE', 'source': 'design_reconstruction', 'kind': 'two_sample',
                    'case_gsms': ['synthetic-TAMPERED'], 'control_gsms': result['control_gsms'],
                    'n_case': 1, 'n_control': 1, 'design_rules': result['rules']}
        self.write(folder, 'ledger.jsonl', ledger)
        self.write(folder, 'comparisons.jsonl', records)
        self.write(folder, 'rows.jsonl', [exported])
        with patch.object(audit.reconstruction, 'load_meta', return_value=meta):
            report = audit.audit_run(folder, 'reconstruction')
        self.assertEqual((report['records'], report['rederived_identically'], report['differing']), (3, 2, 1))
        self.assertEqual(report['differing_examples'][0]['file'], 'rows.jsonl')

    def test_output_is_exclusive_and_inputs_unchanged(self):
        verifier, out = self.root / 'verifier', self.root / 'out'
        folder = verifier / 'adjudication-dnam-001'
        answers = [self.label_answer(), self.label_answer()]
        self.write(folder, 'ledger.jsonl', self.ledger(answers))
        self.write(folder, 'decisions-all.jsonl', [self.label_record(answers)])
        before = {p: p.read_bytes() for p in folder.iterdir()}
        with contextlib.redirect_stdout(io.StringIO()):
            report = audit.run(verifier, out)
        self.assertEqual(report['runs'][0]['rederived_identically'], 1)
        self.assertTrue((out / 'audit.json').exists())
        self.assertTrue((out / 'README.md').exists())
        self.assertEqual(before, {p: p.read_bytes() for p in folder.iterdir()})
        with self.assertRaises(FileExistsError):
            audit.run(verifier, out)
        with self.assertRaises(ValueError):  # never inside an audited run
            audit.run(verifier, folder / 'new-output')
        with contextlib.redirect_stdout(io.StringIO()):  # a final run directory may hold the audit
            audit.run(verifier, verifier / 'final-test' / 'decision-audit')
        self.assertEqual(before, {p: p.read_bytes() for p in folder.iterdir()})


if __name__ == '__main__':
    unittest.main()
