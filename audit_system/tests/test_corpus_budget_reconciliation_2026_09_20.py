"""Software fixtures for immutable aggregate API-budget reconciliation."""
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import types
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from audit_system.agent_workflow_2026_09_10 import write_new
from audit_system.corpus_audit_2026_09_20 import Budget, StudyClient
from src.pipeline_runtime_2026_09_10 import PipelineIntegrityError


class Message:
    def __init__(self, content, kind):
        self.content, self.type = content, kind


class BudgetReconciliationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.output = self.root / 'study'
        self.output.mkdir()
        self.config = {'api': {'model': 'gpt-5.6-luna', 'input_per_million': 0.2,
            'output_per_million': 1.2, 'output_tokens': 10, 'prompt_bytes': 10000,
            'max_calls_per_study': 4, 'max_cost_per_study': 1.0,
            'json_response': True}}
        messages = types.ModuleType('langchain_core.messages')
        messages.SystemMessage = lambda content: Message(content, 'system')
        messages.HumanMessage = lambda content: Message(content, 'human')
        messages.AIMessage = lambda content: Message(content, 'ai')
        package = types.ModuleType('langchain_core')
        package.messages = messages
        self.modules = {'langchain_core': package, 'langchain_core.messages': messages}

    def tearDown(self):
        self.tmp.cleanup()

    def call(self, budget, *, success=True):
        response = {'model': 'gpt-5.6-luna', 'usage': {'prompt_tokens': 10,
            'completion_tokens': 3, 'total_tokens': 13},
            'choices': [{'finish_reason': 'stop', 'message': {'content': '{}'}}]}
        transport = (SimpleNamespace(returncode=0,
                     stdout=json.dumps(response).encode() + b'\n200') if success else
                     SimpleNamespace(returncode=1, stdout=b''))
        with patch.dict(sys.modules, self.modules), patch.dict(
                os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}), patch(
                'src.pipeline_runtime_2026_09_10.subprocess.run', return_value=transport):
            return StudyClient(self.output, self.config, budget).ask(
                'source_design_00', 'software fixture', {'study_id': 'GSE1'})

    def test_rebuilds_completed_usage_from_verified_call_files(self):
        database = self.root / 'budget.sqlite3'
        budget = Budget(database, 1.0, 4)
        self.assertEqual(self.call(budget), '{}')
        charged = (10 * 0.2 + 3 * 1.2) / 1e6
        self.assertAlmostEqual(budget.report()['estimated_or_reserved_usd'], charged)
        database.unlink()
        rebuilt = Budget(database, 1.0, 4)
        result = rebuilt.reconcile(self.root)
        self.assertEqual(result['bound_calls'], 1)
        self.assertEqual(result['usage_estimates'], 1)
        self.assertEqual(result['uncertain_reservations'], 0)
        self.assertAlmostEqual(result['estimated_or_reserved_usd'], charged)
        self.assertAlmostEqual(rebuilt.report()['estimated_or_reserved_usd'], charged)

    def test_uncertain_transport_retains_full_reservation(self):
        database = self.root / 'budget.sqlite3'
        budget = Budget(database, 1.0, 4)
        with self.assertRaisesRegex(RuntimeError, 'OpenAI transport failed'):
            self.call(budget, success=False)
        reserved = budget.report()['estimated_or_reserved_usd']
        self.assertGreater(reserved, 0)
        database.unlink()
        rebuilt = Budget(database, 1.0, 4)
        result = rebuilt.reconcile(self.root)
        self.assertEqual(result['usage_estimates'], 0)
        self.assertEqual(result['uncertain_reservations'], 1)
        self.assertAlmostEqual(result['estimated_or_reserved_usd'], reserved)

    def test_tampered_request_receipt_and_binding_fail_closed(self):
        database = self.root / 'budget.sqlite3'
        budget = Budget(database, 1.0, 4)
        self.call(budget)
        slot = self.output / 'api' / 'source_design_00-000'
        original = (slot / 'request.json').read_text()
        original_receipt = (slot / 'receipt.json').read_text()
        (slot / 'request.json').write_text(original + ' ')
        with self.assertRaisesRegex(PipelineIntegrityError, 'request changed'):
            budget.reconcile(self.root)
        (slot / 'request.json').write_text(original)
        receipt = json.loads((slot / 'receipt.json').read_text())
        receipt['estimated_cost_usd'] = 0.0
        (slot / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
        with self.assertRaisesRegex(PipelineIntegrityError, 'receipt or response changed'):
            budget.reconcile(self.root)
        (slot / 'receipt.json').write_text(original_receipt)
        binding = self.output / 'api' / 'source_design_00-000.aggregate-binding.json'
        value = json.loads(binding.read_text())
        value['reserved_usd'] = math.nextafter(value['reserved_usd'], math.inf)
        binding.write_text(json.dumps(value, indent=2) + '\n')
        with self.assertRaises(PipelineIntegrityError):
            budget.reconcile(self.root)

    def test_legacy_slot_without_binding_refuses_paid_resume(self):
        slot = self.output / 'api' / 'source_design_00-000'
        slot.mkdir(parents=True)
        write_new(slot / 'reservation.json', {'reserved_usd': 0.2,
            'started_unix': 1.0, 'model': 'gpt-5.6-luna',
            'automatic_retries': 0, 'concurrency': 1})
        budget = Budget(self.root / 'budget.sqlite3', 1.0, 4)
        with self.assertRaisesRegex(PipelineIntegrityError, 'lack aggregate bindings'):
            budget.reconcile(self.root)

    def test_carry_forward_keeps_prior_calls_in_the_same_cap_and_preserves_old_ledger(self):
        from audit_system.agent_workflow_2026_09_10 import sha
        prior_root=self.output
        old_path=self.root/'old-budget.sqlite3';old=Budget(old_path,1.0,2)
        self.call(old);old_hash=sha(old_path)
        current_root=self.root/'current';self.output=current_root/'study';self.output.mkdir(parents=True)
        current=Budget(self.root/'new-budget.sqlite3',1.0,2);self.call(current)
        result=current.reconcile(current_root,carry_forward_roots=(prior_root,))
        self.assertEqual(result['bound_calls'],2)
        self.assertEqual(current.report()['call_reservations'],2)
        self.assertEqual(sha(old_path),old_hash)
        with self.assertRaisesRegex(RuntimeError,'AGGREGATE_BUDGET_EXHAUSTED'):
            current.reserve('f'*64,.01,'e'*64)

    def test_carry_forward_retains_uncertain_reservations(self):
        prior_root=self.output
        prior=Budget(self.root/'old-budget.sqlite3',1.0,4)
        with self.assertRaisesRegex(RuntimeError,'OpenAI transport failed'):
            self.call(prior,success=False)
        reserved=prior.report()['estimated_or_reserved_usd']
        current=Budget(self.root/'new-budget.sqlite3',1.0,4)
        result=current.reconcile(self.root/'current',carry_forward_roots=(prior_root,))
        self.assertEqual(result['uncertain_reservations'],1)
        self.assertEqual(current.report()['estimated_or_reserved_usd'],reserved)

    def test_carry_forward_rejects_overlapping_roots_and_changed_old_requests(self):
        prior_root=self.output
        old=Budget(self.root/'old-budget.sqlite3',1.0,4);self.call(old)
        current=Budget(self.root/'new-budget.sqlite3',1.0,4)
        for roots in ((prior_root,prior_root),(self.root,)):
            with self.assertRaisesRegex(PipelineIntegrityError,'overlapping'):
                current.reconcile(self.root/'current',carry_forward_roots=roots)
        request=prior_root/'api/source_design_00-000/request.json'
        request.write_text(request.read_text()+' ')
        with self.assertRaisesRegex(PipelineIntegrityError,'request changed'):
            current.reconcile(self.root/'current',carry_forward_roots=(prior_root,))


if __name__ == '__main__':
    unittest.main()
