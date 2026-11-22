import json
import os
from pathlib import Path
import subprocess
from types import SimpleNamespace
import tempfile
import traceback
import unittest
from unittest.mock import patch

from langchain_core.messages import HumanMessage

from src.pipeline_runtime_2026_09_10 import BoundedOpenAI, PipelineIntegrityError
from audit_system.agent_workflow_2026_09_10 import sha


class TestPipelineRuntimeHTTPDiagnostics(unittest.TestCase):
    @staticmethod
    def successful_response(content='{}'):
        body = {'model': 'gpt-5.6-sol',
                'usage': {'prompt_tokens': 10, 'completion_tokens': 2},
                'choices': [{'finish_reason': 'stop', 'message': {'content': content}}]}
        return SimpleNamespace(returncode=0, stdout=(json.dumps(body) + '\n200').encode())

    def test_json_response_is_opt_in_and_frozen_in_request(self):
        with tempfile.TemporaryDirectory() as folder, \
                patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}, clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                      return_value=self.successful_response()) as transport:
            default = BoundedOpenAI(Path(folder) / 'default', 'fixture', max_cost=1.0)
            default.invoke([HumanMessage(content='software fixture')])
            default_request = json.loads((Path(folder) / 'default/fixture-000/request.json').read_text())
            self.assertNotIn('response_format', default_request)

            enabled = BoundedOpenAI(Path(folder) / 'enabled', 'fixture', max_cost=1.0,
                                    json_response=True)
            enabled.invoke([HumanMessage(content='software fixture')])
            enabled_request = json.loads((Path(folder) / 'enabled/fixture-000/request.json').read_text())
            self.assertEqual(enabled_request['response_format'], {'type': 'json_object'})
            self.assertEqual(transport.call_count, 2)
        for invalid in (None, 0, 1, 'true', []):
            with self.subTest(invalid=invalid), self.assertRaisesRegex(ValueError, 'boolean'):
                BoundedOpenAI('/unused/software-fixture', 'fixture', json_response=invalid)

    def test_http_error_retains_allowlisted_diagnostic_and_never_retries(self):
        key = 'software-fixture-secret-key'
        private_prompt = 'private source phrase 67d1a9a9'
        body = {'error': {
            'type': 'invalid_request_error',
            'code': 'context_length_exceeded',
            'param': 'messages',
            'message': ('Request contained "' + private_prompt + '" at '
                        '/Entropy/kying/private/input.csv with Authorization: Bearer ' + key),
            'internal': {'headers': {'authorization': key}},
        }, 'request_id': 'private-provider-id'}
        response = SimpleNamespace(returncode=0,
                                   stdout=(json.dumps(body) + '\n400').encode())
        with tempfile.TemporaryDirectory() as folder, \
                patch.dict(os.environ, {'OPENAI_API_KEY': key}, clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                      return_value=response) as transport:
            client = BoundedOpenAI(folder, 'fixture', max_calls=1, max_cost=1.0)
            with self.assertRaisesRegex(RuntimeError, r'OpenAI HTTP 400; response suppressed'):
                client.invoke([HumanMessage(content=private_prompt)])
            diagnostic_path = Path(folder) / 'fixture-000/provider-error.json'
            diagnostic = json.loads(diagnostic_path.read_text())
            self.assertEqual(diagnostic, {
                'schema': 'clockbase-openai-error-v1',
                'provider': 'openai',
                'http_status': 400,
                'billing_status': 'uncertain',
                'automatic_retry': False,
                'raw_response_retained': False,
                'error': {
                    'type': 'invalid_request_error',
                    'code': 'context_length_exceeded',
                    'param': 'messages',
                    'message': ('Request contained [REDACTED_REQUEST_TEXT] at '
                                '[REDACTED_LOCAL_PATH] with '
                                'Authorization: Bearer [REDACTED_CREDENTIAL]'),
                },
            })
            serialized = diagnostic_path.read_text()
            for prohibited in (key, private_prompt, '/Entropy/kying', 'private-provider-id',
                               'headers', 'authorization":'):
                self.assertNotIn(prohibited, serialized)
            self.assertFalse((diagnostic_path.parent / 'response.json').exists())
            self.assertFalse((diagnostic_path.parent / 'receipt.json').exists())
            self.assertTrue((diagnostic_path.parent / 'reservation.json').exists())
            self.assertEqual(transport.call_count, 1)

            resumed = BoundedOpenAI(folder, 'fixture', max_calls=1, max_cost=1.0)
            with self.assertRaisesRegex(RuntimeError, 'uncertain billing'):
                resumed.invoke([HumanMessage(content=private_prompt)])
            self.assertEqual(transport.call_count, 1)

    def test_unstructured_http_error_body_is_not_retained(self):
        raw_body = b'upstream included /home/kying/private and software-fixture-secret-key'
        response = SimpleNamespace(returncode=0, stdout=raw_body + b'\n429')
        with tempfile.TemporaryDirectory() as folder, \
                patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-secret-key'}, clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run', return_value=response):
            client = BoundedOpenAI(folder, 'fixture', max_calls=1, max_cost=1.0)
            with self.assertRaisesRegex(RuntimeError, r'OpenAI HTTP 429; response suppressed'):
                client.invoke([HumanMessage(content='software fixture')])
            diagnostic = json.loads((Path(folder) / 'fixture-000/provider-error.json').read_text())
            self.assertEqual(diagnostic['http_status'], 429)
            self.assertEqual(diagnostic['error'], {
                'type': None, 'code': None, 'param': None,
                'message': 'Provider returned an unstructured error body; content was not retained',
            })
            self.assertNotIn(raw_body.decode(), json.dumps(diagnostic))

    def test_timeout_retains_safe_diagnostic_and_uncertain_slot_without_retry(self):
        key = 'software-fixture-timeout-secret'
        private_output = 'software-only-sensitive-output'
        error = subprocess.TimeoutExpired(['software-private-command', key], 190,
            output=(private_output + key).encode(), stderr=b'software-private-stderr')
        with tempfile.TemporaryDirectory() as folder, \
                patch.dict(os.environ, {'OPENAI_API_KEY': key}, clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run', side_effect=error) as transport:
            client = BoundedOpenAI(folder, 'fixture', max_calls=1, max_cost=1.0)
            with self.assertRaisesRegex(RuntimeError, 'OPENAI_TRANSPORT_TIMEOUT') as caught:
                client.invoke([HumanMessage(content='software fixture')])
            slot = Path(folder) / 'fixture-000'
            diagnostic = json.loads((slot / 'transport-error.json').read_text())
            self.assertEqual(set(diagnostic), {'schema', 'provider', 'category', 'error_type',
                'timeout_seconds', 'elapsed_seconds', 'request_sha256', 'billing_status',
                'automatic_retry', 'raw_output_retained'})
            self.assertEqual(diagnostic['schema'], 'clockbase-openai-transport-error-v1')
            self.assertEqual(diagnostic['category'], 'OPENAI_TRANSPORT_TIMEOUT')
            self.assertEqual(diagnostic['error_type'], 'TimeoutExpired')
            self.assertEqual(diagnostic['timeout_seconds'], 190)
            self.assertIsInstance(diagnostic['elapsed_seconds'], float)
            self.assertGreaterEqual(diagnostic['elapsed_seconds'], 0)
            self.assertEqual(diagnostic['request_sha256'], sha(slot / 'request.json'))
            self.assertEqual(diagnostic['billing_status'], 'uncertain')
            self.assertFalse(diagnostic['automatic_retry'])
            self.assertFalse(diagnostic['raw_output_retained'])
            serialized = json.dumps(diagnostic)
            self.assertLess(len(serialized.encode()), 1024)
            rendered_error = ''.join(traceback.format_exception(caught.exception))
            for prohibited in (key, private_output, 'software-private-command', 'software-private-stderr'):
                self.assertNotIn(prohibited, serialized)
                self.assertNotIn(prohibited, rendered_error)
            self.assertTrue(caught.exception.__suppress_context__)
            self.assertTrue((slot / 'reservation.json').is_file())
            self.assertFalse((slot / 'response.json').exists())
            self.assertFalse((slot / 'receipt.json').exists())
            before = {p.name: sha(p) for p in slot.iterdir()}
            resumed = BoundedOpenAI(folder, 'fixture', max_calls=1, max_cost=1.0)
            with self.assertRaisesRegex(RuntimeError, 'uncertain billing'):
                resumed.invoke([HumanMessage(content='software fixture')])
            self.assertEqual(transport.call_count, 1)
            self.assertEqual(before, {p.name: sha(p) for p in slot.iterdir()})

    def test_timeout_handler_does_not_translate_integrity_errors(self):
        with tempfile.TemporaryDirectory() as folder, \
                patch.dict(os.environ, {'OPENAI_API_KEY': 'software-fixture-key'}, clear=True), \
                patch('src.pipeline_runtime_2026_09_10.subprocess.run',
                      side_effect=PipelineIntegrityError('software integrity failure')) as transport:
            with self.assertRaisesRegex(PipelineIntegrityError, 'software integrity failure'):
                BoundedOpenAI(folder, 'fixture').invoke([HumanMessage(content='software fixture')])
            self.assertFalse((Path(folder) / 'fixture-000/transport-error.json').exists())
            self.assertEqual(transport.call_count, 1)


if __name__ == '__main__':
    unittest.main()
