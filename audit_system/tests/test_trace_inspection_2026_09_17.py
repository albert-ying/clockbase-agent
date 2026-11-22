"""Software trace fixtures; no empirical error-rate observations."""
import hashlib
import json
from pathlib import Path
import pickle
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor

from src.detailed_trace_inspection import TraceArchiveLoader, extract_trace, verify_trace_request


class TraceInspectionTests(unittest.TestCase):
    def test_primary_stage_rejects_a_trace_from_another_study_before_output_ingestion(self):
        from unittest.mock import patch
        from langchain_core.messages import HumanMessage, AIMessage
        from src.pipeline_2026_09_10 import primary_analysis
        class WrongStudyAgent:
            def __init__(self, **kwargs): pass
            def configure(self): pass
            def go(self, prompt):
                return [HumanMessage(content='Study identifier: GSE2\nWrong source fixture').pretty_repr(),
                        AIMessage(content='<solution>done</solution>').pretty_repr()], '<solution>done</solution>'
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); source = root / 'predictions.csv'
            source.write_text('sample_id,fixture_clock\nGSM1,1\n')
            with patch('src.pipeline_2026_09_10.CodingAgent', WrongStudyAgent):
                with self.assertRaisesRegex(ValueError, 'source-bound study request'):
                    primary_analysis({'study_id': 'GSE1', 'target': 'software fixture'},
                                     {'path': str(source), 'clocks': ['fixture_clock']}, None, root)
            trace = json.loads((root / 'trace.json').read_text())
            self.assertEqual(trace['request_check']['status'], 'FAIL')
            self.assertEqual(trace['study_id'], 'GSE1')
            self.assertEqual(trace['input_sha256']['clock_predictions.csv'], hashlib.sha256(source.read_bytes()).hexdigest())

    def test_explicit_request_identity_catches_mixed_studies_without_flagging_related_mentions(self):
        from langchain_core.messages import HumanMessage, AIMessage
        prompt = 'Study identifier: GSE1\nSource metadata refers to related study GSE2.'
        messages = [HumanMessage(content=prompt).pretty_repr(), AIMessage(content='I will inspect GSE2 as a reference.').pretty_repr()]
        self.assertEqual(verify_trace_request(messages, 'GSE1', prompt)['status'], 'PASS')
        self.assertEqual(verify_trace_request(messages, 'GSE2', prompt)['status'], 'FAIL')
        self.assertEqual(verify_trace_request(messages, 'GSE1', prompt + '\nChanged target')['status'], 'FAIL')
        self.assertEqual(verify_trace_request(messages * 2, 'GSE1')['status'], 'FAIL')
        legacy = [HumanMessage(content='Analyze the aging effects of conditions from study GSE2.').pretty_repr()]
        self.assertEqual(verify_trace_request(legacy, 'GSE1')['status'], 'FAIL')
        unknown = [HumanMessage(content='Source metadata only; no explicit request ID').pretty_repr()]
        self.assertEqual(verify_trace_request(unknown, 'GSE1')['status'], 'CANNOT VERIFY')

    def source(self, root, entry=None):
        entry = entry or {'study_id': 'GSE1', 'metadata': {'scope': 'SOFTWARE FIXTURE'},
                          'result': (['GSE2: wrong study in software fixture', 'GSE1 result'], 'Entire terminal reply'),
                          'status': 'success', 'execution_time': 2.0}
        path = root / 'source.pkl'; path.write_bytes(pickle.dumps({'GSE1': entry}))
        return path, entry

    def test_full_selected_entry_keeps_mixed_studies_final_reply_and_tuple_provenance(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); source, original = self.source(root)
            result = extract_trace(source, 'GSE1', root / 'output')
            selected = json.loads(Path(result['selected_trace']).read_text())
            self.assertEqual(selected['entry'], json.loads(json.dumps(original)))
            self.assertEqual(selected['serialization']['tuple_paths'], ['/result'])
            self.assertEqual(result['mentioned_studies'], ['GSE1', 'GSE2'])
            self.assertEqual(result['message_count'], 2)
            self.assertFalse(result['trace_code_executed'])
            before = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (root / 'output').rglob('*') if p.is_file()}
            self.assertEqual(extract_trace(source, 'GSE1', root / 'output'), result)
            self.assertEqual(before, {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (root / 'output').rglob('*') if p.is_file()})

    def test_pickle_reducer_is_refused_before_execution(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); marker = root / 'must-not-exist'
            class ExecutableFixture:
                def __reduce__(self):
                    return eval, ("__import__('pathlib').Path(" + repr(str(marker)) + ").touch()",)
            path = root / 'source.pkl'; path.write_bytes(pickle.dumps({'GSE1': ExecutableFixture()}))
            with self.assertRaisesRegex(pickle.UnpicklingError, 'Executable'):
                extract_trace(path, 'GSE1', root / 'output')
            self.assertFalse(marker.exists())
            self.assertFalse(list((root / 'output').rglob('selected-trace.json')))

    def test_changed_source_missing_study_and_credential_like_text_are_not_saved(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); source, entry = self.source(root)
            extract_trace(source, 'GSE1', root / 'output')
            with self.assertRaisesRegex(ValueError, 'absent'):
                extract_trace(source, 'GSE3', root / 'missing')
            entry['result'][0].append('sk-' + 'SOFTWAREFIXTURETOKENONLY0123456789')
            source.write_bytes(pickle.dumps({'GSE1': entry}))
            with self.assertRaisesRegex(ValueError, 'Configuration'):
                extract_trace(source, 'GSE1', root / 'output')
            with self.assertRaisesRegex(ValueError, 'Credential-like'):
                extract_trace(source, 'GSE1', root / 'secret')
            self.assertFalse(list((root / 'secret').rglob('selected-trace.json')))

    def test_bulk_loader_is_hash_bound_cached_sanitized_and_post_design_only(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            prompt = 'Study identifier: GSE1\nInspect the supplied study.'
            messages = [
                '================================ Human Message =================================\n\n' + prompt,
                'Result used /Entropy/kying/private/input.csv; contact_email: person@example.org; '
                'Authorization: Bearer software-fixture-secret-token',
            ]
            terminal = 'Saved to /home/kying/private/output.csv'
            archive = {'GSE1': {'study_id': 'GSE1', 'status': 'success',
                                'result': (messages, terminal)}}
            path = root / 'traces.pkl'; path.write_bytes(pickle.dumps(archive))
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            loader = TraceArchiveLoader({'mouse_rna': {'path': str(path), 'sha256': digest}})
            design_hash = 'a' * 64

            def read():
                return loader.study_evidence('mouse_rna', 'GSE1',
                                             frozen_source_design_sha256=design_hash,
                                             expected_prompt=prompt)
            with ThreadPoolExecutor(max_workers=8) as pool:
                results = list(pool.map(lambda _: read(), range(16)))
            self.assertTrue(all(result == results[0] for result in results))
            result = results[0]
            self.assertEqual(loader.load_counts(), {'mouse_rna': 1})
            self.assertTrue(result['eligible_for_post_design_context'])
            self.assertEqual(result['evidence_count'], 3)
            self.assertEqual(result['source_design_sha256'], design_hash)
            self.assertEqual(result['source_sha256'], digest)
            self.assertFalse(result['raw_archive_object_exposed'])
            self.assertFalse(result['archived_code_executed'])
            serialized = json.dumps(result)
            for private in ('/Entropy/kying', '/home/kying', 'person@example.org',
                            'software-fixture-secret-token'):
                self.assertNotIn(private, serialized)
            self.assertEqual(result['evidence'][1]['original_sha256'],
                             hashlib.sha256(messages[1].encode()).hexdigest())
            self.assertGreaterEqual(len(result['evidence'][1]['redactions']), 1)
            self.assertTrue({'credential', 'private_path', 'email', 'contact_field'} <=
                            set(result['evidence'][1]['redactions'][0]['kinds']))
            self.assertEqual(path.read_bytes(), pickle.dumps(archive))

    def test_bulk_loader_withholds_request_mismatch_and_rejects_unpinned_or_executable_archives(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            wrong = self.source(root)[0]
            digest = hashlib.sha256(wrong.read_bytes()).hexdigest()
            loader = TraceArchiveLoader({'mouse_rna': {'path': wrong, 'sha256': digest}})
            result = loader.study_evidence('mouse_rna', 'GSE1',
                                           frozen_source_design_sha256='b' * 64)
            self.assertFalse(result['eligible_for_post_design_context'])
            self.assertEqual(result['evidence'], [])
            self.assertEqual(result['request_binding']['status'], 'FAIL')
            with self.assertRaisesRegex(ValueError, 'source-design'):
                loader.study_evidence('mouse_rna', 'GSE1', frozen_source_design_sha256='not-a-hash')

            changed = root / 'changed.pkl'; changed.write_bytes(wrong.read_bytes())
            changed_loader = TraceArchiveLoader({'mouse_rna': {'path': changed, 'sha256': '0' * 64}})
            with self.assertRaisesRegex(ValueError, 'pinned declaration'):
                changed_loader.study_evidence('mouse_rna', 'GSE1', frozen_source_design_sha256='c' * 64)

            marker = root / 'must-not-exist'
            class ExecutableFixture:
                def __reduce__(self):
                    return eval, ("__import__('pathlib').Path(" + repr(str(marker)) + ").touch()",)
            executable = root / 'executable.pkl'
            executable.write_bytes(pickle.dumps({'GSE1': ExecutableFixture()}))
            executable_loader = TraceArchiveLoader({'mouse_rna': {'path': executable,
                'sha256': hashlib.sha256(executable.read_bytes()).hexdigest()}})
            with self.assertRaisesRegex(ValueError, 'Executable pickle object refused'):
                executable_loader.study_evidence('mouse_rna', 'GSE1', frozen_source_design_sha256='d' * 64)
            self.assertFalse(marker.exists())

    def test_bulk_loader_supports_three_distinct_archives_and_final_hash_recheck(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); declarations = {}
            for index, modality in enumerate(('mouse_rna', 'human_rna', 'human_dnam'), 1):
                study = 'GSE' + str(index)
                prompt = 'Study identifier: ' + study + '\nInspect the supplied study.'
                message = ('================================ Human Message =================================\n\n' + prompt)
                path = root / (modality + '.pkl')
                path.write_bytes(pickle.dumps({study: {'study_id': study, 'status': 'success',
                                                       'result': ([message], 'done')}}))
                declarations[modality] = {'path': path,
                                          'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
            loader = TraceArchiveLoader(declarations)
            for index, modality in enumerate(declarations, 1):
                result = loader.study_evidence(modality, 'GSE' + str(index),
                                               frozen_source_design_sha256='e' * 64)
                self.assertTrue(result['eligible_for_post_design_context'])
            self.assertEqual(loader.load_counts(), {key: 1 for key in declarations})
            ledger = loader.verify_sources_unchanged()
            self.assertEqual({row['modality'] for row in ledger}, set(declarations))
            self.assertTrue(all(row['loaded_once'] for row in ledger))
            Path(declarations['mouse_rna']['path']).write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'changed after'):
                loader.verify_sources_unchanged()

    def test_bulk_loader_retains_unbound_context_without_upgrading_identity(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            message = ('================================ Human Message =================================\n\n'
                       'Inspect the supplied metadata without an explicit identifier.')
            path = root / 'human.pkl'
            path.write_bytes(pickle.dumps({'GSE9': {'study_id': 'GSE9', 'status': 'success',
                                                    'result': ([message], 'historical result')}}))
            loader = TraceArchiveLoader({'human_rna': {'path': path,
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}})
            result = loader.study_evidence('human_rna', 'GSE9',
                                           frozen_source_design_sha256='f' * 64)
            self.assertEqual(result['request_binding']['status'], 'CANNOT VERIFY')
            self.assertTrue(result['eligible_for_post_design_context'])
            self.assertEqual(result['identity_scope'], 'UNBOUND_CONTEXT')
            self.assertIn('cannot alone support', result['identity_limit'])
            self.assertEqual(len(result['evidence']), 2)


if __name__ == '__main__':
    unittest.main()
