"""Software fixtures for trace binding; these are not empirical study results."""
import pickle
from pathlib import Path
import tempfile
import unittest
import threading
from concurrent.futures import ThreadPoolExecutor

from src.combine_pickle_traces import combine_pickle_files, trace_binding


def record(study, requested=None, status='success'):
    return {'study_id': study, 'status': status, 'result': ([
        '================================ Human Message =================================\n\n'
        'Study identifier: ' + (requested or study) + '\nInspect the supplied study.'
    ], 'Finished')}


class ExecutableObject:
    def __reduce__(self):
        return (str, ('not allowed',))


class TraceCombinationTests(unittest.TestCase):
    def test_wrong_request_and_key_are_failures(self):
        self.assertEqual(trace_binding('GSE1', record('GSE1', 'GSE2'))['status'], 'FAIL')
        self.assertEqual(trace_binding('GSE1', record('GSE2'))['status'], 'FAIL')

    def test_unicode_digits_cannot_pass_as_geo_accessions(self):
        from src.detailed_trace_inspection import verify_trace_request
        for study in ('GSE１２３', 'GSE١٢٣'):
            item = record(study)
            self.assertEqual(trace_binding(study, item)['status'], 'FAIL')
            self.assertEqual(verify_trace_request(item['result'][0], study)['status'], 'FAIL')

    def test_unknown_request_is_retained_as_unresolved(self):
        item = record('GSE1')
        item['result'][0][0] = '================================ Human Message =================================\n\nInspect this study.'
        self.assertEqual(trace_binding('GSE1', item)['status'], 'CANNOT VERIFY')

    def test_merge_withholds_mismatches_and_records_every_input(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            archive = {'GSE1': record('GSE1'), 'GSE2': record('GSE2', 'GSE9')}
            payload = pickle.dumps(archive)
            (root/'one.pkl').write_bytes(payload)
            ledger = []
            merged, files, total = combine_pickle_files(root, ledger)
            self.assertEqual((files, total), (1, 2))
            self.assertEqual(list(merged), ['GSE1'])
            self.assertEqual(ledger[1]['disposition'], 'WITHHELD')
            self.assertEqual((root/'one.pkl').read_bytes(), payload)

    def test_executable_pickle_and_trailing_bytes_never_merge(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            (root/'executable.pkl').write_bytes(pickle.dumps(ExecutableObject()))
            (root/'trailing.pkl').write_bytes(pickle.dumps({'GSE1': record('GSE1')}) + b'x')
            ledger = []
            merged, files, total = combine_pickle_files(root, ledger)
            self.assertEqual((merged, files, total), ({}, 0, 0))
            self.assertTrue(all(item['disposition'] == 'WITHHELD' for item in ledger))

    def test_sorted_sources_and_status_priority_are_deterministic(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            for filename, marker, state in [('z.pkl', 'last', 'success'), ('a.pkl', 'first', 'success'), ('b.pkl', 'error', 'error')]:
                item = record('GSE1', status=state); item['marker'] = marker
                (root/filename).write_bytes(pickle.dumps({'GSE1': item}))
            merged, _, _ = combine_pickle_files(root)
            self.assertEqual(merged['GSE1']['marker'], 'first')

    def test_file_limit_is_checked_before_load(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name); (root/'one.pkl').write_bytes(pickle.dumps({'GSE1': record('GSE1')}))
            ledger = []
            self.assertEqual(combine_pickle_files(root, ledger, max_bytes=1), ({}, 0, 0))
            self.assertEqual(ledger[0]['status'], 'FAIL')


class CodingAgentTraceIsolationTests(unittest.TestCase):
    def make_agent(self, app):
        from src.coding_agent import CodingAgent
        agent = CodingAgent(llm=object(), printout=False)
        agent.app = app
        return agent

    def test_checkpoint_and_returned_logs_are_isolated_per_call(self):
        class App:
            identifiers = []

            def stream(self, inputs, **kwargs):
                self.identifiers.append(kwargs['config']['configurable']['thread_id'])
                yield inputs

        app = App()
        agent = self.make_agent(app)
        first, _ = agent.go('GSE1')
        second, _ = agent.go('GSE2')
        self.assertEqual(len(set(app.identifiers)), 2)
        self.assertIsNot(first, second)
        second.append('caller mutation')
        self.assertNotIn('caller mutation', agent.log)
        self.assertIn('GSE1', first[0])
        self.assertIn('GSE2', second[0])

    def test_shared_instance_serializes_concurrent_trace_collection(self):
        entered, release, second_started = threading.Event(), threading.Event(), threading.Event()

        class App:
            active = 0
            maximum_active = 0

            def stream(self, inputs, **kwargs):
                self.active += 1
                self.maximum_active = max(self.maximum_active, self.active)
                try:
                    if inputs['messages'][0].content == 'GSE1':
                        entered.set()
                        if not release.wait(3):
                            raise TimeoutError('fixture release missing')
                    yield inputs
                finally:
                    self.active -= 1

        app = App()
        agent = self.make_agent(app)

        def second_call():
            second_started.set()
            return agent.go('GSE2')

        with ThreadPoolExecutor(max_workers=2) as pool:
            first_future = pool.submit(agent.go, 'GSE1')
            try:
                self.assertTrue(entered.wait(3))
                second_future = pool.submit(second_call)
                self.assertTrue(second_started.wait(3))
            finally:
                release.set()
            first, second = first_future.result(), second_future.result()
        self.assertEqual(app.maximum_active, 1)
        self.assertEqual(len(first[0]), 1)
        self.assertEqual(len(second[0]), 1)
        self.assertNotIn('GSE2', first[0][0])
        self.assertNotIn('GSE1', second[0][0])

    def test_failed_stream_releases_instance_lock(self):
        class App:
            calls = 0

            def stream(self, inputs, **kwargs):
                self.calls += 1
                if self.calls == 1:
                    raise ValueError('fixture stream failure')
                yield inputs

        agent = self.make_agent(App())
        with self.assertRaisesRegex(ValueError, 'fixture stream failure'):
            agent.go('GSE1')
        self.assertEqual(agent.go('GSE2')[1], 'GSE2')


if __name__ == '__main__':
    unittest.main()
