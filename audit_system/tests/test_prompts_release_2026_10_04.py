"""Source equality and provenance checks for the released prompts and historical scripts."""
import ast
import csv
import hashlib
import pickle
import re
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from prompts import build_prompts
from src import parallel_hs_dnam
from audit_system import paths
from audit_system.tests import require_development_artifact


REPO = Path(__file__).resolve().parents[2]


def source_constant(source, symbol):
    """Evaluate only the prompt assignment or function, without importing its driver."""
    tree = ast.parse((REPO / source).read_text())
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == symbol
                                                for t in node.targets):
            return ast.literal_eval(node.value)
        if isinstance(node, ast.FunctionDef) and node.name == symbol:
            namespace = {'BASE_DIR': '{BASE_DIR}', 'RESULTS_DIR': '{RESULTS_DIR}'}
            isolated = ast.Module(body=[node], type_ignores=[])
            exec(compile(isolated, source, 'exec'), namespace)
            return namespace[symbol]('{study_id}', '{study_design}').replace('{BASE_DIR}/data', '{data_dir}')
    raise AssertionError(f'Source not found: {source}::{symbol}')


class PromptsReleaseTests(unittest.TestCase):
    def test_every_text_file_has_a_source(self):
        files = {path.name for path in (REPO / 'prompts').glob('*.txt')}
        self.assertEqual(files, set(build_prompts.SOURCES))
        self.assertTrue({'task_intention_human_dnam.txt', 'generate_paper_human_dnam.txt',
                         'generate_paper_human_rna.txt', 'generate_paper_mouse_rna.txt'} <= files)

    def test_every_body_equals_its_source(self):
        for name, (source, symbol) in build_prompts.SOURCES.items():
            with self.subTest(prompt=name):
                header, separator, body = (REPO / 'prompts' / name).read_text().partition('\n')
                self.assertEqual(separator, '\n')
                self.assertTrue(header.startswith(f'# Source: {source}::{symbol}'))
                self.assertEqual(body, source_constant(source, symbol))

    def test_dnam_template_equals_all_saved_trace_prompts(self):
        path = paths.TRACES / parallel_hs_dnam.TRACE_SOURCE  # the bundle's agent-traces/, else the repository root
        if not path.is_file():
            path = REPO / parallel_hs_dnam.TRACE_SOURCE
        require_development_artifact(path)
        self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), parallel_hs_dnam.TRACE_SHA256)
        with path.open('rb') as handle:
            traces = pickle.load(handle)
        self.assertEqual(len(traces), 2296)
        for key, entry in traces.items():
            with self.subTest(trace=key):
                wrapper, separator, body = entry['result'][0][0].partition('\n\n')
                self.assertIn('Human Message', wrapper)
                self.assertEqual(separator, '\n\n')
                study_id = re.search(r'from study (GSE\d+)', body).group(1)
                before, design_header, rest = body.partition(
                    '2. Study Design Information (Use this information instead of fetching from GEO):\n   ')
                design, analysis_header, after = rest.partition('\n\n3. Clock Selection and Analysis:')
                self.assertTrue(design_header and analysis_header)
                expected = parallel_hs_dnam.TASK_TEMPLATE.format(
                    study_id=study_id, study_design=design,
                    data_dir='/dfs/user/hanchenw/Agent_PG/clockbase_agent//data',
                    RESULTS_DIR='/dfs/user/hanchenw/Agent_PG/clockbase_agent//results/hs_dnam_sonnet37')
                self.assertEqual(body, expected)

    def test_dnam_driver_uses_the_configured_geo_root(self):
        with mock.patch.object(parallel_hs_dnam.paths, 'GEO', Path('/configured/geo')), \
             mock.patch.object(parallel_hs_dnam, 'CLOCK_FILE', Path('/configured/geo/result/qced_human_dnam_age.csv')):
            task = parallel_hs_dnam.generate_task_intention('GSE211704', {'title': 'example'})
        self.assertIn('`/configured/geo/gse_me_hs.csv`', task)
        self.assertIn('`/configured/geo/GEO_hs_dnam/processed_metadata_preharmonize/GSE211704.csv`', task)
        self.assertIn('`/configured/geo/result/qced_human_dnam_age.csv`', task)
        self.assertIn('aging_analysis_sonnet37.csv', task)
        self.assertNotIn('/dfs/', task)

    def test_check_detects_a_changed_body_and_an_unregistered_file(self):
        with tempfile.TemporaryDirectory() as directory, \
             mock.patch.object(build_prompts, 'PROMPTS', Path(directory)), \
             mock.patch('sys.stdout'), mock.patch('sys.stderr'):
            self.assertEqual(build_prompts.main([]), 0)
            self.assertEqual(build_prompts.main(['--check']), 0)
            path = Path(directory) / 'biological_interpreter_human_dnam.txt'
            path.write_text(path.read_text() + '\nAbridged')
            self.assertEqual(build_prompts.main(['--check']), 1)
            self.assertEqual(build_prompts.main([]), 0)
            (Path(directory) / 'unregistered.txt').write_text('Unregistered prompt')
            self.assertEqual(build_prompts.main(['--check']), 1)

    def test_historical_copies_match_the_manifest(self):
        root = REPO / 'src/historical_geo_pipeline'
        with (root / 'MANIFEST.tsv').open(newline='') as handle:
            records = list(csv.DictReader(handle, delimiter='\t'))
        self.assertTrue(records)
        copied = {path.relative_to(root).as_posix() for path in root.rglob('*')
                  if path.is_file() and path.name not in {'README.md', 'MANIFEST.tsv'}}
        self.assertEqual(copied, {record['relative_path'] for record in records})
        for record in records:
            with self.subTest(script=record['relative_path']):
                path = root / record['relative_path']
                self.assertIn(path.suffix, {'.py', '.R', '.r', '.sh', '.ipynb'})
                content = path.read_bytes()
                self.assertEqual(len(content), int(record['size']))
                self.assertEqual(hashlib.sha256(content).hexdigest(), record['sha256'])


if __name__ == '__main__':
    unittest.main()
