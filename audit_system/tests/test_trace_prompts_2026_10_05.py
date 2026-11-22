"""Exact task-template recovery from the saved human and mouse RNA agent traces."""
import hashlib
import pickle
import re
import subprocess
import sys
import unittest
from collections import Counter
from pathlib import Path
from unittest import mock

from audit_system import paths
from src import parallel_hs_rna, parallel_mm_rna


def first_human_message(entry):
    """Remove the first human wrapper, including its blank-line separator."""
    for index, message in enumerate(entry['result'][0]):
        if not isinstance(message, str):
            continue
        wrapper, separator, body = message.partition('\n\n')
        if separator and re.fullmatch(r'=+ Human Message =+', wrapper):
            return index, body
    raise AssertionError('Saved trace has no first human message in result[0]')


def study_values(body):
    """Extract study-specific text without changing its saved representation."""
    study_id = re.search(r'from study (GSE\d+)', body).group(1)
    header = re.search(r'2\. Study Design Information[^\n]*:\n   ', body).group(0)
    before, _, rest = body.partition(header)
    design, separator, after = rest.partition('\n\n3. Clock Selection and Analysis:')
    if not separator:
        raise AssertionError('Saved task has no section-3 boundary after its study design')
    data_dir = re.search(r'overview from: `(.+)/dataset_level/', before).group(1)
    results_dir = re.search(r'Create output directory if needed: `(.+)/' + study_id + r'/`', before).group(1)
    output_dir = re.search(r'Save to: (.+)/' + study_id + r'/aging_analysis', after).group(1)
    return {'study_id': study_id, 'study_design': design, 'data_dir': data_dir,
            'RESULTS_DIR': results_dir, 'OUTPUT_DIR': output_dir}


class TracePromptsTests(unittest.TestCase):
    def check_saved_prompts(self, driver, expected_counts, expected_offsets):
        path = paths.TRACES / driver.TRACE_SOURCE
        if not path.is_file():
            self.skipTest(f'saved agent traces not present: {path}')
        try:
            with path.open('rb') as handle:
                self.assertEqual(hashlib.file_digest(handle, 'sha256').hexdigest(), driver.TRACE_SHA256)
                handle.seek(0)
                traces = pickle.load(handle)
        except (OSError, pickle.UnpicklingError, EOFError) as exc:
            self.fail(f'Cannot read saved agent traces from {path}: {exc}')
        self.assertEqual(len(traces), sum(expected_counts.values()))
        counts = Counter()
        offsets = Counter()
        for key, entry in traces.items():
            with self.subTest(trace=key):
                index, body = first_human_message(entry)
                values = study_values(body)
                matches = [name for name, template in driver.TASK_TEMPLATES.items()
                           if template.format(**values) == body]
                self.assertEqual(len(matches), 1, f'Expected exactly one released template, found {matches}')
                counts[matches[0]] += 1
                offsets[index] += 1
        self.assertEqual(dict(counts), expected_counts)
        self.assertEqual(dict(counts), driver.TRACE_COUNTS)
        self.assertEqual(dict(offsets), expected_offsets)

    def test_human_templates_equal_every_saved_first_human_message(self):
        self.check_saved_prompts(parallel_hs_rna, {'sonnet37': 9546}, {0: 9546})

    def test_mouse_templates_equal_every_saved_first_human_message(self):
        self.check_saved_prompts(parallel_mm_rna,
                                {'base': 11385, 'all_samples_local': 3578, 'all_samples': 417},
                                {0: 15367, 1: 13})

    def test_human_driver_renders_the_traced_template(self):
        design = {'title': 'Study {with braces}', 'overall_design': 'first line\nsecond line'}
        expected = parallel_hs_rna.TASK_TEMPLATE.format(
            study_id='GSE123', study_design=design, data_dir='/configured/geo',
            RESULTS_DIR='/configured/results')
        expected = expected.replace('/configured/geo/dataset_level/gse_rna_hs.csv',
                                    '/configured/geo/gse_rna_hs.csv')
        expected = expected.replace('/configured/geo/clock_output/human_rna_age_short.csv',
                                    '/configured/geo/result/human_rna_age.csv')
        with mock.patch.object(paths, 'GEO', Path('/configured/geo')), \
             mock.patch.object(parallel_hs_rna, 'RESULTS_DIR', '/configured/results'), \
             mock.patch.object(parallel_hs_rna, 'CLOCK_FILE', Path('/configured/geo/result/human_rna_age.csv')):
            self.assertEqual(parallel_hs_rna.generate_task_intention('GSE123', design), expected)
        self.assertIn('using SONNET37.', expected)
        self.assertIn('model_version: SONNET37', expected)
        self.assertIn('/GSE123/aging_analysis_sonnet37.csv', expected)
        self.assertNotIn('/dfs/', expected)

    def test_mouse_driver_renders_each_traced_template_and_defaults_to_base(self):
        design = {'title': 'Study {with braces}', 'overall_design': 'first line\nsecond line'}
        with mock.patch.object(paths, 'GEO', Path('/configured/geo')), \
             mock.patch.object(parallel_mm_rna, 'RESULTS_DIR', '/configured/results'), \
             mock.patch.object(parallel_mm_rna, 'CLOCK_FILE', Path('/configured/geo/result/mouse_rna_age.csv')):
            for name, template in parallel_mm_rna.TASK_TEMPLATES.items():
                with self.subTest(template=name):
                    expected = template.format(study_id='GSE123', study_design=design,
                                               data_dir='/configured/geo', RESULTS_DIR='/configured/results',
                                               OUTPUT_DIR='/configured/results')
                    expected = expected.replace('/configured/geo/dataset_level/gse_rna_mm.csv',
                                                '/configured/geo/gse_rna_mm.csv')
                    expected = expected.replace('/configured/geo/clock_output/mouse_rna_age_short.csv',
                                                '/configured/geo/result/mouse_rna_age.csv')
                    self.assertEqual(parallel_mm_rna.generate_task_intention('GSE123', design, name), expected)
                    self.assertNotIn('/dfs/', expected)
            self.assertEqual(parallel_mm_rna.generate_task_intention('GSE123', design),
                             parallel_mm_rna.generate_task_intention('GSE123', design, 'base'))
        with self.assertRaisesRegex(ValueError, 'Unknown mouse RNA task template'):
            parallel_mm_rna.generate_task_intention('GSE123', design, 'unknown')

    def test_mouse_templates_preserve_separate_creation_and_save_directories(self):
        for name, template in parallel_mm_rna.TASK_TEMPLATES.items():
            with self.subTest(template=name):
                task = template.format(study_id='GSE123', study_design={}, data_dir='/configured/geo',
                                       RESULTS_DIR='/configured/create', OUTPUT_DIR='/configured/save')
                self.assertIn('Create output directory if needed: `/configured/create/GSE123/`', task)
                self.assertIn('Save to: /configured/save/GSE123/aging_analysis.csv', task)

    def test_drivers_import_without_loading_model_or_network_clients(self):
        script = '''
import builtins
original_import = builtins.__import__
def checked_import(name, *args, **kwargs):
    if name in {'coding_agent', 'src.coding_agent', 'parallel_processing',
                'src.parallel_processing', 'requests', 'openai', 'anthropic'}:
        raise AssertionError('Driver imported a model or network client: ' + name)
    return original_import(name, *args, **kwargs)
builtins.__import__ = checked_import
from src import parallel_hs_rna, parallel_mm_rna
assert parallel_hs_rna.generate_task_intention('GSE123', {})
assert parallel_mm_rna.generate_task_intention('GSE123', {})
'''
        result = subprocess.run([sys.executable, '-B', '-c', script], cwd=paths.REPO,
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
