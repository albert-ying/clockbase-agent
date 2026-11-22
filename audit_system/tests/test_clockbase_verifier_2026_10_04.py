import contextlib
import io
import json
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from audit_system import paths
from audit_system import label_adjudication_2026_09_27 as adjudication
from scripts import clockbase_verifier as verifier
from audit_system import verifier_selftest_2026_10_01 as selftest


class ModelStagesTest(unittest.TestCase):
    def test_row_comparison_retains_all_recomputed_statistics(self):
        # Synthetic record: each field has a distinct value so omissions show.
        fields = ('n_case', 'n_control', 'mean_case', 'mean_control', 'mean_diff',
                  'cohen_d', 't', 'p', 'fdr_family', 'fdr_family_size', 'undefined_reason')
        record = {'row_ordinal': 0, 'final_class': 'DESCRIPTIVE',
                  'recomputed': {k: i for i, k in enumerate(fields)}}
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'repairs.jsonl'
            path.write_text(json.dumps(record) + '\n')
            self.assertEqual(verifier.load_repairs(path)[0], ('DESCRIPTIVE', *range(len(fields))))

    def test_commands_use_the_saved_pairs_and_providers(self):
        stages = verifier.stages(Path('verifier'), 'final', Path('empty'), Path('detached'))
        count = 0
        for stage in stages:
            for command, env, provider in stage.commands:
                if not provider:
                    self.assertNotIn('AUDIT_MODELS', env)
                    continue
                openai = stage.name in ('design-audit', 'design-audit-added', 'design-audit-design-estimand',
                                        'individuals-audit') or 'supported_sample' in command
                self.assertEqual(provider, 'openai' if openai else 'tokenswitch')
                self.assertEqual(env['AUDIT_MODELS'], 'openai:gpt-6-astra,openai:gpt-6-sol' if openai else 'gp-6-as,cl-op-5')
                count += 1
        self.assertEqual(count, 20)

    def test_pass_3_precedes_its_model_stages_and_the_final_chain(self):
        names = [s.name for s in verifier.stages(Path('verifier'), 'final-x', Path('empty'), Path('detached'))]
        order = [names.index(n) for n in ('design-audit-added', 'pass-3', 'design-audit-design-estimand', 'individuals-audit',
                                          'final-chain', 'decision-audit')]
        self.assertEqual(order, sorted(order))
        stage = next(s for s in verifier.stages(Path('verifier'), 'final-x', Path('empty'), Path('detached')) if s.name == 'pass-3')
        self.assertEqual(str(stage.done), 'verifier/prefinal-x/public-export/manifest.json')

    def test_live_aliases_use_the_stage_provider_in_the_launcher(self):
        launcher = ['credential tool', '--provider={provider}', '--']
        with mock.patch.object(paths, 'MODEL_LAUNCHER', launcher):
            for mode in ('live', 'openai', 'neolink'):
                for provider in ('tokenswitch', 'openai'):
                    env, prefix = verifier.model_env(mode, provider)
                    self.assertEqual(prefix, ['credential tool', f'--provider={provider}', '--'])
                    self.assertEqual(env['AUDIT_CACHE_ONLY'], '0')
            self.assertEqual(verifier.model_env('cache-only', 'openai'), ({'AUDIT_CACHE_ONLY': '1'}, []))
        with mock.patch.object(paths, 'MODEL_LAUNCHER', []):
            self.assertEqual(verifier.model_env('live', 'tokenswitch')[1], [])

    def test_cache_only_does_not_send_a_missing_request(self):
        ledger = SimpleNamespace(cache={})
        env, _ = verifier.model_env('cache-only', 'tokenswitch')
        with mock.patch.dict(os.environ, env, clear=True), mock.patch.object(adjudication.urllib.request, 'urlopen') as send:
            for model in ('gp-6-as', 'cl-op-5', 'openai:gpt-6-astra', 'openai:gpt-6-sol'):
                result, cached = adjudication.call(model, 'cache-only check', ledger, {})
                self.assertEqual(result['status'], 'NOT_CACHED')
                self.assertFalse(cached)
            send.assert_not_called()

    def test_mixed_followup_records_models_and_runs_touch_without_a_launcher(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            stage = next(s for s in verifier.stages(root, 'final', root / 'empty', root / 'detached')
                         if s.name == 'adjudication-followup-dnam')
            args = SimpleNamespace(root=root, final='final', stream=root / 'empty', stream_detached=root / 'detached',
                                   models='live', plan_only=False, force=[])
            with mock.patch.object(verifier, 'stages', return_value=[stage]), \
                    mock.patch.object(paths, 'MODEL_LAUNCHER', ['launcher', '{provider}', '--']), \
                    mock.patch.dict(os.environ, {'AUDIT_MODELS': 'unrelated', 'AUDIT_CACHE_ONLY': '1'}), \
                    mock.patch.object(verifier.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run, \
                    contextlib.redirect_stdout(io.StringIO()):
                verifier.run_stages(args)
            calls = run.call_args_list
            for call, provider, models in zip(calls[:3], ('tokenswitch', 'tokenswitch', 'openai'),
                                              ('gp-6-as,cl-op-5', 'gp-6-as,cl-op-5', 'openai:gpt-6-astra,openai:gpt-6-sol')):
                self.assertEqual(call.args[0][:3], ['launcher', provider, '--'])
                self.assertEqual(call.kwargs['env']['AUDIT_MODELS'], models)
                self.assertEqual(call.kwargs['env']['AUDIT_CACHE_ONLY'], '0')
            self.assertEqual(calls[3].args[0], ['touch', str(stage.done)])
            manifest = json.loads(next((root / 'pipeline-runs').glob('*/manifest.json')).read_text())
            recorded = manifest['stages'][0]['model_commands']
            self.assertEqual([c['provider'] for c in recorded], ['tokenswitch', 'tokenswitch', 'openai'])
            self.assertEqual([c['models'] for c in recorded], [['gp-6-as', 'cl-op-5'], ['gp-6-as', 'cl-op-5'],
                                                              ['openai:gpt-6-astra', 'openai:gpt-6-sol']])


class PlanTest(unittest.TestCase):
    def test_plan_does_not_create_an_absent_root_or_run_commands(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder) / 'absent'
            with mock.patch.object(Path, 'mkdir', side_effect=AssertionError('plan attempted mkdir')), \
                    mock.patch.object(verifier.subprocess, 'run') as run, contextlib.redirect_stdout(io.StringIO()):
                verifier.main(['plan', '--root', str(root), '--final', 'final', '--stream', str(root / 'empty'),
                               '--stream-detached', str(root / 'detached')])
            self.assertFalse(root.exists())
            run.assert_not_called()

    def test_cli_defaults_use_paths_and_cache_only(self):
        with mock.patch.object(verifier, 'run_stages') as run:
            verifier.main(['plan', '--final', 'final'])
        args = run.call_args.args[0]
        self.assertEqual(args.root, paths.VERIFIER_ROOT)
        self.assertEqual(args.stream, paths.STREAM_EMPTY)
        self.assertEqual(args.stream_detached, paths.STREAM_DETACHED)
        self.assertEqual(args.models, 'cache-only')
        self.assertTrue(args.plan_only)


class SelftestTests(unittest.TestCase):
    def test_stages_use_saved_predictions_and_final_repair_policy(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            for name in ('classification-001', 'individuals-audit-001'):
                (root / name).mkdir()
                (root / name / 'comparisons.jsonl').write_text('')
            with mock.patch.object(selftest.subprocess, 'run', return_value=SimpleNamespace(returncode=0)) as run:
                selftest.run_stages(root, root, root / 'corrections', 1)
            commands = [call.args[0] for call in run.call_args_list]
            self.assertEqual(len(commands), 15)
            for i, modality in enumerate(selftest.TABLES):
                recovery, repair = commands[i * 5], commands[i * 5 + 4]
                self.assertEqual(recovery[recovery.index('--prediction-replay') + 1],
                                 str(root / f'replay-{modality}-001'))
                for flag in ('--design-estimand', '--clock-output-decisions', '--classification', '--individuals'):
                    self.assertIn(flag, repair)
                self.assertEqual('--repaired-dnam' in repair, modality == 'dnam')
                self.assertNotIn('--adjudication', repair)  # planted inversions must await review
                self.assertNotIn('--design', commands[i * 5 + 2])  # never copy reference arms by row ordinal

    def test_expectations_stay_strict_for_clean_rows_and_duplicate_pairs(self):
        reference = {'mean_diff': 1., 'cohen_d': 0.5, 'p': 0.1, 'n_case': 4, 'n_control': 4}
        entry = {'defect': 'clean', 'reference': reference}
        verdict = {'reasons': [], 'verdict': 'VERIFIED'}
        repair = {'final_class': 'VERIFIED', 'recomputed': reference.copy(), 'actions': []}
        self.assertEqual(selftest.score(entry, verdict, repair), 'EXPECTED')
        self.assertEqual(selftest.score(entry, verdict, {**repair, 'final_class': 'DESCRIPTIVE'}), 'FALSE_ALARM')
        self.assertEqual(selftest.score({**entry, 'defect': 'duplicate'}, verdict, repair, 'DUPLICATE'), 'EXPECTED')
        self.assertEqual(selftest.score({**entry, 'defect': 'duplicate'}, verdict, repair, 'VERIFIED'), 'MISSED')
        self.assertEqual(selftest.score({**entry, 'defect': 'effect_scaled'}, verdict, repair), 'CORRECTED_UNFLAGGED')


if __name__ == '__main__':
    unittest.main()
