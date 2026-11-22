"""Software fixtures for release preparation; no empirical measurements."""
import json
import hashlib
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / 'scripts/release_check_macos_2026_11_12.sh'


def blocks():
    return re.findall(r"<<'PY'[^\n]*\n(.*?)^PY$", SCRIPT.read_text(), re.M | re.S)


class MacReleasePreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='clockbase-macos-software-fixture-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.release, self.bundle, self.work = [self.root / n for n in ('release', 'bundle', 'work')]
        for p in (self.release, self.bundle / 'clockbase_data/inputs', self.work):
            p.mkdir(parents=True)

    def prepare(self):
        return subprocess.run([sys.executable, '-B', '-', str(self.release), str(self.bundle), str(self.work)],
                              input=blocks()[1], text=True, capture_output=True)

    def test_prepare_creates_read_only_hardlinks_and_scoped_profile(self):
        source = self.bundle / 'clockbase_data/inputs/software_fixture.txt'
        source.write_text('SOFTWARE FIXTURE; not an empirical input\n')
        result = self.prepare()
        self.assertEqual(result.returncode, 0, result.stderr)
        copied = self.work / 'clockbase_data/inputs/software_fixture.txt'
        self.assertEqual(source.stat().st_ino, copied.stat().st_ino)
        self.assertEqual(copied.stat().st_mode & 0o222, 0)
        self.assertEqual(source.read_bytes(), copied.read_bytes())
        policy = (self.work / 'sandbox.sb').read_text()
        self.assertIn('(deny default)', policy)
        self.assertIn('(deny network*)', policy)
        paths = json.loads((self.work / 'sandbox_paths.json').read_text())
        self.assertEqual(paths['write_subpaths'], [str(self.work)])
        self.assertNotIn(str(self.bundle), paths['read_subpaths'])
        self.assertNotIn('/Users', paths['read_subpaths'])
        self.assertNotIn('/opt/homebrew', paths['read_subpaths'])
        self.assertIn('/System/Library/Perl', paths['read_subpaths'])
        if sys.platform == 'darwin':
            developer = Path(subprocess.check_output(['/usr/bin/xcode-select', '-p'], text=True).strip()).resolve()
            for relative in ('usr/bin', 'usr/lib', 'usr/libexec/git-core', 'usr/share/git-core'):
                self.assertIn(str(developer / relative), paths['read_subpaths'])
            self.assertIn('/Library/Preferences/com.apple.dt.Xcode.plist', paths['read_literals'])

    def test_prepare_rejects_symlink_back_to_author_data(self):
        (self.bundle / 'clockbase_data/inputs/escape').symlink_to(self.root)
        result = self.prepare()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('symlink', result.stderr)

    def test_path_validation_rejects_work_inside_release(self):
        result = subprocess.run([sys.executable, '-B', '-', str(self.release), str(self.bundle),
                                 str(self.release / 'work'), 'final-reference'],
                                input=blocks()[0], text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('separate directory trees', result.stderr)

    def test_all_embedded_python_compiles(self):
        self.assertEqual(len(blocks()), 4)
        for source in blocks():
            compile(source, str(SCRIPT), 'exec')

    def test_perl_timeout_fails_without_gnu_timeout(self):
        function = re.search(r'^LIMIT\(\) \{\n.*?^\}', SCRIPT.read_text(), re.M | re.S).group()
        result = subprocess.run(['/bin/bash', '-c', function + '\nLIMIT 1 /bin/sleep 10'],
                                env={**os.environ, 'PATH': '/nonexistent'}, capture_output=True, timeout=5)
        self.assertEqual(result.returncode, 124, result.stderr)

    def test_git_initializes_with_private_home_and_empty_templates(self):
        home, template, repo = (self.work / name for name in ('home', 'git-template', 'git-fixture'))
        for folder in (home, template, repo):
            folder.mkdir()
        env = {'HOME': str(home), 'TMPDIR': str(self.work),
               'PATH': '/Library/Developer/CommandLineTools/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin', 'GIT_CONFIG_NOSYSTEM': '1',
               'GIT_CONFIG_GLOBAL': os.devnull, 'GIT_TEMPLATE_DIR': str(template)}
        if os.environ.get('DEVELOPER_DIR'):  # the sandbox points the git shim at the Command Line Tools
            env['DEVELOPER_DIR'] = os.environ['DEVELOPER_DIR']
        # Use the shell's sandbox environment declaration so the test also checks
        # that the release command passes the isolated Git settings to children.
        declaration = re.search(r'^SANDBOX=\(.*?^  AUDIT_CACHE_ONLY=.*?\)$',
                                SCRIPT.read_text(), re.M | re.S).group()
        for setting in ('GIT_CONFIG_NOSYSTEM=1', 'GIT_CONFIG_GLOBAL=/dev/null',
                        'GIT_TEMPLATE_DIR="$WORK/git-template"'):
            self.assertIn(setting, declaration)
        result = subprocess.run(['git', 'init', '-q', str(repo)], env=env,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((repo / '.git/HEAD').is_file())
        self.assertFalse((repo / '.git/hooks').exists())

    def test_lfs_hydration_uses_only_hash_matched_local_objects(self):
        from scripts import build_release_2026_10_04 as builder
        data = b'SOFTWARE FIXTURE; not an empirical model\n'
        digest = hashlib.sha256(data).hexdigest()
        pointer = f'version https://git-lfs.github.com/spec/v1\noid sha256:{digest}\nsize {len(data)}\n'.encode()
        repo, release, local = (self.root / n for n in ('repo', 'archive', 'local-assets'))
        for folder in (repo, release, local):
            folder.mkdir()
        (repo / 'fixture.pkl').write_bytes(pointer)
        (release / 'fixture.pkl').write_bytes(pointer)
        (local / 'fixture.pkl').write_bytes(b'wrong object')
        with mock.patch.object(builder, 'REPO', repo):
            result = builder.hydrate_lfs(release, local)
            self.assertEqual(result['missing'], ['fixture.pkl'])
            self.assertEqual((release / 'fixture.pkl').read_bytes(), pointer)
            (local / 'fixture.pkl').write_bytes(data)
            result = builder.hydrate_lfs(release, local)
            self.assertEqual(result['missing'], [])
            self.assertEqual((release / 'fixture.pkl').read_bytes(), data)
            self.assertEqual((repo / 'fixture.pkl').read_bytes(), pointer)

    def test_exact_comparison_rejects_sub_tolerance_statistic_change(self):
        # This is an infrastructure fixture, not an empirical reproduction.
        (self.work / 'logs').mkdir()
        (self.work / 'logs/tests.log').write_text('Ran 1 tests in 0.1s\n\nOK\n')
        (self.work / 'sandbox_paths.json').write_text('{}')
        roots = [self.work / 'clockbase_data/verifier' / n for n in ('final-release-check', 'final-reference')]
        for root in roots:
            for m in ('dnam', 'hs_rna', 'mm_rna'):
                folder = root / f'repairs-{m}'
                folder.mkdir(parents=True)
                (folder / 'repairs.jsonl').write_text(json.dumps({'row_ordinal': 0, 'final_class': 'VERIFIED',
                    'recomputed': {'n_case': 2, 'n_control': 2, 'p': 0.25, 'fdr_family': 0.25}}) + '\n')
            (root / 'comparison-ranking').mkdir()
            (root / 'comparison-ranking/ranked_comparisons_all.csv').write_text('software_fixture\n1\n')
            (root / 'public-export/studies').mkdir(parents=True)
            (root / 'public-export/studies/fixture.json').write_text('{"software_fixture": true}')
        (roots[0] / 'verify.json').write_text('{"reproduced": true}')
        summary = self.work / 'clockbase_data/verifier/selftest-release-check/summary_SIMULATED.json'
        summary.parent.mkdir()
        summary.write_text(json.dumps({'simulation': True, 'planted_rows': 1, 'caught_any_route': 1, 'caught_as_expected': 1,
            'clean_rows': 1, 'clean_verified_unchanged': 1,
            'detached_simulation': {'by_treatment': {'fixture': {'EXPECTED': 1}}}}))
        report = self.work / 'release_check.json'
        argv = [sys.executable, '-B', '-', str(report), str(REPO), str(self.bundle), 'final-reference',
                '0', '0', '0', str(self.work), '0']
        result = subprocess.run(argv, input=blocks()[3], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        path = roots[0] / 'repairs-dnam/repairs.jsonl'
        row = json.loads(path.read_text())
        row['recomputed']['p'] += 1e-12
        path.write_text(json.dumps(row) + '\n')
        result = subprocess.run(argv, input=blocks()[3], text=True, capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(json.loads(report.read_text())['passed'])

    def test_summary_closure_recovers_historical_atlas_by_hash(self):
        from scripts import build_release_2026_10_04 as builder
        verifier = self.root / 'verifier'
        atlas = verifier / 'final-historical/atlas-v3/atlas_v3_combined.csv'
        atlas.parent.mkdir(parents=True)
        atlas.write_text('SOFTWARE FIXTURE; historical atlas\n')
        wrong = verifier / 'final-other/atlas-v3/atlas_v3_combined.csv'
        wrong.parent.mkdir(parents=True)
        wrong.write_text('SOFTWARE FIXTURE; different atlas\n')
        summary = verifier / 'classification-001/summary.json'
        summary.parent.mkdir()
        summary.write_text(json.dumps({'atlas_sha256': builder.sha256(atlas)}))
        inputs, records = builder.summary_input_closure([summary], {'verifier': verifier})
        self.assertEqual(inputs, {atlas})
        self.assertEqual(records[0]['sha256'], builder.sha256(atlas))

    def test_summary_closure_checks_nested_and_flat_inputs_recursively(self):
        from scripts import build_release_2026_10_04 as builder
        verifier = self.root / 'verifier'
        verifier.mkdir()
        data = verifier / 'software_fixture.csv'
        data.write_text('SOFTWARE FIXTURE; input\n')
        nested = verifier / 'historical_summary.json'
        nested.write_text(json.dumps({'inputs': {'ranked': {
            'path': 'verifier:software_fixture.csv', 'sha256': builder.sha256(data)}}}))
        summary = verifier / 'summary.json'
        summary.write_text(json.dumps({'input': 'verifier:historical_summary.json',
                                      'input_sha256': builder.sha256(nested),
                                      'outputs': {'unused': {'path': 'verifier:absent.csv', 'sha256': 'a' * 64}},
                                      'module_sha256': 'b' * 64}))
        inputs, records = builder.summary_input_closure([summary], {'verifier': verifier})
        self.assertEqual(inputs, {nested, data})
        self.assertEqual(len(records), 2)

    def test_summary_check_cannot_fall_back_to_author_data(self):
        from scripts import build_release_2026_10_04 as builder
        author = self.root / 'author/verifier'
        author.mkdir(parents=True)
        data = author / 'software_fixture.csv'
        data.write_text('SOFTWARE FIXTURE; author copy\n')
        deposit = self.root / 'deposit/verifier'
        deposit.mkdir(parents=True)
        summary = deposit / 'summary.json'
        summary.write_text(json.dumps({'input': str(data), 'input_sha256': builder.sha256(data)}))
        with mock.patch.dict(builder.paths.ROOTS, {'verifier': author}):
            with self.assertRaisesRegex(ValueError, 'Missing hashed input'):
                builder.summary_input_closure([summary], {'verifier': deposit})
        copied = deposit / data.name
        copied.write_bytes(data.read_bytes())
        builder.summary_input_closure([summary], {'verifier': deposit})
        copied.write_text('SOFTWARE FIXTURE; altered bytes\n')
        with self.assertRaisesRegex(ValueError, 'Missing hashed input'):
            builder.summary_input_closure([summary], {'verifier': deposit})

    def test_absent_input_is_recorded_only_when_the_author_data_lacks_it(self):
        from scripts import build_release_2026_10_04 as builder
        author = self.root / 'author-absent/verifier'
        author.mkdir(parents=True)
        deposit = self.root / 'deposit-absent/verifier'
        deposit.mkdir(parents=True)
        missing = author / 'never_written.csv'
        summary = deposit / 'summary.json'
        summary.write_text(json.dumps({'input': str(missing), 'input_sha256': '0' * 64}))
        with mock.patch.dict(builder.paths.ROOTS, {'verifier': author}):
            with self.assertRaisesRegex(ValueError, 'Missing hashed input'):
                builder.summary_input_closure([summary], {'verifier': deposit})
            _, records = builder.summary_input_closure([summary], {'verifier': author}, absent_ok=True)
            absent = {(r['recorded_path'], r['sha256']) for r in records if r.get('unavailable')}
            self.assertEqual(absent, {(str(missing), '0' * 64)})
            _, records = builder.summary_input_closure([summary], {'verifier': deposit}, absent_ok=absent)
            self.assertTrue(records[0]['unavailable'])
            with self.assertRaisesRegex(ValueError, 'Missing hashed input'):
                builder.summary_input_closure([summary], {'verifier': deposit}, absent_ok={('other', '0' * 64)})

    def test_generated_filter_keeps_historical_runs(self):
        from scripts import build_release_2026_10_04 as builder
        def excluded(name):
            return builder.generated_verifier_input(name, 'final-reference',
                {'arbitrary-rebuild-name'}, {'custom-selftest'})
        for name in ('final-old/input.csv', 'prefinal-old/input.csv', 'selftest-old/input.csv',
                     'final-reference/output.csv'):
            self.assertFalse(excluded(name), name)
        for name in ('arbitrary-rebuild-name/output.csv', 'another-release-trace/output.csv',
                     'custom-selftest/output.csv', 'pipeline-runs/log.json', 'pass-2/output.csv'):
            self.assertTrue(excluded(name), name)

    def test_bundler_includes_hashed_historical_input_and_excludes_new_outputs(self):
        from scripts import build_release_2026_10_04 as builder
        from scripts import clockbase_verifier as verifier_module
        from types import SimpleNamespace
        roots = {label: self.root / directory for label, directory in builder.BUNDLE_DIRS.items()}
        for root in roots.values():
            root.mkdir(exist_ok=True)
        for label in ('stream_empty', 'stream_detached'):
            (roots[label] / 'studies').mkdir()
        verifier = roots['verifier']
        reference = verifier / 'final-reference'
        reference.mkdir()
        historical = verifier / 'final-historical/atlas-v3/atlas_v3_combined.csv'
        historical.parent.mkdir(parents=True)
        historical.write_text('SOFTWARE FIXTURE; historical input\n')
        (reference / 'summary.json').write_text(json.dumps({'atlas_sha256': builder.sha256(historical)}))
        generated = verifier / 'arbitrary-rebuild-name/output.csv'
        generated.parent.mkdir()
        generated.write_text('SOFTWARE FIXTURE; generated output\n')
        prefinal = verifier / ('pre' + generated.parent.name) / 'output.csv'
        prefinal.parent.mkdir()
        prefinal.write_text('SOFTWARE FIXTURE; generated preliminary output\n')
        trace = self.root / 'file-access-lists'
        trace.mkdir()
        (trace / 'paths.txt').write_text(str(generated) + '\n' + str(prefinal) + '\n')
        (trace / 'rebuild-trace_run.json').write_text(json.dumps({'final': generated.parent.name, 'selftest': None}))
        args = SimpleNamespace(out=self.root / 'new-bundle', trace_dir=trace, reference=reference.name)
        env = {'VERIFIER_ROOT': verifier, 'STREAM_EMPTY': roots['stream_empty'],
               'STREAM_DETACHED': roots['stream_detached'], 'TRACES': roots['traces'],
               'GEO': roots['geo'], 'METHYLOME_REFERENCE': self.root / 'absent-reference'}
        with mock.patch.multiple(builder.paths, **env), mock.patch.dict(builder.paths.ROOTS, roots), \
             mock.patch.object(builder, 'AGENT_TRACES', []), mock.patch.object(builder, 'FIGURE_INPUTS', []), \
             mock.patch.object(verifier_module, 'stages', return_value=[]), \
             mock.patch.dict(sys.modules, {'clockbase_verifier': verifier_module}):
            builder.bundle(args)
        deposit = args.out / 'clockbase_data/verifier'
        self.assertTrue((deposit / historical.relative_to(verifier)).is_file())
        self.assertFalse((deposit / generated.relative_to(verifier)).exists())
        self.assertFalse((deposit / prefinal.relative_to(verifier)).exists())
        report = json.loads((args.out / 'BUNDLE.json').read_text())
        self.assertEqual(report['summary_input_hashes_checked'], 1)
        self.assertEqual(report['traced_run_outputs_skipped'], 2)
        # A bad copy must fail the bundle step before it writes a success manifest.
        args.out = self.root / 'corrupted-bundle'
        copy = builder.shutil.copy2
        def corrupt(src, dst):
            result = copy(src, dst)
            if Path(src).resolve() == historical.resolve():
                Path(dst).write_text('SOFTWARE FIXTURE; corrupted copy\n')
            return result
        with mock.patch.multiple(builder.paths, **env), mock.patch.dict(builder.paths.ROOTS, roots), \
             mock.patch.object(builder, 'AGENT_TRACES', []), mock.patch.object(builder, 'FIGURE_INPUTS', []), \
             mock.patch.object(verifier_module, 'stages', return_value=[]), \
             mock.patch.dict(sys.modules, {'clockbase_verifier': verifier_module}), \
             mock.patch.object(builder.shutil, 'copy2', side_effect=corrupt):
            with self.assertRaisesRegex(ValueError, 'Missing hashed input'):
                builder.bundle(args)
        self.assertFalse((args.out / 'clockbase_data/BUNDLE_MANIFEST.tsv').exists())
        self.assertFalse((args.out / 'BUNDLE.json').exists())

    def test_staging_links_historical_inputs_without_unrelated_finals(self):
        from scripts import build_release_2026_10_04 as builder
        from scripts import clockbase_verifier as verifier_module
        from types import SimpleNamespace
        verifier = self.root / 'author-verifier'
        atlas = verifier / 'final-historical/atlas-v3/atlas_v3_combined.csv'
        atlas.parent.mkdir(parents=True)
        atlas.write_text('SOFTWARE FIXTURE; prerequisite\n')
        summary = verifier / 'classification-001/summary.json'
        summary.parent.mkdir()
        summary.write_text(json.dumps({'atlas_sha256': builder.sha256(atlas)}))
        reference = verifier / 'final-reference'
        reference.mkdir()
        unrelated = verifier / 'final-unrelated/output.csv'
        unrelated.parent.mkdir()
        unrelated.write_text('SOFTWARE FIXTURE; unused historical output\n')
        target = self.root / 'staged/verifier'
        args = SimpleNamespace(reference=reference.name, out=target)
        with mock.patch.object(builder.paths, 'VERIFIER_ROOT', verifier), \
             mock.patch.dict(builder.paths.ROOTS, {'verifier': verifier}, clear=True), \
             mock.patch.object(builder, 'BUNDLE_DIRS', {'verifier': 'verifier'}), \
             mock.patch.object(verifier_module, 'stages', return_value=[]), \
             mock.patch.dict(sys.modules, {'clockbase_verifier': verifier_module}):
            builder.stage_inputs(args)
        staged_atlas = target / atlas.relative_to(verifier)
        self.assertEqual(staged_atlas.stat().st_ino, atlas.stat().st_ino)
        self.assertTrue((target / summary.relative_to(verifier)).is_file())
        self.assertFalse((target / unrelated.relative_to(verifier)).exists())


if __name__ == '__main__':
    unittest.main()
