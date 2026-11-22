#!/usr/bin/env python3
"""Build the code, rebuild-input and verification-record deposits, and check them.

  code    the curated source release: the analysis system, prompts, verifier, tests and scripts from one commit, with a
          file manifest and a scan for private paths and credentials
  bundle  the rebuild inputs: files a cache-only rebuild opened (traced with scripts/release_trace), hash-recorded
          data inputs of the included run summaries, the
          completion files of every reused stage, the saved agent traces and the reference final run, laid out as
          clockbase_data/ (see audit_system/paths.py), with a hash manifest; also the DNA methylation re-derivation
          outputs of every attempted study and their processed sample metadata
  records the final public export and the recovered source tables, scripts, summaries and receipts cited by the
          manuscript, copied unchanged with a statement index and hash manifest

  python scripts/build_release_2026_10_04.py code   --commit HEAD --out DIR
  python scripts/build_release_2026_10_04.py trace  --final NAME --trace-dir DIR      (runs the cache-only rebuild traced)
  python scripts/build_release_2026_10_04.py stage-inputs --reference NAME --out DIR (hard-links verifier prerequisites)
  python scripts/build_release_2026_10_04.py trace  --final REFERENCE --selftest ROOT --trace-dir DIR
                                                     (runs the planted-defect self-test against REFERENCE traced)
  python scripts/build_release_2026_10_04.py bundle --trace-dir DIR --reference NAME --out DIR
  python scripts/build_release_2026_10_04.py records --final NAME --out DIR

The sandboxed check that a reader can rebuild from these two deposits alone is scripts/release_check_2026_10_04.sh.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from audit_system import paths  # noqa: E402

# Released source: the analysis system as a reader needs it. Peer-review files, manuscripts, data and generated outputs
# (paper/, result*/, figures/, *_sonnet37/, pickles) are not code; the agent traces travel in the bundle instead.
CODE_INCLUDE = ['audit_system/', 'scripts/', 'src/', 'prompts/', 'minipaper_fork/', 'pipeline_assets/', 'notebooks/',
                'PIPELINE_2026_09_10.md', 'README.md', 'setup_pipeline_2026_09_10.sh', 'run_pipeline.sh',
                'human_prompts_for_minipaper.py', 'requirements-verifier.txt', 'pipeline_requirements_2026_09_10.txt', 'LICENSE',
                '.gitattributes', 'questionaire.md', 'minipaper_figures_2026-07-01/pipeline.py',
                'minipaper_figures_2026-07-01/README.md', 'src/historical_geo_pipeline/analysis/',
                'scripts/legacy_sensitivity_analyses_2026_05_14/', 'benchmark/']
CODE_ROOT_SCRIPTS = re.compile(r'^[^/]+\.py$')  # figure and analysis scripts at the repository root
# Also excluded: the tooling that checks the manuscript's statements against the files, and scripts/scheduled_rebuild.py, an unused prototype of a scheduled
# rebuild (its GitHub Actions workflow never existed; ClockBase runs no scheduled updates).
CODE_EXCLUDE = re.compile(r'(^|/)(_archive|__pycache__|\.clockbase-cache|outdated)(/|$)|\.(parquet|csv\.gz|zip|tar\.gz)$|'
                          r'^src/.*\.csv$|^src/Untitled|^audit_system/promises_[^/]*\.json$|^audit_system/(tests/test_)?promise_audit_[^/]*\.py$|'
                          r'^scripts/scheduled_rebuild\.py$|'
                          # The fitted transcriptomic clock models belong to collaborators and are distributed
                          # separately; diffs against earlier figure scripts are working files, not code.
                          r'^pipeline_assets/rna_4_6_2026_09_10/[^/]*\.pkl$|^scripts/figures_2026_11_06_july_exact/[^/]*\.diff$')
# Names that belong to the development checkout (its folder, its manuscript folder, the location of the figure-style
# module on the development machine) are replaced in the released text files by the names the release uses. The list
# is kept beside the manuscript sources, outside the release.
RELEASE_RENAMES_FILE = REPO / 'paper' / 'release_renames_2026_11_22.json'
RELEASE_RENAMES = ([tuple(pair) for pair in json.loads(RELEASE_RENAMES_FILE.read_text())]
                   if RELEASE_RENAMES_FILE.is_file() else [])
PRIVATE_PATH = re.compile(r'/Entropy/|/home/kying|/Users/A\.Y|/dfs/user/|/Lasagna/')
# API keys start with a lowercase 'sk-'; the uppercase 'SK-' of cell-line names (SK-N-SH, SK-MEL-28) is not a key.
# Key prefixes are case-sensitive: case-insensitive 'AKIA' matched inside base64 signatures of Claude Code logs.
CREDENTIAL = re.compile(r'(?-i:sk-[A-Za-z0-9_-]{20,})|Bearer\s+[A-Za-z0-9._-]{20,}|(?-i:AKIA[0-9A-Z]{16})|ghp_[A-Za-z0-9]{30,}|'
                        r'(api[_-]?key|secret|password|token)\s*[:=]\s*["\'][^"\'\s]{12,}["\']', re.I)
NOT_A_CREDENTIAL = re.compile(r'fixture|example|placeholder|REDACTED|token\s*=\s*["\'][-+]?[0-9.eE]+["\']', re.I)
DETECTION_PATTERN = re.compile(r're\.(compile|search|match)\(|PRIVATE_PATH|_RE\b|pattern', re.I)
BUNDLE_DIRS = {'geo': 'GEO_data', 'verifier': 'verifier', 'stream_empty': 'dnam-stream-repair-empty',
               'stream_detached': 'dnam-stream-repair-detached', 'corrections': 'corrections', 'inputs': 'inputs',
               'traces': 'agent-traces'}
AGENT_TRACES = ['successful_traces_human_dnam.pkl', 'successful_traces_human_rna.pkl', 'successful_traces_mouse.pkl',
                'hs_dnam_sonnet37', 'hs_rna_sonnet37', 'result_hs_rna', 'result_mm_rna']
# The DNA methylation re-derivation of every study it attempted, not only the studies the verifier reads: its status,
# sample QC, clock outputs and coverage (which the website serves and the manuscript cites), and the
# processed sample metadata of those studies (recorded age for the comparisons with clock outputs).
REDERIVATION_FILES = {'status.json', 'sample_status.csv', 'sample_map.csv', 'clock_predictions.csv',
                      'clock_coverage.csv', 'log.txt'}
REDERIVATION_SUMMARIES = ['inventory.csv', 'summary.csv', 'summary.json', 'run.json']
# Revision figure inputs outside the reference final run (which bundle() already copies in full).
FIGURE_INPUTS = [
    'verifier:r6-6-blood/structure-002/summary.json',
    'verifier:r6-6-blood/structure-002/cohort.csv',
    'verifier:r6-6-blood/structure-002/correlation_repaired_residual_spearman_historical_horvath_age_input.csv',
    'verifier:r6-6-blood/structure-002/correlation_repaired_residual_spearman_chronological_age_input.csv',
    'verifier:r6-6-blood/grimage-age-input-002/grimage_age_input.csv',
    'verifier:input-integrity-dnam-001/samples.csv',
    'verifier:replay-dnam-001/inputs/GSE73950.csv',
    'geo:GEO_hs_dnam/metadata/GSE73950.csv',
]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_manifest(root, name):
    """Tab-separated path, sha256 and bytes; a path with a tab, newline or quote is CSV-quoted (one agent output
    directory is named after a printed pandas Series), so any CSV reader with delimiter='\\t' parses every row."""
    rows = []
    for path in sorted(p for p in Path(root).rglob('*') if p.is_file()):
        rows.append((path.relative_to(root).as_posix(), sha256(path), path.stat().st_size))
    with open(Path(root) / name, 'w', newline='') as sink:
        writer = csv.writer(sink, delimiter='\t', lineterminator='\n')
        writer.writerow(('path', 'sha256', 'bytes'))
        writer.writerows(rows)
    return rows


def scan(root):
    """Private paths and credential-like strings in text files of the release."""
    findings = {'private_path': [], 'detection_pattern': [], 'credential': []}
    for path in sorted(p for p in Path(root).rglob('*') if p.is_file()):
        try:
            text = path.read_text(errors='strict')
        except (UnicodeDecodeError, OSError):
            continue
        for number, line in enumerate(text.splitlines(), 1):
            where = f'{path.relative_to(root).as_posix()}:{number}: {line.strip()[:160]}'
            if CREDENTIAL.search(line) and not NOT_A_CREDENTIAL.search(line):
                findings['credential'].append(where)
            if PRIVATE_PATH.search(line):
                findings['detection_pattern' if DETECTION_PATTERN.search(line) else 'private_path'].append(where)
    return findings


def hydrate_lfs(out, extra_source=None):
    """Replace LFS pointers with hash-matched local objects, without downloading."""
    replaced, missing = [], []
    for path in sorted(p for p in out.rglob('*') if p.is_file() and p.stat().st_size < 1024):
        text = path.read_bytes()
        if not text.startswith(b'version https://git-lfs.github.com/spec/v1'):
            continue
        oid = re.search(rb'oid sha256:([0-9a-f]{64})', text).group(1).decode()
        candidates = [REPO / path.relative_to(out)]
        if extra_source:
            candidates.append(Path(extra_source) / path.relative_to(out))
        source = next((p for p in candidates if p.is_file() and sha256(p) == oid), None)
        if source is not None:
            subprocess.run(['/bin/cp', '-f', str(source), str(path)], check=True)
            replaced.append(path.relative_to(out).as_posix())
        else:
            missing.append(path.relative_to(out).as_posix())
    return {'replaced': replaced, 'missing': missing}


def build_code(args):
    commit = subprocess.run(['git', 'rev-parse', args.commit], cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()
    listing = ['git', 'ls-files'] if args.worktree else ['git', 'ls-tree', '-r', '--name-only', commit]
    files = subprocess.run(listing, cwd=REPO, capture_output=True, text=True, check=True).stdout.splitlines()
    if args.worktree:  # include the new recovered-source directories before they are staged by the release worker
        for name in ('src/historical_geo_pipeline/analysis', 'scripts/legacy_sensitivity_analyses_2026_05_14',
                     'audit_system/policy'):
            files.extend(p.relative_to(REPO).as_posix() for p in (REPO / name).rglob('*') if p.is_file())
        for name in ('audit_system/final_verification_rules_2026_11_13.py',
                     'audit_system/tests/test_final_verification_rules_2026_11_13.py'):
            if (REPO / name).is_file():
                files.append(name)
        if (REPO / 'scripts/records_README_2026_10_05.md').is_file():
            files.append('scripts/records_README_2026_10_05.md')
        files = sorted(set(files))
    chosen = [f for f in files if (any(f == inc or f.startswith(inc) for inc in CODE_INCLUDE) or CODE_ROOT_SCRIPTS.match(f))
              and not CODE_EXCLUDE.search(f)]
    out = Path(args.out) / (f'clockbase-worktree-{commit[:12]}' if args.worktree else f'clockbase-{commit[:12]}')
    out.mkdir(parents=True, exist_ok=False)
    if args.worktree:  # tracked files as they are in the working tree (for testing uncommitted changes)
        for name in chosen:
            if (REPO / name).is_file():
                (out / name).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(REPO / name, out / name)
    else:
        archive = subprocess.run(['git', 'archive', '--format=tar', commit, '--', *chosen], cwd=REPO, capture_output=True,
                                 env={**os.environ, 'GIT_LFS_SKIP_SMUDGE': '1'}, check=True)
        tar_path = Path(args.out) / f'clockbase-{commit[:12]}.tar'
        tar_path.write_bytes(archive.stdout)
        with tarfile.open(tar_path) as handle:
            handle.extractall(out, filter='data')
        tar_path.unlink()
    lfs = hydrate_lfs(out, args.lfs_source)
    if lfs['missing']:
        raise RuntimeError(f"Missing local LFS objects for {len(lfs['missing'])} release files; no download attempted")
    # The release ships the hydrated files, and the other attribute rules name paths of the development repository,
    # so the attributes file is not part of the release.
    (out / '.gitattributes').unlink(missing_ok=True)
    renamed = {name: 0 for name, _ in RELEASE_RENAMES}
    for path in sorted(p for p in out.rglob('*') if p.is_file()):
        try:
            text = path.read_text(encoding='utf-8')
        except (UnicodeDecodeError, OSError):
            continue
        changed = text
        for name, replacement in RELEASE_RENAMES:
            renamed[name] += changed.count(name)
            changed = changed.replace(name, replacement)
        if changed != text:
            path.write_text(changed, encoding='utf-8')
    rows = write_manifest(out, 'RELEASE_MANIFEST.tsv')
    findings = scan(out)
    report = {'commit': commit, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'files': len(rows),
              'bytes': sum(r[2] for r in rows), 'include': CODE_INCLUDE, 'exclude_pattern': CODE_EXCLUDE.pattern,
              'private_path_lines': len(findings['private_path']), 'detection_pattern_lines': len(findings['detection_pattern']),
              'credential_lines': len(findings['credential']), 'lfs_objects': lfs,
              'renamed_strings': renamed,
              'findings': findings}
    (Path(args.out) / f'clockbase-{commit[:12]}.release.json').write_text(json.dumps(report, indent=2))
    with tarfile.open(Path(args.out) / f'clockbase-{commit[:12]}.tar.gz', 'w:gz') as handle:
        handle.add(out, arcname=out.name)
    print(json.dumps({k: v for k, v in report.items() if k != 'findings'}, indent=1))
    raise SystemExit(1 if findings['credential'] else 0)


def trace(args):
    """Run the cache-only rebuild (a new final directory under the verifier root), or the self-test against a
    reference final run, with the file-access tracer; the bundle ships every file either opened."""
    trace_dir = Path(args.trace_dir)
    trace_dir.mkdir(parents=True, exist_ok=False)
    # Distinguish this invocation's outputs from historical final inputs.
    (trace_dir / 'trace_run.json').write_text(json.dumps({
        'final': args.final, 'selftest': paths.portable(args.selftest) if args.selftest else None,
    }, indent=2) + '\n')
    env = {**os.environ, 'CLOCKBASE_TRACE_DIR': str(trace_dir),
           'PYTHONPATH': os.pathsep.join(filter(None, [str(REPO / 'scripts/release_trace'), os.environ.get('PYTHONPATH')]))}
    command = [sys.executable, '-B', str(REPO / 'scripts/clockbase_verifier.py'), 'run', '--final', args.final,
               '--models', 'cache-only']
    if args.selftest:
        command = [sys.executable, '-B', str(REPO / 'scripts/clockbase_verifier.py'), 'selftest', '--root',
                   str(args.selftest), '--reference', str(paths.VERIFIER_ROOT / args.final), '--per-modality', '160']
    raise SystemExit(subprocess.run(command, cwd=REPO, env=env).returncode)


SERIES_MATRIX = re.compile(r'_series_matrix[^/]*\.txt\.gz$')


def copy_series_matrix_header(source, target):
    """The verifier reads a series matrix only up to !series_matrix_table_begin (sample annotations); copy that part."""
    import gzip
    with gzip.open(source, 'rt', newline='', errors='replace') as reader, gzip.open(target, 'wt', newline='') as writer:
        for line in reader:
            writer.write(line)
            if line.startswith('!series_matrix_table_begin'):
                writer.write('!series_matrix_table_end\n')
                break



def sha256_file(path, chunk=1 << 20):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(chunk), b''):
            digest.update(block)
    return digest.hexdigest()


def generated_verifier_input(relative, reference, trace_finals=(), selftest_roots=()):
    """Exclude the traced outputs, while retaining historical final inputs."""
    relative = Path(relative)
    if not relative.parts or relative.parts[0] == reference:
        return False
    top = relative.parts[0]
    return (top in trace_finals or top.endswith('-release-trace') or top == 'pipeline-runs'
            or top.startswith('pass-')
            or any(relative.is_relative_to(Path(root)) for root in selftest_roots))


# Pathless hashes written by the older summary producers. A filename match
# never suffices without the complete recorded SHA-256.
SUMMARY_INPUT_NAMES = {
    'atlas': ('atlas_v3_combined.csv', 'atlas_repaired.csv'),
    'classification': ('comparisons.jsonl',), 'membership': ('rows.jsonl',),
    'repairs': ('repairs.jsonl',), 'ranking': ('ranked_comparisons_all.csv',),
    'projects': ('projects.json',), 'cohort': ('cohort.csv',),
    'genage': ('genage_models.csv',), 'drugage': ('drugage.csv',),
    'meta': ('intervention_meta.csv',),
}


def summary_inputs(summary):
    """Yield data-input hashes, excluding outputs and source-code identities."""
    try:
        data = json.loads(Path(summary).read_text())
    except (OSError, ValueError) as error:
        raise ValueError(f'Cannot read input provenance {summary}: {error}') from error

    def walk(value, context=()):
        if isinstance(value, list):
            for i, item in enumerate(value):
                yield from walk(item, (*context, str(i)))
        elif isinstance(value, dict):
            for key, item in value.items():
                if key in ('outputs', 'pipeline_receipts'):
                    continue
                if isinstance(item, (dict, list)):
                    yield from walk(item, (*context, key))
                    continue
                if not isinstance(item, str) or not re.fullmatch('[0-9a-f]{64}', item):
                    continue
                recorded, names = None, ()
                if key == 'sha256' and isinstance(value.get('path'), str):
                    recorded = value['path']
                elif key.endswith('_sha256'):
                    stem = key[:-7]
                    if stem == 'output' or any(t in stem for t in ('module', 'ast', 'heuristic')):
                        continue
                    recorded = next((value[k] for k in (stem, stem + '_path', stem + '_source')
                                     if isinstance(value.get(k), str)), None)
                    names = SUMMARY_INPUT_NAMES.get(stem, ())
                    if stem == 'atlas' and isinstance(value.get('dir'), str):
                        recorded = value['dir'].rstrip('/') + '/atlas_repaired.csv'
                    if not recorded and not names:
                        raise ValueError(f'Unknown hashed input {summary}: {".".join((*context, key))}')
                elif 'inputs' in context and ('/' in key or ':' in key):
                    recorded = key  # older blood summaries use {path: sha256}
                else:
                    continue
                if recorded and (recorded.startswith('repo:') or Path(recorded).suffix in ('.py', '.sh')):
                    continue  # source belongs to the separate commit-based deposit
                yield {'field': '.'.join((*context, key)), 'path': recorded,
                       'names': names, 'sha256': item}
    yield from walk(data)


class SummaryInputCatalog:
    """Resolve archived paths and pathless hashes within supplied roots only."""
    def __init__(self, roots):
        self.roots = {label: Path(root) for label, root in roots.items()}
        self.by_name = None
        self.hashes = {}

    def digest(self, path):
        stat = path.stat()
        key = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)
        if key not in self.hashes:
            self.hashes[key] = sha256(path)
        return self.hashes[key]

    def resolve(self, item):
        recorded = item['path']
        if recorded:
            label, sep, relative = recorded.partition(':')
            candidates = []
            if sep and label in self.roots:
                candidates.append(self.roots[label] / relative)
            else:
                original = Path(recorded)
                for label, root in self.roots.items():
                    for base in (root, paths.ROOTS.get(label, root)):
                        try:
                            candidates.append(root / original.relative_to(base))
                        except ValueError:
                            pass
                    for name in (root.name, paths.ROOTS.get(label, root).name, BUNDLE_DIRS.get(label, root.name)):
                        if name in original.parts:
                            index = original.parts.index(name)
                            candidates.append(root.joinpath(*original.parts[index + 1:]))
            for candidate in candidates:
                if candidate.is_file() and self.digest(candidate) == item['sha256']:
                    return candidate
            if sep and label in self.roots:
                raise ValueError(f'Missing or mismatched hashed input {recorded} ({item["sha256"]})')
        # Absolute archive paths can name relocated inputs. Require the hash.
        names = (Path(recorded).name,) if recorded else item['names']
        if self.by_name is None:
            self.by_name = {}
            for root in self.roots.values():
                for path in sorted(root.rglob('*')):
                    if path.is_file():
                        self.by_name.setdefault(path.name, []).append(path)
        for name in names:
            for candidate in self.by_name.get(name, ()):
                if self.digest(candidate) == item['sha256']:
                    return candidate
        raise ValueError(f'Missing hashed input {recorded or names} ({item["sha256"]})')


def summary_input_closure(files, roots, excluded=lambda p: False, absent_ok=None):
    """Resolve summaries recursively; return inputs and provenance records.

    absent_ok decides which unresolved inputs are recorded instead of raised: True for all (staging from the author
    data, where an input of a reused stage can be absent, such as the metadata table behind a cached blood summary),
    or a set of (recorded path, sha256) pairs that the author data also lacks (the deposit check). Any other missing
    or altered input raises."""
    catalog = SummaryInputCatalog(roots)
    pending = [Path(p) for p in files if 'summary' in Path(p).name and str(p).endswith('.json')]
    seen, inputs, records = set(), set(), []
    while pending:
        summary = pending.pop()
        if summary in seen:
            continue
        seen.add(summary)
        for item in summary_inputs(summary):
            try:
                source = catalog.resolve(item)
            except (OSError, ValueError) as error:
                if not (absent_ok is True or (absent_ok and (item['path'], item['sha256']) in absent_ok)):
                    raise ValueError(f'{summary} [{item["field"]}]: {error}') from error
                # Absent from the author data as well; the sandboxed rebuild fails on any input it opens and lacks.
                records.append({'summary': str(summary), 'field': item['field'], 'recorded_path': item['path'],
                                'path': None, 'sha256': item['sha256'], 'unavailable': str(error)})
                continue
            if excluded(source):
                raise ValueError(f'{summary} [{item["field"]}]: input is a generated run output: {source}')
            inputs.add(source)
            records.append({'summary': str(summary), 'field': item['field'],
                            'recorded_path': item['path'], 'path': str(source), 'sha256': item['sha256']})
            if 'summary' in source.name and source.suffix == '.json':
                pending.append(source)
    return inputs, records


def stage_inputs(args):
    """Hard-link reused caches and the hash-matched historical prerequisites."""
    sys.path.insert(0, str(REPO / 'scripts'))
    import clockbase_verifier as cv
    source, target = paths.VERIFIER_ROOT, Path(args.out)
    seeds = [p for p in (source / args.reference).rglob('*summary*.json') if p.is_file()]
    plan = cv.stages(source, args.reference, paths.STREAM_EMPTY, paths.STREAM_DETACHED)
    seeds.extend(Path(s.done) for s in plan if s.done and Path(s.done).is_file())
    # Later standalone runs also reuse earlier final inputs, for example the
    # composite and questionnaire runs' recorded ranking from October 10.
    for summary in source.glob('*/summary.json'):
        if not summary.parent.name.startswith(('final-', 'prefinal-', 'selftest', 'pipeline-runs', 'pass-')):
            seeds.extend(summary.parent.glob('*summary*.json'))
    # The decision audit reads all saved decision runs, including pilots and
    # newer caches that did not exist when the reference audit was first made.
    for folder in source.iterdir():
        if folder.is_dir() and folder.name.startswith(('adjudication-', 'design-', 'classification-', 'individuals-')):
            seeds.extend(folder.glob('*summary*.json'))
    inputs, records = summary_input_closure(seeds, {label: paths.ROOTS[label] for label in BUNDLE_DIRS},
                                            absent_ok=True)
    try:
        target.mkdir(parents=True, exist_ok=False)
        for root, dirs, files in os.walk(source, followlinks=False):
            root = Path(root)
            if root == source:
                dirs[:] = [n for n in dirs if n == args.reference or
                           not n.startswith(('final-', 'prefinal-', 'selftest', 'pipeline-runs', 'pass-'))]
            dest = target / root.relative_to(source)
            dest.mkdir(parents=True, exist_ok=True)
            for name in dirs + files:
                if (root / name).is_symlink():
                    raise ValueError(f'verifier input symlink must be resolved before staging: {root / name}')
            for name in files:
                os.link(root / name, dest / name)
        added = []
        for path in sorted(inputs):
            if not path.is_relative_to(source):
                continue  # other input roots remain supplied by the environment
            relative = path.relative_to(source)
            if generated_verifier_input(relative, args.reference):
                raise ValueError(f'Historical input is under a generated output root: {path}')
            dest = target / relative
            if not dest.exists():
                dest.parent.mkdir(parents=True, exist_ok=True)
                os.link(path, dest)
                added.append(relative.as_posix())
        (target.parent / 'staged_summary_inputs.json').write_text(json.dumps({
            'reference': args.reference, 'historical_files_added': added, 'inputs': records,
        }, indent=2) + '\n')
    except (OSError, ValueError) as error:
        raise ValueError(f'Verifier input staging failed: {error}') from error
    print(f'Staged verifier caches and {len(added)} historical input files; checked {len(records)} input hashes')


def bundle(args):
    sys.path.insert(0, str(REPO / 'scripts'))
    import clockbase_verifier as cv  # the stage list defines which completion files mark reused stages
    out = Path(args.out) / 'clockbase_data'
    out.mkdir(parents=True, exist_ok=False)
    wanted, outside = set(), set()
    for log in Path(args.trace_dir).glob('*.txt'):
        for line in log.read_text().splitlines():
            if not line:
                continue
            label = paths.portable(line).partition(':')[0]
            (wanted if label in BUNDLE_DIRS else outside).add(line)
    plan = cv.stages(paths.VERIFIER_ROOT, args.reference, paths.STREAM_EMPTY, paths.STREAM_DETACHED)
    done_files = [Path(s.done) for s in plan if s.done and Path(s.done).exists() and Path(s.done).is_file()]
    wanted.update(map(str, done_files))
    reference = paths.VERIFIER_ROOT / args.reference
    wanted.update(str(p) for p in reference.rglob('*') if p.is_file())
    rederived = set()
    for root in (paths.STREAM_EMPTY, paths.STREAM_DETACHED):
        rederived.update(str(p) for p in (root / 'studies').glob('*/*') if p.name in REDERIVATION_FILES)
        rederived.update(str(root / name) for name in REDERIVATION_SUMMARIES if (root / name).is_file())
        for study in sorted((root / 'studies').iterdir()):
            metadata = paths.GEO / 'GEO_hs_dnam/processed_metadata' / f'{study.name}.csv'
            if metadata.is_file():
                rederived.add(str(metadata))
    if paths.METHYLOME_REFERENCE.is_file():  # the population-mean methylome of the detached-profile screen
        rederived.add(str(paths.METHYLOME_REFERENCE))
    rederivation_added = len(rederived - wanted)
    wanted.update(rederived)
    wanted.update(str(paths.resolve(name)) for name in FIGURE_INPUTS)
    trace_finals, selftest_roots = set(), set()
    for metadata in Path(args.trace_dir).glob('*trace_run.json'):
        run = json.loads(metadata.read_text())
        if run.get('selftest'):
            label, _, relative = run['selftest'].partition(':')
            if label != 'verifier':
                raise ValueError(f'Trace self-test must be under the verifier root: {metadata}')
            selftest_roots.add(relative)
        else:
            trace_finals.update((run['final'], 'pre' + run['final']))
    copied, headers, run_outputs = 0, 0, 0
    verifier_root = paths.VERIFIER_ROOT.resolve()

    def excluded(source):
        try:
            relative = Path(source).resolve().relative_to(verifier_root)
        except ValueError:
            return False
        return generated_verifier_input(relative, args.reference, trace_finals, selftest_roots)

    inputs, author_records = summary_input_closure(
        [p for p in wanted if not excluded(p)],
        {label: paths.ROOTS[label] for label in BUNDLE_DIRS}, excluded, absent_ok=True)
    absent_from_author_data = {(r['recorded_path'], r['sha256']) for r in author_records if r.get('unavailable')}
    wanted.update(map(str, inputs))
    for source in sorted(wanted):
        label, _, relative = paths.portable(source).partition(':')
        if label not in BUNDLE_DIRS:
            continue
        # Files that the traced runs themselves wrote (their run, self-test and log directories under the verifier
        # root) are outputs, not inputs: the staging step excludes the same directories, and a rebuild recreates them.
        if excluded(source):
            run_outputs += 1
            continue
        target = out / BUNDLE_DIRS[label] / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if label == 'geo' and SERIES_MATRIX.search(relative) and Path(source) not in inputs:
            copy_series_matrix_header(Path(source).resolve(), target)
            headers += 1
        else:
            shutil.copy2(Path(source).resolve(), target)
        copied += 1
    for name in AGENT_TRACES:  # the saved agent traces and per-study outputs (the input of the trace replay)
        source = paths.TRACES / name
        if source.is_dir():
            def keep_or_copy(src, dst):  # files already copied as traced inputs are kept when identical
                if Path(dst).exists():
                    if sha256_file(dst) != sha256_file(src):
                        raise SystemExit(f'agent trace differs from its traced copy: {dst}')
                    return dst
                return shutil.copy2(src, dst)
            shutil.copytree(source.resolve(), out / 'agent-traces' / name, symlinks=False, dirs_exist_ok=True,
                            copy_function=keep_or_copy)
        elif source.exists():
            (out / 'agent-traces').mkdir(exist_ok=True)
            target = out / 'agent-traces' / name
            if target.exists():  # already copied as a traced input (read-only copy): keep it if identical
                if sha256_file(target) != sha256_file(source.resolve()):
                    raise SystemExit(f'agent trace differs from its traced copy: {target}')
                continue
            shutil.copy2(source.resolve(), target)
    # Check the copied files using deposit roots only, with no author-data fallback.
    _, checked_inputs = summary_input_closure(
        out.rglob('*summary*.json'), {label: out / directory for label, directory in BUNDLE_DIRS.items()},
        absent_ok=absent_from_author_data)
    rows = write_manifest(out, 'BUNDLE_MANIFEST.tsv')
    report = {'reference_final': args.reference, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'files': len(rows),
              'bytes': sum(r[2] for r in rows), 'traced_files_copied': copied, 'stage_completion_files': len(done_files),
              'rederivation_files_added': rederivation_added, 'series_matrices_header_only': headers,
              'traced_run_outputs_skipped': run_outputs,
              'summary_input_hashes_checked': sum(1 for r in checked_inputs if not r.get('unavailable')),
              'summary_inputs_not_in_deposit': sorted({r['recorded_path'] or r['sha256'] for r in checked_inputs
                                                       if r.get('unavailable')}),
              'summary_input_files_added': len(inputs),
              'series_matrix_note': 'GEO series matrices are copied up to !series_matrix_table_begin (the sample '
                                    'annotations the verifier reads); the full files are at GEO under the same accession',
              'opened_outside_data_roots': sorted(p for p in outside if not re.search(r'site-packages|/lib/python|/usr/', p))[:200]}
    (Path(args.out) / 'BUNDLE.json').write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != 'opened_outside_data_roots'}, indent=1))
    print('opened outside the data roots (first 20):', report['opened_outside_data_roots'][:20])


def records_inventory(final=None):
    """Only the source files and saved outputs needed by the cited statements; no manuscript or review bundles."""
    evidence = REPO / 'paper/revision_2026-09-07/evidence'
    sources, groups = {}, {}
    benchmark = Path(os.environ.get('CLOCKBASE_BENCHMARK_DATA') or paths.DATA / 'benchmark')

    def add(group, root, names):
        deposited = []
        for name in names:
            target = f'supporting/{group}/{name}'
            sources[target] = root / name
            deposited.append(target)
        groups[group] = deposited

    add('census-sources', paths.GEO / 'result', [
        'qced_human_dnam_age.csv', 'human_rna_age.csv', 'mouse_rna_age.csv',
        'human_dnam_umap_age.csv', 'human_rna_umap_age.csv', 'mouse_rna_umap_age.csv'])
    add('historical-atlas', REPO, ['figures/combined_3modality_plot_data.csv'])
    add('evaluators', REPO, ['validation_results_100_traces.csv', 'validation_results_100_traces_evaluator2.csv',
                            'src/validate_100_traces.py', 'src/validate_100_traces_evaluator2.py'])
    add('figure3', evidence / 'figure3-score-components-2026-09-18', [
        'plan.json', 'execution-002.json', 'execution-002.log', 'saved-bar-fdr-readback.json',
        'run-002/component-summary.csv', 'run-002/display-collisions.json',
        'run-002/saved-bar-source-bindings.json', 'run-002/saved-output-agreement.json',
        'run-002/selection-boundaries.json', 'run-002/source-accounting.csv', 'run-002/summary.json',
        'run-002/table1_pro_aging_pooled_top30-readback.csv', 'run-002/table2_anti_aging_pooled_top30-readback.csv',
        'run-002/legacy_recipe-display_units.csv.gz', 'run-002/legacy_recipe-qualifying_rows.csv.gz',
        'run-002/legacy_recipe-selected_bars.csv.gz', 'run-002/uncapped_p_sensitivity-display_units.csv.gz',
        'run-002/uncapped_p_sensitivity-qualifying_rows.csv.gz', 'run-002/uncapped_p_sensitivity-selected_bars.csv.gz',
        'run-002/uniform_floor100_sensitivity-display_units.csv.gz',
        'run-002/uniform_floor100_sensitivity-qualifying_rows.csv.gz',
        'run-002/uniform_floor100_sensitivity-selected_bars.csv.gz',
        'independent-sol/verify_figure3_stdlib.py', 'independent-sol/figure3-source-receipt-v3.json'])
    add('figure3-script', evidence / 'comparison-identity-master-2026-09-18', [
        'reconstruct_comparison_identity_master_2026_09_18.py'])
    add('saved-ranking-tables', paths.PUBLISHED, ['top_hits_combined_modalities_v2.csv'])
    add('external-db', evidence / 'ed19-source-recovery-2026-09-18', [
        'external_db_crossvalidation_2026-05-14.csv'])
    add('gse179325', evidence / 'gse179325-age-validation-2026-09-18', [
        'input-manifest.json', 'run-004/descriptive-age-summaries.csv', 'run-004/descriptive-age-summaries.json',
        'run-004/joined-source-predictions.csv', 'run-004/receipt.json', 'validate_gse179325_age_2026_09_18.py',
        'independent-source/verify_age_calculations_2026_09_18.py',
        'independent-source/age-calculations-receipt-run004-v2.json',
        'independent-source/age-calculations-root-readback-v2.json'])
    add('predictor-dictionary', evidence / 'predictor-target-followup-2026-09-18', [
        'dictionary-v2/clock_target_unit_dictionary_v2.csv', 'dictionary-v2/clock_target_unit_dictionary_v2.json',
        'dictionary-v2/component-definition-counts.csv', 'dictionary-v2/definition-readback.json',
        'source_receipts.json'])
    add('sample-census', evidence / 'sample-atlas-census-2026-09-18', [
        'plan-v2.json', 'plan-for-replay.json', 'run-001/conditional-dnam-plot-selection.json',
        'run-001/coordinate-membership-differences.csv', 'run-001/duplicate-or-invalid-identities.csv',
        'run-001/numeric-and-identity-anomalies.csv', 'run-001/per-column-coverage.csv',
        'run-001/per-study-identity-counts.csv', 'run-001/source-and-command-receipt.json', 'run-001/summary.json',
        'independent/check_sample_atlas_census_stdlib.py', 'independent/independent-verification-003.json'])
    add('blood-reader', evidence / 'blood-source-audit-2026-09-18', [
        'reanalyze_blood_panel_2026_09_18.py'])
    add('fdr-withholding', evidence / 'testing-family-arithmetic-2026-09-18', [
        'withholding-002/corrections.json', 'withholding-002/application.json', 'withholding-002/summary.json',
        'withholding-002/tables/2026-09-08-DNAm-atlas_WORKING_NOT_VALIDATED.csv',
        'withholding-002/tables/2026-09-08-human-RNA-atlas_WORKING_NOT_VALIDATED.csv',
        'withholding-002/tables/2026-09-08-mouse-RNA-atlas_WORKING_NOT_VALIDATED.csv',
        'withholding-002/tables/2026-09-08-combined-plot_WORKING_NOT_VALIDATED.csv',
        'guard-v2-overlay-replay.json'])
    add('historical-analysis', REPO / 'src/historical_geo_pipeline/analysis', [
        '0_clean_up_data_hs_dnam.R', '1_dnam_analysis.R', '2_clean_up_data_rna.R', '3_rna_analysis.R',
        'summarize_plot.R', 'MANIFEST.tsv', 'README.md'])
    add('legacy-generators', REPO / 'scripts/legacy_sensitivity_analyses_2026_05_14', [
        'external_db_crossvalidation_2026-05-14.py', 'tissue_direction_enrichment_2026-05-14.py',
        'drug_target_class_enrichment_2026-05-14.py', 'README.md'])
    add('global-null', REPO, ['run_fpr_analysis.py'])
    add('questionaire', REPO, ['questionaire.md'])
    add('manuscript-agent', REPO / 'minipaper_figures_2026-07-01', [
        'pipeline.py', 'README.md', 'hits_score_gt100.csv'])
    # Benchmark of the final run (benchmark/published_clock_2026_10_05): the published results extracted from the
    # literature, the benchmark units and their classification, and the per-row records scored for this run. The
    # third-party sample annotations (Tyshkovskiy et al., Zenodo; ComputAgeBench) are not copied.
    add('benchmark', benchmark, [
        'ground_truth_t1_v2_final_2026-10-05.jsonl', 'groupbench/b-v2/linked.jsonl', 'groupbench/prepared-v2/units.jsonl',
        'groupbench/rules-v2/e1_rules.jsonl', 'groupbench/rules-v2/e2_rules.jsonl', 'groupbench/final-v2/outcomes.json',
        'groupbench/final-v2/final_report.json'] + [
        f'groupbench/design-checks/{final}/{name}' for name in (
            'dnam/rows.jsonl', 'hs_rna/rows.jsonl', 'mm_rna/rows.jsonl', 'a_result_stages.jsonl',
            'a_result_stages.summary.json', 'external.json', 'auroc.json', 'in_sample_corrected.json',
            'in_sample_with_readers_corrected.json', 'export_as_checks.json')] + [
        f'agentbench/{scored}/{name}' for scored in ('scored-001', 'scored-claude-001', 'scored-biomni-001')
        for name in ('agent_reruns.json', 'runs.csv', 'records.jsonl')] + [
        f'groupbench/published-released-001/{name}' for name in (
            'published_released_test.json', 'published_released_test.csv', 'individuals.jsonl', 'individuals_summary.json')])
    # Every agent rerun: its task, result files and Codex, Claude Code or Biomni event log. The statement index names
    # the directory; the Claude Code and Biomni runs (protocol items 22 and 24) are published beside the Codex runs.
    for run in sorted(list((benchmark / 'agentbench/runs').glob('codex-*/GSE*'))
                      + list((benchmark / 'agentbench/runs-claude').glob('claude-*/GSE*'))
                      + list((benchmark / 'agentbench/runs-biomni').glob('biomni-*/GSE*'))):
        for name in ('TASK.md', 'results/aging_analysis.csv', 'results/comparison_samples.csv', 'events.jsonl'):
            if (run / name).is_file():
                sources[f'supporting/benchmark/agentbench/runs/{run.parent.name}/{run.name}/{name}'] = run / name
    groups['benchmark-agent-runs'] = ['supporting/benchmark/agentbench/runs/']
    # Agent reruns of the ranked comparisons (protocol item 20): scores and every run's task, results and event log.
    add('benchmark-ranked', benchmark / 'agentbench/scored-ranked-001', ['ranked_reruns.json', 'ranked_reruns.csv', 'runs.csv'])
    add('benchmark-ranked-claude', benchmark / 'agentbench/scored-ranked-claude-001',
        ['ranked_reruns.json', 'ranked_reruns.csv', 'runs.csv'])
    add('benchmark-ranked-biomni', benchmark / 'agentbench/scored-ranked-biomni-001',
        ['ranked_reruns.json', 'ranked_reruns.csv', 'runs.csv'])
    for run in sorted(list((benchmark / 'agentbench/ranked-runs').glob('codex-*/GSE*'))
                      + list((benchmark / 'agentbench/ranked-runs-claude').glob('claude-*/GSE*'))
                      + list((benchmark / 'agentbench/ranked-runs-biomni').glob('biomni-*/GSE*'))):
        for name in ('TASK.md', 'results/aging_analysis.csv', 'results/comparison_samples.csv', 'events.jsonl'):
            if (run / name).is_file():
                sources[f'supporting/benchmark/agentbench/ranked-runs/{run.parent.name}/{run.name}/{name}'] = run / name
    groups['benchmark-ranked-runs'] = ['supporting/benchmark/agentbench/ranked-runs/']
    # Replication across independent projects (audit_system/intervention_meta_2026_09_29.py): the projects built
    # from public GEO series records, the primary analysis and the three sensitivity runs, and every GEO record read.
    meta = paths.VERIFIER_ROOT / 'intervention-meta-2026-10-17'
    add('intervention-meta', meta, ['projects.json', 'geo.log'] + [
        f'{run}/{name}' for run in ('primary', 'by-series', 'all-contexts', 'no-contact-link', 'keep-tumour')
        for name in ('intervention_meta.csv', 'unit_effects.csv', 'summary.json')] + [
        'concordance/summary.json', 'concordance/lifespan_concordance.csv'])
    # The design audit of the meta-analysis's own comparisons (2026-10-17): every verdict that excluded a comparison
    # as tumour material or cultured samples, with both models' answers.
    add('intervention-meta-design-audit', paths.VERIFIER_ROOT / 'design-audit-meta-001',
        ['candidates.csv', 'comparisons.jsonl', 'summary.json'])
    # The judgement version of the ranking (2026-10-18): both models' answers for every ranked comparison.
    add('ranking-judgement', paths.VERIFIER_ROOT / 'followup-judgment-001',
        ['comparisons.jsonl', 'summary.json', 'ranked_with_judgement.csv', 'judgement_age-lowering.csv',
         'judgement_age-raising.csv'])
    # The GenAge and DrugAge files themselves (HAGR) are not copied; concordance/summary.json records their SHA-256.
    for record in sorted((paths.VERIFIER_ROOT / 'geo-series-records-2026-10-14').glob('GSE*.txt')):
        sources[f'supporting/intervention-meta/geo-series-records/{record.name}'] = record
    groups['intervention-meta-geo'] = ['supporting/intervention-meta/geo-series-records/']
    statements = [
        ('res-prediction-tables', ['census-sources', 'sample-census'], 'Prediction records and the saved six-file census.'),
        ('res-fig3-bar-decomp', ['figure3', 'figure3-script', 'saved-ranking-tables', 'historical-atlas'],
         'Saved bar values, component decomposition, selection limits, script and execution receipts.'),
        ('res-257-table', ['external-db', 'legacy-generators'], 'Recovered external-database annotations and generator.'),
        ('res-evaluator-tables', ['evaluators'], 'The two saved study-level evaluator tables and their scripts.'),
        ('m-gse179325-summaries', ['gse179325'], 'All cohort-by-disease summaries, joined source ages and receipts.'),
        ('m-predictor-dictionary', ['predictor-dictionary'], 'Targets, scales, tissue context and implementation evidence.'),
        ('m-cleanup-source', ['historical-analysis'], 'Recovered GrimAge component selection and missing-score masking.'),
        ('m-strict-reader', ['blood-reader', 'sample-census', 'census-sources'],
         'Strict CSV reader, independent stdlib checker, saved census outputs and source tables.'),
        ('m-umap-script', ['historical-analysis'], 'Saved DNAm predictor filtering and coordinate-generation source.'),
        ('m-coordinate-table', ['census-sources', 'sample-census'], 'Recovered coordinate tables and membership counts.'),
        ('m-s11-generator', ['legacy-generators'], 'Historical tissue-direction enrichment generator.'),
        ('m-s12-generator', ['legacy-generators'], 'Historical drug target-class enrichment generator.'),
        ('m-flagged-public', ['historical-atlas'], 'Complete historical result table with outlier_status flags.'),
        ('m-eval-scripts', ['evaluators'], 'Seeded evaluator scripts and their saved outputs.'),
        ('m-eval-tables', ['evaluators'], 'Criterion-level labels and notes in both saved evaluator tables.'),
        ('m-global-null-script', ['global-null', 'historical-atlas'], 'Global-null script and its comparison-record input.'),
        ('m-fdr-manifest', ['fdr-withholding'], 'Exact-row before states, applied overlay tables and replay receipt.'),
        ('res-mouse-dnam-constant', ['historical-analysis'], 'Summary-source constant; no recovered mouse DNAm roster.'),
        ('m-questionaire', ['questionaire'], 'Historical plausibility filter and scoring rubric.'),
        ('m-ms-agent-gate', ['manuscript-agent'], 'Figure-validation source, coverage description and supplied candidate list.'),
        ('m-benchmark', ['benchmark', 'benchmark-agent-runs'], 'Published results extracted from the literature, benchmark units with their '
         'classification, the per-row benchmark records and scores of the final run, and the agent reruns with their '
         'tasks, results, event logs and scores (Supplementary Table S18).'),
        ('m-ranked-reruns', ['benchmark-ranked', 'benchmark-ranked-claude', 'benchmark-ranked-biomni',
                             'benchmark-ranked-runs'], 'Agent reruns of '
         'the 60 ranked comparisons of Tables 1 and 2 with Codex, Claude Code and Biomni: scores, run audit, and every run with '
         'its task, results and event log (Extended Data Figure 23f).'),
        ('m-ranking-judgement', ['ranking-judgement'], 'Judgement version of the ranking: both models\' answers to the '
         'three questions for every ranked comparison and the comparisons each direction keeps (Table 3; Supplementary '
         'Table S17).'),
        ('m-replication', ['intervention-meta', 'intervention-meta-geo', 'intervention-meta-design-audit'],
         'Projects built from the public GEO series records, the intervention meta-analysis by independent project with '
         'its sensitivity runs, the design audit that excluded tumour material and cultured samples, every GEO record '
         'read, and the concordance with GenAge and DrugAge (Extended Data Figure 22; Supplementary Tables S19 and S21).'),
    ]
    return sources, [{'id': identity, 'paths': [p for group in selected for p in groups[group]], 'shows': shows}
                     for identity, selected, shows in statements]


def scan_records_file(path):
    """Stream large source tables; inspect gzip-compressed text as well as ordinary UTF-8 files."""
    findings = {'private_path_lines': 0, 'credential_lines': 0}
    opener = gzip.open if path.suffix == '.gz' else open
    try:
        with opener(path, 'rt', encoding='utf-8', errors='strict') as handle:
            for number, line in enumerate(handle, 1):
                if CREDENTIAL.search(line) and not NOT_A_CREDENTIAL.search(line):
                    # Report the location only; never echo a possible credential into a build log.
                    raise ValueError(f'Credential-like string in {path}:{number}; records build stopped')
                if PRIVATE_PATH.search(line) and not DETECTION_PATTERN.search(line):
                    findings['private_path_lines'] += 1
    except UnicodeDecodeError:
        pass  # binary files do not have text lines
    return findings


def records(args):
    """Copy the final public records and the statement-specific evidence without changing source bytes."""
    destination = Path(args.out)
    if destination.exists():
        raise FileExistsError(f'Refusing existing records destination: {destination}')
    if not args.final or Path(args.final).name != args.final or args.final in {'.', '..'}:
        raise ValueError('--final must name one run directory under the verifier root')
    public = paths.VERIFIER_ROOT / args.final / 'public-export'
    if not public.is_dir():
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {public}')
    sources, statements = records_inventory(args.final)
    for source in sorted(p for p in public.rglob('*') if p.is_file()):
        sources['public-export/' + source.relative_to(public).as_posix()] = source
    if not any(name.startswith('public-export/') for name in sources):
        raise ValueError(f'Final public export is empty: {public}')
    statements.insert(0, {'id': 'per-comparison-verification-records', 'paths': ['public-export/'],
                          'shows': 'Unchanged public-export tree of the named final verifier run.'})
    for target, source in sources.items():
        if not source.is_file():
            raise FileNotFoundError(f'Cannot proceed - required data file not found: {source} ({target})')
    census_receipt = json.loads(sources['supporting/sample-census/run-001/source-and-command-receipt.json'].read_text())
    census_hashes = {Path(row['path']).name: row['sha256'] for row in census_receipt['sources']
                     if row.get('role') in {'predictions', 'coordinates'}}
    template = REPO / 'scripts/records_README_2026_10_05.md'
    text = template.read_text()
    if text.count('{{STATEMENTS}}') != 1:
        raise ValueError(f'Records README must contain one statement-table placeholder: {template}')
    out = destination / 'records'
    provenance, private_paths = {}, {}
    copied_bytes = 0
    try:
        destination.mkdir(parents=True, exist_ok=False)
        out.mkdir()
        for number, (relative, source) in enumerate(sorted(sources.items()), 1):
            target = out / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            before = source.stat()
            digest = sha256(source)
            if relative.startswith('supporting/census-sources/') and digest != census_hashes.get(source.name):
                raise ValueError(f'Census source hash disagrees with the saved receipt: {source}')
            shutil.copy2(source, target)
            if sha256(target) != digest or source.stat().st_mtime_ns != before.st_mtime_ns:
                raise ValueError(f'Source changed or copy hash disagrees: {source}')
            findings = scan_records_file(target)
            if findings['private_path_lines']:
                private_paths[relative] = findings['private_path_lines']
            provenance[relative] = {'source': str(source.resolve()), 'sha256': digest,
                                    'bytes': before.st_size, 'source_mtime_ns': before.st_mtime_ns}
            copied_bytes += before.st_size
            if number % 5000 == 0 or number == len(sources):
                print(f'Records: copied and checked {number}/{len(sources)} files', flush=True)
        table = '\n'.join('| ' + row['id'] + ' | ' + '<br>'.join(f'`{p}`' for p in row['paths']) + ' | ' +
                          row['shows'] + ' |' for row in statements)
        readme = out / 'RECORDS_README.md'
        readme.write_text(text.replace('{{FINAL}}', args.final).replace('{{STATEMENTS}}', table))
        scan_records_file(readme)
        provenance['RECORDS_README.md'] = {'source': str(template), 'generated': True}
        builder = str(Path(__file__).resolve())
        provenance['RECORDS_MANIFEST.tsv'] = {'source': builder, 'generated': True}
        provenance['RECORDS.json'] = {'source': builder, 'generated': True}
        report = {'final_run': args.final, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                  'files': len(sources) + 3, 'copied_files': len(sources), 'copied_bytes': copied_bytes,
                  'bytes': copied_bytes + readme.stat().st_size,
                  'bytes_scope': 'Copied files and README; RECORDS.json and RECORDS_MANIFEST.tsv are excluded.',
                  'credential_lines': 0, 'private_path_lines': sum(private_paths.values()),
                  'private_path_files': private_paths, 'statements': statements, 'sources': provenance}
        metadata = out / 'RECORDS.json'
        metadata.write_text(json.dumps(report, indent=2) + '\n')
        scan_records_file(metadata)
        write_manifest(out, 'RECORDS_MANIFEST.tsv')
        scan_records_file(out / 'RECORDS_MANIFEST.tsv')
    except (OSError, ValueError) as error:
        raise RuntimeError(f'Records build failed at {out}: {error}; partial output is not a completed deposit') from error
    print(json.dumps({k: v for k, v in report.items() if k not in {'sources', 'statements', 'private_path_files'}}, indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = parser.add_subparsers(dest='command', required=True)
    c = sub.add_parser('code')
    c.add_argument('--commit', default='HEAD')
    c.add_argument('--out', type=Path, required=True)
    c.add_argument('--worktree', action='store_true', help='copy tracked files from the working tree instead of a commit')
    c.add_argument('--lfs-source', type=Path, help='additional local code tree with hash-matched LFS assets; never downloads')
    t = sub.add_parser('trace')
    t.add_argument('--final', required=True)
    t.add_argument('--trace-dir', type=Path, required=True)
    t.add_argument('--selftest', type=Path, help='self-test root (new directory); trace the self-test instead')
    s = sub.add_parser('stage-inputs')
    s.add_argument('--reference', required=True)
    s.add_argument('--out', type=Path, required=True, help='new staged verifier directory')
    b = sub.add_parser('bundle')
    b.add_argument('--trace-dir', type=Path, required=True)
    b.add_argument('--reference', required=True, help='final run directory under the verifier root to ship for verify')
    b.add_argument('--out', type=Path, required=True)
    r = sub.add_parser('records')
    r.add_argument('--final', required=True, help='final run directory under the verifier root to deposit')
    r.add_argument('--out', type=Path, required=True, help='new destination; files are written under records/')
    args = parser.parse_args(argv)
    try:
        {'code': build_code, 'trace': trace, 'stage-inputs': stage_inputs,
         'bundle': bundle, 'records': records}[args.command](args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'Release {args.command} failed: {error}\n')


if __name__ == '__main__':
    main()
