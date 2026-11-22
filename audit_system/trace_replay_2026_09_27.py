"""Replay saved CodingAgent traces in a network-less sandbox to recover exactly what each analysis did.

Each saved trace holds the agent's executed Python blocks verbatim. The replay rewrites the
original cluster prefix (/dfs/user/hanchenw/Agent_PG/clockbase_agent/) to a per-trace sandbox
root whose data/ entries link to local, value-identical inputs, runs the blocks in order in one
namespace (as the agent did), and instruments scipy.stats and statsmodels so that every
statistical test records its input vectors, formula/data identifiers and P values. Output files
land in the sandbox and are compared with the saved historical CSV.

  python -m audit_system.trace_replay_2026_09_27 extract --modality dnam --out ROOT
  python -m audit_system.trace_replay_2026_09_27 split-inputs --modality dnam --out ROOT
  python -m audit_system.trace_replay_2026_09_27 run --modality dnam --out ROOT [--workers 24] [--limit N]

No model calls. Replays run under `unshare -rn` (no network) with a per-trace timeout.
GEO inputs and the replay interpreter are selected by audit_system.paths.
"""
import argparse
import hashlib
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from audit_system import paths

REPO = Path(__file__).resolve().parents[1]
GEO = paths.GEO
PREFIX = re.compile(r'/dfs/user/hanchenw/Agent_PG/clockbase_agent/+')
MODALITIES = {
    'dnam': {'traces': 'successful_traces_human_dnam.pkl', 'clock_file': 'human_dnam_age_short.csv',
             'local_clock': GEO / 'result/qced_human_dnam_age.csv', 'geo_folder': 'GEO_hs_dnam',
             'overview': 'gse_me_hs.csv', 'saved_results': ['hs_dnam_sonnet37']},
    'hs_rna': {'traces': 'successful_traces_human_rna.pkl', 'clock_file': 'human_rna_age_short.csv',
               'local_clock': GEO / 'result/human_rna_age.csv', 'geo_folder': 'GEO_hs_rna',
               'overview': 'gse_rna_hs.csv', 'saved_results': ['hs_rna_sonnet37', 'result_hs_rna']},
    'mm_rna': {'traces': 'successful_traces_mouse.pkl', 'clock_file': 'mouse_rna_age_short.csv',
               'local_clock': GEO / 'result/mouse_rna_age.csv', 'geo_folder': 'GEO_mm_rna',
               'overview': 'gse_rna_mm.csv', 'saved_results': ['result_mm_rna']},
}
TIMEOUT_SECONDS = 600
# CLOCKBASE_REPLAY_PYTHON selects the pinned historical replay environment when available.
REPLAY_PYTHON = paths.REPLAY_PYTHON
SCHEMA = 'clockbase-trace-replay-v1'

# Injected before the agent's first block. It must not change any numerical behaviour.
INSTRUMENT = r'''
import json as _json, numpy as _np, re as _re, sys as _sys
_LOG = open(_REPLAY_LOG, 'a')
_GSM = _re.compile(r'^GSM\d+$')

def _vec(x):
    try:
        return _np.asarray(x, dtype=float).ravel().tolist()
    except Exception:
        return None

def _ids(obj):
    try:
        import pandas as _pd
        if isinstance(obj, _pd.DataFrame):
            for c in obj.columns:
                vals = obj[c].astype(str)
                if len(vals) and vals.str.match(r'^GSM\d+$').all():
                    return vals.tolist()
            idx = obj.index.astype(str)
            if len(idx) and idx.str.match(r'^GSM\d+$').all():
                return idx.tolist()
        if isinstance(obj, _pd.Series):
            idx = obj.index.astype(str)
            if len(idx) and idx.str.match(r'^GSM\d+$').all():
                return idx.tolist()
    except Exception:
        pass
    return None

def _emit(record):
    try:
        _LOG.write(_json.dumps(record, default=str) + '\n'); _LOG.flush()
    except Exception:
        pass

import scipy.stats as _st
def _wrap(name):
    original = getattr(_st, name, None)
    if original is None:
        return
    def wrapped(*args, **kwargs):
        result = original(*args, **kwargs)
        try:
            p = float(getattr(result, 'pvalue', result[1]))
        except Exception:
            p = None
        try:
            stat = float(getattr(result, 'statistic', result[0]))
        except Exception:
            stat = None
        _emit({'fn': name, 'vectors': [_vec(a) for a in args[:2]], 'ids': [_ids(a) for a in args[:2]],
               'kwargs': {k: str(v)[:80] for k, v in kwargs.items()}, 'p': p, 'statistic': stat})
        return result
    setattr(_st, name, wrapped)
for _n in ('ttest_ind', 'ttest_rel', 'ttest_1samp', 'mannwhitneyu', 'wilcoxon', 'kruskal', 'f_oneway',
           'pearsonr', 'spearmanr', 'linregress', 'ranksums'):
    _wrap(_n)

def _patch_fit(cls, label):
    original = cls.fit
    def fit(self, *args, **kwargs):
        res = original(self, *args, **kwargs)
        try:
            names = list(getattr(self, 'exog_names', None) or [])
            exog = _np.asarray(self.exog, dtype=float) if getattr(self, 'exog', None) is not None else None
            rows = getattr(getattr(self, 'data', None), 'row_labels', None)
            rows = [str(r) for r in rows] if rows is not None else None
            _emit({'fn': label, 'formula': str(getattr(self, 'formula', '') or ''), 'endog': _vec(self.endog),
                   'exog_names': names, 'exog': exog.tolist() if exog is not None and exog.size <= 200000 else None,
                   'row_labels': rows, 'ids': [rows if rows and all(_GSM.match(r) for r in rows) else None],
                   'pvalues': dict(zip(names, map(float, _np.asarray(res.pvalues, dtype=float).ravel()))) if names else None,
                   'params': dict(zip(names, map(float, _np.asarray(res.params, dtype=float).ravel()))) if names else None,
                   'nobs': float(getattr(res, 'nobs', _np.nan))})
        except Exception as error:
            _emit({'fn': label, 'instrument_error': str(error)[:200]})
        return res
    cls.fit = fit

try:
    from statsmodels.regression.linear_model import OLS as _OLS, WLS as _WLS, GLS as _GLS
    from statsmodels.genmod.generalized_linear_model import GLM as _GLM
    from statsmodels.regression.mixed_linear_model import MixedLM as _MixedLM
    for _cls, _label in ((_OLS, 'sm.OLS'), (_WLS, 'sm.WLS'), (_GLS, 'sm.GLS'), (_GLM, 'sm.GLM'), (_MixedLM, 'sm.MixedLM')):
        _patch_fit(_cls, _label)
except Exception:
    pass
'''

RUNNER = r'''
import json, sys, traceback, io, contextlib, os, resource
# One replay may not take the host down: cap its address space (MemoryError is recorded as a block failure).
resource.setrlimit(resource.RLIMIT_AS, (12 << 30, 12 << 30))
spec = json.load(open(sys.argv[1]))
namespace = {'__name__': '__main__', '_REPLAY_LOG': spec['log']}
os.chdir(spec['cwd'])
exec(compile(spec['instrument'], '<instrument>', 'exec'), namespace)
# Some saved traces begin mid-conversation; define only the standard module names the agent assumed.
exec('import os, re, glob, warnings\nimport numpy as np\nimport pandas as pd\nfrom scipy import stats\n'
     'import statsmodels.api as sm\nimport statsmodels.formula.api as smf', namespace)
status = []
for i, block in enumerate(spec['blocks']):
    out = io.StringIO()
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
            exec(compile(block, f'<block {i}>', 'exec'), namespace)
        status.append({'block': i, 'ok': True})
    except BaseException as error:
        status.append({'block': i, 'ok': False, 'error': f'{type(error).__name__}: {str(error)[:200]}'})
json.dump(status, open(spec['status'], 'w'))
'''


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def python_blocks(messages):
    blocks = []
    for message in messages:
        if 'Ai Message' not in message[:120]:
            continue
        for block in re.findall(r'<execute>(.*?)</execute>', message, re.S):
            text = block.strip('\n')
            first = text.lstrip()
            if first.startswith(('#!R', '#!BASH', 'library(', '!')):
                continue
            blocks.append(text)
    return blocks


def extract(modality, root):
    """One JSON per trace: prompt GSE, prompt paths and the ordered Python blocks."""
    config = MODALITIES[modality]
    target = Path(root) / 'traces'
    target.mkdir(parents=True, exist_ok=True)
    source = paths.TRACES / config['traces']
    traces = pickle.load(open(source, 'rb'))
    index = []
    for key, entry in traces.items():
        messages = entry['result'][0]
        prompt = messages[0]
        match = re.search(r'from study (GSE\d+)', prompt)
        gse = match.group(1) if match else None
        record = {'trace_key': key, 'stored_study_id': entry.get('study_id'), 'prompt_gse': gse,
                  'prompt_paths': re.findall(r'`(/dfs[^`]+)`', prompt), 'blocks': python_blocks(messages),
                  'status': entry.get('status'), 'execution_time': entry.get('execution_time')}
        name = f'{key}.json'.replace('/', '_')
        (target / name).write_text(json.dumps(record))
        index.append({'file': name, 'trace_key': key, 'prompt_gse': gse, 'blocks': len(record['blocks'])})
    (Path(root) / 'trace_index.json').write_text(json.dumps({'source': str(source), 'source_sha256': sha256(source),
                                                               'traces': index}, indent=1))
    print(f'{len(index)} traces extracted')


def split_inputs(modality, root):
    """Per-study copies of the local clock table, so each replay reads only its own study."""
    import pandas as pd
    config = MODALITIES[modality]
    target = Path(root) / 'inputs'
    target.mkdir(parents=True, exist_ok=True)
    columns = None
    for chunk in pd.read_csv(config['local_clock'], chunksize=200_000, low_memory=False):
        columns = list(chunk.columns)
        key = 'gse_id' if 'gse_id' in chunk.columns else 'GSE_ID'
        for gse, part in chunk.groupby(key):
            path = target / f'{gse}.csv'
            part.to_csv(path, mode='a', header=not path.exists(), index=False)
    (Path(root) / 'inputs.json').write_text(json.dumps({'source': str(config['local_clock']),
                                                          'source_sha256': sha256(config['local_clock']),
                                                          'columns': columns}, indent=1))
    print(f'{len(list(target.iterdir()))} study inputs written')


def referenced_studies(record, limit=20):
    """Prompt study plus every GSE the executed code names (mouse traces mix studies)."""
    found = [record['prompt_gse']] if record.get('prompt_gse') else []
    for block in record['blocks']:
        for gse in re.findall(r'GSE\d+', block):
            if gse not in found:
                found.append(gse)
    return found[:limit]


def build_sandbox(modality, root, gse, workdir, studies=None):
    config = MODALITIES[modality]
    base = workdir / 'clockbase_agent'
    (base / 'data/clock_output').mkdir(parents=True)
    (base / 'data/dataset_level').mkdir(parents=True)
    (base / 'data' / config['geo_folder']).mkdir(parents=True)
    inputs = [Path(root) / 'inputs' / f'{g}.csv' for g in (studies or [gse])]
    inputs = [p for p in inputs if p.exists()]
    target = base / 'data/clock_output' / config['clock_file']
    if len(inputs) == 1:
        os.symlink(inputs[0], target)
    elif inputs:
        with open(target, 'w') as sink:
            for i, path in enumerate(inputs):
                with open(path) as source:
                    header = source.readline()
                    if i == 0:
                        sink.write(header)
                    shutil.copyfileobj(source, sink)
    os.symlink(GEO / config['overview'], base / 'data/dataset_level' / config['overview'])
    os.symlink(GEO / config['geo_folder'] / 'processed_metadata_preharmonize',
               base / 'data' / config['geo_folder'] / 'processed_metadata_preharmonize')
    for writable in ('results', 'result'):
        (base / writable).mkdir()
    return base


def replay_one(modality, root, trace_file):
    root = Path(root)
    record = json.loads((root / 'traces' / trace_file).read_text())
    gse = record['prompt_gse']
    out_dir = root / 'replays' / trace_file.removesuffix('.json')
    if (out_dir / 'result.json').exists():
        return json.loads((out_dir / 'result.json').read_text())
    out_dir.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as tmp:
        workdir = Path(tmp)
        studies = referenced_studies(record)
        base = build_sandbox(modality, root, gse, workdir, studies)
        blocks = [PREFIX.sub(str(base) + '/', b) for b in record['blocks']]
        spec = {'blocks': blocks, 'instrument': INSTRUMENT, 'log': str(out_dir / 'tests.jsonl'),
                'status': str(out_dir / 'blocks.json'), 'cwd': str(workdir)}
        (workdir / 'spec.json').write_text(json.dumps(spec))
        (workdir / 'runner.py').write_text(RUNNER)
        env = {'PATH': os.environ['PATH'], 'HOME': str(workdir), 'MPLBACKEND': 'Agg', 'OMP_NUM_THREADS': '1',
               'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1', 'PYTHONHASHSEED': '0'}
        try:
            proc = subprocess.run(['unshare', '-rn', REPLAY_PYTHON, '-B', str(workdir / 'runner.py'),
                                   str(workdir / 'spec.json')], capture_output=True, text=True,
                                  timeout=TIMEOUT_SECONDS, env=env, cwd=workdir)
            outcome, stderr = ('COMPLETED' if proc.returncode == 0 else 'RUNNER_FAILED'), proc.stderr[-2000:]
        except subprocess.TimeoutExpired:
            outcome, stderr = 'TIMEOUT', ''
        outputs = []
        for path in sorted((base / 'results').rglob('*.csv')) + sorted((base / 'result').rglob('*.csv')):
            target = out_dir / 'outputs' / path.relative_to(base)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)
            outputs.append(str(target.relative_to(out_dir)))
    blocks_status = json.loads((out_dir / 'blocks.json').read_text()) if (out_dir / 'blocks.json').exists() else []
    tests = sum(1 for _ in open(out_dir / 'tests.jsonl')) if (out_dir / 'tests.jsonl').exists() else 0
    first_error = next((b.get('error', '') for b in blocks_status if not b['ok']), '')
    incomplete = bool(re.match(r"NameError: name '(?!pd'|np'|os'|stats'|sm'|smf')", first_error)) and bool(blocks_status) and not blocks_status[0]['ok']
    result = {'schema': SCHEMA, 'trace_file': trace_file, 'prompt_gse': gse, 'outcome': outcome,
              'studies_in_sandbox': studies, 'trace_starts_mid_analysis': incomplete,
              'input_available': (root / 'inputs' / f'{gse}.csv').exists(),
              'blocks': len(blocks), 'blocks_failed': sum(not b['ok'] for b in blocks_status),
              'tests_logged': tests, 'outputs': outputs, 'seconds': round(time.time() - started, 1),
              'stderr_tail': stderr}
    (out_dir / 'result.json').write_text(json.dumps(result, indent=1))
    return result


def saved_result_tables(modality, gse):
    frames = []
    for folder in MODALITIES[modality]['saved_results']:
        for path in sorted((paths.TRACES / folder / gse).glob('aging_analysis*.csv')):
            frames.append(path)
    return frames


def _key_rows(frame):
    import pandas as pd
    frame = frame.rename(columns={'pvalue': 'p_value'})
    need = {'comparison_name', 'clock', 'p_value'}
    if not need <= set(frame.columns):
        return None
    frame = frame.copy()
    frame['p_value'] = pd.to_numeric(frame['p_value'], errors='coerce')
    return frame


def _value_index(study_input):
    """Map each exact clock value to its GSM, per clock column."""
    import pandas as pd
    table = pd.read_csv(study_input, low_memory=False)
    gsm_column = 'gsm_id' if 'gsm_id' in table.columns else next(c for c in table.columns if c.lower() in ('gsm', 'gsm_id', 'sample_id'))
    index = {}
    for column in table.columns:
        values = pd.to_numeric(table[column], errors='coerce')
        if values.notna().sum() == 0:
            continue
        mapping = {}
        for value, gsm in zip(values, table[gsm_column]):
            if pd.notna(value):
                mapping.setdefault(float(value), set()).add(str(gsm))
        index[column] = mapping
    return index


def _members(vector, index):
    """GSMs whose clock values equal the vector's values, for the unique clock column that explains all of them."""
    if not vector:
        return None, None
    hits = []
    for column, mapping in index.items():
        found = [mapping.get(float(v)) for v in vector]
        if all(f is not None and len(f) == 1 for f in found):
            hits.append((column, [next(iter(f)) for f in found]))
    if len(hits) == 1:
        return hits[0]
    columns = {c for c, _ in hits}
    members = {tuple(sorted(m)) for _, m in hits}
    return (sorted(columns)[0] if len(members) == 1 and hits else None), (hits[0][1] if len(members) == 1 and hits else None)


def analyze_one(modality, root, trace_file):
    """Compare replayed outputs with the saved CSV and recover GSM membership of every logged test."""
    import pandas as pd
    root = Path(root)
    replay = root / 'replays' / trace_file.removesuffix('.json')
    result = json.loads((replay / 'result.json').read_text())
    gse = result['prompt_gse']
    summary = {'trace_file': trace_file, 'prompt_gse': gse, 'outcome': result['outcome']}
    def read(path):
        try:
            return pd.read_csv(path)
        except (pd.errors.EmptyDataError, pd.errors.ParserError, UnicodeDecodeError):
            return None
    saved = [f for f in map(read, saved_result_tables(modality, gse)) if f is not None]
    replayed = [f for f in (read(replay / o) for o in result['outputs'] if o.endswith('.csv')) if f is not None]
    saved_rows = pd.concat([f for f in map(_key_rows, saved) if f is not None], ignore_index=True) if saved else None
    replay_rows = pd.concat([f for f in map(_key_rows, replayed) if f is not None], ignore_index=True) if replayed and any(_key_rows(f) is not None for f in replayed) else None
    if saved_rows is not None and replay_rows is not None and len(saved_rows):
        s = saved_rows.dropna(subset=['p_value'])
        r = replay_rows.dropna(subset=['p_value'])
        matched = s.apply(lambda row: bool(((r.comparison_name == row.comparison_name) & (r.clock == row.clock)
                                            & (abs(r.p_value - row.p_value) <= 1e-9 * max(abs(row.p_value), 1e-300))).any()), axis=1)
        summary.update(saved_rows=int(len(s)), saved_rows_reproduced=int(matched.sum()))
    else:
        summary.update(saved_rows=0 if saved_rows is None else int(len(saved_rows)), saved_rows_reproduced=0)
    memberships = []
    study_inputs = [root / 'inputs' / f'{g}.csv' for g in (result.get('studies_in_sandbox') or [gse])]
    study_inputs = [p for p in study_inputs if p.exists()]
    tests_path = replay / 'tests.jsonl'
    if study_inputs and tests_path.exists():
        index = {}
        for study_input in study_inputs:
            for column, mapping in _value_index(study_input).items():
                merged = index.setdefault(column, {})
                for value, gsms in mapping.items():
                    merged.setdefault(value, set()).update(gsms)
        for line in open(tests_path):
            test = json.loads(line)
            entry = {'fn': test['fn'], 'p': test.get('p'), 'kwargs': test.get('kwargs')}
            if test['fn'].startswith('sm.'):
                entry.update(formula=test.get('formula'), pvalues=test.get('pvalues'), params=test.get('params'),
                             exog_names=test.get('exog_names'), exog=test.get('exog'),
                             members=(test.get('ids') or [None])[0])
                if entry['members'] is None and test.get('endog'):
                    entry['clock'], entry['members'] = _members(test['endog'], index)
            else:
                vectors = test.get('vectors') or []
                groups = []
                for i, vector in enumerate(vectors):
                    ids = (test.get('ids') or [None, None])[i] if i < len(test.get('ids') or []) else None
                    clock, gsms = _members(vector, index) if ids is None else (None, ids)
                    groups.append({'clock': clock, 'gsms': gsms, 'n': len(vector) if vector else 0})
                entry['groups'] = groups
            memberships.append(entry)
    summary['tests'] = len(memberships)
    summary['tests_with_members'] = sum(1 for m in memberships
                                        if (m.get('groups') and all(g['gsms'] for g in m['groups']))
                                        or m.get('members'))
    (replay / 'memberships.jsonl').write_text(''.join(json.dumps(m, default=str) + '\n' for m in memberships))
    (replay / 'analysis.json').write_text(json.dumps(summary, indent=1))
    return summary


def _analyze_safe(modality, root, trace_file):
    try:
        return analyze_one(modality, root, trace_file)
    except Exception as error:
        return {'trace_file': trace_file, 'analysis_error': f'{type(error).__name__}: {error}'[:300]}


def analyze(modality, root, workers=24):
    from concurrent.futures import ProcessPoolExecutor
    root = Path(root)
    files = [p.name + '.json' for p in sorted((root / 'replays').iterdir()) if (p / 'result.json').exists()]
    rows = []
    with ProcessPoolExecutor(workers) as pool:
        for summary in pool.map(_analyze_safe, [modality] * len(files), [root] * len(files), files, chunksize=8):
            rows.append(summary)
    total = {k: sum(r.get(k, 0) for r in rows) for k in ('saved_rows', 'saved_rows_reproduced', 'tests', 'tests_with_members')}
    total['traces'] = len(rows)
    total['analysis_errors'] = sum(1 for r in rows if 'analysis_error' in r)
    (root / 'analysis_summary.json').write_text(json.dumps({'totals': total, 'traces': rows}, indent=1))
    print(json.dumps(total))


def run(modality, root, workers=24, limit=None, only=None):
    root = Path(root)
    index = json.loads((root / 'trace_index.json').read_text())['traces']
    if only:
        wanted = set(only.split(','))
        index = [t for t in index if t['prompt_gse'] in wanted]
    if limit:
        index = index[:limit]
    counts, started = {}, time.time()
    with ThreadPoolExecutor(workers) as pool:
        for i, result in enumerate(pool.map(lambda t: replay_one(modality, root, t['file']), index)):
            counts[result['outcome']] = counts.get(result['outcome'], 0) + 1
            if i % 200 == 0:
                print(f'{i}/{len(index)} {counts} {time.time() - started:.0f}s', flush=True)
    print(json.dumps(counts))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    commands = parser.add_subparsers(dest='command', required=True)
    for name in ('extract', 'split-inputs', 'run', 'analyze'):
        command = commands.add_parser(name)
        command.add_argument('--modality', choices=sorted(MODALITIES), required=True)
        command.add_argument('--out', type=Path, required=True)
        if name in ('run', 'analyze'):
            command.add_argument('--workers', type=int, default=24)
        if name == 'run':
            command.add_argument('--limit', type=int)
            command.add_argument('--only')
    args = parser.parse_args(argv)
    if args.command == 'extract':
        extract(args.modality, args.out)
    elif args.command == 'split-inputs':
        split_inputs(args.modality, args.out)
    elif args.command == 'analyze':
        analyze(args.modality, args.out, args.workers)
    else:
        run(args.modality, args.out, args.workers, args.limit, args.only)


if __name__ == '__main__':
    main()
