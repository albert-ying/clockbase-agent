"""Prepare the agent reruns of the published-results series (PROTOCOL.md of groupbench, addendum item 16).

For every series and model a run directory is written with the three input files the historical task names, laid out
as in the verifier's trace replay, and TASK.md: a harness note, the historical task text of the series recovered from
its saved trace (paths pointing into the run directory, model label replaced) and one added deliverable, the samples
of each arm.

  python -m benchmark.published_clock_2026_10_05.agentbench.prepare --bench DATA --root RUNS --python PYTHON
"""
import argparse
import csv
import json
import pickle
import re
import shutil
from pathlib import Path

from audit_system import paths

TASK_CLOCKS = ('Horvathv1', 'PhenoAge', 'YingDamAge', 'DunedinPACE', 'GrimAgeV2')
MODELS = ('gpt-6-astra', 'gpt-6-sol', 'gpt-6-luna', 'gpt-5.5')
T1_RUNS = ('t1-openai-002-s4', 't1new-openai-001-s4')
PREFIX = re.compile(r'/dfs/user/hanchenw/Agent_PG/clockbase_agent/+')
HEADER = re.compile(r'^=+ Human Message =+\s*', re.S)
NOTE = """You are running as a command-line coding agent in the directory {run}.
Write Python scripts in this directory and run them with `{python}` (pandas, scipy and statsmodels are installed).
Read only the files under {run}/data. Do not read other locations and do not use the network.
Where the task below says not to execute external commands or create temporary files, it means: no R, no downloads
and no files outside this directory. Running your own Python scripts here is expected.

"""
EXTRA = """
Additional deliverable (required)
    - Save {run}/results/comparison_samples.csv with the columns comparison_name, condition_name, arm, gsm_id.
    - One row per sample of every comparison in the results table; arm is `case` (the condition) or `control`.
    - comparison_name and condition_name must match the results table exactly.
    - For a paired or covariate-adjusted analysis, list the samples that entered the test.
"""


def label(model, framework='codex'):
    return f'{framework.upper()}_' + re.sub(r'[^A-Z0-9]+', '_', model.upper()).strip('_')


def published_records(bench):
    """Published results on the task's clocks whose groups two models selected identically."""
    out = []
    for run in T1_RUNS:
        path = bench / 'runs' / run / 'records.jsonl'
        try:
            lines = path.read_text().splitlines()
        except OSError as error:
            raise SystemExit(f'Cannot proceed - required data file not found: {path} ({error})')
        for line in lines:
            record = json.loads(line)
            if (record['modality'] == 'human DNAm' and record['groups'].get('status') == 'AGREED'
                    and record['clock_clockbase'] in TASK_CLOCKS):
                out.append(record)
    return out


def task_text(trace, gse, run, python, model, framework='codex'):
    """The historical task of one series for a rerun."""
    first = trace['result'][0][0]
    text = HEADER.sub('', first if isinstance(first, str) else first.content).strip('\n')
    text = PREFIX.sub(f'{run}/', text)
    text = text.replace(f'{run}/results/hs_dnam_sonnet37/{gse}/', f'{run}/results/')
    text = text.replace('aging_analysis_sonnet37.csv', 'aging_analysis.csv').replace('SONNET37', label(model, framework))
    if '/dfs/' in text or 'sonnet' in text.lower():
        raise ValueError(f'{gse}: historical path or model label left in the task text')
    return NOTE.format(run=run, python=python) + text + '\n' + EXTRA.format(run=run)


def write_inputs(gse, run):
    overview = paths.GEO / 'gse_me_hs.csv'
    metadata = paths.GEO / 'GEO_hs_dnam/processed_metadata_preharmonize' / f'{gse}.csv'
    clocks = paths.VERIFIER_ROOT / 'replay-dnam-001/inputs' / f'{gse}.csv'
    for path in (overview, metadata, clocks):
        if not path.is_file():
            raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}')
    (run / 'data/dataset_level').mkdir(parents=True)
    (run / 'data/GEO_hs_dnam/processed_metadata_preharmonize').mkdir(parents=True)
    (run / 'data/clock_output').mkdir(parents=True)
    (run / 'results').mkdir()
    with overview.open(newline='') as source, (run / 'data/dataset_level/gse_me_hs.csv').open('w', newline='') as sink:
        reader = csv.DictReader(source)
        writer = csv.DictWriter(sink, fieldnames=reader.fieldnames)
        writer.writeheader()
        rows = [row for row in reader if row.get('accession') == gse]
        writer.writerows(rows)
    if not rows:
        raise ValueError(f'{gse}: no row in the dataset overview')
    shutil.copyfile(metadata, run / 'data/GEO_hs_dnam/processed_metadata_preharmonize' / f'{gse}.csv')
    shutil.copyfile(clocks, run / 'data/clock_output/human_dnam_age_short.csv')
    return {'overview': str(overview), 'metadata': str(metadata), 'clocks': str(clocks)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--bench', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True, help='new directory for the run directories')
    parser.add_argument('--python', required=True, help='interpreter named in the harness note')
    parser.add_argument('--series', nargs='*', help='only these series (pilot)')
    parser.add_argument('--models', nargs='*', default=list(MODELS))
    parser.add_argument('--framework', choices=['codex', 'claude', 'biomni'], default='codex',
                        help='agent framework; run directories are <framework>-<model>/<GSE> (PROTOCOL.md items 22, 24)')
    args = parser.parse_args()
    records = published_records(args.bench)
    with open(paths.TRACES / 'successful_traces_human_dnam.pkl', 'rb') as handle:
        traces = pickle.load(handle)
    published = sorted({r['gse'] for r in records})
    series = [g for g in published if g in traces]
    plan = {'published_results': len(records), 'published_series': len(published), 'series': series,
            'series_without_saved_trace': [g for g in published if g not in traces], 'models': args.models,
            'task_clocks': list(TASK_CLOCKS), 'python': args.python, 'framework': args.framework, 'runs': []}
    selected = [g for g in series if not args.series or g in args.series]
    args.root.mkdir(parents=True, exist_ok=True)
    for model in args.models:
        for gse in selected:
            run = (args.root / f'{args.framework}-{model}' / gse).resolve()
            if run.exists():
                continue
            try:
                sources = write_inputs(gse, run)
                (run / 'TASK.md').write_text(task_text(traces[gse], gse, run, args.python, model, args.framework))
            except (OSError, ValueError) as error:
                shutil.rmtree(run, ignore_errors=True)
                plan['runs'].append({'model': model, 'gse': gse, 'status': 'NOT_PREPARED', 'error': str(error)})
                continue
            plan['runs'].append({'model': model, 'gse': gse, 'status': 'PREPARED', 'run': str(run), 'sources': sources})
    name = 'plan.json' if not args.series else 'plan_pilot.json'
    (args.root / name).write_text(json.dumps(plan, indent=1) + '\n')
    print(len(series), 'series;', sum(r['status'] == 'PREPARED' for r in plan['runs']), 'runs prepared;',
          [r for r in plan['runs'] if r['status'] != 'PREPARED'][:5])


if __name__ == '__main__':
    main()
