"""Score the agent reruns against the published results (PROTOCOL.md of groupbench, addendum item 16).

Each run directory holds the agent's results table and the samples of each arm. A comparison of the run matches a
published result when it is on the same clock and each arm has a Jaccard index of at least 0.5 with the corresponding
published group, in either orientation. Two results enter the published-clock records as systems of their own:

  R:<model>:raw          the direction and P value the agent reported, oriented to the published contrast
  R:<model>:recomputed   the verifier's statistic on the run's arms (design-respecting estimand, the verifier's clock
                         values, overflow rule, at least three samples, individuals or pairs per arm)

The historical agent's row (S1), the released row (S2c, from the checks records of the scored verifier run), ClockBase
with the paper's groups (S3) and the no-data models (S4) are scored on the same records with the same functions as the
published-clock summary. Every run is audited from its Codex or Claude Code event log: commands, web searches, spawned
agents, paths named outside the run directory, token usage.

  python -m benchmark.published_clock_2026_10_05.agentbench.score --bench DATA --root RUNS --checks CHECKS_DIR --out DIR
"""
import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import design_estimand_2026_10_10 as design
from audit_system import paths
from audit_system.row_repairs_2026_09_28 import MIN_ARM, recompute

from ..engine.metrics import paper_key
from ..engine.predictions import PredictionReader, atlas_rows
from ..engine.stats import threshold
from ..groupbench.design_eval import finding_metrics, signed_z
from ..summarize_2026_10_05 import judged, shared, t1_metrics, with_interval
from .prepare import TASK_CLOCKS, published_records

JACCARD = 0.5
KEYS = ('direction_agreement', 'call_agreement', 'call_balanced_accuracy')
CLOCK = {c.lower(): c for c in TASK_CLOCKS}
ALLOWED = ('/usr/', '/bin/', '/tmp/', '/dev/', '/proc/', '/etc/')
ABSOLUTE = re.compile(r'(?<![\w.~-])/(?:[\w.+-]+/)+[\w.+-]*')
NETWORK = re.compile(r'https?://|\b(curl|wget|git clone|pip3? install)\b')


def number(value):
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def read_run(run):
    """The comparisons of one run: arms by (comparison, condition) and one row per comparison and clock."""
    table, samples = run / 'results/aging_analysis.csv', run / 'results/comparison_samples.csv'
    out = {'rows': [], 'arms': {}, 'problems': []}
    if not table.is_file() or not samples.is_file():
        out['problems'].append('NO_RESULTS_TABLE' if not table.is_file() else 'NO_SAMPLES_FILE')
        return out
    try:
        with samples.open(newline='', encoding='utf-8-sig') as handle:
            for row in csv.DictReader(handle):
                key = ((row.get('comparison_name') or '').strip(), (row.get('condition_name') or '').strip())
                arm, gsm = (row.get('arm') or '').strip().lower(), (row.get('gsm_id') or '').strip()
                if arm in ('case', 'control') and gsm:
                    out['arms'].setdefault(key, {'case': [], 'control': []})[arm].append(gsm)
        with table.open(newline='', encoding='utf-8-sig') as handle:
            for row in csv.DictReader(handle):
                key = ((row.get('comparison_name') or '').strip(), (row.get('condition_name') or '').strip())
                clock = CLOCK.get((row.get('clock') or '').strip().lower())
                if clock is None:
                    out['problems'].append('OTHER_CLOCK')
                    continue
                if key not in out['arms']:
                    out['problems'].append('ROW_WITHOUT_SAMPLES')
                    continue
                out['rows'].append({'comparison': key[0], 'condition': key[1], 'clock': clock, 'p': number(row.get('p_value')),
                                    'effect': number(row.get('effect_size')), 'direction': number(row.get('direction')),
                                    'test': (row.get('statistical_test') or '').strip()})
    except (OSError, ValueError, csv.Error) as error:
        out['problems'].append(f'UNREADABLE: {error}')
    return out


def audit(run, python):
    """What the Codex event log shows of one run."""
    out = {'status': None, 'commands': 0, 'web_searches': 0, 'spawned_agents': 0, 'outside_paths': [], 'usage': None,
           'seconds': None, 'network_attempts': 0}
    try:
        out['status'] = (run / 'status').read_text().strip()
        started, finished = (pd.Timestamp((run / n).read_text().strip()) for n in ('started', 'finished'))
        out['seconds'] = float((finished - started).total_seconds())
    except (OSError, ValueError):
        pass
    allowed = (str(run), str(Path(python).parent.parent)) + ALLOWED
    try:
        lines = (run / 'events.jsonl').read_text().splitlines()
    except OSError:
        return out
    for line in lines:
        try:
            event = json.loads(line)
        except ValueError:
            continue
        if event.get('usage'):
            out['usage'] = event['usage']
        item = event.get('item') or {}
        if event.get('type') != 'item.completed':
            continue
        kind = str(item.get('type'))
        if kind == 'command_execution':
            out['commands'] += 1
            command = str(item.get('command', ''))
            out['network_attempts'] += bool(NETWORK.search(command))  # the sandbox gives commands no network
            for path in ABSOLUTE.findall(command):
                # A string such as "/clock_output/x.csv" joined to a variable is not a path: the first component must
                # be a directory at the filesystem root.
                if (not path.startswith(allowed) and Path('/' + path.split('/')[1]).is_dir()
                        and path not in out['outside_paths']):
                    out['outside_paths'].append(path)
        elif 'web_search' in kind:
            out['web_searches'] += 1
        elif 'collab' in kind or 'agent' in kind and kind != 'agent_message':
            out['spawned_agents'] += 1
    return out


FILE_TOOLS = {'Read': 'file_path', 'Write': 'file_path', 'Edit': 'file_path', 'NotebookEdit': 'notebook_path',
              'Glob': 'path', 'Grep': 'path'}
WEB_TOOLS = {'WebSearch', 'WebFetch'}
AGENT_TOOLS = {'Agent', 'Task'}


def run_dir(root, framework, model, gse):
    """Run directory of one series and model: <framework>-<model>/<GSE> (PROTOCOL.md items 16, 20 and 22)."""
    return Path(root) / f'{framework}-{model}' / gse


def audit_claude(run, python):
    """What the Claude Code stream-json log shows of one run (PROTOCOL.md item 22).

    The rule of the Codex audit, extended to Claude Code's file tools: a path outside the run directory, the Python
    environment and system locations counts whether it is named in a shell command or passed to Read, Write, Edit,
    NotebookEdit, Glob or Grep. Usage is taken from the final result event.
    """
    out = {'status': None, 'commands': 0, 'web_searches': 0, 'spawned_agents': 0, 'outside_paths': [], 'usage': None,
           'seconds': None, 'network_attempts': 0}
    try:
        out['status'] = (run / 'status').read_text().strip()
        started, finished = (pd.Timestamp((run / n).read_text().strip()) for n in ('started', 'finished'))
        out['seconds'] = float((finished - started).total_seconds())
    except (OSError, ValueError):
        pass
    allowed = (str(run), str(Path(python).parent.parent)) + ALLOWED

    def note(path):
        path = str(path)
        if not path.startswith('/'):
            return                                     # relative paths resolve inside the run directory
        if (not path.startswith(allowed) and Path('/' + path.split('/')[1]).is_dir()
                and path not in out['outside_paths']):
            out['outside_paths'].append(path)

    try:
        lines = (run / 'events.jsonl').read_text().splitlines()
    except OSError:
        return out
    for line in lines:
        try:
            event = json.loads(line)
        except ValueError:
            continue
        if event.get('type') == 'result' and event.get('usage'):
            usage = event['usage']
            out['usage'] = {'input_tokens': usage.get('input_tokens', 0) + usage.get('cache_creation_input_tokens', 0)
                            + usage.get('cache_read_input_tokens', 0),
                            'cached_input_tokens': usage.get('cache_read_input_tokens', 0),
                            'output_tokens': usage.get('output_tokens', 0),
                            'cost_usd': event.get('total_cost_usd')}
        if event.get('type') != 'assistant':
            continue
        for block in (event.get('message') or {}).get('content') or []:
            if not isinstance(block, dict) or block.get('type') != 'tool_use':
                continue
            name, args = block.get('name'), block.get('input') or {}
            if name == 'Bash':
                out['commands'] += 1
                command = str(args.get('command', ''))
                out['network_attempts'] += bool(NETWORK.search(command))
                for path in ABSOLUTE.findall(command):
                    note(path)
            elif name in FILE_TOOLS:
                if args.get(FILE_TOOLS[name]):
                    note(args[FILE_TOOLS[name]])
            elif name in WEB_TOOLS:
                out['web_searches'] += 1
            elif name in AGENT_TOOLS:
                out['spawned_agents'] += 1
    return out


EXECUTE = re.compile(r'<execute>(.*?)</execute>', re.S | re.I)
SOLUTION = re.compile(r'<solution>', re.I)
BIOMNI_TOOL = re.compile(r'\bbiomni\.tool\.(?!support_tools\b)\w+|\bfrom\s+biomni\.tool\s+import\b')


def executed_block(content):
    """The block Biomni A1 executes for one model message, or None (biomni/agent/a1.py, generate and execute): the
    first <execute> block, closed if unclosed, unless the message holds a <solution>, which ends the run."""
    text = str(content)
    if text.lstrip().startswith('<observation>') or SOLUTION.search(text) or '<execute>' not in text.lower():
        return None
    if '</execute>' not in text.lower():
        text += '</execute>'
    found = EXECUTE.search(text)
    return found.group(1) if found else None


def audit_biomni(run, python):
    """What the Biomni A1 events of one run show (PROTOCOL.md item 24).

    The commands are the executed blocks (Python, R or Bash). The path rule of items 16 and 22 applies to their text.
    A block that matches the network pattern or imports or calls a Biomni tool module other than support_tools counts
    as web access, which excludes the run. Usage is the sum of the recorded model calls.
    """
    out = {'status': None, 'commands': 0, 'web_searches': 0, 'spawned_agents': 0, 'outside_paths': [], 'usage': None,
           'seconds': None, 'network_attempts': 0, 'tool_module_calls': 0}
    try:
        out['status'] = (run / 'status').read_text().strip()
        started, finished = (pd.Timestamp((run / n).read_text().strip()) for n in ('started', 'finished'))
        out['seconds'] = float((finished - started).total_seconds())
    except (OSError, ValueError):
        pass
    allowed = (str(run), str(Path(python).parent.parent)) + ALLOWED
    try:
        lines = (run / 'events.jsonl').read_text().splitlines()
    except OSError:
        return out
    usage = {'input_tokens': 0, 'cached_input_tokens': 0, 'output_tokens': 0}
    for line in lines:
        try:
            event = json.loads(line)
        except ValueError:
            continue
        if event.get('type') == 'usage':
            for key in usage:
                usage[key] += int(event.get(key) or 0)
            out['usage'] = usage
        if event.get('type') != 'message' or event.get('role') != 'ai':
            continue
        block = executed_block(event.get('content'))
        if block is None:
            continue
        out['commands'] += 1
        network, tool = bool(NETWORK.search(block)), bool(BIOMNI_TOOL.search(block))
        out['network_attempts'] += network
        out['tool_module_calls'] += tool
        out['web_searches'] += network or tool
        for path in ABSOLUTE.findall(block):
            if (not path.startswith(allowed) and Path('/' + path.split('/')[1]).is_dir()
                    and path not in out['outside_paths']):
                out['outside_paths'].append(path)
    return out


def audit_run(framework, run, python):
    if framework == 'biomni':
        return audit_biomni(run, python)
    return audit_claude(run, python) if framework == 'claude' else audit(run, python)


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 0.0


def match(record, run):
    """The run's comparison for one published result: (row, swapped, smaller Jaccard index), or None."""
    groups, best = record['groups'], None
    for row in run['rows']:
        if row['clock'] != record['clock_clockbase']:
            continue
        arms = run['arms'][(row['comparison'], row['condition'])]
        for swapped in (False, True):
            case, control = (arms['control'], arms['case']) if swapped else (arms['case'], arms['control'])
            a, b = jaccard(groups['case_gsms'], case), jaccard(groups['control_gsms'], control)
            rank = (min(a, b), a + b, not swapped)
            if min(a, b) >= JACCARD and (best is None or rank > best[0]):
                best = (rank, row, swapped)
    return None if best is None else (best[1], best[2], best[0][0])


def raw_result(row, swapped, overlap):
    """The agent's reported direction and P value in the published orientation."""
    sign = row['direction'] if row['direction'] in (1.0, -1.0) else row['effect']
    if sign is None or sign == 0:
        return {'covered': False, 'status': 'NO_REPORTED_DIRECTION'}
    higher = (sign > 0) != swapped
    return {'covered': True, 'status': 'REPORTED', 'direction': 'higher' if higher else 'lower', 'p': row['p'],
            'swapped': swapped, 'jaccard_min': overlap, 'comparison': row['comparison'], 'test': row['test']}


def recomputed_result(record, row, arms, swapped, overlap, reader, series_cache):
    """The verifier's statistic on the run's arms, oriented to the published contrast."""
    gse, clock = record['gse'], record['clock_clockbase']
    prediction = reader.load(record['modality'], gse, clock, atlas_rows(record['modality'], gse))
    if prediction['status'] != 'OK':
        return {'covered': False, 'status': prediction['status']}
    values = prediction['samples']['clock']
    case, control = (arms['control'], arms['case']) if swapped else (arms['case'], arms['control'])
    dropped = set(design.overflow(values, case + control))
    case, control = [g for g in case if g not in dropped], [g for g in control if g not in dropped]
    if gse not in series_cache:
        series_cache[gse] = design.load_series(paths.GEO, 'GEO_hs_dnam', gse, True)
    estimate = design.estimate(values, case, control, series_cache[gse], 'dnam', clock, (row['comparison'], row['condition']), ())
    if estimate is not None and estimate['recomputed'] is not None:
        result, method = estimate['recomputed'], estimate['design']['method']
    else:
        result, method = recompute(values, case, control), 'WELCH'
    diff = result.get('mean_diff')
    if diff is None:
        return {'covered': False, 'status': 'EMPTY_ARM'}
    tested = min(result['n_case'], result['n_control']) >= MIN_ARM and result.get('p') is not None
    return {'covered': True, 'status': 'TESTED' if tested else 'DESCRIPTIVE', 'method': method,
            'direction': 'higher' if diff > 0 else 'lower' if diff < 0 else 'none', 'p': result.get('p') if tested else None,
            'mean_diff': diff, 'n_case': result['n_case'], 'n_control': result['n_control'], 'swapped': swapped,
            'jaccard_min': overlap, 'overflow_samples_dropped': len(dropped), 'comparison': row['comparison']}


def released_result(record, checks):
    """The released row of the scored verifier run for the historical agent's matched row (S2c)."""
    s2 = record['systems'].get('S2', {})
    row = checks.get((record['gse'], s2.get('row_ordinal'))) if s2.get('covered') is True else None
    if row is None or row['candidate']['status'] != 'OK' or row['candidate'].get('p') is None:
        return {'covered': False, 'status': 'NO_RELEASED_ROW'}
    estimate = -row['candidate']['estimate'] if s2.get('swapped') else row['candidate']['estimate']
    return {'covered': True, 'status': 'RELEASED', 'direction': 'higher' if estimate > 0 else 'lower' if estimate < 0 else 'none',
            'p': row['candidate']['p'], 'mean_diff': estimate, 'method': row['candidate']['method']}


def findings(records, system):
    """Direction AUROC and same-direction counts of one system (item 12)."""
    items = []
    for record in records:
        j = judged(record, system)
        result = record['systems'].get(system, {})
        if not j or not isinstance(result.get('p'), (int, float)) or not isinstance(j['significant'], bool):
            continue
        truth = record['direction'] if record.get('significant') is True else 'none' if record.get('significant') is False else None
        if truth not in ('higher', 'lower', 'none'):
            continue
        sign = {'higher': 1.0, 'lower': -1.0, 'none': 0.0}[j['direction']]
        items.append({'cluster': paper_key(record), 'truth': truth, 'z': signed_z(sign, result['p']),
                      'call': j['direction'] if j['significant'] and sign else 'none'})
    return finding_metrics(items) if items else None


def describe(records, system):
    out = with_interval(records, system, t1_metrics, paper_key, KEYS)
    out['findings'] = findings(records, system)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--bench', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True, help='run directories written by prepare.py')
    parser.add_argument('--checks', type=Path, required=True, help='checks records of the scored verifier run (export_as_checks)')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    plan = json.loads((args.root / 'plan.json').read_text())
    series, models = set(plan['series']), plan['models']
    framework = plan.get('framework', 'codex')
    records = [r for r in published_records(args.bench) if r['gse'] in series]
    checks = {(r['gse'], r['row_ordinal']): r for r in map(json.loads, (args.checks / 'dnam/rows.jsonl').read_text().splitlines())}
    reader, series_cache, runs, audits = PredictionReader(), {}, {}, []
    for model in models:
        for gse in sorted(series):
            run = run_dir(args.root, framework, model, gse)
            runs[(model, gse)] = read_run(run)
            entry = {'model': model, 'gse': gse, **audit_run(framework, run, plan['python']), 'comparisons': len(runs[(model, gse)]['arms']),
                     'rows': len(runs[(model, gse)]['rows']), 'problems': dict(Counter(runs[(model, gse)]['problems']))}
            entry['excluded'] = bool(entry['web_searches'] or entry['spawned_agents'] or entry['outside_paths'])
            audits.append(entry)
            if entry['excluded']:
                runs[(model, gse)] = {'rows': [], 'arms': {}, 'problems': ['EXCLUDED_BY_AUDIT']}
    for record in records:
        record['systems']['S2c'] = released_result(record, checks)
        for model in models:
            run = runs[(model, record['gse'])]
            found = match(record, run)
            if found is None:
                status = 'NO_RESULTS' if not run['rows'] else 'NO_MATCHING_COMPARISON'
                record['systems'][f'R:{model}:raw'] = record['systems'][f'R:{model}:recomputed'] = {'covered': False, 'status': status}
                continue
            row, swapped, overlap = found
            arms = run['arms'][(row['comparison'], row['condition'])]
            record['systems'][f'R:{model}:raw'] = raw_result(row, swapped, overlap)
            record['systems'][f'R:{model}:recomputed'] = recomputed_result(record, row, arms, swapped, overlap, reader, series_cache)
    reruns = [f'R:{m}:{kind}' for m in models for kind in ('raw', 'recomputed')]
    baselines = sorted({s for r in records for s in r['systems'] if s.startswith('S4:')})
    systems = ['S0', 'S1', 'S2c', 'S3'] + baselines + reruns
    every = shared(records, ['S1', 'S2c'] + reruns)
    result = {'records': len(records), 'series': len(series), 'papers': len({paper_key(r) for r in records}), 'models': models,
              'runs': {'planned': len(audits), 'completed': sum(a['status'] == '0' for a in audits),
                       'with_results': sum(a['rows'] > 0 for a in audits), 'excluded_by_audit': sum(a['excluded'] for a in audits),
                       'with_network_attempt_in_a_command': sum(a['network_attempts'] > 0 for a in audits),
                       'input_tokens': sum((a['usage'] or {}).get('input_tokens', 0) for a in audits),
                       'cached_input_tokens': sum((a['usage'] or {}).get('cached_input_tokens', 0) for a in audits),
                       'output_tokens': sum((a['usage'] or {}).get('output_tokens', 0) for a in audits),
                       'median_seconds': float(np.median([a['seconds'] for a in audits if a['seconds'] is not None] or [np.nan]))},
              'coverage_status': {s: dict(Counter(r['systems'].get(s, {}).get('status') for r in records)) for s in reruns},
              'own_coverage': {s: describe(records, s) for s in systems},
              'shared_agent_systems': {'records': len(every), 'papers': len({paper_key(r) for r in every}),
                                       'metrics': {s: describe(every, s) for s in ['S0', 'S1', 'S2c'] + baselines + reruns}},
              'model_with_and_without_data': {}}
    for model in models:
        covered = shared(records, [f'R:{model}:raw', f'R:{model}:recomputed'])
        names = ['S0', f'S4:openai:{model}', f'R:{model}:raw', f'R:{model}:recomputed', 'S1', 'S2c']
        result['model_with_and_without_data'][model] = {
            'records': len(covered), 'papers': len({paper_key(r) for r in covered}),
            'metrics': {s: describe(covered, s) for s in names if s == 'S0' or any(s in r['systems'] for r in covered)}}
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'agent_reruns.json').write_text(json.dumps(result, indent=1, default=float) + '\n')
    with (args.out / 'records.jsonl').open('w') as sink:
        for record in records:
            sink.write(json.dumps({'id': record['id'], 'gse': record['gse'], 'clock': record['clock_clockbase'],
                                   'direction': record.get('direction'), 'significant': record.get('significant'),
                                   'paper': paper_key(record),
                                   'systems': {s: record['systems'].get(s) for s in systems if s != 'S0'}}, default=float) + '\n')
    with (args.out / 'runs.csv').open('w', newline='') as sink:
        fields = ['model', 'gse', 'status', 'seconds', 'commands', 'web_searches', 'spawned_agents', 'network_attempts',
                  'outside_paths', 'excluded',
                  'comparisons', 'rows', 'problems', 'usage']
        writer = csv.DictWriter(sink, fieldnames=fields)
        writer.writeheader()
        for entry in audits:
            writer.writerow({k: json.dumps(entry[k]) if isinstance(entry[k], (dict, list)) else entry[k] for k in fields})
    print(json.dumps(result['runs'], indent=1))
    for name, block in (('own coverage', result['own_coverage']), (f"shared ({len(every)} results)", result['shared_agent_systems']['metrics'])):
        print(name)
        for system, m in block.items():
            f = m.get('findings') or {}
            print(f"  {system:32s} covered {m['covered']:3d}/{m['records']:3d}  call {m['call_agreement'] if m['call_agreement'] is None else round(m['call_agreement'], 3)} "
                  f"(n {m['call_n']})  balanced {m['call_balanced_accuracy'] if m['call_balanced_accuracy'] is None else round(m['call_balanced_accuracy'], 3)}  "
                  f"direction {m['direction_agreement'] if m['direction_agreement'] is None else round(m['direction_agreement'], 3)} (n {m['direction_n']})  "
                  f"AUROC {f.get('direction_auroc') if f.get('direction_auroc') is None else round(f['direction_auroc'], 3)}")


if __name__ == '__main__':
    main()
