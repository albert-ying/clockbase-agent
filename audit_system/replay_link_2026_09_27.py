"""Link every atlas row to the replayed statistical test that produced it.

A logged test explains an atlas row when (1) its P value (two-sample tests) or one of its model
coefficient P values (statsmodels fits) equals the stored P to relative 1e-9, (2) its values come
from the row's clock and (3) every member GSM belongs to the row's study. The case arm is the arm
whose mean difference has the sign of the stored effect (two-sample tests) or the arm coded 1 in the
matched coefficient's design column (model fits). A link to a trace whose prompt named another
study is recorded as misfiled output.

  python -m audit_system.replay_link_2026_09_27 --modality dnam --replay DIR --out DIR
"""
import argparse
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system.membership_recovery_2026_09_27 import CORRECTIONS, MODALITIES, sha256
from audit_system.trace_replay_2026_09_27 import _members, _value_index

SCHEMA = 'clockbase-replay-link-v1'
P_REL = 1e-9


def bucket(p):
    return round(math.log10(p), 6) if p and p > 0 else None


def load_tests(replay_root):
    """All logged tests with their trace, keyed by rounded log10(P) for candidate lookup."""
    by_bucket = defaultdict(list)
    replays = Path(replay_root) / 'replays'
    for folder in sorted(replays.iterdir()):
        path = folder / 'memberships.jsonl'
        if not path.exists():
            continue
        prompt = json.loads((folder / 'result.json').read_text())['prompt_gse']
        for line_no, line in enumerate(open(path)):
            test = json.loads(line)
            test['trace'] = folder.name
            test['prompt_gse'] = prompt
            test['line'] = line_no
            ps = [test.get('p')] if not test['fn'].startswith('sm.') else list((test.get('pvalues') or {}).values())
            for p in ps:
                b = bucket(p)
                if b is not None and math.isfinite(b):
                    by_bucket[b].append(test)
    return by_bucket


def candidates(by_bucket, p):
    b = bucket(p)
    if b is None:
        return []
    out = []
    for key in (b - 1e-6, b, b + 1e-6):
        out.extend(by_bucket.get(round(key, 6), []))
    return out


def explain(row, test, index, study_gsms):
    """Return a link record when this test explains the row, else None."""
    p, effect, clock = row['pvalue'], row['effect_size'], row['clock']
    if test['fn'].startswith('sm.'):
        matched = [name for name, value in (test.get('pvalues') or {}).items()
                   if value is not None and abs(value - p) <= P_REL * p]
        if not matched:
            return None
        endog_clock, members = test.get('clock'), test.get('members')
        if members is None and test.get('endog'):
            endog_clock, members = _members(test['endog'], index)
        clock_in_formula = bool(re.match(rf'\s*Q?\(?[\'"]?{re.escape(clock)}[\'"]?\)?\s*~', test.get('formula') or ''))
        if not members or not set(members) <= study_gsms or not (endog_clock == clock or clock_in_formula):
            return None
        term = matched[0]
        names = test.get('exog_names') or []
        record = {'kind': 'model', 'fn': test['fn'], 'formula': test.get('formula'), 'term': term,
                  'term_param': (test.get('params') or {}).get(term),
                  'members': members, 'n': len(members), 'covariates': [n for n in names if n not in (term, 'Intercept', 'const')]}
        exog = test.get('exog')
        if exog is not None and term in names:
            column = np.asarray(exog)[:, names.index(term)]
            if set(np.unique(column)) <= {0.0, 1.0}:
                record['case_gsms'] = [g for g, x in zip(members, column) if x == 1]
                record['control_gsms'] = [g for g, x in zip(members, column) if x == 0]
            else:
                record['continuous_term'] = True
        return record
    groups = test.get('groups') or []
    if abs((test.get('p') or -1) - p) > P_REL * p:
        return None
    if len(groups) == 2 and all(g.get('gsms') for g in groups):
        clocks = {g.get('clock') for g in groups}
        if clocks != {clock} and None not in clocks:
            return None
        a, b = groups[0]['gsms'], groups[1]['gsms']
        if not (set(a) | set(b)) <= study_gsms:
            return None
        values = index.get(clock, {})
        inverse = {next(iter(v)): k for k, v in values.items() if len(v) == 1}
        mean_a = np.mean([inverse[g] for g in a if g in inverse]) if a else np.nan
        mean_b = np.mean([inverse[g] for g in b if g in inverse]) if b else np.nan
        if np.isfinite(effect) and np.sign(mean_a - mean_b) != np.sign(effect):
            a, b = b, a
        paired = {'paired_case': list(a), 'paired_control': list(b)} if test['fn'] == 'ttest_rel' else {}
        return {**paired, 'kind': 'two_sample', 'fn': test['fn'], 'kwargs': test.get('kwargs'), 'case_gsms': sorted(a),
                'control_gsms': sorted(b), 'n_case': len(a), 'n_control': len(b),
                'first_argument_is_case': groups[0]['gsms'] == a}
    if len(groups) == 1 and groups[0].get('gsms'):
        members = groups[0]['gsms']
        if groups[0].get('clock') not in (clock, None) or not set(members) <= study_gsms:
            return None
        return {'kind': 'one_sample', 'fn': test['fn'], 'members': sorted(members), 'n': len(members)}
    return None


def run(modality, replay_root, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    replay_root = Path(replay_root)
    table = CORRECTIONS / f"2026-09-08-{MODALITIES[modality]['table']}-atlas_WORKING_NOT_VALIDATED.csv"
    atlas = pd.read_csv(table, low_memory=False)
    atlas['pvalue'] = pd.to_numeric(atlas.pvalue, errors='coerce')
    atlas['effect_size'] = pd.to_numeric(atlas.effect_size, errors='coerce')
    by_bucket = load_tests(replay_root)
    indexes, members_of = {}, {}
    records, counts = [], defaultdict(int)
    for ordinal, row in enumerate(atlas.to_dict('records')):
        gse, p = row['gse_id'], row['pvalue']
        record = {'row_ordinal': ordinal, 'gse_id': gse, 'comparison_name': row['comparison_name'],
                  'condition_name': row['condition_name'], 'clock': row['clock'], 'stored_p': p,
                  'stored_effect': row['effect_size'], 'stored_test': row['statistical_test']}
        if not (isinstance(p, float) and 0 < p < 1):
            record['link_status'] = 'STORED_P_UNUSABLE'
        else:
            study_input = replay_root / 'inputs' / f'{gse}.csv'
            if not study_input.exists():
                record['link_status'] = 'STUDY_INPUT_UNAVAILABLE'
            else:
                if gse not in indexes:
                    if len(indexes) > 64:
                        indexes.clear()
                    indexes[gse] = _value_index(study_input)
                    column = next(c for c in pd.read_csv(study_input, nrows=0).columns if c.lower() in ('gsm_id', 'gsm', 'sample_id'))
                    members_of[gse] = set(pd.read_csv(study_input, usecols=[column])[column].astype(str))
                links = []
                for test in candidates(by_bucket, p):
                    linked = explain(row, test, indexes[gse], members_of[gse])
                    if linked:
                        linked.update(trace=test['trace'], trace_prompt_gse=test['prompt_gse'], test_line=test['line'])
                        links.append(linked)
                distinct = {json.dumps({k: v for k, v in l.items() if k in ('case_gsms', 'control_gsms', 'members')},
                                       sort_keys=True) for l in links}
                if not links:
                    record['link_status'] = 'UNLINKED'
                elif len(distinct) == 1:
                    record['link_status'] = 'LINKED'
                    record['link'] = links[0]
                    record['link_traces'] = sorted({l['trace'] for l in links})
                    record['misfiled'] = all(l['trace_prompt_gse'] != gse for l in links)
                else:
                    record['link_status'] = 'AMBIGUOUS'
                    record['link_candidates'] = len(distinct)
        counts[record['link_status']] += 1
        records.append(record)
    with open(out_dir / 'row_links.jsonl', 'w') as sink:
        for record in records:
            sink.write(json.dumps(record, default=str) + '\n')
    linked = [r for r in records if r['link_status'] == 'LINKED']
    manifest = {'schema': SCHEMA, 'modality': modality, 'table': str(table), 'table_sha256': sha256(table),
                'replay_root': str(replay_root), 'module_sha256': sha256(__file__), 'rows': len(records),
                'link_status': dict(counts), 'linked_misfiled': sum(r['misfiled'] for r in linked),
                'linked_kinds': pd.Series([r['link']['kind'] for r in linked]).value_counts().to_dict() if linked else {},
                'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'rows_sha256': sha256(out_dir / 'row_links.jsonl')}
    (out_dir / 'manifest.json').write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k: manifest[k] for k in ('link_status', 'linked_misfiled', 'linked_kinds')}, indent=1))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--modality', choices=sorted(MODALITIES), required=True)
    parser.add_argument('--replay', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    run(args.modality, args.replay, args.out)


if __name__ == '__main__':
    main()
