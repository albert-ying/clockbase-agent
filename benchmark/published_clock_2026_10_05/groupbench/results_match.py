"""Benchmark B of PROTOCOL.md: published results against the agent's result for the same groups.

A published contrast whose groups were assigned to samples (engine groups status AGREED or MANUAL) is linked to the
agent comparisons of its series: the best comparison whose arms overlap the published arms with Jaccard >= MATCH in
both arms, in either orientation, regardless of clock. The primary analysis uses only arms that are the agent's own
(replay or recovery); the sensitivity analysis also admits arms reconstructed by the verifier. For a linked pair the
agent's result is the row of the published clock: the historical row (direction and P as stored by the agent) and
the verified recomputation. Calls are three-way: higher or lower when significant (the paper's stated threshold,
P < 0.05 otherwise), none when not significant.

  python -m benchmark.published_clock_2026_10_05.groupbench.results_match --records RUN/records.jsonl [...]
      --comparisons COMPARISONS.jsonl --ground-truth FINAL.jsonl --out DIR
"""
import argparse
import json
from collections import Counter
from pathlib import Path

from ..engine.stats import threshold

MATCH, EXACT = 0.5, 0.9
CONTAIN = 0.8  # protocol v2 primary match: published arm inside the agent arm, and agent arm from the published group
DEFAULT_ALPHA = 0.05


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 0.0


def call(direction, significant):
    if significant is None or direction not in ('higher', 'lower', 'none'):
        return None
    return direction if significant and direction in ('higher', 'lower') else 'none'


def sign(value):
    if value is None:
        return None
    return 'higher' if value > 0 else 'lower' if value < 0 else 'none'


def link(record, comparisons, own_only, clock=None):
    """The best agent comparison for a published contrast (only comparisons with a row of `clock` when given)."""
    case, control = record['groups']['case_gsms'], record['groups']['control_gsms']
    best = None
    for item in comparisons:
        if own_only and not item['agent_arms']:
            continue
        if clock and not any(c['clock'] == clock for c in item['clocks']):
            continue
        if not item['case_gsms'] or not item['control_gsms']:
            continue
        for swapped in (False, True):
            a_case, a_control = (item['control_gsms'], item['case_gsms']) if swapped else (item['case_gsms'], item['control_gsms'])
            jc, jt = jaccard(case, a_case), jaccard(control, a_control)
            if min(jc, jt) < MATCH:
                continue
            inside = min(len(set(case) & set(a_case)) / len(set(case)), len(set(control) & set(a_control)) / len(set(control)))
            purity = min(len(set(case) & set(a_case)) / len(set(a_case)), len(set(control) & set(a_control)) / len(set(a_control)))
            contained = inside >= CONTAIN and purity >= CONTAIN
            key = (contained, min(jc, jt), jc + jt, not swapped)
            if best is None or key > best[0]:
                best = (key, {'comparison_id': item['comparison_id'], 'swapped': swapped, 'jaccard_case': jc,
                              'jaccard_control': jt, 'published_inside_agent': inside, 'agent_from_published': purity,
                              'contained': contained,
                              'tier': 'contained' if contained else 'jaccard_only',
                              'membership_source': item['membership_source'], 'item': item})
    return best[1] if best else None


def agent_result(record, linked):
    """Historical and verified calls of the linked comparison's row for the published clock."""
    rows = [c for c in linked['item']['clocks'] if c['clock'] == record.get('clock_clockbase')]
    if not rows:
        return {'status': 'CLOCK_NOT_IN_COMPARISON'}
    row = rows[0]
    alpha = threshold(record) or DEFAULT_ALPHA
    flip = {'higher': 'lower', 'lower': 'higher', 'none': 'none', None: None}
    historical = row.get('historical') or {}
    h_direction = sign(historical.get('direction') if historical.get('direction') is not None else historical.get('effect_size'))
    h_p = historical.get('pvalue')
    recomputed = row.get('recomputed') or {}
    v_direction = sign(recomputed.get('mean_diff'))
    v_p = recomputed.get('p')
    if linked['swapped']:
        h_direction, v_direction = flip[h_direction], flip[v_direction]
    return {'status': 'OK', 'row_ordinal': row['row_ordinal'], 'final_class': row.get('final_class'),
            'soft_flags': row.get('soft_flags'), 'alpha': alpha, 'alpha_assumed': threshold(record) is None,
            'historical': {'direction': h_direction, 'p': h_p,
                           'call': call(h_direction, (h_p < alpha) if isinstance(h_p, (int, float)) else None)},
            'verified': {'direction': v_direction, 'p': v_p,
                         'call': call(v_direction, (v_p < alpha) if isinstance(v_p, (int, float)) else None)}}


def run(records, comparisons, keep):
    by_series = {}
    for item in comparisons:
        by_series.setdefault((item['modality'], item['gse']), []).append(item)
    out, counts = [], Counter()
    for record in records:
        if record['id'] not in keep or not record.get('clock_clockbase'):
            continue
        counts['records_with_clockbase_clock'] += 1
        series = by_series.get((record['modality'], record['gse']), [])
        entry = {'id': record['id'], 'gse': record['gse'], 'modality': record['modality'],
                 'clock': record['clock_clockbase'], 'measure': record.get('measure'),
                 'published': {'direction': record.get('direction'), 'significant': record.get('significant'),
                               'p_value': record.get('p_value'), 'threshold': record.get('threshold'),
                               'call': call(record.get('direction'), record.get('significant'))},
                 'series_in_atlas': bool(series), 'groups_status': record['groups'].get('status')}
        if not series:
            entry['status'] = 'SERIES_NOT_IN_ATLAS'
        elif record['groups'].get('status') not in ('AGREED', 'MANUAL'):
            entry['status'] = 'PUBLISHED_GROUPS_UNRESOLVED'
        else:
            for label, own_only in (('primary', True), ('sensitivity', False)):
                linked = link(record, series, own_only)
                if linked:
                    # Rows of different clocks can differ in arms (missing predictions); take the row of the published
                    # clock from the best-overlapping comparison that has it.
                    with_clock = link(record, series, own_only, clock=record['clock_clockbase'])
                    result = agent_result(record, with_clock) if with_clock else {'status': 'CLOCK_NOT_IN_COMPARISON'}
                    entry[label] = {k: v for k, v in linked.items() if k != 'item'} | {
                        'result': result, 'result_comparison_id': with_clock['comparison_id'] if with_clock else None}
            entry['status'] = 'LINKED' if 'primary' in entry else 'LINKED_RECONSTRUCTED_ONLY' if 'sensitivity' in entry else 'NOT_LINKED'
        counts[entry['status']] += 1
        out.append(entry)
    return out, dict(counts)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--records', type=Path, nargs='+', required=True)
    parser.add_argument('--comparisons', type=Path, required=True)
    parser.add_argument('--ground-truth', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    keep = {json.loads(line)['id'] for line in args.ground_truth.read_text().splitlines()}
    records = [json.loads(line) for path in args.records for line in path.read_text().splitlines()]
    comparisons = [json.loads(line) for line in args.comparisons.read_text().splitlines()]
    out, counts = run(records, comparisons, keep)
    (args.out / 'linked.jsonl').write_text(''.join(json.dumps(e) + '\n' for e in out))
    (args.out / 'counts.json').write_text(json.dumps(counts, indent=1) + '\n')
    print(json.dumps(counts, indent=1))


if __name__ == '__main__':
    main()
