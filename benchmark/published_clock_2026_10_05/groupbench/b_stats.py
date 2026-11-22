"""Benchmark B statistics (PROTOCOL.md v2): published results against the agent's result for the same samples.

Calls are three-way (higher, lower, none). The paper's call uses its reported P < 0.05 when a P value is reported and
its stated significance otherwise; the agent's call uses its stored P < 0.05 and stored direction (historical row) or
the verifier's recomputation (released row). Reported: coverage, the 3x3 table with Cohen's kappa, call agreement,
direction agreement among contrasts significant on at least one side, all with a paper-cluster bootstrap; and the list
of discordant pairs for adjudication.

  python -m benchmark.published_clock_2026_10_05.groupbench.b_stats --linked LINKED.jsonl --out DIR
"""
import argparse
import json
from collections import Counter
from pathlib import Path

import numpy as np

SEED = 20261006
CALLS = ('higher', 'lower', 'none')


def paper_of(record_id):
    return record_id.split('-GSE')[0]


def published_call(entry):
    published = entry['published']
    direction = published.get('direction')
    p = published.get('p_value')
    significant = (p < 0.05) if isinstance(p, (int, float)) else published.get('significant')
    if significant is None or direction not in CALLS:
        return None
    return direction if significant and direction in ('higher', 'lower') else 'none'


def agent_call(result, which):
    side = result.get(which) or {}
    direction, p = side.get('direction'), side.get('p')
    if direction not in CALLS or not isinstance(p, (int, float)):
        return None
    return direction if p < 0.05 and direction in ('higher', 'lower') else 'none'


def kappa(pairs):
    if not pairs:
        return None
    n = len(pairs)
    observed = sum(a == b for a, b in pairs) / n
    pa, pb = Counter(a for a, _ in pairs), Counter(b for _, b in pairs)
    expected = sum(pa[c] * pb[c] for c in CALLS) / (n * n)
    return (observed - expected) / (1 - expected) if expected < 1 else None


def metrics(pairs):
    table = {p: {a: sum(1 for x, y in pairs if x == p and y == a) for a in CALLS} for p in CALLS}
    significant = [(x, y) for x, y in pairs if x != 'none' or y != 'none']
    direction = [(x, y) for x, y in pairs if x in ('higher', 'lower') and y in ('higher', 'lower')]
    return {'n': len(pairs), 'table_published_by_agent': table, 'kappa': kappa(pairs),
            'call_agreement': sum(x == y for x, y in pairs) / len(pairs) if pairs else None,
            'agreement_significant_either_side': (sum(x == y for x, y in significant) / len(significant)) if significant else None,
            'n_significant_either_side': len(significant),
            'direction_agreement_both_significant': (sum(x == y for x, y in direction) / len(direction)) if direction else None,
            'n_both_significant': len(direction)}


def cluster_bootstrap(items, statistic, replicates=10000, seed=SEED):
    clusters = {}
    for item in items:
        clusters.setdefault(item['paper'], []).append(item['pair'])
    keys = sorted(clusters)
    if len(keys) < 2:
        return None
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(replicates):
        sample = [p for i in rng.integers(0, len(keys), len(keys)) for p in clusters[keys[i]]]
        value = statistic(sample)
        if value is not None:
            values.append(value)
    return [float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))] if values else None


def run(linked, tier='contained'):
    coverage = Counter(e['status'] for e in linked)
    out = {'coverage': dict(coverage)}
    for which in ('historical', 'verified'):
        items, discordant = [], []
        for entry in linked:
            link = entry.get('primary')
            if not link or (tier and link['tier'] != tier) or link['result']['status'] != 'OK':
                continue
            pub, agent = published_call(entry), agent_call(link['result'], which)
            if pub is None or agent is None:
                continue
            items.append({'paper': paper_of(entry['id']), 'pair': (pub, agent), 'id': entry['id']})
            if pub != agent:
                discordant.append({'id': entry['id'], 'gse': entry['gse'], 'clock': entry['clock'], 'published': pub,
                                   'agent': agent, 'comparison_id': link['result_comparison_id'],
                                   'final_class': link['result'].get('final_class')})
        pairs = [i['pair'] for i in items]
        block = metrics(pairs)
        block['papers'] = len({i['paper'] for i in items})
        block['call_agreement_interval'] = cluster_bootstrap(items, lambda s: sum(x == y for x, y in s) / len(s) if s else None)
        block['kappa_interval'] = cluster_bootstrap(items, kappa)
        block['discordant'] = discordant
        out[which] = block
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--linked', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--tier', default='contained', help="'contained' (primary) or '' for all Jaccard links")
    args = parser.parse_args()
    linked = [json.loads(line) for line in args.linked.read_text().splitlines()]
    result = run(linked, args.tier or None)
    args.out.mkdir(parents=True, exist_ok=True)
    name = f"b_stats_{args.tier or 'all'}.json"
    (args.out / name).write_text(json.dumps(result, indent=1) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'discordant'}, indent=1, default=str)[:3000])


if __name__ == '__main__':
    main()
