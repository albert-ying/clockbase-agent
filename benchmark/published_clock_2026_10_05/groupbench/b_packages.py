"""Blinded review packages for Benchmark B discordances (PROTOCOL.md v2), one file per discordant published contrast.

Each package holds the published record (source quote, group definitions, reported statistics), the published groups as
assigned to samples, the agent's comparison (names, arms, stored test, direction and P) and per-sample GEO metadata.
It holds no verifier output (classes, flags, repairs) and no recomputed values.

  python -m benchmark.published_clock_2026_10_05.groupbench.b_packages --stats B_STATS.json --linked LINKED.jsonl
      --records RUN/records.jsonl [...] --comparisons COMPARISONS.jsonl --out DIR
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from audit_system import paths

GEO = {'human DNAm': paths.GEO / 'GEO_hs_dnam/processed_metadata_preharmonize',
       'human RNA-seq': paths.GEO / 'GEO_hs_rna/processed_metadata_preharmonize',
       'mouse RNA-seq': paths.GEO / 'GEO_mm_rna/processed_metadata_preharmonize'}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--stats', type=Path, required=True)
    parser.add_argument('--linked', type=Path, required=True)
    parser.add_argument('--records', type=Path, nargs='+', required=True)
    parser.add_argument('--comparisons', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    stats = json.loads(args.stats.read_text())
    linked = {e['id']: e for e in map(json.loads, args.linked.read_text().splitlines())}
    records = {r['id']: r for path in args.records for r in map(json.loads, path.read_text().splitlines())}
    comparisons = {c['comparison_id']: c for c in map(json.loads, args.comparisons.read_text().splitlines())}
    args.out.mkdir(parents=True, exist_ok=False)
    index = []
    for item in stats['historical']['discordant']:
        record, entry = records[item['id']], linked[item['id']]
        comparison = comparisons[item['comparison_id']]
        geo_path = GEO[record['modality']] / f"{record['gse']}.csv"
        geo = pd.read_csv(geo_path, dtype=str).fillna('') if geo_path.exists() else pd.DataFrame()
        involved = set(record['groups']['case_gsms']) | set(record['groups']['control_gsms']) | \
            set(comparison['case_gsms']) | set(comparison['control_gsms'])
        samples = [{k: v for k, v in row.items()} for row in geo.to_dict('records') if row.get('GSM') in involved]
        clock = next(c for c in comparison['clocks'] if c['clock'] == record['clock_clockbase'])
        package = {
            'published': {k: record.get(k) for k in ('id', 'gse', 'modality', 'tissue', 'clock_reported', 'clock_clockbase',
                                                     'measure', 'case_label', 'control_label', 'case_definition',
                                                     'control_definition', 'n_case', 'n_control', 'direction', 'significant',
                                                     'threshold', 'p_value', 'effect', 'test', 'covariates')},
            'published_source': {k: (record.get('source') or {}).get(k) for k in ('citation', 'year', 'doi', 'pmid', 'location', 'quote')},
            'published_groups_assigned': {'case_gsms': record['groups']['case_gsms'], 'control_gsms': record['groups']['control_gsms']},
            'calls': {'published': item['published'], 'agent': item['agent']},
            'agent_comparison': {'comparison_name': comparison['comparison_name'], 'condition_name': comparison['condition_name'],
                                 'case_gsms': comparison['case_gsms'], 'control_gsms': comparison['control_gsms'],
                                 'clock': clock['clock'], 'stored_test': (clock.get('historical') or {}).get('statistical_test'),
                                 'stored_direction': (clock.get('historical') or {}).get('direction'),
                                 'stored_p': (clock.get('historical') or {}).get('pvalue'),
                                 'stored_effect': (clock.get('historical') or {}).get('effect_size')},
            'overlap': {k: entry['primary'].get(k) for k in ('jaccard_case', 'jaccard_control', 'published_inside_agent',
                                                             'agent_from_published', 'swapped')},
            'samples': samples,
        }
        name = f"{item['id']}.json"
        (args.out / name).write_text(json.dumps(package, indent=1, default=str) + '\n')
        index.append({'id': item['id'], 'file': name})
    (args.out / 'index.json').write_text(json.dumps(index, indent=1) + '\n')
    print(len(index), 'packages')


if __name__ == '__main__':
    main()
