"""Build the Benchmark A units and the blinded input packages (PROTOCOL.md v2).

Units: within a covered series, the distinct unordered pairs of arm sets among comparisons whose arms are the agent's
own (membership UNIQUE from recovery, replay or both, with effect reproduction PASS for at least one clock). Every other
comparison of a covered series is listed with its exclusion reason.

Packages (no verifier output in any of them):
- names/<set>/<GSE>.json: each unit's name variants, the series' GEO field values and the expert vocabulary; no arms.
- samples/<set>/<GSE>.json: every sample in a unit arm with its GEO fields and expert labels, and which arms hold it.

  python -m benchmark.published_clock_2026_10_05.groupbench.prepare --root GROUPBENCH --out DIR
"""
import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from audit_system import paths

GEO_META = {'A1': paths.GEO / 'GEO_mm_rna/processed_metadata_preharmonize',
            'A2': paths.GEO / 'GEO_hs_dnam/processed_metadata_preharmonize'}
MODALITY = {'A1': 'mouse RNA-seq', 'A2': 'human DNAm'}


def own_arms(item):
    """Arms that are the agent's own: unique membership from recovery or replay, reproduced for at least one clock."""
    if item['membership_status'] != 'UNIQUE' or item['membership_source'] not in ('recovery', 'replay', 'replay+recovery'):
        return False, f"membership {item['membership_status']} from {item['membership_source']}"
    if not item['case_gsms'] or not item['control_gsms']:
        return False, 'empty arm'
    if not any(c.get('effect_reproduction') == 'PASS' for c in item['clocks']):
        return False, 'stored effect not reproduced on these arms'
    return True, None


def expert_labels(root):
    """{set: {gsm: label}} from the A1 mapping plus annotation and the A2 ComputAgeBench labels."""
    annotation = pd.read_excel(root / 'a1/source/Data_annotation_absolute_rodents.xlsx').set_index('Sample')
    a1 = {}
    for path in sorted((root / 'a1/mapping').glob('GSE*.csv')):
        mapping = pd.read_csv(path, dtype=str).fillna('')
        for sample, gsm in zip(mapping['sample'], mapping['gsm']):
            if gsm:
                row = annotation.loc[sample]
                a1[gsm] = {'gse': path.stem, 'sample': sample, **{k: (None if pd.isna(row[k]) else row[k]) for k in (
                    'Strain', 'Sex', 'Tissue', 'Tissue.details', 'Intervention', 'Intervention.type', 'Age.days')}}
    labels = pd.read_csv(root / 'a2/labels.csv', dtype=str).fillna('')
    a2 = defaultdict(dict)
    for _, r in labels.iterrows():
        a2[r.gsm][r.dataset] = {k: r[k] for k in ('condition', 'class', 'age', 'tissue', 'sex')}
    meta = pd.read_csv(root / 'a2/computage_bench_meta_3cf5268.tsv', sep='\t', dtype=str).fillna('')
    cell = {(r.DatasetID, r.SampleID): r.CellType for r in meta.itertuples()}
    for gsm, by_dataset in a2.items():
        for dataset, value in by_dataset.items():
            value['cell_type'] = cell.get((dataset, gsm), '')
    return {'A1': a1, 'A2': dict(a2)}


def geo_table(set_name, gse):
    path = GEO_META[set_name] / f'{gse}.csv'
    return pd.read_csv(path, dtype=str).fillna('') if path.exists() else pd.DataFrame()


def build(root, out):
    comparisons = [json.loads(line) for line in (root / 'agent_comparisons_2026-10-06.jsonl').read_text().splitlines()]
    labels = expert_labels(root)
    a1_series = sorted(p.stem for p in (root / 'a1/mapping').glob('GSE*.csv'))
    a2_ids = set(pd.read_csv(root / 'a2/labels.csv', dtype=str).dataset)
    covered = {'A1': set(a1_series), 'A2': set()}
    for item in comparisons:  # A2: atlas DNAm series named in ComputAgeBench or holding its samples under another ID
        if item['modality'] == MODALITY['A2'] and (item['gse'] in a2_ids or any(
                g in labels['A2'] for g in item['case_gsms'] + item['control_gsms'])):
            covered['A2'].add(item['gse'])
    units, excluded, counts = [], [], Counter()
    for set_name in ('A1', 'A2'):
        for gse in sorted(covered[set_name]):
            series = [c for c in comparisons if c['modality'] == MODALITY[set_name] and c['gse'] == gse]
            pairs = {}
            for item in series:
                ok, reason = own_arms(item)
                if not ok:
                    excluded.append({'set': set_name, 'gse': gse, 'comparison_id': item['comparison_id'],
                                     'comparison_name': item['comparison_name'], 'reason': reason,
                                     'membership_source': item['membership_source']})
                    counts[(set_name, 'excluded')] += 1
                    continue
                key = tuple(sorted([tuple(sorted(item['case_gsms'])), tuple(sorted(item['control_gsms']))]))
                pairs.setdefault(key, []).append(item)
            for index, (key, items) in enumerate(sorted(pairs.items())):
                arm_a, arm_b = list(key[0]), list(key[1])
                variants = [{'comparison_id': it['comparison_id'], 'comparison_name': it['comparison_name'],
                             'condition_name': it['condition_name'], 'membership_source': it['membership_source'],
                             'case_is_arm_a': sorted(it['case_gsms']) == arm_a,
                             'clocks': it['clocks']} for it in items]
                orientations = {v['case_is_arm_a'] for v in variants}
                units.append({'unit_id': f'{set_name}:{gse}:{index}', 'set': set_name, 'gse': gse,
                              'arm_a': arm_a, 'arm_b': arm_b, 'variants': variants,
                              'orientation_conflict': len(orientations) > 1})
                counts[(set_name, 'units')] += 1
                counts[(set_name, 'comparisons_in_units')] += len(items)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'units.jsonl').write_text(''.join(json.dumps(u) + '\n' for u in units))
    (out / 'excluded.jsonl').write_text(''.join(json.dumps(e) + '\n' for e in excluded))
    for set_name in ('A1', 'A2'):
        (out / 'names' / set_name).mkdir(parents=True)
        (out / 'samples' / set_name).mkdir(parents=True)
        for gse in sorted(covered[set_name]):
            series_units = [u for u in units if u['set'] == set_name and u['gse'] == gse]
            if not series_units:
                continue
            geo = geo_table(set_name, gse)
            fields = {c: dict(Counter(geo[c]).most_common(40)) for c in geo.columns if c not in ('GSM',)}
            series_labels = {g: v for g, v in labels[set_name].items()
                             if (v.get('gse') == gse if set_name == 'A1' else True)}
            if set_name == 'A1':
                vocabulary = {k: sorted({str(v[k]) for v in series_labels.values()}) for k in (
                    'Intervention', 'Intervention.type', 'Strain', 'Sex', 'Tissue.details', 'Age.days')}
            else:
                in_series = {g for u in series_units for g in u['arm_a'] + u['arm_b']}
                found = [d for g in in_series for d in labels['A2'].get(g, {}).values()]
                vocabulary = {k: sorted({d[k] for d in found}) for k in ('condition', 'class', 'tissue', 'cell_type', 'sex')}
            names = {'set': set_name, 'gse': gse, 'geo_fields': fields, 'expert_vocabulary': vocabulary,
                     'units': [{'unit_id': u['unit_id'], 'names': [{'comparison_name': v['comparison_name'],
                                                                     'condition_name': v['condition_name']}
                                                                    for v in u['variants']]} for u in series_units]}
            (out / 'names' / set_name / f'{gse}.json').write_text(json.dumps(names, indent=1) + '\n')
            geo_by_gsm = {r['GSM']: {k: v for k, v in r.items() if k != 'GSM'} for r in geo.to_dict('records')} if len(geo) else {}
            samples = {}
            for u in series_units:
                for arm in ('arm_a', 'arm_b'):
                    for g in u[arm]:
                        entry = samples.setdefault(g, {'gsm': g, 'geo': geo_by_gsm.get(g, {}),
                                                       'expert': labels[set_name].get(g), 'in_arms': []})
                        entry['in_arms'].append(f"{u['unit_id']}:{arm}")
            package = {'set': set_name, 'gse': gse, 'all_geo_samples': len(geo_by_gsm),
                       'samples': sorted(samples.values(), key=lambda s: s['gsm'])}
            (out / 'samples' / set_name / f'{gse}.json').write_text(json.dumps(package, indent=1, default=str) + '\n')
    summary = {f'{s}|{k}': n for (s, k), n in sorted(counts.items())}
    summary['covered_series'] = {s: len(v) for s, v in covered.items()}
    summary['series_with_units'] = {s: len({u['gse'] for u in units if u['set'] == s}) for s in ('A1', 'A2')}
    summary['exclusion_reasons'] = dict(Counter(f"{e['set']}|{e['reason']}" for e in excluded))
    (out / 'summary.json').write_text(json.dumps(summary, indent=1) + '\n')
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.root, args.out), indent=1))


if __name__ == '__main__':
    main()
