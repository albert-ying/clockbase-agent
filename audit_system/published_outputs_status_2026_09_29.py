"""Trace every published ranked item (S14, Figure 3 bars, Tables 1/2) to its atlas row and report its
verifier-v3 fate: final class, repair actions, recomputed effect/P/FDR, and whether it remains
significant in the same direction.

Items are matched to repaired atlas rows within the same study and modality by exact stored effect
size (and P where the published file records it). If no exact match remains, an otherwise identical
row with the stored effect sign reversed is reported as SIGN_REVERSED. Unmatched and ambiguous
items are reported, not guessed.

  python -m audit_system.published_outputs_status_2026_09_29 --repairs dnam=DIR --repairs hs_rna=DIR \
      --repairs mm_rna=DIR --out DIR
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system import paths
from audit_system.membership_recovery_2026_09_27 import sha256

MODALITY_LABELS = {'dnam': 'human DNAm', 'hs_rna': 'human RNA-seq', 'mm_rna': 'mouse RNA-seq'}
PUBLISHED = {
    'S14': (paths.PUBLISHED / 'historical_s14_identity.csv', {'effect': 'cb_effect_size', 'p': None, 'claim': 'cb_direction'}),
    'Figure3_bars': (paths.PUBLISHED / 'top_hits_combined_modalities_v2.csv',
                     {'effect': 'log2FoldChange', 'p': 'pvalue', 'claim': 'direction'}),
    'Table1_pro_aging': (paths.PUBLISHED / 'table1_pro_aging_pooled_top30.csv',
                         {'effect': 'effect_size', 'p': 'pvalue', 'claim': 'direction'}),
    'Table2_anti_aging': (paths.PUBLISHED / 'table2_anti_aging_pooled_top30.csv',
                          {'effect': 'effect_size', 'p': 'pvalue', 'claim': 'direction'}),
}
FDR = 0.05


def close(a, b, rel=1e-9):
    return np.isfinite(a) and np.isfinite(b) and abs(a - b) <= rel * max(abs(a), abs(b), 1e-300)


def claimed_direction(value):
    """+1 = claimed to raise the clock (pro-aging), -1 = claimed to lower it (anti-aging)."""
    text = str(value).lower()
    if 'anti' in text or 'decel' in text:
        return -1
    if 'pro' in text or 'accel' in text:
        return 1
    number = pd.to_numeric(value, errors='coerce')
    return int(np.sign(number)) if pd.notna(number) and number != 0 else None


def fate(row, claim):
    cls = row['v3_final_class']
    if cls == 'DUPLICATE':  # the same test as its representative row, with the same recomputed values
        cls = row.get('v3_representative_class')
    actions = str(row.get('v3_actions') or '')
    flipped = 'SWAP_ARMS_CONFIRMED_INVERSION' in actions
    new_sign = np.sign(row['v3_cohen_d']) if pd.notna(row.get('v3_cohen_d')) else np.nan
    significant = pd.notna(row.get('v3_fdr')) and row['v3_fdr'] < FDR
    if cls in ('VERIFIED', 'PROVISIONAL'):
        if significant and claim is not None and new_sign == claim:
            outcome = 'CONFIRMED_CLAIMED_DIRECTION'
        elif significant and claim is not None:
            outcome = 'CONTRADICTS_CLAIMED_DIRECTION'
        elif significant:
            outcome = 'SIGNIFICANT_CLAIM_UNKNOWN'
        else:
            outcome = 'NOT_SIGNIFICANT_AFTER_REPAIR'
    elif cls == 'DESCRIPTIVE':
        outcome = 'DESCRIPTIVE_ONLY_ARM_BELOW_3'
    else:
        outcome = 'UNRESOLVED_' + cls
    return {'outcome': outcome, 'direction_flipped_by_repair': flipped, 'v3_significant': bool(significant)}


def run(repairs, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    atlas = []
    for modality, folder in repairs.items():
        frame = pd.read_csv(Path(folder) / 'atlas_repaired.csv', low_memory=False)
        frame['modality_label'] = MODALITY_LABELS[modality]
        frame['modality_key'] = modality
        frame['atlas_row_ordinal'] = np.arange(len(frame))
        atlas.append(frame)
    atlas = pd.concat(atlas, ignore_index=True)
    for column in ('effect_size', 'pvalue', 'v3_cohen_d', 'v3_fdr', 'v3_p'):
        atlas[column] = pd.to_numeric(atlas[column], errors='coerce')
    if 'v3_representative_class' not in atlas:
        atlas['v3_representative_class'] = np.nan
    effective = atlas.v3_final_class.where(atlas.v3_final_class != 'DUPLICATE', atlas.v3_representative_class)
    by_study = {k: g for k, g in atlas.groupby(['modality_label', 'gse_id'])}
    summary = {}
    for name, (path, keys) in PUBLISHED.items():
        if not path.exists():
            summary[name] = {'status': 'SOURCE_MISSING', 'path': paths.portable(path)}
            continue
        published = pd.read_csv(path, low_memory=False)
        records = []
        for i, item in published.iterrows():
            candidates = by_study.get((item.get('modality'), item['gse_id']))
            effect = pd.to_numeric(item.get(keys['effect']), errors='coerce')
            p = pd.to_numeric(item.get(keys['p']), errors='coerce') if keys['p'] else np.nan
            claim = claimed_direction(item.get(keys['claim']))
            record = {'published_file': name, 'published_index': int(i), 'modality': item.get('modality'),
                      'gse_id': item['gse_id'], 'label': item.get('condition_name', item.get('intervention')),
                      'published_effect': effect, 'published_p': p, 'claimed_direction': claim,
                      'historical_effect_sign_matches_claim': bool(claim is not None and np.sign(effect) == claim)}
            if candidates is None:
                record['match'] = 'STUDY_NOT_IN_ATLAS'
            else:
                hits = candidates[candidates.effect_size.apply(lambda e: close(e, effect))]
                if keys['p'] and len(hits) > 1:
                    hits = hits[hits.pvalue.apply(lambda q: close(q, p))]
                match = 'EXACT'
                if len(hits) == 0:
                    hits = candidates[candidates.effect_size.apply(lambda e: close(e, -effect))]
                    if keys['p']:
                        hits = hits[hits.pvalue.apply(lambda q: close(q, p))]
                    match = 'SIGN_REVERSED'
                if len(hits) == 0:
                    record['match'] = 'NO_EXACT_ROW'
                elif len(hits) > 1 and hits[['v3_cohen_d', 'v3_fdr']].assign(cls=effective[hits.index]).drop_duplicates().shape[0] > 1:
                    record['match'] = 'AMBIGUOUS'
                    record['candidates'] = int(len(hits))
                else:
                    row = hits.iloc[0]
                    record.update(match=match, atlas_modality=row['modality_key'], atlas_row_ordinal=int(row['atlas_row_ordinal']),
                                  comparison_name=row['comparison_name'], clock=row['clock'],
                                  atlas_historical_effect=row['effect_size'], atlas_historical_p=row['pvalue'],
                                  v3_final_class=row['v3_final_class'], v3_actions=row.get('v3_actions'),
                                  v3_reasons=row.get('v3_reasons'), v3_pending=row.get('v3_pending'),
                                  v3_soft_flags=row.get('v3_soft_flags'), v3_cohen_d=row['v3_cohen_d'], v3_p=row['v3_p'],
                                  v3_fdr=row['v3_fdr'], v3_n_case=row.get('v3_n_case'), v3_n_control=row.get('v3_n_control'),
                                  **fate(row, claim))
            records.append(record)
        frame = pd.DataFrame(records)
        frame.to_csv(out_dir / f'{name}_v3_status.csv', index=False)
        matched = frame[frame.match.isin(['EXACT', 'SIGN_REVERSED'])]
        summary[name] = {'items': len(frame), 'match': frame.match.value_counts().to_dict(),
                         'outcomes': matched.outcome.value_counts().to_dict() if len(matched) else {},
                         'direction_flipped_by_repair': int(matched.direction_flipped_by_repair.sum()) if len(matched) else 0,
                         'historical_effect_sign_contradicts_claim': int((~frame.historical_effect_sign_matches_claim).sum()),
                         'source': paths.portable(path), 'source_sha256': sha256(path)}
    summary['created'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    summary['repairs'] = {k: {'dir': paths.portable(v), 'atlas_sha256': sha256(Path(v) / 'atlas_repaired.csv')} for k, v in repairs.items()}
    try:
        (out_dir / 'summary.json').write_text(json.dumps(summary, allow_nan=False, indent=2, default=str))
    except (OSError, ValueError, TypeError) as exc:
        raise RuntimeError(f'Cannot write published status summary {out_dir}: {exc}') from exc
    print(json.dumps({k: v for k, v in summary.items() if k not in ('repairs', 'created')}, indent=1, default=str))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--repairs', action='append', required=True, help='modality=DIR')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    run(dict(item.split('=', 1) for item in args.repairs), args.out)


if __name__ == '__main__':
    main()
