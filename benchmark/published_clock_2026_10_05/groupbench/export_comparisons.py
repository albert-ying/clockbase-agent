"""Collapse the verifier's public-export rows into agent comparisons (one per series, names and arms), all clocks.

Each comparison keeps the membership source of its arms (replay, recovery, both, design_reconstruction or none), the
comparison and condition names, and per clock the row ordinal, final class, soft flags, the agent's stored result and
the verifier's recomputation. Rows without a unique two-sample membership are kept with empty arms so that their
numbers can be reported.

  python -m benchmark.published_clock_2026_10_05.groupbench.export_comparisons --out FILE.jsonl
"""
import argparse
import glob
import hashlib
import json
import os
from collections import Counter
from pathlib import Path

from audit_system import paths

EXPORT = paths.VERIFIER_ROOT / os.environ.get('CLOCKBASE_FINAL_NAME', 'final-2026-10-05') / 'public-export' / 'studies'
MODALITIES = {'mm_rna': 'mouse RNA-seq', 'hs_rna': 'human RNA-seq', 'dnam': 'human DNAm'}
AGENT_SOURCES = {'replay', 'recovery', 'replay+recovery'}


def arm(value):
    return sorted(value) if isinstance(value, list) else []


def comparison_key(gse, row):
    membership = row.get('membership') or {}
    case, control = arm(membership.get('case_gsms')), arm(membership.get('control_gsms'))
    text = json.dumps([gse, row.get('comparison_name'), row.get('condition_name'), case, control])
    return hashlib.sha256(text.encode()).hexdigest()[:16], case, control


def export(out):
    comparisons, counts = {}, Counter()
    for key, modality in MODALITIES.items():
        for path in sorted(glob.glob(str(EXPORT / key / '*.json'))):
            study = json.loads(Path(path).read_text())
            gse = study['gse_id']
            for row in study['rows']:
                identity, case, control = comparison_key(gse, row)
                membership = row.get('membership') or {}
                item = comparisons.setdefault(identity, {
                    'comparison_id': f'{gse}:{identity}', 'modality': modality, 'gse': gse,
                    'comparison_name': row.get('comparison_name'), 'condition_name': row.get('condition_name'),
                    'membership_status': membership.get('status'), 'membership_source': row.get('membership_source'),
                    'agent_arms': row.get('membership_source') in AGENT_SOURCES and membership.get('status') == 'UNIQUE',
                    'case_gsms': case, 'control_gsms': control, 'clocks': []})
                checks = row.get('checks') or {}
                item['clocks'].append({
                    'clock': row.get('clock'), 'row_ordinal': row.get('row_ordinal'), 'final_class': row.get('final_class'),
                    'soft_flags': row.get('soft_flags'), 'historical': row.get('historical'),
                    'recomputed': row.get('recomputed'), 'duplicate_of': row.get('duplicate_of'),
                    'effect_reproduction': (checks.get('effect_reproduction') or {}).get('outcome'),
                    'arm_orientation': (checks.get('arm_orientation') or {}).get('outcome')})
                counts[(modality, row.get('membership_source'))] += 1
    with open(out, 'w') as handle:
        for item in comparisons.values():
            handle.write(json.dumps(item) + '\n')
    return {'comparisons': len(comparisons),
            'agent_arm_comparisons': sum(c['agent_arms'] for c in comparisons.values()),
            'rows_by_modality_and_source': {f'{m}|{s}': n for (m, s), n in sorted(counts.items(), key=str)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f'refusing to overwrite {args.out}')
    print(json.dumps(export(args.out), indent=1))


if __name__ == '__main__':
    main()
