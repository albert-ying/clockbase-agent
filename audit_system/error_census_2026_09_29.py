"""Census of historical defects over every atlas row (no sampling), from verifier-v3 repairs.

Each row is classed by the defects the verifier found in the historical record. Counts are over all
rows; rates use two denominators: all rows, and rows whose historical analysis could be reproduced
(membership recovered from the agent's own replayed test or by exact numerical reproduction), where
material defects are directly observable. Arm inversions count only when GPT-6 Astra and Claude
Opus 5 independently confirmed them; that dual-model agreement is the reference for orientation.

  python -m audit_system.error_census_2026_09_29 --repairs dnam=DIR --repairs hs_rna=DIR --repairs mm_rna=DIR \
      --published DIR --out DIR
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system.membership_recovery_2026_09_27 import CORRECTIONS, MODALITIES, sha256

# Defects that make the historical value wrong or unsupported as reported.
MATERIAL = {
    'CLOCK_OUTPUTS_OUTSIDE_PLAUSIBLE_RANGE': 'clock outputs outside the plausible range',
    'SINGLE_CELLS_AS_INDEPENDENT_SAMPLES': 'single cells tested as independent samples without individual aggregation',
    'ARMS_DEFINED_BY_CLOCK_OUTPUT': 'arms defined by a clock output',
    'ARMS_DO_NOT_MATCH_GEO_SAMPLE_RECORDS': 'arms do not match the GEO sample records',
    'ARMS_INVERTED_CONFIRMED': 'case and control arms swapped relative to the label (dual-model confirmed)',
    'EFFECT_NOT_REPRODUCED': 'stored effect not reproduced by any standard definition on the recovered arms',
    'NON_STANDARD_ASSAY_IN_ARMS': 'enrichment or oxidative-assay samples analysed as bisulfite beta values',
    'EMPTY_OR_ARTIFACT_INPUT': 'clock computed from an empty DNAm matrix or artefactual RNA input',
    'DETACHED_INPUT_PROFILES': 'clock computed from DNAm samples whose values are detached from their probe identities',
    'FDR_BELOW_RAW_P': 'stored FDR smaller than the raw P of the same row',
    'TEST_NOT_A_GROUP_CONTRAST': 'one-sample test or continuous term reported as a group comparison',
    'STORED_P_UNUSABLE': 'stored P value missing, 0 or 1',
}
# Design limitation: computable but underpowered; the declared policy reports these descriptively.
DESIGN = {
    'ARM_BELOW_3_WITH_INFERENCE': 'P value reported although an arm had fewer than 3 samples',
}
# Record-keeping defects: the values may be right, but the record is internally inconsistent or unattributable.
RECORD = {
    'DIRECTION_FIELD_CONTRADICTS_EFFECT': 'direction field and effect sign disagree within the row',
    'OUTPUT_MISFILED': "producing trace stored under another study's key",
    'HISTORICAL_ANALYSIS_NOT_REPRODUCIBLE': 'membership could be recovered neither from the trace nor numerically',
}


HISTORICAL_STATUS_DEFECTS = {'STORED_P_UNUSABLE': 'STORED_P_UNUSABLE', 'CONTINUOUS_TERM': 'TEST_NOT_A_GROUP_CONTRAST',
                             'ONE_SAMPLE_TEST': 'TEST_NOT_A_GROUP_CONTRAST'}


def classify(repair, row, member=None):
    """Defects of the HISTORICAL record. Rows whose arms were newly reconstructed are not reproductions."""
    reasons, actions = set(repair['reasons']), {a['action'] for a in repair['actions']}
    defects = set()
    designed = (member or {}).get('source') == 'design_reconstruction' or 'DESIGN_RECONSTRUCTED' in reasons
    if designed and (member or {}).get('historical_status') in HISTORICAL_STATUS_DEFECTS:
        defects.add(HISTORICAL_STATUS_DEFECTS[member['historical_status']])
    if 'SWAP_ARMS_CONFIRMED_INVERSION' in actions:
        defects.add('ARMS_INVERTED_CONFIRMED')
    if 'EFFECT_NOT_REPRODUCED' in reasons:
        defects.add('EFFECT_NOT_REPRODUCED')
    if 'ARM_TOO_SMALL' in reasons and pd.notna(row.get('pvalue')):
        defects.add('ARM_BELOW_3_WITH_INFERENCE')
    if 'NON_STANDARD_ASSAY' in reasons:
        defects.add('NON_STANDARD_ASSAY_IN_ARMS')
    if reasons & {'INPUT_EMPTY', 'RNA_INPUT_ARTIFACT'}:
        defects.add('EMPTY_OR_ARTIFACT_INPUT')
    if 'INPUT_PROFILE_DETACHED' in reasons:
        defects.add('DETACHED_INPUT_PROFILES')
    if 'FDR_BELOW_RAW_P' in reasons:
        defects.add('FDR_BELOW_RAW_P')
    if reasons & {'ONE_SAMPLE_TEST_NOT_CONTRAST', 'CONTINUOUS_TERM_NOT_GROUP_CONTRAST'}:
        defects.add('TEST_NOT_A_GROUP_CONTRAST')
    defects.update(reasons & {'ARMS_DEFINED_BY_CLOCK_OUTPUT', 'ARMS_DO_NOT_MATCH_GEO_SAMPLE_RECORDS'})
    defects.update(reasons & {'CLOCK_OUTPUTS_OUTSIDE_PLAUSIBLE_RANGE', 'SINGLE_CELLS_AS_INDEPENDENT_SAMPLES'})
    if 'STORED_P_UNUSABLE' in reasons:
        defects.add('STORED_P_UNUSABLE')
    effect, direction = pd.to_numeric(row.get('effect_size'), errors='coerce'), pd.to_numeric(row.get('direction'), errors='coerce')
    if pd.notna(effect) and pd.notna(direction) and direction != 0 and effect != 0 and np.sign(effect) != np.sign(direction):
        defects.add('DIRECTION_FIELD_CONTRADICTS_EFFECT')
    if 'OUTPUT_MISFILED' in reasons:
        defects.add('OUTPUT_MISFILED')
    reproduced = 'MEMBERSHIP_RECOVERED' in reasons and not designed
    if designed:
        # Arm-level findings on new arms describe the repair, not the historical analysis.
        defects -= {'ARM_BELOW_3_WITH_INFERENCE', 'NON_STANDARD_ASSAY_IN_ARMS'}
    if not reproduced:
        defects.add('HISTORICAL_ANALYSIS_NOT_REPRODUCIBLE')
    return defects, reproduced


def rate(numerator, denominator):
    return {'n': int(numerator), 'of': int(denominator), 'fraction': round(numerator / denominator, 4) if denominator else None}


def run(repairs, published_dir, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    rows_out, summary = [], {'material_definitions': MATERIAL, 'design_definitions': DESIGN,
                             'record_definitions': RECORD, 'modalities': {}}
    for modality, folder in repairs.items():
        table = CORRECTIONS / f"2026-09-08-{MODALITIES[modality]['table']}-atlas_WORKING_NOT_VALIDATED.csv"
        atlas = pd.read_csv(table, low_memory=False)
        records = [json.loads(line) for line in open(Path(folder) / 'repairs.jsonl')]
        merged_path = Path(folder).parent / f'merged-{modality}' / 'rows.jsonl'
        members = {r['row_ordinal']: r for r in map(json.loads, open(merged_path))} if merged_path.exists() else {}
        per = {k: 0 for k in list(MATERIAL) + list(DESIGN) + list(RECORD)}
        reproduced_n = material_any = material_any_reproduced = 0
        for repair, row in zip(records, atlas.to_dict('records')):
            defects, reproduced = classify(repair, row, members.get(repair['row_ordinal']))
            for d in defects:
                per[d] += 1
            material = bool(defects & set(MATERIAL))
            reproduced_n += reproduced
            material_any += material
            material_any_reproduced += material and reproduced
            rows_out.append({'modality': modality, 'row_ordinal': repair['row_ordinal'], 'gse_id': repair['gse_id'],
                             'clock': repair['clock'], 'reproduced': reproduced, 'material_defect': material,
                             'defects': ';'.join(sorted(defects)), 'final_class': repair['final_class']})
        n = len(records)
        summary['modalities'][modality] = {
            'rows': n, 'historical_analysis_reproduced': rate(reproduced_n, n),
            'any_material_defect_all_rows': rate(material_any, n),
            'any_material_defect_among_reproduced': rate(material_any_reproduced, reproduced_n),
            'defects': {k: rate(v, n) for k, v in per.items()},
            'table_sha256': sha256(table), 'repairs_sha256': sha256(Path(folder) / 'repairs.jsonl')}
    frame = pd.DataFrame(rows_out)
    frame.to_csv(out_dir / 'row_defects.csv', index=False)
    total = len(frame)
    reproduced = frame.reproduced.sum()
    summary['all'] = {
        'rows': total, 'historical_analysis_reproduced': rate(reproduced, total),
        'any_material_defect_all_rows': rate(frame.material_defect.sum(), total),
        'any_material_defect_among_reproduced': rate((frame.material_defect & frame.reproduced).sum(), reproduced),
        'defects': {k: rate(frame.defects.str.contains(k).sum(), total) for k in list(MATERIAL) + list(DESIGN) + list(RECORD)}}
    if published_dir:
        published = {}
        for path in sorted(Path(published_dir).glob('*_v3_status.csv')):
            p = pd.read_csv(path)
            matched = p[p.match == 'EXACT']
            keyed = frame.set_index(['modality', 'row_ordinal'])
            flags = [keyed.loc[(m, int(o))] for m, o in zip(matched.atlas_modality, matched.atlas_row_ordinal)]
            defects = pd.DataFrame(flags)
            published[path.name.replace('_v3_status.csv', '')] = {
                'items_matched': int(len(matched)),
                'any_material_defect': rate(defects.material_defect.sum(), len(defects)) if len(defects) else None,
                'defects': {k: int(defects.defects.str.contains(k).sum()) for k in list(MATERIAL) + list(DESIGN) + list(RECORD)} if len(defects) else {},
                'outcomes': matched.outcome.value_counts().to_dict()}
        summary['published'] = published
    summary['reference_note'] = ('Arm inversions are counted only when two independent models agreed; other defects are '
                                 'deterministic checks. Counts are a census of all rows, not an estimate from a sample.')
    summary['created'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    summary['module_sha256'] = sha256(__file__)
    (out_dir / 'census.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps({k: summary[k] for k in ('all',)}, indent=1, default=str))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--repairs', action='append', required=True, help='modality=DIR')
    parser.add_argument('--published', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    run(dict(item.split('=', 1) for item in args.repairs), args.published, args.out)


if __name__ == '__main__':
    main()
