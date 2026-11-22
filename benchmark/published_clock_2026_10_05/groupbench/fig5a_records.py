"""Benchmark B ground truth from Tyshkovskiy et al. 2026 Fig. 5a (mouse disease models, mortality tAge).

Inputs: two independent readings of the figure image (b_rna/fig5a_reads.json; the paper has no source-data file). The
readers agreed on series, labels, printed adjusted P, direction and significance for 14 of 15 comparisons; the one
difference is a non-significant rat comparison outside the atlas. Records are written for the comparisons on atlas
series. The definitions carry the dataset total printed in the figure. Where the two-model sample assignment pools
genotypes and so contradicts the printed total (GSE120977, GSE197699), the published groups are not identifiable:
those records are excluded from the primary analysis and, for a sensitivity analysis only, given wild-type groups
(the disease-model background; the genotype is not stated in the paper).

  python -m benchmark.published_clock_2026_10_05.groupbench.fig5a_records --dir B_RNA
"""
import argparse
import json
from pathlib import Path

import pandas as pd

from audit_system import paths

ATLAS = {'GSE142633', 'GSE129586', 'GSE137482', 'GSE189377', 'GSE120977', 'GSE197699'}
MODEL = {'AD': "Alzheimer's disease model (5xFAD)", 'IS': 'ischaemic stroke model', 'CKD': 'chronic kidney disease model',
         'NASH': 'nonalcoholic steatohepatitis model', 'DN': 'diabetic nephropathy model'}
SENSITIVITY_NOTE = ('SENSITIVITY: Fig. 5a prints the dataset total ({n}), which equals one genotype; the genotype is not '
                    'stated, and wild type (the disease-model background) is assumed. Not used in the primary analysis.')


def records(reads):
    reader = reads[0]
    totals = {c['series']: c['n_shown'] for c in reader['comparisons']}
    out = []
    for i, c in enumerate(reader['comparisons']):
        if c['series'] not in ATLAS:
            continue
        age = '' if 'not printed' in c['time_point_or_age'] else c['time_point_or_age']
        where = f"{c['tissue']}{', ' + age if age else ''}"
        suffix = f"; Fig. 5a prints {totals[c['series']]} animals for this dataset (all groups shown)"
        out.append({
            'id': f"TYSH2026-FIG5A-{c['series']}-{i}-v2", 'tier': 'T1', 'gse': c['series'], 'modality': 'mouse RNA-seq',
            'tissue': c['tissue'], 'clock_reported': 'BR rodent multi-tissue mortality clock (Tyshkovskiy et al. 2026)',
            'clock_clockbase': 'Hazard_Mouse_All_All_scaleddiff_EN_tAge', 'measure': 'age',
            'case_label': c['case_label'], 'control_label': 'Control',
            'case_definition': f"{MODEL.get(c['case_label'].split()[0], c['disease_model'])} animals ({where}), as labelled '{c['case_label']}' in Fig. 5a{suffix}",
            'control_definition': f"age-matched control animals ({where}) of the same dataset, labelled 'Control' in Fig. 5a{suffix}",
            'n_case': None, 'n_control': None, 'direction': c['direction'], 'significant': c['significant'],
            'threshold': 'Benjamini-Hochberg-adjusted P < 0.05', 'p_value': None, 'effect': None,
            'test': 'mixed-effects model, Benjamini-Hochberg-adjusted', 'covariates': [],
            'source': {'paper_id': 'DOI:10.1038/s41586-026-10542-3', 'doi': '10.1038/s41586-026-10542-3', 'pmcid': 'PMC13233323',
                       'citation': 'Tyshkovskiy et al., Universal transcriptomic hallmarks of mammalian ageing and mortality, Nature 654:173-188 (2026)',
                       'year': '2026', 'location': 'Fig. 5a',
                       'quote': f"{reader['legend_quote']} [Fig. 5a, read from the image by two independent readers: {c['series']} "
                                f"{c['case_label']} vs Control{', ' + age if age else ''}, BH-adjusted P = {c['p_adj_printed']}, {c['n_shown']}]"},
            'extraction': {'extractors': ['figure-reader-A', 'figure-reader-B (sonnet)'], 'agreement': 'agreed',
                           'notes': ['Values read from the figure image; the paper has no source-data file.',
                                     'Clock differs from ClockBase: BR rodent mortality clock vs ClockBase EN mouse mortality clock (version difference).',
                                     f"Printed adjusted P {c['p_adj_printed']} kept in the quote; significance from it."]}})
    return out


def sensitivity_groups():
    meta = paths.GEO / 'GEO_mm_rna/processed_metadata_preharmonize'
    a = pd.read_csv(meta / 'GSE120977.csv', dtype=str)
    b = pd.read_csv(meta / 'GSE197699.csv', dtype=str)
    return {
        'TYSH2026-FIG5A-GSE120977-13-v2': {
            'case_gsms': a[(a['genotype/variation'] == 'Wild Type') & (a.diet == 'CDAHFD')].GSM.tolist(),
            'control_gsms': a[(a['genotype/variation'] == 'Wild Type') & (a.diet == 'Normal Diet')].GSM.tolist(),
            'provenance': SENSITIVITY_NOTE.format(n='n = 10')},
        'TYSH2026-FIG5A-GSE197699-14-v2': {
            'case_gsms': b[b.title.str.contains('WT diabetic rep')].GSM.tolist(),
            'control_gsms': b[b.title.str.contains('WT non-diabetic')].GSM.tolist(),
            'provenance': SENSITIVITY_NOTE.format(n='n = 6')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--dir', type=Path, required=True)
    args = parser.parse_args()
    reads = json.loads((args.dir / 'fig5a_reads.json').read_text())
    rows = records(reads)
    target = args.dir / 'ground_truth_fig5a_v2_2026-10-06.jsonl'
    rebuilt = ''.join(json.dumps(r) + '\n' for r in rows)
    if target.exists() and target.read_text() != rebuilt:
        raise SystemExit(f'{target} differs from the rebuild')
    target.write_text(rebuilt)
    print(len(rows), 'records;', {k: (len(v['case_gsms']), len(v['control_gsms'])) for k, v in sensitivity_groups().items()})


if __name__ == '__main__':
    main()
