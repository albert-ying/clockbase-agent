"""Export the atlas evidence for the interventions of the experimental follow-up: PLX5622 and other CSF1R inhibitors,
sleep deprivation and cold exposure.

For each group it writes, from one verifier run and its intervention meta-analysis:
  effects.csv          every atlas comparison x clock with the intervention (all sample contexts; verification class)
  ranked.csv           the group's comparisons in the ranking, with both models' judgement answers
  meta_pooled.csv      pooled estimates of the primary meta-analysis and each sensitivity run
  meta_unit_effects.csv  the per-project effects behind the primary pooled estimates
  records/<modality>/<GSE>.json  the verification record of every series involved
  SUMMARY.md           counts and pooled estimates computed from these files
and provenance.json with the inputs' SHA-256.

  python scripts/export_focus_interventions_2026_11_02.py --out DIR \
      [--final final-2026-10-10] [--meta intervention-meta-2026-10-17] [--judgement followup-judgment-001]
"""
import argparse
import hashlib
import json
import shutil
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from audit_system import paths  # noqa: E402
from audit_system.intervention_meta_2026_09_29 import usable_effects  # noqa: E402

# Canonical targets (two-model intervention classification) of each group; names the classification gives these
# interventions in the atlas, chosen by name, not by result.
GROUPS = {
    'plx5622_csf1r': ['plx5622', 'pexidartinib', 'plx3397', 'blz945', 'csf1r inhibitor', 'anti csf1r antibody', 'csf1r'],
    'sleep_deprivation': ['sleep deprivation'],
    'cold_exposure': ['cold exposure', 'cold water immersion'],
}
KEY = ['modality', 'gse_id', 'comparison_name', 'condition_name']
MODALITY_KEY = {'human DNAm': 'dnam', 'human RNA-seq': 'hs_rna', 'mouse RNA-seq': 'mm_rna'}
EFFECT_COLUMNS = ['modality', 'gse_id', 'comparison_name', 'condition_name', 'condition_category', 'canonical_target',
                  'intervention_class', 'manipulation', 'clock', 'v3_context', 'v3_final_class', 'v3_n_case',
                  'v3_n_control', 'v3_mean_diff', 'v3_cohen_d', 'v3_p', 'v3_fdr', 'v3_reasons', 'v3_soft_flags']


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_csv(path):
    try:
        return pd.read_csv(path, low_memory=False)
    except OSError as error:
        raise SystemExit(f'Cannot proceed - required data file not found: {path} ({error})')


def summary(name, effects, ranked, pooled):
    lines = [f'# {name}', '', f'Canonical targets: {", ".join(GROUPS[name])}', '']
    usable = effects[effects.v3_final_class.isin(['VERIFIED', 'PROVISIONAL'])]
    lines.append(f'- Atlas comparison x clock rows: {len(effects)} ({len(usable)} verified or provisional) from '
                 f'{effects.gse_id.nunique()} series; in vivo rows: {int((effects.v3_context == "in_vivo").sum())}.')
    lines.append(f'- Ranked comparisons: {len(ranked)} (with lower clock outputs: '
                 f'{int((ranked.direction_v3 == "age-lowering").sum())}; kept in the judgement version: '
                 f'{int(ranked.judgement_admitted.astype(bool).sum()) if len(ranked) else 0}).')
    primary = pooled[pooled.run == 'primary']
    for row in primary.sort_values('k', ascending=False).itertuples():
        lines.append(f'- Pooled ({row.canonical_target}, {row.condition_category}, primary): {row.k} projects, '
                     f"Hedges' g = {row.estimate:.2f} (95% CI {row.ci_low:.2f} to {row.ci_high:.2f}), FDR = {row.fdr:.3g}.")
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--final', default='final-2026-10-10')
    parser.add_argument('--meta', default='intervention-meta-2026-10-17')
    parser.add_argument('--judgement', default='followup-judgment-001')
    args = parser.parse_args()
    root = paths.VERIFIER_ROOT
    final, meta, judgement = root / args.final, root / args.meta, root / args.judgement
    atlas = final / 'atlas-v3' / 'atlas_v3_combined.csv'
    effects = usable_effects(atlas, root / 'classification-001', None)
    ranked = read_csv(judgement / 'ranked_with_judgement.csv')
    runs = sorted(p.parent.name for p in meta.glob('*/intervention_meta.csv'))
    pooled_all = pd.concat([read_csv(meta / run / 'intervention_meta.csv').assign(run=run) for run in runs])
    units = read_csv(meta / 'primary' / 'unit_effects.csv')
    args.out.mkdir(parents=True, exist_ok=True)
    report = {}
    for name, targets in GROUPS.items():
        folder = args.out / name
        folder.mkdir(exist_ok=True)
        chosen = effects[effects.canonical_target.isin(targets)]
        chosen[[c for c in EFFECT_COLUMNS if c in chosen.columns]].sort_values(['canonical_target', 'gse_id']).to_csv(
            folder / 'effects.csv', index=False)
        keys = chosen[KEY].drop_duplicates()
        in_ranking = ranked.merge(keys, on=KEY, how='inner')
        in_ranking.to_csv(folder / 'ranked.csv', index=False)
        pooled = pooled_all[pooled_all.canonical_target.isin(targets)]
        pooled.to_csv(folder / 'meta_pooled.csv', index=False)
        units[units.canonical_target.isin(targets)].to_csv(folder / 'meta_unit_effects.csv', index=False)
        for modality, gse in chosen[['modality', 'gse_id']].drop_duplicates().itertuples(index=False):
            source = final / 'public-export' / 'studies' / MODALITY_KEY[modality] / f'{gse}.json'
            if source.is_file():
                target = folder / 'records' / MODALITY_KEY[modality]
                target.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, target / source.name)
        (folder / 'SUMMARY.md').write_text(summary(name, chosen, in_ranking, pooled))
        report[name] = {'effects': len(chosen), 'series': int(chosen.gse_id.nunique()), 'ranked': len(in_ranking),
                        'pooled_rows': len(pooled)}
    provenance = {'final': args.final, 'meta': args.meta, 'judgement': args.judgement, 'groups': GROUPS,
                  'inputs_sha256': {'atlas_v3_combined.csv': sha256(atlas),
                                    'ranked_with_judgement.csv': sha256(judgement / 'ranked_with_judgement.csv'),
                                    'primary/intervention_meta.csv': sha256(meta / 'primary' / 'intervention_meta.csv')},
                  'meta_runs': runs, 'counts': report, 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (args.out / 'provenance.json').write_text(json.dumps(provenance, indent=1) + '\n')
    print(json.dumps(report, indent=1))


if __name__ == '__main__':
    main()
