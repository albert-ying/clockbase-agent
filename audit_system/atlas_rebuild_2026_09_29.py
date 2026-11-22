"""Rebuild the ranked atlas outputs from verifier-v3 repaired rows with one statistical ranking.

Every retained comparison uses the declared estimand (Welch test; Cohen's d with pooled SD, case minus
control). A comparison is a robust candidate when its final class is VERIFIED or PROVISIONAL, its
family-wise BH FDR is below 0.05 and the 95% confidence interval of d (Hedges-Olkin standard error)
excludes zero. Candidates are ranked by the confidence bound nearest zero, |bound|, so ranking rewards
effects that are both large and precisely estimated. No author-defined score enters the ranking.
All clocks here increase with age or risk, so d < 0 is age-lowering (anti-aging) and d > 0 age-raising.

Ranked outputs admit only intervention categories (Genetic, Drug, Disease, Environment) whose samples are
mostly in vivo (per-sample vivo_vitro annotations, >= 50% of annotated members). In vitro comparisons are
ranked separately: replicate cultures have very small within-group variance, so their standardized
effects are not comparable with in vivo effects.

Outputs (all rows keep their v3 class, soft flags and repair actions):
  atlas_v3_combined.csv       every atlas row with v3 values, standard error, CI and ranking bound
  category_top15_v3.csv       Figure 3 replacement: top 15 per category (Genetic, Drug, Disease, Environment)
  table1_pro_aging_v3.csv     top 30 age-raising comparisons, pooled across modalities
  table2_anti_aging_v3.csv    top 30 age-lowering comparisons, pooled across modalities
  s14_v3_candidates.csv       all robust age-lowering comparisons (Supplementary Table S14 replacement)
  *_in_vitro_v3.csv           the same rankings for in vitro comparisons

  python -m audit_system.atlas_rebuild_2026_09_29 --repairs dnam=DIR --repairs hs_rna=DIR \
      --repairs mm_rna=DIR --out DIR
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from audit_system.membership_recovery_2026_09_27 import GEO, MODALITIES, sha256

MODALITY_LABELS = {'dnam': 'human DNAm', 'hs_rna': 'human RNA-seq', 'mm_rna': 'mouse RNA-seq'}
CATEGORIES = ['Genetic', 'Drug', 'Disease', 'Environment']
FDR, Z = 0.05, 1.959963984540054


def hedges_olkin_se(d, n1, n2):
    return np.sqrt((n1 + n2) / (n1 * n2) + d ** 2 / (2 * (n1 + n2)))


def annotate(frame):
    d = pd.to_numeric(frame.v3_cohen_d, errors='coerce')
    n1, n2 = pd.to_numeric(frame.v3_n_case, errors='coerce'), pd.to_numeric(frame.v3_n_control, errors='coerce')
    se = hedges_olkin_se(d, n1, n2)
    frame['v3_d_se'] = se
    frame['v3_d_ci_low'], frame['v3_d_ci_high'] = d - Z * se, d + Z * se
    frame['v3_direction'] = np.sign(d)
    usable = frame.v3_final_class.isin(['VERIFIED', 'PROVISIONAL'])
    significant = pd.to_numeric(frame.v3_fdr, errors='coerce') < FDR
    excludes_zero = (frame.v3_d_ci_low > 0) | (frame.v3_d_ci_high < 0)
    frame['v3_robust'] = usable & significant & excludes_zero
    frame['v3_rank_bound'] = np.where(d < 0, -frame.v3_d_ci_high, frame.v3_d_ci_low)  # |bound nearest zero|
    frame['v3_label'] = frame['intervention'].fillna(frame['condition_name']).astype(str).str.slice(0, 80)
    return frame


def sample_context(modality, repairs_path):
    """Per atlas row: share of annotated members that are in vivo, from the GPT-4 per-sample annotations."""
    annotations = pd.read_csv(GEO / 'result' / MODALITIES[modality]['gpt4'], usecols=['GSM_ID', 'vivo_vitro'], dtype=str)
    vivo = annotations.assign(v=annotations.vivo_vitro.str.lower()).drop_duplicates('GSM_ID').set_index('GSM_ID').v
    shares = []
    for line in open(repairs_path):
        record = json.loads(line)
        inputs = record.get('inputs') or {}
        members = list(inputs.get('case_gsms') or []) + list(inputs.get('control_gsms') or [])
        values = vivo.reindex(members).dropna()
        values = values[values.isin(['in_vivo', 'in_vitro'])]
        shares.append(float((values == 'in_vivo').mean()) if len(values) else np.nan)
    return shares


def top(frame, n, keys=('gse_id', 'v3_label')):
    """Best comparison per study and intervention label, then the n highest ranking bounds."""
    ordered = frame.sort_values('v3_rank_bound', ascending=False)
    return ordered.drop_duplicates(list(keys)).head(n)


COLUMNS = ['modality', 'gse_id', 'comparison_name', 'condition_name', 'v3_label', 'condition_category', 'clock',
           'v3_final_class', 'v3_soft_flags', 'v3_actions', 'v3_n_case', 'v3_n_control', 'v3_mean_diff', 'v3_cohen_d',
           'v3_d_ci_low', 'v3_d_ci_high', 'v3_p', 'v3_fdr', 'v3_rank_bound', 'v3_membership_source',
           'effect_size', 'pvalue', 'fdr', 'direction', 'final_score']


def run(repairs, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    frames = []
    for modality, folder in repairs.items():
        frame = pd.read_csv(Path(folder) / 'atlas_repaired.csv', low_memory=False)
        frame['modality'] = MODALITY_LABELS[modality]
        frame['atlas_row_ordinal'] = np.arange(len(frame))
        frame['v3_in_vivo_share'] = sample_context(modality, Path(folder) / 'repairs.jsonl')
        frames.append(frame)
    atlas = annotate(pd.concat(frames, ignore_index=True))
    atlas['v3_context'] = np.where(atlas.v3_in_vivo_share.isna(), 'unknown',
                                   np.where(atlas.v3_in_vivo_share >= 0.5, 'in_vivo', 'in_vitro'))
    atlas['v3_rankable_category'] = atlas.condition_category.isin(CATEGORIES)
    atlas.to_csv(out_dir / 'atlas_v3_combined.csv', index=False)
    robust = atlas[atlas.v3_robust]
    columns = [c for c in COLUMNS if c in atlas.columns] + ['v3_context', 'v3_in_vivo_share']
    for context, suffix in (('in_vivo', ''), ('in_vitro', '_in_vitro')):
        pool = robust[robust.v3_rankable_category & (robust.v3_context == context)]
        categories = [top(pool[pool.condition_category == c], 15)[columns].assign(category=c) for c in CATEGORIES]
        pd.concat(categories).to_csv(out_dir / f'category_top15{suffix}_v3.csv', index=False)
        top(pool[pool.v3_direction > 0], 30)[columns].to_csv(out_dir / f'table1_pro_aging{suffix}_v3.csv', index=False)
        top(pool[pool.v3_direction < 0], 30)[columns].to_csv(out_dir / f'table2_anti_aging{suffix}_v3.csv', index=False)
        candidates = top(pool[pool.v3_direction < 0], len(pool))
        candidates[columns].to_csv(out_dir / f's14{suffix}_v3_candidates.csv', index=False)
        if context == 'in_vivo':
            s14 = candidates
    summary = {
        'rows': len(atlas),
        'final_classes': atlas.v3_final_class.value_counts().to_dict(),
        'robust_comparisons': int(len(robust)),
        'robust_by_modality_direction': {f'{m} {"age-raising" if s > 0 else "age-lowering"}': int(n)
                                         for (m, s), n in robust.groupby(['modality', 'v3_direction']).size().items()},
        'robust_by_category': robust.condition_category.value_counts().to_dict(),
        'robust_by_context': robust.v3_context.value_counts().to_dict(),
        'ranked_pool_in_vivo': int((robust.v3_rankable_category & (robust.v3_context == 'in_vivo')).sum()),
        'ranked_pool_in_vitro': int((robust.v3_rankable_category & (robust.v3_context == 'in_vitro')).sum()),
        'robust_provisional_share': round(float((robust.v3_final_class == 'PROVISIONAL').mean()), 4) if len(robust) else None,
        's14_v3_candidates_deduplicated': int(len(s14)),
        'ranking': 'confidence bound of Cohen d nearest zero (Hedges-Olkin SE, 95%); deduplicated by study and intervention label',
        'inclusion': ('VERIFIED or PROVISIONAL; family BH FDR < 0.05; 95% CI of d excludes 0; ranked outputs: '
                      'category in Genetic/Drug/Disease/Environment and >= 50% of annotated samples in vivo'),
        'repairs': {k: sha256(Path(v) / 'atlas_repaired.csv') for k, v in repairs.items()},
        'module_sha256': sha256(__file__), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=1, default=str))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--repairs', action='append', required=True, help='modality=DIR')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args(argv)
    run(dict(item.split('=', 1) for item in args.repairs), args.out)


if __name__ == '__main__':
    main()
