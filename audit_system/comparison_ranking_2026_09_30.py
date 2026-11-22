"""Comparison-level ranking of verifier-v3 atlas rows (replaces the author-defined composite score).

Studies differ in tissue, dose, duration and model, so comparisons are ranked individually rather than
pooled per intervention. A comparison is eligible when:
  - its final class is VERIFIED or PROVISIONAL, its family BH FDR is below 0.05 and each arm has at
    least 3 samples;
  - both classification models give the same contrast type, and that type matches the displayed
    category: INTERVENTION for Genetic, Drug and Environment, DISEASE_VS_HEALTHY for Disease. Tissue,
    cell-type, developmental and cell-state contrasts therefore never enter the ranking. The models may
    still differ on intervention class, direction or target name (annotations, not eligibility);
  - its clock output is not known to be miscomputed: GrimAgeV2 is excluded because the historical recipe
    and the stream repairs supplied the Horvath prediction as its chronological-age input;
  - no design concern agreed by both label-adjudicating models is attached (ADJUDICATED_* flags);
  - the dual-model design audit (design_audit_2026_09_30) admits it: both models read the samples as
    profiled directly from organisms, neither arm as tumour material, and no confound. Without an audit
    directory the GPT-4 per-sample in vivo label is used instead and the summary says so; that label
    calls cultured primary cells in vivo.
Eligible comparisons are ranked by the bound of the exact confidence interval of Hedges' g (noncentral
t) nearest zero, so small-sample extremes are not favoured. Each displayed item is the best of the k
usable rows that share its study and label (up to five clocks for human DNAm, and repeated tissues or
time points), so its interval is Bonferroni-adjusted for k (level 1 - 0.05/k). The dual-model relevance
score, produced without results, is reported beside each comparison; a sensitivity version keeps
intervention comparisons with relevance >= 5 (the rubric scores non-interventions 0-2, so Disease is
not filtered).

  python -m audit_system.comparison_ranking_2026_09_30 --atlas atlas_v3_combined.csv --classification DIR \
      --out DIR [--design-audit DIR]
"""
import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

from audit_system.membership_recovery_2026_09_27 import sha256

FDR, MIN_ARM, RELEVANCE_THRESHOLD, LEVEL = 0.05, 3, 5, 0.95
CATEGORY_CONTRAST = {'Genetic': 'INTERVENTION', 'Drug': 'INTERVENTION', 'Environment': 'INTERVENTION',
                     'Disease': 'DISEASE_VS_HEALTHY'}
KEY = ['modality', 'gse_id', 'comparison_name', 'condition_name']
# Clocks whose atlas outputs are miscomputed; their rows stay in the atlas (flagged) but are never ranked.
UNRANKED_CLOCKS = {'GrimAgeV2': 'Horvathv1 prediction supplied as the chronological-age input (row_verdicts DNAM_CLOCKS)'}


def contrast_agreed(classes):
    """Comparisons whose two classifications give the same contrast type (status AGREED, or MODELS_DISAGREE on
    other fields only: the record keeps contrast_type when the two types agree)."""
    return classes[classes.status.isin(['AGREED', 'MODELS_DISAGREE']) & classes.contrast_type.notna()]


def design_audit_table(design_audit_dirs):
    """Concatenate design-audit runs (the first run, then runs restricted to comparisons it did not cover)."""
    dirs = [design_audit_dirs] if isinstance(design_audit_dirs, (str, Path)) else list(design_audit_dirs or [])
    frames = [pd.read_json(Path(d) / 'comparisons.jsonl', lines=True) for d in dirs]
    return pd.concat(frames, ignore_index=True).drop_duplicates(KEY, keep='first') if frames else None
COLUMNS = ['modality', 'gse_id', 'comparison_name', 'condition_name', 'condition_category', 'clock', 'v3_final_class',
           'v3_soft_flags', 'v3_actions', 'v3_n_case', 'v3_n_control', 'v3_mean_diff', 'hedges_g', 'g_ci_low', 'g_ci_high',
           'ci_level', 'selection_k', 'rank_bound', 'v3_p', 'v3_fdr', 'relevance_score', 'relevance_score_abs_diff',
           'canonical_target', 'design_audit', 'design_audit_reasons', 'v3_membership_source', 'v3_context',
           'effect_size', 'pvalue', 'fdr', 'direction', 'final_score']


def exact_d_ci(d, n1, n2, level=LEVEL):
    """Exact CI of the standardized mean difference by inverting the noncentral t distribution."""
    df = n1 + n2 - 2
    scale = np.sqrt(1 / n1 + 1 / n2)
    t = d / scale
    alpha = (1 - level) / 2

    def solve(target):
        f = lambda nc: stats.nct.cdf(t, df, nc) - target
        lo, hi = t - 10 - abs(t), t + 10 + abs(t)
        return optimize.brentq(f, lo, hi)

    return solve(1 - alpha) * scale, solve(alpha) * scale


def hedges_j(n1, n2):
    return 1 - 3 / (4 * (n1 + n2 - 2) - 1)


def run(atlas_path, classification_dir, out_dir, design_audit_dir=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)
    atlas = pd.read_csv(atlas_path, low_memory=False)
    classes = pd.read_json(Path(classification_dir) / 'comparisons.jsonl', lines=True)
    classes = contrast_agreed(classes)
    atlas = atlas.merge(classes[KEY + ['contrast_type', 'canonical_target', 'relevance_score', 'relevance_score_abs_diff']],
                        on=KEY, how='left')
    atlas['label'] = atlas['intervention'].fillna(atlas['condition_name']).astype(str).str.slice(0, 80)
    usable = (atlas.v3_final_class.isin(['VERIFIED', 'PROVISIONAL']) & (atlas.v3_n_case >= MIN_ARM)
              & (atlas.v3_n_control >= MIN_ARM) & ~atlas.clock.isin(UNRANKED_CLOCKS))
    selection = atlas[usable].groupby(['gse_id', 'label']).size().rename('selection_k')
    atlas = atlas.merge(selection, left_on=['gse_id', 'label'], right_index=True, how='left')
    base = usable & (pd.to_numeric(atlas.v3_fdr, errors='coerce') < FDR) & atlas.condition_category.isin(CATEGORY_CONTRAST)
    matches = atlas.contrast_type == atlas.condition_category.map(CATEGORY_CONTRAST)
    design_flagged = atlas.v3_soft_flags.fillna('').astype(str).str.contains('ADJUDICATED_')
    audit = design_audit_table(design_audit_dir)
    if audit is not None:
        audit = audit[KEY + ['status', 'reasons']].rename(columns={'status': 'design_audit', 'reasons': 'design_audit_reasons'})
        atlas = atlas.merge(audit, on=KEY, how='left')
        admitted = atlas.design_audit == 'RANKABLE'
        context_source = 'dual-model design audit (ORGANISM_DIRECT, tumour NO, no confound; both models)'
    else:
        admitted = atlas.v3_context == 'in_vivo'
        context_source = 'GPT-4 per-sample vivo_vitro label; known to call cultured primary cells in vivo'
    eligible = atlas[base & matches & ~design_flagged & admitted].copy()
    j = hedges_j(eligible.v3_n_case, eligible.v3_n_control)
    eligible['hedges_g'] = j * eligible.v3_cohen_d
    eligible['ci_level'] = 1 - (1 - LEVEL) / eligible.selection_k
    bounds = [exact_d_ci(d, a, b, level) for d, a, b, level in
              zip(eligible.v3_cohen_d, eligible.v3_n_case, eligible.v3_n_control, eligible.ci_level)]
    eligible['g_ci_low'] = j * np.array([b[0] for b in bounds])
    eligible['g_ci_high'] = j * np.array([b[1] for b in bounds])
    excludes_zero = (eligible.g_ci_low > 0) | (eligible.g_ci_high < 0)
    eligible = eligible[excludes_zero]
    eligible['rank_bound'] = np.where(eligible.hedges_g < 0, -eligible.g_ci_high, eligible.g_ci_low)
    eligible['direction_v3'] = np.where(eligible.hedges_g < 0, 'age-lowering', 'age-raising')
    ordered = eligible.sort_values('rank_bound', ascending=False).drop_duplicates(['gse_id', 'label'])
    columns = [c for c in COLUMNS if c in ordered.columns] + ['direction_v3', 'label']
    ordered[columns].to_csv(out_dir / 'ranked_comparisons_all.csv', index=False)
    outputs = {}
    relevant = (ordered.condition_category == 'Disease') | (ordered.relevance_score >= RELEVANCE_THRESHOLD)
    for suffix, pool in (('', ordered), (f'_relevance{RELEVANCE_THRESHOLD}', ordered[relevant])):
        pd.concat([pool[pool.condition_category == c].head(15).assign(category=c) for c in CATEGORY_CONTRAST])[columns + ['category']] \
            .to_csv(out_dir / f'figure3_category_top15{suffix}.csv', index=False)
        pool[pool.direction_v3 == 'age-raising'].head(30)[columns].to_csv(out_dir / f'table1_pro_aging{suffix}.csv', index=False)
        pool[pool.direction_v3 == 'age-lowering'].head(30)[columns].to_csv(out_dir / f'table2_anti_aging{suffix}.csv', index=False)
        pool[pool.direction_v3 == 'age-lowering'][columns].to_csv(out_dir / f's14_age_lowering{suffix}.csv', index=False)
        outputs[suffix or 'all'] = {'ranked': int(len(pool)), 'age_lowering': int((pool.direction_v3 == 'age-lowering').sum()),
                                    'age_raising': int((pool.direction_v3 == 'age-raising').sum()),
                                    'by_category': pool.condition_category.value_counts().to_dict(),
                                    'by_modality': pool.modality.value_counts().to_dict()}
    step = lambda mask: {'rows': int(mask.sum()), 'comparisons': int(atlas[mask].drop_duplicates(KEY).shape[0])}
    summary = {'context_source': context_source, 'unranked_clocks': UNRANKED_CLOCKS,
               'unranked_clock_rows_verified_or_provisional': int((atlas.v3_final_class.isin(['VERIFIED', 'PROVISIONAL'])
                                                                    & atlas.clock.isin(UNRANKED_CLOCKS)).sum()),
               'contrast_rule': 'both models give the same contrast type, matching the category',
               'funnel': {'usable_significant_n3_rankable_category': step(base),
                          'agreed_matching_contrast_type': step(base & matches),
                          'no_agreed_adjudication_design_concern': step(base & matches & ~design_flagged),
                          'context_admitted': step(base & matches & ~design_flagged & admitted),
                          'adjusted_exact_ci_excludes_zero': {'rows': int(len(eligible))},
                          'after_dedup_by_study_and_label': {'rows': int(len(ordered))}},
               'excluded_by_contrast_type': atlas[base & ~matches].drop_duplicates(KEY).contrast_type
                                            .fillna('NOT_CLASSIFIED_OR_DISAGREED').value_counts().to_dict(),
               'excluded_by_design_audit': (pd.Series([r for rs in atlas[base & matches & ~design_flagged & ~admitted]
                                                       .drop_duplicates(KEY).design_audit_reasons.dropna() for r in rs])
                                            .value_counts().to_dict() if audit is not None else None),
               'outputs': outputs, 'atlas_sha256': sha256(atlas_path),
               'classification_sha256': sha256(Path(classification_dir) / 'comparisons.jsonl'),
               'design_audit_sha256': [sha256(Path(d) / 'comparisons.jsonl') for d in
                                       ([design_audit_dir] if isinstance(design_audit_dir, (str, Path)) else design_audit_dir or [])],
               'module_sha256': sha256(__file__), 'created': time.strftime('%Y-%m-%dT%H:%M:%S%z')}
    (out_dir / 'summary.json').write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary, indent=1, default=str))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--atlas', type=Path, required=True)
    parser.add_argument('--classification', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--design-audit', type=Path, action='append', default=[],
                        help='design-audit run; repeat for runs restricted to comparisons an earlier run did not cover')
    args = parser.parse_args(argv)
    run(args.atlas, args.classification, args.out, args.design_audit)


if __name__ == '__main__':
    main()
