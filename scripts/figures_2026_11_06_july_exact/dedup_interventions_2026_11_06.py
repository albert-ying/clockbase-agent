#!/usr/bin/env python3
"""One displayed entry per intervention in Figure 3 (display only; the ranking itself is unchanged).

Ranked comparisons can name the same intervention several times (time points, doses, backgrounds, two
sub-series, or the same exposure filed under two categories). For display, each intervention keeps its
best-ranked comparison (largest absolute CI bound of Hedges' g nearest zero) and appears in one panel only,
the category of that best comparison. The intervention identity is, in order of preference:
  1. the canonical target that two models assigned for the replication analysis
     (intervention-meta-2026-10-17/primary/meta_input_effects.csv, joined on GEO series, comparison and condition);
  2. the canonical target of the ranking table;
  3. the comparison label with time, dose and schedule qualifiers removed;
then the synonym table below merges names that denote the same manipulation, and each GEO series contributes at
most one displayed entry (its best-ranked comparison). Every merge is written to
evidence/figures-rebuild-2026-11-06/figure3/display_dedup.csv. All ranked comparisons stay in Table S17.

  python dedup_interventions_2026_11_06.py
"""
import os
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
WORK = ROOT / 'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/july_exact'
OUTDIR = ROOT / 'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06/figure3'
META = Path(os.environ.get('CLOCKBASE_META') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'intervention-meta-2026-10-17') / 'primary/meta_input_effects.csv'
JUDGED = Path(os.environ.get('CLOCKBASE_JUDGEMENT') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'followup-judgment-001') / 'judgement_age-lowering.csv'
RANKED = WORK / 'result/volcano_ranked_comparisons_july.csv'

# Names that denote the same manipulation (left: normalized key, right: shared key).
SYNONYMS = {
    'lps': 'lipopolysaccharide', 'lps treatment': 'lipopolysaccharide', 'seb+lps exposure': 'lipopolysaccharide',
    'pexidartinib': 'csf1r inhibitor plx3397', 'plx3397': 'csf1r inhibitor plx3397',
    'csf1r inhibitor plx3397': 'csf1r inhibitor plx3397',
    'lanifibranor 8w': 'lanifibranor', 'cold exposure': 'cold exposure', '4°c': 'cold exposure',
    'heat shock': 'heat stress', 'heat stress in mageb4 ko': 'heat stress',
    'sod1 with mutant tbk1': 'sod1', 'sod1 with wild-type tbk1': 'sod1',
    'mir155': 'mir-155', 'mir-155 overexpression': 'mir-155', 'mir-155ki': 'mir-155',
    'overload': 'mechanical overload', 'synergist ablation (overload)': 'mechanical overload',
    'diabetes with anril ko': 'diabetes', 'type 2 diabetes': 'diabetes',
    '2weeks bleomycin': 'bleomycin', 'progressive eae': 'eae',
    '5xfad': "alzheimer's disease model", 'amyloid pathology': "alzheimer's disease model",
    "alzheimer's disease model (arctic)": "alzheimer's disease model", 'app/ps1': "alzheimer's disease model",
}
QUALIFIERS = re.compile(r'\b(\d+\s*(h|hr|hrs|hour|hours|d|day|days|w|wk|week|weeks|weeks?|mo|month|months)\b|'
                        r'at zt\d+|zt\d+|day \d+|\d+\s*weeks?|\d+w|\d+ug|treatment|exposure)\b')


def normalize(text):
    text = str(text or '').casefold().strip()
    text = re.sub(r'\(.*?\)', ' ', text)
    text = QUALIFIERS.sub(' ', text)
    text = re.sub(r'\s+', ' ', text).strip(' ,;-')
    # disease models of one disease share a key
    if re.search(r'\bcold\b', text) or text in ('4°c', '6ºc'):
        return 'cold exposure'
    if 'huntington' in text:
        return "huntington's disease model"
    if 'alzheimer' in text or re.search(r'\bad\b', text) or text in ('5xfad', 'amyloid pathology', 'app/ps1'):
        return "alzheimer's disease model"
    return SYNONYMS.get(text, text)


def identity(frame, meta):
    keys = ['gse_id', 'comparison_name', 'condition_name']
    m = meta[keys + ['canonical_target']].dropna(subset=['canonical_target']).drop_duplicates(keys)
    out = frame.merge(m.rename(columns={'canonical_target': 'meta_target'}), on=keys, how='left')
    source = out['meta_target'].where(out['meta_target'].notna(), out.get('canonical_target'))
    source = source.where(source.notna(), out['label'])
    out['display_key'] = [normalize(s) for s in source]
    # a key that is still unlike any synonym falls back to the normalized label when the label maps to one
    label_key = [normalize(s) for s in out['label']]
    out['display_key'] = [lk if lk in SYNONYMS.values() and dk not in SYNONYMS.values() else dk
                          for dk, lk in zip(out['display_key'], label_key)]
    return out


def main():
    meta = pd.read_csv(META, low_memory=False)
    ranked = pd.read_csv(RANKED, low_memory=False)
    # idempotent: drop the columns an earlier run added
    ranked = ranked.drop(columns=[c for c in ('meta_target', 'display_key', 'display_category', 'display_best')
                                  if c in ranked.columns])
    ranked = identity(ranked, meta)
    ranked['_abs_bound'] = ranked['rank_bound'].abs()
    best = ranked.sort_values('_abs_bound', ascending=False).drop_duplicates('display_key').drop_duplicates('gse_id')
    ranked['display_category'] = ranked['display_key'].map(dict(zip(best['display_key'], best['condition_category'])))
    ranked['display_best'] = False
    ranked.loc[best.index, 'display_best'] = True
    ranked.drop(columns=['_abs_bound']).to_csv(RANKED, index=False)
    judged = identity(pd.read_csv(JUDGED, low_memory=False), meta)
    judged['_abs_bound'] = judged['rank_bound'].abs()
    keep = judged.sort_values('_abs_bound', ascending=False).drop_duplicates('display_key').drop_duplicates('gse_id').index
    judged['display_best'] = judged.index.isin(keep)
    judged.drop(columns=['_abs_bound']).to_csv(WORK / 'result/judgement_display_2026_11_06.csv', index=False)
    groups = ranked.groupby('display_key').filter(lambda g: len(g) > 1)
    groups = groups.sort_values(['display_key', 'rank_bound'], key=lambda s: s.abs() if s.name == 'rank_bound' else s,
                                ascending=[True, False])
    OUTDIR.mkdir(parents=True, exist_ok=True)
    groups[['display_key', 'display_category', 'display_best', 'label', 'condition_category', 'gse_id',
            'comparison_name', 'hedges_g', 'rank_bound']].to_csv(OUTDIR / 'display_dedup.csv', index=False)
    print(f'{len(ranked)} ranked rows, {ranked.display_key.nunique()} interventions; '
          f'{groups.display_key.nunique()} merged groups; judgement {len(judged)} rows -> {int(judged.display_best.sum())}')


if __name__ == '__main__':
    main()
