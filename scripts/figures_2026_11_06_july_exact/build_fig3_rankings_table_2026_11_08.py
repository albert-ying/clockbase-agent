#!/usr/bin/env python3
"""Figure 3 rankings table: one row per intervention and direction with its rank and score in every ranking version.

The figure rankings and legacy calculation fields are read from saved verifier outputs:
  statistical  final-2026-10-10/comparison-ranking/ranked_comparisons_all.csv (all categories, both directions)
  legacy fields  composite-ranking-001/ranked_with_composite.csv (internal calculations)
  replication  intervention-meta-2026-10-17/primary/ranked_replicated_age-{lowering,raising}.csv
  judgement    followup-judgment-001/judgement_age-{lowering,raising}.csv, with the three agreed design
               answers of followup-judgment-001/ranked_with_judgement.csv
An intervention is the display_key of Figure 3 (identity() of dedup_interventions_2026_11_06.py applied to the ranking
as the figure applies it, that is without the ranking's canonical_target column). Every version is ranked twice:
per comparison in statistical order and per intervention after the figure's display rule (walking down the version's order,
the first comparison of each intervention, then the first of each GEO series). The replication version pools across
series, so only the first unit per intervention is kept. Nothing is judged, pooled or typed here; every value comes
from the source tables. The table is checked against the current Figure 3 records (40 b-e bars, 4 f bars, 5 g bars,
43 panel a labels, the f and g rank columns and the ranks quoted in the text) and the build stops on any mismatch.

  python build_fig3_rankings_table.py
Outputs (beside this script): fig3_rankings.csv, fig3_rankings_columns.csv (definition of every column),
fig3_rankings_provenance.json (inputs with UTC modification times and SHA-256, counts, checks passed).
"""
import sys

sys.dont_write_bytecode = True  # read-only use of the repository modules: no __pycache__ written there

import hashlib
import json
import os
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

LOCAL = Path('/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base/clockbase_revision_local_2026-11-03')
ROOT = LOCAL / 'clockbase-agent'
SCRIPTS = ROOT / 'scripts/figures_2026_11_06_july_exact'
REBUILD = ROOT / 'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06'
EVIDENCE = REBUILD / 'figure3'
JULY = REBUILD / 'july_exact/result'
os.environ.setdefault('CLOCKBASE_VERIFIER_ROOT', str(LOCAL / 'data/clockbase-verifier-v3-2026-09-27'))
VER = Path(os.environ['CLOCKBASE_VERIFIER_ROOT'])
OUT = EVIDENCE  # fig3_rankings.csv and its column and provenance records sit beside the Figure 3 records

SOURCES = {
    'ranked_all': Path(os.environ.get('CLOCKBASE_FINAL') or VER / 'final-2026-10-10') / 'comparison-ranking/ranked_comparisons_all.csv',
    'composite': Path(os.environ.get('CLOCKBASE_COMPOSITE') or VER / 'composite-ranking-001') / 'ranked_with_composite.csv',
    'answers': Path(os.environ.get('CLOCKBASE_JUDGEMENT') or VER / 'followup-judgment-001') / 'ranked_with_judgement.csv',
    'judgement_lower': Path(os.environ.get('CLOCKBASE_JUDGEMENT') or VER / 'followup-judgment-001') / 'judgement_age-lowering.csv',
    'judgement_higher': Path(os.environ.get('CLOCKBASE_JUDGEMENT') or VER / 'followup-judgment-001') / 'judgement_age-raising.csv',
    'admitted_lower': Path(os.environ.get('CLOCKBASE_JUDGEMENT') or VER / 'followup-judgment-001') / 'judgement_age-lowering.csv',
    'replication_lower': Path(os.environ.get('CLOCKBASE_META') or VER / 'intervention-meta-2026-10-17') / 'primary/ranked_replicated_age-lowering.csv',
    'replication_higher': Path(os.environ.get('CLOCKBASE_META') or VER / 'intervention-meta-2026-10-17') / 'primary/ranked_replicated_age-raising.csv',
    'meta_input': Path(os.environ.get('CLOCKBASE_META') or VER / 'intervention-meta-2026-10-17') / 'primary/meta_input_effects.csv',
}
FIGURE_RECORDS = {
    'july_ranked': JULY / 'volcano_ranked_comparisons_july.csv',
    'be_rows': EVIDENCE / 'selected_rankings.csv',
    'be_cells': EVIDENCE / 'condition_grid.csv',
    'f_rows': EVIDENCE / 'panel_f_rows.csv',
    'g_rows': EVIDENCE / 'panel_g_rows.csv',
    'fg_ranks': EVIDENCE / 'statistical_ranks.csv',
    'a_labels': JULY / 'combined_volcano_annotations.csv',
}
MODULES = {'dedup': SCRIPTS / 'dedup_interventions_2026_11_06.py',
           'be': SCRIPTS / 'create_combined_modality_top_hits_v2.py',
           'f': SCRIPTS / 'create_figure3_replication_2026_11_07.py',
           'g': SCRIPTS / 'create_figure3_judgement_2026_11_06.py'}
for _p in list(SOURCES.values()) + list(FIGURE_RECORDS.values()) + list(MODULES.values()):
    if not _p.exists():
        raise SystemExit(f'Cannot proceed - required data file not found: {_p}')

sys.path.insert(0, str(SCRIPTS))
import dedup_interventions_2026_11_06 as dedup                    # noqa: E402  identity(), normalize()
import create_combined_modality_top_hits_v2 as be_fig              # noqa: E402  b-e label, cell and reason rules
import create_figure3_replication_2026_11_07 as f_fig              # noqa: E402  f display names and dagger
import create_figure3_judgement_2026_11_06 as g_fig                # noqa: E402  g label qualifier and dagger

KEYS4 = ['modality', 'gse_id', 'comparison_name', 'condition_name']
KEYS3 = ['gse_id', 'comparison_name', 'condition_name']
DIRECTIONS = ('lower', 'higher')
STATE_WORDS = {'met': 'met', 'not met': 'not met', 'unknown': 'not stated or not agreed'}
PANEL_OF = dict(zip(be_fig.CATEGORIES, 'bcde'))

COLUMNS = {
    'display_key': "Intervention identity of Figure 3: identity() of dedup_interventions_2026_11_06.py on the "
                   "statistical ranking without its canonical_target column (meta-analysis canonical_target joined "
                   "on GEO series, comparison and condition, else the comparison label; time, dose and schedule "
                   "qualifiers removed; synonym table applied). Rows found only in the replication version use "
                   "normalize() of the meta-analysis canonical_target.",
    'direction': "lower = lower clock outputs (Hedges' g < 0, age-lowering); higher = higher clock outputs (g > 0, "
                 "age-raising).",
    'label': "Name of the intervention as Figure 3 prints it, without the dagger: the b-e bar text, else the f bar "
             "text, else the g bar text, else the label of the best statistical comparison cut to 45 characters as "
             "in b-e, else the meta-analysis canonical_target.",
    'label_source': "Where label comes from: panel b-e, panel f, panel g, ranking label or meta-analysis "
                    "canonical_target.",
    'category': "Category (Genetic, Drug, Disease, Environment) of the best statistical comparison of the row; for "
                "rows found only in the replication version, the meta-analysis condition_category (which also uses "
                "Diet and Other).",
    'panels': "Panels of the current Figure 3 that print this intervention in this direction (a; b, c, d or e; f; g).",
    'printed_a': "Label text in panel a for comparisons of this row (the b-e bar text, or the followed-up comparison "
                 "label with the dagger); several are joined by ' | '.",
    'printed_be': "Bar text in panels b-e (July label, cut to 45 characters, qualified by the followed-up "
                  "comparison that shares its canonical_target).",
    'printed_f': "Bar text in panel f, with the dagger where printed.",
    'printed_g': "Bar text in panel g, with the dagger where printed.",
    'followed_up': "True when the row contains one of the three comparisons followed up in aged mice (ouabain "
                   "GSE122080, PLX5622 GSE134151, cold exposure GSE165974; TESTED_KEYS of the b-e generator).",
    'followed_up_gse_id': "GEO series of the followed-up comparison.",
    'followed_up_comparison_name': "Comparison name of the followed-up comparison.",
    'followed_up_statistical_rank_comparison': "Statistical rank of the followed-up comparison itself (of "
                                               "statistical_n_comparisons); differs from statistical_rank_comparison "
                                               "when another comparison of the intervention ranks higher.",
    'followed_up_composite_rank_comparison': "Composite rank of the followed-up comparison itself (of "
                                             "composite_n_comparisons).",
    'n_comparisons': "Ranked comparisons of this intervention in this direction (of the 946 ranked comparisons).",
    'n_series': "Distinct GEO series among those comparisons.",
    'gse_id': "GEO series of the best statistical comparison of the row (largest statistical_score).",
    'comparison_name': "Comparison name of the best statistical comparison.",
    'condition_name': "Condition name of the best statistical comparison.",
    'modality': "Modality of the best statistical comparison (mouse RNA-seq, human RNA-seq, human DNAm).",
    'hedges_g': "Corrected Hedges' g of the best statistical comparison.",
    'g_ci_low': "Lower bound of the confidence interval of that g (ci_level of the ranking, Bonferroni over "
                "selection_k).",
    'g_ci_high': "Upper bound of that interval.",
    'v3_p': "Verified P value of the best statistical comparison.",
    'v3_fdr': "Verified BH-FDR of the best statistical comparison.",
    'statistical_score': "Statistical score: absolute bound of the confidence interval of Hedges' g nearest zero "
                         "(|rank_bound|) of the best statistical comparison; the bar length of b-e.",
    'statistical_rank_comparison': "Comparison-level statistical rank of the best statistical comparison: its place "
                                   "among all ranked comparisons of the direction ordered by statistical_score "
                                   "(statistical_rank of ranked_with_judgement.csv; the rank in the text, e.g. "
                                   "ouabain 16 of 198).",
    'statistical_n_comparisons': "Ranked comparisons in the direction (198 lower, 748 higher).",
    'statistical_display_entry': "True when the best statistical comparison is a displayed entry of Figure 3 "
                                 "(display_best of the figure: walking down statistical_score across both "
                                 "directions, first comparison of each display_key, then first of each GEO series).",
    'statistical_rank_intervention': "Intervention-level statistical rank: place among the displayed entries of the "
                                     "direction ordered by statistical_score; empty when the row has no displayed "
                                     "entry (see dedup_notes).",
    'statistical_n_interventions': "Displayed entries in the direction.",
    'statistical_rank_in_category': "Place among the displayed entries of the same category and direction; b-e show "
                                    "ranks 1-5.",
    'statistical_n_in_category': "Displayed entries in that category and direction.",
    'composite_score': "Composite score of the best composite comparison of the row: min(-log10(max(v3_p, 1e-300)), "
                       "30) x relevance_score / 10 (relevance_score: blinded rubric score, mean of two models, 0-10); "
                       "Genetic, Drug and Environment only.",
    'composite_gse_id': "GEO series of the best composite comparison (smallest composite_rank).",
    'composite_comparison_name': "Comparison name of the best composite comparison.",
    'composite_relevance_score': "relevance_score of that comparison.",
    'composite_neg_log10_p_capped': "min(-log10 P, 30) of that comparison.",
    'composite_rank_comparison': "Comparison-level composite rank (composite_rank of ranked_with_composite.csv): place "
                                 "among the scored comparisons of the direction by statistical rank, preserving the "
                                 "statistical rank (e.g. ouabain 9 of 168).",
    'composite_n_comparisons': "Scored comparisons in the direction (168 lower, 601 higher).",
    'composite_rank_intervention': "Intervention-level composite rank: walking down composite_rank_comparison within "
                                   "the direction, the first comparison of each display_key, then the first of each "
                                   "GEO series; empty when the row's first comparison is dropped by the GEO-series "
                                   "rule.",
    'composite_n_interventions': "Composite entries in the direction after that rule.",
    'replication_canonical_target': "canonical_target of the meta-analysis unit (best unit of this intervention).",
    'replication_category': "condition_category of that unit (the meta-analysis does not pool across categories).",
    'replication_pooled_g': "Pooled Hedges' g (DerSimonian-Laird across projects, primary meta-analysis).",
    'replication_ci_low': "Lower bound of the 95% interval of the pooled g.",
    'replication_ci_high': "Upper bound of that interval.",
    'replication_score': "Replication score: bound of the pooled interval nearest zero (rank_bound of the source); "
                         "the bar length of f.",
    'replication_k': "Independent projects pooled (k).",
    'replication_series': "GEO series pooled (n_series).",
    'replication_fdr': "BH-FDR over interventions of the pooled effect.",
    'replication_rank_unit': "Place of the unit in ranked_replicated_age-{lowering,raising}.csv (robust replicated "
                             "units ordered by replication_score).",
    'replication_n_units': "Robust replicated units in the direction (4 lower, 73 higher).",
    'replication_rank_intervention': "Intervention-level replication rank: units walked in that order, the first unit "
                                     "of each display_key (normalize() of canonical_target); f shows the 4 lower.",
    'replication_n_interventions': "Replicated interventions in the direction after that rule.",
    'replication_best_single_rank': "Best single rank of panel f: smallest statistical_rank_comparison among the "
                                    "ranked comparisons of the direction whose canonical_target equals the unit's; "
                                    "empty (f prints 'not ranked') when there is none.",
    'replication_modalities': "Modalities of the pooled effects.",
    'judgement_admitted_comparisons': "Comparisons of the row admitted by the judgement version (agreed answers; "
                                      "lower: old at start, otherwise healthy, untreated control; higher: otherwise "
                                      "healthy, untreated control), Disease included.",
    'judgement_score': "Absolute adjusted CI bound nearest zero of the row's first admitted comparison in the judgement "
                       "order (the bar length of g).",
    'judgement_gse_id': "GEO series of that comparison.",
    'judgement_comparison_name': "Comparison name of that comparison.",
    'judgement_entry_statistical_rank': "statistical rank of that comparison (of statistical_n_comparisons); the "
                                        "'Statistical rank' column of g.",
    'judgement_rank_comparison': "Comparison-level judgement rank (judgement_rank): place among the "
                                 "admitted comparisons of the direction in statistical order.",
    'judgement_n_comparisons': "Number of admitted comparisons in the direction, including all categories.",
    'judgement_rank_intervention': "Intervention-level judgement rank: walking down judgement_rank_comparison, the "
                                   "first comparison of each display_key, then the first of each GEO series "
                                   "(computed with the same identity and series rule in each direction).",
    'judgement_n_interventions': "Judgement entries in the direction after that rule (5 lower).",
    'answer_age_at_start': "Agreed answer of the two judging models for the best statistical comparison "
                           "(OLD, ADULT, YOUNG, GERMLINE_OR_PRENATAL, NOT_STATED; 'not agreed' when they differed).",
    'answer_background': "Agreed background answer (HEALTHY, DISEASE_MODEL, NOT_STATED, or 'not agreed').",
    'answer_control': "Agreed control answer (UNTREATED, ACTIVE_COMPARATOR, or 'not agreed').",
    'state_old_at_start': "Cell state of b-e: met, not met, or not stated or not agreed.",
    'state_otherwise_healthy': "Cell state of b-e for otherwise healthy (HEALTHY met, DISEASE_MODEL not met).",
    'state_untreated_control': "Cell state of b-e for untreated control (UNTREATED met, ACTIVE_COMPARATOR not met).",
    'reason': "Unmet conditions in cell order (germline, adult, young, disease model, active control), or 'all three "
              "met', from the b-e reason rule; computed for every row with a ranked comparison.",
    'printed_reason': "Reason text printed beside the bar in b-e (lower bars only).",
    'dedup_notes': "Why a version has no intervention-level rank for this row (the GEO-series rule, or the figure's "
                   "display entry of the intervention lies in the other direction).",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def first_per_key_then_series(frame):
    """Rows of frame, already in ranking order, kept by the figure rule: the first of each display_key, then the
    first of each GEO series among those (dedup_interventions_2026_11_06.py and the panel g rule)."""
    first = ~frame['display_key'].duplicated()
    return first & ~(frame['gse_id'].where(first).duplicated() & first)


def series_note(version, frame, kept):
    """For each display_key whose first comparison was dropped by the GEO-series rule, the key that holds the series."""
    firsts = frame[~frame['display_key'].duplicated()]
    holder = dict(zip(frame.loc[kept, 'gse_id'], frame.loc[kept, 'display_key']))
    return {k: f"{version}: GEO series {g} shown for '{holder[g]}'"
            for k, g, keep in zip(firsts['display_key'], firsts['gse_id'], kept[firsts.index]) if not keep}


def load():
    try:
        src = {k: pd.read_csv(p, low_memory=False) for k, p in SOURCES.items()}
        rec = {k: pd.read_csv(p, low_memory=False) for k, p in FIGURE_RECORDS.items()}
    except (OSError, pd.errors.ParserError) as exc:
        raise SystemExit(f'Cannot read an input table: {exc}')
    return src, rec


def comparisons(src, rec, checks):
    """The 946 ranked comparisons with identity, display entry, statistical, composite and judgement fields."""
    ra = src['ranked_all'].rename(columns={'canonical_target': 'ranking_canonical_target'})
    ra['_order'] = np.arange(len(ra))
    ra['direction'] = np.where(ra['hedges_g'] < 0, 'lower', 'higher')
    assert (ra['direction'].map({'lower': 'age-lowering', 'higher': 'age-raising'}) == ra['direction_v3']).all()
    nearest = np.where(ra['g_ci_low'] > 0, ra['g_ci_low'], np.where(ra['g_ci_high'] < 0, -ra['g_ci_high'], 0))
    assert np.allclose(nearest, ra['rank_bound'].abs(), rtol=0, atol=1e-12)
    ra['statistical_score'] = ra['rank_bound'].abs()
    # statistical rank: within direction by statistical_score, stable on the published order
    ra['statistical_rank_comparison'] = 0
    for d in DIRECTIONS:
        idx = ra.index[ra['direction'] == d]
        order = idx[np.argsort(-ra.loc[idx, 'statistical_score'].to_numpy(), kind='stable')]
        ra.loc[order, 'statistical_rank_comparison'] = np.arange(1, len(order) + 1)
    # identity exactly as the figure: the July ranking table carries no canonical_target column
    ident = dedup.identity(ra, src['meta_input'])
    assert len(ident) == len(ra) and (ident['_order'].to_numpy() == ra['_order'].to_numpy()).all()
    ident.index = ra.index
    july = rec['july_ranked']
    assert (july[KEYS4 + ['clock']].to_numpy() == ra[KEYS4 + ['clock']].to_numpy()).all()
    assert (ident['display_key'].to_numpy() == july['display_key'].to_numpy()).all()
    ra['display_key'] = ident['display_key']
    ra['meta_target'] = ident['meta_target']
    # display entries: the three lines of dedup_interventions_2026_11_06.main()
    tmp = ra.assign(_abs_bound=ra['rank_bound'].abs())
    best = tmp.sort_values('_abs_bound', ascending=False).drop_duplicates('display_key').drop_duplicates('gse_id')
    ra['display_best'] = ra.index.isin(best.index)
    assert (ra['display_best'].to_numpy() == july['display_best'].astype(bool).to_numpy()).all()
    checks.append('display_key and display_best equal those of volcano_ranked_comparisons_july.csv (946 rows)')
    # answers and statistical rank of the judgement source
    ans = src['answers']
    assert not ans.duplicated(KEYS4).any()
    ra = ra.merge(ans[KEYS4 + ['statistical_rank', 'canonical_target', 'age_at_start', 'background', 'control',
                               'judgement_admitted']].rename(columns={'canonical_target': 'answers_canonical_target'}),
                  on=KEYS4, how='left', validate='one_to_one')
    assert ra['statistical_rank'].notna().all()
    assert (ra['statistical_rank'] == ra['statistical_rank_comparison']).all()
    assert (ra['answers_canonical_target'].fillna('') == ra['ranking_canonical_target'].fillna('')).all()
    # composite
    comp = src['composite']
    assert not comp.duplicated(KEYS4).any()
    ra = ra.merge(comp[KEYS4 + ['category', 'statistical_rank', 'composite', 'composite_rank', 'relevance_score',
                                'neg_log10_p_capped', 'v3_p']].rename(
        columns={'statistical_rank': '_rank_c', 'relevance_score': '_rel_c', 'v3_p': '_p_c'}),
                  on=KEYS4, how='left', validate='one_to_one')
    assert (ra['category'] == ra['condition_category']).all() and (ra['_rank_c'] == ra['statistical_rank']).all()
    scored = ra['composite'].notna()
    assert set(ra.loc[scored, 'condition_category']) == {'Genetic', 'Drug', 'Environment'}
    assert not ra.loc[ra['condition_category'] != 'Disease', 'composite'].isna().any()
    recomputed = (np.minimum(-np.log10(np.maximum(ra.loc[scored, '_p_c'], 1e-300)), 30)
                  * ra.loc[scored, '_rel_c'] / 10)
    assert np.allclose(recomputed, ra.loc[scored, 'composite'], rtol=1e-12, atol=0)
    for d in DIRECTIONS:
        s = ra[scored & (ra['direction'] == d)].sort_values(['composite', 'statistical_rank'],
                                                             ascending=[False, True], kind='stable')
        assert s['composite_rank'].tolist() == list(range(1, len(s) + 1))
    checks.append('composite recomputed from v3_p and relevance_score; composite_rank recomputed within direction')
    ra = ra.drop(columns=['category', '_rank_c', '_rel_c', '_p_c', 'answers_canonical_target'])
    return ra


def judgement_tables(src, cmp_, checks):
    """Admitted comparisons in statistical order with the figure display rule."""
    out = {}
    for d, key in (('lower', 'judgement_lower'), ('higher', 'judgement_higher')):
        j = src[key].drop(columns=['direction'], errors='ignore').copy().sort_values('statistical_rank', kind='stable')
        assert j['judgement_admitted'].all()
        assert j['judgement_rank'].tolist() == list(range(1, len(j) + 1))
        j = j.merge(cmp_[KEYS4 + ['display_key', 'statistical_rank_comparison', 'direction']],
                    on=KEYS4, how='left', validate='one_to_one')
        assert (j['direction'] == d).all()
        assert (j['statistical_rank'] == j['statistical_rank_comparison']).all()
        assert j['rank_bound'].abs().is_monotonic_decreasing
        j['kept'] = first_per_key_then_series(j)
        j['display_label'] = j['label'].str.replace(r' at ZT\d+', '', regex=True)
        j.loc[j['comparison_name'].eq('Overload vs Sham (Senolytic)'), 'display_label'] = 'Overload in senolytic-treated mice'
        j['judgement_rank_intervention'] = np.where(j['kept'], j['kept'].cumsum(), np.nan)
        checks.append(f'judgement ({d}): {len(j)} admitted comparisons in statistical order, '
                      f'{int(j["kept"].sum())} after one intervention and series')
        out[d] = j.reset_index(drop=True)
    return out


def replication_tables(src, cmp_, checks):
    """Robust replicated units per direction, in source order, with the first unit of each intervention ranked."""
    out = {}
    for d, key in (('lower', 'replication_lower'), ('higher', 'replication_higher')):
        r = src[key].copy().reset_index(drop=True)
        assert r['rank_bound'].is_monotonic_decreasing and (r['tier'] == 'replicated').all()
        bound = np.where(r['ci_low'] > 0, r['ci_low'], np.where(r['ci_high'] < 0, -r['ci_high'], np.nan))
        assert np.allclose(bound, r['rank_bound'], rtol=0, atol=1e-12)
        assert ((r['estimate'] < 0) == (d == 'lower')).all() and (r['fdr'] < 0.05).all() and (r['k'] >= 2).all()
        r['replication_rank_unit'] = np.arange(1, len(r) + 1)
        r['display_key'] = [dedup.normalize(t) for t in r['canonical_target']]
        first = ~r['display_key'].duplicated()
        r['replication_rank_intervention'] = np.where(first, first.cumsum(), np.nan)
        # best single rank of panel f: the best ranked comparison of the direction with this canonical_target
        same = cmp_[cmp_['direction'] == d]
        best = same.groupby('ranking_canonical_target')['statistical_rank_comparison'].min()
        r['replication_best_single_rank'] = r['canonical_target'].map(best)
        out[d] = r
        checks.append(f'replication ({d}): {len(r)} units, {int(first.sum())} interventions after one unit per '
                      'display_key')
    return out


def build(src, rec, checks):
    cmp_ = comparisons(src, rec, checks)
    jt = judgement_tables(src, cmp_, checks)
    rt = replication_tables(src, cmp_, checks)

    # statistical display entries ranked within direction and within category
    ent = cmp_[cmp_['display_best']].sort_values('statistical_rank_comparison', kind='stable').copy()
    ent = ent.sort_values(['direction'], kind='stable')
    ent['statistical_rank_intervention'] = ent.groupby('direction').cumcount() + 1
    ent['statistical_rank_in_category'] = ent.groupby(['direction', 'condition_category']).cumcount() + 1
    n_ent = ent.groupby('direction').size()
    n_cat = ent.groupby(['direction', 'condition_category']).size()

    # composite entries
    comp_rank, comp_notes = {}, {}
    for d in DIRECTIONS:
        s = cmp_[cmp_['composite'].notna() & (cmp_['direction'] == d)].sort_values('composite_rank', kind='stable')
        kept = first_per_key_then_series(s)
        k = s[kept]
        comp_rank.update({(key, d): i + 1 for i, key in enumerate(k['display_key'])})
        comp_notes.update({(key, d): n for key, n in series_note('composite', s, kept).items()})

    judge_notes = {}
    for d, j in jt.items():
        judge_notes.update({(key, d): n for key, n in series_note('judgement', j, j['kept']).items()})

    # b-e picks and their printed text, as the b-e generator selects and labels them
    be = ent[ent['statistical_rank_in_category'] <= be_fig.TOP_N].copy()
    be['intervention'] = be['label']
    be['display'] = be.apply(be_fig.make_label, axis=1)
    picks = {}
    for cat in be_fig.CATEGORIES:
        sub = be[be['condition_category'] == cat]
        sub = pd.concat([sub[sub['direction'] == 'higher'], sub[sub['direction'] == 'lower']])
        sub = sub.sort_values('statistical_score', ascending=False, kind='stable').copy()
        sub['display'], _ = be_fig.qualify_with_tested(sub, src['answers'])
        picks[cat] = sub
    be = pd.concat([p.assign(category=c) for c, p in picks.items()], ignore_index=True)
    printed_be = {(r.display_key, r.direction): r.display for r in be.itertuples()}

    # g text: the label without " at ZT<n>", the senolytic arm named, b-e qualifier, dagger
    gl = jt['lower'][jt['lower']['kept']].copy()
    gl['display'] = gl['label'].str.replace(r' at ZT\d+', '', regex=True)
    gl.loc[gl['comparison_name'].eq('Overload vs Sham (Senolytic)'), 'display'] = 'Overload in senolytic-treated mice'
    assert gl['display'].tolist() == gl['display_label'].tolist()
    other = [len(b) > 0 and not ((b['gse_id'] == r.gse_id) & (b['comparison_name'] == r.comparison_name)).any()
             for b, r in zip([be[be['display'].str.casefold() == x.casefold()] for x in gl['display']],
                             gl.itertuples())]
    gl['display'] = [g_fig.qualifier(r) if o else r['display'] for o, (_, r) in zip(other, gl.iterrows())]
    tested_g = (gl['gse_id'] == g_fig.TESTED_KEY[0]) & (gl['comparison_name'] == g_fig.TESTED_KEY[1])
    gl.loc[tested_g, 'display'] = gl.loc[tested_g, 'display'] + g_fig.TESTED_MARK
    printed_g = {(r.display_key, 'lower'): r.display for r in gl.itertuples()}

    # f text
    fl = rt['lower']
    printed_f = {}
    for r in fl.itertuples():
        text = f_fig.DISPLAY[r.canonical_target] + (f_fig.TESTED_MARK if r.canonical_target in f_fig.TESTED else '')
        printed_f[(r.display_key, 'lower')] = text

    # followed-up comparisons and panel a labels (the 40 b-e texts plus the three followed-up comparisons)
    tested = pd.DataFrame(be_fig.TESTED_KEYS, columns=['gse_id', 'comparison_name']).merge(
        cmp_, on=['gse_id', 'comparison_name'], how='left', validate='one_to_one')
    assert tested['display_key'].notna().all() and len(tested) == 3
    followed = {(r.display_key, r.direction): r for r in tested.itertuples()}
    printed_a = {}
    for r in be.itertuples():
        printed_a.setdefault((r.display_key, r.direction), []).append(r.display)
    for r in tested.itertuples():
        printed_a.setdefault((r.display_key, r.direction), []).append(f'{r.label}{g_fig.TESTED_MARK}')

    # rows: every (display_key, direction) with a ranked comparison or a replicated unit
    rows = []
    groups = {k: g.sort_values('statistical_rank_comparison') for k, g in cmp_.groupby(['display_key', 'direction'])}
    rep_first = {(r.display_key, d): r for d, t in rt.items() for r in t[t['replication_rank_intervention'].notna()]
                 .itertuples()}
    keys = sorted(set(groups) | set(rep_first), key=lambda k: (DIRECTIONS.index(k[1]), k[0]))
    ent_by = {(r.display_key, r.direction): r for r in ent.itertuples()}
    judge_first = {}
    for d, j in jt.items():
        for r in j.itertuples():
            judge_first.setdefault((r.display_key, d), r)
    admitted_n = cmp_[cmp_['judgement_admitted'].astype(bool)].groupby(['display_key', 'direction']).size()
    any_entry = {r.display_key: r.direction for r in ent.itertuples()}
    holder = dict(zip(ent['gse_id'], ent['display_key']))
    n_cmp = cmp_.groupby('direction').size()
    n_comp = cmp_[cmp_['composite'].notna()].groupby('direction').size()

    for key, d in keys:
        row = {'display_key': key, 'direction': d}
        g = groups.get((key, d))
        rep = rep_first.get((key, d))
        notes = []
        if g is not None:
            b = g.iloc[0]
            row.update({'category': b['condition_category'], 'n_comparisons': len(g), 'n_series': g['gse_id'].nunique(),
                        'gse_id': b['gse_id'], 'comparison_name': b['comparison_name'],
                        'condition_name': b['condition_name'], 'modality': b['modality'], 'hedges_g': b['hedges_g'],
                        'g_ci_low': b['g_ci_low'], 'g_ci_high': b['g_ci_high'], 'v3_p': b['v3_p'],
                        'v3_fdr': b['v3_fdr'], 'statistical_score': b['statistical_score'],
                        'statistical_rank_comparison': int(b['statistical_rank_comparison']),
                        'statistical_n_comparisons': int(n_cmp[d])})
            e = ent_by.get((key, d))
            row['statistical_display_entry'] = e is not None
            if e is not None:
                assert e.gse_id == b['gse_id'] and e.comparison_name == b['comparison_name']
                row.update({'statistical_rank_intervention': e.statistical_rank_intervention,
                            'statistical_n_interventions': int(n_ent[d]),
                            'statistical_rank_in_category': e.statistical_rank_in_category,
                            'statistical_n_in_category': int(n_cat[(d, e.condition_category)])})
            elif key in any_entry:
                notes.append(f'statistical: the displayed entry of this intervention has {any_entry[key]} outputs')
            else:
                top = cmp_[cmp_['display_key'] == key].sort_values('statistical_score', ascending=False).iloc[0]
                notes.append(f"statistical: GEO series {top['gse_id']} shown for '{holder[top['gse_id']]}'")
            # condition states of the best statistical comparison (b-e rules)
            states = be_fig.condition_states(b)
            row.update({f'answer_{f}': ('not agreed' if pd.isna(b[f]) else b[f]) for f in be_fig.FIELDS})
            row.update({c: STATE_WORDS[s] for c, s in zip(be_fig.STATE_COLUMNS, states)})
            row['reason'] = be_fig.reason_text(b)
            gc = g[g['composite'].notna()].sort_values('composite_rank')
            if len(gc):
                c = gc.iloc[0]
                row.update({'composite_score': c['composite'], 'composite_gse_id': c['gse_id'],
                            'composite_comparison_name': c['comparison_name'],
                            'composite_relevance_score': c['relevance_score'],
                            'composite_neg_log10_p_capped': c['neg_log10_p_capped'],
                            'composite_rank_comparison': int(c['composite_rank']),
                            'composite_n_comparisons': int(n_comp[d]),
                            'composite_rank_intervention': comp_rank.get((key, d), np.nan),
                            'composite_n_interventions': len([k for k in comp_rank if k[1] == d])})
                if (key, d) in comp_notes:
                    notes.append(comp_notes[(key, d)])
            row['judgement_admitted_comparisons'] = int(admitted_n.get((key, d), 0))
        else:
            row.update({'category': rep.condition_category, 'n_comparisons': 0, 'n_series': 0,
                        'statistical_display_entry': False, 'judgement_admitted_comparisons': 0})
        jf = judge_first.get((key, d))
        if jf is not None:
            row.update({'judgement_score': abs(jf.rank_bound), 'judgement_gse_id': jf.gse_id,
                        'judgement_comparison_name': jf.comparison_name,
                        'judgement_entry_statistical_rank': int(jf.statistical_rank_comparison),
                        'judgement_rank_comparison': int(jf.judgement_rank),
                        'judgement_n_comparisons': len(jt[d]),
                        'judgement_rank_intervention': jf.judgement_rank_intervention,
                        'judgement_n_interventions': int(jt[d]['kept'].sum())})
            if (key, d) in judge_notes:
                notes.append(judge_notes[(key, d)])
        if rep is not None:
            row.update({'replication_canonical_target': rep.canonical_target,
                        'replication_category': rep.condition_category, 'replication_pooled_g': rep.estimate,
                        'replication_ci_low': rep.ci_low, 'replication_ci_high': rep.ci_high,
                        'replication_score': rep.rank_bound, 'replication_k': rep.k, 'replication_series': rep.n_series,
                        'replication_fdr': rep.fdr, 'replication_rank_unit': rep.replication_rank_unit,
                        'replication_n_units': len(rt[d]),
                        'replication_rank_intervention': rep.replication_rank_intervention,
                        'replication_n_interventions': int(rt[d]['replication_rank_intervention'].notna().sum()),
                        'replication_best_single_rank': rep.replication_best_single_rank,
                        'replication_modalities': rep.modalities})
        f_ = followed.get((key, d))
        row['followed_up'] = f_ is not None
        if f_ is not None:
            row.update({'followed_up_gse_id': f_.gse_id, 'followed_up_comparison_name': f_.comparison_name,
                        'followed_up_statistical_rank_comparison': int(f_.statistical_rank_comparison),
                        'followed_up_composite_rank_comparison': f_.composite_rank})
        row['printed_be'] = printed_be.get((key, d), '')
        row['printed_f'] = printed_f.get((key, d), '')
        row['printed_g'] = printed_g.get((key, d), '')
        row['printed_a'] = ' | '.join(printed_a.get((key, d), []))
        row['printed_reason'] = row.get('reason', '') if (row['printed_be'] and d == 'lower') else ''
        dagger = g_fig.TESTED_MARK
        if row['printed_be']:
            row['label'], row['label_source'] = row['printed_be'], 'panel b-e'
        elif row['printed_f']:
            row['label'], row['label_source'] = row['printed_f'].removesuffix(dagger), 'panel f'
        elif row['printed_g']:
            row['label'], row['label_source'] = row['printed_g'].removesuffix(dagger), 'panel g'
        elif g is not None:
            row['label'] = be_fig.make_label(g.iloc[0].to_dict() | {'intervention': g.iloc[0]['label']})
            row['label_source'] = 'ranking label'
        else:
            row['label'], row['label_source'] = rep.canonical_target, 'meta-analysis canonical_target'
        panels = (['a'] if row['printed_a'] else []) + ([PANEL_OF[row['category']]] if row['printed_be'] else []) \
            + (['f'] if row['printed_f'] else []) + (['g'] if row['printed_g'] else [])
        row['panels'] = ';'.join(panels)
        row['dedup_notes'] = '; '.join(notes)
        rows.append(row)

    table = pd.DataFrame(rows).reindex(columns=list(COLUMNS))
    int_cols = [c for c in table.columns if c.endswith(('_rank_comparison', '_rank_intervention', '_rank_in_category',
                                                       '_n_comparisons', '_n_interventions', '_n_in_category',
                                                       '_rank_unit', '_n_units', '_best_single_rank'))
                or c in ('n_comparisons', 'n_series', 'replication_k', 'replication_series',
                         'judgement_admitted_comparisons', 'judgement_entry_statistical_rank')]
    for c in int_cols:
        table[c] = table[c].astype('Int64')
    table = table.sort_values(['direction', 'statistical_rank_comparison', 'replication_rank_unit'],
                              key=lambda s: s.map({'lower': 0, 'higher': 1}) if s.name == 'direction' else s,
                              na_position='last', kind='stable').reset_index(drop=True)
    return table, cmp_, be, gl, fl


def src_ouabain_rank(cmp_):
    admitted = cmp_[(cmp_.direction == 'lower') & cmp_.judgement_admitted].sort_values('statistical_rank')
    return list(admitted.gse_id).index(g_fig.TESTED_KEY[0]) + 1


def verify(table, cmp_, be, gl, fl, rec, checks):
    """Reproduce b-e, f, g, the panel a labels and the ranks quoted in the text from the table."""
    t = table
    # b-e: 40 bars, membership, order, text, values, cells and reasons
    sel = rec['be_rows']
    rows = []
    for cat in be_fig.CATEGORIES:
        s = t[(t['category'] == cat) & (t['statistical_rank_in_category'] <= be_fig.TOP_N)]
        s = pd.concat([s[s['direction'] == 'higher'], s[s['direction'] == 'lower']])
        rows.append(s.sort_values('statistical_score', ascending=False, kind='stable'))
    mine = pd.concat(rows, ignore_index=True)
    assert len(mine) == len(sel)
    for a, b in (('category', 'category'), ('gse_id', 'gse_id'), ('comparison_name', 'comparison_name'),
                 ('condition_name', 'condition_name'), ('printed_be', 'display'), ('modality', 'modality')):
        assert mine[a].tolist() == sel[b].tolist(), a
    assert np.allclose(mine['statistical_score'], sel['composite'], rtol=0, atol=1e-12)
    assert np.allclose(mine['hedges_g'], sel['hedges_g'], rtol=0, atol=1e-12)
    low = mine['direction'] == 'lower'
    for c in be_fig.STATE_COLUMNS:
        assert mine.loc[low, c].tolist() == sel.loc[low, c].map(STATE_WORDS).tolist(), c
    # an empty reason (no unmet condition, at least one not stated or not agreed) is read back from the CSV as NaN
    assert mine.loc[low, 'printed_reason'].tolist() == sel.loc[low, 'reason'].fillna('').tolist()
    grid = rec['be_cells']
    g2 = mine[low].merge(grid, on=['gse_id', 'comparison_name'], suffixes=('', '_grid'), validate='one_to_one')
    assert len(g2) == int(low.sum()) and (g2['printed_reason'] == g2['reason_grid'].fillna('')).all()
    checks.append('b-e: the 40 bars (category, order, comparison, text, modality, bar length, g) and the 20 lower '
                  'bars\' three cells and reason text reproduced from the table')
    # f: 4 bars
    fr = rec['f_rows']
    f = t[(t['direction'] == 'lower') & t['replication_rank_intervention'].notna()].sort_values(
        'replication_rank_intervention')
    assert f['printed_f'].tolist() == fr['display'].tolist()
    assert f['replication_canonical_target'].tolist() == fr['canonical_target'].tolist()
    for a, b in (('replication_pooled_g', 'estimate'), ('replication_ci_low', 'ci_low'),
                 ('replication_ci_high', 'ci_high'), ('replication_score', 'rank_bound'), ('replication_fdr', 'fdr')):
        assert np.allclose(f[a].astype(float), fr[b], rtol=0, atol=1e-12), a
    assert f['replication_k'].tolist() == fr['k'].tolist() and f['replication_series'].tolist() == fr['n_series'].tolist()
    text = [f_fig.NOT_RANKED if pd.isna(x) else str(int(x)) for x in f['replication_best_single_rank']]
    assert text == fr['rank_text'].astype(str).tolist()
    # where the intervention has ranked lower comparisons, the best one is the row's best statistical comparison
    has = f['replication_best_single_rank'].notna()
    assert (f.loc[has, 'replication_best_single_rank'] == f.loc[has, 'statistical_rank_comparison']).all()
    checks.append('f: the 4 bars (order, text, pooled g, CI, bar length, FDR, projects, series, best single rank) '
                  'reproduced from the table')
    # g: 5 bars
    gr = rec['g_rows']
    g = t[(t['direction'] == 'lower') & t['judgement_rank_intervention'].notna()].sort_values(
        'judgement_rank_intervention')
    assert g['printed_g'].tolist() == gr['display'].tolist()
    assert g['judgement_gse_id'].tolist() == gr['gse_id'].tolist()
    assert g['judgement_comparison_name'].tolist() == gr['comparison_name'].tolist()
    assert np.allclose(g['judgement_score'].astype(float), gr['statistical_score'], rtol=0, atol=1e-12)
    assert g['judgement_entry_statistical_rank'].tolist() == gr['statistical_rank'].tolist()
    assert g['judgement_rank_comparison'].tolist() == gr['judgement_rank'].tolist()
    checks.append('g: the displayed bars (order, text, comparison, CI bound, statistical rank and admitted rank) '
                  'reproduced from the table')
    # the rank columns of f and g as printed
    fg = rec['fg_ranks']
    printed = dict(zip(fg['bar'], fg['printed'].astype(str)))
    for r in f.itertuples():
        assert printed[r.printed_f] == (f_fig.NOT_RANKED if pd.isna(r.replication_best_single_rank)
                                        else str(int(r.replication_best_single_rank)))
    for r in g.itertuples():
        assert printed[r.printed_g] == str(int(r.judgement_entry_statistical_rank))
    assert set(fg['of']) == {int(t.loc[t['direction'] == 'lower', 'statistical_n_comparisons'].dropna().iloc[0])}
    # panel a: 43 labels
    ann = rec['a_labels']
    labels = [x for s in t['printed_a'] if s for x in s.split(' | ')]
    assert sorted(labels) == sorted(ann['label'].tolist())
    checks.append('a: the 43 labels (40 b-e texts and the 3 followed-up comparisons) reproduced from the table')
    # ranks quoted in the text
    o = t[(t['display_key'] == 'ouabain') & (t['direction'] == 'lower')].iloc[0]
    source = cmp_[(cmp_.direction == 'lower') & cmp_.gse_id.eq(g_fig.TESTED_KEY[0])
                  & cmp_.comparison_name.eq(g_fig.TESTED_KEY[1])].iloc[0]
    assert o['statistical_rank_comparison'] == source['statistical_rank_comparison']
    assert o['statistical_n_comparisons'] == int(cmp_.direction.eq('lower').sum())
    assert o['judgement_rank_comparison'] == src_ouabain_rank(cmp_) and o['judgement_rank_intervention'] == g.loc[g.display_key.eq('ouabain'), 'judgement_rank_intervention'].iloc[0]
    checks.append('ouabain statistical and judgement positions match the source comparisons')
    assert t.duplicated(['display_key', 'direction']).sum() == 0
    assert t['followed_up'].sum() == 3


def same_inputs_as_figure(checks):
    """The sources shared with the current Figure 3 build are byte-identical to those it recorded."""
    try:
        recorded = json.loads((EVIDENCE / 'input_provenance.json').read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f'Cannot read the Figure 3 input record: {exc}')
    mine = {str(p): p for p in list(SOURCES.values()) + list(FIGURE_RECORDS.values())}
    shared = [r for r in recorded if r['path'] in mine]
    for r in shared:
        assert sha256(mine[r['path']]) == r['sha256'], f"{r['path']} changed since the Figure 3 build"
    checks.append(f'{len(shared)} inputs shared with the Figure 3 build match its recorded SHA-256 '
                  '(evidence/figures-rebuild-2026-11-06/figure3/input_provenance.json)')


def main():
    checks = []
    same_inputs_as_figure(checks)
    src, rec = load()
    table, cmp_, be, gl, fl = build(src, rec, checks)
    verify(table, cmp_, be, gl, fl, rec, checks)
    table.to_csv(OUT / 'fig3_rankings.csv', index=False)
    pd.DataFrame({'column': list(COLUMNS), 'definition': list(COLUMNS.values())}).to_csv(
        OUT / 'fig3_rankings_columns.csv', index=False)
    counts = {d: {'rows': int((table['direction'] == d).sum()),
                  'with_ranked_comparison': int(((table['direction'] == d) & (table['n_comparisons'] > 0)).sum()),
                  'statistical_entries': int(((table['direction'] == d) & table['statistical_rank_intervention']
                                              .notna()).sum()),
                  'composite_scored': int(((table['direction'] == d) & table['composite_rank_comparison']
                                           .notna()).sum()),
                  'composite_entries': int(((table['direction'] == d) & table['composite_rank_intervention']
                                            .notna()).sum()),
                  'replication_entries': int(((table['direction'] == d) & table['replication_rank_intervention']
                                              .notna()).sum()),
                  'replication_only': int(((table['direction'] == d) & (table['n_comparisons'] == 0)).sum()),
                  'judgement_scored': int(((table['direction'] == d) & table['judgement_rank_comparison']
                                           .notna()).sum()),
                  'judgement_entries': int(((table['direction'] == d) & table['judgement_rank_intervention']
                                            .notna()).sum())}
              for d in DIRECTIONS}
    stamp = lambda p: {'path': str(p), 'mtime_utc': datetime.fromtimestamp(p.stat().st_mtime, timezone.utc)
                       .isoformat(), 'sha256': sha256(p)}
    prov = {'script': str(Path(__file__).resolve()), 'run_utc': datetime.now(timezone.utc).isoformat(),
            'repository_head': subprocess.run(['git', '-C', str(ROOT), 'show', '-s', '--format=%H', 'HEAD'], capture_output=True,
                                              text=True, check=True).stdout.strip(),
            'sources': {k: stamp(p) for k, p in SOURCES.items()},
            'figure_records_checked': {k: stamp(p) for k, p in FIGURE_RECORDS.items()},
            'modules_imported': {k: stamp(p) for k, p in MODULES.items()},
            'rows': len(table), 'counts': counts, 'checks_passed': checks}
    (OUT / 'fig3_rankings_provenance.json').write_text(json.dumps(prov, indent=1, ensure_ascii=False))
    print(f'{len(table)} rows -> {OUT / "fig3_rankings.csv"}')
    print(json.dumps(counts, indent=1))
    for c in checks:
        print('ok', c)


if __name__ == '__main__':
    main()
