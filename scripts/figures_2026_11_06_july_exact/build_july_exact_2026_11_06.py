#!/usr/bin/env python3
"""Adapt verified inputs and run the July figure generators.

Style source: the July figure scripts (*.py) and figures (figures/*.pdf).
Geometry, type, colours and annotation placement remain in the copied sources.
All outputs and installed dependencies are isolated under figures-rebuild-2026-11-06.
"""
from pathlib import Path
import argparse
import contextlib
import difflib
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
ORIGINAL = Path('/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base/Aging-AI-Agent')
OUT = ROOT / 'paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-06'
PREVIOUS = OUT.parent / 'figures-rebuild-2026-11-05'
WORK = OUT / 'july_exact'
sys.path.insert(0, str(WORK / 'dependencies'))
os.environ['MPLBACKEND'] = 'Agg'
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

GENERATORS = [
    'create_combined_modality_volcano.py',
    'create_combined_modality_top_hits_v2.py',
    'create_radar_plot_combined.py', 'create_ridge_plot_combined.py',
    'analysis_2_genetic_intervention_types.py',
    'analysis_3_pathway_enrichment.py', 'analysis_4_fda_approved_drugs.py',
    'create_geneage_waffle_plot_nature.py', 'create_drugage_waffle_plot_nature.py',
]
INPUTS = []


def record(path):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}')
    INPUTS.append({'path': str(path.resolve()), 'mtime_utc': datetime.fromtimestamp(
        path.stat().st_mtime, timezone.utc).isoformat(),
        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
    return path


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, default=lambda x: x.item() if hasattr(x, 'item') else str(x)) + '\n')


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def save_copy(name, text, original):
    dest = SCRIPTS / name
    if dest.exists():
        return  # Preserve the customized plotting source of the current package.
    dest.write_text(text)
    delta = ''.join(difflib.unified_diff(original.read_text().splitlines(True),
        text.splitlines(True), fromfile=str(original), tofile=str(dest), n=3))
    dest.with_suffix('.diff').write_text(delta)


def prepare_scripts():
    for name in GENERATORS:
        if (SCRIPTS / name).exists():
            continue
        src = record(ORIGINAL / name)
        text = src.read_text()
        if name == 'create_combined_modality_volcano.py':
            start = text.index('    ann_pool = df[df["fdr"] < 0.05]')
            end = text.index('    # Persist picked annotations for audit', start)
            # Only annotation membership changes. The budgets and placement below stay July's.
            text = text[:start] + '''    ann_pool = pd.read_csv("result/volcano_ranked_comparisons_july.csv", low_memory=False)
    ann_pool["_nlf"] = -np.log10(ann_pool["fdr"].clip(lower=1e-300, upper=1.0))
    pos_picks, neg_picks = [], []
    for cat, budget in ANNOTATION_COUNTS.items():
        cat_rows = ann_pool[ann_pool["condition_category"] == cat]
        pos_picks.append(cat_rows[cat_rows["log2FoldChange"] > 0].head(budget["positive"]))
        neg_picks.append(cat_rows[cat_rows["log2FoldChange"] < 0].head(budget["negative"]))
    ann_pos = pd.concat(pos_picks, ignore_index=True)
    ann_neg = pd.concat(neg_picks, ignore_index=True)
''' + text[end:]
            text = text.replace('ax.set_xlabel("Effect Size")', 'ax.set_xlabel("Cohen\'s d")')
            text = text.replace('    df["_nlf"] = -np.log10(df["fdr"])\n', '    df["_nlf"] = -np.log10(df["fdr"])\n    # Display range of the July panel: Cohen\'s d within +/-30 and -log10 FDR up to 15. The 0.4% of rows beyond\n    # are drawn at the boundary; no row is removed and the header totals count every plotted row.\n    XCAP, YCAP = 30.0, 15.0\n    df["log2FoldChange"] = df["log2FoldChange"].clip(-XCAP, XCAP)\n    df["_nlf"] = df["_nlf"].clip(upper=YCAP)\n')
            text = text.replace('    ann_pool["_nlf"] = -np.log10(ann_pool["fdr"].clip(lower=1e-300, upper=1.0))\n', '    ann_pool["_nlf"] = -np.log10(ann_pool["fdr"].clip(lower=1e-300, upper=1.0))\n    ann_pool["log2FoldChange"] = ann_pool["log2FoldChange"].clip(-30.0, 30.0)\n    ann_pool["_nlf"] = ann_pool["_nlf"].clip(upper=15.0)\n')
            text = text.replace('cbar.set_label("Effect Size",', 'cbar.set_label("Cohen\'s d",')
        elif name == 'create_combined_modality_top_hits_v2.py':
            start = text.index('    # Filter: significant + final_score')
            end = text.index('    # Figure: 2×2 of categories', start)
            text = text[:start] + '''    # Statistical category ranking supplied with verified confidence intervals.
    qual = pd.read_csv("result/figure3_category_top15_july.csv", low_memory=False)
    qual["composite"] = qual["rank_bound"].abs()
    picks = {}
    for cat in CATEGORIES:
        sub = qual[qual["condition_category"] == cat].copy()
        sub["display"] = sub.apply(make_label, axis=1)
        picks[cat] = sub.head(TOP_N)

''' + text[end:]
            text = text.replace('ax.set_xlabel("Score",', 'ax.set_xlabel("Hedges\' g (CI bound nearest zero)",')
            text = text.replace('anti-aging (log2FC < 0)', 'lower clock output')
            text = text.replace('pro-aging (log2FC > 0)', 'higher clock output')
            text = text.replace(', ranked by composite (−log10 p × scaled score)',
                                " ranked by Hedges' g (CI bound nearest zero)")
        elif name == 'create_ridge_plot_combined.py':
            text = text.replace('log2 Fold Change (Age Acceleration)', "Cohen's d")
            text = text.replace('No effect (log2FC = 0)', "No effect (Cohen's d = 0)")
            text = text.replace('Decelerating Aging (left)', 'Lower clock output (left)')
            text = text.replace('Accelerating Aging (right)', 'Higher clock output (right)')
        elif name == 'analysis_2_genetic_intervention_types.py':
            text = text.replace('Enrichment of age-lowering signatures', 'Enrichment of lower clock output')
            text = text.replace('No change (log2FC=0)', "No change (Cohen's d=0)")
            text = text.replace('f"Effect Size"', 'f"Cohen\'s d"')
        elif name == 'analysis_4_fda_approved_drugs.py':
            # Change displayed text only, preserving identifiers, matching, selection and drawing.
            text = text.replace('label="Age-lowering"', 'label="Lower clock output"')
            text = text.replace('label="Age-increasing"', 'label="Higher clock output"')
            text = text.replace('% age-reducing)', '% lower clock output)')
            text = text.replace('"← age-lowering"', '"← lower clock output"')
            text = text.replace('"age-increasing →"', '"higher clock output →"')
            text = text.replace('"No change (log2FC=0)"', '"No change (Cohen\'s d=0)"')
            text = text.replace('"Effect Size"', '"Cohen\'s d"')
            text = text.replace('FDA-approved\\nage-reducing drug signatures',
                                'FDA-approved\\nlower clock output')
        elif name.startswith('create_geneage_waffle'):
            text = text.replace(str(ORIGINAL / 'data/agent_result_genage.csv'), 'data/agent_result_genage.csv')
        elif name.startswith('create_drugage_waffle'):
            text = text.replace(str(ORIGINAL / 'data/agent_result_drugage_significant.csv'), 'data/agent_result_drugage_significant.csv')
        save_copy(name, text, src)
    for num in [3, 4]:
        name = f'assemble_figure{num}_vector.py'
        if (SCRIPTS / name).exists():
            continue
        src = record(ROOT / 'paper/legacy_reference/figure_rev1' / name)
        text = src.read_text()
        text = text.replace('ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))', f'ROOT = {str(WORK)!r}')
        text = text.replace('os.path.join(ROOT, "Aging-AI-Agent", "figures")', 'os.path.join(ROOT, "figures")')
        # Also isolate dormant deployment destinations so copies cannot overwrite manuscript files.
        text = text.replace(f'/tmp/fig{num}_build', str(WORK / f'assembly_build_{num}'))
        save_copy(name, text, src)


def adapt():
    source = record(PREVIOUS / 'data/combined_3modality_plot_data_verified_2026-11-05.csv')
    provenance = record(PREVIOUS / 'data/verified_plotting_row_provenance_2026-11-05.csv')
    df = pd.read_csv(source, low_memory=False)
    rowprov = pd.read_csv(provenance, low_memory=False)
    selected = pd.read_csv(record(Path(os.environ['CLOCKBASE_FINAL']) / 'public-export/atlas_v3_combined.csv'), low_memory=False)
    selected = selected[selected.v3_final_class.isin(['VERIFIED', 'PROVISIONAL'])]
    assert len(df) == len(rowprov) == len(selected)
    assert rowprov.atlas_row_ordinal.tolist() == selected.atlas_row_ordinal.tolist()
    assert np.allclose(rowprov.v3_fdr, selected.v3_fdr)
    assert rowprov.plot_row_index.tolist() == list(range(len(df)))
    assert np.allclose(df.effect_size, rowprov.v3_cohen_d)
    assert df.final_score.isna().all()
    # Normal is an input compatibility flag for an already verified inclusion set.
    # It does not claim that the historical outlier classifier was rerun.
    df['outlier_status'] = 'Normal'
    df.to_csv(WORK / 'figures/combined_3modality_plot_data.csv', index=False)
    mouse = df[df.modality == 'mouse RNA-seq'].copy()
    mouse.to_csv(WORK / 'result/main_results_result_mm_rna_outlier_marked.csv', index=False)
    keyed = df.assign(clock=rowprov.clock.values, verified_row_index=np.arange(len(df)))
    keys = ['modality', 'gse_id', 'comparison_name', 'condition_name', 'clock']
    rankroot = Path(os.environ['CLOCKBASE_FINAL']) / 'comparison-ranking'
    for fn in ['figure3_category_top15.csv', 'ranked_comparisons_all.csv',
               'table1_pro_aging.csv', 'table2_anti_aging.csv']:
        record(rankroot / fn)
    for fn, dest in [('figure3_category_top15.csv', 'figure3_category_top15_july.csv'),
                     ('ranked_comparisons_all.csv', 'volcano_ranked_comparisons_july.csv')]:
        rank = pd.read_csv(rankroot / fn, low_memory=False)
        rank['_source_order'] = np.arange(len(rank))
        # The ranking export retains a historical effect_size alias in some
        # rows. Reconstruct d only for identity checks from its corrected g.
        correction = 1 - 3 / (4 * (rank.v3_n_case + rank.v3_n_control - 2) - 1)
        rank['_rank_source_d'] = rank.hedges_g / correction
        rank['_rank_source_p'] = rank.v3_p
        # Source order is the published statistical ranking. Plot coordinates and
        # FDR stay the verified atlas values, not the ranking export's FDR aliases.
        additions = ['hedges_g', 'g_ci_low', 'g_ci_high', 'rank_bound', 'ci_level',
                     'selection_k', 'label', '_source_order', '_rank_source_d', '_rank_source_p']
        joined = rank[keys + additions].merge(keyed, on=keys, how='left', indicator=True)
        if not joined['_merge'].eq('both').all():
            raise ValueError(f'Ranking rows absent from verified table: {joined.loc[joined._merge != "both", keys].to_dict("records")}')
        joined = joined[np.isclose(joined._rank_source_d, joined.effect_size, rtol=1e-5, atol=1e-9)
                        & np.isclose(joined._rank_source_p, joined.pvalue, rtol=1e-10, atol=1e-300)]
        joined = joined.sort_values(['_source_order','verified_row_index']).drop_duplicates('_source_order')
        assert len(joined) == len(rank), f'Unresolved verified ranking join for {fn}'
        joined = joined.drop(columns='_merge')
        joined['intervention'] = joined['label']
        nearest = np.where(joined.g_ci_low > 0, joined.g_ci_low,
                           np.where(joined.g_ci_high < 0, -joined.g_ci_high, 0))
        assert np.allclose(nearest, joined.rank_bound.abs())
        assert (np.sign(joined.hedges_g) == np.sign(joined.log2FoldChange)).all()
        joined.to_csv(WORK / 'result' / dest, index=False)
    life_path = record(Path(os.environ.get('CLOCKBASE_META') or Path(os.environ['CLOCKBASE_VERIFIER_ROOT']) / 'intervention-meta-2026-10-17') / 'concordance/lifespan_concordance.csv')
    life = pd.read_csv(life_path)
    assert life.species.eq('mouse').all()
    life.to_csv(WORK / 'data/lifespan_concordance_source.csv', index=False)
    for source, prefix, outfile in [('GenAge', 'genage', 'agent_result_genage.csv'),
                                    ('DrugAge', 'drugage', 'agent_result_drugage_significant.csv')]:
        sub = life[life.source == source].copy()
        names = sub.canonical_target.copy()
        # Distinct manipulations can have opposite expected signs. A compound
        # key preserves the verified concordance unit through July's groupby.
        repeated = sub.canonical_target.duplicated(keep=False)
        names.loc[repeated] = names.loc[repeated] + ' ' + sub.loc[repeated, 'manipulation'].str.lower()
        assert not names.duplicated().any()
        compatible = pd.DataFrame({
            # +inf is an explicit pass-through sentinel, NOT a measured trust score.
            'trust_score': np.inf,
            f'{prefix}_annotated': True,
            f'{prefix}_matched_{"genes" if source == "GenAge" else "drugs"}': names,
            f'{prefix}_annotation': np.sign(sub.expected).astype(int),
            'inferred_influence': np.sign(sub.estimate).astype(int),
            'intervention': names,
        })
        assert ((compatible[f'{prefix}_annotation'] == compatible.inferred_influence).values == sub.concordant.values).all()
        compatible.to_csv(WORK / 'data' / outfile, index=False)
    for fn in ['KEGG_2019_Mouse.gmt', 'KEGG_2019_Human.gmt', 'opentargets_annotation.json',
               'panel_f_drug_matches.csv', 'panel_e_all_ora_results.csv', 'panel_e_gene_assignments.csv']:
        src = record(PREVIOUS / 'figure4' / fn)
        shutil.copy2(src, WORK / 'data' / fn)
    for fn in ['gene_reference_mouse.csv']:
        src = record(ROOT / 'result' / fn)
        shutil.copy2(src, WORK / 'result' / fn)
    write_json(WORK / 'data/input_provenance.json', INPUTS)
    write_json(WORK / 'data/adapter_rules.json', {
        'combined_rows': len(df), 'mouse_rows': len(mouse), 'significant_rows': int((df.fdr < .05).sum()),
        'outlier_status': 'Normal is a compatibility inclusion literal for all supplied inference-supporting rows. No July outlier labels are imported.',
        'final_score': 'Remains missing. The unchanged July normalizer uses its default 50 for point size; no empirical score is imputed.',
        'waffle_trust_score': '+inf is a compatibility pass-through sentinel, not an inferred or measured trust score. All supplied mouse concordance units are included.',
        'waffle_unit': 'canonical_target and manipulation; append manipulation only for targets with multiple rows to prevent the July gene-only majority vote from conflating opposite perturbations.',
        'waffle_unmatched': '0 in these matched-only inputs. The supplied concordance table contains no unmatched population; none is invented.',
        'ranking': 'Preserve source-file order, use absolute rank_bound as bar length; verified atlas Cohen d and BH-FDR supply scatter positions and direction.',
        'fdr': 'Input BH-FDR unchanged. July clips only volcano plotting values to [1e-300,1].',
    })


def run_script(name):
    env = os.environ.copy()
    env['PYTHONPATH'] = str(WORK / 'dependencies')
    env['MPLBACKEND'] = 'Agg'
    with (WORK / 'logs' / name.replace('.py', '.log')).open('w') as log:
        print('Running', name, flush=True)
        subprocess.run([sys.executable, str(SCRIPTS / 'run_generator_2026_11_06.py'), name], cwd=WORK,
                       env=env, stdout=log, stderr=subprocess.STDOUT, check=True)


def run_pathway():
    mod = load_module('july_pathway', SCRIPTS / 'analysis_3_pathway_enrichment.py')
    mouse = pd.read_csv(WORK / 'result/main_results_result_mm_rna_outlier_marked.csv')
    ref = pd.read_csv(WORK / 'result/gene_reference_mouse.csv')
    significant, background, gene_df = mod.extract_genes(mouse, set(ref.symbol.str.upper()))
    rows = []
    for line in (WORK / 'data/KEGG_2019_Mouse.gmt').read_text().splitlines():
        term, _, *genes = line.split('\t')
        members = {g.upper() for g in genes} & background
        overlap = significant & members
        if not overlap:
            continue
        a = len(overlap); b = len(significant) - a
        c = len(members) - a; d = len(background) - a - b - c
        odds, p = stats.fisher_exact([[a,b],[c,d]], alternative='greater')
        rows.append({'Pathway': term, 'Database': 'KEGG_2019_Mouse', 'N_Overlap': a,
            'Pathway_Size': len(members), 'Library_Pathway_Size': len(set(genes)),
            'P_value': p, 'Odds_Ratio': odds, 'Genes': ';'.join(sorted(overlap)),
            'a': a, 'b': b, 'c': c, 'd': d})
    all_results = pd.DataFrame(rows)
    all_results['FDR'] = fdrcorrection(all_results.P_value)[1]
    all_results['Neg_log10_FDR'] = -np.log10(all_results.FDR + 1e-300)
    all_results = all_results.sort_values('Neg_log10_FDR', ascending=False)
    all_results.to_csv(WORK / 'result/panel_e_mouse_all_ora_results.csv', index=False)
    # Match July's numeric overlap and pathway-size filter, before its plot filter.
    enrichment = all_results[(all_results.Pathway_Size < 100) & (all_results.N_Overlap >= 15)].copy()
    enrichment.to_csv(WORK / 'result/analysis_3_pathway_enrichment_ora.csv', index=False)
    write_json(WORK / 'result/panel_e_method.json', {
        'background_genes': sorted(background), 'query_genes': sorted(significant),
        'background_n': len(background), 'query_n': len(significant),
        'terms_with_overlap_tested': len(all_results), 'passing_july_size_overlap': len(enrichment),
        'passing_july_fdr': int((enrichment.FDR < .05).sum()),
        'method': 'Cached KEGG_2019_Mouse; official mouse reference; original extract_genes; observed-gene custom background; one-sided Fisher; BH across nonzero-overlap terms; July size<100 and overlap>=15; original create_figure applies FDR<.05 and cancer exclusion.',
    })
    mod.create_figure(enrichment, gene_df)


def run_fda():
    mod = load_module('july_fda', SCRIPTS / 'analysis_4_fda_approved_drugs.py')
    raw = pd.read_csv(WORK / 'data/panel_f_drug_matches.csv', low_memory=False)
    mouse = raw[(raw.modality == 'mouse RNA-seq') & (raw.fdr < .05) & raw.matched].copy()
    mouse['opentargets_info'] = mouse.apply(lambda r: {'name': r.matched_name, 'match_type': r.match_type}, axis=1)
    results = []
    for status in ['FDA_Approved', 'Phase_3', 'Phase_2', 'Phase_1', 'Preclinical']:
        sub = mouse[mouse.clinical_status == status]
        if sub.empty: continue
        results.append({'Clinical_Status':status, 'N_Drugs':len(sub),
            'Mean_Effect':sub.log2FoldChange.mean(), 'Median_Effect':sub.log2FoldChange.median(),
            'Prop_Anti_Aging':(sub.log2FoldChange < 0).mean(),
            'N_Anti_Aging':int((sub.log2FoldChange < 0).sum()),
            'N_Pro_Aging':int((sub.log2FoldChange > 0).sum()),
            'Strongest_Anti_Aging':sub.log2FoldChange.min()})
    results = pd.DataFrame(results).sort_values('Prop_Anti_Aging', ascending=False)
    results.to_csv(WORK / 'result/analysis_4_drug_clinical_status.csv', index=False)
    mouse.drop(columns='opentargets_info').to_csv(WORK / 'result/panel_f_mouse_drug_matches.csv', index=False)
    fda_lower = mouse[(mouse.clinical_status == 'FDA_Approved') & (mouse.log2FoldChange < 0)]
    baseline_rows = raw[(raw.modality == 'mouse RNA-seq') & (raw.fdr < .05)]
    # The original function also makes a density PDF in the isolated workspace.
    # A missing positive-drug list suppresses its unused fourth plot.
    mod.create_figure(mouse, results, pd.DataFrame(), fda_lower, None,
                     (baseline_rows.log2FoldChange < 0).mean())
    top = fda_lower[fda_lower.match_type.isin(['exact','synonym','tradename'])].sort_values('log2FoldChange').drop_duplicates('matched_name').head(15)
    top.drop(columns='opentargets_info').to_csv(WORK / 'result/panel_f_mouse_approved_lower_top15.csv', index=False)


def run_all():
    os.chdir(WORK)
    for name in ['create_combined_modality_top_hits_v2.py', 'create_combined_modality_volcano.py',
                 'create_radar_plot_combined.py', 'create_ridge_plot_combined.py',
                 'analysis_2_genetic_intervention_types.py']:
        run_script(name)
    # Each module executes in a fresh process, preserving July's independent
    # matplotlib/seaborn settings instead of leaking settings across panels.
    for part in ['pathway', 'fda']:
        with (WORK / 'logs' / f'{part}.log').open('w') as log:
            subprocess.run([sys.executable, str(SCRIPTS / 'run_generator_2026_11_06.py'), part],
                env={**os.environ, 'PYTHONPATH':str(WORK / 'dependencies')},
                stdout=log, stderr=subprocess.STDOUT, cwd=WORK, check=True)
    for name in GENERATORS[-2:]: run_script(name)
    figure4_numbers()


def figure4_numbers():
    """Export legend counts and bind every panel table to the current verified inputs."""
    final = Path(os.environ['CLOCKBASE_FINAL']).resolve()
    assert (final / 'reclassification_provenance.json').is_file(), final
    inputs = []

    def source(path):
        path = Path(path).resolve()
        try:
            inputs.append({'path': str(path), 'mtime_ns': path.stat().st_mtime_ns,
                           'sha256': hashlib.sha256(path.read_bytes()).hexdigest()})
        except OSError as exc:
            raise FileNotFoundError(f'Cannot proceed - required data file not found: {path}') from exc
        return path

    def csv(name):
        return pd.read_csv(source(WORK / name), low_memory=False)

    # The adapter receipt must identify this run and its still-identical source files.
    receipt = json.loads(source(WORK / 'data/input_provenance.json').read_text())
    required = [final / 'public-export/atlas_v3_combined.csv',
                final / 'intervention-meta-2026-10-17/concordance/lifespan_concordance.csv']
    for path in required:
        matches = [r for r in receipt if Path(r['path']).resolve() == path]
        assert len(matches) == 1 and matches[0]['sha256'] == hashlib.sha256(source(path).read_bytes()).hexdigest(), path
    frame = csv('figures/combined_3modality_plot_data.csv')
    atlas = pd.read_csv(required[0], low_memory=False)
    usable = atlas[atlas.v3_final_class.isin(['VERIFIED', 'PROVISIONAL'])].copy()
    keys = ['modality', 'gse_id', 'comparison_name', 'condition_name']
    classification = pd.read_json(source(final / 'public-export/classification/comparisons.jsonl'), lines=True)
    usable = usable.merge(classification[keys + ['category']], on=keys, how='left', validate='many_to_one', sort=False)
    # The plotting adapter uses the classification output and maps other categories to Other.
    usable['condition_category'] = usable.category.where(
        usable.category.isin(['Disease', 'Drug', 'Environment', 'Genetic']), 'Other')
    identity = ['modality', 'gse_id', 'comparison_name', 'condition_name', 'condition_category']
    pd.testing.assert_frame_equal(frame[identity].reset_index(drop=True), usable[identity].reset_index(drop=True))
    for plot_field, atlas_field in [('effect_size', 'v3_cohen_d'), ('fdr', 'v3_fdr'), ('pvalue', 'v3_p')]:
        assert np.allclose(frame[plot_field], usable[atlas_field], rtol=1e-12, atol=0, equal_nan=True), plot_field
    categories = csv('result/panel_a_category_counts.csv').set_index('condition_category')
    sig = frame[frame.fdr < .05]
    for category, row in categories.iterrows():
        assert row.all_rows == frame.condition_category.eq(category).sum()
        assert row.significant_rows == sig.condition_category.eq(category).sum()
    ridge = csv('result/panel_b_ridge_statistics.csv')
    ridge_record = json.loads(source(WORK / 'logs/create_ridge_plot_combined_figure_record.json').read_text())
    xlim = ridge_record['axes'][0]['xlim']
    outside = int(((sig.effect_size < xlim[0]) | (sig.effect_size > xlim[1])).sum())
    assert ridge.n.sum() == len(sig) and ridge.outside_display_n.sum() == outside
    genetic = csv('result/analysis_2_genetic_types.csv')
    fold_changes = csv('result/analysis_2_genetic_types_fold_change.csv')
    fisher = csv('result/analysis_2_genetic_pairwise_fold_change.csv')
    method = json.loads(source(WORK / 'result/panel_e_method.json').read_text())
    ora = csv('result/panel_e_mouse_all_ora_results.csv')
    plotted_pathways = csv('result/panel_e_plotted_pathways.csv')
    assert len(ora) == method['terms_with_overlap_tested']
    drugs = csv('result/panel_f_mouse_drug_matches.csv')
    phases = csv('result/analysis_4_drug_clinical_status.csv')
    assert phases.N_Drugs.sum() == len(drugs)
    approved = drugs[drugs.clinical_status == 'FDA_Approved']
    inset = csv('result/panel_f_mouse_approved_lower_top15.csv')
    life = csv('data/lifespan_concordance_source.csv')
    pd.testing.assert_frame_equal(life, pd.read_csv(required[1]))
    genes = life[life.source == 'GenAge']
    gene_names = genes.groupby('canonical_target').manipulation.nunique()
    mouse_genetic = sig[(sig.modality == 'mouse RNA-seq') & (sig.condition_category == 'Genetic')]
    # h-j have independent benchmark builders. Check their tables against the
    # benchmark scored on this final run before exposing their counts to the legend.
    bench = Path(os.environ['CLOCKBASE_BENCH_RESULTS'])
    scored = json.loads(source(bench / 'export_as_checks.json').read_text())
    assert scored['scored'] == final.name
    truth = json.loads(source(bench / 'ground_truth_counts.json').read_text())
    released = json.loads(source(bench / 'published_released_test.json').read_text())
    by_class = json.loads(source(bench / 'three_way_by_class.json').read_text())
    schematic = csv('result/figure4_benchmark_schematic_numbers.csv').set_index('quantity').value.astype(int).to_dict()
    expected_benchmark = {'truth_results': truth['results'], 'truth_papers': truth['papers'],
                          'truth_series': truth['series'], 'records': released['against_models_without_data']['records'],
                          'papers_scored': released['systems']['S3c']['clusters'],
                          'models': len([k for k in released['against_models_without_data']['systems'] if k.startswith('S4:')]),
                          **{name: by_class['systems']['S3c']['recall'][key][1] for name, key in
                             [('decreases', 'lower'), ('increases', 'higher'), ('not_significant', 'none')]}}
    for key, value in expected_benchmark.items():
        assert schematic[key] == value, f'Stale Figure 4 benchmark schematic: {key}'
    benchmark_bars = csv('result/figure4_benchmark_bars.csv')
    for row in benchmark_bars.itertuples():
        for key in ('lower', 'higher', 'none'):
            hits, total = by_class['systems'][row.source_key]['recall'][key]
            assert getattr(row, key + '_hits') == hits and getattr(row, key + '_total') == total
        if row.source_key != 'S0':
            expected_accuracy = released['against_models_without_data']['systems'][row.source_key]['call_balanced_accuracy']
            assert np.isclose(row.balanced_accuracy, expected_accuracy, rtol=1e-12)
    values = {'run': final.name, 'all_rows': len(frame), 'significant_rows': len(sig),
              'categories': categories.to_dict('index'),
              'ridge': {'axis_limit': xlim[1], 'outside_display_n': outside},
              'ridge_statistics': ridge.to_dict('records'),
              'genetic': {'significant': len(mouse_genetic), 'classified': int(genetic.N_Interventions.sum()),
                          'types': dict(zip(genetic.Type, genetic.N_Interventions)), 'type_count': len(genetic),
                          'statistics': fold_changes.to_dict('records'), 'min_pairwise_fisher_fdr': fisher.FDR.min(),
                          'baseline_lower_fraction': float(genetic.N_Anti_Aging.sum() / genetic.N_Interventions.sum())},
              'pathways': {**{k: method[k] for k in ('query_n', 'background_n', 'terms_with_overlap_tested', 'passing_july_fdr')},
                           'plotted': plotted_pathways.to_dict('records')},
              'drugs': {'matched_rows': len(drugs), 'approved_rows': len(approved),
                        'approved_lower': int(approved.log2FoldChange.lt(0).sum()),
                        'approved_higher': int(approved.log2FoldChange.gt(0).sum()),
                        'approved_binomial_p': stats.binomtest(int(approved.log2FoldChange.lt(0).sum()), len(approved), .5).pvalue,
                        'approved_drugs': approved.matched_name.nunique(), 'inset_drugs': len(inset)},
              'genage': {'squares': len(genes), 'genes': len(gene_names),
                         'both_directions': int(gene_names.gt(1).sum())},
              'benchmark': {key: schematic[key] for key in expected_benchmark},
              'sources': inputs}
    write_json(WORK / 'result/figure4_legend_numbers.json', values)
    print(json.dumps({k: v for k, v in values.items() if k != 'sources'}, indent=2))


def assemble():
    for num in [3,4]:
        mod = load_module(f'assembly_{num}', SCRIPTS / f'assemble_figure{num}_vector.py')
        mod.build(str(OUT / f'figure{num}/figure{num}.pdf'))
        mod.FIGSRC = str(ORIGINAL / 'figures')
        if num == 3:
            mod.VOLCANO = str(ORIGINAL / 'figures/volcano_combined_modalities.pdf')
            mod.BARS = str(ORIGINAL / 'figures/top_hits_barplot_combined_modalities_v2.pdf')
        mod.build(str(OUT / f'july_reference/figure{num}_july.pdf'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('step', choices=['prepare','run','pathway','fda','assemble','figure3','figure4','numbers'])
    step = parser.parse_args().step
    try:
        if step == 'prepare': prepare_scripts(); adapt()
        elif step == 'figure3':
            for name in ['create_combined_modality_top_hits_v2.py',
                         'create_combined_modality_volcano.py',
                         'create_figure3_replication_2026_11_07.py',
                         'create_figure3_judgement_2026_11_06.py']:
                run_script(name)
            mod = load_module('assembly_3', SCRIPTS / 'assemble_figure3_vector.py')
            mod.build(str(OUT / 'figure3/figure3.pdf'))
        elif step == 'figure4':  # assembly only: panels a-g from the July generators, h,i from ED23
            mod = load_module('assembly_4', SCRIPTS / 'assemble_figure4_vector.py')
            mod.build(str(OUT / 'figure4/figure4.pdf'))
        elif step == 'run': run_all()
        elif step == 'numbers': figure4_numbers()
        elif step == 'pathway': os.chdir(WORK); run_pathway()
        elif step == 'fda': os.chdir(WORK); run_fda()
        elif step == 'assemble': assemble()
    except Exception as exc:
        print(f'July reproduction failed during {step}: {exc}', file=sys.stderr)
        raise


if __name__ == '__main__': main()
