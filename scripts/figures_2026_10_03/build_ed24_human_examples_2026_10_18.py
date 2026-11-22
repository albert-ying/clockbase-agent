"""Extended Data Figure 24 (revision): human comparisons of the verified ranking.

Five ranked human comparisons, chosen on 2026-10-18 among verified rows for having no published aging-clock study
(checked by literature search; atopic dermatitis, the second-ranked human DNA methylation row, was left out because a
2024 study reports the same cohort): VEXAS syndrome in whole blood, critical limb ischaemia in gastrocnemius muscle,
ataxia-telangiectasia in cerebellar vermis, gonadotropin stimulation in granulosa cells (the one human intervention
among the top-ranked human rows) and varicose veins in the vein wall.

Each panel shows the values the verifier tested: per-individual means where samples repeat an individual, and values
adjusted to the mean age of the comparison where the verifier fitted clock ~ arm + age (the difference between the
adjusted arm means equals the fitted arm coefficient). The builder stops unless the plotted values reproduce the
recorded mean difference of every row.

  python -B scripts/figures_2026_10_03/build_ed24_human_examples_2026_10_18.py --final final-2026-10-10 --out DIR
"""
import ast
import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from nature_figure_utils import (release_arguments, retro_style, trim_axes, add_panel_labels, save_figure,  # noqa: E402
                                 figsize_nature, set_domain, allow_dead_margin, RUST, GRAY, INK, TYPE)
from audit_system import paths  # noqa: E402
from audit_system import design_estimand_2026_10_10 as design  # noqa: E402
from audit_system.row_repairs_2026_09_28 import individuals_by_comparison  # noqa: E402

EXAMPLES = [  # (modality key, GEO folder, series, comparison, panel title, tissue)
    ('hs_rna', 'GEO_hs_rna', 'GSE160268', 'VEXAS vs Healthy (Whole blood)', 'VEXAS syndrome', 'whole blood'),
    ('hs_rna', 'GEO_hs_rna', 'GSE120642', 'Critical limb ischemia vs Healthy Adult', 'Critical limb ischaemia',
     'gastrocnemius'),
    ('hs_rna', 'GEO_hs_rna', 'GSE175776', 'Ataxia Telangiectasia vs Control', 'Ataxia-telangiectasia',
     'cerebellar vermis'),
    ('hs_rna', 'GEO_hs_rna', 'GSE124177', 'Gonadotropins-stimulated vs Natural cycle', 'Gonadotropin stimulation',
     'granulosa cells'),
    ('dnam', 'GEO_hs_dnam', 'GSE68319', 'Varicose Veins vs Normal Veins', 'Varicose veins', 'vein wall'),
]
CONTROL_LABEL = {'GSE124177': 'Natural\ncycle'}
CASE_LABEL = {'GSE160268': 'VEXAS', 'GSE120642': 'CLI', 'GSE175776': 'A-T', 'GSE124177': 'Stimulated',
              'GSE68319': 'Varicose'}
CLOCK_LABEL = {'Horvathv1': 'Horvath clock (years)'}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def as_list(value):
    return ast.literal_eval(value) if isinstance(value, str) else list(value)


def tested_values(final, individuals, modality, folder, gse, comparison, rank):
    """(case values, control values, record row, ranked row, method) of one comparison as the verifier tested it."""
    record = json.loads((final / 'public-export' / 'studies' / modality / f'{gse}.json').read_text())
    ranked = rank[(rank.gse_id == gse) & (rank.comparison_name == comparison)].iloc[0]
    row = [r for r in record['rows'] if r['comparison_name'] == comparison and r['clock'] == ranked.clock][0]
    table = pd.read_csv(paths.resolve(row['predictions']['path']), low_memory=False)
    values = pd.to_numeric(table[ranked.clock], errors='coerce')
    values.index = table['gsm_id'].astype(str)
    values = values[~values.index.duplicated(keep=False)].dropna()
    series = design.load_series(paths.GEO, folder, gse, True)
    method = (row.get('design') or {}).get('method')
    subjects = series['subjects']
    if (row.get('design') or {}).get('subject_source') == 'two-reader individuals audit':
        subjects = individuals[(gse, comparison, row['condition_name'])]
    case = design.unit_means(values, as_list(row['case_gsms']), subjects, series['ages'])
    control = design.unit_means(values, as_list(row['control_gsms']), subjects, series['ages'])
    a = np.array([v for v, _, _ in case.values()])
    b = np.array([v for v, _, _ in control.values()])
    if method == 'AGE_ADJUSTED':
        slope = row['design']['age_slope']
        units = list(case.values()) + list(control.values())
        mean_age = np.mean([age for _, age, _ in units])
        a = np.array([v - slope * (age - mean_age) for v, age, _ in case.values()])
        b = np.array([v - slope * (age - mean_age) for v, age, _ in control.values()])
    recorded = row['recomputed']['mean_diff']
    if abs((a.mean() - b.mean()) - recorded) > 1e-6:
        raise SystemExit(f'{gse}: plotted values give {a.mean() - b.mean():.6f}, the record {recorded:.6f}')
    return a, b, row, ranked, method


def strip(ax, values, y, color, rng):
    jitter = rng.uniform(-0.12, 0.12, len(values))
    ax.boxplot([values], positions=[y], widths=0.5, vert=False, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor='white', edgecolor=color, lw=0.7),
               medianprops=dict(color=color, lw=1.0), whiskerprops=dict(color=color, lw=0.6),
               capprops=dict(color=color, lw=0.6))
    ax.scatter(values, y + jitter, s=9, color=color, edgecolors='none', alpha=0.85, zorder=3, clip_on=False)


def sci(p):
    exponent = int(np.floor(np.log10(p)))
    mantissa = p / 10 ** exponent
    return f'{mantissa:.1f}×10$^{{{exponent}}}$' if exponent < -2 else f'{p:.3f}'


def main():
    final, out = release_arguments(__doc__)
    retro_style()
    rank = pd.read_csv(final / 'comparison-ranking' / 'ranked_comparisons_all.csv')
    individuals = individuals_by_comparison(paths.VERIFIER_ROOT / 'individuals-audit-001')
    fig, axd = plt.subplot_mosaic([['a', 'b'], ['c', 'd'], ['e', '.']], figsize=figsize_nature('double', aspect=0.52),
                                  layout='constrained')
    fig.get_layout_engine().set(w_pad=0.05, h_pad=0.05)
    axes = [axd[k] for k in 'abcde']
    rng = np.random.default_rng(20261018)
    facts = []
    for ax, (modality, folder, gse, comparison, title, tissue) in zip(axes, EXAMPLES):
        a, b, row, ranked, method = tested_values(final, individuals, modality, folder, gse, comparison, rank)
        strip(ax, b, 0, GRAY['500'], rng)
        strip(ax, a, 1, RUST['400'], rng)
        ax.set_yticks([0, 1])
        ax.set_yticklabels([CONTROL_LABEL.get(gse, 'Control'), CASE_LABEL[gse]])
        ax.set_ylim(-0.6, 1.6)
        xlabel = CLOCK_LABEL.get(ranked.clock, 'Clock output')
        if method == 'AGE_ADJUSTED':
            xlabel += ', adjusted for age'
        ax.set_xlabel(xlabel)
        unit = 'individuals' if method in ('UNIT_MEANS', 'AGE_ADJUSTED') else 'samples'
        ax.set_title(f"{title} ({tissue})\n{gse}; g = {ranked.hedges_g:.2f}; FDR = {sci(row['recomputed']['fdr_family'])}; "
                     f"{len(a)} vs {len(b)} {unit}", loc='left')
        facts.append({'gse': gse, 'comparison': comparison, 'tissue': tissue, 'clock': ranked.clock,
                      'modality': ranked.modality, 'method': method or 'WELCH', 'n_case': len(a), 'n_control': len(b),
                      'hedges_g': float(ranked.hedges_g), 'fdr': float(row['recomputed']['fdr_family']),
                      'mean_diff': float(row['recomputed']['mean_diff']), 'final_class': row['final_class'],
                      'unit': unit})
    trim_axes(fig)
    # Panels on the same clock scale share one domain, so they are compared by eye.
    for label in {ax.get_xlabel() for ax in axes}:
        group = [ax for ax in axes if ax.get_xlabel() == label]
        low = min(ax.get_xlim()[0] for ax in group); high = max(ax.get_xlim()[1] for ax in group)
        for ax in group:
            set_domain(ax, xlim=(low, high))
            ax.xaxis.set_major_locator(MaxNLocator(5))
            if len(group) > 1:
                allow_dead_margin(ax, 'x', reason='shares the domain of the other panels on this clock scale')
    add_panel_labels(fig, {k: axd[k] for k in 'abcde'})
    pdf = out / 'ed24_human_examples_2026_10_18.pdf'
    try:
        save_figure(fig, pdf, axes_dict={k: axd[k] for k in 'abcde'},
                    review_png=out / 'ed24_human_examples_2026_10_18_review.png',
                    footer=f'data: {final.name} (verified rows of the comparison ranking)')
    except Exception:
        fig.savefig(out / 'ed24_FAILED_preview.png', dpi=200)
        raise
    (out / 'provenance.json').write_text(json.dumps({
        'final': final.name, 'examples': facts, 'builder_sha256': sha256(__file__), 'output_sha256': sha256(pdf)},
        indent=2) + '\n')
    print(f'wrote {pdf}')


if __name__ == '__main__':
    main()
