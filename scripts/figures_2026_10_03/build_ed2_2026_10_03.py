"""Extended Data Figure 2: cross-clock structure in whole blood.

Empirical inputs (no simulation):
  structure-002   audit_system/blood_clock_structure_2026_09_30.py: whole-blood samples aged 18-100 with all
                  17 clocks, historical predictions for non-empty studies (detached profiles removed) and
                  stream-repaired predictions for empty-matrix studies; within-study age residuals.
  grimage-age-input-002     audit_system/grimage_age_input_2026_09_30.py: GrimAge V1/V2 recomputed on the 51
                  stream-repaired studies (19,005 samples) with the historical inputs (Horvathv1 prediction
                  as age, sex labels read as male), with sex coded, and with recorded age and coded sex.
Panels: a, b residual correlation matrices (historical versus recorded-age GrimAge inputs); c the
Horvathv1-GrimAge residual correlation under the three inputs; d per-sample residuals (samples beyond the
99.9th percentile of absolute residual are not drawn; their count is written to provenance.json).

  MPLBACKEND=Agg python scripts/figures_2026_10_03/build_ed2_2026_10_03.py --final NAME --out DIR
"""
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from audit_system import paths  # noqa: E402
from nature_figure_utils import (release_arguments, retro_style, figsize_nature, trim_axes, add_panel_labels, add_colorbar,  # noqa: E402
                                 save_figure, set_domain, label_panel_group, RUST, BLUE, GRAY, SAGE, INK,
                                 CMAP_DIVERGE_CRISP)

RUN = paths.VERIFIER_ROOT / 'r6-6-blood'
COHORT = RUN / 'structure-002' / 'cohort.csv'
GRIMAGE = RUN / 'grimage-age-input-002' / 'grimage_age_input.csv'
INPUTS = RUN / 'structure-002'
ORDER = ['Horvathv1', 'Horvathv2', 'Hannum', 'Lin', 'Zhang_10', 'PhenoAge', 'GrimAgeV1', 'GrimAgeV2',
         'DunedinPoAm38', 'DunedinPACE']
LABELS = {'Horvathv1': 'Horvath', 'Horvathv2': 'Skin & blood', 'Hannum': 'Hannum', 'Lin': 'Lin',
          'Zhang_10': 'Zhang', 'PhenoAge': 'PhenoAge', 'GrimAgeV1': 'GrimAge', 'GrimAgeV2': 'GrimAge2',
          'DunedinPoAm38': 'DunedinPoAm', 'DunedinPACE': 'DunedinPACE'}
HIST = RUST['400']
CHRONO = BLUE['400']


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def within_study_residuals(frame, columns):
    out = frame[columns].copy()
    for _, index in frame.groupby('gse_id').groups.items():
        age = frame.loc[index, 'age'].to_numpy(float)
        design = np.column_stack([np.ones_like(age), age])
        values = frame.loc[index, columns].to_numpy(float)
        coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
        out.loc[index, columns] = values - design @ coefficients
    return out


def matrix_panel(ax, corr, title):
    """Lower triangle without the diagonal: rows are clocks 2..n, columns clocks 1..n-1."""
    rows, cols = ORDER[1:], ORDER[:-1]
    c = corr.loc[rows, cols].to_numpy()
    n = len(rows)
    masked = np.ma.masked_array(c, mask=np.triu(np.ones_like(c, dtype=bool), k=1))
    cmap = CMAP_DIVERGE_CRISP.copy()
    cmap.set_bad('white')
    mesh = ax.pcolormesh(np.arange(n + 1), np.arange(n + 1), masked[::-1], cmap=cmap, vmin=-1, vmax=1,
                         edgecolors='white', linewidth=0.5)
    # Outline the Horvath-GrimAge cells (rows GrimAgeV1/V2, column Horvathv1) and print their values.
    for row in ('GrimAgeV1', 'GrimAgeV2'):
        i, j = rows.index(row), cols.index('Horvathv1')
        ax.add_patch(plt.Rectangle((j, n - 1 - i), 1, 1, fill=False, edgecolor=INK, linewidth=0.9))
        ax.text(j + 0.5, n - 1 - i + 0.5, f"{c[i, j]:.2f}", ha='center', va='center', fontsize=6,
                color='white' if abs(c[i, j]) > 0.55 else INK)
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_xticklabels([LABELS[k] for k in cols], rotation=60, ha='right', rotation_mode='anchor')
    ax.set_yticks(np.arange(n) + 0.5)
    ax.set_yticklabels([LABELS[k] for k in rows][::-1])
    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    ax.set_title(title, fontsize=8, fontweight='normal', loc='left')
    return mesh


def main():
    retro_style('soft')
    structure = json.loads((INPUTS / 'summary.json').read_text())
    contrast = structure['structures']['grimage_age_input_contrast']
    hist = pd.read_csv(INPUTS / 'correlation_repaired_residual_spearman_historical_horvath_age_input.csv', index_col=0)
    chrono = pd.read_csv(INPUTS / 'correlation_repaired_residual_spearman_chronological_age_input.csv', index_col=0)

    cohort = pd.read_csv(COHORT, dtype={'gsm_id': str})
    grim = pd.read_csv(GRIMAGE, dtype={'gsm_id': str})
    repaired = cohort[cohort.source == 'stream_repaired'].merge(
        grim[['gse_id', 'gsm_id', 'GrimAgeV1_chronological']], on=['gse_id', 'gsm_id'])
    sizes = repaired.groupby('gse_id').size()
    repaired = repaired[repaired.gse_id.isin(sizes.index[sizes >= 20])].reset_index(drop=True)
    resid = within_study_residuals(repaired, ['Horvathv1', 'GrimAgeV1', 'GrimAgeV1_chronological'])

    # Two hand-set grids: the heatmaps' long row labels must not narrow the panels below them.
    fig = plt.figure(figsize=(7.2, 5.7))
    top = fig.add_gridspec(1, 2, left=0.115, right=0.925, top=0.955, bottom=0.515, wspace=0.42)
    bottom = fig.add_gridspec(1, 3, left=0.205, right=0.985, top=0.330, bottom=0.090, wspace=0.30,
                              width_ratios=[1.5, 1.0, 1.0])
    axd = {'a': fig.add_subplot(top[0, 0]), 'b': fig.add_subplot(top[0, 1]), 'c': fig.add_subplot(bottom[0, 0]),
           'd1': fig.add_subplot(bottom[0, 1]), 'd2': fig.add_subplot(bottom[0, 2])}

    matrix_panel(axd['a'], hist, 'GrimAge from Horvath age (historical)')
    mesh = matrix_panel(axd['b'], chrono, 'GrimAge from recorded age')
    add_colorbar(mesh, axd['b'], label='Spearman $\\rho$', ticks=[-1, 0, 1], inset=True)

    # c: Horvath-GrimAge residual correlation under the three inputs (rows), both GrimAge versions.
    ax = axd['c']
    variants = [('historical_horvath_age_input', 'Horvath age, all male'),
                ('sex_corrected_horvath_age_input', 'Horvath age, sex coded'),
                ('chronological_age_input', 'Recorded age, sex coded')]
    y = np.arange(len(variants))[::-1].astype(float)
    colors = [HIST, HIST, CHRONO]
    for version, marker, offset in (('GrimAgeV1', 'o', 0.12), ('GrimAgeV2', 's', -0.12)):
        values = [contrast[v]['horvath_grimage_residual_r'][version] for v, _ in variants]
        ax.plot(values, y + offset, color=GRAY['400'], lw=0.6, zorder=1)
        ax.scatter(values, y + offset, s=22, marker=marker, c=colors, edgecolors=INK, linewidths=0.4,
                   zorder=2, clip_on=False)
    ax.set_yticks(y)
    ax.set_yticklabels([label for _, label in variants])
    ax.set_xlabel('Residual $\\rho$ with Horvath')
    ax.set_title(f"GrimAge inputs ($n$ = {contrast['samples']:,}; {contrast['studies']} studies)", fontsize=7,
                 fontweight='normal', loc='left')
    handles = [plt.Line2D([], [], marker=m, ls='', color=GRAY['300'], markeredgecolor=INK, markeredgewidth=0.4,
                          markersize=4.5, label=l) for m, l in (('o', 'GrimAge'), ('s', 'GrimAge2'))]
    ax.legend(handles=handles, title='Clock', loc='lower right', frameon=False, fontsize=7, title_fontsize=7)

    # d: per-sample within-study residuals, historical versus recorded-age GrimAge (shared limits).
    limit = float(np.nanpercentile(np.abs(resid.to_numpy()), 99.9))
    shown = (resid.abs() <= limit).all(axis=1)
    hidden = int((~shown).sum())
    span = 1.08 * limit
    pooled = {'GrimAgeV1': contrast['historical_horvath_age_input']['horvath_grimage_residual_r']['GrimAgeV1'],
              'GrimAgeV1_chronological': contrast['chronological_age_input']['horvath_grimage_residual_r']['GrimAgeV1']}
    for key, column, color, title in (('d1', 'GrimAgeV1', HIST, 'Horvath age input'),
                                      ('d2', 'GrimAgeV1_chronological', CHRONO, 'Recorded age input')):
        ax = axd[key]
        ax.scatter(resid.loc[shown, 'Horvathv1'], resid.loc[shown, column], s=1.2, c=color, alpha=0.18,
                   edgecolors='none', rasterized=True, clip_on=False)
        # Same statistic as the matrices: within-study Spearman correlations pooled by Fisher z.
        ax.text(0.04, 0.96, f"$\\rho$ = {pooled[column]:.2f}", transform=ax.transAxes, ha='left', va='top', fontsize=7,
                color=INK)
        ax.set_xlabel('Horvath residual (years)')
        ax.set_ylabel('GrimAge residual (years)' if key == 'd1' else '')
        ax.set_title(title, fontsize=7, fontweight='normal', loc='left')
        ax.set_xlim(-span, span)
        ax.set_ylim(-span, span)

    trim_axes(fig, exclude=[axd['a'], axd['b']])
    for key in ('d1', 'd2'):
        set_domain(axd[key], xlim=(-span, span), ylim=(-span, span))
    add_panel_labels(fig, {'a': axd['a'], 'b': axd['b'], 'c': axd['c']})
    label_panel_group(fig, [axd['d1'], axd['d2']], 'd')
    provenance = {'cohort_csv': {'path': paths.portable(COHORT), 'sha256': sha256(COHORT)},
                  'grimage_csv': {'path': paths.portable(GRIMAGE), 'sha256': sha256(GRIMAGE)},
                  'structure_summary_sha256': sha256(INPUTS / 'summary.json'),
                  'script_sha256': sha256(__file__), 'scatter_samples_not_drawn_beyond_99_9th_percentile': hidden,
                  'scatter_limit_years': limit}
    (OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    save_figure(fig, str(OUT / 'ed2_blood_clock_structure_2026_10_03.pdf'), axes_dict=axd,
                footer='data: GEO whole-blood DNAm (ClockBase verifier v3); statistics computed in build_figure_2026_10_03.py',
                review_png=str(OUT / 'ed2_blood_clock_structure_2026_10_03_review.png'))
    plt.close(fig)
    print('saved')


if __name__ == '__main__':
    FINAL, OUT = release_arguments(__doc__)
    main()
