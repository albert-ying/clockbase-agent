"""Figure 4 (revision): intervention effects of the verified atlas across independent experiments and lifespan data.

Reads the intervention meta-analysis of a final run (audit_system/intervention_meta_2026_09_29.py; primary analysis,
series-as-unit run and the GenAge/DrugAge concordance under verifier/intervention-meta-2026-10-17/) and writes
figure4_verified_2026_10_14.pdf, a review PNG and provenance.json. Every value is read from those outputs. The
historical Figure 4 is archived as Supplementary Figure 2.

  python -B scripts/figures_2026_10_03/build_figure4_verified_2026_10_14.py --out "$TMPDIR/figure4"
"""
import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from nature_figure_utils import (release_arguments, retro_style, trim_axes, add_panel_labels, save_figure,  # noqa: E402
                                 set_domain, ref_line, pad_categorical_axis, figsize_nature, smart_annotate,
                                 ERRORBAR, BLUE, RUST, GRAY, INK, TYPE)
from audit_system import paths  # noqa: E402
from build_ed22_replication_2026_10_14 import (META_RUN, KEYS, DISPLAY, label, robust,  # noqa: E402
                                               with_category_suffix)

LOWER, HIGHER = BLUE['400'], RUST['400']
CAP = 12          # -log10 FDR shown up to this value; larger values are drawn at the cap and labelled
XCLIP = 6         # panel d: intervals are drawn within +-XCLIP; a clipped end carries an arrowhead
RENAME = {'Tcpobop': 'TCPOBOP', 'Dextran sodium sulfate': 'DSS', 'Lipopolysaccharide': 'LPS'}
# Panel a names the interventions the Results discuss there, one row each (the row with the most projects).
# Rapamycin is shown project by project in Extended Data Figure 22f and TCPOBOP in its panel d.
NAMED = ['lipopolysaccharide', 'bleomycin', 'sleep deprivation', 'plx5622', 'hyperoxia',
         'dextran sodium sulfate']


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    final, out = release_arguments(__doc__)
    meta = paths.VERIFIER_ROOT / META_RUN
    primary = with_category_suffix(pd.read_csv(meta / 'primary' / 'intervention_meta.csv'))
    by_series = pd.read_csv(meta / 'by-series' / 'intervention_meta.csv')
    lifespan = with_category_suffix(pd.read_csv(meta / 'concordance' / 'lifespan_concordance.csv'))
    summary = json.loads((meta / 'concordance' / 'summary.json').read_text())
    retro_style()
    replicated = primary[primary.k >= 2].copy()
    replicated['robust'] = (replicated.fdr < 0.05) & (replicated.rank_bound > 0)
    fig, axd = plt.subplot_mosaic([['a', 'a', 'd'], ['b', 'c', 'd']], figsize=figsize_nature('double', aspect=0.72),
                                  layout='constrained',
                                  gridspec_kw=dict(width_ratios=[0.95, 1.1, 1.15], height_ratios=[1.3, 0.6]))
    fig.get_layout_engine().set(h_pad=0.05, w_pad=0.05)

    # a: every intervention with results from at least two independent projects
    ax = axd['a']
    y = (-np.log10(replicated.fdr.clip(lower=1e-300))).clip(upper=CAP)
    colors = np.where(~replicated.robust, GRAY['300'], np.where(replicated.estimate < 0, LOWER, HIGHER))
    size = 6 + 5 * np.sqrt(replicated.k)
    order = np.argsort(replicated.robust.values)  # robust points on top
    ax.scatter(replicated.estimate.values[order], y.values[order], s=size.values[order], c=colors[order],
               edgecolors='none', alpha=0.85, clip_on=False, zorder=2)
    ax.set_xlabel("Pooled Hedges' g across projects")
    ax.set_ylabel('$-$log$_{10}$ BH-FDR')
    names = replicated[replicated.canonical_target.isin(NAMED) & (replicated.species == 'mouse')]
    names = names.sort_values('k', ascending=False).drop_duplicates(['species', 'canonical_target'])
    names = with_category_suffix(names.drop(columns='category_suffix'))  # suffix only where repeated in this panel
    points = list(zip(names.estimate, (-np.log10(names.fdr.clip(lower=1e-300))).clip(upper=CAP)))
    labels = []
    for r in names.itertuples():
        text = label(r).replace(' (mm)', '')
        text = RENAME.get(text, text)
        if -np.log10(r.fdr) > CAP:
            text += f' (FDR 10$^{{{int(np.floor(np.log10(r.fdr)))}}}$)'
        labels.append(text)
    n_lower = int((replicated.robust & (replicated.estimate < 0)).sum())
    n_higher = int((replicated.robust & (replicated.estimate > 0)).sum())
    ax.set_title(f'Interventions with $\\geq$2 independent projects ($n$ = {len(replicated):,})', loc='left')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc=HIGHER, mec='none', ms=4.5, label=f'Robust, higher ({n_higher})'),
                       Line2D([], [], marker='o', ls='', mfc=LOWER, mec='none', ms=4.5, label=f'Robust, lower ({n_lower})'),
                       Line2D([], [], marker='o', ls='', mfc=GRAY['300'], mec='none', ms=4.5, label='Not robust')],
              loc='upper right', frameon=False, fontsize=TYPE['legend'], handletextpad=0.3, borderaxespad=0.2)

    # b: what counting GEO series as independent studies does
    ax = axd['b']
    groups = [('2+ units', int((by_series.k >= 2).sum()), int((primary.k >= 2).sum())),
              ('Robust,\nhigher', int((robust(by_series).estimate > 0).sum()), n_higher),
              ('Robust,\nlower', int((robust(by_series).estimate < 0).sum()), n_lower)]
    yy = np.arange(len(groups))[::-1]
    ax.barh(yy + 0.19, [g[1] - 1 for g in groups], height=0.36, left=1, color=GRAY['200'], edgecolor=GRAY['400'],
            lw=0.5, label='Series')
    ax.barh(yy - 0.19, [g[2] - 1 for g in groups], height=0.36, left=1, color=BLUE['200'], edgecolor=BLUE['400'],
            lw=0.5, label='Projects')
    for yi, (_, s_count, p_count) in zip(yy, groups):
        for value, offset in ((s_count, 0.19), (p_count, -0.19)):
            ax.text(value * 1.12, yi + offset, f'{value:,}', ha='left', va='center', fontsize=TYPE['annotation'], color=INK)
    ax.set_xscale('log')
    ax.set_yticks(yy)
    ax.set_yticklabels([g[0] for g in groups])
    ax.set_xlabel('Interventions')
    ax.set_title('Counting unit', loc='left')
    ax.legend(loc='lower right', frameon=False, fontsize=TYPE['legend'], handlelength=1.0, borderaxespad=0.1)

    # c: agreement with lifespan databases, by the direction the lifespan effect predicts
    ax = axd['c']
    post = summary['post_hoc']['all matched']
    rows = [('Shortens\nlifespan', post['concordant_when_expected_higher'], post['expected_higher'], HIGHER),
            ('Extends\nlifespan', post['concordant_when_expected_lower'], post['expected_lower'], LOWER)]
    yy = np.arange(len(rows))[::-1]
    for yi, (name, k, n, color) in zip(yy, rows):
        ax.barh(yi, k / n * 100, height=0.42, color=color, alpha=0.35, edgecolor=color, lw=0.6)
        ax.barh(yi, 100 - k / n * 100, left=k / n * 100, height=0.42, color='white', edgecolor=GRAY['400'], lw=0.5)
        ax.text(k / n * 100 / 2, yi, f'{k} / {n}', ha='center', va='center', fontsize=TYPE['annotation'], color=INK)
    ax.set_yticks(yy)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel('With the expected direction (%)')
    ax.set_title(f"GenAge and DrugAge ($n$ = {summary['sets']['all matched']['both']['n']})", loc='left')

    # d: matched interventions with at least two projects, by expected direction
    ax = axd['d']
    part = lifespan[lifespan.k >= 2].copy()
    part = part.sort_values(['expected', 'estimate'], ascending=[False, True])
    y = np.arange(len(part))[::-1]
    for yi, r in zip(y, part.itertuples()):
        color = HIGHER if r.expected > 0 else LOWER
        low, high = max(r.ci_low, -XCLIP), min(r.ci_high, XCLIP)
        ax.plot([low, high], [yi, yi], color=color, lw=0.8, solid_capstyle='butt', zorder=1)
        for end, clipped, marker in ((low, r.ci_low < -XCLIP, '<'), (high, r.ci_high > XCLIP, '>')):
            if clipped:
                ax.plot([end], [yi], marker=marker, ms=3, color=color, zorder=2, clip_on=False)
        ax.scatter([r.estimate], [yi], s=16, marker='o', facecolor=color if r.concordant else 'white', edgecolor=color,
                   lw=0.8, zorder=3, clip_on=False)
    names = []
    for r in part.itertuples():
        base = label(r).replace(' (mm)', '')
        names.append(f'{base} ({r.k})')
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    for tick, r in zip(ax.get_yticklabels(), part.itertuples()):
        if r.source == 'GenAge':
            tick.set_fontstyle('italic')
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=3)
    ax.set_xlabel("Pooled Hedges' g (95% CI)")
    ax.set_title('Matched, $\\geq$2 projects', loc='left')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', mfc=HIGHER, mec=HIGHER, ms=4, label='Higher'),
                       Line2D([], [], marker='o', ls='', mfc=LOWER, mec=LOWER, ms=4, label='Lower'),
                       Line2D([], [], marker='o', ls='', mfc='white', mec=GRAY['500'], ms=4, label='Opposite')],
              title='Expected', title_fontsize=TYPE['legend'], loc='upper right', frameon=False,
              fontsize=TYPE['legend'], handletextpad=0.3, borderaxespad=0.1)
    ax.xaxis.set_major_locator(MaxNLocator(5, integer=True))
    ax.set_xlim(-XCLIP, XCLIP)

    trim_axes(fig, exclude=[axd['c']])
    placed = smart_annotate(axd['a'], labels, points, mode='near')
    shown = {t.get_text() for item in placed for t in (item if isinstance(item, tuple) else (item,))
             if hasattr(t, 'get_text')}
    missing = [t for t in labels if t not in shown]
    if missing:  # a named intervention must never drop silently from panel a
        raise SystemExit(f'panel a could not place {missing}')
    ref_line(axd['a'], x=0)
    ref_line(axd['d'], x=0)
    set_domain(axd['d'], xlim=(-XCLIP, XCLIP))
    set_domain(axd['a'], ylim=(-0.5, CAP + 0.7))
    pad_categorical_axis(axd['d'], 'y', pad=0.6)
    pad_categorical_axis(axd['b'], 'y', pad=0.45)
    pad_categorical_axis(axd['c'], 'y', pad=0.45)
    set_domain(axd['b'], xlim=(1, 10 ** np.ceil(np.log10(groups[0][1] * 2.5))))
    add_panel_labels(fig, axd)
    pdf = out / 'figure4_verified_2026_10_14.pdf'
    try:
        save_figure(fig, pdf, axes_dict=axd, review_png=out / 'figure4_verified_2026_10_14_review.png',
                    footer=f'data: {final.name} and {META_RUN} (meta-analysis by independent project; GenAge build 21, DrugAge build 5)')
    except Exception:
        fig.savefig(out / 'figure4_FAILED_preview.png', dpi=200)
        raise
    inputs = [meta / 'primary' / 'intervention_meta.csv', meta / 'by-series' / 'intervention_meta.csv',
              meta / 'concordance' / 'lifespan_concordance.csv', meta / 'concordance' / 'summary.json']
    (out / 'provenance.json').write_text(json.dumps({
        'final': final.name, 'meta': META_RUN,
        'inputs_sha256': {str(p.relative_to(paths.VERIFIER_ROOT)): sha256(p) for p in inputs},
        'replicated': len(replicated), 'robust_higher': n_higher, 'robust_lower': n_lower, 'panel_d_rows': len(part),
        'builder_sha256': sha256(__file__), 'output_sha256': sha256(pdf)}, indent=2) + '\n')
    print(f'wrote {pdf}')


if __name__ == '__main__':
    main()
