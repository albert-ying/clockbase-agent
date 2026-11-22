"""Extended Data Figure 22: replication of intervention effects across independent projects.

Reads the intervention meta-analysis of a final run (audit_system/intervention_meta_2026_09_29.py; primary analysis
and three sensitivity runs under verifier/intervention-meta-2026-10-17/) and writes
ed22_replication_2026_10_14.pdf, a review PNG and provenance.json. Every value is read from those outputs.

  python -B scripts/figures_2026_10_03/build_ed22_replication_2026_10_14.py --out "$TMPDIR/ed22"
"""
import hashlib
import json
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Polygon

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from nature_figure_utils import (release_arguments, retro_style, trim_axes, add_panel_labels, save_figure,  # noqa: E402
                                 set_domain, ref_line, pad_categorical_axis, nature_legend, figsize_nature,
                                 allow_dead_margin,
                                 ERRORBAR, BLUE, RUST, GRAY, INK, TYPE)
from audit_system import paths  # noqa: E402

META_RUN = 'intervention-meta-2026-10-17'
RUNS = {'primary': 'Primary', 'by-series': 'Series as unit', 'all-contexts': 'All contexts',
        'no-contact-link': 'No contact link'}
SHORT = {'primary': 'Primary', 'by-series': 'Series', 'all-contexts': 'All contexts', 'no-contact-link': 'No contact',
         'hartung-knapp': 'H–K'}
TOP = 12           # rows shown per direction in the forests
REFERENCE = ('mouse', 'rapamycin', 'EXPOSURE', 'EXPOSURE', 'Drug')   # an established lifespan-extending intervention
KEYS = ['species', 'canonical_target', 'intervention_class', 'manipulation', 'condition_category']  # one row per intervention and atlas category since the 2026-10-15 correction
LARGE_ARM = 120    # arms above this size are usually single cells (the individuals audit does not read them)
# Display names for model-written targets whose text is garbled or long; the data keep the original strings.
DISPLAY = {'3;4 dihydrocollidine;5 diethoxycarbonyl 1': 'DDC', "theiler's murine encephalomyelitis virus": "Theiler's virus",
           'respiratory syncytial virus': 'RSV', 'hepatitis virus': 'Mouse hepatitis virus', 'escherichia coli': 'E. coli',
           'lipopolysaccharide': 'LPS', 'dextran sodium sulfate': 'DSS', 'anti cd40 antibody': 'Anti-CD40',
           'heligmosomoides polygyrus': 'H. polygyrus', 'conduritol b epoxide': 'Conduritol B epoxide'}
LOWER, HIGHER = BLUE['400'], RUST['400']


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def label(row):
    target = str(row.canonical_target)
    if target in DISPLAY:
        name = DISPLAY[target]
    elif row.intervention_class == 'GENE':
        name = target.upper() if row.species == 'human' else target[:1].upper() + target[1:]
        name += {'DECREASE': ' loss', 'INCREASE': ' gain'}.get(row.manipulation, '')
    elif re.match(r'^[a-z]{2,5}[ -]?\d', target):
        name = target.upper()
    else:
        name = target[:1].upper() + target[1:]
    flag = '†' if getattr(row, 'large_arm', False) else ''
    # Since the 2026-10-15 correction one target can hold a pooled row per atlas category; name the category so
    # that two points of the same intervention are not read as one plotted twice.
    category = getattr(row, 'category_suffix', '')
    return f"{name}{flag}{category} ({'mm' if row.species == 'mouse' else 'hs'})"


def with_category_suffix(table, column='canonical_target'):
    """Add `category_suffix` to targets that appear under more than one atlas category *in these rows*.

    Applied to the rows a panel draws, so a label carries the category only where two points of the same
    intervention would otherwise be read as one plotted twice.
    """
    table = table.copy()
    if 'condition_category' not in table.columns or table.empty:
        table['category_suffix'] = ''
        return table
    repeated = table.groupby(['species', column]).condition_category.transform('nunique') > 1
    table['category_suffix'] = np.where(repeated, ', ' + table.condition_category.astype(str).str.lower(), '')
    return table


def robust(table):
    return table[(table.k >= 2) & (table.fdr < 0.05) & (table.rank_bound > 0)]


def load(meta):
    tables = {name: pd.read_csv(meta / name / 'intervention_meta.csv') for name in RUNS}
    units = pd.read_csv(meta / 'primary' / 'unit_effects.csv')
    projects = json.loads((meta / 'projects.json').read_text())
    effects = pd.read_csv(meta / 'primary' / 'meta_input_effects.csv', low_memory=False)
    large = effects[effects[['v3_n_case', 'v3_n_control']].max(axis=1) > LARGE_ARM].groupby(KEYS).size()
    primary = tables['primary'].set_index(KEYS)
    primary['large_arm'] = primary.index.isin(large.index)
    tables['primary'] = primary.reset_index()
    return tables, units, projects


def forest(ax, rows, color):
    rows = with_category_suffix(rows)
    y = np.arange(len(rows))[::-1]
    ax.errorbar(rows.estimate, y, xerr=[rows.estimate - rows.ci_low, rows.ci_high - rows.estimate], fmt='o',
                color=color, ms=3.2, mec=color, mfc=color, ecolor=color, **ERRORBAR['forest'])
    ax.set_yticks(y)
    ax.set_yticklabels([label(r) for r in rows.itertuples()])
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=3)
    ax.set_xlabel("Pooled Hedges' g (95% CI)")
    # The projects (series) column sits inside the axes, at its right edge; the domain is widened for it after the trim.
    for yy, r in zip(y, rows.itertuples()):
        ax.text(0.99, yy, f'{r.k} ({r.n_series})', transform=ax.get_yaxis_transform(), ha='right', va='center',
                fontsize=TYPE['annotation'], color=INK)
    ax.text(0.99, len(rows) - 0.55, 'Projects\n(series)', transform=ax.get_yaxis_transform(), ha='right',
            va='bottom', fontsize=TYPE['annotation'], color=INK, linespacing=1.0)
    ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
    ax._data_x = (float(min(rows.ci_low.min(), 0)), float(max(rows.ci_high.max(), 0)))
    return y


def lab(record):
    """Surname of the GEO contact of a series (the project's laboratory), or '' when the record has none."""
    if not record or not record.get('contact'):
        return ''
    return record['contact'].split(' | ')[0].split()[-1].title()


def per_project(ax, units, table, projects, key, color):
    part = units[(units[KEYS] == pd.Series(dict(zip(KEYS, key)))).all(axis=1)].copy()
    row = table[(table[KEYS] == pd.Series(dict(zip(KEYS, key)))).all(axis=1)].iloc[0]
    part = part.sort_values('g')
    y = np.arange(len(part))[::-1] + 1.4
    half = 1.959964 * np.sqrt(part['var'])
    ax.errorbar(part.g, y, xerr=half, fmt='s', color=color, ms=3.0, ecolor=color, **ERRORBAR['forest'])
    names = []
    for r in part.itertuples():
        series = r.series.split(';')
        who = lab(projects['records'].get(series[0]))
        more = f' +{len(series) - 1}' if len(series) > 1 else ''
        names.append(f"{who + ', ' if who else ''}{series[0]}{more}")
    m, lo, hi = row.estimate, row.ci_low, row.ci_high
    ax.add_patch(Polygon([[lo, 0], [m, 0.28], [hi, 0], [m, -0.28]], closed=True, facecolor=color, edgecolor=INK,
                         lw=0.6, alpha=0.9))
    ax.set_yticks(list(y) + [0])
    ax.set_yticklabels(names + [f'Pooled, {row.k} projects'])
    for t in ax.get_yticklabels():
        t.set_linespacing(1.0)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=3)
    ax.set_xlabel("Hedges' g (95% CI)")
    ax.set_title(f'{label(row)}', loc='left')
    ax._extent = (float(min((part.g - half).min(), lo)), float(max((part.g + half).max(), hi)))
    return part


def main():
    final, out = release_arguments(__doc__)
    meta = paths.VERIFIER_ROOT / META_RUN
    tables, units, projects = load(meta)
    retro_style()
    primary = tables['primary']
    robust_lower = robust(primary)[robust(primary).estimate < 0]
    robust_higher = robust(primary)[robust(primary).estimate > 0]
    lower, higher = robust_lower.head(TOP), robust_higher.head(TOP)

    fig, axd = plt.subplot_mosaic(
        [['a', 'b', 'd', 'd'], ['c', 'c', 'e', 'e'], ['c', 'c', 'f', 'f']],
        figsize=figsize_nature('double', aspect=0.79), layout='constrained',
        gridspec_kw=dict(width_ratios=[2.35, 1.8, 1.55, 1.5], height_ratios=[0.84, 1.45, 1.6]))
    fig.get_layout_engine().set(h_pad=0.04, w_pad=0.04)

    # a, c: robust replicated interventions, lower and higher clock outputs, ranked by the bound nearest zero
    ya = forest(axd['a'], lower, LOWER)
    axd['a'].set_title('Lower outputs', loc='left')
    forest(axd['c'], higher, HIGHER)
    axd['c'].set_title('Higher outputs', loc='left')

    # b: the rows of a under each specification (filled: robust; open: not robust; dash: under two units)
    ax = axd['b']
    columns = list(RUNS) + ['hartung-knapp']
    for j, name in enumerate(columns):
        table = tables['primary' if name == 'hartung-knapp' else name]
        for yy, r in zip(ya, lower.itertuples()):
            hit = table[(table[KEYS] == pd.Series({k: getattr(r, k) for k in KEYS})).all(axis=1)]
            if hit.empty or hit.iloc[0].k < 2:
                ax.plot([j - 0.15, j + 0.15], [yy, yy], color=GRAY['400'], lw=0.8)
                continue
            h = hit.iloc[0]
            ok = (h.ci_high_hk < 0) if name == 'hartung-knapp' else (h.fdr < 0.05 and h.rank_bound > 0)
            ax.scatter([j], [yy], s=24, marker='o', facecolor=LOWER if ok else 'white', edgecolor=LOWER, lw=0.8,
                       clip_on=False, zorder=3)
    ax.set_xticks(range(len(columns)))
    ax.set_xticklabels([SHORT[c] for c in columns], rotation=90, ha='center', va='bottom')
    ax.xaxis.set_ticks_position('top')
    ax.tick_params(axis='x', pad=2)
    ax.set_yticks(ya)
    ax.set_yticklabels([])
    ax.set_xlim(-0.5, len(columns) - 0.5)
    for side in ('left', 'bottom'):
        ax.spines[side].set_visible(False)
    ax.tick_params(length=0)
    handles = [Line2D([], [], marker='o', ls='', mfc=LOWER, mec=LOWER, ms=4.5, label='Robust'),
               Line2D([], [], marker='o', ls='', mfc='white', mec=LOWER, ms=4.5, label='Not robust'),
               Line2D([], [], ls='-', color=GRAY['400'], lw=0.8, label='< 2 units')]
    ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(-0.05, -0.02), frameon=False, ncol=1,
              fontsize=TYPE['legend'], handlelength=1.2, borderaxespad=0.0, labelspacing=0.25)

    # d: what counting series as independent studies does to the counts
    ax = axd['d']
    series_table = tables['by-series']
    groups = [('Two or more units', (series_table.k >= 2).sum(), (primary.k >= 2).sum()),
              ('Robust, lower', len(robust(series_table)[robust(series_table).estimate < 0]), len(robust_lower)),
              ('Robust, higher', len(robust(series_table)[robust(series_table).estimate > 0]), len(robust_higher))]
    y = np.arange(len(groups))[::-1]
    ax.barh(y + 0.19, [g[1] - 1 for g in groups], height=0.36, left=1, color=GRAY['200'], edgecolor=GRAY['400'],
            lw=0.5, label='Series')  # bars run from 1 to the count on the log axis
    ax.barh(y - 0.19, [g[2] - 1 for g in groups], height=0.36, left=1, color=BLUE['200'], edgecolor=BLUE['400'],
            lw=0.5, label='Projects')
    for yi, (_, s_count, p_count) in zip(y, groups):
        for value, offset in ((s_count, 0.19), (p_count, -0.19)):
            ax.text(value * 1.12, yi + offset, f'{value:,}', ha='left', va='center', fontsize=TYPE['annotation'],
                    color=INK)
    ax.set_xscale('log')
    ax.set_yticks(y)
    ax.set_yticklabels([g[0].replace('Two or more units', '2+ units') for g in groups])
    ax.set_xlabel('Interventions')
    ax.legend(loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=2, frameon=False, fontsize=TYPE['legend'],
              handlelength=1.0, columnspacing=0.8, borderaxespad=0.1)

    # e: the robust intervention with lower outputs that has the most projects; f: the reference intervention
    most = robust_lower.sort_values(['k', 'rank_bound'], ascending=False).iloc[0]
    follow = per_project(axd['e'], units, primary, projects, tuple(most[KEYS]), LOWER)
    reference = per_project(axd['f'], units, primary, projects, REFERENCE, LOWER)

    trim_axes(fig, exclude=[axd['b']])
    allow_dead_margin(axd['e'], 'x', reason='shares its domain with panel f so the two are compared by eye')
    for key in ('a', 'c', 'e', 'f'):
        pad_categorical_axis(axd[key], 'y', pad=0.45)
    low = min(axd['e']._extent[0], axd['f']._extent[0])
    high = max(axd['e']._extent[1], axd['f']._extent[1])
    for key in ('e', 'f'):
        set_domain(axd[key], xlim=(np.floor(low), np.ceil(high)))
    for key in ('a', 'c'):  # leave the right 24% of the forest for its projects (series) column
        lo, hi = axd[key]._data_x
        lo, hi = np.floor(lo), np.ceil(hi)
        set_domain(axd[key], xlim=(lo, lo + (hi - lo) / 0.76))
        ticks = [v for v in axd[key].get_xticks() if lo <= v <= hi]
        axd[key].set_xticks(ticks)
        axd[key].spines['bottom'].set_bounds(ticks[0], ticks[-1])
    for key in ('a', 'c', 'e', 'f'):
        ref_line(axd[key], x=0)
    axd['b'].set_ylim(axd['a'].get_ylim())
    set_domain(axd['d'], xlim=(1, 10 ** np.ceil(np.log10(max(g[1] for g in groups) * 2.5))))
    pad_categorical_axis(axd['d'], 'y', pad=0.45)
    add_panel_labels(fig, axd)
    pdf = out / 'ed22_replication_2026_10_14.pdf'
    try:
        save_figure(fig, pdf, axes_dict=axd, review_png=out / 'ed22_replication_2026_10_14_review.png',
                    footer=f'data: {final.name} and {META_RUN} (intervention meta-analysis by independent project)')
    except Exception:
        fig.savefig(out / 'ed22_FAILED_preview.png', dpi=200)  # for the visual review of a draft that failed a gate
        raise
    inputs = [meta / name / 'intervention_meta.csv' for name in RUNS] + [meta / 'primary' / 'unit_effects.csv',
                                                                        meta / 'primary' / 'meta_input_effects.csv',
                                                                        meta / 'projects.json']
    (out / 'provenance.json').write_text(json.dumps({
        'final': final.name, 'meta': META_RUN, 'inputs_sha256': {str(p.relative_to(paths.VERIFIER_ROOT)): sha256(p)
                                                                for p in inputs},
        'rows_lower': len(lower), 'rows_higher': len(higher), 'robust_lower': len(robust_lower),
        'robust_higher': len(robust_higher), 'panel_e': label(most), 'panel_e_units': len(follow),
        'reference_units': len(reference), 'builder_sha256': sha256(__file__), 'output_sha256': sha256(pdf)},
        indent=2) + '\n')
    print(f'wrote {pdf}')


if __name__ == '__main__':
    main()
