"""Extended Data figure: verification of every atlas comparison.

Empirical inputs (no simulation), all from one verifier-v3 chain directory (FINAL) and the DNAm
input-integrity screen:
  repairs-{modality}/manifest.json      final class counts
  error-census/census.json              historical defects per modality
  published-status/summary.json         fate of the published S14, Figure 3 and Table 1/2 items
  input-integrity-dnam-001/samples.csv  per-sample correlation with the GSE40279 reference methylome
  replay-dnam-001/inputs/GSE73950.csv   historical Horvathv1 predictions of the index study

  MPLBACKEND=Agg python scripts/figures_2026_10_03/build_ed21_2026_10_03.py --final NAME --out DIR
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
from audit_system import paths  # noqa: E402
from nature_figure_utils import (release_arguments, retro_style, trim_axes, add_panel_labels, save_figure, set_domain, ref_line,  # noqa: E402
                                 label_panel_group, RUST, BLUE, GRAY, SAGE, GOLD, PLUM, INK)

V = paths.VERIFIER_ROOT
SCREEN = V / 'input-integrity-dnam-001' / 'samples.csv'
INDEX_STUDY = V / 'replay-dnam-001' / 'inputs' / 'GSE73950.csv'
MODALITIES = [('dnam', 'Human DNAm'), ('hs_rna', 'Human RNA-seq'), ('mm_rna', 'Mouse RNA-seq')]
MOD_COLOR = {'dnam': PLUM['400'], 'hs_rna': BLUE['400'], 'mm_rna': GOLD['400']}
CLASSES = [('VERIFIED', 'Verified', SAGE['400']), ('PROVISIONAL', 'Provisional', SAGE['200']),
           ('DESCRIPTIVE', 'Descriptive', GRAY['200']), ('DUPLICATE', 'Duplicate test', GRAY['100']),
           ('PENDING', 'Withheld, pending', GOLD['200']), ('NOT_REPAIRABLE', 'Not repairable', RUST['300'])]
DEFECTS = [('STORED_P_UNUSABLE', 'Stored P missing, 0 or 1'), ('ARMS_INVERTED_CONFIRMED', 'Arms swapped (two models)'),
           ('EFFECT_NOT_REPRODUCED', 'Effect not reproduced'), ('EMPTY_OR_ARTIFACT_INPUT', 'Empty or artefactual input'),
           ('TEST_NOT_A_GROUP_CONTRAST', 'Test not a group contrast'), ('DETACHED_INPUT_PROFILES', 'Detached DNAm profiles'),
           ('NON_STANDARD_ASSAY_IN_ARMS', 'Non-bisulfite assay'), ('FDR_BELOW_RAW_P', 'FDR below raw P')]
OUTPUTS = [('S14', 'Table S14'), ('Figure3_bars', 'Figure 3'), ('Table1_pro_aging', 'Table 1'),
           ('Table2_anti_aging', 'Table 2')]
FATES = [('CONFIRMED_CLAIMED_DIRECTION', 'Confirmed', SAGE['400']), ('NOT_SIGNIFICANT_AFTER_REPAIR', 'Not significant', GRAY['200']),
         ('DESCRIPTIVE_ONLY_ARM_BELOW_3', 'Descriptive', GRAY['100']), ('PENDING', 'Pending', GOLD['200']),
         ('CONTRADICTS_CLAIMED_DIRECTION', 'Reversed', RUST['400']), ('NO_EXACT_ROW', 'No atlas row', 'white')]


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def class_counts():
    counts = {}
    for key, _ in MODALITIES:
        raw = json.loads((FINAL / f'repairs-{key}' / 'manifest.json').read_text())['final_classes']
        merged = {c: 0 for c, _, _ in CLASSES}
        for name, n in raw.items():
            merged['PENDING' if name.startswith('PENDING') else name] += n
        counts[key] = merged
    return counts


def fate_counts():
    status = json.loads((FINAL / 'published-status' / 'summary.json').read_text())
    out = {}
    for key, _ in OUTPUTS:
        merged = {f: 0 for f, _, _ in FATES}
        for name, n in status[key]['outcomes'].items():
            merged['PENDING' if name.startswith('UNRESOLVED') else name] += n
        merged['NO_EXACT_ROW'] = status[key]['match'].get('NO_EXACT_ROW', 0)
        out[key] = merged
    return out


def stacked(ax, rows, labels, keys, spec, total_label=True):
    y = np.arange(len(rows))[::-1]
    for i, (row_key, row) in enumerate(zip(keys, rows)):
        total = sum(row.values())
        left = 0.0
        for name, _, color in spec:
            share = row[name] / total
            if share > 0:
                ax.barh(y[i], share, left=left, height=0.62, color=color, edgecolor=INK, linewidth=0.35)
            left += share
        if total_label:
            ax.text(1.01, y[i], f"$n$ = {total:,}", va='center', ha='left', fontsize=6, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(['0', '25', '50', '75', '100'])


def main():
    retro_style('soft')
    fig = plt.figure(figsize=(7.2, 6.3))
    # Hand-placed axes: categorical panels must stay short-wide and long row labels must not shift neighbours.
    axd = {'a': fig.add_axes([0.135, 0.855, 0.77, 0.085]), 'b': fig.add_axes([0.275, 0.445, 0.355, 0.300]),
           'c': fig.add_axes([0.765, 0.585, 0.215, 0.150]), 'd': fig.add_axes([0.095, 0.085, 0.335, 0.235]),
           'e': fig.add_axes([0.560, 0.085, 0.420, 0.235])}

    # a: final class of every atlas row, per modality.
    counts = class_counts()
    ax = axd['a']
    stacked(ax, [counts[k] for k, _ in MODALITIES], [l for _, l in MODALITIES], [k for k, _ in MODALITIES], CLASSES)
    ax.set_xlabel('Atlas comparisons (%)')
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=INK, linewidth=0.35) for _, _, c in CLASSES]
    ax.legend(handles, [l for _, l, _ in CLASSES], ncol=6, loc='lower left', bbox_to_anchor=(0.0, 1.02), frameon=False,
              fontsize=7, handlelength=1.0, columnspacing=1.0, borderaxespad=0.0)

    # b: material defects of the historical record, per modality (percentage of that modality's rows).
    census = json.loads((FINAL / 'error-census' / 'census.json').read_text())['modalities']
    ax = axd['b']
    y = np.arange(len(DEFECTS))[::-1]
    for offset, (key, label) in zip((0.2, 0.0, -0.2), MODALITIES):
        share = [100 * census[key]['defects'][d]['fraction'] for d, _ in DEFECTS]
        ax.scatter(share, y + offset, s=14, c=MOD_COLOR[key], edgecolors=INK, linewidths=0.35, zorder=2, clip_on=False,
                   label=label)
    ax.set_yticks(y)
    ax.set_yticklabels([l for _, l in DEFECTS])
    ax.set_xlabel('Rows with defect (%)')
    ax.legend(title='Modality', loc='lower right', frameon=False, fontsize=7, title_fontsize=7, handletextpad=0.2)

    # c: fate of the published items after repair, judged against their claimed direction.
    fates = fate_counts()
    ax = axd['c']
    stacked(ax, [fates[k] for k, _ in OUTPUTS], [l for _, l in OUTPUTS], [k for k, _ in OUTPUTS], FATES)
    ax.set_xlabel('Published items (%)')
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=c, edgecolor=INK, linewidth=0.35) for _, _, c in FATES]
    ax.legend(handles, [l for _, l, _ in FATES], ncol=2, loc='upper left', bbox_to_anchor=(-0.55, -0.42), frameon=False,
              fontsize=7, handlelength=0.9, columnspacing=0.8, borderaxespad=0.0, labelspacing=0.3)

    # d: the index study of detached profiles, GSE73950: input correlation against historical Horvath age.
    screen = pd.read_csv(SCREEN)
    study = screen[screen.gse_id == 'GSE73950'].set_index('gsm')
    meta = pd.read_csv(paths.GEO / 'GEO_hs_dnam/metadata/GSE73950.csv', low_memory=False).set_index('geo_accession')
    arm = meta.reindex(study.index).title.str[0].map({'E': 'Endometriosis', 'H': 'Control'})
    historical = pd.read_csv(INDEX_STUDY).set_index('gsm_id')['Horvathv1'].reindex(study.index)
    ax = axd['d']
    for name, color in (('Endometriosis', BLUE['300']), ('Control', RUST['300'])):
        mask = (arm == name).to_numpy()
        ax.scatter(study.r_reference.to_numpy()[mask], historical.to_numpy()[mask], s=12, c=color, edgecolors=INK,
                   linewidths=0.3, clip_on=False, zorder=2, label=f"{name} ($n$ = {int(mask.sum())})")
    ax.set_xlabel('Correlation with reference methylome ($r$)')
    ax.set_ylabel('Historical Horvath age (years)')
    ax.set_title('GSE73950', fontsize=7, fontweight='normal', loc='left')
    ax.legend(loc='upper left', bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=7, handletextpad=0.2)
    controls = arm == 'Control'
    ax.annotate(f"{int(controls.sum())} controls overlap\n($r$ = {study.r_reference[controls].median():.3f}; "
                f"{historical[controls].median():.1f} years)", xy=(study.r_reference[controls].median(), historical[controls].median()),
                xytext=(0.22, 4), fontsize=6, color=RUST['400'], ha='left', va='bottom',
                arrowprops=dict(arrowstyle='-', color=RUST['400'], lw=0.5))

    # e: all screened samples, with the detached (hard) and atypical (soft) thresholds.
    ax = axd['e']
    r = screen.r_reference.dropna().to_numpy()
    bins = np.linspace(-0.45, 1.0, 59)
    ax.hist(r, bins=bins, color=GRAY['200'], edgecolor=INK, linewidth=0.3, log=True)
    ax.hist(r[r < 0.10], bins=bins, color=RUST['300'], edgecolor=INK, linewidth=0.3, log=True)
    ax.set_xlabel('Correlation with reference methylome ($r$)')
    ax.set_ylabel('Samples')
    ax.set_title(f"$n$ = {len(r):,} samples in {screen.dropna(subset=['r_reference']).gse_id.nunique():,} studies",
                 fontsize=7, fontweight='normal', loc='left')

    trim_axes(fig, exclude=[axd['a'], axd['c'], axd['e']])
    for key in ('d', 'e'):  # same quantity on both x axes
        set_domain(axd[key], xlim=(-0.5, 1.05))
    yb = axd['b'].get_ylim()
    set_domain(axd['b'], ylim=(min(yb[0], -0.55), max(yb[1], len(DEFECTS) - 0.45)))
    ref_line(axd['e'], x=0.10)
    ref_line(axd['e'], x=0.50)
    axd['e'].text(0.085, 0.96, f"{int((r < 0.10).sum()):,} below 0.10", transform=axd['e'].get_xaxis_transform(),
                  ha='right', va='top', fontsize=6, color=RUST['400'])
    add_panel_labels(fig, axd)
    provenance = {'final_dir': paths.portable(FINAL), 'screen_sha256': sha256(SCREEN), 'index_study_sha256': sha256(INDEX_STUDY),
                  'census_sha256': sha256(FINAL / 'error-census' / 'census.json'),
                  'published_status_sha256': sha256(FINAL / 'published-status' / 'summary.json'),
                  'script_sha256': sha256(__file__)}
    (OUT / 'provenance.json').write_text(json.dumps(provenance, indent=2))
    save_figure(fig, str(OUT / 'ed_verification_census_2026_10_03.pdf'), axes_dict=axd,
                footer=f'data: ClockBase atlas, verifier v3 ({FINAL.name}); statistics computed in build_figure_2026_10_03.py',
                review_png=str(OUT / 'ed_verification_census_2026_10_03_review.png'))
    plt.close(fig)
    print('saved')


if __name__ == '__main__':
    FINAL, OUT = release_arguments(__doc__)
    main()
