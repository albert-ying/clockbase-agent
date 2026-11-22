"""Replacement Figure 3, Tables 1 and 2, Supplementary Table S17 and the historical Tables 1 and 2 with their verified fate.

Empirical input (no simulation): FINAL/comparison-ranking from audit_system/comparison_ranking_2026_09_30.py,
run with the dual-model design audits, and FINAL/public-export (historical tables with their verified fate). Figure 3
shows, per category, the 15 comparisons with the largest bound of the exact confidence interval of Hedges' g nearest
zero; Tables 1 and 2 list the top 30 comparisons with higher and lower clock outputs; Supplementary Table S17 lists
every ranked comparison; the historical Tables 1 and 2 (Supplementary Table [[HIST-T]]) carry each row's verified fate.

  MPLBACKEND=Agg python scripts/figures_2026_10_03/build_figure3_tables_2026_10_03.py --final NAME --out DIR
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
                                 RUST, BLUE, INK, ERRORBAR)

CATEGORIES = ['Genetic', 'Drug', 'Environment', 'Disease']
MOD_SHORT = {'mouse RNA-seq': 'mm RNA', 'human RNA-seq': 'hs RNA', 'human DNAm': 'hs DNAm'}
CLOCK_SHORT = {'Hazard_Mouse_All_All_scaleddiff_EN_tAge': 'Mouse mortality tAge',
               'Hazard2_All_All_All_scaleddiff_EN_tAge': 'Multi-species mortality tAge'}
COLUMNS = {'rank': 'Rank', 'category': 'Category', 'modality': 'Modality', 'gse_id': 'GSE', 'label': 'Comparison',
           'comparison_name': 'Comparison (atlas label)', 'clock_name': 'Clock', 'v3_n_case': 'n case',
           'v3_n_control': 'n control', 'v3_mean_diff': 'Mean difference (clock units)', 'hedges_g': "Hedges' g",
           'g_ci_low': 'CI low', 'g_ci_high': 'CI high', 'ci_level': 'CI level', 'selection_k': 'Rows compared (k)',
           'v3_fdr': 'FDR', 'relevance_score': 'Relevance (0-10)', 'v3_final_class': 'Verification class',
           'v3_soft_flags': 'Soft flags'}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def tidy(frame):
    out = frame.copy().reset_index(drop=True)
    out['rank'] = np.arange(1, len(out) + 1)
    out['clock_name'] = out.clock.map(CLOCK_SHORT).fillna(out.clock)
    out['category'] = out.get('category', out.condition_category)
    return out[[c for c in COLUMNS if c in out.columns]].rename(columns=COLUMNS)


def tables():
    t1 = pd.read_csv(RANK / 'table1_pro_aging.csv')
    t2 = pd.read_csv(RANK / 'table2_anti_aging.csv')
    s17 = pd.read_csv(RANK / 'ranked_comparisons_all.csv')
    out = {'Table 1 (higher outputs)': tidy(t1), 'Table 2 (lower outputs)': tidy(t2), 'Table S17 (all ranked)': tidy(s17)}
    for name, path in (('Historical Table 1 (fate)', 'Table1_pro_aging_v3_status.csv'),
                       ('Historical Table 2 (fate)', 'Table2_anti_aging_v3_status.csv')):
        out[name] = pd.read_csv(FINAL / 'public-export' / path)
    with pd.ExcelWriter(OUT / 'ranking_tables_2026_10_03.xlsx') as writer:
        for name, frame in out.items():
            frame.to_excel(writer, sheet_name=name[:31], index=False)
            frame.to_csv(OUT / (name.split(' (')[0].replace(' ', '_').lower() + '_2026_10_03.csv'), index=False)
    return out


def figure(summary):
    data = pd.read_csv(RANK / 'figure3_category_top15.csv')
    retro_style('soft')
    fig = plt.figure(figsize=(7.2, 9.6))
    axd = {'a': fig.add_axes([0.40, 0.875, 0.50, 0.095]), 'b': fig.add_axes([0.46, 0.055, 0.50, 0.745])}

    # a: eligibility sequence of the ranking.
    funnel = summary['funnel']
    steps = [('FDR < 0.05, n ≥ 3, ranked category', funnel['usable_significant_n3_rankable_category']['comparisons']),
             ('Contrast type agreed (two models)', funnel['agreed_matching_contrast_type']['comparisons']),
             ('No agreed design concern', funnel['no_agreed_adjudication_design_concern']['comparisons']),
             ('Design audit admitted (two models)', funnel['context_admitted']['comparisons']),
             ('Ranked (CI excludes 0; best row per study and label)', funnel['after_dedup_by_study_and_label']['rows'])]
    ax = axd['a']
    y = np.arange(len(steps))[::-1]
    ax.barh(y, [n for _, n in steps], height=0.6, color=BLUE['200'], edgecolor=BLUE['400'], linewidth=0.5)
    for yy, (_, n) in zip(y, steps):
        ax.text(n + 80, yy, f"{n:,}", va='center', ha='left', fontsize=6, color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([l for l, _ in steps])
    ax.set_xlabel('Comparisons')

    # b: per category, the 15 comparisons with the largest confidence bound nearest zero.
    ax = axd['b']
    rows, positions, headers, yy = [], [], [], 0.0
    for category in CATEGORIES:
        part = data[data.category == category].reset_index(drop=True)
        headers.append((yy, f"{category} (top {len(part)} of {int(summary['outputs']['all']['by_category'][category])})"))
        yy += 1.0
        for _, r in part.iterrows():
            rows.append(r)
            positions.append(yy)
            yy += 1.0
        yy += 0.6
    frame = pd.DataFrame(rows).reset_index(drop=True)
    ypos = -np.array(positions)
    colors = np.where(frame.hedges_g < 0, BLUE['400'], RUST['400'])
    ax.errorbar(frame.hedges_g, ypos, xerr=[frame.hedges_g - frame.g_ci_low, frame.g_ci_high - frame.hedges_g],
                fmt='none', ecolor=INK, **ERRORBAR['forest'])
    ax.scatter(frame.hedges_g, ypos, s=14, c=colors, edgecolors=INK, linewidths=0.35, zorder=3, clip_on=False)
    labels = [f"{l[:28]} ({g}, {MOD_SHORT.get(m, m)}, {int(a)}/{int(b)})" for l, g, m, a, b in
              zip(frame.label, frame.gse_id, frame.modality, frame.v3_n_case, frame.v3_n_control)]
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=7)
    for y0, text in headers:
        ax.text(-0.02, -y0, text, transform=ax.get_yaxis_transform(), ha='right', va='center', fontsize=7,
                fontweight='bold', color=INK)
    ax.set_xlabel("Hedges' g with adjusted 95% CI (blue, lower clock output; rust, higher)")
    ax.set_ylim(-(yy - 0.1), 0.9)
    trim_axes(fig)
    set_domain(ax, ylim=(-(yy - 0.1), 0.9))  # room around the first and last rows
    ref_line(ax, x=0)
    add_panel_labels(fig, axd)
    save_figure(fig, str(OUT / 'figure3_ranking_2026_10_03.pdf'), axes_dict=axd,
                footer=f'data: ClockBase atlas, verifier v3 ({FINAL.name}); ranking computed in comparison_ranking_2026_09_30.py',
                review_png=str(OUT / 'figure3_ranking_2026_10_03_review.png'))
    plt.close(fig)


if __name__ == '__main__':
    FINAL, OUT = release_arguments(__doc__)
    RANK = FINAL / 'comparison-ranking'
    summary = json.loads((RANK / 'summary.json').read_text())
    tables()
    figure(summary)
    (OUT / 'provenance.json').write_text(json.dumps(
        {'ranking_dir': paths.portable(RANK), 'summary_sha256': sha256(RANK / 'summary.json'),
         'figure3_source_sha256': sha256(RANK / 'figure3_category_top15.csv'), 'script_sha256': sha256(__file__)}, indent=2))
    print('saved')
