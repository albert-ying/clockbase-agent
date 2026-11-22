"""Extended Data Figure 23: ClockEffectBench and the agent reruns.

Reads the benchmark results of a final run that the code release holds (benchmark/published_clock_2026_10_05:
groupbench/results_<RUN>/ground_truth_funnel.json and published_released_test.json; agentbench/results_2026_10_11/
agent_reruns.json and results_2026_10_14/ranked_reruns.json for Codex, results_2026_10_18 for Claude Code and
results_2026_10_19 for Biomni) and writes ed23_benchmark_2026_10_14.pdf, a review PNG and provenance.json. Every value
is read from those files.

  python -B scripts/figures_2026_10_03/build_ed23_benchmark_2026_10_14.py --out "$TMPDIR/ed23"
"""
import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from nature_figure_utils import (release_arguments, retro_style, trim_axes, add_panel_labels, save_figure,  # noqa: E402
                                 set_domain, ref_line, pad_categorical_axis, figsize_nature, ERRORBAR,
                                 BLUE, RUST, GRAY, SAGE, INK, TYPE)

REPO = Path(__file__).resolve().parents[2]
BENCH = REPO / 'benchmark/published_clock_2026_10_05'
MODEL = {'gpt-6-astra': 'GPT-6 Astra', 'gpt-6-sol': 'GPT-6 Sol', 'gpt-6-luna': 'GPT-6 Luna', 'gpt-5.5': 'GPT-5.5',
         'gpt-4.1': 'GPT-4.1', 'claude-opus-5-5': 'Claude Opus 5.5'}
FUNNEL = [("DNAm: scored with the paper's groups", 'Scored'), ('DNAm: NO_CLOCK_MAPPING', 'Clock not in ClockBase'),
          ('DNAm: UNRESOLVED_GROUPS', 'Groups not resolved'), ('DNAm: no published direction or significance',
                                                               'No direction or P'),
          ('RNA', 'RNA-seq clock'), ('DNAm: INSUFFICIENT_RECORDED_AGE', 'Ages not recorded'),
          ('DNAm: UNSUPPORTED_MEASURE', 'Measure not formable'), ('DNAm: EMPTY_ARM', 'Empty group'),
          ('DNAm: covered, no statistic', 'No statistic')]
OUTCOMES = [('significant, same direction', 'Significant, same', BLUE['400']),
            ('not significant, same direction', 'Same, not significant', BLUE['200']),
            ('not significant, opposite direction', 'Opposite, not significant', RUST['200']),
            ('significant, opposite direction', 'Significant, opposite', RUST['400']),
            ('not testable', 'Not testable', GRAY['200'])]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    final, out = release_arguments(__doc__)
    stamp = final.name.replace('final-', '').replace('-', '_')
    files = {'funnel': BENCH / f'groupbench/results_{stamp}/ground_truth_funnel.json',
             'released': BENCH / f'groupbench/results_{stamp}/published_released_test.json',
             'reruns': BENCH / 'agentbench/results_2026_10_11/agent_reruns.json',
             'ranked': BENCH / 'agentbench/results_2026_10_14/ranked_reruns.json',
             'claude_reruns': BENCH / 'agentbench/results_2026_10_18/agent_reruns.json',
             'claude_ranked': BENCH / 'agentbench/results_2026_10_18/ranked_reruns.json',
             'biomni_reruns': BENCH / 'agentbench/results_2026_10_19/agent_reruns.json',
             'biomni_ranked': BENCH / 'agentbench/results_2026_10_19/ranked_reruns.json'}
    for path in files.values():
        if not path.is_file():
            raise SystemExit(f'Cannot proceed - required data file not found: {path}')
    funnel = json.loads(files['funnel'].read_text())
    released = json.loads(files['released'].read_text())
    reruns = json.loads(files['reruns'].read_text())
    ranked = json.loads(files['ranked'].read_text())
    claude_reruns = json.loads(files['claude_reruns'].read_text())
    claude_ranked = json.loads(files['claude_ranked'].read_text())
    biomni_reruns = json.loads(files['biomni_reruns'].read_text())
    biomni_ranked = json.loads(files['biomni_ranked'].read_text())
    retro_style()
    fig, axd = plt.subplot_mosaic([['a', 'b'], ['c', 'd'], ['e', 'f']],
                                  figsize=figsize_nature('double', aspect=0.98), layout='constrained',
                                  gridspec_kw=dict(width_ratios=[1.0, 1.0], height_ratios=[0.8, 0.85, 1.0]))
    fig.get_layout_engine().set(h_pad=0.14, w_pad=0.05)

    # a: the ground truth and why each published effect is or is not scored with the paper's groups
    ax = axd['a']
    reasons = dict(funnel['reasons'])
    reasons['RNA'] = sum(v for k, v in reasons.items() if not k.startswith('DNAm'))
    values = [reasons.get(key, 0) for key, _ in FUNNEL]
    y = np.arange(len(FUNNEL))[::-1]
    colors = [BLUE['300'] if i == 0 else GRAY['300'] for i in range(len(FUNNEL))]
    ax.barh(y, values, height=0.62, color=colors, edgecolor='none')
    for yi, value in zip(y, values):
        ax.text(value + 3, yi, f'{value:,}', va='center', ha='left', fontsize=TYPE['annotation'], color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([name for _, name in FUNNEL])
    ax.set_xlabel('Published effects')
    ax.set_title(f"Ground truth ({funnel['ground_truth']:,} effects)", loc='left')
    ax.set_xticks([0, 100, 200])

    # b: published significant effects with the paper's groups: Welch test for every effect and the released test
    ax = axd['b']
    rows = [('Welch test', released['published_significant']['S3']), ('Released test', released['published_significant']['S3c'])]
    y = np.arange(len(rows))[::-1]
    for yi, (_, counts) in zip(y, rows):
        left = 0
        for key, name, color in OUTCOMES:
            value = counts.get(key, 0)
            if value:
                ax.barh(yi, value, left=left, height=0.55, color=color, edgecolor='white', lw=0.5, label=name)
                if value >= 3:
                    ax.text(left + value / 2, yi, str(value), ha='center', va='center', fontsize=TYPE['annotation'],
                            color='white' if color in (BLUE['400'], RUST['400']) else INK)
                left += value
    ax.set_yticks(y)
    ax.set_yticklabels([name for name, _ in rows])
    ax.set_xlabel('Published significant effects')
    ax.set_title("Paper's groups, ClockBase clock outputs", loc='left')
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc='upper center', bbox_to_anchor=(0.45, -0.38), ncol=3, frameon=False,
              fontsize=TYPE['legend'], handlelength=1.0, columnspacing=0.8, borderaxespad=0.1)

    # c, d: against language models given the study description without the data (same rows)
    systems = released['against_models_without_data']['systems']
    order = [('S3c', 'ClockBase')] + [(f'S4:openai:{m}', MODEL[m]) for m in
                                      ('gpt-6-astra', 'gpt-6-sol', 'gpt-6-luna', 'gpt-5.5', 'gpt-4.1')] + [('S0', "Always 'higher'")]
    total = released['against_models_without_data']['published_significant']
    y = np.arange(len(order))[::-1]
    ax = axd['c']
    for yi, (key, _) in zip(y, order):
        lower = systems[key]['published_lower_identified']
        higher = systems[key]['published_higher_identified']
        ax.plot([higher[0] / higher[1] * 100], [yi], 'o', ms=4.5, color=RUST['400'], clip_on=False, zorder=3)
        ax.plot([lower[0] / lower[1] * 100], [yi], 'o', ms=4.5, color=BLUE['400'], clip_on=False, zorder=3)
        ax.plot([lower[0] / lower[1] * 100, higher[0] / higher[1] * 100], [yi, yi], color=GRAY['300'], lw=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([name for _, name in order])
    ax.set_xlabel('Published effects identified (%)')
    ax.set_title('Without the data', loc='left')
    ax.set_xticks([0, 50, 100])
    ax.plot([], [], 'o', ms=4.5, color=BLUE['400'], label=f"Decreases ({total['lower']})")
    ax.plot([], [], 'o', ms=4.5, color=RUST['400'], label=f"Increases ({total['higher']})")
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.32), ncol=2, frameon=False, fontsize=TYPE['legend'],
              handletextpad=0.2, borderaxespad=0.1)
    ax = axd['d']
    for yi, (key, _) in zip(y, order):
        value = systems[key]['call_balanced_accuracy']
        low, high = systems[key]['call_balanced_accuracy_interval']
        color = BLUE['400'] if key == 'S3c' else GRAY['400']
        ax.errorbar([value], [yi], xerr=[[value - low], [high - value]], fmt='o', ms=4, color=color, ecolor=color,
                    **ERRORBAR['forest'])
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.set_xlabel('Balanced accuracy')
    ax.set_title('Three-way call', loc='left')
    ax.set_xticks([0.3, 0.6, 0.9])

    # e: the discovery task rerun on the benchmark series: agreement with the published call
    # Biomni runs GPT-6 Astra, as Codex does, so its rows are keyed by framework.
    own = {**reruns['own_coverage'], **{k: v for k, v in claude_reruns['own_coverage'].items() if k.startswith('R:')},
           **{'biomni:' + k: v for k, v in biomni_reruns['own_coverage'].items() if k.startswith('R:')}}
    rows = ([('S1', 'Historical agent', SAGE['400']), ('S2c', 'Released rows', BLUE['400'])]
            + [(f'R:{m}:raw', f'Codex + {MODEL[m]}', RUST['400']) for m in reruns['models']]
            + [(f'R:{m}:raw', f"Claude Code + {MODEL[m].replace('Claude ', '')}", RUST['300'])
               for m in claude_reruns['models']]
            + [(f'biomni:R:{m}:raw', f'Biomni + {MODEL[m]}', RUST['300']) for m in biomni_reruns['models']]
            + [(f'S4:openai:{m}', f'{MODEL[m]}, no data', GRAY['400']) for m in reruns['models']]
            + [('S0', "Always 'higher'", GRAY['300'])])
    ax = axd['e']
    y = np.arange(len(rows))[::-1]
    for yi, (key, _, color) in zip(y, rows):
        block = own[key]
        value = block['call_agreement']
        low, high = block.get('call_agreement_interval') or (value, value)
        ax.errorbar([value * 100], [yi], xerr=[[value * 100 - low * 100], [high * 100 - value * 100]], fmt='o', ms=4,
                    color=color, ecolor=color, **ERRORBAR['forest'])
        ax.text(1.02, yi, f"{block['call_n']}", transform=ax.get_yaxis_transform(), ha='left', va='center',
                fontsize=TYPE['annotation'], color=INK)
    ax.text(1.02, len(rows) - 0.75, '$n$', transform=ax.get_yaxis_transform(), ha='left', va='bottom',
            fontsize=TYPE['annotation'], color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([name for _, name, _ in rows])
    ax.set_xlabel('Agreement with published call (%)')
    ax.set_title(f"Benchmark series ({reruns['series']})", loc='left')
    ax.set_xticks([0, 50, 100])

    # f: the ranked comparisons of Tables 1 and 2 reanalysed by the reruns
    ax = axd['f']
    configurations = ([('Codex', m, ranked['by_model'][m]) for m in ranked['models']]
                      + [('Claude Code', m, claude_ranked['by_model'][m]) for m in claude_ranked['models']]
                      + [('Biomni', m, biomni_ranked['by_model'][m]) for m in biomni_ranked['models']])
    # Horizontal, with the configuration labels of panel e: six two-line tick labels do not fit under a half-width panel.
    y = np.arange(len(configurations))[::-1]
    height = 0.19
    for j, (table, color) in enumerate((('higher', RUST), ('lower', BLUE))):
        found = [b[table]['found'] for _, _, b in configurations]
        recurs = [b[table]['recurs'] for _, _, b in configurations]
        offset = (0.5 - j) * 2 * height
        ax.barh(y + offset + height / 2, found, height=height, color=color['200'], edgecolor=color['400'], lw=0.5,
                label=f'{table.capitalize()}: found')
        ax.barh(y + offset - height / 2, recurs, height=height, color=color['400'], edgecolor=color['400'], lw=0.5,
                label=f'{table.capitalize()}: recurs')
    ax.set_yticks(y)
    ax.set_yticklabels([f"{framework} + {MODEL[m].replace('Claude ', '') if framework == 'Claude Code' else MODEL[m]}"
                        for framework, m, _ in configurations])
    ax.set_xlabel('Ranked comparisons (of 30)')
    ax.set_title('Tables 1 and 2 reanalysed', loc='left')
    ax.legend(loc='upper center', bbox_to_anchor=(0.45, -0.3), ncol=2, frameon=False, fontsize=TYPE['legend'],
              handlelength=1.0, columnspacing=0.8, borderaxespad=0.1)

    trim_axes(fig, exclude=[axd['b']])
    for key in ('a', 'c', 'd', 'e'):
        pad_categorical_axis(axd[key], 'y', pad=0.5)
    pad_categorical_axis(axd['b'], 'y', pad=0.5)
    set_domain(axd['c'], xlim=(0, 100))
    set_domain(axd['a'], xlim=(0, 200))
    axd['a'].set_xticks([0, 100, 200])
    axd['a'].spines['bottom'].set_bounds(0, 200)
    set_domain(axd['e'], xlim=(0, 100))
    set_domain(axd['d'], xlim=(0.3, 0.9))
    pad_categorical_axis(axd['f'], 'y', pad=0.5)
    set_domain(axd['f'], xlim=(0, 30))
    axd['f'].set_xticks([0, 10, 20, 30])
    ref_line(axd['d'], x=1 / 3)
    add_panel_labels(fig, axd)
    pdf = out / 'ed23_benchmark_2026_10_14.pdf'
    try:
        save_figure(fig, pdf, axes_dict=axd, review_png=out / 'ed23_benchmark_2026_10_14_review.png',
                    footer=f'data: ClockEffectBench on {final.name}; agent reruns 2026-10-11, 2026-10-14, 2026-10-18 and 2026-10-19')
    except Exception:
        fig.savefig(out / 'ed23_FAILED_preview.png', dpi=200)
        raise
    (out / 'provenance.json').write_text(json.dumps({
        'final': final.name, 'inputs_sha256': {str(p.relative_to(REPO)): sha256(p) for p in files.values()},
        'builder_sha256': sha256(__file__), 'output_sha256': sha256(pdf)}, indent=2) + '\n')
    print(f'wrote {pdf}')


if __name__ == '__main__':
    main()
