"""Extended Data Figure 23: ClockEffectBench and the agent reruns.

Reads the benchmark results of a final run that the code release holds (benchmark/published_clock_2026_10_05:
groupbench/results_<RUN>/ground_truth_funnel.json and published_released_test.json; agentbench/results_2026_10_11/
agent_reruns.json and results_2026_10_14/ranked_reruns.json for Codex, results_2026_10_18 for Claude Code and
results_2026_10_19 for Biomni) and writes ed23_benchmark_2026_10_14.pdf, a review PNG and provenance.json. Every value
is read from those files.

  python -B scripts/figures_2026_10_03/build_ed23_benchmark_2026_10_14.py --out "$TMPDIR/ed23"
"""
import os
import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'figures_2026_10_03'))
from nature_figure_utils import (release_arguments, retro_style, trim_axes, add_panel_labels, save_figure,  # noqa: E402
                                 set_domain, ref_line, pad_categorical_axis, figsize_nature, ERRORBAR,
                                 BLUE, RUST, GRAY, SAGE, INK, TYPE)

from ed21_24_restore_common_2026_11_05 import (style, axes, grid, save, readme, table_md,
    WIDTH, FINAL, CLASSIC, HIGHER, LOWER)
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
BLUE = {'400': LOWER, '200': '#9fc2df', '300': CLASSIC[0]}
RUST = {'400': HIGHER, '200': '#f2c49b', '300': CLASSIC[1]}
SAGE = {'400': CLASSIC[2]}
GRAY = {'200': '#cccccc', '300': '#aaaaaa', '400': '#777777'}
TYPE['annotation'] = 6.5
OUTCOMES = [('significant, same direction', 'Sig., same', BLUE['400']),
            ('not significant, same direction', 'NS, same', BLUE['200']),
            ('not significant, opposite direction', 'NS, opposite', RUST['200']),
            ('significant, opposite direction', 'Sig., opposite', RUST['400']),
            ('not testable', 'Not testable', GRAY['200'])]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    final = FINAL
    stamp = final.name.replace('final-', '').replace('-', '_')
    results = Path(os.environ.get('CLOCKBASE_BENCH_RESULTS') or BENCH / f'groupbench/results_{stamp}')
    files = {'funnel': results / 'ground_truth_funnel.json',
             'released': results / 'published_released_test.json',
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
    style()
    fig = plt.figure(figsize=(WIDTH, 9.7))
    axd = {
        'a': axes(fig, [.235,.85,.255,.11], letter='a', label_x=.015),
        'b': axes(fig, [.715,.90,.255,.05], letter='b', label_x=.565),
        'c': axes(fig, [.215,.61,.275,.14], letter='c', label_x=.015),
        'd': axes(fig, [.735,.61,.235,.14], letter='d', label_x=.565),
        'e': axes(fig, [.36,.385,.57,.16], letter='e', label_x=.015),
        'f': axes(fig, [.36,.075,.61,.245], letter='f', label_x=.015)}
    printed = {}

    # a: the ground truth and why each published effect is or is not scored with the paper's groups
    ax = axd['a']
    reasons = dict(funnel['reasons'])
    reasons['RNA'] = sum(v for k, v in reasons.items() if not k.startswith('DNAm'))
    values = [reasons.get(key, 0) for key, _ in FUNNEL]
    y = np.arange(len(FUNNEL))[::-1]
    colors = [CLASSIC[2], CLASSIC[0], CLASSIC[3], CLASSIC[1], '#C44E52',
              CLASSIC[4], CLASSIC[5], '#999999', '#cccccc']
    ax.barh(y, values, height=0.62, color=colors, edgecolor='none')
    for yi, value in zip(y, values):
        ax.text(value + 3, yi, f'{value:,}', va='center', ha='left', fontsize=TYPE['annotation'], color=INK)
    ax.set_yticks(y)
    ax.set_yticklabels([name for _, name in FUNNEL])
    ax.set_xlabel('Published effects')
    ax.set_title(f"Ground truth ({funnel['ground_truth']:,} effects)", loc='left')
    ax.set_xticks([0, 100, 200])
    printed['a'] = dict(zip([n for _, n in FUNNEL], values))

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
    ax.set_title("Paper’s groups\nClockBase clock outputs", loc='left', pad=13)
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), loc='upper left', bbox_to_anchor=(-.04, -1.22), ncol=2, frameon=False,
              fontsize=7, handlelength=.8, columnspacing=.7, borderaxespad=0.1, labelspacing=.15)

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
    ax.set_title('Without the data', loc='left', pad=18)
    ax.set_xticks([0, 50, 100])
    ax.plot([], [], 'o', ms=4.5, color=BLUE['400'], label=f"Lower ({total['lower']})")
    ax.plot([], [], 'o', ms=4.5, color=RUST['400'], label=f"Higher ({total['higher']})")
    ax.legend(loc='upper center', bbox_to_anchor=(.5, 1.04), ncol=2, frameon=False, fontsize=7,
              handletextpad=0.2, columnspacing=.8, borderaxespad=0.1)
    ax = axd['d']
    for yi, (key, _) in zip(y, order):
        value = systems[key]['call_balanced_accuracy']
        low, high = systems[key]['call_balanced_accuracy_interval']
        color = BLUE['400'] if key == 'S3c' else GRAY['400']
        ax.errorbar([value], [yi], xerr=[[value - low], [high - value]], fmt='o', ms=4, color=color, ecolor=color,
                    **ERRORBAR['forest'])
    ax.set_yticks(y)
    ax.set_yticklabels([name for _, name in order])
    ax.set_xlabel('Balanced accuracy')
    ax.set_title('Three-way call', loc='left', pad=18)
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
    ax.text(1.02, len(rows) - 0.10, '$n$', transform=ax.get_yaxis_transform(), ha='left', va='bottom',
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
    height = 0.23
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
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(.65, .005), ncol=4, frameon=False, fontsize=7,
               handlelength=.8, columnspacing=.7, handletextpad=.3, borderaxespad=0)

    # Restore July outlined bars, substantial markers, dotted value-axis grids.
    for key, ax in axd.items():
        ax.set_ylim(-.65, len(ax.get_yticks())-.35)
        ax.tick_params(axis='y', length=0)
        grid(ax)
        for patch in ax.patches:
            patch.set_edgecolor('#333333'); patch.set_linewidth(.4)
    axd['a'].set_xlim(0,205); axd['a'].set_xticks([0,100,200])
    axd['b'].set_xlim(0,55); axd['b'].set_xticks([0,25,50])
    axd['c'].set_xlim(-2,102); axd['c'].set_xticks([0,50,100])
    axd['c'].set_ylim(-.65,7.3)
    axd['d'].set_ylim(-.65,7.3)
    axd['d'].set_xlim(.3,.9); axd['d'].set_xticks([.3,.6,.9])
    axd['d'].axvline(1/3, color='#777777', lw=.7, ls='--')
    axd['e'].set_xlim(0,100); axd['e'].set_xticks([0,25,50,75,100])
    axd['f'].set_xlim(0,32); axd['f'].set_xticks([0,10,20,30])
    # Values on the grouped bars follow July ED7; each of the four bars keeps its value.
    for patch in axd['f'].patches:
        width = patch.get_width()
        if width:
            axd['f'].text(width+.2, patch.get_y()+patch.get_height()/2, f'{width:g}',
                          va='center', fontsize=6)
    printed['ground_truth_total'] = funnel['ground_truth']
    printed['b'] = released['published_significant']
    printed['c_d'] = {key: systems[key] for key,_ in order}
    printed['c_direction_denominators'] = total
    printed['e_series'] = reruns['series']
    printed['e'] = [{'key':key,'label':label,**own[key]} for key,label,_ in rows]
    printed['f'] = [{'framework':framework,'model':model,**block}
                    for framework,model,block in configurations]
    out = save(fig, 'ed23', axd, 13, sources=list(files.values()), facts=printed, script=__file__)
    import pandas as pd
    pd.DataFrame([{'reason':label,'n':n} for label,n in printed['a'].items()]).to_csv(out/'panel_a_funnel.csv',index=False)
    edata = pd.DataFrame([{'system':label,'agreement':own[key]['call_agreement'],
                          'n':own[key]['call_n'],'interval':own[key].get('call_agreement_interval')}
                         for key,label,_ in rows])
    edata.to_csv(out/'panel_e_agreement.csv',index=False)
    fdata = pd.DataFrame([{'framework':framework,'model':model,
                         **{f'{direction}_{metric}':block[direction][metric]
                            for direction in ['higher','lower'] for metric in ['found','recurs']}}
                        for framework,model,block in configurations])
    fdata.to_csv(out/'panel_f_ranked_counts.csv',index=False)
    cdata = pd.DataFrame([{'system':label,'lower_n':systems[key]['published_lower_identified'][0],
                         'lower_total':systems[key]['published_lower_identified'][1],
                         'higher_n':systems[key]['published_higher_identified'][0],
                         'higher_total':systems[key]['published_higher_identified'][1],
                         'balanced_accuracy':systems[key]['call_balanced_accuracy'],
                         'interval':systems[key]['call_balanced_accuracy_interval']}
                        for key,label in order])
    cdata.to_csv(out/'panel_c_d_models.csv',index=False)
    readme(out, f"""# Benchmark performance

ED23 retains all six October panels. July ED7 supplies outlined category bars
and bar-value labels; July Figure 3 supplies the orange/blue directional key.
July ED19 supplies compact horizontal rows and bold panel titles. The July
comparison uses page 13 at the same physical scale.

## Panel map

| October panel | Restored panel | Content |
| --- | --- | --- |
| a | a | All nine ground-truth coverage categories |
| b | b | The same released-test and Welch-test outcome counts |
| c | c | Direction-specific identification for the same seven systems |
| d | d | Three-way balanced accuracy and source intervals for those systems |
| e | e | All benchmark-series agent and no-data model results with their intervals and n |
| f | f | All six rerun configurations and all four found/recurs count series |

## Data and numbers

The eight input JSON files are listed by absolute path, modification time and
SHA-256 in `provenance.json`. They are the same benchmark outputs read by the
October builder; no new model run is performed. `numbers.json` retains every
plotted value and printed numerical annotation. Panel CSV files provide the
coverage, identification, accuracy, agreement, sample-size and ranked counts.
Axis ticks are in the stated units; the 1/3 reference is three-class chance.

Panel a uses `groupbench/results_2026_10_10/ground_truth_funnel.json`.
Ground-truth total = {funnel['ground_truth']}; nine bar counts are:

{table_md(pd.DataFrame([{'reason':k,'n':v} for k,v in printed['a'].items()]))}

Panel b uses `published_released_test.json`, `published_significant.S3` and
`.S3c`. Counts by outcome are retained exactly, including unlabelled narrow
segments; all counts are present in `numbers.json`. Same-direction categories
are blue, opposite-direction categories orange, and untestable effects grey.

Panels c/d use `against_models_without_data.systems` in that same file.
Panel c divides each identified count by its published-direction denominator
and multiplies by 100. The published significant effects comprise
{total['lower']} decreases and {total['higher']} increases. Panel d reads the
balanced-accuracy estimate and interval without changing their calculation.
The 95% intervals in d and e are percentile bootstrap intervals over papers,
using 1,000 replicates and seed 20261005 in
`benchmark/published_clock_2026_10_05/summarize_2026_10_05.py`. The published
test comparison and agent rerun scorer call its `with_interval` function.

{table_md(cdata)}

Panel e reads the `own_coverage` blocks from agent reruns dated October 11
(Codex), October 18 (Claude Code) and October 19 (Biomni). There are
{reruns['series']} benchmark series. Agreement and supplied intervals are
multiplied by 100; the right-hand n is the system's call_n, not the number
of series. Coverage differs by system, so this panel does not imply a shared
set of scored calls. Missing intervals retain the source point-only display.

{table_md(edata)}

Panel f reads ranked reruns dated October 14 (Codex), October 18 (Claude Code)
and October 19 (Biomni), fields `by_model[model][higher/lower].found/recurs`.
Each direction has 30 ranked comparisons.

{table_md(fdata)}

## Draft legend

**a,** Coverage of {funnel['ground_truth']} published clock effects in the
benchmark. Bar labels give effects in each coverage or exclusion category.
**b,** Outcomes for published significant effects when the paper's groups are
tested using ClockBase clock outputs with a Welch test or the released test.
Dark and pale blue indicate significant and non-significant same-direction
results; pale and dark orange indicate non-significant and significant
opposite-direction results; grey indicates untestable effects. **c,**
Percentage of published decreases (n = {total['lower']}, blue) and increases
(n = {total['higher']}, orange) identified by ClockBase, five language models
given the study description without data, and an always-higher baseline.
Connected dots pair the two directions within a system. **d,** Three-way
call balanced accuracy for the same systems, with 95% percentile bootstrap
intervals over papers (1,000 replicates). The dashed line marks 1/3 chance
accuracy. **e,** Agreement with
the published call for historical and released rows, Codex, Claude Code and
Biomni reruns, and no-data language-model baselines across
{reruns['series']} benchmark series. Points and whiskers show agreement and
95% percentile bootstrap intervals over papers (1,000 replicates); n at right is the number of calls evaluated for that
system. **f,** Numbers of higher- and lower-output ranked comparisons found
and recurring in each rerun configuration, out of 30 per direction.
Pale bars show found comparisons and dark bars recurring effects; orange
and blue denote higher and lower output. Performance summaries and intervals
are reproduced from the benchmark; no new statistical tests are added.

## Checks

`checks.json` records text-floor, overlap, clipping, exact PDF size and vector
checks. The page is 183 × 246.38 mm. Review PNGs use 200 dpi. The same source
system order and every outcome category, configuration and uncertainty
interval are retained. The 200-dpi render and same-scale July comparison are
inspected before delivery. The original builders and combined package remain
unchanged.
""")


if __name__ == '__main__':
    main()
