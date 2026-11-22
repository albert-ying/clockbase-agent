# Biomni reruns (PROTOCOL.md items 24, 25 and 26, 2026-10-18)

Third agent framework: Biomni A1. Public code from
github.com/snap-stanford/Biomni at commit 400c1f366b96a35ca253e13c9b06c5076af41d65 (version 0.0.8, Apache-2.0), with
`gpt-6-astra` through the OpenAI API at medium reasoning effort. This is the model of the first Codex configuration,
so the two configurations differ only in the framework. Same tasks, inputs and scoring as the Codex reruns (items 16
and 20) and the Claude Code reruns (items 22 and 23).

| Set | Runs | Excluded | Result |
|---|---:|---:|---|
| Published results (61 series) | 61 | 0 | call agreement 0.81 (n 68; recomputed 0.78, n 67), direction agreement 0.92 (n 26), direction AUROC 0.977, 82 of 158 results covered |
| Ranked comparisons (56 series) | 56 | 0 | found 56 of 60, 55 in the ranked direction; 54 recur (27 higher, 27 lower) |

Across Codex (GPT-6 Astra, GPT-6 Sol, GPT-6 Luna, GPT-5.5), Claude Code (Claude Opus 5.5) and Biomni (GPT-6 Astra):
56 of the 60 ranked comparisons recur in at least two configurations and 43 in all six
(`ranked_reruns_six_configurations.json`).

One comparison is found with the opposite sign: ouabain (GSE122080), in the Biomni run only.
- Biomni saved a control arm that pooled the six old saline-treated mice with six young ones (Jaccard index with the
  verified arm 0.5, the matching threshold) and fitted an age-stratified model.
- On the pooled arms the verifier's Welch statistic is +0.15 (P = 0.68, not significant).
- Biomni's own age-stratified estimate has the ranked direction (P = 9 × 10⁻⁶).
- The five other configurations found it with the verified arms (Jaccard 1.0) and P = 9 × 10⁻⁶.

On the 70 published results covered by both the historical agent and the Biomni rerun, the calls agreed in 82%
(historical agent) and 83% (Biomni).

Runs:
- All 117 runs exited 0 and wrote a results table. No run was excluded: there was no path outside the run directory,
  no network command and no Biomni tool module.
- 8 benchmark runs reported no testable comparison: no rows, only "not estimable" rows, or rows without sample lists.
- Median 95 and 73.5 seconds per run.
- Tokens, benchmark set: 2.95 M input (0.76 M cached), 1.27 M output.
- Tokens, ranked set: 2.34 M input (0.52 M cached), 0.26 M output.

First attempt, unscored (item 25): all 117 runs exited 0, but 30 ended without a results table.
- The Responses API does not apply A1's stop sequences. In 37 runs, one model message held both an executable block
  and a `<solution>` tag, and A1 ended before running the block.
- The runner now cuts each reply at the first `</execute>` or `</solution>`, and all 117 runs were repeated.
- The first attempt is kept in `/Entropy/kying/clock-agent-runs-biomni-2026-10-18` and
  `/Entropy/kying/clock-agent-runs-ranked-biomni-2026-10-18`.
- Its tokens: 4.77 M input (0.98 M cached), 0.73 M output.

`ranked combine` now treats a configuration as scored directory plus model, because gpt-6-astra runs under two
frameworks. The five-configuration result of item 23 is unchanged (55 and 45).

Files:
- Published results: `agent_reruns.json`, `runs.csv`, `records.jsonl`.
- Ranked comparisons: `ranked_reruns.json`, `ranked_reruns.csv`, `ranked_runs.csv`.
- Combination: `ranked_reruns_six_configurations.json`.
- Example task: `TASK_example_GSE122080_biomni-gpt-6-astra.md`.
- Biomni environment: `biomni_pip_freeze.txt`.
- Run roots: `/Entropy/kying/clock-agent-runs-biomni-r2-2026-10-18` and
  `/Entropy/kying/clock-agent-runs-ranked-biomni-r2-2026-10-18`, linked as `DATA/agentbench/runs-biomni` and
  `ranked-runs-biomni`.

Commands (DATA is `/Entropy/kying/clockbase-benchmark-2026-10-05`; GEO is `/Entropy/kying/clock-agent-runs-ranked-2026-10-14/geo`):

    git clone https://github.com/snap-stanford/Biomni.git src && git -C src checkout 400c1f366b96a35ca253e13c9b06c5076af41d65
    python -m venv venv && venv/bin/pip install ./src langgraph langchain-openai pandas numpy scipy statsmodels matplotlib tqdm requests
    python -m benchmark.published_clock_2026_10_05.agentbench.prepare --bench DATA --root RUNS --python venv/bin/python --models gpt-6-astra --framework biomni
    python -m benchmark.published_clock_2026_10_05.agentbench.ranked prepare --final final-2026-10-10 --root RANKED --python venv/bin/python --geo-cache GEO --models gpt-6-astra --framework biomni
    bash benchmark/published_clock_2026_10_05/agentbench/run_biomni.sh RUNS venv 4
    python -m benchmark.published_clock_2026_10_05.agentbench.score --bench DATA --root RUNS --checks DATA/groupbench/design-checks/final-2026-10-10 --out DATA/agentbench/scored-biomni-001
    python -m benchmark.published_clock_2026_10_05.agentbench.ranked score --final final-2026-10-10 --root RANKED --out DATA/agentbench/scored-ranked-biomni-001
    python -m benchmark.published_clock_2026_10_05.agentbench.ranked combine --scored DATA/agentbench/scored-ranked-001 DATA/agentbench/scored-ranked-claude-001 DATA/agentbench/scored-ranked-biomni-001 --out results_2026_10_19/ranked_reruns_six_configurations.json
