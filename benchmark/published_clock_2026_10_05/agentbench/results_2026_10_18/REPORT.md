# Claude Code reruns (PROTOCOL.md items 22 and 23, 2026-10-18)

Second agent framework: Claude Code CLI 2.1.302 with Claude Opus 5.5, in an isolated
configuration directory (credentials only) and a bubblewrap sandbox with the filesystem read-only outside the run
directory. Same tasks, inputs and scoring as the Codex reruns of items 16 and 20.

| Set | Runs | Excluded | Result |
|---|---:|---:|---|
| Published results (61 series) | 61 | 0 | call agreement 0.83 (n 87), direction agreement 0.95 (n 37), 107 of 158 results covered |
| Ranked comparisons (56 series) | 56 | 0 | found 52 of 60, all in the ranked direction; 51 recur |

Across Codex (GPT-6 Astra, GPT-6 Sol, GPT-6 Luna, GPT-5.5) and Claude Code (Claude Opus 5.5): 55 of the 60 ranked
comparisons recur in at least two configurations and 45 in all five; none is reversed
(`ranked_reruns_five_configurations.json`).

Files: `agent_reruns.json`, `runs.csv`, `records.jsonl` (published results); `ranked_reruns.json`,
`ranked_reruns.csv`, `ranked_runs.csv` (ranked comparisons); `TASK_example_GSE122080_claude-opus-5-5.md`.

Commands:

    python -m benchmark.published_clock_2026_10_05.agentbench.prepare --bench DATA --root RUNS --python PYTHON --models claude-opus-5-5 --framework claude
    python -m benchmark.published_clock_2026_10_05.agentbench.ranked prepare --final final-2026-10-10 --root RANKED --python PYTHON --geo-cache GEO --models claude-opus-5-5 --framework claude
    bash benchmark/published_clock_2026_10_05/agentbench/run_claude.sh RUNS CLAUDE_HOME 2
    python -m benchmark.published_clock_2026_10_05.agentbench.score --bench DATA --root RUNS --checks DATA/groupbench/design-checks/final-2026-10-10 --out DATA/agentbench/scored-claude-001
    python -m benchmark.published_clock_2026_10_05.agentbench.ranked score --final final-2026-10-10 --root RANKED --out DATA/agentbench/scored-ranked-claude-001
    python -m benchmark.published_clock_2026_10_05.agentbench.ranked combine --scored DATA/agentbench/scored-ranked-001 DATA/agentbench/scored-ranked-claude-001 --out results_2026_10_18/ranked_reruns_five_configurations.json
