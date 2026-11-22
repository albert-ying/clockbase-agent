# ClockEffectBench: agent reruns with other language models (2026-10-11)

Question: does the discovery task give results of similar accuracy when another agent framework and other language
models perform it? Rules: `../../groupbench/PROTOCOL.md`, addendum item 16, committed before any rerun was scored.
All values are computed by `../score.py` from the run outputs; nothing is simulated.

- Series and results: 61 human DNA methylation series with a published result on one of the five clocks of
  the task, groups assigned identically by two models, and a saved trace of the historical agent; 158
  published results from 47 papers.
- Agents: Codex CLI 0.161.0 through the OpenAI API with GPT-6 Astra, GPT-6 Sol, GPT-6 Luna, GPT-5.5, reasoning effort medium,
  in a Codex home without user instructions, plugins, memories or web search (`codex_home_config.toml`). One run per
  series and model: 244 runs.
- Task: the historical task text of each series, recovered from its saved trace, with one added deliverable (the
  samples of each arm); an example is `TASK_example_GSE25395_gpt-6-astra.md`.
- Run audit (`runs.csv`): every run wrote its result files; 1 reached the 30-minute limit after doing
  so and is scored; 1 is excluded because a command named a path outside its run directory (a Codex skills
  folder); no run made a web search, spawned a sub-agent or attempted a network command; 13 runs reported no
  testable comparison. Token use: 32.2 million input tokens (27.0 million cached) and
  1.41 million output tokens; median run time 63 seconds.

## Each system on the results it covers (158 published results)

| System | Results covered | Call agrees (n) | Direction agrees, published significant (n) | Direction AUROC (95% CI) |
|---|---|---|---|---|
| Historical agent (CodingAgent, Claude 3.7 Sonnet), stored result | 88 | 82% (74) | 87% (30) | 0.97 (0.87-1.00) |
| Released row of the historical agent | 93 | 85% (78) | 94% (33) | 0.98 (0.90-1.00) |
| Codex with GPT-6 Astra, reported result | 84 | 83% (70) | 93% (27) | 0.98 (0.88-1.00) |
| Codex with GPT-6 Sol, reported result | 86 | 90% (72) | 93% (27) | 0.99 (0.96-1.00) |
| Codex with GPT-6 Luna, reported result | 92 | 77% (75) | 90% (30) | 0.98 (0.89-1.00) |
| Codex with GPT-5.5, reported result | 100 | 79% (82) | 94% (34) | 0.98 (0.91-1.00) |
| Codex with GPT-6 Astra, verifier recomputation | 84 | 80% (70) | 93% (27) | 0.97 (0.86-1.00) |
| Codex with GPT-6 Sol, verifier recomputation | 86 | 86% (70) | 93% (27) | 0.99 (0.96-1.00) |
| Codex with GPT-6 Luna, verifier recomputation | 93 | 82% (73) | 94% (31) | 0.98 (0.89-1.00) |
| Codex with GPT-5.5, verifier recomputation | 101 | 82% (79) | 94% (34) | 0.98 (0.91-1.00) |
| GPT-6 Astra without data | 143 | 70% (114) | 60% (43) | n/a |
| GPT-6 Sol without data | 143 | 68% (114) | 60% (43) | n/a |
| GPT-6 Luna without data | 143 | 63% (114) | 53% (43) | n/a |
| GPT-5.5 without data | 143 | 67% (114) | 70% (43) | n/a |
| GPT-4.1 without data | 143 | 62% (114) | 72% (43) | n/a |
| Paper's groups, ClockBase clock outputs | 99 | 74% (73) | 94% (35) | 0.95 (0.89-1.00) |
| Constant rule (higher, significant) | 158 | 32% (128) | 80% (51) | n/a |

## The 68 results covered by the historical agent, its released row and every rerun (28 papers)

| System | Results covered | Call agrees (n) | Direction agrees, published significant (n) | Direction AUROC (95% CI) |
|---|---|---|---|---|
| Historical agent (CodingAgent, Claude 3.7 Sonnet), stored result | 68 | 83% (59) | 86% (21) | 0.98 (0.88-1.00) |
| Released row of the historical agent | 68 | 85% (59) | 95% (21) | 1.00 (1.00-1.00) |
| Codex with GPT-6 Astra, reported result | 68 | 85% (59) | 95% (21) | 1.00 (1.00-1.00) |
| Codex with GPT-6 Sol, reported result | 68 | 92% (59) | 95% (21) | 1.00 (1.00-1.00) |
| Codex with GPT-6 Luna, reported result | 68 | 81% (59) | 90% (21) | 1.00 (1.00-1.00) |
| Codex with GPT-5.5, reported result | 68 | 83% (59) | 95% (21) | 1.00 (1.00-1.00) |
| Codex with GPT-6 Astra, verifier recomputation | 68 | 83% (59) | 95% (21) | 1.00 (1.00-1.00) |
| Codex with GPT-6 Sol, verifier recomputation | 68 | 86% (59) | 95% (21) | 1.00 (1.00-1.00) |
| Codex with GPT-6 Luna, verifier recomputation | 68 | 83% (59) | 95% (21) | 1.00 (1.00-1.00) |
| Codex with GPT-5.5, verifier recomputation | 68 | 83% (59) | 95% (21) | 1.00 (1.00-1.00) |
| GPT-6 Astra without data | 64 | 73% (55) | 63% (19) | n/a |
| GPT-6 Sol without data | 64 | 71% (55) | 68% (19) | n/a |
| GPT-6 Luna without data | 64 | 62% (55) | 58% (19) | n/a |
| GPT-5.5 without data | 64 | 67% (55) | 68% (19) | n/a |
| GPT-4.1 without data | 64 | 58% (55) | 68% (19) | n/a |
| Constant rule (higher, significant) | 68 | 25% (59) | 71% (21) | n/a |

## Reading

- Under another agent framework and four other language models, the task agrees with published results about as
  often as the historical agent's stored results did, in calls and in direction.
- The same models without the data are at or below the constant rule in direction.
- The verifier's recomputation on a rerun's saved arms changes its agreement little, in either direction.
- Intervals are wide: results cluster within papers.

## Limits and notes

- Published-results part only: DNA methylation series, five clocks. The annotation part and RNA-seq were not rerun,
  and the historical CodingAgent itself was not run with another model.
- A rerun's comparison is matched to a published result by sample overlap alone (Jaccard index of at least 0.5 for
  each arm); the historical agent's row keeps the match of the published-clock benchmark.
- The recomputation uses the verifier's clock values and design-respecting estimand on the saved arms; membership
  recovery, orientation adjudication and input checks of the full verification are not applied.
- Audit detail fixed after the first runs and before scoring was read: a string such as `/clock_output/x.csv` that an
  agent joins to a variable is not counted as a path outside the run directory; a path counts when its first
  component is a directory at the filesystem root.
- Agreement metrics use the paper bootstrap of the published-clock summary (1,000 replicates, seed 20261005); the
  AUROC uses that of item 12 (2,000 replicates, seed 20261009).
