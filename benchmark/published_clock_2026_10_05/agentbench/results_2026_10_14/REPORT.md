# Agent reruns of the ranked comparisons (protocol item 20)

The 60 comparisons of Tables 1 and 2 of verifier run final-2026-10-10 (30 with higher and 30 with lower clock outputs; 56 mouse RNA-seq series) were reanalysed by Codex with four language models, one run per series and model, with the released mouse task template and the inputs the historical task names. Every number below is computed from the run files by `ranked.py score`; nothing is simulated.

Runs: 224 planned, 224 exited 0, 0 reached the time limit, 0 excluded by the audit, 0 without a comparison, 0 web searches; 26,282,223 input tokens (20,867,484 cached) and 1,077,252 output tokens; median 68 s per run. Recomputation check: the statistic on the verified arms equals the released row for 60 of 60 comparisons. Codex CLI 0.162.1 (item 16 used 0.161.0; the CLI updated itself).

| Model | Higher: found | same direction | recurs | Lower: found | same direction | recurs |
|---|---|---|---|---|---|---|
| Codex + GPT-6 Astra | 28 | 28 | 28 | 29 | 29 | 28 |
| Codex + GPT-6 Sol | 27 | 27 | 27 | 27 | 27 | 25 |
| Codex + GPT-6 Luna | 25 | 25 | 25 | 27 | 27 | 26 |
| Codex + GPT-5.5 | 27 | 27 | 27 | 27 | 27 | 24 |

Across the four models: of the 30 comparisons with higher outputs, 28 were found by at least one, 28 recurred in at least two and 24 in all four; of the 30 with lower outputs, 29, 27 and 23. No found comparison had the opposite direction.

Comparisons that did not recur in all four models (found, recurs, of four):

| Table | Rank | Series | Comparison | Found | Recurs |
|---|---|---|---|---|---|
| higher | 4 | GSE198631 | Bleomycin vs Saline (all) | 3 | 3 |
| higher | 13 | GSE117868 | AD model with PU.1 KO vs Control (4.5 months) | 0 | 0 |
| higher | 17 | GSE146141 | TBK1 mut SOD_es vs TBK1 WT | 0 | 0 |
| higher | 20 | GSE137914 | Tet2/Tet3-deficient_CD93-CD23-CD21 Io vs control_CD93-CD23-CD21 Io | 3 | 3 |
| higher | 23 | GSE159135 | Lung: Alb-Cre, Bmal fl/fl (BLKO) vs Alb-Cre (AL feeding) | 3 | 3 |
| higher | 29 | GSE141729 | IWC 7d + infection 24 hpi vs Control | 2 | 2 |
| lower | 3 | GSE180384 | YF vs 3M-052 (Day D1, DCs) | 0 | 0 |
| lower | 5 | GSE195707 | Overload vs Sham (Vehicle) | 3 | 3 |
| lower | 14 | GSE155711 | ITE_MOG_L vs Vehicle in Astrocytes | 1 | 1 |
| lower | 24 | GSE100888 | Atg7 cKO vs Atg7 cWT (SOD1G93A background) | 4 | 3 |
| lower | 25 | GSE173510 | shBrd4 vs shCD19 | 4 | 0 |
| lower | 27 | GSE167360 | RUXO vs untreated | 4 | 2 |
| lower | 30 | GSE134151 | 5xFAD Hippocampus PLX5622 vs Vehicle | 2 | 2 |

Found: a comparison of the run whose arms overlap the verified arms with Jaccard indices of at least 0.5. Recurs: the verifier's statistic on the run's arms has the ranked direction and P < 0.05. Files: `ranked_reruns.json` (summary), `ranked_reruns.csv` (one row per model and ranked comparison), `runs.csv` (run audit), `TASK_example_GSE122080_gpt-6-astra.md`.
