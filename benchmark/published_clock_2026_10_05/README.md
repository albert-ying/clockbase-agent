# Published aging-clock benchmark, 2026-10-05 (ClockEffectBench)

The manuscript reports this benchmark and the group-choice benchmark in `groupbench/` together as **ClockEffectBench**:
a benchmark of whether an analysis system recovers published aging-clock effects from public data. Its systems S3 and
S4 below are reported there as "ClockBase with the paper's groups" and "language models without data".

Ground truth: results that published papers report for aging clocks applied to group comparisons in GEO series that
the ClockBase atlas also analysed. Systems compared on each published contrast:

| System | What it is |
|---|---|
| S0 | constant rule: every contrast higher in the case and significant (the base rate) |
| S1 | the historical ClockBase agent row for the same series, clock and contrast |
| S2 | the same row after verification (recomputed statistics, final class) |
| S3 | the paper's groups, assigned from the GEO sample metadata by two models that must select identical samples, tested with ClockBase's per-sample clock outputs (Welch test) |
| S4 | language models without data: GEO description, group definitions, a neutral clock description and the measure; each predicts direction, significance, confidence and whether it recognizes the publication |

Data root: `/Entropy/kying/clockbase-benchmark-2026-10-05/` (`DATA` below). Design: `DATA/DESIGN.md`; record schema:
`DATA/GROUND_TRUTH_SCHEMA.md`; extraction rules: `DATA/EXTRACTION_INSTRUCTIONS.md`.

## Pipeline

Run from the repository root (`PY=.venv-pipeline/bin/python -B`, `TMPDIR=/Entropy/kying/tmp`). Model calls go through
the verifier's audited client (`audit_system/label_adjudication_2026_09_27.call`): every request hash, raw response,
usage and failure is kept in a `ledger.jsonl`, and a rerun replays answered requests instead of re-sampling them.

1. Retrieval (below): Europe PMC search, atlas intersection and triage of paper-series pairs into LIKELY, POSSIBLE and
   UNLIKELY (`retrieve.py`, `validate.py`). Curated sources: `curated/` (T2 expected directions from ComputAgeBench).
2. Extraction, LIKELY pairs (19 papers): two independent Claude extractions and a reconciliation (2026-10-05 session;
   `DATA/extraction/`), giving `DATA/ground_truth_t1_2026-10-05.jsonl` (168 records) and the T2 file (66 records).
3. Extraction, POSSIBLE pairs (379 papers with open full text; 96 pairs without it are listed in
   `DATA/extraction/possible/skipped.jsonl`):
   - `extraction.possible_inputs` writes one input per paper.
   - `extraction/run_codex_extract.sh` runs one read-only Codex extraction per paper (first extractor gpt-6-astra;
     second extractor gpt-6.1-sol on every paper with first-pass records plus 40 sampled papers without, seed 20261005).
   - `extraction.reconcile compare` checks every record against the cached text (quotes as verbatim pieces of the
     document, every copied statistic present in its quote) and accepts records identical in series, clock, measure,
     direction, P value, group sizes and effect; the rest go to `review/`.
   - The same driver with `PROMPT=codex_review_prompt.md` runs a Codex reviewer (gpt-6.1-sol, high effort) on each
     review file; `extraction.reconcile collect` keeps reviewer records that pass the same checks.
   - `extraction.reconcile assemble` adds agreed and reviewed records to the base ground truth, only for series listed
     for the paper, not already in the base and in the atlas modality of the series.
4. Audit: a stratified random sample of 40 new records (seed 20261005) checked against the text by Claude subagents
   (`DATA/audit/`); records with a major error are excluded (`DATA/audit/exclusions_after_audit_2026-10-05.json`),
   giving `DATA/ground_truth_t1_v2_final_2026-10-05.jsonl`.
5. `curated.expand_t2` expands each T2 expectation to one record per ClockBase clock of its series.
6. S1-S3: `engine.cli run --ground-truth GT --out RUN [--max-requests N]` with `BENCH_GROUP_MODELS` (two models).
   Group prompt v2 treats published sample counts as descriptive; v1 runs (`*-openai-001`) are kept.
7. S4: `llm_baseline.cli run --ground-truth GT --out DIR --models M` (one directory per model).
8. Scoring with S4 attached, without new requests: `AUDIT_CACHE_ONLY=1 BENCH_GROUP_LEDGER=RUN/ledger.jsonl engine.cli
   run ... --cache-only --s4-predictions MERGED`.
9. `summarize_2026_10_05 --t1 RUN... --t2 RUN --ground-truth FINAL --out DIR`: per-system tables, records shared by
   systems, the three-class call balanced accuracy, S4 by publication recognition, S3 against each model on
   publications it did not recognize, and T2 expected directions.

## Group-choice benchmark and the released run

`groupbench/` holds the second benchmark: the agent's choice of groups and its results against expert sample
annotation (Tyshkovskiy et al. 2026 for mouse RNA-seq, ComputAgeBench for human DNA methylation) and against the
published results above. Its protocol with every amendment is `groupbench/PROTOCOL.md`; the reports are
`groupbench/results_2026_10_06/REPORT.md` (census of design errors), `groupbench/results_2026_10_07/DESIGN_CHECKS.md`
(design checks and their validation) and `groupbench/results_2026_10_10/REPORT.md` (the released verifier run
final-2026-10-10). The stage "verifier, Welch test only" is the earlier run final-2026-10-05, which the rebuild bundle
does not ship; the released records keep the Welch result of every row tested otherwise (`recomputed_declared`).
The benchmark needs its data root (`DATA`), which is not part of the rebuild bundle; the records deposit holds the
extracted published results, the benchmark units and the per-row records. To score a verifier run:

```sh
python -m benchmark.published_clock_2026_10_05.groupbench.export_as_checks --scored final-2026-10-10 --out DATA/groupbench/design-checks/final-2026-10-10
python -m benchmark.published_clock_2026_10_05.groupbench.design_eval --root DATA/groupbench --checks DATA/groupbench/design-checks/final-2026-10-10 --out DATA/groupbench/design-checks/final-2026-10-10 --corrected   # also --external
python -m benchmark.published_clock_2026_10_05.groupbench.a_result --root DATA/groupbench --prepared DATA/groupbench/prepared-v2 --frozen DATA/groupbench/frozen-v2 --out DATA/groupbench/design-checks/final-2026-10-10/a_result_stages.jsonl --stages DATA/groupbench/design-checks/final-2026-10-10
python -m benchmark.published_clock_2026_10_05.groupbench.design_eval --root DATA/groupbench --checks DATA/groupbench/design-checks/final-2026-10-10 --out DATA/groupbench/design-checks/final-2026-10-10 --auroc
CLOCKBASE_BENCH_SCORED=final-2026-10-10 python -m benchmark.published_clock_2026_10_05.groupbench.figure_benchmark_2026_10_07 --bench DATA --out RESULTS --modality-table   # also --table, --tidy, --scorecard
```

## Agent reruns with other language models

`agentbench/` reruns the discovery task on the series with published results, with Codex as the agent framework and
four GPT models, and scores the runs against the published results (`groupbench/PROTOCOL.md`, addendum item 16;
report `agentbench/results_2026_10_11/REPORT.md`):

```sh
python -m benchmark.published_clock_2026_10_05.agentbench.prepare --bench DATA --root RUNS --python PYTHON
CLOCKBASE_MODEL_LAUNCHER="..." bash benchmark/published_clock_2026_10_05/agentbench/run_codex.sh RUNS CODEX_HOME 8
python -m benchmark.published_clock_2026_10_05.agentbench.score --bench DATA --root RUNS --checks DATA/groupbench/design-checks/final-2026-10-10 --out DATA/agentbench/scored-001
```

The same tasks were run with Claude Code and Claude Opus 5.5 as a second agent framework (`run_claude.sh`, a
bubblewrap sandbox with an isolated configuration directory; addendum items 22 and 23). The commands and results,
including the combination of the ranked-comparison reruns over the five configurations, are in
`agentbench/results_2026_10_18/REPORT.md`.

## Deviations from the design

- NeoLink returned HTTP 402 (payment required) for every model on 2026-10-05 (`DATA/runs/*-live-001`). All model calls
  then used the OpenAI API directly: group rules by gpt-6-astra and gpt-6-sol; S4 by gpt-6-astra, gpt-6-sol,
  gpt-6-luna, gpt-5.5 and gpt-4.1. The design's Claude S4 model and the S5 pilot (CodingAgent with other models) were
  not run. The direct OpenAI branch sends no temperature (reasoning models accept none); answers are sampled once and
  replayed from the ledger.
- The POSSIBLE extractions were made by two OpenAI models with an OpenAI reviewer; the cross-family check is the
  Claude audit of a sample.

## Retrieval

This package finds papers that mention aging clocks and GEO series represented in the
2026-10-05 ClockBase atlas. It screens paper/series pairs for a published group comparison;
it does not extract statistics or establish matching contrasts.

The assigned output folder is
`/Entropy/kying/clockbase-benchmark-2026-10-05/retrieval/`. Its generated `README.md`
contains the counts, method, limits and 30 most promising LIKELY pairs. `queries.json`
records exact queries and cursor pages; `papers.jsonl` records all retrieved papers;
`candidates.jsonl` records the complete atlas intersection with passages and locations;
`summary.json` records counts and source failures. Cached bytes and SHA256 metadata
remain under `cache/`; normalized prose remains under `texts/`. GEO accession
annotations and core records are also preserved separately. Duplicate MED/PMC annotation
records for an article are combined, so an empty record cannot erase another record’s GEO links. GPL/GSM mentions do not
enter the series intersection.

Curated database linkage uses NCBI ELink from PubMed to GEO and ESummary for the
linked GEO records. Repeated PMID parameters preserve separate per-paper link sets;
linked GDS records retain their originating GSE field. Europe PMC's databaseLinks
endpoint rejects GEO as a database filter, as documented by cached probes. Core
database-name annotations are retained alongside the NCBI cross-reference results.

Run from the repository root with Python 3.12:

```sh
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m benchmark.published_clock_2026_10_05.retrieve
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m unittest benchmark.published_clock_2026_10_05.test_retrieve
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m benchmark.published_clock_2026_10_05.validate
```

Add `--offline` to retrieval for a replay that makes no network requests. Both commands
accept `--atlas` and `--output`. Cached responses, including HTTP errors, are reused.
The client limits request starts to 2.5 per second and retries transient failures with
backoff. Search paging fails if it cannot reach the reported hit count. Per-paper
source failures are reported in the summary rather than treated as absent evidence.

Automatic LIKELY requires a clock-result sentence with comparison language and a
reported statistic or direction, plus the GSE, application language and groups
within one adjacent sentence in the same prose paragraph.
Training context, reviews, methods, introductions, supplementary descriptions,
data-availability sections and whole-table text cannot
produce automatic LIKELY. Nearby/unbound comparisons are POSSIBLE. A manual review
in `retrieval/manual_reviews.jsonl` can override the class with a reason and passages;
each passage must match the cached text verbatim, including its section and label.
Reviewed candidates retain their automatic classification. No paywalled full texts,
supplementary files or figure pixels are read by this step.

`validate` independently reads the cached search pages, verifies checksums, traces
every accession and candidate passage to its fetched source, and recomputes atlas
row/clock counts and the complete intersection. The unit tests use synthetic prose
and mocked HTTP responses to check failure handling and classification boundaries;
these fixtures do not enter the benchmark outputs.

Final validation on 2026-10-05: 19 unit tests passed; independent source readback passed 68,876 checks over 1,725 papers, 1,978 candidates and 30 ranked pairs. The final offline replay made zero network requests. All 55 LIKELY pairs were manually reviewed; 65 pairs have manual review records.
