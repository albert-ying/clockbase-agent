# Published-clock comparison engine

Run from the verifier repository root with Python 3.12:

```bash
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m \
  benchmark.published_clock_2026_10_05.engine.cli run \
  --ground-truth ground_truth.jsonl --out engine-results \
  --groups-manual manual-groups.json --cache-only
```

The output directory must be new. Inputs and verifier data are read-only. Every data
root comes from `audit_system.paths`, including relocated rebuild bundles. The engine
adds no dependencies beyond the verifier's numpy, pandas and scipy environment.

## Inputs

Ground-truth JSONL follows `GROUND_TRUTH_SCHEMA.md` in
`/Entropy/kying/clockbase-benchmark-2026-10-05/`. Records need a unique `id`, quoted
source, tier, GSE, modality, mapped clock (possibly null), measure and both group
definitions. Extraction agreement remains the upstream collector's responsibility;
the engine does not certify a single extraction or supply missing published statistics.

An optional manual-groups JSON object maps record IDs to lists:

```json
{
  "paper2020-GSE123-Horvath-1": {
    "case_gsms": ["GSM1", "GSM2"],
    "control_gsms": ["GSM3", "GSM4"],
    "provenance": "Source section and literal metadata rule used to make these lists"
  }
}
```

Manual assignments are marked `MANUAL`; they do not claim dual-model agreement.
Overlapping arms, duplicate IDs, empty arms and IDs absent from metadata are unresolved.

S4 is optional. Pass `--s4-predictions FILE`, or place `s4_predictions.jsonl` beside
ground truth. Each line has `ground_truth_id`, `model`, `status`, `direction`,
`significant` and optional `confidence`, as exported by `llm_baseline/scorer.py`.
`id` and `model_requested` are accepted aliases. Nested `prediction` or `parsed`
objects are also accepted. Each model remains a separate `S4:<model>` system;
failed predictions remain uncovered. Duplicate record/model entries are rejected.

## Prediction sources and provenance

The dated public export under `paths.VERIFIER_ROOT/final-2026-10-05/public-export`
contains study JSON files and the prediction receipts used by the verifier. An
unambiguous exported source takes precedence. Replay inputs use the imported
`row_repairs.study_values`; repaired DNAm uses the imported `repaired_dnam_values`,
including `REPAIRED` status and `RETAINED` QC filtering. Missing sources or conflicting
receipts are unresolved; a stored hash mismatch prevents S3 recomputation.

Without an exported source, repaired DNAm is checked in empty-then-detached root
order. Historical predictions delegate clock-column selection to the imported
`membership_recovery._prediction_sources`: DNAm raw and QCed tables, and RNA
per-study `clock_result` files (including `_mean` columns). Recovered raw/QCed
source labels are honored; otherwise QCed DNAm precedes raw, as in replay inputs.
GrimAge V1/V2 names follow membership recovery. Duplicate prediction IDs are dropped
in full, rather than selecting an arbitrary duplicate.

Recorded age and sex come from processed metadata, preferring `harmonized.age` and
`harmonized.sex`. Ages must be numeric; strings or missing ages are not inferred.
Prediction, metadata, export, QC, status, ground-truth, manual-group and code receipts
include SHA256 and file modification time.

## Grouping and matching

Two distinct models from `BENCH_GROUP_MODELS` (default `gp-6-as,cl-op-5`) independently
propose literal field-to-allowed-values rules for each arm. Fields combine with AND;
values within a field combine with OR. Unknown fields or values fail validation.
Different rules may be accepted only if their selected case and control GSM sets
are identical. Both selections and raw answers are retained on disagreement.

The client imports the verifier's `call` and `Ledger`, preserving request hashing,
temperature, raw responses, token usage, budget and `AUDIT_CACHE_ONLY=1` behavior.
The default request budget is zero. `--cache-only` prevents transport regardless of
budget. `BENCH_GROUP_LEDGER` names an existing read-only ledger for cache replay.
Positive `--max-requests` without cache-only enables audited transport and needs
the verifier's configured credentials; none of the tests or smoke runs use it.

S1/S2 candidates require the same GSE and clock. Both orientations are scored by
Jaccard overlap; each arm must reach 0.5. Ranking uses minimum overlap, then summed
overlap, then prefers a representative over a tied duplicate and uses row ordinal
for a stable tie-break. Every candidate and both orientation scores are retained.
With unresolved S3 groups, two models must agree on SAME or SWAPPED from textual
descriptions. Their judgments never substitute a sample assignment for S3.

## Systems and statistics

- **S1:** the closest matched historical atlas row. Coverage means a row matched;
  missing stored statistics still remain missing. Historical pooled Cohen's d is
  compared with a published d only when its recovered effect definition confirms it.
- **S2:** the same row's recomputed values when its class is `VERIFIED`, `PROVISIONAL`
  or `DESCRIPTIVE`. Pending or withheld rows are uncovered. Descriptive rows have
  no inferential P value. Duplicate rows are not counted as independent results.
- **S3:** ClockBase outputs tested using published sample groups. `age` and `rate`
  use the output directly; `acceleration_difference` subtracts recorded age;
  `acceleration_residual` fits OLS with an intercept on all series samples with finite
  predictions and recorded ages. Fit coefficients and sample IDs are recorded.
- **S4:** each supplied language model's no-data predictions, kept separately.

S3 imports the verifier's Welch two-sided t-test and pooled-SD Cohen's d routine.
Effects and mean differences are case minus control. An arm with fewer than three
usable samples is descriptive: means and eligible d remain, P/significance do not.
The tested sample counts and excluded GSMs are explicit. Significance uses a numeric
published threshold (`P < 0.05`, `p ≤ .05`, `5e-2`); an absent or unparseable threshold
has no significance decision. FDR thresholds are not interpreted as raw-P thresholds.

S1/S2 usually compare clock age; a published acceleration measure receives an
`estimand_mismatch` flag. Matching does not repair that mismatch. S3 controls only
chronological age when residuals are requested: pairing, family effects, cell
composition and other published adjustments are not reproduced. Control-only
residual fits and rank tests in papers can therefore differ from S3.

## Outputs and metrics

- `records.jsonl`: published records, source receipts, model answers and selections,
  all atlas candidates, oriented system results and S3 residual-fit provenance.
- `per-record.csv`: one contrast per row, with published and all system results.
- `summary.json`: separate T1/T2 summaries overall and by modality, clock family
  and source. Unknown source collections use the publication citation as the key.
- `report.md`: coverage and agreement overview, denominators, intervals and limits.
- `manifest.json`: inputs, code hashes, configuration and requests sent.
- `ledger.jsonl`: raw live responses or cache-miss receipts when model work is attempted.

Direction agreement is among published-significant contrasts with an available
system direction. Significance agreement requires both decisions; non-significant
agreement is also reported. Descriptive results can provide a direction. Effect
Spearman correlations require at least three pairs, matching measures and compatible
effect definitions/units. Betas and unstandardized values are not relabeled as d.
Different unit groups have separate correlations; no pooled rank correlation mixes
them. Confidence Brier score treats confidence as probability of a correct direction;
mean confidence minus accuracy is also reported. A `post_cutoff=true` record is
scored separately, but the engine does not infer training cutoffs from model names.

95% percentile intervals resample papers, identified by DOI, PMID, PMCID, citation
or URL, with all their contrasts together. Seed defaults to 20261005; use `--seed`
and `--bootstrap` to override. Fewer than two papers or undefined bootstrap metrics
give null interval bounds. JSON uses null for undefined statistics, never NaN.

## Validation

```bash
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m unittest \
  benchmark.published_clock_2026_10_05.engine.tests.test_engine
```

Tests use labeled synthetic constants for software verification, and a fake HTTP
transport to verify the imported client/ledger. They make no network/model calls.
The empirical smoke inputs, source receipt, deterministic group-list preparation
and results are under `/Entropy/kying/clockbase-benchmark-2026-10-05/engine-smoke/`.
Smoke records are single-extractor checks, separate from the collector's agreed
ground-truth corpus. A two-series smoke run does not establish benchmark accuracy.
