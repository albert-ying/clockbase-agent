# S4 language models without data

S4 predicts a published contrast from a GEO study description, group definitions,
sample counts, a neutral clock description, and the measure. It reads no clock
outputs or sample measurements. It uses the ground-truth schema in
`/Entropy/kying/clockbase-benchmark-2026-10-05/GROUND_TRUTH_SCHEMA.md`.

```bash
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m \
  benchmark.published_clock_2026_10_05.llm_baseline.cli run \
  --ground-truth FILE --out DIR --limit 10
```

Run live requests later through the credential launcher. This implementation was
checked with a fake HTTP transport only. `BENCH_S4_MODELS` defaults to
`gp-6-as,cl-op-5,gp-4o,deepseek-v4-pro,qwen3.8-max,ge-3.1-pr-p`;
`--models` accepts a comma-separated override. The audited client's gateway branch
sends temperature 0. Direct `openai:` names (for example `openai:gpt-6-astra`) go
to the OpenAI API without a temperature, because OpenAI reasoning models accept
none, and with `max_completion_tokens` 16000 including reasoning; each answer is
sampled once and replayed from the ledger. `manifest.json` records both settings
under `sampling`. The prompt requests JSON; the client has no separate JSON
response-format option.

`BENCH_S4_MAX_REQUESTS` defaults to 3000 and `--max-requests` overrides it.
The budget counts logical client requests cumulatively across this output
directory, including failures; cache hits consume none. The shared client may
retry a transport failure up to three times within one logical request.
`AUDIT_CACHE_ONLY=1` prevents all transport, including retries, and produces
`NOT_CACHED` rows for misses. CLI exit codes are 0 for a complete run, 2 for
incomplete predictions, and 1 for an input or file error.

## Inputs and blinding

`audit_system.paths.GEO` supplies `gse_me_hs.csv`, `gse_rna_hs.csv`, or
`gse_rna_mm.csv`. The loader selects `accession`, `title`, `description`, and
`samples`. Missing descriptions or series counts trigger a local series-matrix
header fallback (`!Series_title`, `!Series_summary`, `!Series_overall_design`,
and unique `!Sample_geo_accession` values); reading stops before the matrix.
No network retrieval occurs. A missing required study description stops the run
with `Cannot proceed - required data file not found: ...`. Unreported counts stay
unreported; case and control counts are not treated as the total series count.

Only explicit input fields reach the prompt. Record IDs, citations, publication
metadata, result fields, and source quotes are omitted. Reported clock names map
to the small neutral dictionary in `prompts.py`; unknown clocks or measures
produce `PROMPT_FAILED` rows without model requests. Supported measures are clock
output, residual on chronological age, clock minus chronological age, and pace.

GEO sentences containing source text, known paper titles, quoted statistics,
thresholds, or result wording are withheld and flagged. Exact numeric result
values that collide with a sample count cause that count to be withheld.
`check_blinding()` checks the study-specific stimulus. The fixed task necessarily
includes the answer options and the target P < 0.05, so it is excluded from the
check; these are not the record's published answer or threshold.

Case and control definitions containing directional or significance wording
remain in the prompt with `case_control_answer_leakage`. This flags an unavoidable
exception to stimulus blinding, including phenotype definitions such as
accelerated aging. An exact source quote or numeric result in a group definition
produces a flagged `PROMPT_FAILED` row instead of sending the prompt. The record
stays in the output. This lexical check cannot certify the absence of indirect
clues or a model's recall of a publication.

## Outputs and replay

The output directory contains:

| File | Contents |
| --- | --- |
| `manifest.json` | Ground-truth path, SHA256, modification timestamp, schema, completion limit |
| `requests.jsonl` | Canonical request body and hash, prompt flags, source paths/hashes/timestamps, token estimate |
| `ledger.jsonl` | Audited raw response, requested/returned model, usage, parse result, status, attempts |
| `predictions.jsonl` | Append-only record/model attempts and normalized predictions |
| `cost_summary.json` | Tokens, requests, retries, failures and missing usage per requested model |
| `s4_predictions.jsonl` | Latest record/model attempt in the engine's input format |

Request hashes match `audit_system.label_adjudication_2026_09_27.call` exactly.
`S4Ledger` subclasses its `Ledger`, including caching of returned unparseable
answers so a restart cannot purchase a replacement answer. Transport failures
can be retried by a later run subject to the cumulative budget. Raw-response
references are `ledger.jsonl#request_sha256=<hash>`; use the last matching ledger
entry when a transport failure was followed by a later response. Responses from
cache do not increase usage totals. Dollar costs remain null because no provider
price schedule was supplied.

`predictions.jsonl` records `record_id`, `model`, `status`, `prediction`
(`higher`, `lower`, `none`), `significance` (boolean), `confidence` (0–1),
`recognizes_publication` (boolean), `publication`, `raw_response_reference`,
`request_hash`, `prompt_flags`, and `cached`. Unavailable predictions have null
answer fields. Recognition requires an identified publication. Invalid JSON,
duplicate keys, missing answers, out-of-range confidence, and contradictory
answers are retained as `UNPARSEABLE`; no values are imputed.

The scorer exports flat objects with `ground_truth_id`, `model`, `system: "S4"`,
`status`, `direction`, `significant`, `confidence`, `recognizes_publication`,
`publication`, `request_hash`, `raw_response_reference`, and `prompt_flags`.
This matches `engine/README.md` and `engine.cli.load_s4`, checked directly in the tests. Failed rows
are retained and the engine marks them uncovered. Pass
`--s4-predictions DIR/s4_predictions.jsonl` to the engine; its default discovery
looks beside the ground-truth file.

The same output directory can resume the same ground-truth file. Requests,
ledger entries and prediction attempts append; summaries and the engine export
are regenerated from those files. Existing outputs without an S4 manifest and
changed ground-truth identity are refused. `--limit` counts records, not requests;
each selected record is paired with each model. A smaller later selection does
not remove earlier rows from the export.

## Checks and example

```bash
TMPDIR=/Entropy/kying/tmp .venv-pipeline/bin/python -B -m unittest \
  benchmark.published_clock_2026_10_05.llm_baseline.test_baseline
```

Tests use synthetic text fixtures and a fake transport, including prompt
blinding, clock mapping, matrix-header fallback, JSON/fenced JSON parsing,
malformed responses, ledger caching, budget stops, cache-only hits and misses,
failure accounting, output preservation, and the engine reader. Test files are
created and removed inside this package. No benchmark results are simulated.

The following prompt uses `curated-00001` in
`/Entropy/kying/clockbase-benchmark-2026-10-05/curated/records.jsonl` (GSE25395).
The extraction table uses `record_id`, `gse_id`, `clock_source`, and `quote`;
this example maps them to the ground-truth schema's `id`, `gse`, `clock_reported`,
and `source.quote`. Its note explicitly identifies the residual measure. This
mapping is only for the example check; the run CLI expects the ground-truth
schema. Series count 60 comes from the GEO index; the published case count is
35 and the control count is unreported. The published subset's total is not
substituted for the GEO series count.

The GEO description source is `geo:gse_me_hs.csv`, SHA256
`fc49595b4e13e3e18b4a12d82bc0234b24b5a2c776c175fb7311009056b4eefa`,
modification time `1746437099773362339` nanoseconds since the Unix epoch.

```text
GEO series: GSE25395
Modality: human DNAm
Title: Altered DNA Methylation in Leukocytes with Trisomy 21
Summary: The primary abnormality in Down syndrome (DS), trisomy 21, is well known, but how this chromosomal gain produces the complex DS phenotype, including immune system defects, is not well understood. We profiled DNA methylation in total peripheral blood leukocytes (PBL) and T-lymphocytes from adults with DS and normal controls and found gene-specific abnormalities of CpG methylation in DS, with many of the differentially methylated genes having known or predicted roles in lymphocyte development and function.<br><br>NOTE 1: the analysis in Kerkel et al., PLOS Genet, 2010 was carried out after removing probes for genes on the X and Y chromosomes. NOTE 2: these GEO data have been analyzed with this latest version of BeadStudio; the data in the supplementary tables of Kerkel et al., PLOS Genetics, 2010, were from these same Bead chips, but extracted using an earlier version of this software. This re-extraction does not affect the conclusions in Kerkel et al., PLOS Genetics, 2010.
Overall design: [not supplied]
Clock: Horvath pan-tissue age. Estimates chronological age from DNA methylation across tissues.
Measure: Residual of clock output regressed on chronological age.
Series samples: 60
Case samples: 35
Control samples: not reported
Case definition: Down syndrome (trisomy 21)
Control definition: Non-Down-syndrome samples in this dataset

Predict the case group relative to the control group for the specified clock measure in this series.
Use the study description and your prior knowledge only. Do not access sample data, tools or external sources.
Predict higher, lower, or none (no difference), and whether the difference would be significant at P < 0.05.
Confidence is your probability (0 to 1) that the direction prediction is correct.
Say whether you recognize a publication reporting this clock result for this GEO series; if yes, identify it.
Return JSON only:
{"prediction": "higher|lower|none", "significance": true,
 "confidence": 0.5, "recognizes_publication": false, "publication": null}
Use a boolean for significance and recognition. If recognition is true, publication must identify the paper.
Choose each value independently; the example values do not describe this study.
```

This prompt is approximately **579 input tokens**, estimated as the ceiling of
characters divided by four. This is a planning estimate, not a tokenizer count
or billed usage. A compact answer will often be about 80–150 output tokens;
the completion cap is 400. Longer GEO descriptions increase the input estimate,
which is saved for every request. Ledger usage supplies the reported token counts
after a live run.
