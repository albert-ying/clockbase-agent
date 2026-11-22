You are one of two independent extractors for a benchmark of published aging-clock results on GEO series. Do not look
for any other extractor's output. Do not edit any file; you only read files and run read-only commands.

Read the instructions {DATA}/EXTRACTION_INSTRUCTIONS.md and the record schema {DATA}/GROUND_TRUTH_SCHEMA.md, then the
input {INPUT}: the paper, the path of its cached full text (JSON with `blocks`: text, section, label, element) and the
GEO series to extract with the ClockBase atlas clocks of each.

Work table-first: read every table and its caption, then the Results, then the Methods. Use python to search the cached
text for each series accession, its cohort or dataset name and the clock results. Use the cached text only; do not use
memory of the paper. Copy quotes verbatim from the text blocks. Many papers will have no per-series group comparison
for the listed series (for example, the series served only as training or validation data for a clock); then return
no records and say why under `excluded`.

Your final message must be one JSON object and nothing else:
{"paper_id": "<paper_id from the input>",
 "records": [{"tier": "T1", "gse": ..., "modality": ..., "tissue": ..., "clock_reported": ..., "clock_clockbase": ...,
              "measure": ..., "case_label": ..., "control_label": ..., "case_definition": ..., "control_definition": ...,
              "n_case": ..., "n_control": ..., "direction": ..., "significant": ..., "threshold": ..., "p_value": ...,
              "effect": {"value": ..., "unit": ..., "type": ...} or null, "test": ..., "covariates": [...],
              "source": {"quote": "<verbatim>", "location": "<table/figure/section>"}}],
 "excluded": [{"gse": ..., "what": ..., "reason": ...}]}
Field meanings and allowed values are those of GROUND_TRUTH_SCHEMA.md. Every non-null statistic must appear in its quote.
