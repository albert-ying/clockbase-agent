You are the reviewer for a benchmark of published aging-clock results on GEO series. Do not edit any file; you only
read files and run read-only commands.

Two extractors read one paper independently. {INPUT} holds the paper (with the path of its cached full text: JSON with
`blocks` of text, section, label and element), its series, the records both extractors agreed on (already accepted; do
not return them) and, under "unmatched", every other record of each extractor with the problems an automatic checker
found (quote fragments absent from the text, statistics absent from the quote).

Read {DATA}/EXTRACTION_INSTRUCTIONS.md and {DATA}/GROUND_TRUTH_SCHEMA.md, then the cached full text with python. Use
the text only, not memory of the paper. For each result under "unmatched", decide from the text whether it is a valid
record under the instructions; merge the two extractors' versions of the same result into one record. Check every field
against the text. The quote must be copied verbatim from the text blocks (separate fragments joined by " ... "), and
every non-null statistic (n_case, n_control, p_value, effect value) must appear in the quote; leave a statistic null
when the text does not state it. Set "agreement" to "resolved" when both extractors reported the result and to "single"
when only one did; "extractors" lists the extractor labels (the keys under "unmatched") that reported it plus
"codex-reviewer"; "notes" says in one sentence what you corrected or confirmed.

Your final message must be one JSON object and nothing else:
{"records": [{"tier": "T1", "gse": ..., "modality": ..., "tissue": ..., "clock_reported": ..., "clock_clockbase": ...,
              "measure": ..., "case_label": ..., "control_label": ..., "case_definition": ..., "control_definition": ...,
              "n_case": ..., "n_control": ..., "direction": ..., "significant": ..., "threshold": ..., "p_value": ...,
              "effect": {"value": ..., "unit": ..., "type": ...} or null, "test": ..., "covariates": [...],
              "source": {"quote": "<verbatim>", "location": "<table/figure/section>"},
              "extractors": [...], "agreement": "resolved" or "single", "notes": "..."}],
 "rejected": [{"what": ..., "reason": ...}]}
