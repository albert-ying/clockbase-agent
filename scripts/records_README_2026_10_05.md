# ClockBase verification records and recovered source tables

This deposit contains the unchanged `public-export/` tree from verifier run `{{FINAL}}` and the saved files cited by the manuscript. `supporting/` contains source tables, historical scripts, summaries and execution receipts. These are recovered records; the deposit builder does not run clocks, recompute statistics, repair rows or infer missing historical results.

`RECORDS_MANIFEST.tsv` records each file's relative path, SHA256 and byte count, including this README and `RECORDS.json`. It excludes its own hash. Read it with a CSV reader using a tab delimiter: paths containing tabs, newlines or quotes are CSV-quoted. `RECORDS.json` records the final run, creation time, copied-file count and size, source path and timestamp of each copied file, and the generated files' source. Its total file count includes both inventories and this README; its byte total covers the copied files and README only.

Every copied file is preserved byte for byte and checked against its source hash. All copied UTF-8 text, including compressed text tables, is scanned with the code-release builder's credential patterns; a match stops the build. The builder also reports local filesystem path references. Historical source paths inside dictionaries, scripts and receipts are retained for provenance. `RECORDS.json` maps them to the deposited files; readers must relocate historical script inputs before executing them. The copies do not establish the historical model version, acquisition roster, testing family or biological validity.

`questionaire.md` is preserved as used. Two headings understate their options: Q2.2 reads '(7 points)' and Q3.1 reads '(5 points)', but both offer options up to 10, matching the section maxima (25 and 10), the saved scores and Methods/Extended Data Figure 8 (0–10).

The index below names the manuscript statements and the deposited files that support them. The code release also contains `questionaire.md`, `minipaper_figures_2026-07-01/pipeline.py` and its README, `src/historical_geo_pipeline/analysis/`, and `scripts/legacy_sensitivity_analyses_2026_05_14/`.

| Manuscript statement id | Deposited path(s), relative to this directory | What the files show |
|---|---|---|
{{STATEMENTS}}

The six census source tables describe prediction and coordinate records. Coordinate records must not be added to prediction records as independent samples. The historical mouse DNAm count of 1,749 remains a source-code constant without a recovered roster. The Figure 3 decomposition preserves tied-record discrepancies and selection limits. The evaluator tables contain automated labels and notes; the builder does not convert them into independent reviewer judgments. The withholding overlay is labeled `WORKING_NOT_VALIDATED` in its source filenames and preserves that status here.

The inventory holds the cited numerical summaries and scripts, with their source receipts. No missing empirical data are replaced with simulations. The records builder refuses missing required inputs and an existing destination.

Build from the revision repository with:

```sh
python scripts/build_release_2026_10_04.py records --final FINAL_RUN --out NEW_DIRECTORY
```

The result is `NEW_DIRECTORY/records/`. To check it, parse `RECORDS_MANIFEST.tsv` with Python's `csv.DictReader(delimiter='\t')`, check each file's byte count and SHA256, and compare the copied-file hashes against `RECORDS.json`. Publication or upload of this local deposit is a separate step.
