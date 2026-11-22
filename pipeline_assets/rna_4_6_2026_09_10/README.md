# ClockBase RNA producer assets

These are preserved copies of the inspected GEO producer and its fitted 4.6
models. `manifest.json` records source paths, timestamps, checksums and the GEO
producer commit. Runtime code never follows those server-specific paths.

The new entry point uses the existing `CodingAgent` and the primary analysis
clocks named in `src/parallel_mm_rna.py` and `src/parallel_hs_rna.py`: two mouse
Hazard clocks and human Hazard2, with scaled differential preprocessing. The
additional All Hazard asset is retained for explicit comparison, not selected
as a silent substitute. Other historical clocks and YuGene are outside this
entry point's validated clock profile.

The R transform and producer Python file are unchanged copies. The pipeline
calls their selected functions, checks exact sample identity and observation
coverage, and independently verifies both preprocessing and prediction. The
historical `1_preprocess_rna.R` wrapper is provenance only; it is never executed
as a global corpus job.

The four fitted models (`EN_*.pkl`; names, sizes and checksums in `manifest.json`) are
distributed separately from the code release. Place them in this folder before
running the RNA clock stage. Only load model files that match the manifest:
pickle loading executes code.
