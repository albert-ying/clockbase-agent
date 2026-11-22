# Historical GEO pipeline source snapshot

Copied on 2026-10-04 from `/Entropy/kying/GEO_data/src`, preserving relative paths and bytes for every `.py`, `.R`, `.r`, `.sh` and `.ipynb` file present there. `MANIFEST.tsv` records each original path, SHA256, size and source modification time (UTC and nanoseconds). All 19 copies match their source hashes. No data, outputs or bytecode caches are included.

This is a snapshot of the source tree available on the copy date. Some producers and the two guard-test files include September 2026 corrections; the snapshot does not claim to recover the exact versions used for every historical prediction. Historical paths, dependencies, incomplete stubs and interactive prototypes are retained. Run the revision pipeline with its documented inputs and dependencies; these copies alone are not a self-contained rebuild.

`src/dnam_stream_repair_2026_09_26.py` and `src/pipeline_inputs_2026_09_10.py` are the revision's re-implementations. They remain the revision entry points for DNAm recovery and input acquisition/adapters; the copied producers preserve the prior GEO workflow for inspection.

Before copying, all source text was scanned for API-key, secret, token and password terms, recognizable key formats and private-key markers. No credential values were found, so no values were replaced with `REDACTED`. The `token` matches in the DNAm guard test are local numeric parser-test variables.

| Script | Role read from the code |
|---|---|
| `0_get_meta_data.R` | Retrieved GEO dataset descriptions and designs, selected human/mouse RNA and DNAm accessions, downloaded series matrices and extracted sample metadata. |
| `1.2_extract_human_matrix.R` | Extracted human DNAm feature matrices from GEO series matrices into compressed per-study CSVs. |
| `1_process_meta_data.R` | Extracted and harmonized sample fields for all four modalities, retaining pre-harmonization metadata alongside harmonized CSVs. |
| `2_download_supplementary_data.R` | Downloaded and unpacked supplementary RAW archives for human and mouse RNA, with mouse DNAm download code retained as comments. |
| `3.1_make_list_for_succeed_files.R` | Combined accession lists from ARCHS4 and supplementary RNA processing into human and mouse success lists. |
| `3_supplementary_file_compile.R` | Parsed supplementary RNA counts, mapped gene identifiers and wrote scaled and normalized expression matrices for mouse and human studies. |
| `4_compute_age.R` | Applied ClockBasis human DNAm clocks to series matrices and wrote per-study sample-age CSVs. |
| `4_compute_age.py` | Computed biolearn DNAm clock predictions after GSE40279-reference imputation and exported predictions, QC and source-identity records. |
| `4_compute_missing.py` | Calculated clock CpG missingness and summarized imputation, prediction, metadata and intervention artifacts with source logs. |
| `4c_mortality_trans_clock.R` | Applied multi-tissue and blood mortality-clock coefficient sets to human DNAm series matrices. |
| `5_create_meta_annot_for_train.R` | Prepared human DNAm descriptions, extracted metadata and target annotation tables for metadata-model training. |
| `6_intervention_analysis_human_dnam.R` | Defined the harmonized-output folder in a three-line intervention-analysis stub; no analysis is implemented. |
| `old_4a_compute_causal_age.R` | Applied the archived CausAge, DamAge and AdaptAge coefficient sets to human DNAm matrices. |
| `old_4b_compute_grimage.R` | Computed archived GrimAge v1/v2 predictions and missingness from DNAm protein models and cohort-calibration parameters. |
| `tAge.py` | Applied 18 mouse and two human transcriptomic clock models to scaled RNA expression and exported per-study ages. |
| `temp.py` | Prototyped per-clock methylation-site missingness summaries using biolearn QC reports; depends on interactive variables. |
| `temp_mv_file.R` | Combined per-sample mouse DNAm files into per-study feature tables. |
| `tests/test_clock_input_guards_2026_09_10.py` | Tested later DNAm empty-input, GEO parsing, source-identity and R input guards without running clocks. |
| `tests/test_rna_identity_guards_2026_09_10.py` | Tested later RNA feature/sample identity and clock-input validation without applying trained models. |
