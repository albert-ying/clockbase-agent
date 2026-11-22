# Recovered atlas analysis source

These five R scripts are byte-identical copies of `GEO_data/src_analysis/*.R`. Tracked files are read from `git show HEAD:<path>`; files absent from Git are copied from the working tree and labeled as untracked snapshots. `MANIFEST.tsv` records the origin, SHA256, size, source timestamp and each script's role; the parent manifest contains the same rows. The copies preserve historical paths and dependencies. Source recovery does not establish the producing runtime or validate the outputs.

`0_clean_up_data_hs_dnam.R` selects and masks the GrimAge components, `1_dnam_analysis.R` filters DNAm predictors and builds coordinates and blood-panel summaries, `2_clean_up_data_rna.R` cleans RNA inputs, and `3_rna_analysis.R` builds RNA coordinates and plots. `summarize_plot.R` contains the mouse DNAm count of 1,749 as a constant; a source roster for that count remains unavailable.
