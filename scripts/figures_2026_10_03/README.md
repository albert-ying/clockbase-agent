# Revision figure builders (2026-10-03)

These are the preserved evidence builders with portable inputs and a new output directory. Plotting code, labels,
selection rules and the shared plotting helper are unchanged; the helper adds only argument handling.

- `build_figure3_tables_2026_10_03.py` makes `figure3_ranking_2026_10_03.pdf`, Tables 1 and 2, Supplementary Table S17,
  the historical Tables 1 and 2 with verified fate (five CSVs), and `ranking_tables_2026_10_03.xlsx` (five sheets).
- `build_ed2_2026_10_03.py` makes `ed2_blood_clock_structure_2026_10_03.pdf`, the whole-blood cross-clock figure.
- `build_ed21_2026_10_03.py` makes `ed_verification_census_2026_10_03.pdf`, the verification census.
- `build_figure4_verified_2026_10_14.py` (added 2026-10-14) makes `figure4_verified_2026_10_14.pdf`, the revision's
  Figure 4 (verified intervention effects across independent projects and against GenAge and DrugAge), from
  `verifier/intervention-meta-2026-10-17/` including `concordance/`.
- `build_ed23_benchmark_2026_10_14.py` (added 2026-10-14) makes `ed23_benchmark_2026_10_14.pdf` (ClockEffectBench and
  the agent reruns) from the benchmark result files of the code release (`benchmark/published_clock_2026_10_05/`).
- `build_ed22_replication_2026_10_14.py` (added 2026-10-14) makes `ed22_replication_2026_10_14.pdf`, the replication of
  intervention effects across independent projects, from `verifier/intervention-meta-2026-10-17/` (the primary
  analysis, its three sensitivity runs and `projects.json`; see audit_system/README.md). These inputs are in the
  records deposit under `supporting/intervention-meta/`, not in the 2026-10-05 or 2026-10-10 bundle.

Each also writes a review PNG and `provenance.json`. Use the release's Python environment (NumPy, pandas,
Matplotlib and its Excel-writing dependency). From the code-release root, set `CLOCKBASE_DATA` to the unpacked
bundle's `clockbase_data/` directory, then run:

```sh
export CLOCKBASE_DATA=/path/to/bundle/clockbase_data
export MPLBACKEND=Agg
python -B scripts/figures_2026_10_03/build_figure3_tables_2026_10_03.py --out "$TMPDIR/w21-figure3"
python -B scripts/figures_2026_10_03/build_ed2_2026_10_03.py --out "$TMPDIR/w21-ed2"
python -B scripts/figures_2026_10_03/build_ed21_2026_10_03.py --out "$TMPDIR/w21-ed21"
```

Every `--out` directory must be new and outside the code and data roots. `--final NAME` selects a run under
`audit_system.paths.VERIFIER_ROOT`; when omitted, the reference is read from `BUNDLE.json` in `CLOCKBASE_DATA`
or its parent (the deposited layout). If neither manifest exists, `--final` is required. ED2 uses the fixed
`r6-6-blood/structure-002` and `grimage-age-input-002` analyses; its final argument only validates the selected run.
Existing `audit_system.paths` environment overrides still apply. Choose a writable `TMPDIR` outside the release.

All empirical inputs opened by the three builders are listed below, relative to `clockbase_data/`. `<FINAL>` is
the selected run. All 20 inputs are present in the 2026-10-05 release bundle, whose `BUNDLE.json` names
`final-2026-10-05` as `reference_final`; the five `FIGURE_INPUTS` additions are included. Presence was checked
against `clockbase_data/BUNDLE_MANIFEST.tsv` with a CSV reader using a tab delimiter.

| Builder | Bundle input | 2026-10-05 release bundle |
| --- | --- | --- |
| Figure 3 / tables | `verifier/<FINAL>/comparison-ranking/figure3_category_top15.csv` | Present |
| Figure 3 / tables | `verifier/<FINAL>/comparison-ranking/ranked_comparisons_all.csv` | Present |
| Figure 3 / tables | `verifier/<FINAL>/comparison-ranking/summary.json` | Present |
| Figure 3 / tables | `verifier/<FINAL>/comparison-ranking/table1_pro_aging.csv` | Present |
| Figure 3 / tables | `verifier/<FINAL>/comparison-ranking/table2_anti_aging.csv` | Present |
| Figure 3 / tables | `verifier/<FINAL>/public-export/Table1_pro_aging_v3_status.csv` | Present |
| Figure 3 / tables | `verifier/<FINAL>/public-export/Table2_anti_aging_v3_status.csv` | Present |
| ED2 | `verifier/r6-6-blood/grimage-age-input-002/grimage_age_input.csv` | Present (FIGURE_INPUTS) |
| ED2 | `verifier/r6-6-blood/structure-002/cohort.csv` | Present (FIGURE_INPUTS) |
| ED2 | `verifier/r6-6-blood/structure-002/correlation_repaired_residual_spearman_chronological_age_input.csv` | Present (FIGURE_INPUTS) |
| ED2 | `verifier/r6-6-blood/structure-002/correlation_repaired_residual_spearman_historical_horvath_age_input.csv` | Present (FIGURE_INPUTS) |
| ED2 | `verifier/r6-6-blood/structure-002/summary.json` | Present |
| ED21 | `GEO_data/GEO_hs_dnam/metadata/GSE73950.csv` | Present (FIGURE_INPUTS) |
| ED21 | `verifier/<FINAL>/error-census/census.json` | Present |
| ED21 | `verifier/<FINAL>/published-status/summary.json` | Present |
| ED21 | `verifier/<FINAL>/repairs-dnam/manifest.json` | Present |
| ED21 | `verifier/<FINAL>/repairs-hs_rna/manifest.json` | Present |
| ED21 | `verifier/<FINAL>/repairs-mm_rna/manifest.json` | Present |
| ED21 | `verifier/input-integrity-dnam-001/samples.csv` | Present |
| ED21 | `verifier/replay-dnam-001/inputs/GSE73950.csv` | Present |

ED2's paper-only `inputs/structure_summary.json` is replaced by `verifier/r6-6-blood/structure-002/summary.json`.
Its two paper-only `inputs/correlation_repaired_residual_spearman_*.csv` files are replaced by the corresponding
historical-Horvath-age and chronological-age matrices under that same verifier directory. All three replacements
were verified byte-identical to the evidence snapshots. No paper file is required or moved.

The shared helper probes `references/published_corpus_metrics.json` for optional style diagnostics. That file is
absent in both the evidence builder folders and these released builders; the original documented fallback remains
in use. It supplies style diagnostics, not empirical figure values. Runtime libraries and fonts are environment
inputs rather than bundle data. Preserve the Arial font stack for comparable rendering.

The earlier 2026-10-05 reproduction against the 2026-10-04 run used Python 3.12, Matplotlib 3.10.6, NumPy 1.26.4, pandas 2.2.3, openpyxl 3.1.5,
Pillow 12.3.0 and `--final final-2026-10-04`. All three builders exited 0. Read tracing confirmed exactly 7, 5 and 8
empirical inputs respectively, with no uncovered dependencies. A standalone copy of these 20 inputs (15 from the
existing bundle, five from the added list), without data symlinks, also rebuilt all three with exit 0 using the
parent `BUNDLE.json` default. All three review PNGs and five CSVs were byte-identical between the local-data and
standalone-bundle runs. Standalone-bundle PDFs also matched the local-data rebuilds with zero pixel difference at
150 dpi. All 15 CLI checks rejected existing outputs, repository/data outputs, a missing final without a manifest,
and a traversal final name (exit 2). AST parsing of all four released Python files passed.

All five CSVs are byte-identical to the evidence. All cells, sheet names and dimensions match in the Excel workbook:
Tables 1 and 2 are each 31 by 19 including headers, S17 is 1001 by 19, and both historical tables are 31 by 27.

PDFs were rendered using `pdftoppm -r 150 -png -singlefile` for both evidence and rebuilt files. Pixel differences
were computed on RGB arrays; values below use the 0-255 channel range. Figure 3 needs one column of white padding
on the evidence image's right edge for comparison (top-left alignment, no scaling).

| PDF | Evidence / rebuild raster (width x height) | Changed pixels / total | Mean absolute channel difference |
| --- | --- | --- | --- |
| Figure 3 | 994 x 1439 / 995 x 1439 | 181897 / 1431805 | 14.7041603663 |
| ED2 | 1069 x 879 / 1069 x 879 | 44205 / 939651 | 1.1233511165 |
| ED21 | 1111 x 955 / 1111 x 955 | 42417 / 1061005 | 1.6849662977 |

These PDFs are not pixel-identical. The evidence PDF metadata records Matplotlib 3.9.2, whereas the required
pipeline environment uses 3.10.6; embedded fonts remain Arial. With the plotting code and plotted values preserved,
the changed text spacing, small crop shifts and rasterization are consistent with that renderer-version change.
Figure 3 and ED21 also retain the original dynamic footer and therefore show `final-2026-10-04` instead of the
evidence's `final-2026-10-03`. JSON input hash differences reflect timestamps, source paths and upstream repair
hashes; the plotted counts, fractions, rankings and blood inputs agree. Provenance paths now use portable root
labels, and builder script hashes change with this release adaptation. Plotting code was not altered to force equality.

The saved review PNGs were visually inspected. Figure preflights had no blockers; the preserved figures still report
layout/axis-margin warnings (Figure 3: one warning and one advisory; ED2: three warnings; ED21: one warning).
The 2026-10-05 release bundle contains all 20 inputs and names `final-2026-10-05` as its reference run, which the
builders use by default. Its `final-2026-10-04` directory holds only `public-export/`, so the earlier reproduction
command with that reference cannot be repeated from this bundle. The revision's Figure 3, ED2 and ED21 in figure
package v8 were rebuilt with `--final final-2026-10-05` (`evidence/figures-rebuild-2026-10-05`); Figure 3 and ED21
footers show `final-2026-10-05`. The comparison above remains a dated record of the earlier run.

The 2026-10-10 release bundle also contains all 20 inputs (checked against its `clockbase_data/BUNDLE_MANIFEST.tsv`
in the same way) and names `final-2026-10-10` as its reference run, the run with the design-respecting estimand.
Figure 3, ED2 and ED21 of figure package v9 were rebuilt with `--final final-2026-10-10`
(`evidence/figures-rebuild-2026-10-10`); Figure 3 and ED21 footers show `final-2026-10-10`, and their values differ
from the 2026-10-05 pages because the run differs.
