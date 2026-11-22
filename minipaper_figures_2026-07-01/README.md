# Mini-papers + Fig-5-style figures for FDR<0.05 & score>100 hits (2026-07-01)

Generated a sample-level Fig-5-style figure (and, where one did not exist, a
Claude-SDK-written mini-paper) for every ClockBase hit with FDR<0.05 and
final_score>100, on the reaudited canonical data.

## Coverage (final)

| | count |
|---|---|
| Total FDR<0.05 & score>100 hits | 420 |
| Mapped to a case_id | 412 (8 unmapped: intervention-label mismatch) |
| **Mini-paper present** | **412 / 412 mapped (100%)** |
| **Embedded validated figure** | **316 / 412 (77%)** |

By modality (figure / mini-paper / mapped): mouse RNA-seq 271/356/356 ·
human RNA-seq 23/31/31 · human DNAm 22/25/25.

Three mini-papers were newly written (Opus 4.7 with Paperclip PMC evidence); the rest
already existed and only received the figure.

## What each figure is

LEFT = experimental-design schematic (cohort -> two arms -> assay -> clock ->
biological age). RIGHT = REAL per-sample ClockBase clock predictions for the two
arms with a Welch t-test. No simulated data.

## Data-integrity gate (why 96 hits have a mini-paper but no figure)

For each hit two arm assignments are proposed: an exact-metadata heuristic (maps
the "CASE vs CONTROL" comparison name onto a harmonized metadata field, e.g.
disease=='diabetic retinopathy' vs 'normal') and a Claude-SDK group choice. The
reliability `engine.py` recomputes the statistic; a figure renders ONLY if the
reconstructed arms reproduce the published **direction** AND remain **significant
(p<0.05)**. Otherwise the hit keeps its mini-paper but is flagged (no figure) --
never fabricated. Flagged reasons: complex multi-arm designs, degenerate
comparisons (e.g. "HIV-COPD vs HIV-COPD"), covariate-dependent effects, or arms
with <2 reconstructable samples. DNAm uses the published-significant clock among
the 5 (validation-driven).

## Layout

- `pipeline.py` — the generator (arm derivation, engine validation, figure, prose)
- `hits_score_gt100.csv` — the 420 hits + case_id + coverage flags
- `figures/` — central collection of every rendered figure (png @600dpi + pdf)
- `arms_cache/` — per-hit SDK arm-derivation cache (resumable)
- Deliverables live in each case folder:
  `data/report{,_hs_rna,_hs_dnam}/<case_id>/<case_id>_figure.png|.pdf` and the
  figure is referenced at the bottom of that case's `parsed_report.md`.

## Re-run

```
./venv/bin/python minipaper_figures_2026-07-01/pipeline.py --build-hits
./venv/bin/python minipaper_figures_2026-07-01/pipeline.py --mode full   # resumable
```

## 8 still-unmapped hits (intervention label mismatch, close by hand if needed)

GSE148794 Serpina3n; GSE147608 Inactive SJIA; GSE130306 OIS D6vsD0; GSE135965
Ambisome; GSE97158 post-malaria recovery; GSE163356 HIIT; GSE129748 nasal allergen
challenge; GSE140344 Down syndrome. Each has FDR<0.05 & score>100 but the combined
CSV intervention label does not exactly match `scores_summary`, so no case folder
was resolved. (Several, e.g. Ambisome, were deliberately dropped from the
manuscript in earlier audits.)
