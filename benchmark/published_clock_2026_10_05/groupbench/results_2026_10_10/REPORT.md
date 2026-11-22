# ClockEffectBench: verifier run final-2026-10-10 (2026-10-10)

The verifier now tests each row with a Welch t-test or, where the sample records call for it, with means per
individual, a paired t-test or age as a covariate; it drops clock outputs above 1,000 in absolute value and reads
ranking candidates for shared individuals with two models. This folder scores that run with the benchmark of
`../PROTOCOL.md` (addendum item 14). All values are computed from the run's public export and the benchmark data;
nothing is simulated. The design checks that led to the change, and their validation, are in
`../results_2026_10_07/DESIGN_CHECKS.md`.

| File | Content |
|---|---|
| `Table_benchmark_by_modality_2026-10-10.md` | direction, significance and effect size by modality and stage |
| `Table_benchmark_2026-10-10.md` | the same benchmark by stage, with intervals |
| `Table_benchmark_2026-10-10.csv` | one row per reference, modality, measure and stage (Supplementary Table S18) |
| `Figure_benchmark_stages_2026-10-10.pdf` | the stages drawn |
| `auroc.json`, `external.json`, `conclusion_stages.json`, `in_sample_corrected.json`, `in_sample_with_readers_corrected.json` | the scored results |
| `published_released_test.json`, `published_released_test.csv` | the paper's groups with the Welch test and with the released test, one row per published result (protocol item 19) |
| `export_as_checks.json` | rows compared between final-2026-10-05 and final-2026-10-10, by test |
| `ground_truth_counts.json` | size of the published ground truth |

Stages: **raw agent** (the stored result), **verifier, Welch test only** (final-2026-10-05), **released**
(final-2026-10-10).

## Results

- Published DNA methylation clock results, the paper's groups with ClockBase clock outputs and the released test
  (110 results, 32 papers; protocol item 19): direction AUROC 1.00
  (0.98-1.00); of 52 published significant results,
  51 have the published direction, 40 are significant in that direction and none is significant in the
  opposite direction; 51 of 58 published null results are null. With a Welch test for every
  result (the first scoring): AUROC 0.97, 33 of 54 significant, 33 of 34 in the same direction when both
  are significant. Reading individuals from the sample records with two models where the series has no individual
  field (13 comparisons with agreed shared individuals) leaves 40 of 52. Effects in
  years (30 results, 7 papers): Spearman 0.61 (0.27-0.93), 0.33 before four results with clock outputs above 1,000
  years were recomputed without them.
- Against language models without the data, on the 110 results every system answers (36 published increases,
  16 decreases): ClockBase identifies 16 of the decreases and 35 of the increases; the five
  models identify 6-6 decreases and 20-31 increases. Balanced accuracy of the three-way call:
  0.78 for ClockBase, 0.46-0.52 for the models, 0.33 for a constant rule.
- Atlas rows against the same published results: direction AUROC 0.95 raw agent, 0.92 Welch only, 0.98 released;
  same direction 21/22, 19/21, 25/26; three-way call agreement on the 93 scored calls of the 115 results with a row
  at both verifier stages 68.8% Welch only and 81.7% released (difference 12.9 points, 95% CI 5.3 to 21.6, bootstrap
  over 43 papers; 13 of the 15 changed calls move to the paper's call and 1 away).
- Conclusions against expert-defined groups: mouse RNA-seq 93/98 raw agent, 88/89 Welch only, 88/89 released; human
  DNAm, rows with confirmed orientation, 114/130 Welch only and 126/128 released; all released human DNAm rows 144/180
  and 164/178. No row of either verifier stage reaches a significant conclusion opposite to the expert-group analysis.
- Strict count (any design flaw): mouse 20/98 raw agent, 5/89 Welch only, 5/89 released; human DNAm 39/44, 34/43,
  13/43; with the two readers applied to every comparison 3/89 and 8/43.
- In the atlas, the released statistic differs from the Welch test in 7,157 of 64,482 rows (11.1%); 456 of 19,103
  significant rows lose or reverse their call, 493 rows become significant and 429 are no longer released.

## Why the Welch test recovered fewer published significant results

S3 tested the paper's groups with an unpaired Welch test. Of its 21 unrecovered published significant results, 19
were in the published direction, and the papers behind most of them had used paired tests, mixed models with an
individual term, or regressions with age. The released test pairs or averages by individual where a GEO field
records individuals and enters age where it is recorded: it recovers 8 more and loses 1, and the one result that was
significant in the opposite direction under the Welch test (GSE52588, cases compared with older relatives) reverses
to the published direction. No rule by sample size, clock or published P value separates the remaining results
(with at least ten samples per group the released test recovers 21 of 28), so none is excluded.

Published significant results the released test does not recover (with read individuals, S3r):

| Series | Clock | Published | Paper's test | ClockBase |
|---|---|---|---|---|
| GSE25395 | Horvathv1 | higher (P = 0.0086) | multivariate linear regression of DNAm age on chronological age and DS | not significant, same direction; P = 0.06; n = 29 and 21; welch |
| GSE59457 | Horvathv1 | higher (P = 0.036) | Kruskal-Wallis (nonparametric) group comparison test | not significant, same direction; P = 0.08; n = 6 and 4; age adjusted |
| GSE75008 | Horvathv1 | lower (P = not given) | test not stated | not significant, same direction; P = 0.29; n = 40 and 40; paired |
| GSE142439 | Horvathv2 | lower (P = 0.042) | one-sided linear mixed-effects model (lmer, lmerTest) | not significant, same direction; P = 0.08; n = 8 and 8; paired |
| GSE100249 | PhenoAge | higher (P = 0.005) | linear regression | not testable; n = 1 and 3; welch |
| GSE100249 | Horvathv2 | higher (P = 0.048) | linear regression | not testable; n = 2 and 2; welch |
| GSE45727 | Horvathv1 | higher (P = 0.013) | regression coefficient t-test | not significant, opposite direction; P = 0.47; n = 3 and 6; welch |
| GSE162560 | Horvathv2 | lower (P = 0.003) | paired t-test | not significant, same direction; P = 0.46; n = 32 and 32; paired |
| GSE193836 | Horvathv1 | higher (P = 0.036) | test not stated | not significant, same direction; P = 0.83; n = 7 and 8; welch |
| GSE185973 | Horvathv2 | lower (P = not given) | One-way analysis of variance followed by Dunnett’s test | not significant, same direction; P = 0.09; n = 6 and 3; welch |
| GSE210245 | GrimAgeV1 | lower (P = 0.0093) | Two-sided paired Student's t-test, uncorrected | not significant, same direction; P = 0.40; n = 18 and 18; paired |
| GSE106648 | Horvathv1 | lower (P = 0.0001) | MANCOVA, clinical group term | not significant, same direction; P = 0.10; n = 140 and 139; age adjusted |
| GSE53740 | Hannum | higher (P = 0.002) | Pairwise Wilcoxon test with Benjamini–Hochberg multiple testing correc | not significant, same direction; P = 0.33; n = 127 and 193; age adjusted |
| GSE53740 | Hannum | higher (P = 0.03) | Pairwise Wilcoxon test with Benjamini–Hochberg multiple testing correc | not significant, same direction; P = 0.33; n = 127 and 193; age adjusted |

1 of these papers used a one-sided test and matches ClockBase's one-sided P value. Of the rest, the papers adjusted for
covariates the deposits do not hold, used sample sets that differ from the deposit, or compared groups of fewer than
ten samples; in three paired designs the readers found the pairs and the paired test still shows a smaller
difference in ClockBase clock outputs than the paper reports.

## Limits

- The estimand was designed from the expert-annotated benchmark, so the expert-group rows are in-sample; the
  published results were scored after the estimand was fixed.
- The conclusion-level comparison and its restriction to rows with confirmed orientation are post hoc (item 11).
- Results cluster within series and papers; design-effect-adjusted lower bounds are 87% (mouse) and 83% (human DNAm,
  confirmed orientation).
- The individuals audit read 1,092 ranking candidates; at least one reader could not decide in 456 (431 mouse). No
  benchmark comparison changed through the audit.
- Human RNA-seq has no reference and is not benchmarked.
