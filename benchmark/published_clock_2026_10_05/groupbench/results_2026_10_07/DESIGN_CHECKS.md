# Design checks for repeated subjects, pairing and age imbalance (2026-10-07)

Follow-up to `../results_2026_10_06/REPORT.md`. Protocol: `../PROTOCOL.md`, addendum of 2026-10-07. Code:
`../design_checks.py`, `../design_eval.py`, `../design_validate.py`. Data:
`/Entropy/kying/clockbase-benchmark-2026-10-05/groupbench/design-checks/` (`v3` checks, `validation-v2` and
`headline-v2` reads). The checks are a prototype outside the released verifier; nothing in the release changed.

## Question

The group-choice benchmark found three kinds of error in released human DNA methylation rows that the verifier does
not test for: one person contributing several samples to an arm, paired arms tested as independent, and arms of
different age compared on a clock that tracks age. This report asks how much of the released atlas they touch, whether
checks on recorded fields catch them, and whether the ranked tables are affected.

## Method

1. **Checks on recorded fields** (deterministic). Subjects from a GEO sample field named as a subject identifier; ages
   from the harmonized age, else for human series a GEO age-at-sampling field. Findings on the released arms:
   REPEATED_SUBJECT, PAIRED_AS_INDEPENDENT (at least half of the subjects in both arms), AGE_IMBALANCE (human series;
   median gap above 5 years or standardized difference above 0.5; pace clocks and declared age or time contrasts
   exempt).
2. **Recomputation.** Each flagged row is recomputed on subject means, with a paired t-test for paired arms and age as
   a covariate under age imbalance, and its call (direction, two-sided P < 0.05) is compared with the released call.
   The released P value was reproduced for every flagged row before recomputing (0 of 5,343 not reproduced).
3. **Two model readers** (openai:gpt-6-astra, openai:gpt-6-sol) read the GEO record of a comparison (series text and
   the fields of the arms' samples; no check output, no clock value) and state which samples come from the same
   individual, each sample's age, and whether age or time is the contrast. Findings are derived from each answer with
   the same rule functions; a finding counts when both readers give it. On the 464 comparisons of the validation and
   benchmark samples the two readers gave the same finding in 461 (repeated subject), 459 (pairing) and 462 (age).
4. **Samples read:** 285 released comparisons from series outside the benchmark (per modality 20 flagged for subjects,
   15 flagged for age, 60 unflagged; one per series; seed 20261007), the 179 released benchmark comparisons, and all
   257 comparisons of Table 1, Table 2, Figure 3 and Table S14 whose arms hold at most 120 samples.

The first version of the checks (v2) also used the harmonized subject field and applied the age check to mouse series.
The readers confirmed 0 of 2 sampled flags from the harmonized subject field and 9 of 15 from mouse ages (five wrong
harmonized ages, one developmental contrast), so the reported version (v3) leaves both out. The rule changes made after
seeing results are listed in the protocol addendum.

## Results

### 1. Released atlas, checks on recorded fields

| | Human DNAm | Human RNA-seq | Mouse RNA-seq | All |
|---|---|---|---|---|
| Released rows | 13,022 | 24,379 | 27,081 | 64,482 |
| Rows with a recorded subject field | 1,419 (10.9%) | 2,682 (11.0%) | 635 (2.3%) | 4,736 (7.3%) |
| Rows with recorded age | 4,216 (32.4%) | 1,610 (6.6%) | not used | |
| **Flagged rows** | **2,337 (17.9%)** | **2,705 (11.1%)** | **301 (1.1%)** | **5,343 (8.3%)** |
| Repeated subject | 427 | 730 | 183 | 1,340 |
| Paired arms | 951 | 1,912 | 164 | 3,027 |
| Age imbalance | 1,233 | 453 | not applied | 1,686 |
| Call changes on recomputation | 463 | 331 | 36 | 830 (1.3% of released) |
| Not testable on recomputation | 134 | 224 | 32 | 390 |

Of the 19,103 released rows with P < 0.05, 1,418 are flagged; 407 (2.1% of the significant rows) lose significance or
reverse direction on recomputation, and 131 cannot be tested (one subject in an arm, fewer than three pairs, or age
separating the arms).

The three findings act differently. Paired arms tested as independent mostly lose power: of 3,027 rows, 2,517 keep
their call, 309 become significant and 81 lose significance. Repeated subjects inflate significance: of 1,340 rows,
147 lose significance, 46 gain it, and 265 cannot be tested on subject means (151 have a single subject in an arm).
Age imbalance moves calls both ways: of 1,686 rows, 266 lose significance, 105 gain it and 4 reverse.

### 2. The verifier's recomputation removes the agent's own pairing and age adjustment

Every released row is recomputed as a Welch two-sample test. The agent's own test was a paired test or a mixed model
for 1,264 released rows and a model with age as a covariate for 1,994. Where the checks flag such a row and can
recompute it:

| Agent's own analysis | Rows | Released call = agent's call | Design-respecting call = agent's call |
|---|---|---|---|
| Model with age, arms age-imbalanced | 603 | 435 (72%) | 562 (93%) |
| Paired test or mixed model, arms paired on a recorded subject field | 498 | 427 (86%) | 484 (97%) |

In 104 of the 603 rows the agent's age-adjusted result was not significant and the released row is; in 62 of the 498
the agent's paired result was significant and the released row is not. Part of the released human-cohort errors
therefore comes from the verifier's uniform estimand, where the agent had analysed the design correctly.

Correction to the 2026-10-06 report: its released-atlas estimate for human DNAm (33 of 44 units) counted only error
types confirmed for the agent's own analysis. Two more released units (GSE62867, tissues of the same six patients)
were tested as paired by the agent and are released unpaired; counting them, the estimate is 35 of 44 (80%;
cluster-adjusted 95% CI 46-97%). The mouse estimate (5 of 89) is unchanged.

### 3. The benchmark's released errors (in-sample)

| Human DNAm (44 released units) | Error units flagged | Fully addressed | Errors left if flagged rows are recomputed |
|---|---|---|---|
| As released, with the correction above | | | 35 / 44 (95% CI 46-97%) |
| Checks on recorded fields | 29 of 35 | 22 | 13 / 44 (95% CI 8-62%) |
| Recorded fields plus both readers | 35 of 35 | 27 | 8 / 44 (95% CI 3-49%) |

Neither layer flags any of the 9 correct units. The 8 errors left are of other kinds: arms not carrying the named
levels (3), another factor separating the arms (4), a mixed arm (1). Recorded fields miss the repeated subjects of
GSE43976 and GSE118468 and the pairing of GSE62867, whose subject identifiers exist only in sample titles or in a
misspelt field; the readers find them. (Against the 33 errors of the 2026-10-06 report, the same rows read 29 of 33
and 11 / 44, then 33 of 33 and 8 / 44.)

Mouse (89 released units, 5 errors): the recorded fields flag none; the readers find the two repeated-animal units
(GSE101657). Three errors remain (two inversions, one pairing the readers did not find).

The checks were designed from these errors, so this section shows coverage, not validation.

### 4. Validation on series outside the benchmark

Flags confirmed by both readers (comparisons where both readers could decide):

| | Human DNAm | Human RNA-seq | Mouse RNA-seq | All |
|---|---|---|---|---|
| Repeated subject | 3 / 4 | 6 / 6 | 10 / 10 | 19 / 20 |
| Paired arms | 13 / 14 | 14 / 15 | 11 / 12 | 38 / 41 |
| Age imbalance | 17 / 17 | 15 / 16 | not applied | 32 / 33 |

89 of 94 flags are confirmed (95%). The five others: identifiers numbered within each group, read as shared subjects
(2); donor-derived cultures the readers classed as cell lines (2); a harmonized age absent from the GEO record (1).

Unflagged comparisons in which both readers find shared individuals:

| | Human DNAm | Human RNA-seq | Mouse RNA-seq |
|---|---|---|---|
| Both readers could decide | 50 of 60 | 51 of 60 | 32 of 60 |
| Repeated subject or pairing found | 12 (24%; 95% CI 13-38%) | 8 (16%; 7-29%) | 5 (16%; 5-33%) |
| Age imbalance found (of 60) | 0 | 0 | 0 |

The recorded fields are precise and incomplete. In most series the individual is identified only in sample titles or
the series text (patient codes at two visits, tumour and adjacent tissue, one donor's cells under two treatments,
several tissues of one mouse, one library on several sequencing lanes). Of the 25 comparisons found, 17 are paired
designs and 10 contain repeated subjects (2 both). No unflagged comparison had an age imbalance the readers could see.

### 5. Ranked tables

Checks on recorded fields: none of the 30 Table 1 rows, 30 Table 2 rows or 60 Figure 3 rows is flagged; 30 of the
1,000 ranked comparisons are (20 repeated subject, 5 paired, 9 age imbalance), of which 9 lose significance on
recomputation and 1 cannot be tested.

Readers (257 comparisons read; 14 Table S14 comparisons above 120 samples were not read):

| Table | Read | Both readers could decide on individuals | Repeated subject (both readers / one) | Paired | Age imbalance |
|---|---|---|---|---|---|
| Table 1 | 30 | 14 | 0 / 0 | 0 | 0 |
| Table 2 | 30 | 19 | 2 / 1 | 0 | 0 |
| Figure 3 | 60 | 33 | 1 / 0 | 0 | 0 |
| Table S14 | 201 | 124 | 16 / 4 | 2 | 3 |

Table 2 entries with a finding:

- **GSE133505, Myf6 knockout vs wild type (10 months):** the 8 vs 8 samples are four technical replicates of two mice
  per genotype. On mouse means (2 vs 2) the direction holds (P = 0.0046; released P = 7.5e-08). Two mice per arm is
  below the ranking's minimum of three.
- **GSE162613, miR-155KI vs WT (also Figure 3):** the 6 vs 6 samples are two cell populations from three mice per
  genotype. On mouse means (3 vs 3) the result holds (P = 1.8e-05; released P = 2.7e-10).
- **GSE151160, CD36 KO vs wild type (one reader):** the titles show one library per genotype on four sequencing lanes
  (`12_S28_L001` to `L004`, `8_S25_L001` to `L004`), which is 1 vs 1; the other reader could not tell.

Table S14: of the 16 comparisons with repeated subjects by both readers, 8 keep their call on subject means and 8 lose
significance (GSE191200, GSE154816, GSE159014 heart and liver, GSE178768, GSE79408, GSE79661 two contrasts). The two
paired comparisons keep their call. The three age-imbalanced comparisons (GSE113018, GSE50759, GSE191200) lose
significance with age as a covariate. Per-comparison values: `headline.json`, `ranked.json`.

### 6. Why the human DNAm rate is high

`Figure_benchmark_2026-10-07.pdf` (built by `../figure_benchmark_2026_10_07.py`; values in
`figure_benchmark_data.json`) shows the benchmark in seven panels; its intervals are Clopper-Pearson on the
design-effect-adjusted sample size for every row.

- **The sample is a hard case.** The annotated human series are disease cohorts, many longitudinal. Four series
  (GSE185389, GSE49904, GSE130491, GSE151355) hold 24 of the 44 comparisons and 23 of the 35 errors; 11 of the 17
  series have an error.
- **Most errors are of three kinds of design statistics.** 27 of the 35 error units carry only age imbalance (15
  units), pairing (12) or repeated individuals (6); 8 carry another kind.
- **An error here is a design flaw, whether or not the conclusion changes.** Of the 29 error units the recorded fields
  flag, the design-respecting test changes a call in 15, changes none in 12, and cannot be computed for 2.
- **Three mechanisms.** The clocks are compared as raw clock age, so any age difference between arms enters the
  result. Individuals are identified in sample titles or series text, where neither the agent nor the verifier reads
  them systematically. And the verifier recomputes every row as an unpaired Welch test, which removes pairing and age
  adjustment where the agent had them (section 2).
- **It is not the atlas-wide rate.** The same checks flag 18% of released human DNAm rows, 11% of human RNA-seq and
  1% of mouse RNA-seq, and 2.1% of significant released rows lose their call.

### 7. The age instruction

The agent's task prompt (`prompts/task_intention_human_dnam.txt`) says: "DO NOT compare young and old populations
directly", "Determine relevant covariates (age, sex, tissue type, batch effects)", "Focus on comparing treatment vs.
control within the same age group", "Use appropriate tests based on data structure (t-test, ANOVA, or linear models)",
and asks for a `covariate` column, "Variables controlled for in analysis". It names age as a covariate to consider; it
does not require age as a term of the model.

| Released human rows | DNAm | RNA-seq |
|---|---|---|
| All | 13,022 | 24,379 |
| Age recorded for both arms | 4,216 | 1,610 |
| of these: age is a term of the agent's model | 1,498 (36%) | 199 (12%) |
| of these: age declared in the `covariate` column, test without it | 1,603 (38%) | 916 (57%) |
| of these: age neither declared nor modelled | 1,115 (26%) | 495 (31%) |
| Released statistic carries age | 0 | 0 |

Three things went wrong. The prompt left the choice of test open, and the agent ran a two-sample test in most rows.
In 1,603 DNAm and 916 RNA-seq rows it wrote age into the `covariate` column of a result whose test did not contain
age, so the table says more than the analysis did. And where the agent did model age, the verifier's uniform Welch
recomputation removed it: over the 1,142 such DNAm rows with a stored agent call, the released call equals the agent's
in 922 (81%), and a recomputation with age equals it in 1,092 (96%).

It can be corrected without re-running the agent, for the rows that record age. With clock ~ arm + age as the
estimand for every released human row with recorded ages (age contrasts and the pace clock left as they are):

| | DNAm | RNA-seq |
|---|---|---|
| Rows recomputed | 3,088 | 1,451 |
| Call unchanged | 2,258 | 1,072 |
| Significance lost | 278 | 39 |
| Significance gained | 167 | 33 |
| Direction reversed | 9 | 0 |
| Age term not estimable (one age for all samples) | 310 | 252 |
| Age separates the arms (cannot be adjusted) | 66 | 55 |

It cannot be corrected where GEO records no age: 8,806 of the 13,022 released DNAm rows and 22,769 of the 24,379
RNA-seq rows. Those rows can only be labelled as unadjusted. An earlier version of this report counted 2,906 released
rows with age in the agent's model; that count matched the letters "age" anywhere in the formula, including clock
names such as PhenoAge and terms such as agent and passage. Read term by term, it is 1,994 (`design_checks.formula_has_age`).

### 8. The candidate estimand against published results (2026-10-09)

The candidate release estimand was written into the protocol (addendum item 10) and committed before scoring: subject
means and a paired test where a subject field shows repeats or pairing; least squares on arm and age for human rows
with recorded ages on a clock that tracks age; the released Welch test otherwise. One variant, scored with the
existing metrics (`design_eval.py --external`, `external.json`).

| Agreement with published results | Raw agent | Released | Candidate |
|---|---|---|---|
| Three-way call, 105 published results on released rows | 79% | 71% | 80% |
| Direction, 35 published significant results | 86% | 89% | 94% |
| Balanced accuracy of the call | 0.67 | 0.60 | 0.67 |
| Three-way call, 55 same-sample contrasts | 85.5% | 80.0% | 85.5% |
| Cohen's kappa, same contrasts | 0.74 | 0.64 | 0.74 |

On all 116 published results with a released row, call agreement rises from 68.1% to 80.9%: 12.8 points (95% CI 5.2
to 21.3, bootstrap over 43 papers). Of the 15 calls that change, 13 move to the paper's call, 1 moves away and 1
disagrees both ways; in the same-sample contrasts 4 of 5 move to it and 1 away. The verifier's Welch recomputation had
lowered agreement with published results below the agent's own; the candidate restores it while keeping the
verifier's repairs.

In the released atlas the candidate differs from the Welch test in 7,215 of 64,482 rows (11.2%): 458 of the 19,103
significant rows lose or reverse their call (2.4%), 494 rows become significant, and 398 have no candidate call (one
subject in an arm, fewer than three pairs, or age separating the arms). `Table_benchmark_2026-10-09.md` gives the
benchmark by stage and `Figure_benchmark_stages_2026-10-09.pdf` draws it.

### 9. Conclusion-level accuracy against expert-defined groups (2026-10-09, post hoc)

The benchmark above counts a comparison as an error when its design has any flaw. This section asks the narrower
question of protocol amendment 3: does the comparison reach the conclusion (direction, P < 0.05) that the complete
expert-defined groups give? It is scored at three stages, against the declared groups analysed with the
design-respecting estimand (`a_result.py --stages`, `conclusion_stages.json`; protocol addendum item 11).

| | Raw agent | Released | Candidate |
|---|---|---|---|
| Mouse RNA-seq (29 series) | 93 / 98 (94.9%) | 88 / 89 (98.9%) | 88 / 89 (98.9%) |
| Human DNAm, labels confirmed (14 series) | | 114 / 130 (87.7%) | 128 / 130 (98.5%) |
| Human DNAm, all released (16 series) | 147 / 192 (76.6%) | 144 / 180 (80.0%) | 166 / 176 (94.3%) |

- No released or candidate row reaches a significant conclusion opposite to the expert-group analysis; the raw agent
  does in 1 mouse and 10 human DNAm rows.
- "Labels confirmed" leaves out released rows the verifier marked LABEL_UNCONFIRMED (2,571 of 64,482 released rows;
  1,280 of the 13,022 human DNAm rows). Eight of the ten remaining human DNAm disagreements are in such rows
  (GSE118468, GSE32148). The other stratifications examined: VERIFIED rows, candidate 80 / 82; rows without an agreed
  CASE_NOT_CONDITION issue, candidate 129 / 129.
- Against the Welch reference of amendment 3 the released rows agree in 88 / 89 (mouse), 133 / 134 (human DNAm,
  labels confirmed) and 166 / 184 (human DNAm, all).
- Intervals are wide because results cluster in series: design-effect-adjusted lower bounds are 87% (mouse) and 83%
  (human DNAm, labels confirmed, candidate); unadjusted 94% and 95%. The point estimates exceed 95%; a 95% lower bound
  above 95% would need more annotated series.
- Published results (section 8): the candidate contradicts a published significant result in 1 of 94 calls and 1 of 55
  same-sample contrasts, and has the same direction in 25 of 26 and 21 of 22 results significant on both sides.
  Three-way call agreement stays at 80-86% because the paper is significant where ClockBase is not in 14 of the 18
  disagreeing calls.

### 10. The benchmark on significant findings, by modality (2026-10-09)

At K.A.Y.'s request the benchmark is also summarised on findings that are significant in the reference, with a
threshold-free direction AUROC, separately for DNA methylation and RNA-seq clocks (protocol addendum item 12, defined
before computing; `design_eval.py --auroc`, `auroc.json`, `Table_benchmark_by_modality_2026-10-09.md`).

| | Direction AUROC | Same direction, both significant |
|---|---|---|
| DNAm, published results, published groups recomputed | 0.97 (0.91-1.00) | 33 / 34 |
| DNAm, published results, released row / candidate | 0.92 (0.79-0.99) / 0.98 (0.91-1.00) | 19 / 21 / 25 / 26 |
| DNAm, expert-defined groups, raw agent / released / candidate | 0.71 / 1.00 / 1.00 | 24 / 34 / 28 / 28 / 35 / 35 |
| Mouse RNA-seq, expert-defined groups, raw agent / released | 0.96 / 1.00 | 35 / 36 / 34 / 34 |

The metric is sound for the question "when the reference finds an effect, does ClockBase find the same direction".
It does not see findings that only ClockBase calls significant (6 of 59 against published results with the published
groups; 27 of 143 released human DNAm results against the expert-defined groups, 8 of 139 with the candidate) or
published effects ClockBase does not reach (33 of 54 replicated), which the table reports beside it. Against the
expert-defined groups the samples and clock values are ClockBase's own, so 1.00 there means the right groups in the
right orientation, not an independent replication. Published results exist for DNA methylation only; mouse RNA-seq
rests on the expert-defined groups and three published contrasts; human RNA-seq has no reference.

### 11. Direction, significance and effect size against the papers (2026-10-09)

The benchmark's question, as K.A.Y. put it: when ClockBase compares an intervention with its control, does it find
the effect the paper reports - same direction, same significance, similar size? `Table_benchmark_by_modality_2026-10-09.md`
answers in those three terms for each modality (protocol addendum items 12 and 13).

- **Direction** agrees: against published DNA methylation results the direction AUROC is 0.97 for the paper's groups
  and 0.98 for the candidate on atlas rows (0.92 as released); the same direction in 33 of 34 and 25 of 26 results
  significant on both sides.
- **Significance** agrees less: of 54 published significant results, ClockBase finds 33 significant in the same
  direction (61%); of 59 not significant, 53 are not significant in ClockBase (90%); overall 87 of 113 (77%). On atlas
  rows the candidate reaches 77 of 94 (82%), the released rows 66 of 94 (70%).
- **Effect size** can be compared for 30 results from 7 papers that report an effect in years: Spearman 0.61
  (0.27-0.93), same sign in 25, median absolute difference 1.5 years against a median published effect of 1.9 years.
- **Overflow values.** Two series in this comparison carry clock values of the order of 1e39 years; without a guard
  the Spearman correlation is 0.33. In the released atlas, 46 human DNAm rows have an arm mean above 1,000 years and a
  further 270 an arm mean above 200 years on a clock in years. The released verifier has no physical-range check.
- Mouse RNA-seq: 3 of 3 published contrasts on the same samples agree in direction and call; against the
  expert-defined groups the released results agree in 88 of 89.

## What this supports

- Checks on recorded subject and age fields can be added to the verifier as they are: 95% of their flags were
  confirmed on series outside the benchmark, and they flag 8.3% of released rows (18% of human DNAm, 1% of mouse).
- The larger share of non-independent samples is not in a subject field. Readers found shared individuals in 16-24% of
  unflagged comparisons where they could decide. Covering these needs the two-reader step, which agreed 99% of the
  time and found every repeated-subject error of the human benchmark.
- Recomputing with the agent's own design where its test was paired or age-adjusted would restore the agent's call in
  most of those rows.
- In the human DNAm benchmark, the three checks with both layers would leave 8 of 44 released units in error (from 33
  as reported, 35 with the correction), all of other kinds.
- Table 1 and Figure 3 are unaffected apart from GSE162613, whose result holds. Table 2 has two entries (GSE133505,
  GSE151160) that rest on fewer than three animals per arm.

## Limits

- The reference standard is two language models. Their answers were inspected for the Table 2 entries and for the
  unconfirmed flags and missed comparisons listed above; no person outside this analysis has reviewed them.
- A comparison where a reader could not decide is counted as unresolved, most often in mouse series that do not
  identify animals (28 of 60 unflagged mouse comparisons). Unshared individuals cannot be shown where GEO does not
  record them.
- The validation samples are small (intervals above) and exclude comparisons with more than 120 samples in the arms.
- Losing significance on subject means partly reflects fewer units; it marks a result the design does not support at
  P < 0.05, not a result shown to be false.
- The readers' individuals were used for recomputation only in the ranked tables. Atlas-wide counts of affected rows
  beyond the recorded fields were not computed.

## Decisions for the authors

1. Adopt the recorded-field checks and the agent-design rule in the verifier. This changes released statistics, needs
   a new final run, and changes numbers in the manuscript and the release.
2. Add the two-reader step for individuals: for the 1,000 ranked comparisons, or for all released comparisons (about
   54,400 distinct arm pairs; the 721 comparisons read here used 2.25 million prompt and 0.82 million completion
   tokens).
3. Table 2: replace or annotate GSE133505 and GSE151160; keep GSE162613 with the statistic on mouse means.
4. Table S14: the eight comparisons that lose significance on subject means and the three age-imbalanced ones.
