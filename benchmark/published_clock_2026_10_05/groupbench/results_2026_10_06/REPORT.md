# Group-choice and result benchmark of the ClockBase agent (2026-10-06)

Protocol: `../PROTOCOL.md` (v2 with amendments 1-4; v2 committed before any classification, inputs frozen with hashes
in `../PREREGISTRATION.json`). Data: `/Entropy/kying/clockbase-benchmark-2026-10-05/groupbench/`.

## What was tested

For each comparison the agent made (its own case and control samples, recovered from its saved code or as the unique
split reproducing its stored P value), whether the arms are the groups the comparison names, whether anything else
separates them, whether the statistics respect the design (repeated subjects, pairing, arm size), and whether the
reported direction agrees with the named groups. Ground truth for the groups is expert per-sample annotation made for
other studies:

- **A1, mouse RNA-seq:** Tyshkovskiy et al., *Nature* 654:173-188 (2026), rodent meta-dataset annotation (Zenodo
  10.5281/zenodo.18763485); 1,253 samples in 42 atlas series mapped one-to-one to GSMs; 98 agent units in 32 series.
- **A2, human DNA methylation:** ComputAgeBench per-sample labels; 47 agent units in 18 series (44 assessable).

Each unit was classified by fixed rules and then reviewed by three independent model reviewers (design, metadata,
skeptic briefs) who saw names, arms, per-sample labels and GEO metadata, but no rule outcome and no verifier output
(Fleiss' kappa 0.81). A rule error stood unless all three reviewers gave the same checkable reason; a rule-correct unit
became an error when at least two reviewers found one. Before scoring, the rules detected six kinds of SIMULATED
planted errors in 100% of applicable cases (`SIMULATED_detection.json`).

**E1** is the raw agent (its arms and reported direction). **E2** is the released atlas (post-repair arms and the
verifier's recomputed direction, rows released as VERIFIED or PROVISIONAL).

## Results

### A1, mouse RNA-seq (98 units, 32 series)

| Estimand | Errors / units | Cluster-adjusted 95% CI | Series with an error |
|---|---|---|---|
| E1 raw agent, census as run | 20 / 98 (20%) | 17-36% (model; marginal 33%) | 11 / 32 |
| E1, adjudication-material artifact reclassified | 14 / 98 (14%) | 3-38% | 8 / 32 |
| E1 design errors only (no statistical) | 13 / 98 (13%) | 3-34% | 7 / 32 |
| E1 any error, counting effect-sign record errors | 31 / 98 (32%) | 12-58% | 13 / 32 |
| **E2 released atlas** | **5 / 89 (5.6%)** | **0.3-23%** | **3 / 29** |

Error types (E1): P value from an arm under 3 samples 8; reported direction contradicting the named groups 6; arms
not carrying the named levels 5; one animal contributing liver and kidney to an arm 2; paired cells tested as
independent 1. Separately, 17 units carry an effect size whose sign contradicts the agent's own labels while its
direction field is right (GSE190939 and five other series); the verifier swapped these arms.

The verifier withheld 8 and repaired 7 of the 20 E1 errors and released 5 (two repeated-subject and one pairing error,
two inversions); it withheld 1 of 78 correct units. Errors by contrast type: age contrasts 7 / 12, cell-type contrasts
3 / 5, interventions 9 / 54, genotypes 1 / 27. Replayed arms 0 / 18 errors, recovery-only arms 12 / 46.

### A2, human DNA methylation (44 assessable units, 17 series)

| Estimand | Errors / units | Cluster-adjusted 95% CI | Series with an error |
|---|---|---|---|
| E1 raw agent, census as run | 39 / 44 (89%) | 69-94% | 14 / 17 |
| E1, age confounding only when median ages differ > 5 years | 31 / 44 (70%) | 43-90% | 13 / 17 |
| E1 design errors only | 29 / 44 (66%) | 51-80% | 13 / 17 |
| **E2 released atlas** | **33 / 44 (75%)** | **41-95%** | **10 / 17** |

Error types (E1): raw clock age compared between arms of different age 15; repeated measures of the same people
analysed as independent 10; one person contributing several samples to an arm 6; reported direction contradicting the
named groups 7; arms not carrying the named levels 4; another factor separating the arms 4; impure arm 1. Ten units
come from one longitudinal HIV study (GSE185389), four each from GSE130491 (DMF, same patients before and after) and
GSE151355. The verifier repaired 6 of the 39 errors and released 33: it does not check pairing, repeated subjects or
age imbalance between arms.

### B, agreement with published results for the same samples

56 published contrasts (26 papers) whose groups were assigned to samples and contained in an agent comparison on the
same clock: three-way call agreement 84% (paper-clustered 95% CI 75-93%), Cohen's kappa 0.71 (0.55-0.88), same
direction in 20 of 21 contrasts significant on both sides. Of the 9 disagreements, three reviewers attributed 3 to the
agent (5.4% of matched contrasts): an unpaired test of repeated visits in GSE176394 (two contrasts, whose arms were also
built one sample off), and raw clock age compared between young and old patients in GSE100850. The other 6: age
adjustment or sample subsets the GEO metadata do not record (3), clock or preprocessing differences (2), power (1).
From Tyshkovskiy et al. 2026 Fig. 5a, three mouse disease contrasts were matched and all three agree.

## What this supports

- The released mouse atlas, in expert-annotated series: 5.6% of comparisons with an error (95% CI 0.3-23%). The
  pre-stated threshold for calling it low (upper bound below 10%) is not met with 89 units in 29 series.
- The raw agent makes errors at a rate the verifier must catch: 14-20% in mouse (32% counting effect-sign record
  errors); the verifier caught 15 of 20.
- Human DNA methylation comparisons in these disease cohorts have design-statistics errors the verifier does not test
  for. Of the 33 released errors, 25 involve only repeated subjects, paired samples or age imbalance under raw clock
  age, 6 involve one of these together with another error (another separating factor, label or purity), and 2 are label
  mismatches. Checks for those three would remove the 25 and flag the 6; that change would need evaluation on series
  outside this benchmark.
- Where the agent made the same comparison as a published paper, its conclusion agreed in 84% of contrasts, and 5% of
  matched contrasts had an agent error.

## Limits

- Reviewers are language models; a stratified 40-unit workbook for blind review by two domain experts is prepared
  (`human-review/`), and an external audit of a random atlas sample (A0) is proposed in the protocol.
- Units are comparisons with recoverable agent arms in expert-annotated series. Excluded comparisons (A1 44, A2 58)
  bound the A1 E1 rate between 14% and 45% and the A2 rate between 38% and 95%.
- A1 labels come from the ClockBase authors' laboratory (clock training set); A2 series are disease cohorts, many
  longitudinal, and are not a random sample of the DNA methylation atlas.
- In clean mouse designs a naive designer splitting one GEO field passes the design rules (all 177 assessable splits of
  200; none flagged), so A1 design checks bound errors but do not show skill at choosing comparisons.
- Amendments made after results (4) are reported side by side with the analyses as run.

## Correction (2026-10-07)

The released-atlas estimate for human DNA methylation (33 of 44 units) kept only error types the census confirmed for
the agent's own analysis. Two released units of GSE62867 (tissues of the same six patients) were tested as paired by
the agent and are released as unpaired Welch tests; counting them, the estimate is 35 of 44 (80%; cluster-adjusted
95% CI 46-97%). The mouse estimate is unchanged. The sentence "Checks for those three would remove the 25 and flag the
6" was a projection; the measured result is in `../results_2026_10_07/DESIGN_CHECKS.md` (against the corrected 35, checks on
recorded fields leave 13 of 44 units in error, with two model readers 8 of 44).
