# Group-choice and result benchmark of the ClockBase agent: protocol v2 (2026-10-06)

Question: how often does the ClockBase agent compare the wrong groups, or report a result that the right groups do not
give? The clock models are not the target; the agent's design (which samples form each arm, and which arm is the case)
and its reported result are.

Version history: draft v1 (commit 7b84aed37, revised in 0289f3883) was reviewed by three independent critics before
any classification (`DATA/groupbench/acquire_workflow_result_2026-10-06.json`). This v2 replaces it. No benchmark
classification had been run when v2 was written. Information already seen at that point is listed at the end.
Changes after classification starts go to a dated amendment section with the reason, and both analyses are reported.

## Estimands

- **E1, raw agent error rate:** design or result errors in the agent's comparisons as the agent made them (its own arms
  and its reported direction, before any verifier repair).
- **E2, released-atlas error rate:** the same check on the comparisons as released (arms and direction after repair),
  restricted to rows released as VERIFIED or PROVISIONAL.
- **Verifier sensitivity and false alarms:** of E1 errors, the share the verifier repaired, withheld or flagged; of E1
  correct comparisons, the share it repaired or withheld; agreement of its inversion calls with the expert-label
  orientation.

A1 and A2 are reported separately; no pooled rate. Claims are limited to two-sample comparisons in series with expert
sample labels, by contrast type. No atlas-wide or human RNA-seq rate is claimed from this benchmark (see Benchmark A0).

## Units

The unit is a distinct unordered pair of arm sets within a series. Names, orientations, clocks and test variants are
attributes; a pair seen in both orientations is assessed once, with an orientation-conflict flag. Counts are reported at
unit and row level.

Agent arms are the arms in the membership record (before repair) when they are the agent's own: `recovery` (the single
split of the study's samples reproducing the stored P to relative 1e-6), `replay` with effect reproduction PASS for at
least one clock, or both with identical ordered arms. Other comparisons (design reconstruction, replay without
reproduction, no unique membership, other test families) are counted per series and enter the exclusion bounds.
Recovery searches only splits defined by metadata fields, so recovery-only units cannot contain irregular arms; they
are reported as a separate stratum.

## Ground truth (Benchmark A)

- **A1, mouse RNA-seq:** Tyshkovskiy et al., *Nature* 654:173-188 (2026), rodent meta-dataset annotation (Zenodo
  10.5281/zenodo.18763485, `Data_annotation_absolute_rodents.xlsx`, SHA-256 d6d6cf25...); 42 of its 96 series are in
  the atlas. Sample names were mapped to GSMs one-to-one, with mapped counts equal to Supplementary Table 1B and each
  mapped sample's tissue, age, sex and intervention checked against its GEO record (`a1/mapping/`).
- **A2, human DNA methylation:** ComputAgeBench (`computage/computage_bench`, revision 3cf5268), keyed by (dataset,
  GSM); 34 datasets are in the atlas.

A1 labels come from the ClockBase authors' laboratory and the tAge clock training set; a provenance table states for
each A1 series whether ClockBase authors generated or curated it and whether it reached the agent. A1 is reported with
and without those series.

### Coding table (frozen before classification)

- A1: tissue = `Tissue.details`; intervention parsed into components (genotype, compound, diet, exercise, hypoxia,
  each with dose, schedule and start age) by two independent reads, reconciled; background strain = `Strain` with
  "X control" mapped to the background of X, the genotype moving to the intervention components; age group = the
  series' GEO age field when present, otherwise ages clustered by single linkage with a 25% relative gap; sex
  "Unknown" = missing.
- A2: condition, refined by the GEO disease field where ComputAgeBench pools subtypes (for example Crohn's disease and
  ulcerative colitis); tissue x cell type; sex "U" or conflicting = missing; age in years (continuous).
- Every label is checked against the sample's GEO record by two reads; a conflict sets that field to missing and is
  listed.

### Declared contrast (frozen before classification)

For each unit, two independent reads of comparison_name and condition_name (not of the arms) record the declared
factor (an intervention component, tissue or cell type, sex, background, age group, condition, or a named GEO field),
the declared case level and the declared reference level. condition_name names the case level (the label under which
the atlas reports the effect); comparison_name word order is used when condition_name names no level. Disagreements or
contradictory names are resolved by a third read; a unit that stays ambiguous is classified **ambiguous label** and
reported, outside the inversion count.

### Unlabelled samples and subjects

For each unlabelled sample in an arm, from GEO: (a) same GEO group as labelled samples of that arm (takes their label,
flagged curator-excluded); (b) repeats a subject already in the comparison (subject rule); (c) another GEO group (a
distinct level). Subject identifiers (individual, patient, animal, twin pair) come from GEO fields or title tokens by
two reads. A unit is not assessable when more than 20% of an arm cannot be assigned by (a)-(c).

## Classification (deterministic rules, then adjudication)

1. **Not assessable:** arms unassignable (above); or the declared factor is recorded neither by the annotation nor by a
   named GEO field (unencoded factor). Both counted in the flow diagram.
2. **Degenerate test** (error): identical arms, an empty arm, or stored P = 1 with effect 0.
3. **Wrong factor** (error): the arms separate completely on a recorded factor other than the declared one.
4. **Impure arm / label mismatch** (error): less than 90% of an arm carries its declared level (case arm the declared
   case level, reference arm the declared reference level); for intervention contrasts the reference arm is the arm
   lacking the declared component, which handles nested and factorial references.
5. **Confounded** (error): a nuisance factor differs between arms, meaning total variation distance > 0.2, Fisher
   (categorical) or Wilcoxon (age) P < 0.05, and at least two samples would have to move to remove the imbalance;
   unless the arms equal the investigators' complete GEO groups (recorded as study-design imbalance, not an agent
   error) or the agent's recorded model adjusts for the factor. A2 age: arm medians differ by more than 5 years or the
   standardized mean difference exceeds 0.5, and the compared clock output is not age-adjusted.
6. **Repeated subject** (error): a subject contributes more than one sample to one arm. **Paired as independent**
   (statistical error): at least half the subjects are in both arms and the recorded test is unpaired.
7. **Small-arm inference** (statistical error): a P value from an arm of fewer than 3 samples.
8. **Reported-result inversion** (error): orient the arms by the declared levels; the sign of mean(declared case) -
   mean(declared reference) of the atlas clock values disagrees with the agent's reported direction for that clock.
   When the direction agrees but the stored effect-size sign does not, the unit gets an **effect-sign inconsistency**
   (a record error reported separately, not a design error).
9. **Correct** otherwise. Flags: incomplete arm (below 80% of the expert stratum), pooled variants of one component,
   pooled ages within non-overlapping age groups.

Before the agent is scored, the rules are run on simulated planted errors (labelled SIMULATED; seed 20261006; code and
inputs recorded): arms swapped with reported direction kept, two strata merged, samples of another tissue or age group
added, the declared label replaced, half an arm dropped, a subject duplicated. Each error type must be detected at 90%
or more; otherwise the rules are revised (and the revision recorded) before scoring. A naive designer (every
single-field split of the GEO characteristics) is scored with the same rules as a reference point.

### Adjudication (census)

Every assessable unit is reviewed, in a shuffled queue that shows neither the rule class nor any verifier output.
Three model reviewers with distinct briefs (design and groups; GEO metadata and annotation; skeptic arguing for error)
each record a judgment from the expert labels, GEO metadata and names before seeing the agent's recorded design text,
then a final judgment; both are kept. A rule error is overturned only when all three reviewers give the same specific,
checkable reason; a rule Correct becomes an error when at least two reviewers find one. Outcomes: confirmed error,
overturned, unresolved, annotation defect, ambiguous label. Primary = confirmed / assessable; upper estimate =
(confirmed + unresolved) / assessable. Fleiss' kappa and the overturn reasons are reported. A stratified random 40
units (10 per rule class, seed 20261006) are packaged for review by two domain experts outside the ClockBase team who
do not see model calls (Benchmark A0 below); until they report, the model adjudication is labelled as such.

## A-result endpoint

For each expert single-factor contrast (each non-control level against the level lacking only that component, within
the same tissue, background and age group; time courses: oldest against youngest; sexes pooled unless the GEO groups are
sex-specific; listed before matching) targeted by an agent unit with the same declared factor and stratum, a
two-sided Welch test on the atlas clock values of the expert's groups is compared with the agent's reported result:
direction discordance among contrasts significant (P < 0.05) on either side, and discordance of the significance call.

## Statistics and reporting

- Primary: confirmed errors / assessable units for E1 and E2, from a random-intercept logistic model with series as
  cluster (profile-likelihood 95% CI); with fewer than 5 errors or non-convergence, Clopper-Pearson on effective n =
  n / DEFF (DEFF = 1 + (m - 1) ICC, ANOVA ICC floored at 0). Co-primary: the proportion of series with at least one
  confirmed error (Clopper-Pearson).
- Secondary: unadjusted Clopper-Pearson; series bootstrap (BCa, 10,000, seed 20261006); leave-one-series-out range
  (at least without GSE131754 and GSE190939); rates by error type, contrast type (age in controls, lifespan
  intervention, genotype, disease condition, other) and membership stratum; rule-based (pre-adjudication) rates;
  sensitivity at total-variation thresholds 0.1 and 0.3 and label coverage 50% and 100%.
- Exclusion bounds: the primary rate recomputed counting every excluded comparison of the covered series as an error,
  and as correct.
- Flow diagram per set (all atlas comparisons in the covered series to final calls), per-series table, provenance
  table, simulated-error detection table, naive-designer table, reviewer kappa.
- Wording: the agent's error rate is called low only if the one-sided 95% upper bound of the cluster-adjusted E1 rate
  is below 10%; the released atlas only if that bound for E2 is below 10%. Otherwise the numbers are reported without
  that word.

## Benchmark A0: external audit sample (proposed; requires people outside the team)

A random sample of about 350 units, at most one per series, stratified by modality, membership source (including
design reconstruction and unrecoverable rows) and verifier class, curated blind from GEO and the source paper by two
curators outside the ClockBase team and the Gladyshev laboratory, a third resolving disagreements; weighted back to
the atlas composition. This is the only route to an atlas-wide or human RNA-seq error rate. The sample and the
curation packages are prepared; running it is an author decision.

## Benchmark B: statistical results against published results

Ground truth: `ground_truth_t1_v2_final_2026-10-05.jsonl` and the per-dataset tAge comparisons of Tyshkovskiy et al.
2026 Fig. 5a and Extended Data Fig. 8a (two independent extractions from the figure images and legends, reconciled).
The unit is the published contrast, clocks nested; intervals cluster by paper.

- Coverage endpoint: published contrasts in atlas series with any matching agent unit; unmatched contrasts are scored
  by reason (no agent comparison; an agent comparison with the same declared contrast but different samples, which is
  adjudicated as a candidate group error; reconstruction only).
- Primary match: published groups assigned to samples by two models selecting identical samples; at least 80% of each
  published arm inside the agent arm and at least 80% of the agent arm from the published group. Jaccard >= 0.5 and
  name-only matches are reported separately.
- Calls: two-sided P < 0.05 on the agent's stored P and on the paper's reported P or stated significance; FDR as a
  sensitivity. Report the 3x3 table (higher, lower, no difference) with Cohen's kappa, direction agreement among
  contrasts significant on at least one side, and effect concordance where the paper reports effects.
- Every discordant pair is adjudicated as in Benchmark A into: agent group error, agent label or reported-result
  inversion, agent statistical error, clock-value or preprocessing difference, estimand difference, sample-subset
  difference, power or test difference (same direction), published result or extraction error. An estimand or subset
  difference counts as an agent error when the agent's choice is invalid for the design. Clocks trained on the same
  samples are flagged (in-sample agreement).

## Information already seen when v2 was written

Verifier outputs (final classes, orientation outcomes, flags, repairs) for atlas rows, including those in the benchmark
series, as recorded in the public export; the published-clock benchmark results of 2026-10-05 (`results_2026_10_05`);
the critics' read-only scratch tallies quoted in their reports (for example 16 of 99 A1 agent-arm comparisons with an
orientation FAIL from the verifier, ten of them in GSE190939; a literal run of the draft rules on A2); the mapping
agents' notes on which groups each annotation labels.

## Data integrity

Ground truth comes only from the cited sources; nothing is filled in. Every model read, reviewer judgment and simulated
input is kept with its inputs and labelled. S3 and S4 of the published-clock benchmark are context, not endpoints.

## Amendment 1 (2026-10-06, before any classification; found while writing and testing the rules)

1. Rodent age as a nuisance factor: the arms differ in age when their median ages differ by more than 25% of the larger
   median. Single-linkage clustering of exact ages was dropped because it splits cohorts whose ages are jittered
   (GSE165409: 62-95 days in one young group). Declared age contrasts still use the declared-contrast rules.
2. Nuisance factors are of two kinds. Cohort factors (sex; and A2 age, with its own rule) can be imbalanced by
   recruitment or chance, so a partial imbalance counts only with a permutation P < 0.05 (2,000 permutations of the
   total variation distance, seed 20261006). Design factors (tissue or cell type, background, intervention components,
   condition, rodent age group) are set by the investigators, so a partial imbalance counts when the total variation
   distance exceeds 0.2 and at least two samples would have to move. Complete separation counts for every factor.
   Reason: with three or four samples per arm no permutation P reaches 0.05, so a case arm mixing liver and kidney would
   otherwise pass.
3. Unencoded factor: a declared contrast is unencoded when fewer than 80% of an arm's samples have a value in any field
   the declared rules read. Arms whose samples record the field but carry neither declared level are a label mismatch.

## Amendment 2 (2026-10-06, after the first simulated run, before the agent was scored)

The first simulated run (simulated-v2) exchanged the stored arms and recomputed the reported direction from them. For 18
host units whose stored arm order was already opposite to their names (stored direction and effect sign disagree), the
exchange restored the declared orientation, so no error was planted and detection read 83% (88 of 106). This was a
construction fault of the simulation, not of the rules. The planted inversion now sets the agent's case arm to the arm
holding the declared reference level and recomputes the reported direction from the agent's arms. Both runs are kept.

## Amendment 3 (2026-10-06, before the census results were read)

A-result endpoint, as implemented: the expert groups of a contrast are the complete groups that the frozen declared
rules select in the series (the declared levels), not a separately enumerated list of expert single-factor contrasts.
The declared rules already express the levels in expert or GEO fields, and every unit has them; a separate enumeration
would repeat that work for the same contrasts. The endpoint is reported for units the census finds correct.

## Amendment 4 (2026-10-06, after the census; reported alongside the analyses as run)

1. Interval validity: the profile-likelihood interval of the random-intercept model is used only when it contains the
   model estimate and has positive width; otherwise the protocol's fallback (Clopper-Pearson on the design-effect
   adjusted n) is reported. With few or many events the profile was inverted or collapsed in some analyses.
2. Adjudication-material defect: the census packages described the reported direction as relative to the agent's
   stored case arm. For units whose stored effect-size sign is opposite to their names (the record error the protocol
   reports separately), that description is false, and reviewers called some of them label mismatches or inversions.
   Results are reported as the census ran and with those calls reclassified as effect-sign record errors, together
   with the combined rate counting every effect-sign record error.
3. Age confounding (A2): reviewers applied the protocol's criterion (median ages differ by more than 5 years or
   standardized mean difference above 0.5), using GEO ages where ComputAgeBench gives none, and without the
   study-design exemption that the rules applied. Results are also reported counting age confounding only when the
   arms' median GEO ages differ by more than 5 years.

## Addendum (2026-10-07, after the benchmark results): prototype design checks

The benchmark found three kinds of error in released rows that the verifier does not test for: repeated subjects
within an arm, paired arms tested as independent, and arms of different age compared on a clock that tracks age.
`design_checks.py` is a prototype of checks for them; it is not part of the released verifier, and the checks were
designed from the benchmark's findings, so the benchmark cannot validate them.

1. Inputs are recorded fields only. Subjects: a GEO sample field named as a subject identifier, else the harmonized
   subject field (rules in the module). Ages: the harmonized age, else for human series a GEO age-at-sampling field in
   years. A row is flagged on its released (post-repair) arms.
2. Each flagged row is recomputed with subject means, a paired t-test when at least half of the subjects are in both
   arms (three or more pairs), and age as a covariate under age imbalance; the released call (direction, two-sided
   P < 0.05) is compared with the recomputed call.
3. Rule change after the first run. The first version exempted a row from the age check when its names contained any
   age or time term; scored against the census it missed comparisons whose names carry a time qualifier but contrast
   something else. The second version exempts a row when either of the verifier's two contrast classifications says
   AGE or DEVELOPMENT_OR_TIME or the names carry an age-group word, and adds the GEO age fallback. Both runs are kept
   (`design-checks/v1`, `design-checks/v2`); results are reported for v2.
4. In-sample scoring (`design_eval.py`): census-confirmed types against check findings for the released benchmark
   rows.
5. Validation outside the benchmark (`design_validate.py`, seed 20261007): released units (distinct arm pairs, one per
   series, arms of at most 120 samples) from series not annotated for A1 or A2, per modality 20 flagged for subjects,
   15 flagged for age only and 60 without a finding, plus every released benchmark unit. Two model readers
   (openai:gpt-6-astra, openai:gpt-6-sol) see the series text, the comparison names and the GEO fields of the arms'
   samples, and no check output, clock value or harmonized field. Each answer is turned into findings by the same rule
   functions. A finding is the reference when both readers give it, absent when neither does, unresolved otherwise.
   Reported: the share of flagged units the readers confirm, and the share of unflagged units in which both readers
   find a design issue the recorded fields do not show.
6. Every comparison of the ranked tables (Table 1, Table 2, Figure 3, Table S14) with arms of at most 120 samples is
   read the same way.
7. The readers are language models; their agreement is reported, and no person has reviewed their answers.
8. After validation (v3). The readers confirmed 0 of 2 sampled flags that came from the harmonized subject field and
   9 of 15 flags of the mouse age check (five wrong harmonized ages, one developmental contrast). The default variant
   therefore uses GEO subject fields only and applies the age check to human series only; `--variant v2` reproduces
   the earlier behaviour. The validation strata were drawn on v2 flags; v3 flags are a subset of them, and the
   confirmed shares for v3 are computed on the same reads.
9. Released rows whose own agent test was paired, a mixed model or a model with age are compared three ways (agent's
   stored call, released call, design-respecting call) in `design_eval.py --atlas`.
10. Candidate release estimand (2026-10-09; stated and committed before it was scored against published results).
    For every released row, on its released arms: (a) a row flagged for repeated subjects or pairing uses subject
    means, with a paired t-test when the arms are paired and otherwise a Welch test on subject means, or least squares
    on arm and age when (b) applies; (b) a human row on a clock that tracks age, with ages recorded for both arms and
    no declared age or time contrast, uses least squares of the clock on arm and age, and keeps the released Welch
    test when the age term is not estimable; (c) every other row keeps the released Welch test. A row whose candidate
    statistic cannot be computed (one subject in an arm, fewer than three pairs, age separating the arms) has no
    candidate call and counts as not covered. Scoring uses the existing metrics and nothing else: the published-clock
    benchmark T1 with systems S1 (historical row), S2 (released row) and S2c (candidate on the released row), and
    Benchmark B with the agent's, the released and the candidate call. No other variant of the estimand is tried.
11. Conclusion-level endpoint by stage (2026-10-09, after the results above; post hoc). The A-result endpoint of
    amendment 3 (the reported call against a test on the complete declared groups) is extended in two ways: the same
    row is scored at three stages (agent's reported call, released call, candidate call, each oriented to the declared
    case), and a second reference analyses the declared groups with the candidate estimand instead of the Welch test,
    so that a paired or age-imbalanced design is not scored against an unpaired, unadjusted reference. Three
    stratifications by fields the verifier already carries were examined and all are reported: final class (VERIFIED,
    PROVISIONAL), the LABEL_UNCONFIRMED soft flag, and an agreed CASE_NOT_CONDITION adjudication issue. Intervals:
    Clopper-Pearson on the design-effect-adjusted sample size, with the unadjusted interval beside it.
12. Significant-finding metrics by modality (2026-10-09, requested by K.A.Y.; defined before computing). The benchmark
    is summarised on findings that are significant in the reference, separately for DNA methylation clocks and for
    RNA-seq clocks. Score of a ClockBase result: the signed z value of its two-sided P (the sign of its estimate times
    the normal quantile of 1 - P/2). (a) Direction AUROC: among reference results significant as higher or lower, the
    probability that a "higher" result has a larger signed z than a "lower" one (ties one half). (b) Same direction
    among results significant on both sides. (c) Detection AUROC: the absolute z for reference-significant against
    reference-not-significant results. (d) Replicated: reference-significant results that ClockBase also finds
    significant in the same direction. (e) Significant in ClockBase where the reference is not. References: published
    results (T1; paper-cluster bootstrap) and the expert-defined groups analysed with the design-respecting estimand
    (item 11; series-cluster bootstrap); 2,000 replicates, seed 20261009. Stages: raw agent, released, candidate, and
    for published results the published groups recomputed (S3). The stricter estimands above stay in the report.
13. Effect agreement and overflow guard (2026-10-09). (a) Effect agreement with papers: for published results that
    report an effect in years (case minus control, as reported; records whose extraction notes call the sign
    convention uncertain or whose sign contradicts the stated direction are left out), the published effect is set
    against the ClockBase mean difference for the published groups in the paper's measure (S3): Spearman and Pearson
    correlation, share with the same sign, median absolute difference; paper-cluster bootstrap for the correlations.
    (b) Overflow guard (post hoc; found while computing (a)): some series carry clock values of the order of 1e39
    years. As a labelled sensitivity, clock values with an absolute value above 1,000 years are treated as missing and
    the affected published-group results are recomputed; results are reported with and without the guard. The released
    atlas has no such guard.
14. Scoring the released run (2026-10-10, after the results above). The candidate estimand of item 10, the overflow
    guard of item 13(b) and a two-reader individuals audit of ranking candidates (the reader prompt validated in items
    1 to 9) entered the verifier; its run final-2026-10-10 releases the result. The benchmark scores that run with the
    functions and metrics above and nothing new: `export_as_checks.py` writes, for every row released in
    final-2026-10-05, the statistic the row has in final-2026-10-10, and the existing functions read it in place of the
    candidate. Stage names in the tables of this run: raw agent (S1), verifier with a Welch test for every row
    (final-2026-10-05; S2), released (final-2026-10-10; S2c). Three rules, fixed here: (a) a design flaw counts as
    corrected only when the released statistic respects it (individual means, a paired test, or age as a covariate);
    a soft flag without a corrected statistic (age separating the arms) does not count; (b) rows that
    final-2026-10-10 no longer releases (fewer than three individuals or pairs in an arm, overflow values leaving an
    arm empty) leave the denominators and are counted beside them; (c) the column "released + two readers" applies the
    validated readers to every benchmark comparison and is hypothetical, because the run reads only ranking
    candidates. Differences from the prototype of item 10: individuals come from a GEO sample field named as an
    individual identifier or from the audit, not from processed metadata; when age separates the arms the row keeps
    its Welch statistic with a flag instead of having no call. No metric, threshold or stratum was added or changed.
15. Name (2026-10-10). The manuscript reports the published-clock benchmark (`../README.md`) and this group-choice
    benchmark together as ClockBase-Bench. Systems S3 and S4 of the published-clock benchmark are reported as ClockBase
    with the paper's groups and as language models without data; the tidy table of a scored run lists them with the
    constant rule S0, on the results scored for every system. No analysis was added or changed for the name.
16. Agent reruns with other language models (2026-10-11, requested by K.A.Y.; written and committed before any rerun
    was scored). Question: does the discovery task give results of similar accuracy when another agent framework and
    other language models perform it? Code: `../agentbench/`.
    (a) Agents. Codex CLI 0.161.0 (`codex exec`, workspace-write sandbox without network for commands) with GPT-6
    Astra, GPT-6 Sol, GPT-6 Luna and GPT-5.5 through the OpenAI API at reasoning effort medium: the models of the
    no-data baseline (S4) that this Codex version runs (it has no model metadata for GPT-4.1). The Codex home of the
    runs holds no user or project instructions, plugins, memories or web search. A run is excluded and reported when
    its event log shows a web search, a spawned sub-agent, or a command that names a path outside its run directory,
    the Python environment and system locations.
    (b) Series. Human DNA methylation series with at least one published result on one of the five clocks of the task
    (Horvathv1, PhenoAge, YingDamAge, DunedinPACE, GrimAgeV2) whose groups the two group-assignment models selected
    identically (status AGREED), and with a saved trace of the historical agent: 61 series (164 results in 65 series;
    the 4 series without a saved trace are left out).
    (c) Task. The historical task text of the series, recovered from its saved trace, with the data and result paths
    pointing into the run directory and the model label replaced. A harness note precedes it (how to run Python in the
    run directory; that running its own scripts is expected), and one deliverable is added: `comparison_samples.csv`,
    the sample identifiers of the case and control arms of every comparison, so that a run can be scored without
    membership recovery. Inputs are the three files the historical task names, laid out as in the verifier's trace
    replay: the series' row of the dataset overview, its sample metadata, and its clock outputs
    (`replay-dnam-001/inputs`). One run per series and model, time limit 30 minutes; a run that ends without a
    results table is repeated once, and runs that fail twice are reported.
    (d) Scoring. Against the published results of (b), matched by sample overlap only: a comparison of the run matches
    a published result when it is on the same clock and each arm has a Jaccard index of at least 0.5 with the
    corresponding published group, in either orientation; the match with the largest smaller index is taken. Each
    run is scored twice: raw (the direction and P value the agent reported, oriented to the published contrast) and
    recomputed (the verifier's statistic on the run's arms: the design-respecting estimand with the verifier's clock
    values, the overflow rule, at least three samples, individuals or pairs per arm). The membership recovery,
    orientation adjudication and input checks of the full verification are not applied to reruns, which save their
    arms. Significance uses the paper's threshold, 0.05 when none is stated.
    (e) Metrics, as in the published-clock summary and item 12: coverage; three-way call agreement and balanced
    accuracy; direction agreement on published significant results; direction AUROC; same direction when both are
    significant. Reported on each system's covered results and on the results covered by every compared system
    (the historical agent's row S1, the released row, ClockBase with the paper's groups S3, the no-data models S4
    and the reruns), with the paper-cluster bootstraps of those summaries.
    (f) Not covered: the annotation part (it needs the declared-contrast coding for every rerun), RNA-seq, and the
    historical CodingAgent itself run with another model.
17. Name changed (2026-10-11, K.A.Y.): the benchmark is reported as ClockEffectBench; the name of item 15 is no longer
    used. No analysis changed.
18. Outcome of item 16 (2026-10-11). 244 runs; report `../agentbench/results_2026_10_11/REPORT.md`. Two details
    were fixed while the runs were in progress and before any score was read: the run audit counts a path as outside
    the run directory only when its first component is a directory at the filesystem root (strings such as
    `/clock_output/x.csv` that an agent joins to a variable are not paths), and the bootstraps are those of the
    existing summaries (1,000 replicates with seed 20261005 for agreement, 2,000 with seed 20261009 for the AUROC).
    The runs went through the OpenAI API with a key taken from the environment, because the Codex home of the
    experiment holds no login.
19. The paper's groups with the released test (2026-10-12, post hoc, requested by K.A.Y. after reading the results:
    "figure out why" 33 of 54 published significant results were recovered). Finding that prompted it: S3 tests the
    paper's groups with an unpaired Welch test. Of its 21 misses, 19 are in the published direction, and the papers
    behind most of them used paired tests, mixed models with an individual term, or regressions with age. The
    released verifier does not use the Welch test for such comparisons. Two systems are therefore added beside S3,
    with the rules of the released verifier and nothing tuned on this benchmark:
    S3c: the design-respecting estimand of `audit_system/design_estimand_2026_10_10.py` on the paper's groups
    (individuals from a GEO sample field, age as a covariate under its conditions, overflow values dropped, at least
    three samples, individuals or pairs per arm; a result with fewer has no call and leaves the denominators).
    S3r: the same, with individuals taken from the two-reader individuals audit where the series has no GEO individual
    field: the verifier's reader prompt, GPT-6 Astra and GPT-6 Sol, arms of at most 120 samples, individuals used only
    when both readers report the same shared individuals (the rule the verifier applies to ranking candidates).
    S3c was computed once in exploration before this item was written (40 of 52); S3r was specified here before any
    read. S3 stays in every table. Code: `published_released_test.py`. Also reported, without changing a call: the
    published results whose paper states a one-sided test, and for every published significant result that S3r does
    not recover, its direction, P value, sample sizes and the paper's test.
20. Agent reruns of the ranked comparisons (2026-10-14, requested by K.A.Y. so that
    top-ranked interventions and candidate rankings are checked under other agent frameworks and models; written and
    committed before any run). Item 16 measured agreement with published results on human DNAm series; the ranked
    comparisons are mouse RNA-seq, which has no published-result reference, so this item measures whether the ranked
    comparisons recur when another configuration analyses the same series.
    (a) Comparisons: the 60 comparisons of Tables 1 and 2 of verifier run final-2026-10-10 (30 with higher and 30 with
    lower clock outputs; `comparison-ranking/table1_pro_aging.csv`, `table2_anti_aging.csv`), with their verified arms
    (`case_gsms`, `control_gsms` of the per-study records in `public-export/studies/mm_rna/`).
    (b) Runs: every series of these comparisons with the four models and the Codex configuration of item 16 (same
    CODEX_HOME, medium reasoning effort, workspace sandbox, 30-minute limit). Task text: the released mouse template
    `all_samples_local` (`src/parallel_mm_rna.py`; it forbids online retrieval) with the series' title, summary and
    overall design from its public GEO record read on 2026-10-14, the paths pointed into the run directory, the harness
    note and the samples deliverable of item 16. Inputs: the series' row of `gse_rna_mm.csv`, its sample metadata
    (`GEO_mm_rna/processed_metadata_preharmonize/`) and its clock outputs (`replay-mm_rna-001/inputs/`). A single
    template is used for every series; the historical runs used three template versions.
    (c) Scoring: a ranked comparison is found by a run when one of the run's comparisons has arms whose Jaccard indices
    with the verified arms are both at least 0.5, in either orientation (the best overlap is kept). For a found
    comparison the verifier's statistic is recomputed on the run's arms with the ranked comparison's clock
    (Hazard_Mouse_All_All_scaleddiff_EN_tAge), the design-respecting estimand and the overflow rule, oriented to the
    ranked comparison; it recurs when its mean difference has the ranked direction and P < 0.05. Also reported: the
    direction and P value the run reported for that clock, when it did. As a check of the recomputation, the
    statistic on the verified arms must equal the released row.
    (d) Reported per model and for the four models together: of the 30 higher and 30 lower comparisons, the numbers
    found, found in the ranked direction, and recurring; and the number of comparisons found by at least two and by
    all four models. Runs are audited as in item 16; a run that names a path outside its run directory, or whose
    results cannot be read, counts as finding nothing and is reported.
21. Outcome of item 20 (2026-10-14). 224 runs, all exited 0, none excluded, no web search; Codex CLI 0.162.1 (it had
    updated itself since item 16, which used 0.161.0; configuration otherwise unchanged). Two details of the scoring
    code were fixed after item 20 was committed and before any run was prepared: when a series' record holds several rows with the label of
    a ranked comparison, the row whose released statistic is the table row's is used (GSE185761); and where GEO has no
    individual field, the two readers' individuals of the ranked comparison enter the recomputation, as in the verifier,
    so that the statistic on the verified arms equals the released row for all 60. Report:
    `../agentbench/results_2026_10_14/REPORT.md`.
22. Claude Code as a second agent framework (2026-10-18, requested by K.A.Y. to test another
    agent framework, Claude Code; Biomni was not run because its source is not available to us; written and
    committed before any run). The tasks, inputs, series, scoring and reported metrics are those of item 16 (61 human
    DNAm series against published results) and item 20 (the 56 series of the 60 ranked comparisons); only the agent
    framework and model change.
    (a) Agent. Claude Code CLI 2.1.302 in print mode (`claude -p`, stream-json events) with Claude Opus 5.5
    (`claude-opus-5-5`), signed in with the laboratory's Claude account. The run's configuration directory
    (`CLAUDE_CONFIG_DIR`, also HOME) holds only the sign-in credentials: no user or project instructions (CLAUDE.md),
    memories, skills, plugins, hooks or MCP servers (`--strict-mcp-config` without a configuration). Tools are limited
    to the built-in coding set (`--tools Bash Read Write Edit Glob Grep NotebookEdit TodoWrite`): no web search or
    fetch, no sub-agents. Permissions are not prompted (`--dangerously-skip-permissions`); instead the process runs in
    bubblewrap with the filesystem read-only except its run directory, its configuration directory and a private
    /tmp, which matches the workspace-write sandbox of Codex. Unlike that sandbox, commands keep network access,
    because the agent process itself needs the API; network use by commands is counted, as in item 16.
    (b) Audit. A run is excluded and reported when its events show a web tool, a sub-agent, or a command or a file-tool
    argument (Read, Write, Edit, Glob, Grep, NotebookEdit) that names a path outside its run directory, the Python
    environment and system locations: the rule of item 16 extended to Claude Code's file tools, which reach files
    without a shell command.
    (c) Runs. One run per series, 30-minute limit, the same run directories' inputs and task text as items 16 and 20
    (prepared again under `claude-claude-opus-5-5/`); a run that ends without a results table is repeated once.
23. Outcome of item 22 (2026-10-18). 117 runs (61 published-result series and 56 ranked series), all exited 0, none
    excluded: no web tool, sub-agent, network command or path outside the run directory; median 32 and 30.5 seconds
    per run; Claude Code CLI 2.1.302 with Claude Opus 5.5. Published results (item 16's scoring): covered 107 of 158;
    call agreement 0.83 (n = 87) raw and 0.86 recomputed; direction agreement 0.95 (n = 37); within the range of the
    four Codex configurations (call 0.77-0.90, direction 0.90-0.94; historical agent 0.82 and 0.87). Ranked
    comparisons (item 20's scoring): found 52 of 60, all in the ranked direction, 51 recurring (26 higher, 25 lower);
    the recomputation check equals the released row for all 60. Across the five configurations of the two frameworks
    (`agentbench.ranked combine`), 55 of the 60 recur in at least two and 45 in all five; none is reversed; the
    four Codex configurations alone give 55 and 47, as reported in item 21. Report:
    `../agentbench/results_2026_10_18/REPORT.md`.
24. Biomni as a third agent framework (2026-10-18, requested by K.A.Y.; written and committed before any run,
    pilot included). Item 22 said Biomni's source was not available to us; it is public (github.com/snap-stanford/Biomni,
    Apache-2.0) and is run here as a biomedical agent framework. Tasks, inputs, series, matching (Jaccard 0.5),
    scoring and reported metrics are those of items 16 and 20; only the agent framework changes.
    (a) Agent. Biomni A1 at commit 400c1f366b96a35ca253e13c9b06c5076af41d65 (package version 0.0.8), installed from
    that clone into its own virtual environment with only the packages A1 imports and pandas, numpy, scipy,
    statsmodels and matplotlib (the large bioinformatics environment of Biomni is not installed; versions in the
    results folder). Model `gpt-6-astra` through the OpenAI API at medium reasoning effort: the model of the first Codex
    configuration, so this configuration differs from it by the framework. Biomni's defaults are otherwise kept: the
    resource retriever (the agent's model selects tools, libraries and know-how for the prompt), 600 seconds per
    executed block, no self-critic, recursion limit 500.
    (b) Adaptations, made in our runner (`agentbench/biomni_runner.py`), never in Biomni's files.
    (1) Model client: Biomni sends `gpt-5` models through the Responses API and drops the stop and temperature
    parameters those models reject; the same handling is applied to `gpt-6-astra`, with medium reasoning effort and a
    callback that records token usage.
    (2) Data lake: Biomni's own option `expected_data_lake_files=[]` skips the download of its data lake and benchmark
    files. The task needs only its three input files, and the prompt lists no data-lake items.
    (3) Trajectory: the runner drives A1's compiled graph with the inputs and configuration of `A1.go`, after the
    retrieval step of `A1.go`, and writes each message to events.jsonl as it is produced, so that a run cut by the limit
    keeps its trajectory.
    (4) Credentials: the API key is removed from the process environment once the client is built, because executed
    code shares A1's process.
    (c) Sandbox. Bubblewrap with the filesystem read-only except the run directory and a private /tmp, as in item 22.
    Network stays available, because A1 calls the API from the process that executes the code.
    (d) Task text. The run directories are prepared under `biomni-gpt-6-astra/` with the task text of items 16 and 20.
    The harness note names the interpreter of Biomni's environment, because A1 executes Python in its own process.
    (e) Audit. The commands of a run are its executed blocks (Python, R or Bash). A run is excluded and reported when an
    executed block does any of the following:
    - names a path outside its run directory, Biomni's environment and system locations;
    - matches the network pattern of item 16 (a URL, curl, wget, git clone or pip install);
    - imports or calls a Biomni tool module (`biomni.tool.*` other than `support_tools`). These modules query external
      databases and web services or run software outside the task's inputs.
    (f) Runs. One run per series with a 30-minute limit; a run that ends without a results table is repeated once.
    A pilot of one run per set checks the output format. Pilot runs count as the runs of their series unless the
    harness has to change, in which case they are rerun and the change is reported here.
25. Amendment of item 24 after its first attempt and before any scoring (2026-10-18; written and committed before
    the second attempt). The first attempt prepared and ran all 117 runs. All exited 0, but 30 of them (22 benchmark,
    8 ranked) ended without a results table. In 37 of the 117 runs, one model message held both an executed block and
    a `<solution>` tag; A1 treats the tag as the end, so the block that would have saved the results never ran. A1
    builds its client with the stop sequences `</execute>` and `</solution>`, which end a reply at its first block. The
    Responses API, which Biomni uses for gpt-5 models and item 24 (b1) extended to gpt-6, does not accept stop
    sequences, so they were dropped. Amendment: the client cuts each reply before the first stop sequence, as an API
    that applies them returns it (`truncate_at_stop` in `agentbench/biomni_runner.py`, tested in `test_score.py`).
    Everything else in item 24 is unchanged. All 117 runs are prepared again under new roots and rerun. The rerun is
    scored. The first attempt is kept unscored (`/Entropy/kying/clock-agent-runs-biomni-2026-10-18` and
    `...-ranked-biomni-2026-10-18`) and reported with its failure count.
26. Outcome of items 24 and 25 (2026-10-18). Second attempt, which is the scored one.
    - Runs: 117 (61 published-result series and 56 ranked series). All exited 0 with a results table; none was
      excluded (no path outside the run directory, no network command, no Biomni tool module); none was repeated.
    - Biomni A1 at 400c1f3 with gpt-6-astra. Median 95 and 73.5 seconds per run.
    - Published results (item 16's scoring): covered 82 of 158. Call agreement 0.81 (n = 68) raw and 0.78 recomputed;
      direction agreement 0.92 (n = 26); direction AUROC 0.977. Eight runs reported no testable comparison.
    - On the 70 results that the historical agent and the Biomni rerun both covered, the calls agreed in 0.82
      (historical agent) and 0.83 (Biomni).
    - Ranked comparisons (item 20's scoring): found 56 of 60, 55 in the ranked direction, 54 recurring (27 higher,
      27 lower). The recomputation check equals the released row for all 60.
    - The one comparison found with the opposite sign is ouabain (GSE122080). Biomni's control arm pooled the six old
      saline-treated mice with six young ones (Jaccard 0.5, the threshold). On those arms the verifier's Welch
      statistic is +0.15 (P = 0.68). Biomni's own age-stratified model reported the ranked direction (P = 9 x 10^-6).
    - Across the six configurations of the three frameworks: 56 of the 60 recur in at least two and 43 in all six;
      one is reversed in the sense of item 23 (found with the opposite sign), the ouabain case above.
    - `ranked combine` now counts a configuration as a scored directory and a model, since gpt-6-astra runs in two
      frameworks; the five-configuration result of item 23 is unchanged (55 and 45).
    - First attempt (item 25), unscored: 30 of 117 runs ended without a results table.
    - Report: `../agentbench/results_2026_10_19/` (REPORT.md written by the main session).
27. Repeated runs of one configuration (2026-10-22, requested by K.A.Y.; written and committed before any run).
    Items 16-26 vary the agent framework and the model, with one run per series. This item asks whether independent
    end-to-end runs of one configuration reach the same conclusions. Series, inputs, matching (Jaccard 0.5) and
    scoring are those of items 16 and 20.
    (a) Agent. The online agent of the ClockBase website (repository clockbase_agent_website, branch
    online-agent-2026-10-17, commit 1d2ff92b2): a tool loop over the OpenAI Responses API with two tools,
    `run_python` (a fresh Python interpreter in bubblewrap, no network, no file access outside the run's workspace,
    read-only inputs, 120 s per execution) and `finish`. Model `gpt-6-sol`, medium reasoning effort, at most 14 model
    calls of at most 16,000 output tokens, 900 s per run.
    (b) Task text. The website's template for the modality: the human DNAm template recovered from the saved traces
    (`src/parallel_hs_dnam.py`) for the 61 series of item 16, and the released mouse template `all_samples_local`
    for the 56 series of item 20. The design block is the series' title, summary and overall design from its public
    GEO record, read at run time; item 16 used the historical task text of each series, which differs only in that
    block. Paths point into the sandbox. The harness note of the website and the samples deliverable of item 16 are
    added; the note asks for arm `other` for samples that enter a model only as covariate or nuisance strata, which
    scoring ignores, as it ignores any arm other than case and control.
    (c) Inputs. As items 16 and 20, except that 5 DNAm series whose historical matrix was empty or detached
    (GSE199700, GSE72774, GSE72777, GSE72778, GSE75704) receive the re-derived predictions the website serves.
    (d) Runs. Three independent replicates per series: separate model sessions with no shared state, submitted
    through the website's run manager with the shared-key limits raised for this experiment (per run at most
    USD 1.50 for the three replicates together; at most USD 30 in total). One attempt per replicate, no repeats. A
    replicate that ends without results is counted as finding nothing and is reported.
    (e) Scoring of each replicate. Item 16's functions for the 158 published results (raw and recomputed; call
    agreement, direction agreement, direction AUROC), and item 20's for the 60 ranked comparisons (found, found in the
    ranked direction, recurring). Website runs are not audited from an event log; their code ran without network or
    access outside the workspace.
    (f) Stability endpoints.
    - Published results. Each replicate's status per result is: not covered, or the recomputed call (higher and
      significant, lower and significant, or not significant). Primary: among the results covered by all three
      replicates, the fraction with the same recomputed call in all three. Secondary: the fraction of the 158 results
      with the same status in all three replicates, coverage included; and the range of the per-replicate metrics.
    - Ranked comparisons. For each of the 60: the number of replicates (0-3) in which it recurs. Reported: recurring in
      all three, in at least two, in at least one; and the number found with the opposite sign in any replicate.
    - The website's own label per series (significant comparisons reproduced by every replicate, from the website's
      recheck on each replicate's arms), reported as the distribution of labels.
    (g) Every result is reported, including failed replicates and low-stability series. No threshold of this item is
    changed after the runs are seen.
