# ClockBase independent audit system

## Verifier (v3)

Start with [Verify the atlas](../README.md#verify-the-atlas) for the code release plus rebuild bundle.
Run commands from the release root with the installed verifier environment active. Locations and overrides are
in [paths.py](paths.py). `M` below is `dnam`, `hs_rna` or `mm_rna`; `V` is the verifier root
(default `clockbase_data/verifier`), and `NAME` is a new final directory name under `V`.

```bash
python scripts/clockbase_verifier.py plan --final NAME
python scripts/clockbase_verifier.py run --final NAME --models cache-only
```

To execute a stage in the table, substitute its stage name into the full controller command
`python scripts/clockbase_verifier.py run --final NAME --models cache-only --force STAGE`.
Replace `M` with the modality; use `--models live` for a fresh model request. This still runs the normal plan,
including other incomplete stages. Source-dependent stages require the additional inputs described below.

The controller executes the stages below in order. It reuses stages with completion files and skips preliminary
passes when all their readers are already complete. The bundle supplies the upstream and model-stage files the rebuild reads and stage completion files. It does
not contain every intermediate file of those stages. Fresh membership, replay and integrity runs require the complete source predictions, molecular inputs
and traces they read; these are not all in the bundle. Fresh stream-repair inventory classification also needs
omitted development census inputs. The release supports reuse of those saved stages, rather than promising their
complete regeneration from the bundle. Source acquisition requires network access; live model stages require the
credentials described in the root README. `plan` lists each stage's reuse status.

| Stage | Module or script | Inputs | Outputs under `V` unless stated otherwise |
|---|---|---|---|
| `membership-M` | `audit_system.membership_recovery_2026_09_27` | Correction tables, metadata, per-sample clock predictions, saved study outputs | `membership-M-001/` |
| `replay-M` | `audit_system.trace_replay_2026_09_27` (`extract`, `split-inputs`, `run`, `analyze`) | Saved traces and per-sample prediction inputs; replay interpreter | `replay-M-001/`, including split inputs and replay summary |
| `link-M` | `audit_system.replay_link_2026_09_27` | Correction rows and `replay-M-001/` | `link-M-001/` |
| `integrity-dnam` | `audit_system.input_integrity_2026_09_30` (`dnam`) | DNAm source matrices, sample metadata and clock outputs | `input-integrity-dnam-001/samples.csv` and screen records |
| `stream-repair-empty` | `src/dnam_source_inventory_2026_09_30/geo_inventory.py`, `classify.py`; `src.dnam_stream_repair_2026_09_26` | Empty-matrix census, GEO supplements, inventory classifier census inputs and methylome reference | `empty-repair-inventory/`; status, QC and clock outputs under `CLOCKBASE_STREAM_EMPTY` |
| `stream-repair-detached` | `audit_system.input_integrity_2026_09_30` (`detached-list`); the same inventory and stream-repair code | Integrity screen, GEO supplements, classifier census inputs and methylome reference | `detached-repair-inventory/`; status, QC and clocks under `CLOCKBASE_STREAM_DETACHED` |
| `pass-1` | `audit_system.membership_merge_2026_09_27`, `audit_system.row_verdicts_2026_09_27`, `audit_system.row_repairs_2026_09_28` | Recovery, links, replay, DNAm integrity | `pass-1/merged-M`, `verdicts-M`, `repairs-M` |
| `design-M` | `audit_system.design_reconstruction_2026_09_28` | Pass-1 unrecovered rows and source metadata | `design-M-001/` and model ledgers |
| `adjudication-decisive-M` | `audit_system.label_adjudication_2026_09_27` | Pass-1 membership and decisive label-dependent verdicts | `adjudication-M-001/` and model ledgers |
| `pass-2` | The same merge, verdict and repair modules as pass 1 | Recovery, links, replay, design and adjudication | `pass-2/merged-M`, `verdicts-M`, `repairs-M` |
| `adjudication-followup-M` | `audit_system.label_adjudication_2026_09_27` | Pass-2 pending/unconfirmed rows and supported validation sample | Additional `adjudication-M-001/` decisions and `.followup-done` |
| `pass-2-atlas` | `audit_system.atlas_rebuild_2026_09_29` | All three pass-2 repair directories | `pass-2/atlas-v3/` |
| `classification` | `audit_system.intervention_classification_2026_09_29` | Pass-2 atlas and repair records | `classification-001/` and model ledgers |
| `design-audit`, `design-audit-added` | `audit_system.design_audit_2026_09_30` | Pass-2 atlas, classification, repairs; the added run excludes the first run's candidates | `design-audit-001/`, `design-audit-002/` and model ledgers |
| `pass-3` | The final-chain script with the design-respecting estimand, before the individuals audit | Saved upstream/model stages and DNAm re-derived outputs | `preNAME/` (a complete chain; its atlas names the ranking candidates of the two stages below) |
| `design-audit-design-estimand` | `audit_system.design_audit_2026_09_30` | Pass-3 atlas, classification, repairs; excludes candidates of the first two design audits | `design-audit-003/` and model ledgers |
| `individuals-audit` | `audit_system.individuals_audit_2026_10_10` | Pass-3 atlas and repairs, classification, the three design audits, GEO sample records | `individuals-audit-001/` and model ledgers |
| `final-chain` | [Final-chain script](../scripts/verifier_v3_final_chain_2026_09_30.sh); modules below | Saved upstream/model stages and DNAm re-derived outputs | `NAME/` verdicts, repairs, published fates, census, atlas, ranking, export |
| `decision-audit` | `audit_system.model_decision_audit_2026_10_04` | Saved decisions and raw-response ledgers | `NAME/decision-audit/audit.json` |
| `blood-structure-cohort` | `audit_system.blood_clock_structure_2026_09_30` | DNAm integrity screen and cohort source data | `r6-6-blood/cohort-002/` |
| `grimage-input` | `audit_system.grimage_age_input_2026_09_30` | Blood cohort and source clock inputs | `r6-6-blood/grimage-age-input-002/` |
| `blood-structure` | `audit_system.blood_clock_structure_2026_09_30` | Integrity screen and GrimAge input run | `r6-6-blood/structure-002/` |
| `validation-frame` | `audit_system.validation_frame_2026_10_05` | Saved verifier selection and adjudication records | `validation-frame-001/` |
| `tests` | Standard-library `unittest` (the controller's `VERIFIER_TESTS`) | Released test fixtures; bundle traces where available | Test results in the controller log |

The final chain runs merge, verdicts and repairs for each modality. The repairs recompute each row with a Welch
t-test or, where the sample records call for it, with the design-respecting estimand of
`audit_system.design_estimand_2026_10_10` (means per individual, a paired t-test, or age as a covariate; clock
outputs above 1,000 in absolute value are dropped); individuals come from a GEO sample field or from the two-reader
individuals audit, and the Welch result stays in the row's record (`recomputed_declared`). The repairs are followed by
`audit_system.published_outputs_status_2026_09_29` (repair records and published-item inputs → published fates),
`audit_system.error_census_2026_09_29` (repairs and published fates → census),
`audit_system.atlas_rebuild_2026_09_29` (repairs → atlas),
`audit_system.comparison_ranking_2026_09_30` (atlas, classification and design audits → ranking), and
`audit_system.v3_public_export_2026_09_30` (final chain and verifier stages → per-study records and hash manifest).
These commands run through the controller above; the shell script also accepts `V NAME` directly:

```bash
CLOCKBASE_PYTHON=python bash scripts/verifier_v3_final_chain_2026_09_30.sh clockbase_data/verifier NAME
python -m audit_system.model_decision_audit_2026_10_04 --out NEW_DECISION_AUDIT
```

The direct chain refuses an existing final directory. The decision audit makes no requests; failed model calls
have no response to re-derive. For a particular stage, use the controller's full command
`python scripts/clockbase_verifier.py run --final NAME --models cache-only --force STAGE`.
This forces that stage within the normal dependency-ordered run; it is not an isolated-stage selector.
It does not supply omitted source inputs or clear an output directory whose module refuses reuse.
The exact stage commands, arguments, environment overrides and completion files are defined in
[the controller's `stages()`](../scripts/clockbase_verifier.py).

## Additional released entry points

[The pipeline guide](../PIPELINE_2026_09_10.md) covers new source ingestion and its separate scientific limits.
The standalone source-design workflow uses provider-specific child processes: `--provider openai` needs
`OPENAI_API_KEY`, `--provider tokenswitch` needs `NEOLINK_API_KEY`, and `--provider cascade` needs both
(or an appropriate `CLOCKBASE_MODEL_LAUNCHER` prefix). No private operator skill is required.
Use metadata you acquired for the study, a declared target and a new output directory:

```bash
python -m audit_system.agent_workflow_2026_09_10 --metadata SOURCE_METADATA.csv --study GSE... --target 'SOURCE-SUPPORTED COMPARISON' --provider cascade --output NEW_DESIGN_RUN
python -m audit_system.source_packet prepare --metadata SOURCE_METADATA.csv --study GSE... --output NEW_PACKET
python -m audit_system.evaluation_2026_09_10 --replay-comparison OWN_COMPARISON_ARCHIVE --output NEW_REPLAY_DIRECTORY
python -m audit_system.evaluation_2026_09_10 --freeze-reference-unavailable OWN_COMPARISON_ARCHIVE --study GSE... --output NEW_REFERENCE_DIRECTORY
```

`SOURCE_METADATA.csv` is a reader-supplied complete study table, not a promised bundle example.
The evaluation modes require a development comparison archive produced by `--compare-plan`, including its pinned
sources and completed run records. They cannot replay arbitrary verifier directories or missing development
archives. No example comparison plan/archive is shipped. `--freeze-reference-unavailable` supports only the
archive's declared shared-prerequisite failure; absence of review after successful execution is unsupported.

### Replication across independent projects

`intervention_meta_2026_09_29.py` pools the verified effects of each intervention across studies (Extended Data
Figure 22, Supplementary Table S19). `geo` reads the public GEO record of every series in an intervention with
results in more than one series and joins series that share a SuperSeries, a PubMed ID or a contact name and
institute into projects; `run` pools by project (or by series without `--projects`). The records deposit holds
the GEO records read on 2026-10-14 (`supporting/intervention-meta/geo-series-records/`), so `geo --offline` with
that directory as `--cache` rebuilds `projects.json` without network access:

```bash
python -m audit_system.intervention_meta_2026_09_29 geo --atlas V/NAME/atlas-v3/atlas_v3_combined.csv --classification V/NAME/public-export/classification --cache GEO_RECORDS --out NEW/projects.json --offline
python -m audit_system.intervention_meta_2026_09_29 run --atlas V/NAME/atlas-v3/atlas_v3_combined.csv --classification V/NAME/public-export/classification --projects NEW/projects.json --context in_vivo --out NEW/primary
```

The sensitivity runs drop `--projects` (series as the unit), drop `--context` (all contexts) or add
`--no-contact-link`. `concordance --meta NEW/primary/intervention_meta.csv --genage genage_models.csv --drugage
drugage.csv --out NEW/concordance` matches the pooled mouse effects to GenAge (build 21) and DrugAge (build 5); the HAGR
files are not redistributed, and `concordance/summary.json` records their SHA-256. `scripts/figures_2026_10_03/build_ed22_replication_2026_10_14.py` draws the figure from
`V/intervention-meta-2026-10-14/`.

## New-round design review

`source_packet.build_packet(study_id, samples, comparisons)` prepares a
study-design review packet without generating hashes. It keeps every supplied
sample separate, rejects duplicate or mixed-study sample IDs, and withholds
clock predictions, P values, effect sizes, scores and previous model labels.
Raw metadata and harmonized fields remain explicitly named. The caller supplies
the full study sample set; the function does not retrieve missing data.

`source_packet.parse_review(text)` requires complete JSON with `status`,
`findings` and `unresolved`. It does not salvage truncated responses or turn
model judgments into scientific approval. The caller must also reject responses
whose completion reason indicates truncation. This module does not execute
numerical reconstruction or independently verify a model's findings.

## Existing W0 index (legacy reference)

This package implements the W0 evidence index and deterministic claim decision
for the ClockBase independent audit. It records existing artifacts and audit
statuses without changing the source reports.

> **Note:** A successful ledger build confirms that selected report bytes were
> read and recorded. It does not approve a scientific result or manuscript claim.

## Interpret the records

The measurement adapter preserves `FAIL`, `BLOCKED`, `NEEDS_REVIEW`, and
`NOT_EVALUATED` as distinct evidence states. It does not convert a completed
software check into scientific approval.

The protected-repair adapter verifies every artifact digest declared by the
manifest. Repairer identity and independent verification authority remain
`BLOCKED` until separate actor and authority records establish them.

The claim decision returns `APPROVED` only when all required evidence and
checks pass the recorded authority policy. Failed checks, quarantined evidence,
simulation used as empirical evidence, unresolved material findings, source
identity changes, and self-verification prevent approval.

## Declared covariate designs

`source_packet.design_matrix_diagnostic(design, column_names, sample_ids, contrasts)`
checks a supplied finite design matrix and named coefficient contrasts. It returns
the SVD rank, residual degrees of freedom and each contrast's estimability, with
explicit tolerances. Redundant nuisance covariates can coexist with an estimable
condition contrast; full matrix rank is therefore reported separately. Missing
biological independence or measurement applicability remains `CANNOT VERIFY`.
This utility does not choose covariates, infer animal identity or approve a repair.

## Run the tests

Activate the environment installed in [README.md](../README.md#install), then run:

```bash
python -B -m unittest discover -s audit_system/tests -t .
```

Tests use software fixtures and skip named cases when their development artifacts are absent. They do not
estimate analytical error rates. The current release receipt is the `release_check.json` generated by the
[root release check](../README.md#release-check), rather than a dated development test count.
