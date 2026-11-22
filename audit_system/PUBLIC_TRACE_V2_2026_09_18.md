# Public pipeline trace version 2

Version 1 source-review traces remain unchanged. Version 2 adds saved pipeline execution and numerical readback. The exporter verifies every completed stage receipt and artifact hash before reading results. It does not rerun molecular preprocessing or clocks, endorse a scientific claim, or infer an effect-size definition from an unlabeled value. This initial exporter requires all eight completed stages; incomplete pipeline runs need their own terminal receipt before export.

The root schema is `clockbase-audit-trace-v2`. All version 1 fields remain required with their existing meanings. These additional root objects are required:

```text
execution:
  stage_count: integer >= 0
  completed_stage_count: integer >= 0
  clocks_rerun_during_export: false
  stages: [{stage: string, status: PASS|FAIL|CANNOT VERIFY,
            elapsed_seconds: finite number >= 0, artifact_count: integer >= 0}]
data:
  sample_ids: unique GSM strings
  sample_count: integer == length(sample_ids)
  molecular_feature_count: integer >= 0
  molecular_sample_count: integer >= 0
  clock_names: unique strings
  preprocessing: [{name: string, value: string}]
  checks: [{name: string, status: PASS|FAIL|CANNOT VERIFY, scope: string, reason: string}]
analysis:
  scope: string
  multiplicity:
    method: string|null
    attempted_tests: integer >= 0|null
    denominator_status: PASS|FAIL|CANNOT VERIFY
  results:
    - result_id: unique string
      comparison_id: string
      clock: string
      case_ids: unique GSM strings
      control_ids: unique GSM strings, disjoint from case_ids
      original: numerical_record
      proposal: numerical_record|null
      proposal_groups: {case_ids: unique GSM strings, control_ids: unique GSM strings}|null
      selection:
        output: original|candidate|none|unrecorded
        status: PASS|FAIL|CANNOT VERIFY
        accepted_result: boolean
        repair_applied: boolean|null
        reason: string
      checks:
        - name: string
          status: PASS|FAIL|CANNOT VERIFY
          scope: string
          reported: finite number|null
          recomputed: finite number|null
          absolute_difference: finite number >= 0|null
          reason: string
```

`numerical_record` contains exactly `state`, `statistical_test`, `p_value`, `mean_difference`, `effect_size`, `effect_size_definition`, and `reason`. State is `COMPUTED`, `ABSTAINED`, `INCONSISTENT`, or `NOT_REPORTED`. The test and reason are strings. Numeric fields are finite numbers or null; P values must lie in [0,1]. The effect definition is a string or null. Null means unavailable, and zero must remain zero. Empty or missing estimates are never fabricated. `ABSTAINED` is distinct from a failed numerical result.

`proposal_groups` is null exactly when `proposal` is null. Otherwise its case and control identities are disjoint and drawn from the source roster. Original and proposed memberships are displayed separately, including when a repair changes grouping.

Individual source or arithmetic checks can pass while overall scientific status remains `CANNOT VERIFY`. The interface must label original statistics, repair proposals and selected output separately. Native clock-score differences and standardized effects require separate labels. An absent effect definition must display as unspecified. Software stage completion and response validity must not appear as scientific acceptance.

All downloadable files use the existing registered-artifact contract and SHA256 verification. API envelopes, credentials, contact metadata and absolute server paths are excluded. The public source snapshot, sample predictions, bounded molecular input, numerical checks and stage/hash manifest can be downloaded when explicitly registered. Large molecular inputs can be represented by their hash and source record without copying the matrix into the website bundle. The exporter records whether each check was repeated during export or read from a saved stage.
