#!/bin/bash
# Verifier v3 final chain: rebuild every atlas verdict, repair and ranked output from the saved
# recovery, replay, adjudication and design-reconstruction artefacts.
#
#   bash scripts/verifier_v3_final_chain_2026_09_30.sh V_ROOT OUT_NAME [--classify] [--design-audit]
#
# V_ROOT    verifier run root (audit_system.paths.VERIFIER_ROOT in the rebuild bundle)
# OUT_NAME  new directory name under V_ROOT for this chain (refuses to overwrite)
# --classify      also run the incremental dual-model classification (paid; needs TokenSwitch credentials)
# --design-audit  also run the dual-model design audit of ranking candidates (paid) into V_ROOT/design-audit-001
# The DNAm input-integrity screen (V_ROOT/input-integrity-dnam-001) and an existing design audit are used
# when present; the ranking summary records which context source it used.
set -euo pipefail
V=$1
OUT=$V/$2
shift 2
CLASSIFY=; AUDIT=
for flag in "$@"; do
  case $flag in
    --classify) CLASSIFY=1 ;;
    --design-audit) AUDIT=1 ;;
    *) echo "unknown flag $flag"; exit 1 ;;
  esac
done
INTEGRITY=$V/input-integrity-dnam-001
REPO=$(cd "$(dirname "$0")/.." && pwd)
VENV_PY="${CLOCKBASE_PYTHON:-$REPO/.venv-pipeline/bin/python}"
# Checked clock-output arm policy travels with the code deposit. An explicit
# override is retained; an unset or empty override selects the deposited policy.
CLOCK_OUTPUT_DECISIONS=${CLOCKBASE_CLOCK_OUTPUT_DECISIONS:-$REPO/audit_system/policy/clock_output_arm_decisions_2026-11-12.json}
[ -f "$CLOCK_OUTPUT_DECISIONS" ] || { echo "missing clock-output decisions: $CLOCK_OUTPUT_DECISIONS" >&2; exit 1; }
PY="$VENV_PY -B -m"
cd "$REPO"
STREAM=${STREAM_ROOT:-${CLOCKBASE_STREAM_EMPTY:-$("$VENV_PY" -B -c 'from audit_system import paths; print(paths.STREAM_EMPTY)')}}
STREAM_DETACHED=${STREAM_DETACHED_ROOT:-${CLOCKBASE_STREAM_DETACHED:-$("$VENV_PY" -B -c 'from audit_system import paths; print(paths.STREAM_DETACHED)')}}
# Optional credential prefix; shlex parsing preserves quoted arguments and paths with spaces.
API=()
if [ -n "$CLASSIFY$AUDIT" ]; then
  while IFS= read -r -d '' token; do API+=("$token"); done < <(
    "$VENV_PY" -B -c 'from audit_system import paths; import sys; [sys.stdout.write(token + "\0") for token in paths.MODEL_LAUNCHER]'
  )
fi
model_run() {
  local provider=$1
  shift
  local launcher=() token
  for token in ${API[@]+"${API[@]}"}; do launcher+=("${token//\{provider\}/$provider}"); done
  ${launcher[@]+"${launcher[@]}"} "$VENV_PY" -B -m "$@"
}
[ -e "$OUT" ] && { echo "refusing to overwrite $OUT"; exit 1; }
mkdir -p "$OUT"
for m in dnam hs_rna mm_rna; do
  $PY audit_system.membership_merge_2026_09_27 --modality $m --recovery $V/membership-$m-001 \
      --link $V/link-$m-001 --replay $V/replay-$m-001 --design $V/design-$m-001 --out $OUT/merged-$m
  checks=()
  [ $m = dnam ] && [ -f $INTEGRITY/samples.csv ] && checks=(--integrity $INTEGRITY)
  $PY audit_system.row_verdicts_2026_09_27 --modality $m --membership $OUT/merged-$m ${checks[@]+"${checks[@]}"} --out $OUT/verdicts-$m
  extra=()
  [ $m = dnam ] && extra=(--repaired-dnam $STREAM)
  [ $m = dnam ] && [ -d $STREAM_DETACHED/studies ] && extra+=(--repaired-dnam $STREAM_DETACHED)
  # The final chain uses the design-respecting estimand; age and time contrasts and audited individuals are read
  # from the classification and individuals-audit runs when they exist.
  design=(--design-estimand)
  design+=(--clock-output-decisions "$CLOCK_OUTPUT_DECISIONS")
  [ -f $V/classification-001/comparisons.jsonl ] && design+=(--classification $V/classification-001)
  [ -f $V/individuals-audit-001/comparisons.jsonl ] && design+=(--individuals $V/individuals-audit-001)
  $PY audit_system.row_repairs_2026_09_28 --modality $m --membership $OUT/merged-$m --verdicts $OUT/verdicts-$m \
      --replay $V/replay-$m-001 --adjudication $V/adjudication-$m-001 ${extra[@]+"${extra[@]}"} ${design[@]+"${design[@]}"} --out $OUT/repairs-$m
done
REP=(--repairs dnam=$OUT/repairs-dnam --repairs hs_rna=$OUT/repairs-hs_rna --repairs mm_rna=$OUT/repairs-mm_rna)
$PY audit_system.published_outputs_status_2026_09_29 ${REP[@]+"${REP[@]}"} --out $OUT/published-status
$PY audit_system.error_census_2026_09_29 ${REP[@]+"${REP[@]}"} --published $OUT/published-status --out $OUT/error-census
$PY audit_system.atlas_rebuild_2026_09_29 ${REP[@]+"${REP[@]}"} --out $OUT/atlas-v3
if [ -n "$CLASSIFY" ]; then
  AUDIT_MODELS=gp-6-as,cl-op-5 model_run tokenswitch audit_system.intervention_classification_2026_09_29 \
      --atlas $OUT/atlas-v3/atlas_v3_combined.csv --out $V/classification-001 ${REP[@]+"${REP[@]}"} --workers 12 --max-requests 7000
fi
if [ -n "$AUDIT" ]; then
  AUDIT_MODELS=openai:gpt-6-astra,openai:gpt-6-sol model_run openai audit_system.design_audit_2026_09_30 \
      --atlas $OUT/atlas-v3/atlas_v3_combined.csv --classification $V/classification-001 ${REP[@]+"${REP[@]}"} \
      --out $V/design-audit-001 --workers 12
fi
if [ -f $V/classification-001/comparisons.jsonl ]; then
  audit=()
  for d in $V/design-audit-0*; do [ -f $d/comparisons.jsonl ] && audit+=(--design-audit $d); done  # 001, then added runs
  $PY audit_system.comparison_ranking_2026_09_30 --atlas $OUT/atlas-v3/atlas_v3_combined.csv \
      --classification $V/classification-001 ${audit[@]+"${audit[@]}"} --out $OUT/comparison-ranking
fi
$PY audit_system.v3_public_export_2026_09_30 --final $OUT --verifier $V --out $OUT/public-export
echo "chain complete: $OUT"
