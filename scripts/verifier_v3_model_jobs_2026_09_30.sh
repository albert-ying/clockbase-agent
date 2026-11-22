#!/bin/bash
# Paid dual-model jobs that wait for NeoLink quota, then a chain rebuild that uses their results.
#
#   bash scripts/verifier_v3_model_jobs_2026_09_30.sh V_ROOT FINAL_NAME NEW_FINAL_NAME
#
# V_ROOT          verifier run root (e.g. clockbase_data/verifier)
# FINAL_NAME      chain directory whose atlas and repairs define the candidates (e.g. final-2026-09-30c)
# NEW_FINAL_NAME  chain directory to build afterwards (refuses to overwrite)
# Every request goes through the ledger of its output directory, so a rerun re-sends only failed calls.
set -euo pipefail
V=$1
F=$V/$2
NEW=$3
REPO=$(cd "$(dirname "$0")/.." && pwd)
PY=("$REPO/.venv-pipeline/bin/python" -B -m)
REP=(--repairs dnam=$F/repairs-dnam --repairs hs_rna=$F/repairs-hs_rna --repairs mm_rna=$F/repairs-mm_rna)
cd "$REPO"
# CLOCKBASE_MODEL_LAUNCHER is an optional credential command prefix; {provider} selects the provider.
mapfile -d '' -t API < <("${PY[0]}" -B -c 'import sys; from audit_system import paths; sys.stdout.buffer.write(b"".join(token.replace("{provider}", "tokenswitch").encode() + b"\0" for token in paths.MODEL_LAUNCHER))')
# 1. Comparisons whose classification failed (402) or was truncated: smaller batches, larger token limit.
"${API[@]}" "${PY[@]}" audit_system.intervention_classification_2026_09_29 --atlas $F/atlas-v3/atlas_v3_combined.csv \
    --out $V/classification-001 "${REP[@]}" --workers 12 --batch 8 --max-tokens 10000
# 2. Design audit of every ranking candidate (orientation, specimen, tumour, confounds).
"${API[@]}" "${PY[@]}" audit_system.design_audit_2026_09_30 --atlas $F/atlas-v3/atlas_v3_combined.csv \
    --classification $V/classification-001 "${REP[@]}" --out $V/design-audit-001 --workers 12
# 3. Orientation: a random sample of heuristic-only rows (validation), then the rows waiting for adjudication.
for m in dnam hs_rna mm_rna; do
  "${API[@]}" "${PY[@]}" audit_system.label_adjudication_2026_09_27 --modality $m --membership $F/merged-$m \
      --out $V/adjudication-$m-001 --select supported_sample --limit 300 --verdicts $F/repairs-$m --workers 8
  "${API[@]}" "${PY[@]}" audit_system.label_adjudication_2026_09_27 --modality $m --membership $F/merged-$m \
      --out $V/adjudication-$m-001 --select pending --verdicts $F/repairs-$m --workers 8
done
# 4. Rebuild every verdict, repair, ranking and export with the new decisions and the design audit.
bash "$REPO/scripts/verifier_v3_final_chain_2026_09_30.sh" $V $NEW
echo "model jobs and chain complete: $V/$NEW"
