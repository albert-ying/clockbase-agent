#!/bin/bash
# Highest-yield model jobs through the OpenAI API, then a chain rebuild.
# Both judges are OpenAI models (GPT-6 Astra and GPT-6 Sol); only agreement decides, as before.
#
#   bash scripts/verifier_v3_openai_jobs_2026_10_01.sh V_ROOT FINAL_NAME NEW_FINAL_NAME
set -euo pipefail
V=$1
F=$V/$2
NEW=$3
REPO=$(cd "$(dirname "$0")/.." && pwd)
PY=("$REPO/.venv-pipeline/bin/python" -B -m)
export AUDIT_MODELS=openai:gpt-6-astra,openai:gpt-6-sol
REP=(--repairs dnam=$F/repairs-dnam --repairs hs_rna=$F/repairs-hs_rna --repairs mm_rna=$F/repairs-mm_rna)
cd "$REPO"
# CLOCKBASE_MODEL_LAUNCHER is an optional credential command prefix; {provider} selects the provider.
mapfile -d '' -t API < <("${PY[0]}" -B -c 'import sys; from audit_system import paths; sys.stdout.buffer.write(b"".join(token.replace("{provider}", "openai").encode() + b"\0" for token in paths.MODEL_LAUNCHER))')
# 1. Design audit of every ranking candidate (orientation, specimen, tumour material, confounds).
"${API[@]}" "${PY[@]}" audit_system.design_audit_2026_09_30 --atlas $F/atlas-v3/atlas_v3_combined.csv \
    --classification $V/classification-001 "${REP[@]}" --out $V/design-audit-001 --workers 12
# 2. Validation sample: 300 usable rows per modality whose orientation rests on the label-token heuristic alone.
for m in dnam hs_rna mm_rna; do
  "${API[@]}" "${PY[@]}" audit_system.label_adjudication_2026_09_27 --modality $m --membership $F/merged-$m \
      --out $V/adjudication-$m-001 --select supported_sample --limit 300 --verdicts $F/repairs-$m --workers 8
done
# 3. Rebuild with the audit (ranking) and the validation decisions (swaps where both judges confirm an inversion).
bash "$REPO/scripts/verifier_v3_final_chain_2026_09_30.sh" $V $NEW
echo "openai jobs and chain complete: $V/$NEW"
