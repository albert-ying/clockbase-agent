#!/usr/bin/env bash
# Prepare the portable source-to-audit profile without accessing credentials.
set -euo pipefail
TASK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$TASK_ROOT"
for dependency in curl Rscript bwrap; do
  command -v "$dependency" >/dev/null || { echo "Missing system dependency: $dependency (Debian/Ubuntu: apt install curl r-base-core bubblewrap)" >&2; exit 1; }
done
TASK_PYTHON="${CLOCKBASE_PYTHON:-python3.12}"
if command -v uv >/dev/null; then
  if [[ ! -d .venv-pipeline ]]; then uv venv --python "$TASK_PYTHON" .venv-pipeline; fi
  uv pip install --python .venv-pipeline/bin/python -r pipeline_requirements_2026_09_10.txt
else
  if [[ ! -d .venv-pipeline ]]; then "$TASK_PYTHON" -m venv .venv-pipeline; fi
  .venv-pipeline/bin/python -m pip install -r pipeline_requirements_2026_09_10.txt
fi
.venv-pipeline/bin/python -c 'from src.pipeline_inputs_2026_09_10 import checked_assets; checked_assets(); print("Model and preprocessing payloads verified")'
echo 'Setup complete. Set OPENAI_API_KEY, then use PYTHON=.venv-pipeline/bin/python bash run_pipeline.sh --help'
