#!/usr/bin/env bash
# Run one Codex agent per prepared run directory that has no status yet (PROTOCOL.md of groupbench, item 16).
#   bash run_codex.sh RUNS CODEX_HOME JOBS [RUN_DIR ...]
# RUNS is the root written by prepare.py (RUNS/codex-<model>/<GSE>/TASK.md); CODEX_HOME is the isolated Codex home
# (no user instructions, plugins, memories or web search; the API key comes from the launcher's environment).
# Each run writes events.jsonl (Codex JSON events with commands and token usage), final_message.txt, stderr.log and
# status (the exit code). CLOCKBASE_MODEL_LAUNCHER supplies the command prefix that sets OPENAI_API_KEY.
set -uo pipefail
RUNS=$(realpath "$1"); HOME_DIR=$(realpath "$2"); JOBS=$3; shift 3
LIMIT=${AGENTBENCH_SECONDS:-1800}
export RUNS HOME_DIR LIMIT

one() {
    local dir=$1 model
    [[ -s "$dir/status" ]] && return 0
    model=$(basename "$(dirname "$dir")"); model=${model#codex-}
    date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/started"
    ${CLOCKBASE_MODEL_LAUNCHER:-} env CODEX_HOME="$HOME_DIR" timeout --kill-after=60 "$LIMIT" \
        codex exec --skip-git-repo-check --ephemeral --strict-config -s workspace-write -C "$dir" -m "$model" --json \
        -o "$dir/final_message.txt" - < "$dir/TASK.md" > "$dir/events.jsonl" 2> "$dir/stderr.log"
    echo $? > "$dir/status"
    date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/finished"
}
export -f one

if (($#)); then printf '%s\n' "$@"; else ls -d "$RUNS"/codex-*/GSE*; fi | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
