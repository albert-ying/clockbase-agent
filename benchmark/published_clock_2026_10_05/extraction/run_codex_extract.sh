#!/usr/bin/env bash
# Run one Codex extraction per paper input that has no output yet (resumable).
#   [EXTRACTOR=codex] [PROMPT=codex_prompt.md] bash run_codex_extract.sh DATA POSSIBLE_DIR JOBS MODEL EFFORT [INPUT ...]
# Writes POSSIBLE_DIR/$EXTRACTOR/<paper>.json (the final message), <paper>.log and <paper>.status (exit code;
# 'blocked' when the provider refused the content). Inputs default to every POSSIBLE_DIR/inputs/*.json.
set -uo pipefail
DATA=$1; POSSIBLE=$2; JOBS=$3; MODEL=$4; EFFORT=$5; shift 5
HERE=$(cd "$(dirname "$0")" && pwd)
OUT="$POSSIBLE/${EXTRACTOR:-codex}"; mkdir -p "$OUT"
PROMPT=${PROMPT:-codex_prompt.md}
export DATA OUT MODEL EFFORT HERE PROMPT

one() {
    local input=$1 name; name=$(basename "$input" .json)
    [[ -s "$OUT/$name.status" ]] && return 0
    local prompt; prompt=$(sed -e "s|{DATA}|$DATA|g" -e "s|{INPUT}|$input|g" "$HERE/$PROMPT")
    timeout --kill-after=60 3600 codex exec --skip-git-repo-check -s read-only -C "$DATA" -m "$MODEL" \
        -c model_reasoning_effort="\"$EFFORT\"" -o "$OUT/$name.json" - <<<"$prompt" > "$OUT/$name.log" 2>&1
    local code=$?
    if grep -q 'flagged for possible biological risk' "$OUT/$name.log"; then echo blocked > "$OUT/$name.status"
    else echo "$code" > "$OUT/$name.status"; fi
}
export -f one

if (($#)); then printf '%s\n' "$@"; else ls "$POSSIBLE"/inputs/*.json; fi | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
