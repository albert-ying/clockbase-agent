#!/usr/bin/env bash
# Run one Biomni A1 session per prepared run directory that has no status yet (PROTOCOL.md of groupbench, item 24).
#   bash run_biomni.sh RUNS BIOMNI_VENV JOBS [RUN_DIR ...]
# RUNS is the root written by prepare.py or ranked.py with --framework biomni (RUNS/biomni-<model>/<GSE>/TASK.md);
# BIOMNI_VENV is the virtual environment with Biomni installed from its public repository. The credential launcher
# (CLOCKBASE_MODEL_LAUNCHER) supplies OPENAI_API_KEY. Each run is confined
# by bubblewrap: the filesystem is read-only except the run directory and a private /tmp. Each run writes events.jsonl
# (biomni_runner.py: messages, executed blocks, token usage), stdout.log, stderr.log and status (the exit code).
set -uo pipefail
RUNS=$(realpath "$1"); VENV=$(realpath "$2"); JOBS=$3; shift 3
LIMIT=${AGENTBENCH_SECONDS:-1800}
RUNNER=$(realpath "$(dirname "$0")/biomni_runner.py")
SCRATCH=${AGENTBENCH_TMP:-/Entropy/kying/tmp/biomni-bench-tmp}
LAUNCHER=${CLOCKBASE_MODEL_LAUNCHER:?set CLOCKBASE_MODEL_LAUNCHER to the command prefix that supplies OPENAI_API_KEY}
export RUNS VENV LIMIT RUNNER SCRATCH LAUNCHER

one() {
    local dir=$1 model tmp
    [[ -s "$dir/status" ]] && return 0
    model=$(basename "$(dirname "$dir")"); model=${model#biomni-}
    tmp="$SCRATCH/$model/$(basename "$dir")"; mkdir -p "$tmp/home"
    date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/started"
    # shellcheck disable=SC2086  # the launcher is a command prefix
    $LAUNCHER timeout --kill-after=60 "$LIMIT" \
        bwrap --ro-bind / / --bind "$dir" "$dir" --bind "$tmp" /tmp --dev /dev --proc /proc --chdir "$dir" -- \
        env HOME=/tmp/home TMPDIR=/tmp MPLCONFIGDIR=/tmp/home/mpl \
        "$VENV/bin/python" -I -B "$RUNNER" "$model" > "$dir/stdout.log" 2> "$dir/stderr.log"
    echo $? > "$dir/status"
    date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/finished"
}
export -f one

if (($#)); then printf '%s\n' "$@"; else ls -d "$RUNS"/biomni-*/GSE*; fi | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
