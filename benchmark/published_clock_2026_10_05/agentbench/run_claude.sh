#!/usr/bin/env bash
# Run one Claude Code agent per prepared run directory that has no status yet (PROTOCOL.md of groupbench, item 22).
#   bash run_claude.sh RUNS CLAUDE_HOME JOBS [RUN_DIR ...]
# RUNS is the root written by prepare.py or ranked.py with --framework claude (RUNS/claude-<model>/<GSE>/TASK.md);
# CLAUDE_HOME is the isolated Claude Code configuration directory, holding only the sign-in credentials (no
# CLAUDE.md, memories, skills, plugins, hooks or MCP servers). Each run is confined by bubblewrap: the filesystem is
# read-only except the run directory, CLAUDE_HOME and a private /tmp, as Codex's workspace-write sandbox. Each run
# writes events.jsonl (stream-json events with tool calls and usage), stderr.log and status (the exit code).
set -uo pipefail
RUNS=$(realpath "$1"); HOME_DIR=$(realpath "$2"); JOBS=$3; shift 3
LIMIT=${AGENTBENCH_SECONDS:-1800}
CLAUDE=${CLAUDE_BIN:-claude}
SCRATCH=${AGENTBENCH_TMP:-/Entropy/kying/tmp/claude-bench-tmp}
export RUNS HOME_DIR LIMIT CLAUDE SCRATCH

one() {
    local dir=$1 model tmp
    [[ -s "$dir/status" ]] && return 0
    model=$(basename "$(dirname "$dir")"); model=${model#claude-}
    tmp="$SCRATCH/$model/$(basename "$dir")"; mkdir -p "$tmp"
    date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/started"
    timeout --kill-after=60 "$LIMIT" \
        bwrap --ro-bind / / --bind "$dir" "$dir" --bind "$HOME_DIR" "$HOME_DIR" --bind "$tmp" /tmp \
              --dev /dev --proc /proc --chdir "$dir" -- \
        env HOME="$HOME_DIR" CLAUDE_CONFIG_DIR="$HOME_DIR" TMPDIR=/tmp \
        "$CLAUDE" -p --model "$model" --output-format stream-json --verbose \
              --tools Bash Read Write Edit Glob Grep NotebookEdit TodoWrite \
              --strict-mcp-config --dangerously-skip-permissions --no-session-persistence \
        < "$dir/TASK.md" > "$dir/events.jsonl" 2> "$dir/stderr.log"
    echo $? > "$dir/status"
    date -u +%Y-%m-%dT%H:%M:%SZ > "$dir/finished"
}
export -f one

if (($#)); then printf '%s\n' "$@"; else ls -d "$RUNS"/claude-*/GSE*; fi | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {}
