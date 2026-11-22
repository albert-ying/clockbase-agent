#!/usr/bin/env bash
# Rebuild the verifier's final run from the two deposits alone (the code release and the rebuild bundle) and compare it
# row by row with the reference run shipped in the bundle; then run the planted-defect self-test (a simulation) against
# the reference and the unit tests. Everything runs in a bubblewrap sandbox without network and without this server's
# data trees, so it fails on any input that the deposits do not contain.
#
#   bash scripts/release_check_2026_10_04.sh RELEASE_DIR BUNDLE_DIR WORK_DIR REFERENCE_FINAL
#
# RELEASE_DIR  the curated source tree (scripts/build_release_2026_10_04.py code)
# BUNDLE_DIR   directory containing clockbase_data/ (scripts/build_release_2026_10_04.py bundle)
# WORK_DIR     new directory for the environment, a writable copy of the data (hard links, inputs read-only) and the report
# REFERENCE_FINAL  name of the reference final run inside clockbase_data/verifier/
set -euo pipefail
RELEASE=$(realpath "$1"); BUNDLE=$(realpath "$2"); WORK=$(realpath -m "$3"); REFERENCE=$4
[ -e "$WORK" ] && { echo "refusing to reuse $WORK"; exit 1; }
mkdir -p "$WORK/home" "$WORK/logs"
REPORT=$WORK/release_check.json
step() { echo "== $1"; }

step "environment from $RELEASE/requirements-verifier.txt"
uv venv --python 3.12 "$WORK/venv" > "$WORK/logs/venv.log" 2>&1
UV_CACHE_DIR=${UV_CACHE_DIR:-$WORK/uv-cache} uv pip install --python "$WORK/venv/bin/python" -r "$RELEASE/requirements-verifier.txt" \
    --extra-index-url https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match >> "$WORK/logs/venv.log" 2>&1

step "writable data root: hard links to the bundle, every input read-only"
cp -al "$BUNDLE/clockbase_data" "$WORK/clockbase_data"
find "$WORK/clockbase_data" -type f -exec chmod a-w {} +

SANDBOX=(bwrap --unshare-all --die-with-parent --new-session
         --ro-bind /usr /usr --ro-bind /etc /etc --symlink usr/bin /bin --symlink usr/lib /lib --symlink usr/lib64 /lib64
         --proc /proc --dev /dev --tmpfs /tmp
         --ro-bind "$RELEASE" "$RELEASE" --bind "$WORK" "$WORK"
         --setenv HOME "$WORK/home" --setenv TMPDIR /tmp --setenv PATH /usr/bin
         --setenv CLOCKBASE_DATA "$WORK/clockbase_data" --setenv CLOCKBASE_PYTHON "$WORK/venv/bin/python"
         --chdir "$RELEASE")
PY="$WORK/venv/bin/python"
# Every step has a time limit, so a hung step fails the check instead of blocking it.
LIMIT() { local limit=$1; shift; timeout --kill-after=60 "$limit" "$@"; }

status=0
step "rebuild (cache-only, sandboxed)"
LIMIT 4h "${SANDBOX[@]}" "$PY" -B scripts/clockbase_verifier.py run --final final-release-check --models cache-only \
    > "$WORK/logs/run.log" 2>&1 || status=$?
run_status=$status

step "verify against $REFERENCE"
verify_status=1
if [ $run_status -eq 0 ]; then
  LIMIT 1h "${SANDBOX[@]}" "$PY" -B scripts/clockbase_verifier.py verify --run "$WORK/clockbase_data/verifier/final-release-check" \
      --reference "$WORK/clockbase_data/verifier/$REFERENCE" > "$WORK/logs/verify.log" 2>&1 && verify_status=0 || verify_status=$?
fi

step "planted-defect self-test against $REFERENCE (simulation, sandboxed)"
selftest_status=0
LIMIT 1h "${SANDBOX[@]}" "$PY" -B scripts/clockbase_verifier.py selftest --root "$WORK/clockbase_data/verifier/selftest-release-check" \
    --reference "$WORK/clockbase_data/verifier/$REFERENCE" --per-modality 160 > "$WORK/logs/selftest.log" 2>&1 || selftest_status=$?

step "unit tests (sandboxed)"
test_status=0
LIMIT 2h "${SANDBOX[@]}" "$PY" -B -m unittest discover -s audit_system/tests -t . > "$WORK/logs/tests.log" 2>&1 || test_status=$?

"$PY" - "$REPORT" "$RELEASE" "$BUNDLE" "$REFERENCE" "$run_status" "$verify_status" "$test_status" "$WORK" "$selftest_status" <<'EOF'
import json, re, sys, time
from pathlib import Path
report, release, bundle, reference, run_status, verify_status, test_status, work, selftest_status = sys.argv[1:]
tests = Path(work, 'logs/tests.log').read_text(errors='ignore')
ran = re.search(r'Ran (\d+) tests', tests)
verify = Path(work, 'clockbase_data/verifier/final-release-check/verify.json')
selftest = Path(work, 'clockbase_data/verifier/selftest-release-check/summary_SIMULATED.json')
out = {'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'release': release, 'bundle': bundle, 'reference': reference,
       'sandbox': 'bubblewrap: no network; /usr and /etc read-only; release read-only; only the work directory writable; '
                  'no other filesystem visible',
       'rebuild_exit': int(run_status), 'verify_exit': int(verify_status), 'tests_exit': int(test_status),
       'tests_ran': int(ran.group(1)) if ran else None, 'tests_tail': tests.strip().splitlines()[-3:],
       'verify': json.loads(verify.read_text()) if verify.exists() else None,
       'selftest_exit': int(selftest_status),
       'selftest': json.loads(selftest.read_text()) if selftest.exists() else None}
out['passed'] = (out['rebuild_exit'] == 0 and out['verify_exit'] == 0 and out['tests_exit'] == 0
                 and out['selftest_exit'] == 0)
Path(report).write_text(json.dumps(out, indent=2, default=str))
print(json.dumps({k: v for k, v in out.items() if k not in ('verify', 'selftest')}, indent=1))
EOF
