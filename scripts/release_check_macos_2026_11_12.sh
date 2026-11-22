#!/usr/bin/env bash
# Reproduce the deposited final run on macOS. Environment installation happens
# before isolation; rebuild, comparison, simulated self-test and tests use Seatbelt.
# Usage: bash SCRIPT RELEASE_DIR BUNDLE_DIR NEW_WORK_DIR REFERENCE_FINAL
# --prepare-only as a fifth argument prepares inputs/profile without running a
# sandbox or producing release_check.json. It is for preparation checks only.
set -euo pipefail
umask 077
[[ $# -eq 4 || ( $# -eq 5 && $5 == --prepare-only ) ]] || {
  echo "usage: $0 RELEASE BUNDLE NEW_WORK REFERENCE [--prepare-only]" >&2; exit 2;
}
[[ $(uname -s) == Darwin ]] || { echo "macOS is required" >&2; exit 1; }
[[ -x /usr/bin/sandbox-exec ]] || { echo "sandbox-exec is required" >&2; exit 1; }
command -v uv >/dev/null || { echo "uv is required" >&2; exit 1; }
# uv supplies a native Python 3.12 interpreter. No caller's data/credential
# environment is inherited by any sandboxed command.
BOOTSTRAP=$(uv python find --offline 3.12)
DEVELOPER_DIR=$(/usr/bin/xcode-select -p)
eval "$("$BOOTSTRAP" -B - "$1" "$2" "$3" "$4" <<'PY'
import shlex, sys
from pathlib import Path
release, bundle, work = (Path(p).resolve() for p in sys.argv[1:4])
reference = sys.argv[4]
if reference in ('', '.', '..') or Path(reference).name != reference:
    raise SystemExit('REFERENCE must be one directory name')
for a, b in ((release, work), (bundle, work), (release, bundle)):
    if a == b or a in b.parents or b in a.parents:
        raise SystemExit('release, bundle and work must be separate directory trees')
for name, value in zip(('RELEASE', 'BUNDLE', 'WORK', 'REFERENCE'), (release, bundle, work, reference)):
    print(f'{name}={shlex.quote(str(value))}')
PY
)"
[[ ! -e $WORK ]] || { echo "refusing to reuse $WORK" >&2; exit 1; }
[[ -f $RELEASE/requirements-verifier.txt && -f $RELEASE/scripts/clockbase_verifier.py ]] || {
  echo "incomplete code release: $RELEASE" >&2; exit 1;
}
[[ -f $BUNDLE/clockbase_data/verifier/$REFERENCE/public-export/manifest.json ]] || {
  echo "missing reference public export: $REFERENCE" >&2; exit 1;
}
mkdir -p "$WORK/home" "$WORK/logs" "$WORK/tmp" "$WORK/home/matplotlib" "$WORK/git-template"
REPORT=$WORK/release_check.json
trap 'echo "release check failed at line $LINENO; logs: $WORK/logs" >&2' ERR
step() { echo "== $1"; }

# Kill the complete process group on timeout, including subprocess workers.
# Perl is available on stock macOS when GNU timeout is absent.
LIMIT() {
  local seconds=$1; shift
  if command -v gtimeout >/dev/null; then
    gtimeout --kill-after=60 "$seconds" "$@"
  elif command -v timeout >/dev/null; then
    timeout --kill-after=60 "$seconds" "$@"
  else
    /usr/bin/perl -e '
      use POSIX qw(setpgid WNOHANG);
      my $limit = shift @ARGV;
      my $pid = fork(); defined $pid or die "fork: $!";
      if (!$pid) { setpgid(0, 0); exec @ARGV; die "exec: $!"; }
      eval { setpgid($pid, $pid) };
      $SIG{ALRM} = sub {
        kill "TERM", -$pid;
        $SIG{ALRM} = sub { kill "KILL", -$pid; waitpid($pid, 0); exit 124; };
        alarm 60; waitpid($pid, 0); exit 124;
      };
      alarm $limit; waitpid($pid, 0); alarm 0;
      exit(($? & 127) ? 128 + ($? & 127) : $? >> 8);
    ' "$seconds" "$@"
  fi
}

step "environment from release requirements-verifier.txt"
UV_CACHE_DIR="$WORK/uv-cache" uv venv --python "$BOOTSTRAP" "$WORK/venv" > "$WORK/logs/venv.log" 2>&1
UV_CACHE_DIR="$WORK/uv-cache" uv pip install --python "$WORK/venv/bin/python" -r "$RELEASE/requirements-verifier.txt" \
  --extra-index-url https://download.pytorch.org/whl/cpu --index-strategy unsafe-best-match >> "$WORK/logs/venv.log" 2>&1
PY=$WORK/venv/bin/python
step "hard-linked read-only data copy and explicit runtime profile"
"$PY" -B - "$RELEASE" "$BUNDLE" "$WORK" <<'PY'
import json, os, stat, subprocess, sys
from pathlib import Path
release, bundle, work = map(Path, sys.argv[1:])
source, target = bundle / 'clockbase_data', work / 'clockbase_data'
try:
    # No symlink may lead back to the author's data. A hard-link failure (for
    # example a different filesystem) stops rather than silently copying data.
    for folder, dirs, files in os.walk(source, followlinks=False):
        folder = Path(folder)
        dest = target / folder.relative_to(source)
        dest.mkdir(parents=True, exist_ok=True)
        for name in dirs + files:
            if (folder / name).is_symlink():
                raise ValueError(f'bundle contains a symlink: {folder / name}')
        for name in files:
            src = folder / name
            if not src.is_file():
                raise ValueError(f'not a regular bundle file: {src}')
            os.link(src, dest / name)
            # Permissions belong to the shared inode; the bundle becomes
            # read-only too. Its file content and timestamps are preserved.
            os.chmod(dest / name, stat.S_IMODE(src.stat().st_mode) & ~0o222)
    runtime = Path(sys.base_prefix).resolve()
    if runtime == Path('/') or runtime == Path.home():
        raise ValueError(f'unsafe Python runtime root: {runtime}')
    subpaths = [release, work, runtime, Path('/bin'), Path('/usr/bin'),
                Path('/usr/lib'), Path('/usr/libexec'),
                Path('/System/Library/Frameworks'), Path('/System/Library/PrivateFrameworks'),
                Path('/System/Library/Perl'),
                Path('/usr/share/locale'), Path('/usr/share/zoneinfo'),
                Path('/Library/Developer/CommandLineTools/usr/bin'),
                Path('/Library/Developer/CommandLineTools/usr/lib'),
                Path('/Library/Developer/CommandLineTools/usr/libexec/git-core'),
                Path('/Library/Developer/CommandLineTools/usr/share/git-core'),
                Path('/private/var/db/timezone')]
    # /usr/bin/git is an Apple shim. Its selected installation may be Xcode
    # rather than Command Line Tools; allow only that installation's Git runtime.
    developer = None
    if sys.platform == 'darwin':
        developer = Path(subprocess.check_output(['/usr/bin/xcode-select', '-p'], text=True).strip()).resolve()
        subpaths += [developer / 'usr/bin', developer / 'usr/lib',
                     developer / 'usr/libexec/git-core', developer / 'usr/share/git-core']
    literals = ['/dev/null', '/dev/urandom', '/dev/random', '/dev/fd',
                '/private/etc/localtime', '/private/etc/passwd', '/private/etc/group',
                '/private/etc/hosts', '/private/var/run/utmpx',
                '/private/var/db/xcode_select_link',
                '/Library/Preferences/com.apple.dt.Xcode.plist',
                '/System/Library/CoreServices/SystemVersion.plist',
                # dyld reads the root directory itself; /bin/sh resolves through /private/var/select/sh.
                '/', '/private/var/select', '/private/var/select/sh']
    if developer and developer.parent.name == 'Contents':
        literals += [str(developer.parent / 'Info.plist'), str(developer.parent / 'version.plist')]
    # The environment's python link can pass through uv's minor-version alias directory (a symlink to the
    # versioned runtime); resolving the link needs read access to each link on the way, not to their siblings.
    link = work / 'venv' / 'bin' / 'python'
    hops = 0
    while link.is_symlink() and hops < 8:
        target = Path(os.readlink(link))
        target = target if target.is_absolute() else (link.parent / target)
        for part in [target] + list(target.parents):
            if part.is_symlink():
                literals.append(str(part))
        link, hops = target, hops + 1
    # Metadata access to ancestors permits traversal without reading siblings.
    ancestors = {str(p) for root in subpaths for p in root.parents}
    quote = lambda s: json.dumps(str(s))
    profile = ['(version 1)', '(deny default)', '(deny network*)',
               '(allow process*)', '(allow signal (target same-sandbox))',
               '(allow sysctl-read)', '(allow ipc-posix*)',
               '(allow mach-lookup (global-name "com.apple.system.logger") '
               '(global-name "com.apple.cfprefsd.agent") (global-name "com.apple.cfprefsd.daemon") '
               '(global-name "com.apple.system.opendirectoryd.libinfo") '
               '(global-name "com.apple.FontObjectsServer") (global-name "com.apple.fonts"))']
    profile += [f'(allow file-read* (subpath {quote(p)}))' for p in subpaths]
    profile += [f'(allow file-read* (literal {quote(p)}))' for p in literals]
    profile += [f'(allow file-read-metadata (literal {quote(p)}))' for p in sorted(ancestors)]
    # Loading compiled extension modules and the runtime library maps them executable.
    executable_roots = [runtime, work, release, Path('/usr/lib'), Path('/System'),
                        Path('/Library/Developer/CommandLineTools/usr/lib')]
    if developer:
        executable_roots += [developer / 'usr/bin', developer / 'usr/lib', developer / 'usr/libexec/git-core']
    profile += [f'(allow file-map-executable (subpath {quote(p)}))' for p in executable_roots]
    # /dev/null is a character-device sink used by subprocess.DEVNULL in the
    # existing tests. It cannot store data outside WORK.
    profile += [f'(allow file-write* (subpath {quote(work)}))']
    profile += ['(allow file-write-data (literal "/dev/null"))']
    (work / 'sandbox.sb').write_text('\n'.join(profile) + '\n')
    (work / 'sandbox_paths.json').write_text(json.dumps({
        'read_subpaths': list(map(str, subpaths)), 'read_literals': literals,
        'executable_subpaths': list(map(str, executable_roots)),
        'metadata_ancestors': sorted(ancestors), 'write_subpaths': [str(work)],
        'private_tmp': str(work / 'tmp'), 'write_device_sink': '/dev/null',
        'network': 'deny network*'}, indent=2))
except (OSError, ValueError) as error:
    raise SystemExit(f'Cannot prepare isolated data/profile: {error}') from error
PY

if [[ ${5:-} == --prepare-only ]]; then
  echo "Prepared $WORK/sandbox.sb; sandbox untested; no release_check.json produced"
  exit 0
fi
SANDBOX=(/usr/bin/sandbox-exec -f "$WORK/sandbox.sb" /usr/bin/env -i
  HOME="$WORK/home" TMPDIR="$WORK/tmp/" PATH=/Library/Developer/CommandLineTools/usr/bin:/usr/bin:/bin:/usr/sbin:/sbin
  LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 PYTHONDONTWRITEBYTECODE=1
  DEVELOPER_DIR="$DEVELOPER_DIR" GIT_CONFIG_NOSYSTEM=1 GIT_CONFIG_GLOBAL=/dev/null
  GIT_TEMPLATE_DIR="$WORK/git-template"
  MPLCONFIGDIR="$WORK/home/matplotlib" XDG_CACHE_HOME="$WORK/home/.cache"
  CLOCKBASE_DATA="$WORK/clockbase_data" CLOCKBASE_PYTHON="$PY"
  AUDIT_CACHE_ONLY=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
  DEVELOPER_DIR=/Library/Developer/CommandLineTools)
cd "$RELEASE"
step "sandbox network and filesystem probes"
# Only EPERM/EACCES establishes sandbox enforcement. DNS failure, a refused
# connection or a timeout does not count as a successful isolation check.
LIMIT 30 "${SANDBOX[@]}" "$PY" -B - "$WORK" "$RELEASE" "$BUNDLE" "$REFERENCE" <<'PY' > "$WORK/logs/sandbox-probe.log" 2>&1
import errno, json, os, socket, sys
from pathlib import Path
work, release, bundle = map(Path, sys.argv[1:4])
reference = sys.argv[4]
denied = {errno.EPERM, errno.EACCES}
with socket.socket() as server:
    server.settimeout(2)
    try:
        server.bind(('127.0.0.1', 0))
    except OSError as error:
        if error.errno not in denied:
            raise
    else:
        raise SystemExit('network probe succeeded: local bind was allowed')
with socket.socket() as client:
    client.settimeout(2)
    try:
        client.connect(('1.1.1.1', 443))
    except OSError as error:
        if error.errno not in denied:
            raise
    else:
        raise SystemExit('network probe succeeded: outbound connection was allowed')
# This file was checked before isolation, so ENOENT cannot masquerade as a
# successful denial. Read the original bundle, outside its hard-linked copy.
for path in (bundle / 'clockbase_data/verifier' / reference / 'public-export/manifest.json',):
    try:
        path.read_bytes()
    except OSError as error:
        if error.errno not in denied:
            raise
    else:
        raise SystemExit(f'outside read probe succeeded: {path}')
try:
    (release / '.clockbase-write-probe').write_text('must be denied')
except OSError as error:
    if error.errno not in denied:
        raise
else:
    raise SystemExit('release write probe succeeded')
(work / 'tmp/probe').write_text('private temporary directory works')
print(json.dumps({'network_denied': True, 'outside_reads_denied': True, 'release_writes_denied': True}))
PY

run_status=0
step "rebuild (cache-only, sandboxed)"
LIMIT 14400 "${SANDBOX[@]}" "$PY" -B scripts/clockbase_verifier.py run --final final-release-check --models cache-only \
  > "$WORK/logs/run.log" 2>&1 || run_status=$?
verify_status=1
step "verify against $REFERENCE"
if [[ $run_status -eq 0 ]]; then
  LIMIT 3600 "${SANDBOX[@]}" "$PY" -B scripts/clockbase_verifier.py verify --run "$WORK/clockbase_data/verifier/final-release-check" \
    --reference "$WORK/clockbase_data/verifier/$REFERENCE" > "$WORK/logs/verify.log" 2>&1 && verify_status=0 || verify_status=$?
fi
selftest_status=0
step "planted-defect self-test (simulation, sandboxed)"
LIMIT 3600 "${SANDBOX[@]}" "$PY" -B scripts/clockbase_verifier.py selftest --root "$WORK/clockbase_data/verifier/selftest-release-check" \
  --reference "$WORK/clockbase_data/verifier/$REFERENCE" --per-modality 160 > "$WORK/logs/selftest.log" 2>&1 || selftest_status=$?
test_status=0
step "unit tests (sandboxed)"
LIMIT 7200 "${SANDBOX[@]}" "$PY" -B -m unittest discover -s audit_system/tests -t . > "$WORK/logs/tests.log" 2>&1 || test_status=$?

"$PY" -B - "$REPORT" "$RELEASE" "$BUNDLE" "$REFERENCE" "$run_status" "$verify_status" "$test_status" "$WORK" "$selftest_status" <<'PY'
import csv, hashlib, json, re, sys, time
from pathlib import Path
report, release, bundle, reference, run_status, verify_status, test_status, work, selftest_status = sys.argv[1:]
work = Path(work)
tests = (work / 'logs/tests.log').read_text(errors='replace')
ran = re.search(r'Ran (\d+) tests', tests)
new, ref = (work / 'clockbase_data/verifier' / name for name in ('final-release-check', reference))
verify = new / 'verify.json'
selftest = work / 'clockbase_data/verifier/selftest-release-check/summary_SIMULATED.json'
out = {'created': time.strftime('%Y-%m-%dT%H:%M:%S%z'), 'release': release, 'bundle': bundle, 'reference': reference,
       'sandbox': 'sandbox-exec: deny network*; default-deny reads except release, work and explicitly listed '
                  'system/Python runtime locations in sandbox_paths.json; stored file writes only in work, including private tmp; '
                  '/dev/null character-device sink allowed for subprocesses; '
                  'clean environment; network and filesystem denial probes passed',
       'rebuild_exit': int(run_status), 'verify_exit': int(verify_status), 'tests_exit': int(test_status),
       'tests_ran': int(ran.group(1)) if ran else None, 'tests_tail': tests.strip().splitlines()[-3:],
       'verify': json.loads(verify.read_text()) if verify.exists() else None,
       'selftest_exit': int(selftest_status), 'selftest': json.loads(selftest.read_text()) if selftest.exists() else None}
# The manuscript says identical values and byte-identical study records.
# Also check exact row values and complete ranking CSV content rather than
# relying on the verifier's 1e-9 comparison tolerance and rounded ranking keys.
exact = {}
try:
    def exact_rows(path):
        with path.open() as handle:
            rows = [json.loads(line) for line in handle]
        return {r['row_ordinal']: (r['final_class'], r.get('recomputed')) for r in rows}
    for modality in ('dnam', 'hs_rna', 'mm_rna'):
        name = f'repairs-{modality}/repairs.jsonl'
        exact[f'rows-{modality}'] = exact_rows(new / name) == exact_rows(ref / name)
    name = 'comparison-ranking/ranked_comparisons_all.csv'
    with (new / name).open() as a, (ref / name).open() as b:
        exact['ranking_all_columns'] = list(csv.DictReader(a)) == list(csv.DictReader(b))
    def study_hashes(root):
        return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (root / 'public-export/studies').rglob('*.json')}
    a, b = study_hashes(new), study_hashes(ref)
    exact['study_bytes'] = bool(b) and a == b
except (OSError, ValueError, KeyError) as error:
    exact['error'] = str(error)
out['exact_content'] = exact
out['sandbox_paths'] = json.loads((work / 'sandbox_paths.json').read_text())
simulation = out['selftest'] or {}
out['selftest_checks_passed'] = (
    simulation.get('simulation') is True and simulation.get('planted_rows', 0) > 0
    and simulation.get('planted_rows') == simulation.get('caught_any_route')
    and simulation.get('clean_rows', 0) > 0
    and simulation.get('clean_rows') == simulation.get('clean_verified_unchanged')
    and bool(simulation.get('detached_simulation', {}).get('by_treatment'))
    and all(set(counts) == {'EXPECTED'} for counts in
            simulation.get('detached_simulation', {}).get('by_treatment', {}).values()))
out['passed'] = (all(out[k] == 0 for k in ('rebuild_exit', 'verify_exit', 'tests_exit', 'selftest_exit'))
                 and out['verify'] is not None and out['verify'].get('reproduced') is True
                 and out['selftest_checks_passed']
                 and bool(out['tests_ran']) and len(exact) == 5 and all(v is True for v in exact.values()))
Path(report).write_text(json.dumps(out, indent=2, default=str) + '\n')
print(json.dumps({k: v for k, v in out.items() if k not in ('verify', 'selftest', 'sandbox_paths')}, indent=1))
raise SystemExit(0 if out['passed'] else 1)
PY
