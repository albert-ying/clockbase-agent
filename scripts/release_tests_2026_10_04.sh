#!/usr/bin/env bash
# Run the unit tests of a code release the way the release check does: release tree read-only, no network, no other
# data, TMPDIR on a private tmpfs. Without unittest arguments, discovers every test under audit_system/tests.
#   bash scripts/release_tests_2026_10_04.sh RELEASE_DIR VENV_DIR [unittest arguments...]
set -euo pipefail
RELEASE=$(realpath "$1"); VENV=$(realpath "$2"); shift 2
ARGS=("$@")
[ ${#ARGS[@]} -eq 0 ] && ARGS=(discover -s audit_system/tests -t .)
HOME_DIR=$(mktemp -d)
bwrap --unshare-all --die-with-parent --new-session \
      --ro-bind /usr /usr --ro-bind /etc /etc --symlink usr/bin /bin --symlink usr/lib /lib --symlink usr/lib64 /lib64 \
      --proc /proc --dev /dev --tmpfs /tmp \
      --ro-bind "$RELEASE" "$RELEASE" --ro-bind "$VENV" "$VENV" --bind "$HOME_DIR" "$HOME_DIR" \
      --setenv HOME "$HOME_DIR" --setenv TMPDIR /tmp --setenv PATH /usr/bin --chdir "$RELEASE" \
      "$VENV/bin/python" -B -m unittest "${ARGS[@]}"
