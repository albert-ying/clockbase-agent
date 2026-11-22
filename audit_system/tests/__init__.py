"""Tests for the ClockBase independent audit system."""

from pathlib import Path
import unittest

from audit_system.paths import REPO


def require_development_artifact(path):
    """Skip a development regression only when its source artifact is absent."""
    path = Path(path)
    if not path.exists():
        try:
            name = path.relative_to(REPO).as_posix()
        except ValueError:
            name = str(path)
        raise unittest.SkipTest('development artifact not in the code release: ' + name)
    return path
