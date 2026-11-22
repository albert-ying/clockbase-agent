"""Adjudication decision files: one file per run, read back in run order so later runs override earlier ones."""
import re
import time
from pathlib import Path


def decision_files(folder):
    """decisions-*.jsonl sorted by the run timestamp suffix (YYYYmmddTHHMMSS), else by modification time."""
    def key(path):
        found = re.search(r'-(\d{8}T\d{6})\.jsonl$', path.name)
        return found.group(1) if found else time.strftime('%Y%m%dT%H%M%S', time.localtime(path.stat().st_mtime))
    return sorted(Path(folder).glob('decisions-*.jsonl'), key=key)
