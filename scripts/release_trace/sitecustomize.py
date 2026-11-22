"""Record every file a Python process opens for reading (rebuild-bundle tracing).

Active only when CLOCKBASE_TRACE_DIR is set and this directory is on PYTHONPATH; each process appends the absolute
paths it opens for reading to CLOCKBASE_TRACE_DIR/<pid>.txt as it opens them. Paths are written at once, not at exit,
because multiprocessing workers end with os._exit and never run exit handlers. Used by
scripts/build_release_2026_10_04.py to collect exactly the inputs a cache-only rebuild reads.
"""
import os
import sys

_TRACE_DIR = os.environ.get('CLOCKBASE_TRACE_DIR')
if _TRACE_DIR:
    _opened = set()
    _writing = False

    def _hook(event, args):
        global _writing
        if event != 'open' or not args or _writing:
            return
        path, mode = args[0], args[1] if len(args) > 1 else 'r'
        if isinstance(path, int) or path is None:
            return
        if isinstance(mode, str) and any(flag in mode for flag in ('w', 'a', 'x', '+')):
            return
        try:
            path = os.path.abspath(os.fsdecode(path))
        except (TypeError, ValueError):
            return
        if path in _opened:
            return
        _opened.add(path)
        _writing = True  # the trace file's own open must not be traced
        try:
            with open(os.path.join(_TRACE_DIR, f'{os.getpid()}.txt'), 'a') as handle:
                handle.write(path + '\n')
        finally:
            _writing = False

    sys.addaudithook(_hook)
