"""Where the verifier reads its inputs and writes its runs.

Every location can be set by an environment variable. Unset locations resolve under CLOCKBASE_DATA (default
<repository>/clockbase_data), which has the layout of the deposited rebuild bundle:

  clockbase_data/
    GEO_data/                      per-sample clock predictions, sample metadata, DNAm matrices and QC tables
    verifier/                      verifier runs (membership-*, replay-*, link-*, design-*, adjudication-*, ...)
    dnam-stream-repair-empty/      studies re-derived from GEO because their historical DNAm matrix was empty
    dnam-stream-repair-detached/   studies re-derived because their historical profiles were detached from probes
    corrections/                   the three atlas correction tables (the verified record)
    inputs/                        the empty-matrix census and the published items checked against the atlas
    agent-traces/                  the saved CodingAgent traces and per-study outputs (input of the trace replay)

  CLOCKBASE_DATA, CLOCKBASE_GEO, CLOCKBASE_VERIFIER_ROOT, CLOCKBASE_STREAM_EMPTY, CLOCKBASE_STREAM_DETACHED,
  CLOCKBASE_CORRECTIONS, CLOCKBASE_INPUTS, CLOCKBASE_TRACES, CLOCKBASE_METHYLOME_REFERENCE,
  CLOCKBASE_REPLAY_PYTHON   interpreter of the sandboxed trace replay (default: the running interpreter),
  CLOCKBASE_MODEL_LAUNCHER  optional command prefix that supplies API credentials to model stages
                            (default: none; verifier model stages read OPENAI_API_KEY for openai: models, or
                            TOKENSWITCH_BASE_URL and TOKENSWITCH_API_KEY for the NeoLink pair; the standalone
                            agent workflow reads OPENAI_API_KEY or NEOLINK_API_KEY).

`portable(path)` writes a path under one of these roots as 'root:relative/path', so exported records do not
depend on where a reader unpacks the bundle.
"""
import os
import shlex
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _env(name, default):
    value = os.environ.get(name)
    return Path(value) if value else Path(default)


DATA = _env('CLOCKBASE_DATA', REPO / 'clockbase_data')
GEO = _env('CLOCKBASE_GEO', DATA / 'GEO_data')
VERIFIER_ROOT = _env('CLOCKBASE_VERIFIER_ROOT', DATA / 'verifier')
STREAM_EMPTY = _env('CLOCKBASE_STREAM_EMPTY', DATA / 'dnam-stream-repair-empty')
STREAM_DETACHED = _env('CLOCKBASE_STREAM_DETACHED', DATA / 'dnam-stream-repair-detached')
CORRECTIONS = _env('CLOCKBASE_CORRECTIONS', DATA / 'corrections')
INPUTS = _env('CLOCKBASE_INPUTS', DATA / 'inputs')
TRACES = _env('CLOCKBASE_TRACES', DATA / 'agent-traces')
METHYLOME_REFERENCE = _env('CLOCKBASE_METHYLOME_REFERENCE', STREAM_EMPTY / 'reference-cache' /
                           'GSE40279-means-aa10085d9e53a17a-geo_stream_round_trip_v2.csv')
EMPTY_MATRIX_LIST = INPUTS / 'empty-matrix-verification.jsonl'
PUBLISHED = INPUTS / 'published'
REPLAY_PYTHON = os.environ.get('CLOCKBASE_REPLAY_PYTHON') or sys.executable
MODEL_LAUNCHER = shlex.split(os.environ.get('CLOCKBASE_MODEL_LAUNCHER', ''))

ROOTS = {'geo': GEO, 'verifier': VERIFIER_ROOT, 'stream_empty': STREAM_EMPTY, 'stream_detached': STREAM_DETACHED,
         'corrections': CORRECTIONS, 'inputs': INPUTS, 'traces': TRACES, 'repo': REPO}


def portable(path):
    """'root:relative/path' for a path under a known root (as given or with symlinks resolved); else the path."""
    path = Path(path)
    candidates = [path]
    try:
        candidates.append(path.resolve())
    except OSError:
        pass
    for label, root in ROOTS.items():
        for base in {root, root.resolve() if root.exists() else root}:
            for candidate in candidates:
                try:
                    return f'{label}:{candidate.relative_to(base).as_posix()}'
                except ValueError:
                    continue
    return str(path)


def resolve(portable_path):
    """Inverse of portable(): the local path of 'root:relative/path'."""
    label, _, relative = str(portable_path).partition(':')
    if label in ROOTS and relative:
        return ROOTS[label] / relative
    original = Path(portable_path)
    if original.exists():
        return original
    # Archived records can predate portable prefixes. Resolve a known bundle
    # directory name under the explicitly configured local root after transfer.
    for root in sorted(ROOTS.values(), key=lambda p: len(p.parts), reverse=True):
        if root.name in original.parts:
            position = original.parts.index(root.name)
            candidate = root.joinpath(*original.parts[position + 1:])
            if candidate.exists():
                return candidate
    return original
