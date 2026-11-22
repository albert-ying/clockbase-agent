"""
Autonomous post-audit driver for ClockBase Agent results.

Loads the FDR-significant rows produced by the primary Sonnet-3.7 analysis
layer (the `main_results_*_with_scores.csv` files in `result/`) and runs an
independent second-pass biological-plausibility audit using Anthropic
Claude Opus 4.7 at temperature 0.

Two operating modes (selected via `--mode`):

  plausibility   re-audit the displayed effect direction and emit a
                 q0_adjustment in {-50, 0, +30}. Per-case verdicts are
                 written to `result/case_audit_results/{case_id}.json` and
                 aggregated into `result/dir_audit_reaudit.json`.

  curate         emit short_name, display_name, control_is_healthy,
                 comparison_inverted, and an optional naming_correction
                 per case. Outputs are chunked into
                 `result/case_audit_results/chunk_NN.json` mimicking the
                 layout already present from the original interactive
                 audit run.

This script is the canonical CLI entry point for the autonomous Opus 4.7
post-audit pass that the manuscript Methods section references. It uses
the same Anthropic API pattern as `src/coding_agent.py`, but overrides the
default temperature (0.7) to 0 because this is an audit, not generation.

The prompts are kept in `prompts/post_audit_opus47.txt` and
`prompts/case_curator_opus47.txt`; see those files for the audit rubric.

Usage
-----

    # plausibility re-audit of all FDR-significant mouse RNA rows
    export ANTHROPIC_API_KEY=...
    python src/autonomous_post_audit.py \\
        --mode plausibility \\
        --input result/main_results_result_mm_rna_with_scores.csv \\
        --output-dir result/case_audit_results \\
        --aggregate result/dir_audit_reaudit.json

    # curate a single chunk of 100 cases for a smoke test
    python src/autonomous_post_audit.py \\
        --mode curate \\
        --input result/dir_audit_chunks/chunk_01.csv \\
        --output-dir result/case_audit_results \\
        --chunk-name chunk_01 \\
        --limit 5

Requirements
------------

The script only needs `anthropic`, `pandas`, and the standard library. The
existing `src/coding_agent.py` already pins `langchain_anthropic`, which
pulls in the lower-level `anthropic` client transitively. If running this
script in a slimmer environment, install:

    pip install anthropic pandas
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

try:
    # The official Anthropic Python SDK. This is the same client that
    # langchain_anthropic wraps in coding_agent.py.
    import anthropic
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "anthropic package is required. Install with `pip install anthropic`."
    ) from exc


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parent.parent
PROMPTS_DIR = REPO_ROOT / "prompts"

# Model identifier for the post-audit pass. Manuscript Methods commits to
# Claude Opus 4.7 here, distinct from the primary `claude-3-7-sonnet-...`
# coding agent.
DEFAULT_MODEL = "claude-opus-4-7"

# Temperature override. The primary coding agent runs at 0.7 (for plan-and-
# act exploration), but plausibility audit is an adjudication task and runs
# at 0 for reproducibility.
DEFAULT_TEMPERATURE = 0.0

DEFAULT_MAX_TOKENS = 1024

# How many seconds to back off between Anthropic API errors.
RETRY_BASE_DELAY_SEC = 4.0
MAX_RETRIES = 5


# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("autonomous_post_audit")


# ---------------------------------------------------------------------------
# Prompt loading
# ---------------------------------------------------------------------------

def load_prompt(name: str) -> str:
    """Load a prompt template from `prompts/{name}.txt`.

    Strips leading `#`-prefixed comment lines so the actual prompt body is
    what gets sent to the model.
    """
    path = PROMPTS_DIR / f"{name}.txt"
    if not path.exists():
        raise FileNotFoundError(
            f"Prompt template not found: {path}. "
            f"The `prompts/` directory is part of the repository — "
            f"did you delete it?"
        )
    raw = path.read_text(encoding="utf-8").splitlines()

    # Strip leading metadata block (consecutive lines starting with '#').
    body_start = 0
    for i, line in enumerate(raw):
        if line.startswith("#") or line.strip() == "":
            body_start = i + 1
            continue
        break

    return "\n".join(raw[body_start:]).strip()


# ---------------------------------------------------------------------------
# Anthropic client
# ---------------------------------------------------------------------------

def get_client(api_key_env: str) -> anthropic.Anthropic:
    """Construct an Anthropic client using an API key from the environment."""
    key = os.environ.get(api_key_env)
    if not key:
        raise SystemExit(
            f"Environment variable {api_key_env!r} is not set. "
            f"Export it (e.g., `export {api_key_env}=sk-ant-...`) before "
            f"running this script."
        )
    return anthropic.Anthropic(api_key=key)


def call_model(
    client: anthropic.Anthropic,
    model: str,
    prompt: str,
    temperature: float,
    max_tokens: int,
) -> str:
    """Call the Anthropic Messages API with retry-on-error.

    The post-audit prompts are designed to return STRICT JSON; we leave
    JSON parsing to the caller so retries are decoupled from validation.
    """
    last_exc: Exception | None = None
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = client.messages.create(
                model=model,
                max_tokens=max_tokens,
                temperature=temperature,
                messages=[{"role": "user", "content": prompt}],
            )
            # The Anthropic SDK returns a list of content blocks; we want
            # the concatenated text.
            parts = []
            for block in response.content:
                text = getattr(block, "text", None)
                if text:
                    parts.append(text)
            return "".join(parts).strip()
        except Exception as exc:  # pragma: no cover - network path
            last_exc = exc
            delay = RETRY_BASE_DELAY_SEC * attempt
            log.warning(
                "Anthropic call failed (attempt %d/%d): %s. Sleeping %.1fs.",
                attempt, MAX_RETRIES, exc, delay,
            )
            time.sleep(delay)
    assert last_exc is not None
    raise last_exc


# ---------------------------------------------------------------------------
# JSON parsing helpers
# ---------------------------------------------------------------------------

def parse_strict_json(raw: str) -> dict[str, Any]:
    """Parse a model response that is supposed to be a single JSON object.

    Tolerates a leading/trailing code fence (```json ... ```) but errors out
    on anything else — we deliberately do not let the model wander off the
    schema.
    """
    text = raw.strip()
    if text.startswith("```"):
        # strip first line and trailing fence
        lines = text.splitlines()
        if lines[0].startswith("```"):
            lines = lines[1:]
        if lines and lines[-1].startswith("```"):
            lines = lines[:-1]
        text = "\n".join(lines).strip()
    return json.loads(text)


def safe_str(v: Any, default: str = "") -> str:
    """Coerce a pandas cell value into a safe-for-prompt string."""
    if v is None:
        return default
    try:
        if pd.isna(v):
            return default
    except (TypeError, ValueError):
        pass
    return str(v)


# ---------------------------------------------------------------------------
# Modes
# ---------------------------------------------------------------------------

def _fill_plausibility_prompt(template: str, row: pd.Series) -> str:
    return template.format(
        case_id=safe_str(row.get("case_id"), default="NA"),
        gse_id=safe_str(row.get("gse_id"), default="NA"),
        study_title=safe_str(row.get("title"), default=""),
        study_summary=safe_str(row.get("summary"), default="")[:1500],
        comparison_name=safe_str(row.get("comparison_name"), default=""),
        condition_name=safe_str(row.get("condition_name"), default=""),
        condition_category=safe_str(row.get("condition_category"), default=""),
        clock=safe_str(row.get("clock"), default=""),
        log2fc=safe_str(row.get("log2fc") or row.get("effect_size"), default=""),
        p_value=safe_str(row.get("p_value"), default=""),
        trust_score=safe_str(row.get("trust_score"), default=""),
        primary_q0_rationale=safe_str(row.get("q0_rationale"), default="")[:1000],
    )


def _fill_curator_prompt(template: str, row: pd.Series) -> str:
    return template.format(
        case_id=safe_str(row.get("case_id"), default="NA"),
        gse_id=safe_str(row.get("gse_id"), default="NA"),
        study_title=safe_str(row.get("title"), default=""),
        study_summary=safe_str(row.get("summary"), default="")[:1500],
        comparison_name=safe_str(row.get("comparison_name"), default=""),
        condition_name=safe_str(row.get("condition_name"), default=""),
        condition_category=safe_str(row.get("condition_category"), default=""),
    )


def run_plausibility(
    df: pd.DataFrame,
    client: anthropic.Anthropic,
    model: str,
    temperature: float,
    output_dir: Path,
    aggregate_path: Path | None,
    limit: int | None,
    overwrite: bool,
) -> list[dict[str, Any]]:
    template = load_prompt("post_audit_opus47")
    output_dir.mkdir(parents=True, exist_ok=True)

    aggregated: list[dict[str, Any]] = []
    rows: Iterable[pd.Series] = (row for _, row in df.iterrows())
    n_done = 0

    for row in rows:
        if limit is not None and n_done >= limit:
            break

        case_id = safe_str(row.get("case_id"), default="").strip()
        if not case_id:
            log.warning("Row missing case_id; skipping. row=%s", dict(row))
            continue

        per_case_path = output_dir / f"{case_id}.json"
        if per_case_path.exists() and not overwrite:
            log.info("[%s] cached — loading existing verdict", case_id)
            try:
                aggregated.append(json.loads(per_case_path.read_text()))
            except Exception as exc:
                log.warning("[%s] cached file unparseable (%s); re-running",
                            case_id, exc)
            else:
                n_done += 1
                continue

        prompt = _fill_plausibility_prompt(template, row)
        try:
            raw = call_model(client, model, prompt, temperature, DEFAULT_MAX_TOKENS)
            verdict = parse_strict_json(raw)
        except Exception as exc:
            log.error("[%s] audit failed: %s", case_id, exc)
            verdict = {
                "case_id": case_id,
                "q0_adjustment": 0,
                "action": "revert_to_zero",
                "rationale": f"AUDIT_ERROR: {exc}",
                "_error": True,
            }

        if verdict.get("q0_adjustment") not in (-50, 0, 30):
            log.warning("[%s] invalid q0_adjustment=%r; coercing to 0",
                        case_id, verdict.get("q0_adjustment"))
            verdict["q0_adjustment"] = 0
            verdict["action"] = "revert_to_zero"

        verdict.setdefault("case_id", case_id)
        per_case_path.write_text(json.dumps(verdict, indent=2))
        aggregated.append(verdict)
        n_done += 1
        log.info("[%s] q0_adjustment=%s action=%s",
                 case_id, verdict.get("q0_adjustment"), verdict.get("action"))

    if aggregate_path is not None:
        aggregate_path.parent.mkdir(parents=True, exist_ok=True)
        # Dump the q0_adjustment-only summary in the legacy
        # dir_audit_reaudit.json schema: list of {case_id, action, rationale}.
        legacy = [
            {
                "case_id": v["case_id"],
                "action": v.get("action", "revert_to_zero"),
                "rationale": v.get("rationale", ""),
            }
            for v in aggregated
        ]
        aggregate_path.write_text(json.dumps(legacy, indent=2))
        log.info("aggregated %d verdicts -> %s", len(legacy), aggregate_path)

    return aggregated


def run_curate(
    df: pd.DataFrame,
    client: anthropic.Anthropic,
    model: str,
    temperature: float,
    output_dir: Path,
    chunk_name: str | None,
    limit: int | None,
    overwrite: bool,
) -> list[dict[str, Any]]:
    template = load_prompt("case_curator_opus47")
    output_dir.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, Any]] = []
    n_done = 0

    for _, row in df.iterrows():
        if limit is not None and n_done >= limit:
            break

        case_id = safe_str(row.get("case_id"), default="").strip()
        if not case_id:
            log.warning("Row missing case_id; skipping.")
            continue

        prompt = _fill_curator_prompt(template, row)
        try:
            raw = call_model(client, model, prompt, temperature, DEFAULT_MAX_TOKENS)
            record = parse_strict_json(raw)
        except Exception as exc:
            log.error("[%s] curation failed: %s", case_id, exc)
            record = {
                "case_id": case_id,
                "short_name": "",
                "display_name": "",
                "control_is_healthy": None,
                "comparison_inverted": None,
                "naming_correction": None,
                "_error": str(exc),
            }

        record.setdefault("case_id", case_id)
        records.append(record)
        n_done += 1
        log.info("[%s] short_name=%r control_is_healthy=%s",
                 case_id, record.get("short_name"),
                 record.get("control_is_healthy"))

    chunk_name = chunk_name or "curated"
    chunk_path = output_dir / f"{chunk_name}.json"
    if chunk_path.exists() and not overwrite:
        backup = chunk_path.with_suffix(chunk_path.suffix + ".precurate_bak")
        chunk_path.rename(backup)
        log.info("backed up existing %s -> %s", chunk_path, backup)
    chunk_path.write_text(json.dumps(records, indent=2))
    log.info("wrote %d curated records -> %s", len(records), chunk_path)

    return records


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Autonomous post-audit driver (Claude Opus 4.7).",
    )
    p.add_argument(
        "--mode",
        choices=["plausibility", "curate"],
        required=True,
        help="plausibility: q0_adjustment audit. curate: short/display name + direction flags.",
    )
    p.add_argument(
        "--input",
        type=Path,
        required=True,
        help="Path to the CSV holding the cases to audit. Expected columns include "
             "case_id, gse_id, condition_name, comparison_name, etc.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "result" / "case_audit_results",
        help="Directory to write per-case JSONs into.",
    )
    p.add_argument(
        "--aggregate",
        type=Path,
        default=None,
        help="(plausibility mode) Path for the aggregated dir_audit_reaudit.json file.",
    )
    p.add_argument(
        "--chunk-name",
        default=None,
        help="(curate mode) Basename (without .json) for the chunk output file.",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Stop after this many rows (useful for smoke tests).",
    )
    p.add_argument(
        "--model",
        default=DEFAULT_MODEL,
        help=f"Anthropic model id (default: {DEFAULT_MODEL}).",
    )
    p.add_argument(
        "--temperature",
        type=float,
        default=DEFAULT_TEMPERATURE,
        help=f"Temperature (default: {DEFAULT_TEMPERATURE}).",
    )
    p.add_argument(
        "--api-key-env",
        default="ANTHROPIC_API_KEY",
        help="Environment variable to read the API key from.",
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-run cases even if a per-case JSON already exists.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.input.exists():
        log.error("input CSV not found: %s", args.input)
        return 2

    log.info("loading %s", args.input)
    df = pd.read_csv(args.input)
    log.info("loaded %d rows", len(df))

    client = get_client(args.api_key_env)

    if args.mode == "plausibility":
        run_plausibility(
            df=df,
            client=client,
            model=args.model,
            temperature=args.temperature,
            output_dir=args.output_dir,
            aggregate_path=args.aggregate,
            limit=args.limit,
            overwrite=args.overwrite,
        )
    elif args.mode == "curate":
        run_curate(
            df=df,
            client=client,
            model=args.model,
            temperature=args.temperature,
            output_dir=args.output_dir,
            chunk_name=args.chunk_name,
            limit=args.limit,
            overwrite=args.overwrite,
        )
    else:  # pragma: no cover - argparse enforces the choices
        raise ValueError(f"unknown mode {args.mode!r}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
