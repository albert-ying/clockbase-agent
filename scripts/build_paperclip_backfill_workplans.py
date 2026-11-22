#!/usr/bin/env python3
"""
Build workplan JSONs for the paperclip backfill — every parsed_report.md
across the three corpora that lacks `citations.gxl.ai` references.

Mirrors the 2026-05-14 _tier80 schema so the runner can treat backfill chunks
the same as the original Round-2 chunks.

Output: data/minipaper_workplans_paperclip_backfill/
  <modality>__<side>__chunk<i>of<N>.json
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PAPERCLIP_SIG = re.compile(r"citations\.gxl\.ai|ncbi\.nlm\.nih\.gov/pmc", re.IGNORECASE)

CORPORA = [
    dict(
        modality="mouse RNA-seq",
        report_dir="data/report",
        artifact_base_dir="data/output",
        output_base_dir="data/report",
        prompt_template_name="GENERATE_PAPER_2",
    ),
    dict(
        modality="human RNA-seq",
        report_dir="data/report_hs_rna",
        artifact_base_dir="minipaper_fork/output",
        output_base_dir="data/report_hs_rna",
        prompt_template_name="GENERATE_PAPER_2_HUMAN_RNA",
    ),
    dict(
        modality="human DNAm",
        report_dir="data/report_hs_dnam",
        artifact_base_dir="minipaper_fork/output_hs_dnam",
        output_base_dir="data/report_hs_dnam",
        prompt_template_name="GENERATE_PAPER_2_HUMAN_DNAM",
    ),
]


def side_of(direction) -> str:
    try:
        return "pro" if float(direction) > 0 else "anti"
    except (TypeError, ValueError):
        return "unknown"


def slug(modality: str) -> str:
    return modality.replace(" ", "_").replace("-", "_")


def build_case(case_id: str, artifact_dir: Path, output_path: Path, rank: int) -> dict | None:
    """Build a workplan case from the artifact dir's sample.json."""
    sample_p = artifact_dir / "sample.json"
    if not sample_p.exists():
        return None
    try:
        s = json.loads(sample_p.read_text())
    except json.JSONDecodeError:
        return None
    return {
        "rank": rank,
        "case_id": case_id,
        "gse_id": s.get("gse_id"),
        "comparison_name": s.get("comparison_name"),
        "condition": s.get("condition_name"),
        "condition_category": s.get("condition_category"),
        "intervention": s.get("condition_name"),  # mirror Round-2 fallback
        "effect_size": s.get("effect_size"),
        "effect_direction": "positive"
        if (s.get("direction") or 0) > 0
        else ("negative" if (s.get("direction") or 0) < 0 else "unknown"),
        "effect_pvalue": s.get("pvalue"),
        "fdr": None,
        "final_score": None,
        "composite_score": None,
        "artifact_dir": str(artifact_dir.relative_to(ROOT)),
        "output_path": str(output_path.relative_to(ROOT)),
    }


def find_missing_cases(corpus: dict) -> list[dict]:
    report_dir = ROOT / corpus["report_dir"]
    artifact_root = ROOT / corpus["artifact_base_dir"]
    missing = []
    for sub in sorted(p for p in report_dir.iterdir() if p.is_dir()):
        case_id = sub.name
        md = sub / "parsed_report.md"
        has_paperclip = md.exists() and bool(
            PAPERCLIP_SIG.search(md.read_text(errors="ignore"))
        )
        if has_paperclip:
            continue
        artifact = artifact_root / case_id
        if not artifact.exists():
            # try matching by gse_id prefix
            cands = list(artifact_root.glob(f"{case_id.split('-')[0]}-*"))
            if cands:
                artifact = cands[0]
        if not artifact.exists():
            continue
        c = build_case(case_id, artifact, md, rank=0)
        if c is None:
            continue
        missing.append(c)
    return missing


def chunk(items: list, size: int) -> list[list]:
    return [items[i : i + size] for i in range(0, len(items), size)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunk-size", type=int, default=50)
    ap.add_argument(
        "--out",
        default="data/minipaper_workplans_paperclip_backfill",
    )
    args = ap.parse_args()

    out_root = ROOT / args.out
    out_root.mkdir(parents=True, exist_ok=True)

    grand_total = 0
    for corpus in CORPORA:
        missing = find_missing_cases(corpus)
        # Split by pro/anti (using direction)
        by_side: dict[str, list[dict]] = {"pro": [], "anti": [], "unknown": []}
        for c in missing:
            side = "pro" if c["effect_direction"] == "positive" else (
                "anti" if c["effect_direction"] == "negative" else "unknown"
            )
            by_side[side].append(c)

        for side, cases in by_side.items():
            if not cases:
                continue
            # rank by |effect_size| descending, then by 1/pvalue
            cases.sort(
                key=lambda c: (
                    -abs(c["effect_size"] or 0),
                    c["effect_pvalue"] if c["effect_pvalue"] is not None else 1.0,
                )
            )
            for i, c in enumerate(cases, 1):
                c["rank"] = i

            chunks = chunk(cases, args.chunk_size)
            n_chunks = len(chunks)
            for idx, ch in enumerate(chunks, 1):
                fname = (
                    f"{slug(corpus['modality'])}__{side}__"
                    f"chunk{idx}of{n_chunks}.json"
                )
                payload = {
                    "modality": corpus["modality"],
                    "side": side,
                    "chunk": f"{idx}of{n_chunks}",
                    "tier": "paperclip_backfill",
                    "prompt_template_name": corpus["prompt_template_name"],
                    "artifact_base_dir": corpus["artifact_base_dir"],
                    "output_base_dir": corpus["output_base_dir"],
                    "n_cases": len(ch),
                    "cases": ch,
                }
                (out_root / fname).write_text(json.dumps(payload, indent=2))
            print(
                f"  {corpus['modality']:15s} | {side:7s} | "
                f"{len(cases):4d} cases -> {n_chunks} chunks"
            )
            grand_total += len(cases)

    print(f"\nTOTAL backfill cases written: {grand_total}")
    print(f"Workplan dir: {out_root.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
