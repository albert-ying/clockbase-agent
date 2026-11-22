#!/usr/bin/env python3
"""Identity-only score reconciliation. Legacy production merge is disabled.

No scores are transferred. A match does not validate science, scoring or FDR.
run_pipeline.sh intentionally stops here until score publication is repaired.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
TEXT_KEYS = ("gse_id", "comparison_name", "condition_name", "condition_category")
NUMBER_KEYS = ("effect_size", "pvalue", "direction")
PVALUE_RTOL = 1e-10


def resolve_identity(master, sample, primary_clock):
    """Return unique historical tuple matches without validating their science."""
    required = (*TEXT_KEYS, "clock", *NUMBER_KEYS)
    if not set(required).issubset(master.columns) or not master.index.is_unique:
        raise ValueError("Master requires identity columns and unique row indices")
    if any(not isinstance(sample.get(k), str) or not sample[k] for k in TEXT_KEYS):
        return {"status": "invalid_scoring_identity", "reason": "Missing text identity"}
    try:
        numbers = {k: float(sample[k]) for k in NUMBER_KEYS}
    except (KeyError, TypeError, ValueError):
        return {"status": "invalid_scoring_identity", "reason": "Missing numeric identity"}
    if (not all(np.isfinite(v) for v in numbers.values())
            or not 0 <= numbers["pvalue"] <= 1
            or numbers["direction"] not in (-1, 1)
            or not isinstance(primary_clock, str) or not primary_clock
            or sample.get("clock", primary_clock) != primary_clock):
        return {"status": "invalid_scoring_identity", "reason": "Invalid statistics or clock"}
    candidates, matches = [], []
    for index, row in master[master.gse_id == sample["gse_id"]].iterrows():
        reasons = [k for k in TEXT_KEYS if row[k] != sample[k]]
        if row["clock"] != primary_clock:
            reasons.append("clock")
        for key in NUMBER_KEYS:
            try:
                value = float(row[key])
                equal = (np.isclose(value, numbers[key], rtol=PVALUE_RTOL, atol=0)
                         if key == "pvalue" else value == numbers[key])
                if not np.isfinite(value) or not equal:
                    reasons.append(key)
            except (ValueError, TypeError):
                reasons.append(key)
        candidates.append({"master_index": int(index), "csv_line": int(index) + 2,
                           "mismatched_fields": reasons})
        if not reasons:
            matches.append(int(index))
    result = {"status": "unique_identity_match" if len(matches) == 1 else
                        "ambiguous_identity" if matches else "no_identity_match",
              "matching_indices": matches, "candidates": candidates,
              "scientific_status": "not_validated"}
    if len(matches) == 1:
        result["master_index"] = matches[0]
    if sample["gse_id"] == "GSE52588":
        result["scientific_status"] = "historical_analysis_invalid_clock_derived_groups"
    return result


def evidence(path):
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "mtime_ns": path.stat().st_mtime_ns}


def audit(modality, primary_clock):
    master_path = ROOT / f"result/main_results_hs_{modality}_outlier_marked.csv"
    scores_path = ROOT / f"data/scores_summary_hs_{modality}.csv"
    master = pd.read_csv(master_path, low_memory=False)
    scores = pd.read_csv(scores_path)
    if not {"dataset", "gse_id"}.issubset(scores.columns):
        raise ValueError("Score summary lacks case identity")
    if scores.dataset.isna().any() or scores.dataset.duplicated().any():
        raise ValueError("Missing or duplicate scored case identity")
    inputs = [evidence(master_path), evidence(scores_path), evidence(Path(__file__))]
    if modality == "dnam":
        inputs.append(evidence(ROOT / "s02_prepare_scoring_input_human.py"))
    cases = []
    for _, score in scores.iterrows():
        case = score.dataset
        if not isinstance(case, str) or not re.fullmatch(r"GSE[0-9]+-[0-9]+", case):
            raise ValueError("Invalid case-directory identity")
        sample_path = ROOT / f"minipaper_fork/output_hs_{modality}" / case / "sample.json"
        record = {"case_id": case, "sample_path": str(sample_path)}
        if not sample_path.is_file():
            record.update(status="missing_scoring_input")
        else:
            inputs.append(evidence(sample_path))
            sample = json.loads(sample_path.read_text())
            if sample.get("gse_id") != score.gse_id or case.split("-")[0] != score.gse_id:
                record.update(status="case_gse_conflict")
            else:
                record.update(resolve_identity(master, sample, primary_clock))
                record["scoring_identity"] = {k: sample.get(k) for k in (*TEXT_KEYS, *NUMBER_KEYS)}
                record["primary_clock"] = primary_clock
        cases.append(record)
    return {"status": "identity_only; scores and science not validated; no production writes",
            "primary_clock": primary_clock,
            "primary_clock_basis": "s02 DNAM_PRIMARY_CLOCK contract" if modality == "dnam" else "explicit CLI argument",
            "numeric_matching": {"pvalue_rtol": PVALUE_RTOL, "pvalue_atol": 0,
                                 "effect_size": "exact signed parsed-float equality", "direction": "literal equality"},
            "warning": "Historical matches do not validate groups, models, predictions, scores, FDR or ranks. GSE52588 used clock-derived groups.",
            "inputs": inputs, "counts": dict(Counter(c["status"] for c in cases)), "cases": cases}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modality", choices=("dnam", "rna"), required=True)
    parser.add_argument("--primary-clock")
    parser.add_argument("--audit-output", type=Path, required=True,
                        help="New JSON only; production CSV merge disabled")
    args = parser.parse_args()
    if args.audit_output.suffix.lower() != ".json" or args.audit_output.exists():
        parser.error("Specify a new .json report; overwriting files is forbidden")
    if args.modality == "rna" and not args.primary_clock:
        parser.error("RNA requires an explicit --primary-clock")
    if args.modality == "dnam" and args.primary_clock not in (None, "GrimAgeV2"):
        parser.error("DNAm scorer primary clock is GrimAgeV2")
    try:
        report = audit(args.modality, args.primary_clock or "GrimAgeV2")
        serialized = json.dumps(report, indent=2, allow_nan=False) + "\n"
        with args.audit_output.open("x") as handle:
            handle.write(serialized)
        print(json.dumps(report["counts"]))
        print("Identity audit saved. No scores or production tables changed.")
    except (OSError, ValueError, KeyError, TypeError) as exc:
        parser.exit(1, f"Identity audit failed: {exc}\n")


if __name__ == "__main__":
    main()
