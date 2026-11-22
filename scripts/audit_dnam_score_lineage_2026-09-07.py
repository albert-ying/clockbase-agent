#!/usr/bin/env python3
"""Read existing DNAm scoring inputs and report comparison/clock mismatches.

This diagnostic does not change scientific source tables or rerun predictions.
The JSON report records physical CSV rows (header is row 1) and input hashes.
"""

import argparse
import csv
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import math
import re
import pickle
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CASES = ("GSE140344-40", "GSE52588-95")
FIELDS = (
    "case_id", "gse_id", "comparison_name", "condition_name", "clock",
    "pvalue", "p_value", "log2FoldChange", "effect_size", "direction",
    "fdr", "final_score", "intervention", "covariate", "statistical_test",
)


def evidence(path):
    data = path.read_bytes()
    return {
        "path": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(data).hexdigest(),
        "modified_utc": datetime.fromtimestamp(
            path.stat().st_mtime, timezone.utc).isoformat(),
    }


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [(i, row) for i, row in enumerate(csv.DictReader(handle), 2)]


def selected(row):
    return {key: row[key] for key in FIELDS if key in row}


def same_number(a, b):
    try:
        return math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=0)
    except (TypeError, ValueError):
        return False


def audit():
    master = ROOT / "result/main_results_hs_dnam_with_scores.csv"
    report = {"scope": list(CASES), "inputs": [evidence(master)], "cases": []}
    master_rows = rows(master)
    for case in CASES:
        gse = case.split("-")[0]
        sample_path = ROOT / "minipaper_fork/output_hs_dnam" / case / "sample.json"
        raw_path = ROOT / "hs_dnam_sonnet37" / gse / "aging_analysis_sonnet37.csv"
        report["inputs"].extend([evidence(sample_path), evidence(raw_path)])
        sample = json.loads(sample_path.read_text())
        raw = rows(raw_path)
        matched = [(i, row) for i, row in master_rows if row["case_id"] == case]
        if len(matched) != 1:
            raise ValueError(f"Expected one scored master row for {case}; got {len(matched)}")
        index, row = matched[0]
        same_comparison = [
            (i, r) for i, r in raw
            if r["comparison_name"].strip() == row["comparison_name"].strip()
        ]
        p_matches = [
            {"csv_row": i, **selected(r)} for i, r in same_comparison
            if same_number(r["p_value"], row["pvalue"])
        ]
        same_clock_p = [r for r in p_matches if r["clock"] == row["clock"]]
        components = {key: float(value) for key, value in row.items()
                      if key.startswith(("q0_", "q1_", "q2_", "q3_", "q4_", "q5_"))}
        report["cases"].append({
            "case_id": case,
            "scoring_input": selected(sample),
            "scored_master": {"csv_row": index, **selected(row)},
            "score_comparison_matches": (
                sample["comparison_name"].strip() == row["comparison_name"].strip()),
            "master_p_matches_source_clock": bool(same_clock_p),
            "source_rows_matching_master_p": p_matches,
            "score_components": components,
            "component_sum": sum(components.values()),
            "stored_score_equals_component_sum": same_number(
                row["final_score"], sum(components.values())),
            "disposition": "Unresolved: no source-table correction applied",
        })
    return report


def inventory():
    """Inventory saved identities without choosing replacement rows or scores."""
    master = ROOT / "result/main_results_hs_dnam_with_scores.csv"
    inputs = {str(master): evidence(master)}
    records = []
    seen = set()
    for index, row in rows(master):
        case = row["case_id"]
        if not case or case in seen:
            raise ValueError(f"Missing or duplicate scored case identity: {case!r}")
        seen.add(case)
        sample_path = ROOT / "minipaper_fork/output_hs_dnam" / case / "sample.json"
        source_path = ROOT / "hs_dnam_sonnet37" / row["gse_id"] / "aging_analysis_sonnet37.csv"
        record = {"case_id": case, "master_csv_row": index,
                  "master_comparison": row["comparison_name"],
                  "sample_path": str(sample_path.relative_to(ROOT)),
                  "source_path": str(source_path.relative_to(ROOT))}
        if not sample_path.exists():
            record["status"] = "missing_scoring_input"
            records.append(record)
            continue
        inputs[str(sample_path)] = evidence(sample_path)
        sample = json.loads(sample_path.read_text())
        if sample.get("gse_id") != row["gse_id"]:
            raise ValueError(f"Scoring-input GSE mismatch for {case}")
        comparison = sample.get("comparison_name")
        if not isinstance(comparison, str) or not comparison.strip():
            raise ValueError(f"Missing scoring-input comparison for {case}")
        record["scoring_input"] = selected(sample)
        record["status"] = ("comparison_name_match" if comparison.strip() ==
                            row["comparison_name"].strip() else "comparison_name_differs")
        if not source_path.exists():
            record["source_status"] = "missing_analysis_file"
            records.append(record)
            continue
        inputs[str(source_path)] = evidence(source_path)
        candidates = [(i, r) for i, r in rows(source_path)
                      if r["comparison_name"].strip() == comparison.strip()]
        record["source_status"] = "comparison_found" if candidates else "comparison_not_found"
        record["source_comparison_rows"] = [{"csv_row": i, **selected(r)} for i, r in candidates]
        # Matching the signed effect and p value avoids silently reversing a
        # contrast. Multiple hits remain ambiguous; none is selected here.
        record["primary_tuple_matches"] = [i for i, r in candidates
            if same_number(sample.get("pvalue"), r["p_value"])
            and same_number(sample.get("effect_size"), r["effect_size"])]
        records.append(record)
    return {"scope": "All current scored DNAm rows; exact-name inventory only",
            "limitations": "Name agreement and tuple agreement do not validate model choice, score correctness, or biological interpretation. No repair applied.",
            "counts": dict(Counter(r["status"] for r in records)),
            "source_counts": dict(Counter(r.get("source_status", "not_checked_missing_input") for r in records)),
            "total_rows": len(records), "inputs": list(inputs.values()), "records": records}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Write a new JSON evidence report")
    parser.add_argument("--inventory", action="store_true", help="Check every scored DNAm identity without repairing source tables")
    parser.add_argument("--scores", action="store_true", help="Audit stored scoring artifacts and literal cross-clock inputs")
    parser.add_argument("--design", action="store_true", help="Summarize the two Down syndrome studies from saved sample metadata")
    parser.add_argument("--trace", action="store_true", help="Extract retained Down syndrome execution messages without executing their code")
    args = parser.parse_args()
    try:
        if sum((args.inventory, args.scores, args.design, args.trace)) > 1:
            raise ValueError("Choose only one audit mode")
        report = execution_traces() if args.trace else study_design() if args.design else score_artifacts() if args.scores else inventory() if args.inventory else audit()
        content = json.dumps(report, indent=2, allow_nan=False) + "\n"
        if args.output:
            with args.output.open("x") as handle:
                handle.write(content)
        else:
            print(content, end="")
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"DNAm lineage audit failed: {exc}\n")


def score_artifacts():
    """Check saved arithmetic and prompt-input consistency, not scientific validity."""
    master = ROOT / "result/main_results_hs_dnam_with_scores.csv"
    inputs = [evidence(master)]
    records = []
    prefixes = ("q0_", "q1_", "q2_", "q3_", "q4_", "q5_")
    clocks = ("PhenoAge", "Horvathv1", "YingDamAge", "DunedinPACE")

    def components(data):
        result = {}
        for key, value in data.items():
            if key.startswith(prefixes) and isinstance(value, dict):
                field = "penalty_or_bonus_applied" if key.startswith("q0_") else "points"
                if field not in value:
                    raise ValueError(f"Missing {field} in {key}")
                result[key] = float(value[field])
        if len(result) != 8 or not all(math.isfinite(v) for v in result.values()):
            raise ValueError("Expected eight finite scored components, excluding Q7")
        return result

    for index, row in rows(master):
        case = row["case_id"]
        folder = ROOT / "minipaper_fork/output_hs_dnam" / case
        paths = {name: folder / name for name in (
            "sample.json", "identified_relevance.json",
            "identify_relevance_prompt.txt", "identify_relevance_response.txt")}
        record = {"case_id": case, "master_csv_row": index,
                  "missing_artifacts": [name for name, path in paths.items() if not path.exists()]}
        if record["missing_artifacts"]:
            records.append(record)
            continue
        inputs.extend(evidence(path) for path in paths.values())
        sample = json.loads(paths["sample.json"].read_text())
        parsed = json.loads(paths["identified_relevance.json"].read_text())
        response_text = paths["identify_relevance_response.txt"].read_text()
        start = response_text.find("{")
        try:
            response, _ = json.JSONDecoder().raw_decode(response_text[start:])
        except ValueError as exc:
            record["response_parse_error"] = str(exc)
            response = None
        parsed_components = components(parsed)
        response_components = components(response) if response is not None else None
        prompt = paths["identify_relevance_prompt.txt"].read_text()
        cross = {clock: {"effect": sample.get(clock + "_effect_size"),
                         "pvalue": sample.get(clock + "_pvalue")} for clock in clocks}
        complete = all(isinstance(v, (int, float)) and math.isfinite(v)
                       for values in cross.values() for v in values.values())
        primary = float(sample["effect_size"])
        literal_support = sum(v["pvalue"] < 0.05 and v["effect"] * primary > 0
                              for v in cross.values()) if complete else None
        nonfinite = []
        for clock, values in cross.items():
            for key, value in values.items():
                if isinstance(value, float) and not math.isfinite(value):
                    nonfinite.append(f"{clock}.{key}={value}")
                    values[key] = None
        record.update({
            "comparison": sample["comparison_name"],
            "primary_effect": primary, "direction_field": sample["direction"],
            "prompt_names_grimage_primary": "Primary clock — GrimAgeV2" in prompt,
            "prompt_excludes_q7": "Additional Context Questions (Not scored)" in prompt,
            "parsed_components": parsed_components,
            "parsed_component_sum": sum(parsed_components.values()),
            "parsed_final_score": parsed["final_score"],
            "parsed_arithmetic_matches": same_number(sum(parsed_components.values()), parsed["final_score"]),
            "response_components": response_components,
            "response_final_score": response["final_score"] if response is not None else None,
            "response_arithmetic_matches": same_number(sum(response_components.values()), response["final_score"]) if response is not None else None,
            "parsed_matches_response_scoring": (parsed_components == response_components and same_number(parsed["final_score"], response["final_score"])) if response is not None else None,
            "master_final_score": row["final_score"],
            "master_matches_parsed_final": same_number(row["final_score"], parsed["final_score"]),
            "cross_clock_inputs": cross, "stored_n_support": sample.get("n_support"),
            "nonfinite_cross_clock_inputs": nonfinite,
            "literal_signed_support_count": literal_support,
            "literal_support_matches_stored": literal_support == sample.get("n_support") if complete else None,
            "effect_sign_matches_direction": primary * float(sample["direction"]) > 0,
        })
        records.append(record)
    checked = [r for r in records if "parsed_arithmetic_matches" in r]
    flags = ("parsed_arithmetic_matches", "response_arithmetic_matches",
             "parsed_matches_response_scoring", "master_matches_parsed_final",
             "literal_support_matches_stored", "effect_sign_matches_direction",
             "prompt_names_grimage_primary", "prompt_excludes_q7")
    return {"scope": "Stored DNAm score artifacts; no rescore or source-table repair",
            "limitations": "Raw responses use strict JSON parsing; failure disables only response-dependent checks. Each flag has its own available denominator. Literal support compares displayed signed effects only; it does not validate biological orientation, statistical models, source literature, or justify restoring excluded scores. Correction provenance requires separate review.",
            "total_rows": len(records), "checked_rows": len(checked),
            "parseable_response_rows": sum(r["response_arithmetic_matches"] is not None for r in checked),
            "available_denominators": {flag: sum(r[flag] is not None for r in checked) for flag in flags},
            "flag_counts": {flag: dict(Counter(str(r[flag]) for r in checked)) for flag in flags},
            "inputs": inputs, "records": records}


def study_design():
    """Describe saved sample metadata; do not fit models or recompute clocks."""
    path = ROOT.parent / "_downloads_archive/2026-06-18/revision-and-figures/clockbase_export_2026-06-15/clockbase_persample_human_dnam_2026-06-15.csv"
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    studies = {g: [] for g in ("GSE52588", "GSE140344")}
    with path.open(newline="") as handle:
        for index, row in enumerate(csv.DictReader(handle), 2):
            if row["gse_id"] not in studies:
                continue
            raw = row["raw_characteristics"]
            def field(name):
                match = re.search(r"(?:=|\|\|\s*)" + re.escape(name) + r":\s*([^|]+)", raw)
                return match.group(1).strip() if match else None
            if row["gse_id"] == "GSE52588":
                title = row["title"].lower()
                arm = "mother" if "mother" in title else "sibling" if "sibling" in title else "DS" if "subject affected" in title else "unresolved"
                group = (field("group") or "").lower()
                raw_arm = "mother" if "mother" in group else "sibling" if "sibling" in group else "DS" if "subject affected" in group else "unresolved"
                if arm != raw_arm:
                    raise ValueError(f"Title and raw group disagree at CSV row {index}")
            else:
                disease, cognition = field("disease state"), field("cognition state")
                arm = "control" if disease == "control" else "T21_" + str(cognition) if disease == "trisomy 21" else "unresolved"
            if arm == "unresolved" or field("age") is None:
                raise ValueError(f"Unresolved arm or missing raw age at CSV row {index}")
            age = float(field("age"))
            if not math.isfinite(age):
                raise ValueError(f"Nonfinite raw age at CSV row {index}")
            studies[row["gse_id"]].append({
                "csv_row": index, "gsm_id": row["gsm_id"], "arm": arm,
                "raw_age": age, "raw_family": field("family"),
                "raw_gender": field("gender"), "raw_cognition": field("cognition state"),
                "raw_characteristics": raw,
            })
    results = {}
    for gse, samples in studies.items():
        if not samples or len({s["gsm_id"] for s in samples}) != len(samples):
            raise ValueError(f"Missing samples or duplicate GSMs for {gse}")
        groups = {}
        for arm in sorted({s["arm"] for s in samples}):
            subset = [s for s in samples if s["arm"] == arm]
            ages = [s["raw_age"] for s in subset]
            groups[arm] = {"n": len(subset), "age_min": min(ages), "age_max": max(ages),
                           "age_mean": sum(ages) / len(ages),
                           "family_counts": dict(Counter(s["raw_family"] for s in subset if s["raw_family"] is not None))}
        if gse == "GSE52588":
            if set(groups) != {"DS", "mother", "sibling"}:
                raise ValueError("GSE52588 expected three observed family arms")
            families = [set(group["family_counts"]) for group in groups.values()]
            if not all(family == families[0] for family in families) or not all(
                len(group["family_counts"]) == group["n"] and all(n == 1 for n in group["family_counts"].values())
                for group in groups.values()
            ):
                raise ValueError("GSE52588 family identities do not form complete triplets")
        source = ROOT / "hs_dnam_sonnet37" / gse / "aging_analysis_sonnet37.csv"
        results[gse] = {"n_samples": len(samples), "groups": groups, "samples": samples,
                        "analysis_source": evidence(source),
                        "grimage_source_rows": [{"csv_row": i, **selected(r)} for i, r in rows(source) if r["clock"] == "GrimAgeV2"]}
    return {"scope": "Raw exported sample metadata and existing source model labels; no prediction or model rerun",
            "input": {"path": str(path), "sha256": digest.hexdigest(), "modified_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()},
            "limitations": "Metadata identifies observed sample arms and age distributions. It does not establish the executed regression formula, clock implementation, or scientific validity of stored model results.",
            "studies": results}


def execution_traces():
    """Deserialize only plain stored data; never execute stored analysis code."""
    class DataOnlyUnpickler(pickle.Unpickler):
        def find_class(self, module, name):
            raise ValueError(f"Unexpected pickle global: {module}.{name}")

        def persistent_load(self, pid):
            raise ValueError("Unexpected persistent pickle reference")

    path = ROOT / "successful_traces_human_dnam.pkl"
    with path.open("rb") as handle:
        data = DataOnlyUnpickler(handle).load()
    studies = {}
    for gse in ("GSE52588", "GSE140344"):
        result = data[gse]["result"]
        messages = result[0]
        if not isinstance(messages, list) or not all(isinstance(m, str) for m in messages):
            raise ValueError(f"Unexpected message representation for {gse}")
        studies[gse] = {
            "stored_status": data[gse]["status"],
            "messages": [{"message_index_zero_based": i, "text": message}
                         for i, message in enumerate(messages)],
            "stored_summary": result[1],
        }
    return {"scope": "Retained execution messages for two studies; extraction only",
            "input": evidence(path),
            "limitations": "Stored success status and generated summaries do not establish scientific validity. No stored code was executed and no model was rerun.",
            "studies": studies}


if __name__ == "__main__":
    main()
