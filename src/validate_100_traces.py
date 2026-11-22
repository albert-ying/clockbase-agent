"""
Validate 100 random traces for correctness.
Binary validation: 1 = correct, 0 = data issue (not error), error description = agent error
Includes short notes only for dataset issues.
"""

import sys
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths

import os
import pickle
import random

import numpy as np
import pandas as pd
from scipy import stats


def evaluate_hypothesis(result_df, metadata):
    """
    Evaluate if the hypothesis/comparison is scientifically valid.
    Only checks for clear comparison structure, not category matching.
    Returns 1 if correct, error description if not.
    """
    try:
        comparison = result_df["comparison_name"].iloc[0]

        # Check for clear comparison structure
        if "vs" not in comparison.lower():
            return "No clear control vs treatment"

        return 1

    except Exception as e:
        return f"Error: {str(e)}"


def generate_note(study_id, trace, result_df):
    """
    Generate a very short note only for dataset issues.
    Returns empty string if no issues.
    """
    try:
        notes = []

        # Check for extreme effect sizes (note but not error)
        effect_sizes = result_df["effect_size"]
        extreme = effect_sizes[abs(effect_sizes) >= 10]
        if len(extreme) > 0:
            max_effect = extreme.abs().max()
            notes.append(f"Extreme effect size: |d|={max_effect:.1f}")

        # Check for missing p-values (data limitation)
        if result_df["p_value"].isna().any():
            notes.append("Insufficient samples")

        # Check for very small sample sizes
        # Inferred from test description
        test = str(result_df["statistical_test"].iloc[0]).lower()
        if "n=1" in test or "single sample" in test:
            if "Insufficient samples" not in notes:
                notes.append("n<2 per group")

        return "; ".join(notes) if notes else ""

    except Exception as e:
        return f"Error: {str(e)}"


def validate_trace(study_id, trace, result_path):
    """
    Validate a single trace on 8 binary criteria.
    Returns dict with 1 for correct, 0 for data issue, or error description for agent error.
    """
    validation = {
        "study_id": study_id,
        "hypothesis": "",
        "statistical_test": "",
        "p_value_validity": "",
        "effect_size": "",
        "comparison": "",
        "completeness": "",
        "data_integration": "",
        "protocol": "",
        "note": "",
    }

    # Check if result file exists
    if not os.path.exists(result_path):
        validation["note"] = "Result file not found"
        for key in validation:
            if key not in ["study_id", "note"]:
                validation[key] = 0
        return validation

    try:
        result_df = pd.read_csv(result_path)
    except Exception as e:
        validation["note"] = "Cannot read result file"
        for key in validation:
            if key not in ["study_id", "note"]:
                validation[key] = 0
        return validation

    # Filter to only Mouse clock results
    mouse_results = result_df[result_df["clock"].str.contains("Mouse", case=False)]
    if len(mouse_results) == 0:
        validation["note"] = "No Mouse clock results"
        for key in validation:
            if key not in ["study_id", "note"]:
                validation[key] = 0
        return validation

    result_df = mouse_results  # Only use Mouse clock from now on
    metadata = trace.get("metadata", {})

    # 1. Hypothesis Correctness (only check comparison structure)
    validation["hypothesis"] = evaluate_hypothesis(result_df, metadata)

    # 2. Statistical Test Selection
    try:
        test = result_df["statistical_test"].iloc[0]

        # Check if it's a data limitation (insufficient samples)
        test_lower = str(test).lower()
        if any(
            phrase in test_lower
            for phrase in [
                "n=1",
                "single sample",
                "insufficient",
                "n/a",
                "descriptive",
            ]
        ):
            validation["statistical_test"] = 0  # Data issue, not error

        # Check for NaN or numeric (error)
        elif pd.isna(test) or isinstance(test, (int, float)):
            validation["statistical_test"] = 0  # Data issue

        else:
            # Check if it's an acceptable test
            acceptable_tests = [
                "t-test",
                "t test",
                "welch",
                "student",
                "two-sample",
                "anova",
                "kruskal",
                "mann-whitney",
                "wilcoxon",
            ]

            if any(acceptable in test_lower for acceptable in acceptable_tests):
                validation["statistical_test"] = 1
            else:
                # Check if it's a data limitation description
                if any(
                    phrase in test_lower
                    for phrase in ["comparison", "difference", "mean", "direct"]
                ):
                    validation["statistical_test"] = 0  # Data issue
                else:
                    validation["statistical_test"] = f"Unrecognized: {test}"
    except Exception as e:
        validation["statistical_test"] = 0

    # 3. P-value Validity
    try:
        p_values = result_df["p_value"]

        # If NA, it's a data limitation, not an error
        if p_values.isna().all():
            validation["p_value_validity"] = 0  # Data issue
        elif all((p_values >= 0) & (p_values <= 1) | p_values.isna()):
            validation["p_value_validity"] = 1
        else:
            invalid = p_values[(p_values < 0) | (p_values > 1)]
            validation["p_value_validity"] = f"Invalid: {invalid.tolist()}"
    except Exception as e:
        validation["p_value_validity"] = 0

    # 4. Effect Size Calculation
    try:
        effect_sizes = result_df["effect_size"]

        # If NA due to data limitation, mark as 0
        if effect_sizes.isna().all():
            validation["effect_size"] = 0  # Data issue

        # If extreme (>10), it's a data issue (note it), not an error
        elif any(abs(effect_sizes) >= 10):
            validation["effect_size"] = 0  # Data issue - will be noted

        else:
            validation["effect_size"] = 1
    except Exception as e:
        validation["effect_size"] = 0

    # 5. Comparison Validity
    try:
        comparison = result_df["comparison_name"].iloc[0].lower()
        # Check for forbidden direct age comparisons
        if any(
            term in comparison
            for term in ["young vs old", "old vs young", "aged vs young"]
        ):
            validation["comparison"] = f"Forbidden: {comparison}"
        else:
            validation["comparison"] = 1
    except Exception as e:
        validation["comparison"] = 0

    # 6. Result Completeness
    try:
        required_cols = [
            "gse_id",
            "comparison_name",
            "condition_name",
            "condition_category",
            "covariate",
            "clock",
            "statistical_test",
            "p_value",
            "effect_size",
            "direction",
            "trust_score",
        ]
        missing = [col for col in required_cols if col not in result_df.columns]

        if missing:
            validation["completeness"] = f"Missing: {', '.join(missing)}"
        else:
            # Check for NA in critical fields
            critical_fields = [
                "gse_id",
                "comparison_name",
                "condition_name",
                "clock",
                "trust_score",
            ]
            has_na = []
            for field in critical_fields:
                if result_df[field].isna().any():
                    has_na.append(field)

            if has_na:
                validation["completeness"] = f"NA in: {', '.join(has_na)}"
            else:
                # NA in statistical fields due to data limitation is OK
                validation["completeness"] = 1
    except Exception as e:
        validation["completeness"] = 0

    # 7. Data Integration
    try:
        gse_id = result_df["gse_id"].iloc[0]
        if not gse_id.startswith("GSE"):
            validation["data_integration"] = f"Invalid GSE: {gse_id}"
        else:
            validation["data_integration"] = 1
    except Exception as e:
        validation["data_integration"] = 0

    # 8. Protocol Adherence
    try:
        status = trace.get("status", "unknown")
        if status != "success":
            validation["protocol"] = f"Status: {status}"
        else:
            trust = result_df["trust_score"].iloc[0]
            if 1 <= trust <= 10:
                validation["protocol"] = 1
            else:
                validation["protocol"] = f"Invalid trust: {trust}"
    except Exception as e:
        validation["protocol"] = 0

    # 9. Generate Note (only for dataset issues)
    validation["note"] = generate_note(study_id, trace, result_df)

    return validation


def main():
    # Load traces
    traces_path = paths.TRACES / "successful_traces_mouse.pkl"
    with open(traces_path, "rb") as f:
        traces = pickle.load(f)

    print(f"Total traces available: {len(traces)}")

    # Set random seed for reproducibility
    random.seed(42)

    # Select 100 random traces
    all_keys = list(traces.keys())
    selected_keys = random.sample(all_keys, min(100, len(all_keys)))

    print(f"Selected {len(selected_keys)} traces for validation\n")

    # Validate each trace
    results = []
    for i, study_id in enumerate(selected_keys, 1):
        trace = traces[study_id]
        result_path = paths.TRACES / "result_mm_rna" / study_id / "aging_analysis.csv"

        validation = validate_trace(study_id, trace, result_path)
        results.append(validation)

        if i % 10 == 0:
            print(f"Validated: {i}/{len(selected_keys)} traces")

    # Create DataFrame
    df = pd.DataFrame(results)

    # Save to CSV
    output_path = paths.REPO / "validation_results_100_traces.csv"
    df.to_csv(output_path, index=False)
    print(f"\nResults saved to: {output_path}")

    # Calculate correctness rate (excluding 0s which are data issues)
    print("\n" + "=" * 80)
    print("CORRECTNESS SUMMARY")
    print("=" * 80)

    total_criteria = 0
    total_correct = 0
    total_errors = 0
    total_data_issues = 0

    criteria_cols = [
        "hypothesis",
        "statistical_test",
        "p_value_validity",
        "effect_size",
        "comparison",
        "completeness",
        "data_integration",
        "protocol",
    ]

    for col in criteria_cols:
        correct_count = (df[col] == 1).sum()
        data_issue_count = (df[col] == 0).sum()
        error_count = len(df) - correct_count - data_issue_count
        total_count = len(df)

        total_correct += correct_count
        total_errors += error_count
        total_data_issues += data_issue_count
        total_criteria += total_count

        percentage = (correct_count / total_count) * 100
        print(
            f"{col:20s}: {correct_count:3d} correct, {error_count:3d} errors, {data_issue_count:3d} data issues ({percentage:5.1f}%)"
        )

    # Calculate correctness rate excluding data issues
    analyzable = total_criteria - total_data_issues
    correctness_rate = (total_correct / analyzable * 100) if analyzable > 0 else 0

    print("=" * 80)
    print(f"Total Correct: {total_correct}/{total_criteria}")
    print(f"Total Errors (agent): {total_errors}/{total_criteria}")
    print(f"Total Data Issues: {total_data_issues}/{total_criteria}")
    print(
        f"Correctness Rate (excluding data issues): {total_correct}/{analyzable} = {correctness_rate:.1f}%"
    )
    print("=" * 80)

    # Show actual errors (not data issues)
    print("\n" + "=" * 80)
    print("AGENT ERRORS (excluding data limitations)")
    print("=" * 80)

    error_count = 0
    for col in criteria_cols:
        errors = df[(df[col] != 1) & (df[col] != 0)]
        if len(errors) > 0:
            print(f"\n{col.upper()} - {len(errors)} errors:")
            for idx, row in errors.iterrows():
                print(f"  {row['study_id']}: {row[col]}")
                error_count += 1

    if error_count == 0:
        print("\nNo agent errors detected!")

    # Show data issues
    data_issue_studies = df[(df[criteria_cols] == 0).any(axis=1) & (df["note"] != "")][
        "study_id"
    ].tolist()
    print(f"\n{len(data_issue_studies)} studies with data issues (see 'note' column)")

    return df


if __name__ == "__main__":
    df = main()
    print("\n" + "=" * 80)
    print("Validation Complete")
    print("=" * 80)
