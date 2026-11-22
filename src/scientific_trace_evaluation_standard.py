"""
Scientific Trace Evaluation Standard for PhD-Level Scientists

This module provides a comprehensive framework for evaluating aging clock analysis traces
with a focus on scientific rigor, statistical validity, and biological plausibility.

Scoring: Each criterion is scored from 0-10, with weighted importance.
Final score is normalized to 0-100.

Author: AI Agent
Date: November 12, 2025
"""

import sys
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from audit_system import paths

import os
import pickle
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats


@dataclass
class EvaluationCriterion:
    """Individual evaluation criterion with score and justification"""

    name: str
    score: float  # 0-10
    weight: float  # Importance weight
    justification: str
    issues: List[str]
    recommendations: List[str]


class ScientificTraceEvaluator:
    """Comprehensive evaluator for aging clock analysis traces"""

    def __init__(self):
        self.criteria_weights = {
            "hypothesis_formulation": 1.5,
            "experimental_design": 1.5,
            "statistical_methodology": 2.0,
            "effect_size_interpretation": 1.5,
            "sample_size_power": 1.5,
            "confounding_control": 1.5,
            "biological_plausibility": 1.0,
            "reproducibility": 1.0,
            "multiple_testing": 1.0,
            "data_quality": 1.0,
        }

    def evaluate_trace(
        self, study_id: str, trace: Dict, result_df: pd.DataFrame
    ) -> Dict:
        """Comprehensive evaluation of a single trace"""

        evaluation = {
            "study_id": study_id,
            "criteria": {},
            "overall_score": 0.0,
            "grade": "",
            "summary": "",
            "critical_issues": [],
            "strengths": [],
            "recommendations": [],
        }

        # 1. Hypothesis Formulation
        evaluation["criteria"]["hypothesis_formulation"] = (
            self._evaluate_hypothesis_formulation(trace, result_df)
        )

        # 2. Experimental Design Understanding
        evaluation["criteria"]["experimental_design"] = (
            self._evaluate_experimental_design(trace, result_df)
        )

        # 3. Statistical Methodology
        evaluation["criteria"]["statistical_methodology"] = (
            self._evaluate_statistical_methodology(trace, result_df)
        )

        # 4. Effect Size Interpretation
        evaluation["criteria"]["effect_size_interpretation"] = (
            self._evaluate_effect_size_interpretation(trace, result_df)
        )

        # 5. Sample Size and Statistical Power
        evaluation["criteria"]["sample_size_power"] = self._evaluate_sample_size_power(
            trace, result_df
        )

        # 6. Confounding Variable Control
        evaluation["criteria"]["confounding_control"] = (
            self._evaluate_confounding_control(trace, result_df)
        )

        # 7. Biological Plausibility
        evaluation["criteria"]["biological_plausibility"] = (
            self._evaluate_biological_plausibility(trace, result_df)
        )

        # 8. Reproducibility and Transparency
        evaluation["criteria"]["reproducibility"] = self._evaluate_reproducibility(
            trace, result_df
        )

        # 9. Multiple Testing Considerations
        evaluation["criteria"]["multiple_testing"] = self._evaluate_multiple_testing(
            trace, result_df
        )

        # 10. Data Quality Assessment
        evaluation["criteria"]["data_quality"] = self._evaluate_data_quality(
            trace, result_df
        )

        # Calculate overall score
        evaluation["overall_score"] = self._calculate_overall_score(
            evaluation["criteria"]
        )
        evaluation["grade"] = self._assign_grade(evaluation["overall_score"])

        # Compile summary
        evaluation["summary"] = self._generate_summary(evaluation)
        evaluation["critical_issues"] = self._extract_critical_issues(evaluation)
        evaluation["strengths"] = self._extract_strengths(evaluation)
        evaluation["recommendations"] = self._compile_recommendations(evaluation)

        return evaluation

    def _evaluate_hypothesis_formulation(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 1: Hypothesis Formulation and Scientific Question

        Evaluates:
        - Is there a clear, testable hypothesis?
        - Is the comparison biologically meaningful?
        - Are the aging clocks appropriate for the question?
        """
        score = 10.0
        issues = []
        recommendations = []
        justification = ""

        # Check if comparison makes sense
        for _, row in result_df.iterrows():
            comparison = row["comparison_name"]
            condition = row["condition_name"]
            category = row["condition_category"]

            # Check for meaningful hypothesis
            if "vs" not in comparison.lower():
                score -= 2
                issues.append(f"Comparison '{comparison}' lacks clear control group")

            # Check if aging analysis is appropriate
            if category in ["Other", "Unknown"]:
                score -= 1
                issues.append(f"Unclear condition category: {category}")
                recommendations.append(
                    "Clarify condition categorization for biological interpretation"
                )

            # Check for implicit hypothesis testing
            if condition.lower() in ["young", "old", "aging", "aged"]:
                score -= 3
                issues.append(
                    "CRITICAL: Direct age comparison detected - violates protocol"
                )
                recommendations.append(
                    "Do not compare young vs old directly; use clock deviations within age groups"
                )

        justification = (
            f"Hypothesis formulation assessment based on comparison clarity, "
            f"biological relevance, and protocol adherence. "
            f"Comparison: {result_df['comparison_name'].iloc[0]}"
        )

        return EvaluationCriterion(
            name="Hypothesis Formulation",
            score=max(0, score),
            weight=self.criteria_weights["hypothesis_formulation"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_experimental_design(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 2: Experimental Design Understanding

        Evaluates:
        - Correct identification of control vs treatment groups
        - Understanding of experimental conditions
        - Appropriate handling of study design
        """
        score = 10.0
        issues = []
        recommendations = []

        metadata = trace.get("metadata", {})
        title = metadata.get("title", "")
        summary = metadata.get("summary", "")
        design = metadata.get("overall_design", "")

        # Check if design is understood
        if not design or len(design) < 50:
            score -= 2
            issues.append("Limited experimental design information available")

        # Check trust score appropriateness
        trust_score = result_df["trust_score"].iloc[0]

        # Higher trust should correlate with better design
        has_covariates = not pd.isna(result_df["covariate"].iloc[0])

        if has_covariates and trust_score < 5:
            score -= 1
            issues.append(
                "Trust score may be too low for study with controlled covariates"
            )

        if not has_covariates and trust_score > 7:
            score -= 1
            issues.append(
                "Trust score may be too high for study without covariate control"
            )

        # Check condition categorization
        category = result_df["condition_category"].iloc[0]
        condition = result_df["condition_name"].iloc[0].lower()

        category_matches = {
            "Diet": ["diet", "food", "calori", "nutrient", "fast", "restrict"],
            "Drug": ["drug", "compound", "treatment", "inhibitor", "agonist"],
            "Genetic": ["knockout", "mutation", "transgenic", "-/-", "overexpression"],
            "Disease": ["disease", "cancer", "tumor", "pathology", "infection"],
            "Environment": ["stress", "temperature", "hypoxia", "radiation"],
        }

        if category in category_matches:
            if not any(keyword in condition for keyword in category_matches[category]):
                score -= 1.5
                issues.append(
                    f"Condition categorization questionable: '{condition}' as '{category}'"
                )
                recommendations.append("Review condition categorization logic")

        justification = (
            f"Experimental design understanding based on control group identification, "
            f"trust score appropriateness, and condition categorization accuracy."
        )

        return EvaluationCriterion(
            name="Experimental Design Understanding",
            score=max(0, score),
            weight=self.criteria_weights["experimental_design"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_statistical_methodology(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 3: Statistical Methodology

        Evaluates:
        - Appropriateness of statistical test
        - Assumptions met (normality, variance, independence)
        - Correct test selection for data structure
        - Effect size calculation correctness
        """
        score = 10.0
        issues = []
        recommendations = []

        for idx, row in result_df.iterrows():
            test = row["statistical_test"]
            p_value = row["p_value"]
            effect_size = row["effect_size"]
            direction = row["direction"]

            # Check test appropriateness
            if test == "Welch t-test":
                score += 0  # Excellent choice - doesn't assume equal variances
            elif test == "t-test":
                score -= 0.5
                recommendations.append(
                    "Consider Welch t-test instead of standard t-test for robustness"
                )
            elif "ANOVA" in test:
                # ANOVA appropriate for 3+ groups
                score += 0
            else:
                score -= 2
                issues.append(f"Unusual or inappropriate test: {test}")

            # Check p-value validity
            if not (0 <= p_value <= 1):
                score -= 5
                issues.append(f"CRITICAL: Invalid p-value: {p_value}")

            # Check effect size and direction consistency
            expected_direction = 1 if effect_size > 0 else -1
            if direction != expected_direction:
                score -= 5
                issues.append(
                    f"CRITICAL: Direction ({direction}) inconsistent with effect size ({effect_size})"
                )

            # Check for suspicious perfect p-values
            if p_value in [0.0, 1.0]:
                score -= 1
                issues.append(
                    f"Suspicious extreme p-value: {p_value} - check calculation"
                )

            # Check effect size magnitude reasonableness
            if abs(effect_size) > 10:
                score -= 2
                issues.append(
                    f"Unusually large effect size: {effect_size} - potential outlier or error"
                )
                recommendations.append("Investigate extreme effect size values")

        justification = (
            f"Statistical methodology evaluated based on test selection, "
            f"p-value validity, and effect size calculation consistency. "
            f"Test used: {result_df['statistical_test'].iloc[0]}"
        )

        return EvaluationCriterion(
            name="Statistical Methodology",
            score=max(0, score),
            weight=self.criteria_weights["statistical_methodology"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_effect_size_interpretation(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 4: Effect Size Interpretation

        Evaluates:
        - Is effect size reported (not just p-value)?
        - Is effect size magnitude appropriate?
        - Is biological significance considered alongside statistical significance?
        """
        score = 10.0
        issues = []
        recommendations = []

        for idx, row in result_df.iterrows():
            effect_size = row["effect_size"]
            p_value = row["p_value"]

            abs_effect = abs(effect_size)

            # Good: Effect size is reported
            # This gets full marks by default

            # Check for discrepancy between effect size and significance
            if abs_effect > 0.8 and p_value > 0.05:
                issues.append(
                    f"Large effect size ({effect_size:.2f}) but not significant (p={p_value:.3f}) - "
                    f"likely underpowered"
                )
                recommendations.append(
                    "Consider biological significance despite lack of statistical significance"
                )
                score -= 0.5  # Minor deduction, not necessarily wrong

            if abs_effect < 0.2 and p_value < 0.05:
                issues.append(
                    f"Small effect size ({effect_size:.2f}) but statistically significant - "
                    f"may lack biological relevance"
                )
                recommendations.append(
                    "Interpret statistical significance with caution"
                )
                score -= 0.5

            # Check for reporting precision
            # Effect size should be Cohen's d with reasonable precision
            if abs_effect < 0.01 and abs_effect != 0:
                recommendations.append(
                    "Very small effect size - verify biological relevance"
                )

            # Multiple clocks should show consistency
            if idx > 0:
                prev_effect = result_df.iloc[idx - 1]["effect_size"]
                prev_direction = result_df.iloc[idx - 1]["direction"]

                # Check if clocks show opposite directions
                if np.sign(effect_size) != np.sign(prev_effect):
                    score -= 1.5
                    issues.append(
                        f"Clocks show opposite directions: "
                        f"{prev_effect:.2f} vs {effect_size:.2f} - "
                        f"interpret with caution"
                    )
                    recommendations.append(
                        "Investigate biological reasons for clock discordance"
                    )

        justification = (
            f"Effect size interpretation based on magnitude appropriateness, "
            f"consistency with p-values, and biological relevance considerations."
        )

        return EvaluationCriterion(
            name="Effect Size Interpretation",
            score=max(0, score),
            weight=self.criteria_weights["effect_size_interpretation"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_sample_size_power(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 5: Sample Size and Statistical Power

        Evaluates:
        - Is sample size adequate for detecting effects?
        - Are small sample size limitations acknowledged?
        - Is statistical power considered?
        """
        score = 10.0
        issues = []
        recommendations = []

        # Note: We don't have direct sample size in result_df
        # Infer from trust score and context
        trust_score = result_df["trust_score"].iloc[0]

        # Lower trust scores might indicate sample size issues
        if trust_score <= 3:
            score -= 1
            issues.append("Low trust score may indicate sample size limitations")
            recommendations.append("Verify adequate sample size for statistical power")

        # Check for large effect sizes with non-significant p-values
        # This often indicates underpowered studies
        for _, row in result_df.iterrows():
            effect_size = row["effect_size"]
            p_value = row["p_value"]

            if abs(effect_size) > 0.8 and p_value > 0.1:
                score -= 2
                issues.append(
                    f"Large effect ({effect_size:.2f}) but p={p_value:.3f} suggests "
                    f"study is underpowered"
                )
                recommendations.append(
                    "Sample size may be insufficient - consider replication with larger n"
                )

            # Very small p-values with small effects might indicate large sample
            if abs(effect_size) < 0.3 and p_value < 0.001:
                recommendations.append(
                    "Small effect is significant due to large sample size - "
                    "consider biological relevance"
                )

        justification = (
            f"Sample size and power assessment based on effect size/p-value relationships "
            f"and trust score indicators."
        )

        return EvaluationCriterion(
            name="Sample Size and Statistical Power",
            score=max(0, score),
            weight=self.criteria_weights["sample_size_power"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_confounding_control(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 6: Confounding Variable Control

        Evaluates:
        - Are confounding variables identified?
        - Are covariates appropriately controlled?
        - Is age properly handled (critical for aging studies)?
        """
        score = 10.0
        issues = []
        recommendations = []

        has_covariates = not pd.isna(result_df["covariate"].iloc[0])

        if has_covariates:
            covariates = str(result_df["covariate"].iloc[0]).split(",")
            covariates = [c.strip() for c in covariates]

            # Good: Covariates are identified
            score += 0

            # Check if age is included (very important for aging studies)
            if "age" in covariates:
                score += 0  # Excellent
            else:
                score -= 2
                issues.append(
                    "Age not listed as covariate - critical for aging studies"
                )
                recommendations.append(
                    "Ensure samples are age-matched or age is controlled"
                )

            # Check for common confounders
            recommended_covariates = ["age", "sex", "batch", "tissue", "treatment"]
            for cov in recommended_covariates:
                if cov in str(result_df["covariate"].iloc[0]).lower():
                    justification_text = f"Good: {cov} identified as covariate"

        else:
            score -= 3
            issues.append("No covariates identified - potential confounding bias")
            recommendations.append(
                "Identify and control for confounding variables (age, sex, batch, etc.)"
            )

        # Check condition category for inherent confounding
        category = result_df["condition_category"].iloc[0]
        condition = result_df["condition_name"].iloc[0]

        if category == "Disease":
            if not has_covariates or "age" not in str(result_df["covariate"].iloc[0]):
                score -= 1.5
                issues.append(
                    "Disease studies require careful age matching - not clearly controlled"
                )
                recommendations.append(
                    "Verify age matching in disease vs control groups"
                )

        justification = (
            f"Confounding control evaluated based on covariate identification and "
            f"appropriateness for aging studies. "
            f"Covariates: {result_df['covariate'].iloc[0] if has_covariates else 'None'}"
        )

        return EvaluationCriterion(
            name="Confounding Variable Control",
            score=max(0, score),
            weight=self.criteria_weights["confounding_control"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_biological_plausibility(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 7: Biological Plausibility

        Evaluates:
        - Do results make biological sense?
        - Are interpretations consistent with known biology?
        - Are mechanistic considerations appropriate?
        """
        score = 10.0
        issues = []
        recommendations = []

        metadata = trace.get("metadata", {})
        summary = metadata.get("summary", "").lower()

        category = result_df["condition_category"].iloc[0]
        condition = result_df["condition_name"].iloc[0].lower()
        direction = result_df["direction"].iloc[0]

        # Known aging modulators and expected directions
        pro_aging = ["high fat", "western diet", "stress", "radiation", "toxin"]
        anti_aging = [
            "caloric restriction",
            "rapamycin",
            "metformin",
            "exercise",
            "resveratrol",
        ]

        # Check plausibility
        if any(term in condition for term in pro_aging):
            if direction == -1:  # Decelerated aging
                score -= 1.5
                issues.append(
                    f"Unexpected: {condition} shows decelerated aging - "
                    f"typically expected to accelerate aging"
                )
                recommendations.append("Verify result and consider biological context")

        if any(term in condition for term in anti_aging):
            if direction == 1:  # Accelerated aging
                score -= 1.5
                issues.append(
                    f"Unexpected: {condition} shows accelerated aging - "
                    f"typically expected to decelerate aging"
                )
                recommendations.append("Verify result and consider biological context")

        # Check for tissue-specific considerations
        if "brain" in summary or "neuron" in summary:
            recommendations.append(
                "Brain tissue may show different aging patterns - interpret accordingly"
            )

        if "liver" in summary:
            recommendations.append(
                "Liver has high regenerative capacity - aging signals may differ"
            )

        # Check clock consistency with tissue
        clocks_used = result_df["clock"].unique()
        if len(clocks_used) >= 2:
            # Good: Multiple clocks provide validation
            score += 0

        justification = (
            f"Biological plausibility assessed based on known aging modulators, "
            f"tissue-specific considerations, and mechanistic reasoning. "
            f"Condition: {condition}, Direction: {'accelerated' if direction == 1 else 'decelerated'} aging"
        )

        return EvaluationCriterion(
            name="Biological Plausibility",
            score=max(0, score),
            weight=self.criteria_weights["biological_plausibility"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_reproducibility(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 8: Reproducibility and Transparency

        Evaluates:
        - Are methods clearly documented?
        - Are results reproducible from provided information?
        - Is transparency maintained?
        """
        score = 10.0
        issues = []
        recommendations = []

        # Check if all required fields are present
        required_fields = [
            "gse_id",
            "comparison_name",
            "condition_name",
            "condition_category",
            "clock",
            "statistical_test",
            "p_value",
            "effect_size",
            "direction",
            "trust_score",
        ]

        for field in required_fields:
            if field not in result_df.columns or pd.isna(result_df[field].iloc[0]):
                score -= 1
                issues.append(f"Missing or NA value in required field: {field}")

        # Check if study ID is valid
        gse_id = result_df["gse_id"].iloc[0]
        if not gse_id.startswith("GSE"):
            score -= 2
            issues.append(f"Invalid GEO accession format: {gse_id}")

        # Check if methods are specified
        test_method = result_df["statistical_test"].iloc[0]
        if test_method in ["Unknown", "NA", ""]:
            score -= 3
            issues.append("Statistical method not specified")

        # Execution time check
        execution_time = trace.get("execution_time", 0)
        if execution_time < 10:
            score -= 1
            issues.append(
                f"Very short execution time ({execution_time:.1f}s) - "
                f"analysis may be incomplete"
            )
        elif execution_time > 600:  # 10 minutes
            recommendations.append(
                f"Long execution time ({execution_time/60:.1f} min) - "
                f"consider optimization"
            )

        justification = (
            f"Reproducibility assessed based on documentation completeness, "
            f"method specification, and result transparency."
        )

        return EvaluationCriterion(
            name="Reproducibility and Transparency",
            score=max(0, score),
            weight=self.criteria_weights["reproducibility"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_multiple_testing(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 9: Multiple Testing Considerations

        Evaluates:
        - Are multiple comparisons accounted for?
        - Should corrections be applied?
        - Are family-wise error rates considered?
        """
        score = 10.0
        issues = []
        recommendations = []

        num_tests = len(result_df)

        # If testing multiple clocks on same data, this IS multiple testing
        if num_tests > 1:
            # Check if any corrections are applied
            # (We can't directly see this, but can infer from p-values)

            min_p = result_df["p_value"].min()

            if num_tests >= 2:
                score -= 1
                recommendations.append(
                    f"Testing {num_tests} clocks - consider Bonferroni correction "
                    f"(p < {0.05/num_tests:.4f}) or FDR control"
                )

            # If multiple tests and significant results
            if min_p < 0.05:
                bonferroni_threshold = 0.05 / num_tests
                if min_p > bonferroni_threshold:
                    issues.append(
                        f"Significant result (p={min_p:.4f}) does not survive "
                        f"Bonferroni correction (p < {bonferroni_threshold:.4f})"
                    )
                    recommendations.append(
                        "Apply multiple testing correction or interpret with caution"
                    )
                    score -= 1.5

        # Check for different comparisons in same study
        unique_comparisons = result_df["comparison_name"].nunique()
        if unique_comparisons > 1:
            score -= 2
            issues.append(
                f"{unique_comparisons} comparisons made - multiple testing correction required"
            )
            recommendations.append(
                "Apply Benjamini-Hochberg FDR or Bonferroni correction"
            )

        justification = (
            f"Multiple testing considerations based on number of tests ({num_tests}) "
            f"and comparisons ({unique_comparisons if 'unique_comparisons' in locals() else 1})."
        )

        return EvaluationCriterion(
            name="Multiple Testing Considerations",
            score=max(0, score),
            weight=self.criteria_weights["multiple_testing"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _evaluate_data_quality(
        self, trace: Dict, result_df: pd.DataFrame
    ) -> EvaluationCriterion:
        """
        Criterion 10: Data Quality Assessment

        Evaluates:
        - Are there data quality indicators?
        - Are outliers considered?
        - Is data preprocessing appropriate?
        """
        score = 10.0
        issues = []
        recommendations = []

        # Check trust score as proxy for data quality
        trust_score = result_df["trust_score"].iloc[0]

        if trust_score <= 3:
            score -= 2
            issues.append(
                f"Low trust score ({trust_score}) indicates potential data quality concerns"
            )
            recommendations.append(
                "Review data quality metrics and preprocessing steps"
            )
        elif trust_score >= 8:
            score += 0  # High quality, no change needed

        # Check for extreme values that might indicate issues
        for _, row in result_df.iterrows():
            effect_size = row["effect_size"]
            p_value = row["p_value"]

            if abs(effect_size) > 5:
                score -= 1
                issues.append(
                    f"Extreme effect size ({effect_size:.2f}) may indicate outliers or data issues"
                )
                recommendations.append(
                    "Check for outliers and data preprocessing quality"
                )

            if p_value == 0.0 or p_value == 1.0:
                score -= 0.5
                issues.append(f"Extreme p-value ({p_value}) - verify calculation")

        # Check metadata completeness
        metadata = trace.get("metadata", {})
        if not metadata.get("summary") or len(metadata.get("summary", "")) < 100:
            score -= 1
            issues.append("Limited study metadata - may affect interpretation")

        # Check status
        if trace.get("status") != "success":
            score -= 5
            issues.append(f"Analysis status: {trace.get('status')} - not successful")

        justification = (
            f"Data quality assessed based on trust score, extreme value detection, "
            f"and metadata completeness."
        )

        return EvaluationCriterion(
            name="Data Quality Assessment",
            score=max(0, score),
            weight=self.criteria_weights["data_quality"],
            justification=justification,
            issues=issues,
            recommendations=recommendations,
        )

    def _calculate_overall_score(
        self, criteria: Dict[str, EvaluationCriterion]
    ) -> float:
        """Calculate weighted overall score (0-100)"""
        total_weighted_score = 0
        total_weight = 0

        for criterion in criteria.values():
            total_weighted_score += criterion.score * criterion.weight
            total_weight += criterion.weight

        return (total_weighted_score / total_weight) * 10  # Scale to 0-100

    def _assign_grade(self, score: float) -> str:
        """Assign letter grade based on score"""
        if score >= 90:
            return "A (Excellent)"
        elif score >= 80:
            return "B (Good)"
        elif score >= 70:
            return "C (Acceptable)"
        elif score >= 60:
            return "D (Poor)"
        else:
            return "F (Failing)"

    def _generate_summary(self, evaluation: Dict) -> str:
        """Generate summary of evaluation"""
        score = evaluation["overall_score"]
        grade = evaluation["grade"]

        num_issues = sum(len(c.issues) for c in evaluation["criteria"].values())
        num_recommendations = sum(
            len(c.recommendations) for c in evaluation["criteria"].values()
        )

        summary = f"""
Scientific Evaluation Summary for {evaluation['study_id']}

Overall Score: {score:.1f}/100 ({grade})

This analysis received a {grade.split('(')[1].rstrip(')')} rating based on 
comprehensive evaluation of 10 scientific criteria including hypothesis formulation,
statistical methodology, and biological plausibility.

Issues Identified: {num_issues}
Recommendations: {num_recommendations}

Key Strengths: {len(self._extract_strengths(evaluation))}
Critical Issues: {len(self._extract_critical_issues(evaluation))}
"""
        return summary

    def _extract_critical_issues(self, evaluation: Dict) -> List[str]:
        """Extract critical issues (score < 5 on important criteria)"""
        critical = []

        important_criteria = [
            "hypothesis_formulation",
            "statistical_methodology",
            "confounding_control",
        ]

        for criterion_name in important_criteria:
            criterion = evaluation["criteria"].get(criterion_name)
            if criterion and criterion.score < 5:
                critical.extend(criterion.issues)

        # Also add any issue with "CRITICAL" in it
        for criterion in evaluation["criteria"].values():
            for issue in criterion.issues:
                if "CRITICAL" in issue and issue not in critical:
                    critical.append(issue)

        return critical

    def _extract_strengths(self, evaluation: Dict) -> List[str]:
        """Extract strengths (score >= 8)"""
        strengths = []

        for name, criterion in evaluation["criteria"].items():
            if criterion.score >= 8:
                strengths.append(f"{criterion.name}: {criterion.justification}")

        return strengths

    def _compile_recommendations(self, evaluation: Dict) -> List[str]:
        """Compile all recommendations"""
        all_recommendations = []

        for criterion in evaluation["criteria"].values():
            all_recommendations.extend(criterion.recommendations)

        # Remove duplicates while preserving order
        seen = set()
        unique_recommendations = []
        for rec in all_recommendations:
            if rec not in seen:
                seen.add(rec)
                unique_recommendations.append(rec)

        return unique_recommendations

    def print_detailed_report(self, evaluation: Dict):
        """Print detailed evaluation report"""
        print("\n" + "=" * 90)
        print(f"SCIENTIFIC EVALUATION REPORT: {evaluation['study_id']}")
        print("=" * 90)

        print(f"\n{'='*70}")
        print(f"OVERALL ASSESSMENT")
        print(f"{'='*70}")
        print(f"Score: {evaluation['overall_score']:.1f}/100")
        print(f"Grade: {evaluation['grade']}")

        print(f"\n{'='*70}")
        print(f"DETAILED CRITERIA SCORES")
        print(f"{'='*70}\n")

        for name, criterion in evaluation["criteria"].items():
            print(f"{'─'*70}")
            print(f"{criterion.name}")
            print(f"Score: {criterion.score:.1f}/10 (weight: {criterion.weight})")
            print(f"{'─'*70}")
            print(f"Justification: {criterion.justification}\n")

            if criterion.issues:
                print("⚠️  Issues:")
                for issue in criterion.issues:
                    print(f"   • {issue}")
                print()

            if criterion.recommendations:
                print("💡 Recommendations:")
                for rec in criterion.recommendations:
                    print(f"   • {rec}")
                print()

        if evaluation["critical_issues"]:
            print(f"\n{'='*70}")
            print(f"🔴 CRITICAL ISSUES REQUIRING ATTENTION")
            print(f"{'='*70}\n")
            for issue in evaluation["critical_issues"]:
                print(f"• {issue}")

        if evaluation["strengths"]:
            print(f"\n{'='*70}")
            print(f"✅ STRENGTHS")
            print(f"{'='*70}\n")
            for strength in evaluation["strengths"]:
                print(f"• {strength}")

        print(f"\n{'='*70}")
        print(f"📋 SUMMARY RECOMMENDATIONS")
        print(f"{'='*70}\n")
        for rec in evaluation["recommendations"][:10]:  # Top 10
            print(f"• {rec}")


def main():
    """Main evaluation function"""
    # Load traces
    traces_path = paths.TRACES / "successful_traces_mouse.pkl"
    with open(traces_path, "rb") as f:
        traces = pickle.load(f)

    # Set random seed
    random.seed(42)

    # Select same two traces
    all_keys = list(traces.keys())
    selected_keys = random.sample(all_keys, 2)

    # Initialize evaluator
    evaluator = ScientificTraceEvaluator()

    print("=" * 90)
    print("PhD-LEVEL SCIENTIFIC EVALUATION OF AGING CLOCK ANALYSIS TRACES")
    print("=" * 90)
    print(f"\nEvaluating {len(selected_keys)} randomly selected traces")
    print(f"Evaluation Framework: 10 criteria with weighted scoring\n")

    results_summary = []

    for study_id in selected_keys:
        trace = traces[study_id]

        # Load result file
        result_path = paths.TRACES / "result_mm_rna" / study_id / "aging_analysis.csv"

        if not os.path.exists(result_path):
            print(f"\n⚠️  Skipping {study_id}: Result file not found")
            continue

        result_df = pd.read_csv(result_path)

        # Evaluate
        evaluation = evaluator.evaluate_trace(study_id, trace, result_df)

        # Print report
        evaluator.print_detailed_report(evaluation)

        # Store summary
        results_summary.append(
            {
                "study_id": study_id,
                "score": evaluation["overall_score"],
                "grade": evaluation["grade"],
                "critical_issues": len(evaluation["critical_issues"]),
            }
        )

    # Final summary
    print("\n\n" + "=" * 90)
    print("EVALUATION SUMMARY")
    print("=" * 90)
    print("\nStudy Scores:")
    for result in results_summary:
        print(
            f"  {result['study_id']}: {result['score']:.1f}/100 ({result['grade']}) - "
            f"{result['critical_issues']} critical issues"
        )

    avg_score = np.mean([r["score"] for r in results_summary])
    print(f"\nAverage Score: {avg_score:.1f}/100")

    if avg_score >= 80:
        overall = "EXCELLENT - High scientific rigor"
    elif avg_score >= 70:
        overall = "GOOD - Acceptable scientific standards"
    elif avg_score >= 60:
        overall = "ACCEPTABLE - Some concerns present"
    else:
        overall = "NEEDS IMPROVEMENT - Significant issues identified"

    print(f"Overall Assessment: {overall}")


if __name__ == "__main__":
    main()
