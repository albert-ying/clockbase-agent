#!/usr/bin/env python3
"""
ANALYSIS 4 (combined modalities): FDA-Approved Drug Clinical Validation
========================================================================

Pools mouse RNA-seq + human RNA-seq + human DNAm into a single significant-drug
table, applies the same per-modality SCORE_FLOOR as the combined volcano:
    mouse RNA-seq : 100
    human RNA-seq : 100
    human DNAm    :  80

Then runs the full OpenTargets matching + clinical-stage classification +
FDA anti-aging / pro-aging / clinical-stage panels, writing combined outputs.

This script imports the heavy machinery (matching, classification, plotting)
from `analysis_4_fda_approved_drugs.py` so we do not duplicate the ~1800-line
matcher. We override only:
  - `load_data`            : pulls from the combined 3-modality CSV
  - figure save paths      : monkey-patch matplotlib.pyplot.savefig wrapper
  - csv save paths         : redirect via PATH_REWRITE table

Outputs:
  figures/analysis_4_fda_approved_drugs_combined.png   (panel: clinical-stage barplot)
  figures/analysis_4_fda_approved_drugs_combined.pdf
  figures/analysis_4_fda_antiaging_drugs_combined.{png,pdf}
  figures/analysis_4_fda_proaging_drugs_combined.{png,pdf}
  figures/analysis_4_drug_clinical_density_combined.{png,pdf}
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Import heavy machinery from the original mouse-only analysis_4
import analysis_4_fda_approved_drugs as a4

SCORE_FLOOR = {"mouse RNA-seq": 100, "human RNA-seq": 100, "human DNAm": 80}

# Map original output filenames -> combined filenames (used by both savefig
# and to_csv monkey-patches).
PATH_REWRITE = {
    # CSV outputs
    "result/analysis_4_drug_clinical_status.csv":              "result/analysis_4_drug_clinical_status_combined.csv",
    "result/analysis_4_fda_antiaging_drugs.csv":               "result/analysis_4_fda_antiaging_drugs_combined.csv",
    "result/analysis_4_fda_proaging_drugs.csv":                "result/analysis_4_fda_proaging_drugs_combined.csv",
    "result/analysis_4_drug_pairwise_tests.csv":               "result/analysis_4_drug_pairwise_tests_combined.csv",
    "result/analysis_4_drug_pairwise_fold_change.csv":         "result/analysis_4_drug_pairwise_fold_change_combined.csv",
    "result/analysis_4_categories_vs_background_fisher.csv":   "result/analysis_4_categories_vs_background_fisher_combined.csv",
    "result/analysis_4_fuzzy_matches.csv":                     "result/analysis_4_fuzzy_matches_combined.csv",
    # Figure outputs (one entry per stem; both .pdf and .png handled below)
    "figures/analysis_4_drug_clinical_barplot.pdf":            "figures/analysis_4_fda_approved_drugs_combined.pdf",
    "figures/analysis_4_drug_clinical_barplot.png":            "figures/analysis_4_fda_approved_drugs_combined.png",
    "figures/analysis_4_drug_clinical_density.pdf":            "figures/analysis_4_drug_clinical_density_combined.pdf",
    "figures/analysis_4_drug_clinical_density.png":            "figures/analysis_4_drug_clinical_density_combined.png",
    "figures/analysis_4_fda_antiaging_drugs.pdf":              "figures/analysis_4_fda_antiaging_drugs_combined.pdf",
    "figures/analysis_4_fda_antiaging_drugs.png":              "figures/analysis_4_fda_antiaging_drugs_combined.png",
    "figures/analysis_4_fda_proaging_drugs.pdf":               "figures/analysis_4_fda_proaging_drugs_combined.pdf",
    "figures/analysis_4_fda_proaging_drugs.png":               "figures/analysis_4_fda_proaging_drugs_combined.png",
}


def combined_load_data():
    """Replacement for a4.load_data: pool 3 modalities. NO SCORE_FLOOR (matches mouse-only behavior; FDA matching downstream uses fdr<0.05 + Drug filter)."""
    print("Loading combined 3-modality data...")
    df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)
    if "outlier_status" in df.columns:
        df = df[df["outlier_status"] == "Normal"].copy()
    print(f"  rows after outlier filter: {len(df):,}")
    print(df["modality"].value_counts().to_string())

    # The downstream matcher/plotter uses these columns: condition_name,
    # condition_category, fdr, log2FoldChange. Ensure dtype + presence.
    needed = {"condition_name", "condition_category", "fdr", "log2FoldChange"}
    missing = needed - set(df.columns)
    if missing:
        raise RuntimeError(f"missing required columns in combined CSV: {missing}")

    # Try to load DrugAge (used only for descriptive printout)
    try:
        drugage_df = pd.read_csv("data/agent_result_drugage_significant.csv")
        print(f"Loaded DrugAge database: {len(drugage_df)} compounds")
    except Exception:
        drugage_df = None

    # Load OpenTargets drug molecule data
    try:
        opentargets_df = pd.read_parquet(
            "data/drug_molecule/part-00000-35323264-2241-4256-b643-3a92cd0230e5-c000.snappy.parquet"
        )
        print(f"Loaded OpenTargets drug data: {len(opentargets_df):,} molecules")
    except Exception as e:
        print(f"Warning: could not load OpenTargets data: {e}")
        opentargets_df = None

    print(f"Loaded {len(df):,} combined comparisons")
    return df, drugage_df, opentargets_df


# ---- monkey-patch savefig & to_csv to redirect to combined paths ----
_orig_savefig = plt.savefig
_orig_figure_savefig = plt.Figure.savefig
_orig_df_to_csv = pd.DataFrame.to_csv


def _rewrite_path(p):
    s = str(p)
    return PATH_REWRITE.get(s, s)


def _patched_savefig(fname, *args, **kwargs):
    return _orig_savefig(_rewrite_path(fname), *args, **kwargs)


def _patched_figure_savefig(self, fname, *args, **kwargs):
    return _orig_figure_savefig(self, _rewrite_path(fname), *args, **kwargs)


def _patched_to_csv(self, path_or_buf=None, *args, **kwargs):
    if isinstance(path_or_buf, (str, Path)):
        path_or_buf = _rewrite_path(path_or_buf)
    return _orig_df_to_csv(self, path_or_buf, *args, **kwargs)


def main():
    print("\n" + "=" * 80)
    print("FDA-APPROVED DRUG CLINICAL VALIDATION (combined modalities)")
    print("=" * 80)

    # Patch matplotlib + pandas IO
    plt.savefig = _patched_savefig
    plt.Figure.savefig = _patched_figure_savefig
    pd.DataFrame.to_csv = _patched_to_csv

    # Patch the loader on the imported module
    a4.load_data = combined_load_data

    # Run the original pipeline (now sourced from combined data, writing to
    # combined output paths)
    main_df, drugage_df, opentargets_df = a4.load_data()
    drug_df, results_df, overall_baseline, lookup = a4.classify_drugs(
        main_df, drugage_df, opentargets_df
    )
    print(results_df)
    fda_anti, fda_pro = a4.identify_fda_drugs(drug_df)
    chi2, p_chi2, comp_df = a4.perform_statistical_tests(drug_df, results_df)
    a4.create_pairwise_fold_change_table(results_df)
    a4.test_categories_vs_background(main_df, opentargets_df, lookup)
    a4.create_figure(drug_df, results_df, comp_df, fda_anti, fda_pro, overall_baseline)
    a4.print_summary(results_df, chi2, p_chi2, fda_anti)
    print("\nDone.")


if __name__ == "__main__":
    main()
