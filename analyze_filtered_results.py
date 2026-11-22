#!/usr/bin/env python3
"""
Analyze agent results WITH FILTERING to generate accurate statistics for manuscript
Following the exact filtering process from s01_visualize_results.py
"""

import pandas as pd
import numpy as np
import glob
import os
from statsmodels.stats.multitest import fdrcorrection

# Load main data
main_df = pd.read_csv('data/main.csv')
scores_df = pd.read_csv('data/scores_summary.csv')

print("="*60)
print("AGENT ANALYSIS STATISTICS WITH FILTERING")
print("="*60)

# 1. BEFORE FILTERING
print("\n1. BEFORE OUTLIER FILTERING")
print("-"*40)
print(f"Total RNA-seq comparisons: {len(main_df):,}")
print(f"Unique datasets (GSE IDs): {main_df['gse_id'].nunique():,}")

# Calculate FDR before filtering
def calculate_fdr(df):
    """Calculate FDR-adjusted p-values"""
    df = df.copy()
    if "pvalue" in df.columns:
        # Replace 0 p-values with smallest non-zero p-value
        mask = ~df["pvalue"].isna()
        min_nonzero_p = df.loc[mask & (df["pvalue"] > 0), "pvalue"].min()
        if pd.isna(min_nonzero_p):
            min_nonzero_p = 1e-16
        df.loc[mask & (df["pvalue"] == 0), "pvalue"] = min_nonzero_p

        # Calculate FDR
        _, fdr = fdrcorrection(df.loc[mask, "pvalue"])
        df.loc[mask, "fdr"] = fdr

        # Replace 0 FDR with smallest non-zero FDR
        min_nonzero_fdr = df.loc[mask & (df["fdr"] > 0), "fdr"].min()
        if pd.isna(min_nonzero_fdr):
            min_nonzero_fdr = 1e-16
        df.loc[mask & (df["fdr"] == 0), "fdr"] = min_nonzero_fdr
    return df

# Apply FDR calculation
combined_df = calculate_fdr(main_df)

# Before filtering statistics
sig_before = (combined_df['fdr'] < 0.05).sum()
print(f"Significant (FDR < 0.05): {sig_before:,} ({sig_before/len(combined_df)*100:.1f}%)")

# Category breakdown before filtering
cat_before = combined_df['condition_category'].value_counts()
print("\nTop categories before filtering:")
for cat, count in cat_before.head(5).items():
    print(f"  {cat}: {count:,} ({count/len(combined_df)*100:.1f}%)")

# 2. APPLY EXACT FILTERING FROM s01_visualize_results.py
print("\n2. OUTLIER FILTERING PROCESS")
print("-"*40)

# Define outlier thresholds for log2FoldChange (using IQR method)
q1_log2fc = combined_df["log2FoldChange"].quantile(0.25)
q3_log2fc = combined_df["log2FoldChange"].quantile(0.75)
iqr_log2fc = q3_log2fc - q1_log2fc

lower_bound_log2fc = q1_log2fc - 10 * iqr_log2fc
upper_bound_log2fc = q3_log2fc + 10 * iqr_log2fc

# Define outlier thresholds for p-values
combined_df["neg_log10_pvalue"] = -np.log10(combined_df["pvalue"])

q1_pval = combined_df["neg_log10_pvalue"].quantile(0.25)
q3_pval = combined_df["neg_log10_pvalue"].quantile(0.75)
iqr_pval = q3_pval - q1_pval

lower_bound_neg_log10_pval = q1_pval - 10 * iqr_pval
upper_bound_neg_log10_pval = q3_pval + 10 * iqr_pval

# Convert back to p-values for filtering
upper_bound_pval = 10 ** (-lower_bound_neg_log10_pval)
lower_bound_pval = 10 ** (-upper_bound_neg_log10_pval)

# Ensure bounds are within valid p-value range
lower_bound_pval = max(0, lower_bound_pval)
upper_bound_pval = min(1, upper_bound_pval)

print(f"Log2FC bounds: [{lower_bound_log2fc:.2f}, {upper_bound_log2fc:.2f}]")
print(f"P-value bounds: [{lower_bound_pval:.2e}, {upper_bound_pval:.2e}]")

# Filter out both p-value and effect size outliers
filtered_df = combined_df[
    (combined_df["log2FoldChange"] >= lower_bound_log2fc) &
    (combined_df["log2FoldChange"] <= upper_bound_log2fc) &
    (combined_df["pvalue"] >= lower_bound_pval) &
    (combined_df["pvalue"] <= upper_bound_pval)
]

print(f"\nData points before filtering: {len(combined_df):,}")
print(f"Data points after filtering: {len(filtered_df):,}")
print(f"Outliers removed: {len(combined_df) - len(filtered_df):,} ({(len(combined_df) - len(filtered_df))/len(combined_df)*100:.1f}%)")

# 3. AFTER FILTERING STATISTICS (FOR MANUSCRIPT)
print("\n3. AFTER FILTERING (MANUSCRIPT NUMBERS)")
print("-"*40)
print(f"Total comparisons analyzed: {len(filtered_df):,}")
print(f"Unique datasets: {filtered_df['gse_id'].nunique():,}")

# Significant findings after filtering
sig_mask = filtered_df['fdr'] < 0.05
print(f"Significant age effects (FDR < 0.05): {sig_mask.sum():,} ({sig_mask.mean()*100:.1f}%)")

# Category breakdown after filtering
print("\n4. INTERVENTION CATEGORIES (AFTER FILTERING)")
print("-"*40)
category_counts = filtered_df['condition_category'].value_counts()
for cat, count in category_counts.head(5).items():
    pct = count/len(filtered_df)*100
    sig_count = (filtered_df[filtered_df['condition_category'] == cat]['fdr'] < 0.05).sum()
    sig_pct = sig_count/count*100 if count > 0 else 0
    print(f"{cat}: {count:,} ({pct:.1f}%) - {sig_count} significant ({sig_pct:.1f}%)")

# Effect size distribution after filtering
print("\n5. EFFECT SIZE DISTRIBUTION (AFTER FILTERING)")
print("-"*40)
effect_sizes = filtered_df['log2FoldChange'].dropna()
print(f"Mean log2 effect size: {effect_sizes.mean():.3f}")
print(f"Median log2 effect size: {effect_sizes.median():.3f}")
print(f"Max positive effect: {effect_sizes.max():.3f}")
print(f"Max negative effect: {effect_sizes.min():.3f}")
print(f"Effect size IQR: [{effect_sizes.quantile(0.25):.3f}, {effect_sizes.quantile(0.75):.3f}]")

# Top drug candidates FROM FILTERED DATA
print("\n6. TOP DRUG CANDIDATES (FROM FILTERED DATA)")
print("-"*40)

# Get drugs from filtered data
drug_df = filtered_df[filtered_df['condition_category'] == 'Drug'].copy()
print(f"Total drug comparisons (after filtering): {len(drug_df):,}")

sig_drugs = drug_df[drug_df['fdr'] < 0.05].copy()
print(f"Significant drug effects: {len(sig_drugs):,} ({len(sig_drugs)/len(drug_df)*100:.1f}%)")

# Find specific drugs in filtered data
drugs_to_find = ['ouabain', 'BMH21', '3M-052', 'SLR14', 'rapamycin', 'metformin']

for drug in drugs_to_find:
    drug_results = filtered_df[filtered_df['condition_name'].str.contains(drug, case=False, na=False)]
    if not drug_results.empty:
        best = drug_results.nsmallest(1, 'fdr').iloc[0]
        print(f"\n{drug.upper()}:")
        print(f"  Log2FC: {best['log2FoldChange']:.3f}")
        print(f"  P-value: {best['pvalue']:.2e}")
        print(f"  FDR: {best['fdr']:.2e}")
        print(f"  GSE: {best['gse_id']}")
    else:
        print(f"\n{drug.upper()}: FILTERED OUT AS OUTLIER")

# Genetic interventions from filtered data
print("\n7. GENETIC INTERVENTIONS (FROM FILTERED DATA)")
print("-"*40)
genetic_df = filtered_df[filtered_df['condition_category'] == 'Genetic'].copy()
print(f"Genetic interventions (after filtering): {len(genetic_df):,}")
sig_genetic = (genetic_df['fdr'] < 0.05).sum()
print(f"Significant genetic effects: {sig_genetic} ({sig_genetic/len(genetic_df)*100:.1f}%)")

# Top genetic by effect size from filtered data
if 'abs_log2FoldChange' not in genetic_df.columns:
    genetic_df['abs_log2FoldChange'] = genetic_df['log2FoldChange'].abs()

top_genetic = genetic_df.nlargest(5, 'abs_log2FoldChange')
print("\nTop 5 genetic interventions by effect size:")
for _, row in top_genetic.iterrows():
    print(f"  {row['condition_name'][:30]:30} | Log2FC: {row['log2FoldChange']:6.2f} | FDR: {row['fdr']:.2e}")

# Summary for manuscript
print("\n" + "="*60)
print("SUMMARY FOR MANUSCRIPT")
print("="*60)
print(f"Initial dataset: {len(combined_df):,} comparisons")
print(f"After filtering (10× IQR): {len(filtered_df):,} comparisons")
print(f"Outliers removed: {len(combined_df) - len(filtered_df):,}")
print(f"Final significant (FDR<0.05): {sig_mask.sum():,} ({sig_mask.mean()*100:.1f}%)")
print(f"\nCategory distribution (filtered):")
for cat, count in category_counts.head(5).items():
    print(f"  {cat}: {count:,} ({count/len(filtered_df)*100:.1f}%)")

# Save filtered statistics
with open('result/filtered_manuscript_statistics.txt', 'w') as f:
    f.write("FILTERED DATA STATISTICS FOR MANUSCRIPT\n")
    f.write("="*60 + "\n\n")
    f.write(f"Initial comparisons: {len(combined_df):,}\n")
    f.write(f"After outlier filtering: {len(filtered_df):,}\n")
    f.write(f"Significant (FDR<0.05): {sig_mask.sum():,} ({sig_mask.mean()*100:.1f}%)\n")
    f.write(f"\nTop categories (after filtering):\n")
    for cat, count in category_counts.head(5).items():
        f.write(f"  - {cat}: {count:,}\n")
    f.write(f"\nEffect size range: {effect_sizes.min():.2f} to {effect_sizes.max():.2f}\n")

print("\nFiltered statistics saved to result/filtered_manuscript_statistics.txt")