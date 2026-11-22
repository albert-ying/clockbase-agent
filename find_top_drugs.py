#!/usr/bin/env python3
"""
Find top drug candidates from the data
"""

import pandas as pd
import numpy as np

# Load data files
main_df = pd.read_csv('data/main.csv')
scores_df = pd.read_csv('data/scores_summary.csv')

print("="*60)
print("TOP DRUG CANDIDATES ANALYSIS")
print("="*60)

# Filter for drugs
drug_df = main_df[main_df['condition_category'] == 'Drug'].copy()
print(f"\nTotal drug comparisons: {len(drug_df):,}")

# Significant drugs
if 'fdr' in drug_df.columns:
    sig_drugs = drug_df[drug_df['fdr'] < 0.05].copy()
else:
    sig_drugs = drug_df[drug_df['pvalue'] < 0.05].copy()

print(f"Significant drug effects: {len(sig_drugs):,}")

# Extract base GSE for matching with scores
scores_df['base_gse'] = scores_df['dataset'].str.extract(r'(GSE\d+)')

# Find drugs with high scores
drug_scores = scores_df[scores_df['intervention'].notna()].copy()

# Filter for actual drug interventions (based on the text patterns)
drug_keywords = ['ouabain', 'BMH21', '3M-052', 'SLR14', 'rapamycin', 'metformin',
                'resveratrol', 'doxycycline', 'insulin', 'tamoxifen', 'dexamethasone',
                'imatinib', 'tretinoin', 'isotretinoin', 'acetaminophen', 'aspirin']

# Look for drugs in intervention names
drug_mask = drug_scores['intervention'].str.lower().str.contains('|'.join(drug_keywords), na=False)
drug_interventions = drug_scores[drug_mask].copy()

print(f"\nDrug interventions with scores: {len(drug_interventions)}")

# Get top drugs by score
if len(drug_interventions) > 0:
    print("\nTOP DRUGS BY FINAL SCORE:")
    print("-"*40)
    top_drugs = drug_interventions.nlargest(20, 'final_score')
    for idx, row in top_drugs.iterrows():
        print(f"{row['intervention'][:30]:30} | Score: {row['final_score']:3.0f} | Dataset: {row['dataset']}")

# Now find the actual effect sizes for top drugs
print("\n\nDRUG EFFECT SIZES FROM MAIN DATA:")
print("-"*40)

# Search for specific drugs in condition_name
drugs_to_find = ['ouabain', 'BMH21', '3M-052', 'SLR14', 'rapamycin', 'metformin']

for drug in drugs_to_find:
    drug_results = main_df[main_df['condition_name'].str.contains(drug, case=False, na=False)]
    if not drug_results.empty:
        # Get the most significant result
        if 'fdr' in drug_results.columns:
            best = drug_results.nsmallest(1, 'fdr').iloc[0]
        else:
            best = drug_results.nsmallest(1, 'pvalue').iloc[0]

        print(f"\n{drug.upper()}:")
        print(f"  GSE ID: {best['gse_id']}")
        print(f"  Condition: {best['condition_name']}")
        print(f"  Log2FC: {best['log2FoldChange']:.3f}")
        print(f"  P-value: {best['pvalue']:.2e}")
        if 'fdr' in best:
            print(f"  FDR: {best['fdr']:.2e}")

        # Try to find the score
        score_match = scores_df[scores_df['intervention'].str.contains(drug, case=False, na=False)]
        if not score_match.empty:
            print(f"  Final Score: {score_match.iloc[0]['final_score']:.0f}")

# Look at all significant drug results
print("\n\nALL SIGNIFICANT DRUGS (TOP 20 BY EFFECT SIZE):")
print("-"*40)

if len(sig_drugs) > 0:
    sig_drugs['abs_effect'] = sig_drugs['log2FoldChange'].abs()
    top_sig = sig_drugs.nlargest(20, 'abs_effect')

    for idx, row in top_sig.iterrows():
        condition = row['condition_name'][:30] if pd.notna(row['condition_name']) else 'Unknown'
        print(f"{condition:30} | Log2FC: {row['log2FoldChange']:7.3f} | P: {row['pvalue']:.2e} | GSE: {row['gse_id']}")

# Look for age-reducing drugs
print("\n\nTOP AGE-REDUCING DRUGS (negative log2FC):")
print("-"*40)

age_reducing = sig_drugs[sig_drugs['log2FoldChange'] < 0].nsmallest(20, 'log2FoldChange')
for idx, row in age_reducing.iterrows():
    condition = row['condition_name'][:30] if pd.notna(row['condition_name']) else 'Unknown'
    print(f"{condition:30} | Log2FC: {row['log2FoldChange']:7.3f} | P: {row['pvalue']:.2e} | GSE: {row['gse_id']}")

# Look for age-accelerating drugs
print("\n\nTOP AGE-ACCELERATING DRUGS (positive log2FC):")
print("-"*40)

age_accelerating = sig_drugs[sig_drugs['log2FoldChange'] > 0].nlargest(20, 'log2FoldChange')
for idx, row in age_accelerating.iterrows():
    condition = row['condition_name'][:30] if pd.notna(row['condition_name']) else 'Unknown'
    print(f"{condition:30} | Log2FC: {row['log2FoldChange']:7.3f} | P: {row['pvalue']:.2e} | GSE: {row['gse_id']}")