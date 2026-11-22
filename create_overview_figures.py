#!/usr/bin/env python3
"""
Create Nature publication-quality overview figures for intervention categories
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from aquarel import load_theme
from scipy import stats
from adjustText import adjust_text
import matplotlib.patches as mpatches
import matplotlib.patheffects as pe

# Load and apply theme
theme = (
    load_theme("scientific")
    .set_grid(draw=False)
    .set_font(size=15)
    .set_ticks(direction="out")
    .set_axis_labels(pad=10)
)
theme.apply()

# Load data
df = pd.read_csv('figures/combined_volcano_plot_data.csv')

# Filter for significant results only
df = df[df['is_significant'] == True].copy()

# Rename column for consistency
df['category'] = df['condition_category']

# Map direction values to text
df['direction'] = df['direction'].map({1.0: 'positive', -1.0: 'negative', 0.0: 'neutral'})

# Define category colors
category_colors = {
    'Drug': '#2E86AB',
    'Disease': '#A23B72',
    'Genetic': '#F18F01',
    'Environment': '#C73E1D',
    'Other': '#6B7A8F'
}

# Figure 1: Bar plot showing number of interventions per category
fig, ax = plt.subplots(figsize=(10, 8))

# Count interventions per category
category_counts = df['category'].value_counts()
categories = category_counts.index.tolist()
counts = category_counts.values.tolist()

# Create bars with custom colors
bars = ax.bar(categories, counts,
               color=[category_colors.get(cat, '#6B7A8F') for cat in categories],
               edgecolor='black', linewidth=2, alpha=0.9)

# Add value labels on bars
for bar, count in zip(bars, counts):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
            f'{int(count)}',
            ha='center', va='bottom', fontsize=14, fontweight='bold')

ax.set_xlabel('Intervention Category', fontsize=16, fontweight='bold')
ax.set_ylabel('Number of Interventions', fontsize=16, fontweight='bold')
ax.set_title('Distribution of Interventions Across Categories', fontsize=18, fontweight='bold', pad=20)

# Customize axes
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0, max(counts) * 1.1)

# Apply theme transformations and save
theme.apply_transforms()
plt.tight_layout()
plt.savefig('figures/intervention_counts_by_category.pdf', bbox_inches='tight')
plt.savefig('figures/intervention_counts_by_category.png', dpi=600, bbox_inches='tight')
plt.close()

# Figure 2: Aging direction analysis with statistical test
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Count accelerating vs decelerating per category
aging_direction_data = []
for category in df['category'].unique():
    cat_df = df[df['category'] == category]
    accelerating = len(cat_df[cat_df['direction'] == 'positive'])
    decelerating = len(cat_df[cat_df['direction'] == 'negative'])
    total = len(cat_df)
    aging_direction_data.append({
        'Category': category,
        'Accelerating': accelerating,
        'Decelerating': decelerating,
        'Total': total,
        'Proportion Accelerating': accelerating/total if total > 0 else 0,
        'Proportion Decelerating': decelerating/total if total > 0 else 0
    })

direction_df = pd.DataFrame(aging_direction_data).sort_values('Total', ascending=False)

# Stacked bar plot
categories = direction_df['Category'].tolist()
accelerating = direction_df['Accelerating'].tolist()
decelerating = direction_df['Decelerating'].tolist()

x = np.arange(len(categories))
width = 0.6

# Create stacked bars with outline
bars1 = ax1.bar(x, decelerating, width, label='Decelerating Aging',
                color='#4A90E2', edgecolor='black', linewidth=2, alpha=0.9)
bars2 = ax1.bar(x, accelerating, width, bottom=decelerating, label='Accelerating Aging',
                color='#E94B3C', edgecolor='black', linewidth=2, alpha=0.9)

# Add percentage labels
for i, (dec, acc, total) in enumerate(zip(decelerating, accelerating, direction_df['Total'].tolist())):
    # Decelerating percentage
    if dec > 0:
        ax1.text(i, dec/2, f'{dec/total*100:.1f}%',
                ha='center', va='center', fontweight='bold', color='white', fontsize=12)
    # Accelerating percentage
    if acc > 0:
        ax1.text(i, dec + acc/2, f'{acc/total*100:.1f}%',
                ha='center', va='center', fontweight='bold', color='white', fontsize=12)

ax1.set_xlabel('Intervention Category', fontsize=16, fontweight='bold')
ax1.set_ylabel('Number of Interventions', fontsize=16, fontweight='bold')
ax1.set_title('Aging Direction by Category', fontsize=18, fontweight='bold', pad=20)
ax1.set_xticks(x)
ax1.set_xticklabels(categories)
ax1.legend(loc='upper right', frameon=False, fontsize=12)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Statistical analysis - Chi-square test
contingency_table = direction_df[['Accelerating', 'Decelerating']].values
chi2, p_value, dof, expected = stats.chi2_contingency(contingency_table)

# Proportion plot with statistical annotation
props_acc = direction_df['Proportion Accelerating'].tolist()
props_dec = direction_df['Proportion Decelerating'].tolist()

x2 = np.arange(len(categories))
width2 = 0.35

bars3 = ax2.bar(x2 - width2/2, props_dec, width2, label='Decelerating',
                color='#4A90E2', edgecolor='black', linewidth=2, alpha=0.9)
bars4 = ax2.bar(x2 + width2/2, props_acc, width2, label='Accelerating',
                color='#E94B3C', edgecolor='black', linewidth=2, alpha=0.9)

# Add significance indicators
ax2.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5, linewidth=2)

ax2.set_xlabel('Intervention Category', fontsize=16, fontweight='bold')
ax2.set_ylabel('Proportion of Interventions', fontsize=16, fontweight='bold')
ax2.set_title(f'Aging Direction Proportions\n(χ² = {chi2:.2f}, p = {p_value:.3f})',
              fontsize=18, fontweight='bold', pad=20)
ax2.set_xticks(x2)
ax2.set_xticklabels(categories)
ax2.set_ylim(0, 1.0)
ax2.legend(loc='upper right', frameon=False, fontsize=12)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Apply theme transformations
theme.apply_transforms()
plt.tight_layout()
plt.savefig('figures/aging_direction_analysis.pdf', bbox_inches='tight')
plt.savefig('figures/aging_direction_analysis.png', dpi=600, bbox_inches='tight')
plt.close()

# Figure 3: Effect size distributions by category
fig, axes = plt.subplots(2, 2, figsize=(16, 14))
axes = axes.flatten()

for idx, category in enumerate(df['category'].unique()[:4]):
    ax = axes[idx]
    cat_data = df[df['category'] == category]['log2FoldChange']

    # Create histogram with KDE overlay
    counts, bins, patches = ax.hist(cat_data, bins=20, alpha=0.7,
                                    color=category_colors.get(category, '#6B7A8F'),
                                    edgecolor='black', linewidth=1.5)

    # Add KDE curve
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(cat_data)
    x_range = np.linspace(cat_data.min(), cat_data.max(), 200)
    kde_values = kde(x_range)

    # Scale KDE to match histogram
    kde_scaled = kde_values * len(cat_data) * (bins[1] - bins[0])
    ax.plot(x_range, kde_scaled, color='black', linewidth=2.5, alpha=0.8)

    # Add vertical line at zero
    ax.axvline(x=0, color='red', linestyle='--', linewidth=2, alpha=0.7)

    # Add median line
    median_val = cat_data.median()
    ax.axvline(x=median_val, color='navy', linestyle='-', linewidth=2, alpha=0.7,
              label=f'Median: {median_val:.2f}')

    # Statistics annotation
    mean_val = cat_data.mean()
    std_val = cat_data.std()
    skew_val = cat_data.skew()

    stats_text = (f'n = {len(cat_data)}\n'
                 f'Mean = {mean_val:.2f}\n'
                 f'SD = {std_val:.2f}\n'
                 f'Skew = {skew_val:.2f}')

    ax.text(0.95, 0.95, stats_text, transform=ax.transAxes,
           fontsize=11, verticalalignment='top', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='black'))

    ax.set_xlabel('log₂ Fold Change', fontsize=14, fontweight='bold')
    ax.set_ylabel('Frequency', fontsize=14, fontweight='bold')
    ax.set_title(f'{category} Interventions', fontsize=16, fontweight='bold')
    ax.legend(loc='upper left', frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add shaded regions for effect direction
    y_max = ax.get_ylim()[1]
    ax.fill_between([cat_data.min(), 0], 0, y_max, alpha=0.1, color='blue', label='Decelerating')
    ax.fill_between([0, cat_data.max()], 0, y_max, alpha=0.1, color='red', label='Accelerating')

# Overall title
fig.suptitle('Distribution of Effect Sizes Across Intervention Categories',
            fontsize=20, fontweight='bold', y=1.02)

# Apply theme transformations
theme.apply_transforms()
plt.tight_layout()
plt.savefig('figures/effect_size_distributions.pdf', bbox_inches='tight')
plt.savefig('figures/effect_size_distributions.png', dpi=600, bbox_inches='tight')
plt.close()

# Figure 4: Violin plot comparing effect sizes across categories
fig, ax = plt.subplots(figsize=(14, 10))

# Prepare data for violin plot
categories_order = df.groupby('category')['log2FoldChange'].median().sort_values().index.tolist()

# Create violin plot with scatter overlay
parts = ax.violinplot([df[df['category'] == cat]['log2FoldChange'].values
                       for cat in categories_order],
                      positions=range(len(categories_order)),
                      widths=0.7,
                      showmeans=False,
                      showmedians=False,
                      showextrema=False)

# Customize violin colors
for i, (pc, cat) in enumerate(zip(parts['bodies'], categories_order)):
    pc.set_facecolor(category_colors.get(cat, '#6B7A8F'))
    pc.set_edgecolor('black')
    pc.set_linewidth(2)
    pc.set_alpha(0.7)

# Add scatter points with outline
for i, cat in enumerate(categories_order):
    cat_data = df[df['category'] == cat]['log2FoldChange'].values
    y = np.random.normal(i, 0.04, size=len(cat_data))

    # Outline
    ax.scatter(cat_data, y, s=80, c='black', alpha=0.8, zorder=1)
    # Data points
    ax.scatter(cat_data, y, s=50, c=category_colors.get(cat, '#6B7A8F'),
              alpha=0.6, zorder=2)

# Add median lines
for i, cat in enumerate(categories_order):
    median_val = df[df['category'] == cat]['log2FoldChange'].median()
    ax.plot([median_val, median_val], [i-0.35, i+0.35],
           color='white', linewidth=3, zorder=3)
    ax.plot([median_val, median_val], [i-0.35, i+0.35],
           color='black', linewidth=2, zorder=4)

# Customize axes
ax.axvline(x=0, color='red', linestyle='--', linewidth=2.5, alpha=0.7, zorder=0)
ax.set_yticks(range(len(categories_order)))
ax.set_yticklabels(categories_order, fontsize=14)
ax.set_xlabel('log₂ Fold Change (Age Acceleration)', fontsize=16, fontweight='bold')
ax.set_ylabel('Intervention Category', fontsize=16, fontweight='bold')
ax.set_title('Effect Size Distribution Across Intervention Categories',
            fontsize=18, fontweight='bold', pad=20)

# Add annotations for aging directions
ax.text(0.02, 0.98, '← Decelerating Aging', transform=ax.transAxes,
       fontsize=12, color='blue', fontweight='bold',
       verticalalignment='top', horizontalalignment='left')
ax.text(0.98, 0.98, 'Accelerating Aging →', transform=ax.transAxes,
       fontsize=12, color='red', fontweight='bold',
       verticalalignment='top', horizontalalignment='right')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3, linestyle=':', linewidth=1)

# Apply theme transformations
theme.apply_transforms()
plt.tight_layout()
plt.savefig('figures/effect_size_violin_plot.pdf', bbox_inches='tight')
plt.savefig('figures/effect_size_violin_plot.png', dpi=600, bbox_inches='tight')
plt.close()

# Create summary statistics table
summary_stats = []
for category in df['category'].unique():
    cat_df = df[df['category'] == category]
    summary_stats.append({
        'Category': category,
        'Total Interventions': len(cat_df),
        'Accelerating': len(cat_df[cat_df['direction'] == 'positive']),
        'Decelerating': len(cat_df[cat_df['direction'] == 'negative']),
        'Mean Effect Size': cat_df['log2FoldChange'].mean(),
        'Median Effect Size': cat_df['log2FoldChange'].median(),
        'SD Effect Size': cat_df['log2FoldChange'].std(),
        'Min Effect Size': cat_df['log2FoldChange'].min(),
        'Max Effect Size': cat_df['log2FoldChange'].max()
    })

summary_df = pd.DataFrame(summary_stats)
summary_df = summary_df.round(2)
summary_df.to_csv('figures/intervention_summary_statistics.csv', index=False)

print("Figures created successfully!")
print("\nSummary Statistics:")
print(summary_df.to_string())
print(f"\nChi-square test for category differences:")
print(f"χ² = {chi2:.3f}, p-value = {p_value:.4f}")

# Perform ANOVA to test if effect sizes differ between categories
categories_list = []
effects_list = []
for category in df['category'].unique():
    cat_effects = df[df['category'] == category]['log2FoldChange'].values
    categories_list.append([category] * len(cat_effects))
    effects_list.append(cat_effects)

all_categories = np.concatenate(categories_list)
all_effects = np.concatenate(effects_list)

# One-way ANOVA
f_stat, anova_p = stats.f_oneway(*[df[df['category'] == cat]['log2FoldChange'].values
                                   for cat in df['category'].unique()])

print(f"\nOne-way ANOVA for effect size differences:")
print(f"F-statistic = {f_stat:.3f}, p-value = {anova_p:.4f}")

# Post-hoc pairwise t-tests if significant
if anova_p < 0.05:
    print("\nPost-hoc pairwise t-tests (Bonferroni corrected):")
    categories = df['category'].unique()
    n_comparisons = len(categories) * (len(categories) - 1) / 2
    alpha_corrected = 0.05 / n_comparisons

    for i, cat1 in enumerate(categories):
        for cat2 in categories[i+1:]:
            data1 = df[df['category'] == cat1]['log2FoldChange'].values
            data2 = df[df['category'] == cat2]['log2FoldChange'].values
            t_stat, p_val = stats.ttest_ind(data1, data2)
            sig = '*' if p_val < alpha_corrected else ''
            print(f"  {cat1} vs {cat2}: t = {t_stat:.3f}, p = {p_val:.4f} {sig}")