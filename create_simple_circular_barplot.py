#!/usr/bin/env python3
"""
Create Nature publication-quality circular barplot showing intervention counts by category
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle
from textwrap import wrap

# Load data
df = pd.read_csv('figures/combined_volcano_plot_data.csv')

# Filter for significant results only
df = df[df['is_significant'] == True].copy()

# Count interventions per category
category_counts = df['condition_category'].value_counts()

# Sort categories for better visualization
categories_order = ['Genetic', 'Other', 'Drug', 'Disease', 'Environment']
categories = []
counts = []
for cat in categories_order:
    if cat in category_counts.index:
        categories.append(cat)
        counts.append(category_counts[cat])

# Set up the data for circular plot
N = len(categories)
ANGLES = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
counts = counts
categories = categories

# Store original values for plotting
counts_plot = counts.copy()
categories_plot = categories.copy()

# Make the plot close properly for polar coordinates
counts = counts + [counts[0]]
ANGLES = ANGLES + [ANGLES[0]]

# Define colors for each category
category_colors = {
    'Drug': '#2E86AB',
    'Disease': '#A23B72',
    'Genetic': '#F18F01',
    'Environment': '#C73E1D',
    'Other': '#6B7A8F'
}

# Get colors for bars
colors = [category_colors[cat] for cat in categories_plot]

# Create the figure
fig, ax = plt.subplots(figsize=(12, 12), subplot_kw=dict(projection='polar'))

# Set background
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

# Plot the bars
bars = ax.bar(ANGLES[:-1], counts_plot, width=2*np.pi/N - 0.01,
               bottom=0, color=colors, edgecolor='black', linewidth=2.5,
               alpha=0.85)

# Set the labels
ax.set_xticks(ANGLES[:-1])
ax.set_xticklabels(categories_plot, size=16, fontweight='bold')

# Add value labels on bars
for angle, count, cat in zip(ANGLES[:-1], counts_plot, categories_plot):
    # Calculate label position
    label_angle = angle
    if label_angle > np.pi:
        label_angle = label_angle - np.pi
        ha = 'right'
    else:
        ha = 'left'

    # Add count value
    ax.text(angle, count + 100, str(int(count)),
            ha='center', va='bottom', size=14, fontweight='bold')

# Customize the plot
ax.set_ylim(0, max(counts_plot) * 1.15)
ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)

# Add gridlines
ax.set_rticks([500, 1000, 1500, 2000])
ax.set_rlabel_position(0)
ax.tick_params(axis='y', labelsize=11, colors='gray')
ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.8)

# Remove radial axis label
ax.set_ylabel('')

# Add title
plt.title('Intervention Categories Distribution\n5,763 Significant Interventions',
          pad=30, fontsize=20, fontweight='bold')

# Add a subtle center circle
center_circle = Circle((0, 0), 0.25, transform=ax.transData._b,
                       color='white', edgecolor='black', linewidth=2, zorder=10)
ax.add_patch(center_circle)

# Save the figure
plt.tight_layout()
plt.savefig('figures/simple_circular_barplot_categories.pdf', dpi=300, bbox_inches='tight')
plt.savefig('figures/simple_circular_barplot_categories.png', dpi=600, bbox_inches='tight')
plt.show()

print("Simple circular barplot created successfully!")
print("\nCategory counts:")
for cat, count in zip(categories_plot, counts_plot):
    print(f"{cat}: {count}")

# Create an alternative cleaner version
fig2, ax2 = plt.subplots(figsize=(10, 10), subplot_kw=dict(projection='polar'))

# Set background
fig2.patch.set_facecolor("white")
ax2.set_facecolor("white")

# Adjust angles to add gaps between bars
width = 2 * np.pi / N - 0.15  # Smaller width for gaps
offset = 0.075  # Half of the gap

# Plot bars with gaps
for i, (angle, count, cat) in enumerate(zip(ANGLES[:-1], counts_plot, categories_plot)):
    ax2.bar(angle, count, width=width, bottom=200,
            color=category_colors[cat], edgecolor='white', linewidth=3,
            alpha=0.9)

    # Add category labels
    ax2.text(angle, 100, cat.upper(), ha='center', va='center',
             fontsize=13, fontweight='bold', color='black')

    # Add count values
    ax2.text(angle, count + 350, f'{count:,}', ha='center', va='bottom',
             fontsize=15, fontweight='bold')

# Clean up the plot
ax2.set_ylim(0, max(counts_plot) + 500)
ax2.set_theta_offset(np.pi / 2)
ax2.set_theta_direction(-1)

# Remove all grids and ticks
ax2.set_rticks([])
ax2.set_xticks([])
ax2.grid(False)
ax2.spines['polar'].set_visible(False)
ax2.set_facecolor('#F8F8F8')

# Add concentric circles for reference
for r in [500, 1000, 1500, 2000]:
    if r <= max(counts_plot):
        circle = plt.Circle((0, 0), r, transform=ax2.transData._b,
                           fill=False, edgecolor='gray', linewidth=0.5,
                           alpha=0.3, linestyle='--')
        ax2.add_artist(circle)
        # Add value labels for circles
        if r > 0:
            ax2.text(np.pi/2, r, str(r), ha='center', va='center',
                    fontsize=9, color='gray', alpha=0.7)

# Add title with larger font
ax2.text(0, 1.15, 'Aging Intervention Landscape',
         transform=ax2.transAxes, ha='center', va='center',
         fontsize=24, fontweight='bold')

ax2.text(0, 1.08, f'Distribution of {sum(counts_plot):,} Significant Interventions',
         transform=ax2.transAxes, ha='center', va='center',
         fontsize=14, style='italic', color='#444444')

# Add center decoration
center_circle2 = Circle((0, 0), 180, transform=ax2.transData._b,
                        color='white', edgecolor='none', zorder=10)
ax2.add_patch(center_circle2)

# Save the cleaner version
plt.tight_layout()
plt.savefig('figures/clean_circular_barplot_categories.pdf', dpi=300, bbox_inches='tight')
plt.savefig('figures/clean_circular_barplot_categories.png', dpi=600, bbox_inches='tight')
plt.show()

print("\nClean version created successfully!")