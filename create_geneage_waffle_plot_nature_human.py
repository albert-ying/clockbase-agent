#!/usr/bin/env python3
"""
Create Nature publication-quality waffle plot showing replicability of genetic interventions
with the GenAge database - high-quality interventions (trust score > 80)
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats

warnings.filterwarnings("ignore")

# Configuration flag: Set to True to merge overlapping entries using majority voting
USE_MERGED = True  # Set to True to merge duplicate genes by majority vote

# Apply consistent theme
from aquarel import load_theme

# ---- HUMAN CONFIG --------------------------------------------------
MODALITY = "rna"  # "rna" or "dnam"  — change and rerun for the other modality
METHOD_NAME = f"hs_{MODALITY}"
SUFFIX = f"_{METHOD_NAME}"
# --------------------------------------------------------------------


theme = (
    load_theme("scientific")
    .set_grid(draw=False)
    .set_font(size=14)
    .set_ticks(direction="out")
    .set_axis_labels(pad=10)
)
theme.apply()

# Set white background for Nature publication
plt.rcParams["savefig.transparent"] = False
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"
plt.rcParams["font.family"] = "Arial"

# Load data
df = pd.read_csv(f"data/agent_result_genage_{METHOD_NAME}.csv")

print("\n" + "=" * 60)
print("Data Preprocessing")
print("=" * 60)
print(f"Original dataset size: {len(df)}")

# Remove exact duplicates
df = df.drop_duplicates()
print(f"After removing duplicates: {len(df)}")

# Filter by trust score > 80
df = df[df["trust_score"] > 80]
print(f"After filtering trust_score > 80: {len(df)}")

# Analyze consistency
print("\n" + "=" * 60)
print("GenAge Database Validation Analysis")
print("=" * 60)

# Count matched genes (genage_annotated = True)
matched_genes = df[df["genage_annotated"] == True].copy()
print(f"\nTotal interventions in dataset: {len(df)}")
print(f"Matched with GenAge: {len(matched_genes)}")

# Among matched genes, check agreement
# Agreement: both genage_annotation and inferred_influence are the same
matched_genes["agreement"] = (
    matched_genes["genage_annotation"] == matched_genes["inferred_influence"]
)

# Handle merging if USE_MERGED is True
if USE_MERGED:
    print("\n=== Using MERGED mode (majority voting) ===")
    # Group by gene and use majority voting for inferred_influence
    gene_groups = matched_genes.groupby("genage_matched_genes")

    merged_results = []
    for gene, group in gene_groups:
        # Count votes for each direction
        vote_counts = group["inferred_influence"].value_counts()
        majority_influence = vote_counts.idxmax()

        # Get the GenAge annotation (should be same for all entries of this gene)
        genage_annotation = group["genage_annotation"].iloc[0]

        # Determine agreement based on majority vote
        agreement = genage_annotation == majority_influence

        merged_results.append(
            {
                "genage_matched_genes": gene,
                "inferred_influence": majority_influence,
                "genage_annotation": genage_annotation,
                "agreement": agreement,
                "num_interventions": len(group),
            }
        )

    matched_genes_display = pd.DataFrame(merged_results)
    agreement_count = matched_genes_display["agreement"].sum()
    disagreement_count = len(matched_genes_display) - agreement_count
    not_in_genage = len(df) - len(matched_genes)

    print(
        f"\nMerged to {len(matched_genes_display)} unique genes (from {len(matched_genes)} interventions)"
    )
    print(f"  Agreement: {agreement_count}")
    print(f"  Disagreement: {disagreement_count}")
    print(f"Not in GenAge: {not_in_genage}")
else:
    print("\n=== Using UNMERGED mode (all interventions) ===")
    matched_genes_display = matched_genes
    agreement_count = matched_genes["agreement"].sum()
    disagreement_count = len(matched_genes) - agreement_count
    not_in_genage = len(df) - len(matched_genes)

    # Count unique genes for reference
    unique_genes = matched_genes["genage_matched_genes"].nunique()

    print(f"\nAmong matched interventions:")
    print(f"  Agreement: {agreement_count}")
    print(f"  Disagreement: {disagreement_count}")
    print(f"  Unique genes: {unique_genes}")
    print(f"Not in GenAge: {not_in_genage}")

# Print disagreement details
if disagreement_count > 0:
    print("\n" + "=" * 60)
    print("Disagreement Details")
    print("=" * 60)
    disagreement_genes = matched_genes_display[
        matched_genes_display["agreement"] == False
    ]
    for idx, row in disagreement_genes.iterrows():
        gene_name = row["genage_matched_genes"]
        genage_ann = row["genage_annotation"]
        inferred_ann = row["inferred_influence"]
        if USE_MERGED:
            num_interventions = row["num_interventions"]
            print(
                f"{gene_name}: GenAge={genage_ann}, Inferred={inferred_ann} (from {num_interventions} interventions)"
            )
        else:
            trust_score = row["trust_score"]
            intervention = row["intervention"]
            print(
                f"{gene_name} ({intervention}): GenAge={genage_ann}, Inferred={inferred_ann}, Trust Score={trust_score}"
            )
    print("=" * 60)

# Get matched gene names for annotation (keep all instances, not unique)
matched_gene_names = matched_genes_display["genage_matched_genes"].tolist()
unique_genes = matched_genes_display["genage_matched_genes"].nunique()
if USE_MERGED:
    print(f"\nMatched unique genes: {len(matched_gene_names)} total")
else:
    print(
        f"\nMatched gene interventions: {len(matched_gene_names)} total ({unique_genes} unique genes)"
    )

# Calculate statistics
if len(matched_genes_display) > 0:
    agreement_rate = agreement_count / len(matched_genes_display)
    print(f"\nAgreement rate: {agreement_rate:.1%}")

    # Binomial test
    from scipy.stats import binomtest

    result = binomtest(
        agreement_count, len(matched_genes_display), p=0.5, alternative="greater"
    )
    p_value = result.pvalue
    print(f"Binomial test p-value: {p_value:.4e}")
else:
    agreement_rate = 0
    p_value = 1.0
    print("\nNo matched genes found for statistical analysis")

print("=" * 60)

# Prepare data for waffle plot
if USE_MERGED:
    TOTAL_INTERVENTIONS = len(matched_genes_display) + not_in_genage
    MATCHED_GENES = len(matched_genes_display)
else:
    TOTAL_INTERVENTIONS = len(df)
    MATCHED_GENES = len(matched_genes)

AGREEMENT_COUNT = agreement_count
DISAGREEMENT_COUNT = disagreement_count
NOT_IN_GENAGE = not_in_genage

# Nature color palette
AGREE_COLOR = "#3B8EA5"  # Nature blue - agreement
DISAGREE_COLOR = "#E64B35"  # Nature red - disagreement
NOT_IN_GENAGE_COLOR = "#D3D3D3"  # Gray for not in GenAge
BORDER_COLOR = "#FFFFFF"  # White border between squares

# Create figure - increased height for vertical labels
fig = plt.figure(figsize=(24, 8))

# Create main axis for waffle plot
ax = fig.add_subplot(111)

# Waffle plot configuration - adjust based on total number
# Let's use 10 rows x appropriate columns to fit all interventions
n_cols = 100
n_rows = int(np.ceil(TOTAL_INTERVENTIONS / n_cols))

print(
    f"\nWaffle configuration: {n_rows} rows x {n_cols} cols = {n_rows * n_cols} squares"
)

# Calculate square positions
square_size = 0.18
gap = 0.003

# Create waffle data - preserve original order so labels match
# First, add the matched genes in their original order
waffle_data = []
for idx, row in matched_genes_display.iterrows():
    if row["agreement"]:
        waffle_data.append("agree")
    else:
        waffle_data.append("disagree")

# Then add the non-matched ones
waffle_data.extend(["not_in_genage"] * NOT_IN_GENAGE)

# Ensure we have the right number
assert (
    len(waffle_data) == TOTAL_INTERVENTIONS
), f"Waffle data length mismatch: {len(waffle_data)} vs {TOTAL_INTERVENTIONS}"

# Draw waffle squares and collect positions for matched genes
square_idx = 0
matched_positions = []  # Store positions of matched genes for annotation

for row in range(n_rows - 1, -1, -1):  # Start from top row
    for col in range(n_cols):
        if square_idx < len(waffle_data):
            # Determine color based on category
            if waffle_data[square_idx] == "agree":
                color = AGREE_COLOR
                edge_color = BORDER_COLOR
                edge_width = 0.2
                alpha = 1.0
                # Store position for annotation
                x = col * (square_size + gap) + square_size / 2
                y = row * (square_size + gap) + square_size / 2
                matched_positions.append((x, y, "agree"))
            elif waffle_data[square_idx] == "disagree":
                color = DISAGREE_COLOR
                edge_color = BORDER_COLOR
                edge_width = 0.2
                alpha = 1.0
                # Store position for annotation
                x = col * (square_size + gap) + square_size / 2
                y = row * (square_size + gap) + square_size / 2
                matched_positions.append((x, y, "disagree"))
            else:  # not_in_genage
                color = NOT_IN_GENAGE_COLOR
                edge_color = BORDER_COLOR
                edge_width = 0.15
                alpha = 0.6

            # Draw square
            rect = Rectangle(
                (col * (square_size + gap), row * (square_size + gap)),
                square_size,
                square_size,
                facecolor=color,
                edgecolor=edge_color,
                linewidth=edge_width,
                alpha=alpha,
            )
            ax.add_patch(rect)

            square_idx += 1

# Add gene name labels alternating above and below boxes with 45-degree rotation and connecting lines
label_padding = 0.15  # Vertical padding between square and label
line_padding = 0.05  # Small gap between square and line start

for idx, (x, y, category) in enumerate(matched_positions):
    if idx < len(matched_gene_names):
        gene_name = matched_gene_names[idx]
        if pd.notna(gene_name) and gene_name:
            # Set color based on agreement/disagreement
            text_color = AGREE_COLOR if category == "agree" else DISAGREE_COLOR

            # Alternate between above (even idx) and below (odd idx)
            if idx % 2 == 0:  # Even: place above
                square_top_y = y + square_size / 2
                line_start_y = square_top_y + line_padding
                label_y = square_top_y + label_padding
                va_align = "bottom"
                ha_align = "left"
            else:  # Odd: place below
                square_bottom_y = y - square_size / 2
                line_start_y = square_bottom_y - line_padding
                label_y = square_bottom_y - label_padding
                va_align = "top"
                ha_align = "right"

            # Draw connecting line from square to label
            if idx % 2 == 0:
                ax.plot(
                    [x, x],
                    [y + square_size / 2, label_y],
                    color=text_color,
                    linewidth=2,
                    alpha=0.8,
                    zorder=10,
                )
            else:
                ax.plot(
                    [x, x],
                    [y - square_size / 2, label_y],
                    color=text_color,
                    linewidth=2,
                    alpha=0.8,
                    zorder=10,
                )

            # Add text label at 45 degrees - larger font since we have more space
            ax.text(
                x,
                label_y,
                gene_name,
                fontsize=16,
                ha=ha_align,
                va=va_align,
                color=text_color,
                weight="bold",
                rotation=45,
                zorder=11,
            )

# Configure axes - extend y-axis to accommodate labels above
ax.set_xlim(-0.5, n_cols * (square_size + gap))
ax.set_ylim(
    -0.8, n_rows * (square_size + gap) + 1.5
)  # Extended for labels, reduced bottom spacing
ax.set_aspect("equal")
ax.axis("off")

# Add title and p-value
title_y = n_rows * (square_size + gap) + 0.8  # Further reduced spacing for labels
ax.text(
    n_cols * (square_size + gap) / 2,
    title_y,
    f"GenAge Database Validation (n = {TOTAL_INTERVENTIONS})",
    fontsize=24,
    fontweight="bold",
    ha="center",
    color="#000000",
)

# Display exact p-value
if len(matched_genes) > 0:
    p_text = f"Binomial test: p = {p_value:.4e}"
    ax.text(
        n_cols * (square_size + gap) / 2,
        title_y - 0.35,
        p_text,
        fontsize=18,
        fontweight="normal",
        ha="center",
        color="#333333",
    )

# Add legend below with larger font
legend_elements = [
    Line2D(
        [0],
        [0],
        marker="s",
        color="w",
        markerfacecolor=AGREE_COLOR,
        markersize=14,
        label=f"Agreement with GenAge ({AGREEMENT_COUNT})",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        color="w",
        markerfacecolor=DISAGREE_COLOR,
        markersize=14,
        label=f"Disagreement with GenAge ({DISAGREEMENT_COUNT})",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        color="w",
        markerfacecolor=NOT_IN_GENAGE_COLOR,
        alpha=0.6,
        markersize=14,
        label=f"Not in GenAge database ({NOT_IN_GENAGE})",
    ),
]

legend = ax.legend(
    handles=legend_elements,
    loc="upper center",
    bbox_to_anchor=(0.5, 1),  # Further reduced spacing between legend and figure
    ncol=3,
    frameon=False,
    fontsize=16,
    columnspacing=3.0,
)
for text in legend.get_texts():
    text.set_fontweight("normal")

# Apply theme transformations
theme.apply_transforms()

# Fine-tune layout
plt.tight_layout()

# Save figure with high quality for Nature submission
plt.savefig(
    f"figures/geneage_waffle_nature{SUFFIX}.pdf",
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)
plt.savefig(
    f"figures/geneage_waffle_nature{SUFFIX}.png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)
plt.savefig(
    f"figures/geneage_waffle_nature{SUFFIX}.svg",
    format="svg",
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)

plt.show()

# Print detailed summary for manuscript
print("\n" + "=" * 60)
print("GenAge Database Validation Summary (For Manuscript)")
print("=" * 60)
print(f"Total genetic interventions analyzed: {TOTAL_INTERVENTIONS}")
print(
    f"Matched with GenAge database: {MATCHED_GENES} interventions ({MATCHED_GENES/TOTAL_INTERVENTIONS:.1%})"
)
print(f"  Representing {unique_genes} unique genes")
if MATCHED_GENES > 0:
    print(f"Agreement: {AGREEMENT_COUNT}/{MATCHED_GENES} ({agreement_rate:.1%})")
    print(
        f"Disagreement: {DISAGREEMENT_COUNT}/{MATCHED_GENES} ({DISAGREEMENT_COUNT/MATCHED_GENES:.1%})"
    )
    print(f"Binomial test p-value: {p_value:.4e}")
print(f"Not in GenAge: {NOT_IN_GENAGE} ({NOT_IN_GENAGE/TOTAL_INTERVENTIONS:.1%})")
print("=" * 60)

# Print matched gene details
if len(matched_genes) > 0:
    print("\n" + "=" * 60)
    print("Matched Gene Intervention Details")
    print("=" * 60)
    for idx, row in matched_genes.iterrows():
        gene_name = row["genage_matched_genes"]
        genage_ann = row["genage_annotation"]
        inferred_ann = row["inferred_influence"]
        intervention = row["intervention"]
        agree_status = "AGREE" if row["agreement"] else "DISAGREE"
        print(
            f"{gene_name} ({intervention}): GenAge={genage_ann}, Inferred={inferred_ann} [{agree_status}]"
        )
    print("=" * 60)
