#!/usr/bin/env python3
"""
Create Nature publication-quality waffle plot showing replicability of drug interventions
with the DrugAge database - high-quality interventions (trust score > 80)
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
USE_MERGED = True  # Set to True to merge duplicate drugs by majority vote

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
df = pd.read_csv(
    f"data/agent_result_drugage_{METHOD_NAME}.csv"
)

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
print("DrugAge Database Validation Analysis")
print("=" * 60)

# Count matched drugs (drugage_annotated = True)
matched_drugs = df[df["drugage_annotated"] == True].copy()
print(f"\nTotal interventions in dataset: {len(df)}")
print(f"Matched with DrugAge: {len(matched_drugs)}")

# Among matched drugs, check agreement
# Agreement: both drugage_annotation and inferred_influence are the same
matched_drugs["agreement"] = (
    matched_drugs["drugage_annotation"] == matched_drugs["inferred_influence"]
)

# Handle merging if USE_MERGED is True
if USE_MERGED:
    print("\n=== Using MERGED mode (majority voting) ===")
    # Group by drug and use majority voting for inferred_influence
    drug_groups = matched_drugs.groupby("drugage_matched_drugs")

    merged_results = []
    for drug, group in drug_groups:
        # Count votes for each direction
        vote_counts = group["inferred_influence"].value_counts()
        majority_influence = vote_counts.idxmax()

        # Get the DrugAge annotation (should be same for all entries of this drug)
        drugage_annotation = group["drugage_annotation"].iloc[0]

        # Determine agreement based on majority vote
        agreement = drugage_annotation == majority_influence

        merged_results.append(
            {
                "drugage_matched_drugs": drug,
                "inferred_influence": majority_influence,
                "drugage_annotation": drugage_annotation,
                "agreement": agreement,
                "num_interventions": len(group),
            }
        )

    matched_drugs_display = pd.DataFrame(merged_results)
    agreement_count = matched_drugs_display["agreement"].sum()
    disagreement_count = len(matched_drugs_display) - agreement_count
    not_in_drugage = len(df) - len(matched_drugs)

    print(
        f"\nMerged to {len(matched_drugs_display)} unique drugs (from {len(matched_drugs)} interventions)"
    )
    print(f"  Agreement: {agreement_count}")
    print(f"  Disagreement: {disagreement_count}")
    print(f"Not in DrugAge: {not_in_drugage}")
else:
    print("\n=== Using UNMERGED mode (all interventions) ===")
    matched_drugs_display = matched_drugs
    agreement_count = matched_drugs["agreement"].sum()
    disagreement_count = len(matched_drugs) - agreement_count
    not_in_drugage = len(df) - len(matched_drugs)

    print(f"\nAmong matched drugs:")
    print(f"  Agreement: {agreement_count}")
    print(f"  Disagreement: {disagreement_count}")
    print(f"Not in DrugAge: {not_in_drugage}")

# Print disagreement details
if disagreement_count > 0:
    print("\n" + "=" * 60)
    print("Disagreement Details")
    print("=" * 60)
    disagreement_drugs = matched_drugs_display[
        matched_drugs_display["agreement"] == False
    ]
    for idx, row in disagreement_drugs.iterrows():
        drug_name = row["drugage_matched_drugs"]
        drugage_ann = row["drugage_annotation"]
        inferred_ann = row["inferred_influence"]
        if USE_MERGED:
            num_interventions = row["num_interventions"]
            print(
                f"{drug_name}: DrugAge={drugage_ann}, Inferred={inferred_ann} (from {num_interventions} interventions)"
            )
        else:
            trust_score = row["trust_score"]
            drugage_annotated = row["drugage_annotated"]
            print(
                f"{drug_name}: DrugAge={drugage_ann}, Inferred={inferred_ann}, Trust Score={trust_score}, DrugAge Annotated={drugage_annotated}"
            )
    print("=" * 60)

# Get matched drug names for annotation (keep all instances, not unique)
matched_drug_names = matched_drugs_display["drugage_matched_drugs"].tolist()
if USE_MERGED:
    print(f"\nMatched unique drugs: {len(matched_drug_names)} total")
else:
    print(f"\nMatched drugs: {matched_drug_names}")

# Calculate statistics
if len(matched_drugs_display) > 0:
    agreement_rate = agreement_count / len(matched_drugs_display)
    print(f"\nAgreement rate: {agreement_rate:.1%}")

    # Binomial test
    from scipy.stats import binomtest

    result = binomtest(
        agreement_count, len(matched_drugs_display), p=0.5, alternative="greater"
    )
    p_value = result.pvalue
    print(f"Binomial test p-value: {p_value:.4e}")
else:
    agreement_rate = 0
    p_value = 1.0
    print("\nNo matched drugs found for statistical analysis")

print("=" * 60)

# Prepare data for waffle plot
if USE_MERGED:
    TOTAL_INTERVENTIONS = len(matched_drugs_display) + not_in_drugage
    MATCHED_DRUGS = len(matched_drugs_display)
else:
    TOTAL_INTERVENTIONS = len(df)
    MATCHED_DRUGS = len(matched_drugs)

AGREEMENT_COUNT = agreement_count
DISAGREEMENT_COUNT = disagreement_count
NOT_IN_DRUGAGE = not_in_drugage

# Nature color palette
AGREE_COLOR = "#3B8EA5"  # Nature blue - agreement
DISAGREE_COLOR = "#E64B35"  # Nature red - disagreement
NOT_IN_DRUGAGE_COLOR = "#D3D3D3"  # Gray for not in DrugAge
BORDER_COLOR = "#FFFFFF"  # White border between squares

# Create figure - increased height for vertical labels
fig = plt.figure(figsize=(24, 6))

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
# First, add the matched drugs in their original order
waffle_data = []
for idx, row in matched_drugs_display.iterrows():
    if row["agreement"]:
        waffle_data.append("agree")
    else:
        waffle_data.append("disagree")

# Then add the non-matched ones
waffle_data.extend(["not_in_drugage"] * NOT_IN_DRUGAGE)

# Ensure we have the right number
assert (
    len(waffle_data) == TOTAL_INTERVENTIONS
), f"Waffle data length mismatch: {len(waffle_data)} vs {TOTAL_INTERVENTIONS}"

# Draw waffle squares and collect positions for matched drugs
square_idx = 0
matched_positions = []  # Store positions of matched drugs for annotation

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
            else:  # not_in_drugage
                color = NOT_IN_DRUGAGE_COLOR
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

# Add drug name labels alternating above and below boxes with 45-degree rotation and connecting lines
label_padding = 0.15  # Vertical padding between square and label
line_padding = 0.05  # Small gap between square and line start

for idx, (x, y, category) in enumerate(matched_positions):
    if idx < len(matched_drug_names):
        drug_name = matched_drug_names[idx]
        if pd.notna(drug_name) and drug_name:
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
                drug_name,
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
    f"DrugAge Database Validation (n = {TOTAL_INTERVENTIONS}) ",
    fontsize=24,
    fontweight="bold",
    ha="center",
    color="#000000",
)

# Display exact p-value
if len(matched_drugs) > 0:
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
        label=f"Agreement with DrugAge ({AGREEMENT_COUNT})",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        color="w",
        markerfacecolor=DISAGREE_COLOR,
        markersize=14,
        label=f"Disagreement with DrugAge ({DISAGREEMENT_COUNT})",
    ),
    Line2D(
        [0],
        [0],
        marker="s",
        color="w",
        markerfacecolor=NOT_IN_DRUGAGE_COLOR,
        alpha=0.6,
        markersize=14,
        label=f"Not in DrugAge database ({NOT_IN_DRUGAGE})",
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
    f"figures/drugage_waffle_nature{SUFFIX}.pdf",
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)
plt.savefig(
    f"figures/drugage_waffle_nature{SUFFIX}.png",
    dpi=600,
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)
plt.savefig(
    f"figures/drugage_waffle_nature{SUFFIX}.svg",
    format="svg",
    bbox_inches="tight",
    facecolor="white",
    edgecolor="none",
)

plt.show()

# Print detailed summary for manuscript
print("\n" + "=" * 60)
print("DrugAge Database Validation Summary (For Manuscript)")
print("=" * 60)
print(f"Total drug interventions analyzed: {TOTAL_INTERVENTIONS}")
print(
    f"Matched with DrugAge database: {MATCHED_DRUGS} ({MATCHED_DRUGS/TOTAL_INTERVENTIONS:.1%})"
)
if MATCHED_DRUGS > 0:
    print(f"Agreement: {AGREEMENT_COUNT}/{MATCHED_DRUGS} ({agreement_rate:.1%})")
    print(
        f"Disagreement: {DISAGREEMENT_COUNT}/{MATCHED_DRUGS} ({DISAGREEMENT_COUNT/MATCHED_DRUGS:.1%})"
    )
    print(f"Binomial test p-value: {p_value:.4e}")
print(f"Not in DrugAge: {NOT_IN_DRUGAGE} ({NOT_IN_DRUGAGE/TOTAL_INTERVENTIONS:.1%})")
print("=" * 60)

# Print matched drug details
if len(matched_drugs) > 0:
    print("\n" + "=" * 60)
    print("Matched Drug Details")
    print("=" * 60)
    for idx, row in matched_drugs.iterrows():
        drug_name = row["drugage_matched_drugs"]
        drugage_ann = row["drugage_annotation"]
        inferred_ann = row["inferred_influence"]
        agree_status = "AGREE" if row["agreement"] else "DISAGREE"
        print(
            f"{drug_name}: DrugAge={drugage_ann}, Inferred={inferred_ann} [{agree_status}]"
        )
    print("=" * 60)
