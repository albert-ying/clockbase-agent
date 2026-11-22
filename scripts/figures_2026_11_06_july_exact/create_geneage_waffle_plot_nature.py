#!/usr/bin/env python3
"""
Create Nature publication-quality waffle plot showing replicability of genetic interventions
with the GenAge database - high-quality interventions (trust score > 80)
"""

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy import stats

warnings.filterwarnings("ignore")

# Configuration flag: Set to True to merge overlapping entries using majority voting
USE_MERGED = True  # Set to True to merge duplicate genes by majority vote

# Apply consistent theme
from aquarel import load_theme

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
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["mathtext.fontset"] = "custom"
plt.rcParams["mathtext.rm"] = "Arial"
plt.rcParams["mathtext.it"] = "Arial:italic"
plt.rcParams["mathtext.bf"] = "Arial:bold"
plt.rcParams["pdf.fonttype"] = 42

# Load data
df = pd.read_csv(
    "data/agent_result_genage.csv"
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
    print("\n=== Using MERGED mode (majority vote per gene; one entry per unique gene) ===")
    # Group by gene only: same gene across multiple studies counts as one observation
    # (avoids inflating N due to multi-study repeats of the same gene)
    gene_groups = matched_genes.groupby("genage_matched_genes")

    merged_results = []
    for gene, group in gene_groups:
        vote_counts = group["inferred_influence"].value_counts()
        majority_influence = vote_counts.idxmax()
        genage_annotation = group["genage_annotation"].iloc[0]
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

# ---------------------------------------------------------------------------
# Print-size layout (2026-11-06): the panel is drawn at its final size,
# 183 mm x 28 mm (Nature full width), and placed at 100% scale. One axes
# covers the whole figure with data units in millimetres, so positions below
# are physical distances and font sizes are final point sizes.
# ---------------------------------------------------------------------------
FIG_W_MM, FIG_H_MM = 183.0, 28.0
MM_PER_PT = 25.4 / 72.0

GENE_FONT_PT = 6.0  # gene labels, italic (gene symbols)
LEGEND_FONT_PT = 6.5
TITLE_FONT_PT = 7.0
LEADER_LW_PT = 0.6
SQUARE_EDGE_LW_PT = 0.2
SQUARE_GAP_MM = 0.12
EDGE_MARGIN_MM = 0.5  # clearance between ink and the page edge
CORNER_FREE_MM = 4.0  # top-left corner reserved for the panel letter
HEADER_GAP_MM = 0.5  # clearance between header line and the gene labels
# Direction arrows after the gene symbol: Arial's arrow glyphs have a shaft
# of about 0.23 pt at 6 pt, so a stroke in the label colour thickens them to
# about 0.5 pt; the ink gap between symbol and arrow is 0.15 em
ARROW_STROKE_PT = 0.25
ARROW_GAP_MM = 0.15 * GENE_FONT_PT * MM_PER_PT


def format_p_value(p):
    """P value to two significant figures; mathtext for small values."""
    if p <= 0:
        return "0"
    exponent = int(np.floor(np.log10(p)))
    if exponent >= -3:
        return f"{p:.{max(0, 1 - exponent)}f}"
    mantissa = p / 10**exponent
    if round(mantissa, 1) >= 10:
        mantissa /= 10
        exponent += 1
    return rf"${mantissa:.1f}\times10^{{{exponent}}}$"


fig = plt.figure(figsize=(FIG_W_MM / 25.4, FIG_H_MM / 25.4))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, FIG_W_MM)
ax.set_ylim(0, FIG_H_MM)
ax.set_xticks([])
ax.set_yticks([])
ax.axis("off")

# Waffle plot configuration: one row holding every square
n_cols = TOTAL_INTERVENTIONS
n_rows = 1
print(
    f"\nWaffle configuration: {n_rows} rows x {n_cols} cols = {n_rows * n_cols} squares"
)

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

# Gene labels: one per square, alternating above and below. A gene
# manipulated in both directions has two adjacent squares ("<gene> decrease",
# "<gene> increase"); their labels carry a down or up arrow.
DIRECTION_MARKS = {"decrease": "↓", "increase": "↑"}
labels = []  # (square index, gene symbol, arrow or None, colour, above?)
matched_idx = 0
for square_idx, category in enumerate(waffle_data):
    if category == "not_in_genage":
        continue
    gene_name = matched_gene_names[matched_idx]
    matched_idx += 1
    if pd.notna(gene_name) and gene_name:
        parts = gene_name.split()
        symbol = parts[0].capitalize()
        mark = DIRECTION_MARKS[parts[1]] if len(parts) > 1 else None
        text_color = AGREE_COLOR if category == "agree" else DISAGREE_COLOR
        labels.append((square_idx, symbol, mark, text_color, len(labels) % 2 == 0))
n_pairs = sum(mark is not None for _, _, mark, _, _ in labels) // 2
print(f"Labels: {len(labels)} squares, {len(labels) - n_pairs} genes, "
      f"{n_pairs} genes manipulated in both directions")

# Squares (positions are set by place_elements below)
squares = []
for category in waffle_data:
    if category == "agree":
        color, edge_width, alpha = AGREE_COLOR, SQUARE_EDGE_LW_PT, 1.0
    elif category == "disagree":
        color, edge_width, alpha = DISAGREE_COLOR, SQUARE_EDGE_LW_PT, 1.0
    else:  # not_in_genage
        color, edge_width, alpha = NOT_IN_GENAGE_COLOR, SQUARE_EDGE_LW_PT * 0.75, 0.6
    rect = Rectangle(
        (0, 0), 1, 1,
        facecolor=color,
        edgecolor=BORDER_COLOR,
        linewidth=edge_width,
        alpha=alpha,
    )
    ax.add_patch(rect)
    squares.append(rect)

# Gene label texts and leader ticks; labels run at 45 degrees from the
# tick end (up-right above the strip, down-left below it). The direction
# arrow is a separate upright (non-italic) text set after the symbol.
label_artists = []
for square_idx, gene_name, mark, text_color, above in labels:
    (leader,) = ax.plot(
        [0, 0], [0, 0],
        color=text_color,
        linewidth=LEADER_LW_PT,
        alpha=0.8,
        solid_capstyle="butt",
        zorder=10,
    )
    label_style = dict(
        fontsize=GENE_FONT_PT,
        ha="left" if above else "right",
        va="center",
        rotation=45,
        rotation_mode="anchor",
        color=text_color,
        zorder=11,
    )
    parts = [ax.text(0, 0, gene_name, fontstyle="italic", **label_style)]
    if mark is not None:
        parts.append(
            ax.text(
                0, 0,
                mark,
                path_effects=[
                    patheffects.withStroke(linewidth=ARROW_STROKE_PT, foreground=text_color)
                ],
                **label_style,
            )
        )
    label_artists.append((square_idx, above, leader, parts))

# Header line: title, legend and test result
header_title = ax.text(
    0, 0,
    f"GenAge gene manipulations (n = {TOTAL_INTERVENTIONS})",
    fontsize=TITLE_FONT_PT,
    fontweight="bold",
    ha="left",
    va="baseline",
    color="#000000",
)
legend_items = [
    (AGREE_COLOR, 1.0, f"Agreement with GenAge ({AGREEMENT_COUNT})", AGREEMENT_COUNT),
    (DISAGREE_COLOR, 1.0, f"Disagreement with GenAge ({DISAGREEMENT_COUNT})", DISAGREEMENT_COUNT),
    (NOT_IN_GENAGE_COLOR, 0.6, f"Not in GenAge database ({NOT_IN_GENAGE})", NOT_IN_GENAGE),
]
# Legend entries with zero members are not drawn
legend_artists = []
for color, alpha, label, count in legend_items:
    if count == 0:
        continue
    swatch = Rectangle((0, 0), 1, 1, facecolor=color, edgecolor="none", alpha=alpha)
    ax.add_patch(swatch)
    text = ax.text(0, 0, label, fontsize=LEGEND_FONT_PT, ha="left", va="baseline",
                   color="#000000")
    legend_artists.append((swatch, text))
# Test result as plain-text pieces (italic P set separately): a mathtext
# string sits about 0.3 mm above the plain-text baseline of the title and
# legend, so mathtext is used only for a value below 0.001.
p_artists = []
if len(matched_genes) > 0:
    for piece, style in (
        ("One-sided binomial", "normal"),
        ("P", "italic"),
        ("=", "normal"),
        (format_p_value(p_value), "normal"),
    ):
        p_artists.append(
            ax.text(
                0, 0,
                piece,
                fontsize=LEGEND_FONT_PT,
                fontstyle=style,
                ha="right",
                va="baseline",
                color="#333333",
            )
        )
SPACE_MM = 0.278 * LEGEND_FONT_PT * MM_PER_PT  # Arial word space (569/2048 em)

renderer = fig.canvas.get_renderer()
to_data = ax.transData.inverted()


def extent_mm(artist):
    bb = artist.get_window_extent(renderer)
    (x0, y0), (x1, y1) = to_data.transform([[bb.x0, bb.y0], [bb.x1, bb.y1]])
    return x0, y0, x1, y1


def baseline_width_mm(text):
    """Width of a rotated text along its own baseline."""
    text.set_rotation(0)
    x0, _, x1, _ = extent_mm(text)
    text.set_rotation(45)
    return x1 - x0


BASELINE_DIR = np.array([np.cos(np.pi / 4), np.sin(np.pi / 4)])
label_widths = [[baseline_width_mm(t) for t in parts] for _, _, _, parts in label_artists]


def place_label(parts, widths, above, x, y):
    """Set a label (symbol, then arrow if any) whose tick-side end is at (x, y)."""
    anchor = np.array([x, y])
    if len(parts) == 1:
        parts[0].set_position((x, y))
    elif above:  # reads up-right from the tick end; arrow after the symbol
        parts[0].set_position((x, y))
        parts[1].set_position(tuple(anchor + BASELINE_DIR * (widths[0] + ARROW_GAP_MM)))
    else:  # reads down-left to the tick end; arrow at the tick end
        parts[1].set_position((x, y))
        parts[0].set_position(tuple(anchor - BASELINE_DIR * (widths[1] + ARROW_GAP_MM)))


def union_extent(artists):
    es = [extent_mm(a) for a in artists]
    return (min(e[0] for e in es), min(e[1] for e in es),
            max(e[2] for e in es), max(e[3] for e in es))


# Rotated label extents relative to their anchor (independent of position)
label_offsets = []
for (_, above, _, parts), widths in zip(label_artists, label_widths):
    place_label(parts, widths, above, 0, 0)
    label_offsets.append(union_extent(parts))
below_depth = max(-off[1] for (_, above, _, _), off in zip(label_artists, label_offsets) if not above)
above_height = max(off[3] for (_, above, _, _), off in zip(label_artists, label_offsets) if above)

header_title.set_position((0, 0))
title_ext = extent_mm(header_title)
header_ascent = title_ext[3]
header_descent = -title_ext[1]
for _, text in legend_artists:
    text.set_position((0, 0))
    e = extent_mm(text)
    header_ascent = max(header_ascent, e[3])
    header_descent = max(header_descent, -e[1])
for p_artist in p_artists:
    p_artist.set_position((0, 0))
    e = extent_mm(p_artist)
    header_ascent = max(header_ascent, e[3])
    header_descent = max(header_descent, -e[1])


def place_elements(x_start, pitch, leader_len):
    square_size = pitch - SQUARE_GAP_MM
    y_below_anchor = EDGE_MARGIN_MM + below_depth
    y_square = y_below_anchor + leader_len
    y_above_anchor = y_square + square_size + leader_len
    for i, rect in enumerate(squares):
        rect.set_xy((x_start + i * pitch + SQUARE_GAP_MM / 2, y_square))
        rect.set_width(square_size)
        rect.set_height(square_size)
    for (square_idx, above, leader, parts), widths in zip(label_artists, label_widths):
        x = x_start + square_idx * pitch + pitch / 2
        if above:
            leader.set_data([x, x], [y_square + square_size, y_above_anchor])
            place_label(parts, widths, above, x, y_above_anchor)
        else:
            leader.set_data([x, x], [y_square, y_below_anchor])
            place_label(parts, widths, above, x, y_below_anchor)
    return y_above_anchor + above_height


# Solve the horizontal layout so the strip and its labels span the page width
x_start, pitch = 2.0, (FIG_W_MM - 8.0) / n_cols
for _ in range(8):
    place_elements(x_start, pitch, 1.0)
    xs = [extent_mm(r) for r in squares] + [union_extent(p) for _, _, _, p in label_artists]
    ink_left = min(e[0] for e in xs)
    ink_right = max(e[2] for e in xs)
    target_span = FIG_W_MM - 2 * EDGE_MARGIN_MM
    strip_left, strip_right = x_start, x_start + n_cols * pitch
    lead = strip_left - ink_left
    tail = ink_right - strip_right
    pitch = (target_span - lead - tail) / n_cols
    x_start = EDGE_MARGIN_MM + lead

# Vertical layout: leader ticks take the space left after the labels,
# squares and header line
square_size = pitch - SQUARE_GAP_MM
# Header baseline level with the baseline of the 10-pt panel letter that
# assemble_figure4_vector.py sets 10 pt x 0.78 below the panel top
header_baseline = FIG_H_MM - max(10 * 0.78 * MM_PER_PT, EDGE_MARGIN_MM + header_ascent)
header_bottom = header_baseline - header_descent
fixed_height = EDGE_MARGIN_MM + below_depth + square_size + above_height + HEADER_GAP_MM
leader_len = (header_bottom - fixed_height) / 2
labels_top = place_elements(x_start, pitch, leader_len)
print(
    f"Layout: pitch {pitch:.2f} mm, square {square_size:.2f} mm, "
    f"leader {leader_len:.2f} mm, labels top {labels_top:.2f} mm, "
    f"header baseline {header_baseline:.2f} mm"
)

# Header: title at left (clear of the panel-letter corner), legend centred,
# test result right-aligned to the end of the strip
swatch_mm = LEGEND_FONT_PT * MM_PER_PT * 0.75
swatch_gap_mm = 1.0
entry_gap_mm = 4.0
header_title.set_position((CORNER_FREE_MM + 0.6, header_baseline))
legend_width = 0.0
entry_widths = []
for swatch, text in legend_artists:
    text.set_position((0, header_baseline))
    e = extent_mm(text)
    entry_widths.append(e[2] - e[0])
legend_width = sum(swatch_mm + swatch_gap_mm + w for w in entry_widths) + entry_gap_mm * (
    len(entry_widths) - 1
)
x_cursor = FIG_W_MM / 2 - legend_width / 2
for (swatch, text), w in zip(legend_artists, entry_widths):
    swatch.set_xy((x_cursor, header_baseline - 0.05))
    swatch.set_width(swatch_mm)
    swatch.set_height(swatch_mm)
    text.set_position((x_cursor + swatch_mm + swatch_gap_mm, header_baseline))
    x_cursor += swatch_mm + swatch_gap_mm + w + entry_gap_mm
x_right = min(x_start + n_cols * pitch - SQUARE_GAP_MM / 2, FIG_W_MM - EDGE_MARGIN_MM)
for p_artist in reversed(p_artists):
    p_artist.set_position((x_right, header_baseline))
    x_right = extent_mm(p_artist)[0] - SPACE_MM

# Apply theme transformations
theme.apply_transforms()

# Save at the exact print size (no bbox_inches='tight')
plt.savefig(
    "figures/geneage_waffle_nature.pdf",
    dpi=600,
    facecolor="white",
    edgecolor="none",
)
plt.savefig(
    "figures/geneage_waffle_nature.png",
    dpi=600,
    facecolor="white",
    edgecolor="none",
)
plt.savefig(
    "figures/geneage_waffle_nature.svg",
    format="svg",
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
