import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

# Set Nature paper style
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 8
plt.rcParams["axes.linewidth"] = 0

# Define color palette
colors = {
    "plausibility": "#FFE6E6",  # Light red
    "model": "#E3F2FD",  # Light blue
    "aging": "#E8F5E9",  # Light green
    "rigor": "#F3E5F5",  # Light purple
    "translational": "#FFF3E0",  # Light orange
    "novelty": "#FCE4EC",  # Light pink
    "context": "#F5F5F5",  # Light gray
}

colors_dark = {
    "plausibility": "#E74C3C",
    "model": "#2196F3",
    "aging": "#4CAF50",
    "rigor": "#9C27B0",
    "translational": "#FF9800",
    "novelty": "#E91E63",
    "context": "#757575",
}

# Create figure with better proportions
fig = plt.figure(figsize=(8.5, 11))
ax = fig.add_subplot(111)
ax.set_xlim(0, 8.5)
ax.set_ylim(0, 11)
ax.axis("off")

y_pos = 10.5  # Start from top


def add_title(ax, y, text, fontsize=11, bold=True):
    """Add a centered title"""
    weight = "bold" if bold else "normal"
    ax.text(
        4.25,
        y,
        text,
        ha="center",
        va="top",
        fontsize=fontsize,
        fontweight=weight,
        color="#2C3E50",
    )
    return y - 0.3


def add_section_header(ax, y, section_num, title, points, color, color_dark):
    """Add a colored section header"""
    # Background rectangle
    rect = Rectangle(
        (0.3, y - 0.35), 7.9, 0.35, facecolor=color, edgecolor=color_dark, linewidth=2
    )
    ax.add_patch(rect)

    # Section title
    ax.text(
        0.5,
        y - 0.18,
        f"Section {section_num}: {title}",
        ha="left",
        va="center",
        fontsize=9,
        fontweight="bold",
        color=color_dark,
    )

    # Points
    ax.text(
        8.0,
        y - 0.18,
        points,
        ha="right",
        va="center",
        fontsize=9,
        fontweight="bold",
        color=color_dark,
    )

    return y - 0.45


def add_question(ax, y, question_text, color, max_width=7.4):
    """Add a question with colored background"""
    # Background rectangle
    rect = Rectangle((0.5, y - 0.25), 7.5, 0.25, facecolor=color, edgecolor="none")
    ax.add_patch(rect)

    # Question text
    ax.text(
        0.6,
        y - 0.125,
        question_text,
        ha="left",
        va="center",
        fontsize=8,
        fontweight="bold",
        color="#2C3E50",
        wrap=True,
    )

    return y - 0.3


def add_option(ax, y, text, points, color_dark):
    """Add an option with points"""
    ax.text(0.8, y, f"• {text}", ha="left", va="top", fontsize=7.5, color="#424242")
    ax.text(
        7.8,
        y,
        f"→ {points} pts",
        ha="right",
        va="top",
        fontsize=7.5,
        fontweight="bold",
        color=color_dark,
    )
    return y - 0.2


def add_text(ax, y, text, indent=0.5, fontsize=7.5, italic=False):
    """Add regular text"""
    style = "italic" if italic else "normal"
    ax.text(
        indent,
        y,
        text,
        ha="left",
        va="top",
        fontsize=fontsize,
        color="#424242",
        style=style,
    )
    return y - 0.18


# Main title
y_pos = add_title(
    ax, y_pos, "Aging Intervention Evidence Assessment Questionnaire", fontsize=12
)
y_pos = add_title(
    ax, y_pos, "For Systematic Evaluation of Aging Research", fontsize=9, bold=False
)
y_pos -= 0.15

# Section 0: Plausibility Filter
y_pos = add_section_header(
    ax,
    y_pos,
    "0",
    "Plausibility Filter (CRITICAL)",
    "±50 pts",
    colors["plausibility"],
    colors_dark["plausibility"],
)
y_pos = add_question(
    ax, y_pos, "Q0. Biological Plausibility Assessment", colors["plausibility"]
)
y_pos = add_option(
    ax,
    y_pos,
    "Major concerns (errors, implausible)",
    "-50",
    colors_dark["plausibility"],
)
y_pos = add_option(
    ax, y_pos, "Highly plausible (strong mechanism)", "+30", colors_dark["plausibility"]
)
y_pos -= 0.1

# Section 1: Experimental Model
y_pos = add_section_header(
    ax,
    y_pos,
    "1",
    "Experimental Model Quality",
    "Max 25 pts",
    colors["model"],
    colors_dark["model"],
)
y_pos = add_question(
    ax, y_pos, "Q1. Highest-quality experimental model used", colors["model"]
)
y_pos = add_option(
    ax, y_pos, "Immortalized/cancer cell lines", "0", colors_dark["model"]
)
y_pos = add_option(ax, y_pos, "Primary cell cultures", "15", colors_dark["model"])
y_pos = add_option(
    ax, y_pos, "Organoids/3D culture systems", "18", colors_dark["model"]
)
y_pos = add_option(ax, y_pos, "Ex vivo tissue models", "20", colors_dark["model"])
y_pos = add_option(
    ax, y_pos, "In vivo studies (mice, rats, humans)", "25", colors_dark["model"]
)
y_pos -= 0.1

# Section 2: Aging Biology Relevance
y_pos = add_section_header(
    ax,
    y_pos,
    "2",
    "Aging Biology Relevance",
    "Max 25 pts",
    colors["aging"],
    colors_dark["aging"],
)
y_pos = add_question(
    ax, y_pos, "Q2.1. Targets/affects known aging pathways? (15 pts)", colors["aging"]
)
y_pos = add_option(ax, y_pos, "No mention of aging pathways", "0", colors_dark["aging"])
y_pos = add_option(
    ax, y_pos, "Indirect effect on aging pathways", "5", colors_dark["aging"]
)
y_pos = add_option(
    ax, y_pos, "Direct modulation of one pathway", "10", colors_dark["aging"]
)
y_pos = add_option(
    ax, y_pos, "Affects multiple aging pathways", "15", colors_dark["aging"]
)

y_pos = add_question(
    ax, y_pos, "Q2.2. Links to age-related diseases? (10 pts)", colors["aging"]
)
y_pos = add_option(
    ax, y_pos, "No connection to age-related diseases", "0", colors_dark["aging"]
)
y_pos = add_option(ax, y_pos, "Theoretical discussion only", "5", colors_dark["aging"])
y_pos = add_option(
    ax, y_pos, "Experimental evidence in disease models", "8", colors_dark["aging"]
)
y_pos = add_option(
    ax, y_pos, "Multiple disease impacts with evidence", "10", colors_dark["aging"]
)
y_pos -= 0.1

# Section 3: Experimental Rigor
y_pos = add_section_header(
    ax,
    y_pos,
    "3",
    "Experimental Rigor & Quality",
    "Max 10 pts",
    colors["rigor"],
    colors_dark["rigor"],
)
y_pos = add_question(
    ax, y_pos, "Q3.1. Experimental design rigor (10 pts)", colors["rigor"]
)
y_pos = add_option(
    ax, y_pos, "Poor design with major confounds", "0", colors_dark["rigor"]
)
y_pos = add_option(
    ax, y_pos, "Rigorous design with proper controls", "10", colors_dark["rigor"]
)
y_pos -= 0.1

# Section 4: Translational Potential
y_pos = add_section_header(
    ax,
    y_pos,
    "4",
    "Translational Potential",
    "Max 5 pts",
    colors["translational"],
    colors_dark["translational"],
)
y_pos = add_question(
    ax,
    y_pos,
    "Q4.1. Human relevance and safety profile (5 pts)",
    colors["translational"],
)
y_pos = add_option(
    ax, y_pos, "No human data available", "0", colors_dark["translational"]
)
y_pos = add_option(
    ax,
    y_pos,
    "Limited human data or safety concerns",
    "2",
    colors_dark["translational"],
)
y_pos = add_option(
    ax,
    y_pos,
    "Some human use with acceptable safety",
    "5",
    colors_dark["translational"],
)
y_pos -= 0.1

# Section 5: Novelty and Innovation
y_pos = add_section_header(
    ax,
    y_pos,
    "5",
    "Novelty and Innovation",
    "Max 35 pts",
    colors["novelty"],
    colors_dark["novelty"],
)
y_pos = add_question(
    ax, y_pos, "Q5.1. Literature saturation assessment (15 pts)", colors["novelty"]
)
y_pos = add_option(
    ax, y_pos, "Well-established (50+ aging studies)", "0", colors_dark["novelty"]
)
y_pos = add_option(
    ax, y_pos, "Moderately studied (10-50 publications)", "3", colors_dark["novelty"]
)
y_pos = add_option(
    ax, y_pos, "Emerging intervention (3-10 studies)", "10", colors_dark["novelty"]
)
y_pos = add_option(
    ax, y_pos, "Novel/repurposed (<3 aging studies)", "15", colors_dark["novelty"]
)

y_pos = add_question(
    ax, y_pos, "Q5.2. Mechanistic novelty and insight (20 pts)", colors["novelty"]
)
y_pos = add_option(
    ax, y_pos, "Well-characterized pathways only", "5", colors_dark["novelty"]
)
y_pos = add_option(
    ax, y_pos, "New combination of known pathways", "10", colors_dark["novelty"]
)
y_pos = add_option(
    ax, y_pos, "Novel aging mechanism discovered", "15", colors_dark["novelty"]
)
y_pos = add_option(
    ax, y_pos, "Paradigm-shifting discovery", "20", colors_dark["novelty"]
)
y_pos -= 0.1

# Additional Context (not scored)
y_pos = add_section_header(
    ax,
    y_pos,
    "6",
    "Additional Context (Not Scored)",
    "",
    colors["context"],
    colors_dark["context"],
)
y_pos = add_question(
    ax, y_pos, "Q6. Intervention name (concise label)", colors["context"]
)
y_pos -= 0.05
y_pos = add_question(
    ax, y_pos, "Q7. Experimental model/species/tissue", colors["context"]
)
y_pos -= 0.15

# Final scoring box
rect = Rectangle(
    (0.3, y_pos - 0.5), 7.9, 0.5, facecolor="#34495E", edgecolor="#2C3E50", linewidth=2
)
ax.add_patch(rect)
ax.text(
    4.25,
    y_pos - 0.15,
    "FINAL SCORE = Base Score (0-100) + Plausibility Modifier",
    ha="center",
    va="center",
    fontsize=9,
    fontweight="bold",
    color="white",
)
ax.text(
    4.25,
    y_pos - 0.35,
    "Range: 0-130 points",
    ha="center",
    va="center",
    fontsize=8,
    color="white",
)

plt.tight_layout()

# Create figure folder if it doesn't exist
import os

os.makedirs("figure", exist_ok=True)

# Save figure
plt.savefig(
    "figure/questionnaire_framework.png",
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)
plt.savefig(
    "figure/questionnaire_framework.pdf", bbox_inches="tight", facecolor="white"
)

print(
    "✓ Figure saved as 'figure/questionnaire_framework.png' (300 DPI) and 'figure/questionnaire_framework.pdf'"
)
print("✓ Questionnaire formatted with color-coded sections (no legend)")

plt.show()
