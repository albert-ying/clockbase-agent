#!/usr/bin/env python3
"""
Publication-quality visualizations for score distributions.
Creates three key figures: overall score, subscore distributions, and correlation network.
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
import networkx as nx
from scipy.stats import gaussian_kde

# Try to import aquarel for scientific theme
try:
    from aquarel import load_theme

    USE_AQUAREL = True
except ImportError:
    USE_AQUAREL = False
    print("Note: aquarel not found, using standard matplotlib theme")


def setup_theme():
    """Set up publication-quality theme."""
    global USE_AQUAREL
    if USE_AQUAREL:
        try:
            theme = (
                load_theme("scientific")
                .set_grid(draw=False)
                .set_font(size=14)
                .set_ticks(direction="out")
                .set_axis_labels(pad=10)
            )
            theme.apply()
            return theme
        except:
            USE_AQUAREL = False

    # Fallback to standard scientific settings
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 12,
            "axes.labelsize": 14,
            "axes.titlesize": 16,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 11,
            "figure.titlesize": 18,
            "axes.grid": False,
            "axes.linewidth": 1.5,
            "axes.edgecolor": "#333333",
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "figure.dpi": 100,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.1,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )
    return None


# Define beautiful color palette
COLORS = {
    "q0_plausibility": "#FF6B6B",  # Coral red
    "q1_experimental_model": "#4ECDC4",  # Turquoise
    "q2_1_aging_pathways": "#45B7D1",  # Sky blue
    "q2_2_age_diseases": "#96CEB4",  # Sage green
    "q3_1_experimental_rigor": "#FFEAA7",  # Pale yellow
    "q4_1_translational_potential": "#DDA0DD",  # Plum
    "q5_1_literature_saturation": "#98D8C8",  # Mint
    "q5_2_mechanistic_novelty": "#F7B731",  # Golden
    "final_score": "#2C3E50",  # Dark blue-gray
}

# Define better names for publication
NAMES = {
    "q0_plausibility": "Biological Plausibility",
    "q1_experimental_model": "Model Relevance",
    "q2_1_aging_pathways": "Aging Pathways",
    "q2_2_age_diseases": "Age-Related Diseases",
    "q3_1_experimental_rigor": "Experimental Rigor",
    "q4_1_translational_potential": "Translational Potential",
    "q5_1_literature_saturation": "Literature Coverage",
    "q5_2_mechanistic_novelty": "Mechanistic Novelty",
    "final_score": "Composite Score",
}


def create_overall_score_distribution(df):
    """Create a beautiful histogram for the overall final score."""

    fig, ax = plt.subplots(figsize=(10, 6))

    data = df["final_score"].dropna()

    # Create histogram
    n, bins, patches = ax.hist(
        data,
        bins=40,
        color=COLORS["final_score"],
        alpha=0.7,
        edgecolor="white",
        linewidth=1.5,
    )

    # Color gradient for bars
    cm = plt.cm.viridis
    norm = plt.Normalize(vmin=bins.min(), vmax=bins.max())
    for i, patch in enumerate(patches):
        patch.set_facecolor(cm(norm(bins[i])))

    # Add KDE overlay
    kde = gaussian_kde(data, bw_method="scott")
    x_range = np.linspace(data.min(), data.max(), 300)
    kde_values = kde(x_range)
    kde_scaled = kde_values * len(data) * (bins[1] - bins[0])

    ax2 = ax.twinx()
    ax2.fill_between(x_range, 0, kde_scaled, color="#FF6B6B", alpha=0.3)
    ax2.plot(x_range, kde_scaled, color="#FF6B6B", linewidth=2.5, alpha=0.8)
    ax2.set_ylabel("Density", fontsize=14, fontweight="bold")
    ax2.grid(False)
    ax2.spines["top"].set_visible(False)

    # Statistics
    mean_val = data.mean()
    median_val = data.median()
    std_val = data.std()

    # Add vertical lines
    ax.axvline(
        mean_val,
        color="#E74C3C",
        linestyle="--",
        linewidth=2.5,
        alpha=0.8,
        label=f"Mean: {mean_val:.1f}",
    )
    ax.axvline(
        median_val,
        color="#27AE60",
        linestyle="--",
        linewidth=2.5,
        alpha=0.8,
        label=f"Median: {median_val:.1f}",
    )

    # Add percentile shading
    q25 = data.quantile(0.25)
    q75 = data.quantile(0.75)
    ax.axvspan(q25, q75, alpha=0.1, color="gray", label="IQR")

    # Labels and title
    ax.set_xlabel("Composite Score", fontsize=16, fontweight="bold")
    ax.set_ylabel("Frequency", fontsize=16, fontweight="bold")
    ax.set_title(
        "Distribution of Aging Research Relevance Scores",
        fontsize=18,
        fontweight="bold",
        pad=20,
    )

    # Statistics box
    stats_text = f"n = {len(data):,}\nμ = {mean_val:.1f}\nσ = {std_val:.1f}\nIQR: [{q25:.0f}, {q75:.0f}]"
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        fontsize=13,
        verticalalignment="top",
        bbox=dict(
            boxstyle="round,pad=0.8", facecolor="white", edgecolor="gray", alpha=0.95
        ),
    )

    # Legend
    ax.legend(loc="upper right", fontsize=12, frameon=True, fancybox=True, shadow=True)

    # Grid
    ax.grid(True, alpha=0.2, linestyle="--", axis="y")
    ax.set_axisbelow(True)

    plt.tight_layout()
    return fig


def create_subscore_distributions(df, score_columns):
    """Create overlapping density plots for all subscores."""

    fig, ax = plt.subplots(figsize=(12, 7))

    # Normalize scores to 0-1 range for comparison
    for col in score_columns:
        data = df[col].dropna()

        # Get the actual range for this score
        if col == "q0_plausibility":
            min_val, max_val = -50, 30
        elif col == "q1_experimental_model":
            min_val, max_val = 0, 25
        elif col in ["q2_1_aging_pathways", "q5_1_literature_saturation"]:
            min_val, max_val = 0, 15
        elif col in ["q2_2_age_diseases", "q3_1_experimental_rigor"]:
            min_val, max_val = 0, 10
        elif col == "q4_1_translational_potential":
            min_val, max_val = 0, 5
        elif col == "q5_2_mechanistic_novelty":
            min_val, max_val = 0, 20
        else:
            continue

        # Normalize to 0-1
        normalized = (data - min_val) / (max_val - min_val)

        # Create KDE
        if len(normalized) > 1:
            kde = gaussian_kde(normalized, bw_method="scott")
            x_range = np.linspace(0, 1, 300)
            kde_values = kde(x_range)

            # Plot with custom styling
            ax.plot(
                x_range,
                kde_values,
                label=NAMES[col],
                color=COLORS[col],
                linewidth=2.5,
                alpha=0.85,
            )

            # Fill under curve
            ax.fill_between(x_range, 0, kde_values, color=COLORS[col], alpha=0.15)

    # Styling
    ax.set_xlabel("Normalized Score (0-1)", fontsize=16, fontweight="bold")
    ax.set_ylabel("Density", fontsize=16, fontweight="bold")
    ax.set_title(
        "Distribution of Score Components", fontsize=18, fontweight="bold", pad=20
    )

    # Legend with custom styling
    legend = ax.legend(
        loc="upper right",
        fontsize=11,
        title="Score Components",
        title_fontsize=12,
        frameon=True,
        fancybox=True,
        shadow=True,
        ncol=2,
        columnspacing=1,
        handlelength=2.5,
    )
    legend.get_frame().set_alpha(0.95)

    # Grid
    ax.grid(True, alpha=0.2, linestyle="--")
    ax.set_axisbelow(True)

    # Set x-axis limits
    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    return fig


def get_text_color(bg_color):
    """Determine optimal text color (black or white) based on background color."""
    # Convert hex to RGB
    if isinstance(bg_color, str) and bg_color.startswith("#"):
        hex_color = bg_color.lstrip("#")
        rgb = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    else:
        return "white"  # Default to white for unknown formats

    # Calculate luminance
    luminance = (0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]) / 255

    # Return black for light backgrounds, white for dark
    return "black" if luminance > 0.5 else "white"


def create_correlation_network(df, score_columns):
    """Create an enhanced network graph showing correlations between scores."""

    fig, ax = plt.subplots(figsize=(12, 10), facecolor="#F8F9FA")
    ax.set_facecolor("#F8F9FA")

    # Calculate correlation matrix
    corr_matrix = df[score_columns + ["final_score"]].corr()

    # Create network graph
    G = nx.Graph()

    # Add nodes with correlation to final score as attribute
    for col in score_columns + ["final_score"]:
        if col != "final_score":
            final_corr = corr_matrix.loc[col, "final_score"]
        else:
            final_corr = 1.0
        G.add_node(col, final_correlation=final_corr)

    # Add edges for significant correlations
    threshold = 0.25  # Lowered threshold for more connections
    edge_data = []
    for i, col1 in enumerate(score_columns + ["final_score"]):
        for j, col2 in enumerate(score_columns + ["final_score"]):
            if i < j:  # Upper triangle only
                corr_val = corr_matrix.loc[col1, col2]
                if abs(corr_val) > threshold and not (
                    col1 == "final_score" or col2 == "final_score"
                ):
                    # Don't show edges to final_score to reduce clutter
                    edge_data.append((col1, col2, corr_val))
                elif (col1 == "final_score" or col2 == "final_score") and abs(
                    corr_val
                ) > 0.1:
                    # Show all edges to final_score with lower threshold
                    edge_data.append((col1, col2, corr_val))

    for col1, col2, corr_val in edge_data:
        G.add_edge(col1, col2, weight=abs(corr_val), correlation=corr_val)

    # Create circular layout with final_score at center
    # Arrange other nodes in a circle based on their correlation with final_score
    pos = {}
    if "final_score" in G.nodes():
        # Center position for final_score
        pos["final_score"] = np.array([0, 0])

        # Sort other nodes by correlation with final_score
        other_nodes = [n for n in G.nodes() if n != "final_score"]
        node_correlations = [(n, G.nodes[n]["final_correlation"]) for n in other_nodes]
        node_correlations.sort(key=lambda x: x[1], reverse=True)

        # Arrange in circle
        n_nodes = len(other_nodes)
        radius = 0.8
        for i, (node, _) in enumerate(node_correlations):
            angle = 2 * np.pi * i / n_nodes - np.pi / 2  # Start from top
            pos[node] = np.array([radius * np.cos(angle), radius * np.sin(angle)])
    else:
        pos = nx.circular_layout(G)

    # Draw edges with bezier curves for better aesthetics
    for edge in G.edges(data=True):
        x1, y1 = pos[edge[0]]
        x2, y2 = pos[edge[1]]

        weight = edge[2]["weight"]
        correlation = edge[2]["correlation"]

        # Determine if edge connects to final_score
        is_final_edge = "final_score" in [edge[0], edge[1]]

        # Color and style based on correlation
        if correlation > 0:
            edge_color = "#2ECC71" if not is_final_edge else "#27AE60"
        else:
            edge_color = "#E74C3C" if not is_final_edge else "#C0392B"

        # Draw curved edge for better visibility
        if is_final_edge:
            # Straight lines for final_score connections
            linewidth = 1 + weight * 4
            alpha = 0.3 + weight * 0.5
            linestyle = "-"
        else:
            # Curved lines for inter-component connections
            linewidth = 0.5 + weight * 3
            alpha = 0.2 + weight * 0.3
            linestyle = "--" if weight < 0.5 else "-"

            # Add slight curve
            mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
            offset = 0.1
            mid_x += offset * (y2 - y1)
            mid_y -= offset * (x2 - x1)

            # Draw curved path using quadratic bezier
            t = np.linspace(0, 1, 100)
            bezier_x = (1 - t) ** 2 * x1 + 2 * (1 - t) * t * mid_x + t**2 * x2
            bezier_y = (1 - t) ** 2 * y1 + 2 * (1 - t) * t * mid_y + t**2 * y2
            ax.plot(
                bezier_x,
                bezier_y,
                color=edge_color,
                linewidth=linewidth,
                alpha=alpha,
                linestyle=linestyle,
                zorder=1,
            )
            continue

        # Draw straight edge for final_score connections
        ax.plot(
            [x1, x2],
            [y1, y2],
            color=edge_color,
            linewidth=linewidth,
            alpha=alpha,
            linestyle=linestyle,
            zorder=1,
        )

        # Add correlation value for strong correlations
        if abs(correlation) > 0.5 and not is_final_edge:
            mid_x, mid_y = (x1 + x2) / 2, (y1 + y2) / 2
            ax.text(
                mid_x,
                mid_y,
                f"{correlation:.2f}",
                fontsize=8,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor=edge_color,
                    alpha=0.9,
                ),
                ha="center",
                zorder=4,
            )

    # Draw nodes with enhanced styling
    for node in G.nodes():
        x, y = pos[node]

        # Node size based on correlation with final score
        if node == "final_score":
            node_radius = 0.15
            node_color = "#2C3E50"
            edge_color = "#F39C12"
            edge_width = 4
            text_color = get_text_color(node_color)
            font_weight = "bold"
            font_size = 12
        else:
            final_corr = abs(G.nodes[node]["final_correlation"])
            node_radius = 0.08 + 0.05 * final_corr
            node_color = COLORS.get(node, "#95A5A6")
            edge_color = "white"
            edge_width = 2.5
            text_color = get_text_color(node_color)
            font_weight = "normal"
            font_size = 9 + 2 * final_corr

        # Draw node with gradient effect
        circle = plt.Circle(
            (x, y),
            node_radius,
            color=node_color,
            ec=edge_color,
            linewidth=edge_width,
            zorder=2,
        )
        ax.add_patch(circle)

        # Add subtle shadow
        shadow = plt.Circle(
            (x + 0.01, y - 0.01), node_radius, color="gray", alpha=0.2, zorder=1.5
        )
        ax.add_patch(shadow)

        # Add label with better formatting
        label = NAMES.get(node, node)

        # Smart label wrapping
        if len(label) > 12 and " " in label:
            words = label.split()
            if len(words) == 2:
                label = "\n".join(words)
            elif len(words) > 2:
                mid = len(words) // 2
                label = "\n".join([" ".join(words[:mid]), " ".join(words[mid:])])

        ax.text(
            x,
            y,
            label,
            fontsize=font_size,
            fontweight=font_weight,
            ha="center",
            va="center",
            color=text_color,
            zorder=3,
        )

        # Add correlation value below node (except for final_score)
        if node != "final_score":
            corr_with_final = G.nodes[node]["final_correlation"]
            ax.text(
                x,
                y - node_radius - 0.05,
                f"r={corr_with_final:.2f}",
                fontsize=7,
                ha="center",
                va="top",
                color="#555555",
                style="italic",
                zorder=3,
            )

    # Enhanced legend
    legend_elements = [
        plt.Line2D(
            [0],
            [0],
            color="#27AE60",
            linewidth=3,
            label="Positive correlation",
            alpha=0.8,
        ),
        plt.Line2D(
            [0],
            [0],
            color="#C0392B",
            linewidth=3,
            label="Negative correlation",
            alpha=0.8,
        ),
        plt.Line2D(
            [0],
            [0],
            color="gray",
            linewidth=1,
            linestyle="--",
            label=f"Weak (|r| < 0.5)",
            alpha=0.4,
        ),
        plt.Line2D(
            [0], [0], color="gray", linewidth=3, label="Strong (|r| ≥ 0.5)", alpha=0.8
        ),
        plt.Circle(
            (0, 0),
            0.1,
            facecolor="#F39C12",
            edgecolor="#F39C12",
            label="Composite Score",
            linewidth=2,
        ),
    ]

    legend = ax.legend(
        handles=legend_elements,
        loc="upper left",
        fontsize=10,
        frameon=True,
        fancybox=True,
        shadow=True,
        framealpha=0.95,
        title="Legend",
        title_fontsize=11,
    )
    legend.get_frame().set_facecolor("white")

    # Set limits and remove axes
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.15, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")

    # Add subtle grid in background
    for i in np.linspace(-1.0, 1.0, 9):
        ax.axhline(y=i, color="gray", linewidth=0.1, alpha=0.2, zorder=0)
        ax.axvline(x=i, color="gray", linewidth=0.1, alpha=0.2, zorder=0)

    plt.tight_layout()
    return fig


def main():
    """Main execution function."""

    # Load data
    print("Loading score data...")
    df = pd.read_csv("data/scores_summary.csv")
    print(f"Loaded {len(df)} records")

    # Get score columns (excluding q7)
    score_columns = [
        col
        for col in df.columns
        if col.startswith("q") and col != "q7_experimental_model"
    ]
    print(f"Processing {len(score_columns)} score components (q7 excluded)")

    # Setup theme
    theme = setup_theme()

    # Create output directory
    output_dir = Path("figures")
    output_dir.mkdir(exist_ok=True)

    # 1. Overall score distribution
    print("\nCreating overall score distribution...")
    fig1 = create_overall_score_distribution(df)
    if theme and USE_AQUAREL:
        theme.apply_transforms()

    fig1.savefig(output_dir / "score_distribution_overall.pdf", bbox_inches="tight")
    fig1.savefig(
        output_dir / "score_distribution_overall.png", dpi=600, bbox_inches="tight"
    )
    print(f"Saved: {output_dir}/score_distribution_overall.pdf/png")
    plt.close()

    # 2. Subscore distributions
    print("\nCreating subscore distributions...")
    fig2 = create_subscore_distributions(df, score_columns)
    if theme and USE_AQUAREL:
        theme.apply_transforms()

    fig2.savefig(output_dir / "score_distribution_components.pdf", bbox_inches="tight")
    fig2.savefig(
        output_dir / "score_distribution_components.png", dpi=600, bbox_inches="tight"
    )
    print(f"Saved: {output_dir}/score_distribution_components.pdf/png")
    plt.close()

    # 3. Correlation network
    print("\nCreating correlation network...")
    fig3 = create_correlation_network(df, score_columns)
    if theme and USE_AQUAREL:
        theme.apply_transforms()

    fig3.savefig(output_dir / "score_correlation_network.pdf", bbox_inches="tight")
    fig3.savefig(
        output_dir / "score_correlation_network.png", dpi=600, bbox_inches="tight"
    )
    print(f"Saved: {output_dir}/score_correlation_network.pdf/png")
    plt.close()

    print("\n" + "=" * 60)
    print("All visualizations created successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
