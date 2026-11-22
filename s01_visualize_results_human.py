# HUMAN VERSION of s01_visualize_results.py
# Aggregates per-GSE aging_analysis_sonnet37.csv files from hs_{MODALITY}_sonnet37/
# into a master outlier-marked table used by analysis_1..4_human scripts.
# Mouse version untouched; flip MODALITY below and run twice (rna + dnam).
import os
from pathlib import Path

# ---- HUMAN CONFIG --------------------------------------------------
MODALITY = "rna"  # "rna" or "dnam"  — change and rerun for the other modality
METHOD_DIR = f"hs_{MODALITY}_sonnet37"  # folder to walk for per-GSE CSVs
METHOD_NAME = f"hs_{MODALITY}"  # used in output filenames
INPUT_FILENAME = "aging_analysis_sonnet37.csv"  # per-GSE file pattern
# Clocks to retain. RNA uses a single transcriptomic clock; DNAm uses 5 canonical clocks.
CLOCK_FILTER = {
    "rna": ["Hazard2_All_All_All_scaleddiff_EN_tAge"],
    "dnam": ["PhenoAge", "Horvathv1", "YingDamAge", "DunedinPACE", "GrimAgeV2"],
}[MODALITY]
# --------------------------------------------------------------------

import matplotlib.patheffects as pe

# Data visualization and analysis libraries
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from aquarel import load_theme
from statsmodels.stats.multitest import fdrcorrection

theme = (
    load_theme("scientific")
    .set_grid(draw=False)
    .set_font(size=15)
    .set_ticks(direction="out")
    .set_axis_labels(pad=10)
)
theme.apply()

# Configuration for annotations
ANNOTATIONS_PER_SIDE = 15  # Number of annotations for positive and negative effects
SHOW_GSE_ID = False  # Set to False to hide GSE IDs in annotations
FORCE_REFRESH_ANNOTATIONS = (
    False  # Set to True to override existing manual adjustments and start fresh
)


def normalize_trust_scores(df):
    """Normalize final scores for point size and alpha"""
    df = df.copy()

    # Always use final_score (we've already set it with default value of 50)
    df["final_score"] = df["final_score"].fillna(50)

    # Normalize scores for size and alpha - using smaller dots
    max_score = df["final_score"].max()
    min_score = df["final_score"].min()

    if max_score == min_score:
        df["size"] = 10  # Smaller default size
        df["alpha"] = 0.6
    else:
        # Cubic scale for score sizing
        normalized_score = (df["final_score"] - min_score) / (max_score - min_score)
        df["size"] = normalized_score**3 * 120  # Cubic scale for stronger distinction
        df["alpha"] = (
            0.3 + (df["final_score"] - min_score) / (max_score - min_score) * 0.7
        )

    return df


def add_size_legend(ax, min_score=0, max_score=150, score_type="Final Score"):
    """Add a legend showing the relationship between score and point size"""
    import numpy as np
    from matplotlib.lines import Line2D

    # Adjust score values based on the scale
    if max_score > 100:
        score_values = [50, 70, 90, 110, 130]
    else:
        score_values = [2, 4, 6, 8, 10]

    legend_elements = []

    # Create custom legend elements for each score
    for score in score_values:
        if score >= min_score and score <= max_score:
            normalized = (score - min_score) / (max_score - min_score)
            size = (
                normalized**3 * 120
            )  # Match the actual size calculation from normalize_trust_scores
            # Create a Line2D object that displays as a point
            # markersize needs to be scaled appropriately for display
            # Using sqrt(size/pi) to convert area to radius, then scale for display
            marker_size = np.sqrt(size / np.pi) * 2
            legend_elements.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="w",
                    markerfacecolor="gray",
                    markersize=marker_size,
                    alpha=1,
                    linestyle="None",
                    label=f"{score_type}: {score}",
                )
            )

    # Create the size legend
    size_legend = ax.legend(
        handles=legend_elements,
        title=f"{score_type} (Size)",
        loc="center left",
        bbox_to_anchor=(0.85, 0.9),  # Lower position to avoid overlap
        frameon=True,
        fancybox=False,
        shadow=False,
        fontsize=10,
    )

    # Add the legend as an artist to ensure it persists

    return size_legend


def lump_categories(df, min_freq=3, n_top=5):
    """Lump all categories except the top n into 'Other'"""
    category_counts = df["condition_category"].value_counts()
    top_categories = category_counts.nlargest(n_top).index
    df = df.copy()
    df.loc[~df["condition_category"].isin(top_categories), "condition_category"] = (
        "Other"
    )
    return df


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


def create_volcano_plot(df, title, output_file):
    if df.empty:
        print(f"No data available for volcano plot: {title}")
        return

    # Create figure with more space for legend
    plt.figure(figsize=(15, 8))

    # Calculate FDR and transform effect size to log scale
    df = calculate_fdr(df)
    df["log_effect"] = np.log2(np.abs(df["log2FoldChange"]) + 1) * np.sign(
        df["log2FoldChange"]
    )

    # Lump rare categories
    df = lump_categories(df)

    # Define significance threshold using FDR
    significant = df["fdr"] < 0.05

    # Get unique intervention types and sort them
    interventions = sorted(df["condition_category"].unique())
    colors = sns.color_palette("husl", len(interventions))
    color_dict = dict(zip(interventions, colors))

    # Normalize trust scores
    df = normalize_trust_scores(df)

    # Create main scatter plot with two-layer approach
    for intervention in interventions:
        mask = df["condition_category"] == intervention
        if mask.any():
            # Layer 1: Black outline (larger size)
            plt.scatter(
                df.loc[mask, "log_effect"],
                -np.log10(df.loc[mask, "fdr"]),
                s=df.loc[mask, "size"] * 1.5,
                c="black",
                alpha=1,
                zorder=1,
            )
            # Layer 2: Colored data points
            plt.scatter(
                df.loc[mask, "log_effect"],
                -np.log10(df.loc[mask, "fdr"]),
                s=df.loc[mask, "size"],
                c=color_dict[intervention],
                alpha=df.loc[mask, "alpha"],
                label=intervention,
                zorder=2,
            )

    # Add threshold lines
    plt.axhline(
        y=-np.log10(0.05), color="gray", linestyle="--", alpha=0.5, label="FDR = 0.05"
    )
    plt.axvline(
        x=np.log2(2), color="gray", linestyle="--", alpha=0.5, label="2-fold change"
    )
    plt.axvline(x=-np.log2(2), color="gray", linestyle="--", alpha=0.5)

    # Add labels and title
    plt.xlabel("Effect Size")
    plt.ylabel("-log10(FDR)")
    plt.title(title)

    # Add annotations for top results with text adjustment
    df["significance"] = -np.log10(df["fdr"])
    df["significance"] = df["significance"].replace(
        [np.inf, -np.inf], df["significance"].replace([np.inf, -np.inf], np.nan).max()
    )

    # Select points to annotate:
    # 1. Most significant points
    # 2. Points with largest effect sizes (both positive and negative)
    # 3. Points with high trust scores
    # 4. Points that are both significant and have large effect sizes
    # Calculate significance * abs(log_effect) score for ranking
    df["score"] = df["significance"] * np.abs(df["log_effect"])

    # Filter by final_score > 80 and take top 20 by combined score
    to_annotate = df[df["final_score"] > 80].nlargest(20, "score").drop_duplicates()

    # If we don't have enough points with high score, lower threshold
    if len(to_annotate) < 5:
        print(f"  Only {len(to_annotate)} points with final_score > 80, trying > 70")
        to_annotate = df[df["final_score"] > 70].nlargest(20, "score").drop_duplicates()
    if len(to_annotate) < 5:
        print(
            f"  Only {len(to_annotate)} points with final_score > 70, using top scores"
        )
        to_annotate = df.nlargest(20, "final_score").drop_duplicates()

    texts = []
    for _, row in to_annotate.iterrows():
        # Create more informative labels
        # Always use intervention name (we've loaded it from scores_summary)
        label_text = (
            row["intervention"][:40]
            if pd.notna(row.get("intervention", "")) and row.get("intervention", "")
            else "Unknown"
        )

        if SHOW_GSE_ID:
            label = f"{row['gse_id']}\n{label_text}"
        else:
            label = label_text

        texts.append(
            plt.text(
                row["log_effect"],
                row["significance"],
                label,
                fontsize=10,
                ha="center",
                va="bottom",
                path_effects=[
                    pe.withStroke(linewidth=3, foreground="white"),
                    pe.Normal(),
                ],
            )
        )

    # Adjust text positions to avoid overlap with more iterations and force
    if texts:
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle="->", color="gray", alpha=0.5, lw=0.5),
            expand_points=(2, 2),  # Push text further from points
            force_points=(0.8, 0.8),
            lim=500,  # More iterations for better positioning
        )

    # Create a more prominent legend for intervention types
    # Move legend outside the plot
    ax = plt.gca()

    # First create and add the main legend
    main_legend = ax.legend(
        title="Intervention Types",
        bbox_to_anchor=(1.25, 1),
        loc="upper left",
        borderaxespad=0,
        frameon=True,
        fancybox=True,
        shadow=True,
        fontsize=10,
    )

    # Add the main legend as an artist so it persists
    ax.add_artist(main_legend)

    # Now add the size legend (this becomes the 'current' legend)
    size_legend = add_size_legend(ax)

    # Add statistics text box with more information
    stats_text = (
        f"Total points: {len(df)}\n"
        f"Significant (FDR<0.05): {sum(significant)} ({sum(significant)/len(df)*100:.1f}%)\n"
        f"Final score range: {df['final_score'].min():.1f}-{df['final_score'].max():.1f}\n"
        f"Effect size range: {df['log2FoldChange'].min():.2f} to {df['log2FoldChange'].max():.2f}\n"
        f"Categories: {len(interventions)}"
    )

    plt.text(
        0.95,
        0.95,
        stats_text,
        transform=plt.gca().transAxes,
        verticalalignment="bottom",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="white", edgecolor="black", alpha=0.8),
        fontsize=9,
    )

    # Adjust layout to prevent text cutoff
    plt.tight_layout()
    theme.apply_transforms()
    # Save in both PDF and PNG formats
    plt.savefig(output_file.replace(".png", ".pdf"), bbox_inches="tight")
    plt.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close()


def load_data(file_path):
    """Load and process aging analysis results"""
    try:
        df = pd.read_csv(file_path)

        required_columns = [
            "gse_id",
            "comparison_name",
            "condition_name",
            "condition_category",
            "clock",
            "effect_size",
            "p_value",
            "trust_score",
            "direction",
        ]

        if not all(col in df.columns for col in required_columns):
            print(f"Missing required columns in file: {file_path}")
            return None

        # Convert numeric columns
        df["effect_size"] = pd.to_numeric(df["effect_size"], errors="coerce")
        df["p_value"] = pd.to_numeric(df["p_value"], errors="coerce")
        df["trust_score"] = pd.to_numeric(df["trust_score"], errors="coerce")
        df["direction"] = pd.to_numeric(df["direction"], errors="coerce")

        # Calculate log2FoldChange correctly: abs(effect) * direction
        df["log2FoldChange"] = np.abs(df["effect_size"]) * df["direction"]

        # Rename p_value column
        df = df.rename(columns={"p_value": "pvalue"})

        # Drop rows with missing critical values
        df = df.dropna(subset=["log2FoldChange", "pvalue"])

        # Ensure valid p-values
        df = df[df["pvalue"].between(0, 1)]

        if len(df) == 0:
            print(f"No valid data in file: {file_path}")
            return None

        return df

    except Exception as e:
        print(f"Error reading file {file_path}: {str(e)}")
        return None


def process_directory(base_dir):
    """Process all per-GSE analysis CSVs in directory (filename set by INPUT_FILENAME).

    CONTAMINATION FILTER: ~10-13% of the agent's per-GSE CSV files carry rows
    whose `gse_id` column disagrees with the parent folder name (evidence of
    multi-trace-rerun accidents in the agent output). These rows systematically
    ship wrong-direction effect sizes. Drop any row where `gse_id` != folder name.
    """
    results = []
    gse_ids = set()
    total_rows = 0
    dropped_rows = 0
    contaminated_files = 0

    for root, _, files in os.walk(base_dir):
        if INPUT_FILENAME in files:
            file_path = os.path.join(root, INPUT_FILENAME)
            df = load_data(file_path)

            if df is not None:
                folder_gse = os.path.basename(root)
                before = len(df)
                keep_mask = df["gse_id"].astype(str).str.strip() == folder_gse
                dropped_this_file = int((~keep_mask).sum())
                if dropped_this_file > 0:
                    contaminated_files += 1
                df = df[keep_mask].copy()
                total_rows += before
                dropped_rows += dropped_this_file

                if len(df) > 0:
                    results.append(df)
                    gse_ids.update(df["gse_id"].unique())

    if total_rows > 0:
        print(
            f"CONTAMINATION FILTER: dropped {dropped_rows:,}/{total_rows:,} rows "
            f"({100*dropped_rows/total_rows:.1f}%) across {contaminated_files:,} "
            f"files where gse_id != folder name"
        )
    return results, gse_ids


def load_scores_and_interventions(scores_file="data/scores_summary.csv"):
    """Load scores summary with interventions and final scores"""
    try:
        scores_df = pd.read_csv(scores_file)
        # Create a mapping dictionary using dataset ID as key
        # dataset ID format: GSE{number}-{row_index_in_main_csv}
        scores_mapping = {}

        for _, row in scores_df.iterrows():
            dataset_id = row["dataset"]  # e.g., GSE100035-0
            scores_mapping[dataset_id] = {
                "intervention": row.get("intervention", ""),
                "final_score": row.get("final_score", 0),
            }

        print(f"Loaded {len(scores_mapping)} dataset entries from scores file")
        return scores_mapping
    except Exception as e:
        print(f"Error loading scores file: {e}")
        return {}


## Run
# Create output directory
os.makedirs("figures", exist_ok=True)

# Load scores and interventions mapping
scores_mapping = load_scores_and_interventions(scores_file=f"data/scores_summary_{METHOD_NAME}.csv")
print(f"Loaded {len(scores_mapping)} dataset scores and interventions")

# Process each method's directory

method = METHOD_NAME  # used for combined_df / main_results output filenames

if os.path.exists(f"./result/combined_df_{method}.csv"):
    combined_df = pd.read_csv(f"./result/combined_df_{method}.csv")
    combined_df.iloc[0]

    # Load case_mapping.csv which contains matched dataset_id (mouse-only; skip if absent)
    if os.path.exists("data/case_mapping.csv"):
        mapping_df = pd.read_csv("data/case_mapping.csv")
    else:
        print("No data/case_mapping.csv — skipping dataset_id merge (human pipeline)")
        mapping_df = pd.DataFrame()
    # Try to merge in dataset_id directly using available keys
    preferred_keys = ["gse_id", "comparison_name", "condition_name"]
    merge_keys = [
        k
        for k in preferred_keys
        if k in combined_df.columns and k in mapping_df.columns
    ]
    if "case_id" in mapping_df.columns and len(merge_keys) > 0:
        combined_df = combined_df.merge(
            mapping_df[merge_keys + ["case_id"]].drop_duplicates(),
            on=merge_keys,
            how="left",
        )
        # Assign dataset_id from case_id
        if "case_id" in combined_df.columns:
            combined_df["dataset_id"] = combined_df["case_id"]
    else:
        # Fallback: if mapping_df has row_index, construct dataset_id from it
        if "row_index" in mapping_df.columns and len(merge_keys) > 0:
            combined_df = combined_df.merge(
                mapping_df[merge_keys + ["row_index"]].drop_duplicates(),
                on=merge_keys,
                how="left",
            )
            combined_df["dataset_id"] = combined_df.apply(
                lambda row: (
                    f"{row['gse_id']}-{int(row['row_index'])}"
                    if pd.notna(row.get("row_index"))
                    else ""
                ),
                axis=1,
            )
        else:
            # If no mapping available, ensure dataset_id exists as empty
            if "dataset_id" not in combined_df.columns:
                combined_df["dataset_id"] = ""

    # Create dataset_id for display/tracking purposes
    combined_df.iloc[0]

    # Function to get intervention and score with fallback logic
    def get_intervention_and_score(row):
        # First try exact dataset_id match
        dataset_id = row["dataset_id"]
        if dataset_id in scores_mapping:
            return pd.Series(
                {
                    "intervention": scores_mapping[dataset_id]["intervention"],
                    "final_score": scores_mapping[dataset_id]["final_score"],
                }
            )

        return pd.Series({"intervention": "", "final_score": 50})

    # Apply mapping to get interventions and final scores
    combined_df[["intervention", "final_score"]] = combined_df.apply(
        get_intervention_and_score, axis=1
    )

    # Use final_score, default to 50 if not available
    combined_df["original_trust_score"] = combined_df.get("trust_score", 50).copy()
    combined_df["final_score"] = combined_df["final_score"].fillna(50)
    combined_df["trust_score"] = combined_df["final_score"]

    # Drop the temporary row_index column
    combined_df = combined_df.drop(columns=["row_index"], errors="ignore")
else:
    results, gse_ids = process_directory(METHOD_DIR)

    print(len(results))

    combined_df = pd.concat(results, ignore_index=True)

    # Keep only the clocks defined in CLOCK_FILTER (human-specific: 1 for RNA, 5 for DNAm)
    combined_df = combined_df[combined_df["clock"].isin(CLOCK_FILTER)]
    print(f"After clock filter ({CLOCK_FILTER}): {len(combined_df):,} rows")

    # Load case_mapping.csv to obtain dataset_id directly (mouse-only; skip if absent)
    if os.path.exists("data/case_mapping.csv"):
        mapping_df = pd.read_csv("data/case_mapping.csv")
    else:
        print("No data/case_mapping.csv — skipping dataset_id merge (human pipeline)")
        mapping_df = pd.DataFrame()

    # Merge to get row indices
    preferred_keys = ["gse_id", "comparison_name", "condition_name"]
    merge_keys = [
        k
        for k in preferred_keys
        if k in combined_df.columns and k in mapping_df.columns
    ]
    if "case_id" in mapping_df.columns and len(merge_keys) > 0:
        combined_df = combined_df.merge(
            mapping_df[merge_keys + ["case_id"]].drop_duplicates(),
            on=merge_keys,
            how="left",
        )
        if "case_id" in combined_df.columns:
            combined_df["dataset_id"] = combined_df["case_id"]
    else:
        if "row_index" in mapping_df.columns and len(merge_keys) > 0:
            combined_df = combined_df.merge(
                mapping_df[merge_keys + ["row_index"]].drop_duplicates(),
                on=merge_keys,
                how="left",
            )
            combined_df["dataset_id"] = combined_df.apply(
                lambda row: (
                    f"{row['gse_id']}-{int(row['row_index'])}"
                    if pd.notna(row.get("row_index"))
                    else ""
                ),
                axis=1,
            )
        else:
            combined_df["dataset_id"] = ""

    # Now map interventions and final scores using the dataset_id
    def get_intervention_and_score(dataset_id):
        if dataset_id in scores_mapping:
            return pd.Series(
                {
                    "intervention": scores_mapping[dataset_id]["intervention"],
                    "final_score": scores_mapping[dataset_id]["final_score"],
                }
            )
        return pd.Series({"intervention": "", "final_score": np.nan})

    # Apply mapping to get interventions and final scores
    combined_df[["intervention", "final_score"]] = combined_df["dataset_id"].apply(
        get_intervention_and_score
    )

    # Use final_score, default to 50 if not available
    combined_df["original_trust_score"] = combined_df["trust_score"].copy()
    combined_df["final_score"] = combined_df["final_score"].fillna(50)
    combined_df["trust_score"] = combined_df["final_score"]

    # Drop the temporary row_index column before saving
    combined_df = combined_df.drop(columns=["row_index"], errors="ignore")

    combined_df.to_csv(f"./result/combined_df_{method}.csv", index=False)

# HUMAN PATCH (runs in BOTH branches): overwrite intervention + final_score using
# the merged-scores file produced by s03_merge_scores_human.py. Matches by
# (gse_id, comparison, condition), robust to enumeration differences between
# s01 and the scorer. Without this, every human row gets final_score=50 default
# and Figure 3 point sizing + annotations degenerate.
_scored_path = f"./result/main_results_{method}_with_scores.csv"
if os.path.exists(_scored_path):
    _scored = pd.read_csv(_scored_path, low_memory=False,
                          usecols=["gse_id", "comparison_name", "condition_name",
                                   "final_score", "intervention"])
    _scored = _scored.drop_duplicates(subset=["gse_id", "comparison_name", "condition_name"])
    _merge_keys = ["gse_id", "comparison_name", "condition_name"]
    combined_df = combined_df.drop(columns=["final_score", "intervention"], errors="ignore")
    combined_df = combined_df.merge(_scored, on=_merge_keys, how="left")
    combined_df["final_score"] = combined_df["final_score"].fillna(50)
    combined_df["intervention"] = combined_df["intervention"].fillna("")
    combined_df["trust_score"] = combined_df["final_score"]
    if "original_trust_score" not in combined_df.columns:
        combined_df["original_trust_score"] = combined_df["trust_score"].copy()
    print(f"HUMAN PATCH: merged final_score for {combined_df['final_score'].notna().sum():,} rows; "
          f"{(combined_df['final_score']!=50).sum():,} with real scores")

# Filter out inappropriate Environment entries
print("\nFiltering out inappropriate Environment category entries...")
# These are tissue/cell types, not environmental treatments
# inappropriate_environment_conditions = [
#     "plasma cells",
#     "cerebellum",
#     "Blood",
#     "Bone marrow plasma cells",
#     "Colonic Biopsy",
#     "Colon",
#     "Colon mucosa",
#     "Brain Region",
#     "Brainstem Astrocytes",
#     "Cortex Astrocytes",
#     "Cell Culture",
#     "Cell culture",
#     "B cell culture",
#     "B cell differentiation",
#     "Spleen",
#     "Lymph nodes",
#     "Thymus",
#     "Bone marrow",
#     "Liver",
#     "Kidney",
#     "Heart",
#     "Lung",
#     "Muscle",
#     "plasmablasts",
#     "Ileum macrophages day 12",
#     "macrophages",
#     "neurons",
#     "astrocytes",
#     "microglia",
#     "T cells",
#     "B cells",
# ]
# before_filter = len(combined_df)
# combined_df = combined_df[
#     ~(
#         (combined_df["condition_category"] == "Environment")
#         & (combined_df["condition_name"].isin(inappropriate_environment_conditions))
#     )
# ]
# after_filter = len(combined_df)
# print(f"Removed {before_filter - after_filter} inappropriate Environment entries")

##
# Remove outliers
print(f"\nRemoving outliers for {method}...")

# Before outlier removal
print(f"Data points before outlier removal: {len(combined_df)}")

# Define outlier thresholds for log2FoldChange (using IQR method)
q1_log2fc = combined_df["log2FoldChange"].quantile(0.25)
q3_log2fc = combined_df["log2FoldChange"].quantile(0.75)
iqr_log2fc = q3_log2fc - q1_log2fc

lower_bound_log2fc = q1_log2fc - 10 * iqr_log2fc
upper_bound_log2fc = q3_log2fc + 10 * iqr_log2fc

# Define outlier thresholds for p-values
# For p-values, we'll use -log10(pvalue) to better identify outliers
combined_df["neg_log10_pvalue"] = -np.log10(combined_df["pvalue"])
combined_df = calculate_fdr(combined_df)

q1_pval = combined_df["neg_log10_pvalue"].quantile(0.25)
q3_pval = combined_df["neg_log10_pvalue"].quantile(0.75)
iqr_pval = q3_pval - q1_pval

lower_bound_neg_log10_pval = q1_pval - 10 * iqr_pval
upper_bound_neg_log10_pval = q3_pval + 10 * iqr_pval

# Convert back to p-values for filtering
# Note: higher -log10(p) means smaller p-value
upper_bound_pval = 10 ** (-lower_bound_neg_log10_pval)
lower_bound_pval = 10 ** (-upper_bound_neg_log10_pval)

# Ensure bounds are within valid p-value range
lower_bound_pval = max(0, lower_bound_pval)
upper_bound_pval = min(1, upper_bound_pval)

## Filter out both p-value and effect size outliers
filtered_df = combined_df[
    (combined_df["log2FoldChange"] >= lower_bound_log2fc)
    & (combined_df["log2FoldChange"] <= upper_bound_log2fc)
    & (combined_df["pvalue"] >= lower_bound_pval)
    & (combined_df["pvalue"] <= upper_bound_pval)
]

# After outlier removal
print(f"Data points after outlier removal: {len(filtered_df)}")
print(f"Removed {len(combined_df) - len(filtered_df)} outliers")
print(f"Log2FC bounds: [{lower_bound_log2fc:.2f}, {upper_bound_log2fc:.2f}]")
print(f"P-value bounds: [{lower_bound_pval:.2e}, {upper_bound_pval:.2e}]")

# Create main result table with outlier flags
print("\nCreating main result table with outlier flags...")

# Create main result table
result_df = combined_df.copy()


# Determine outlier status for each entry
def get_outlier_status(row):
    effect_outlier = (row["log2FoldChange"] < lower_bound_log2fc) | (
        row["log2FoldChange"] > upper_bound_log2fc
    )
    pvalue_outlier = (row["pvalue"] < lower_bound_pval) | (
        row["pvalue"] > upper_bound_pval
    )

    if effect_outlier and pvalue_outlier:
        return "Both Effect Size & P-value Outlier"
    elif effect_outlier:
        return "Effect Size Outlier"
    elif pvalue_outlier:
        return "P-value Outlier"
    else:
        return "Normal"


result_df["outlier_status"] = result_df.apply(get_outlier_status, axis=1)

# Add useful analytical columns
result_df["abs_log2FoldChange"] = np.abs(result_df["log2FoldChange"])
result_df["neg_log10_pvalue"] = -np.log10(result_df["pvalue"])

# Sort by outlier status, then by significance
result_df = result_df.sort_values(
    ["outlier_status", "neg_log10_pvalue", "abs_log2FoldChange"],
    ascending=[True, False, False],
)

# Save the main result table (includes ALL data with outlier flags)
result_table_path = f"result/main_results_{method}_outlier_marked.csv"
result_df.to_csv(result_table_path, index=False)

# Print summary statistics
print(f"\nResult Table Summary:")
print(f"Total entries: {len(result_df)}")
outlier_counts = result_df["outlier_status"].value_counts()
for status, count in outlier_counts.items():
    print(f"  {status}: {count} ({count/len(result_df)*100:.1f}%)")

print(
    f"\nSaved main results table (all data with outlier flags) to: {result_table_path}"
)

# Create volcano plot with filtered data
## Create category-specific volcano plots for top categories
print("\nCreating category-specific volcano plots...")

# Get top 4 categories by count
top_categories = (
    filtered_df["condition_category"].value_counts().nlargest(5).index.tolist()
)
print(f"Top 5 categories: {', '.join(top_categories)}")

# Create a custom colormap for trust scores
trust_cmap = plt.cm.YlOrRd

##
# Create separate plots for each category
for category in top_categories:

    category_df = filtered_df[filtered_df["condition_category"] == category]
    print(f"Category '{category}': {len(category_df)} data points")

    # Highlight significant points - ALL with FDR < 0.05
    significant_df = category_df[category_df["fdr"] < 0.05]

    # For annotations only - filter by original trust_score > 6 AND final_score >= 70
    annotate_candidates = significant_df[
        (significant_df.get("trust_score", 0) > 6)
        & (significant_df["final_score"] >= 70)
    ].copy()
    if len(annotate_candidates) > 0:
        annotate_candidates["score"] = -annotate_candidates[
            "log2FoldChange"
        ] * np.log10(annotate_candidates["fdr"])
        annotate_candidates["abs_score"] = np.abs(annotate_candidates["score"])
        annotate_candidates["score"] = annotate_candidates["score"].replace(
            [np.inf, -np.inf], 0
        )

        annotate_df = annotate_candidates.sort_values("abs_score", ascending=False)
        annotate_df = annotate_df.drop_duplicates(subset=["gse_id"], keep="first")
        annotate_df = annotate_df.sort_values("score", ascending=False)
        top_n = annotate_df.head(ANNOTATIONS_PER_SIDE)
        if category == "Disease":
            bottom_n = annotate_df.tail(0)
        else:
            bottom_n = annotate_df.tail(ANNOTATIONS_PER_SIDE)
        annotate_df = pd.concat([top_n, bottom_n])
    else:
        annotate_df = pd.DataFrame()

    # Save significant results for reference
    if len(significant_df) > 0:
        significant_df["score"] = -significant_df["log2FoldChange"] * np.log10(
            significant_df["fdr"]
        )
        significant_df["abs_score"] = np.abs(significant_df["score"])
        significant_df["score"] = significant_df["score"].replace([np.inf, -np.inf], 0)
        significant_df.to_csv(f"result/significant_df_{category}.csv", index=False)

    # Create a new figure for this category
    plt.figure(figsize=(16, 8))
    ax = plt.gca()

    # Normalize trust scores for this category
    category_df = normalize_trust_scores(category_df)

    # Plot all points in gray first (no outline for background points)
    ax.scatter(
        category_df["log2FoldChange"],
        -np.log10(category_df["fdr"]),
        s=5,
        color="lightgray",
        alpha=0.3,
        label="Not significant",
        zorder=1,
    )

    # Normalize trust scores for significant points to get proper sizing
    significant_df = normalize_trust_scores(significant_df)

    # Color points by direction and trust score with two-layer approach
    if len(significant_df) > 0:
        # Layer 1: Black outline for significant points
        ax.scatter(
            significant_df["log2FoldChange"],
            -np.log10(significant_df["fdr"]),
            s=significant_df["size"] * 1.5,
            c="black",
            alpha=1,
            zorder=2,
        )
        # Layer 2: Colored significant points
        points = ax.scatter(
            significant_df["log2FoldChange"],
            -np.log10(significant_df["fdr"]),
            s=significant_df["size"],
            c=significant_df["log2FoldChange"],
            cmap="RdBu_r",
            alpha=0.9,
            edgecolor="none",
            vmin=-np.abs(significant_df["log2FoldChange"]).max(),
            vmax=np.abs(significant_df["log2FoldChange"]).max(),
            zorder=3,
        )

    # Add labels for top significant results
    top_results = annotate_df
    texts = []
    for _, result in top_results.iterrows():
        # Always use intervention name
        label_text = (
            result["intervention"][:20]
            if pd.notna(result.get("intervention", ""))
            and result.get("intervention", "")
            else "Unknown"
        )

        if SHOW_GSE_ID:
            label = f"{result['gse_id']}\n{label_text}"
        else:
            label = label_text
        texts.append(
            ax.text(
                result["log2FoldChange"],
                -np.log10(result["fdr"]),
                label,
                fontsize=12,
                ha="center",
                va="bottom",
                path_effects=[
                    pe.withStroke(linewidth=2, foreground="white"),
                    pe.Normal(),
                ],
            )
        )

    # Adjust text to avoid overlaps
    if texts:
        adjust_text(
            texts,
            arrowprops=dict(arrowstyle="->", color="black", alpha=0.5, lw=0.5),
            expand_points=(2, 2),
            force_points=(0.8, 0.8),
            lim=500,
        )

    # Add threshold lines
    ax.axhline(
        y=-np.log10(0.05), color="gray", linestyle="--", alpha=0.5, label="FDR = 0.05"
    )
    ax.axvline(
        x=np.log2(2), color="gray", linestyle="--", alpha=0.5, label="2-fold change"
    )
    ax.axvline(x=-np.log2(2), color="gray", linestyle="--", alpha=0.5)

    # Set axis labels and title
    ax.set_xlabel("Effect Size")
    ax.set_ylabel("-log10(FDR)")
    plt.title(f"{category} (n={len(category_df)})")

    # Add colorbar for trust scores
    if len(significant_df) > 0:
        cbar = plt.colorbar(points, ax=ax)
        cbar.set_label("Effect Size")

    # Add size legend for final scores
    min_score = category_df["final_score"].min()
    max_score = category_df["final_score"].max()
    size_legend = add_size_legend(
        ax, min_score=min_score, max_score=max_score, score_type="Final Score"
    )

    ax.add_artist(size_legend)

    # Adjust layout to prevent text cutoff
    plt.tight_layout()

    theme.apply_transforms()

    # Save the plot in both PDF and PNG formats
    plt.savefig(
        f"figures/volcano_plot_{method}_{category.lower().replace(' ', '_')}.pdf",
        bbox_inches="tight",
        facecolor="white",
    )
    plt.savefig(
        f"figures/volcano_plot_{method}_{category.lower().replace(' ', '_')}.png",
        dpi=600,
        bbox_inches="tight",
        facecolor="white",
    )
    plt.close()

    print(f"  Saved plot for {category}")


print("\nAll category-specific volcano plots have been saved to the figures directory.")

## Create combined multi-panel figure for Nature publication
print("\nCreating combined multi-panel figure for Nature publication...")

# Determine layout based on number of categories
n_categories = len(top_categories)
if n_categories <= 4:
    nrows, ncols = 2, 2
elif n_categories <= 6:
    nrows, ncols = 2, 3
else:
    nrows, ncols = 3, 3

# Set Arial font for Nature publication standards
import matplotlib as mpl

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["font.sans-serif"] = ["Arial"]

# Create figure with appropriate size (Nature double column width is ~180mm)
fig_width = 180 / 12  # Convert mm to inches
fig_height = fig_width * (nrows / ncols) * 0.85  # Maintain good aspect ratio
fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))
axes = axes.flatten()

# Panel labels
panel_labels = ["a", "b", "c", "d", "e", "f", "g", "h", "i"]

# Create plots for each category
for idx, category in enumerate(top_categories):
    if idx >= len(axes):
        break

    ax = axes[idx]

    category_df = filtered_df[filtered_df["condition_category"] == category]
    print(f"  Panel {panel_labels[idx]}: {category} (n={len(category_df)})")

    # Highlight significant points
    significant_df = category_df[category_df["fdr"] < 0.05]

    # For annotations - filter by trust_score > 6 AND final_score >= 70
    annotate_candidates = significant_df[
        (significant_df.get("trust_score", 0) > 6)
        & (significant_df["final_score"] >= 70)
    ].copy()

    if len(annotate_candidates) > 0:
        annotate_candidates["score"] = -annotate_candidates[
            "log2FoldChange"
        ] * np.log10(annotate_candidates["fdr"])
        annotate_candidates["abs_score"] = np.abs(annotate_candidates["score"])
        annotate_candidates["score"] = annotate_candidates["score"].replace(
            [np.inf, -np.inf], 0
        )
        annotate_df = annotate_candidates.sort_values("abs_score", ascending=False)
        annotate_df = annotate_df.drop_duplicates(subset=["gse_id"], keep="first")
        annotate_df = annotate_df.sort_values("score", ascending=False)

        # Reduce annotations for combined figure
        top_n = annotate_df.head(5)
        if category == "Disease":
            bottom_n = annotate_df.tail(0)
        else:
            bottom_n = annotate_df.tail(5)
        annotate_df = pd.concat([top_n, bottom_n])
    else:
        annotate_df = pd.DataFrame()

    # Normalize trust scores
    category_df = normalize_trust_scores(category_df)
    significant_df = normalize_trust_scores(significant_df)

    # Plot non-significant points
    ax.scatter(
        category_df["log2FoldChange"],
        -np.log10(category_df["fdr"]),
        s=3,
        color="lightgray",
        alpha=0.3,
        zorder=1,
    )

    # Plot significant points with outline
    if len(significant_df) > 0:
        # Black outline
        ax.scatter(
            significant_df["log2FoldChange"],
            -np.log10(significant_df["fdr"]),
            s=significant_df["size"] * 1.2,
            c="black",
            alpha=1,
            zorder=2,
        )
        # Colored points
        points = ax.scatter(
            significant_df["log2FoldChange"],
            -np.log10(significant_df["fdr"]),
            s=significant_df["size"] * 0.8,
            c=significant_df["log2FoldChange"],
            cmap="RdBu_r",
            alpha=0.9,
            edgecolor="none",
            vmin=-np.abs(significant_df["log2FoldChange"]).max(),
            vmax=np.abs(significant_df["log2FoldChange"]).max(),
            zorder=3,
        )

    # Add annotations
    texts = []
    for _, result in annotate_df.iterrows():
        label_text = (
            result["intervention"][:15]
            if pd.notna(result.get("intervention", ""))
            and result.get("intervention", "")
            else "Unknown"
        )

        texts.append(
            ax.text(
                result["log2FoldChange"],
                -np.log10(result["fdr"]),
                label_text,
                fontsize=7,
                ha="center",
                va="bottom",
                path_effects=[
                    pe.withStroke(linewidth=1, foreground="white"),
                    pe.Normal(),
                ],
            )
        )

    # Adjust text to avoid overlaps
    if texts:
        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="->", color="black", alpha=0.5, lw=0.3),
            expand_points=(1.5, 1.5),
            force_points=(0.5, 0.5),
            lim=300,
        )

    # Add threshold lines
    ax.axhline(
        y=-np.log10(0.05), color="gray", linestyle="--", alpha=0.5, linewidth=0.5
    )
    ax.axvline(x=np.log2(2), color="gray", linestyle="--", alpha=0.5, linewidth=0.5)
    ax.axvline(x=-np.log2(2), color="gray", linestyle="--", alpha=0.5, linewidth=0.5)

    # Set labels with appropriate font sizes for Nature
    ax.set_xlabel("Effect Size", fontsize=8)
    ax.set_ylabel("-log10(FDR)", fontsize=8)
    ax.set_title(f"{category} (n={len(category_df)})", fontsize=9, pad=5)
    ax.tick_params(labelsize=7)

    # Add bold panel label in Arial
    ax.text(
        -0.15,
        1.05,
        panel_labels[idx],
        transform=ax.transAxes,
        fontsize=18,
        fontweight="bold",
        va="top",
        ha="right",
        family="Arial",
    )

# Remove extra subplots if any
for idx in range(n_categories, len(axes)):
    fig.delaxes(axes[idx])

# Adjust layout
plt.tight_layout()

# Save in both PDF and PNG formats with high resolution
pdf_path = f"figures/volcano_plot_combined_categories_nature_{method}.pdf"
png_path = f"figures/volcano_plot_combined_categories_nature_{method}.png"

plt.savefig(
    pdf_path,
    bbox_inches="tight",
    facecolor="white",
    dpi=600,
)
plt.savefig(
    png_path,
    bbox_inches="tight",
    facecolor="white",
    dpi=600,
)
plt.close()

print(f"  Saved combined figure to:")
print(f"    PDF: {pdf_path}")
print(f"    PNG: {png_path}")

# Reset font to default
mpl.rcParams.update(mpl.rcParamsDefault)
theme.apply()

print("Combined multi-panel figure created successfully!")

## Add a plot for all categories only show annotation for Drug, Disease, and Genetics (Use different annotation colors to annotate each and add a legend for that.)

print("\nCreating combined plot with selective annotations...")

# Import pypalettes for consistent coloring
from pypalettes import load_cmap

# Create figure with adjusted size for better spacing
# Halve height for DNAm where the data span is narrower
_combined_fig_height = 5 if MODALITY == "dnam" else 10
plt.figure(figsize=(18, _combined_fig_height))
ax = plt.gca()

# Normalize trust scores and lump categories (top 5 to ensure Disease is included)
plot_df = normalize_trust_scores(filtered_df.copy())
plot_df = lump_categories(plot_df, n_top=5)  # Changed to top 5 to include Disease

# Save the complete data table for combined volcano plot reproduction
# This includes all data points that will be plotted
combined_plot_data = plot_df.copy()

# Add additional computed columns for plot reproduction
combined_plot_data["neg_log10_fdr"] = -np.log10(combined_plot_data["fdr"])
combined_plot_data["is_significant"] = combined_plot_data["fdr"] < 0.05
combined_plot_data["annotation_score"] = combined_plot_data["neg_log10_fdr"] * np.abs(
    combined_plot_data["log2FoldChange"]
)

# Sort by significance and score for easy reference
combined_plot_data = combined_plot_data.sort_values(
    ["is_significant", "annotation_score", "final_score"],
    ascending=[False, False, False],
)

# Select relevant columns for saving
columns_to_save = [
    "gse_id",  # Original GSE ID
    "dataset_id",  # Full dataset ID from scores_summary (e.g., GSE100035-0)
    "comparison_name",  # Original comparison
    "condition_name",  # Original condition
    "intervention",  # New intervention name from scores_summary
    "condition_category",  # Category (Drug, Disease, etc.)
    "log2FoldChange",  # Effect size
    "pvalue",  # Original p-value
    "fdr",  # FDR-adjusted p-value
    "neg_log10_fdr",  # -log10(FDR) for y-axis
    "final_score",  # Final score from scoring system
    "original_trust_score",  # Original trust score for reference
    "is_significant",  # Whether FDR < 0.05
    "annotation_score",  # Score used for annotation ranking
    "size",  # Point size for plotting
    "alpha",  # Point transparency for plotting
    "clock",  # Clock used
    "direction",  # Direction of effect
    "outlier_status",  # Outlier status from filtering
]

# Only include columns that exist
columns_to_save = [col for col in columns_to_save if col in combined_plot_data.columns]

# Save the complete data table
combined_plot_file = "figures/combined_volcano_plot_data.csv"
combined_plot_data[columns_to_save].to_csv(combined_plot_file, index=False)
print(
    f"\nSaved complete combined plot data ({len(combined_plot_data)} points) to: {combined_plot_file}"
)

# Also save a summary of the data
summary_stats = {
    "Total points": len(combined_plot_data),
    "Significant (FDR<0.05)": combined_plot_data["is_significant"].sum(),
    "Categories": combined_plot_data["condition_category"].nunique(),
    "Unique interventions": combined_plot_data["intervention"].nunique(),
    "Final score range": f"{combined_plot_data['final_score'].min():.1f}-{combined_plot_data['final_score'].max():.1f}",
    "Log2FC range": f"{combined_plot_data['log2FoldChange'].min():.2f} to {combined_plot_data['log2FoldChange'].max():.2f}",
}

with open("figures/combined_volcano_plot_summary.txt", "w") as f:
    f.write("Combined Volcano Plot Data Summary\n")
    f.write("=" * 50 + "\n\n")
    for key, value in summary_stats.items():
        f.write(f"{key}: {value}\n")
    f.write("\n" + "=" * 50 + "\n")
    f.write("\nCategory breakdown:\n")
    for cat, count in combined_plot_data["condition_category"].value_counts().items():
        f.write(f"  {cat}: {count} points\n")
        sig_count = combined_plot_data[
            (combined_plot_data["condition_category"] == cat)
            & (combined_plot_data["is_significant"])
        ].shape[0]
        f.write(f"    - Significant: {sig_count}\n")

print("Saved combined plot summary to: figures/combined_volcano_plot_summary.txt")

# Create detailed scores table with final ranking score and drop sparse columns
try:
    detailed_df = combined_plot_data.copy()
    # Compute final ranking score: -log10(pvalue) * abs(log2FoldChange) * final_score
    if (
        "pvalue" in detailed_df.columns
        and "log2FoldChange" in detailed_df.columns
        and "final_score" in detailed_df.columns
    ):
        # Replace zero/invalid p-values for log calculation
        safe_p = detailed_df["pvalue"].replace(0, np.nan)
        min_nonzero = safe_p[safe_p > 0].min()
        if pd.isna(min_nonzero):
            min_nonzero = 1e-300
        safe_p = safe_p.fillna(min_nonzero)
        detailed_df["final_ranking_score"] = (
            -np.log10(safe_p)
            * np.abs(detailed_df["log2FoldChange"])
            * detailed_df["final_score"].fillna(50)
        )
    else:
        detailed_df["final_ranking_score"] = np.nan

    # Drop columns with >99% missingness
    threshold = 0.99
    missing_frac = detailed_df.isna().mean()
    keep_cols = missing_frac[missing_frac <= threshold].index.tolist()
    detailed_df = detailed_df[keep_cols]

    # Remove redundant plotting-only columns if present
    redundant_plot_cols = [
        "size",
        "alpha",
        "annotation_score",
        "is_significant",
        "neg_log10_fdr",
    ]
    drop_cols = [c for c in redundant_plot_cols if c in detailed_df.columns]
    if drop_cols:
        detailed_df = detailed_df.drop(columns=drop_cols, errors="ignore")

    # Remove exact-duplicate columns (identical content)
    # Build a map of column content signature to first column name seen
    seen_signatures = {}
    duplicate_cols = []
    for col in detailed_df.columns:
        try:
            sig = tuple(detailed_df[col].astype(str).fillna("<NA>").values)
        except Exception:
            # Fallback without astype if mixed types cause issues
            sig = tuple(
                detailed_df[col]
                .apply(lambda v: "<NA>" if pd.isna(v) else str(v))
                .values
            )
        if sig in seen_signatures:
            duplicate_cols.append(col)
        else:
            seen_signatures[sig] = col
    if duplicate_cols:
        detailed_df = detailed_df.drop(columns=duplicate_cols, errors="ignore")

    # Sort by final_ranking_score desc for convenience
    if "final_ranking_score" in detailed_df.columns:
        detailed_df = detailed_df.sort_values("final_ranking_score", ascending=False)

    out_path = f"result/detailed_scores_{method}.csv"
    detailed_df.to_csv(out_path, index=False)
    print(f"Saved detailed scores table to: {out_path}")
except Exception as e:
    print(f"WARNING: Failed to save detailed scores table: {e}")

# Define significance threshold using FDR
significant = plot_df["fdr"] < 0.05

# Get unique categories and create color/marker mappings
# After lumping, we should have exactly 5 categories (top 5)
category_counts = plot_df["condition_category"].value_counts()
unique_categories = category_counts.index.tolist()

# Ensure "Other" is last if it exists
if "Other" in unique_categories:
    unique_categories.remove("Other")
    unique_categories.append("Other")

# Define thematic colors for each category
category_colors = {
    "Disease": "#DC143C",  # Crimson red for disease
    "Drug": "#4B0082",  # Indigo for pharmaceutical
    "Genetic": "#228B22",  # Forest green for genetic
    "Environment": "#FF8C00",  # Dark orange for environment
    "Other": "#808080",  # Gray for other
}

# Define distinct markers for each category
markers = ["s", "D", "p", "h", "^"]
category_markers = {
    category: markers[i % len(markers)] for i, category in enumerate(unique_categories)
}

# Create a list to track which categories have been added to legend
legend_added = set()

# Prepare data for density calculation (will plot after dots)
x_data = plot_df["log2FoldChange"].values
y_data = -np.log10(plot_df["fdr"]).values

# Remove infinite values
mask = np.isfinite(x_data) & np.isfinite(y_data)
x_data = x_data[mask]
y_data = y_data[mask]

# Plot non-significant points first (all in gray)
non_sig_mask = plot_df["fdr"] > 0.05
ax.scatter(
    plot_df.loc[non_sig_mask, "log2FoldChange"],
    -np.log10(plot_df.loc[non_sig_mask, "fdr"]),
    s=5,
    c="lightgray",
    alpha=0.3,
    zorder=1,
)

# Plot significant points by category with effect size coloring
sig_mask = plot_df["fdr"] < 0.05
if sig_mask.any():
    significant_df = plot_df[sig_mask].copy()

    # Calculate color limits for consistent coloring
    vmax = np.abs(significant_df["log2FoldChange"]).max()
    vmin = -vmax

    # Plot by category to apply different markers
    for category in unique_categories:
        category_mask = significant_df["condition_category"] == category
        if category_mask.any():
            category_data = significant_df[category_mask]

            # Layer 1: Black outline (with same marker shape)
            ax.scatter(
                category_data["log2FoldChange"],
                -np.log10(category_data["fdr"]),
                s=category_data["size"] * 1.5,
                c="black",
                alpha=1,
                marker=category_markers[category],
                zorder=2,
            )

            # Layer 2: Colored points by effect size
            scatter = ax.scatter(
                category_data["log2FoldChange"],
                -np.log10(category_data["fdr"]),
                s=category_data["size"],
                c=category_data["log2FoldChange"],
                cmap="RdBu_r",
                alpha=0.9,
                marker=category_markers[category],
                label=category if category not in legend_added else "",
                vmin=vmin,
                vmax=vmax,
                zorder=3,
            )
            legend_added.add(category)

# Add density contour plot on top of dots with more transparency
if len(x_data) > 10:  # Only plot if we have enough data
    # Create a temporary dataframe for seaborn
    import pandas as pd

    kde_df = pd.DataFrame({"log2FoldChange": x_data, "neg_log10_fdr": y_data})

    # Add contour lines with fill on top of dots - stronger colors
    contour = sns.kdeplot(
        data=kde_df,
        x="log2FoldChange",
        y="neg_log10_fdr",
        fill=True,
        cmap="rocket_r",
        alpha=0.4,  # More visible
        levels=6,
        linewidths=1.0,  # Thicker contour lines
        ax=ax,
        zorder=4,  # On top of dots
    )

    # Add contour lines only for better visibility
    sns.kdeplot(
        data=kde_df,
        x="log2FoldChange",
        y="neg_log10_fdr",
        fill=False,
        color="black",  # Dark contour lines
        alpha=0.8,
        levels=6,
        linewidths=1.0,
        ax=ax,
        zorder=5,  # On top of filled contours
    )

# Add threshold lines
plt.axhline(
    y=-np.log10(0.05), color="gray", linestyle="--", alpha=0.5, label="P = 0.05"
)
plt.axvline(x=0, color="gray", linestyle="--", alpha=0.5)
plt.axvline(
    x=np.log2(2), color="gray", linestyle="--", alpha=0.5, label="2-fold change"
)
plt.axvline(x=-np.log2(2), color="gray", linestyle="--", alpha=0.5)

# Define annotation counts per category - configurable per side.
# DNAm's scored candidate pool is much smaller and heavily biased toward
# disease pro-aging; rebalance so the anti-aging side (reprogramming, DNMT
# inhibitors, senolytics) gets more shelf space on the volcano.
if MODALITY == "dnam":
    ANNOTATION_COUNTS = {
        "Disease":     {"positive": 3,  "negative": 5},
        "Drug":        {"positive": 3,  "negative": 12},
        "Genetic":     {"positive": 3,  "negative": 12},
        "Environment": {"positive": 3,  "negative": 5},
        "Other":       {"positive": 0,  "negative": 0},
    }
else:
    ANNOTATION_COUNTS = {
        "Disease":     {"positive": 10, "negative": 3},
        "Drug":        {"positive": 10, "negative": 20},
        "Genetic":     {"positive": 8,  "negative": 8},
        "Environment": {"positive": 8,  "negative": 5},
        "Other":       {"positive": 0,  "negative": 0},
    }

# Categories are now properly configured

# Collect all annotations first
all_positive_annotations = []
all_negative_annotations = []
# Split seen-sets by direction so a positive hit from GSE_X doesn't suppress an
# anti-aging hit from the same GSE (common in DNAm where 5 clocks per
# comparison produce both signs).
seen_gse_ids_pos = set()
seen_gse_ids_neg = set()
seen_interventions_pos = set()
seen_interventions_neg = set()

# Create list to store annotation data for review
annotation_data_for_review = []

# ---- PLAUSIBILITY GATE -----------------------------------------------------
# Post-hoc filter: drop "anti-aging" annotations that are biologically
# implausible (disease states, tumors, acute infections) or reflect cohort-age
# framings rather than true interventions. Applied BEFORE top-N selection so
# the limited annotation budget is spent on interventions that make sense.
_IMPLAUSIBLE_KEYWORDS = (
    # Neoplastic / dedifferentiation
    "cancer", "carcinoma", "metastasis", "metastatic", "tumor", "tumour",
    "glioma", "glioblastoma", "melanoma", "lymphoma", "leukemia", "leukaemia",
    "sarcoma", "adenoma", "neoplas",
    # Senescence models (pro-aging)
    "senescen", "passage",
    # Age-related diseases
    "Alzheimer", "Parkinson", "COPD", "AMD", "macular degeneration",
    "cirrhosis", "fibrosis", "atherosclero", "colitis", "inflammat",
    # Acute pathologies
    "infection", "viral", "influenza", "COVID", "sepsis", "fatal",
    # Mutational / oncogenic drivers
    "mutation", "mutant", " R[0-9]", "oncogen",
    # Imprinting / uniparental disomy / Lynch-syndrome criteria
    "UPD", "imprinting", "Paternal UPD", "Maternal UPD",
    "Amsterdam", "Lynch", "Silver-Russell", "SRS",
    # Non-intervention framings
    "vs (Age)", "vs (Survival)", "(Age)", "(Survival)", "Proxy",
    "clinical criteria", "clinical SRS",
    # Disorders
    "syndrome", "resistance", "disorder",
)
_IMPLAUSIBLE_RE = "|".join(_IMPLAUSIBLE_KEYWORDS)


def _is_implausible_antiaging(df):
    """Boolean mask: anti-aging (log2FC<0) rows that match implausibility keywords.
    Checks condition_name, comparison_name, AND intervention (LLM label)."""
    anti = df["log2FoldChange"] < 0
    keyword_hit = False
    for col in ["condition_name", "comparison_name", "intervention"]:
        if col in df.columns:
            hit = df[col].astype(str).str.contains(
                _IMPLAUSIBLE_RE, case=False, na=False, regex=True
            )
            keyword_hit = hit if keyword_hit is False else (keyword_hit | hit)
    return anti & keyword_hit
# ----------------------------------------------------------------------------

for category in unique_categories:
    if category not in plot_df["condition_category"].values:
        continue

    # Skip if no annotations configured for this category
    if category not in ANNOTATION_COUNTS or (
        ANNOTATION_COUNTS[category]["positive"] == 0
        and ANNOTATION_COUNTS[category]["negative"] == 0
    ):
        continue

    # Get data for this category
    # Filter by: FDR < 0.05, final_score floor (modality-aware), AND plausibility.
    # DNAm scores skew higher (median 84.5) than RNA (60), so a >=100 bar there
    # starves the anti-aging side of candidates. Use a lower DNAm floor.
    _score_floor = 80 if MODALITY == "dnam" else 100
    _cat_mask = (
        (plot_df["condition_category"] == category)
        & (plot_df["fdr"] < 0.05)
        & (plot_df["final_score"] >= _score_floor)
        & ~_is_implausible_antiaging(plot_df)
    )
    category_df = plot_df[_cat_mask].copy()

    if len(category_df) == 0:
        continue

    # Calculate score for ranking (significance * effect size)
    category_df["score"] = -np.log10(category_df["fdr"]) * np.abs(
        category_df["log2FoldChange"]
    )

    # Deduplicate by intervention name first, then by GSE ID
    # Keep the entry with the highest score for each intervention
    category_df = category_df.sort_values("score", ascending=False)
    category_df = category_df.drop_duplicates(subset=["intervention"], keep="first")

    # Separate positive and negative effects
    positive_df = category_df[category_df["log2FoldChange"] > 0].sort_values(
        "score", ascending=False
    )
    negative_df = category_df[category_df["log2FoldChange"] < 0].sort_values(
        "score", ascending=False
    )

    # Get configured number of annotations
    n_pos = ANNOTATION_COUNTS[category]["positive"]
    n_neg = ANNOTATION_COUNTS[category]["negative"]

    # Save twice as many to table for LLM review flexibility
    n_pos_save = n_pos * 2 if n_pos > 0 else 0
    n_neg_save = n_neg * 2 if n_neg > 0 else 0

    # Select top results for saving (2x for review)
    top_pos_save = positive_df.head(n_pos_save) if n_pos_save > 0 else pd.DataFrame()
    top_neg_save = negative_df.head(n_neg_save) if n_neg_save > 0 else pd.DataFrame()

    # Select top results for plotting (configured amount)
    top_pos = positive_df.head(n_pos) if n_pos > 0 else pd.DataFrame()
    top_neg = negative_df.head(n_neg) if n_neg > 0 else pd.DataFrame()

    # First, save ALL candidates to review table (2x amount)
    for _, row in top_pos_save.iterrows():
        annotation_data_for_review.append(
            {
                "gse_id": row["gse_id"],
                "category": category,
                "condition_name": row["condition_name"],
                "intervention": row.get("intervention", "Unknown"),
                "current_annotation": row.get("intervention", "Unknown")[:25],
                "log2FoldChange": row["log2FoldChange"],
                "fdr": row["fdr"],
                "pvalue": row["pvalue"],
                "trust_score": row.get(
                    "original_trust_score", row.get("trust_score", 0)
                ),
                "final_score": row.get("final_score", np.nan),
                "direction": "positive",
                "comparison": row.get("comparison", ""),
                "sample_info": row.get("sample_info", ""),
                "condition_category": row["condition_category"],
            }
        )

    for _, row in top_neg_save.iterrows():
        annotation_data_for_review.append(
            {
                "gse_id": row["gse_id"],
                "category": category,
                "condition_name": row["condition_name"],
                "intervention": row.get("intervention", "Unknown"),
                "current_annotation": row.get("intervention", "Unknown")[:25],
                "log2FoldChange": row["log2FoldChange"],
                "fdr": row["fdr"],
                "pvalue": row["pvalue"],
                "trust_score": row.get(
                    "original_trust_score", row.get("trust_score", 0)
                ),
                "final_score": row.get("final_score", np.nan),
                "direction": "negative",
                "comparison": row.get("comparison", ""),
                "sample_info": row.get("sample_info", ""),
                "condition_category": row["condition_category"],
            }
        )

    # Then add to annotation lists for PLOTTING (skip duplicates by intervention name)
    for _, row in top_pos.iterrows():
        intervention = row.get("intervention", "Unknown")
        if intervention not in seen_interventions_pos and row["gse_id"] not in seen_gse_ids_pos:
            all_positive_annotations.append(
                {
                    "x": row["log2FoldChange"],
                    "y": -np.log10(row["fdr"]),
                    "label": intervention[:25],
                    "category": category,
                    "color": category_colors[category],
                    "gse_id": row["gse_id"],
                    "intervention": intervention,
                }
            )
            seen_interventions_pos.add(intervention)
            seen_gse_ids_pos.add(row["gse_id"])

    for _, row in top_neg.iterrows():
        intervention = row.get("intervention", "Unknown")
        if intervention not in seen_interventions_neg and row["gse_id"] not in seen_gse_ids_neg:
            all_negative_annotations.append(
                {
                    "x": row["log2FoldChange"],
                    "y": -np.log10(row["fdr"]),
                    "label": intervention[:25],
                    "category": category,
                    "color": category_colors[category],
                    "gse_id": row["gse_id"],
                    "intervention": intervention,
                }
            )
            seen_interventions_neg.add(intervention)
            seen_gse_ids_neg.add(row["gse_id"])

# Save annotation data for review
if annotation_data_for_review:
    annotation_df = pd.DataFrame(annotation_data_for_review)
    annotation_df = annotation_df.sort_values(
        ["category", "direction", "log2FoldChange"], ascending=[True, True, False]
    )

    # Always save the base annotation table
    base_annotation_file = f"figures/volcano_plot_annotations_{method}.csv"
    annotation_df.to_csv(base_annotation_file, index=False)
    print(f"\nSaved {len(annotation_df)} annotations to {base_annotation_file}")

    # Filter to only include annotations that are actually plotted
    # (union of positive and negative seen-sets after the direction-aware split)
    _plotted_gse = seen_gse_ids_pos | seen_gse_ids_neg
    plotted_annotation_df = annotation_df[
        annotation_df["gse_id"].isin(_plotted_gse)
    ].copy()

    # Use current annotations directly; skip refined table
    improved_annotations = {}
    for _, row in plotted_annotation_df.iterrows():
        gse_id = row["gse_id"]
        improved_annotations[gse_id] = row["current_annotation"]

    # Add final annotations to dataframe for reference
    plotted_annotation_df["final_annotation"] = plotted_annotation_df["gse_id"].map(
        improved_annotations
    )
    plotted_annotation_df.to_csv(
        f"figures/volcano_plot_annotations_{method}_with_final.csv", index=False
    )

# ============================================================================
# SAVE ALL ANNOTATIONS FOR MANUAL ADJUSTMENT BEFORE PLOTTING
# ============================================================================
print("\n" + "=" * 80)
print("ANNOTATION MANAGEMENT - CHECKING FOR MANUAL ADJUSTMENTS")
print("=" * 80)

manual_adjustment_file = (
    f"figures/volcano_plot_annotations_FOR_MANUAL_ADJUSTMENT_{method}.csv"
)

# Check if force refresh is enabled
if FORCE_REFRESH_ANNOTATIONS and os.path.exists(manual_adjustment_file):
    print(f"\n⚠️  FORCE_REFRESH_ANNOTATIONS = True")
    print(f"   Deleting existing manual adjustment file to start fresh...")
    os.remove(manual_adjustment_file)
    print(f"   ✓ Deleted: {manual_adjustment_file}")

# First, check if manual adjustment file already exists
if os.path.exists(manual_adjustment_file):
    print(f"\n✓ Found existing manual adjustment file!")
    print(f"  Loading: {manual_adjustment_file}")

    try:
        existing_manual_df = pd.read_csv(manual_adjustment_file)
        print(f"  Existing file has {len(existing_manual_df)} entries")

        # Get currently selected annotations (from the collection process)
        current_gse_ids = set()
        for ann in all_positive_annotations + all_negative_annotations:
            current_gse_ids.add(ann["gse_id"])

        # Get ALL existing GSE IDs from manual file (both included and excluded)
        # This ensures we don't suggest annotations that were already tried and excluded
        all_existing_gse_ids = set(existing_manual_df["gse_id"])

        # Get existing GSE IDs from manual file (that are not excluded)
        existing_included = existing_manual_df[
            existing_manual_df.get("plot_annotation", "yes").str.lower() != "no"
        ]
        existing_gse_ids = set(existing_included["gse_id"])

        # Count excluded entries
        excluded_count = len(existing_manual_df) - len(existing_included)

        print(f"  Manual file status:")
        print(f"    - {len(existing_included)} entries marked to plot")
        print(f"    - {excluded_count} entries excluded (plot_annotation = 'no')")
        print(
            f"    - {len(all_existing_gse_ids)} total unique GSE IDs in file "
            f"(won't be re-suggested)"
        )

        # If some entries were excluded, we need to add replacement entries
        if excluded_count > 0:
            print(f"\n  Need to add {excluded_count} replacement annotations...")

            # Get additional candidates from category_df
            replacement_annotations = []
            for category in unique_categories:
                if category not in plot_df["condition_category"].values:
                    continue

                # Get configured annotation counts
                if category not in ANNOTATION_COUNTS:
                    continue

                n_pos = ANNOTATION_COUNTS[category]["positive"]
                n_neg = ANNOTATION_COUNTS[category]["negative"]

                # Get candidate data
                category_df_for_replacements = plot_df[
                    (plot_df["condition_category"] == category)
                    & (plot_df["fdr"] < 0.05)
                    & (plot_df["final_score"] >= 100)
                ].copy()

                if len(category_df_for_replacements) == 0:
                    continue

                # Calculate ranking score
                category_df_for_replacements["score"] = -np.log10(
                    category_df_for_replacements["fdr"]
                ) * np.abs(category_df_for_replacements["log2FoldChange"])

                # Deduplicate by intervention
                category_df_for_replacements = category_df_for_replacements.sort_values(
                    "score", ascending=False
                )
                category_df_for_replacements = (
                    category_df_for_replacements.drop_duplicates(
                        subset=["intervention"], keep="first"
                    )
                )

                # Get positive and negative separately
                pos_candidates = category_df_for_replacements[
                    category_df_for_replacements["log2FoldChange"] > 0
                ].head(
                    n_pos * 3
                )  # Get 3x more candidates

                neg_candidates = category_df_for_replacements[
                    category_df_for_replacements["log2FoldChange"] < 0
                ].head(n_neg * 3)

                # Add to replacement pool
                # Check against all_existing_gse_ids to avoid suggesting previously excluded annotations
                for _, row in pos_candidates.iterrows():
                    if row["gse_id"] not in all_existing_gse_ids:
                        replacement_annotations.append(
                            {
                                "gse_id": row["gse_id"],
                                "category": category,
                                "direction": "positive",
                                "row": row,
                            }
                        )

                for _, row in neg_candidates.iterrows():
                    if row["gse_id"] not in all_existing_gse_ids:
                        replacement_annotations.append(
                            {
                                "gse_id": row["gse_id"],
                                "category": category,
                                "direction": "negative",
                                "row": row,
                            }
                        )

            # Sort replacement annotations by score and take top ones needed
            print(
                f"  Found {len(replacement_annotations)} candidate replacements "
                f"(excluding already-tried annotations)"
            )
            replacement_annotations = sorted(
                replacement_annotations, key=lambda x: x["row"]["score"], reverse=True
            )[:excluded_count]

            # Add new entries to the existing file
            new_entries = []
            for rep in replacement_annotations:
                row = rep["row"]
                new_entries.append(
                    {
                        "gse_id": row["gse_id"],
                        "category": rep["category"],
                        "direction": rep["direction"],
                        "x_position": row["log2FoldChange"],
                        "y_position": -np.log10(row["fdr"]),
                        "current_label": row.get("intervention", "Unknown")[:25],
                        "full_intervention": row.get("intervention", "Unknown"),
                        "color": category_colors[rep["category"]],
                        "adjusted_label": row.get("intervention", "Unknown")[:25],
                        "plot_annotation": "yes",
                        "refilled": "YES - NEW REPLACEMENT",  # Mark as refilled
                    }
                )

            if new_entries:
                new_df = pd.DataFrame(new_entries)
                # Ensure existing_manual_df has 'refilled' column
                if "refilled" not in existing_manual_df.columns:
                    existing_manual_df["refilled"] = ""
                updated_df = pd.concat([existing_manual_df, new_df], ignore_index=True)
                updated_df = updated_df.sort_values(
                    ["category", "direction", "y_position"],
                    ascending=[True, True, False],
                )
                updated_df.to_csv(manual_adjustment_file, index=False)
                print(f"  ✓ Added {len(new_entries)} replacement annotations to file")
                print(
                    f"     (marked with 'refilled' = 'YES - NEW REPLACEMENT' for easy identification)"
                )
            else:
                print(f"  No suitable replacement annotations found")
                updated_df = existing_manual_df
        else:
            updated_df = existing_manual_df
            print(f"  No replacements needed - using existing annotations as-is")

        # Now apply the manual adjustments to the annotation lists
        print(f"\n  Applying manual adjustments to plot...")
        manual_adjustments = {}
        exclude_gse_ids = set()

        for _, row in updated_df.iterrows():
            gse_id = row["gse_id"]
            if row.get("plot_annotation", "yes").lower() == "no":
                exclude_gse_ids.add(gse_id)
            else:
                manual_adjustments[gse_id] = {
                    "label": row.get("adjusted_label", row["current_label"]),
                    "category": row.get(
                        "category", "Other"
                    ),  # Respect category changes
                }

        # Filter annotation lists to only include entries from manual file
        all_positive_annotations = [
            ann
            for ann in all_positive_annotations
            if ann["gse_id"] in manual_adjustments
        ]
        all_negative_annotations = [
            ann
            for ann in all_negative_annotations
            if ann["gse_id"] in manual_adjustments
        ]

        # Update labels and categories with values from CSV file
        # This ensures your manual edits in 'adjusted_label' and 'category' columns are used
        adjusted_count = 0
        category_changed_count = 0

        for ann in all_positive_annotations:
            if ann["gse_id"] in manual_adjustments:
                adj = manual_adjustments[ann["gse_id"]]
                ann["label"] = adj["label"]

                # Check if category was changed
                if ann["category"] != adj["category"]:
                    ann["category"] = adj["category"]
                    # Update color to match new category
                    if adj["category"] in category_colors:
                        ann["color"] = category_colors[adj["category"]]
                    category_changed_count += 1

                adjusted_count += 1

        for ann in all_negative_annotations:
            if ann["gse_id"] in manual_adjustments:
                adj = manual_adjustments[ann["gse_id"]]
                ann["label"] = adj["label"]

                # Check if category was changed
                if ann["category"] != adj["category"]:
                    ann["category"] = adj["category"]
                    # Update color to match new category
                    if adj["category"] in category_colors:
                        ann["color"] = category_colors[adj["category"]]
                    category_changed_count += 1

                adjusted_count += 1

        print(f"  ✓ Using {len(manual_adjustments)} annotations for plotting")
        print(f"  ✓ Applied {adjusted_count} adjusted labels from CSV file")
        if category_changed_count > 0:
            print(
                f"  ✓ Applied {category_changed_count} category changes from CSV file"
            )
        print(f"  ✓ Excluded {len(exclude_gse_ids)} annotations")

    except Exception as e:
        print(f"\n⚠ Error processing manual adjustments: {e}")
        print("  Creating fresh annotation file...")
        # Fall through to create new file below
        existing_manual_df = None

else:
    print(f"\n📝 No existing manual adjustment file found.")
    print(f"   Creating new file: {manual_adjustment_file}")
    existing_manual_df = None

# If no existing file, create a new one from current annotations
if not os.path.exists(manual_adjustment_file) or existing_manual_df is None:
    all_annotations_to_plot = []

    for ann in all_positive_annotations:
        all_annotations_to_plot.append(
            {
                "gse_id": ann["gse_id"],
                "category": ann["category"],
                "direction": "positive",
                "x_position": ann["x"],
                "y_position": ann["y"],
                "current_label": ann["label"],
                "full_intervention": ann["intervention"],
                "color": ann["color"],
                "adjusted_label": ann["label"],
                "plot_annotation": "yes",
                "refilled": "",  # Empty for original annotations
            }
        )

    for ann in all_negative_annotations:
        all_annotations_to_plot.append(
            {
                "gse_id": ann["gse_id"],
                "category": ann["category"],
                "direction": "negative",
                "x_position": ann["x"],
                "y_position": ann["y"],
                "current_label": ann["label"],
                "full_intervention": ann["intervention"],
                "color": ann["color"],
                "adjusted_label": ann["label"],
                "plot_annotation": "yes",
                "refilled": "",  # Empty for original annotations
            }
        )

    if all_annotations_to_plot:
        annotations_df = pd.DataFrame(all_annotations_to_plot)
        annotations_df = annotations_df.sort_values(
            ["category", "direction", "y_position"], ascending=[True, True, False]
        )
        annotations_df.to_csv(manual_adjustment_file, index=False)

        print(f"\n   ✓ Created file with {len(annotations_df)} annotations")
        print(f"\n   Annotation counts by category and direction:")
        for category in annotations_df["category"].unique():
            cat_df = annotations_df[annotations_df["category"] == category]
            pos_count = len(cat_df[cat_df["direction"] == "positive"])
            neg_count = len(cat_df[cat_df["direction"] == "negative"])
            print(
                f"     {category:15s}: {pos_count:2d} positive, {neg_count:2d} negative"
            )

print("\n" + "=" * 80)
print("MANUAL ADJUSTMENT INSTRUCTIONS:")
print("=" * 80)
print("\nFile Columns:")
print("  • gse_id: Study identifier")
print(
    "  • category: ✏️  EDIT to change intervention category (Drug, Disease, Genetic, Environment, Other)"
)
print("  • direction: positive or negative effect")
print("  • current_label: Original auto-generated label")
print("  • adjusted_label: ✏️  EDIT THIS to change annotation text shown on plot")
print("  • plot_annotation: ✏️  Set to 'no' to exclude from plot")
print("  • refilled: Shows 'YES - NEW REPLACEMENT' for auto-filled replacements")
print("\nWorkflow:")
print("  1. Edit 'adjusted_label' column to customize annotation text")
print("  2. Edit 'category' column to reclassify interventions (affects color on plot)")
print("  3. Set 'plot_annotation' to 'no' to exclude specific annotations")
print("  4. Save the file and re-run this script")
print("  5. Excluded annotations will be automatically replaced with new ones")
print("  6. Review the new replacements (marked in 'refilled' column) and repeat")
print(
    f"\nTo start completely fresh, set FORCE_REFRESH_ANNOTATIONS = True (currently: {FORCE_REFRESH_ANNOTATIONS})"
)
print("=" * 80 + "\n")

# Manual positioning similar to category-specific plots
texts = []

# Get current axis limits for positioning
xlim = ax.get_xlim()
ylim = ax.get_ylim()

# Extend x-axis limits to make room for annotations (reduced extension)
x_range = xlim[1] - xlim[0]
new_xlim = (
    xlim[0] - x_range * 0.15,
    xlim[1] + x_range * 0.15,
)  # Reduced from 0.3 to 0.15
ax.set_xlim(new_xlim)

# Position positive annotations on the right
if all_positive_annotations:
    # Sort by y-coordinate
    all_positive_annotations = sorted(all_positive_annotations, key=lambda x: x["y"])

    # Find the maximum log2 effect size in SIGNIFICANT data only
    significant_data = plot_df[plot_df["fdr"] < 0.05]
    if len(significant_data) > 0:
        max_log2_effect = significant_data["log2FoldChange"].max()
    else:
        max_log2_effect = 0

    # Position text closer to plot
    x_start = (
        max_log2_effect + x_range * 0.01
    )  # Start horizontal line just after max significant data point
    x_end = new_xlim[1] - x_range * 0.02  # Reduced from 0.05 to 0.02

    # Spread annotations more vertically
    y_min = min([a["y"] for a in all_positive_annotations])
    y_max = max([a["y"] for a in all_positive_annotations])
    y_range = y_max - y_min
    y_positions = np.linspace(
        y_min - y_range * 0.2, y_max + y_range * 0.2, len(all_positive_annotations)
    )

    for ann, y_pos in zip(all_positive_annotations, y_positions):
        # Use the label from annotation (which has been updated from adjusted_label CSV column)
        label_to_use = ann["label"]

        # Skip if label is empty
        if not label_to_use or not label_to_use.strip():
            continue

        text = plt.text(
            x_end,
            y_pos,
            label_to_use,
            ha="right",
            va="center",  # Right align for right side
            fontsize=9,
            color=ann["color"],
            path_effects=[pe.withStroke(linewidth=3, foreground="white")],
            weight="bold",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.8
            ),
        )
        texts.append(text)

        # Draw connection line: horizontal-slanted-horizontal (3 segments)
        # Segment 1: Short horizontal line from dot
        x_h1_end = ann["x"] + x_range * 0.01  # Very short horizontal

        # Segment 2: Slanted line (keep it short in x-direction)
        x_slant_end = x_h1_end + x_range * 0.025  # Short slanted segment

        # Segment 3: Horizontal line to text
        x_text_actual = x_end - x_range * 0.005  # Leave small gap before text

        plt.plot(
            [ann["x"], x_h1_end, x_slant_end, x_text_actual],
            [ann["y"], ann["y"], y_pos, y_pos],
            color="black",
            ls=":",
            alpha=0.4,
            lw=1.0,
        )

# Position negative annotations on the left
if all_negative_annotations:
    # Sort by y-coordinate
    all_negative_annotations = sorted(all_negative_annotations, key=lambda x: x["y"])

    # Find the minimum log2 effect size in SIGNIFICANT data only
    significant_data = plot_df[plot_df["fdr"] < 0.05]
    if len(significant_data) > 0:
        min_log2_effect = significant_data["log2FoldChange"].min()
    else:
        min_log2_effect = 0

    # Position text closer to plot
    x_start = (
        min_log2_effect - x_range * 0.01
    )  # Start horizontal line just after min significant data point
    x_end = new_xlim[0] + x_range * 0.02  # Reduced from 0.05 to 0.02

    # Spread annotations more vertically
    y_min = min([a["y"] for a in all_negative_annotations])
    y_max = max([a["y"] for a in all_negative_annotations])
    y_range = y_max - y_min
    y_positions = np.linspace(
        y_min - y_range * 0.2, y_max + y_range * 0.1, len(all_negative_annotations)
    )

    for ann, y_pos in zip(all_negative_annotations, y_positions):
        # Use the label from annotation (which has been updated from adjusted_label CSV column)
        label_to_use = ann["label"]

        # Skip if label is empty
        if not label_to_use or not label_to_use.strip():
            continue

        text = plt.text(
            x_end,
            y_pos,
            label_to_use,
            ha="left",
            va="center",  # Left align for left side
            fontsize=9,
            color=ann["color"],
            path_effects=[pe.withStroke(linewidth=3, foreground="white")],
            weight="bold",
            bbox=dict(
                boxstyle="round,pad=0.3", facecolor="white", edgecolor="none", alpha=0.8
            ),
        )
        texts.append(text)

        # Draw connection line: horizontal-slanted-horizontal (3 segments)
        # Segment 1: Short horizontal line from dot
        x_h1_end = ann["x"] - x_range * 0.01  # Very short horizontal (leftward)

        # Segment 2: Slanted line (keep it short in x-direction)
        x_slant_end = x_h1_end - x_range * 0.025  # Short slanted segment

        # Segment 3: Horizontal line to text
        x_text_actual = x_end + x_range * 0.005  # Leave small gap before text

        plt.plot(
            [ann["x"], x_h1_end, x_slant_end, x_text_actual],
            [ann["y"], ann["y"], y_pos, y_pos],
            color="black",
            ls=":",
            alpha=0.4,
            lw=1.0,
        )

# Add labels and title
plt.xlabel("Effect Size")
plt.ylabel("-log10(Adjusted P-value)")

# Get statistics for title
significant_count = sum(plot_df["fdr"] < 0.05)
total_count = len(plot_df)

# Create comprehensive title with stats - smaller font
_MODALITY_LABEL = {"rna": "RNA-seq", "dnam": "DNA methylation"}[MODALITY]
plt.title(
    f"Human {_MODALITY_LABEL} Aging Interventions\n"
    + f"Total comparisons: {total_count:,} | Significant (FDR<0.05): {significant_count:,} ({significant_count/total_count*100:.1f}%)",
    fontsize=12,
)

# Add colorbar for effect size
if sig_mask.any() and "scatter" in locals():
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label("Effect Size", rotation=270, labelpad=20)

# Create two separate legends - one for categories, one for sizes
from matplotlib.lines import Line2D

# Create category legend elements with markers and colors
category_handles = []
for category in unique_categories:
    if category in category_markers and category in category_colors:
        category_handles.append(
            Line2D(
                [0],
                [0],
                marker=category_markers[category],
                color="w",
                markerfacecolor=category_colors[category],
                markersize=10,
                label=category,
                linestyle="None",
            )
        )

# Create size legend elements for final scores
min_score = plot_df["final_score"].min()
max_score = plot_df["final_score"].max()

# Always use the 50-130 scale for final scores
score_values = [50, 70, 90, 110, 130]

size_handles = []
if pd.notna(min_score) and pd.notna(max_score):
    if max_score == min_score:
        # Single-value range: show one representative size
        representative_scores = [min_score]
        denom = 1.0
    else:
        representative_scores = [
            s for s in score_values if s >= min_score and s <= max_score
        ]
        denom = max_score - min_score

    for score in representative_scores:
        normalized = (score - min_score) / denom if denom != 0 else 1.0
        size = normalized**3 * 120
        marker_size = np.sqrt(size / np.pi) * 2
        size_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="gray",
                markersize=marker_size,
                alpha=1,
                linestyle="None",
                label=f"{int(score) if not np.isnan(score) else ''}",
            )
        )

# HUMAN: skip the category and size legends — colorbar alone is sufficient.
# (Mouse script still draws both; we just suppress them here.)

# Stats have been moved to the title

# Apply theme and save
theme.apply_transforms()
plt.tight_layout()

# Save in both PDF and PNG formats
plt.savefig(
    f"figures/volcano_plot_{method}_all_categories_new.pdf",
    bbox_inches="tight",
)
plt.savefig(
    f"figures/volcano_plot_{method}_all_categories_new.png",
    dpi=600,
    bbox_inches="tight",
)
plt.close()

print("Combined plot with new coloring approach saved.")
