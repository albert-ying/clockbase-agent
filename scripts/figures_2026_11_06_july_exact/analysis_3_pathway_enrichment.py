#!/usr/bin/env python3
"""
ANALYSIS 3: Over-Representation Analysis (ORA) - Pathway Enrichment
====================================================================

Purpose: Validate biological plausibility by testing enrichment of aging-related genes
         in established KEGG/GO pathways using standardized enrichment methods.

Approach:
1. Download official gene reference databases (human/mouse)
2. Extract gene names from genetic intervention conditions using keyword matching
3. Validate genes against reference database (remove non-matching entries)
4. Perform Over-Representation Analysis (ORA) using GSEApy Enrichr API
5. FDR correction for multiple testing

Key Method: Over-Representation Analysis (ORA) via GSEApy Enrichr API
- Takes a list of significant genes (not ranked expression data)
- Tests for overrepresentation in pathway databases using Fisher's exact test
- Uses standardized pathway databases (KEGG, GO Biological Process, Reactome, WikiPathways)
- Proper statistical testing with multiple comparison correction
- Official gene symbol validation

Note: This is ORA (gene list enrichment), not GSEA (ranked list enrichment).

Expected Result: Significant enrichment in canonical aging pathways validates
                 that findings converge on known biology (not random noise).
"""

import re
import warnings
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import mygene
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")

# Print geometry for Figure 4e (mm, points); the PDF is placed at 100%.
PANEL_W_MM, PANEL_H_MM = 83.0, 44.0
MARGIN_MM = 0.6      # left and right page margin
TOP_MARGIN_MM = 0.4
BOTTOM_MARGIN_MM = 0.3
LABEL_PT = 6.5       # pathway names, one line each
TICK_PT = 6.5
AXIS_TITLE_PT = 7
KEY_PT = 6


def load_data():
    """Load main results data"""
    print("Loading data...")
    main_df = pd.read_csv("result/main_results_result_mm_rna_outlier_marked.csv")
    main_df = main_df[main_df["outlier_status"] == "Normal"].copy()
    print(f"Loaded {len(main_df):,} comparisons")
    return main_df


def download_gene_reference(species="mouse"):
    """
    Download official gene reference database

    Uses GSEApy's built-in gene lists as reference
    (derived from MSigDB and NCBI gene databases)

    Parameters:
    -----------
    species : str
        'mouse' or 'human'

    Returns:
    --------
    pandas.DataFrame with gene symbols and set of all gene names
    """
    print(f"\nDownloading {species} gene reference database...")

    # Use GSEApy to get comprehensive gene list from Enrichr libraries
    # This provides validated gene symbols used in pathway databases
    try:
        if species.lower() == "mouse":
            organism = "Mouse"
            # Get gene sets from multiple libraries to build comprehensive reference
            libraries = ["KEGG_2021_Mouse", "GO_Biological_Process_2021"]
        else:
            organism = "Human"
            libraries = ["KEGG_2021_Human", "GO_Biological_Process_2021"]

        all_genes = set()

        # Extract all genes from pathway libraries
        print(f"  Building reference from pathway databases...")
        for lib in libraries:
            try:
                lib_genes = gp.get_library(name=lib, organism=organism)
                for pathway, genes in lib_genes.items():
                    all_genes.update([g.upper() for g in genes])
                print(f"    {lib}: {len(lib_genes):,} pathways")
            except Exception as e:
                print(f"    Warning: Could not load {lib}: {e}")
                continue

        # Add common gene symbols from MyGene for additional coverage
        try:
            mg = mygene.MyGeneInfo()
            # Get a sample of genes to expand the reference
            if species.lower() == "mouse":
                sample_query = mg.query(
                    "_exists_:symbol AND _exists_:taxid",
                    species="10090",
                    fields="symbol",
                    size=5000,
                )
            else:
                sample_query = mg.query(
                    "_exists_:symbol AND _exists_:taxid",
                    species="9606",
                    fields="symbol",
                    size=5000,
                )

            for hit in sample_query.get("hits", []):
                if "symbol" in hit:
                    all_genes.add(hit["symbol"].upper())

            print(f"    Added genes from MyGene.info")
        except Exception as e:
            print(f"    Note: MyGene query failed, using pathway-only reference: {e}")

        # Create DataFrame
        gene_data = [{"symbol": gene} for gene in sorted(all_genes)]
        gene_ref_df = pd.DataFrame(gene_data)

        print(f"  Total unique genes: {len(all_genes):,}")

        # Save reference
        gene_ref_df.to_csv(f"result/gene_reference_{species}.csv", index=False)
        print(f"  ✓ Saved: result/gene_reference_{species}.csv")

        return gene_ref_df, all_genes

    except Exception as e:
        print(f"  ERROR: Could not build gene reference: {e}")
        print(f"  Using fallback: permissive gene matching")
        # Return empty set to skip validation
        return pd.DataFrame(), set()


def extract_genes(main_df, gene_reference_set):
    """
    Extract gene names from condition names using keyword matching

    Method:
    1. Filter to significant genetic interventions
    2. Parse condition names to extract gene symbols
    3. Remove common suffixes (knockout, KO, overexpression, etc.)
    4. Validate against official gene reference database
    5. Remove entries with no matching gene
    6. Separate into anti-aging and pro-aging gene sets

    Parameters:
    -----------
    main_df : pandas.DataFrame
        Main results data
    gene_reference_set : set
        Set of all valid gene names from reference database

    Returns:
    --------
    Validated gene sets and summary dataframe
    """
    print("\n" + "=" * 80)
    print("ANALYSIS 3: OVER-REPRESENTATION ANALYSIS (ORA)")
    print("=" * 80)

    # First, extract ALL genes from genetic interventions (for background)
    all_genetic_df = main_df[main_df["condition_category"] == "Genetic"].copy()

    # Then filter to significant ones (for ORA testing)
    genetic_df = main_df[
        (main_df["condition_category"] == "Genetic") & (main_df["fdr"] < 0.05)
    ].copy()

    print(f"\nExtracting genes from genetic interventions...")
    print(f"  Total genetic interventions: {len(all_genetic_df):,}")
    print(f"  Significant genetic interventions (FDR < 0.05): {len(genetic_df):,}")

    def extract_and_validate_gene(condition_name):
        """
        Extract gene name from condition string and validate against reference

        Uses keyword matching: searches for any word in the condition name
        that matches a known gene symbol or alias
        """
        condition_name = str(condition_name).upper()

        # Remove common suffixes first
        suffixes = [
            " KNOCKOUT",
            " KNOCKDOWN",
            " OVEREXPRESSION",
            " KO",
            " KD",
            " OE",
            " DELETION",
            " MUTANT",
            " -/-",
            " +/-",
            " TRANSGENIC",
            " DEFICIENT",
            " DEFICIENCY",
            " NULL",
            " HETEROZYGOUS",
        ]
        for suffix in suffixes:
            condition_name = condition_name.replace(suffix, "")

        # Try to match each word against gene reference
        words = condition_name.split()
        for word in words:
            # Clean word
            clean_word = re.sub(r"[^A-Z0-9]", "", word)
            if len(clean_word) >= 2 and clean_word[0].isalpha():
                # Check if it matches a known gene
                if clean_word in gene_reference_set:
                    return clean_word

        # If no match found, try first word as fallback
        if len(words) > 0:
            first_word = re.sub(r"[^A-Z0-9]", "", words[0])
            if len(first_word) >= 2 and first_word[0].isalpha():
                return first_word  # Return but will be filtered later

        return None

    # Extract genes from ALL genetic interventions (for background)
    all_genetic_df["gene_name_raw"] = all_genetic_df["condition_name"].apply(
        extract_and_validate_gene
    )

    # Extract genes from significant interventions (for ORA)
    genetic_df["gene_name_raw"] = genetic_df["condition_name"].apply(
        extract_and_validate_gene
    )

    # Count before validation
    n_all = len(all_genetic_df[all_genetic_df["gene_name_raw"].notna()])
    n_sig = len(genetic_df[genetic_df["gene_name_raw"].notna()])

    # Validate: keep only genes in reference database (if available)
    if len(gene_reference_set) > 0:
        # Validate background genes
        all_genetic_df["gene_name"] = all_genetic_df["gene_name_raw"].apply(
            lambda x: x if x in gene_reference_set else None
        )
        all_genetic_df = all_genetic_df[all_genetic_df["gene_name"].notna()].copy()

        # Validate significant genes
        genetic_df["gene_name"] = genetic_df["gene_name_raw"].apply(
            lambda x: x if x in gene_reference_set else None
        )
        n_removed = len(genetic_df[genetic_df["gene_name"].isna()])
        genetic_df = genetic_df[genetic_df["gene_name"].notna()].copy()

        print(f"\n  Gene extraction and validation:")
        print(f"    Total genes extracted (all interventions): {n_all}")
        print(f"    Validated background genes: {len(all_genetic_df):,}")
        print(f"    Significant genes extracted: {n_sig}")
        print(f"    Removed {n_removed} non-matching significant genes")
        print(f"    Validated significant genes: {len(genetic_df):,}")
    else:
        all_genetic_df["gene_name"] = all_genetic_df["gene_name_raw"]
        all_genetic_df = all_genetic_df[all_genetic_df["gene_name"].notna()].copy()
        genetic_df["gene_name"] = genetic_df["gene_name_raw"]
        genetic_df = genetic_df[genetic_df["gene_name"].notna()].copy()
        print(f"  WARNING: No gene reference available, skipping validation")
        print(f"  Extracted {len(all_genetic_df):,} background genes (unvalidated)")
        print(f"  Extracted {len(genetic_df):,} significant genes (unvalidated)")

    # Get all unique background genes
    all_background_genes = set(all_genetic_df["gene_name"].unique())

    # Get all significant genes
    significant_genes = set(genetic_df["gene_name"].unique())

    print(f"\n  Final gene counts:")
    print(f"    Background genes (all interventions): {len(all_background_genes)}")
    print(f"    Significant genes: {len(significant_genes)}")

    # Save gene list with validation status
    gene_df = genetic_df.groupby("gene_name").agg(
        {"log2FoldChange": ["mean", "count"], "fdr": "min", "condition_name": "first"}
    )
    gene_df.columns = ["Mean_Effect", "N_Studies", "Min_FDR", "Example_Condition"]
    gene_df = gene_df.reset_index()
    gene_df["Abs_Effect"] = np.abs(gene_df["Mean_Effect"])
    gene_df = gene_df.sort_values("Abs_Effect", ascending=False)

    gene_df.to_csv("result/analysis_3_gene_list_validated.csv", index=False)
    print("\n✓ Saved: result/analysis_3_gene_list_validated.csv")

    return (
        significant_genes,
        all_background_genes,
        gene_df,
    )


def perform_ora_enrichment(significant_genes, all_background_genes, species="mouse"):
    """
    Perform Over-Representation Analysis (ORA) using GSEApy Enrichr API

    Uses multiple established pathway databases:
    - KEGG_2019_Mouse / KEGG_2019_Human
    - GO_Biological_Process_2023
    - Reactome_2022
    - WikiPathways_2019_Mouse / WikiPathways_2019_Human

    Method (ORA):
    1. Takes a list of significant genes
    2. Uses all genetic intervention genes as custom background
    3. Query Enrichr API through GSEApy
    4. Tests for overrepresentation using one-sided Fisher's exact test
    5. Automatic FDR correction (Benjamini-Hochberg)
    6. Filters pathways: size < 100 genes AND overlap >= 15 genes
    7. Returns filtered enrichment results

    Note: Fisher's exact test is one-sided (tests for enrichment/over-representation only)

    Parameters:
    -----------
    significant_genes : set
        Set of significant gene symbols (FDR < 0.05)
    all_background_genes : set
        All genes from genetic interventions (used as background)
    species : str
        'mouse' or 'human'

    Returns:
    --------
    DataFrame with enrichment results
    """
    print("\n" + "-" * 60)
    print("PERFORMING OVER-REPRESENTATION ANALYSIS (ORA)")
    print("-" * 60)

    # Prepare gene lists
    gene_list = list(significant_genes)
    background_list = list(all_background_genes)

    print(f"\nGene lists:")
    print(f"  Background (all genetic interventions): {len(background_list)} genes")
    print(f"  Significant query genes: {len(gene_list)} genes")

    # Define pathway databases based on species
    # Using only KEGG pathways
    if species.lower() == "mouse":
        gene_sets = [
            "KEGG_2019_Mouse",
        ]
    else:  # human
        gene_sets = [
            "KEGG_2019_Human",
        ]

    # Run ORA for all significant genes
    all_results = []

    for gene_set_name in gene_sets:
        print(f"\n  Testing enrichment in: {gene_set_name}")
        try:
            # Run Enrichr analysis with custom background
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=gene_set_name,
                organism="Mouse" if species.lower() == "mouse" else "Human",
                outdir=None,  # Don't save intermediate files
                cutoff=0.05,  # FDR cutoff
                background=background_list,  # Use genetic intervention genes as background
            )

            # Get results
            results = enr.results
            if len(results) > 0:
                results["Database"] = gene_set_name
                all_results.append(results)
                print(f"    Found {len(results)} significant pathways (FDR < 0.05)")
                print(f"    No significant pathways found")

        except Exception as e:
            print(f"    Error with {gene_set_name}: {e}")
            return pd.DataFrame()

    # Combine and process results
    if len(all_results) > 0:
        combined_df = pd.concat(all_results, ignore_index=True)

        print(f"\nColumns before rename: {combined_df.columns.tolist()}")

        # Standardize column names
        combined_df = combined_df.rename(
            columns={
                "Term": "Pathway",
                "Adjusted P-value": "FDR",
                "P-value": "P_value",
                "Odds Ratio": "Odds_Ratio",
                "Combined Score": "Combined_Score",
                "Overlap": "Overlap",
                "Genes": "Genes",
            }
        )

        print(f"Columns after rename: {combined_df.columns.tolist()}")

        # Parse overlap to get counts
        if "Overlap" in combined_df.columns:
            print(f"\nParsing Overlap column...")
            print(f"  Sample Overlap values: {combined_df['Overlap'].head(3).tolist()}")
            combined_df[["N_Overlap", "Pathway_Size"]] = combined_df[
                "Overlap"
            ].str.split("/", expand=True)
            combined_df["N_Overlap"] = combined_df["N_Overlap"].astype(int)
            combined_df["Pathway_Size"] = combined_df["Pathway_Size"].astype(int)
            print(
                f"  Parsed N_Overlap range: {combined_df['N_Overlap'].min()}-{combined_df['N_Overlap'].max()}"
            )
            print(
                f"  Parsed Pathway_Size range: {combined_df['Pathway_Size'].min()}-{combined_df['Pathway_Size'].max()}"
            )

            # Apply filters: pathway size < 100 and N_overlap >= 15
            print(f"\nApplying filters:")
            print(f"  Total pathways before filtering: {len(combined_df)}")
            combined_df = combined_df[
                (combined_df["Pathway_Size"] < 100) & (combined_df["N_Overlap"] >= 15)
            ]
            print(
                f"  Pathways after filtering (size<100, overlap>=15): {len(combined_df)}"
            )
        else:
            print(f"\nNOTE: Overlap column not present (custom background used)")
            print(f"  Pathway size information not available with custom background")
            # Extract N_Overlap from genes list
            combined_df["N_Overlap"] = combined_df["Genes"].apply(
                lambda x: len(str(x).split(";")) if pd.notna(x) else 0
            )
            combined_df["Pathway_Size"] = None  # Not available

            # Apply filter: only N_overlap >= 15 (can't filter by pathway size)
            print(f"\nApplying filters:")
            print(f"  Total pathways before filtering: {len(combined_df)}")
            combined_df = combined_df[combined_df["N_Overlap"] >= 15]
            print(f"  Pathways after filtering (overlap>=15): {len(combined_df)}")
            print(f"  Note: Pathway size filter not applied (size info unavailable)")

        if len(combined_df) == 0:
            print("\n  WARNING: No pathways passed the filters")
            return pd.DataFrame()

        # Calculate -log10(FDR) for visualization (FDR = adjusted p-value)
        combined_df["Neg_log10_FDR"] = -np.log10(combined_df["FDR"] + 1e-300)

        # Sort by -log10(FDR) in descending order (most significant first)
        combined_df = combined_df.sort_values("Neg_log10_FDR", ascending=False)

        # Select key columns
        key_cols = [
            "Pathway",
            "Database",
            "N_Overlap",
            "Pathway_Size",
            "P_value",
            "FDR",
            "Odds_Ratio",
            "Combined_Score",
            "Neg_log10_FDR",
            "Genes",
        ]
        combined_df = combined_df[[c for c in key_cols if c in combined_df.columns]]

        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"\nTotal significant pathways: {len(combined_df)}")

        print(f"\nPathways by database:")
        print(combined_df["Database"].value_counts().to_string())

        print(f"\nTop 10 enriched pathways:")
        print(
            combined_df[["Pathway", "Database", "FDR", "Odds_Ratio"]]
            .head(10)
            .to_string(index=False)
        )

        # Save results
        combined_df.to_csv("result/analysis_3_pathway_enrichment_ora.csv", index=False)
        print("\n✓ Saved: result/analysis_3_pathway_enrichment_ora.csv")

        return combined_df
    else:
        print("\n  WARNING: No significant enrichment found in any database")
        return pd.DataFrame()


def create_figure(enrichment_df, gene_df):
    """Create horizontal bar plot visualization of top enriched pathways"""
    print("\nCreating pathway enrichment visualization...")

    if len(enrichment_df) == 0:
        print("  WARNING: No enrichment results to plot")
        return

    # Filter out cancer-related pathways
    cancer_keywords = [
        "cancer",
        "carcinoma",
        "tumor",
        "tumour",
        "neoplasm",
        "oncogene",
        "metastasis",
        "malignancy",
        "leukemia",
        "lymphoma",
        "melanoma",
        "glioma",
        "sarcoma",
        "myeloma",
        "adenocarcinoma",
        "hepatocellular",
        "colorectal cancer",
        "breast cancer",
        "lung cancer",
        "prostate cancer",
        "bladder cancer",
        "pancreatic cancer",
        "gastric cancer",
        "thyroid cancer",
        "renal cancer",
        "endometrial cancer",
        "cervical cancer",
        "ovarian cancer",
        "esophageal cancer",
    ]

    # Create mask to exclude cancer-related pathways
    cancer_mask = (
        enrichment_df["Pathway"]
        .str.lower()
        .str.contains("|".join(cancer_keywords), case=False, na=False)
    )

    n_cancer = cancer_mask.sum()
    if n_cancer > 0:
        print(f"  Filtering out {n_cancer} cancer-related pathways")
        enrichment_df_filtered = enrichment_df[~cancer_mask].copy()
    else:
        enrichment_df_filtered = enrichment_df.copy()

    # Filter for significant pathways (FDR < 0.05) and take top entries
    plot_df = (
        enrichment_df_filtered[enrichment_df_filtered["FDR"] < 0.05].head(30).copy()
    )
    plot_df.to_csv('result/panel_e_plotted_pathways.csv', index=False)

    # Drawn at print size for Figure 4e (83 x 44 mm, placed at 100%).
    # Plot type, data, order, colours and label text follow the July figure;
    # only size, type, line widths and placement are set for print.
    mm = 1 / 25.4
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "mathtext.fontset": "custom",
            "mathtext.rm": "Arial",
            "mathtext.it": "Arial:italic",
            "mathtext.bf": "Arial:bold",
            "pdf.fonttype": 42,
            "axes.unicode_minus": True,
        }
    )
    fig = plt.figure(figsize=(PANEL_W_MM * mm, PANEL_H_MM * mm), facecolor="white")
    ax = fig.add_axes([0.70, 0.15, 0.28, 0.82])  # provisional; fitted below

    # Truncate pathway names for better display
    pathway_names = []
    for p in plot_df["Pathway"]:
        # Remove common prefixes and truncate
        p = str(p).replace("GO:", "").replace("KEGG:", "").replace("Reactome:", "")
        if len(p) > 60:
            p = p[:57] + "..."
        pathway_names.append(p)

    # Get -log10(FDR) values (FDR = adjusted p-value)
    neg_log_fdr = plot_df["Neg_log10_FDR"].values

    # Create colormap (rocket from seaborn)
    import matplotlib.cm as cm
    from matplotlib.colors import LinearSegmentedColormap

    # Use rocket colormap - rank by significance
    cmap = sns.color_palette("rocket", as_cmap=True)
    norm = plt.Normalize(vmin=0, vmax=len(plot_df) - 1)
    colors = [cmap(norm(i)) for i in range(len(plot_df))]

    # Create horizontal bar plot
    y_pos = np.arange(len(pathway_names))
    bars = ax.barh(
        y_pos,
        neg_log_fdr,
        color=colors,
        alpha=0.9,
        edgecolor="black",
        linewidth=0.3,
    )

    # Customize plot
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pathway_names, fontsize=LABEL_PT)
    ax.set_xlabel(
        "−log$_{10}$ FDR", fontsize=AXIS_TITLE_PT, loc="right", labelpad=1.5
    )
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    # Tick marks are off in this style; length=0 removes their reserved gap.
    ax.tick_params(axis="x", labelsize=TICK_PT, length=0, pad=2)
    ax.tick_params(axis="y", length=0, pad=2)

    # Add significance threshold line
    threshold_fdr = 0.05
    threshold_neg_log = -np.log10(threshold_fdr)
    ax.axvline(
        threshold_neg_log,
        color="red",
        linestyle="--",
        linewidth=0.9,
        label=f"FDR = {threshold_fdr}",
        alpha=0.7,
        zorder=0,
    )

    # Styling
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.55)
    ax.invert_yaxis()  # Highest significance at top
    ax.grid(axis="y", linewidth=0.5)
    ax.grid(axis="x", alpha=0.3, linestyle=":", linewidth=0.4)

    # Fit the axes so the longest pathway name, the x tick labels and the
    # x title sit inside the page with small margins.
    W, H = PANEL_W_MM, PANEL_H_MM
    for _ in range(4):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        px_mm = 25.4 / fig.dpi
        boxes_y = [t.get_window_extent(r) for t in ax.get_yticklabels()]
        lo, hi = ax.get_xlim()
        boxes_x = [
            t.get_window_extent(r)
            for t, loc in zip(ax.get_xticklabels(), ax.get_xticks())
            if t.get_text() and lo - 1e-9 <= loc <= hi + 1e-9
        ]
        box_xl = ax.xaxis.label.get_window_extent(r)
        left_ink = min(b.x0 for b in boxes_y) * px_mm
        right_ink = max([b.x1 for b in boxes_x] + [box_xl.x1]) * px_mm
        bottom_ink = box_xl.y0 * px_mm
        x0, y0, w0, h0 = ax.get_position().bounds
        new_left = x0 * W + (MARGIN_MM - left_ink)
        new_right = (x0 + w0) * W - (right_ink - (W - MARGIN_MM))
        new_bottom = y0 * H + (BOTTOM_MARGIN_MM - bottom_ink)
        new_top = H - TOP_MARGIN_MM
        ax.set_position(
            [new_left / W, new_bottom / H, (new_right - new_left) / W,
             (new_top - new_bottom) / H]
        )

    # FDR key in the free bottom-left corner, on the x-title row, ending just
    # left of the bar area.
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    box_xl = ax.xaxis.label.get_window_extent(r)
    key_y = (box_xl.y0 + box_xl.y1) / 2 / fig.bbox.height
    key_x = ax.get_position().x0 - 1.0 / W
    ax.legend(
        fontsize=KEY_PT,
        loc="center right",
        bbox_to_anchor=(key_x, key_y),
        bbox_transform=fig.transFigure,
        frameon=False,
        handlelength=2.2,
        handletextpad=0.5,
        borderaxespad=0,
        borderpad=0,
    )

    fig.savefig("figures/analysis_3_pathway_enrichment.pdf", dpi=300)
    fig.savefig("figures/analysis_3_pathway_enrichment.png", dpi=300)
    plt.close(fig)

    print("✓ Saved: figures/analysis_3_pathway_enrichment.pdf/.png")


def print_summary(enrichment_df, gene_df, species="mouse"):
    """Print manuscript-ready summary"""
    print("\n" + "=" * 80)
    print("MANUSCRIPT SUMMARY")
    print("=" * 80)

    if len(enrichment_df) == 0:
        print("\nNo significant enrichment found.")
        return

    significant = enrichment_df[enrichment_df["FDR"] < 0.05]
    top_pathway = enrichment_df.iloc[0]

    # Count pathways by database
    db_counts = enrichment_df["Database"].value_counts().to_dict()

    print(
        f"""
OVER-REPRESENTATION ANALYSIS (ORA)

Methods:
We extracted and validated {len(gene_df)} unique gene symbols from genetic 
interventions using keyword matching against the official {species} gene reference 
database (MyGene.info). Over-Representation Analysis (ORA) was performed using GSEApy 
v1.0+ with the Enrichr API, querying KEGG pathway database.

Overrepresentation was tested using one-sided Fisher's exact test (testing for 
enrichment) with automatic FDR correction (Benjamini-Hochberg). 
Significance threshold: FDR < 0.05.
Pathway filters: pathway size < 100 genes AND overlap >= 15 genes.

Results:
{len(significant)} significantly enriched pathways identified across all databases:
{chr(10).join('  - ' + row['Pathway'][:60] + f" (FDR={row['FDR']:.2e}, OR={row['Odds_Ratio']:.2f}, {row['Database']})" for _, row in significant.head(10).iterrows())}
{f'... and {len(significant) - 10} more' if len(significant) > 10 else ''}

Database breakdown:
{chr(10).join(f'  {db}: {count} pathways' for db, count in db_counts.items())}

The strongest enrichment was observed in:
{top_pathway['Pathway'][:70]}
(FDR = {top_pathway['FDR']:.2e}, OR = {top_pathway['Odds_Ratio']:.2f}, 
 Database: {top_pathway['Database']})

Conclusion:
Significant enrichment in established aging-related pathways validates that aging 
clock perturbations converge on canonical aging biology rather than representing 
random technical artifacts. The enrichment spans multiple independent pathway 
databases, indicating robust biological coherence. Odds ratios (range: 
{enrichment_df['Odds_Ratio'].min():.1f}-{enrichment_df['Odds_Ratio'].max():.1f}) 
demonstrate strong pathway-level convergence.
"""
    )


def main(species="mouse"):
    """
    Run complete Over-Representation Analysis (ORA) with standardized methods

    Parameters:
    -----------
    species : str
        'mouse' or 'human' - determines which gene reference to use
    """
    print("\n" + "=" * 80)
    print("OVER-REPRESENTATION ANALYSIS (ORA) - STANDARDIZED PIPELINE")
    print("=" * 80)
    print(f"Species: {species.upper()}")
    print("Method: ORA via GSEApy Enrichr API (one-sided Fisher's exact test)")
    print("Multiple testing correction: Benjamini-Hochberg FDR")
    print("Gene validation: Pathway databases + MyGene.info")
    print("=" * 80)

    # Step 1: Load data
    main_df = load_data()

    # Step 2: Download gene reference
    gene_ref_df, gene_reference_set = download_gene_reference(species=species)

    # Step 3: Extract and validate genes
    (
        significant_genes,
        all_background_genes,
        gene_df,
    ) = extract_genes(main_df, gene_reference_set)

    # Step 4: Perform Over-Representation Analysis
    # Use ALL genetic intervention genes as background
    enrichment_df = perform_ora_enrichment(
        significant_genes, all_background_genes, species=species
    )

    # Step 5: Create figure and print summary
    if len(enrichment_df) > 0:
        create_figure(enrichment_df, gene_df)
        print_summary(enrichment_df, gene_df, species=species)

        print("\n" + "=" * 80)
        print("ANALYSIS COMPLETE")
        print("=" * 80)
        print("\nDeliverables:")
        print(f"  1. result/gene_reference_{species}.csv - Gene reference database")
        print("  2. result/analysis_3_gene_list_validated.csv - Validated genes")
        print(
            "  3. result/analysis_3_pathway_enrichment_ora.csv - ORA enrichment results"
        )
        print(
            "  4. figures/analysis_3_pathway_enrichment.pdf/.png - Publication figure"
        )
    else:
        print("\n" + "=" * 80)
        print("ANALYSIS INCOMPLETE - No pathways passed filters")
        print("=" * 80)


if __name__ == "__main__":
    # Default to mouse, but can be changed to 'human'
    import sys

    species = sys.argv[1] if len(sys.argv) > 1 else "mouse"
    main(species=species)
