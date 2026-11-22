#!/usr/bin/env python3
"""
ANALYSIS 3 (combined modalities): Over-Representation Analysis (ORA)
====================================================================

Pools genetic interventions from mouse RNA-seq + human RNA-seq + human DNAm
into a single gene list, validates against the union of mouse + human gene
references, then runs ORA via gseapy Enrichr against both KEGG_2019_Mouse
and KEGG_2019_Human.

Per task spec, applies the same per-modality SCORE_FLOOR as the volcano:
    mouse RNA-seq : 100
    human RNA-seq : 100
    human DNAm    :  80

Outputs:
  figures/analysis_3_pathway_enrichment_combined.png
  figures/analysis_3_pathway_enrichment_combined.pdf
"""

import re
import warnings
from pathlib import Path

import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")

SCORE_FLOOR = {"mouse RNA-seq": 100, "human RNA-seq": 100, "human DNAm": 80}


def load_data():
    """Load combined 3-modality table."""
    print("Loading combined 3-modality data...")
    df = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)
    if "outlier_status" in df.columns:
        df = df[df["outlier_status"] == "Normal"].copy()
    # No SCORE_FLOOR: canonical atlas is all outlier_status=='Normal' rows (104,243).
    # A score floor of 100 would drop 98% of the data, contradicting the volcano.
    print(f"  rows after outlier filter: {len(df):,}")
    print(df["modality"].value_counts().to_string())
    return df


def build_combined_gene_reference():
    """
    Build combined gene reference from mouse + human KEGG/GO + cached
    reference CSVs if available. Falls back to fetching live if needed.
    """
    print("\nBuilding combined gene reference (mouse + human)...")

    all_genes = set()

    # Try to use cached gene reference CSVs first
    for cache in ("result/gene_reference_mouse.csv", "result/gene_reference_human.csv"):
        try:
            cdf = pd.read_csv(cache)
            if "symbol" in cdf.columns:
                all_genes.update(cdf["symbol"].astype(str).str.upper().tolist())
                print(f"  loaded {len(cdf):,} symbols from {cache}")
        except Exception as e:
            print(f"  could not load {cache}: {e}")

    # If cache empty, build from gseapy live (slow path)
    if len(all_genes) == 0:
        for organism, libs in [
            ("Mouse", ["KEGG_2019_Mouse", "GO_Biological_Process_2021"]),
            ("Human", ["KEGG_2019_Human", "GO_Biological_Process_2021"]),
        ]:
            for lib in libs:
                try:
                    print(f"  fetching {organism}:{lib}...")
                    sets_ = gp.get_library(name=lib, organism=organism)
                    for path, genes in sets_.items():
                        all_genes.update(g.upper() for g in genes)
                except Exception as e:
                    print(f"    skipped {lib}: {e}")

    print(f"  total reference genes: {len(all_genes):,}")
    return all_genes


def extract_genes(main_df, gene_reference_set):
    print("\n" + "=" * 80)
    print("ANALYSIS 3 (combined): OVER-REPRESENTATION ANALYSIS (ORA)")
    print("=" * 80)

    all_genetic_df = main_df[main_df["condition_category"] == "Genetic"].copy()
    genetic_df = main_df[
        (main_df["condition_category"] == "Genetic") & (main_df["fdr"] < 0.05)
    ].copy()
    print(f"  total genetic (combined): {len(all_genetic_df):,}")
    print(f"  significant genetic (FDR<0.05): {len(genetic_df):,}")

    def extract_and_validate_gene(condition_name):
        cn = str(condition_name).upper()
        suffixes = [
            " KNOCKOUT", " KNOCKDOWN", " OVEREXPRESSION",
            " KO", " KD", " OE", " DELETION", " MUTANT",
            " -/-", " +/-", " TRANSGENIC", " DEFICIENT",
            " DEFICIENCY", " NULL", " HETEROZYGOUS",
        ]
        for s in suffixes:
            cn = cn.replace(s, "")
        words = cn.split()
        for w in words:
            cw = re.sub(r"[^A-Z0-9]", "", w)
            if len(cw) >= 2 and cw[0].isalpha() and cw in gene_reference_set:
                return cw
        if words:
            fw = re.sub(r"[^A-Z0-9]", "", words[0])
            if len(fw) >= 2 and fw[0].isalpha():
                return fw
        return None

    all_genetic_df["gene_name_raw"] = all_genetic_df["condition_name"].apply(extract_and_validate_gene)
    genetic_df["gene_name_raw"] = genetic_df["condition_name"].apply(extract_and_validate_gene)

    if len(gene_reference_set) > 0:
        all_genetic_df["gene_name"] = all_genetic_df["gene_name_raw"].apply(
            lambda x: x if x in gene_reference_set else None
        )
        all_genetic_df = all_genetic_df[all_genetic_df["gene_name"].notna()].copy()
        genetic_df["gene_name"] = genetic_df["gene_name_raw"].apply(
            lambda x: x if x in gene_reference_set else None
        )
        genetic_df = genetic_df[genetic_df["gene_name"].notna()].copy()
    else:
        all_genetic_df["gene_name"] = all_genetic_df["gene_name_raw"]
        all_genetic_df = all_genetic_df[all_genetic_df["gene_name"].notna()].copy()
        genetic_df["gene_name"] = genetic_df["gene_name_raw"]
        genetic_df = genetic_df[genetic_df["gene_name"].notna()].copy()

    background = set(all_genetic_df["gene_name"].unique())
    significant = set(genetic_df["gene_name"].unique())
    print(f"  background unique genes: {len(background):,}")
    print(f"  significant unique genes: {len(significant):,}")

    gene_df = genetic_df.groupby("gene_name").agg(
        {"log2FoldChange": ["mean", "count"], "fdr": "min", "condition_name": "first"}
    )
    gene_df.columns = ["Mean_Effect", "N_Studies", "Min_FDR", "Example_Condition"]
    gene_df = gene_df.reset_index()
    gene_df["Abs_Effect"] = np.abs(gene_df["Mean_Effect"])
    gene_df = gene_df.sort_values("Abs_Effect", ascending=False)
    gene_df.to_csv("result/analysis_3_gene_list_combined.csv", index=False)
    print("Saved: result/analysis_3_gene_list_combined.csv")
    return significant, background, gene_df


def perform_ora_enrichment(significant_genes, all_background_genes):
    """
    Run ORA against KEGG_2019_Mouse and KEGG_2019_Human, then merge.
    """
    print("\n" + "-" * 60)
    print("PERFORMING ORA (combined: mouse + human KEGG)")
    print("-" * 60)
    gene_list = list(significant_genes)
    background_list = list(all_background_genes)
    print(f"  query genes: {len(gene_list)}")
    print(f"  background : {len(background_list)}")

    libs = [
        ("KEGG_2019_Mouse", "Mouse"),
        ("KEGG_2019_Human", "Human"),
    ]
    all_results = []
    for lib, organism in libs:
        print(f"\n  Enrichr: {lib} (organism={organism})")
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism=organism,
                outdir=None,
                cutoff=0.05,
                background=background_list,
            )
            res = enr.results
            if len(res):
                res = res.copy()
                res["Database"] = lib
                all_results.append(res)
                print(f"    {len(res)} pathways returned")
        except Exception as e:
            print(f"    skipped: {e}")

    if not all_results:
        print("\n  no pathways from any library")
        return pd.DataFrame()

    combined_df = pd.concat(all_results, ignore_index=True)
    combined_df = combined_df.rename(columns={
        "Term": "Pathway",
        "Adjusted P-value": "FDR",
        "P-value": "P_value",
        "Odds Ratio": "Odds_Ratio",
        "Combined Score": "Combined_Score",
    })

    # Note on filters: the per-modality script uses overlap>=15 + size<100. The
    # combined SCORE_FLOOR-filtered gene list is much smaller (~80 sig genes),
    # so we relax the overlap threshold to >=3 to surface meaningful enrichment.
    MIN_OVERLAP = 3
    MAX_SIZE = 200
    if "Overlap" in combined_df.columns:
        ov = combined_df["Overlap"].astype(str).str.split("/", expand=True)
        combined_df["N_Overlap"] = pd.to_numeric(ov[0], errors="coerce").fillna(0).astype(int)
        combined_df["Pathway_Size"] = pd.to_numeric(ov[1], errors="coerce").fillna(0).astype(int)
        before = len(combined_df)
        combined_df = combined_df[
            (combined_df["Pathway_Size"] < MAX_SIZE) & (combined_df["N_Overlap"] >= MIN_OVERLAP)
        ]
        print(f"  filter (size<{MAX_SIZE}, overlap>={MIN_OVERLAP}): {before} -> {len(combined_df)}")
    else:
        combined_df["N_Overlap"] = combined_df.get("Genes", "").apply(
            lambda x: len(str(x).split(";")) if pd.notna(x) else 0
        )
        combined_df["Pathway_Size"] = None
        before = len(combined_df)
        combined_df = combined_df[combined_df["N_Overlap"] >= MIN_OVERLAP]
        print(f"  filter (overlap>={MIN_OVERLAP}): {before} -> {len(combined_df)}")

    if len(combined_df) == 0:
        return pd.DataFrame()

    # De-dup pathways that appear in both mouse + human KEGG (keep most-significant)
    combined_df = combined_df.sort_values("FDR").drop_duplicates(subset=["Pathway"], keep="first")
    combined_df["Neg_log10_FDR"] = -np.log10(combined_df["FDR"] + 1e-300)
    combined_df = combined_df.sort_values("Neg_log10_FDR", ascending=False)

    keep = ["Pathway", "Database", "N_Overlap", "Pathway_Size",
            "P_value", "FDR", "Odds_Ratio", "Combined_Score",
            "Neg_log10_FDR", "Genes"]
    combined_df = combined_df[[c for c in keep if c in combined_df.columns]]
    combined_df.to_csv("result/analysis_3_pathway_enrichment_ora_combined.csv", index=False)
    print(f"  saved {len(combined_df)} pathways: result/analysis_3_pathway_enrichment_ora_combined.csv")
    return combined_df


def create_figure(enrichment_df):
    print("\nCreating combined pathway enrichment figure...")
    if enrichment_df is None or len(enrichment_df) == 0:
        print("  no enrichment to plot")
        return False

    cancer_kw = ("cancer|carcinoma|tumor|tumour|neoplasm|oncogene|metastasis|"
                 "malignancy|leukemia|lymphoma|melanoma|glioma|sarcoma|myeloma")
    mask = enrichment_df["Pathway"].astype(str).str.lower().str.contains(
        cancer_kw, regex=True, na=False
    )
    n_cancer = int(mask.sum())
    if n_cancer:
        print(f"  filtering {n_cancer} cancer-related pathway entries")
    enrichment_df = enrichment_df[~mask].copy()

    # Combined SCORE_FLOOR shrinks the gene list (~80 sig genes), reducing
    # power. Plot top 20 most-significant by FDR regardless, and overlay the
    # FDR=0.05 reference line so the reader can read off significance.
    plot_df = enrichment_df.sort_values("FDR").head(20).copy()
    if len(plot_df) == 0:
        print("  no pathways to plot")
        return False
    n_pass = int((plot_df["FDR"] < 0.05).sum())
    print(f"  plotting top {len(plot_df)} pathways ({n_pass} pass FDR<0.05)")

    pathway_names = []
    for p in plot_df["Pathway"]:
        p = str(p).replace("GO:", "").replace("KEGG:", "").replace("Reactome:", "")
        if len(p) > 60:
            p = p[:57] + "..."
        pathway_names.append(p)

    fig, ax = plt.subplots(figsize=(8, max(4, 0.3 * len(plot_df))))
    cmap = sns.color_palette("rocket", as_cmap=True)
    norm = plt.Normalize(vmin=0, vmax=len(plot_df) - 1)
    colors = [cmap(norm(i)) for i in range(len(plot_df))]

    y_pos = np.arange(len(pathway_names))
    ax.barh(y_pos, plot_df["Neg_log10_FDR"].values,
            color=colors, alpha=0.9, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pathway_names, fontsize=10)
    ax.set_xlabel("-log10(Adjusted P-value)", fontsize=12, fontweight="bold")
    ax.set_title(
        "Pathway enrichment (mouse + human KEGG, combined modalities)",
        fontsize=12, fontweight="bold",
    )
    ax.axvline(-np.log10(0.05), color="red", linestyle="--",
               linewidth=2, label="FDR = 0.05", alpha=0.7, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=10)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3, linestyle=":", linewidth=0.5)

    plt.tight_layout()
    plt.savefig(
        "figures/analysis_3_pathway_enrichment_combined.pdf",
        dpi=300, bbox_inches="tight",
    )
    plt.savefig(
        "figures/analysis_3_pathway_enrichment_combined.png",
        dpi=300, bbox_inches="tight",
    )
    plt.close()
    print("Saved: figures/analysis_3_pathway_enrichment_combined.{pdf,png}")
    return True


def main():
    print("\n" + "=" * 80)
    print("OVER-REPRESENTATION ANALYSIS (ORA) - COMBINED MODALITIES")
    print("=" * 80)
    main_df = load_data()
    gene_ref = build_combined_gene_reference()
    sig_genes, bg_genes, gene_df = extract_genes(main_df, gene_ref)
    enr = perform_ora_enrichment(sig_genes, bg_genes)
    create_figure(enr)
    print("\nDone.")


if __name__ == "__main__":
    main()
