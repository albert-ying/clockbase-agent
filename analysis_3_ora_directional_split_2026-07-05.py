#!/usr/bin/env python3
"""
ANALYSIS 3 (directional split): KEGG ORA for pro-aging vs anti-aging genetic interventions
==========================================================================================

Reproduces the manuscript Figure 4e KEGG Over-Representation Analysis (ORA) exactly, then
SPLITS the significant genetic gene set by direction to answer the question:
  "Is pathway enrichment pronounced for both pro- and anti-aging interventions, and is
   there a difference (e.g. pro-aging more in senescence, anti-aging more in insulin)?"

Method is IDENTICAL to analysis_3_pathway_enrichment.py:
  - Input: result/main_results_result_mm_rna_outlier_marked.csv, outlier_status=='Normal'
  - Significant set: condition_category=='Genetic' & fdr<0.05
  - Background: ALL Genetic interventions (Normal)
  - Gene symbol extracted from condition_name by keyword matching + validation against
    the official mouse gene reference (result/gene_reference_mouse.csv, the manuscript's).
  - ORA: gseapy Enrichr, KEGG_2019_Mouse, organism=Mouse, custom background = all genetic
    intervention genes. One-sided Fisher, BH FDR. Filter overlap N>=15.

Direction convention (from analysis_2): log2FoldChange > 0 = age-INCREASING (pro-aging),
                                        log2FoldChange < 0 = age-LOWERING (anti-aging).

A gene may appear in BOTH lists if it is significant in one study in the pro-aging direction
and in another study in the anti-aging direction. We keep it in both (per-comparison
assignment) and report the overlap explicitly. NO fabrication: only real extracted gene
symbols and real Enrichr results are used.
"""

import re
import warnings
from datetime import date

import gseapy as gp
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

DATA = "result/main_results_result_mm_rna_outlier_marked.csv"
GENE_REF = "result/gene_reference_mouse.csv"
KEGG_LIB = "KEGG_2019_Mouse"
TODAY = "2026-07-05"

SUFFIXES = [
    " KNOCKOUT", " KNOCKDOWN", " OVEREXPRESSION", " KO", " KD", " OE",
    " DELETION", " MUTANT", " -/-", " +/-", " TRANSGENIC", " DEFICIENT",
    " DEFICIENCY", " NULL", " HETEROZYGOUS",
]


def extract_and_validate_gene(condition_name, gene_reference_set):
    """EXACT copy of the manuscript extractor in analysis_3_pathway_enrichment.py."""
    condition_name = str(condition_name).upper()
    for suffix in SUFFIXES:
        condition_name = condition_name.replace(suffix, "")
    words = condition_name.split()
    for word in words:
        clean_word = re.sub(r"[^A-Z0-9]", "", word)
        if len(clean_word) >= 2 and clean_word[0].isalpha():
            if clean_word in gene_reference_set:
                return clean_word
    if len(words) > 0:
        first_word = re.sub(r"[^A-Z0-9]", "", words[0])
        if len(first_word) >= 2 and first_word[0].isalpha():
            return first_word
    return None


def run_ora(gene_list, background_list, tag):
    """Run KEGG_2019_Mouse Enrichr ORA with custom background; return cleaned df."""
    print(f"\n--- ORA [{tag}]: {len(gene_list)} query genes, {len(background_list)} background ---")
    enr = gp.enrichr(
        gene_list=list(gene_list),
        gene_sets=KEGG_LIB,
        organism="mouse",
        outdir=None,
        cutoff=0.05,
        background=list(background_list),
    )
    res = enr.results.copy()
    res["Database"] = KEGG_LIB
    res = res.rename(columns={
        "Term": "Pathway", "Adjusted P-value": "FDR", "P-value": "P_value",
        "Odds Ratio": "Odds_Ratio", "Combined Score": "Combined_Score",
    })
    # Custom background -> no numeric Overlap denominator; derive N_Overlap from Genes
    if "Overlap" in res.columns and res["Overlap"].astype(str).str.contains("/").all():
        parts = res["Overlap"].str.split("/", expand=True)
        res["N_Overlap"] = parts[0].astype(int)
        res["Pathway_Size"] = parts[1].astype(int)
    else:
        res["N_Overlap"] = res["Genes"].apply(
            lambda x: len(str(x).split(";")) if pd.notna(x) else 0
        )
        res["Pathway_Size"] = None
    res = res[res["N_Overlap"] >= 15].copy()
    res["Neg_log10_FDR"] = -np.log10(res["FDR"] + 1e-300)
    res = res.sort_values("Neg_log10_FDR", ascending=False)
    keep = ["Pathway", "Database", "N_Overlap", "Pathway_Size", "P_value", "FDR",
            "Odds_Ratio", "Combined_Score", "Neg_log10_FDR", "Genes"]
    res = res[[c for c in keep if c in res.columns]]
    print(f"    pathways with overlap>=15: {len(res)} | FDR<0.05: {(res['FDR']<0.05).sum()}")
    return res


def main():
    print("Loading data and gene reference...")
    df = pd.read_csv(DATA, low_memory=False)
    df = df[df["outlier_status"] == "Normal"].copy()

    gene_ref = set(pd.read_csv(GENE_REF)["symbol"].astype(str).str.upper())
    print(f"  gene reference symbols: {len(gene_ref):,}")

    all_genetic = df[df["condition_category"] == "Genetic"].copy()
    sig = df[(df["condition_category"] == "Genetic") & (df["fdr"] < 0.05)].copy()
    print(f"  all genetic (Normal): {len(all_genetic):,}")
    print(f"  significant genetic (fdr<0.05): {len(sig):,}")

    # Extract + validate gene per row
    for frame in (all_genetic, sig):
        frame["gene_raw"] = frame["condition_name"].apply(
            lambda c: extract_and_validate_gene(c, gene_ref))
        frame["gene_name"] = frame["gene_raw"].apply(
            lambda x: x if x in gene_ref else None)
    all_genetic = all_genetic[all_genetic["gene_name"].notna()].copy()
    sig = sig[sig["gene_name"].notna()].copy()

    background = set(all_genetic["gene_name"].unique())
    sig_genes = set(sig["gene_name"].unique())
    print(f"\n  validated background genes: {len(background)}")
    print(f"  validated significant genes: {len(sig_genes)}")

    # Directional split (per-comparison assignment)
    pro_rows = sig[sig["log2FoldChange"] > 0]
    anti_rows = sig[sig["log2FoldChange"] < 0]
    pro_genes = set(pro_rows["gene_name"].unique())
    anti_genes = set(anti_rows["gene_name"].unique())
    both = pro_genes & anti_genes
    print(f"\n  pro-aging rows (log2FC>0): {len(pro_rows)} -> {len(pro_genes)} unique genes")
    print(f"  anti-aging rows (log2FC<0): {len(anti_rows)} -> {len(anti_genes)} unique genes")
    print(f"  genes in BOTH directions: {len(both)}")

    # 1) Verify combined ORA reproduces manuscript key pathways
    combined = run_ora(sig_genes, background, "COMBINED (verify vs Fig4e)")
    key = ["Cellular senescence", "Longevity regulating pathway", "Th17 cell differentiation",
           "Insulin signaling pathway", "Autophagy"]
    print("\n  Verification of manuscript key pathways (combined set):")
    for k in key:
        row = combined[combined["Pathway"].str.contains(k, case=False, na=False)]
        if len(row):
            r = row.iloc[0]
            print(f"    {k:40s} FDR={r['FDR']:.3g}  overlap={r['N_Overlap']}  OR={r['Odds_Ratio']:.2f}")
        else:
            print(f"    {k:40s} NOT in overlap>=15 set")

    # 2) Directional ORAs
    pro = run_ora(pro_genes, background, "PRO-AGING (age-increasing)")
    anti = run_ora(anti_genes, background, "ANTI-AGING (age-lowering)")

    pro_path = f"result/analysis_3_ora_pro_aging_KEGG_{TODAY}.csv"
    anti_path = f"result/analysis_3_ora_anti_aging_KEGG_{TODAY}.csv"
    pro.to_csv(pro_path, index=False)
    anti.to_csv(anti_path, index=False)
    combined.to_csv(f"result/analysis_3_ora_combined_KEGG_{TODAY}.csv", index=False)
    print(f"\nSaved:\n  {pro_path}\n  {anti_path}")

    # 3) Comparison report
    def top(df_, n=15):
        return df_[df_["FDR"] < 0.05].head(n)

    print("\n" + "=" * 90)
    print(f"PRO-AGING (age-increasing): {len(pro_genes)} genes | sig KEGG pathways FDR<0.05: {(pro['FDR']<0.05).sum()}")
    print("=" * 90)
    for _, r in top(pro).iterrows():
        print(f"  {r['Pathway'][:48]:48s} FDR={r['FDR']:.2e} n={r['N_Overlap']:>2} OR={r['Odds_Ratio']:.2f}")

    print("\n" + "=" * 90)
    print(f"ANTI-AGING (age-lowering): {len(anti_genes)} genes | sig KEGG pathways FDR<0.05: {(anti['FDR']<0.05).sum()}")
    print("=" * 90)
    for _, r in top(anti).iterrows():
        print(f"  {r['Pathway'][:48]:48s} FDR={r['FDR']:.2e} n={r['N_Overlap']:>2} OR={r['Odds_Ratio']:.2f}")

    pro_sig = set(pro[pro["FDR"] < 0.05]["Pathway"])
    anti_sig = set(anti[anti["FDR"] < 0.05]["Pathway"])
    print("\n" + "=" * 90)
    print("SHARED (FDR<0.05 in both directions):")
    for p in sorted(pro_sig & anti_sig):
        print(f"  {p}")
    print("\nPRO-AGING ONLY (FDR<0.05):")
    for p in sorted(pro_sig - anti_sig):
        print(f"  {p}")
    print("\nANTI-AGING ONLY (FDR<0.05):")
    for p in sorted(anti_sig - pro_sig):
        print(f"  {p}")

    # Targeted pathways for Alex's hypothesis (report FDR in each direction regardless of sig)
    print("\n" + "=" * 90)
    print("ALEX'S TARGET PATHWAYS (FDR in each direction):")
    targets = ["Cellular senescence", "Longevity regulating pathway", "Insulin signaling pathway",
               "Autophagy", "Th17 cell differentiation", "FoxO signaling pathway",
               "mTOR signaling pathway", "AMPK signaling pathway", "NF-kappa B signaling pathway",
               "PI3K-Akt signaling pathway", "JAK-STAT signaling pathway"]
    print(f"  {'Pathway':42s} {'PRO FDR':>10s} {'PRO n':>6s} {'ANTI FDR':>10s} {'ANTI n':>6s}")
    for t in targets:
        pr = pro[pro["Pathway"].str.contains(t, case=False, na=False)]
        an = anti[anti["Pathway"].str.contains(t, case=False, na=False)]
        pf = f"{pr.iloc[0]['FDR']:.2e}" if len(pr) else "  -  "
        pn = f"{pr.iloc[0]['N_Overlap']}" if len(pr) else "-"
        af = f"{an.iloc[0]['FDR']:.2e}" if len(an) else "  -  "
        an_ = f"{an.iloc[0]['N_Overlap']}" if len(an) else "-"
        print(f"  {t:42s} {pf:>10s} {pn:>6s} {af:>10s} {an_:>6s}")


if __name__ == "__main__":
    main()
