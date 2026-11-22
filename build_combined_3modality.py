#!/usr/bin/env python3
"""Refresh figures/combined_3modality_plot_data.csv final_score from the CURRENT
compiled scores (including plausibility q0 corrections), preserving the comparison
structure (104,243 Normal rows, multi-clock DNAm, directions, outlier flags).

Join keys differ by modality because the combined file lacks human case_id:
  mouse RNA-seq : by case_id  -> data/scores_summary.csv (dataset)
  human RNA-seq : by (gse_id, comparison_name, condition_name) -> with_scores
  human DNAm    : by (gse_id, comparison_name, condition_name) -> with_scores
Unscored comparisons get NaN final_score (as before).
"""
import numpy as np
import pandas as pd

DEST = "figures/combined_3modality_plot_data.csv"
KEY = ["gse_id", "comparison_name", "condition_name"]


def main():
    df = pd.read_csv(DEST, low_memory=False)
    fs = pd.Series(np.nan, index=df.index)

    # mouse: by case_id
    mm = pd.read_csv("data/scores_summary.csv")[["dataset", "final_score"]].dropna(subset=["dataset"])
    mmap = dict(zip(mm["dataset"].astype(str), mm["final_score"]))
    mask = df["modality"] == "mouse RNA-seq"
    fs.loc[mask] = df.loc[mask, "case_id"].astype(str).map(mmap)
    print(f"mouse: {fs.loc[mask].notna().sum():,} scored")

    # human: by (gse_id, comparison_name, condition_name)
    for mod, wf in [("human RNA-seq", "result/main_results_hs_rna_with_scores.csv"),
                    ("human DNAm", "result/main_results_hs_dnam_with_scores.csv")]:
        w = pd.read_csv(wf, low_memory=False)[KEY + ["final_score"]].dropna(subset=KEY)
        w = w.drop_duplicates(subset=KEY)
        mask = df["modality"] == mod
        merged = df.loc[mask, KEY].merge(w, on=KEY, how="left")
        fs.loc[mask] = merged["final_score"].values
        print(f"{mod}: {merged['final_score'].notna().sum():,} scored")

    df["final_score"] = fs
    df["neg_log10_fdr"] = -np.log10(pd.to_numeric(df["fdr"], errors="coerce").clip(lower=1e-300))
    df.to_csv(DEST, index=False)
    n = df[df.outlier_status == "Normal"]
    print(f"written {DEST}: {len(df):,} rows | Normal {len(n):,} | FDR<0.05 (Normal) {int((n.fdr<0.05).sum()):,}")


if __name__ == "__main__":
    main()
