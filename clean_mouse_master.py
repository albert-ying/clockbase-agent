#!/usr/bin/env python3
"""
Regenerate the mouse RNA master table with the CONTAMINATION FILTER applied.

The agent's per-GSE `aging_analysis.csv` files occasionally carry rows whose
`gse_id` column disagrees with the parent folder name. These "bad" rows ship
wrong-direction effect sizes and pollute every downstream analysis. ~10% of
mouse RNA rows are affected.

This script:
  1. Walks result_mm_rna/*/aging_analysis.csv
  2. Drops rows where gse_id != parent folder name
  3. Concatenates, filters to mouse's published clock, computes log2FC +
     outlier flags the same way mouse s01 does
  4. Writes result/main_results_result_mm_rna_outlier_marked_clean.csv

Does NOT modify s01_visualize_results.py. If the mouse downstream analyses (analysis_1..4, etc.) need to use
this cleaner master, point their `pd.read_csv(...)` to the _clean.csv path.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

BASE_DIR = "result_mm_rna"
INPUT_FILENAME = "aging_analysis.csv"
OUTPUT_PATH = "result/main_results_result_mm_rna_outlier_marked_clean.csv"
# The single clock that mouse s01 filters to (from s01_visualize_results.py:512)
MOUSE_CLOCK = "Hazard_Mouse_All_All_scaleddiff_EN_tAge"


def load_one(fp):
    required = [
        "gse_id", "comparison_name", "condition_name", "condition_category",
        "clock", "effect_size", "p_value", "trust_score", "direction",
    ]
    try:
        df = pd.read_csv(fp)
        if not all(c in df.columns for c in required):
            return None
        df["effect_size"] = pd.to_numeric(df["effect_size"], errors="coerce")
        df["p_value"]    = pd.to_numeric(df["p_value"], errors="coerce")
        df["trust_score"] = pd.to_numeric(df["trust_score"], errors="coerce")
        df["direction"]  = pd.to_numeric(df["direction"], errors="coerce")
        df["log2FoldChange"] = np.abs(df["effect_size"]) * df["direction"]
        df = df.rename(columns={"p_value": "pvalue"})
        df = df.dropna(subset=["log2FoldChange", "pvalue"])
        df = df[df["pvalue"].between(0, 1)]
        if len(df) == 0:
            return None
        return df
    except Exception as e:
        return None


def main():
    print(f"walking {BASE_DIR}/*/**{INPUT_FILENAME} (this takes a minute for 17k dirs)...")
    frames = []
    total = 0
    dropped = 0
    bad_files = 0
    for d in sorted(Path(BASE_DIR).iterdir()):
        if not d.is_dir():
            continue
        f = d / INPUT_FILENAME
        if not f.exists():
            continue
        df = load_one(f)
        if df is None:
            continue
        # Contamination filter: enforce gse_id == folder name
        before = len(df)
        keep = df["gse_id"].astype(str).str.strip() == d.name
        df = df[keep].copy()
        dropped_this = before - len(df)
        if dropped_this > 0:
            bad_files += 1
        total += before
        dropped += dropped_this
        if len(df) > 0:
            frames.append(df)

    if total > 0:
        print(f"  contamination dropped: {dropped:,}/{total:,} rows "
              f"({100*dropped/total:.1f}%) across {bad_files:,} files")
    print(f"  clean rows: {sum(len(f) for f in frames):,}")

    combined = pd.concat(frames, ignore_index=True)

    # Single-clock filter (matches mouse s01 line 512)
    combined = combined[combined["clock"] == MOUSE_CLOCK].copy()
    print(f"  after clock filter ({MOUSE_CLOCK}): {len(combined):,}")

    # FDR correction
    combined["neg_log10_pvalue"] = -np.log10(combined["pvalue"].clip(lower=1e-300))
    mask = combined["pvalue"].notna()
    min_nonzero = combined.loc[mask & (combined["pvalue"] > 0), "pvalue"].min()
    combined.loc[mask & (combined["pvalue"] == 0), "pvalue"] = min_nonzero
    _, fdr = fdrcorrection(combined.loc[mask, "pvalue"].values)
    combined.loc[mask, "fdr"] = fdr

    # IQR-based outlier flags on log2FC and on −log10(p), 10× fence
    def iqr_bounds(series, mult=10):
        q1, q3 = series.quantile(0.25), series.quantile(0.75)
        iqr = q3 - q1
        return q1 - mult * iqr, q3 + mult * iqr

    lf, uf = iqr_bounds(combined["log2FoldChange"])
    lp, up = iqr_bounds(combined["neg_log10_pvalue"])

    def status(row):
        ef = (row["log2FoldChange"] < lf) or (row["log2FoldChange"] > uf)
        pf = (row["neg_log10_pvalue"] < lp) or (row["neg_log10_pvalue"] > up)
        if ef and pf: return "Both Effect Size & P-value Outlier"
        if ef:        return "Effect Size Outlier"
        if pf:        return "P-value Outlier"
        return "Normal"

    combined["outlier_status"] = combined.apply(status, axis=1)
    combined["abs_log2FoldChange"] = np.abs(combined["log2FoldChange"])

    Path(OUTPUT_PATH).parent.mkdir(exist_ok=True)
    combined.to_csv(OUTPUT_PATH, index=False)

    # Summary
    print(f"\n✓ wrote {OUTPUT_PATH}  ({len(combined):,} rows)")
    print(combined["outlier_status"].value_counts())


if __name__ == "__main__":
    main()
