#!/usr/bin/env python3
"""
s02_prepare_scoring_input_human.py
===================================

Filters the master human results table down to significant non-outlier
comparisons (FDR < 0.05 AND outlier_status == "Normal") and writes a
case_mapping-style CSV ready for the LLM scoring agent.

Mirrors the mouse scoring workflow (see data/case_mapping.csv +
data/scores_summary.csv). Flip MODALITY to process DNAm.

Output:
  result/for_scoring_hs_{MODALITY}.csv  — one row per (gse_id, comparison,
      condition, clock) that needs scoring, with case_id = {gse_id}-{idx}.
"""

# ---- HUMAN CONFIG --------------------------------------------------
MODALITY = "rna"  # "rna" or "dnam"
METHOD_NAME = f"hs_{MODALITY}"
FDR_THRESHOLD = 0.05

# Only these four categories go to the scoring agent. Everything else is
# lumped into an "Other" bucket that is not scored (saved separately for
# record-keeping but not included in the scoring list).
SCORING_CATEGORIES = ["Disease", "Environment", "Drug", "Genetic"]

# DNAm-only filtering. GrimAgeV2 is the primary clock. A comparison enters
# the scoring list iff its GrimAgeV2 row has FDR < FDR_THRESHOLD AND at least
# DNAM_MIN_SUPPORT of the other 4 clocks have pvalue < DNAM_SUPPORT_PVALUE.
DNAM_PRIMARY_CLOCK = "GrimAgeV2"
DNAM_OTHER_CLOCKS = ["PhenoAge", "Horvathv1", "YingDamAge", "DunedinPACE"]
DNAM_SUPPORT_PVALUE = 0.05
DNAM_MIN_SUPPORT = 2
# --------------------------------------------------------------------

from pathlib import Path

import pandas as pd


def main():
    master_path = f"result/main_results_{METHOD_NAME}_outlier_marked.csv"
    out_path = f"result/for_scoring_{METHOD_NAME}.csv"

    print(f"Loading {master_path}...")
    df = pd.read_csv(master_path, low_memory=False)
    print(f"  {len(df):,} total rows")

    # Restrict to scoreable categories (everything else is lumped into "Other"
    # and not sent to the scorer).
    before = len(df)
    df = df[df["condition_category"].isin(SCORING_CATEGORIES)].copy()
    print(
        f"  after category ∈ {SCORING_CATEGORIES}: {len(df):,} "
        f"({before - len(df):,} dropped as non-scoreable 'Other')"
    )

    if MODALITY == "rna":
        # RNA: single clock, just require FDR < threshold. Outliers kept.
        before = len(df)
        df = df[df["fdr"] < FDR_THRESHOLD].copy()
        print(f"  after fdr<{FDR_THRESHOLD}: {len(df):,} ({before - len(df):,} dropped)")
    elif MODALITY == "dnam":
        # DNAm: GrimAgeV2 primary + cross-clock support + directionality check.
        # 1) GrimAgeV2 FDR < threshold
        # 2) In at least DNAM_MIN_SUPPORT of the 4 other clocks:
        #       pvalue < DNAM_SUPPORT_PVALUE  AND  direction matches GrimAgeV2
        # 3) Output only the GrimAgeV2 row per passing comparison, with extra
        #    columns carrying each other clock's pvalue + effect_size.
        import numpy as np
        group_keys = ["gse_id", "comparison_name", "condition_name"]
        all_clocks = [DNAM_PRIMARY_CLOCK] + DNAM_OTHER_CLOCKS

        # Pivot to wide: one row per comparison, clock metrics as columns
        sub = df[df["clock"].isin(all_clocks)].copy()
        # Use first() to collapse any accidental duplicates per (comparison, clock)
        wide_parts = {}
        for metric in ["pvalue", "effect_size", "direction", "log2FoldChange", "fdr"]:
            w = (
                sub.pivot_table(
                    index=group_keys,
                    columns="clock",
                    values=metric,
                    aggfunc="first",
                )
            )
            wide_parts[metric] = w

        # Flatten into a single wide DataFrame with f"{clock}_{metric}" columns
        wide = wide_parts["pvalue"].iloc[:, :0].copy()  # just the index
        for metric, w in wide_parts.items():
            for clk in w.columns:
                wide[f"{clk}_{metric}"] = w[clk]
        wide = wide.reset_index()

        # Bring along per-comparison metadata from the GrimAgeV2 row
        meta_cols = [
            c
            for c in [
                "condition_category",
                "covariate",
                "statistical_test",
                "trust_score",
                "neg_log10_pvalue",
                "outlier_status",
            ]
            if c in sub.columns
        ]
        meta = (
            sub[sub["clock"] == DNAM_PRIMARY_CLOCK][group_keys + meta_cols]
            .drop_duplicates(subset=group_keys)
        )
        wide = wide.merge(meta, on=group_keys, how="left")

        print(f"  comparisons with GrimAgeV2 data: {wide[f'{DNAM_PRIMARY_CLOCK}_fdr'].notna().sum():,}")

        # 1) GrimAgeV2 FDR < threshold
        before = len(wide)
        wide = wide[wide[f"{DNAM_PRIMARY_CLOCK}_fdr"] < FDR_THRESHOLD].copy()
        print(f"  after {DNAM_PRIMARY_CLOCK} FDR<{FDR_THRESHOLD}: {len(wide):,}")

        # 2) Cross-clock support with directionality agreement
        grim_sign = np.sign(wide[f"{DNAM_PRIMARY_CLOCK}_log2FoldChange"])
        support = np.zeros(len(wide), dtype=int)
        for clk in DNAM_OTHER_CLOCKS:
            pcol = f"{clk}_pvalue"
            lcol = f"{clk}_log2FoldChange"
            if pcol not in wide.columns or lcol not in wide.columns:
                continue
            this_sig = wide[pcol] < DNAM_SUPPORT_PVALUE
            this_match = np.sign(wide[lcol]) == grim_sign
            support = support + (this_sig.fillna(False) & this_match.fillna(False)).astype(int).to_numpy()
        wide["n_support"] = support
        before = len(wide)
        wide = wide[wide["n_support"] >= DNAM_MIN_SUPPORT].copy()
        print(
            f"  after ≥{DNAM_MIN_SUPPORT} of {DNAM_OTHER_CLOCKS} with p<{DNAM_SUPPORT_PVALUE} "
            f"AND same direction as {DNAM_PRIMARY_CLOCK}: {len(wide):,} ({before - len(wide):,} dropped)"
        )

        # Promote GrimAgeV2 columns to the "primary" schema the scorer expects
        rename_map = {
            f"{DNAM_PRIMARY_CLOCK}_pvalue": "pvalue",
            f"{DNAM_PRIMARY_CLOCK}_effect_size": "effect_size",
            f"{DNAM_PRIMARY_CLOCK}_direction": "direction",
            f"{DNAM_PRIMARY_CLOCK}_log2FoldChange": "log2FoldChange",
            f"{DNAM_PRIMARY_CLOCK}_fdr": "fdr",
        }
        wide = wide.rename(columns=rename_map)
        wide["clock"] = DNAM_PRIMARY_CLOCK
        df = wide
    else:
        raise ValueError(f"Unknown MODALITY: {MODALITY}")

    # Assign a unique case_id = {gse_id}-{row_index}. Stable across runs
    # as long as the input order is stable; re-running s01 will re-assign.
    df = df.reset_index(drop=True)
    df["case_id"] = df.apply(
        lambda r: f"{r['gse_id']}-{r.name}", axis=1
    )

    # Columns that the mouse scoring agent's sample.json consumes.
    # Primary clock metrics (pvalue / effect_size / ... without prefix) come
    # from the single RNA clock or from GrimAgeV2 for DNAm.
    keep_cols = [
        "case_id",
        "gse_id",
        "comparison_name",
        "condition_name",
        "condition_category",
        "covariate",
        "clock",
        "statistical_test",
        "pvalue",
        "effect_size",
        "direction",
        "log2FoldChange",
        "fdr",
        "neg_log10_pvalue",
        "outlier_status",
        "trust_score",
        "n_support",  # DNAm: count of supporting non-primary clocks (direction-matched, p<DNAM_SUPPORT_PVALUE)
    ]
    # DNAm: append per-other-clock pvalue + effect_size columns for the scorer
    if MODALITY == "dnam":
        for clk in DNAM_OTHER_CLOCKS:
            keep_cols.extend([f"{clk}_pvalue", f"{clk}_effect_size"])
    keep_cols = [c for c in keep_cols if c in df.columns]
    out = df[keep_cols].copy()

    Path("result").mkdir(exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"\n✓ Saved: {out_path}  ({len(out):,} rows, {len(keep_cols)} columns)")

    # Summary stats for a sanity check before sending to scorer
    print("\n--- Summary ---")
    print(f"Unique GSE studies: {out['gse_id'].nunique():,}")
    print(f"Unique condition_name: {out['condition_name'].nunique():,}")
    print(f"\nBy condition_category:")
    print(out["condition_category"].value_counts().head(10))
    anti = (out["log2FoldChange"] < 0).sum()
    pro = (out["log2FoldChange"] > 0).sum()
    print(
        f"\nDirection: anti-aging (log2FC<0) = {anti:,}   pro-aging (log2FC>0) = {pro:,}"
    )
    print(f"Top 10 anti-aging hits (lowest log2FoldChange):")
    top_anti = out.nsmallest(10, "log2FoldChange")[
        ["case_id", "condition_name", "condition_category", "log2FoldChange", "fdr"]
    ]
    print(top_anti.to_string(index=False))


if __name__ == "__main__":
    main()
