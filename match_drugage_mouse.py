#!/usr/bin/env python3
"""
Match mouse RNA-seq Drug interventions against DrugAge and write
data/agent_result_drugage_significant.csv for the combined waffle-plot step.

Mouse counterpart of match_drugage_human.py. Closes a reproducibility gap:
data/agent_result_drugage_significant.csv previously had no builder script.
Reuses the (bug-fixed) word-boundary drug matcher and reads the canonical
outlier-removed mouse analytic table (detailed_scores_result_mm_rna.csv).
"""

from pathlib import Path

import pandas as pd

from match_drugage_human import (
    build_drugage_index,
    match_drug,
    infer_drug_influence,
    DRUGAGE_CSV,
)

INPUT = "result/detailed_scores_result_mm_rna.csv"
OUT = "data/agent_result_drugage_significant.csv"


def main():
    print(f"loading DrugAge: {DRUGAGE_CSV}")
    da = pd.read_csv(DRUGAGE_CSV)
    da_index = build_drugage_index(da)
    print(f"  {len(da_index):,} unique compounds in DrugAge index")

    print(f"loading {INPUT}")
    df = pd.read_csv(INPUT, low_memory=False)
    df = df[df["condition_category"] == "Drug"].copy()
    print(f"  {len(df):,} Drug rows (canonical outlier-removed mouse set)")
    # Manuscript DrugAge set is FDR<0.05 (score floor applied downstream by the
    # waffle script via trust_score).
    df = df[df["fdr"] < 0.05].copy()
    print(f"  {len(df):,} FDR-significant Drug rows")

    matches = df.apply(
        lambda r: match_drug(
            f"{r.get('intervention','')}  {r.get('condition_name','')}", da_index),
        axis=1,
    )
    df["drugage_matched_drugs"] = matches.apply(lambda x: x[0])
    df["drugage_annotation"] = matches.apply(lambda x: x[1])
    df["drugage_annotated"] = df["drugage_matched_drugs"].astype(bool)
    df["inferred_influence"] = df["log2FoldChange"].apply(infer_drug_influence)

    Path("data").mkdir(exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\n✓ saved {OUT}  ({len(df):,} rows)")
    print(f"  matched in DrugAge: {df['drugage_annotated'].sum():,} "
          f"({100*df['drugage_annotated'].mean():.1f}%)")


if __name__ == "__main__":
    main()
