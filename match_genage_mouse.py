#!/usr/bin/env python3
"""
Match mouse RNA-seq Genetic interventions against GenAge and write
data/agent_result_genage.csv for the combined waffle-plot step.

This is the mouse counterpart of match_genage_human.py. It exists to close a
reproducibility gap: data/agent_result_genage.csv previously had no builder
script in the repository. It reuses the (bug-fixed) gene-extraction and
longevity-influence logic from match_genage_human.py and reads the canonical
outlier-removed mouse analytic table (detailed_scores_result_mm_rna.csv).
"""

from pathlib import Path

import pandas as pd

from match_genage_human import (
    extract_gene_symbol,
    infer_influence,
    build_direction_lookup,
    GENAGE_HUMAN,
)

INPUT = "result/detailed_scores_result_mm_rna.csv"
OUT = "data/agent_result_genage.csv"


def main():
    print(f"loading GenAge reference: {GENAGE_HUMAN}")
    ga = pd.read_csv(GENAGE_HUMAN)
    gene_set = set(ga["symbol"].dropna().astype(str).str.upper())
    direction_map = build_direction_lookup()
    gene_set = gene_set | set(direction_map.keys())
    print(f"  combined reference: {len(gene_set):,} gene symbols")

    print(f"loading {INPUT}")
    df = pd.read_csv(INPUT, low_memory=False)
    df = df[df["condition_category"] == "Genetic"].copy()
    print(f"  {len(df):,} Genetic rows (canonical outlier-removed mouse set)")
    # Manuscript GenAge set is FDR<0.05 (the score floor is applied downstream
    # by the waffle script via trust_score).
    df = df[df["fdr"] < 0.05].copy()
    print(f"  {len(df):,} FDR-significant Genetic rows")

    df["genage_matched_genes"] = df.apply(
        lambda r: extract_gene_symbol(r.get("condition_name", ""),
                                      r.get("intervention", ""), gene_set),
        axis=1,
    )
    df["genage_annotated"] = df["genage_matched_genes"].astype(bool)
    df["genage_annotation"] = df["genage_matched_genes"].apply(
        lambda g: direction_map.get(g.upper() if isinstance(g, str) else "", "")
    )
    df["inferred_influence"] = df.apply(
        lambda r: infer_influence(r.get("condition_name", ""),
                                  r.get("intervention", ""),
                                  r.get("log2FoldChange")),
        axis=1,
    )

    Path("data").mkdir(exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\n✓ saved {OUT}  ({len(df):,} rows)")
    print(f"  matched in GenAge: {df['genage_annotated'].sum():,} "
          f"({100*df['genage_annotated'].mean():.1f}%)")


if __name__ == "__main__":
    main()
