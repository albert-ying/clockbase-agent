#!/usr/bin/env python3
"""Regenerate manuscript Tables 1 (pro-aging) and 2 (anti-aging) from the CURRENT
combined scores, using the documented composite metric:

  qualifying: FDR<0.05 AND final_score >= modality floor (mouse RNA 100,
              human RNA 100, human DNAm 80)
  scaled    : (final_score - min)/(max - min), min-max across the whole
              qualifying subset (global; matches Methods + Figure 3 bar plot)
  composite : -log10(pvalue) * scaled
  dedup     : by (gse_id, intervention), keep max composite
  headline  : intervention categories Genetic, Drug, Environment only
              (Disease/Other excluded); extreme-tail rows retained.

Writes result/table1_pro_aging_pooled_top20.csv, table2_..., and prints the
markdown tables (9 columns) for pasting into the manuscript.
"""
import numpy as np
import pandas as pd

FLOOR = {"mouse RNA-seq": 100, "human RNA-seq": 100, "human DNAm": 80}
MODSHORT = {"mouse RNA-seq": "muRNA", "human RNA-seq": "huRNA", "human DNAm": "huDNAm"}
# Cap on -log10(p) in the composite so extreme-significance outliers (e.g. sorted
# immune-cell comparisons with p < 1e-100) do not dominate the ranking. Matches the
# bar plot (create_combined_modality_top_hits_v2.py).
CAP_NEGLOG10P = 30


def sci(x):
    if pd.isna(x):
        return "NA"
    if x == 0:
        return "0"
    e = int(np.floor(np.log10(abs(x))))
    m = x / 10**e
    if round(m, 2) >= 10:  # rounding pushed mantissa to 10.00
        m /= 10
        e += 1
    return f"{m:.2f}×10^{e}^"


def main():
    d = pd.read_csv("figures/combined_3modality_plot_data.csv", low_memory=False)
    d = d[d["outlier_status"].notna()].copy()
    d["floor"] = d["modality"].map(FLOOR)
    q = d[(d["fdr"] < 0.05) & (d["final_score"] >= d["floor"])].copy()

    # min-max normalise the final score across the whole qualifying subset
    # (global, matching the Methods definition and the Figure 3 bar plot, which
    # both scale "within the qualifying subset" rather than per modality).
    fs_min = q["final_score"].min()
    fs_max = q["final_score"].max()
    span = fs_max - fs_min
    q["scaled"] = (q["final_score"] - fs_min) / span if span > 0 else 1.0
    q["composite"] = np.minimum(-np.log10(q["pvalue"].clip(lower=1e-300)), CAP_NEGLOG10P) * q["scaled"]
    q = q.sort_values("composite", ascending=False).drop_duplicates(["gse_id", "intervention"])

    TOP_N = 30
    gde = q[q["condition_category"].isin(["Genetic", "Drug", "Environment"])].copy()
    t1 = gde[gde["log2FoldChange"] > 0].nlargest(TOP_N, "composite")
    t2 = gde[gde["log2FoldChange"] < 0].nlargest(TOP_N, "composite")
    t1.to_csv(f"result/table1_pro_aging_pooled_top{TOP_N}.csv", index=False)
    t2.to_csv(f"result/table2_anti_aging_pooled_top{TOP_N}.csv", index=False)

    def md(df):
        rows = ["| GSE ID | Intervention | Category | Modality | β | FDR | Score | P-value | Composite |",
                "|---|---|---|---|---:|---:|---:|---:|---:|"]
        for _, r in df.iterrows():
            b = f"{r.log2FoldChange:+.2f}"
            rows.append(f"| {r['gse_id']} | {r['intervention']} | {r['condition_category']} | "
                        f"{r['modality']} | {b} | {sci(r['fdr'])} | "
                        f"{int(r['final_score'])} | {sci(r['pvalue'])} | {r['composite']:.2f} |")
        return "\n".join(rows)

    print("### TABLE 1 (pro-aging)\n"); print(md(t1))
    print("\n### TABLE 2 (anti-aging)\n"); print(md(t2))


if __name__ == "__main__":
    main()
