#!/usr/bin/env python3
"""
Human-version score summarizer.
Walks minipaper_fork/output{_hs_dnam}/GSE*-i/identified_relevance.json and
produces data/scores_summary_hs_{modality}.csv with one row per scored case.
Flip MODALITY at the top and rerun for each modality.
"""

# ---- HUMAN CONFIG --------------------------------------------------
MODALITY = "rna"  # "rna" or "dnam"
# Where the scorer wrote its outputs. For RNA the orchestrator left output/
# in place after the run ended; for DNAm it was renamed to output_hs_dnam/.
OUTPUT_DIRS = {
    "rna":  "minipaper_fork/output",
    "dnam": "minipaper_fork/output_hs_dnam",
}
# --------------------------------------------------------------------

import json
from pathlib import Path

import pandas as pd


def extract_scores_from_json(json_path):
    try:
        with open(json_path) as f:
            data = json.load(f)
        scores = {}
        for key, value in data.items():
            if isinstance(value, dict) and "points" in value:
                scores[key] = value["points"]
            elif key == "q0_plausibility" and isinstance(value, dict) and "penalty_or_bonus_applied" in value:
                scores["q0_plausibility"] = value["penalty_or_bonus_applied"]
        if "final_score" in data:
            scores["final_score"] = data["final_score"]
        if "q6_intervention" in data:
            scores["intervention"] = data["q6_intervention"]
        return scores
    except Exception as e:
        print(f"  error reading {json_path}: {e}")
        return None


def main():
    output_dir = Path(OUTPUT_DIRS[MODALITY])
    if not output_dir.exists():
        raise SystemExit(f"output dir does not exist: {output_dir}")
    print(f"walking {output_dir} for identified_relevance.json files...")

    rows = []
    for sub in sorted(output_dir.iterdir()):
        if not sub.is_dir():
            continue
        j = sub / "identified_relevance.json"
        if not j.exists():
            continue
        s = extract_scores_from_json(j)
        if s is None:
            continue
        # dataset name is the case_id, e.g. GSE12345-67
        s["dataset"] = sub.name
        rows.append(s)

    df = pd.DataFrame(rows)
    print(f"  collected {len(df)} scored cases")

    # Stable column order: dataset, q*, intervention, final_score
    q_cols = sorted([c for c in df.columns if c.startswith("q")])
    order = ["dataset"] + q_cols
    if "intervention" in df.columns:
        order.append("intervention")
    if "final_score" in df.columns:
        order.append("final_score")
    order = [c for c in order if c in df.columns]
    df = df[order].sort_values("dataset").reset_index(drop=True)

    # Derive gse_id + row_index from dataset string for easy merging later
    df["gse_id"] = df["dataset"].str.rsplit("-", n=1).str[0]
    df["row_index"] = df["dataset"].str.rsplit("-", n=1).str[1]

    out_path = Path(f"data/scores_summary_hs_{MODALITY}.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"\n✓ saved {out_path}  ({len(df)} rows, {len(df.columns)} cols)")

    if "final_score" in df.columns:
        print(f"\nfinal_score stats: mean={df['final_score'].mean():.1f}, "
              f"median={df['final_score'].median():.1f}, "
              f"min={df['final_score'].min()}, max={df['final_score'].max()}")
    print(f"unique GSE studies: {df['gse_id'].nunique()}")
    if "q6_intervention" in df.columns or "intervention" in df.columns:
        top = df.nlargest(10, "final_score")[["dataset", "intervention", "final_score"]]
        print("\ntop 10 by final_score:")
        print(top.to_string(index=False))


if __name__ == "__main__":
    main()
