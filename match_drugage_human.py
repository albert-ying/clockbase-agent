#!/usr/bin/env python3
"""
Match human significant Drug interventions against the DrugAge compound database
and write data/agent_result_drugage_hs_{MODALITY}.csv for the waffle-plot step.

DrugAge entries carry avg_lifespan_change_percent which yields a per-compound
longevity direction:
  positive → Pro-Longevity
  negative → Anti-Longevity
So the human drug waffle supports full 4-bucket agreement analysis.

Flip MODALITY and rerun for each modality.
"""

# ---- CONFIG --------------------------------------------------------
MODALITY = "rna"       # "rna" or "dnam"
METHOD_NAME = f"hs_{MODALITY}"
DRUGAGE_CSV = "data/refdata/drugage.csv"
OUT = f"data/agent_result_drugage_hs_{MODALITY}.csv"
# --------------------------------------------------------------------

import re
from pathlib import Path

import pandas as pd
import numpy as np


def normalise_drug(name: str) -> str:
    """Lowercase + alnum-only for fuzzy matching."""
    if not isinstance(name, str):
        return ""
    return re.sub(r"[^a-z0-9]", "", name.lower())


def build_drugage_index(drugage_df: pd.DataFrame) -> dict:
    """Collapse DrugAge rows to one entry per compound with majority-direction label.

    Returns:  { canonical_name_lower: (display_name, 'Pro-Longevity'|'Anti-Longevity'|'Unclear') }
    """
    # Keep only compounds with a numeric lifespan-change value
    d = drugage_df.copy()
    d["avg_lifespan_change_percent"] = pd.to_numeric(
        d["avg_lifespan_change_percent"], errors="coerce"
    )
    d = d.dropna(subset=["compound_name", "avg_lifespan_change_percent"])
    # Majority vote per compound
    out = {}
    for name, g in d.groupby("compound_name"):
        pos = (g["avg_lifespan_change_percent"] > 0).sum()
        neg = (g["avg_lifespan_change_percent"] < 0).sum()
        if pos > neg:
            label = "Pro-Longevity"
        elif neg > pos:
            label = "Anti-Longevity"
        else:
            label = "Unclear"
        out[normalise_drug(name)] = (name, label)
    return out


def _text_ngrams(text: str, max_n: int = 4) -> set:
    """Word-boundary n-grams of the text, each normalised to alnum-only.

    The text is split into words on any non-alphanumeric character, so a
    compound name only matches when it lines up with whole words. This avoids
    spurious substring hits such as 'cytidine' inside '5-aza-2-deoxycytidine'
    while still catching multi-word compounds ('nicotinamide riboside').
    """
    words = [w for w in re.split(r"[^A-Za-z0-9]+", text.lower()) if w]
    grams = set()
    for n in range(1, max_n + 1):
        for i in range(len(words) - n + 1):
            grams.add("".join(words[i:i + n]))
    return grams


def match_drug(text: str, da_index: dict):
    """Find a DrugAge entry whose canonical name appears as whole word(s) in text."""
    if not isinstance(text, str) or not text.strip():
        return "", ""
    grams = _text_ngrams(text)
    if not grams:
        return "", ""
    # Word-boundary match; prefer the longest matched compound name.
    best_hit = None
    for key, (name, label) in da_index.items():
        if len(key) < 4:
            continue  # skip very short names (ambiguous)
        if key in grams:
            if best_hit is None or len(key) > len(best_hit[0]):
                best_hit = (key, name, label)
    if best_hit:
        return best_hit[1], best_hit[2]
    return "", ""


def infer_drug_influence(log2fc) -> str:
    """Pro/Anti based on aging-clock direction: anti-aging = Pro-Longevity."""
    if log2fc is None or pd.isna(log2fc):
        return ""
    return "Pro-Longevity" if log2fc < 0 else "Anti-Longevity"


def main():
    print(f"loading DrugAge: {DRUGAGE_CSV}")
    da = pd.read_csv(DRUGAGE_CSV)
    da_index = build_drugage_index(da)
    print(f"  {len(da_index):,} unique compounds in DrugAge index")

    master = f"result/main_results_{METHOD_NAME}_with_scores.csv"
    print(f"loading {master}")
    df = pd.read_csv(master, low_memory=False)
    df = df[df["condition_category"] == "Drug"].copy()
    print(f"  {len(df):,} Drug scored rows")

    print("matching against DrugAge...")
    matches = df.apply(
        lambda r: match_drug(
            f"{r.get('intervention','')}  {r.get('condition_name','')}",
            da_index,
        ),
        axis=1,
    )
    df["drugage_matched_drugs"] = matches.apply(lambda x: x[0])
    df["drugage_annotation"]    = matches.apply(lambda x: x[1])
    df["drugage_annotated"]     = df["drugage_matched_drugs"].astype(bool)
    df["inferred_influence"]    = df["log2FoldChange"].apply(infer_drug_influence)

    Path("data").mkdir(exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"\n✓ saved {OUT}  ({len(df):,} rows)")
    print(f"  matched in DrugAge: {df['drugage_annotated'].sum():,} "
          f"({100*df['drugage_annotated'].mean():.1f}%)")
    ann = df[df["drugage_annotated"]].copy()
    ann["agreement"] = ann["drugage_annotation"] == ann["inferred_influence"]
    print(f"  of matched:  agreement={ann['agreement'].sum():,}  "
          f"disagreement={(~ann['agreement']).sum():,}")

    # Show a few matched examples
    print("\nexample DrugAge matches:")
    for _, r in ann.head(10).iterrows():
        print(f"  {r['drugage_matched_drugs']:<30} ← {str(r.get('intervention',''))[:35]:<37}  "
              f"DrugAge={r['drugage_annotation']:<14}  inferred={r['inferred_influence']}")


if __name__ == "__main__":
    main()
