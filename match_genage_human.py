#!/usr/bin/env python3
"""
Match human significant genetic interventions against the GenAge Human gene list
and write data/agent_result_genage_hs_{MODALITY}.csv for the waffle-plot step.

GenAge Human has no per-gene Pro/Anti-longevity direction column — only a
curated list of ~300 human genes implicated in aging. The human waffle is
therefore a 2-bucket partition:
  - annotated in GenAge (blue)
  - not in GenAge      (gray)

Flip MODALITY and rerun for each modality.
"""

# ---- CONFIG --------------------------------------------------------
MODALITY = "rna"      # "rna" or "dnam"
METHOD_NAME = f"hs_{MODALITY}"
GENAGE_HUMAN  = "data/refdata/genage_human.csv"
GENAGE_MODELS = "data/refdata/genage_models.csv"     # has longevity_influence across species
OUT = f"data/agent_result_genage_hs_{MODALITY}.csv"
# Priority order when resolving ambiguous multi-species annotations
SPECIES_PRIORITY = ["Mus musculus", "Homo sapiens", "Rattus norvegicus",
                    "Drosophila melanogaster", "Caenorhabditis elegans",
                    "Saccharomyces cerevisiae"]
# --------------------------------------------------------------------

import re
from pathlib import Path

import pandas as pd


def extract_gene_symbol(condition_name: str, intervention: str, gene_set: set) -> str:
    """Pull a gene symbol out of condition_name / intervention text.

    Strategy: tokenise both strings, strip suffixes like '-KO', '+/+', '(miRNA)',
    then look up each token in the GenAge symbol set (case-insensitive match
    with the canonical symbol preserved). Returns the matched symbol or ''.
    """
    suffixes = [
        "knockout", "knockdown", "overexpression", "deletion", "transgenic",
        "mutant", "mutation", "variant", "deficiency", "deficient",
        "ko", "kd", "oe", "-/-", "+/+", "+/-", "hom", "het",
    ]
    gene_set_upper = {g.upper() for g in gene_set}
    for text in (condition_name or "", intervention or ""):
        text = re.sub(r"[()\[\]{}]", " ", str(text))
        # Remove suffixes
        for suf in suffixes:
            text = re.sub(rf"\b{re.escape(suf)}\b", " ", text, flags=re.IGNORECASE)
        # Tokenise on non-alnum (keeping intragenic hyphens like "BRCA1-2")
        tokens = re.findall(r"[A-Za-z][A-Za-z0-9\-]*[A-Za-z0-9]", text)
        for tok in tokens:
            T = tok.upper()
            # Require >=3 chars: 2-letter tokens (TO, ON, OR) collide with English
            # words and short cross-species symbols and are too ambiguous.
            if T in gene_set_upper and len(T) >= 3 and not T.isdigit():
                # Biology abbreviations / English words that collide with gene
                # symbols but are virtually never the intended gene in this text
                # (e.g. ESC = embryonic stem cell, not the Drosophila esc gene).
                if T in {"AND", "FOR", "THE", "OFF", "ESC", "ESCS", "IPSC",
                         "IPSCS", "IPS", "NSC", "HSC", "MSC", "EPC"}:
                    continue
                return T
    return ""


def infer_influence(condition_name: str, intervention: str, log2fc) -> str:
    """Guess Pro/Anti-Longevity for the target gene, given intervention type + direction.

    KO / KD + log2FC > 0 (loss drives aging)   → gene is Pro-Longevity
    KO / KD + log2FC < 0 (loss slows aging)    → gene is Anti-Longevity
    OE       + log2FC > 0 (gain drives aging)  → gene is Anti-Longevity
    OE       + log2FC < 0 (gain slows aging)   → gene is Pro-Longevity
    """
    text = f"{condition_name or ''} {intervention or ''}".lower()
    is_loss = any(kw in text for kw in [
        "knockout", "knockdown", "-/-", "deletion", "deficient", "deficiency",
        " ko", " kd", "shrna", "sirna", "crispr", "mutation", "mutant", "null",
    ])
    is_gain = any(kw in text for kw in [
        "overexpression", " oe", "transgenic", " tg", "knockin", "gain",
    ])
    if not is_loss and not is_gain:
        is_loss = True  # default assumption when nothing matches
    if log2fc is None or pd.isna(log2fc):
        return ""
    if is_loss and log2fc > 0:   return "Pro-Longevity"
    if is_loss and log2fc < 0:   return "Anti-Longevity"
    if is_gain and log2fc > 0:   return "Anti-Longevity"
    if is_gain and log2fc < 0:   return "Pro-Longevity"
    return ""


def build_direction_lookup() -> dict:
    """symbol→longevity_influence map built from GenAge Models (Mus first)."""
    mdl = pd.read_csv(GENAGE_MODELS)
    mdl["symbol_u"] = mdl["symbol"].astype(str).str.upper()
    mdl["longevity_influence"] = mdl["longevity influence"].astype(str)
    # Assign a priority rank by species; best row per symbol wins
    rank = {s: i for i, s in enumerate(SPECIES_PRIORITY)}
    mdl["_rank"] = mdl["organism"].map(lambda o: rank.get(o, 99))
    mdl = mdl.sort_values(["symbol_u", "_rank"])
    dedup = mdl.drop_duplicates(subset=["symbol_u"], keep="first")
    return dict(zip(dedup["symbol_u"], dedup["longevity_influence"]))


def main():
    print(f"loading GenAge Human: {GENAGE_HUMAN}")
    ga = pd.read_csv(GENAGE_HUMAN)
    gene_set = set(ga["symbol"].dropna().astype(str).str.upper())
    print(f"  {len(gene_set):,} unique human genes listed in GenAge Human")

    direction_map = build_direction_lookup()
    print(f"  {len(direction_map):,} gene→longevity_influence entries from GenAge Models")
    # Union the lookup keys into gene_set so we can match conserved genes that
    # are in mouse but not in the human-specific list
    gene_set = gene_set | set(direction_map.keys())
    print(f"  combined reference: {len(gene_set):,} gene symbols")

    master = f"result/main_results_{METHOD_NAME}_with_scores.csv"
    print(f"loading {master}")
    df = pd.read_csv(master, low_memory=False)
    df = df[df["condition_category"] == "Genetic"].copy()
    print(f"  {len(df):,} Genetic scored rows")

    print("extracting gene symbols + inferring longevity influence...")
    df["genage_matched_genes"] = df.apply(
        lambda r: extract_gene_symbol(r.get("condition_name", ""),
                                       r.get("intervention", ""),
                                       gene_set),
        axis=1,
    )
    df["genage_annotated"] = df["genage_matched_genes"].astype(bool)
    # Direction from GenAge Models lookup (any species, Mus prioritised)
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
    print(f"  matched in GenAge Human: {df['genage_annotated'].sum():,} "
          f"({100*df['genage_annotated'].mean():.1f}%)")
    print(f"  inferred Pro-Longevity:  {(df['inferred_influence']=='Pro-Longevity').sum():,}")
    print(f"  inferred Anti-Longevity: {(df['inferred_influence']=='Anti-Longevity').sum():,}")

    matches = df[df["genage_annotated"]].head(10)
    if len(matches):
        print("\nexample GenAge matches:")
        for _, r in matches.iterrows():
            print(f"  {r['genage_matched_genes']:<10} ← {str(r.get('intervention',''))[:40]:<42}  log2FC={r['log2FoldChange']:+.2f}  inferred={r['inferred_influence']}")


if __name__ == "__main__":
    main()
