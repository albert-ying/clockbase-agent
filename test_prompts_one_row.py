#!/usr/bin/env python3
"""
Dry-run: take one row from each for_scoring_hs_*.csv, build the same ir_kw /
paper_kw dict that minipaper/src/run.py would build, and render both prompts.

Purpose is to catch missing placeholders or wrong column names BEFORE we send
anything through the LLM pipeline. GEO metadata + retrieved_documents are
stubbed with realistic-looking placeholders because we don't have the
external services running here.
"""

import re
import sys
from pathlib import Path

import pandas as pd
import requests

# Make the prompt module importable
sys.path.insert(0, str(Path(__file__).parent))
import human_prompts_for_minipaper as HP


def fetch_geo_metadata(gse_id: str) -> dict:
    """Mirror of get_data.get_geo_metadata() — direct GEO XML API, no proxy."""
    url = (
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        f"acc={gse_id}&targ=self&form=xml&view=quick"
    )
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    text = r.text

    def first(tag):
        m = re.search(rf"<{tag}>(.+?)</{tag}>", text, re.S)
        return m.group(1).strip() if m else ""

    return {
        "title": first("Title"),
        "summary": first("Summary"),
        "overall_design": first("Overall-Design"),
        "citations": [],  # PubMed abstracts not needed for prompt-format smoke test
    }


# RAG output is still stubbed — that's a separate service, not fetchable here
MOCK_RETRIEVED_DOCUMENTS = "<RAG RETRIEVAL RESULTS — fetched from the RAG endpoint at runtime>"
MOCK_INTERVENTION_SOURCE_LITERATURE = "<PubMed abstracts from GEO citations — fetched at runtime>"


def build_ir_kw_rna(row, meta):
    """Mirror the kwargs run.py passes into IDENTIFY_RELEVANCE.format()."""
    return {
        "gse_id": row["gse_id"],
        "title": meta["title"],
        "summary": meta["summary"],
        "study_design": meta["overall_design"],
        "condition": row["condition_name"],
        "comparison": row["comparison_name"],
        "condition_cate": row["condition_category"],
        "effect_size": row["effect_size"],
        "effect_direction": row["direction"],
        "effect_pvalue": row["pvalue"],
        "extreme_value": row["outlier_status"],
        "retrieved_documents": MOCK_RETRIEVED_DOCUMENTS,
        "intervention_source_literature": MOCK_INTERVENTION_SOURCE_LITERATURE,
    }


def build_ir_kw_dnam(row, meta):
    """Same as RNA plus the 9 cross-clock fields for the multi-clock prompt."""
    kw = build_ir_kw_rna(row, meta)
    kw.update(
        {
            "phenoage_effect_size": row.get("PhenoAge_effect_size"),
            "phenoage_pvalue": row.get("PhenoAge_pvalue"),
            "horvathv1_effect_size": row.get("Horvathv1_effect_size"),
            "horvathv1_pvalue": row.get("Horvathv1_pvalue"),
            "yingdamage_effect_size": row.get("YingDamAge_effect_size"),
            "yingdamage_pvalue": row.get("YingDamAge_pvalue"),
            "dunedinpace_effect_size": row.get("DunedinPACE_effect_size"),
            "dunedinpace_pvalue": row.get("DunedinPACE_pvalue"),
            "n_support": row.get("n_support"),
        }
    )
    return kw


def build_paper_kw_rna(row, meta):
    return {
        k: build_ir_kw_rna(row, meta)[k]
        for k in ["gse_id", "title", "condition", "comparison", "effect_size", "effect_direction", "effect_pvalue"]
    }


def build_paper_kw_dnam(row, meta):
    kw = build_paper_kw_rna(row, meta)
    kw.update(
        {
            "phenoage_effect_size": row.get("PhenoAge_effect_size"),
            "phenoage_pvalue": row.get("PhenoAge_pvalue"),
            "horvathv1_effect_size": row.get("Horvathv1_effect_size"),
            "horvathv1_pvalue": row.get("Horvathv1_pvalue"),
            "yingdamage_effect_size": row.get("YingDamAge_effect_size"),
            "yingdamage_pvalue": row.get("YingDamAge_pvalue"),
            "dunedinpace_effect_size": row.get("DunedinPACE_effect_size"),
            "dunedinpace_pvalue": row.get("DunedinPACE_pvalue"),
            "n_support": row.get("n_support"),
        }
    )
    return kw


def pick_representative_row(df):
    """Pick a top-effect-size row to exercise the prompt meaningfully."""
    return df.nsmallest(1, "log2FoldChange").iloc[0]


def dryrun(csv_path, modality):
    print("=" * 80)
    print(f"  {modality.upper()}:  {csv_path}")
    print("=" * 80)
    df = pd.read_csv(csv_path, low_memory=False)
    print(f"  Loaded {len(df):,} rows, columns: {list(df.columns)[:8]}...")
    row = pick_representative_row(df)
    print(f"\n  Test row:")
    print(f"    case_id           = {row['case_id']}")
    print(f"    gse_id            = {row['gse_id']}")
    print(f"    condition_name    = {row['condition_name']}")
    print(f"    condition_category= {row['condition_category']}")
    print(f"    effect_size       = {row['effect_size']}")
    print(f"    direction         = {row['direction']}")
    print(f"    pvalue            = {row['pvalue']}")
    print(f"    fdr               = {row['fdr']}")
    print(f"    outlier_status    = {row['outlier_status']}")
    if modality == "dnam":
        print(f"    n_support         = {row['n_support']}")
        for clk in ["PhenoAge", "Horvathv1", "YingDamAge", "DunedinPACE"]:
            print(
                f"    {clk:<13}   effect_size={row[f'{clk}_effect_size']}, pvalue={row[f'{clk}_pvalue']}"
            )

    print(f"\n  Fetching GEO metadata for {row['gse_id']}...")
    try:
        meta = fetch_geo_metadata(row["gse_id"])
        print(f"    ✓ title           ({len(meta['title']):,} chars)")
        print(f"    ✓ summary         ({len(meta['summary']):,} chars)")
        print(f"    ✓ overall_design  ({len(meta['overall_design']):,} chars)")
    except Exception as e:
        print(f"    ✗ GEO fetch failed: {e}")
        return

    if modality == "rna":
        ir_kw = build_ir_kw_rna(row, meta)
        paper_kw = build_paper_kw_rna(row, meta)
        ir_template = HP.IDENTIFY_RELEVANCE_HUMAN_RNA
        paper_template = HP.GENERATE_PAPER_2_HUMAN_RNA
    else:
        ir_kw = build_ir_kw_dnam(row, meta)
        paper_kw = build_paper_kw_dnam(row, meta)
        ir_template = HP.IDENTIFY_RELEVANCE_HUMAN_DNAM
        paper_template = HP.GENERATE_PAPER_2_HUMAN_DNAM

    print(f"\n  --- Formatting IDENTIFY_RELEVANCE ---")
    try:
        ir_prompt = ir_template.format(**ir_kw)
        print(f"  ✓ ir_prompt rendered ({len(ir_prompt):,} chars)")
    except KeyError as e:
        print(f"  ✗ MISSING KEY in ir_kw: {e}")
        return

    print(f"\n  --- Formatting GENERATE_PAPER_2 ---")
    try:
        paper_prompt = paper_template.format(**paper_kw)
        print(f"  ✓ paper_prompt rendered ({len(paper_prompt):,} chars)")
    except KeyError as e:
        print(f"  ✗ MISSING KEY in paper_kw: {e}")
        return

    # Write to tmp for visual inspection
    out_ir = Path(f"/tmp/test_ir_{modality}.txt")
    out_paper = Path(f"/tmp/test_paper_{modality}.txt")
    out_ir.write_text(ir_prompt)
    out_paper.write_text(paper_prompt)
    print(f"\n  Wrote {out_ir} and {out_paper} for inspection.")

    # Show the opening recap block (first 40 lines)
    print(f"\n  --- IDENTIFY_RELEVANCE opening (first 40 lines of {len(ir_prompt.splitlines()):,}) ---")
    for line in ir_prompt.splitlines()[:40]:
        print(f"    {line}")
    print(f"\n  --- GENERATE_PAPER_2 opening (first 20 lines of {len(paper_prompt.splitlines()):,}) ---")
    for line in paper_prompt.splitlines()[:20]:
        print(f"    {line}")
    print()


if __name__ == "__main__":
    dryrun("result/for_scoring_hs_rna.csv", "rna")
    dryrun("result/for_scoring_hs_dnam.csv", "dnam")
