#!/usr/bin/env python3
"""
Quality audit of manuscript-agent drafts.

For each draft in data/report*/<case_id>/parsed_report.md scores:
  (1) format compliance: Introduction / Results / Discussion / References
  (2) summary-statistic presence: p-value/FDR token, effect-size/log2FC token
  (3) word count
  (4) reference / citation token count

Outputs revision_r1_2026-05-18/manuscript_agent_quality_S11.csv with per-draft
scores and prints the summary statistics.

Run from the Aging-AI-Agent directory:  python3 audit_manuscript_agent_quality.py
"""
import glob
import random
import re
import statistics as st
from pathlib import Path

import pandas as pd

OUTDIR = Path("../revision_r1_2026-05-18")
SAMPLE_N = 100
SEED = 42

SECTIONS = {
    "Introduction": r"(?i)^#+\s*introduction",
    "Results": r"(?i)^#+\s*results",
    "Discussion": r"(?i)^#+\s*discussion",
    "References": r"(?i)^#+\s*references",
}


def score_draft(path: str) -> dict:
    txt = Path(path).read_text(errors="ignore")
    sections = {k: bool(re.search(p, txt, re.M)) for k, p in SECTIONS.items()}
    return {
        "path": path,
        "case_id": Path(path).parent.name,
        "word_count": len(re.findall(r"\b\w+\b", txt)),
        **sections,
        "all_sections": all(sections.values()),
        "has_pvalue": bool(re.search(r"p[-\s]?[=<]\s*\d|FDR\s*[<=]\s*\d", txt, re.I)),
        "has_effect_size": bool(re.search(r"(effect\s*size|log2[Ff]old|β\s*=)", txt)),
        "unique_gse_mentions": len(set(re.findall(r"GSE\d+", txt))),
        "citation_tokens": len(re.findall(r"(?:PMC\d+|PMID:?\s*\d+|doi:\s*\S+|https?://\S+)", txt)),
    }


def main():
    random.seed(SEED)
    all_files = sorted(glob.glob("data/report*/*/parsed_report.md"))
    sample = random.sample(all_files, min(SAMPLE_N, len(all_files)))
    print(f"available drafts: {len(all_files)} ; sampling {len(sample)} (seed {SEED})")

    rows = [score_draft(f) for f in sample]
    df = pd.DataFrame(rows)
    OUTDIR.mkdir(exist_ok=True)
    out_csv = OUTDIR / "manuscript_agent_quality_S11.csv"
    df.to_csv(out_csv, index=False)
    print(f"saved {out_csv}  ({len(df)} drafts)")

    n = len(df)
    print(f"\n=== summary (random sample of {n}, seed {SEED}) ===")
    for sec in SECTIONS:
        c = int(df[sec].sum())
        print(f"  has {sec}: {c}/{n} ({100*c/n:.0f}%)")
    print(f"  all four sections: {int(df['all_sections'].sum())}/{n}")
    print(f"  p-value/FDR present: {int(df['has_pvalue'].sum())}/{n} "
          f"({100*df['has_pvalue'].mean():.0f}%)")
    print(f"  effect-size present: {int(df['has_effect_size'].sum())}/{n} "
          f"({100*df['has_effect_size'].mean():.0f}%)")
    wc = df["word_count"].tolist()
    qs = st.quantiles(wc, n=4)
    print(f"  word count: median {int(st.median(wc))}, IQR [{int(qs[0])}-{int(qs[2])}]")
    refs = df["citation_tokens"].tolist()
    print(f"  citation tokens: median {int(st.median(refs))} (range {min(refs)}-{max(refs)})")
    print(f"  drafts with >=5 citations: {int((df['citation_tokens']>=5).sum())}/{n}")


if __name__ == "__main__":
    main()
