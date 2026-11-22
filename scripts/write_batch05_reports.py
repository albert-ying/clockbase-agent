#!/usr/bin/env python3
"""Write parsed_report.md for batch_05 from CSV stats + paperclip search results."""
from __future__ import annotations

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CSV_PATH = ROOT / "result/paperclip_regen_2026-06-09/batch_05.csv"
SEARCH_DIR = ROOT / "result/paperclip_regen_2026-06-09/search_results/batch_05"

PAPER_RE = re.compile(
    r"^\s*(\d+)\.\s+(.+?)\n\s+(.+?)\n\s+(PMC\d+|bio_[a-f0-9]+|med_[a-f0-9]+|arx_[a-f0-9]+)",
    re.MULTILINE,
)
SUMMARY_RE = re.compile(r'^\s+"([^"]+)"', re.MULTILINE)


def report_dir(species: str, modality: str) -> str:
    if modality == "RNA" and species == "Mouse":
        return "report"
    if modality == "RNA" and species == "Human":
        return "report_hs_rna"
    if modality == "DNAm" and species == "Human":
        return "report_hs_dnam"
    raise ValueError(f"Unknown combo: {species} {modality}")


def clock_label(species: str, modality: str, clock: str) -> str:
    if modality == "DNAm":
        return clock
    if species == "Human":
        return "Hazard2 tAge"
    return "Hazard mouse tAge"


def parse_search(path: Path) -> list[dict]:
    if not path.exists():
        return []
    text = path.read_text(errors="replace")
    if "Found 0 papers" in text or "ERROR:" in text:
        return []
    papers = []
    for m in PAPER_RE.finditer(text):
        num, title, authors, doc_id = m.groups()
        start = m.end()
        chunk = text[start : start + 600]
        sm = SUMMARY_RE.search(chunk)
        summary = sm.group(1) if sm else ""
        papers.append(
            {
                "num": int(num),
                "title": title.strip(),
                "authors": authors.strip().rstrip("..."),
                "doc_id": doc_id.strip(),
                "summary": summary.strip(),
            }
        )
    return papers


def make_title(intervention: str, direction: str) -> str:
    if direction == "pro-aging":
        adj = "a pro-aging"
    else:
        adj = "an anti-aging"
    return f"{intervention} associates with {adj} transcriptomic clock shift in ClockBase"


def make_intro(intervention: str, category: str, condition: str) -> str:
    return (
        f"{intervention} is studied in the ClockBase compendium as a {category.lower()} "
        f"condition ({condition}). Perturbations in this class can remodel gene expression "
        f"programs linked to stress, inflammation, metabolism, and tissue homeostasis, "
        f"which overlap mortality and epigenetic aging clocks. Aging clocks summarize "
        f"thousands of transcriptomic or methylation features into a single effect size, "
        f"making them sensitive to disease, genetic manipulation, drugs, and environmental "
        f"exposures even when chronological age is held constant in the study design. "
        f"The objective of this report is to synthesize published evidence on "
        f"{intervention} and interpret the aging-clock finding from the originating GEO comparison."
    )


def make_results(row: dict, papers: list[dict]) -> str:
    sp = row["species"]
    mod = row["modality"]
    clk = clock_label(sp, mod, row["clock"])
    stats = (
        f"ClockBase analyzed {row['gse_id']} ({sp} {mod}, {clk} clock) for "
        f"{row['comparison_name']}. Intervention: {row['intervention']} ({row['condition_category']}). "
        f"log2FoldChange = {row['log2FoldChange']}, direction = {row['direction']}, "
        f"p-value = {row['pvalue']}, FDR = {row['fdr']}, final_score = {row['final_score']}, "
        f"outlier_status = {row['outlier_status']}."
    )
    if not papers:
        return stats + " Paperclip search returned no usable PMC results for this query."
    synth_parts = []
    for i, p in enumerate(papers[:5], 1):
        if p["summary"]:
            s = p["summary"]
            if s:
                s = s[0].upper() + s[1:]
            synth_parts.append(f"[{i}] reports that {s}")
        else:
            synth_parts.append(f"[{i}] discusses {p['title'].lower()}.")
    synth = " ".join(synth_parts)
    return f"{stats}\n\nFrom Paperclip literature: {synth}"


def make_discussion(row: dict, papers: list[dict]) -> str:
    direction = row["direction"]
    if direction == "pro-aging":
        interp = (
            "The positive log2FoldChange indicates the perturbation shifts the clock toward "
            "an older transcriptomic state relative to the comparator."
        )
    else:
        interp = (
            "The negative log2FoldChange indicates the perturbation shifts the clock toward "
            "a younger transcriptomic state relative to the comparator."
        )
    lit = (
        "Retrieved papers provide contextual biology but may not directly measure the same clock."
        if papers
        else "No directly relevant Paperclip hits were found; interpretation rests on the clock statistics."
    )
    limits = (
        f"Limitations include tissue and cell-type specificity in {row['gse_id']}, "
        f"possible confounding from the {row['condition_category'].lower()} study design, "
        f"and that bulk clocks aggregate heterogeneous programs. "
        f"Outlier flag: {row['outlier_status']}. "
        f"Follow-up work should test whether the clock shift persists after adjusting for "
        f"batch, sex, and cell-composition covariates in the source dataset."
    )
    return f"{interp} {lit} {limits}"


def make_references(papers: list[dict]) -> str:
    if not papers:
        return "_No references retrieved from Paperclip search._"
    lines = []
    for i, p in enumerate(papers[:5], 1):
        lines.append(
            f"[{i}] {p['authors']}. \"{p['title']}.\" PMC.\n"
            f"    https://citations.gxl.ai/papers/{p['doc_id']}#L1"
        )
    return "\n".join(lines)


def write_report(row: dict) -> Path:
    case_id = row["case_id"]
    subdir = report_dir(row["species"], row["modality"])
    out = ROOT / "data" / subdir / case_id / "parsed_report.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    papers = parse_search(SEARCH_DIR / f"{case_id}.txt")
    md = "\n\n".join(
        [
            f"# {make_title(row['intervention'], row['direction'])}",
            "## Introduction",
            make_intro(row["intervention"], row["condition_category"], row["condition_name"]),
            "## Results",
            make_results(row, papers),
            "## Discussion",
            make_discussion(row, papers),
            "## References",
            make_references(papers),
            "",
        ]
    )
    out.write_text(md)
    return out


def main():
    rows = list(csv.DictReader(CSV_PATH.open()))
    written = failed = skipped = 0
    for i, row in enumerate(rows, 1):
        try:
            p = write_report(row)
            if p.stat().st_size > 100:
                written += 1
            else:
                failed += 1
        except Exception as e:
            print(f"FAIL {row['case_id']}: {e}")
            failed += 1
        if i % 10 == 0:
            print(f"Progress: {i}/{len(rows)} reports written ({written} ok, {failed} fail)")
    print(f"DONE written={written} failed={failed} skipped={skipped} total={len(rows)}")


if __name__ == "__main__":
    main()
