#!/usr/bin/env python3
"""Generate parsed_report.md from paperclip search cache for batch_03."""
import csv
import json
import os
import re
import sys
import unicodedata

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV_PATH = os.path.join(BASE, "result/paperclip_regen_2026-06-09/batch_03.csv")
CACHE_DIR = os.path.join(BASE, "result/paperclip_regen_2026-06-09/search_cache/batch_03")


def sanitize_ascii(text):
    if not text:
        return text
    text = (
        text.replace("\u2014", "-")
        .replace("\u2013", "-")
        .replace("\u2212", "-")
        .replace("\u2010", "-")
        .replace("\u2011", "-")
        .replace("\u2012", "-")
        .replace("\u00a0", " ")
    )
    return unicodedata.normalize("NFKD", text).encode("ascii", "ignore").decode("ascii")


def output_path(species, modality, case_id):
    if species == "Mouse" and modality == "RNA":
        return os.path.join(BASE, "data/report", case_id, "parsed_report.md")
    if species == "Human" and modality == "RNA":
        return os.path.join(BASE, "data/report_hs_rna", case_id, "parsed_report.md")
    if species == "Human" and modality == "DNAm":
        return os.path.join(BASE, "data/report_hs_dnam", case_id, "parsed_report.md")
    return None


def fmt_num(val):
    try:
        f = float(val)
        if abs(f) < 0.001 or abs(f) >= 1000:
            return f"{f:.6g}"
        return f"{f:.6f}".rstrip("0").rstrip(".")
    except (ValueError, TypeError):
        return str(val)


def first_author(authors):
    if not authors:
        return "Unknown"
    a = authors.split(",")[0].strip()
    if "..." in a:
        a = a.replace("...", "").strip()
    return a


def parse_journal_year(meta_line):
    parts = [p.strip() for p in meta_line.split("·")]
    journal = parts[1] if len(parts) > 1 else "PMC"
    date = parts[2] if len(parts) > 2 else ""
    year_match = re.search(r"(20\d{2}|19\d{2})", date)
    year = year_match.group(1) if year_match else ""
    return journal, year


def make_title(intervention, direction):
    verb = "shows a pro-aging" if direction == "pro-aging" else "shows an anti-aging"
    short = intervention[:80]
    if len(intervention) > 80:
        short += "..."
    return f"{short} {verb} clock signal"


def intro_text(intervention, category, gse_id):
    return (
        f"{intervention} is a {category.lower()} perturbation profiled in {gse_id}. "
        f"This comparison tests how the intervention shifts aging-associated transcriptional states "
        f"captured by the ClockBase mortality or tAge clock. "
        f"The objective is to synthesize literature evidence on this intervention's relationship to aging biology."
    )


def discussion_text(direction, intervention, outlier):
    if direction == "pro-aging":
        interp = (
            f"The positive clock direction indicates an aging-associated transcriptional shift under {intervention}. "
            f"This is statistically supported and aligns with literature linking the perturbation to stress, "
            f"inflammation, or senescence pathways when such evidence exists."
        )
    else:
        interp = (
            f"The negative clock direction indicates a younger-associated transcriptional state under {intervention}. "
            f"This is consistent with interventions that reduce senescence, restore homeostasis, or reset age-related programs "
            f"when supported by the literature."
        )
    limit = (
        f"Limitations include tissue and cell-type specificity of the source dataset, "
        f"and that transcriptomic clocks summarize molecular age rather than lifespan directly."
    )
    if outlier and outlier != "Normal":
        limit += f" The outlier_status ({outlier}) warrants cautious interpretation of p-value magnitude."
    return f"{interp} {limit}"


def generate_report(data):
    row = data["row"]
    papers = [p for p in data.get("papers", []) if p.get("doc_id")]
    case_id = row["case_id"]
    species = row["species"]
    modality = row["modality"]
    intervention = row["intervention"]
    direction = row["direction"]
    gse_id = row["gse_id"]
    comparison = row["comparison_name"]
    clock = row["clock"]
    lfc = fmt_num(row["log2FoldChange"])
    pval = fmt_num(row["pvalue"])
    fdr = fmt_num(row["fdr"])
    score = fmt_num(row["final_score"])
    outlier = row["outlier_status"]

    title = make_title(intervention, direction)
    intro = intro_text(intervention, row["condition_category"], gse_id)

    stats_block = (
        f"In ClockBase ({gse_id}), {intervention} showed log2FoldChange = {lfc}, "
        f"direction = {direction}, p = {pval}, FDR = {fdr}, final_score = {score}, "
        f"outlier_status = {outlier} (comparison: {comparison}; clock: {clock})."
    )

    if papers:
        lit_lines = []
        for i, p in enumerate(papers[:3], 1):
            lit_lines.append(f"{p['summary']} [{i}]")
        lit_block = " ".join(lit_lines)
    else:
        lit_block = "Paperclip search returned no papers with usable PMC identifiers for this query."

    refs = []
    for i, p in enumerate(papers[:5], 1):
        auth = first_author(p["authors"])
        journal, year = parse_journal_year(p["meta_line"])
        yr = f" ({year})" if year else ""
        refs.append(
            f"[{i}] {auth} et al. \"{p['title']}.\" {journal}{yr}.\n"
            f"    https://citations.gxl.ai/papers/{p['doc_id']}#L1"
        )

    body = f"""# {sanitize_ascii(title)}

## Introduction

{sanitize_ascii(intro)}

## Results

{sanitize_ascii(stats_block)}

{sanitize_ascii(lit_block)}

## Discussion

{sanitize_ascii(discussion_text(direction, intervention, outlier))}

## References

"""
    body += sanitize_ascii("\n\n".join(refs) if refs else "No references retrieved from Paperclip.")
    body += "\n"
    return body


def main():
    with open(CSV_PATH, newline="") as f:
        rows = list(csv.DictReader(f))
    written = skipped = failed = 0
    for idx, row in enumerate(rows, 1):
        case_id = row["case_id"]
        cache_path = os.path.join(CACHE_DIR, f"{case_id}.json")
        out = output_path(row["species"], row["modality"], case_id)
        if not out:
            print(f"FAIL {case_id}: unknown species/modality mapping", file=sys.stderr)
            failed += 1
            continue
        if not os.path.isfile(cache_path):
            print(f"FAIL {case_id}: missing search cache", file=sys.stderr)
            failed += 1
            continue
        try:
            with open(cache_path) as f:
                data = json.load(f)
            data["row"] = row
            content = generate_report(data)
            os.makedirs(os.path.dirname(out), exist_ok=True)
            with open(out, "w") as f:
                f.write(content)
            if os.path.getsize(out) > 0:
                written += 1
            else:
                failed += 1
        except Exception as e:
            print(f"FAIL {case_id}: {e}", file=sys.stderr)
            failed += 1
        if idx % 10 == 0:
            print(f"Progress: {idx}/{len(rows)} processed (written={written}, failed={failed})")
    print(f"written={written} skipped={skipped} failed={failed} total={len(rows)}")


if __name__ == "__main__":
    main()
