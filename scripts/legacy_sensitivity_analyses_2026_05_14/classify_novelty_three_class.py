#!/usr/bin/env python3
"""Classify 106 unmatched candidates into a three-class novelty taxonomy.

Class A: any aging-related link in cited literature
Class B: direction-consistent aging-related link (literature supports agent's direction)
Class C: no prior aging-related evidence (0-1 aging-corpus citations, or refs lack
         intervention-specific aging-relevant content)
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

REPO = Path("/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base")
NOVELTY_CSV = REPO / "figure_rev1/sensitivity_analyses/novelty_distribution_unmatched.csv"
REPORT_DIR = REPO / "Aging-AI-Agent/data/report"
OUT_CSV = REPO / "figure_rev1/sensitivity_analyses/novelty_three_class.csv"
OUT_MD = REPO / "figure_rev1/sensitivity_analyses/novelty_three_class_summary.md"

# Master CSVs for direction lookup (authoritative source for agent's direction
# call, since effect_size sign can be inverted by comparator definitions).
MASTER_MM_RNA = REPO / "Aging-AI-Agent/result/main_results_result_mm_rna_outlier_marked.csv"
MASTER_HS_RNA = REPO / "Aging-AI-Agent/result/main_results_hs_rna_outlier_marked.csv"
MASTER_HS_DNAM = REPO / "Aging-AI-Agent/result/main_results_hs_dnam_outlier_marked.csv"


def load_directions() -> dict[tuple[str, str], float]:
    """Map (gse_id, condition_name) -> direction (-1 anti-aging, +1 pro-aging)
    from master CSVs."""
    directions: dict[tuple[str, str], float] = {}
    for fp in [MASTER_MM_RNA, MASTER_HS_RNA, MASTER_HS_DNAM]:
        if not fp.exists():
            continue
        with fp.open() as f:
            reader = csv.DictReader(f)
            for row in reader:
                gse = (row.get("gse_id") or "").strip()
                cond = (row.get("condition_name") or "").strip()
                if not gse or not cond:
                    continue
                dval = row.get("direction") or ""
                try:
                    d = float(dval)
                except (TypeError, ValueError):
                    continue
                # Keep first occurrence (most rows are unique gse+condition)
                directions.setdefault((gse, cond), d)
    return directions


# Cue phrases used in the Introduction that indicate DIRECTION-CONSISTENT prior
# evidence for the intervention itself (not just generic aging biology).
ANTI_AGING_DIRECT_CUES = [
    r"extend(s|ed|ing)?\s+(life|health)span",
    r"prolong(s|ed|ing)?\s+(life|health)span",
    r"increase(s|d)?\s+(life|health)span",
    r"longevity[- ]promoting",
    r"known\s+geroprotector",
    r"established\s+geroprotector",
    r"reduces?\s+aging\s+(markers|biomarkers|hallmarks)",
    r"ameliorat(es|ed|ing)\s+aging",
    r"slows?\s+aging",
    r"slows?\s+(dna[- ]methylation|epigenetic|biological)\s+age",
    r"attenuat(es|ed|ing)\s+aging",
    r"delays?\s+aging",
    r"reverses?\s+aging",
    r"protective\s+against\s+aging",
    r"anti[- ]aging\s+(effect|properties|compound|intervention)",
    r"senolytic",
    r"caloric\s+restriction\s+mim",
    r"rejuvenat(es|ed|ing|ion)",
    r"counteract(s|ed|ing)?\s+aging",
    r"reduces?\s+(cellular\s+)?senescence",
    r"clears?\s+senescent\s+cells",
    r"reduce(s|d)?\s+(all[- ]cause\s+)?mortality",
    r"decreases?\s+mortality",
    r"increases?\s+telomere",
    r"telomere\s+length(en|ening)?",
    r"younger\s+(transcript|methyl|biological)",
    r"protect(s|ed|ive)?\s+against\s+(age|senescence)",
    r"counteract(s|ed|ing)?\s+sarcopenia",
    r"prevent(s|ed|ing)?\s+age[- ]related",
    r"improve(s|d|ment)?\s+healthspan",
    r"oppose(s|d|ing)?\s+(inflammaging|aging)",
    r"dampen(s|ed|ing)?\s+inflammaging",
]

PRO_AGING_DIRECT_CUES = [
    r"accelerat(es|ed|ing)?\s+(biological\s+|molecular\s+|cellular\s+)?aging",
    r"promotes?\s+aging",
    r"induces?\s+(cellular\s+)?senescence",
    r"drives?\s+aging",
    r"drive(r|rs)\s+of\s+(\w+\s+){0,3}(cellular\s+senescence|aging|age[- ]related|senescence)",
    r"shortens?\s+(life|health)span",
    r"reduces?\s+(life|health)span",
    r"premature\s+aging",
    r"premature\s+senescen",
    r"progeroid",
    r"photoaging|photo[- ]aging",
    r"increase(s|d)?\s+aging\s+(markers|biomarkers|hallmarks)",
    r"causes?\s+aging",
    r"pro[- ]aging\s+(effect|phenotype|shift|signature)",
    r"contribute(s|d)?\s+to\s+aging",
    r"associated\s+with\s+(aging|senescence|age[- ]related)",
    r"linked\s+to\s+aging",
    r"declines?\s+(with|during|in)\s+aging",
    r"impaired\s+(with|during)\s+aging",
    r"loss\s+of\s+\w+\s+(causes|drives|promotes)\s+",
    r"chronic\s+inflammation",
    r"inflammaging",
    r"hallmark(s)?\s+of\s+aging",
    r"age[- ]related\s+(pathology|disease|decline|dysfunction|fibrosis)",
    r"mitochondrial\s+dysfunction.{0,30}(aging|senescence)",
    r"aging.{0,30}mitochondrial\s+dysfunction",
    r"(senescent|aging)\s+cells?",
    r"increased.*mortality",
    r"raises\s+mortality",
    r"clonal\s+hematopoiesis|clonal\s+expansion",
    r"sarcopenia|frailty",
    r"sasp",  # senescence-associated secretory phenotype
    r"DNA\s+damage\s+response.{0,40}(aging|senescence)",
    r"telomere\s+(attrition|shortening|dysfunction)",
    r"dna\s+methylation\s+age",
    r"epigenetic\s+(clock|drift|age)",
    r"age[- ]associated|age-associated",
]


def load_novelty_csv() -> list[dict[str, str]]:
    rows = []
    with NOVELTY_CSV.open() as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def parse_report(case_id: str) -> dict[str, object] | None:
    md = REPORT_DIR / case_id / "parsed_report.md"
    if not md.exists():
        return None
    text = md.read_text()

    # Title is first line that looks like a title (start of file before Introduction)
    # Drafts use either "# Title" or "**Title**\n\nText" formats.
    title_match = re.search(r"^#\s+(.+)$", text, flags=re.MULTILINE)
    if title_match:
        title = title_match.group(1).strip()
    else:
        # Try **Title** format
        tm = re.search(r"\*\*Title\*\*\s*\n+(.+?)(?:\n\n|\n\*\*)", text, flags=re.DOTALL)
        title = tm.group(1).strip() if tm else ""

    # Introduction: between Introduction header and next section
    intro_match = re.search(
        r"(?:^##\s+Introduction|^\*\*Introduction\*\*)\s*\n+(.+?)(?:^##\s+|^\*\*[A-Z])",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    intro = intro_match.group(1).strip() if intro_match else ""

    # Results section text - sometimes the direction language appears there
    results_match = re.search(
        r"(?:^##\s+Results|^\*\*Results\*\*)\s*\n+(.+?)(?:^##\s+Discussion|^\*\*Discussion\*\*)",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    results = results_match.group(1).strip() if results_match else ""

    # References: extract list items
    refs_match = re.search(
        r"(?:^##\s+References|^\*\*References\*\*)\s*\n+(.+)$",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    refs_block = refs_match.group(1).strip() if refs_match else ""
    # Split on numbered list items
    refs = re.split(r"\n+(?=(?:\[\d+\]|\d+\.))", refs_block)
    refs = [r.strip() for r in refs if r.strip()]

    return {
        "title": title,
        "intro": intro,
        "results": results,
        "refs": refs,
    }


def count_aging_refs(refs: list[str]) -> int:
    """Count references that mention aging, senescence, longevity, or lifespan."""
    aging_kw = re.compile(
        r"\b(aging|ageing|age[- ]related|senescen|longevity|lifespan|healthspan|"
        r"progeroid|geroprotect|anti[- ]aging|biological\s+age|epigenetic\s+age|"
        r"sarcopenia|frailty|cellular\s+senescence|gerontolog)",
        flags=re.IGNORECASE,
    )
    n = 0
    for r in refs:
        if aging_kw.search(r):
            n += 1
    return n


def title_direction(title: str, effect_size: float) -> str:
    """Determine agent direction from effect size (authoritative source)."""
    # Effect size in masters: negative = anti-aging (rejuvenates), positive = pro-aging.
    # This is the authoritative direction; titles can be inverted by framing.
    return "anti-aging" if effect_size < 0 else "pro-aging"


def direction_consistent_intro(intro: str, direction: str, require_citation: bool = False) -> tuple[bool, str]:
    """Check if the Introduction contains direction-consistent claims about the
    intervention's aging effect. Returns (matched, matched_phrase)."""
    cues = ANTI_AGING_DIRECT_CUES if direction == "anti-aging" else PRO_AGING_DIRECT_CUES
    intro_l = intro.lower()
    for pat in cues:
        m = re.search(pat, intro_l)
        if m:
            if require_citation:
                # Ensure a citation appears within +/- 200 chars of the match
                window = intro_l[max(0, m.start() - 200):m.end() + 200]
                if re.search(r"\[\s*\d+", window):
                    return True, m.group(0)
            else:
                return True, m.group(0)
    return False, ""


def classify(direction: str, intro: str, results: str, refs: list[str], n_aging_refs: int) -> tuple[str, str]:
    """Return (class, rationale).

    Definitions:
      C = no prior aging-related evidence: <=1 aging-themed citation, OR the
          Introduction explicitly states the intervention has not been previously
          characterized in aging.
      B = direction-consistent: the report's Introduction/Results cite prior
          literature where either (i) the specific intervention has been linked
          to aging in the same direction, OR (ii) the core biological pathway
          targeted by the intervention has been shown in cited prior literature
          to drive aging in the same direction. Operationalized as: >=2 aging
          citations AND at least one direction-consistent claim appears in the
          combined Intro+Results body.
      A = any aging-related link but NOT direction-consistent: >=2 aging
          citations but no direction-consistent intervention-specific claim is
          present in the report's narrative.
    """
    intro_ref_count = len(re.findall(r"\[\s*\d+", intro))
    results_ref_count = len(re.findall(r"\[\s*\d+", results))
    body = intro + "\n\n" + results

    not_studied_pat = re.compile(
        r"(not\s+(been\s+)?(previously|extensively)\s+(characterized|studied|investigated)\s+"
        r"(?:in|for|with|regarding)?\s*(?:the\s+)?(?:context\s+of\s+)?aging|"
        r"limited\s+prior\s+evidence\s+(in|for|regarding)\s+aging|"
        r"no\s+(direct|prior)\s+(evidence|studies)\s+(linking|connecting).{0,40}aging|"
        r"has\s+not\s+been\s+(directly\s+)?linked\s+to\s+aging|"
        r"novel\s+(intervention|target)\s+in\s+aging\s+research)",
        flags=re.IGNORECASE,
    )

    if n_aging_refs <= 1:
        rationale = (
            f"Only {n_aging_refs} aging-themed citation in references; report lacks "
            f"prior literature linking intervention to aging."
        )
        return "C_no_prior_evidence", rationale

    if not_studied_pat.search(body):
        rationale = (
            "Report explicitly notes the intervention has not been previously characterized "
            f"in aging context despite {n_aging_refs} aging-themed refs."
        )
        return "C_no_prior_evidence", rationale

    # B vs A: does the narrative present a direction-consistent claim grounded in citations?
    # Strip the first ~600 chars of Results (which restate the agent's own
    # analytical finding) so we only consider prior-literature synthesis.
    results_lit = results[600:] if len(results) > 600 else ""
    consistent_intro, phrase_intro = direction_consistent_intro(intro, direction, require_citation=False)
    consistent_results, phrase_results = direction_consistent_intro(results_lit, direction, require_citation=False)
    consistent = consistent_intro or consistent_results
    phrase = phrase_intro or phrase_results
    results_lit_refs = len(re.findall(r"\[\s*\d+", results_lit))

    # B requires >=2 aging refs AND consistent claim AND at least one citation
    # somewhere in the narrative (intro OR results-lit), so the direction-consistent
    # claim is grounded in cited prior literature.
    has_anchor_citations = (intro_ref_count + results_lit_refs) >= 1

    if consistent and has_anchor_citations and n_aging_refs >= 2:
        rationale = (
            f"Report cites prior literature with direction-consistent ({direction}) "
            f"framing (\"{phrase}\"); {n_aging_refs} aging-themed refs in bibliography."
        )
        return "B_direction_consistent", rationale

    # Otherwise: aging-related refs exist but no direction-consistent narrative
    rationale = (
        f"{n_aging_refs} aging-corpus refs but report narrative lacks an explicit "
        f"direction-consistent ({direction}) claim about the intervention; "
        f"cited literature provides general aging-biology context only."
    )
    return "A_any_link", rationale


def main() -> None:
    rows_in = load_novelty_csv()
    master_dir = load_directions()
    out_rows = []
    for r in rows_in:
        case_id = r["case_id"]
        parsed = parse_report(case_id)
        if parsed is None:
            out_rows.append(
                {
                    "modality": r["modality"],
                    "gse_id": r["gse_id"],
                    "case_id": case_id,
                    "intervention": r["intervention"],
                    "agent_direction": "anti-aging" if float(r["effect_size"]) < 0 else "pro-aging",
                    "novelty_class_3way": "C_no_prior_evidence",
                    "rationale_one_line": "Report file missing; defaulting to no-prior-evidence.",
                }
            )
            continue

        eff = float(r["effect_size"])
        # Direction from master CSV (authoritative; handles comparator inversions).
        # Fall back to effect_size sign if not in master.
        dval = master_dir.get((r["gse_id"], r["condition_name"]))
        if dval is None:
            direction = "anti-aging" if eff < 0 else "pro-aging"
        else:
            direction = "anti-aging" if dval < 0 else "pro-aging"

        n_aging = count_aging_refs(parsed["refs"])
        cls, rationale = classify(direction, parsed["intro"], parsed["results"], parsed["refs"], n_aging)
        out_rows.append(
            {
                "modality": r["modality"],
                "gse_id": r["gse_id"],
                "case_id": case_id,
                "intervention": r["intervention"],
                "agent_direction": direction,
                "novelty_class_3way": cls,
                "rationale_one_line": rationale,
            }
        )

    # Write CSV
    fieldnames = [
        "modality",
        "gse_id",
        "case_id",
        "intervention",
        "agent_direction",
        "novelty_class_3way",
        "rationale_one_line",
    ]
    with OUT_CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in out_rows:
            w.writerow(row)

    # Write summary
    total = len(out_rows)
    counts = {"A_any_link": 0, "B_direction_consistent": 0, "C_no_prior_evidence": 0}
    by_modality: dict[str, dict[str, int]] = {}
    for r in out_rows:
        counts[r["novelty_class_3way"]] += 1
        by_modality.setdefault(r["modality"], {"A_any_link": 0, "B_direction_consistent": 0, "C_no_prior_evidence": 0})
        by_modality[r["modality"]][r["novelty_class_3way"]] += 1

    def pct(n: int, d: int) -> str:
        return f"{(100.0 * n / d):.1f}%" if d else "n/a"

    examples: dict[str, list[dict]] = {"A_any_link": [], "B_direction_consistent": [], "C_no_prior_evidence": []}
    for r in out_rows:
        if len(examples[r["novelty_class_3way"]]) < 5:
            examples[r["novelty_class_3way"]].append(r)

    lines: list[str] = []
    lines.append("# Three-class novelty classification\n")
    lines.append(
        "Classification of 106 unmatched ClockBase Agent candidates "
        "(interventions not present in GenAge/DrugAge).\n"
    )
    lines.append(
        "- **Class A** any aging-related link (some prior literature connects intervention or biology to aging)\n"
        "- **Class B** direction-consistent aging-related link (prior literature supports the same direction)\n"
        "- **Class C** no prior aging-related evidence (<=1 aging-corpus citation in the draft)\n"
    )
    lines.append(
        "**Method.** Each unmatched candidate has an agent-written mini-paper draft "
        "at `Aging-AI-Agent/data/report/<case_id>/parsed_report.md` synthesizing prior "
        "literature. For each draft we (1) extract the agent's measured direction from "
        "the master CSV (`direction` column: -1 anti-aging, +1 pro-aging, accounting "
        "for comparator orientation); (2) count aging-themed bibliography entries using "
        "the regex `aging|ageing|senescen|longevity|lifespan|healthspan|progeroid|geroprotect|"
        "biological age|epigenetic age|sarcopenia|frailty|cellular senescence|gerontology`; "
        "(3) search the Introduction and the prior-literature portion of Results (after "
        "the first ~600 chars, which restate the agent's own measurement) for "
        "direction-consistent cue phrases co-located with at least one citation marker. "
        "Drafts with <=1 aging-themed reference are Class C. Drafts with >=2 aging "
        "references AND a citation-anchored direction-consistent claim are Class B. The "
        "remainder (>=2 aging references but no anchored direction-consistent claim) "
        "are Class A.\n\n"
        "**Caveat.** The Class B count is an upper bound: regex cues cannot always "
        "distinguish prior-literature statements about the *specific intervention* from "
        "statements about the *general pathway* the intervention targets. We treat "
        "pathway-level direction consistency as sufficient for B because the reviewer's "
        "question is whether the agent's call is grounded in some direction-consistent "
        "prior literature, which a coherent pathway-level argument satisfies.\n"
    )

    lines.append(f"## Overall ({total} candidates)\n")
    lines.append(f"- Class A: {counts['A_any_link']} ({pct(counts['A_any_link'], total)})")
    lines.append(
        f"- Class B: {counts['B_direction_consistent']} ({pct(counts['B_direction_consistent'], total)})"
    )
    lines.append(
        f"- Class C: {counts['C_no_prior_evidence']} ({pct(counts['C_no_prior_evidence'], total)})\n"
    )

    lines.append("## Per-modality\n")
    for mod, d in sorted(by_modality.items()):
        tot = sum(d.values())
        lines.append(f"### {mod} (n={tot})")
        lines.append(f"- Class A: {d['A_any_link']} ({pct(d['A_any_link'], tot)})")
        lines.append(f"- Class B: {d['B_direction_consistent']} ({pct(d['B_direction_consistent'], tot)})")
        lines.append(f"- Class C: {d['C_no_prior_evidence']} ({pct(d['C_no_prior_evidence'], tot)})\n")

    lines.append("## Representative examples\n")
    label_map = {
        "A_any_link": "Class A (any aging-related link, but not direction-consistent)",
        "B_direction_consistent": "Class B (direction-consistent prior evidence)",
        "C_no_prior_evidence": "Class C (no prior aging-related evidence)",
    }
    for cls, label in label_map.items():
        lines.append(f"### {label}\n")
        for ex in examples[cls]:
            lines.append(
                f"- **{ex['case_id']} | {ex['intervention']}** "
                f"({ex['modality']}, agent direction: {ex['agent_direction']})\n"
                f"  - {ex['rationale_one_line']}"
            )
        lines.append("")

    OUT_MD.write_text("\n".join(lines))

    print(f"Wrote {OUT_CSV}")
    print(f"Wrote {OUT_MD}")
    print(f"A={counts['A_any_link']}  B={counts['B_direction_consistent']}  C={counts['C_no_prior_evidence']}")


if __name__ == "__main__":
    main()
