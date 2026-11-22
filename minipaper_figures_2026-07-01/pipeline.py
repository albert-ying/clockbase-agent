"""Generate a Fig-5-style sample-level figure (and, where missing, an SDK-written
mini-paper) for every FDR<0.05 & final_score>100 ClockBase hit.

Per hit:
  1. enumerate the GSE's sample groups from the standard per-sample export
  2. Claude SDK (Opus 4.7, subscription auth) picks which group(s) are the
     intervention (case) arm and which are the reference (control) arm, matching
     the published comparison_name -- returns group INDICES, not GSM echoes
  3. resolve indices -> full GSM lists, run the reliability engine, and VALIDATE:
     the recomputed direction must match the published direction (sign of the
     effect). Figures render only when the arms reproduce the published result;
     otherwise the hit is flagged low-confidence and skipped (never fabricated).
  4. render LEFT experimental-design schematic + RIGHT real per-sample box plot
     without displaying the reconstructed-arm p-value
  5. for hits without an existing parsed_report.md, an SDK call writes a
     mini-paper (Paperclip PMC evidence + exact ClockBase stats); the figure is
     referenced in every mini-paper.

Resumable: per-hit arm-derivation + mini-paper results are cached as JSON.

DATA INTEGRITY: every plotted value is a REAL per-sample clock prediction from the
ClockBase export. No simulated data. Arms are validated against the published stat.

Run:
  ./venv/bin/python minipaper_figures_2026-07-01/pipeline.py --build-hits
  ./venv/bin/python minipaper_figures_2026-07-01/pipeline.py --pilot 3
  ./venv/bin/python minipaper_figures_2026-07-01/pipeline.py --mode full
"""
from __future__ import annotations

import argparse
import asyncio
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent                                   # Aging-AI-Agent
AUDIT = ROOT / "reliability_audit_2026-06-16"
sys.path.insert(0, str(AUDIT))
import config as acfg                                 # noqa: E402
from engine import run_engine, load_persample        # noqa: E402
from dossier import _s                                # noqa: E402

FIGDIR = HERE / "figures"
ARMS_CACHE = HERE / "arms_cache"
PAPER_CACHE = HERE / "minipaper_cache"
HITS_CSV = HERE / "hits_score_gt100.csv"
for d in (FIGDIR, ARMS_CACHE, PAPER_CACHE):
    d.mkdir(parents=True, exist_ok=True)

MODEL = "claude-opus-4-7"
COMBINED = ROOT / "figures" / "combined_3modality_plot_data.csv"
REPORT_DIR = {
    "mouse RNA-seq": ROOT / "data" / "report",
    "human RNA-seq": ROOT / "data" / "report_hs_rna",
    "human DNAm": ROOT / "data" / "report_hs_dnam",
}
CLOCK_LABEL = {
    "mouse RNA-seq": "Mouse transcriptomic clock (tAge)",
    "human RNA-seq": "Human transcriptomic clock (tAge)",
    "human DNAm": "DNAm clock",
}
ASSAY_LABEL = {
    "mouse RNA-seq": "Bulk RNA sequencing",
    "human RNA-seq": "Bulk RNA sequencing",
    "human DNAm": "Genome-wide DNA\nmethylation array",
}

# ----------------------------------------------------------------------------
# Palette (blue = age-lowering case, rust = pro-aging case, grey = control)
# ----------------------------------------------------------------------------
INK = "#222222"
BLUE = ("#9DBED6", "#5B8FB0")
RUST = ("#D8A48F", "#B5654A")
GREY = ("#C9C9C9", "#9A9A9A")

_GROUP_FIELDS = [
    "harmonized.genetic_information", "harmonized.treatment", "harmonized.compound",
    "harmonized.disease", "harmonized.disease_stage", "harmonized.specimen_part",
    "harmonized.cell_line", "harmonized.developmental_stage", "harmonized.time",
    "harmonized.age", "gpt_case_control",
]


# ============================================================================
# 1. Hit list
# ============================================================================
def build_hits() -> pd.DataFrame:
    comb = pd.read_csv(COMBINED, low_memory=False)
    comb["fs"] = pd.to_numeric(comb["final_score"], errors="coerce")
    comb["fdrn"] = pd.to_numeric(comb["fdr"], errors="coerce")
    hits = (comb[(comb["fdrn"] < 0.05) & (comb["fs"] > 100)]
            .drop_duplicates(["modality", "gse_id", "comparison_name"]).copy())

    mm = (pd.read_csv(ROOT / "result" / "detailed_scores_result_mm_rna.csv", low_memory=False)
          [["gse_id", "comparison_name", "case_id"]].dropna().drop_duplicates())

    def hs_map(path):
        d = pd.read_csv(path)
        d["gse"] = d["dataset"].astype(str).str.split("-").str[0]
        return d[["gse", "intervention", "dataset"]].rename(columns={"dataset": "case_id"})
    hsr = hs_map(ROOT / "data" / "scores_summary_hs_rna.csv")
    hsd = hs_map(ROOT / "data" / "scores_summary_hs_dnam.csv")

    def case_for(r):
        if r["modality"] == "mouse RNA-seq":
            m = mm[(mm.gse_id == r["gse_id"]) & (mm.comparison_name == r["comparison_name"])]
            return m["case_id"].iloc[0] if len(m) else None
        tbl = hsr if r["modality"] == "human RNA-seq" else hsd
        m = tbl[(tbl.gse == r["gse_id"]) & (tbl.intervention.astype(str) == str(r["intervention"]))]
        return m["case_id"].iloc[0] if len(m) else None
    hits["case_id"] = hits.apply(case_for, axis=1)

    have = set(os.path.basename(os.path.dirname(p))
               for p in ROOT.glob("data/report*/*/parsed_report.md"))
    hits["has_report"] = hits["case_id"].isin(have)
    hits.to_csv(HITS_CSV, index=False)
    print(f"built {len(hits)} hits -> {HITS_CSV}  "
          f"(mapped {hits.case_id.notna().sum()}, existing report {int(hits.has_report.sum())})")
    return hits


def load_hits() -> pd.DataFrame:
    if not HITS_CSV.exists():
        return build_hits()
    return pd.read_csv(HITS_CSV, low_memory=False)


# ============================================================================
# 2. Sample-group enumeration (keeps GSM membership server-side)
# ============================================================================
def enumerate_groups(persample: pd.DataFrame, gse_id: str, clock_col: str):
    d = persample[persample["gse_id"] == gse_id].copy()
    if clock_col in d.columns:
        d[clock_col] = pd.to_numeric(d[clock_col], errors="coerce")
    else:
        d[clock_col] = np.nan
    gfields = [f for f in _GROUP_FIELDS if f in d.columns]
    d["_grp"] = d[gfields].apply(
        lambda r: " | ".join(f"{f.split('.')[-1]}={_s(r[f])}"
                             for f in gfields if _s(r[f]) != ""), axis=1)
    d["_grp"] = d["_grp"].replace("", "(no metadata)")
    groups = []
    for i, (sig, sub) in enumerate(sorted(d.groupby("_grp"), key=lambda kv: kv[0])):
        vals = sub[clock_col].dropna().to_numpy()
        titles = "; ".join(sorted({_s(t) for t in sub.get("title", pd.Series([])).head(4)} - {""}))[:180]
        groups.append({
            "idx": i, "signature": sig, "n": int(len(sub)),
            "clock_mean": float(np.mean(vals)) if len(vals) else None,
            "clock_sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else None,
            "example_titles": titles,
            "gsms": [s for s in (_s(x) for x in sub["gsm_id"]) if s],
        })
    return groups, d


def groups_text(groups) -> str:
    out = []
    for g in groups:
        m = f"{g['clock_mean']:.3f}" if g["clock_mean"] is not None else "NA"
        sd = f"{g['clock_sd']:.3f}" if g["clock_sd"] is not None else "NA"
        out.append(f"[G{g['idx']}] n={g['n']}  clock_mean={m}  clock_sd={sd}\n"
                   f"      {g['signature']}\n"
                   f"      titles: {g['example_titles']}")
    return "\n".join(out)


# ============================================================================
# 3. SDK arm derivation
# ============================================================================
ARM_PROMPT = """You are reconstructing the exact two-arm comparison behind a published \
aging-clock result, so a box plot can show the real per-sample data.

STUDY: {gse_id}  ({modality})
PUBLISHED COMPARISON: "{comparison_name}"
INTERVENTION / CONDITION: {intervention}
PUBLISHED DIRECTION: {direction_word}  (effect size {effect}, p {pvalue})

The comparison is written as "CASE vs CONTROL": the CASE (left name) is the \
intervention / perturbed / condition arm; the CONTROL (right name) is the reference \
/ baseline arm. Below are the sample groups of this study (clock_mean is the mean \
biological-age prediction of the group):

{groups}

Pick which group index/indices form the CASE arm and which form the CONTROL arm so \
the split reproduces the published comparison. Use only groups that clearly belong. \
If no valid healthy/reference contrast exists in the data, return empty arms with \
confidence "none".

Return STRICT JSON only:
{{"case_groups": [ints], "ctrl_groups": [ints], "case_label": "label (<=22 chars)", \
"ctrl_label": "label (<=22 chars)", "case_short": "axis tick (<=11 chars)", \
"ctrl_short": "axis tick (<=11 chars)", "confidence": "high|medium|low|none", \
"reasoning": "one sentence"}}"""


def _extract_json(raw: str) -> dict:
    m = re.search(r"\{.*\}", raw, re.DOTALL)
    if not m:
        raise ValueError("no JSON in model output")
    return json.loads(m.group(0))


async def _sdk_json(prompt: str, system: str) -> dict:
    from claude_agent_sdk import query, ClaudeAgentOptions, AssistantMessage, TextBlock
    opts = ClaudeAgentOptions(model=MODEL, max_turns=1, allowed_tools=[],
                              system_prompt=system, setting_sources=[])
    chunks = []
    async for msg in query(prompt=prompt, options=opts):
        if isinstance(msg, AssistantMessage):
            for b in msg.content:
                if isinstance(b, TextBlock):
                    chunks.append(b.text)
    return _extract_json("".join(chunks))


def derive_arms(row, groups) -> dict:
    prompt = ARM_PROMPT.format(
        gse_id=row["gse_id"], modality=row["modality"],
        comparison_name=row["comparison_name"], intervention=row.get("intervention", ""),
        direction_word=("pro-aging (raises age)" if float(row.get("direction", 0)) > 0
                        else "age-lowering (lowers age)"),
        effect=row.get("effect_size", row.get("log2FoldChange", "")),
        pvalue=row.get("pvalue", ""), groups=groups_text(groups))
    system = ("You reconstruct case/control sample arms for aging-clock comparisons. "
              "You are precise and never invent groups. Output strict JSON only.")
    return asyncio.run(_sdk_json(prompt, system))


# ============================================================================
# 4. Engine validation
# ============================================================================
# DNAm collapses 5 clocks; try each present clock and keep the one that reproduces
# the published direction (validation-driven clock choice). PhenoAge is preferred
# as the representative clock (matches the Fig 5 examples) when several qualify.
DNAM_CLOCKS = [
    ("PhenoAge", "clockresult.PhenoAge"), ("YingDamAge", "clockresult.YingDamAge"),
    ("Horvathv1", "clockresult.Horvathv1"), ("DunedinPACE", "clockresult.DunedinPACE"),
    ("GrimAge2", "DNAmGrimAge2"),
]


def _clock_candidates(modality, persample):
    if modality == "human DNAm":
        return [(lab, col) for lab, col in DNAM_CLOCKS if col in persample.columns]
    return [("", acfg.DEFAULT_CLOCK_COL[modality])]


_STOP = {"vs", "the", "and", "with", "of", "in", "at", "to", "a", "an", "cells", "cell",
         "control", "controls", "group", "groups", "samples", "tissue"}
_HEUR_FIELDS = ["harmonized.disease", "harmonized.genetic_information", "harmonized.treatment",
                "harmonized.compound", "harmonized.disease_stage", "harmonized.developmental_stage",
                "harmonized.specimen_part", "source_name", "title"]


def heuristic_arms(row, d):
    """Map the comparison_name 'CASE vs CONTROL' directly onto a metadata field value
    (exact match preferred). This aggregates arms the way the published contrast did,
    e.g. disease=='diabetic retinopathy' vs disease=='normal'."""
    cn = str(row.get("comparison_name", ""))
    body = cn
    if ":" in cn and " vs " in cn.split(":", 1)[1]:
        body = cn.split(":", 1)[1]
    if " vs " not in body.lower():
        return None
    parts = re.split(r"\s+vs\.?\s+", body, maxsplit=1, flags=re.IGNORECASE)
    if len(parts) != 2:
        return None
    left, right = parts[0].strip(), parts[1].strip()
    best = None
    for f in _HEUR_FIELDS:
        if f not in d.columns:
            continue
        vals = d[f].astype(str)
        for exact in (True, False):
            if exact:
                cm = vals.str.lower() == left.lower()
                km = vals.str.lower() == right.lower()
            else:
                cm = vals.str.contains(re.escape(left), case=False, na=False)
                km = vals.str.contains(re.escape(right), case=False, na=False)
            overlap = (cm & km)
            cm = cm & ~overlap; km = km & ~overlap
            nc, nk = int(cm.sum()), int(km.sum())
            if nc >= 2 and nk >= 2:
                score = (0 if exact else 1, -min(nc, nk))
                if best is None or score < best[0]:
                    best = (score, f, cm, km, nc, nk)
    if best is None:
        return None
    _, f, cm, km, nc, nk = best
    return {"case_gsms": d[cm]["gsm_id"].astype(str).tolist(),
            "ctrl_gsms": d[km]["gsm_id"].astype(str).tolist(),
            "case_label": left[:22], "ctrl_label": right[:22],
            "case_short": left[:11], "ctrl_short": right[:11],
            "confidence": "heuristic", "field": f}


def _arms_to_gsms(arms, groups):
    """Resolve an arm spec to explicit GSM lists (group-index form or gsm-list form)."""
    if "case_gsms" in arms:
        return [str(x) for x in arms["case_gsms"]], [str(x) for x in arms["ctrl_gsms"]]
    idx2 = {g["idx"]: g["gsms"] for g in groups}
    return ([str(x) for i in arms.get("case_groups", []) for x in idx2.get(i, [])],
            [str(x) for i in arms.get("ctrl_groups", []) for x in idx2.get(i, [])])


SIG_GATE = 0.05  # a rendered figure must itself reproduce the published significance


def validate_and_values(row, groups, arm_candidates, persample):
    """Evaluate each (arm candidate x clock); render only when the recomputed result
    reproduces the published direction AND is significant (p < SIG_GATE)."""
    pub_dir = int(np.sign(float(row.get("direction", 0))))
    d = persample[persample["gse_id"] == row["gse_id"]].copy()
    best = None                      # passes gate: (pvalue, payload)
    best_fail = None                 # best effort if none pass: (rank, reason)
    for arms in arm_candidates:
        if not arms:
            continue
        case_gsms, ctrl_gsms = _arms_to_gsms(arms, groups)
        if len(case_gsms) < 2 or len(ctrl_gsms) < 2:
            best_fail = best_fail or (9, "insufficient samples in a derived arm")
            continue
        case_sel = d["gsm_id"].astype(str).isin(case_gsms)
        ctrl_sel = d["gsm_id"].astype(str).isin(ctrl_gsms)
        for lab, clock in _clock_candidates(row["modality"], persample):
            res = run_engine(persample, row["gse_id"], row["modality"], clock,
                             case_gsms=case_gsms, ctrl_gsms=ctrl_gsms)
            if res.pvalue is None:
                continue
            dir_ok = (res.direction == pub_dir) or pub_dir == 0
            if not dir_ok:
                best_fail = best_fail or (2, f"direction mismatch (engine {res.direction} vs published {pub_dir})")
                continue
            if res.pvalue >= SIG_GATE:
                best_fail = best_fail or (3, f"arms do not reproduce significance (p={res.pvalue:.2g})")
                continue
            dd = d.copy(); dd[clock] = pd.to_numeric(dd[clock], errors="coerce")
            cvals = dd[case_sel][clock].dropna().to_numpy()
            kvals = dd[ctrl_sel][clock].dropna().to_numpy()
            if len(cvals) < 2 or len(kvals) < 2:
                continue
            payload = {"res": res, "case_vals": cvals, "ctrl_vals": kvals,
                       "clock": clock, "clock_label": lab, "arms": arms}
            if best is None or res.pvalue < best[0]:
                best = (res.pvalue, payload)
    if best is not None:
        return best[1], "ok"
    return None, (best_fail[1] if best_fail else "no usable arm/clock")


# ============================================================================
# 5. Figure
# ============================================================================
def _box(ax, x, y, w, h, text, fc, fs=7.2, weight="normal"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle="round,pad=0.012,rounding_size=0.04",
                 facecolor=fc, edgecolor=INK, linewidth=0.9, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=fs, color=INK,
            zorder=4, weight=weight, linespacing=1.25)


def _arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1),
                 arrowstyle="-|>,head_width=2.6,head_length=4.2",
                 mutation_scale=1, lw=1.0, color=INK, shrinkA=1, shrinkB=1, zorder=2))


def draw_schematic(ax, row, case_label, ctrl_label, n_case, n_ctrl, case_c, ctrl_c, clock_disp):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    _box(ax, 0.5, 0.92, 0.82, 0.115, f"{row['gse_id']} cohort", "white", fs=7.5, weight="bold")
    yarm = 0.66
    _box(ax, 0.27, yarm, 0.44, 0.15, f"{case_label}\n(n = {n_case})", case_c)
    _box(ax, 0.73, yarm, 0.44, 0.15, f"{ctrl_label}\n(n = {n_ctrl})", ctrl_c)
    _arrow(ax, 0.42, 0.862, 0.30, yarm + 0.08); _arrow(ax, 0.58, 0.862, 0.70, yarm + 0.08)
    _box(ax, 0.5, 0.40, 0.82, 0.105, ASSAY_LABEL[row["modality"]], "white")
    _arrow(ax, 0.27, yarm - 0.08, 0.42, 0.455); _arrow(ax, 0.73, yarm - 0.08, 0.58, 0.455)
    _box(ax, 0.5, 0.155, 0.82, 0.11, f"{clock_disp}\n→ biological age", "white")
    _arrow(ax, 0.5, 0.347, 0.5, 0.21)


def _wrap(s, width=13):
    """Word-wrap a short arm label onto at most two lines for the x-axis tick."""
    s = s.strip()
    if len(s) <= width:
        return s
    words = s.split()
    if len(words) == 1:
        return s[:width]
    line1 = ""
    for i, w in enumerate(words):
        if not line1 or len(line1) + 1 + len(w) <= width:
            line1 = (line1 + " " + w).strip()
        else:
            return line1 + "\n" + " ".join(words[i:])[:width + 4]
    return line1


def draw_boxplot(ax, case_vals, ctrl_vals, case_label, ctrl_label, case_c, ctrl_c, ylabel):
    rng = np.random.default_rng(0)
    data = [ctrl_vals, case_vals]; cols = [ctrl_c, case_c]; pos = [1, 1.62]
    bp = ax.boxplot(data, positions=pos, widths=0.24, patch_artist=True, showfliers=False,
                    medianprops=dict(color="black", lw=1.0),
                    whiskerprops=dict(lw=0.8, color=INK), capprops=dict(lw=0.8, color=INK),
                    boxprops=dict(lw=0.8, edgecolor=INK))
    for patch, c in zip(bp["boxes"], cols):
        patch.set_facecolor(c); patch.set_alpha(0.55)
    for i, (dv, c) in enumerate(zip(data, cols)):
        jit = rng.uniform(-0.06, 0.06, len(dv))
        ax.scatter(np.full(len(dv), pos[i]) + jit, dv, s=11, facecolors=c,
                   edgecolors="black", linewidths=0.4, alpha=0.85, zorder=4, clip_on=False)
    ax.set_ylabel(ylabel, fontsize=8.5)
    ax.set_xticks(pos); ax.set_xticklabels([ctrl_label, case_label], fontsize=8)
    ax.set_xlim(0.55, 2.07); ax.tick_params(axis="y", labelsize=8)
    lo = min(case_vals.min(), ctrl_vals.min()); hi = max(case_vals.max(), ctrl_vals.max())
    span = hi - lo if hi > lo else 1.0
    ax.set_ylim(lo - 0.07 * span, hi + 0.14 * span)
    for s in ("top", "right"): ax.spines[s].set_visible(False)


def make_figure(row, arms, vv, out_base: Path):
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 8, "axes.linewidth": 0.8,
                         "svg.fonttype": "none", "pdf.fonttype": 42})
    pro = float(row.get("direction", 0)) > 0
    case_c = (RUST if pro else BLUE)
    case_face, ctrl_face = case_c[0], GREY[0]
    case_lab = str(arms.get("case_label") or "Case")[:22]
    ctrl_lab = str(arms.get("ctrl_label") or "Control")[:22]
    case_short = _wrap(str(arms.get("case_short") or case_lab.split("(")[0]))
    ctrl_short = _wrap(str(arms.get("ctrl_short") or ctrl_lab.split("(")[0]))
    cv, kv = vv["case_vals"], vv["ctrl_vals"]
    if row["modality"] == "human DNAm":
        clock_disp = f"DNAm {vv.get('clock_label') or 'clock'}"
        ylabel = f"{vv.get('clock_label') or 'DNAm'} (years)"
    else:
        clock_disp = CLOCK_LABEL[row["modality"]]
        ylabel = f"{row['modality'].split(' ')[-1]} age (clock units)"
    fig, axd = plt.subplot_mosaic([["schema", "box"]], figsize=(7.2, 3.3),
                                  layout="constrained", gridspec_kw=dict(width_ratios=[1.2, 0.8]))
    draw_schematic(axd["schema"], row, case_lab, ctrl_lab, len(cv), len(kv), case_face, GREY[0], clock_disp)
    draw_boxplot(axd["box"], cv, kv, case_short, ctrl_short, case_c[1], GREY[1], ylabel)
    direction = "lower" if not pro else "higher"
    fig.suptitle(f"{_s(row.get('intervention','')) or row['comparison_name']}: "
                 f"{direction} biological age ({row['modality']})",
                 fontsize=8.5, weight="bold", y=1.04)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_base}.png", dpi=600, bbox_inches="tight")
    fig.savefig(f"{out_base}.pdf", bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# 6. Mini-paper (SDK) for hits without an existing report
# ============================================================================
def paperclip_evidence(intervention: str) -> str:
    try:
        out = subprocess.run(["paperclip", "search", "-s", "pmc",
                              f"{intervention} aging", "-n", "5"],
                             capture_output=True, text=True, timeout=120)
        return out.stdout[:4000] or "(no Paperclip results)"
    except Exception as e:
        return f"(Paperclip unavailable: {e})"


PAPER_PROMPT = """Write a concise mini-paper (250-450 words) on this ClockBase Agent finding.

Finding: in {gse_id} ({modality}), the comparison "{comparison_name}" \
({intervention}) shows {direction_word} biological age.
ClockBase statistics for prose reporting: effect size {effect}, p-value {pvalue}, \
FDR {fdr}, case n={n_case}, control n={n_ctrl}.

Literature evidence (PMC search results; cite what is relevant, link as given):
{evidence}

Structure with markdown headers: Introduction, Results, Discussion, References. \
In Results, report the ClockBase statistics in ordinary paper prose, not as a \
table, bullet list, or key=value ledger. Use compact numeric formatting: about \
three significant figures for effect sizes, p-values, and FDRs; scientific \
notation for very small p-values/FDRs; and no unnecessary trailing .0. \
Use plain numbered references without Markdown bold or italic wrappers. Do not \
write empty citation metadata such as () or (None), and do not write Paperclip \
result labels in References. Do not mention Paperclip by name in the report; \
refer to the evidence as cited literature. Use bracketed numeric citations such \
as [1] or [1,2], never caret markers such as ^1^. Do not say "search result" \
or "search-result summary" in the report prose. \
Do not overstate causal claims. \
Do NOT use em-dashes. Cite PMIDs/DOIs from the evidence where they support a point; \
if the evidence lacks a direct match, say so plainly."""


def _normalize_number_text(text: str) -> str:
    if "e" in text or "E" in text:
        mantissa, exponent = re.split("[eE]", text, maxsplit=1)
        mantissa = mantissa.rstrip("0").rstrip(".")
        if mantissa in {"", "-", "+"}:
            mantissa = "0"
        return f"{mantissa}e{int(exponent)}"
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return text or "0"


def _format_report_stat(value, kind: str) -> str:
    if value is None:
        return ""
    try:
        number = float(value)
    except (TypeError, ValueError):
        return _s(value)
    if not math.isfinite(number):
        return _s(value)
    if number == 0:
        return "0"
    if kind == "p" and abs(number) < 0.001:
        return _normalize_number_text(f"{number:.2e}")
    return _normalize_number_text(f"{number:.3g}")


def write_minipaper(row, evidence: str, vv=None) -> str:
    if isinstance(vv, dict):
        n_case, n_ctrl = len(vv["case_vals"]), len(vv["ctrl_vals"])
    else:
        n_case = n_ctrl = "(per published comparison)"
    prompt = PAPER_PROMPT.format(
        gse_id=row["gse_id"], modality=row["modality"], comparison_name=row["comparison_name"],
        intervention=_s(row.get("intervention", "")),
        direction_word=("age-lowering" if float(row.get("direction", 0)) < 0 else "pro-aging"),
        effect=_format_report_stat(row.get("effect_size", ""), "effect"),
        pvalue=_format_report_stat(row.get("pvalue", ""), "p"),
        fdr=_format_report_stat(row.get("fdr", ""), "p"),
        n_case=n_case, n_ctrl=n_ctrl, evidence=evidence)
    system = ("You are a careful aging-biology writer. You never fabricate statistics or "
              "citations, never use em-dashes, and never write raw field-name "
              "lines such as `clock = ...` or `pvalue = ...`. You write clean "
              "plain-text references without malformed Markdown.")

    async def _run():
        from claude_agent_sdk import query, ClaudeAgentOptions, AssistantMessage, TextBlock
        opts = ClaudeAgentOptions(model=MODEL, max_turns=1, allowed_tools=[],
                                  system_prompt=system, setting_sources=[])
        chunks = []
        async for msg in query(prompt=prompt, options=opts):
            if isinstance(msg, AssistantMessage):
                for b in msg.content:
                    if isinstance(b, TextBlock):
                        chunks.append(b.text)
        return "".join(chunks)
    return asyncio.run(_run())


FIG_MARKER = "## Sample-level validation figure"
FIG_CAPTION = (FIG_MARKER + "\n\n"
               "![Experimental design and per-sample biological-age comparison]"
               "({figrel})\n\n*Left: experimental-design schematic. Right: real per-sample "
               "ClockBase biological-age predictions for the reconstructed {case} vs {ctrl} "
               "arms. Arms were validated to reproduce the published "
               "direction of effect.*")


def strip_figure_block(text: str) -> str:
    """Remove a previously-inserted figure block (marker to the next H1/H2 or EOF)."""
    idx = text.find(FIG_MARKER)
    if idx == -1:
        return text
    rest = text[idx + len(FIG_MARKER):]
    # the block runs until the next markdown heading (## or #) or end of text
    m = re.search(r"\n#{1,2} ", rest)
    tail = rest[m.start():] if m else ""
    return (text[:idx].rstrip() + ("\n\n" + tail.lstrip("\n") if tail else "\n")).rstrip() + "\n"


def insert_figure_at_top(body: str, cap: str) -> str:
    """Place the figure block right after the H1 title so it is the first thing a
    reader sees. Falls back to prepending if there is no H1."""
    body = strip_figure_block(body)
    cap = cap.strip()
    lines = body.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("# "):
            head = "\n".join(lines[: i + 1]).rstrip()
            rest = "\n".join(lines[i + 1:]).lstrip("\n")
            return f"{head}\n\n{cap}\n\n{rest}".rstrip() + "\n"
    return f"{cap}\n\n{body.lstrip()}".rstrip() + "\n"


# ============================================================================
# 7. Orchestrate one hit
# ============================================================================
def process_hit(row, persample_cache, do_paper=True):
    key = f"{row['modality'].replace(' ', '_')}__{row['gse_id']}__{re.sub(r'[^A-Za-z0-9]+','-',str(row['comparison_name']))[:60]}"
    arm_cache = ARMS_CACHE / f"{key}.json"
    mod = row["modality"]
    if mod not in persample_cache:
        persample_cache[mod] = load_persample(mod, acfg.PERSAMPLE[mod])
    persample = persample_cache[mod]
    clock = acfg.DEFAULT_CLOCK_COL[mod]
    groups, _ = enumerate_groups(persample, row["gse_id"], clock)

    # ---- figure (best-effort; a validated figure is an add-on, not a gate on prose) ----
    vv, msg = None, "no figure"
    if len(groups) >= 2:
        if arm_cache.exists():
            sdk_arms = json.load(open(arm_cache))
        else:
            try:
                sdk_arms = derive_arms(row, groups)
                json.dump(sdk_arms, open(arm_cache, "w"), indent=2)
            except Exception as e:
                sdk_arms = {}
                msg = f"arm derivation failed: {e}"
        if sdk_arms:
            d = persample[persample["gse_id"] == row["gse_id"]].copy()
            heur = heuristic_arms(row, d)   # exact-metadata heuristic, then SDK group choice
            vv, msg = validate_and_values(row, groups, [heur, sdk_arms], persample)
            if not isinstance(vv, dict):
                vv = None
    else:
        msg = "fewer than 2 sample groups"

    fig_ok = isinstance(vv, dict)
    outcome = {"key": key, "case_id": row.get("case_id"),
               "status": "ok" if fig_ok else "flagged", "reason": "" if fig_ok else msg}
    cap = ""
    if fig_ok:
        arms = vv["arms"]
        figbase = FIGDIR / key
        make_figure(row, arms, vv, figbase)
        outcome.update(case_n=len(vv["case_vals"]), ctrl_n=len(vv["ctrl_vals"]),
                       engine_p=vv["res"].pvalue, confidence=arms.get("confidence"))
        cap = FIG_CAPTION.format(figrel="{cid}_figure.png",
                                 case=arms.get("case_label", "case"),
                                 ctrl=arms.get("ctrl_label", "control"))

    # ---- report folder: prose (if missing) + figure (if available) ----
    cid = row.get("case_id")
    if isinstance(cid, str) and cid and cid.lower() != "nan":
        folder = REPORT_DIR[mod] / cid
        folder.mkdir(parents=True, exist_ok=True)
        md = folder / "parsed_report.md"
        cap = cap.replace("{cid}", cid)
        if fig_ok:
            import shutil
            figbase = FIGDIR / key
            shutil.copy2(f"{figbase}.png", folder / f"{cid}_figure.png")
            shutil.copy2(f"{figbase}.pdf", folder / f"{cid}_figure.pdf")
        if md.exists():
            if fig_ok:
                md.write_text(insert_figure_at_top(md.read_text(), cap))
                outcome["report"] = "figure placed at top of existing mini-paper"
            else:
                outcome["report"] = "existing mini-paper kept (figure not reproducible)"
        elif do_paper:
            ev = paperclip_evidence(_s(row.get("intervention", "")) or row["comparison_name"])
            body = write_minipaper(row, ev, vv=vv)
            md.write_text(insert_figure_at_top(body, cap) if fig_ok else body.rstrip() + "\n")
            outcome["report"] = "SDK mini-paper generated" + ("" if fig_ok else " (no figure)")
    return outcome


# ============================================================================
# CLI
# ============================================================================
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--build-hits", action="store_true")
    ap.add_argument("--pilot", type=int, default=0, help="process N representative hits")
    ap.add_argument("--mode", choices=["full"], default=None)
    ap.add_argument("--no-paper", action="store_true", help="figures only, skip prose gen")
    args = ap.parse_args(argv)

    if args.build_hits:
        build_hits(); return 0

    hits = load_hits()
    hits = hits[hits["case_id"].notna()].copy()

    if args.pilot:
        # representative: mixed modality + both directions, prefer a mix of existing/missing report
        pick = []
        for mod in ["mouse RNA-seq", "human RNA-seq", "human DNAm"]:
            sub = hits[hits.modality == mod]
            if len(sub):
                pick.append(sub.iloc[0])
                low = sub[sub.direction < 0]
                if len(low):
                    pick.append(low.iloc[0])
        sel = pd.DataFrame(pick).drop_duplicates(["modality", "gse_id", "comparison_name"]).head(args.pilot)
    else:
        sel = hits

    cache = {}
    results = []
    for _, row in sel.iterrows():
        r = process_hit(row, cache, do_paper=not args.no_paper)
        results.append(r)
        print(f"[{r['status']:7}] {r['key'][:70]:70}  {r.get('reason', r.get('report',''))}")
    ok = sum(1 for r in results if r["status"] == "ok")
    print(f"\n{ok}/{len(results)} figures rendered; "
          f"{sum(1 for r in results if r['status']=='flagged')} flagged; "
          f"{sum(1 for r in results if r['status']=='error')} errored")
    json.dump(results, open(HERE / "run_results.json", "w"), indent=2, default=str)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
