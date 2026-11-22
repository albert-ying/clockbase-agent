#!/usr/bin/env python3
"""
GenAge / DrugAge comparability audit.

Implements revision_r1_2026-05-18/comparability_audit_protocol.md: a
directional-concordance match is retained only when the GEO experimental
setting is comparable to the GenAge/DrugAge basis. Gates are defined purely
from experimental-setting features and applied blind to agreement. Every
matched case is written out with its comparability verdict and (if excluded)
the first gate it failed.

Run from the Aging-AI-Agent directory:  python3 comparability_audit.py
"""

import re
from pathlib import Path

import pandas as pd
from scipy import stats

OUTDIR = Path("../revision_r1_2026-05-18")

# ---- context / tooling vocabularies --------------------------------------
CONTEXT_EXCLUDE = [
    "tumor", "tumour", "cancer", "carcinoma", "glioma", "glioblastoma",
    "leukemia", "leukaemia", "lymphoma", "melanoma", "sarcoma", "metastas",
    "hela", "ipsc", " esc", "esc ", "embryonic stem", "organoid",
    "embryonic", "fetal", "foetal", "neonatal", "cell line", "immortalized",
    "immortalised",
]
LOF_KW = ["knockout", "knock-out", "knock out", "-/-", " ko", "ko ", "deletion",
          "deleted", "deficiency", "deficient", "null", "knockdown",
          "knock-down", " kd", "kd ", "shrna", "sirna", "rnai", "silencing"]
GOF_KW = ["overexpression", "overexpress", " oe", "oe ", "transgenic",
          "knock-in", "knockin", "knock in", "gain-of-function", "gain of function"]
MULTI_SEP = ["/", "+", "&", " and ", " x ", ";", " plus "]
TOOL_DRUGS = {"doxycycline", "tetracycline", "tamoxifen", "4-hydroxytamoxifen",
              "4-oht", "hydroxytamoxifen", "puromycin", "geneticin", "g418",
              "blasticidin", "hygromycin", "iptg", "dmso", "cre"}


def _t(row, *fields):
    return " ".join(str(row.get(f, "") or "") for f in fields).lower()


def genage_gate(row):
    """Return '' if comparable, else the first failed gate id."""
    cond = _t(row, "condition_name", "comparison_name", "intervention")
    raw = str(row.get("condition_name", "") or "") + " " + str(row.get("intervention", "") or "")
    # G1 single-target
    if any(sep in raw for sep in MULTI_SEP):
        return "G1_multi_target"
    # G2 defined LOF/GOF
    if not (any(k in cond for k in LOF_KW) or any(k in cond for k in GOF_KW)):
        return "G2_no_LOF_GOF"
    # G3 physiological context
    if any(k in cond for k in CONTEXT_EXCLUDE):
        return "G3_context"
    return ""


def drugage_gate(row):
    cond = _t(row, "condition_name", "comparison_name", "intervention")
    raw = str(row.get("condition_name", "") or "") + " " + str(row.get("intervention", "") or "")
    drug = str(row.get("drugage_matched_drugs", "") or "").strip().lower()
    # D1 tool reagent
    if drug in TOOL_DRUGS:
        return "D1_tool_reagent"
    # D2 combination
    if any(sep in raw for sep in MULTI_SEP):
        return "D2_combination"
    # D3 context
    if any(k in cond for k in CONTEXT_EXCLUDE):
        return "D3_context"
    # D4 direction defined
    ann = str(row.get("drugage_annotation", "") or "")
    if ann in ("", "Unclear", "nan"):
        return "D4_unclear"
    return ""


def score_col(df):
    return "final_score" if "final_score" in df.columns else "trust_score"


def audit(df, kind):
    """kind = 'genage' or 'drugage'. Returns (comparable_df, all_audited_df)."""
    sc = score_col(df)
    df = df[pd.to_numeric(df[sc], errors="coerce") > 80].copy()
    matched_col = "genage_annotated" if kind == "genage" else "drugage_annotated"
    m = df[df[matched_col] == True].copy()
    gate = genage_gate if kind == "genage" else drugage_gate
    m["failed_gate"] = m.apply(gate, axis=1)
    m["comparable"] = m["failed_gate"] == ""
    return m[m["comparable"]].copy(), m


def concordance(comp, kind, label):
    gene_col = "genage_matched_genes" if kind == "genage" else "drugage_matched_drugs"
    ann_col = "genage_annotation" if kind == "genage" else "drugage_annotation"
    rows = []
    for g, grp in comp.groupby(gene_col):
        inf = grp["inferred_influence"]
        maj = inf.mode().iloc[0] if len(inf.mode()) else ""
        ann = grp[ann_col].iloc[0]
        if ann in ("", "Unclear", "nan") or maj == "":
            continue
        rows.append((g, ann, maj, ann == maj))
    res = pd.DataFrame(rows, columns=[gene_col, "annotation", "inferred", "agree"])
    n = len(res)
    a = int(res["agree"].sum())
    p = stats.binomtest(a, n, 0.5).pvalue if n else float("nan")
    print(f"  {label}: comparable unique {kind[:-3] if False else kind} entries = {n}; "
          f"agreement {a}/{n} = {100*a/n:.1f}%  (binomial p = {p:.3g})" if n
          else f"  {label}: 0 comparable entries")
    return res


def run(kind, mouse_file, human_files):
    print(f"\n{'='*66}\n{kind.upper()} comparability audit\n{'='*66}")
    # Mouse-only
    mdf = pd.read_csv(mouse_file, low_memory=False)
    mcomp, mall = audit(mdf, kind)
    print(f"\n[mouse-only]  matched (score>80): {len(mall)}  "
          f"-> comparable: {len(mcomp)}  excluded: {len(mall)-len(mcomp)}")
    print("  exclusions by gate:", mall[~mall["comparable"]]["failed_gate"].value_counts().to_dict())
    mres = concordance(mcomp, kind, "mouse-only")
    # Human combined
    hdf = pd.concat([pd.read_csv(f, low_memory=False) for f in human_files],
                    ignore_index=True)
    hcomp, hall = audit(hdf, kind)
    print(f"\n[human RNA+DNAm combined]  matched (score>80): {len(hall)}  "
          f"-> comparable: {len(hcomp)}  excluded: {len(hall)-len(hcomp)}")
    print("  exclusions by gate:", hall[~hall["comparable"]]["failed_gate"].value_counts().to_dict())
    hres = concordance(hcomp, kind, "human-combined")
    # write audit trail
    mall["cohort"] = "mouse"; hall["cohort"] = "human"
    OUTDIR.mkdir(exist_ok=True)
    keepcols = [c for c in ["cohort", "gse_id", "condition_name", "intervention",
                "genage_matched_genes", "drugage_matched_drugs",
                "genage_annotation", "drugage_annotation", "inferred_influence",
                "comparable", "failed_gate"] if c in mall.columns or c in hall.columns]
    out = pd.concat([mall, hall], ignore_index=True)
    out[[c for c in keepcols if c in out.columns]].to_csv(
        OUTDIR / f"comparability_audit_{kind}.csv", index=False)
    print(f"\n  audit trail -> revision_r1_2026-05-18/comparability_audit_{kind}.csv")
    return mres, hres


if __name__ == "__main__":
    run("genage", "data/agent_result_genage.csv",
        ["data/agent_result_genage_hs_rna.csv", "data/agent_result_genage_hs_dnam.csv"])
    run("drugage", "data/agent_result_drugage_significant.csv",
        ["data/agent_result_drugage_hs_rna.csv", "data/agent_result_drugage_hs_dnam.csv"])
