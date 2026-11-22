#!/usr/bin/env python3
"""Restore archived ED11 summaries using July horizontal bar styling.

Copied from the September ED11 source-summary correction builder. The
original July generator was not recovered. Outputs remain in the ED11
restoration folder; source CSV and prior figures are unchanged.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
TERMS = (
    ("mean_0.5", "n_missing_mean", "0.5", "worstcase_bias_yr_(err=0.5)"),
    ("p95_0.5", "n_missing_p95", "0.5", "worstcase_bias_yr_p95_(err=0.5)"),
    ("mean_0.05", "n_missing_mean", "0.05", "realistic_bias_yr_(err=0.05)"),
)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1048576), b""):
            h.update(b)
    return h.hexdigest()


def write_json(path, data):
    with path.open("x", encoding="utf-8") as f:
        json.dump(data, f, indent=2)
        f.write("\n")


def interval(lexeme):
    value = Decimal(lexeme)
    half_step = Decimal(10) ** value.as_tuple().exponent / 2
    return max(Decimal(0), value - half_step), value + half_step


def source_products(rows):
    out = []
    for r in rows:
        w = Decimal(r["mean_abs_weight"])
        wlo, whi = interval(r["mean_abs_weight"])
        for name, nkey, eps, storedkey in TERMS:
            n = Decimal(r[nkey])
            nlo, nhi = interval(r[nkey])
            slo, shi = interval(r[storedkey])
            computed = n * w * Decimal(eps)
            lo, hi = nlo * wlo * Decimal(eps), nhi * whi * Decimal(eps)
            stored = Decimal(r[storedkey])
            out.append({
                "clock": r["clock"], "term": name,
                "missing_count_lexeme": r[nkey], "mean_weight_lexeme": r["mean_abs_weight"],
                "epsilon_assumption": eps, "source_column": storedkey,
                "stored_lexeme": r[storedkey], "rounded_input_product": str(computed),
                "product_minus_stored": str(computed - stored),
                "exact_equality": computed == stored,
                "product_interval_inclusive": [str(lo), str(hi)],
                "stored_rounding_interval_inclusive": [str(slo), str(shi)],
                "rounding_intervals_overlap": max(lo, slo) <= min(hi, shi),
            })
    return out


# Reference tokens measured from July PDF page 17: horizontal grouped bars,
# seaborn deep blue/red/green, sans-serif labels, two adjacent panels.
# The July generator is absent from the searched repository; this is a copy
# of the archived correction builder with the plotting routine reconstructed.
def build():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import fitz
    from PIL import Image, ImageDraw
    source = ROOT / "paper/legacy_reference/figure_rev1/sensitivity_analyses/imputation_sensitivity_results.csv"
    out = ROOT / "paper/revision_2026-09-07/evidence/figures-rebuild-2026-11-05/ed11"
    out.mkdir(parents=True, exist_ok=True)
    if not source.is_file():
        raise FileNotFoundError(f"Cannot proceed - required data file not found: {source}")
    with source.open() as f:
        rows = list(csv.DictReader(f))
    order = ["DunedinPoAm", "HannumAge", "DunedinPACE", "HorvathAge", "PedBE", "ZhangAge", "PhenoAge"]
    lookup = {r["clock"]: r for r in rows}
    rows = [lookup[c] for c in order]
    assert len(rows) == 7
    plt.rcParams.update({"font.family":"sans-serif", "font.sans-serif":["Arial", "DejaVu Sans"],
                        "font.size":7, "axes.labelsize":7, "xtick.labelsize":7,
                        "ytick.labelsize":7, "legend.fontsize":7,
                        "pdf.fonttype":42, "ps.fonttype":42, "axes.linewidth":.6})
    fig, (a,b) = plt.subplots(1,2, figsize=(183/25.4, 68/25.4))
    fig.subplots_adjust(left=.125, right=.99, bottom=.18, top=.85, wspace=.27)
    y=np.arange(7); h=.36
    med=np.array([float(r["median_missing_frac"]) for r in rows])
    p95=np.array([float(r["p95_missing_frac"]) for r in rows])
    terms=np.array([float(r["realistic_bias_yr_(err=0.05)"]) for r in rows])
    a.barh(y-h/2, med, height=h, color="#4c72b0", label="Median")
    a.barh(y+h/2, p95, height=h, color="#c44e52", label="95th percentile")
    a.axvline(.1, color="#bdbdbd", linestyle="--", linewidth=.7, label="10% missing CpGs")
    a.set_xlim(0,.30); a.set_xticks([0,.05,.10,.15,.20,.25,.30])
    a.set_xlabel("Fraction of clock CpGs missing")
    b.barh(y, terms, height=.76, color="#55a868")
    b.set_xlim(0,2.4); b.set_xticks([0,.5,1,1.5,2])
    b.set_xlabel("Scenario term (clock-specific scale)")
    for ax in (a,b):
        ax.set_yticks(y, order); ax.set_ylim(6.6,-.6)
        ax.spines[["top","right"]].set_visible(False)
        ax.tick_params(axis="y",length=0,pad=3)
        ax.tick_params(axis="x",width=.6,length=3,pad=2)
    a.set_title("Per-clock missingness", fontsize=8,pad=10)
    b.set_title("Imputation scenario term (ε = 0.05)", fontsize=8,pad=10)
    for ax, label in zip((a,b),('a','b')):
        ax.text(-.19,1.10,label,transform=ax.transAxes,fontsize=10,fontweight="bold",va="bottom")
    # The July lower-right legend position is clear of both short bottom bars.
    a.legend(loc="lower right",frameon=True,edgecolor="#dddddd",facecolor="white",
             borderpad=.3, handlelength=1.5,labelspacing=.2,handletextpad=.4)
    fig.canvas.draw()
    textcheck=[]
    for ax in fig.axes:
        for t in [*ax.get_xticklabels(),*ax.get_yticklabels(),ax.xaxis.label,ax.yaxis.label,ax.title]:
            if t.get_text():
                textcheck.append({"text":t.get_text(),"pt":t.get_fontsize()})
    assert min(t["pt"] for t in textcheck)>=7
    # Full skill checker is recorded without applying its layout/style recipe.
    sys.path.insert(0,str(Path("pipeline_assets/figure_style_2026_09_10")))
    try:
        from nature_figure_utils import figure_preflight
        rep=figure_preflight(fig,{"a":a,"b":b})
        (out/"skill_preflight.txt").write_text(str(rep))
    except Exception as exc:
        (out/"skill_preflight.txt").write_text("Checker error: "+str(exc))
    fig.savefig(out/"ed11.pdf")
    fig.savefig(out/"ed11_review.png",dpi=200)
    plt.close(fig)
    doc=fitz.open(out/"ed11.pdf"); page=doc[0]
    assert abs(page.rect.width-183/25.4*72)<.01 and page.rect.height<=247/25.4*72
    assert len(page.get_images())==0
    spans=[s for block in page.get_text("dict")["blocks"] if "lines" in block for line in block["lines"] for s in line["spans"]]
    assert all(s["size"]>=6.99 for s in spans)
    original=fitz.open(ROOT/"paper/revision_2026-09-07/baseline/03_Figures_2026-07-06.pdf")
    jp=original[16]
    # Crop the July artwork to exclude package header and page whitespace.
    # The full artwork widths are normalized to 183 mm for direct comparison.
    target_width=1441
    july_clip=fitz.Rect(80,101,550,260.5)
    jp_pix=jp.get_pixmap(matrix=fitz.Matrix(3,3),clip=july_clip)
    july_im=Image.frombytes("RGB",[jp_pix.width,jp_pix.height],jp_pix.samples)
    pp=page.get_pixmap(matrix=fitz.Matrix(3,3))
    new_im=Image.frombytes("RGB",[pp.width,pp.height],pp.samples)
    normalized=[]
    for im in (july_im,new_im):
        normalized.append(im.resize((target_width,round(im.height*target_width/im.width)),Image.Resampling.LANCZOS))
    margin=24; label_height=38; gap=28
    joined=Image.new("RGB",(2*target_width+2*margin+gap,max(im.height for im in normalized)+label_height+2*margin),"white")
    draw=ImageDraw.Draw(joined)
    for i,(im,label) in enumerate(zip(normalized,("July reference", "Restored"))):
        x=margin+i*(target_width+gap)
        draw.text((x,margin),label,fill="black")
        joined.paste(im,(x,margin+label_height))
    joined.save(out/"ed11_vs_july.png")
    with (out/"plotted_values.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    report={"source":str(source),"source_sha256":sha(source),"source_mtime_utc":datetime.fromtimestamp(source.stat().st_mtime,timezone.utc).isoformat(),
            "order":order,"median":med.tolist(),"p95":p95.tolist(),"scenario_terms":terms.tolist(),
            "width_mm":183,"height_mm":68,"embedded_images":0,"pdf_text_min_pt":min(s["size"] for s in spans),"text_checks":textcheck,
            "preprocessing":"Source rows reordered to July order only; no recalculation or filtering.",
            "scenario_parameters":{"epsilon":.05,"formula":"n_missing_mean * mean_abs_weight * epsilon; stored rounded source values plotted"}}
    (out/"checks.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))

if __name__ == "__main__":
    try:
        build()
    except (OSError,ValueError,KeyError,RuntimeError,AssertionError) as exc:
        raise SystemExit(f"ED11 restoration stopped: {exc}") from exc
