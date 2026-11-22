#!/usr/bin/env python3
"""
Combined-modality GenAge waffle plot.

Pools genetic interventions from all 3 modalities:
  - data/agent_result_genage.csv          (mouse RNA-seq, trust_score)
  - data/agent_result_genage_hs_rna.csv   (human RNA-seq, final_score)
  - data/agent_result_genage_hs_dnam.csv  (human DNAm, final_score)

Per-modality SCORE_FLOOR (matches volcano):
    mouse RNA-seq : 100
    human RNA-seq : 100
    human DNAm    :  80

Outputs:
  figures/create_geneage_waffle_plot_nature_combined.png
  figures/create_geneage_waffle_plot_nature_combined.pdf
"""

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.stats import binomtest

warnings.filterwarnings("ignore")

USE_MERGED = True

try:
    from aquarel import load_theme
    theme = (
        load_theme("scientific")
        .set_grid(draw=False)
        .set_font(size=14)
        .set_ticks(direction="out")
        .set_axis_labels(pad=10)
    )
    theme.apply()
except Exception:
    theme = None

plt.rcParams["savefig.transparent"] = False
plt.rcParams["figure.facecolor"] = "white"
plt.rcParams["axes.facecolor"] = "white"
plt.rcParams["savefig.facecolor"] = "white"
plt.rcParams["font.family"] = "Arial"

# Use trust_score > 80 for all modalities, matching the original mouse-only
# waffle (create_geneage_waffle_plot_nature.py) which also uses trust_score > 80.
# A floor of 100 applied to the human files would keep < 5% of rows.
SCORE_FLOOR = {"mouse RNA-seq": 80, "human RNA-seq": 80, "human DNAm": 80}

MODALITY_FILES = [
    ("mouse RNA-seq", "data/agent_result_genage.csv",         "trust_score"),
    ("human RNA-seq", "data/agent_result_genage_hs_rna.csv",  "trust_score"),
    ("human DNAm",    "data/agent_result_genage_hs_dnam.csv", "trust_score"),
]


def load_combined():
    frames = []
    for modality, path, score_col in MODALITY_FILES:
        if not Path(path).exists():
            print(f"  WARN: {path} not found, skipping {modality}")
            continue
        d = pd.read_csv(path, low_memory=False)
        d = d.drop_duplicates()
        # Apply per-modality SCORE_FLOOR (use the modality's score column)
        if score_col not in d.columns:
            print(f"  WARN: {path} missing {score_col}; skipping {modality}")
            continue
        floor = SCORE_FLOOR[modality]
        before = len(d)
        d = d[d[score_col] > floor].copy()
        d["modality"] = modality
        # Normalise columns we depend on
        if "genage_annotated" not in d.columns:
            print(f"  WARN: {path} missing genage_annotated; skipping {modality}")
            continue
        d["genage_annotated"] = d["genage_annotated"].astype(str).str.upper().isin(["TRUE", "1"])
        keep = ["modality", "genage_matched_genes", "genage_annotated",
                "genage_annotation", "inferred_influence", "intervention"]
        keep = [c for c in keep if c in d.columns]
        frames.append(d[keep])
        print(f"  {modality:14s}: {len(d):,} rows after {score_col} > {floor} (was {before:,})")
    if not frames:
        raise RuntimeError("no GenAge tables loaded")
    return pd.concat(frames, ignore_index=True)


df = load_combined()

print("\n" + "=" * 60)
print("Combined GenAge Database Validation Analysis")
print("=" * 60)
print(f"Total combined rows: {len(df):,}")

matched_genes = df[df["genage_annotated"] == True].copy()
print(f"Matched with GenAge: {len(matched_genes):,}")

matched_genes["agreement"] = (
    matched_genes["genage_annotation"] == matched_genes["inferred_influence"]
)

if USE_MERGED:
    print("\n=== MERGED mode (majority voting per gene across modalities) ===")
    merged = []
    for gene, grp in matched_genes.groupby("genage_matched_genes"):
        if pd.isna(gene) or not gene:
            continue
        votes = grp["inferred_influence"].value_counts()
        majority = votes.idxmax() if not votes.empty else None
        annotation = grp["genage_annotation"].iloc[0]
        agree = (annotation == majority) if majority is not None else False
        merged.append({
            "genage_matched_genes": gene,
            "inferred_influence": majority,
            "genage_annotation": annotation,
            "agreement": agree,
            "num_interventions": len(grp),
            "modalities": ",".join(sorted(grp["modality"].unique())),
        })
    matched_genes_display = pd.DataFrame(merged)
    agreement_count = int(matched_genes_display["agreement"].sum())
    disagreement_count = len(matched_genes_display) - agreement_count
    not_in_genage = len(df) - len(matched_genes)
    print(f"Merged to {len(matched_genes_display):,} unique genes (from {len(matched_genes):,} interventions)")
    print(f"  Agreement: {agreement_count}")
    print(f"  Disagreement: {disagreement_count}")
    print(f"Not in GenAge: {not_in_genage}")
else:
    matched_genes_display = matched_genes
    agreement_count = int(matched_genes["agreement"].sum())
    disagreement_count = len(matched_genes) - agreement_count
    not_in_genage = len(df) - len(matched_genes)

# Statistics
if len(matched_genes_display) > 0:
    agreement_rate = agreement_count / len(matched_genes_display)
    p_value = binomtest(agreement_count, len(matched_genes_display),
                        p=0.5, alternative="greater").pvalue
else:
    agreement_rate = 0.0
    p_value = 1.0

print(f"\nAgreement rate: {agreement_rate:.1%}")
print(f"Binomial test p-value: {p_value:.4e}")

# Waffle data
TOTAL = len(matched_genes_display) + not_in_genage if USE_MERGED else len(df)
MATCHED = len(matched_genes_display)
matched_gene_names = matched_genes_display["genage_matched_genes"].tolist()

# Colors (Nature palette)
AGREE_COLOR = "#3B8EA5"
DISAGREE_COLOR = "#E64B35"
NOT_IN_GENAGE_COLOR = "#D3D3D3"
BORDER_COLOR = "#FFFFFF"

# Choose grid
n_cols = 100
n_rows = int(np.ceil(TOTAL / n_cols))
square_size = 0.18
gap = 0.003

# Build waffle data: matched first (in deterministic order), then non-matched
waffle_data = ["agree" if r["agreement"] else "disagree"
               for _, r in matched_genes_display.iterrows()]
waffle_data.extend(["not_in_genage"] * not_in_genage)
assert len(waffle_data) == TOTAL, (len(waffle_data), TOTAL)

# Figure
fig = plt.figure(figsize=(24, max(8, n_rows * 0.5 + 2.5)))
ax = fig.add_subplot(111)

square_idx = 0
matched_positions = []
for row in range(n_rows - 1, -1, -1):
    for col in range(n_cols):
        if square_idx >= len(waffle_data):
            break
        kind = waffle_data[square_idx]
        if kind == "agree":
            color = AGREE_COLOR; alpha = 1.0; ew = 0.2
        elif kind == "disagree":
            color = DISAGREE_COLOR; alpha = 1.0; ew = 0.2
        else:
            color = NOT_IN_GENAGE_COLOR; alpha = 0.6; ew = 0.15
        x_lo = col * (square_size + gap)
        y_lo = row * (square_size + gap)
        ax.add_patch(Rectangle(
            (x_lo, y_lo), square_size, square_size,
            facecolor=color, edgecolor=BORDER_COLOR,
            linewidth=ew, alpha=alpha,
        ))
        if kind in ("agree", "disagree"):
            matched_positions.append((x_lo + square_size / 2,
                                      y_lo + square_size / 2, kind))
        square_idx += 1

# Gene labels alternating above/below
label_padding = 0.15
line_padding = 0.05
for idx, (x, y, kind) in enumerate(matched_positions):
    if idx >= len(matched_gene_names):
        break
    name = matched_gene_names[idx]
    if pd.isna(name) or not name:
        continue
    text_color = AGREE_COLOR if kind == "agree" else DISAGREE_COLOR
    if idx % 2 == 0:
        sq_top = y + square_size / 2
        label_y = sq_top + label_padding
        va_align = "bottom"; ha_align = "left"
        ax.plot([x, x], [sq_top, label_y], color=text_color,
                linewidth=2, alpha=0.8, zorder=10)
    else:
        sq_bot = y - square_size / 2
        label_y = sq_bot - label_padding
        va_align = "top"; ha_align = "right"
        ax.plot([x, x], [sq_bot, label_y], color=text_color,
                linewidth=2, alpha=0.8, zorder=10)
    ax.text(x, label_y, name, fontsize=14, ha=ha_align, va=va_align,
            color=text_color, weight="bold", rotation=45, zorder=11)

ax.set_xlim(-0.5, n_cols * (square_size + gap))
ax.set_ylim(-0.8, n_rows * (square_size + gap) + 1.5)
ax.set_aspect("equal")
ax.axis("off")

title_y = n_rows * (square_size + gap) + 0.8
ax.text(n_cols * (square_size + gap) / 2, title_y,
        f"GenAge Database Validation, all modalities combined (n = {TOTAL})",
        fontsize=22, fontweight="bold", ha="center", color="#000000")
if MATCHED > 0:
    ax.text(n_cols * (square_size + gap) / 2, title_y - 0.35,
            f"Binomial test: p = {p_value:.4e}",
            fontsize=16, ha="center", color="#333333")

legend_elements = [
    Line2D([0], [0], marker="s", color="w", markerfacecolor=AGREE_COLOR,
           markersize=14, label=f"Agreement with GenAge ({agreement_count})"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=DISAGREE_COLOR,
           markersize=14, label=f"Disagreement with GenAge ({disagreement_count})"),
    Line2D([0], [0], marker="s", color="w", markerfacecolor=NOT_IN_GENAGE_COLOR,
           alpha=0.6, markersize=14,
           label=f"Not in GenAge database ({not_in_genage})"),
]
ax.legend(handles=legend_elements, loc="upper center",
          bbox_to_anchor=(0.5, 1), ncol=3, frameon=False,
          fontsize=14, columnspacing=3.0)

if theme is not None:
    try:
        theme.apply_transforms()
    except Exception:
        pass

plt.tight_layout()
plt.savefig("figures/create_geneage_waffle_plot_nature_combined.pdf",
            dpi=300, bbox_inches="tight", facecolor="white", edgecolor="none")
plt.savefig("figures/create_geneage_waffle_plot_nature_combined.png",
            dpi=600, bbox_inches="tight", facecolor="white", edgecolor="none")
print("Saved figures/create_geneage_waffle_plot_nature_combined.{pdf,png}")

print("\n" + "=" * 60)
print("Combined GenAge Validation Summary")
print("=" * 60)
print(f"Total interventions analysed: {TOTAL:,}")
print(f"Matched with GenAge (unique): {MATCHED:,}")
if MATCHED:
    print(f"Agreement: {agreement_count}/{MATCHED} ({agreement_rate:.1%})")
    print(f"Disagreement: {disagreement_count}/{MATCHED}")
    print(f"Binomial test p-value: {p_value:.4e}")
print(f"Not in GenAge: {not_in_genage}")
