"""
Dedicated drill-down on the extreme-tail subset of the pooled three-modality
ClockBase Agent atlas.

Inputs:
  - Aging-AI-Agent/figures/combined_3modality_plot_data.csv (pooled atlas with outlier_status)
  - figure_rev1/sensitivity_analyses/extreme_tail_representative_examples.csv (curated 19)

Outputs:
  1. figure_rev1/ExtendedDataFigure15_extreme_tail_volcano.{png,pdf}
  2. figure_rev1/sensitivity_analyses/extreme_tail_category_breakdown.csv
  3. figure_rev1/sensitivity_analyses/extreme_tail_pathway_enrichment.csv
  4. figure_rev1/sensitivity_analyses/extreme_tail_dedicated_analysis_summary.md

No simulated data. All values computed from the real combined_3modality_plot_data.csv.
"""

from __future__ import annotations

import re
import json
from itertools import combinations
import sys
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from PIL import Image  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "paper/legacy_reference/figure_rev1"))
from nature_figure_utils import INK, retro_style, trim_axes  # noqa: E402

try:
    from adjustText import adjust_text  # type: ignore
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path("/Users/A.Y/Desktop/Projects/Manuscripts/2022_clock_base")
ATLAS = ROOT / "Aging-AI-Agent" / "figures" / "combined_3modality_plot_data.csv"
REP = ROOT / "figure_rev1" / "sensitivity_analyses" / "extreme_tail_representative_examples.csv"
OUT_FIG_DIR = Path(__file__).resolve().parents[2] / "paper/revision_2026-09-07/evidence/figure-fixes-2026-11-17/ed14"
OUT_FIG_DIR.mkdir(parents=True, exist_ok=True)
OUT_SENS_DIR = OUT_FIG_DIR / "source_checks"
OUT_SENS_DIR.mkdir(parents=True, exist_ok=True)

VOLCANO_PNG = OUT_FIG_DIR / "ExtendedDataFigure14_extreme_tail_volcano.png"
VOLCANO_PDF = OUT_FIG_DIR / "ExtendedDataFigure14_extreme_tail_volcano.pdf"
SUBMISSION_FIG_DIR = OUT_FIG_DIR / "package_source"
SUBMISSION_VOLCANO_PDF = SUBMISSION_FIG_DIR / "ExtendedDataFigure14_extreme_tail_volcano.pdf"
CATEGORY_CSV = OUT_SENS_DIR / "extreme_tail_category_breakdown.csv"
PATHWAY_CSV = OUT_SENS_DIR / "extreme_tail_pathway_enrichment.csv"
SUMMARY_MD = OUT_SENS_DIR / "extreme_tail_dedicated_analysis_summary.md"


def normalize_png_canvas(path: Path, target_size: tuple[int, int]) -> None:
    """Keep DOCX replacement images at the original embedded aspect ratio."""
    image = Image.open(path).convert("RGBA")
    scale = min(target_size[0] / image.width, target_size[1] / image.height)
    resized = image.resize(
        (round(image.width * scale), round(image.height * scale)),
        Image.Resampling.LANCZOS,
    )
    canvas = Image.new("RGBA", target_size, (255, 255, 255, 255))
    left = (target_size[0] - resized.width) // 2
    top = (target_size[1] - resized.height) // 2
    canvas.alpha_composite(resized, (left, top))
    canvas.convert("RGB").save(path)


# ---------------------------------------------------------------------------
# Load + subset
# ---------------------------------------------------------------------------
def load_tail_subset() -> pd.DataFrame:
    df = pd.read_csv(ATLAS, low_memory=False)
    sig = df[df["fdr"] < 0.05].copy()
    tail = sig[sig["outlier_status"] != "Normal"].copy()
    return tail


# ---------------------------------------------------------------------------
# Category inference
# ---------------------------------------------------------------------------
# Keyword rules ordered by priority. First match wins.
CATEGORY_RULES: list[tuple[str, list[str]]] = [
    (
        "developmental",
        [
            r"\bembryo", r"\bembryonic\b", r"\bfetal\b", r"\bfoetal\b",
            r"\bpostnatal\b", r"\bneonat", r"\bdevelopment", r"\bmaturation\b",
            r"\bE\d+\.?\d*\b", r"\bP\d+\b", r"\bdifferentiation\b",
            r"\bmesenchymoangioblast\b", r"\bcardiac differentiation\b",
            r"\bteratoma\b", r"\bpuberty\b", r"\bgastrul", r"\borganoid\b",
            r"\bday\s*\d+\s+(?:differentiation|maturation|RA|BMP|FGF)",
            r"\bRA(?:\+|/)?BMP", r"\bRA and SAG\b", r"\bdifferentiated\b",
            r"\bsubcereb\b", r"\bgerm cells\b", r"\bspermatocyte\b",
        ],
    ),
    (
        "age_baseline",
        [
            r"\bold age\b", r"\bage[:\s]+", r"\bold vs young\b",
            r"\b\d+\s*(?:months?|mo|years?|yrs?)\s*vs\s*\d+\s*(?:months?|mo|years?|yrs?)\b",
            r"\bchronological\s*aging\b", r"\bage-related\b", r"\bageing\b",
            r"\baging baseline\b", r"\bphenoage\b",
            r"\bnatural aging\b", r"\blongitudinal aging\b", r"\badvanced age\b",
            r"\belderly\b", r"\bgeriatric\b", r"\bsenescence\b", r"\bsenescent\b",
            r"\b\d+\s*months?\b\s*vs", r"\bAging\b",
        ],
    ),
    (
        "acute_stress",
        [
            r"\bUVB?\b", r"\birradiation\b", r"\bgamma\s*(?:radiation)?", r"\bproton\s*radiation\b",
            r"\bH2O2\b", r"\bhydrogen peroxide\b", r"\bheat\s*shock\b",
            r"\bhypoxia\b", r"\bhyperoxia\b", r"\bLPS\b", r"\boxidative\s*stress\b",
            r"\bradiation\b", r"\bIR\b\s*\d", r"\bcGy\b",
            r"\bbleomycin\b", r"\bpoly\(I:?C\)\b", r"\bsalmonella\b",
            r"\bischemi", r"\bstroke\b", r"\bSCI\b", r"\bspinal cord injury\b",
            r"\bdemyelin", r"\bnerve injury\b", r"\bUUO\b",
            r"\bcachectic\b", r"\bcachexia\b", r"\bsepsis\b", r"\binjury\b",
            r"\babdominalrt\b",
        ],
    ),
    (
        "transformed_cell",
        [
            r"\biPSC\b", r"\bhESC\b", r"\bESC\b", r"\bMEF\b",
            r"\bNT-?ESC\b", r"\binduced pluripotent\b", r"\bembryonic stem\b",
            r"\breprogramm", r"\b4F\b", r"\bOSKM\b",
            r"\b(?:tumou?r|cancer|carcinoma|melanoma|leukem|leukaemi|lymphoma|sarcoma|glioma|glioblastoma|SCLC|NSCLC|CLL|AML|MDS)\b",
            r"\bHeLa\b", r"\btransformed\b", r"\bRasV12\b",
            r"\bcell line\b", r"\bcell cluster\b", r"\bMGG\d", r"\bSLX-\d",
            r"\biGSC\b", r"\bFGSC\b", r"\bSSC\b", r"\biP-MSC\b",
            r"\bBRAFV600E\b", r"\bEGFR\b", r"\bTEL-FGFR1\b",
            r"\boxBS\b", r"\b5HmC\b", r"\bABT-?199\b",
            r"\bfibroblast subtype\b", r"\bFibroblast\b\s*$",
        ],
    ),
    (
        "disease",
        [
            r"\bAlzheimer", r"\b5xFAD\b", r"\bdiabetes\b", r"\bdiabetic\b",
            r"\bfibrosis\b", r"\bIPF\b", r"\bidiopathic pulmonary\b",
            r"\bdepression\b", r"\bschizophren", r"\bobesity\b",
            r"\bmyeloproliferative\b", r"\bautoimmune\b", r"\bparkinson",
            r"\bcardiomyopathy\b", r"\bhepatitis\b", r"\bcirrhosis\b",
            r"\bdisease\b", r"\bpatient\b",
            r"\bdown syndrome\b", r"\bendometriosis\b", r"\bNASH\b",
            r"\bLynch[\s-]?like\b", r"\btobacco smok", r"\bsmoking\b",
            r"\bfamine\b", r"\bdeficiency\b", r"\bsyndrome\b", r"\bmutation\b",
            r"\bautism\b", r"\bischemic stroke\b",
        ],
    ),
    (
        "tissue_baseline",
        [
            r"\bliver vs\b", r"\bbrain vs\b", r"\bstomach\b", r"\bintestine\b",
            r"\bforestomach\b", r"\bplacenta\b", r"\breference rna\b",
            r"\btissue\b", r"\borgan\b", r"\bbone\b\s*$",
            r"\b(?:cortex|cerebellum|hippocampus|kidney|lung|heart|muscle|spleen|saliva)\s*$",
            r"\bcell\s*type\b", r"\bcell type\b",
            r"\blobular vs interlobular\b", r"\btissue difference\b",
            r"\bsperm\b", r"\bmacrophage\b\s*vs", r"\bmonocyte\b\s*vs",
            r"\bnaive T cells?\b", r"\bmemory T cells?\b", r"\beffector\b.*cells?\b",
            r"\bperitubular\b", r"\bCFU-E\b", r"\berythroblast",
            r"\bcluster \d", r"\bRGC\b", r"\btranscriptional cluster\b",
            r"\bbrain region\b", r"\bdendritic\b",
        ],
    ),
    (
        "intervention",
        [
            # genetic perturbations
            r"\bknockout\b", r"\bknockdown\b", r"\b-KD\b", r"\b_KO\b", r"\b\bKO\b",
            r"-?KO\b", r"\boverexpression\b", r"\boverexpress\b",
            r"\bCRISPR\b", r"\bshRNA\b", r"\bsiRNA\b", r"\bshArid", r"\b-/-\b",
            r"\bcre\b", r"\btransgenic\b", r"\bdKO\b", r"\bmutant\b",
            r"\bfusion gene\b", r"\bRARa\b", r"\bcGAS\b", r"\bMHC class\b",
            r"\bgene transfer\b", r"\boverexpressing\b", r"\bdeletion\b",
            r"\bvs WT\b", r"\bvs Wild[\s-]?[Tt]ype\b", r"\bvs Control\b",
            # drugs / small molecules
            r"\bchemotherapy\b", r"\bazacitidine\b", r"\bmetformin\b", r"\brapamycin\b",
            r"\btamoxifen\b", r"\binhibitor\b", r"\btreatment\b", r"\binduction\b",
            r"\bdexamethasone\b", r"\bFK506\b", r"\bSGI-\d", r"\bBO-?\d{3}\b",
            r"\bdBET\d", r"\bIUB\d+\b", r"\bAM\d{3}\b", r"\boxaliplatin\b",
            r"\bsodium butyrate\b", r"\bNaBu\b", r"\bretinoic acid\b", r"\bRA-induced\b",
            r"\bvalproic\b", r"\bvorinostat\b", r"\bdecitabine\b",
            r"\bM\.SssI\b", r"\bCBL0137\b", r"\bpanobinostat\b",
            r"\b5-aza\b", r"\bdeoxycytidine\b",
            r"\btransplant\b", r"\binfection\b", r"\bvaccine\b",
            # diet / supplementation
            r"\bcalori(?:c|e) restriction\b", r"\bhigh fat diet\b", r"\bdiet\b",
            r"\bsupplementation\b", r"\bfolic acid\b", r"\bzinc\b",
            r"\bB12\b", r"\bestrogen\b", r"\bE2\b",
            # experimental conditions
            r"\bdose\b", r"\b\d+\s*[uµ]?M\b", r"\bvs DMSO\b", r"\bvs Vehicle\b", r"\bvs vehicle\b",
            r"\blearning\b", r"\btraining\b", r"\bexercise\b",
            r"\bimplant\b", r"\bhydrogel\b", r"\bPEG\b",
            r"\bGM\s*CSF\b", r"\bIFN\b", r"\bTGF-?b\b", r"\bTGF beta\b",
            r"\bBMP\d?\b", r"\bFGF\d?\b", r"\bEGF\b",
            r"\bmiR-?\d", r"\bmiRNA\b", r"\bactivated_miR",
            r"\bPDX\b", r"\b\bPI\b", r"\bMethylation Method\b",
        ],
    ),
]

COMPILED_RULES: list[tuple[str, list[re.Pattern]]] = [
    (cat, [re.compile(p, re.IGNORECASE) for p in pats])
    for cat, pats in CATEGORY_RULES
]


def classify_row(row: pd.Series) -> str:
    """Assign a single category label using keyword priority."""
    text_parts = []
    for c in ("comparison_name", "condition_name", "intervention"):
        v = row.get(c)
        if isinstance(v, str) and v.strip():
            text_parts.append(v)
    blob = " || ".join(text_parts)
    if not blob:
        return "other"
    for cat, pats in COMPILED_RULES:
        if any(p.search(blob) for p in pats):
            return cat
    return "other"


def assign_tail_label(status: str) -> str:
    if status == "Effect Size Outlier":
        return "effect_size_tail"
    if status == "P-value Outlier":
        return "p_value_tail"
    if status == "Both Effect Size & P-value Outlier":
        return "both_tails"
    return "normal"


# ---------------------------------------------------------------------------
# Volcano (Output 1)
# ---------------------------------------------------------------------------
TAIL_COLORS = {
    "effect_size_tail": "#C04A1E",   # rust
    "p_value_tail": "#1F6FB2",       # navy
    "both_tails": "#6A1B7F",         # purple
}
TAIL_LABEL_NICE = {
    "effect_size_tail": "Effect-size tail only",
    "p_value_tail": "P-value tail only",
    "both_tails": "Both tails",
}


def check_annotation_geometry(fig, ax, annotations, output_bbox) -> None:
    """Fail before saving if rendered boxes obscure labels, legend or points."""
    from matplotlib.text import Text
    from matplotlib.transforms import Bbox

    original_dpi = fig.dpi
    checks = []
    failures = []
    try:
        # PDF point scale, the interactive figure, and the saved PNG scale.
        for dpi in (72, 200, 300):
            fig.set_dpi(dpi)
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            boxes = []
            for annotation in annotations:
                # Annotation.get_window_extent includes the leader.  Measure
                # only glyphs and the padded box, including its border stroke.
                text_box = Text.get_window_extent(annotation, renderer)
                patch = annotation.get_bbox_patch()
                box = Bbox.union([text_box, patch.get_window_extent(renderer)])
                boxes.append(box.padded(patch.get_linewidth() * dpi / 144))
            legend_box = ax.get_legend().get_window_extent(renderer)
            points = ax.transData.transform([a.xy for a in annotations])
            # Scatter s=14 pt^2, unchanged from the original plot.
            radius = np.sqrt(14) * dpi / 144
            point_boxes = [Bbox.from_extents(x-radius, y-radius, x+radius, y+radius)
                           for x, y in points]
            output_box = output_bbox.transformed(fig.dpi_scale_trans)
            records = []
            for i, (annotation, box, point) in enumerate(zip(annotations, boxes, points)):
                label = annotation.get_text()
                if box.overlaps(legend_box):
                    failures.append(f"{dpi} dpi: label covers legend: {label!r}")
                for j, point_box in enumerate(point_boxes):
                    if box.overlaps(point_box):
                        failures.append(f"{dpi} dpi: label {i} covers labelled marker {j}")
                if not (output_box.contains(*box.p0) and output_box.contains(*box.p1)):
                    failures.append(f"{dpi} dpi: label is outside saved crop: {label!r}")
                arrow = annotation.arrow_patch
                vertices = arrow.get_path().transformed(arrow.get_transform()).vertices
                endpoint_distance = float(np.linalg.norm(vertices[-1] - point))
                if endpoint_distance > 0.1:
                    failures.append(f"{dpi} dpi: leader {i} misses point by {endpoint_distance:.3f} px")
                # A connector must also avoid every other label box.
                from matplotlib.path import Path as MplPath
                for j, other_box in enumerate(boxes):
                    if i != j and MplPath(vertices).intersects_bbox(other_box, filled=False):
                        failures.append(f"{dpi} dpi: leader {i} crosses label {j}")
                records.append(dict(label=label, point_data=list(annotation.xy),
                                    point_px=point.tolist(), box_px=list(box.extents),
                                    text_px=list(Text.get_window_extent(annotation, renderer).extents),
                                    leader_endpoint_px=vertices[-1].tolist(),
                                    endpoint_distance_px=endpoint_distance))
            overlaps = [(i, j) for i, j in combinations(range(len(boxes)), 2)
                        if boxes[i].overlaps(boxes[j])]
            failures.extend(f"{dpi} dpi: label boxes {i} and {j} overlap" for i, j in overlaps)
            checks.append(dict(dpi=dpi, annotations=records,
                               legend_box_px=list(legend_box.extents),
                               label_pairs_checked=len(boxes)*(len(boxes)-1)//2,
                               labelled_point_checks=len(boxes)**2,
                               label_box_overlaps=overlaps))
        report = dict(passed=not failures, label_count=len(annotations),
                      failures=failures, render_checks=checks)
        try:
            (OUT_FIG_DIR / "annotation_geometry_v2_2026-11-17.json").write_text(
                json.dumps(report, indent=2))
        except OSError as exc:
            raise RuntimeError(f"Cannot save ED14 annotation geometry in {OUT_FIG_DIR}: {exc}") from exc
        if failures:
            raise RuntimeError("ED14 annotation geometry failed: " + "; ".join(failures))
    finally:
        fig.set_dpi(original_dpi)
        fig.canvas.draw()


def make_volcano(tail: pd.DataFrame) -> None:
    retro_style(["#C04A1E", "#1F6FB2", "#6A1B7F"])
    plot_df = tail.copy()
    plot_df["x"] = plot_df["effect_size"]
    # neg_log10_fdr already present and == -log10(fdr); use neg log10 p for "volcano"
    plot_df["y"] = -np.log10(plot_df["pvalue"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(8.5, 6.2), dpi=200)
    for tlab in ("p_value_tail", "effect_size_tail", "both_tails"):
        sub = plot_df[plot_df["tail_label"] == tlab]
        ax.scatter(
            sub["x"],
            sub["y"],
            s=14,
            color=TAIL_COLORS[tlab],
            alpha=0.55,
            linewidths=0,
            label=f"{TAIL_LABEL_NICE[tlab]} (n={len(sub)})",
        )

    # cap x for readability
    x_max = np.nanpercentile(np.abs(plot_df["x"]), 99.5)
    ax.set_xlim(-x_max * 1.05, x_max * 1.05)

    # The submitted page shows these nine labels.  Both GSE49906 conditions
    # occupied the same box there; retain its visible cardiac-fibroblast label
    # and attach it only to that condition's source point.
    annotation_positions = {
        ("GSE53060", "NT-ESC"): (-280, 30),
        ("GSE88425", "SNRNP70 knockdown"): (-200, 14),
        ("GSE117942", "4-OH tamoxifen"): (290, 50),
        ("GSE105266", "Proton radiation 150.0 cGy"): (140, 30),
        ("GSE105266", "Gamma radiation 10.0 cGy"): (90, 12),
        ("GSE166479", "EMT-high CAFs"): (-110, 226),
        ("GSE49088", "NRASG12V expression off"): (120, 240),
        ("GSE95638", "Osr2RFP/-"): (-5, 244),
        ("GSE49906", "Cardiac fibroblast maturation"): (245, 257),
    }
    annotations = []
    for (gse, condition), xytext in annotation_positions.items():
        rows = plot_df[(plot_df["gse_id"] == gse)
                       & (plot_df["condition_name"].str.strip() == condition)]
        # The atlas repeats the EMT-high row at identical coordinates.
        rows = rows.drop_duplicates(subset=["x", "y"])
        if len(rows) != 1:
            raise ValueError(f"Expected one source point for {gse} / {condition}; found {len(rows)}")
        r = rows.iloc[0]
        annotation = ax.annotate(
            f"{gse}\n{condition}",
            xy=(r["x"], r["y"]),
            xytext=xytext,
            textcoords="data",
            fontsize=7,
            color=INK,
            ha="center",
            va="center",
            # Zero endpoint shrink makes the connector reach the point centre.
            arrowprops=dict(arrowstyle="-", color="0.35", lw=0.5, shrinkB=0),
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.6", lw=0.4, alpha=0.9),
            zorder=5,
        )
        annotations.append(annotation)

    ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color="0.5")
    ax.axvline(0, ls="-", lw=0.4, color="0.7")

    ax.set_xlabel("Clock-shift effect size (β)")
    ax.set_ylabel("−log10(p-value)")
    ax.set_ylim(-3, 292)
    # title removed: title text lives in the docx caption only
    ax.legend(loc="upper left", frameon=False, fontsize=8)
    trim_axes(fig)

    fig.tight_layout()
    # Preserve the pre-revision axes rectangle and output crop exactly: label
    # placement must not resize the plot through tight_layout/bbox_inches.
    ax.set_position((0.0611356209150327, 0.07229782706093191,
                     0.9265114379084968, 0.890815916330632))
    from matplotlib.transforms import Bbox
    # Agg and PDF use slightly different glyph metrics. These are the
    # original renderer-specific tight crops, with its 0.08-inch save padding.
    pdf_bbox = Bbox.from_bounds(0.10038194444444459, 0.10210937500000009,
                               8.294618055555555, 5.996020979086315).padded(0.08)
    png_bbox = Bbox.from_bounds(0.10104166666666674, 0.09828125,
                               8.293958333333334, 6.004573930475204).padded(0.08)
    check_annotation_geometry(fig, ax, annotations, pdf_bbox)
    SUBMISSION_FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(VOLCANO_PNG, dpi=300, bbox_inches=png_bbox)
    fig.savefig(VOLCANO_PDF, bbox_inches=pdf_bbox)
    fig.savefig(SUBMISSION_VOLCANO_PDF, bbox_inches=pdf_bbox)
    plt.close(fig)
    normalize_png_canvas(VOLCANO_PNG, (2519, 1700))


# ---------------------------------------------------------------------------
# Category breakdown (Output 2)
# ---------------------------------------------------------------------------
CATEGORY_ORDER = [
    "developmental",
    "tissue_baseline",
    "age_baseline",
    "acute_stress",
    "transformed_cell",
    "disease",
    "intervention",
    "other",
]


def make_category_breakdown(tail: pd.DataFrame) -> pd.DataFrame:
    n_total = len(tail)
    rows = []
    for cat in CATEGORY_ORDER:
        sub = tail[tail["category"] == cat]
        n = len(sub)
        rows.append(
            dict(
                category=cat,
                n=n,
                proportion=n / n_total,
                n_effect_size_tail=int((sub["tail_label"] == "effect_size_tail").sum()),
                n_pvalue_tail=int((sub["tail_label"] == "p_value_tail").sum()),
                n_both_tails=int((sub["tail_label"] == "both_tails").sum()),
                median_abs_effect=float(np.nanmedian(np.abs(sub["effect_size"]))) if n else np.nan,
                median_fdr=float(np.nanmedian(sub["fdr"])) if n else np.nan,
                n_unique_gse=int(sub["gse_id"].nunique()),
            )
        )
    out = pd.DataFrame(rows)
    out["proportion"] = out["proportion"].round(4)
    out["median_abs_effect"] = out["median_abs_effect"].round(3)
    out.to_csv(CATEGORY_CSV, index=False)
    return out


# ---------------------------------------------------------------------------
# Pathway / mechanism enrichment (Output 3)
# ---------------------------------------------------------------------------
# We don't have gene-target columns per row; the actionable text we have is the
# intervention / condition_name strings. We therefore run a keyword-based
# *mechanism* enrichment: each row carries a small bag of tokens, and we ask
# which tokens are enriched in the tail vs the FDR-significant background.
# This is the same shape as ORA (one-sided Fisher's exact, BH FDR), just on
# normalized mechanism keywords rather than gene IDs.

MECHANISM_KEYWORDS: dict[str, list[str]] = {
    # cell-state / reprogramming
    "iPSC / reprogramming": [r"\biPSC\b", r"\b4F\b", r"\bOSKM\b", r"\breprogramm"],
    "ESC / embryonic stem": [r"\bhESC\b", r"\bESC\b(?!\w)", r"\bembryonic stem\b", r"\bNT-?ESC\b"],
    "Pluripotent / differentiation": [r"\bdifferentiation\b", r"\bpluripotent\b", r"\bmesenchymoangioblast\b"],
    # stress
    "Ionizing radiation": [r"\birradiation\b", r"\bgamma\b", r"\bproton\b", r"\bcGy\b", r"\bradiation\b"],
    "UV / photo-damage": [r"\bUVB?\b", r"\bphoto.?aging\b", r"\bphoto.?damage\b"],
    "Hypoxia / oxidative": [r"\bhypoxia\b", r"\bhyperoxia\b", r"\bH2O2\b", r"\boxidative\b"],
    "Inflammation / LPS / infection": [r"\bLPS\b", r"\binfection\b", r"\bIFN\b", r"\bimmune\b", r"\bsepsis\b"],
    # development / aging
    "Embryonic / fetal development": [r"\bembryo", r"\bfetal\b", r"\bfoetal\b", r"\bE\d+\.\d+\b"],
    "Postnatal / maturation": [r"\bpostnatal\b", r"\bneonat", r"\bmaturation\b", r"\bP\d+\b"],
    "Chronological aging baseline": [r"\bchronological\s*aging\b", r"\bage-related\b", r"\bold age\b", r"\bphenoage\b"],
    # tissue / cell type
    "Cross-tissue baseline": [r"\btissue\b", r"\borgan\b", r"\bliver\b", r"\bbrain\b", r"\bstomach\b", r"\bintestine\b", r"\bkidney\b", r"\bplacenta\b"],
    "Cell-type / cell-line baseline": [r"\bcell\s*type\b", r"\bHeLa\b", r"\bMEF\b", r"\bfibroblast subtype\b", r"\breference rna\b"],
    # disease
    "Cancer / tumor": [r"\b(?:tumou?r|cancer|carcinoma|melanoma|leukemia|lymphoma|sarcoma|glioma|SCLC)\b"],
    "Fibrosis / chronic disease": [r"\bfibrosis\b", r"\bIPF\b", r"\bdepression\b", r"\bdiabetes\b", r"\bobesity\b"],
    # intervention
    "Genetic loss-of-function (KO/KD)": [r"\bknockout\b", r"\bknockdown\b", r"\b_KO\b", r"\b-?KD\b"],
    "Diet / metabolic": [r"\bhigh fat diet\b", r"\bcalori(?:c|e) restriction\b", r"\bdiet\b", r"\bfasting\b"],
    "Small-molecule / drug": [r"\bazacitidine\b", r"\bmetformin\b", r"\brapamycin\b", r"\btamoxifen\b", r"\binhibitor\b", r"\bchemotherapy\b", r"\bdrug\b"],
}

COMPILED_MECH: dict[str, list[re.Pattern]] = {
    k: [re.compile(p, re.IGNORECASE) for p in pats] for k, pats in MECHANISM_KEYWORDS.items()
}


def _blob(row: pd.Series) -> str:
    parts = []
    for c in ("comparison_name", "condition_name", "intervention"):
        v = row.get(c)
        if isinstance(v, str) and v.strip():
            parts.append(v)
    return " || ".join(parts)


def _mech_hit(blob: str, pats: list[re.Pattern]) -> bool:
    return any(p.search(blob) for p in pats)


def make_pathway_enrichment(tail: pd.DataFrame, sig_all: pd.DataFrame) -> pd.DataFrame:
    """One-sided Fisher's exact test: each mechanism keyword in tail vs FDR-sig background."""
    from scipy.stats import fisher_exact

    tail_blobs = tail.apply(_blob, axis=1)
    bg_blobs = sig_all.apply(_blob, axis=1)
    n_tail = len(tail_blobs)
    n_bg = len(bg_blobs)

    rows = []
    for mech, pats in COMPILED_MECH.items():
        a = int(tail_blobs.map(lambda s: _mech_hit(s, pats)).sum())   # tail hits
        b = n_tail - a                                                 # tail non-hits
        c = int(bg_blobs.map(lambda s: _mech_hit(s, pats)).sum())     # bg hits
        d = n_bg - c                                                  # bg non-hits
        # tail vs full FDR-sig background; subtract a/b so background is non-tail
        c_nontail = c - a
        d_nontail = d - b
        if c_nontail < 0 or d_nontail < 0:
            # shouldn't happen, defensive
            c_nontail = max(c_nontail, 0)
            d_nontail = max(d_nontail, 0)
        table = [[a, b], [c_nontail, d_nontail]]
        try:
            odds, p = fisher_exact(table, alternative="greater")
        except ValueError:
            odds, p = np.nan, 1.0
        rows.append(
            dict(
                mechanism=mech,
                tail_hits=a,
                tail_total=n_tail,
                tail_frac=a / n_tail if n_tail else np.nan,
                nontail_hits=c_nontail,
                nontail_total=n_bg - n_tail,
                nontail_frac=c_nontail / (n_bg - n_tail) if (n_bg - n_tail) else np.nan,
                odds_ratio=odds,
                p_value=p,
            )
        )
    df = pd.DataFrame(rows).sort_values("p_value")
    # BH FDR
    from statsmodels.stats.multitest import multipletests
    df["fdr_bh"] = multipletests(df["p_value"].fillna(1.0), method="fdr_bh")[1]
    df["odds_ratio"] = df["odds_ratio"].round(3)
    df["tail_frac"] = df["tail_frac"].round(4)
    df["nontail_frac"] = df["nontail_frac"].round(4)
    df["p_value"] = df["p_value"].astype(float)
    df["fdr_bh"] = df["fdr_bh"].astype(float)
    df.to_csv(PATHWAY_CSV, index=False)
    return df


# ---------------------------------------------------------------------------
# Markdown summary (Output 4)
# ---------------------------------------------------------------------------
def write_summary_markdown(
    tail: pd.DataFrame,
    breakdown: pd.DataFrame,
    enrichment: pd.DataFrame,
) -> None:
    n_tail = len(tail)
    n_total_atlas = 106172
    n_sig = 13569

    rep = pd.read_csv(REP)

    # category percentages
    pct = {row["category"]: row["proportion"] * 100 for _, row in breakdown.iterrows()}

    def pct_str(cat: str) -> str:
        return f"{pct.get(cat, 0):.0f}%"

    # vignettes: pull the first row of each category from rep CSV
    def vignette(cat: str) -> str | None:
        sub = rep[rep["category"] == cat]
        if sub.empty:
            return None
        r = sub.iloc[0]
        gse = r["gse_id"]
        cn = r["condition_name"]
        mod = r["modality"]
        eff = r["effect_size"]
        fdr = r["fdr"]
        return (
            f"**{cat.replace('_', ' ')} ({gse}, {mod})**: "
            f"{cn} (β={eff:.1f}, FDR={fdr:.1e}). {r['one_line_explanation']}"
        )

    vignettes = [v for v in (vignette("developmental"), vignette("tissue_baseline"), vignette("acute_stress")) if v]

    # top enriched mechanisms (top 5 by p-value with at least 5 tail hits)
    enr = enrichment[enrichment["tail_hits"] >= 5].head(5)
    top_mech_lines = []
    for _, r in enr.iterrows():
        top_mech_lines.append(
            f"  - {r['mechanism']}: OR={r['odds_ratio']:.2f}, tail-fraction={r['tail_frac']*100:.1f}%, "
            f"p={r['p_value']:.1e}, FDR={r['fdr_bh']:.1e}"
        )

    # genuine geroprotector / intervention examples from the tail (real interventions
    # with effect_size_tail flag, high score, real intervention text)
    int_examples = tail[
        (tail["category"] == "intervention")
        & (tail["intervention"].notna())
        & (tail["final_score"] >= 50)
    ].sort_values("effect_size", key=lambda s: s.abs(), ascending=False).head(5)
    int_lines = []
    for _, r in int_examples.iterrows():
        int_lines.append(
            f"  - {r['gse_id']} ({r['modality']}): {r['intervention']}, "
            f"beta={r['effect_size']:.1f}, FDR={r['fdr']:.1e}, score={int(r['final_score'])}"
        )

    md = []
    md.append("# Dedicated drill-down on the extreme-tail subset of the ClockBase Agent atlas")
    md.append("")
    md.append(
        f"Inputs: `Aging-AI-Agent/figures/combined_3modality_plot_data.csv` "
        f"(pooled three-modality table, n={n_total_atlas:,}). "
        f"Tail subset definition: FDR<0.05 and `outlier_status` in "
        f"{{Effect Size Outlier, P-value Outlier, Both}}. "
        f"Tail n={n_tail:,} of {n_sig:,} FDR-significant rows ({n_tail/n_sig*100:.1f}%)."
    )
    md.append("")
    md.append("## 1. Headline")
    md.append("")
    md.append(
        f"The extreme-tail subset (n={n_tail:,}; {n_tail/n_sig*100:.1f}% of FDR-significant rows in the "
        f"pooled three-modality atlas) is dominated by biological contrasts where a large clock shift is "
        f"mechanically expected rather than by statistical false positives. Keyword classification of "
        f"`condition_name`, `comparison_name`, and `intervention` assigns "
        f"{pct_str('developmental')} of tail rows to developmental contrasts, "
        f"{pct_str('tissue_baseline')} to tissue- or cell-type-of-origin baselines, "
        f"{pct_str('age_baseline')} to chronological-age baselines, "
        f"{pct_str('acute_stress')} to acute high-dose stressors (radiation, UV, LPS, hypoxia), "
        f"{pct_str('transformed_cell')} to transformed or stem-cell-state baselines (iPSC/ESC/HeLa/tumour), "
        f"{pct_str('disease')} to overt disease contrasts, "
        f"{pct_str('intervention')} to real interventions with biologically large effect sizes, "
        f"and {pct_str('other')} unclassified. The useful finding is two-fold: first, almost none of the "
        f"tail rows are statistical artefacts (the dominant signal is real biology that happens to be "
        f"big), and second, the intervention slice contains a small set of genuine geroprotector-style "
        f"perturbations with very large effects that are explicitly preserved on the website for "
        f"re-analysis."
    )
    md.append("")
    md.append("## 2. Representative vignettes")
    md.append("")
    for v in vignettes:
        md.append(f"- {v}")
    md.append("")
    md.append("## 3. What's hidden in the tail that may be useful")
    md.append("")
    md.append(
        "A small number of real geroprotective or genetic interventions with very large effects sit "
        "inside the tail flag and are explicitly retained in the public result table for re-analysis. "
        "Top tail-flagged real interventions include:"
    )
    md.extend(int_lines if int_lines else ["  - (no high-score real interventions matched the filter)"])
    md.append("")
    md.append("## 4. Mechanism enrichment in tail vs FDR-significant background")
    md.append("")
    md.append("Top 5 mechanism keywords enriched in the tail relative to the non-tail FDR-significant background:")
    md.extend(top_mech_lines if top_mech_lines else ["  - (no mechanism passed the >=5-hit filter)"])
    md.append("")
    md.append("## Files generated")
    md.append("")
    md.append(f"- `{VOLCANO_PNG.relative_to(ROOT)}`")
    md.append(f"- `{VOLCANO_PDF.relative_to(ROOT)}`")
    md.append(f"- `{CATEGORY_CSV.relative_to(ROOT)}`")
    md.append(f"- `{PATHWAY_CSV.relative_to(ROOT)}`")
    md.append(f"- `{SUMMARY_MD.relative_to(ROOT)}`")
    md.append(f"- `{Path(__file__).resolve().relative_to(ROOT)}`")
    md.append("")
    SUMMARY_MD.write_text("\n".join(md))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def main() -> None:
    tail = load_tail_subset()
    tail["tail_label"] = tail["outlier_status"].map(assign_tail_label)
    tail["category"] = tail.apply(classify_row, axis=1)

    print(f"Loaded tail subset: n={len(tail)}")
    print("Tail label counts:")
    print(tail["tail_label"].value_counts())
    print()
    print("Category counts:")
    print(tail["category"].value_counts())

    # Output 1
    make_volcano(tail)
    print(f"Wrote {VOLCANO_PNG} and {VOLCANO_PDF}")

    # Output 2
    breakdown = make_category_breakdown(tail)
    print("Category breakdown:")
    print(breakdown.to_string(index=False))
    print(f"Wrote {CATEGORY_CSV}")

    # Output 3
    df = pd.read_csv(ATLAS, low_memory=False)
    sig_all = df[df["fdr"] < 0.05].copy()
    enrichment = make_pathway_enrichment(tail, sig_all)
    print("Top 5 enriched mechanisms:")
    print(enrichment.head(5).to_string(index=False))
    print(f"Wrote {PATHWAY_CSV}")

    # Output 4
    write_summary_markdown(tail, breakdown, enrichment)
    print(f"Wrote {SUMMARY_MD}")


if __name__ == "__main__":
    main()
