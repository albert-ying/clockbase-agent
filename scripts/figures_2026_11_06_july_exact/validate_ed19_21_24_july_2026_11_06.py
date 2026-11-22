"""Check the delivered Extended Data pages against their supplied sources."""
from __future__ import annotations

import os
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pymupdf as fitz
from PIL import Image

REPO = Path(__file__).resolve().parents[2]
EVIDENCE = REPO / "paper/revision_2026-09-07/evidence"
OLD = EVIDENCE / "figures-rebuild-2026-11-05"
NEW = EVIDENCE / "figures-rebuild-2026-11-06"


def digest(path):
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        raise RuntimeError(f"Cannot read required file {path}: {exc}") from exc


def validate_page(number):
    name = f"ed{number}"
    folder = NEW / name
    for filename in (f"{name}.pdf", f"{name}_review.png",
                     f"{name}_vs_july.png", "README.md"):
        path = folder / filename
        assert path.is_file() and path.stat().st_size, str(path)
    with fitz.open(folder / f"{name}.pdf") as doc:
        assert len(doc) == 1, name
        page = doc[0]
        width_mm = page.rect.width * 25.4 / 72
        # User permits either 183 mm or the July assembled page width.
        assert min(abs(width_mm - w) for w in (183, 210.142667, 215.9)) < .02, width_mm
        assert len(page.get_drawings()) > 20, name
        assert not page.get_images(), f"{name} contains raster plot images"
        spans = [s for b in page.get_text("dict")["blocks"]
                 for line in b.get("lines", []) for s in line["spans"]
                 if s["text"].strip()]
        assert spans, name
        assert min(s["size"] for s in spans) >= 6 - .015, (
            name, [(s["text"], s["size"]) for s in spans if s["size"] < 5.985])
        clipped = [s["text"] for s in spans if not
                   (page.rect + (-.5, -.5, .5, .5)).contains(fitz.Rect(s["bbox"]))]
        assert not clipped, (name, clipped)
        text = page.get_text()
        assert "Extended Data Figure" not in text, name
        assert "computed in" not in text.lower() and "data:" not in text.lower(), name
        info = {
            "width_mm": width_mm,
            "height_mm": page.rect.height * 25.4 / 72,
            "vector_paths": len(page.get_drawings()),
            "raster_images": len(page.get_images()),
            "minimum_pdf_span_size_pt": min(s["size"] for s in spans),
            "fonts": sorted({s["font"] for s in spans}),
            "pdf_sha256": digest(folder / f"{name}.pdf"),
        }
    with Image.open(folder / f"{name}_review.png") as png:
        assert abs(png.info.get("dpi", (0,))[0] - 200) < .1
    # Conserved tables remain byte-identical; ED21b is recomputed from the current census.
    if number != 19:
        unchanged = []
        for source in sorted((OLD / name).glob("*.csv")):
            target = folder / source.name
            assert target.is_file(), f"Missing conservation table: {target}"
            if number == 21 and source.name == "panel_b_defects.csv":
                census_path = Path(os.environ["CLOCKBASE_FINAL"]) / "error-census/census.json"
                census = json.loads(census_path.read_text())["modalities"]
                provenance = json.loads((folder / "provenance.json").read_text())
                record = next(r for r in provenance if Path(r["path"]).resolve() == census_path.resolve())
                assert record["sha256"] == digest(census_path), "ED21 census provenance is stale"
                defects = pd.read_csv(target)
                old_codes = set(pd.read_csv(source).defect)
                expected_codes = old_codes | {"ARMS_DEFINED_BY_CLOCK_OUTPUT", "ARMS_DO_NOT_MATCH_GEO_SAMPLE_RECORDS", "CLOCK_OUTPUTS_OUTSIDE_PLAUSIBLE_RANGE", "SINGLE_CELLS_AS_INDEPENDENT_SAMPLES"}
                expected_keys = {(m, code) for m in ("dnam", "hs_rna", "mm_rna") for code in expected_codes}
                assert set(zip(defects.modality, defects.defect)) == expected_keys
                assert len(defects) == len(expected_keys), "Duplicate ED21 defect rows"
                for row in defects.itertuples():
                    expected = census[row.modality]["defects"][row.defect]
                    assert row.n == expected["n"] and row.of == expected["of"]
                    assert np.isclose(row.fraction, expected["fraction"], rtol=0, atol=1e-12)
                    assert np.isclose(row.percent, 100 * expected["n"] / expected["of"], rtol=0, atol=1e-12)
                info["current_census_verified_table"] = source.name
                continue
            assert digest(source) == digest(target), f"Changed source values: {target}"
            unchanged.append(source.name)
        info["byte_identical_original_tables"] = unchanged
    return info


def main():
    summary = {f"ed{n}": validate_page(n) for n in (19, 21, 22, 23, 24)}
    funnel = pd.read_csv(OLD / "figure3/panel_f_values.csv")
    final = Path(os.environ['CLOCKBASE_FINAL'])
    ranking = json.loads((final / 'comparison-ranking/summary.json').read_text())
    stages = ('usable_significant_n3_rankable_category', 'agreed_matching_contrast_type',
              'no_agreed_adjudication_design_concern', 'context_admitted')
    assert funnel.comparisons.tolist() == [ranking['funnel'][k]['comparisons'] for k in stages] + [ranking['outputs']['all']['ranked']]
    replication = pd.read_csv(OLD / "figure4/panel_i_replication_points.csv")
    robust = replication[replication.robust.astype(str).str.lower().eq("true")]
    meta = Path(os.environ.get('CLOCKBASE_META') or final.parent / 'intervention-meta-2026-10-17')
    primary = pd.read_csv(meta / 'primary/intervention_meta.csv')
    expected = primary[primary.tier.eq('replicated')]
    assert len(replication) == len(expected)
    expected_robust = expected[expected.fdr.lt(.05) & expected.rank_bound.gt(0)]
    assert (robust.estimate > 0).sum() == (expected_robust.estimate > 0).sum()
    assert (robust.estimate < 0).sum() == (expected_robust.estimate < 0).sum()
    named = pd.read_csv(OLD / "figure4/panel_i_named_points.csv")
    assert set(named.canonical_target) == {
        "lipopolysaccharide", "sleep deprivation", "plx5622", "bleomycin"}
    # The copied existing example values retain their tested mean differences.
    stats = pd.read_csv(NEW / "ed24/example_statistics.csv")
    values = pd.read_csv(NEW / "ed24/plotted_sample_values.csv")
    for row in stats.itertuples():
        data = values[values.panel == row.panel]
        case = data[data.group == "case"].clock_value
        control = data[data.group == "control"].clock_value
        assert len(case) == row.n_case and len(control) == row.n_control
        assert np.isclose(case.mean() - control.mean(), row.mean_diff, atol=1e-6)
    summary["ed21"]["added_funnel_values"] = funnel.comparisons.tolist()
    summary["ed22"]["added_volcano_counts"] = {
        "interventions": len(replication), "robust_higher": int((robust.estimate > 0).sum()), "robust_lower": int((robust.estimate < 0).sum())}
    # Keep independent verification with each requested page.
    for name, report in summary.items():
        (NEW / name / "independent_validation.json").write_text(
            json.dumps({"passed": True, **report}, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        raise SystemExit(f"Extended Data validation failed: {exc}") from exc
