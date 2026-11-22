#!/usr/bin/env python3
"""Targeted model reconstruction from saved predictions and observed metadata.

Outputs are internal reconstruction results, not replacements for atlas ranks.
No clock predictions, methylation processing, API calls, or source edits occur.
Run with the repository venv/bin/python.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess

import numpy as np
import pandas as pd
import scipy
from scipy.stats import ttest_rel, shapiro, chi2
import statsmodels
from statsmodels.regression.linear_model import OLS
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
CLOCKS = {
    "Horvathv1": "clockresult.Horvathv1",
    "PhenoAge": "clockresult.PhenoAge",
    "YingDamAge": "clockresult.YingDamAge",
    "DunedinPACE": "clockresult.DunedinPACE",
    "GrimAgeV2": "clockresult.DNAmGrimAge.1",
}
MISSING = {k: k + "__missing_frac" for k in CLOCKS}


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fit(y, design, term):
    values = design.to_numpy(dtype=float)
    if not np.isfinite(values).all() or not np.isfinite(y).all():
        raise ValueError("Nonfinite model input")
    rank = int(np.linalg.matrix_rank(values))
    if rank != design.shape[1] or len(y) <= rank:
        raise ValueError("Rank-deficient design or no residual degrees of freedom")
    model = OLS(np.asarray(y, dtype=float), design).fit(cov_type="nonrobust")
    ci = model.conf_int(alpha=0.05).loc[term].tolist()
    return {
        "coefficient": float(model.params[term]), "se": float(model.bse[term]),
        "ci95": ci, "pvalue": float(model.pvalues[term]),
        "n": int(model.nobs), "rank": rank, "residual_df": int(model.df_resid),
        "condition_number": float(np.linalg.cond(values)),
        "residual_sum_squares": float(model.ssr),
        "max_leverage": float(model.get_influence().hat_matrix_diag.max()),
        "design_columns": list(design),
    }


def run():
    design_path = ROOT / "docs/2026-09-07-down-syndrome-design.json"
    trace_path = ROOT / "docs/2026-09-07-down-syndrome-execution-traces.json"
    source_path = ROOT / "hs_dnam_sonnet37/GSE140344/aging_analysis_sonnet37.csv"
    design = json.loads(design_path.read_text())
    export = Path(design["input"]["path"])
    if sha(export) != design["input"]["sha256"]:
        raise ValueError("Sample export changed since raw metadata review")
    wanted = {"gse_id", "gsm_id", "DNAmGrimAge2", "has_prediction", "is_imputed",
              *CLOCKS.values(), *MISSING.values()}
    chunks = []
    for chunk in pd.read_csv(export, chunksize=50000, usecols=lambda c: c in wanted):
        chunks.append(chunk[chunk.gse_id.isin(design["studies"])])
    predictions = pd.concat(chunks, ignore_index=True)
    if predictions.duplicated(["gse_id", "gsm_id"]).any():
        raise ValueError("Duplicate prediction identity")

    # Direct check against the five sample values printed in the retained trace.
    trace = json.loads(trace_path.read_text())
    lines = trace["studies"]["GSE140344"]["messages"][34]["text"].splitlines()
    header = next(line.split() for line in lines if line.strip().startswith("gse_id  Horvathv1"))
    preview = [line.split() for line in lines if line.startswith("GSM41595")]
    if len(preview) != 5:
        raise ValueError("Expected five trace sample previews")
    mapping_checks = []
    for clock, column in CLOCKS.items():
        delta = []
        for tokens in preview:
            recorded = float(tokens[header.index(clock) + 1])
            matches = predictions[(predictions.gse_id == "GSE140344") & (predictions.gsm_id == tokens[0])]
            if len(matches) != 1:
                raise ValueError("Trace sample does not map uniquely")
            delta.append(abs(float(matches.iloc[0][column]) - recorded))
        if max(delta) > 0.00000051:
            raise ValueError(f"{clock} does not match six-decimal trace output")
        mapping_checks.append({"clock": clock, "export_column": column,
                               "sample_gsms": [tokens[0] for tokens in preview],
                               "max_absolute_difference": max(delta)})

    data = {}
    qc = {}
    for gse, study in design["studies"].items():
        metadata = pd.DataFrame(study["samples"])
        frame = metadata.merge(predictions[predictions.gse_id == gse], on="gsm_id", validate="one_to_one")
        if len(frame) != study["n_samples"] or not np.isfinite(frame[list(CLOCKS.values())].to_numpy()).all():
            raise ValueError(f"Missing samples or nonfinite predictions for {gse}")
        frame["female"] = frame.raw_gender.str.lower().map({"female": 1.0, "male": 0.0})
        if frame.female.isna().any():
            raise ValueError(f"Missing sex coding for {gse}")
        data[gse] = frame
        qc[gse] = {"n": len(frame), "unique_gsms": frame.gsm_id.nunique(),
                   "missing_feature_fractions": {clock: {
                       "min": float(frame[column].min()), "max": float(frame[column].max()),
                       "nonmissing": int(frame[column].notna().sum())}
                       for clock, column in MISSING.items()},
                   "has_prediction_counts": {str(k): int(v) for k, v in frame.has_prediction.value_counts(dropna=False).items()},
                   "is_imputed_counts": {str(k): int(v) for k, v in frame.is_imputed.value_counts(dropna=False).items()}}

    ds = data["GSE52588"].query("arm == 'DS'").set_index("raw_family").sort_index()
    sib = data["GSE52588"].query("arm == 'sibling'").set_index("raw_family").sort_index()
    if len(ds) != 29 or not ds.index.is_unique or not sib.index.is_unique or not ds.index.equals(sib.index):
        raise ValueError("Expected 29 complete DS-sibling pairs")
    difference_design = pd.DataFrame({"const": 1.0, "delta_age": ds.raw_age-sib.raw_age,
                                      "delta_female": ds.female-sib.female}, index=ds.index)
    combined = pd.concat([ds.assign(DS=1.0), sib.assign(DS=0.0)]).reset_index()
    fixed_design = pd.DataFrame({"const": 1.0, "DS": combined.DS,
                                "age_centered": combined.raw_age-combined.raw_age.mean(),
                                "female": combined.female})
    fixed_design = pd.concat([fixed_design, pd.get_dummies(combined.raw_family, prefix="family", drop_first=True, dtype=float)], axis=1)
    results = []
    for clock, column in CLOCKS.items():
        primary = fit((ds[column]-sib[column]).to_numpy(), difference_design, "const")
        differences = (ds[column]-sib[column]).to_numpy()
        model = OLS(differences, difference_design).fit()
        robust = model.get_robustcov_results(cov_type="HC3", use_t=True)
        influence = model.get_influence()
        auxiliary = OLS(np.square(model.resid), difference_design).fit()
        lm = float(len(ds) * auxiliary.rsquared)
        diagnostics = {
            "role": "Post-primary diagnostics and secondary sensitivity; primary inference unchanged",
            "shapiro_residual_pvalue": float(shapiro(model.resid).pvalue),
            "koenker_breusch_pagan_lm": lm,
            "koenker_breusch_pagan_df": 2,
            "koenker_breusch_pagan_pvalue": float(chi2.sf(lm, 2)),
            "hc3_secondary": {"coefficient": float(robust.params[0]),
                              "se": float(robust.bse[0]), "ci95": robust.conf_int()[0].tolist(),
                              "pvalue": float(robust.pvalues[0]), "inference": "two-sided t, residual df 26"},
            "family_influence": [{"family": family, "residual": float(model.resid.iloc[i]),
                                  "leverage": float(influence.hat_matrix_diag[i]),
                                  "cooks_distance": float(influence.cooks_distance[0][i]),
                                  "external_studentized_residual": float(influence.resid_studentized_external[i])}
                                 for i, family in enumerate(ds.index)],
            "leave_one_family_out": [{"family": family,
                                       "model": fit(np.delete(differences, i), difference_design.drop(index=family), "const")}
                                      for i, family in enumerate(ds.index)],
        }
        crosscheck = fit(combined[column].to_numpy(), fixed_design, "DS")
        for key in ("coefficient", "se", "pvalue"):
            if not np.isclose(primary[key], crosscheck[key], rtol=1e-8, atol=1e-10):
                raise ValueError(f"Paired and family-fixed-effect models disagree: {clock} {key}")
        descriptive = ttest_rel(ds[column], sib[column])
        results.append({"study": "GSE52588", "contrast": "DS_vs_sibling", "clock": clock,
                        "units": "pace units" if clock == "DunedinPACE" else "clock years",
                        "primary": primary, "diagnostics": diagnostics, "family_fixed_effect_crosscheck": crosscheck,
                        "unadjusted_paired_secondary": {"mean_difference": float((ds[column]-sib[column]).mean()),
                                                       "pvalue": float(descriptive.pvalue)}})
    q = multipletests([r["primary"]["pvalue"] for r in results], method="fdr_bh")[1]
    for row, value in zip(results, q):
        row["local_five_clock_bh_q"] = float(value)

    replays = []
    for name, arms in (("all_T21_vs_control", ["T21_IQ+", "T21_IQ-", "control"]),
                       ("IQplus_vs_control", ["T21_IQ+", "control"])):
        frame = data["GSE140344"][data["GSE140344"].arm.isin(arms)].copy().reset_index(drop=True)
        matrix = pd.DataFrame({"const": 1.0, "DS": (frame.arm != "control").astype(float),
                               "age_centered": frame.raw_age-frame.raw_age.mean(), "female": frame.female})
        for clock, column in CLOCKS.items():
            replays.append({"study": "GSE140344", "contrast": name, "clock": clock,
                            "units": "pace units" if clock == "DunedinPACE" else "clock years",
                            "model": fit(frame[column].to_numpy(), matrix, "DS"),
                            "gsm_ids": frame.gsm_id.tolist()})

    source = pd.read_csv(source_path)
    comparisons = []
    names = {"all_T21_vs_control": "Trisomy21_vs_Control", "IQplus_vs_control": "Trisomy21_IQplus_vs_Control"}
    for replay in replays:
        matched = source[(source.comparison_name == names[replay["contrast"]]) & (source.clock == replay["clock"])]
        if len(matched) != 1:
            raise ValueError("Source comparison and clock must map uniquely")
        pvalue = float(matched.iloc[0].p_value)
        reconstructed = replay["model"]["pvalue"]
        matches = bool(np.isclose(pvalue, reconstructed, rtol=1e-10, atol=0))
        if not matches:
            raise ValueError("Source p-value replay differs beyond floating-point tolerance")
        comparisons.append({"comparison": names[replay["contrast"]], "clock": replay["clock"],
                            "source_csv_line": int(matched.index[0]) + 2,
                            "source_pvalue": pvalue, "reconstructed_pvalue": reconstructed,
                            "absolute_difference": abs(pvalue-reconstructed), "matches": matches})

    return {"status": "Internal reconstruction; not incorporated into atlas or submission",
            "inference": "Conventional homoskedastic OLS, two-sided t inference, 95% confidence intervals, specified before fitting",
            "primary_formula": "DS-minus-sibling clock difference ~ intercept + age difference + female difference",
            "multiplicity": "BH across exactly five GSE52588 primary tests only; not atlas FDR. GSE140344 results are unadjusted-p diagnostic replays.",
            "limitations": "Existing finite predictions have nonzero clock-feature missing fractions and is_imputed=0; these fields do not establish upstream missing-feature handling. Prediction validity, scores and atlas ranks remain unresolved. Mothers excluded from inference. Diagnostics have limited power with 29 families and do not prove assumptions.",
            "model_assumptions": ["Linear age adjustment", "Additive sex effect", "Independent families", "Correct observed family identifiers", "Primary OLS homoskedastic errors; HC3 is secondary"],
            "effect_schema": "Adjusted signed clock-unit coefficients are not the historical absolute unadjusted standardized effect_size; no automatic substitution.",
            "sources": [{"path": str(p), "sha256": sha(p), "mtime_ns": p.stat().st_mtime_ns} for p in (export, design_path, trace_path, source_path)],
            "execution_provenance": {"script_sha256": sha(Path(__file__)), "checkout_head_at_run": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(), "script_sha_is_authoritative_for_uncommitted_edits": True},
            "source_pvalue_comparisons": comparisons,
            "software": {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__, "statsmodels": statsmodels.__version__},
            "clock_mapping_preview_checks": mapping_checks, "prediction_qc": qc,
            "family_pairs": [{"family": family, "DS_gsm": ds.loc[family, "gsm_id"], "sibling_gsm": sib.loc[family, "gsm_id"]} for family in ds.index],
            "GSE52588_results": results, "GSE140344_replays": replays}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = run()
        content = json.dumps(result, indent=2, allow_nan=False) + "\n"
        with args.output.open("x") as handle:
            handle.write(content)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(1, f"Targeted reconstruction failed: {exc}\n")
