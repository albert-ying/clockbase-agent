"""Independent base-R verification of source-bound conditional calculations.

No producer, pandas, NumPy or SciPy code is imported. Source sample/clock values
are read directly; pairing and covariate values are treated as declared designs,
whose biological/source validity requires separate review. A numerical PASS is
never scientific acceptance. Constructed software tests belong in the test file.
"""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import time


SCHEMA = "clockbase-independent-statistics-v1"
POOLED_DEFINITION = "case-minus-control divided by pooled within-group sample SD"
PAIRED_DEFINITION = "mean within-unit case-minus-control divided by SD of within-unit differences"
OLS_DEFINITION = "unstandardized case coefficient; standardized effect not defined"

# Fixed reference code. User-supplied strings are data, never R expressions.
R_REFERENCE = r'''
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("expected job table and output directory")
jobs <- read.csv(args[[1L]], check.names=FALSE, stringsAsFactors=FALSE,
                 na.strings=character())
outdir <- args[[2L]]
writeLines(capture.output(sessionInfo()), file.path(outdir, "r-session.txt"))
outcomes <- list()
for (j in seq_len(nrow(jobs))) {
  job <- jobs[j,,drop=FALSE]
  tryCatch({
    d <- read.csv(job$data_file, check.names=FALSE, stringsAsFactors=FALSE,
                  colClasses="character", na.strings=character())
    d$value <- as.numeric(d$value)
    for (column in grep("^cov_[0-9]+$",names(d),value=TRUE)) d[[column]] <- as.numeric(d[[column]])
    if (anyDuplicated(d$sample_id) || any(!is.finite(d$value))) stop("invalid sample identities or values")
    confidence <- as.numeric(job$confidence)
    if (!is.finite(confidence) || confidence <= 0 || confidence >= 1) stop("invalid confidence level")
    a <- d$value[d$arm == "case"]
    b <- d$value[d$arm == "control"]
    if (length(a) == 0L || length(b) == 0L || length(a)+length(b) != nrow(d)) stop("invalid arms")
    answer <- list(case_mean=mean(a), control_mean=mean(b), mean_difference=mean(a)-mean(b),
                   n_case=length(a), n_control=length(b))
    mode <- job$mode
    if (mode == "descriptive") {
      # Assay means only. No standardized effect, P value or inferential interval.
    } else if (mode %in% c("paired", "welch", "student")) {
      if (min(length(a),length(b)) < 2L) stop("too few observations for t inference")
      pv <- ((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2)
      if (!is.finite(pv) || pv <= 0) stop("degenerate pooled variance")
      pooled_effect <- (mean(a)-mean(b))/sqrt(pv)
      if (mode == "paired") {
        if (any(!nzchar(d$unit_id)) || anyDuplicated(paste(d$arm,d$unit_id,sep="\r"))) stop("duplicate or blank pair units")
        units <- sort(unique(d$unit_id))
        aa <- d[d$arm == "case",,drop=FALSE]; bb <- d[d$arm == "control",,drop=FALSE]
        if (!identical(sort(aa$unit_id),units) || !identical(sort(bb$unit_id),units)) stop("incomplete pairs")
        a <- aa$value[match(units,aa$unit_id)]; b <- bb$value[match(units,bb$unit_id)]
        fit <- t.test(a,b,paired=TRUE,alternative="two.sided",conf.level=confidence)
        delta <- a-b
        se <- sd(delta)/sqrt(length(delta))
        answer$mean_difference <- mean(delta)
        answer$effect_size <- mean(delta)/sd(delta)
        answer$pooled_standardized_difference <- mean(delta)/sqrt(pv)
        answer$n_biological_units <- length(delta)
        answer$within_unit_difference <- setNames(delta,units)
      } else {
        fit <- t.test(a,b,paired=FALSE,var.equal=(mode == "student"),
                      alternative="two.sided",conf.level=confidence)
        se <- if (mode == "student") sqrt(pv*(1/length(a)+1/length(b))) else sqrt(var(a)/length(a)+var(b)/length(b))
        answer$effect_size <- pooled_effect
        answer$case_variance <- var(a); answer$control_variance <- var(b)
      }
      answer$p_value <- unname(fit$p.value)
      answer$standard_error <- se
      answer$statistic <- unname(fit$statistic)
      answer$degrees_of_freedom <- unname(fit$parameter)
      answer$ci_low <- fit$conf.int[[1L]]; answer$ci_high <- fit$conf.int[[2L]]
    } else if (mode %in% c("ols", "ols_hc3")) {
      covnames <- grep("^cov_[0-9]+$",names(d),value=TRUE)
      if (length(covnames) == 0L) stop("adjusted model requires declared covariates")
      x <- cbind(intercept=1, case=as.numeric(d$arm == "case"), as.matrix(d[,covnames,drop=FALSE]))
      storage.mode(x) <- "double"
      if (any(!is.finite(x)) || nrow(x) <= ncol(x) || qr(x,tol=1e-10)$rank != ncol(x)) stop("nonfinite, saturated or rank-deficient design")
      if (kappa(x,exact=TRUE) > 1e10) stop("ill-conditioned design")
      fit <- lm.fit(x,d$value,tol=1e-10)
      df <- nrow(x)-ncol(x)
      bread <- solve(crossprod(x))
      leverage <- rowSums((x %*% bread)*x)
      if (mode == "ols_hc3") {
        if (any(1-leverage <= 1e-10)) stop("unit leverage prevents HC3 covariance")
        covariance <- bread %*% crossprod(x,x*(fit$residuals/(1-leverage))^2) %*% bread
      } else covariance <- sum(fit$residuals^2)/df*bread
      beta <- unname(fit$coefficients[[2L]])
      se <- sqrt(covariance[2L,2L])
      if (!is.finite(se) || se <= 0) stop("degenerate adjusted variance")
      critical <- qt((1+confidence)/2,df)
      answer$adjusted_mean_difference <- beta
      answer$standard_error <- se
      answer$statistic <- beta/se
      answer$degrees_of_freedom <- df
      answer$residual_degrees_of_freedom <- df
      answer$p_value <- 2*pt(-abs(beta/se),df)
      answer$ci_low <- beta-critical*se; answer$ci_high <- beta+critical*se
      answer$coefficients <- unname(fit$coefficients)
      answer$residuals <- unname(fit$residuals)
      answer$leverage <- unname(leverage)
    } else stop("unsupported independent method")
    records <- list()
    for (field in names(answer)) {
      values <- answer[[field]]
      if (any(!is.finite(values))) stop(paste("nonfinite independent result",field))
      indexes <- names(values)
      if (is.null(indexes)) indexes <- as.character(seq_along(values))
      records[[length(records)+1L]] <- data.frame(field=field,index=indexes,
          value=sprintf("%.17g",as.numeric(values)),stringsAsFactors=FALSE)
    }
    write.csv(do.call(rbind,records),file.path(outdir,paste0(job$job_id,"-reference.csv")),row.names=FALSE,quote=TRUE)
    outcomes[[length(outcomes)+1L]] <- data.frame(job_id=job$job_id,status="PASS",message="",stringsAsFactors=FALSE)
  },error=function(e) {
    outcomes[[length(outcomes)+1L]] <<- data.frame(job_id=job$job_id,status="CANNOT VERIFY",message=conditionMessage(e),stringsAsFactors=FALSE)
  })
}
write.csv(do.call(rbind,outcomes),file.path(outdir,"r-outcomes.csv"),row.names=FALSE,quote=TRUE)
'''


def _sha(path):
    digest = hashlib.sha256()
    try:
        with Path(path).open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
    except OSError as error:
        raise OSError(f"Cannot hash independent-verification file {path}: {error}") from error
    return digest.hexdigest()


def _json_bytes(value):
    return (json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, indent=2) + "\n").encode()


def _write_json(path, value):
    try:
        with Path(path).open("xb") as handle:
            handle.write(_json_bytes(value))
    except OSError as error:
        raise OSError(f"Cannot preserve independent-verification output {path}: {error}") from error


def _number(value):
    if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value):
        raise ValueError("Expected a finite numerical result")
    return float(value)


def _csv(path):
    try:
        with Path(path).open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.reader(handle, strict=True))
    except (OSError, UnicodeError, csv.Error) as error:
        raise ValueError(f"Cannot parse prediction CSV: {error}") from error
    if not rows or not rows[0] or rows[0][0].startswith("version https://git-lfs.github.com/spec/v1"):
        raise ValueError("Prediction CSV is empty or contains an unresolved LFS pointer")
    if not rows[0] or any(len(row) != len(rows[0]) for row in rows[1:]):
        raise ValueError("Prediction CSV has an empty header or ragged records")
    header, records = rows[0], rows[1:]
    positions = [index for index, name in enumerate(header) if name in {"sample_id", "gsm_id", "GSM", "geo_accession"}]
    binding = None
    if header[0] == "":
        positions.insert(0, 0)
        binding = {"index": 0, "header": "", "format": "GSM accession",
                   "scope": "explicit full-column canonical-GSM/uniqueness check; no row-order identity inference"}
    if not positions:
        raise ValueError("Prediction identity column is absent")
    by_id = {}
    for record in records:
        identifiers = {record[index] for index in positions}
        if len(identifiers) != 1:
            raise ValueError("Prediction sample identity columns conflict")
        identifier = next(iter(identifiers))
        if not re.fullmatch(r"GSM[0-9]+", identifier, flags=re.ASCII):
            raise ValueError("Prediction identities must be canonical ASCII GSM accessions")
        if identifier in by_id:
            raise ValueError("Prediction sample identity is duplicated")
        by_id[identifier] = record
    repeated = [{"finding": "repeated_identity_header", "header": name,
                 "positions": [index for index in positions if header[index] == name],
                 "position_index_base": 0, "provenance_status": "UNRESOLVED",
                 "scope": "All identity cells agree as canonical ASCII GSMs in every source row; sample rows are unique. Repeated labels remain a schema finding; their origin and scientific validity are not established."}
                for name, count in Counter(header[index] for index in positions).items() if count > 1]
    return header, by_id, {"columns": positions, "rows": len(records), "blank_header_binding": binding,
                          "schema_findings": repeated,
                          "scope": "Full-file identity checks only; every source header and cell is retained without renaming, dropping columns or rewriting the prediction source. Not a scientific repair."}


def _prepared(record, header, predictions):
    case, control = record.get("case_ids"), record.get("control_ids")
    if not isinstance(case, list) or not isinstance(control, list) or not case or not control:
        raise ValueError("Declared arms are absent or empty")
    members = case + control
    if any(not isinstance(sid, str) or not re.fullmatch(r"GSM[0-9]+", sid, flags=re.ASCII) for sid in members) or len(set(members)) != len(members):
        raise ValueError("Declared sample membership overlaps, duplicates or lacks canonical GSM identity")
    clock = record.get("clock")
    if not isinstance(clock, str) or header.count(clock) != 1:
        raise ValueError("Requested source clock is absent or ambiguous")
    missing = sorted(set(members) - set(predictions))
    if missing:
        raise ValueError("Missing source-defined predictions: " + ",".join(missing))
    column = header.index(clock)
    values = [float(predictions[sid][column]) for sid in members]
    if not all(math.isfinite(value) for value in values):
        raise ValueError("Requested source predictions contain nonfinite values")
    rows = [{"sample_id": sid, "arm": "case" if index < len(case) else "control", "value": value, "unit_id": ""}
            for index, (sid, value) in enumerate(zip(members, values))]
    return case, control, rows


def _design(record, case, control, rows):
    candidate = record["calculation"]
    if not isinstance(candidate, dict):
        raise ValueError("Calculation is not an object")
    for name in ("clock", "case_ids", "control_ids"):
        if name in candidate and candidate[name] != record[name]:
            raise ValueError("Calculation contradicts its declared source binding: " + name)
    if "n_case_assays" in candidate or "n_control_assays" in candidate:
        if any(candidate.get(key) is not None for key in ("p_value", "effect_size", "confidence_interval", "statistic", "degrees_of_freedom")):
            raise ValueError("Descriptive assay output contains inferential values")
        return "descriptive", 0.95, []
    method = candidate.get("statistical_test")
    methods = {"Welch t-test": "welch", "Student t-test": "student", "Paired t-test": "paired", "OLS adjusted contrast": "ols", "OLS HC3 adjusted contrast": "ols_hc3"}
    if method not in methods:
        raise NotImplementedError("Declared calculation method has no independent R adapter")
    mode = methods[method]
    confidence = _number(candidate.get("confidence_level"))
    if not 0 < confidence < 1 or candidate.get("alternative") != "two-sided":
        raise ValueError("Confidence level or alternative hypothesis is absent or unsupported")
    if mode == "paired":
        pairs = candidate.get("pairs")
        if not isinstance(pairs, list) or len(pairs) != len(case) or len(case) != len(control):
            raise ValueError("Explicit one-to-one pairing is required")
        if (any(not isinstance(pair, dict) or not isinstance(pair.get("unit_id"), str) or not pair["unit_id"].strip() or "\r" in pair["unit_id"] or "\n" in pair["unit_id"] for pair in pairs)
                or len({pair["unit_id"] for pair in pairs}) != len(pairs)
                or {pair.get("case_id") for pair in pairs} != set(case) or {pair.get("control_id") for pair in pairs} != set(control)):
            raise ValueError("Pair units or exact GSM membership are invalid")
        units = {pair[arm]: pair["unit_id"] for pair in pairs for arm in ("case_id", "control_id")}
        for row in rows:
            row["unit_id"] = units[row["sample_id"]]
    names = []
    if mode.startswith("ols"):
        estimator = candidate.get("variance_estimator")
        if estimator == "HC3 with residual-df t reference":
            mode = "ols_hc3"
        elif estimator != "classical homoskedastic OLS" or mode == "ols_hc3":
            raise NotImplementedError("OLS variance estimator/reference distribution is not declared exactly")
        covariates = candidate.get("covariates")
        if not isinstance(covariates, list) or not covariates:
            raise ValueError("OLS requires explicit sample-bound covariate values")
        for index, covariate in enumerate(covariates):
            if not isinstance(covariate, dict):
                raise ValueError("OLS covariate is not an object")
            name, values = covariate.get("name"), covariate.get("values")
            if not isinstance(name, str) or not name or name in ["intercept", "case", *names] or not isinstance(values, dict) or set(values) != set(case + control):
                raise ValueError("OLS covariate names or sample bindings are incomplete")
            names.append(name)
            for row in rows:
                row[f"cov_{index}"] = _number(values[row["sample_id"]])
    return mode, confidence, names


def _compare(record, mode, covariate_names, rows, reference, relative, absolute):
    candidate = record["calculation"]
    checks = []
    def check(passed, name, observed=None, expected=None):
        checks.append({"name": name, "status": "PASS" if passed else "FAIL", "observed": observed, "expected": expected})
    def scalar(field):
        values = reference.get(field, [])
        if len(values) != 1:
            raise ValueError("Independent R result lacks one scalar: " + field)
        return values[0][1]
    def close(a, b):
        try:
            return math.isclose(_number(a), _number(b), rel_tol=relative, abs_tol=absolute)
        except (ValueError, TypeError):
            return False
    def numeric(key, ref=None):
        expected = scalar(ref or key)
        observed = candidate.get(key)
        if key in {"p_value", "independent_p_value"}:
            # An absolute tolerance can silently accept large relative errors
            # in small P values. Zero agrees only with independent underflow.
            try:
                passed = 0 <= _number(observed) <= 1 and math.isclose(observed, expected, rel_tol=relative, abs_tol=0)
            except (ValueError, TypeError):
                passed = False
        else:
            passed = close(observed, expected)
        check(passed, key, observed, expected)
    for key in ("case_mean", "control_mean", "mean_difference"):
        numeric(key)
    if mode == "descriptive":
        numeric("n_case_assays", "n_case"); numeric("n_control_assays", "n_control")
        return checks
    for key in ("n_case", "n_control", "p_value", "standard_error", "statistic", "degrees_of_freedom"):
        numeric(key)
    if "independent_p_value" in candidate:
        numeric("independent_p_value", "p_value")
    ci = candidate.get("confidence_interval")
    expected_ci = [scalar("ci_low"), scalar("ci_high")]
    check(isinstance(ci, list) and len(ci) == 2 and all(close(a, b) for a, b in zip(ci, expected_ci)), "confidence_interval", ci, expected_ci)
    expected_estimand = "adjusted case coefficient" if mode.startswith("ols") else "case-minus-control mean"
    check(candidate.get("confidence_interval_estimand") == expected_estimand, "confidence_interval_estimand", candidate.get("confidence_interval_estimand"), expected_estimand)
    definition = OLS_DEFINITION if mode.startswith("ols") else (PAIRED_DEFINITION if mode == "paired" else POOLED_DEFINITION)
    check(candidate.get("effect_size_definition") == definition, "effect_size_definition", candidate.get("effect_size_definition"), definition)
    if mode.startswith("ols"):
        check(candidate.get("effect_size") is None, "OLS standardized effect remains undefined", candidate.get("effect_size"), None)
        numeric("adjusted_mean_difference"); numeric("residual_degrees_of_freedom")
        for key in ("coefficients", "residuals", "leverage"):
            expected = [value for _, value in reference.get(key, [])]
            observed = candidate.get(key)
            check(isinstance(observed, list) and len(observed) == len(expected) and all(close(a, b) for a, b in zip(observed, expected)), key, observed, expected)
        expected_columns = ["intercept", "case", *covariate_names]
        check(candidate.get("design_columns") == expected_columns, "design_columns", candidate.get("design_columns"), expected_columns)
        matrix = [[1., float(row["arm"] == "case"), *[row[f"cov_{index}"] for index in range(len(covariate_names))]] for row in rows]
        observed = candidate.get("design_matrix")
        check(isinstance(observed, list) and len(observed) == len(matrix) and all(isinstance(a, list) and len(a) == len(b) and all(close(x, y) for x, y in zip(a, b)) for a, b in zip(observed, matrix)), "design_matrix binds declared sample covariates", observed, matrix)
    else:
        numeric("effect_size")
        if mode == "paired":
            numeric("n_biological_units"); numeric("pooled_standardized_difference")
            check(candidate.get("pooled_standardized_difference_definition") == POOLED_DEFINITION, "pooled effect definition", candidate.get("pooled_standardized_difference_definition"), POOLED_DEFINITION)
            by_unit = dict(reference["within_unit_difference"])
            expected = [by_unit[pair["unit_id"]] for pair in candidate["pairs"]]
            observed = candidate.get("within_unit_differences")
            check(isinstance(observed, list) and len(observed) == len(expected) and all(close(a, b) for a, b in zip(observed, expected)), "within_unit_differences", observed, expected)
        else:
            numeric("case_variance"); numeric("control_variance")
    return checks


def model_verification_view(receipt):
    """Return a compact, allowlisted numerical view for the separate reviewer.

    Raw reasons, R diagnostics, filesystem paths, commands, pair-unit labels and
    arbitrary source text stay in the local receipt. Only bounded identifiers,
    status/check counts, hashes and finite numerical reference metrics leave it.
    Unusual identifiers are represented by their hash and input index, without
    silently truncating or rewriting the source label.
    """
    if not isinstance(receipt, dict) or receipt.get("schema") != SCHEMA:
        raise ValueError("Expected an independent-statistics receipt")
    statuses = ("PASS", "FAIL", "CANNOT VERIFY")
    methods = {"descriptive", "welch", "student", "paired", "ols", "ols_hc3"}
    scalar_names = {"case_mean", "control_mean", "mean_difference", "n_case", "n_control", "p_value",
                    "standard_error", "statistic", "degrees_of_freedom", "effect_size", "case_variance",
                    "control_variance", "n_biological_units", "pooled_standardized_difference",
                    "adjusted_mean_difference", "residual_degrees_of_freedom"}
    check_names = scalar_names | {"n_case_assays", "n_control_assays", "independent_p_value", "confidence_interval",
        "confidence_interval_estimand", "effect_size_definition", "OLS standardized effect remains undefined",
        "coefficients", "residuals", "leverage", "design_columns", "design_matrix binds declared sample covariates",
        "pooled effect definition", "within_unit_differences"}

    def status(value):
        return value if value in statuses else "CANNOT VERIFY"

    def digest(value):
        return value if isinstance(value, str) and re.fullmatch(r"[a-f0-9]{64}", value) else None

    def finite(value):
        try:
            return _number(value)
        except (ValueError, TypeError):
            return None

    def bounded_label(value):
        if isinstance(value, str) and len(value) <= 160 and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9 _.,:()+%='\[\]-]*", value):
            return value
        return None

    results = []
    for row in receipt.get("results", []):
        item = {name: status(row.get(name)) for name in ("status", "numerical_status", "stop_verification")}
        item["scientific_status"] = "CANNOT VERIFY"
        item["input_index"] = row.get("input_index") if type(row.get("input_index")) is int and row["input_index"] >= 0 else None
        for name in ("comparison_id", "clock"):
            value = row.get(name)
            item[name] = bounded_label(value)
            if item[name] is None:
                item[name + "_sha256"] = hashlib.sha256(_json_bytes(value)).hexdigest()
                item[name + "_label_status"] = "WITHHELD_UNSAFE_OR_UNSUPPORTED_LABEL"
        item["method"] = row.get("method") if row.get("method") in methods else None
        checks = row.get("checks", [])
        item["check_counts"] = {key: sum(status(check.get("status")) == key for check in checks) for key in statuses}
        item["failed_checks"] = [check["name"] for check in checks if check.get("status") == "FAIL" and check.get("name") in check_names]
        item["unlisted_check_count"] = sum(check.get("name") not in check_names for check in checks)
        reference = row.get("reference_values", {})
        metrics = {}
        for name in sorted(scalar_names):
            values = reference.get(name)
            if isinstance(values, list) and len(values) == 1 and isinstance(values[0], list) and len(values[0]) == 2:
                value = finite(values[0][1])
                if value is not None:
                    metrics[name] = value
        interval = [reference.get(key) for key in ("ci_low", "ci_high")]
        if all(isinstance(values, list) and len(values) == 1 and isinstance(values[0], list) and len(values[0]) == 2
               and finite(values[0][1]) is not None for values in interval):
            metrics["confidence_interval"] = [finite(values[0][1]) for values in interval]
        if item["method"] == "descriptive":
            for arm in ("case", "control"):
                if "n_" + arm in metrics:
                    metrics["n_" + arm + "_assays"] = metrics.pop("n_" + arm)
        item["reference_metrics"] = metrics
        results.append(item)
    return {"schema": "clockbase-independent-statistics-model-view-v1",
            "status": status(receipt.get("status")), "scientific_status": "CANNOT VERIFY",
            "scope": "Independent arithmetic under declared source groups/design; numerical agreement is not scientific acceptance",
            "request_sha256": digest(receipt.get("request_sha256")),
            "predictions_sha256": digest(receipt.get("predictions_sha256_after")),
            "source_hash_status": status(receipt.get("source_hash_status")),
            "counts": {key: sum(row["status"] == key for row in results) for key in statuses},
            "results": results}


def verify_calculations(calculations, predictions_path, output_dir, *, rscript="/usr/bin/Rscript",
                        timeout_seconds=60, relative_tolerance=1e-8, absolute_tolerance=1e-10):
    """Verify calculations independently; write immutable receipts and return one.

    Exact GSM/clock identity is checked against the prediction CSV. Blank first
    headers require an explicitly recorded full-column GSM-format check. Repeated
    allowlisted identity labels require exact full-file agreement and retain a
    schema finding with unresolved provenance; no columns are removed. Output
    status is numerical only. A source stop may be independently confirmed while
    numerical status remains CANNOT VERIFY. Existing completed outputs resume
    only when all request, source, implementation and output hashes agree.
    """
    started = time.monotonic()
    if not isinstance(calculations, list) or any(not isinstance(value, dict) for value in calculations):
        raise ValueError("Calculations must be a list of records")
    if not 0 < timeout_seconds <= 600 or not math.isfinite(relative_tolerance) or relative_tolerance < 0 or not math.isfinite(absolute_tolerance) or absolute_tolerance < 0:
        raise ValueError("Invalid bounded timeout or numerical tolerances")
    source = Path(predictions_path).resolve()
    executable = shutil.which(str(rscript))
    executable = Path(executable).resolve() if executable else None
    source_hash = _sha(source) if source.is_file() else None
    request = {"schema": SCHEMA, "calculations": calculations, "predictions": {"path": str(source), "sha256": source_hash},
               "rscript": {"path": str(executable) if executable else str(rscript), "sha256": _sha(executable) if executable else None},
               "implementation_sha256": _sha(Path(__file__)), "r_reference_sha256": hashlib.sha256(R_REFERENCE.encode()).hexdigest(),
               "timeout_seconds": timeout_seconds, "relative_tolerance": relative_tolerance, "absolute_tolerance": absolute_tolerance,
               "p_value_comparison": "relative tolerance only; zero requires agreement on numerical underflow"}
    request_bytes = _json_bytes(request)
    request_sha = hashlib.sha256(request_bytes).hexdigest()
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    receipt_path = output / "receipt.json"
    if receipt_path.exists():
        checksum_path = output / "receipt.sha256"
        if not checksum_path.is_file() or _sha(receipt_path) != checksum_path.read_text().strip():
            raise ValueError("Independent verification receipt changed on resume")
        prior = json.loads(receipt_path.read_text())
        if prior["request_sha256"] != request_sha:
            raise ValueError("Independent verification inputs or implementation changed; use a new output directory")
        for relative, expected in prior["outputs_sha256"].items():
            path = (output / relative).resolve()
            path.relative_to(output)
            if _sha(path) != expected:
                raise ValueError("Independent verification output changed on resume: " + relative)
        return prior
    if any(output.iterdir()):
        raise FileExistsError("Incomplete independent-verification output is preserved; use a new directory")
    with (output / "request.json").open("xb") as handle:
        handle.write(request_bytes)
    (output / "independent-reference.R").write_text(R_REFERENCE)
    source_error, binding = None, None
    try:
        if source_hash is None:
            raise FileNotFoundError("Cannot proceed - required data file not found: " + str(source))
        header, predictions, binding = _csv(source)
        _write_json(output / "source-identity.json", binding)
    except (ValueError, OSError) as error:
        header, predictions = [], {}
        source_error = str(error)
    results, jobs, pending = [], [], {}
    for index, record in enumerate(calculations):
        identifier = f"comparison-{index:06d}"
        result = {"comparison_id": record.get("comparison_id"), "clock": record.get("clock"), "input_index": index,
                  "status": "CANNOT VERIFY", "numerical_status": "CANNOT VERIFY", "stop_verification": "CANNOT VERIFY",
                  "scientific_status": "CANNOT VERIFY", "checks": [], "scope": "conditional arithmetic under declared groups/design only"}
        results.append(result)
        try:
            if source_error:
                raise ValueError(source_error)
            case, control, rows = _prepared(record, header, predictions)
        except (ValueError, TypeError, KeyError) as error:
            result["reason"] = str(error)
            result["stop_verification"] = "PASS" if record.get("status") == "CANNOT VERIFY" and not record.get("calculation") else "CANNOT VERIFY"
            continue
        if record.get("status") != "PASS":
            result["reason"] = "Producer retained an unresolved or failed calculation; scientific/design reasons are not adjudicated by a numerical verifier."
            if record.get("calculation") is not None or record.get("status") not in {"FAIL", "CANNOT VERIFY"}:
                result.update(status="FAIL", numerical_status="FAIL", reason="Calculation or status contradicts the unresolved producer record")
            continue
        try:
            mode, confidence, names = _design(record, case, control, rows)
        except NotImplementedError as error:
            result["reason"] = str(error)
            continue
        except (ValueError, TypeError, KeyError) as error:
            result.update(status="FAIL", numerical_status="FAIL", reason=str(error))
            continue
        data_file = output / f"{identifier}-input.csv"
        with data_file.open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
        jobs.append({"job_id": identifier, "data_file": str(data_file), "mode": mode, "confidence": confidence})
        pending[identifier] = (index, mode, names, rows)
    process_record = {"executed": False, "returncode": None, "timed_out": False}
    if jobs and executable is not None:
        jobs_file = output / "jobs.csv"
        with jobs_file.open("x", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(jobs[0]), lineterminator="\n")
            writer.writeheader(); writer.writerows(jobs)
        command = [str(executable), "--vanilla", str(output / "independent-reference.R"), str(jobs_file), str(output)]
        environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
        process_record.update(executed=True, command=command)
        try:
            completed = subprocess.run(command, capture_output=True, text=True, check=False, timeout=timeout_seconds, env=environment)
            process_record.update(returncode=completed.returncode, stdout=completed.stdout, stderr=completed.stderr)
        except subprocess.TimeoutExpired as error:
            def text(value):
                return value.decode(errors="replace") if isinstance(value, bytes) else (value or "")
            process_record.update(timed_out=True, stdout=text(error.stdout), stderr=text(error.stderr))
        except OSError as error:
            process_record["error"] = str(error)
    elif jobs:
        process_record["error"] = "Rscript is unavailable; no independent calculation executed"
    outcomes = {}
    if process_record["returncode"] == 0 and not process_record["timed_out"]:
        try:
            with (output / "r-outcomes.csv").open(newline="") as handle:
                records = list(csv.DictReader(handle))
            outcomes = {row["job_id"]: row for row in records}
            if len(outcomes) != len(records) or set(outcomes) != set(pending):
                raise ValueError("R outcomes omitted, duplicated or invented a job")
        except (OSError, ValueError, KeyError) as error:
            outcomes = {}
            process_record["result_error"] = str(error)
    _write_json(output / "r-process.json", process_record)
    for identifier, (index, mode, names, rows) in pending.items():
        result = results[index]
        result["method"] = mode
        outcome = outcomes.get(identifier)
        if outcome is None:
            result["reason"] = "Independent R execution unavailable, timed out or failed; no numerical verdict"
            continue
        if outcome["status"] != "PASS":
            result["reason"] = "Independent R stopped: " + outcome.get("message", "unspecified error")
            continue
        try:
            with (output / f"{identifier}-reference.csv").open(newline="") as handle:
                values = list(csv.DictReader(handle))
            reference = {}
            seen = set()
            for value in values:
                key = (value["field"], value["index"])
                if key in seen:
                    raise ValueError("Duplicated independent result field/index")
                seen.add(key)
                number = float(value["value"])
                if not math.isfinite(number):
                    raise ValueError("Nonfinite independent result")
                reference.setdefault(value["field"], []).append([value["index"], number])
            checks = _compare(calculations[index], mode, names, rows, reference, relative_tolerance, absolute_tolerance)
            passed = all(check["status"] == "PASS" for check in checks)
            result.update(status="PASS" if passed else "FAIL", numerical_status="PASS" if passed else "FAIL", checks=checks,
                          reference_values=reference, reason="Independent source-value arithmetic agrees" if passed else "Independent calculation or definition disagrees")
        except (OSError, ValueError, TypeError, KeyError) as error:
            result["reason"] = "Independent result unavailable or malformed: " + str(error)
    final_source_hash = _sha(source) if source.is_file() else None
    stable = final_source_hash == source_hash
    if not stable:
        for result in results:
            result.update(status="CANNOT VERIFY", numerical_status="CANNOT VERIFY", reason="Prediction source changed during independent verification")
    counts = dict(Counter(result["status"] for result in results))
    status = "FAIL" if not stable or counts.get("FAIL") else ("CANNOT VERIFY" if not results or counts.get("CANNOT VERIFY") else "PASS")
    outputs = {str(path.relative_to(output)): _sha(path) for path in sorted(output.rglob("*")) if path.is_file()}
    receipt = {"schema": SCHEMA, "status": status, "scientific_status": "CANNOT VERIFY", "accepted_scientific_repairs": 0,
               "scope": "Independent base-R numerical verification; biological/source design and scientific acceptance are separate",
               "recorded_at_utc": datetime.now(timezone.utc).isoformat(), "elapsed_seconds": time.monotonic()-started,
               "request_sha256": request_sha, "predictions_sha256_before": source_hash, "predictions_sha256_after": final_source_hash,
               "source_identity": binding,
               "source_hash_status": "PASS" if stable else "FAIL", "counts": {key: counts.get(key, 0) for key in ("PASS", "FAIL", "CANNOT VERIFY")},
               "verified_input_stops": sum(row["stop_verification"] == "PASS" for row in results),
               "results": results, "outputs_sha256": outputs, "r_process": process_record,
               "model_calls": 0, "api_cost_usd": 0, "clock_runs": 0,
               "limitations": ["Groups, pair identities and covariate values are supplied declarations; this function does not independently reconstruct them from source metadata.",
                               "PASS confirms calculations and effect/interval definitions under the stated design, not independence, distributional assumptions, measurement applicability, causality or biological aging.",
                               "No multiplicity adjustment, scientific error rate, historical input-lineage certification, accepted repair or original-file mutation is produced."]}
    _write_json(receipt_path, receipt)
    with (output / "receipt.sha256").open("x") as handle:
        handle.write(_sha(receipt_path) + "\n")
    return receipt
