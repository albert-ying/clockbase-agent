"""Deterministic, read-only measurement evidence for completed pipeline runs."""

from __future__ import annotations

import hashlib
import json
import csv
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import pandas as pd


SCHEMA = "clockbase-pipeline-measurement-evidence-v1"
STAGES = ("02_molecular_input", "04_preprocessing", "05_clock_predictions")
REPO_ROOT = Path(__file__).resolve().parents[1]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _checked(path: str | Path, expected: str, what: str) -> Path:
    resolved = Path(path).resolve()
    if not resolved.is_file():
        raise ValueError(f"{what} is not a file: {resolved}")
    observed = _sha256(resolved)
    if observed != expected:
        raise ValueError(f"{what} SHA256 mismatch: expected {expected}, observed {observed}")
    return resolved


def _json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _stage(root: Path, name: str) -> tuple[dict[str, Any], dict[Path, str], dict[str, str]]:
    directory = (root / name).resolve()
    receipt = directory / "complete.json"
    sidecar = directory / "complete.sha256"
    if not sidecar.is_file():
        raise ValueError(f"missing completed-stage sidecar: {sidecar}")
    expected = sidecar.read_text(encoding="ascii").strip()
    _checked(receipt, expected, f"{name} completion receipt")
    document = _json(receipt)
    outputs = document.get("outputs")
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError(f"{name} has no declared output members")
    members: dict[Path, str] = {}
    for relative, digest in outputs.items():
        if not isinstance(relative, str) or not isinstance(digest, str):
            raise ValueError(f"invalid {name} output declaration")
        member = (directory / relative).resolve()
        if directory not in member.parents:
            raise ValueError(f"{name} output escapes stage directory: {relative}")
        _checked(member, digest, f"{name} output {relative}")
        members[member] = digest
    if not isinstance(document.get("result"), dict):
        raise ValueError(f"{name} result is not an object")
    return document["result"], members, {"path": str(receipt), "sha256": expected}


def _samples(source_index: Mapping[str, Any], study_id: str) -> list[str]:
    if source_index.get("study_id", study_id) != study_id:
        raise ValueError("source_index study_id mismatch")
    values = source_index.get("sample_ids")
    if not isinstance(values, list) or not values or not all(isinstance(x, str) and x for x in values):
        raise ValueError("source_index.sample_ids must be a nonempty string list")
    if len(values) != len(set(values)):
        raise ValueError("source_index.sample_ids contains duplicates")
    return values


def _record_file(files: dict[str, dict[str, str]], role: str, path: Path, digest: str) -> None:
    _checked(path, digest, role)
    declaration = {"path": str(path), "sha256": digest}
    previous = files.get(role)
    if previous is not None and previous != declaration:
        raise ValueError(f"conflicting source file role: {role}")
    files[role] = declaration


def _member(
    files: dict[str, dict[str, str]], role: str, path: str | Path, digest: str,
    all_members: Mapping[Path, str],
) -> Path:
    resolved = Path(path).resolve()
    if all_members.get(resolved) != digest:
        raise ValueError(f"{role} is not an exact completed-stage member")
    _record_file(files, role, resolved, digest)
    return resolved


def _frame(path: Path, *, beta: bool) -> pd.DataFrame:
    # The production strict reader performs logical GEO parsing, width and ID
    # validation, domain checks, and before/after hashing. Importing lazily
    # keeps this readback helper out of pipeline module import cycles.
    from src.pipeline_inputs_2026_09_10 import read_matrix

    return read_matrix(path, beta=beta)


def _table(path: Path, *, index_expected: list[str]) -> pd.DataFrame:
    before = _sha256(path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        try:
            header = next(csv.reader(handle))
        except StopIteration as exc:
            raise ValueError(f"empty table: {path}") from exc
    if len(header) < 2 or len(header[1:]) != len(set(header[1:])):
        raise ValueError(f"duplicate or absent literal data columns: {path}")
    frame = pd.read_csv(path, index_col=0, float_precision="round_trip")
    after = _sha256(path)
    if before != after:
        raise ValueError(f"file changed while being read: {path}")
    if frame.index.astype(str).tolist() != index_expected:
        raise ValueError(f"sample identity/order mismatch: {path}")
    if frame.index.has_duplicates or frame.columns.has_duplicates:
        raise ValueError(f"duplicate identity in table: {path}")
    numeric = frame.apply(pd.to_numeric, errors="raise")
    array = numeric.to_numpy(dtype=float)
    if np.isinf(array).any():
        raise ValueError(f"infinite value in table: {path}")
    return numeric


def _compact_sources(
    files: dict[str, dict[str, str]], sources: Mapping[str, Any],
    all_members: Mapping[Path, str], prefix: str, run_root: Path,
    sealed_sources: Mapping[str, Any],
    trusted_external: Mapping[Path, tuple[str, dict[str, str] | None]],
) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for name, declaration in sorted(sources.items()):
        if not isinstance(declaration, Mapping):
            raise ValueError(f"invalid measurement source {name}")
        path, digest = declaration.get("path"), declaration.get("sha256")
        if not isinstance(path, str) or not isinstance(digest, str):
            raise ValueError(f"unbound measurement source {name}")
        resolved = Path(path).resolve()
        if sealed_sources.get(name) != declaration:
            raise ValueError(f"measurement source is absent from sealed asset manifest: {name}")
        # No file beneath the run may be promoted to evidence unless the
        # completed-stage output ledger binds it exactly.
        if run_root == resolved or run_root in resolved.parents:
            _member(files, f"{prefix}:{name}", resolved, digest, all_members)
        else:
            trusted = trusted_external.get(resolved)
            if trusted is None or trusted[0] != digest:
                raise ValueError(f"external measurement source lacks a sealed binding: {name}")
            if trusted[1] is not None:
                _record_file(files, f"asset_manifest:{resolved.parent.name}",
                             Path(trusted[1]["path"]), trusted[1]["sha256"])
            _record_file(files, f"{prefix}:{name}", resolved, digest)
        hashes[name] = digest
    return hashes


def _trusted_external_sources(
    run_root: Path, source_bindings: Mapping[str, Any],
) -> dict[Path, tuple[str, dict[str, str] | None]]:
    trusted: dict[Path, tuple[str, dict[str, str] | None]] = {}
    for declaration in source_bindings.values():
        if isinstance(declaration, Mapping):
            path, digest = declaration.get("path"), declaration.get("sha256")
            if isinstance(path, str) and isinstance(digest, str):
                trusted[Path(path).resolve()] = (digest, None)
    asset_root = REPO_ROOT / "pipeline_assets"
    implementation: Mapping[str, Any] = {}
    version_path = run_root.parent / "version.json"
    if version_path.is_file():
        version = _json(version_path)
        value = version.get("implementation", {})
        if not isinstance(value, Mapping):
            raise ValueError("version implementation binding is not an object")
        implementation = value
    for manifest_path in sorted(asset_root.glob("*/manifest.json")):
        manifest_digest = _sha256(manifest_path)
        relative_manifest = str(manifest_path.relative_to(REPO_ROOT))
        if implementation.get(relative_manifest) != manifest_digest:
            continue
        manifest = _json(manifest_path)
        entries = manifest.get("files")
        if not isinstance(entries, Mapping):
            raise ValueError(f"asset manifest lacks files mapping: {manifest_path}")
        declaration = {"path": str(manifest_path.resolve()), "sha256": manifest_digest}
        for relative, entry in entries.items():
            if not isinstance(relative, str) or not isinstance(entry, Mapping):
                raise ValueError(f"invalid asset manifest entry: {manifest_path}")
            digest = entry.get("sha256")
            if not isinstance(digest, str):
                raise ValueError(f"asset manifest entry lacks SHA256: {relative}")
            asset_path = (manifest_path.parent / relative).resolve()
            if manifest_path.parent.resolve() not in asset_path.parents:
                raise ValueError(f"asset manifest path escapes directory: {relative}")
            trusted[asset_path] = (digest, declaration)
    return trusted


def _recovered_dnam_source(
    run_root: Path, study_id: str, molecular: Mapping[str, Any],
    source_bindings: Mapping[str, Any], source_sample_ids: list[str],
    members: dict[Path, str], source_files: dict[str, dict[str, str]],
    stage_receipts: dict[str, dict[str, str]],
) -> dict[str, Any]:
    """Bind the original empty input, supplement and independently checked recovery."""
    from src.pipeline_inputs_2026_09_10 import propose_dnam_recovery
    from audit_system.report_2026_09_10 import verify_dnam_recovery

    documents = {}
    stage02 = (run_root / "02_molecular_input").resolve()
    for filename in ("molecular-input.json", "recovery-proposal.json",
                     "reconstruction.json", "recovery-verification.json"):
        paths = [path for path in members if path.name == filename and stage02 in path.parents]
        if len(paths) != 1:
            raise ValueError("Recovered DNAm input lacks one sealed " + filename)
        path = paths[0]
        _member(source_files, "stage02:recovery:" + filename, path, members[path], members)
        documents[filename] = _json(path)
    if documents["molecular-input.json"] != molecular:
        raise ValueError("Recovered molecular result differs from its sealed manifest")
    proposal = documents["recovery-proposal.json"]
    reconstruction = documents["reconstruction.json"]
    if (proposal.get("study_id") != study_id or proposal.get("original_feature_count") != 0
            or proposal.get("proposed_action") != "reconstruct_beta_input_from_linked_supplement"
            or molecular.get("repair_status") != "PASS"):
        raise ValueError("Recovered DNAm proposal does not establish this empty-input repair")
    if any(molecular.get(key) != value for key, value in reconstruction.items()
           if key != "repair_status"):
        raise ValueError("Recovered matrix differs from its sealed reconstruction")
    metadata, metadata_members, metadata_receipt = _stage(run_root, "01_metadata")
    members.update(metadata_members)
    stage_receipts["01_metadata"] = metadata_receipt
    sources = proposal.get("source_files", {})
    if not {"failed_matrix", "supplement", "metadata", "series"} <= set(sources):
        raise ValueError("Recovered DNAm proposal lacks source roles")
    if set(sources) - {"failed_matrix", "supplement", "metadata", "series", "sample_records"}:
        raise ValueError("Recovered DNAm proposal contains undeclared source roles")

    def bound(role, declaration, path=None):
        candidate = sources[role]
        if (not isinstance(declaration, Mapping)
                or candidate.get("sha256") != declaration.get("sha256")
                or Path(str(candidate.get("path"))).resolve()
                != Path(str(path or declaration.get("path"))).resolve()):
            raise ValueError("Recovered DNAm source role differs from sealed binding: " + role)
        _record_file(source_files, "stage02:recovery_source:" + role,
                     Path(candidate["path"]), candidate["sha256"])

    bound("failed_matrix", source_bindings.get("matrix"))
    bound("supplement", source_bindings.get("dnam_supplement"))
    if (Path(str(molecular.get("repaired_original"))).resolve()
            != Path(sources["failed_matrix"]["path"]).resolve()
            or Path(str(molecular.get("source"))).resolve()
            != Path(sources["supplement"]["path"]).resolve()
            or molecular.get("source_sha256") != sources["supplement"]["sha256"]):
        raise ValueError("Recovered molecular input loses its original or supplementary source")
    for role, field, version_role in (("metadata", "metadata", None),
                                     ("series", "series_source", "series"),
                                     ("sample_records", "source_sample_records", "samples")):
        if role not in sources:
            continue
        path = Path(str(metadata.get(field))).resolve()
        digest = metadata_members.get(path)
        if digest is None:
            raise ValueError("Recovered DNAm metadata is not a sealed metadata-stage member")
        bound(role, {"path": str(path), "sha256": digest})
        if version_role and source_bindings.get(version_role, {}).get("sha256") != digest:
            raise ValueError("Recovered DNAm primary metadata differs from version source")
    recovered_ids = [row["gsm_id"] for row in proposal.get("crosswalk", [])]
    if len(recovered_ids) != len(set(recovered_ids)) or set(recovered_ids) != set(source_sample_ids):
        raise ValueError("Recovered DNAm crosswalk differs from source sample membership")
    repeated_proposal = propose_dnam_recovery(
        study_id, sources["failed_matrix"]["path"], sources["supplement"]["path"],
        sources["metadata"]["path"], sources["series"]["path"],
        sample_source=sources.get("sample_records", {}).get("path"))
    for key in ("study_id", "source_url", "crosswalk", "source_layout", "original_feature_count"):
        if repeated_proposal[key] != proposal.get(key):
            raise ValueError("Recovered DNAm source proposal no longer reproduces: " + key)
    verified = verify_dnam_recovery(proposal, reconstruction)
    if verified != documents["recovery-verification.json"] or verified.get("status") != "PASS":
        raise ValueError("Recovered DNAm independent verification does not reproduce")
    return {"original_feature_count": 0,
            "original_matrix_sha256": sources["failed_matrix"]["sha256"],
            "supplement_sha256": sources["supplement"]["sha256"],
            "reconstructed_matrix_sha256": molecular["sha256"],
            "source_column_to_gsm": [{"source_column": r["matrix_sample_id"],
                                      "gsm_id": r["gsm_id"], "mapping_rule": r["mapping_rule"]}
                                     for r in proposal["crosswalk"]],
            "independent_input_verification": verified,
            "scope": "Source recovery and every-value verification; no preprocessing, clocks or biological validity inferred"}


def build_measurement_evidence(
    root: str | Path,
    study_id: str,
    modality: str,
    source_index: Mapping[str, Any],
    source_bindings: Mapping[str, Any],
) -> dict[str, Any]:
    """Build a compact measurement bundle without executing preprocessing/clocks."""
    if modality not in {"human_rna", "mouse_rna", "human_dnam"}:
        raise ValueError(f"unsupported measurement modality: {modality}")
    run_root = Path(root).resolve()
    source_sample_ids = _samples(source_index, study_id)
    source_files: dict[str, dict[str, str]] = {}
    trusted_external = _trusted_external_sources(run_root, source_bindings)

    # Recheck all sealed version sources. These bind molecular input and the
    # source index, even when some files are not exposed as M evidence.
    for role, declaration in sorted(source_bindings.items()):
        if not isinstance(declaration, Mapping):
            raise ValueError(f"invalid source binding: {role}")
        path, digest = declaration.get("path"), declaration.get("sha256")
        if not isinstance(path, str) or not isinstance(digest, str):
            raise ValueError(f"unbound source binding: {role}")
        _record_file(source_files, f"version:{role}", Path(path).resolve(), digest)

    stage_results: dict[str, dict[str, Any]] = {}
    members: dict[Path, str] = {}
    stage_receipts: dict[str, dict[str, str]] = {}
    for stage in STAGES:
        result, stage_members, receipt = _stage(run_root, stage)
        stage_results[stage] = result
        members.update(stage_members)
        stage_receipts[stage] = receipt

    # Stages 04/05 must bind their result objects through explicit JSON output
    # members, rather than relying only on fields embedded in complete.json.
    sealed_manifests: dict[str, dict[str, Any]] = {}
    for stage, filename in (("04_preprocessing", "preprocessing.json"),
                            ("05_clock_predictions", "predictions.json")):
        candidates = [path for path in members if path.name == filename and stage in path.parts]
        if len(candidates) != 1:
            raise ValueError(f"{stage} lacks one sealed {filename} manifest")
        manifest = _json(candidates[0])
        if manifest != stage_results[stage]:
            raise ValueError(f"{stage} result differs from sealed {filename}")
        sealed_manifests[stage] = manifest
        _member(source_files, f"{stage}:manifest", candidates[0], members[candidates[0]], members)

    molecular = stage_results["02_molecular_input"]
    raw_path = Path(str(molecular.get("path"))).resolve()
    raw_digest = molecular.get("sha256")
    if not isinstance(raw_digest, str):
        raise ValueError("molecular input lacks SHA256")
    _member(source_files, "stage02:matrix", raw_path, raw_digest, members)
    matrix_binding = source_bindings.get("matrix")
    archs4_binding = source_bindings.get("archs4")
    recovered_dnam = None
    if molecular.get("repaired_original") is not None:
        if modality != "human_dnam":
            raise ValueError("Supplement recovery is only supported for the declared DNAm route")
        recovered_dnam = _recovered_dnam_source(
            run_root, study_id, molecular, source_bindings, source_sample_ids,
            members, source_files, stage_receipts)
    elif isinstance(matrix_binding, Mapping):
        origin_digest = molecular.get("source_sha256", raw_digest)
        if origin_digest != matrix_binding.get("sha256"):
            raise ValueError("molecular input does not bind source_bindings.matrix")
    elif isinstance(archs4_binding, Mapping) and modality in {"human_rna", "mouse_rna"}:
        extraction_candidates = [path for path in members if path.name == "extraction.json"]
        if len(extraction_candidates) != 1 or _json(extraction_candidates[0]) != molecular:
            raise ValueError("ARCHS4 extraction lacks a matching sealed extraction receipt")
        if (Path(str(molecular.get("source_hdf5"))).resolve()
                != Path(str(archs4_binding.get("path"))).resolve()):
            raise ValueError("ARCHS4 extraction source path differs from sealed binding")
        memberships = molecular.get("all_study_memberships")
        if (not isinstance(memberships, Mapping)
                or set(memberships) != set(source_sample_ids)
                or any(not isinstance(values, list) or study_id not in values
                       for values in memberships.values())):
            raise ValueError("ARCHS4 extraction membership differs from source index/study")
        if molecular.get("source_size") != Path(str(archs4_binding.get("path"))).stat().st_size:
            raise ValueError("ARCHS4 extraction source size differs from sealed HDF5")
        _member(source_files, "stage02:extraction_manifest", extraction_candidates[0],
                members[extraction_candidates[0]], members)
    else:
        raise ValueError("molecular input lacks a supported matrix or ARCHS4 source binding")

    raw = _frame(raw_path, beta=modality == "human_dnam")
    sample_ids = raw.columns.astype(str).tolist()
    if set(sample_ids) != set(source_sample_ids):
        raise ValueError("raw molecular matrix sample identity membership mismatch")
    raw_missing = raw.isna().sum(axis=0).astype(int).to_dict()
    if molecular.get("feature_count", raw.shape[0]) != raw.shape[0] or molecular.get(
            "sample_count", raw.shape[1]) != raw.shape[1]:
        raise ValueError("molecular extraction dimensions differ from matrix")
    derived_matrix = None
    if "matrix_provenance" in source_bindings:
        # This function is called only after the source-only design is frozen.
        # Reuse the source graph validator; never expose these receipts during
        # source design or rehash its historical HDF5 origin in every child.
        from src.update_pipeline_2026_09_20 import read_matrix_provenance

        declared = source_bindings["matrix_provenance"]
        derived = read_matrix_provenance(declared["path"], declared["sha256"], source_bindings, study_id, modality)
        expected_roles = {"matrix_provenance" if role == "wrapper" else "matrix_provenance:" + role
                          for role in derived["files"]}
        declared_roles = {role for role in source_bindings if role == "matrix_provenance" or role.startswith("matrix_provenance:")}
        if expected_roles != declared_roles:
            raise ValueError("derived matrix nested sources differ from sealed version inventory")
        for role, source in derived["files"].items():
            key = "matrix_provenance" if role == "wrapper" else "matrix_provenance:" + role
            sealed = source_bindings[key]
            if Path(sealed["path"]).resolve() != Path(source["path"]) or sealed["sha256"] != source["sha256"]:
                raise ValueError("derived matrix file lacks exact sealed version binding: " + role)
        derived_matrix = derived["measurement"]
        if (derived_matrix["feature_count"] != raw.shape[0] or derived_matrix["sample_count"] != raw.shape[1]
                or derived_matrix["sample_ids"] != sample_ids):
            raise ValueError("derived matrix receipt dimensions/order differ from completed input")

    idat_derived = None
    idat_record = None
    idat_roles = {role for role in source_bindings
                  if role == "dnam_idat_provenance" or role.startswith("dnam_idat_provenance:")}
    if idat_roles:
        if (modality != "human_dnam" or "dnam_idat_provenance" not in source_bindings
                or recovered_dnam is not None):
            raise ValueError("DNAm IDAT provenance requires a separately derived DNAm matrix")
        from src.update_pipeline_2026_09_20 import read_dnam_idat_provenance

        declared = source_bindings["dnam_idat_provenance"]
        derived = read_dnam_idat_provenance(declared["path"], declared["sha256"], study_id, sources=source_bindings)
        expected_roles = {"dnam_idat_provenance" if role == "wrapper" else "dnam_idat_provenance:" + role
                          for role in derived["files"]}
        if idat_roles != expected_roles:
            raise ValueError("DNAm IDAT nested sources differ from sealed version inventory")
        for role, source in derived["files"].items():
            key = "dnam_idat_provenance" if role == "wrapper" else "dnam_idat_provenance:" + role
            sealed = source_bindings[key]
            if Path(sealed["path"]).resolve() != Path(source["path"]) or sealed["sha256"] != source["sha256"]:
                raise ValueError("DNAm IDAT file lacks exact sealed version binding: " + role)
        idat_record = {"path": str(Path(declared["path"]).resolve()), "sha256": declared["sha256"]}
        if molecular.get("idat_provenance") != idat_record or derived["sample_ids"] != sample_ids:
            raise ValueError("DNAm IDAT provenance does not bind stage 02 sample order and receipt")
        original_beta = _frame(Path(source_bindings["matrix"]["path"]), beta=True)
        if (not raw.index.equals(original_beta.index) or not raw.columns.equals(original_beta.columns)
                or not np.array_equal(raw.to_numpy(dtype=float), original_beta.to_numpy(dtype=float), equal_nan=True)):
            raise ValueError("DNAm IDAT stage 02 values or identities differ from the bound matrix")
        idat_derived = {**derived["measurement"], "scientific_status": derived["scientific_status"]}
    elif molecular.get("idat_provenance") is not None:
        raise ValueError("DNAm IDAT stage 02 receipt lacks a sealed provenance role")

    prep = stage_results["04_preprocessing"]
    if idat_record is not None:
        if prep.get("idat_provenance") != idat_record or prep.get("idat_processing") is not True:
            raise ValueError("DNAm IDAT preprocessing lost or changed the upstream provenance")
    elif prep.get("idat_provenance") is not None or prep.get("idat_processing") is True:
        raise ValueError("IDAT processing claim lacks a sealed upstream provenance")
    prediction = stage_results["05_clock_predictions"]
    evidence: list[dict[str, Any]] = []
    evidence.append({
        "id": "M00000",
        "kind": "molecular_matrix",
        "value": {
            "feature_count": int(raw.shape[0]),
            "sample_count": int(raw.shape[1]),
            "missing_values_total": int(raw.isna().sum().sum()),
            "missing_values_by_sample": {sample: int(raw_missing[sample]) for sample in sample_ids},
            "numeric_domain": "nonnegative finite counts" if modality in {"human_rna", "mouse_rna"} else "finite beta values in [0,1], with missing values allowed",
            "source_index_sample_order_sha256": hashlib.sha256(
                json.dumps(source_sample_ids, ensure_ascii=False, separators=(",", ":")).encode()
            ).hexdigest(),
            "raw_sample_order_sha256": hashlib.sha256(
                json.dumps(sample_ids, ensure_ascii=False, separators=(",", ":")).encode()
            ).hexdigest(),
            "raw_order_matches_source_index": sample_ids == source_sample_ids,
        },
        "sample_ids": sample_ids,
        "source": {"role": "stage02:matrix", "sha256": raw_digest},
    })

    if modality in {"human_rna", "mouse_rna"}:
        transformed = Path(str(prep.get("scaled_diff"))).resolve()
        transformed_digest = members.get(transformed)
        if transformed_digest is None:
            raise ValueError("RNA transformed matrix is not a completed-stage member")
        _member(source_files, "stage04:transformed", transformed, transformed_digest, members)
        transformed_frame = _table(transformed, index_expected=sample_ids)
        if prep.get("output_sha256", {}).get("scaled_diff.csv") != transformed_digest:
            raise ValueError("RNA transformed matrix differs from preprocessing hash field")
        prep_measurement_sources = prep.get("measurement_sources", {})
        if prep_measurement_sources.get("raw_counts") != {"path": str(raw_path), "sha256": raw_digest}:
            raise ValueError("RNA preprocessing raw input does not bind stage 02")
        prep_sources = _compact_sources(
            source_files, prep_measurement_sources, members, "stage04:source", run_root,
            sealed_manifests["04_preprocessing"].get("measurement_sources", {}),
            trusted_external,
        )
        evidence.append({
            "id": "M00001", "kind": "preprocessing", "sample_ids": sample_ids,
            "value": {
                "profile": prep.get("normalization_profile"),
                "log_transform": prep.get("log_transform"),
                "filter_order": prep.get("filter_order"),
                "filter_count_minimum": prep.get("filter_count_minimum"),
                "filter_sample_fraction": prep.get("filter_sample_fraction"),
                "mapped_filtered_genes": prep.get("mapped_filtered_genes"),
                "transformed_feature_count": int(transformed_frame.shape[1]),
                "transformed_missing_by_sample": {
                    sample: int(transformed_frame.loc[sample].isna().sum()) for sample in sample_ids
                },
                "historical_producing_version_verified": prep.get("historical_producing_version_verified"),
                "source_stratum": prep.get("source_stratum"),
            },
            "source": {"role": "stage04:transformed", "sha256": transformed_digest,
                       "input_hashes": prep_sources},
        })
    else:
        imputed = Path(str(prep.get("path"))).resolve()
        imputed_digest = members.get(imputed)
        if imputed_digest is None:
            raise ValueError("DNAm imputed matrix is not a completed-stage member")
        _member(source_files, "stage04:imputed_beta", imputed, imputed_digest, members)
        imputed_frame = _frame(imputed, beta=True)
        if imputed_frame.columns.astype(str).tolist() != sample_ids:
            raise ValueError("imputed beta matrix sample identity/order mismatch")
        reference = prep.get("reference")
        if not isinstance(reference, Mapping):
            raise ValueError("DNAm preprocessing lacks reference declaration")
        ref_source = _checked(reference.get("source", ""), reference.get("source_sha256", ""), "DNAm reference source")
        _record_file(source_files, "stage04:reference_source", ref_source, reference["source_sha256"])
        if prep.get("raw") != str(raw_path):
            raise ValueError("DNAm preprocessing raw input does not bind stage 02")
        dnam_binding = source_bindings.get("dnam_reference")
        if (not isinstance(dnam_binding, Mapping)
                or Path(str(dnam_binding.get("path"))).resolve() != ref_source
                or dnam_binding.get("sha256") != reference.get("source_sha256")):
            raise ValueError("DNAm reference does not bind source_bindings.dnam_reference")
        evidence.append({
            "id": "M00001", "kind": "preprocessing", "sample_ids": sample_ids,
            "value": {
                "profile": "reference-mean imputation",
                "imputed_feature_count": int(imputed_frame.shape[0]),
                "imputed_missing_values_total": int(imputed_frame.isna().sum().sum()),
                "raw_methylation_preprocessing_recomputed": prep.get("raw_methylation_preprocessing_recomputed"),
                "independent_imputation_verification": prep.get("independent_imputation_verification"),
                "idat_processing": prep.get("idat_processing"),
                "reference_study": reference.get("reference_study"),
                "reference_reader_contract": reference.get("reader_contract"),
                "reference_means_sha256": reference.get("means_sha256"),
            },
            "source": {"role": "stage04:imputed_beta", "sha256": imputed_digest,
                       "reference_source_sha256": reference.get("source_sha256")},
        })

    predictions_path = Path(str(prediction.get("path"))).resolve()
    predictions_digest = members.get(predictions_path)
    if predictions_digest is None:
        raise ValueError("prediction table is not a completed-stage member")
    _member(source_files, "stage05:predictions", predictions_path, predictions_digest, members)
    predictions = _table(predictions_path, index_expected=sample_ids)
    clocks = prediction.get("clocks")
    if (not isinstance(clocks, list) or len(clocks) != len(set(clocks))
            or predictions.columns.astype(str).tolist() != clocks):
        raise ValueError("clock inventory does not match prediction columns")
    if prediction.get("prediction_sha256", predictions_digest) != predictions_digest:
        raise ValueError("prediction table differs from prediction hash field")
    evidence.append({
        "id": "M00002", "kind": "prediction_inventory", "sample_ids": sample_ids,
        "value": {
            "clocks": clocks, "sample_count": len(sample_ids),
            "clock_predictions_recomputed": prediction.get("clock_predictions_recomputed"),
            "scientific_status": prediction.get("scientific_status"),
            "measurement_units": "CANNOT VERIFY",
            "calibration_applicability": "CANNOT VERIFY",
        },
        "source": {"role": "stage05:predictions", "sha256": predictions_digest},
    })

    checks = prediction.get("checks")
    if not isinstance(checks, list) or len(checks) != len(clocks):
        raise ValueError("prediction check inventory does not match clocks")
    expected_check_clocks = [clock[:-5] if clock.endswith("_tAge") else clock for clock in clocks]
    check_clocks = [check.get("clock") if isinstance(check, Mapping) else None for check in checks]
    if check_clocks != expected_check_clocks or len(check_clocks) != len(set(check_clocks)):
        raise ValueError("prediction checks do not correspond uniquely to prediction clocks")
    for number, check in enumerate(checks, start=3):
        if not isinstance(check, Mapping):
            raise ValueError("invalid prediction check")
        measurement_sources = check.get("measurement_sources", {})
        if not isinstance(measurement_sources, Mapping):
            raise ValueError(f"invalid measurement source ledger for {check.get('clock')}")
        raw_role = "raw_beta" if modality == "human_dnam" else "raw_counts"
        if measurement_sources.get(raw_role) != {"path": str(raw_path), "sha256": raw_digest}:
            raise ValueError(f"{check.get('clock')} raw measurement source does not bind stage 02")
        source_hashes = _compact_sources(
            source_files, measurement_sources, members,
            f"stage05:{check.get('clock')}:source", run_root,
            sealed_manifests["05_clock_predictions"]["checks"][number - 3].get("measurement_sources", {}),
            trusted_external,
        )
        required = int(check.get("required_features"))
        clock = str(check.get("clock"))
        if modality in {"human_rna", "mouse_rna"}:
            coverage_path = next((p for p in members if p.name == f"{clock}-feature-coverage.json"), None)
            if coverage_path is None:
                raise ValueError(f"missing feature coverage member for {clock}")
            coverage_digest = members[coverage_path]
            if measurement_sources.get("feature_coverage") != {
                    "path": str(coverage_path), "sha256": coverage_digest}:
                raise ValueError(f"feature coverage source does not bind stage 05 for {clock}")
            if measurement_sources.get("transformed_input") != {
                    "path": str(transformed), "sha256": transformed_digest}:
                raise ValueError(f"transformed input does not bind stage 04 for {clock}")
            if check.get("model_sha256") != measurement_sources.get("model", {}).get("sha256"):
                raise ValueError(f"model hash field differs from measurement source for {clock}")
            _member(source_files, f"stage05:{clock}:coverage", coverage_path, coverage_digest, members)
            coverage = _json(coverage_path)
            feature_ids = coverage.get("feature_ids")
            retained_ids = coverage.get("training_retained_ids")
            dropped_ids = coverage.get("training_dropped_ids")
            selected_ids = coverage.get("selected_ids")
            nonzero_ids = coverage.get("nonzero_coefficient_ids")
            missing_map = coverage.get("missing_ids_by_sample")
            imputed_map = coverage.get("imputed_ids_by_sample")
            if (not isinstance(feature_ids, list) or len(feature_ids) != required
                    or len(feature_ids) != len(set(feature_ids))
                    or not all(isinstance(values, list) for values in
                               (retained_ids, dropped_ids, selected_ids, nonzero_ids))
                    or not isinstance(missing_map, Mapping) or not isinstance(imputed_map, Mapping)):
                raise ValueError(f"invalid RNA feature coverage for {clock}")
            retained_set, dropped_set = set(retained_ids), set(dropped_ids)
            if (retained_set & dropped_set or retained_set | dropped_set != set(feature_ids)
                    or retained_ids != [feature for feature in feature_ids if feature in retained_set]
                    or dropped_ids != [feature for feature in feature_ids if feature in dropped_set]
                    or len(retained_ids) != len(retained_set) or len(dropped_ids) != len(dropped_set)
                    or len(selected_ids) != len(set(selected_ids)) or not set(selected_ids) <= retained_set
                    or len(nonzero_ids) != len(set(nonzero_ids)) or not set(nonzero_ids) <= set(selected_ids)):
                raise ValueError(f"RNA training/selection feature partition is invalid for {clock}")
            selected = transformed_frame.reindex(columns=feature_ids)
            derived_missing = {
                sample: [feature for feature in feature_ids if not pd.notna(selected.at[sample, feature])]
                for sample in sample_ids
            }
            if any(list(missing_map[sample]) != derived_missing[sample] for sample in sample_ids):
                raise ValueError(f"RNA missing-feature ledger differs from transformed matrix for {clock}")
            observed = {sample: required - len(derived_missing[sample]) for sample in sample_ids}
            missing = {sample: len(derived_missing[sample]) for sample in sample_ids}
            derived_imputed = {
                sample: [feature for feature in derived_missing[sample] if feature in retained_set]
                for sample in sample_ids
            }
            if any(list(imputed_map[sample]) != derived_imputed[sample] for sample in sample_ids):
                raise ValueError(f"RNA imputed-feature ledger differs from missing retained features for {clock}")
            imputed = {sample: len(derived_imputed[sample]) for sample in sample_ids}
            missing_dropped = {
                sample: len([feature for feature in derived_missing[sample] if feature in dropped_set])
                for sample in sample_ids
            }
            if any(observed[s] + imputed[s] + missing_dropped[s] != required for s in sample_ids):
                raise ValueError(f"RNA feature coverage partition does not sum to required count for {clock}")
            declared_selected = coverage.get("observed_selected_by_sample")
            declared_nonzero = coverage.get("observed_nonzero_by_sample")
            observed_selected = transformed_frame.reindex(columns=selected_ids).notna().sum(axis=1).astype(int).to_dict()
            observed_nonzero = transformed_frame.reindex(columns=nonzero_ids).notna().sum(axis=1).astype(int).to_dict()
            if (not isinstance(declared_selected, Mapping) or not isinstance(declared_nonzero, Mapping)
                    or any(int(declared_selected[sample]) != observed_selected[sample] for sample in sample_ids)
                    or any(int(declared_nonzero[sample]) != observed_nonzero[sample] for sample in sample_ids)):
                raise ValueError(f"RNA selected/nonzero observed summaries differ from recomputation for {clock}")
            strategy = coverage.get("imputation_strategy")
            source_hashes["coverage"] = coverage_digest
        else:
            coefficient_path = next(
                (p for p in members if p.name == f"{clock}-coefficients.csv"), None
            )
            if coefficient_path is None:
                raise ValueError(f"missing coefficients for {clock}")
            coefficient_digest = members[coefficient_path]
            if measurement_sources.get("coefficients") != {
                    "path": str(coefficient_path), "sha256": coefficient_digest}:
                raise ValueError(f"coefficient source does not bind stage 05 for {clock}")
            _member(source_files, f"stage05:{clock}:coefficients", coefficient_path, coefficient_digest, members)
            coefficient_frame = pd.read_csv(coefficient_path, dtype={"CpGmarker": str})
            markers = [x for x in coefficient_frame["CpGmarker"].tolist() if x.lower() != "intercept"]
            if len(markers) != required or len(markers) != len(set(markers)):
                raise ValueError(f"coefficient feature inventory mismatch for {clock}")
            selected = raw.reindex(markers)
            observed = {sample: int(selected[sample].notna().sum()) for sample in sample_ids}
            imputed = {sample: required - observed[sample] for sample in sample_ids}
            missing = dict(imputed)
            missing_dropped = {sample: 0 for sample in sample_ids}
            retained_ids, dropped_ids, selected_ids = markers, [], markers
            coefficients = dict(zip(coefficient_frame["CpGmarker"].astype(str),
                                    pd.to_numeric(coefficient_frame["CoefficientTraining"], errors="raise")))
            nonzero_ids = [marker for marker in markers if coefficients[marker] != 0]
            observed_selected = observed
            observed_nonzero = {
                sample: int(raw.reindex(nonzero_ids)[sample].notna().sum()) for sample in sample_ids
            }
            strategy = "reference-mean imputation"
            source_hashes["coefficients"] = coefficient_digest
        if min(observed.values()) != int(check.get("observed_features_min")):
            raise ValueError(f"observed feature minimum mismatch for {clock}")
        evidence.append({
            "id": f"M{number:05d}", "kind": "clock_feature_coverage", "sample_ids": sample_ids,
            "value": {
                "clock": clock, "required_features": required,
                "observed_features_by_sample": observed,
                "missing_features_by_sample": missing,
                "imputed_features_by_sample": imputed,
                "missing_training_dropped_features_by_sample": missing_dropped,
                "training_retained_feature_count": len(retained_ids),
                "training_dropped_feature_count": len(dropped_ids),
                "selected_feature_count": len(selected_ids),
                "nonzero_coefficient_feature_count": len(nonzero_ids),
                "observed_selected_features_by_sample": {
                    sample: int(observed_selected[sample]) for sample in sample_ids
                },
                "observed_nonzero_features_by_sample": {
                    sample: int(observed_nonzero[sample]) for sample in sample_ids
                },
                "imputation_strategy": strategy,
                "numerical_status": check.get("numerical_status"),
                "max_abs_independent_difference": check.get("max_abs_independent_difference"),
                "measurement_applicability": check.get("measurement_applicability"),
                "measurement_units": "CANNOT VERIFY",
                "calibration_applicability": "CANNOT VERIFY",
            },
            "source": {"input_hashes": source_hashes},
        })

    if derived_matrix is not None:
        evidence.append({
            "id": f"M{len(evidence):05d}", "kind": "historical_archs4_matrix_derivation", "sample_ids": sample_ids,
            "value": derived_matrix,
            "source": {"role": "version:matrix_provenance", "sha256": source_bindings["matrix_provenance"]["sha256"]},
        })

    if recovered_dnam is not None:
        evidence.append({
            "id": f"M{len(evidence):05d}", "kind": "source_linked_dnam_input_recovery",
            "sample_ids": sample_ids, "value": recovered_dnam,
            "source": {"role": "stage02:recovery:recovery-verification.json",
                       "sha256": source_files["stage02:recovery:recovery-verification.json"]["sha256"]},
        })

    if idat_derived is not None:
        evidence.append({
            "id": f"M{len(evidence):05d}", "kind": "source_bound_dnam_idat_preprocessing",
            "sample_ids": sample_ids, "value": idat_derived,
            "source": {"role": "version:dnam_idat_provenance", "sha256": idat_record["sha256"]},
        })

    # Detect mutation across the full readback interval.
    for role, declaration in source_files.items():
        _checked(declaration["path"], declaration["sha256"], f"final recheck {role}")
    for stage, declaration in stage_receipts.items():
        _checked(declaration["path"], declaration["sha256"], f"final recheck {stage}")

    evidence_map = {
        item["id"]: {key: value for key, value in item.items() if key != "id"}
        for item in evidence
    }
    try:
        encoded_evidence = json.dumps(evidence_map, sort_keys=True, allow_nan=False)
    except ValueError as exc:
        raise ValueError(f"non-finite value in compact measurement evidence: {exc}") from exc
    if str(run_root) in encoded_evidence or "/Entropy/" in encoded_evidence:
        raise ValueError("absolute path leaked into compact measurement evidence")
    return {
        "schema": SCHEMA,
        "study_id": study_id,
        "modality": modality,
        "evidence": evidence_map,
        "source_files": source_files,
        "stage_receipts": stage_receipts,
        "scientific_status": "CANNOT VERIFY",
        "limitations": [
            "Measurement evidence does not establish biological sample units, grouping, pairing, or covariate control.",
            "Clock units, calibration, biological applicability, and scientific acceptance remain CANNOT VERIFY unless separately sourced.",
        ],
    }
